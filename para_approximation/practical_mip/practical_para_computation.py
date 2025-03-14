import os
import time
import json
import numpy as np
from pyomo.opt import TerminationCondition
from definitions import *
from para_approximation.practical_mip.para_inexact_modeler import ParaboloidModel
from para_approximation.exact_mip.test_para_computation import ExactParaboloidHandler
from para_relaxation.pyomo_solver import PyomoSolver
from para_approximation.utilities.result_checker import check_function_violation, check_approximation_violation
from para_approximation.utilities.zigzag_utilities import check_zigzag_violation
from para_approximation.utilities.test_utilities import *


class InexactParaboloidHandler(ExactParaboloidHandler):
    def __init__(self, function_string, function_setting_path, model_setting_path, n_paraboloids, approx_below):
        # initialize the same parameters as in the exact setting
        super().__init__(function_string, function_setting_path, model_setting_path, n_paraboloids, approx_below)

        ## additional initializations
        # initialize function and model settings
        assert os.path.exists(function_setting_path), f"path {function_setting_path} does not exist"
        assert os.path.exists(model_setting_path), f"path {model_setting_path} does not exist"
        function_setting = json.load(open(function_setting_path))
        model_setting = json.load(open(model_setting_path))

        # initialize function and model specific parameters
        self.func_lb = float(function_setting["lb"])
        self.func_ub = float(function_setting["ub"])
        self.func_int_length = self.func_ub - self.func_lb
        self.func_maxL = float(function_setting["L-constant"])
        self.model_eps = float(model_setting["eps"])

        # custom values for the discretization in the inexact MIP approach and increase parameters
        self.n_delta_t_increase = np.ceil(self.func_int_length / self.model_eps * self.func_maxL)
        self.n_delta_d_increase = np.ceil(self.func_int_length / self.model_eps * self.func_maxL)
        self.n_max_delta_t = np.ceil(self.func_int_length / self.model_eps * self.func_maxL / 10)
        self.n_max_delta_d = np.ceil(self.func_int_length / self.model_eps * self.func_maxL / 10)

        # values to set for respecting solutions and to increase the number of paraboloids when solutions are unsatisfying
        self.objective_tolerance = 10.0
        self.n_function_violations = 0
        self.n_approximation_violations = 0

    def solution_step(self, quads, lins, cons, max_time_limit, print_output=False):
        # initialize the start time for modelling and booleans to control the search loop
        end_of_loop = False
        next_loop = False
        start_modelling_time = time.time()

        # initialize paraboloid finding model
        self.para_modeler = ParaboloidModel(self.function_setting_path, self.model_setting_path, self.n_paraboloids,
                                            self.approx_below, self.n_max_delta_t, self.n_max_delta_d)
        # compute the number of binaries and continuous variables to potentially exit for too large problems
        terminate_modelling, end_of_loop = self._check_model_size(print_output)
        if terminate_modelling:
            return end_of_loop, not end_of_loop, None

        results = self._model_and_solve(quads, lins, cons, start_modelling_time, max_time_limit, print_output)
        # if results is None, a memory error was encountered, thus the current no of paraboloids is an upper bound
        if results is None:
            # update the upper bound to be the current number of paraboloids
            self.n_paraboloids_ub = self.n_paraboloids
            # make an update step as solution counts as not valid
            end_of_loop = self.update_search_bounds(valid_solution=False, update_lower_bound=False)
            # update the next loop variable
            next_loop = not end_of_loop

        return end_of_loop, next_loop, results

    def validation_step(self, results, print_output=False):
        # track time for validation
        start_validation_time = time.time()

        # initialize an end of loop parameter
        end_of_loop = False
        # create parameter for validity of solution
        valid_solution = False

        # initialize paraboloid coefficients as empty lists
        quads, lins, cons = [], [], []

        # dependent on the solver status (objective too large, infeasible problem, time limit w/o feasible solution),
        # increase the number of paraboloids immediately:
        time_limit_reached = results.solver.termination_condition == TerminationCondition.maxTimeLimit
        infeasible_problem = results.solver.termination_condition == TerminationCondition.infeasible
        # cannot explicitly catch memory limit, but the error message suffices
        memory_limit = "Solver quit with a problem" in str(results.solver.message)

        # increase the amount of paraboloids if they do not satisfy the requirements
        if (results.problem.lower_bound > self.objective_tolerance or infeasible_problem or
                ((time_limit_reached or memory_limit) and
                 (np.isnan(results.problem.upper_bound) or np.isinf(results.problem.upper_bound)))):
            if infeasible_problem:
                if print_output:
                    print(f"Increase of number of paraboloids due to infeasibility!")
            elif results.problem.lower_bound > self.objective_tolerance:
                if print_output:
                    print(f"Increase of number of paraboloids due to objective tolerance!")
            else:
                if print_output:
                    print(f"Increase of number of paraboloids due to time limit!")

            # reset violation/approximation counter
            self.n_function_violations = 0
            self.n_approximation_violations = 0

            return end_of_loop, valid_solution, quads, lins, cons

        # extract the parameters
        quads, lins, cons = self.para_modeler.extract_results(self.model)

        # extract function string and domain
        domain = [[self.para_modeler.lb], [self.para_modeler.ub]]

        # check if a paraboloid exceeds the function
        if self.function_string == "zigzag":
            function_violation, f_p_distances = check_zigzag_violation(quads, lins, cons, domain, "function",
                                                                       eps=self.para_modeler.eps,
                                                                       dim=self.para_modeler.dim,
                                                                       print_result=print_output,
                                                                       approx_below=self.approx_below)
        else:
            function_violation, f_p_distances = check_function_violation(self.function_string, quads, lins, cons, domain,
                                                                         dim=self.para_modeler.dim,
                                                                         print_result=print_output,
                                                                         approx_below=self.approx_below)

        # compute manipulated constant coefficients of the paraboloids which result from shifting the current ones
        # exactly below function f
        manip_cons = shift_cons_for_violation(function_violation, f_p_distances, cons)

        # sanity check:
        if self.function_string == "zigzag":
            shifted_function_violation, f_p_distances = check_zigzag_violation(quads, lins, manip_cons, domain, "function",
                                                                               eps=self.para_modeler.eps,
                                                                               dim=self.para_modeler.dim,
                                                                               print_result=False,
                                                                               approx_below=self.approx_below)
        else:
            shifted_function_violation, f_p_distances = check_function_violation(self.function_string, quads, lins, manip_cons,
                                                                                 domain, dim=self.para_modeler.dim,
                                                                                 print_result=False,
                                                                                 approx_below=self.approx_below)
        # -100_000 signals infeasibility and, thus, violation in any case
        if min(f_p_distances) > -99_999:
            assert not shifted_function_violation, f"function violation can not be false if shifted paras are considered"

        # if a paraboloid exceeds the function, increase the amount of discretization points
        if function_violation:
            # increase the amount of discretization points linearly on base value
            self.n_max_delta_d += self.n_delta_d_increase
            self.n_function_violations += 1
            # update end of loop, as the number of paraboloids is not to be increased
            end_of_loop = True
            if self.n_function_violations == 20:
                self.n_function_violations = 0
                self.n_paraboloids += 1
                #return end_of_loop, valid_solution, quads, lins, cons

        # check if the eps approximation is violated with the shifted constants
        if self.function_string == "zigzag":
            approximation_violation, approximation_violation_distance = check_zigzag_violation(
                quads, lins, manip_cons, domain, "approximation", eps=self.para_modeler.eps, dim=self.para_modeler.dim,
                print_result=print_output, approx_below=self.approx_below)
        else:
            approximation_violation, approximation_violation_distance = check_approximation_violation(
                self.function_string, quads, lins, manip_cons, domain, eps=self.para_modeler.eps, dim=self.para_modeler.dim,
                print_result=print_output, approx_below=self.approx_below)

        # if epsilon approximation is violated increase the amount of segments
        if approximation_violation:
            self.n_max_delta_t += self.n_delta_t_increase
            self.n_approximation_violations += 1
            # update end of loop, as the number of paraboloids is not to be increased
            end_of_loop = True
            if self.n_approximation_violations == 20:
                self.n_approximation_violations = 0
                self.n_paraboloids += 1

        # compute the validity of the solution
        valid_solution = not function_violation and not approximation_violation

        # update time for validation
        self.times[2] += time.time() - start_validation_time

        return end_of_loop, valid_solution, quads, lins, manip_cons

