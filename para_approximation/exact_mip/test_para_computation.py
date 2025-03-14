import os
import time
import numpy as np
import signal
from pyomo.opt import TerminationCondition
from definitions import *
from para_approximation.exact_mip.para_modeler import ParaboloidModel
from para_relaxation.pyomo_solver import PyomoSolver
from para_approximation.utilities.result_checker import check_function_violation, check_approximation_violation
from para_approximation.utilities.zigzag_utilities import check_zigzag_violation
from para_approximation.utilities.test_utilities import *


class ExactParaboloidHandler(object):
    def __init__(self, function_string, function_setting_path, model_setting_path, n_paraboloids, approx_below):
        # string of function to approximate
        self.function_string = function_string

        # initialize the paths for the function and model settings
        self.function_setting_path = function_setting_path
        self.model_setting_path = model_setting_path

        # current/starting number of paraboloids, lower and upper bound
        self.n_paraboloids = n_paraboloids
        self.n_paraboloids_lb = self.n_paraboloids
        self.n_paraboloids_ub = np.inf

        # approximation from below (True) or above (False)
        self.approx_below = approx_below

        # modeller and model
        self.para_modeler = None
        self.model = None

        # TODO: internal time tracking
        self.times = [0, 0, 0] #modelling, solving, evaluation

    def solution_step(self, quads, lins, cons, max_time_limit, print_output=False):
        # initialize the start time for modelling and booleans to control the search loop
        end_of_loop = False
        next_loop = False
        start_modelling_time = time.time()

        # print output for debugging
        if print_output:
            print(f"Starting Modelling")
        # initialize the modeler for the MIP model
        self.para_modeler = ParaboloidModel(self.function_setting_path, self.model_setting_path, self.n_paraboloids,
                                            self.approx_below)

        # assert for too large values right away by checking for empty d_indices
        if len(self.para_modeler.d_indices) == 0:
            self.n_paraboloids_ub = self.n_paraboloids
            end_of_loop = self.update_search_bounds(valid_solution=False, update_lower_bound=False)
            return end_of_loop, not end_of_loop, None

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

        # return end of loop indication and results
        return end_of_loop, next_loop, results

    def validation_step(self, results, print_output=False):
        # track time for validation
        start_validation_time = time.time()

        # initialize paraboloid coefficients as empty lists
        quads, lins, cons = [], [], []

        # dependent on the solver status, decrease/increase the number of paraboloids:
        # - optimal solution -> decrease
        # - feasible solution with positive objective -> manual test and decision based on this
        # - infeasible/time limit without solution -> increase

        # initialize solver status variables
        optimal_solu_found = results.solver.termination_condition == TerminationCondition.optimal
        time_limit_reached = results.solver.termination_condition == TerminationCondition.maxTimeLimit
        # cannot explicitly catch memory limit, but the error message suffices
        memory_limit = "Solver quit with a problem" in str(results.solver.message)

        # create parameter for validity of solution
        valid_solution = False
        # if the MIP was solved to optimality, store the validity and the result
        if optimal_solu_found:
            valid_solution = True
            quads, lins, cons = self.para_modeler.extract_results(self.model)
        # if the time limit was reached, check if a valid solution exists
        elif time_limit_reached or memory_limit:
            # an upper bound means that there is a feasible solution
            if not (np.isnan(results.problem.upper_bound) or np.isinf(results.problem.upper_bound)):
                # manually check for validity
                # extract the parameters of the paraboloids
                quads, lins, cons = self.para_modeler.extract_results(self.model)

                # extract necessary values to check for the violation
                domain = [[self.para_modeler.lb], [self.para_modeler.ub]]

                # check for a function violation and shift the constant paraboloid parameters accordingly
                if self.function_string == "zigzag":
                    function_violation, f_p_distances = check_zigzag_violation(quads, lins, cons, domain, "function",
                                                                               eps=self.para_modeler.eps,
                                                                               dim=self.para_modeler.dim,
                                                                               print_result=True,
                                                                               approx_below=self.approx_below)
                else:
                    function_violation, f_p_distances = check_function_violation(self.function_string, quads, lins, cons,
                                                                                 domain, dim=self.para_modeler.dim,
                                                                                 print_result=True,
                                                                                 approx_below=self.approx_below)

                shifted_cons = shift_cons_for_violation(function_violation, f_p_distances, cons)

                # check if the eps approximation is violated with the shifted constants
                if self.function_string == "zigzag":
                    approximation_violation, approximation_violation_distance = check_zigzag_violation(
                        quads, lins, shifted_cons, domain, "approximation", eps=self.para_modeler.eps, dim=1,
                        print_result=True, approx_below=self.approx_below)
                else:
                    approximation_violation, approximation_violation_distance = check_approximation_violation(
                        self.function_string, quads, lins, shifted_cons, domain, eps=self.para_modeler.eps, dim=1,
                        print_result=True, approx_below=self.approx_below)

                # set flag to valid solution if there is no approximation violation after shifting
                if not approximation_violation:
                    valid_solution = True
                    cons = shifted_cons

        # update time for validation
        self.times[2] += time.time() - start_validation_time

        # first parameter is end of loop, used in inexact model
        return False, valid_solution, quads, lins, cons

    def update_search_bounds(self, valid_solution, update_lower_bound=True):
        # initialize end of loop boolean
        end_of_loop = False
        # if a valid solution exists, extract the parameters and decrease the search number
        if valid_solution:
            # extract the parameters
            quads, lins, cons = self.para_modeler.extract_results(self.model)

            # set the current number of paraboloids as upper bound
            self.n_paraboloids_ub = self.n_paraboloids

            # exit the loop if lower and upper bound only differ by one
            if self.n_paraboloids_ub <= self.n_paraboloids_lb + 1:
                end_of_loop = True

            # choose the new number of paraboloids as the middle, but at least lb + 1]
            increase = int((self.n_paraboloids_ub - self.n_paraboloids_lb) / 2)
            self.n_paraboloids = max(increase + self.n_paraboloids_lb, self.n_paraboloids_lb + 1)

        # if there is no valid solution, increase the number of paraboloids
        else:
            # safe the current number as lower bound if not stated otherwise
            if update_lower_bound:
                self.n_paraboloids_lb = self.n_paraboloids

            # already exit the loop if lower and upper bound only differ by one
            if self.n_paraboloids_ub <= self.n_paraboloids_lb + 1:
                end_of_loop = True
            # if the upper bound is None, just double
            if np.isinf(self.n_paraboloids_ub):
                self.n_paraboloids *= 2
            else:
                # choose the new number of paraboloids as the middle, but at least lb + 1
                increase = int((self.n_paraboloids_ub - self.n_paraboloids_lb) / 2)
                self.n_paraboloids = max(increase + self.n_paraboloids_lb, self.n_paraboloids_lb + 1)

        return end_of_loop

    def _check_model_size(self, print_output):
        """ check for potential model size and signal termination if necessary """
        # initialize boolean for terminating modelling step
        terminate_modelling = False
        end_of_loop = False

        # initialize the number of relevant indices
        n_para_indices = len(self.para_modeler.para_indices)
        n_t_indices = len(self.para_modeler.t_indices)
        n_d_indices = len(self.para_modeler.d_indices)

        # compute the number of binaries and continuous variables to potentially exit for too large problems
        n_binaries = n_para_indices * n_t_indices
        n_continuous = 3 * n_para_indices + n_para_indices * n_d_indices
        if n_continuous >= n_max_continuous and n_binaries >= n_max_binaries:
            # update the termination parameter
            terminate_modelling = True
            # update the upper bound to be the current number of paraboloids
            self.n_paraboloids_ub = self.n_paraboloids
            # make an update step as solution counts as not valid
            end_of_loop = self.update_search_bounds(valid_solution=False, update_lower_bound=False)

        # print output if necessary
        if terminate_modelling & print_output:
            print(f"Terminating due to model size; #cont = {n_continuous} & #bin = {n_binaries}")

        # return end of loop and termination value, also number of variables
        return terminate_modelling, end_of_loop

    def _model_and_solve(self, quads, lins, cons, start_modelling_time, max_time_limit, print_output):
        """ models the necessary MIP and solves it """
        # set up the MIP model
        self.model = self.para_modeler.setup_model(False, [quads, lins, cons])

        # compute the modelling time
        start_solving_time = time.time()
        assert start_solving_time > start_modelling_time, f"solving must start after modelling"
        self.times[0] += start_solving_time - start_modelling_time

        # print output for debugging
        if print_output:
            print(f"Finished Modelling; Starting solving")

        # solve model
        solver = PyomoSolver("gurobi-direct", time_limit=max_time_limit)
        objective, results = solver.solve_model(self.model, tee=print_output)

        # give output for debugging
        if print_output:
            print(f"Finished solving")

        # compute solving time
        self.times[1] += time.time() - start_solving_time

        return results
