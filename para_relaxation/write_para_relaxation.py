import os

from definitions import *
from osil_parser.osil_parser import OSILParser
from osil_parser.osil_1Dreformulation import reformulate_osil_parser_to_1d, single_reformulation
from osil_parser.osil_expressions import *
from osil_to_pyomo import OsilPyomoConverter
from instance_extractor import attributes, get_instances_by_attribute

import numpy as np
import time
import copy
import json

# test set name
test_set = "minlplib"
# data format
data_suffix = "osil"

# epsilon for approximation guarantee
eps = 1e-2

# choice of attributes
selected_attributes = {
    "name" : ['batchdes'],#, 'contvar', 'eg_all_s', 'eg_disc2_s', 'eg_disc_s', 'eg_int_s', 'ex1222', 'ex14_1_3', 'ex14_1_4', 'ex3pb', 'ex8_1_1', 'ex8_1_2', 'ex8_2_1b', 'ex8_2_4b', 'ex8_4_6', 'ex8_4_7', 'feedtray', 'ghg_1veh', 'ghg_2veh', 'ghg_3veh', 'inscribedsquare01', 'inscribedsquare02', 'inscribedsquare03', 'kriging_peaks-red020', 'lnts100', 'lnts200', 'lnts400', 'lnts50', 'mathopt6', 'polygon100', 'polygon25', 'polygon50', 'polygon75', 'pooling_epa1', 'pooling_epa2', 'pooling_epa3', 'synthes2', 'synthes3', 't1000'],
    #"adddate" : ["2000-01-30", None],
    #"convex" : False,
    "dualbound" : [-np.inf, np.inf], # allows for deleting dualbound is None
    "ncons" : [0, 10000],
    "nvars" : [0, 10000],
    #"probtype" : ["INLP"]
}

solver_output = True
delete_orig = False

# REMINDER: ex14_1_9 needs really much accuracy

instances_all_subed = {}


def substitute_paraboloids(orig_parser, n_new_vars, delete_orig=True):
    re_parser = copy.deepcopy(orig_parser)

    instances_all_subed[re_parser.name] = True

    n_para_below = 0
    n_para_above = 0

    nl_indices_to_remove = []
    nl_constraints_to_readd = []

    # load the parabola parameters
    para_paramters_path = os.path.join(out_path, "para_parameters.json")
    json_file = open(para_paramters_path, "r")
    para_params = json.load(json_file)

    # define a parameter to track constraints which cannot be substituted
    n_cons_not_subable = 0

    # iterate over the nonlinear constraints
    for nl_index, nl in re_parser.nl_constraints.items():
        # substitute sine/cosine or exp functions
        if isinstance(nl, (OSILSine, OSILCosine, OSILExp)):
            # define the function string dependent on the osil expression type
            if isinstance(nl, OSILSine):
                function_string = "sin"
            elif isinstance(nl, OSILCosine):
                function_string = "cos"
            elif isinstance(nl, OSILExp):
                function_string = "exp"
            else:
                print("EXITING DUE TO ERROR")
                exit()

            # extract the non-linearity's coefficient, the variable index of its argument and the argument's bounds
            argument_coefficient = nl.coefficient
            argument_index = nl.expression
            assert isinstance(argument_index, (int,)), f"argument must be variable index in expression tree"
            argument_lb = re_parser.variables[argument_index].lb
            argument_ub = re_parser.variables[argument_index].ub
            argument_lb = -np.inf if argument_lb is None else argument_lb
            argument_ub = np.inf if argument_ub is None else argument_ub

            # # recompute the bounds due to the coefficient, i.e., (x'=) 5x for x in [-1, 1] -> x' in [-5, 5]
            temp_argument_lb = argument_lb * argument_coefficient
            temp_argument_ub = argument_ub * argument_coefficient
            search_argument_lb = min(temp_argument_lb, temp_argument_ub)
            search_argument_ub = max(temp_argument_lb, temp_argument_ub)

            # extract the parabola parameters dependent on the function string and the bounds
            no_valid_parabolas, para_params_until_approx = get_para_params(para_params, search_argument_lb,
                                                                           search_argument_ub, function_string)
            # continue if the arguments where too wide
            if no_valid_parabolas:
                if n_cons_not_subable == 0:
                    print(f"Found nonlinear constraint which can not be substituted")
                    instances_all_subed[re_parser.name] = False
                n_cons_not_subable += 1
                continue

            if search_argument_lb == search_argument_ub:
                continue

            # remove the nl from the parser and add variable with respective bounds instead

            # save current nl index for removal
            nl_indices_to_remove.append(nl_index)

            # per default, we want to add a variable y (= f(x)) and then bound y by the paraboloids from both sides,
            # i.e., y >= p_i and y <= p_i
            # for the case, when we already explicitly detect an inequality y >= f(x) or y <= f(x), we only use one side
            # that is the case when there is only one linear coefficient no quadratic one
            if len(re_parser.lin_coeffs[nl_index]) == 1 and len(re_parser.quad_coeffs[nl_index]) == 0:
                # initialize an indication of using an existing variable
                use_existing_var = True
                # get the variable index associated to y and its coefficient (the latter for using it below)
                var_index = re_parser.lin_coeffs[nl_index][0][0]
                var_coeff = re_parser.lin_coeffs[nl_index][0][1]

                # delete the (in)equality constraint, as we bound the variables by paraboloids below anyways
                if delete_orig:
                    re_parser.lin_coeffs[nl_index] = []
                else:
                    nl_indices_to_remove.pop()
            # else create a new variable and replace the occurrence of (co)sine/exp with it
            else:
                # initialize an indication of using an existing variable
                use_existing_var = False
                # create a new variable
                n_new_vars += 1
                new_variable_name = f"aux{n_new_vars}"
                new_variable = OSILVariable(new_variable_name, nl.lower_bound, nl.upper_bound)
                # save the new variable's index
                var_index = len(re_parser.variables)
                # add the variable instead of the non-linearity
                re_parser.variables.append(new_variable)
                new_lin_coef = (var_index, 1.0)
                re_parser.lin_coeffs[nl_index].append(new_lin_coef)

                # we can add paraboloids in addition to the original constraint
                if not delete_orig:
                    # to additionally create the constraint y = f(x) afterwards, save nonlinearity and variable index
                    nl_constraints_to_readd += [(nl_index, nl, var_index)]

            # if we use an existing variable, we only may need approximations from above or below
            if use_existing_var:
                approx_directions = []
                if re_parser.constraint_infos[nl_index][2] is not None:
                    approx_directions += ["below"]
                if re_parser.constraint_infos[nl_index][1] is not None:
                    approx_directions += ["above"]
            else:
                approx_directions = ["below", "above"]

            # for a new variable y, approximate it by the paraboloids, i.e., y >= p_i and y <= p'_j
            # for an existing variable z, we change f(x) + a * z <= rhs to p_i(x) + a * z <= rhs and
            # lhs <= f(x) + a * z to lhs <= p'_j(x) + a * z
            for approx_string in approx_directions:
                # extract suitable paraboloid parameters
                quads, lins, cons = para_params_until_approx[approx_string]

                # define the constant rhs dependent on whether to use an existing variable (extract bound) or not (0)
                if use_existing_var:
                    bound_index = 2 if approx_string == "below" else 1
                    constant = re_parser.constraint_infos[nl_index][bound_index]
                else:
                    constant = 0.0

                # add quadratic and linear coefficients for the arguments variable index
                # the constant parameter is added as rhs/lhs, -constant >= paraboloid - constant - y respective
                # -constant <= paraboloid - constant - y for upper bounds
                for quad, lin, con in zip(quads, lins, cons):
                    # number of current constraints gives new constraint index
                    n_constraints = len(re_parser.constraint_infos)
                    # depending on constraint type, approximate from below or above
                    if approx_string == "below": #bound_index == 1:
                        constraint_info = (f"para_below{n_para_below}", None, constant - con)
                        n_para_below += 1
                    else:
                        constraint_info = (f"para_above{n_para_above}", constant - con, None)
                        n_para_above += 1

                    # for the argument, say x, and the auxiliary variable, say y, include the linear coefficient
                    # and -1, respectively, as linear coefficients
                    if use_existing_var:
                        new_var_coef = (var_index, var_coeff)
                    else:
                        new_var_coef = (var_index, -1.0)
                    # manipulate linear adn quadratic coefficients by the argument's coeffcient, i.e., f(a * x) <= y
                    # iff f(z) <= y -> (relax) p(z) <= y iff p(a * x) <= y
                    lin_coefs = [(argument_index, lin * argument_coefficient), new_var_coef]
                    # only the argument needs its quadratic coefficient
                    quad_coefs = [(argument_index, argument_index, quad * argument_coefficient ** 2)]

                    # add the linear and quadratic coefficients to the parser as well as the constraint infos,
                    # creating a new constraint
                    re_parser.lin_coeffs[n_constraints] = lin_coefs
                    re_parser.quad_coeffs[n_constraints] = quad_coefs
                    re_parser.constraint_infos.append(constraint_info)

    # update the nonlinear constraints by deleting the original constraint
    for nl_index in nl_indices_to_remove:
        del re_parser.nl_constraints[nl_index]

    # if we keep the original constraints, re-add them as y=f(x)
    if not delete_orig:
        for (nl_index, nl, var_index) in nl_constraints_to_readd:
            # add new constraints on new index
            n_constraints = len(re_parser.constraint_infos)
            # initialize linear coefficient
            new_var_coef = (var_index, -1.0)
            # add to linear and nonlinear constraints
            re_parser.lin_coeffs[n_constraints] = new_var_coef
            re_parser.nl_constraints[n_constraints] = nl
            # add constraint infos, signalling the equality with rhs and lhs zero
            re_parser.constraint_infos.append((f"keeping nl_{nl_index}", 0, 0))

    return re_parser


def get_para_params(params, arg_lb, arg_ub, func_string):
    """ compute whether there exist valid parabolas and return a dict if so """
    # initialize the indication of an existing approximation and the approx direction dict
    invalid_approx = True
    params_until_approx = {}

    # set-up the parameter dictionary for the respective function string
    func_string_for_dict = "sin" if func_string == "cos" else func_string
    func_dict = params[str(eps)][func_string_for_dict]

    # make a case distinction by the function string
    if func_string == "exp":
        # test the upper bound to be below -2.0 which is the smallest lower bound for exp
        if arg_ub <= -2.0:
            # check that the lower bound is above -5.0
            if arg_lb >= -5.0:
                # extract the parameters belonging to the interval [-5, -2] and update the invalid approx
                invalid_approx = False
                params_until_approx = func_dict[str(-5.0)][str(-2.0)]
        # check whether the upper bound is at least valid for our case
        elif arg_ub <= 2.0:
            # find the according lower bound
            if arg_lb >= -2.0:
                # extract the parameters belonging to the interval [-2, 2]
                invalid_approx = False
                params_until_approx = func_dict[str(-2.0)][str(2.0)]
            elif arg_lb >= -5.0:
                # extract the parameters belonging to the interval [-5, -2]
                invalid_approx = False
                params_until_approx = func_dict[str(-5.0)][str(2.0)]
    # for sine and cosine we have to check that the bounds do not exceed 2*pi
    elif func_string in ["sin", "cos"]:
        # for (co)sine, we compute the bounds mod 2
        mod_arg_lb = arg_lb % (2 * np.pi)
        arg_length = arg_ub - arg_lb
        mod_arg_ub = mod_arg_lb + arg_length
        # compute the shift for 2 pi later an
        two_pi_shift = arg_lb // (2 * np.pi)

        # initialize the value of flipping as false
        flip = 0

        # check for validity in [k * pi/2, 2pi + k * pi/2], while respecting necessary re-computations by identities
        # here: case sine on [0, 2pi] and cosine on [pi/2, 5pi/2]
        if ((func_string == "sin" and mod_arg_lb >= 0 and mod_arg_ub <= 2 * np.pi + 1e-9) or
                (func_string == "cos" and mod_arg_lb >= np.pi / 2 - 1e-9 and mod_arg_ub <= 5 * np.pi / 2 + 1e-9)):
            # for sine, just take the corresponding bounds, for cosine shift for  - pi/2
            x_shift = 0 if func_string == "sin" else -np.pi / 2
            # flip if cosine is involved
            if func_string == "cos":
                flip = 1

            # initialize the bound keys and the medium bound
            keys = [str(0.0), str(np.round(np.pi, 6)), str(np.round(2 * np.pi, 6))]
            medium_bound = np.pi if func_string == "sin" else 3 * np.pi / 2

        # here: case sine on [3pi/2, 7pi/2] or cosine on [0, 2pi]
        elif ((func_string == "sin" and mod_arg_lb >= 3 * np.pi / 2 - 1e-9 and mod_arg_ub <= 7 * np.pi / 2 + 1e-9) or
                (func_string == "cos" and mod_arg_lb >= 0 and mod_arg_ub <= 2 * np.pi + 1e-9)):
            # for sine, just take the corresponding bounds, for cosine shift for  - pi/2
            x_shift = -2 * np.pi if func_string == "sin" else -np.pi / 2
            # due to a trigonometric identity, we have to flip the result for cosine
            if func_string == "cos":
                flip = 1

            # initialize the bound keys and the medium bound
            keys = [str(np.round(-np.pi/2, 6)), str(np.round(np.pi/2, 6)), str(np.round(3 * np.pi / 2, 6))]
            medium_bound = 5 * np.pi / 2 if func_string == "sin" else np.pi

        # here: case sine on [pi, 3pi] or cosine on [3pi/2, 7pi/2]
        elif ((func_string == "sin" and mod_arg_lb >= np.pi - 1e-9 and mod_arg_ub <= 3 * np.pi + 1e-9) or
                (func_string == "cos" and mod_arg_lb >= 3 * np.pi / 2 - 1e-9 and mod_arg_ub <= 7 * np.pi / 2 + 1e-9)):
            # for sine, we have a shift of -pi, for cosine, -3pi/2
            x_shift = -np.pi if func_string == "sin" else -3 * np.pi / 2

            # due to a trigonometric identity, flip parabolas for sine
            if func_string == "sin":
                flip = 1

            # initialize the bound keys and the medium bound
            keys = [str(0.0), str(np.round(np.pi, 6)), str(np.round(2 * np.pi, 6))]
            medium_bound = 2 * np.pi if func_string == "sin" else 5 * np.pi / 2

        # here: case sine on [pi/2, 5pi/2] or cosine on [pi, 3pi]
        elif ((func_string == "sin" and mod_arg_lb >= np.pi / 2 - 1e-9 and mod_arg_ub <= 5 * np.pi / 2 + 1e-9) or
              (func_string == "cos" and mod_arg_lb >= 0 and mod_arg_ub <= 2 * np.pi + 1e-9)):
            # for sine, we have a shift of -pi, for cosine, -3pi/2
            x_shift = -np.pi if func_string == "sin" else -3 * np.pi / 2

            # due to a trigonometric identity, flip parabolas for sine
            if func_string == "sin":
                flip = 1

            # initialize the bound keys and the medium bound
            keys = [str(np.round(-np.pi / 2, 6)), str(np.round(np.pi / 2, 6)), str(np.round(3 * np.pi / 2, 6))]
            medium_bound = 3 * np.pi / 2 if func_string == "sin" else np.pi

        # bounds are out of range for our approximation
        else:
            return invalid_approx, params_until_approx

        # check for smaller bounds in half of the interval
        if mod_arg_ub <= medium_bound + 1e-9:
            params_until_approx = func_dict[keys[0]][keys[1]]
        elif mod_arg_lb >= medium_bound - 1e-9:
            params_until_approx = func_dict[keys[1]][keys[2]]
        else:
            params_until_approx = func_dict[keys[0]][keys[2]]

        # update x shift according to modulo operation
        x_shift -= two_pi_shift * 2 * np.pi

        # update the parameters of the parabolas according to the bounds
        params_until_approx = shift_and_flip_parameters(params_until_approx, x_shift, flip)
        invalid_approx = False

    else:
        print(f"UNKNOWN FUNCTION STRING: {func_string}")
        exit()

    return invalid_approx, params_until_approx


def shift_and_flip_parameters(approx_dict, x_shift, flip):
    """ adjust the parameters in approx dict such that they incorporate an x shift and a flip (0, 1)"""
    # initialize the dictionary to return
    temp_dict = {}

    approx_directions = ["below", "above"]
    # iterate the two approximation directions and adjust the parameters
    for approx_index, approx_dir in enumerate(approx_directions):
        # initialize the new quadratic, linear and constant coefficients
        quads_new, lins_new, cons_new = [], [], []
        # extract the original coefficients
        quads_orig, lins_orig, cons_orig = approx_dict[approx_dir]
        # iterate coefficients and manipulate
        for quad, lin, con in zip(quads_orig, lins_orig, cons_orig):
            if flip == 1:
                quad *= -1
                lin *= -1
                con *= -1
            quads_new.append(quad)
            lins_new.append(2 * quad * x_shift + lin)
            cons_new.append(quad * x_shift ** 2 + lin * x_shift + con)
        # if flip is not desired (0), safe in temp dict with identical approximation direction, else other
        final_approx_index = (approx_index + flip) % 2
        temp_dict[approx_directions[final_approx_index]] = [quads_new, lins_new, cons_new]

    return temp_dict


if __name__ == "__main__":
    # extract list of instances as defined by the chosen attributes
    instances, primal_bounds, dual_bounds = get_instances_by_attribute(test_set, selected_attributes)

    primal_values = []
    deviating_instances = []

    # for each instance, parse, reformulate_1d and add/substitute paraboloids
    for index, inst in enumerate(instances):
        # print the current instance and then number of instances checked
        print(f"Adding paraboloids to {inst}\t\t | {index + 1}/{len(instances)}")

        # set up parser
        instance_path = os.path.join(ROOT_DIR, test_set, data_suffix, inst + f".{data_suffix}")
        parser = OSILParser(instance_path)
        # parse
        parser.parse()

        # convert parser to parser with 1-d reformulation
        n_new_variables, parser_1d = reformulate_osil_parser_to_1d(parser)

        # exchange possible functions with parabolas
        parser_para = substitute_paraboloids(parser_1d, n_new_variables, delete_orig=delete_orig)

        # initialize osil to pyomo
        converter_1d = OsilPyomoConverter(parser_1d)
        # set up model
        model_1d = converter_1d.setup_model()

        # same for paraboloid parser
        converter_para = OsilPyomoConverter(parser_para)
        model_para = converter_para.setup_model()
        if delete_orig:
            model_para_path = os.path.join(instances_path, f"{inst}_para.gms")
        else:
            model_para_path = os.path.join(instances_path, f"{inst}_both.gms")
        model_para.write(model_para_path)



