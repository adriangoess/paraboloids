import os
from bisect import bisect_right
import numpy as np
import pyomo.environ as pyo
from pyomo.opt import SolverStatus, TerminationCondition

from definitions import *
from para_relaxation.pyomo_solver import PyomoSolver
from para_approximation.utilities.result_checker import check_coefficient_validity

zigzag_params_path = os.path.join(ROOT_DIR, "para_approximation", "utilities", "zigzag_parameters.csv")
assert os.path.exists(zigzag_params_path)

file = open(zigzag_params_path, "r")
lines = file.readlines()
assert len(lines) >= 2
x_values = lines[0]
y_values = lines[1]
file.close()

x_values = [float(val) for val in x_values.split(",")]
y_values = [float(val) for val in y_values.split(",")]


def zigzag(x):
    """ zigzagging function based on zigzag_parameters.csv """
    # manual serialization of function by checking isinstance ndarray
    if isinstance(x, (np.ndarray, )):
        for val in x:
            assert min(x_values) <= val <= max(x_values)
        iteration_values = x
    else:
        assert min(x_values) <= x <= max(x_values)
        iteration_values = [x]

    return_values = []
    for val in iteration_values:
        # compute the index above the current value
        up_index = bisect_right(x_values, val)

        # check whether the insertion index is the length of the list, then the upper bound is evaluated
        if up_index < len(x_values):
            # compute a linear interpolation between the left and the right index
            slope = (y_values[up_index] - y_values[up_index - 1]) / (x_values[up_index] - x_values[up_index - 1])
            return_values.append(y_values[up_index - 1] + (val - x_values[up_index - 1]) * slope)
        else:
            # return the corresponding value to the upper bound
            return_values.append(y_values[-1])

    # manual serialization
    if isinstance(x, (np.ndarray, )):
        return np.array(return_values)
    else:
        return return_values[0]


def Zigzag(lb, ub):
    """ the integral of the zigzagging function from lb to ub """
    assert min(x_values) - 1e-10 <= lb <= ub <= max(x_values) + 1e-10

    # initialize the integral value
    integral = 0

    # compute the indices of the x_values left and right of [lb, ub]
    left_index = bisect_right(x_values, lb) - 1
    right_index = min(bisect_right(x_values, ub), len(x_values) - 1)

    # iterate consecutive index pairs, compute the integral by the area of the trapezoidal shape
    iter_index = left_index
    while iter_index + 1 <= right_index:
        # left bound is right most value of lb and current index and right bound accordingly
        left_x = max(x_values[iter_index], lb)
        right_x = min(x_values[iter_index + 1], ub)

        # update the integral value with the trapezoid area
        integral += (y_values[iter_index] + y_values[iter_index + 1]) * (right_x - left_x) * 0.5

        iter_index += 1

    return integral


def check_zigzag_violation(quads, lins, conss, domain, viol_type, eps, dim=1, print_result=True, approx_below=True):
    """ compute the maximal function or approximation violation of the zigzagging function by checking each interval """
    # initialize the distances of zigzag to p (min f - p) and the violation flag
    minimal_f_p_distances = [np.inf] * len(quads) if viol_type == "function" else []
    violation = False

    assert viol_type in ["function", "approximation"], "can compute violation w.r.t. function or approximation only"
    assert dim == 1, f"zigzagging is a 1 dimensional function"
    assert check_coefficient_validity(quads, lins, conss, dim)
    assert len(domain) == 2, f"domain must have two entries, lb and ub"
    assert len(domain[0]) == 1, f"domain entries must have dimension 1"
    assert len(domain[1]) == 1, f"domain entries must have dimension 1"
    assert domain[0] < domain[1], f"ub has to exceed lb"

    # compute the left and right most index to be considered for domain
    left_index = bisect_right(x_values, domain[0][0]) - 1
    right_index = min(bisect_right(x_values, domain[1][0]), len(x_values) - 1)
    assert left_index < right_index, "there needs to be at least one nonempty interval"

    # iterate the x intervals by iterating the left index as long as it is smaller than the right one
    iter_index = left_index
    while iter_index < right_index:
        # compute the lower and upper bound for the current check
        lb = max(x_values[iter_index], domain[0][0])
        ub = min(x_values[iter_index + 1], domain[1][0])

        # compute the slope and the constant for the current linear function
        slope = (y_values[iter_index + 1] - y_values[iter_index]) / (x_values[iter_index + 1] - x_values[iter_index])
        constant = y_values[iter_index] - x_values[iter_index] * slope

        # establish the model to check whether the paraboloids violate the current linear function or its approxmation
        if viol_type == "function":
            tmp_viol, tmp_dists = _function_violation_model(quads, lins, conss, lb, ub, slope, constant, approx_below)
            for para_index in range(len(quads)):
                minimal_f_p_distances[para_index] = min(minimal_f_p_distances[para_index], tmp_dists[para_index])
        else:
            tmp_viol, tmp_dist = _approx_violation_model(quads, lins, conss, lb, ub, slope, constant, eps, approx_below)
            minimal_f_p_distances.append(tmp_dist)

        # update the tracking measure
        violation = violation or tmp_viol

        iter_index += 1

    if min(minimal_f_p_distances) < 0 and print_result:
        print(f"{viol_type.capitalize()} violation: {min(minimal_f_p_distances)}")
    elif print_result:
        print(f"No {viol_type} violation!")

    if viol_type == "function":
        return violation, minimal_f_p_distances
    else:
        return violation, min(minimal_f_p_distances)


def _function_violation_model(quads, lins, conss, lb, ub, slope, constant, approx_below=True):
    """
        setup and solve MIPs which checks the violation of the paraboloids defined by (quads, lins, conss)
        on [lb, ub] of the linear function slope * x + constant
    """
    # initialize the f-p distance (= min f - p) and the violation
    minimal_f_p_distances = []
    violation = False

    # define the linear function
    func = lambda x: (-1) ** (1 + approx_below) * (slope * x + constant)

    # iterate the paraboloids and compute the maximal violation by computing min( lin - p_i ) over lb, ub
    for i, (quadratics, linears, constants) in enumerate(zip(quads, lins, conss)):
        # initialize the model
        m = pyo.ConcreteModel()

        # initialize the variable
        m.x = pyo.Var([0], within=pyo.Reals, bounds=[lb, ub])

        # initialize the paraboloid function; directly model linearly to avoid GAMS errors if necessary
        if quadratics[0] == 0.0:
            paraboloid = linears[0] * m.x[0] + constants[0]
            # use gurobi here, as GAMS fails to solve LPs with SCIP
            solver = PyomoSolver("gurobi-direct")
        else:
            paraboloid = quadratics[0] * m.x[0]**2 + linears[0] * m.x[0] + constants[0]
            solver = PyomoSolver("scip")

        # setup the objective
        m.obj = pyo.Objective(expr=func(m.x[0]) - paraboloid, sense=pyo.minimize)

        # solve the model
        obj, results = solver.solve_model(m, False)

        minimal_f_p_distances.append(obj)
        if obj < 0:
            violation = True

    return violation, minimal_f_p_distances


def _approx_violation_model(quads, lins, conss, lb, ub, slope, constant, eps, approx_below=True):
    """
        setup and solve a MIP which checks whether the approximation by paraboloids defined by (quads, lins, conss)
        on [lb, ub] fulfills the epsilon approximation guarantee w.r.t. the linear function slope * x + constant.
        So, max_i p_i >= linear - epsilon on [lb, ub] is wanted. For this, one computes the minimum of
        (y - linear + eps) where y >= p_i for each i.
    """
    # define the linear function
    sign_of_f = -1 if approx_below else 1
    func = lambda x: sign_of_f * (slope * x + constant)

    # initialize the model
    m = pyo.ConcreteModel()

    # initialize the variables
    m.x = pyo.Var([0], within=pyo.Reals, bounds=[lb, ub])
    m.y = pyo.Var("y", within=pyo.Reals)

    # for each paraboloid, model y >= p_i over [lb, ub]
    m.constraints = pyo.ConstraintList()
    for i, (quadratics, linears, constants) in enumerate(zip(quads, lins, conss)):
        paraboloid = quadratics[0] * m.x[0]**2 + linears[0] * m.x[0] + constants[0]
        m.constraints.add(expr=m.y["y"] >= paraboloid)

    # model the objective
    m.obj = pyo.Objective(expr=m.y["y"] + func(m.x[0]) + eps, sense=pyo.minimize)

    # solve the model
    solver = PyomoSolver("scip")
    solver.solve_model(m, False)
    if m.obj() < 0:
        violation = True
    else:
        violation = False

    return violation, m.obj()

