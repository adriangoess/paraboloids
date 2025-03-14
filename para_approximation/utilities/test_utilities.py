import os
import sys
import json
import datetime
from definitions import *

run_settings_arguments = ["exact_mip", "function_setting_names", "model_setting_names", "approximation_directions",
                          "max_time_limit", "max_iterations", "log_search", "write_solution", "visual_output"]


def input_handler(first_argument):
    """ this function reads the first argument and tries to initialize the options from json """
    # first argument must end with .json
    assert first_argument.endswith(".json"), f"the argument {first_argument} must be of json type"
    # concatenate the argument and test for existence
    option_path = os.path.join(run_settings_path, first_argument)
    assert os.path.exists(option_path), f"path {option_path} does not exist"
    # load the json and chck for its structure
    temp_options = json.load(open(option_path))
    assert len(temp_options.keys()) == len(run_settings_arguments)
    for temp_key in run_settings_arguments:
        assert temp_key in temp_options.keys(), f"key {temp_key} must be present in {first_argument}"
        # additionally check for the type of the provided data
        value = temp_options[temp_key]
        if temp_key in ["exact_mip", "log_search", "write_solution", "visual_output"]:
            assert value in [0, 1], f"argument of {temp_key} must be 0 or 1"
        elif temp_key in ["max_time_limit", "max_iterations"]:
            assert isinstance(value, (int, float)) and temp_options[temp_key] >= 0, \
                f"argument of {temp_key} must be non-negative integer or float"
        elif temp_key in ["function_setting_names", "model_setting_names", "approximation_directions"]:
            assert isinstance(value, (list,)), f"argument of {temp_key} must be list"
            for element in value:
                if temp_key == "approximation_directions":
                    assert element in [0, 1], f"element of list in argument {temp_key} must be 0 or 1"
                else:
                    assert isinstance(element, (str,)), f"element of list in argument {temp_key} must be of type string"

    return temp_options


def get_function_string(setting_name):
    """extract the function string based on the function setting name"""
    assert "sine" in setting_name or "exp" in setting_name or "x3" in setting_name or "zigzag" in setting_name, \
        f"so far only implemented for sine, cosine, x^3, exp and zigzag"
    if setting_name.startswith("sine"):
        function_string = "sin"
    elif setting_name.startswith("zigzag"):
        function_string = "zigzag"
    elif "cosine" in setting_name:
        function_string = "cos"
    elif "x3" in setting_name:
        function_string = "x^3"
    else:
        function_string = "exp"
    return function_string


def shift_cons_for_violation(function_violation, f_p_distances, cons):
    """ given whether a function violation occurs in the parabolic approximation, shift the constant values """
    # initialize the shifted cons
    shifted_cons = cons.copy()

    # compute the distances but respect numerical issues
    distances = [max([0, dist - 1e-6]) if not function_violation else dist - 1e-6 for dist in f_p_distances]
    for k, para_constants in enumerate(cons):
        for l, entry in enumerate(para_constants):
            shifted_cons[k][l] += distances[k]

    return shifted_cons


def get_out_params(quads, lins, cons, approx_below, para_modeler, visional_output=False):
    """ flip the coefficients if necessary (approx below) and plot if wanted """
    out_quads = [(-1) ** (1 + approx_below) * quad[0] for quad in quads]
    out_lins = [(-1) ** (1 + approx_below) * lin[0] for lin in lins]
    out_cons = [(-1) ** (1 + approx_below) * con[0] for con in cons]

    if visional_output:
        function_to_plot = lambda x: (-1) ** (1 + approx_below) * para_modeler.f(x)
        # add plotting function

    return out_quads, out_lins, out_cons


def write_params(function_string, eps, func_setting_name, approx_below, exact, quads, lins, cons, times):
    """ write the parameters of the paraboloids and the times in separate json files """
    # initialize the data path as the out path only
    data_path = out_path

    # update the function string for path creation to not include the '^' sign
    temp_func_string = function_string
    if "^" in temp_func_string:
        temp_func_string = temp_func_string.replace("^", "")
    # update the data path, create if it does not exist
    data_path = os.path.join(data_path, temp_func_string)
    if not os.path.exists(data_path):
        os.mkdir(data_path)

    # create sub directory for epsilon
    temp_eps_string = str(eps)
    temp_eps_string = temp_eps_string.replace(".", "_")
    data_path = os.path.join(data_path, temp_eps_string)
    if not os.path.exists(data_path):
        os.mkdir(data_path)

    # create file prefix
    file_prefix = func_setting_name.replace(".json", "")
    if approx_below:
        file_prefix += "_below"
    else:
        file_prefix += "_above"
    if exact:
        file_prefix = "exact_" + file_prefix
    else:
        file_prefix = "inexact_" + file_prefix

    # label a successful run if quads is nonempty, initialize updated flag
    successful = len(quads) > 0
    updated = False

    if successful:
        # create parameters file name and check for existence (as well as number of entries)
        parameters_file_name = file_prefix + "_parameters.csv"
        parameters_path = os.path.join(data_path, parameters_file_name)
        # save whether this path does not exist
        not_existing = not os.path.exists(parameters_path)
        if not_existing:
            access_mode = "w+"
        else:
            access_mode = "r+"
        with open(parameters_path, access_mode) as parameters_file:
            # read all parameter lines
            all_lines = parameters_file.readlines()
            n_entries = len(all_lines)
            # only overwrite parameters if none there or more parabolas before than now
            if not_existing or len(quads) < n_entries:
                # re-set updated flag
                updated = True
                # write one line for one paraboloid
                for i in range(len(quads)):
                    # concatenate line and write it
                    line = f"{quads[i]},{lins[i]},{cons[i]}\n"
                    parameters_file.write(line)

    # if the run was unsuccessful or there was no update, add the time entry, otherwise write a new file
    times_file_name = file_prefix + "_times.json"
    times_path = os.path.join(data_path, times_file_name)
    existing_times = os.path.exists(times_path)

    if updated or not existing_times:
        access_mode = "w+"
    else:
        access_mode = "r+"

    # open/create the times file
    with open(times_path, access_mode) as times_file:
        if updated or not existing_times:
            times_dict = {}
        else:
            times_dict = json.load(times_file)

        # create the current time stamp as key
        time_stamp = datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
        times_dict[time_stamp] = times

        times_file.seek(0)
        json.dump(times_dict, times_file, indent=4)
        times_file.truncate()

