import os
import sys
import json
import numpy as np
from definitions import *
from para_approximation.utilities.test_utilities import input_handler, get_function_string, get_out_params, write_params

from para_approximation.exact_mip.test_para_computation import ExactParaboloidHandler
from para_approximation.practical_mip.practical_para_computation import InexactParaboloidHandler


if __name__ == "__main__":
    # main file requires the name of a run setting to be executed, so at least one argument
    assert len(sys.argv) > 1, f"main.py requires at least one argument, that is, the run setting"
    options = input_handler(sys.argv[1])

    if options["exact_mip"]:
        para_handler_class = ExactParaboloidHandler
    else:
        para_handler_class = InexactParaboloidHandler

    # loop over the approximation directions, True = "from below", False = "from above"
    for approx_below in options["approximation_directions"]:
        # loop over the function settings and extract possibly necessary information
        for function_setting_name in options["function_setting_names"]:
            # create function settings path, test for existence, initialize function settings
            function_setting_path = os.path.join(func_settings_path, function_setting_name)
            assert os.path.exists(function_setting_path), f"path {function_setting_path} does not exist"
            function_setting = json.load(open(function_setting_path))

            # extract function string
            function_string = get_function_string(function_setting_name)
            n_paraboloids = 1

            # loop over model settings, they are assumed to be sorted in decreasing eps
            for model_setting_name in options["model_setting_names"]:
                # initialize writing and plotting value
                visual_output = bool(options["visual_output"])
                write_output = bool(options["write_solution"])

                if visual_output:
                    print("-"*100)

                # create model settings path
                model_setting_path = os.path.join(model_settings_path, model_setting_name)
                assert os.path.exists(model_setting_path), f"path {model_setting_path} does not exist"
                model_setting = json.load(open(model_setting_path))

                # initialize the iteration counter, bounds for the number of paraboloids, the coefficients
                n_iterations = 0
                n_paraboloids_lb = n_paraboloids
                n_paraboloids_ub = np.inf
                quads, lins, cons = [], [], []
                out_quads, out_lins, out_cons = [], [], []

                # initialize the paraboloid handler
                para_handler = para_handler_class(function_string, function_setting_path, model_setting_path, n_paraboloids, approx_below)

                # iterate the search procedure
                while n_iterations < options["max_iterations"]:
                    # make a solution step dependent on which approach to take
                    end_of_loop, next_loop, results = para_handler.solution_step(quads, lins, cons,
                                                                                 options["max_time_limit"],
                                                                                 visual_output)
                    # terminate loop if detected in solution step
                    if end_of_loop:
                        break
                    elif next_loop:
                        n_iterations += 1
                        continue

                    # make a validation step dependent on which approach to take
                    next_loop, valid_solution, quads, lins, cons = para_handler.validation_step(results)

                    # continue if end of loop
                    if next_loop:
                        n_iterations += 1
                        continue

                    # update the search bounds dependent on strategy
                    end_of_loop = para_handler.update_search_bounds(valid_solution)

                    # save the parameters if valid
                    if len(quads) > 0 and valid_solution:
                        if len(out_quads) == 0 or len(quads) < len(out_quads):
                            out_quads, out_lins, out_cons = get_out_params(
                                quads, lins, cons, approx_below, para_handler.para_modeler, visual_output
                            )

                    # end of loop signals an end of the binary search in this case
                    if end_of_loop:
                        break
                    n_iterations += 1

                # save the parameters for current model + function if new best
                write_params(function_string, model_setting["eps"], function_setting_name, approx_below,
                             options["exact_mip"], out_quads, out_lins, out_cons, para_handler.times)

