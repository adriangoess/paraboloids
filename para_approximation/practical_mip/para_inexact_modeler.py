import os
import json
import numpy as np
import pyomo.environ as pyo

from para_approximation.utilities.templates import ParaboloidModelTemplate


class ParaboloidModel(ParaboloidModelTemplate):
    def __init__(self, func_settings_path, model_settings_path, number_paras, approx_below=True,
                 custom_max_n_delta_t=np.inf, custom_max_n_delta_d=np.inf):
        super().__init__(func_settings_path, model_settings_path, number_paras, approx_below,
                         custom_max_n_delta_t, custom_max_n_delta_d)

        # update delta and gamma (set to zero) in comparison to exact MIP approach
        self.delta = 0
        self.gamma = 0.5 * (self.delta/self.eps)

        # take the custom number of discretization points for the approximation as is
        max_n_delta_t = int(custom_max_n_delta_t)
        self.delta_t = (self.ub - self.lb) / max_n_delta_t
        self.t_indices = list(range(max_n_delta_t + 1))

        # take the custom number of discretization points for the function bound as is
        max_n_delta_d = int(custom_max_n_delta_d)
        self.delta_d = (self.ub - self.lb)/max_n_delta_d
        self.d_indices = list(range(max_n_delta_d + 1))

        # [unused] override custom values for the upper and lower para bound constraints
        self.custom_lower_para_delta = 0.9 * self.eps
        self.custom_upper_para_gamma = 0.25

    # overwritten method from superclass, as constraint (5d) is left out
    def _initialize_constraints(self, model):
        """ initializing the constraints of the model """
        # (5b): model lower bound of parabolas as big M formulation; p^l(t) >= f(t) - del - M1 * (1-s^l_t)
        model.approx_bound = pyo.Constraint(self.para_indices, self.t_indices, rule=self._lower_parabola_bound)
        # (5c): enforcing one of the containment binaries to be 1; sum_l s^l_t >= 1
        model.containment = pyo.Constraint(self.t_indices, rule=self._containment_sum)
        # (5d) is left out

        # (5e): model upper bound of parabolas; p^l(t) <= f(d) - gamma * eps
        model.upper_bound = pyo.Constraint(self.para_indices, self.d_indices, rule=self._upper_parabola_bound)
        # (5f): track violation of integral between parabola and function; v^l_d >= integral p - (f - gamma eps)
        model.violation_tracking = pyo.Constraint(self.para_indices, self.d_indices[:-1],
                                                  rule=self._violation_tracking)

        return model
