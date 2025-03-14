import os

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

para_approximation_path = os.path.join(ROOT_DIR, "para_approximation")
out_path = os.path.join(para_approximation_path, "out")
settings_path = os.path.join(para_approximation_path, "settings")
func_settings_path = os.path.join(settings_path, "func_settings")
model_settings_path = os.path.join(settings_path, "model_settings")
run_settings_path = os.path.join(settings_path, "run_settings")

para_relaxation_path = os.path.join(ROOT_DIR, "para_relaxation")
instances_path = os.path.join(para_relaxation_path, "instances")
instances_out_path = os.path.join(para_relaxation_path, "out_files")

minlplib_path = os.path.join(ROOT_DIR, "minlplib")

# maximal number of binary and continuous variables
n_max_binaries = 1_000
n_max_continuous = 200_000