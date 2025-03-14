# Parabolic Approximation and Relaxation for MINLP

This repository includes the code used in the article *Parabolic Approximation and Relaxation for MINLP* by A. Göß, R. Burlacu, and A. Martin (2025+). 
If you use use this code in your work, please cite the article according to the corresponding section below.

A preprint will be available soon.

## Set-up

We encourage a set-up via an [Anaconda][1] environment.
Python version we used: Python 3.11.7

For code usage install the following packages:
- numpy (tested v1.26.4) 
- pyomo (tested v6.7.3)
- pandas (tested v2.2.2)
- gurobi (tested v11.0.3, only for catching errors, can be commented out)

Numpy and pyomo are relevant for the paraboloid computation, see directory para_approximation. 
For the paraboloid relaxation step, pandas is additionally required in the instance extractor, see directory para_relaxation.

For the solution of the created problems, we leverage [GAMS][2]. For this, you need a valid license. As we used [Gurobi][3], you may also require a license for this. There are academic licenses for both.

As a test set we considered the MINLPLib instances that are available in OSIL format, see [instances][4]. If the you want to run the code properly, download the instances and store them in directory minlplib/osil/. 
The provided CSV-file with information about the MINLPLib instances may need to be updated.

## Usage

Note: You may need to export your PYTHONPATH variable to the repository in order to ensure proper functionality of the relative imports, e.g., `export PYTHONPATH=/path/to/repository/paraboloids`.
This can solve problems like `ModuleNotFoundError: No module named ...`.

In order to run the approximations by paraboloids activate the environment if one is created (`conda activate myenv`) and run `python main.py my_run_setting.json` from within the para_approximation directory.
The program searches for run settings in `para_approximation/settings/run_settings/`, so create one or use `test_exact.json` or `test_inexact.json` for testing.
You can find numerous examples in the `run_setting` directory.

For creating instances with para and both type with parabolic approximations, also activate the environment as above and run `python write_para_relaxation.py` from within the para_relaxation directory. 
This writes the adjusted problem in GAMS format which can be solved by gams with any available solver.
We put the resulting files (except the `eg_...` ones due to their size) in the `instances/` directory for completeness.

The type of problem can be adjusted with the `delete_orig` flag in header of the file.
The instance names as well as other attributes can also be controlled in the file header.
The paraboloids are taken from a summary file in `para_approximation/out/para_parameters.json`.



## How to cite

As mentioned in the introducing comment, please consider the following article if you leverage the present code:

> A. Göß, R. Burlacu, and A. Martin (2025+).
> Parabolic Approximation and Relaxation for MINLP
> (link available soon)

You may also use the BibTeX entry: 
```bibtex
@article{Goess_Paraboloid_2024,
    author={G{\"o}{\ss}, Adrian and Burlacu, Robert and Martin, Alexander},
    title={Parabolic Approximation and Relaxation for MINLP},
    month={March},
    year={2025+},
    eprint={2407.06143},
    archivePrefix={arXiv},
    primaryClass={math.OC},
    url={https://arxiv.org/abs/2407.06143}
}
```

This project is licensed under the terms of the MIT license.

[0]: https://arxiv.org/abs/2407.06143

[1]: https://docs.anaconda.com/anaconda/install/

[2]: https://www.gams.com

[3]: https://www.gurobi.com

[4]: https://www.minlplib.org/index.html
