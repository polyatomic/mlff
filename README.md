# Body-ordered Kernel Force Field.

* [mathlib](mathlib): utility library of mathematical algorithms.
* [chemlib](chemlib): cheminformatics utility library.
* [util](util): general utility library.
* [bokfit](bokfit): tool for body-ordered kernel machine learning force field training.
* [bokpredict](bokpredict): tool for prediction.
* [examples](examples): example data for model development.

## General idea

The method works as a typical kernel method such as kernel ridge regression (KRR) but kernel function centers are selected from independently formed grid instead
of putting them directly onto training points. This gives advantage in many-body expansion force field generation by allowing to select
centers independently for different many-body terms.

## Prerequisites

* Intel oneAPI 2021 for LAPACK. Other LAPACK implementations may be used as well but have not been tested.
* Python 3.9 or newer for Python wrappers (currently in early development stage).

## Compile instructions

1. Open solution file `bokfit.sln` with Visual Studio 2019 or newer.
1. Build the projects (Ctrl-B). This is going to produce `bokfit` and `bokpredict` executables in the subdirectories of x64.
To install chemutil module for Python, modify lists `sfc_include_dirs` and `sfc_libraries` in `setup.py` if
needed to point to correct locations of oneAPI headers and libraries in your system. Module can then be installed
with the command:

`pip install <path name of folder containing pyproject.toml>`

## Example usage

Force field can be generated from the file [benzene.xyz](http://www.quantum-machine.org/gdml/data/xyz/md17_benzene2017.zip) and example input files in [examples/training](examples/training) by using the following steps:

1) Preliminary selection of grid points based on the closeness to training data points.  
`bokfit -p p1.txt > res1.log`

2) Calculation of descriptor matrix.  
`bokfit -p p2.txt > res2.log`

3) Ordering of training points according to their calculated leverages in different descriptor blocks (12 blocks and 12 different orderings).  
`bokfit -p p3.txt > res3.log`

4) Creating of the sketch of the whole descriptor matrix, i.e. eliminating redundant descriptors from the matrix.  
`bokfit -p p4.txt > res4.log`

5) Building the model as the linear combination of descriptors by stepwise selection of the highest leverage points into the training set.  
`bokfit -p p5.txt > res5.log`

Steps 1-4 are performed in unsupervised manner, i.e. they don't use any energy data and are used for the generation of sufficiently diverse training set.
Training energies in `e.txt` are given in meV. The model is trained with 21292 selected geometries and the mean absolute error for the whole geometry set of size 627983 will
be around 1 meV (see `res5.log`). Only 745 grid points (kernel function centers) from different many-body terms are used for prediction which is far less than the size of the
training set.
Trained force field data files are given in [examples/prediction](examples/prediction) and prediction for the geometries in `benzene.xyz` can be run as:

`bokpredict -p predict.txt > res.log`

The resulting file `res.log` contains predicted energy values and predicted-observed comparison statistics.
