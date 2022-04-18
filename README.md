<span style="font-size:larger;">BASIC Manual</span>
========

<!--- [//]: # (Badges) 
[![GitHub Actions Build Status](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/actgpaw/workflows/CI/badge.svg)](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/actgpaw/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/ACTgpaw/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/ACTgpaw/branch/master) -->

# Table of Contents

- [About BASIC](#about-basic)
    - [Highlights in v0.3](#highlights-in-basic-v0.3)
    - [Available convergence test](#available-convergence-test)
    - [Workflow of BASIC](#workflow-of-basic)
    - [License and credits](#license-and-credits)
- [Download and Install](#download-and-install)
- [Use BASIC](#use-basic)
- [Code Structure](#code-structure)
- [Author Info](#author-info)

# About BASIC

BASIC (**B**ulk, **A**dsorption, **S**urface **C**alculator with automat**I**c **C**onvergence test) is a python package designed to minimize the effort to compute bulk, surface and adsorption model using DFT by **streamlining script preparation**, **robotizing convergence test** and **automating data storage process**. 

## Highlights in BASIC v0.3

* Support [GPAW](https://wiki.fysik.dtu.dk/gpaw/#) DFT code
* Renewed computation logic:
    * Modulized computation process. Support single computation, multiple calculator parameters convergence test, slab layer and area convergence test.
    * Automated directory creation process. Streamline the preparation process further. 
    * Calculator parameters convergence process is more generalize, making way for other DFT code.
* Support bottom fix and center fix slab relaxation.

## Available Convergence Test

* Calculator parameters (available for bulk, surface and adsorption model)
    * grid spacing (h)
    * k Points (kpts)
    * kpts density (kpts_density)
    * smearing width (occupation_width)
* Slab super cell size
    * Number of layer
    * Area of slab

## Workflow of BASIC

## License and credits

# Download and install

## Github (for developers)

If you would like to contribute to the development of this software,
BASIC can be installed via a clone from Github. First, you'll need to clone the
github repo to your local machine (or wherever you'd like to use BASIC) using
`git clone`. Once the repo has been cloned, you can install BASIC as an editable
package by changing into the created directory (the one with `setup.py`) and installing
via: 
```
pip install -e .
```

# Use BASIC




# Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.5.
