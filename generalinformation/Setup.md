<div align="center">
  <p style="font-size:32px;">MA32070 setup instructions</p>
</div>

----

*&#169; Eike Mueller, University of Bath 2025. These notes are copyright of Eike Mueller, University of Bath. They are provided exclusively for educational purposes at the University and are to be downloaded or copied for your private study only. Further distribution, e.g. by upload to external repositories, is prohibited. html generated with [pandoc](https://pandoc.org/) using [easy-pandoc-templates](https://github.com/ryangrose/easy-pandoc-templates) under the [GPL-3.0.1 license](https://github.com/ryangrose/easy-pandoc-templates?tab=GPL-3.0-1-ov-file#readme)*

----

# General setup

## Using the command line
Scientific software is usually developed in a command line environment: the user interacts with the system by typing commands at a command line prompt, also known as a "shell". This allows the user to do things like creating files and directories and executing Python scripts. Linux and Mac computers provide a Unix-like environment which you can access by launching a terminal session. On Windows systems, you will have to first install the [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/) (see below), which gives you access to a Linux-like system.

If you have never used the command line, familiarise yourself with the basics by working through sections 1 to 3 of the [Software Carpentry tutorial on the Unix shell](https://swcarpentry.github.io/shell-novice/).

## Directory structure
Make sure that you organise your code in a sensible way. It is recommended that you create a single directory, called `ma32070` for this course. Inside this directory, create a subdirectory `ma32070/workspace` for any code you develop.

# Required software

## Notable
Open a terminal window and create a new conda environment as follows:

```
conda create --name ${HOME}/petsc_sandbox -y
conda activate ${HOME}/petsc_sandbox
conda install -c conda-forge petsc -y
conda install -c conda-forge petsc4py -y
conda install -c conda-forge pytest -y
```

Create a directory in which you want to install the finite elements package, for example `git_workspace`. Change to this directory and then run the following commands:

```
git clone https://github.com/eikehmueller/finiteelements.git
cd finiteelements
python3 -m pip install --editable .
python -m pytest -v
```

## Windows computers
If you are working on a Windows computer, you should install the [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/) by following the [WSL installation instructions](https://learn.microsoft.com/en-us/windows/wsl/install), see also the [troubleshooting section](https://learn.microsoft.com/en-us/windows/wsl/troubleshooting?source=recommendations) if you run into problems. Use a suitable Linux distribution, such as Ubuntu-24.04. Make sure you activate the WSL after installation (you will need to create a user account and password for this). Inside the WSL Linux system, you might then have to install additional packages, see below.

## Python
We will use Python (version 3) for software development. Most computers nowadays come with a working Python3 installation, check that you have a recent version (anything from 3.12 upwards should be fine) by running

```
python3 --version
```

in the command line prompt.

## Virtual environment
It is strongly recommended to carry out all work for this course inside a dedicated Python [virtual environment (venv)](https://docs.python.org/3/library/venv.html): this will allow you to install Python packages in a controlled way, without affecting the system-wide Python installation. If you have a working Python installation, you should be able to create a virtual environment with

```
python3 -m venv PATH
```

where `PATH` is the name of the directory in which to create the venv. You probably want to use a subdirectory of your `ma32070` directory, e.g. `ma32070/venv`. Once created, you can activate the virtual environment with

```
source PATH/bin/activate
```

Hence, whenever you work on this course, you should do this at the beginning of each session:

1. Open a terminal to get a command line prompt
2. If you are using Windows, launch the WSL by typing `wsl` in the prompt
3. Navigate to the `ma32070` course directory
4. Activate the virtual environment

## Required Python packages
We will require several Python packages, which can usually be installed with

```
python3 -m pip install PACKAGE
```

where `PACKAGE` is the name of the package to be installed:

* [numpy](https://numpy.org/) for numerical linear algebra
* [matplotlib](https://matplotlib.org/) for plotting and visualisation
* [pytest](https://docs.pytest.org/en/stable/) for testing
* [petsc4py](https://petsc.org/release/petsc4py/) for linear solvers and preconditioners later in this course.

Install these packages inside the virtual environment you created. If this is not possible (this might be the case if you are using the WSL), you will likely have to pass the additional flag `--break-system-packages` to `python3 -m pip install`.

### Installing petsc4py
The installation of petsc4py can be tricky since it requires building the PETSc library itself. For this, you might have to install additional (non-Python) packages such as `libblas-dev` and `liblapack-dev`. We will only need petsc4py towards at the end of the course, so you do not have to get it working in the first week.

## Additional software
We will use the following additional software:

* [git](https://git-scm.com/) for version control. This is a command line tool which should be installed on most systems. When writing your own code, you are welcome to use graphical interfaces such as [Sourcetree](https://www.sourcetreeapp.com/) or [GitKraken](https://www.gitkraken.com/).
* [Visual Studio (VS) Code](https://code.visualstudio.com/) for editing source code.
* [Paraview](https://www.paraview.org/) for visualising results.

Installation should be straightforward with the provided download links.

# Installing the finite element library
You should clone the provided finite element library by running

```
git clone git@github.com:eikehmueller/finiteelements.git
```

in your `ma32070` directory. Once you have done this, change to `ma32070/finitelements` and run

```
python3 -m pip install --editable .
```

A good test of verifying that you have correctly installed all required software is to run

```
pytest -v
```

inside this directory.

If everything worked correctly, you can use functionality from this library in your own code. For example, you could write the following Python script (in any directory):

```Python
from fem.utilitymeshes import RectangleMesh

mesh = RectangleMesh(Lx=1.0, Ly=2.0, nref=2)
print (f"number of cells    = {mesh.ncells}")
print (f"number of vertices = {mesh.nvertices}")
```

to create a rectangle mesh and print out the number of cells and vertices.

# Final directory layout
If you followed all steps above, your `ma32070` directory should now contain the following subfolders:

* `ma32070/venv` for the virtual environment
* `ma32070/finiteelements` for the finite element library
* `ma32070/workspace` for your own work

# Troubleshooting
If you run into problems, please first try to resolve these yourself: in many cases it helps to google the relevant error messages or suitable keywords that summarise the problem. Searching pages such as [stackoverflow](https://stackoverflow.com/) can also be helpful since it is likely that someone else will have encountered the same problem.

If you do get stuck and need help from your course tutor, please make sure that you provide sufficient information with your request. For this, try to answer the following questions:

* What exactly is the problem? Which code failed? Be as specific as possible.
* What exactly did you do to trigger the error? Describe any previous steps you might have gone through.
* What is the exact error message you get? Please copy and paste raw text instead of sending screenshots.
* What did you expect to happen? What did you actually observe?
* Under which environment (operating system, used Python version, ...) did the error occur?
* Which other software have you installed successfully?
* Which steps have you already taken to resolve the issue?
* What do you think the problem is? Is this consistent with what you actually observe?

You might want to review the [Appendix 1. Getting Help](https://object-oriented-python.github.io/a1_help.html) from [David Ham's excellent book](https://object-oriented-python.github.io/index.html) on how to ask good questions.