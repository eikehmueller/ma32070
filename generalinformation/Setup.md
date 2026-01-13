----

*&#169; Eike Mueller, University of Bath 2025. These notes are copyright of Eike Mueller, University of Bath. They are provided exclusively for educational purposes at the University and are to be downloaded or copied for your private study only. Further distribution, e.g. by upload to external repositories, is prohibited. html generated with [pandoc](https://pandoc.org/) using [easy-pandoc-templates](https://github.com/ryangrose/easy-pandoc-templates) under the [GPL-3.0.1 license](https://github.com/ryangrose/easy-pandoc-templates?tab=GPL-3.0-1-ov-file#readme)*

----

# General setup
The purpose of this notebook is to describe the setup of the software environment required to develop code for MA32070.

The crucial ingredients are

* a **terminal** to access the command line
* a working **Python** environment
* several **Python packages** for numerical computations in particular the [petsc4py](https://petsc.org/release/petsc4py/) numerical linear algebra library
* a modern code **editor** such as [VSCode](https://code.visualstudio.com/)
* the MA32070 [**finite element library**](https://github.bath.ac.uk/em459/finiteelements)

At the beginning of the course, please work through the two steps below to get set up. Once you have done this, review the section on the directory layout for the course material.

## Noteable
The default setup is to use the **Noteable** environment, which every student on this course has access to. This provides a suitable Python installation with all required packages. You can access Noteable through the link on moodle. You should choose the "Standard Python 3 with VS-Code editor" notebook, which allows you to edit files with the built-in VScode editor.

## Working on your own computer
You are welcome to set up a suitable Python environment on your own computer and work with this if you prefer. This is relatively straightforward on a Linux or Mac computer, but more challenging on Windows. Below are some instructions that might help you with this. **Note that if you choose to use your own computer it is your responsibility to set it up correctly, and we can only provide limited support for this. You are also responsible for backing up any code you write. Problems with your own computer are not a valid reason for requesting a coursework extension.**

## Using the command line
Scientific software is usually developed in a command line environment: the user interacts with the system by typing commands at a command line prompt, also known as a "shell". This allows the user to do things like creating files and directories and executing Python scripts. Linux and Mac computers provide a Unix-like environment which you can access by launching a terminal session. On Windows systems, you will have to first install the [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/) (see below), which gives you access to a Linux-like system in which you can install the Python packages required for this unit.

If you have never used the command line, familiarise yourself with the basics by working through sections 1 to 3 of the [Software Carpentry tutorial on the Unix shell](https://swcarpentry.github.io/shell-novice/).

### Accessing a terminal
On Noteable, you can launch a terminal by choosing "File" $\rightarrow$ "New" $\rightarrow$ "Terminal" from the menu bar or by clicking the square "Terminal" button in the "Launcher" tab on the right.

In the VSCode environment (both on Noteable and in a standalone installation on your own computer) you can access the command line through a terminal in one of the tabs at the bottom of the window.

## Directory layout
It is strongly recommended that you carry out all work related to this course in a dedicated directory (or folder), which will be referred to as `ma32070/` in the following. Depending on your setup, this folder will contain additional subdirectories, see below for more details on how to organise your code.

# Step 1: Setup of the required Python environment

## Option 1: Noteable
Access Noteable through the link on moodle and choose the "Standard Python 3 with VS-Code editor" notebook.

That's it.

You can now edit files with the VS-Code editor and run them either in the VS-Code terminal or by launching a separate terminal session.

## Option 2: Using your own computer
Please note that the following instructions are rough guidelines, and they will likely have to be adapted to your particular computer. Setting up the Python programming environment is more straightforward on Linux and Mac computer which provide access to a command line environment with an up-to-date Python version. On Windows computers, you will first have to install the WSL, see below. 


Download and install [VSCode](https://code.visualstudio.com/) (or use your own favourite editor).

To access the command line environment, launch a terminal and install the additional required Python packages, as described in the following.

### Python
Most computers nowadays come with a working Python (version 3) installation. Double check that you have a recent version (anything from 3.12 upwards should be fine) by running

```
python3 --version
```

in the command line prompt.

To install the additional Python packages it is strongly recommended that you either do this via `pip install` inside a virtual environment (Option 2a) or use [Anaconda](https://www.anaconda.com/) (Option 2b).

### Option 2a: Virtual environment
Using a dedicated Python [virtual environment (venv)](https://docs.python.org/3/library/venv.html) allows you to install Python packages in a controlled way, without affecting the system-wide Python installation. If you have a working Python installation, you should be able to create a virtual environment with

```
python3 -m venv PATH
```

where `PATH` is the name of the directory in which to create the venv. You probably want to use a subdirectory of your `ma32070/` directory, e.g. `ma32070/venv`. Once created, you can activate the virtual environment with

```
source PATH/bin/activate
```

Hence, whenever you work on this course, you should activate the virtual environment as soon as you have opened a terminal session.

#### Installation of required Python packages
Inside the virtual environment, the required Python packages can be installed with

```
python3 -m pip install PACKAGE
```

where `PACKAGE` is the name of the package to be installed:

* [numpy](https://numpy.org/) for numerical linear algebra
* [matplotlib](https://matplotlib.org/) for plotting and visualisation
* [pytest](https://docs.pytest.org/en/stable/) for testing
* [petsc4py](https://petsc.org/release/petsc4py/) for linear solvers and preconditioners later in this course.

Install these packages inside the virtual environment you created (make sure that you activated the environment first).

#### Installing petsc4py
The installation of petsc4py can be tricky since it requires building the PETSc library itself. For this, you might have to install additional (non-Python) packages such as `libblas-dev` and `liblapack-dev`.

### Option 2b: Anaconda
Follow the [instructions](https://www.anaconda.com/docs/getting-started/anaconda/install) to install and set up Anaconda on your computer. 

#### Installation of required Python packages
Launch a terminal and make sure that Anaconda has been activated: you should see `(base)` in the command line prompt.

Use the `conda create` command to set up a dedicated environment for this course. When doing so, make sure that it includes the packages `petsc petsc4py pytest mpi4py` from the `conda-forge` channel, usually this means including `-c conda-forge petsc petsc4py pytest mpi4py` at the end of the command.

For example, to create a local environment called `petsc_sandbox` on a Linux machine, run the following command:

```
conda create -y --prefix MA32070_DIR/petsc_sandbox -c conda-forge petsc petsc4py pytest mpi4py
```

where `MA32070_DIR` needs to be replaced by the name of your `ma32070/` directory. You will need to look at the [documentation of `conda create`]((https://docs.conda.io/projects/conda/en/stable/commands/create.html) to work out the exact command on your computer. To use the created environment, it first needs to be activated with 

```
conda activate MA32070_DIR/petsc_sandbox
```

every time you work on this course (this will happen automatically if you append this line to your `${HOME}/.bashrc` file on Linux and Mac computers). You can check whether the environment has been activated successfully by looking at the command prompt: this should contain the name of the environment, i.e. `petsc_sandbox` in the above example.

## Windows computers
If you are working on a Windows computer, you should install the [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/) by following the [WSL installation instructions](https://learn.microsoft.com/en-us/windows/wsl/install), see also the [troubleshooting section](https://learn.microsoft.com/en-us/windows/wsl/troubleshooting?source=recommendations) if you run into problems. Use a suitable Linux distribution, such as Ubuntu-24.04. Make sure you activate the WSL after installation (you will need to create a user account and password for this). Inside the WSL Linux system, you will then have to set up git and Python and install the necessary packages, as described above.

# Step 2: Installing the finite element library
The finite element library provides some Python code that we will use throughout the semester.

Check that the `git` [version control](https://git-scm.com/) tool is available by typing

```
git --version
```

in the command line. This should produce output like this (any reasonably recent version of git will be ok):

```
eike@eike-linux:~$ git --version
git version 2.43.0
```

Change to the `ma32070/` directory and then run the following commands:

```
git clone https://github.bath.ac.uk/em459/finiteelements.git
cd finiteelements
python3 -m pip install --editable .
python -m pytest -v
```

This will create a subdirectory `ma32070/finitelements`.

If everything worked correctly, you can use functionality from this library in your own code. For example, you could write the following Python script (in any directory):

```Python
from fem.utilitymeshes import rectangle_mesh

mesh = rectangle_mesh(Lx=1.0, Ly=2.0, nref=2)
print (f"number of cells    = {mesh.ncells}")
print (f"number of vertices = {mesh.nvertices}")
```

to create a rectangle mesh and print out the number of cells and vertices.

Please **never** edit any code in `ma32070/finiteelements` or add any other code to the directory. If you do this by accident, delete the directory and install the finite element library again as described above.

### Updating the finite element library
During the semester, additional functionality and model solutions will be added to the finite element library. To access this code, you will need to change to the `ma32070/finiteelements` directory and run

```
git pull
```

Don't forget to change to the directory with your own code afterwards.

# Directory layout
If you worked through the steps above, your `ma32070/`  directory should contain the following subdirectories:

* `ma32070/finiteelements` for the finite element library
* a subdirectory with the virtual environment or Anaconda environment:
  - `ma32070/venv` (if you installed the Python packages in a virtual environment) or 
  - `ma32070/petsc_sandbox` (if you use Anaconda)

It is suggested that you create an additional directory
* `ma32070/workspace` for your own work

In summary, the `ma32070/` folder should contain exactly three subdirectories and you should only ever edit or add code in `ma32070/workspace`.

