# modali
Structural dynamics solver based on modal integration  
Adrien Crovato and Huseyin Guner
ULiege, 2019

## Description
modali allows to compute the response of a structure subjected to an excitation based on its modal charasteristics.

modali solves the following equation:  
**M_q** q_ddot + **C_q** q_dot + **K_q** q = f_q,  
where **M_q**, **C_q** and **K_q** are the modal mass, damping and stiffness matrices, q is the modal displacements vector and f_q is the modal forces vector.
The modal displacements are then multiplied by the mode shape matrix to get the physical displacements.

## Features
* Static solver: solves **K_q** q = f_q for a diagonal stifness matrix
* Dynamic integrator: solves the full equation with Runge-Kutta order 4-5

## Usage
modali has been designed to be used with [CUPyDO](https://github.com/ulgltas/CUPyDO) but also works in standalone.  
Example/test cases are given under [tests](/tests), they can be run with
```bash
python run.py tests/testname.py
```
For coupled cases, plese refer to CUPyDO repository.

### Input
The following inputs are required to create a case:
* Initial modal displacements, velocities and loads (can be 0)
* Modal matrices formatted as numpy matrices
* Mode shapes in a .csv file. Each line of the file should contain an index, the x, y and z coordinates of a physical node, and the x, y, z coordinates of the m first mode shapes.
* Extractor index list (cannot be empty!)