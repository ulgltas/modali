# Modal Solver
Structural dynamics solver based on modal decomposition
Huseyin Guner and Adrien Crovato
ULiege, 2019

## Description
Modal Solver allows to compute the response of a structure subjected to an excitation based on its modal charasteristics.

Modal Solver solves the following equation: **M_q** q_ddot + **C_q** q_dot + **K_q** q = f_q,
where **M_q**, **C_q** and **K_q** are the modal mass, damping and stiffness matrices, q is the modal displacements vector and f_q is the modal forces vector.
The modal displacements are then multiplied by the mode shape matrix to get the physical displacements

## Features
* Static solver: solves **K_q** q = f_q for a diagonal stifness matrix
* Dynamic solver: solves the full equation with Runge-Kutta order 4-5

## Usage
Modal Solver has been designed to be used with [CUPyDO](https://github.com/ulgltas/CUPyDO) but can be easily improved to work in standalone.

### Configuration script
```python
# Config file
    fconfig = os.path.join(os.path.abspath(os.path.dirname(__file__)),'modes.csv')
    # Paramters
    initialModalDisp = np.array([1., 0.]) # initial displacement
    initialModalVel = np.zeros(2, dtype=float) # intitial velocities
    modalMass = np.diag([1., 2.]) # mass matrix
    modalDamping = np.zeros((2, 2), dtype=float) # damping matrix
    modalStiffness = np.diag([3., 4.]) # stiffness matrix
    nModes = initialModalDisp.shape[0] # number of modes

    # Solver
    solver = ms.ModalSolver(nModes) # intialize solver
    solver.setMatrices(modalMass, modalDamping, modalStiffness) # initialize matrices
    solver.readModes(fconfig) # initial mode shape matrix
    solver.setInitial(initialModalDisp, initialModalVel) # set initial condition on displacements and velocities
    solver.setExtractor([1]) # set a list of extractor by providing the global indeces of the (physical) node 
```

### Input
* Modal matrices formatted as numpy matrices
* Mode shapes in a .csv file. Each line of the file should contain an index, the x, y and z coordinates of a physical node, and the x, y, z coordinates of the m first mode shapes.
```csv
"Global_Index", "x_coord", "y_coord", "z_coord", "dX_mode1", "dY_mode1", "dZ_mode1", "dX_mode2", "dY_mode2", "dZ_mode2"
1, 1.e-01, 1.e-01, 1.-01, 1e-01, 1e-01, 1e-01, 1e-01, 1e-01, 1e-01
```