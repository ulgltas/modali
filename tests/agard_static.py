#! /usr/bin/env python
# -*- coding: utf-8; -*-

import numpy as np
import modali as ms
import os.path

def main():
    # Config file
    fconfig = os.path.join(os.path.abspath(os.path.dirname(__file__)),'models/agard_modes.csv')
    # Paramters
    initialModalDisp = np.zeros(4, dtype=float)
    initialModalVel = np.zeros(4, dtype=float)
    initialModalLoads = np.array([0.406085, -0.338318, -0.013586, -0.101239])
    modalMass = np.diag([2.9107e-04, 8.3181e-05, 1.7447e-04, 3.4281e-05])
    modalDamping = np.zeros((4, 4), dtype=float)
    modalStiffness = np.diag([1.0468, 5.3468, 17.3717, 12.9114])
    nModes = initialModalDisp.shape[0]

    # Initialize solver
    solver = ms.modali(nModes)
    solver.setMatrices(modalMass, modalDamping, modalStiffness)
    solver.readModes(fconfig)
    solver.setInitial(initialModalDisp, initialModalVel, initialModalLoads)
    solver.setExtractor([16, 13808])
    print''

    # Run solver
    solver.runStatic()

    # Save data to disk
    solver.write('agard_sol')
    print''

    # Print reference solution for check (obtained with CUPyDO)
    print 'Check:'
    print 'y0:', solver.y0[0], '?~', 0.387929
    print 'y1:', solver.y0[1], '?~', -0.063275
    print 'y2:', solver.y0[2], '?~', -0.000782
    print 'y3:', solver.y0[3], '?~', -0.007841

    # eof
    print ''

if __name__ == '__main__':
    main()