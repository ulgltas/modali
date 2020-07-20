#!/usr/bin/env python3
# -*- coding: utf-8; -*-

''' 
Copyright 2019 University of Liège

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 

modali.py
Python Modal solver
Huseyin Guner, Adrien Crovato 
'''

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import numpy as np
from scipy import integrate
from numpy.linalg import inv

# ----------------------------------------------------------------------
#  Modal solver
# ----------------------------------------------------------------------

class modali():
    """
    Modal solver
    """
    def __init__(self, m):
        # Say hi!
        print('Hi! I am a modal integrator!')
        print('Adrien Crovato and Huseyin Guner')
        print('ULiege, 2018-2019\n')

        # Get number of modes
        self.nModes = m
        print('Number of modes:', self.nModes)

    def setMatrices(self, _Mq, _Cq, _Kq):
        """Set the modal matrices and the number of modes
        """
        self.Mq = _Mq
        self.invMq = inv(self.Mq)
        self.Cq = _Cq
        self.Kq = _Kq
        print('Initialized modal matrices.')

    def readModes(self, fname):
        """Read the modes
        """
        # Read file
        print('Reading file:', fname)
        fl = open(fname, 'r')
        label = next(fl).split(',')
        fl.close()
        data = np.loadtxt(fname, delimiter=',', skiprows=1)
        # Store data
        self.nNodes = data.shape[0]
        self.nodalGlobalIndex = (data[:,0]).astype(int)
        self.nodalCoord_X = data[:,1]
        self.nodalCoord_Y = data[:,2]
        self.nodalCoord_Z = data[:,3]
        nodalMod_X = np.zeros((self.nNodes, self.nModes))
        nodalMod_Y = np.zeros((self.nNodes, self.nModes))
        nodalMod_Z = np.zeros((self.nNodes, self.nModes))
        for i in range(0, self.nModes):
            nodalMod_X[:,i] = data[:,4+3*i]
            nodalMod_Y[:,i] = data[:,5+3*i]
            nodalMod_Z[:,i] = data[:,6+3*i]
        print('Number of nodes:', self.nNodes)
        # Initialize modal matrix
        self.Phi = np.concatenate((nodalMod_X, nodalMod_Y, nodalMod_Z))    
        self.PhiT = self.Phi.transpose()
        print('Initialized mode shape matrix.')

    def setInitial(self, _xi, _vi, _fi):
        """Set the initial conditions (displacement, velocity and forces)
        """
        self.y0 = np.concatenate((_xi, _vi))
        self.dispX, self.dispY, self.dispZ = self.__getPhysicalDisp(self.y0[0:self.nModes])
        self.fq = _fi
        print('Set initial displacements:', self.y0[0:self.nModes])
        print('Set initial velocities:', self.y0[self.nModes:-1])
        print('Set initial forces:', self.fq)

    def setExtractor(self, _list):
        """Set an extractor list
        """
        self.extractor = {} # dictionnay mapping global to local index
        for gidx in _list:
            lidx = np.argwhere(self.nodalGlobalIndex == gidx)
            self.extractor[gidx] = lidx[0,0]
        print('Initialized extractor list with indices:', self.extractor)

    def updateLoads(self, _fx, _fy, _fz):
        """Set the load before the computation
        """
        f = np.concatenate((_fx, _fy, _fz)) # physical force vector
        self.fq = self.__getModalForce(f) # modal force vector

    def runStatic(self):
        """Run the static modal solver
        """
        print('Running static modal solver...')
        # Solve
        y = np.zeros((2, len(self.y0)))
        y[0, :] = self.y0 # store initial state
        for i in range(0, self.nModes):
            y[1, i] = self.fq[i] / self.Kq[i,i]
        self.y0 = y[1, :] # update initial state
        # Get physical physical displacements
        self.dispX, self.dispY, self.dispZ = self.__getPhysicalDisp(self.y0[0:self.nModes])
        # Printout
        print('{0:>5s}   {1:>12s}   {2:>12s}'.format('Dof', 'y_i', 'y_f'))
        for i in range(0, self.nModes):
            print('{0:5d}   {1:12.6f}   {2:12.6f}'.format(i, y[0, i], y[1, i]))
        print('')

    def runDynamic(self, t1, t2):
        """Run the dynamic modal sovler (time integration)
        """
        def f(t, y, self):   
            return np.concatenate([y[self.nModes:2*self.nModes], np.dot(self.invMq, (-np.dot(self.Cq, y[self.nModes:2*self.nModes]) - np.dot(self.Kq, y[0:self.nModes]) + self.fq))]) # equations of motion in modal coordinates

        print('Running dynamic modal solver...')
        # Sanity check
        if t2 <= 0 or t2 <= t1:
            raise Exception('final time ({0:f}) is either negative or leq. than initial time ({1:f})!\n'.format(t2, t1))
        # Solve
        t = np.array([t1, t2])
        y = np.zeros((len(t), len(self.y0)))
        y[0, :] = self.y0
            
        r = integrate.ode(f).set_integrator("dopri5") # explicit runge-kutta method of order (4)5 due to Dormand & Prince
        r.set_initial_value(self.y0, t1).set_f_params(self)
        for i in range(1, len(t)):
           y[i, :] = r.integrate(t[i])
           if not r.successful():
               raise RuntimeError("Could not integrate!\n")
        self.y0 = y[1, :]              
        # Get physical physical displacements
        self.dispX, self.dispY, self.dispZ = self.__getPhysicalDisp(self.y0[0:self.nModes])
        # Printout
        print('{0:>5s}   {1:>12s}   {2:>12s}   {3:>12s}   {4:>12s}'.format('Dof', 'y_i', 'y_f', 'y_i_dot', 'y_f_dot'))
        for i in range(0, self.nModes):
            print('{0:5d}   {1:12.6f}   {2:12.6f}   {3:12.6f}   {4:12.6f}'.format(i, y[0, i], y[1, i], y[0, i+self.nModes], y[1, i+self.nModes]))
        print('')

    def write(self, fname):
        """Write physical coordinates and modal data to disk
        """
        print('Writing data file:', fname+'.csv')
        file = open(fname+'.csv', 'w')
        file.write('index, x_coord, y_coord, z_coord, ')
        for j in range(0, self.nModes-1):
            file.write('dX_mode{0:d}, dY_mode{0:d}, dZ_mode{0:d}, '.format(j+1))
        file.write('dX_mode{0:d}, dY_mode{0:d}, dZ_mode{0:d}\n'.format(self.nModes))
        for i in range(0, self.nNodes):
            file.write('{0:d}, {1:f}, {2:f}, {3:f}, '.format(self.nodalGlobalIndex[i], self.nodalCoord_X[i]+self.dispX[i], self.nodalCoord_Y[i]+self.dispY[i], self.nodalCoord_Z[i]+self.dispZ[i]))
            for j in range(0, self.nModes-1):
                file.write('{0:f}, {1:f}, {2:f}, '.format(self.Phi[i,j], self.Phi[i+self.nNodes,j], self.Phi[i+2*self.nNodes,j]))
            file.write('{0:f}, {1:f}, {2:f}\n'.format(self.Phi[i,self.nModes-1], self.Phi[i+self.nNodes,self.nModes-1], self.Phi[i+2*self.nNodes,self.nModes-1]))
        file.close()

    def __getModalForce(self, f):
        """Transform a force vector to the modal space
        """
        return np.dot(self.PhiT, f)

    def __getPhysicalDisp(self, d):
        """Transform a displacement vector to the physical space
        """
        d = np.dot(self.Phi, d)
        dX = d[0:self.nNodes]
        dY = d[self.nNodes:2*self.nNodes]
        dZ = d[2*self.nNodes:3*self.nNodes]
        return dX, dY, dZ
