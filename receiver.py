import sys
import os
import re
import numpy as np
from sympy import Symbol, diff, lambdify
import numdifftools as nd
import select
import time
from scipy.optimize import fsolve, newton_krylov, broyden1, broyden2, anderson, newton
import datetime as dt

'''initiate constants'''
def constants():
    master = {}
    vals, labels = [], []
    for line in open('data.dat', 'r'):
        line = line.strip()
        val, label = [tempVar.strip() for tempVar in line.split('/= ')]
        vals.append(float(val)), labels.append(label)

    return vals[0], vals[1], vals[2], vals[3]

pi, c, R, s = constants()

class Receiver:
    def __init__(self, currLocation, solveType):
        #format input in a more usable way
        self.currLocation = [[float(tempVar2.strip()) for tempVar2 in tempVar3.strip().split(' ')] for tempVar3 in currLocation]
        self.solveType = solveType  #select solving method to optimize runtime

        self.initialize_X()
        self.output = []

    def f(self, X):    #create the function f that we want to minimize the gradient of
        f = 0
        for sat in np.arange(1, len(self.currLocation)):
            # sum(||X_s+1 - Xv|| - ||X_s - XV|| - c * (X_s - X_s+1)
            f += (np.linalg.norm(np.array(self.currLocation[sat][2:]) - X, ord=2) -
                   np.linalg.norm(np.array(self.currLocation[sat - 1][2:]) - X, ord=2) -
                   c * (self.currLocation[sat - 1][1] - self.currLocation[sat][1])) ** 2
        return f

    def gradient_fNewton(self, X):
        #nd.Gradient is generates a gradient matrix for the function f
        gradF = nd.Gradient(self.f)
        self.gradF = gradF(X)
        #performing this again creates our jacobian for gradient(f)
        self.grad2F = nd.Gradient(nd.Gradient(self.f))(X)
        return [self.gradF, self.grad2F]

    def gradient_fBroyden(self, X):
        self.gradF = nd.Gradient(self.f)(X)
        return self.gradF

    def solve_s (self, s, arg):
        J, F = arg #jacobian and gradient(f)
        return np.dot(J, s) + F #this is just the residual form of the function sJ = -F that will be fed to a solver

    def newton_solve(self,Xguess):
        #initialize
        X0 = Xguess
        S = 10
        nIter = 0.
        #loop until the newton step S is less than .01 (1 cm)
        tic = dt.datetime.today()
        while np.linalg.norm(S) > .01:
            #calculate F and J, first and second gradients
            F, J = self.gradient_fNewton(X0)
            #solve for S using built in solver for ease
            S = fsolve(self.solve_s, X0, args = [J, F])
            #calculate Xn+1 and reset Xn (which is X0 in my case)
            XN = X0 + S
            X0 = XN
            nIter += 1
            if  (dt.datetime.today() - tic).seconds > 30: #runtime error handling
                return nIter

        return XN

    def solve_X(self, Xguess):
        #initiate solver based on Xguess
        if self.solveType == 'newton':
            self.X = self.newton_solve(Xguess)
            try:    #runtime error handling
                len(self.X)
            except:
                return self.X
        elif self.solveType == 'broyden':
            self.X = broyden1(self.gradient_fBroyden, Xguess)
        else:
            print('solver selection error')
        print(self.X)
        #calculate t based on solution for X
        self.t = np.linalg.norm(np.array(self.currLocation[0][2:]) - self.X, ord=2) / c + self.currLocation[0][1]
        return [self.t, np.dot(self.R3back(), self.X), self.X]

    def initialize_X(self):
        #this just creates the cartesian coordinates for lamppost b12 to use as an initial guess
        self.tv, latVec, NS, longVec, EW, h = 12123.0, np.array([40, 45, 55]).astype(float), \
                                              float(1), np.array([111, 50, 58]).astype(float), float(-1), float(1372)
        lat = 2 * pi * NS * (latVec[0] / 360 + latVec[1] / (360 * 60) + latVec[2] / (360 * 60 ** 2))
        long = 2 * pi * EW * (longVec[0] / 360 + longVec[1] / (360 * 60) + longVec[2] / (360 * 60 ** 2))
        self.X0 = np.array([(R + h) * np.cos(lat) * np.cos(long), (R + h) * np.cos(lat) * np.sin(long), (R + h) * np.sin(lat)])
        return np.dot(self.R3(), self.X0)

    def R3(self): #time correction matrix
        alpha = pi*self.tv/s
        self.r3 = np.array([[np.cos(alpha), -np.sin(alpha), 0],
                            [np.sin(alpha), np.cos(alpha), 0],
                            [0, 0, 1]])
        return self.r3

    def R3back(self): #time correction matrix but backward to t=0
        alpha = -2*pi*self.t/s
        self.r3back = np.array([[np.cos(alpha), -np.sin(alpha), 0],
                            [np.sin(alpha), np.cos(alpha), 0],
                            [0, 0, 1]])
        return self.r3back

loop = 0
output = []
currLocation = []
init = True
solveTypeInit = True
#start log
log = open('receiver.log', 'w')
log.write('RECEIVER INPUT\n\n')

#loop through lines of standard input
for line in sys.stdin.readlines():
    log.write(line)
    if loop == 0:   #if its the first line, initiate the satellite ID
        satID = int(line.split(' ')[0])
        currLocation.append(line)
        loop = -1
        continue
    newsatID = int(line.split(' ')[0]) #set most recent satellite ID
    if newsatID < satID: #check if the new satellite ID is less than the old one, this indicates a new "stop" in the trip is reached
        if solveTypeInit:
            tempCurrLoc = [[float(tempVar2.strip()) for tempVar2 in tempVar3.strip().split(' ')] for tempVar3 in
                             currLocation]

            solveType = 'newton' if round(tempCurrLoc[0][1], 2) == (152.21 or 102122.92) else 'broyden'
            solveTypeInit = False
        else:
            pass

        rec = Receiver(currLocation, solveType) #initiate the receiver object
        if init:
            Xguess = rec.initialize_X()   #if its the first stop in the trip, initiate Xguess as lamppost b12
            init = False
        else:
            Xguess = output[-1][-1] #if it ins't the first stop, set Xguess as the previous location

        # collect outout from solve_X and reinitiate variables for next stop in the trip
        newOutput = rec.solve_X(Xguess)
        if type(newOutput) == float:
            output.append(newOutput)
            break
        else:
            output.append(newOutput)
        currLocation = [line]
        satID = newsatID
    else: #add the next line of satellite data if it isn't a new stop and reset satID variable
        currLocation.append(line)
        satID = newsatID

#calculate final stop's location

if type(output[-1]) == float:    #runtime error handling
    pass

else:
    rec = Receiver(currLocation, solveType)
    Xguess = output[-1][-1]
    output.append(rec.solve_X(Xguess))


def rad_to_deg(radians): #function for converting from radians to degrees
    radians = -(radians - pi) if radians > pi else radians
    angle = (radians * 180 / np.pi) % 360
    angle = angle - 180 if angle > 180 else angle
    degrees, t =  '{:f}'.format(angle).split('.')
    t = float('.' + t)
    t *= 60
    minutes, seconds = '{:f}'.format(t).split('.')
    seconds = round(float('.' + seconds) * 60, 2)
    return [degrees, minutes, seconds]


fullOutput = []

for out in output: #loop through output from receiver class calculations to convert to lat long and altitude
    if type(out) == float:
        fullOutput.append('Timeout Error: Convergence not reached in less than 30 seconds with ' + str(round(out)) + ' iterations.')
        break
    t = out[0]
    x, y, z = out[1]
    h = round(np.sqrt(x**2 + y**2 + z**2) - R, 2)
    NS = 1 if z >= 0 else -1

    if x**2 + y**2 != 0:    #all of this logic represents piecewise functions that are found in HW 1 to convert from x,y,z to lat long
        lat = np.arctan(z / np.sqrt(x**2 + y**2))
    else:
        if z > 0:
            lat = pi / 2
        else:
            lat = -pi / 2
    if x > 0 and y > 0:
        long = np.arctan(y / x)
    elif x < 0:
        long = pi + np.arctan(y / x)
    elif x > 0 and y < 0:
        long = 2 * pi + np.arctan( y / x)
    else:
        sys.stdout.write('error' + '\n')
    EW = 1 if long < pi else -1
    latDeg = rad_to_deg(lat)
    longDeg = rad_to_deg(long)
    #add the line of output in the correct format
    fullOutput.append(str(t) + ' ' + str(latDeg[0]) + ' ' + str(latDeg[1]) + ' ' + str(latDeg[2]) + ' ' + str(NS) \
                  + ' ' + str(longDeg[0]) + ' ' + str(longDeg[1]) + ' ' + str(longDeg[2]) + ' ' + str(EW) + ' ' + str(h))

log.write('\nRECEIVER OUTPUT\n\n')
for out in fullOutput:  #write the output as standard output
    log.write(out +'\n')
    sys.stdout.write(out + '\n')

log.close()