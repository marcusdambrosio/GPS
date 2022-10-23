import sys
import os
import re
import numpy as np


def constants():
    master = {}
    vals, labels = [], []
    for line in open('data.dat', 'r'):
        line = line.strip()
        val, label = [tempVar.strip() for tempVar in line.split('/= ')]
        vals.append(float(val)), labels.append(label)
    master['pi'], master['c'], master['R'], master['s']  = vals[0], vals[1], vals[2], vals[3]
    satID = 0
    countReset = 9
    currData = []
    counter = 1
    for value in vals[4:]:
        currData.append(value)
        if counter % countReset == 0:
            master[satID] = currData
            satID += 1
            currData = []
        counter += 1

    for key, val in master.items():
        if type(val) != float:
            master[key] = [np.array(val[0:3]), np.array(val[3:6]), val[6], val[7], val[8]]

    return master

class SingleSatellite:
    def __init__(self, satID, datafile, tv):
        self.aboveHorizon = False
        self.ID = satID
        self.u, self.v, self.p, self.h, self.theta = datafile[satID]
        self.X = (R + self.h) * (self.u * np.cos(2 * pi * tv / self.p + self.theta) + self.v * np.sin(2 * pi * tv / self.p + self.theta))

    def position(self, t):
        self.X = (R + self.h) * (self.u * np.cos(2 * pi * t / self.p + self.theta) + self.v * np.sin(2 * pi * t / self.p + self.theta))
        return self.X


class Satellite:
    def __init__(self, fullInput):
        '''VEHICLE DATA'''
        self.fullInput = fullInput
        self.tv = fullInput[0]
        self.R3()
        self.cart_at_t0()
        self.cart_at_t()

    def satellite_info(self):
        self.satDict = {}
        self.output = []
        for key, val in datafile.items():
            if type(val) != float:
                currSat = SingleSatellite(key, datafile, self.tv)
                ts, Xs = self.satellite_output(currSat)
                if np.dot(self.Xv.T, (Xs - self.Xv)) > 0:
                    self.output.append([key, ts, Xs])

        return self.output
    
    def rad_to_deg(self):
        #conversion and modulo 360 to prevent continually counting rotations
        angle = (radians * 180 / np.pi ) % 360
        #split into whole degrees and the decimal portion
        degrees, t = str(angle).split('.')
        degrees = int(degrees)
        t = float('.' + t)
        #convert to minutes and split off the seconds like above
        t *= 60
        minutes, seconds = str(t).split('.')
        minutes = int(minutes)
        seconds = round(float('.' + seconds) * 60)
        return [degrees, minutes, seconds]

    def cart_at_t0(self):
        tv, latVec, NS, longVec, EW, h  = self.fullInput[0], self.fullInput[1:4], \
                                          self.fullInput[4], self.fullInput[5:8], self.fullInput[8], self.fullInput[9]
        lat = 2 * pi * NS * (latVec[0] / 360 + latVec[1] / (360 * 60) + latVec[2] / (360 * 60 ** 2))
        long = 2 * pi * EW * (longVec[0] / 360 + longVec[1] / (360 * 60) + longVec[2] / (360 * 60 ** 2))
        self.X0 = np.array([(R + h) * np.cos(lat) * np.cos(long), (R + h) * np.cos(lat) * np.sin(long), (R + h) * np.sin(lat)])

    def R3(self):
        alpha = 2*pi*self.fullInput[0]/s
        self.r3 = np.array([[np.cos(alpha), -np.sin(alpha), 0],
                            [np.sin(alpha), np.cos(alpha), 0],
                            [0, 0, 1]])

    def cart_at_t(self):
        self.Xv = np.dot(self.r3, self.X0)

    def satellite_output(self, sat):
        tk = self.tv
        tk1 = self.tv + 1
        diff = np.abs(tk - tk1)
        while diff > .01/c:
            tk1 = self.tv - np.sqrt(np.sum((sat.position(tk) - self.Xv) ** 2)) / c
            diff = np.abs(tk - tk1)
            tk = tk1
        return tk, sat.position(tk)


    def satellite_output_newtons(self, sat):
        # tk = self.tv - np.linalg.norm(sat.position(self.tv) - self.Xv, ord = 2)
        tk = self.tv
        tk1 = self.tv + 1
        diff = np.abs(tk - tk1)
        while diff > (.01 / c):
            f = np.dot((sat.position(tk) - self.Xv) , (sat.position(tk) - self.Xv)) - c ** 2 * (self.tv - tk) ** 2
            fprime = 4*pi * (R + sat.h) / sat.p * np.dot((sat.position(tk) - self.Xv) ,
                                                         (-sat.u * np.sin(2*pi*tk/sat.p + sat.theta) + sat.v * np.cos(2*pi*tk/sat.p + sat.theta))) + 2 * c**2 * (self.tv - tk)

            tk1 = tk - f / fprime
            diff = (tk1 - tk)
            tk = tk1

        sys.exit()
        return tk, sat.position(tk)



'''constants'''
datafile = constants()
R, pi, s, c = datafile['R'], datafile['pi'], datafile['s'], datafile['c']

log = open('satellite.log', 'w')
log.write('CONSTANTS\n\n')
for key, val in datafile.items():
    log.write(str(key) + ' : ' + str(val) + '\n')
log.write('\n')

newFile = True


for line in sys.stdin.readlines():
    if len(line)<10:
        continue
    # inputFile = open(line.strip(), 'r')
    # inputData = []
    # for dataLine in inputFile:
    #     if len(dataLine)<20:
    #         continue
    #     inputData.append([float(c.strip()) for c in dataLine.split(' ')])

    file = [float(tempVar2.strip()) for tempVar2 in line.strip().split(' ')]

    log.write('\nInput = ' + str(file) + '\n\n')


    S = Satellite(file)
    for out in S.satellite_info():
        newFile = False
        out += out[2].tolist()
        del out[2]
        output = ' '.join([str(tempVar3) for tempVar3 in out])
        print(output)
        sys.stdout.write(output + '\n')
        log.write('Satellite ' + str(out[0]) + ' output = ' + str(output) + '\n')

log.close()
sys.exit()

