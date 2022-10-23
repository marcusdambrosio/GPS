import sys
import os
import re
import numpy as np

def fixed_data(method = 'by sat'):
    master = {}
    vals, labels = [], []
    for line in open('data.dat', 'r'):
        line = line.strip()
        val, label = [c.strip() for c in line.split('/= ')]
        vals.append(float(val)), labels.append(label)
    master['pi'] = vals[0]
    master['c'] = vals[1]
    master['R'] = vals[2]
    master['s'] = vals[3]

    labels = labels[4:]

    if method == 'by sat':
        satID = 0
        countReset = 9
        currData = []
        counter = 1
        for value in vals[4:]:
            currData.append(value)
            if counter%countReset == 0:
                master[satID] = currData
                satID += 1
                currData = []
            counter += 1
        
        for key, val in master.items():
            if type(val) != float:
                master[key] = [np.array(val[0:3]), np.array(val[3:6]), val[6], val[7], val[8]]
            
    else:
        for i, value, in enumerate(vals[4:]):
            if labels[i][0] == 'u' or labels[i][0] == 'v':
                if labels[i][1] == '1':
                    label = labels[i][0] + 'x'
                elif labels[i][1] == '2':
                    label = labels[i][0] + 'y'
                else:
                    label = labels[i][0] + 'z'
                master[f'{label}_{labels[i][-1]}'] = value
            else:
                satID = re.findall('\d+' , labels[i])[0]
                master[f'{labels[i].split(" ")[0]}_{satID}'] = value

    return master
