import numpy as np
import matplotlib.pyplot as plt
import json
import math
import copy

from numpy.matrixlib.defmatrix import matrix

import model
import modelMM
import solver
import solverMM


with open('layersMM.json', 'rt', encoding='UTF8') as json_file:
    layers = json.load(json_file)

properties = {"pile": {"E": 210*1e3, "I": 226000000, "A": 19100}}

# Hem = 6.9
modelObj = modelMM.makeModelObj(layers, 15)
with open("modelObjMM.json", "w", encoding='UTF8') as json_file:
    json.dump(modelObj, json_file)


y = []
for node in modelObj["node"]:
    y.append(node["y"])
xStan = np.zeros(len(y))
dispExcel = np.ones(len(y)) * -1.75
matrixObj = solverMM.makeMatrix(modelObj, properties)



disp = matrixObj["dispVec"][0]
q = 1*matrixObj["qVec"][0]

matrixJson = {"dispVec":disp.tolist(), "qVec": q.tolist()}
with open("matrixObjMM.json", "w", encoding='UTF8') as json_file:
    json.dump(matrixJson, json_file)


qLim = max(abs(max(q)), abs(min(q)))*1.3
dispLim = max(abs(max(disp)), abs(min(disp)))*1.3
dispLim = 1.75*1.5 if dispLim <= 1.75 else dispLim


fig, ax = plt.subplots(2)
ax[0].plot(xStan, y, label='$pile$')
ax[0].plot(q, y, label='$pressure$')
ax[0].set_xlabel('Pressure ($N/mm^2$)')
ax[0].set_ylabel('Depth ($mm$)')
ax[0].legend()
# ax[0].set_xlim(qLim, -qLim)
ax[0].set_xlim(38, -38)
ax[0].set_ylim(-8000, 2000)
ax[0].grid(True)

ax[1].plot(xStan, y, label='$Pile$')
ax[1].plot(disp, y, label='$disp$')
ax[1].plot(dispExcel, y, label='$excel$')
ax[1].set_xlabel('Displacement ($mm$)')
ax[1].set_ylabel('Depth ($mm$)')
ax[1].legend()
ax[1].set_xlim(-dispLim, dispLim)
# ax[1].set_xlim(-0.0025, 0.0025)
ax[1].set_ylim(-8000, 2000)
ax[1].grid(True)

plt.show()
