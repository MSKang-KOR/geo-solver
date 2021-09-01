import numpy as np
import matplotlib.pyplot as plt
import json
import math
import copy

# from numpy.matrixlib.defmatrix import matrix

import model
import modelAP
import solver
import solver2
import solver3
import solver4
import solverAP
# import solverMM


with open('layers2.json', 'rt', encoding='UTF8') as json_file:
    layers = json.load(json_file)

properties = {"pile": {"E": 210*1e3, "I": 226000000, "A": 19100}}

# Hem = 6.9
modelObj = model.makeModelObj(layers, 15)
with open("modelObjM.json", "w", encoding='UTF8') as json_file:
    json.dump(modelObj, json_file)
y = []
for node in modelObj["node"]:
    y.append(node["y"])
xStan = np.zeros(len(y))
dispExcel = np.ones(len(y)) * -1.75*1e-3
matrixObj = solver4.solver(modelObj, properties)
matrixObj2 = solver4.passive_active(modelObj, properties)


disp = matrixObj["dispVec"][0]
beforeDisp = matrixObj["beforeDispVec"][0]
q = 1*matrixObj["qVec"][0]
q_0 = 1*matrixObj2["qVec"][0]
q_0_excav = 1*matrixObj2["qVec2"][0]
q_active = 1*matrixObj2["qVec3"][0]
q_active_excav = 1*matrixObj2["qVec4"][0]
q_passive = 1*matrixObj2["qVec5"][0]
q_passive_excav = 1*matrixObj2["qVec6"][0]

matrixJson = {"dispVec": disp.tolist(), "qVec": q.tolist()}
with open("matrixObjM.json", "w", encoding='UTF8') as json_file:
    json.dump(matrixJson, json_file)


qLim = max(abs(max(q)), abs(min(q)))*1.3
dispLim = max(abs(max(disp)), abs(min(disp)))*1.3
dispLim = 0.00175*1.5 if dispLim <= 0.00175 else dispLim


fig, ax = plt.subplots(2)
ax[0].plot(xStan, y, label='$pile$')
# ax[0].plot(q_0, y, label='$p_0$')
# ax[0].plot(q_active, y, label='$p_a$')
# ax[0].plot(q_passive, y, label='$p_p$')
ax[0].plot(q, y, label='$pressure$')
# ax[0].plot(q_0_excav, y, label='$p_0excav$')
# ax[0].plot(q_active_excav, y, label='$p_a excav$')
# ax[0].plot(q_passive_excav, y, label='$p_p excav$')
ax[0].set_xlabel('Pressure ($kN/m^2$)')
ax[0].set_ylabel('Depth ($m$)')
ax[0].legend()
# ax[0].set_xlim(-qLim, qLim)
# ax[0].set_xlim(38, -38)
ax[0].set_ylim(-8, 2)
ax[0].grid(True)

ax[1].plot(xStan, y, label='$Pile$')
ax[1].plot(disp, y, label='$disp$')
ax[1].plot(beforeDisp, y, label='$disp_before$')
ax[1].plot(dispExcel, y, label='$excel$')
ax[1].set_xlabel('Displacement ($m$)')
ax[1].set_ylabel('Depth ($m$)')
ax[1].legend()
ax[1].set_xlim(dispLim, -dispLim)
# ax[1].set_xlim(-0.0025, 0.0025)
ax[1].set_ylim(-8, 2)
ax[1].grid(True)

plt.show()
