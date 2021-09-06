import numpy as np
import matplotlib.pyplot as plt
import json
import math
import copy

import model
import solver


with open('layers.json', 'rt', encoding='UTF8') as json_file:
    layers = json.load(json_file)

properties = {"pile": {"E":210*1e3, "I":226000000, "A":19100}}

# Hem = 6.9
modelObj = model.makeModelObj(layers,15)
pyCurve = modelObj["pyCurve"][60]
print(pyCurve)
print(pyCurve["coeff"])
print(len(pyCurve["coeff"])==1)
# print(data)
# coeff = data["coeff"]
# xLim = data["xLim"]
# pLim = data["pLim"]
# print(pLim)


dx = 0.00001
xlim = 0.05
xList = np.arange(-xlim,xlim+dx,dx)
p = []
for x in xList:
    p.append(solver.calcForceFromPYCurve(x, pyCurve))




fig, ax = plt.subplots()

ax.plot(xList, p, label='$Before$')
ax.set_xlabel('Displacement ($m$)')
ax.set_ylabel('Presssure ($kN/m^2/m$)')
ax.legend()
# plt.ylim(min(p)*2, max(p)*2)
plt.xlim(-xlim, xlim)
plt.grid(True)
plt.show()