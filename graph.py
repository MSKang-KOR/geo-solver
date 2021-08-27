import numpy as np
import matplotlib.pyplot as plt
import json
import math
import copy

import model

# import pressure as p


with open('layers.json', 'rt', encoding='UTF8') as json_file:
    layers = json.load(json_file)

Hem = 6.9
obj = model.makeModelObj(layers,15)
y = []
for node in obj["node"]:
    y.append(node["y"])
p01 = obj["pressure"]["backSide"]["rest"]
pa1 = obj["pressure"]["backSide"]["active"]
pq1 = obj["pressure"]["backSide"]["passive"]
p02 = obj["pressure"]["excavationSide"]["rest"]
pa2 = obj["pressure"]["excavationSide"]["active"]
pq2 = obj["pressure"]["excavationSide"]["passive"]

print(len(obj["node"]))
print(len(obj["pressure"]["backSide"]["active"]))
print(len(obj["pressure"]["backSide"]["active"]))

xlim = 400
# xlim = max(max(PbA), max(PbP), max(PeA), max(PeA))

fig, ax = plt.subplots()
# # 배면측
ax.plot(p01, y, label='$P_01$')
# ax.plot(pa1, y, label='$P_a1$')
# ax.plot(pq1, y, label='$P_p1$')
# # 굴착측
ax.plot(p02, y, label='$P_02$')
# ax.plot(pa2, y, label='$P_a2$')
# ax.plot(pq2, y, label='$P_p2$')

ax.set_xlabel('Pressure ($kN/m^2$)')
ax.set_ylabel('Depth ($m$)')
ax.legend()
plt.ylim(min(y), max(y))
# plt.xlim(-xlim, xlim)
plt.xlim(-xlim, xlim)
plt.grid(True)
plt.show()
    
