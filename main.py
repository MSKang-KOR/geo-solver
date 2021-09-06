import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import json
import math
import copy

# from numpy.matrixlib.defmatrix import matrix

import model
import model2
import modelAP
import solver
import solver2
import solver3
import solver4
import solver5
import solverAP


properties = {"pile": {"E": 210*1e3, "I": 226000000, "A": 19100}}
with open('input.json', 'rt', encoding='UTF8') as json_file:
    input = json.load(json_file)

# ------------------------------------------------------------

modelList = model2.modelGenerator(input)
modelObj = modelList[0]["model"]
matrixObj2 = solver5.passive_active(modelObj, properties)
q_0 = 1*matrixObj2["qVec"][0]

y = []
for node in modelObj["node"]:
    y.append(node["y"])
xStan = np.zeros(len(y))

result = []
hWall = input["wall"]["h"]
hWater = input["water"]["h"]

struts = []
for i in range(len(input["step"])):
    step = input["step"][i]
    if len(step["strut"]):
        for name in step["strut"]:
            struts.append(name)
hStrutList = [input["strut"][name]["h"] for name in struts]
isStrut = False

fig, ax = plt.subplots(len(modelList), 4)
fig.set_size_inches(15, 8)
for i in range(len(modelList)):
    step = input["step"][i]
    hex = step["excavation"] if step["excavation"] != 0 else hex
    isStrut = True if len(step["strut"])!=0 else isStrut

    modelObj = modelList[i]["model"]
    isSolve = modelList[i]["solve"]

    if i == 0:
        matrixObj = solver5.solver(input, modelObj, solve=isSolve)
    else:
        qVec0 = copy.deepcopy(matrixObj["qVec"])
        kList0 = copy.deepcopy(matrixObj["kList"])
        dispVec0 = copy.deepcopy(matrixObj["dispVec"])
        matrixObj = solver5.solver(
            input, modelObj, qVec0, kList0, dispVec0, solve=isSolve)
    result.append(matrixObj)

    step = result[i]
    q = step["qVec"][0]
    disp = step["dispVec"][0]
    V = step["vVec"][0]
    M = step["mVec"][0]

    qLim = max([abs(min(q)), abs(max(q))]) * 1.3
    ax[i, 0].plot(xStan, y, label='$pile$', color='g', linewidth=2)
    ax[i, 0].plot(q, y, label='$pressure$', color='r')
    ax[i, 0].set_ylabel('Depth ($m$)')
    ax[i, 0].set_xlim(-qLim, qLim)
    ax[i, 0].set_ylim(-hWall, 0)
    ax[i, 0].grid(True)

    dispLim = max([abs(min(disp)), abs(max(disp))]) * 1.3
    ax[i, 1].plot(xStan, y, label='$pile$', color='g', linewidth=2)
    ax[i, 1].plot(disp, y, label='$Displacement$', color='r')
    ax[i, 1].set_xlim(dispLim, -dispLim)
    ax[i, 1].set_ylim(-hWall, 0)
    ax[i, 1].grid(True)

    VLim = max([abs(min(V)), abs(max(V))]) * 1.3
    ax[i, 2].plot(xStan, y, label='$pile$', color='g', linewidth=2)
    ax[i, 2].plot(V, y, label='$ShearForce$', color='r')
    ax[i, 2].set_xlim(-VLim, VLim)
    ax[i, 2].set_ylim(-hWall, 0)
    ax[i, 2].grid(True)

    MLim = max([abs(min(M)), abs(max(M))]) * 1.3
    ax[i, 3].plot(xStan, y, label='$pile$', color='g', linewidth=2)
    ax[i, 3].plot(M, y, label='$ShearForce$', color='r')
    ax[i, 3].set_xlim(-MLim, MLim)
    ax[i, 3].set_ylim(-hWall, 0)
    ax[i, 3].grid(True)

    # print('-----------------')
    # print('i:',i)
    xLimList = [qLim, dispLim, VLim, MLim]
    for j in range(4):
        xLim = xLimList[j]
        # print('j:',j, 'xLim:',xLim)
        # ------- Excavation line --------
        if j != 1:
            ax[i, j].hlines(y=-hex, xmin=-xLim, xmax=0,
                            color='k', linewidth=1.0)
        else:
            ax[i, j].hlines(y=-hex, xmin=0, xmax=xLim,
                            color='k', linewidth=1.0)

        # ------- Water line --------
        yWaterBack = hWater if hex < hWater else hex
        if j != 1:
            ax[i, j].hlines(y=-hWater, xmin=0, xmax=xLim,
                            color='b', linestyle='--', linewidth=1.0)
            ax[i, j].hlines(y=-yWaterBack, xmin=-xLim, xmax=0,
                            color='b', linestyle='--', linewidth=1.0)
        else:
            ax[i, j].hlines(y=-hWater, xmin=-xLim, xmax=0,
                            color='b', linestyle='--', linewidth=1.0)
            ax[i, j].hlines(y=-hWater, xmin=0, xmax=xLim,
                            color='b', linestyle='--', linewidth=1.0)

        # ------- Sturt --------
        if isStrut:
            for hStrut in hStrutList:
                if j != 1:
                    ax[i, j].hlines(y=-hStrut, xmin=-xLim, xmax=0, color='g', linestyle='-', linewidth=2)
                else:
                    ax[i, j].hlines(y=-hStrut, xmin=0, xmax=xLim, color='g', linestyle='-', linewidth=2)

        # ------- Soil line --------
        for soilNum in range(len(input["soil"])):
            hSoil = input["soil"][soilNum]["h"]
            if j != 1:
                ax[i, j].hlines(y=-hSoil, xmin=0, xmax=xLim,
                                color='k', linewidth=1.0)
            else:
                ax[i, j].hlines(y=-hSoil, xmin=-xLim, xmax=0,
                                color='k', linewidth=1.0)

ax[0, 0].set_title('Pressure ($kN/m^2/m$)')
ax[0, 1].set_title('Displacement ($m/m$)')
ax[0, 2].set_title('ShearForce ($kN/m$)')
ax[0, 3].set_title('Moment ($kN \\bullet m/m$)')

plt.show()


# --------------------------------------------------------

# with open('layers2.json', 'rt', encoding='UTF8') as json_file:
#     layers = json.load(json_file)

# # Hem = 6.9
# modelObj = model.makeModelObj(layers, 15)
# with open("modelObjM.json", "w", encoding='UTF8') as json_file:
#     json.dump(modelObj, json_file)
# y = []
# for node in modelObj["node"]:
#     y.append(node["y"])
# xStan = np.zeros(len(y))
# dispExcel = np.ones(len(y)) * -1.75*1e-3
# matrixObj = solver4.solver(modelObj, properties)
# matrixObj2 = solver4.passive_active(modelObj, properties)


# disp = matrixObj["dispVec"][0]
# beforeDisp = matrixObj["beforeDispVec"][0]
# q = 1*matrixObj["qVec"][0]
# q_0 = 1*matrixObj2["qVec"][0]
# q_0_excav = 1*matrixObj2["qVec2"][0]
# q_active = 1*matrixObj2["qVec3"][0]
# q_active_excav = 1*matrixObj2["qVec4"][0]
# q_passive = 1*matrixObj2["qVec5"][0]
# q_passive_excav = 1*matrixObj2["qVec6"][0]

# matrixJson = {"dispVec": disp.tolist(), "qVec": q.tolist()}
# with open("matrixObjM.json", "w", encoding='UTF8') as json_file:
#     json.dump(matrixJson, json_file)


# # qLim = max(abs(max(q)), abs(min(q)))*1.3
# # dispLim = max(abs(max(disp)), abs(min(disp)))*1.3
# # dispLim = 0.00175*1.5 if dispLim <= 0.00175 else dispLim


# # fig, ax = plt.subplots(2)
# ax[1,0].plot(xStan, y, label='$pile$')
# # ax[0].plot(q_0, y, label='$p_0$')
# # ax[0].plot(q_active, y, label='$p_a$')
# # ax[0].plot(q_passive, y, label='$p_p$')
# ax[1,0].plot(q, y, label='$pressure$')
# # ax[0].plot(q_0_excav, y, label='$p_0excav$')
# # ax[0].plot(q_active_excav, y, label='$p_a excav$')
# # ax[0].plot(q_passive_excav, y, label='$p_p excav$')
# ax[1,0].set_xlabel('Pressure ($kN/m^2$)')
# ax[1,0].set_ylabel('Depth ($m$)')
# ax[1,0].legend()
# ax[1,0].set_xlim(-qLim, qLim)
# # ax[0].set_xlim(38, -38)
# ax[1,0].set_ylim(-8, 2)
# ax[1,0].grid(True)

# ax[1,1].plot(xStan, y, label='$Pile$')
# ax[1,1].plot(disp, y, label='$disp$')
# # ax[1].plot(beforeDisp, y, label='$disp_before$')
# # ax[1].plot(dispExcel, y, label='$excel$')
# ax[1,1].set_xlabel('Displacement ($m$)')
# ax[1,1].set_ylabel('Depth ($m$)')
# ax[1,1].legend()
# ax[1,1].set_xlim(dispLim, -dispLim)
# # ax[1].set_xlim(-0.0025, 0.0025)
# ax[1,1].set_ylim(-8, 2)
# ax[1,1].grid(True)

# plt.show()
