import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import json
import math
import copy

import model
import solver


properties = {"pile": {"E": 210*1e3, "I": 226000000, "A": 19100}}
with open('input.json', 'rt', encoding='UTF8') as json_file:
    input = json.load(json_file)

# ------------------------------------------------------------

modelList = model.modelGenerator(input)
modelObj = modelList[0]["model"]
matrixObj2 = solver.passive_active(modelObj, properties)
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
fig.set_size_inches(18, 9)
ax[0, 0].set_title('Pressure ($kN/m^2/m$)', fontsize=10)
ax[0, 1].set_title('Displacement ($m/m$)', fontsize=10)
ax[0, 2].set_title('ShearForce ($kN/m$)', fontsize=10)
ax[0, 3].set_title('Moment ($kN \\bullet m/m$)', fontsize=10)
for i in range(len(modelList)):
    step = input["step"][i]
    hex = step["excavation"] if step["excavation"] != 0 else hex
    isStrut = True if len(step["strut"]) != 0 else isStrut

    modelObj = modelList[i]["model"]
    isSolve = modelList[i]["solve"]

    if i == 0:
        matrixObj = solver.solver(input, modelObj, solve=isSolve)
    else:
        qVec0 = copy.deepcopy(matrixObj["qVec"])
        kList0 = copy.deepcopy(matrixObj["kList"])
        dispVec0 = copy.deepcopy(matrixObj["dispVec"])
        matrixObj = solver.solver(
            input, modelObj, qVec0, kList0, dispVec0, solve=isSolve)
    result.append(matrixObj)

    step = result[i]
    q = step["qVec"][0]
    disp = step["dispVec"][0]
    v = step["vVec"][0]
    m = step["mVec"][0]

    qMax = max([abs(min(q[2:len(q)-3])), abs(max(q[2:len(q)-3]))])
    indexMaxQTuple = np.where(q==qMax) if qMax in q else np.where(q==-qMax)
    indexMaxQ = indexMaxQTuple[0][0]
    xMaxQ = round(q[indexMaxQ],3)
    yMaxQ = y[indexMaxQ]
    qLim = qMax * 2
    ax[i, 0].plot(xStan, y, label='$pile$', color='g', linewidth=2)
    ax[i, 0].plot(q, y, label='$pressure$', color='r')
    ax[i, 0].set_ylabel('Depth ($m$)')
    ax[i, 0].set_xlim(-qLim, qLim)
    ax[i, 0].set_ylim(-hWall, 0)
    ax[i, 0].annotate('Max: {0}'.format(xMaxQ), xy=(xMaxQ, yMaxQ), xytext=(xMaxQ - abs(xMaxQ)/2, yMaxQ), fontsize=7)
    ax[i, 0].grid(True)

    dispMax = max([abs(min(disp[2:len(disp)-3])), abs(max(disp[2:len(disp)-3]))])
    indexMaxDispTuple = np.where(disp==dispMax) if dispMax in disp else np.where(disp==-dispMax)
    indexMaxDisp = indexMaxDispTuple[0][0]
    xMaxDisp = round(disp[indexMaxDisp],5)
    yMaxDisp = y[indexMaxDisp]
    xTextDisp = xMaxDisp + abs(xMaxDisp)/2
    yTextDisp = yMaxDisp if np.round(yMaxDisp,3)!=.0 else yMaxDisp - 0.35
    dispLim = dispMax * 4
    ax[i, 1].plot(xStan, y, label='$pile$', color='g', linewidth=2)
    ax[i, 1].plot(disp, y, label='$Displacement$', color='r')
    ax[i, 1].set_xlim(dispLim, -dispLim)
    ax[i, 1].set_ylim(-hWall, 0)
    ax[i, 1].annotate('Max: {0}'.format(xMaxDisp), xy=(xMaxDisp, yMaxDisp), xytext=(xTextDisp, yTextDisp), fontsize=7)
    ax[i, 1].ticklabel_format(axis="x", style="sci", scilimits=(0,0))
    ax[i, 1].grid(True)

    vMax = max([abs(min(v[2:len(v)-3])), abs(max(v[2:len(v)-3]))])
    indexMaxVTuple = np.where(v==vMax) if vMax in v else np.where(v==-vMax)
    indexMaxV = indexMaxVTuple[0][0]
    xMaxV = round(v[indexMaxV],3)
    yMaxV = y[indexMaxV]
    vLim = vMax * 2
    ax[i, 2].plot(xStan, y, label='$pile$', color='g', linewidth=2)
    ax[i, 2].plot(v, y, label='$ShearForce$', color='r')
    ax[i, 2].set_xlim(-vLim, vLim)
    ax[i, 2].set_ylim(-hWall, 0)
    ax[i, 2].annotate('Max: {0}'.format(xMaxV), xy=(xMaxV, yMaxV), xytext=(xMaxV - abs(xMaxV)/2, yMaxV), fontsize=7)
    ax[i, 2].grid(True)

    mMax = max([abs(min(m[2:len(m)-3])), abs(max(m[2:len(m)-3]))])
    indexMaxMTuple = np.where(m==mMax) if mMax in m else np.where(m==-mMax)
    indexMaxM = indexMaxMTuple[0][0]
    xMaxM = round(m[indexMaxM],3)
    yMaxM = y[indexMaxM]
    mLim = mMax * 2
    mLim = max([abs(min(m)), abs(max(m))]) * 2
    ax[i, 3].plot(xStan, y, label='$pile$', color='g', linewidth=2)
    ax[i, 3].plot(m, y, label='$ShearForce$', color='r')
    ax[i, 3].set_xlim(-mLim, mLim)
    ax[i, 3].set_ylim(-hWall, 0)
    ax[i, 3].annotate('Max: {0}'.format(xMaxM), xy=(xMaxM, yMaxM), xytext=(xMaxM - abs(xMaxM)/2, yMaxM), fontsize=7)
    ax[i, 3].grid(True)

    # print('-----------------')
    # print('i:',i)
    xLimList = [qLim, dispLim, vLim, mLim]
    # ySurfBack = np.linspace(-hWall,0,100)
    ySurfBack = np.ones(100)*(-hWall)

    for j in range(4):
        xLim = xLimList[j]
        xSurfBack = np.linspace(0, xLim, 100)
        xSurfEx = np.linspace(-xLim, 0, 100)
        ySurfEx = np.ones(100)*(-hex)
        # ------- Soil surface --------
        if j != 1:
            ax[i, j].fill_between(xSurfBack, 0, ySurfBack,
                                  color='sandybrown', alpha=0.3)
            ax[i, j].fill_between(xSurfEx, ySurfEx, -hWall,
                                  color='sandybrown', alpha=0.3)
        else:
            ax[i, j].fill_between(xSurfEx, 0, ySurfBack,
                                  color='sandybrown', alpha=0.3)
            ax[i, j].fill_between(xSurfBack, ySurfEx, -
                                  hWall, color='sandybrown', alpha=0.3)

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
            ax[i, j].hlines(y=-yWaterBack, xmin=0, xmax=xLim,
                            color='b', linestyle='--', linewidth=1.0)

        # ------- Sturt --------
        if isStrut:
            for hStrut in hStrutList:
                if j != 1:
                    ax[i, j].hlines(y=-hStrut, xmin=-xLim, xmax=0,
                                    color='g', linestyle='-', linewidth=2)
                else:
                    ax[i, j].hlines(y=-hStrut, xmin=0, xmax=xLim,
                                    color='g', linestyle='-', linewidth=2)

        # ------- Soil line --------
        for soilNum in range(len(input["soil"])):
            hSoil = input["soil"][soilNum]["h"]
            if j != 1:
                ax[i, j].hlines(y=-hSoil, xmin=0, xmax=xLim,
                                color='k', linewidth=1.0)
            else:
                ax[i, j].hlines(y=-hSoil, xmin=-xLim, xmax=0,
                                color='k', linewidth=1.0)

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
