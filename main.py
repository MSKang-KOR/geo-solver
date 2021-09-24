import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import json
import math
import copy

import model
import solver


properties = {"pile": {"E": 210*1e3, "I": 226000000, "A": 19100}}
with open('input2.json', 'rt', encoding='UTF8') as json_file:
    input = json.load(json_file)

# ------------------------------------------------------------

modelList = model.modelGenerator(input)
modelObj = modelList[0]["model"]
matrixObj2 = solver.passive_active(modelList[2]["model"], properties)
qTest = 1*matrixObj2["qVec3"][0]

y = []
for node in modelObj["node"]:
    y.append(node["y"])
xStan = np.zeros(len(y))

result = []
hWall = input["wall"]["h"]
hWater = input["water"]["h"]

struts = []
hStrutList = []

fig, ax = plt.subplots(len(modelList), 4)
fig.set_size_inches(18, 9)

if not len(modelList) in [0, 1]:
    ax[0, 0].set_title('Pressure ($kN/m^2/m$)', fontsize=10)
    ax[0, 1].set_title('Displacement ($m/m$)', fontsize=10)
    ax[0, 2].set_title('ShearForce ($kN/m$)', fontsize=10)
    ax[0, 3].set_title('Moment ($kN \\bullet m/m$)', fontsize=10)
else:
    ax[0].set_title('Pressure ($kN/m^2/m$)', fontsize=10)
    ax[1].set_title('Displacement ($m/m$)', fontsize=10)
    ax[2].set_title('ShearForce ($kN/m$)', fontsize=10)
    ax[3].set_title('Moment ($kN \\bullet m/m$)', fontsize=10)

maxObjList = []
for i in range(len(modelList)):
    step = input["step"][i]
    hex = step["excavation"] if step["excavation"] != 0 else hex

    nowStruts = []
    if len(step["strut"]):
        for name in step["strut"]:
            nowStruts.append(name)

    for m in range(len(nowStruts)):
        hStrutList.append(input["strut"][nowStruts[m]]["h"])

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

    graphFig1 = ax[i, 0] if len(modelList) != 1 else ax[0]
    graphFig2 = ax[i, 1] if len(modelList) != 1 else ax[1]
    graphFig3 = ax[i, 2] if len(modelList) != 1 else ax[2]
    graphFig4 = ax[i, 3] if len(modelList) != 1 else ax[3]

    qMax = max([abs(min(q[2:len(q)-3])), abs(max(q[2:len(q)-3]))])
    indexMaxQTuple = np.where(q == qMax) if qMax in q else np.where(q == -qMax)
    indexMaxQ = indexMaxQTuple[0][0]
    xMaxQ = round(q[indexMaxQ], 3)
    yMaxQ = y[indexMaxQ]
    qLim = qMax * 2
    graphFig1.plot(xStan, y, label='$pile$', color='g', linewidth=2)
    graphFig1.plot(q, y, label='$pressure$', color='r')
    graphFig1.set_ylabel('Depth ($m$)')
    graphFig1.set_xlim(-qLim, qLim)
    graphFig1.set_ylim(-hWall, 0)
    graphFig1.annotate('Max: {0}'.format(xMaxQ), xy=(
        xMaxQ, yMaxQ), xytext=(xMaxQ - abs(xMaxQ)/2, yMaxQ), fontsize=7)
    graphFig1.grid(True)

    dispMax = max([abs(min(disp[2:len(disp)-3])),
                  abs(max(disp[2:len(disp)-3]))])
    indexMaxDispTuple = np.where(
        disp == dispMax) if dispMax in disp else np.where(disp == -dispMax)
    indexMaxDisp = indexMaxDispTuple[0][0]
    xMaxDisp = round(disp[indexMaxDisp], 5)
    yMaxDisp = y[indexMaxDisp]
    xTextDisp = xMaxDisp + abs(xMaxDisp)/2
    yTextDisp = yMaxDisp if np.round(yMaxDisp, 3) != .0 else yMaxDisp - 0.35
    dispLim = dispMax * 4
    graphFig2.plot(xStan, y, label='$pile$', color='g', linewidth=2)
    graphFig2.plot(disp, y, label='$Displacement$', color='r')
    graphFig2.set_xlim(dispLim, -dispLim)
    graphFig2.set_ylim(-hWall, 0)
    graphFig2.annotate('Max: {0}'.format(xMaxDisp), xy=(
        xMaxDisp, yMaxDisp), xytext=(xTextDisp, yTextDisp), fontsize=7)
    graphFig2.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    graphFig2.grid(True)

    vMax = max([abs(min(v[2:len(v)-3])), abs(max(v[2:len(v)-3]))])
    indexMaxVTuple = np.where(v == vMax) if vMax in v else np.where(v == -vMax)
    indexMaxV = indexMaxVTuple[0][0]
    xMaxV = round(v[indexMaxV], 3)
    yMaxV = y[indexMaxV]
    vLim = vMax * 2
    graphFig3.plot(xStan, y, label='$pile$', color='g', linewidth=2)
    graphFig3.plot(v, y, label='$ShearForce$', color='r')
    graphFig3.set_xlim(-vLim, vLim)
    graphFig3.set_ylim(-hWall, 0)
    graphFig3.annotate('Max: {0}'.format(xMaxV), xy=(
        xMaxV, yMaxV), xytext=(xMaxV - abs(xMaxV)/2, yMaxV), fontsize=7)
    graphFig3.grid(True)

    mMax = max([abs(min(m[2:len(m)-3])), abs(max(m[2:len(m)-3]))])
    indexMaxMTuple = np.where(m == mMax) if mMax in m else np.where(m == -mMax)
    indexMaxM = indexMaxMTuple[0][0]
    xMaxM = round(m[indexMaxM], 3)
    yMaxM = y[indexMaxM]
    mLim = mMax * 2
    mLim = max([abs(min(m)), abs(max(m))]) * 2
    graphFig4.plot(xStan, y, label='$pile$', color='g', linewidth=2)
    graphFig4.plot(m, y, label='$ShearForce$', color='r')
    graphFig4.set_xlim(-mLim, mLim)
    graphFig4.set_ylim(-hWall, 0)
    graphFig4.annotate('Max: {0}'.format(xMaxM), xy=(
        xMaxM, yMaxM), xytext=(xMaxM - abs(xMaxM)/2, yMaxM), fontsize=7)
    graphFig4.grid(True)

    xLimList = [qLim, dispLim, vLim, mLim]
    ySurfBack = np.ones(100)*(-hWall)

    for j in range(4):
        xLim = xLimList[j]
        xSurfBack = np.linspace(0, xLim, 100)
        xSurfEx = np.linspace(-xLim, 0, 100)
        ySurfEx = np.ones(100)*(-hex)
        figForDrawing = ax[i, j] if len(modelList) != 1 else ax[j]
        # ------- Soil surface --------
        if j != 1:
            figForDrawing.fill_between(xSurfBack, 0, ySurfBack,
                                       color='sandybrown', alpha=0.3)
            figForDrawing.fill_between(xSurfEx, ySurfEx, -hWall,
                                       color='sandybrown', alpha=0.3)
        else:
            figForDrawing.fill_between(xSurfEx, 0, ySurfBack,
                                       color='sandybrown', alpha=0.3)
            figForDrawing.fill_between(xSurfBack, ySurfEx, -
                                       hWall, color='sandybrown', alpha=0.3)

        # ------- Excavation line --------
        if j != 1:
            figForDrawing.hlines(y=-hex, xmin=-xLim, xmax=0,
                                 color='k', linewidth=1.0)
        else:
            figForDrawing.hlines(y=-hex, xmin=0, xmax=xLim,
                                 color='k', linewidth=1.0)

        # ------- Water line --------
        yWaterBack = hWater if hex < hWater else hex
        if j != 1:
            figForDrawing.hlines(y=-hWater, xmin=0, xmax=xLim,
                                 color='b', linestyle='--', linewidth=1.0)
            figForDrawing.hlines(y=-yWaterBack, xmin=-xLim, xmax=0,
                                 color='b', linestyle='--', linewidth=1.0)
        else:
            figForDrawing.hlines(y=-hWater, xmin=-xLim, xmax=0,
                                 color='b', linestyle='--', linewidth=1.0)
            figForDrawing.hlines(y=-yWaterBack, xmin=0, xmax=xLim,
                                 color='b', linestyle='--', linewidth=1.0)

        # ------- Sturt --------
        # if isStrut:
        for hStrut in hStrutList:
            if j != 1:
                figForDrawing.hlines(y=-hStrut, xmin=-xLim, xmax=0,
                                        color='g', linestyle='-', linewidth=2)
            else:
                figForDrawing.hlines(y=-hStrut, xmin=0, xmax=xLim,
                                        color='g', linestyle='-', linewidth=2)

        # ------- Soil line --------
        for soilNum in range(len(input["soil"])):
            hSoil = input["soil"][soilNum]["h"]
            if j != 1:
                figForDrawing.hlines(y=-hSoil, xmin=0, xmax=xLim,
                                     color='k', linewidth=1.0)
            else:
                figForDrawing.hlines(y=-hSoil, xmin=-xLim, xmax=0,
                                     color='k', linewidth=1.0)
    
    maxObj = {"q": xMaxQ, "u": xMaxDisp, "v": xMaxV, "m": xMaxM}
    maxObjList.append(maxObj)

plt.show()
