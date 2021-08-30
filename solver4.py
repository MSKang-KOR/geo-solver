
import numpy as np
import matplotlib.pyplot as plt
import json
import math
import copy


def calcDispFromPYCurve(p, pyCurve):
    coeff = pyCurve["coeff"]
    pLim = pyCurve["pLim"]
    if len(coeff) == 1:
        if p >= pLim[0]:
            p = pLim[0]
        elif p <= pLim[1]:
            p = pLim[1]
        x = (p - coeff[0][1]) / coeff[0][0]
    else:
        if p >= pLim[0]:
            index = 0
            p = pLim[0]
        elif pLim[1] <= p and p < pLim[0]:
            index = 0
        elif pLim[2] <= p and p < pLim[1]:
            index = 1
        elif pLim[3] <= p and p < pLim[2]:
            index = 2
        else:
            index = 2
            p = pLim[3]
        x = (p - coeff[index][1]) / coeff[index][0]
    return x


def calcForceFromPYCurve(x, pyCurve):
    coeff = pyCurve["coeff"]
    xLim = pyCurve["xLim"]
    if len(coeff) == 1:
        if x <= xLim[0]:
            x = xLim[0]
        elif xLim[1] <= x:
            x = xLim[1]
        p = (coeff[0][0]*x + coeff[0][1])
    else:
        if x <= xLim[0]:
            index = 0
            x = xLim[0]
        elif xLim[0] < x and x <= xLim[1]:
            index = 0
        elif xLim[1] < x and x <= xLim[2]:
            index = 1
        elif xLim[2] < x and x <= xLim[3]:
            index = 2
        elif xLim[3] <= x:
            index = 2
            x = xLim[3]
        p = (coeff[index][0]*x + coeff[index][1])
    return p


def KhFromPYCurve(p, pyCurve, i):
    coeff = pyCurve["coeff"]
    xLim = pyCurve["xLim"]
    pLim = pyCurve["pLim"]

    if len(coeff) == 1:
        Kh = 0 if p >= pLim[0] or pLim[1] >= p else -coeff[0][0]
        testNum = 0
    else:
        if p <= pLim[3]:
            Kh = 0
            testNum = 1
        elif pLim[3] < p and p <= pLim[2]:
            Kh = -coeff[2][0]
            testNum = 2
        elif pLim[2] < p and p <= pLim[1]:
            Kh = -coeff[1][0]
            testNum = 3
        elif pLim[1] < p and p <= pLim[0]:
            Kh = -coeff[0][0]
            testNum = 4
        else:
            Kh = 0
            testNum = 5
    return Kh


def stiffnessMatrix(nodeInfo, EI):
    length = len(nodeInfo)
    kMatrix = np.zeros([length, length])
    for i in range(length):
        h = nodeInfo[i]["dh"]
        Kh = nodeInfo[i]["Kh"]
        if i == 0 or i == 1:
            a = [-1, 2, 0, -2, 1] if i == 0 else [0, 1, -2, 1, 0]
            b = np.zeros(length-len(a))
            Ki = np.append(a, b)
            q = 0
        elif i == length-1 or i == length-2:
            a = [-1, 2, 0, -2, 1] if i == length-2 else [0, 1, -2, 1, 0]
            b = np.zeros(length-len(a))
            Ki = np.append(b, a)
            q = 0
        else:
            a = np.zeros(i-2)
            alpha = 6 + Kh * (h**4)/EI
            b = np.array([1, -4, alpha, -4, 1])
            c = np.append(a, b)
            d = np.zeros(length - len(a) - len(b))
            Ki = np.append(c, d)

        kMatrix[i, :] = Ki * EI
    return kMatrix


def renewMatrix(model, EI, numTotalNode, errMax, nodeForLoopBefore, dispVecBefore, qVecBefore, coeffVec):
    # KhList = np.zeros(numTotalNode-4)
    qVecIN = copy.deepcopy(qVecBefore)
    qVecNew = copy.deepcopy(qVecBefore)
    nodeForLoopNew = copy.deepcopy(nodeForLoopBefore)
    checkVec = np.zeros(numTotalNode-4)

    for i in range(2, numTotalNode-2):
        checkNum = i-2
        isExcavation = model["node"][i]["isExcavation"]
        pyCurve = model["pyCurve"][i]
        xBefore = dispVecBefore[i][0]
        pForCheck = calcForceFromPYCurve(xBefore, pyCurve)
        KhForCheck = KhFromPYCurve(pForCheck, pyCurve, i)
        nodeForLoopNew[i]["Kh"] = KhForCheck
        qVecIN[i] = pForCheck

    kMatrix = stiffnessMatrix(nodeForLoopNew, EI)
    kInverse = np.linalg.inv(kMatrix)
    forceVec = np.multiply(qVecIN, coeffVec).reshape(numTotalNode, 1)
    dispVecNew = np.matmul(kInverse, forceVec)

    for i in range(2, numTotalNode-2):
        checkNum = i-2
        isExcavation = model["node"][i]["isExcavation"]
        pyCurve = model["pyCurve"][i]
        err =  dispVecBefore[i][0] - dispVecNew[i][0]
        check = abs(err) <= errMax
        checkVec[checkNum] = check
        if i==8:
            print("-------------------")
            print("xBefore: ", dispVecBefore[i][0], " ", "xAfter: ", dispVecNew[i][0])
            print("err: ", err)
            print("chekc: ", check)

        # if isExcavation:
        #     if nodeForLoopBefore[i]["Kh"] == 0:
        #         pNew = qVecNew[i]
        #     else:
        #         qVecNew[i] = calcForceFromPYCurve(xBefore, pyCurve)
        # else:
        #     if nodeForLoopBefore[i]["Kh"] != 0:
        #         if dispVecBefore[i][0]*qVecBefore[i] < 0:
        #             qVecNew[i] = qVecBefore[i] - nodeForLoopBefore[i]["Kh"]*dispVecBefore[i][0]
        #         else:
        #             qVecNew[i] = qVecBefore[i] + nodeForLoopBefore[i]["Kh"]*dispVecBefore[i][0]
        #     else:
        #         qVecNew[i] = calcForceFromPYCurve(xBefore, pyCurve)

        # if not check:
        qVecNew[i] = calcForceFromPYCurve(dispVecBefore[i][0], pyCurve)
        # qVecNew[i] = qVecBefore[i] - nodeForLoopBefore[i]["Kh"]*dispVecBefore[i][0]
        if qVecNew[i] >= pyCurve["pLim"][0]:
            qVecNew[i] = pyCurve["pLim"][0]
        elif qVecNew[i] <= pyCurve["pLim"][1]:
            qVecNew[i] = pyCurve["pLim"][1]
        KhNew = KhFromPYCurve(qVecNew[i], pyCurve, i)
        nodeForLoopNew[i]["Kh"] = KhNew
        # # KhList[i] = KhNew
    return {"checkVec": checkVec, "qVec": qVecNew, "nodeForLoop": nodeForLoopNew}


def solver(model, properties):
    numTotalNode = len(model["node"])
    numRealNode = len(model["node"]) - 4

    qVec = np.zeros(numTotalNode)
    nodeForLoop = copy.deepcopy(model["node"])
    coeffVec = np.zeros(numTotalNode)
    forceVec = np.zeros(numTotalNode)
    dispVec = np.zeros(numTotalNode)
    KhList = np.zeros(numRealNode)

    E = properties["pile"]["E"] * 1e3   # unit: kPa
    I = properties["pile"]["I"] * 1e-12  # unit: m^4
    EI = E*I

    kMatrix = stiffnessMatrix(model["node"], EI)

    for i in range(numTotalNode):
        h = model["node"][i]["dh"]
        qVec[i] = - model["pressure"]["backSide"]["rest"][i] + \
            model["pressure"]["excavationSide"]["rest"][i]
        # qVec[i] = model["node"][i]["qTest"]
        coeffVec[i] = (h**4)
        forceVec[i] = qVec[i] * coeffVec[i]

    forceVec = forceVec.reshape(numTotalNode, 1)
    kInverse = np.linalg.inv(kMatrix)
    dispVec = np.matmul(kInverse, forceVec)

    errMax = 0.0001
    obj = renewMatrix(model, EI, numTotalNode, errMax,
                      nodeForLoop, dispVec, qVec, coeffVec)
    checkVec = obj["checkVec"]

    test = 0
    while False in checkVec and test != 100:
        nodeForLoop = obj["nodeForLoop"]
        qVec = obj["qVec"]
        forceVec = np.multiply(qVec,  coeffVec).reshape(numTotalNode, 1)

        kMatrix = stiffnessMatrix(nodeForLoop, EI)
        kInverse = np.linalg.inv(kMatrix)
        dispVec = np.matmul(kInverse, forceVec)

        obj = renewMatrix(model, EI, numTotalNode, errMax,
                      nodeForLoop, dispVec, qVec, coeffVec)
        checkVec = obj["checkVec"]
        test += 1
        print(test)

    # print(test)
    print(checkVec)

    return {"qVec": qVec.reshape(1, numTotalNode), "dispVec": dispVec.reshape(1, numTotalNode)}
