
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


def renewMatrix(model, numTotalNode, errMax, nodeForLoop, dispVec, checkVec, qVec):
    KhList = np.zeros(numTotalNode-4)
    qVecNew = copy.deepcopy(qVec)
    nodeForLoopNew = copy.deepcopy(nodeForLoop)
    checkVecNew = copy.deepcopy(checkVec)
    # print(dispVec)
    for i in range(2, numTotalNode-2):
        checkNum = i-2
        isExcavation = model["node"][i]["isExcavation"]
        pyCurve = model["pyCurve"][i]
        KhBefore = nodeForLoopNew[i]["Kh"]
        xBefore = dispVec[i][0]
        pForCheck = calcForceFromPYCurve(xBefore, pyCurve)
        deff = qVecNew[i] / pForCheck
        err = 1 - abs(deff)
        check = err <= errMax
        checkVecNew[checkNum] = check

        # if check:
        if isExcavation:
            # if KhBefore == 0:
            #     pNew = qVecNew[i]
            # else:
            pNew = calcForceFromPYCurve(xBefore, pyCurve)
        else:
            if KhBefore != 0:
                if xBefore*qVecNew[i] < 0:
                    pNew = qVecNew[i] + KhBefore*xBefore
                else:
                    pNew = qVecNew[i] - KhBefore*xBefore
            else:
                pNew = calcForceFromPYCurve(xBefore, pyCurve)

        qVecNew[i] = pNew
        KhNew = KhFromPYCurve(pNew, pyCurve, i)
        KhList[checkNum] = KhNew
        nodeForLoopNew[i]["Kh"] = KhNew
    return {"qVec": qVecNew, "checkVec": checkVecNew, "nodeForLoop": nodeForLoopNew, "KhList": KhList}


def solver(model, properties):
    numTotalNode = len(model["node"])
    numRealNode = len(model["node"]) - 4

    qVec = np.zeros(numTotalNode)
    coeffVec = np.zeros(numTotalNode)
    forceVec = np.zeros(numTotalNode)
    dispVec = np.zeros(numTotalNode)

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

    checkVec = []
    for i in range(numRealNode):
        checkVec.append(False)

    errMax = 0.1
    KhList = np.zeros(numRealNode)
    nodeForLoop = copy.deepcopy(model["node"])

    obj = renewMatrix(model, numTotalNode, errMax,
                      nodeForLoop, dispVec, checkVec, qVec)
    checkVec = obj["checkVec"]
    KhList = obj["KhList"]
    
    
    test = 0
    while False in checkVec and test != 0:
        nodeForLoop = obj["nodeForLoop"]
        qVec = obj["qVec"]
        forceVec = np.multiply(qVec,  coeffVec).reshape(numTotalNode, 1)
        
        kMatrix = stiffnessMatrix(nodeForLoop, EI)
        kInverse = np.linalg.inv(kMatrix)
        dispVec = np.matmul(kInverse, forceVec)

        obj = renewMatrix(model, numTotalNode, errMax,
                      nodeForLoop, dispVec, checkVec, qVec)
        checkVec = obj["checkVec"]
        KhList = obj["KhList"]
        
        test += 1

    print(test)
    print(checkVec[0:len(checkVec)-1])
    # print(KhList)

    return {"qVec": qVec.reshape(1, numTotalNode), "dispVec": dispVec.reshape(1, numTotalNode)}
