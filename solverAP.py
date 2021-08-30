
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


def renewMatrix(model, numTotalNode, errMax, nodeForLoop, dispVec, checkVecB, checkVecE, p0BVec, p0EVec, qVec):
    KhList = np.zeros(numTotalNode-4)
    p0BVecNew = copy.deepcopy(p0BVec)
    p0EVecNew = copy.deepcopy(p0EVec)
    qVecNew = copy.deepcopy(qVec)
    nodeForLoopNew = copy.deepcopy(nodeForLoop)
    checkVecBNew = copy.deepcopy(checkVecB)
    checkVecENew = copy.deepcopy(checkVecE)
    # print('here')
    # print(p0BVecNew)
    # print(p0EVecNew)
    # print(dispVec)
    for i in range(2, numTotalNode-2):
        checkNum = i-2
        isExcavation = model["node"][i]["isExcavation"]
        pyCurveB = model["pyCurve"][i]["back"]
        pyCurveE = model["pyCurve"][i]["excavation"]
        KhBefore = nodeForLoopNew[i]["Kh"]
        xBefore = dispVec[i][0]
        pForCheckB = calcForceFromPYCurve(xBefore, pyCurveB)
        pForCheckE = calcForceFromPYCurve(xBefore, pyCurveE)
        deffB = p0BVecNew[i] / pForCheckB
        deffE = p0EVecNew[i] / pForCheckE
        errB = 1 - abs(deffB)
        errE = 1 - abs(deffE)
        checkB = errB <= errMax
        checkE = errE <= errMax
        checkVecBNew[checkNum] = checkB
        checkVecENew[checkNum] = checkE

        pBNew = -p0BVecNew[i] + KhBefore*xBefore
        pENew = p0EVecNew[i] + KhBefore*xBefore
        if pBNew >= pyCurveB["pLim"][0]:
            pBNew = pyCurveB["pLim"][0]
        elif pBNew <= pyCurveB["pLim"][1]:
            pBNew = pyCurveB["pLim"][1]
        
        if pENew >= pyCurveE["pLim"][0]:
            pENew = pyCurveE["pLim"][0]
        elif pENew <= pyCurveE["pLim"][1]:
            pENew = pyCurveE["pLim"][1]

        p0BVecNew[i] = pBNew
        p0EVecNew[i] = pENew

        qVecNew[i] = -pBNew + pENew
        KhB = KhFromPYCurve(-pBNew, pyCurveB, i)
        KhE = KhFromPYCurve(-pENew, pyCurveE, i)
        
        KhNew = KhB + KhE
        KhList[checkNum] = KhNew
        nodeForLoopNew[i]["Kh"] = KhNew
    return {"p0BVec": p0BVecNew, "p0EVec": p0EVecNew, "qVec": qVecNew, "checkVecB": checkVecBNew, "checkVecE": checkVecENew, "nodeForLoop": nodeForLoopNew, "KhList": KhList}


def solver(model, properties):
    numTotalNode = len(model["node"])
    numRealNode = len(model["node"]) - 4

    p0BVec = np.zeros(numTotalNode)
    p0EVec = np.zeros(numTotalNode)
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
        p0BVec[i] = model["pressure"]["backSide"]["rest"][i]
        p0EVec[i] = model["pressure"]["excavationSide"]["rest"][i]
        qVec[i] = - p0BVec[i] + p0EVec[i]
        # qVec[i] = model["node"][i]["qTest"]
        coeffVec[i] = (h**4)
        forceVec[i] = qVec[i] * coeffVec[i]

    forceVec = forceVec.reshape(numTotalNode, 1)
    kInverse = np.linalg.inv(kMatrix)
    dispVec = np.matmul(kInverse, forceVec)

    checkVecB = []
    checkVecE = []
    for i in range(numRealNode):
        checkVecB.append(False)
        checkVecE.append(False)

    errMax = 0.01
    KhList = np.zeros(numRealNode)
    nodeForLoop = copy.deepcopy(model["node"])

    obj = renewMatrix(model, numTotalNode, errMax,
                      nodeForLoop, dispVec, checkVecB, checkVecE, p0BVec, p0EVec, qVec)
    checkVecB = obj["checkVecB"]
    checkVecE = obj["checkVecE"]
    KhList = obj["KhList"]

    test = 0
    # print(False in checkVecB)
    # print(False in checkVecE)
    # print((False in checkVecB) and (False in checkVecE) and test != 51)
    # while (False in checkVecB) and (False in checkVecE) and test != 51:
    while test != 51:
        nodeForLoop = obj["nodeForLoop"]
        qVec = obj["qVec"]
        forceVec = np.multiply(qVec,  coeffVec).reshape(numTotalNode, 1)

        kMatrix = stiffnessMatrix(nodeForLoop, EI)
        kInverse = np.linalg.inv(kMatrix)
        dispVec = np.matmul(kInverse, forceVec)

        obj = renewMatrix(model, numTotalNode, errMax,
                      nodeForLoop, dispVec, checkVecB, checkVecE, p0BVec, p0EVec, qVec)
        checkVecB = obj["checkVecB"]
        checkVecE = obj["checkVecE"]
        KhList = obj["KhList"]

        test += 1

    # print(test)
    print(checkVecB)
    print(checkVecE)
    # print(KhList)

    return {"qVec": qVec.reshape(1, numTotalNode), "dispVec": dispVec.reshape(1, numTotalNode)}

