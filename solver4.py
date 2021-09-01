
import numpy as np
import matplotlib.pyplot as plt
import json
import math
import copy


# def calcDispFromPYCurve(p, pyCurve):
#     coeff = pyCurve["coeff"]
#     pLim = pyCurve["pLim"]
#     if len(coeff) == 1:
#         if p >= pLim[0]:
#             p = pLim[0]
#         elif p <= pLim[1]:
#             p = pLim[1]
#         x = (p - coeff[0][1]) / coeff[0][0]
#     else:
#         if p >= pLim[0]:
#             index = 0
#             p = pLim[0]
#         elif pLim[1] <= p and p < pLim[0]:
#             index = 0
#         elif pLim[2] <= p and p < pLim[1]:
#             index = 1
#         elif pLim[3] <= p and p < pLim[2]:
#             index = 2
#         else:
#             index = 2
#             p = pLim[3]
#         x = (p - coeff[index][1]) / coeff[index][0]
#     return x


def calcForceFromPYCurve(x, pyCurve):
    coeff = pyCurve["coeff"]
    xLim = pyCurve["xLim"]
    pLim = pyCurve["pLim"]
    if len(coeff) == 1:
        if x < xLim[0]:
            p = pLim[0]
            Kh = 0
        elif x > xLim[1]:
            p = pLim[1]
            Kh = 0
        else:
            p = (coeff[0][0]*x + coeff[0][1])
            Kh = coeff[0][0]
        # print(coeff)
    else:
        if x < xLim[0]:
            p = pLim[0]
            Kh = 0
        elif xLim[0] <= x and x < xLim[1]:
            index = 0
            p = (coeff[index][0]*x + coeff[index][1])
            Kh = coeff[index][0]
        elif xLim[1] < x and x < xLim[2]:
            index = 1
            p = (coeff[index][0]*x + coeff[index][1])
            Kh = coeff[index][0]
        elif xLim[2] < x and x <= xLim[3]:
            index = 2
            p = (coeff[index][0]*x + coeff[index][1])
            Kh = coeff[index][0]
        elif x > xLim[3]:
            p = pLim[3]
            Kh = 0
        else:
            print("error")
            raise IndexError()

    return p, abs(Kh)

# 동일 할때는 kh가 0이 아님. limit을 넘어갈때만 0?


# def KhFromPYCurve(p, pyCurve, i):
#     coeff = pyCurve["coeff"]
#     xLim = pyCurve["xLim"]
#     pLim = pyCurve["pLim"]

#     if len(coeff) == 1:
#         Kh = 0 if p > pLim[0] or pLim[1] > p else coeff[0][0]
#         testNum = 0
#     else:
#         if p < pLim[3] or p > pLim[0]:
#             Kh = 0
#             testNum = 1
#         elif pLim[3] <= p and p < pLim[2]:
#             Kh = coeff[2][0]
#             testNum = 2
#         elif pLim[2] < p and p < pLim[1]:
#             Kh = coeff[1][0]
#             testNum = 3
#         elif pLim[1] < p and p <= pLim[0]:
#             Kh = coeff[0][0]
#             testNum = 4
#         else:
#             print("error")
#             raise IndexError()
#     return Kh


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
            alpha = 6 + (Kh) * (h**4)/EI
            b = np.array([1, -4, alpha, -4, 1])
            c = np.append(a, b)
            d = np.zeros(length - len(a) - len(b))
            Ki = np.append(c, d)

        kMatrix[i, :] = Ki * EI
    return kMatrix


def calcStiffnessMat(modelProperties, forceMatrix, beforeDispMatrix, checkMatrix, node_len, EI, qVec):
    errMax = 0.0001

    kMatrix = stiffnessMatrix(modelProperties["node"], EI)
    forceVec = forceMatrix.reshape(node_len, 1)
    kInverse = np.linalg.inv(kMatrix)

    dispVec = np.matmul(kInverse, forceVec)

    breakPoint = False

    for i in range(2, node_len-2):
        if not checkMatrix[i-2] or breakPoint:
            # pyCurve = modelProperties["pyCurve"][i]
            # P_prime, _ = calcForceFromPYCurve(dispVec[i], pyCurve)
            # err = abs(P_prime-qVec[i])
            err = abs(beforeDispMatrix[i]-dispVec[i])
            checkMatrix[i-2] = (err < errMax)
            breakPoint = True

    return dispVec, checkMatrix


def modProperties(modelProperties, numTotalNode, dispVec, coeffVec, checkVec, qVecIN, useP0=False):
    # qVecIN = np.zeros(numTotalNode)
    breakPoint = False
    for i in range(2, numTotalNode-2):
        if not checkVec[i-2] or breakPoint:
            pyCurve = modelProperties["pyCurve"][i]
            xBefore = dispVec[i][0]
            P_prime, Kh = calcForceFromPYCurve(xBefore, pyCurve)
            modelProperties["node"][i]["Kh"] = Kh
            qVecIN[i] = P_prime
            breakPoint = True
    forceVec = np.multiply(qVecIN,  coeffVec).reshape(numTotalNode, 1)
    return modelProperties, forceVec, qVecIN


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
        pForCheck, KhForCheck = calcForceFromPYCurve(xBefore, pyCurve)
        # KhForCheck = KhFromPYCurve(pForCheck, pyCurve, i)
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
        err = dispVecBefore[i][0] - dispVecNew[i][0]
        check = abs(err) <= errMax
        checkVec[checkNum] = check
        qVecNew[i], KhNew = calcForceFromPYCurve(dispVecBefore[i][0], pyCurve)

        nodeForLoopNew[i]["Kh"] = KhNew

    return {"checkVec": checkVec, "qVec": qVecNew, "nodeForLoop": nodeForLoopNew}


def solver(model_immutable, properties):
    numTotalNode = len(model_immutable["node"])
    numRealNode = len(model_immutable["node"]) - 4

    qVec = np.zeros(numTotalNode)
    model = copy.deepcopy(model_immutable)
    coeffVec = np.zeros(numTotalNode)
    forceVec = np.zeros(numTotalNode)
    dispVec = np.zeros(numTotalNode)
    KhList = np.zeros(numRealNode)
    checkVec = np.zeros(numRealNode)

    errMax = 0.001

    E = properties["pile"]["E"] * 1e3   # unit: kPa
    I = properties["pile"]["I"] * 1e-12  # unit: m^4
    EI = E*I

    # setup qVec, forceVec and coeffVec
    for i in range(numTotalNode):
        h = model["node"][i]["dh"]
        qVec[i] = model["pressure"]["backSide"]["rest"][i] + \
            model["pressure"]["excavationSide"]["rest"][i]
        # qVec[i] = model["node"][i]["qTest"]
        coeffVec[i] = (h**4)
        forceVec[i] = qVec[i] * coeffVec[i]

    beforeDispVec = copy.deepcopy(dispVec)
    forceVec = forceVec.reshape(numTotalNode, 1)
    dispVec, checkVec = calcStiffnessMat(
        model, forceVec, dispVec, checkVec, numTotalNode, EI, qVec)

    test = 0
    while False in checkVec and test !=25:
        # print("here")
        model, forceVec, qVec = modProperties(
            model, numTotalNode, dispVec, coeffVec, checkVec, copy.deepcopy(qVec))

        beforeDispVec = copy.deepcopy(dispVec)
        dispVec, checkVec = calcStiffnessMat(
            model, forceVec, dispVec, checkVec, numTotalNode, EI, qVec)

        test += 1
        # print(test)

    # print(test)
    print(checkVec)

    return {"qVec": qVec.reshape(1, numTotalNode), "dispVec": dispVec.reshape(1, numTotalNode), "beforeDispVec": beforeDispVec.reshape(1, numTotalNode)}


def passive_active(model, properties):
    numTotalNode = len(model["node"])

    qVec = np.zeros(numTotalNode)
    qVec2 = np.zeros(numTotalNode)
    qVec3 = np.zeros(numTotalNode)
    qVec4 = np.zeros(numTotalNode)
    qVec5 = np.zeros(numTotalNode)
    qVec6 = np.zeros(numTotalNode)

    for i in range(numTotalNode):
        qVec[i] = model["pressure"]["backSide"]["rest"][i] + \
            model["pressure"]["excavationSide"]["rest"][i]
        qVec2[i] = model["pressure"]["excavationSide"]["rest"][i]
        qVec3[i] = model["pressure"]["backSide"]["active"][i] + \
            model["pressure"]["excavationSide"]["passive"][i]
        qVec4[i] = model["pressure"]["excavationSide"]["active"][i]
        qVec5[i] = model["pressure"]["backSide"]["passive"][i] + \
            model["pressure"]["excavationSide"]["active"][i]
        qVec6[i] = model["pressure"]["excavationSide"]["passive"][i]

    return {
        "qVec": qVec.reshape(1, numTotalNode),
        "qVec2": qVec2.reshape(1, numTotalNode),
        "qVec3": qVec3.reshape(1, numTotalNode),
        "qVec4": qVec4.reshape(1, numTotalNode),
        "qVec5": qVec5.reshape(1, numTotalNode),
        "qVec6": qVec6.reshape(1, numTotalNode)
    }
