
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
            Kh = 100
        elif x > xLim[1]:
            p = pLim[1]
            Kh = 100
        else:
            p = (coeff[0][0]*x + coeff[0][1])
            Kh = coeff[0][0]
    else:
        if x < xLim[0]:
            p = pLim[0]
            Kh = 100
        elif xLim[0] < x and x < xLim[1]:
            index = 0
            p = (coeff[index][0]*x + coeff[index][1])
            Kh = coeff[index][0]
        elif xLim[1] < x and x < xLim[2]:
            index = 1
            p = (coeff[index][0]*x + coeff[index][1])
            Kh = coeff[index][0]
        elif xLim[2] < x and x < xLim[3]:
            index = 2
            p = (coeff[index][0]*x + coeff[index][1])
            Kh = coeff[index][0]
        elif x > xLim[3]:
            p = pLim[3]
            Kh = 100
        else:
            print("error")
            raise IndexError()

    return p, abs(Kh)

# 동일 할때는 kh가 0이 아님. limit을 넘어갈때만 0?


def KhFromPYCurve(p, pyCurve, i):
    coeff = pyCurve["coeff"]
    xLim = pyCurve["xLim"]
    pLim = pyCurve["pLim"]

    if len(coeff) == 1:
        Kh = 0 if p > pLim[0] or pLim[1] > p else coeff[0][0]
        testNum = 0
    else:
        if p < pLim[3] or p > pLim[0]:
            Kh = 100
            testNum = 1
        elif pLim[3] < p and p < pLim[2]:
            Kh = coeff[2][0]
            testNum = 2
        elif pLim[2] < p and p < pLim[1]:
            Kh = coeff[1][0]
            testNum = 3
        elif pLim[1] < p and p < pLim[0]:
            Kh = coeff[0][0]
            testNum = 4
        else:
            print("error")
            raise IndexError()
    return abs(Kh)


def stiffnessMatrix(nodeInfo, kList, EI):
    length = len(nodeInfo)
    kMatrix = np.zeros([length, length])
    for i in range(length):
        h = nodeInfo[i]["dh"]
        Kh = kList[i]
        Ks = nodeInfo[i]["Ks"]
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
            alpha = 6 + (Kh + Ks) * (h**4)/EI
            b = np.array([1, -4, alpha, -4, 1])
            c = np.append(a, b)
            d = np.zeros(length - len(a) - len(b))
            Ki = np.append(c, d)

        kMatrix[i, :] = Ki * EI
    return kMatrix


def calcStiffnessMat(modelProperties, kList, forceMatrix, beforeDispMatrix, checkMatrix, node_len, EI, qVec, beforeqVec):
    errMax = 0.00000001

    kMatrix = stiffnessMatrix(modelProperties["node"], kList, EI)
    forceVec = forceMatrix.reshape(node_len, 1)
    kInverse = np.linalg.inv(kMatrix)

    dispVec = np.matmul(kInverse, forceVec)

    for i in range(2, node_len-2):
        # pyCurve = modelProperties["pyCurve"][i]
        # P_prime, _ = calcForceFromPYCurve(dispVec[i], pyCurve)
        # err = abs(beforeqVec[i]-qVec[i])
        err = abs(beforeDispMatrix[i]-dispVec[i])
        checkMatrix[i-2] = (err < errMax)

    return dispVec, checkMatrix


def modProperties(modelProperties, kList, numTotalNode, dispVec, coeffVec, checkVec, qVecIN, p0Vec, beforeqVec, useP0=False):
    # qVecIN = np.zeros(numTotalNode)

    for i in range(2, numTotalNode-2):

        pyCurve = modelProperties["pyCurve"][i]
        xBefore = dispVec[i][0]
        _, Kh = calcForceFromPYCurve(xBefore, pyCurve)
        # Kh = KhFromPYCurve(beforeqVec[i], pyCurve, i)
        # if xBefore > 0:
        #     qVecIN[i] -= Kh*abs(xBefore)
        # else:
        #     qVecIN[i] += Kh*abs(xBefore)
        qVecIN[i] -= Kh*xBefore
        Kh = KhFromPYCurve(qVecIN[i], pyCurve, i)
        if len(pyCurve["pLim"]) == 2:
            if qVecIN[i] > pyCurve["pLim"][0]:
                qVecIN[i] = pyCurve["pLim"][0]
                # Kh = 0
            elif qVecIN[i] < pyCurve["pLim"][1]:
                qVecIN[i] = pyCurve["pLim"][1]
                # Kh = 0
        else:
            if qVecIN[i] > pyCurve["pLim"][0]:
                qVecIN[i] = pyCurve["pLim"][0]
                # Kh = 0
            elif qVecIN[i] < pyCurve["pLim"][3]:
                qVecIN[i] = pyCurve["pLim"][3]
                # Kh = 0
        # modelProperties["node"][i]["Kh"] = Kh
        kList[i] = Kh

    forceVec = np.multiply(qVecIN,  coeffVec).reshape(numTotalNode, 1)
    return modelProperties, kList, forceVec, qVecIN


def calcShearForce(model, dispVec, numTotalNode, EI):
    nodeInfo = model["node"]
    length = len(nodeInfo)
    vMatrix = np.zeros([length, length])
    for i in range(length):
        h = nodeInfo[i]["dh"]
        if not i in [0, 1, numTotalNode-2, numTotalNode-1]:
            a = np.zeros(i-2)
            b = [-1, 2, 0, -2, 1]
            c = np.append(a, b)
            d = np.zeros(length-len(c))
            Ki = np.append(c, d)
            vMatrix[i, :] = Ki
        else:
            vMatrix[i] = np.zeros(numTotalNode)
    vVec = np.matmul(vMatrix, dispVec.reshape(numTotalNode, 1)) * EI / (2*h**3)
    return vVec


def calcMoment(model, dispVec, numTotalNode, EI):
    nodeInfo = model["node"]
    length = len(nodeInfo)
    mMatrix = np.zeros([length, length])
    for i in range(length):
        h = nodeInfo[i]["dh"]
        if not i in [0, 1, numTotalNode-2, numTotalNode-1]:
            a = np.zeros(i-1)
            b = [1, -2, 1]
            c = np.append(a, b)
            d = np.zeros(length-len(c))
            Ki = np.append(c, d)
            mMatrix[i, :] = Ki
        else:
            mMatrix[i] = np.zeros(numTotalNode)
    mVec = np.matmul(mMatrix, dispVec.reshape(numTotalNode, 1)) * EI / (h**2)
    return mVec


def solver(input, model_immutable, qVec0=np.array([]), kList0=np.array([]), dispVec0=np.array([]), **kwargs):
    numTotalNode = len(model_immutable["node"])
    numRealNode = len(model_immutable["node"]) - 4

    qVec = np.zeros(numTotalNode)
    model = copy.deepcopy(model_immutable)
    coeffVec = np.zeros(numTotalNode)
    forceVec = np.zeros(numTotalNode)
    dispVec = np.zeros(numTotalNode)
    kList = np.zeros(numTotalNode)
    checkVec = np.zeros(numRealNode)

    E = input["wall"]["E"]   # unit: kPa
    I = input["wall"]["I"]   # unit: m^4
    EI = E*I
    beforeqVec = copy.deepcopy(qVec)
    # setup qVec, forceVec and coeffVec
    for i in range(numTotalNode):
        h = model["node"][i]["dh"]
        coeffVec[i] = (h**4)

        if len(qVec0) == 0:
            qVec[i] = model["pressure"]["backSide"]["rest"][i] + \
                model["pressure"]["excavationSide"]["rest"][i]
            kList[i] = model["node"][i]["Kh"]
        else:
            qVec[i] = qVec0[0][i]
            kList[i] = kList0[i]

        forceVec[i] = qVec[i] * coeffVec[i]

    isSolve = kwargs["solve"]
    if isSolve:
        beforeDispVec = copy.deepcopy(dispVec)
        forceVec = forceVec.reshape(numTotalNode, 1)
        dispVec, checkVec = calcStiffnessMat(
            model, kList, forceVec, dispVec, checkVec, numTotalNode, EI, qVec, beforeqVec)
        p0Vec = copy.deepcopy(qVec)
        numIter = 0
        while False in checkVec and numIter != 3000:
            beforeqVec = copy.deepcopy(qVec)
            model, kList, forceVec, qVec = modProperties(
                model, kList, numTotalNode, dispVec, coeffVec, checkVec, copy.deepcopy(qVec), p0Vec, beforeqVec)

            beforeDispVec = copy.deepcopy(dispVec)
            dispVec, checkVec = calcStiffnessMat(
                model, kList, forceVec, beforeDispVec, checkVec, numTotalNode, EI, qVec, beforeqVec)

            numIter += 1
        print("Convergence:", "Success" if not False in checkVec else "Fail",
              ", ", "Number of iteration:", numIter)
    else:
        dispVec = dispVec0

    vVec = calcShearForce(model, dispVec, numTotalNode, EI)
    mVec = calcMoment(model, dispVec, numTotalNode, EI)

    # return {"kList": kList, "qVec": qVec.reshape(1, numTotalNode), "dispVec": dispVec.reshape(1, numTotalNode), "beforeDispVec": beforeDispVec.reshape(1, numTotalNode)}
    return {"kList": kList, "qVec": qVec.reshape(1, numTotalNode), "dispVec": dispVec.reshape(1, numTotalNode), "vVec": vVec.reshape(1, numTotalNode), "mVec": mVec.reshape(1, numTotalNode)}


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
