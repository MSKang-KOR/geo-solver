import json
import math
from matplotlib.pyplot import xlim
import numpy as np

with open('layers.json', 'rt', encoding='UTF8') as json_file:
    layers = json.load(json_file)


def makeModelObj(layers, q=0, gammaW=10):
    result = {"node": [], "pressure": {"backSide": {"active": [],
                                                    "passive": [], "rest": []}, "excavationSide": {"active": [], "passive": [], "rest": []}}, "pyCurve": []}

    qUpperSoil1 = 0
    qUpperSoil2 = 0
    nodeNum = 3
    dhInit = 0.01
    dhInit = 0.1
    qTestList = [[-8, -10.5], [-10.5, -13],
                 [-13, -1], [19, 1], [1, -2], [-2, 0]]
    # qTestList = [[0,-3.61],[7.6,2.0],[7.6,0],[0,-2],[-2,1.5],[1.5,-3.2]]
    for i in range(len(layers)):
        layer = layers[i]
        h1 = layer["h1"]/1e3
        h2 = layer["h2"]/1e3
        hw = 0
        Hlayer = abs(h2 - h1)   # unit: m
        mTest = (qTestList[i][1] - qTestList[i][0])/(h2 - h1)
        coeffTest = [mTest, qTestList[i][0] - mTest*h1]

        gammaSAT = layer["gammaSAT"]
        gamma = layer["gamma"] if not layer["isWater"] else gammaSAT - gammaW
        c = layer["c"]
        phi = layer["phi"] * math.pi/180
        isExcavation = layer["isExcavation"]

        K0 = 1 - math.sin(phi)
        Ka = (1 - math.sin(phi)) / (1 + math.sin(phi))
        Kp = 1 / Ka

        numEle = round(Hlayer/dhInit)
        dh = Hlayer/numEle
        Kh = layer["Kh"]    # unit: kN/m3
        # kh = Kh if not layer["isExcavation"] else 0  # unit: kN/m2
        kh = Kh*dh

        if i == 0:
            result["node"].append(
                {"number": 1, "y": h1 + dh*2, "kh": kh, "Kh": Kh, "dh": dh, "isExcavation": isExcavation, "qTest": 0})
            result["node"].append(
                {"number": 2, "y": h1 + dh*1, "kh": kh, "Kh": Kh, "dh": dh, "isExcavation": isExcavation, "qTest": 0})
            for j in range(2):
                result["pyCurve"].append(
                    {"coeff": [[0, 0]], "xLim": [0, 0], "pLim": [0, 0]})
                for loca in ["backSide", "excavationSide"]:
                    for type in ["rest", "active", "passive"]:
                        result["pressure"][loca][type].append(0)

        for num in range(numEle+1):
            if i != 0 and num == 0:
                result["node"].pop(-1)
                result["pyCurve"].pop(-1)
                for loca in ["backSide", "excavationSide"]:
                    for type in ["rest", "active", "passive"]:
                        result["pressure"][loca][type].pop(-1)

            y = h1 - dh * num
            h = abs(h1-y)
            qTest = coeffTest[0] * y + coeffTest[1]
            result["node"].append(
                {"number": nodeNum, "y": y, "kh": kh, "Kh": Kh, "dh": dh, "isExcavation": isExcavation, "qTest": qTest})

            # 배면측 정지, 주동, 수동 토압
            p01 = K0 * (q + qUpperSoil1 + gamma*h) - 2 * c * math.sqrt(K0)
            pa1 = Ka * (q + qUpperSoil1 + gamma*h) - 2 * c * \
                math.sqrt(Ka)+(gammaW * hw)
            pa1 = 0 if pa1 < 0 else pa1
            pp1 = Kp * (q + qUpperSoil1 + gamma*h) + 2 * c * math.sqrt(Kp)

            # 굴착측 정지, 주동, 수동 토압
            p02 = -(K0 * (qUpperSoil2 + gamma * h) - 2 * c *
                    math.sqrt(K0))if not layer["isExcavation"] else 0
            # p02 = 0
            # p02 = 0 if p02 < 0 else p02

            pa2 = -(Ka * (qUpperSoil2 + gamma*h) - 2 * c *
                    math.sqrt(Ka)+(gammaW * hw)) if not layer["isExcavation"] else 0

            pa2 = 0 if pa2 > 0 else pa2
            pp2 = -(Kp * (qUpperSoil2 + gamma*h) + 2 * c *
                    math.sqrt(Kp)) if not layer["isExcavation"] else 0

            result["pressure"]["backSide"]["rest"].append(p01)
            result["pressure"]["backSide"]["active"].append(pa1)
            result["pressure"]["backSide"]["passive"].append(pp1)
            result["pressure"]["excavationSide"]["rest"].append(p02)
            result["pressure"]["excavationSide"]["active"].append(pa2)
            result["pressure"]["excavationSide"]["passive"].append(pp2)

            if layer["isExcavation"]:
                xa = -(pa1-p01)/Kh
                xp = -(pp1-p01)/Kh
                coeff = [[-Kh, p01]]  # p-y curve 1차 방정식 계수
                xLim = [xp, xa]
                pLim = [pp1, pa1]
            else:
                pl = pp1 + pa2
                pr = pa1 + pp2
                xp1 = -(pp1-p01)/Kh
                xa1 = -(pa1-p01)/Kh
                xa2 = -(pa2-p02)/Kh
                xp2 = -(pp2-p02)/Kh
                if xa2 == 0:
                    coeff = [[-Kh, p01]]  # p-y curve 1차 방정식 계수
                    xLim = [xp1, xa1]
                    pLim = [pp1, pa1]
                else:
                    _points = [(xp1, pp1+(pa2 if -Kh*xp1+p02 > pa2 else -Kh*xp1+p02)),
                               (xa1, pa1+(pp2 if -Kh*xa1+p02 < pp2 else -Kh*xa1+p02)),
                               (xa2, pa2+(pp1 if -Kh*xa2+p01 > pp1 else -Kh*xa2+p01)),
                               (xp2, pp2+(pa1 if -Kh*xp2+p01 < pa1 else -Kh*xp2+p01))]
                    _points.sort(key=lambda x: x[0])
                    xLim = [item[0] for item in _points]
                    pLim = [item[1] for item in _points]
                    coeff = [
                        [(_points[1][1]-_points[0][1])/(_points[1][0]-_points[0][0]),
                         _points[0][1]-(_points[1][1]-_points[0][1])/(_points[1][0]-_points[0][0])*_points[0][0]],
                        [(_points[2][1]-_points[1][1])/(_points[2][0]-_points[1][0]),
                         _points[1][1]-(_points[2][1]-_points[1][1])/(_points[2][0]-_points[1][0])*_points[1][0]],
                        [(_points[3][1]-_points[2][1])/(_points[3][0]-_points[2][0]),
                         _points[2][1]-(_points[3][1]-_points[2][1])/(_points[3][0]-_points[2][0])*_points[2][0]]
                    ]

            result["pyCurve"].append(
                {"coeff": coeff, "xLim": xLim, "pLim": pLim})

            nodeNum += num
        qUpperSoil1 += gamma * Hlayer
        qUpperSoil2 += gamma * Hlayer if not layer["isExcavation"] else 0
        hw += Hlayer if layer["isWater"] else 0

        if i == len(layers)-1:
            result["node"].append(
                {"number": nodeNum, "y": h2 - dh*1, "kh": kh, "Kh": Kh, "dh": dh, "isExcavation": isExcavation, "qTest": qTest})
            result["node"].append(
                {"number": nodeNum+1, "y": h2 - dh*2, "kh": kh, "Kh": Kh, "dh": dh, "isExcavation": isExcavation, "qTest": qTest})
            for j in range(2):
                result["pyCurve"].append(
                    {"coeff": [[0, 0]], "xLim": [0, 0], "pLim": [0, 0]})
                result["pressure"]["backSide"]["rest"].append(0)
                result["pressure"]["backSide"]["active"].append(0)
                result["pressure"]["backSide"]["passive"].append(0)
                result["pressure"]["excavationSide"]["rest"].append(0)
                result["pressure"]["excavationSide"]["active"].append(0)
                result["pressure"]["excavationSide"]["passive"].append(0)

    return result
