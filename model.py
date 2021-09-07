import json
import math
import copy
from matplotlib.pyplot import xlim
import numpy as np


def modelGenerator(input):
    result = []
    struts = []
    stepIndex = 0
    hex = 0
    for i in range(len(input["step"])):
        step = input["step"][i]

        if len(step["strut"]):
            for name in step["strut"]:
                struts.append(name)

        # if step["excavation"]!=0:
        stepIndex = i
        isSolve = True if step["excavation"] != 0 else False
        hex = step["excavation"] if step["excavation"] != 0 else hex
        layers = copy.deepcopy(input["layers"])
        model = makeModelObj(input, layers, struts, hex, stepIndex)
        result.append({"solve": isSolve, "model": model})
        print("Step:", i+1, ", ", "isSolve:", isSolve)
        # else:
        #     el = copy.deepcopy(result[stepIndex])
        #     el["solve"] = False

        #     result.append(el)

    return result


def makeModelObj(input, layers, struts, hex, stepIndex):
    result = {"node": [], "pressure": {"backSide": {"active": [],
                                                    "passive": [], "rest": []}, "excavationSide": {"active": [], "passive": [], "rest": []}}, "pyCurve": []}
    qUpperSoil1 = 0
    qUpperSoil2 = 0
    pw1 = 0
    pw2 = 0
    nodeNum = 3
    dhInit = 0.1
    hw = input["water"]["h"]       # 수위
    rhoW = input["water"]["rho"]    # 물의 단위중량
    # hex = input["step"][stepIndex]["excavation"]  # 굴착깊이
    q = input["step"][0]["surcharge"]   # 상재하중
    hStrutList = [input["strut"][name]["h"] for name in struts]

    for i in range(len(layers)):
        layer = layers[i]
        h1 = layer["h1"]
        h2 = layer["h2"]
        Hlayer = abs(h2 - h1)   # unit: m

        isWater = True if hw <= abs(h1) else False
        # print("hw:", hw, "h1:", abs(h1), "h2:", abs(h2), "isWater:", isWater)
        isExcavation = True if hex >= abs(h2) else False
        # print("hex:", hex, "h2:", abs(h2), "isExcavation:", isExcavation)

        rhoSAT = layer["rhoSAT"]
        rho = layer["rho"] if not layer["isWater"] else rhoSAT - rhoW
        c = layer["c"]
        phi = layer["phi"] * math.pi/180

        K0 = 1 - math.sin(phi)
        Ka = (1 - math.sin(phi)) / (1 + math.sin(phi))
        Kp = 1 / Ka

        numEle = round(Hlayer/dhInit)
        dh = Hlayer/numEle
        Kh = layer["Kh"]    # unit: kN/m3
        # kh = Kh if not isExcavation else 0  # unit: kN/m2
        # kh = Kh*dh

        if i == 0:
            result["node"].append(
                {"number": 1, "y": h1 + dh*2, "Kh": Kh, "dh": dh, "isExcavation": isExcavation, "Ks": 0})
            result["node"].append(
                {"number": 2, "y": h1 + dh*1, "Kh": Kh, "dh": dh, "isExcavation": isExcavation, "Ks": 0})
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

            isStrut = True if num == 0 and np.round(
                abs(h1), 3) in np.round(hStrutList, 3) else False

            strutIndex = 0
            Ks = 0
            for m in range(len(hStrutList)):
                if num == 0 and np.round(abs(h1), 3) in np.round(hStrutList, 3):
                    strutIndex = m
                    name = struts[strutIndex]
                    Es = input["strut"][name]["E"]
                    As = input["strut"][name]["A"]
                    Ls = input["strut"][name]["L"]
                    Ks = Es*As/Ls
                    break

            # if isStrut:
            #     print(isStrut, Ks)

            result["node"].append(
                {"number": nodeNum, "y": y, "Kh": Kh, "dh": dh, "isExcavation": isExcavation, "Ks": Ks})

            # 수압
            if isWater:
                if num != numEle:
                    pw1 += rhoW * dh
                if abs(hw) >= abs(h1):
                    pw2 += rhoW * dh
            # 배면측 정지, 주동, 수동 토압
            p01 = K0 * (q + qUpperSoil1 + rho*h) - \
                2 * c * math.sqrt(K0) + pw1
            pa1 = Ka * (q + qUpperSoil1 + rho*h) - \
                2 * c * math.sqrt(Ka) + pw1
            pa1 = 0 if pa1 < 0 else pa1
            pp1 = Kp * (q + qUpperSoil1 + rho*h) + \
                2 * c * math.sqrt(Kp) + pw1

            # 굴착측 정지, 주동, 수동 토압
            p02 = -(K0 * (qUpperSoil2 + rho * h) - 2 * c *
                    math.sqrt(K0) + pw2)if not isExcavation else 0
            # p02 = 0
            # p02 = 0 if p02 < 0 else p02
            pa2 = -(Ka * (qUpperSoil2 + rho*h) - 2 * c *
                    math.sqrt(Ka) + pw2) if not isExcavation else 0
            pa2 = 0 if pa2 > 0 else pa2
            pp2 = -(Kp * (qUpperSoil2 + rho*h) + 2 * c *
                    math.sqrt(Kp) + pw2) if not isExcavation else 0

            result["pressure"]["backSide"]["rest"].append(p01)
            result["pressure"]["backSide"]["active"].append(pa1)
            result["pressure"]["backSide"]["passive"].append(pp1)
            result["pressure"]["excavationSide"]["rest"].append(p02)
            result["pressure"]["excavationSide"]["active"].append(pa2)
            result["pressure"]["excavationSide"]["passive"].append(pp2)

            if isExcavation:
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
        qUpperSoil1 += rho * Hlayer
        qUpperSoil2 += rho * Hlayer if not isExcavation else 0

        if i == len(layers)-1:
            result["node"].append(
                {"number": nodeNum, "y": h2 - dh*1, "Kh": Kh, "dh": dh, "isExcavation": isExcavation, "Ks": 0})
            result["node"].append(
                {"number": nodeNum+1, "y": h2 - dh*2, "Kh": Kh, "dh": dh, "isExcavation": isExcavation, "Ks": 0})
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
