import numpy as np
import matplotlib.pyplot as plt


def polynomialFitForDisp(Hex, step, measured, numerics, order):
    measured = np.array(measured)
    numeric = numerics[step-1]

    i = 0
    check = False
    while not check:
        if abs(numeric[i][1]) >= Hex:
            check =True
            index = i
            break
        i+=1

    dispCalc = numeric[index:,0]
    yCalc = numeric[index:,1]
    dispMeasured = measured[:,0] + dispCalc[0]
    yMeasured = measured[:,1]

    x = np.concatenate((dispCalc.reshape(1,len(dispCalc))[0], dispMeasured.reshape(1,len(dispMeasured))[0]),axis=0)
    y = np.concatenate((yCalc.reshape(1,len(yCalc))[0], yMeasured.reshape(1,len(yMeasured))[0]),axis=0)
    # print(np.shape(x), np.shape(y))


    yp = np.linspace(-6.9,0,100)
    xp = np.poly1d(np.polyfit(y,x,4))

    fig2, ax2 = plt.subplots(1,1)
    ax2.plot(x,y,'.', xp(yp),yp,'-')
    ax2.grid(True)
    plt.show()

    return 'ok'

