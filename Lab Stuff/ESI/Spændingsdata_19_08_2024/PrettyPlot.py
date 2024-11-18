import numpy as np
import matplotlib.pyplot as plt


def makePrettyPlot(DataSet,sqrtN = False):
    # Arguments:
    # DataSet: 2D np.array with 1st column x axis, 2nd column y axis
    # sqrtN: include sqrt N uncertainty on y data? Defaults to False
    figure, axes = plt.subplots()
    xdata = DataSet[:, 0]
    ydata = DataSet[:, 1]
    axes.scatter(xdata, ydata)
    if sqrtN:
        axes.errorbar(xdata, ydata, np.sqrt(ydata), ls='')
    #plt.show()
    return figure, axes

