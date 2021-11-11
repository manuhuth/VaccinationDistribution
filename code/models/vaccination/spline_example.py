import pickle
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

from scipy.interpolate import CubicHermiteSpline as cbs


def finite_differences(xx, yy):
    dd = []

    fd = onesidedFD(yy[0], yy[1], xx[1] - xx[0])
    dd.append(fd)

    for i in range(1, len(xx) - 1):
        dd.append(
            centeredFD(
                yy[i - 1], yy[i], yy[i + 1], xx[i] - xx[i - 1], xx[i + 1] - xx[i]
            )
        )

    fd = onesidedFD(yy[-2], yy[-1], xx[-1] - xx[-2])
    dd.append(fd)

    return np.asarray(dd)


def onesidedFD(y0, y1, h):

    return (y1 - y0) * 1 / h


def centeredFD(ym1, y0, yp1, hm, hp):
    if hm == hp:
        return (yp1 - ym1) / (2 * hm)
    else:
        return ((yp1 - y0) / hp + (y0 - ym1) / hm) / 2


def get_spline(
    array,
    periods,
    length,
    total_length,
    grid_points=6000,
    transform=True,
    x=None,
    transform_y=True,
    grid=None,
):

    if transform_y is True:
        y = np.log(array / (1 - array))

    if x is None:
        x = np.linspace(0, periods * length, periods + 1)

    fd = finite_differences(x, y)

    spline = cbs(x, y, fd)
    if grid is None:
        grid = np.linspace(0, total_length, grid_points)

    spline_vals = spline(grid)
    if transform is True:
        out = 1 / (1 + np.exp(-spline_vals))
    else:
        out = spline_vals

    return out


periods = 10
length = 2
total_length = 20
np.random.seed(123)
array = np.random.uniform(0, 1, periods + 1)
time = grid = np.linspace(0, total_length, 6000)
spline = get_spline(
    array, periods, length, total_length, grid_points=6000, transform=True
)
spline_ = get_spline(
    array, periods, length, total_length, grid_points=6000, transform=False
)

plt.plot(time, spline, label=r"$\rho^{(i)}(t)$")
plt.scatter(np.linspace(0, 20, 11), array, label=r"$\rho^{(i)}_{j}$")
plt.xlabel("$t$ in Days")
plt.ylabel("$R_0$ multiplier in country i")
plt.xticks(np.linspace(0, 20, 11))
plt.legend()

plt.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/transformed_spline_r_value"
)

plt.plot(time, spline_)
plt.scatter(
    np.linspace(0, 20, 11), np.log(array / (1 - array)), label="Spline parameters"
)
plt.xlabel("Days")
plt.ylabel("Spline values")
plt.xticks(np.linspace(0, 20, 11))
plt.legend()

plt.savefig("/home/manuel/Documents/VaccinationDistribution/paper/images/spline")


#############################################################################
# x = spline_xx.values()
# y_str = [x for x in par_R.keys() if "countryC" in x]
# y = []
# for i in y_str:
#    y.append(par_R[i])
# y = np.array(y)
# s = get_spline(1 / (1+np.exp(-y)), len(y)-1, int(420/9) , 420, grid_points=6000, transform=False)
# time = np.linspace(0, 420, 6000)
# plt.plot(time, s)
