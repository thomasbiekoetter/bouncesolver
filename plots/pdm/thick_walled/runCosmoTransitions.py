import numpy as np
import pandas as pd
from scipy.optimize import fmin
from cosmoTransitions import pathDeformation
from cosmoTransitions.helper_functions import gradientFunction
from ableiter.first import First



a = 1.8
b = 0.2
c = 0.0


def V(X):
    global c
    x0, x1 = X[...,0], X[...,1]
    y = (x0 ** 2 + x1 ** 2) * (
      a * (x0 - 1.0) ** 2 +
      b * (x1 - 1.0) ** 2 - c)
    return y


gradV = gradientFunction(V, 1e-6, 2, 4)


def calc_bounce_with_cosmo(path):
    Y = pathDeformation.fullTunneling(
                    path,
                    V,
                    gradV,
                    tunneling_init_params={'alpha': 2},
                    deformation_deform_params={'verbose': 0},
                    verbose=True)
    return Y


def main():
    global c
    dcs = []
    for c in np.linspace(0.3, 0.5, 101):
        x_false = fmin(V, x0=[0, 0])
        x_true = fmin(V, x0=[1, 1])
        path = np.array([x_true, x_false])
        try:
            co = calc_bounce_with_cosmo(path)
            dcs.append({'c': c, 'S': co.action})
            print(c, co.action)
        except ValueError:
            print("cosmo failed.")
    df = pd.DataFrame(dcs)
    df.to_csv("data_CosmoTransitions.csv")







main()
