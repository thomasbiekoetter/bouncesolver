import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("relax.csv")
dg = pd.read_csv("x_of_rho.csv")

fig, ax = plt.subplots()

ax.plot(
    df['rho'],
    df['xso'])

ax.plot(
    df['rho'],
    df['xdo'])

# ax.plot(
#     dg['rho '],
#     dg['x   '],
#     color='magenta',
#     alpha=0.5)

ax.plot(
    df['rho'],
    df['xin'],
    ls=":")

ax.set_ylim(-0.5, 2.0)

plt.savefig("relax.pdf")

