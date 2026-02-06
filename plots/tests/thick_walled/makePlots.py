import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate


# ===============================================


df = pd.read_csv("pot_of_phi1_phi2.csv")
fig, ax = plt.subplots()

# xmins = [
#     [0, 1.057445],
#     [0, 1.553]]

number_points = 200
xmax = np.max(df["phi1"])
xmin = np.min(df["phi1"])
xscale = xmax - xmin
ymax = np.max(df["phi2"])
ymin = np.min(df["phi2"])
yscale = ymax - ymin
scale = np.array([xscale, yscale])
x = np.linspace(xmin, xmax, number_points)
y = np.linspace(ymin, ymax, number_points)
mat = np.array([
    df["phi1"].to_numpy(),
    df["phi2"].to_numpy(),
    df["pot "].to_numpy()]).T
X, Y = np.meshgrid(x, y)
Z = interpolate.griddata(
    mat[:,0:2] / scale,
    mat[:,2],
    (X / xscale, Y / yscale),
    method="linear")
Zmin = np.min(Z)
Zmax = np.max(Z)
lvls = []
for i in range(0, 20):
    lvls.append(Zmin + i * 0.01 * (Zmax - Zmin))
cs = ax.contour(
    X, Y, Z,
    levels=lvls,
)

dg = pd.read_csv("phi_of_rho.csv")
ax.plot(
    dg['phi1'],
    dg['phi2'])

# ax.plot(
#     xmins[0],
#     xmins[1])

plt.savefig("potential.pdf")


# ===============================================


df = pd.read_csv("x_of_rho.csv")

fig, ax = plt.subplots()

ax.plot(
    df['rho '],
    df['x   '])

ax.plot(
    df['rho '],
    df['xdot'])

plt.savefig("x_of_rho.pdf")


# ===============================================


df = pd.read_csv("pot_of_rho.csv")

fig, ax = plt.subplots()

ax.plot(
    df['rho'],
    df['pot'])

plt.savefig("pot_of_rho.pdf")


# ===============================================


df = pd.read_csv("pot_of_x.csv")

fig, ax = plt.subplots()

ax.scatter(
    df['x'],
    df['p'],
    s=6)

plt.savefig("pot_of_x.pdf")


# ===============================================


df = pd.read_csv("phi_of_x.csv")

fig, ax = plt.subplots()

ax.plot(
    df['x   '],
    df['phi1'],
    color="C0")

ax.plot(
    df['x   '],
    df['phi2'],
    color="C0")

plt.savefig("phi_of_x.pdf")


# ===============================================


df = pd.read_csv("forces_of_x.csv")

fig, ax = plt.subplots()

ax.plot(
    df['x '],
    df['N1'],
    color="C0")

ax.plot(
    df['x '],
    df['N2'],
    color="C0")

ax.plot(
    df['x '],
    df['P1'],
    color="C1")

ax.plot(
    df['x '],
    df['P2'],
    color="C1")

plt.savefig("forces_of_x.pdf")

