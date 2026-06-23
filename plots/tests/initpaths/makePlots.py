import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate


df = pd.read_csv("pot_of_phi1_phi2.csv")

fig, ax = plt.subplots(
    figsize=(4.2, 4))

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
for i in range(0, 30):
    lvls.append(Zmin + i * 0.003 * (Zmax - Zmin))
cs = ax.contour(
    X, Y, Z,
    levels=lvls,
    alpha=0.3,
    cmap='plasma'
)

dg = pd.read_csv("phi_of_rho_1_apprx.csv")
ax.plot(
    dg['phi1'],
    dg['phi2'],
    label='init_path_mode=1')

dg = pd.read_csv("phi_of_rho_2_apprx.csv")
ax.plot(
    dg['phi1'],
    dg['phi2'],
    label='init_path_mode=2')

dg = pd.read_csv("phi_of_rho_3_apprx.csv")
ax.plot(
    dg['phi1'],
    dg['phi2'],
    label='init_path_mode=3')

ax.set_xlim(-2, 40)
ax.set_ylim(-2, 40)

ax.text(
    0, 37,
    r'Initial paths')

ax.legend(frameon=False)

plt.savefig("initial_paths.pdf")





fig, ax = plt.subplots(
    figsize=(4.2, 4))

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
for i in range(0, 30):
    lvls.append(Zmin + i * 0.003 * (Zmax - Zmin))
cs = ax.contour(
    X, Y, Z,
    levels=lvls,
    alpha=0.3,
    cmap='plasma'
)

dg = pd.read_csv("phi_of_rho_1_apprx.csv")
ax.plot(
    dg['phi1'],
    dg['phi2'],
    color='C0',
    ls='--',
    label='init_path_mode=1, initial')

dg = pd.read_csv("phi_of_rho_1_deform.csv")
ax.plot(
    dg['phi1'],
    dg['phi2'],
    color='C0',
    ls='-',
    label='init_path_mode=1, deformed')

dg = pd.read_csv("phi_of_rho_3_apprx.csv")
ax.plot(
    dg['phi1'],
    dg['phi2'],
    color='C1',
    ls='--',
    alpha=0.6,
    label='init_path_mode=3, initial')

dg = pd.read_csv("phi_of_rho_3_deform.csv")
ax.plot(
    dg['phi1'],
    dg['phi2'],
    color='C1',
    ls='-',
    alpha=0.6,
    label='init_path_mode=3, deformed')

ax.set_xlim(-2, 35)
ax.set_ylim(-2, 35)

ax.legend(frameon=False)

plt.savefig("deformed_paths.pdf")



fig, ax = plt.subplots(
    figsize=(4.6, 3.6))

dg = pd.read_csv("phi_of_rho_1_apprx.csv")
ax.plot(
    dg['rho '],
    dg['phi1'],
    ls='--',
    alpha=0.4,
    label='init_path_mode=1, init',
    color='C0')
ax.plot(
    dg['rho '],
    dg['phi2'],
    ls='--',
    alpha=0.4,
    color='C0')

dg = pd.read_csv("phi_of_rho_1_deform.csv")
ax.plot(
    dg['rho '],
    dg['phi1'],
    alpha=0.4,
    label='init_path_mode=1, deformed',
    color='C0')
ax.plot(
    dg['rho '],
    dg['phi2'],
    alpha=0.4,
    color='C0')

dg = pd.read_csv("phi_of_rho_3_apprx.csv")
ax.plot(
    dg['rho '],
    dg['phi1'],
    ls='--',
    alpha=0.4,
    label='init_path_mode=3, init',
    color='C1')
ax.plot(
    dg['rho '],
    dg['phi2'],
    ls='--',
    alpha=0.4,
    color='C1')

dg = pd.read_csv("phi_of_rho_3_deform.csv")
ax.plot(
    dg['rho '],
    dg['phi1'],
    alpha=0.4,
    label='init_path_mode=3, deformed',
    color='C1')
ax.plot(
    dg['rho '],
    dg['phi2'],
    alpha=0.4,
    color='C1')

ax.legend()

ax.set_xlim(0, 0.55)

plt.savefig("phi_of_rho.pdf")
