import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import interpolate


df = pd.read_csv("pot_of_phi1_phi2.csv")
fig, ax = plt.subplots()
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
dg = pd.read_csv("1_phi_of_rho.csv")
ax.plot(
    dg['phi1'],
    dg['phi2'],
    label='maxiter = 1')
dg = pd.read_csv("40_phi_of_rho.csv")
ax.plot(
    dg['phi1'],
    dg['phi2'],
    label='maxiter = 40')
ax.legend()
plt.savefig("potential.pdf")




fig, ax = plt.subplots()
df = pd.read_csv("1_pot_of_rho.csv")
ax.plot(
    df['rho'],
    df['pot'],
    label='maxiter = 1')
df = pd.read_csv("40_pot_of_rho.csv")
ax.plot(
    df['rho'],
    df['pot'],
    label='maxiter = 40')
ax.legend()
plt.savefig("pot_of_rho.pdf")





fig, ax = plt.subplots()
df = pd.read_csv("1_pot_of_x.csv")
ax.plot(
    df['x'],
    df['p'],
    label='maxiter = 1')
df = pd.read_csv("40_pot_of_x.csv")
ax.plot(
    df['x'],
    df['p'],
    label='maxiter = 40')
ax.legend()
plt.savefig("pot_of_x.pdf")





fig, ax = plt.subplots()
df = pd.read_csv("1_x_of_rho.csv")
ax.plot(
    df['rho '],
    df['x   '],
    color='C0',
    label='x, maxiter = 1')
ax.plot(
    df['rho '],
    df['xdot'],
    ls=':',
    color='C0',
    label='xdot, maxiter = 1')
df = pd.read_csv("40_x_of_rho.csv")
ax.plot(
    df['rho '],
    df['x   '],
    color='C1',
    label='x, maxiter = 40')
ax.plot(
    df['rho '],
    df['xdot'],
    ls=':',
    color='C1',
    label='xdot, maxiter = 40')
ax.legend()
plt.savefig("x_of_rho.pdf")





fig, ax = plt.subplots()
df = pd.read_csv("1_phi_of_x.csv")
ax.plot(
    df['x   '],
    df['phi1'],
    label='phi1, maxiter = 1',
    color="C0")
ax.plot(
    df['x   '],
    df['phi2'],
    label='phi2, maxiter = 1',
    ls='--',
    color="C0")
df = pd.read_csv("40_phi_of_x.csv")
ax.plot(
    df['x   '],
    df['phi1'],
    label='phi1, maxiter = 40',
    color="C1")
ax.plot(
    df['x   '],
    df['phi2'],
    ls='--',
    label='phi2, maxiter = 40',
    color="C1")
ax.legend()
plt.savefig("phi_of_x.pdf")






fig, ax = plt.subplots()
df = pd.read_csv("1_forces_of_x.csv")
ax.plot(
    df['x '],
    df['N1'],
    label='N1, maxiter = 1',
    color="C0")
ax.plot(
    df['x '],
    df['N2'],
    ls='--',
    label='N2, maxiter = 1',
    color="C0")
df = pd.read_csv("40_forces_of_x.csv")
ax.plot(
    df['x '],
    df['N1'],
    label='N1, maxiter = 40',
    color="C1")
ax.plot(
    df['x '],
    df['N2'],
    ls='--',
    label='N2, maxiter = 40',
    color="C1")
ax.legend()
plt.savefig("normalforces_of_x.pdf")





fig, ax = plt.subplots()
df = pd.read_csv("1_forces_of_x.csv")
ax.plot(
    df['x '],
    df['P1'],
    label='P1, maxiter = 1',
    color="C0")
ax.plot(
    df['x '],
    df['P2'],
    ls='--',
    label='P2, maxiter = 1',
    color="C0")
df = pd.read_csv("40_forces_of_x.csv")
ax.plot(
    df['x '],
    df['P1'],
    label='P1, maxiter = 40',
    color="C1")
ax.plot(
    df['x '],
    df['P2'],
    ls='--',
    label='P2, maxiter = 40',
    color="C1")
ax.legend()
plt.savefig("orthoforces_of_x.pdf")





