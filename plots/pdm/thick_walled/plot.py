import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_csv("data.csv")
dg = pd.read_csv("data_CosmoTransitions.csv")

fig, ax = plt.subplots(
    figsize=(4, 4))

plt.plot(
    df['c '],
    df['S '],
    lw=2,
    label="bouncesolver")

df = df[df['S0'] < 100]
plt.plot(
    df['c '],
    df['S0'],
    lw=2,
    color='C1',
    label="bouncesolver -- straight path")

plt.plot(
    dg['c'],
    dg['S'],
    ls='--',
    lw=2,
    color='magenta',
    label="CosmoTransitions")

plt.xlabel("c")
plt.ylabel("S")

plt.legend()

plt.savefig("thick_walled.pdf")

