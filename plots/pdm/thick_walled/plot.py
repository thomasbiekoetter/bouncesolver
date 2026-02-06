import os
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


mpl.use('pgf')
pgf_dc = {
    'text.usetex': True,
    'pgf.rcfonts': False
}
mpl.rcParams.update(pgf_dc)


df = pd.read_csv("data.csv")
dg = pd.read_csv("data_CosmoTransitions.csv")

fig, ax = plt.subplots(
    2, 1,
    sharex=True,
    gridspec_kw={'height_ratios': [2, 1]},
    figsize=(4.4, 3.2))

fig.subplots_adjust(hspace=0.05)

ax[0].plot(
    df['c '],
    df['S0'],
    lw=2,
    color='black',
    alpha=0.3,
    label=r"$\textrm{Straight-path appr.}$")

ax[0].plot(
    df['c '],
    df['S '],
    lw=2,
    label=r"$\textrm{\texttt{bouncesolver}}$")

ax[0].plot(
    dg['c'],
    dg['S'],
    ls='--',
    label=r"$\textrm{\texttt{CosmoTransitions}}$")

x = 1e2 * (df['S '].to_numpy() - dg['S'].to_numpy()) / df['S '].to_numpy()
ax[1].plot(
    df['c '],
    x,
    label=r"$\textrm{Rel.~difference in \%}$"
)

ax[1].axhline(
    0,
    ls=":",
    color="black")

ax[0].set_yscale('log')
ax[0].set_ylabel(r"$S_E$")

ax[0].set_xlabel(r"$c$")

ax[0].set_xlim(df['c '].min(), df['c '].max())

ax[1].set_xlabel(r"$c$")
majorLocator = MultipleLocator(0.1)
minorLocator = MultipleLocator(0.01)
ax[1].xaxis.set_major_locator(majorLocator)
ax[1].xaxis.set_minor_locator(minorLocator)
ax[1].tick_params(
    axis='x',
    direction='in',
    which='both',
    top=True)

majorLocator = MultipleLocator(1)
minorLocator = MultipleLocator(0.2)
ax[1].set_ylim(-1, 5)
ax[1].yaxis.set_major_locator(majorLocator)
ax[1].yaxis.set_minor_locator(minorLocator)
ax[1].tick_params(
    axis='y',
    direction='in',
    which='both',
    right=True)

ax[0].legend(
    frameon=False)

ax[1].legend(
#   loc=(0.5, 0.02),
    frameon=False)

ax[0].text(
    0.104,
    6.5,
    r"$V(\vec{\phi}) = |\vec{\phi}|^2 [1.8(\phi_1 - 1)^2 + 0.2(\phi_2 - 1)^2 - c]$",
    fontsize="small"
)

plt.savefig("thick_walled.pdf")

os.system("pdfcrop thick_walled.pdf thick_walled.pdf")
