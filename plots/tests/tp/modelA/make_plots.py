import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


df = pd.read_csv("phi.csv")
fig, ax = plt.subplots()

ax.plot(
    df['phi1'],
    df['phi2'])

plt.savefig("phi.pdf")

