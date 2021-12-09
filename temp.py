import numpy as np
import scipy
from scipy import signal
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.pyplot import gca
from matplotlib.ticker import StrMethodFormatter

font = FontProperties()
font.set_family('serif')
font.set_name('Times New Roman')
font.set_style('normal')
font.set_size(14)

Fs = 1000.0
Dt = 1/Fs
T = 1
t = np.arange(0,T,Dt)
#print(t[:2])
fr = 1
Amp = .2

y = np.sin(2 * np.pi * fr * t)*Amp


#%%%%

fig, ax = plt.subplots(figsize=(5, 3))
fig.subplots_adjust(bottom=0.2, left=0.2)
ax.plot(t,y)
ax.set_xlabel('Time (s)', fontproperties=font)
ax.set_ylabel('Signal (m)', fontproperties=font)
ax.set_xticklabels(ax.get_xticks(), fontproperties=font)
ax.set_yticklabels(ax.get_yticks(), fontproperties=font)
plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}')) # 2 decimal places
plt.gca().xaxis.set_major_formatter(StrMethodFormatter('{x:,.2f}')) # 2 decimal places
plt.savefig("D:\PhD\People\Muhammad Waqar\Python\sample.jpg", dpi=300)








