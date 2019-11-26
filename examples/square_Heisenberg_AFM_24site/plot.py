#This code is part of DIETA
#
#Authored by Kirill Alpin

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

f = open('spec_output.txt')
lines = f.readlines()
f.close()

max_omega_step = int(lines[0].split('\t')[0])
k_num_steps = int(lines[0].split('\t')[1])
num_sym = int(lines[0].split('\t')[2])



stuff = []
kvals = []
i = 0
nowstuff = np.zeros((k_num_steps * (num_sym - 1)))
for l in lines[1:]:
    spl = l.split('\t')
    nowstuff[i] = float(spl[4])
    if len(kvals) < k_num_steps * (num_sym - 1):
        kvals.append(np.array([float(spl[1]),float(spl[2]),float(spl[3])]))
    i+=1
    if i >= k_num_steps * (num_sym - 1):
        i = 0
        stuff.append(nowstuff)
        nowstuff = np.zeros((k_num_steps * (num_sym - 1)))
    
a = np.vstack(stuff)[::-1,:]

a = np.where(a > 100, 100, a)
a = np.where(a < 0.0, 0.0, a)

ax = plt.gca()
im = ax.matshow(a, extent=[0, 240.0/40.0, 0.0, 3.0])
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.05)
cb = plt.colorbar(im, cax=cax)
cb.set_ticks([])
ax.plot((60/40.0, 60/40.0), (0, 120/40.0), 'k-')
ax.plot((120/40.0, 120/40.0), (0, 120/40.0), 'k-')
ax.text(6.0-0.05, -0.15, "K'")
ax.text(4.5-0.05, -0.15, "$\Gamma$\'")
ax.text(3.0-0.05, -0.15, "K")
ax.text(1.5-0.05, -0.15, "M")
ax.text(0.0-0.05, -0.15, r'$\Gamma$')
ax.set_xticks([])
ax.plot((180/40.0, 180/40.0), (0, 120/40.0), 'k-')

plt.savefig('spec_output.pdf', bbox_inches='tight')
plt.show()
