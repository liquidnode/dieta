#This code is part of DIETA
#
#Authored by Kirill Alpin

import numpy as np
import matplotlib.pyplot as plt

with open("sweep2.txt") as f:
    content = f.readlines()
    
x =[]
ls = []
for l in content:
    st = l.split('\t')
    x.append(np.abs(float(st[0])))
    ls.append(np.abs(float(st[-1])))
    
lt = '-'
plt.plot(x, ls, lt)
plt.xlabel("$t$")
plt.ylabel("Superfluid density")
plt.savefig('sweep2.pdf')
plt.show()
