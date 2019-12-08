# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mp
import os
import struct

def load_dat(fn, nx, nz, mod="f"):
    f = open(fn, "rb")
    pic = np.zeros((nz, nx))
    for i in range(nx):
        for j in range(nz):
            data = f.read(4)
            elem = struct.unpack(mod, data)[0]
            pic[j][i] = elem
    f.close()
    return pic




fn = "sei"
nx = 200
nz = 8000
iters = 200
seis_data  = load_dat("record",nx,nz)
g_best_location  = load_dat("g_best_location",2,iters,'i')

seis_mul = load_dat("seis_mul",100,1)



# plt.figure(1)
# plt.figure(figsize=(10, 2))
# plt.figure(dpi=100)

plt.figure(1)

plt.imshow(seis_data,aspect="auto")
for irow in range(0,iters):
    plt.plot(g_best_location[irow,0],g_best_location[irow,1],'r*')
    pass


plt.show()


# plt.figure(2)

# for irow in range(0,99):
    # plt.plot(irow,seis_mul[0,irow],'r*')
    # pass

# plt.show()


