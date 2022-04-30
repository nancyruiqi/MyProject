#!/bin/env python

import os
import numpy as np
import time

N = 5

t = time.process_time()

# Strong coupling #####################
for b in np.linspace(0.2,2,10):
    cmd = "mpirun -np 32 -machinefile procs32 pureSUN -nc=" + str(N) + \
        " -beta=" + str(b) + \
        " -L=8x8x8x8 -cold -warms=100 -trajecs=20 -meas=2 >> su"+str(N)+"bloop.out"
    print(cmd)
    os.system(cmd)

# Weak coupling #######################
# Change range for each N 

for b in np.linspace(14,20,13):
    cmd = "mpirun -np 32 -machinefile procs32 pureSUN -nc=" + str(N) + \
        " -beta=" + str(b) + \
        " -L=8x8x8x8 -cold -warms=100 -trajecs=20 -meas=2 >> su"+str(N)+"bloop.out"
    print(cmd)
    os.system(cmd)


#elapsed_time = time.process_time() - t
print("Elapsed Time: {}".format(time.process_time() - t))
