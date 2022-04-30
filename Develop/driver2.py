#!/bin/env python

import os
import numpy as np
import time


#######  SU(2) ######################################
N = 2
t = time.process_time()

# Strong coupling #####################
#for b in np.linspace(0.2,2,10):
#    cmd = "mpirun -np 32 -machinefile procs32 pureSUN -nc=" + str(N) + \
#        " -beta=" + str(b) + \
#        " -L=8x8x8x8 -cold -warms=100 -trajecs=20 -meas=2 >> su"+str(N)+"bloop.out"
#    print(cmd)
#    os.system(cmd)

# Weak coupling #######################
# Change range for each N 

for b in np.arange(4,20,0.5):
    cmd = "mpirun -np 32 -machinefile procs32 pureSUN -nc=" + str(N) + \
        " -beta=" + str(b) + \
        " -L=8x8x8x8 -cold -warms=100 -trajecs=20 -meas=2 >> su"+str(N)+"bloop.out"
    print(cmd)
    os.system(cmd)

#elapsed_time = time.process_time() - t
print("Elapsed Time: {}".format(time.process_time() - t))



########## SU(3) #################################
N = 3
t = time.process_time()

# Strong coupling #####################
for b in np.arange(4.5,6.5,0.1):
    cmd = "mpirun -np 32 -machinefile procs32 pureSUN -nc=" + str(N) + \
        " -beta=" + str(b) + \
        " -L=8x8x8x8 -cold -warms=100 -trajecs=20 -meas=2 >> su"+str(N)+"bloop.out"
    print(cmd)
    os.system(cmd)

# Weak coupling #######################
# Change range for each N 

for b in np.arange(8,20,0.5):
    cmd = "mpirun -np 32 -machinefile procs32 pureSUN -nc=" + str(N) + \
        " -beta=" + str(b) + \
        " -L=8x8x8x8 -cold -warms=100 -trajecs=20 -meas=2 >> su"+str(N)+"bloop.out"
    print(cmd)
    os.system(cmd)

#elapsed_time = time.process_time() - t
print("Elapsed Time: {}".format(time.process_time() - t))




########## SU(4) ###################################
N = 4
t = time.process_time()

# Strong coupling #####################
for b in np.arange(8,11,0.1):
    cmd = "mpirun -np 32 -machinefile procs32 pureSUN -nc=" + str(N) + \
        " -beta=" + str(b) + \
        " -L=8x8x8x8 -cold -warms=100 -trajecs=20 -meas=2 >> su"+str(N)+"bloop.out"
    print(cmd)
    os.system(cmd)

# Weak coupling #######################
# Change range for each N 

for b in np.arange(13,25,0.5):
    cmd = "mpirun -np 32 -machinefile procs32 pureSUN -nc=" + str(N) + \
        " -beta=" + str(b) + \
        " -L=8x8x8x8 -cold -warms=100 -trajecs=20 -meas=2 >> su"+str(N)+"bloop.out"
    print(cmd)
    os.system(cmd)

#elapsed_time = time.process_time() - t
print("Elapsed Time: {}".format(time.process_time() - t))



########## SU(5) #########################################
N = 5
t = time.process_time()

# Strong coupling #####################
#for b in np.linspace(0.2,2,10):
#    cmd = "mpirun -np 32 -machinefile procs32 pureSUN -nc=" + str(N) + \
#        " -beta=" + str(b) + \
#        " -L=8x8x8x8 -cold -warms=100 -trajecs=20 -meas=2 >> su"+str(N)+"bloop.out"
#    print(cmd)
#    os.system(cmd)

# Weak coupling #######################
# Change range for each N 

for b in np.arange(20,30,0.5):
    cmd = "mpirun -np 32 -machinefile procs32 pureSUN -nc=" + str(N) + \
        " -beta=" + str(b) + \
        " -L=8x8x8x8 -cold -warms=100 -trajecs=20 -meas=2 >> su"+str(N)+"bloop.out"
    print(cmd)
    os.system(cmd)

#elapsed_time = time.process_time() - t
print("Elapsed Time: {}".format(time.process_time() - t))



########## SU(10) #########################################
N = 10
t = time.process_time()

# Strong coupling #####################
for b in np.arange(1,30,1):
    cmd = "mpirun -np 32 -machinefile procs32 pureSUN -nc=" + str(N) + \
        " -beta=" + str(b) + \
        " -L=8x8x8x8 -cold -warms=100 -trajecs=20 -meas=2 >> su"+str(N)+"bloop.out"
    print(cmd)
    os.system(cmd)

# Weak coupling #######################
# Change range for each N 

#for b in np.linspace(14,20,13):
#    cmd = "mpirun -np 32 -machinefile procs32 pureSUN -nc=" + str(N) + \
#        " -beta=" + str(b) + \
#        " -L=8x8x8x8 -cold -warms=100 -trajecs=20 -meas=2 >> su"+str(N)+"bloop.out"
#    print(cmd)
#    os.system(cmd)

#elapsed_time = time.process_time() - t
print("Elapsed Time: {}".format(time.process_time() - t))


