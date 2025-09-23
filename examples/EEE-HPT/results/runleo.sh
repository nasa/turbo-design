#!/bin/sh
mpirun -np 1 leo_acc80 EEE-HPT-0.leo 1 : -np 1 leo_acc80 EEE-HPT-1.leo 1 : -np 1 leo_acc80 EEE-HPT-2.leo 1 : -np 1 leo_acc80 EEE-HPT-3.leo 1

