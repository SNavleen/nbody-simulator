#!/bin/bash
make
mpirun -np 1 ./x.project 100 100 100 900 50 0.5 1920 1080 simulation
# mpirun -np 1 ./x.project 6000 3000 1000 900 50 0.5 1920 1080 simulation
# mpirun -np 1 ./x.project 6000 3000 1000 900 50 0.5 1920 1080 simulation
# mpirun -np 1 ./x.project 6000 3000 1000 900 50 0.5 1920 1080 simulation
# mpirun -np 1 ./x.project 6000 3000 1000 900 50 0.5 1920 1080 simulation
# mpirun -np 1 ./x.project 6000 3000 1000 900 50 0.5 1920 1080 simulation


#run command to make the video, 60 fps
# mpirun -np 1 ./x.project 6000 6000 6000 14400 10 0.016 2550 1440 video
