#!/bin/bash
make

# mpirun -np numberOfProcessors ./x.project numParticlesLight numParticleMedium numParticleHeavy numSteps subSteps timeSubStep imageWidth imageHeight imageFilenamePrex
mpirun -np 1 ./x.project 100 100 100 10 10 1 1920 1080 simulation
# mpirun -np 1 ./x.project 6000 3000 1000 900 50 0.5 1920 1080 simulation
# mpirun -np 1 ./x.project 6000 3000 1000 900 50 0.5 1920 1080 simulation
# mpirun -np 1 ./x.project 6000 3000 1000 900 50 0.5 1920 1080 simulation
# mpirun -np 1 ./x.project 6000 3000 1000 900 50 0.5 1920 1080 simulation
# mpirun -np 1 ./x.project 6000 3000 1000 900 50 0.5 1920 1080 simulation

# mpirun -np 1 ./x.project 6000 6000 6000 14400 10 0.016 2550 1440 video
