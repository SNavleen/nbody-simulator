#!/bin/bash
make

# mpirun -np numberOfProcessors ./x.project numParticlesLight numParticleMedium numParticleHeavy numSteps subSteps timeSubStep imageWidth imageHeight imageFilenamePrex

# mpirun -np 6 ./x.project 6000 6000 6000 14400 10 0.016 2550 1440 video
mpirun -np 1 ./x.project 250 100 150 14400 10 0.016 2550 1440 video
mpirun -np 2 ./x.project 1000 500 500 14400 10 0.016 2550 1440 video
mpirun -np 4 ./x.project 2000 1000 1000 14400 10 0.016 2550 1440 video
mpirun -np 8 ./x.project 1000 4000 3000 14400 10 0.016 2550 1440 video
mpirun -np 16 ./x.project 6000 4000 5000 14400 10 0.016 2550 1440 video
mpirun -np 32 ./x.project 12000 10000 10000 14400 10 0.016 2550 1440 video
