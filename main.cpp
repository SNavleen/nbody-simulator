#include <stdio.h>
#include <stdlib.h>

#include <string.h>

#include <mpi.h>
#include <math.h>
#include <time.h> /* time */

#include "vector3d.h"
#include "savebmp.h"
#include "properties.h"

#define epsilon 0.000000000000000222
#define particleRadius 20

using namespace std;


typedef struct {
  double const velocityMin;
  double const velocityMax;
  double const massMin;
  double const massMax;
  vec3 const   colour;
} ParticleProperties;

typedef struct {
  vec3   position;
  vec3   velocity;
  vec3   force;
  vec3   colour;
  double mass;
  int    radius;
} Particle;


void particlesInit(int                startIndex,
                   int                endIndex,
                   int                width,
                   int                height,
                   Particle          *particles,
                   ParticleProperties particleProperties);
void updateParticles(int       numParticles,
                     int       h,
                     Particle *oldParticles,
                     Particle *newParticles);
void imgInit(int            h,
             int            w,
             unsigned char *img);
void imgUpdate(Particle const *particles,
               int             numParticles,
               unsigned char  *img,
               int             height,
               int             width);
void putpixel(int            x,
              int            y,
              unsigned char *img,
              int            width,
              vec3           colour);

const float G = 6.67300E-11;

// MPI variables
int p, my_rank;
// Status of a MPI_Recv
MPI_Status status;
// Request of a MPI_Isend
MPI_Request request;
/* New Datatype for vec3 */
MPI_Datatype mpi_vec3, mpi_Particle, mpiold_Particle[4];
MPI_Aint     offsetsParticle[4], addrParticle[4];

// Input variables
int numParticlesLight, numParticlesMedium, numParticlesHeavy, numSteps, subSteps,
    width, height;
double timeSubStep;

// Particle Informations
int totalParticles;

ParticleProperties lightProperties =
{ velocityLightMin, velocityLightMax, massLightMin, massLightMin, colourLight };
ParticleProperties mediumProperties =
{ velocityMediumMin, velocityMediumMax, massMediumMin, massMediumMin,
  colourMedium };
ParticleProperties heavyProperties =
{ velocityHeavyMin, velocityHeavyMax, massHeavyMin, massHeavyMin, colourHeavy };

Particle *particlesOld, *particlesNew, *particlesNewLocal;

int main(int argc, char *argv[]) {
  if (argc != 10) {
    printf(
      "Usage: %s numParticlesLight numParticleMedium numParticleHeavy numSteps subSteps timeSubStep imageWidth imageHeight imageFilenamePrex\n",
      argv[0]);
  }
  srand48(time(NULL));


  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  // New data type for vec3
  MPI_Type_contiguous(3, MPI_DOUBLE, &mpi_vec3);
  MPI_Type_commit(&mpi_vec3);

  // variables
  numParticlesLight  = atoi(argv[1]);
  numParticlesMedium = atoi(argv[2]);
  numParticlesHeavy  = atoi(argv[3]);

  numSteps    = atoi(argv[4]);
  subSteps    = atoi(argv[5]);
  timeSubStep = atof(argv[6]);

  width  = atoi(argv[7]);
  height = atoi(argv[8]);

  unsigned char *image = (unsigned char *)malloc(height * width * 3);

  totalParticles = numParticlesLight + numParticlesMedium + numParticlesHeavy;

  particlesOld      = (Particle *)malloc(totalParticles * sizeof(Particle));
  particlesNew      = (Particle *)malloc(totalParticles * sizeof(Particle));
  particlesNewLocal = (Particle *)malloc((totalParticles / p) * sizeof(Particle));

  // New data type for Particle struct
  int blockLenParticle[4];
  MPI_Get_address(particlesOld,          &addrParticle[0]);
  MPI_Get_address(&particlesOld->colour, &addrParticle[1]);
  MPI_Get_address(&particlesOld->mass,   &addrParticle[2]);
  MPI_Get_address(&particlesOld->radius, &addrParticle[3]);
  blockLenParticle[0] = 4;
  blockLenParticle[1] = 1;
  blockLenParticle[2] = 1;
  blockLenParticle[3] = 1;
  offsetsParticle[0]  = 0;
  offsetsParticle[1]  = addrParticle[2] - addrParticle[0];
  offsetsParticle[2]  = addrParticle[3] - addrParticle[0];
  offsetsParticle[3]  = sizeof(Particle);
  mpiold_Particle[0]  = mpi_vec3;
  mpiold_Particle[1]  = MPI_DOUBLE;
  mpiold_Particle[2]  = MPI_INT;
  mpiold_Particle[3]  = MPI_UB;
  MPI_Type_struct(4,
                  blockLenParticle,
                  offsetsParticle,
                  mpiold_Particle,
                  &mpi_Particle);
  MPI_Type_commit(&mpi_Particle);


  // root node stuff goes here
  if (my_rank == 0) {
    // printf("Rank %d - totalParticles: %d\n", my_rank, totalParticles);
    particlesInit(0, numParticlesLight,
                  width,
                  height,
                  particlesOld,
                  lightProperties);
    particlesInit(numParticlesLight, numParticlesLight + numParticlesMedium,
                  width,
                  height,
                  particlesOld,
                  mediumProperties);
    particlesInit(numParticlesLight + numParticlesMedium, totalParticles,
                  width,
                  height,
                  particlesOld,
                  heavyProperties);
    for (int dest = 0; dest < p; dest++) {
      MPI_Isend(&totalParticles,
                1,
                MPI_INT,
                dest,
                0,
                MPI_COMM_WORLD,
                &request);
    }
  } else {
    MPI_Recv(&totalParticles, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
  }
  for (int i = 0; i < numSteps * subSteps; i++) {
    MPI_Bcast(&particlesOld[0], totalParticles, mpi_Particle, 0, MPI_COMM_WORLD);
    //   printf("After Rank %d (particlesOld position) - x: %f, y: %f, z: %f\n",
    //          my_rank,
    //          particlesOld[0].position.x,
    //          particlesOld[0].position.y,
    //          particlesOld[0].position.z);
    //  printf("After Rank %d (particlesOld velocity) - x: %f, y: %f, z: %f\n",
    //         my_rank,
    //         particlesOld[0].velocity.x,
    //         particlesOld[0].velocity.y,
    //         particlesOld[0].velocity.z);
    //   printf("After Rank %d (particlesOld force) - x: %f, y: %f, z: %f\n",
    //          my_rank,
    //          particlesOld[0].force.x,
    //          particlesOld[0].force.y,
    //          particlesOld[0].force.z);
    //  printf("After Rank %d (particlesOld colour) - x: %f, y: %f, z: %f\n",
    //         my_rank,
    //         particlesOld[0].colour.x,
    //         particlesOld[0].colour.y,
    //         particlesOld[0].colour.z);
    // printf("After Rank %d (particlesOld) - mass: %f, radius: %d\n",
    //                   my_rank,
    //                   particlesOld[0].mass,
    //                   particlesOld[0].radius);
    updateParticles(totalParticles,
                    timeSubStep,
                    particlesOld,
                    particlesNewLocal);
    MPI_Gather(&particlesNewLocal[0], totalParticles / p, mpi_Particle,
               &particlesNew[0], totalParticles / p, mpi_Particle, 0,
               MPI_COMM_WORLD);
    if (my_rank == 0) {
      swap(*particlesOld, *particlesNew);
      if (i % numSteps == 0) {
        // printf("Rank %d - i: %d\n", my_rank, i);
        imgInit(height, width, image);
        imgUpdate(particlesOld, totalParticles, image, height, width);
        // int imageNumber  = i / numSteps;
        string imageName = argv[9];
        string extention = ".bmp";
        // imageName += "_"
        // imageName += to_string(imageNumber);
        imageName += extention;
        saveBMP(imageName.c_str(), image, width, height);
      }
    }
  }

  // free(image);

  MPI_Finalize();
  return 0;
}

void particlesInit(int                startIndex,
                   int                endIndex,
                   int                width,
                   int                height,
                   Particle          *particles,
                   ParticleProperties particleProperties) {
  double range = particleProperties.velocityMax -
                 particleProperties.velocityMin;
  double massRange = particleProperties.massMax - particleProperties.massMax;

  for (int i = startIndex; i < endIndex; i++) {
    double x = drand48() * width;
    double y = drand48() * height;
    // double z = drand48();
    double z = 0;
    particles[i].position = vec3(x, y, z);

    double vx = (drand48() * range) + particleProperties.velocityMin;
    double vy = (drand48() * range) + particleProperties.velocityMin;
    // double vz = drand48();
    double vz = 0;
    particles[i].velocity = vec3(vx, vy, vz);

    particles[i].mass = (drand48() * massRange) +
                        particleProperties.massMax;

    particles[i].colour = particleProperties.colour;

    particles[i].radius = particleRadius;
    // printf("Rank %d - i: %d, x: %f, y: %f, z: %f\n",
    //        my_rank,
    //        i,
    //        particles[i].position.x,
    //        particles[i].position.y,
    //        particles[i].position.z);
  }
}

void updateParticles(int       numParticles,
                     int       h,
                     Particle *oldParticles,
                     Particle *newParticles) {
  int startIndex = my_rank * (numParticles / p);
  int endIndex   = startIndex + (numParticles / p);
  int i          = 0;

  // Basic alogirthm
  for (int q = startIndex; q < endIndex; q++) {
    // for each other particle that is not q in the same N particles
    for (int k = 0; k < numParticles; k++) {
      if (k != q) {
        vec3 diff = oldParticles[q].position - oldParticles[k].position;
        oldParticles[q].force += diff *
                                 (oldParticles[k].mass / diff.Magnitude());
      }
    }
    oldParticles[q].force = oldParticles[q].force * (-G * oldParticles[q].mass);

    newParticles[i].velocity = oldParticles[q].velocity +
                               ((oldParticles[q].force * h) /
                                oldParticles[q].mass);

    newParticles[i].position = oldParticles[q].position +
                               (oldParticles[q].velocity * h);

    newParticles[i].colour = oldParticles[q].colour;

    newParticles[i].radius = oldParticles[q].radius;
    i++;
  }
}

void imgInit(int            h,
             int            w,
             unsigned char *img) {
  for (int j = 0; j < h; j++) {
    for (int i = 0; i < w; i++) {
      img[(j * w + i) * 3 + 0] = 0;
      img[(j * w + i) * 3 + 1] = 0;
      img[(j * w + i) * 3 + 2] = 0;
    }
  }
}

void imgUpdate(Particle const *particles,
               int             numParticles,
               unsigned char  *img,
               int             height,
               int             width) {
  for (int q = 0; q < numParticles; q++) {
    int y    =  particles[q].position.y - particles[q].radius;
    int yEnd =  particles[q].position.y + particles[q].radius;
    int x    =  particles[q].position.x - particles[q].radius;
    int xEnd =  particles[q].position.x + particles[q].radius;
    int r    = particles[q].radius;

    long rIndex, gIndex, bIndex;
    for (int h = y; h <= yEnd; h++) {
      for (int w = x; w <= xEnd; w++) {
        if ((h < 0) || (h >= height)) continue;
        if ((w < 0) || (w >= width)) continue;
        // putpixel(w, h, img, width, particles[q].colour);
        if (sqrt(w * w + h * h) > r * r) {
          // Get the RGB indexes
          rIndex = (h * width + w) * 3 + 0;
          gIndex = (h * width + w) * 3 + 1;
          bIndex = (h * width + w) * 3 + 2;

          img[rIndex] = 255;
          img[gIndex] = 255;
          img[bIndex] = 255;
        } else {
          // Get the RGB indexes
          rIndex = (h * width + w) * 3 + 0;
          gIndex = (h * width + w) * 3 + 1;
          bIndex = (h * width + w) * 3 + 2;

          img[rIndex] = particles[q].colour.z * 255;
          img[gIndex] = particles[q].colour.y * 255;
          img[bIndex] = particles[q].colour.x * 255;
        }
      }
    }
  }
}

void putpixel(int x, int y, unsigned char  *img, int width, vec3 colour) {
  long rIndex, gIndex, bIndex;

  // Get the RGB indexes
  rIndex = (y * width + x) * 3 + 0;
  gIndex = (y * width + x) * 3 + 1;
  bIndex = (y * width + x) * 3 + 2;

  img[rIndex] = colour.z * 255;
  img[gIndex] = colour.y * 255;
  img[bIndex] = colour.x * 255;
}
