#include <stdio.h>
#include <stdlib.h>

#include <string.h>

#include <mpi.h>
#include <math.h>

#include "vector3d.h"
#include "savebmp.h"
#include "properties.h"

#define epsilon 0.000000000000000222

typedef struct {
	vec3 position;
	// double velocity;
	// double acceleration;
	// double mass;
} particleVect;

int main(int argc, char *argv[]) {
  if (argc != 10) {
    printf(
      "Usage: %s numParticlesLight numParticleMedium numParticleHeavy numSteps subSteps timeSubStep imageWidth imageHeight imageFilenamePrex\n",
      argv[0]);
  }

  MPI_Init(&argc, &argv);

  int p, my_rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  // variables
  int numParticlesLight = 5;
  int numParticleMedium = 0;
  int numParticleHeavy  = 0;

  int numSteps = 0;
  int subSteps = 0;
  double timeSubStep;

  int width, height;

  unsigned char *image;

  // root node stuff goes here
  if (my_rank == 0) {
    // almost done, just save the image
    //saveBMP(argv[9], image, width, height);
		particleVect lightParticles[numParticlesLight];
		double range = velocityLightMax-velocityLightMin;
		double massRange = massLightMax-massLightMin;
		for(int i = 0; i <numParticlesLight; i++){
			double position[3];
			double x = (drand48()*range)+velocityLightMin;
			double y = (drand48()*range)+velocityLightMin;
			double z = drand48();
			lightParticles[i].position =  vec3(x,y,z);
			// lightParticles[i].position.x = (drand48()*range)+velocityLightMin;
			// lightParticles[i].position.y = (drand48()*range)+velocityLightMin;
			// lightParticles[i].position.z = drand48();
			// lightParticles[i].velocity.x = drand48();
			// lightParticles[i].velocity.y = drand48();
			// lightParticles[i].velocity.z = drand48();
			// lightParticles[i].acceleration.x = 0;
			// lightParticles[i].acceleration.y = 0;
			// lightParticles[i].acceleration.z = 0;
			// lightParticles[i].mass = (drand48()*massRange)+massLightMin;
			printf("i: %d, x: %f, y: %f, z: %f\n", i, lightParticles[i].position.x,lightParticles[i].position.y,lightParticles[i].position.z );
		}
  }
  // all other nodes do this
  else {}

  // free(image);

  MPI_Finalize();
  return 0;
}
