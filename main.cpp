#include <stdio.h>
#include <stdlib.h>

#include <string.h>

#include <mpi.h>
#include <math.h>
#include <time.h>       /* time */


#include "vector3d.h"
#include "savebmp.h"
#include "properties.h"

#define epsilon 0.000000000000000222

typedef struct {
	vec3 position;
	vec3 velocity;
	vec3 acceleration;
	double mass;
} particleVect;

void lightInit(int num, int width, int height){
	particleVect lightParticles[num];
	double range = velocityLightMax-velocityLightMin;
	double massRange = massLightMax-massLightMin;
	for(int i = 0; i <num; i++){
		double x = drand48()*width;
		double y = drand48()*height;
		double z = drand48();
		lightParticles[i].position = vec3(x,y,z);
		double vx = (drand48()*range)+velocityLightMin;
		double vy = (drand48()*range)+velocityLightMin;
		double vz = drand48();
		lightParticles[i].velocity = vec3(vx,vy,vz);
		lightParticles[i].acceleration = vec3(0,0,0);
		lightParticles[i].mass = (drand48()*massRange)+massLightMin;
		//printf("i: %d, x: %f, y: %f, z: %f\n", i, lightParticles[i].position.x,lightParticles[i].position.y,lightParticles[i].position.z );
	}
}
void mediumInit(int num, int width, int height){
	particleVect mediumParticles[num];
	double range = velocityMediumMax-velocityMediumMin;
	double massRange = massMediumMax-massMediumMin;
	for(int i = 0; i <num; i++){
		double x = drand48()*width;
		double y = drand48()*height;
		double z = drand48();
		mediumParticles[i].position = vec3(x,y,z);
		double vx = (drand48()*range)+velocityMediumMin;
		double vy = (drand48()*range)+velocityMediumMin;
		double vz = drand48();
		mediumParticles[i].velocity = vec3(vx,vy,vz);
		mediumParticles[i].acceleration = vec3(0,0,0);
		mediumParticles[i].mass = (drand48()*massRange)+massMediumMin;
		//printf("i: %d, x: %f, y: %f, z: %f\n", i, mediumParticles[i].position.x,mediumParticles[i].position.y,mediumParticles[i].position.z );
	}
}
void heavyInit(int num, int width, int height){
	particleVect heavyParticles[num];
	double range = velocityHeavyMax-velocityHeavyMin;
	double massRange = massHeavyMax-massHeavyMin;
	for(int i = 0; i <num; i++){
		double x = drand48()*width;
		double y = drand48()*height;
		double z = drand48();
		heavyParticles[i].position = vec3(x,y,z);
		double vx = (drand48()*range)+velocityHeavyMin;
		double vy = (drand48()*range)+velocityHeavyMin;
		double vz = drand48();
		heavyParticles[i].velocity = vec3(vx,vy,vz);
		heavyParticles[i].acceleration = vec3(0,0,0);
		heavyParticles[i].mass = (drand48()*massRange)+massHeavyMin;
		//printf("i: %d, x: %f, y: %f, z: %f\n", i, heavyParticles[i].position.x,heavyParticles[i].position.y,heavyParticles[i].position.z );
	}
}


int main(int argc, char *argv[]) {
  if (argc != 10) {
    printf(
      "Usage: %s numParticlesLight numParticleMedium numParticleHeavy numSteps subSteps timeSubStep imageWidth imageHeight imageFilenamePrex\n",
      argv[0]);
  }
	srand48(time(NULL));


  MPI_Init(&argc, &argv);

  int p, my_rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  // variables
  int numParticlesLight = atoi(argv[1]);
  int numParticlesMedium = atoi(argv[2]);
  int numParticlesHeavy  = atoi(argv[3]);

  int numSteps = atoi(argv[4]);
  int subSteps = atoi(argv[5]);
  double timeSubStep = atof(argv[6]);

  int width = atoi(argv[7]);
	int height = atoi(argv[8]);

  unsigned char *image;

	lightInit(numParticlesLight,width,height);
	mediumInit(numParticlesMedium,width,height);
	heavyInit(numParticlesHeavy,width,height);

  // root node stuff goes here
  if (my_rank == 0) {
    // almost done, just save the image
    //saveBMP(argv[9], image, width, height);

  }
  // all other nodes do this
  else {}

  // free(image);

  MPI_Finalize();
  return 0;
}
