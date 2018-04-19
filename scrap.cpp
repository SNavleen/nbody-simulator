/*********************************************************
THIS FILE IS JUST FOR NOTES AND PSEUDO CODE - DO NOT SUBMIT THIS FILE
*********************************************************/

// from slides:
// G = 6.673 × 10−11m3/(kg · s2) is the gravitational constant
// from google:
const float G  = 6.67300E-11;

// Number of particles
int N = 10;

// h = subtimesteps
h = subtimesteps

/*********************************************************/

// on a time step from t[i] to t[i+1]
/* particle q properties */
q.speed[i+1] = q.speed[i] + h * q.velocity[i];
q.velocity[i+1] = q.velocity[i] + h * q.F[i] / q.mass; // velocity = speed`
q.acceleration = ?; // acceleration = speed``

// other particle k exerts force on current particle q
q.force[i] = (-G) * (q.mass * k.mass) / (distance between q and k) * (q.speed - k.speed);

// for each particle q in N particles
// basic alogirthm
for (x=0; x<N; x++){
  // q = x-th particle

  q.force[i] =0;

  // for each other particle that is not q in the same N particles
  for(y=0; y<N; y++){
    // k = y-th particle

    if (k!=q){
      // compute f(qk,i) -- see slide because i'm confused
      q.force[i] = q.force[i] + f(qk,i)

    }
  }
}
// reduced alogirthm
for (q=0; q<N; q++){
  p.force[i] =0;

  //
}
