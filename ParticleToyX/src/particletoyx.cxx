#include "particles.hxx"

#include <vect.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cstdlib>
#include <vector>

namespace ParticleToyX {
using namespace Arith;
using namespace ParticlesX;

extern "C" void ParticleToyX_init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ParticleToyX_init;
  DECLARE_CCTK_PARAMETERS;

  const vect<CCTK_REAL, dim> xmin{-1, -1, -1};
  const vect<CCTK_REAL, dim> xmax{+1, +1, +1};

  Particles &particles = get_particles();
  const int nparticles = particles.size();

  // This is a random number
  srand(0xdfcd7e49);

  for (int n = 0; n < nparticles; ++n) {
    const vect<CCTK_REAL, dim> r{
        rand() / CCTK_REAL(RAND_MAX),
        rand() / CCTK_REAL(RAND_MAX),
        rand() / CCTK_REAL(RAND_MAX),
    };
    const vect<CCTK_REAL, dim> pos = xmin + r * (xmax - xmin);
    const Particle p{CCTK_REAL(n), 0, pos, {0, 0, 0}, {0, 0, 0}};
    particles.set(p, n);
  }
}

extern "C" void ParticleToyX_accel(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ParticleToyX_accel;
  DECLARE_CCTK_PARAMETERS;

  Particles &particles = get_particles();
  const int nparticles = particles.size();

  for (int n = 0; n < nparticles; ++n) {
    Particle p = particles.get(n);
    p.acc = -p.pos;
    particles.set(p, n);
  }
}

extern "C" void ParticleToyX_evolve(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ParticleToyX_evolve;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dt = CCTK_DELTA_TIME;

  Particles &particles = get_particles();
  const int nparticles = particles.size();

  for (int n = 0; n < nparticles; ++n) {
    Particle p = particles.get(n);
    p.tau += dt;
    p.pos += dt * p.vel;
    p.vel += dt * p.acc;
    particles.set(p, n);
  }
}

} // namespace ParticleToyX
