#include "puncture.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>

#include <mpi.h>

#include <cassert>
#include <cmath>
#include <cstdio>
#include <array>
#include <ctype.h>

namespace PunctureTracker {

static PunctureContainer *g_punctures = nullptr;

const int max_num_tracked = 10;

extern "C" void PunctureTracker_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_PunctureTracker_Setup;
  DECLARE_CCTK_PARAMETERS;

  if (verbose) {
    CCTK_INFO("Initializing PunctureTracker");
  }

  for (int n = 0; n < max_num_tracked; ++n) {
    if (track[n]) {
      pt_loc_t[n] = cctk_time;
      pt_loc_x[n] = initial_x[n];
      pt_loc_y[n] = initial_y[n];
      pt_loc_z[n] = initial_z[n];
      pt_vel_t[n] = cctk_time;
      pt_vel_x[n] = 0.0;
      pt_vel_y[n] = 0.0;
      pt_vel_z[n] = 0.0;
    } else {
      pt_loc_t[n] = 0.0;
      pt_loc_x[n] = 0.0;
      pt_loc_y[n] = 0.0;
      pt_loc_z[n] = 0.0;
      pt_vel_t[n] = 0.0;
      pt_vel_x[n] = 0.0;
      pt_vel_y[n] = 0.0;
      pt_vel_z[n] = 0.0;
    }
  }

  // Initialize PunctureContainer
  if (g_punctures == nullptr) {
    g_punctures = new PunctureContainer();
  }

  for (int n = 0; n < max_num_tracked; ++n) {
    if (track[n]) {
      g_punctures->getTime().push_back(cctk_time);
      g_punctures->getLocation()[0].push_back(initial_x[n]);
      g_punctures->getLocation()[1].push_back(initial_y[n]);
      g_punctures->getLocation()[2].push_back(initial_z[n]);
      g_punctures->getVelocity()[0].push_back(0.0);
      g_punctures->getVelocity()[1].push_back(0.0);
      g_punctures->getVelocity()[2].push_back(0.0);
    }
  }
  const int nPunctures = g_punctures->getLocation()[0].size();
  g_punctures->getPreviousTime().resize(nPunctures);
  g_punctures->getBeta()[0].resize(nPunctures);
  g_punctures->getBeta()[1].resize(nPunctures);
  g_punctures->getBeta()[2].resize(nPunctures);
  g_punctures->setNumPunctures();
  assert(g_punctures->getNumPunctures() == nPunctures);

  // enabled if refinement regions should follow the punctures
  if (track_boxes) {
    const std::array<std::vector<CCTK_REAL>, Loop::dim> &location =
        g_punctures->getLocation();
    for (size_t i = 0; i < location[0].size(); ++i) {
      CCTK_VINFO("Writing punc coords to box %zu.", i);
      position_x[i] = location[0][i];
      position_y[i] = location[1][i];
      position_z[i] = location[2][i];
    }
  }
}

extern "C" void PunctureTracker_Finalize(CCTK_ARGUMENTS) {
  delete g_punctures;
  g_punctures = nullptr;
}

extern "C" void PunctureTracker_Track(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_PunctureTracker_Track;
  DECLARE_CCTK_PARAMETERS;

  // Do not track while setting up initial data; time interpolation may fail
  if (cctk_iteration == 0) {
    return;
  }

  // Some output
  if (verbose) {
    CCTK_INFO("Tracking punctures...");
  }

  const std::array<std::vector<CCTK_REAL>, Loop::dim> &location =
      g_punctures->getLocation();

  if (verbose) {
    for (size_t n = 0; n < location[0].size(); ++n) {
      if (track[n]) {
        CCTK_VINFO("Puncture #%zu is at (%g,%g,%g)", n, double(location[0][n]),
                   double(location[1][n]), double(location[2][n]));
      }
    }
  }

  // Manual time level cycling
  g_punctures->updatePreviousTime(CCTK_PASS_CTOC);

  // Interpolate
  g_punctures->interpolate(CCTK_PASS_CTOC);

  if (CCTK_MyProc(cctkGH) == 0) {
    const std::array<std::vector<CCTK_REAL>, Loop::dim> &beta =
        g_punctures->getBeta();

    // More output
    if (verbose) {
      for (size_t n = 0; n < beta[0].size(); ++n) {
        CCTK_VINFO("Shift at puncture #%zu is at (%g,%g,%g)", n,
                   double(beta[0][n]), double(beta[1][n]), double(beta[2][n]));
      }
    }

    // Check for NaNs and large shift components
    for (size_t n = 0; n < beta[0].size(); ++n) {
      CCTK_REAL norm =
          sqrt(pow(beta[0][n], 2) + pow(beta[1][n], 2) + pow(beta[2][n], 2));

      if (!CCTK_isfinite(norm) || norm > shift_limit) {
        CCTK_VERROR("Shift at puncture #%zu is (%g,%g,%g).  This likely "
                    "indicates an error in the simulation.",
                    n, double(beta[0][n]), double(beta[1][n]),
                    double(beta[2][n]));
      }
    }
  }

  // Time evolution
  g_punctures->evolve(CCTK_PASS_CTOC);

  // Broadcast result: 3 components for location, 3 components for velocity
  g_punctures->broadcast(CCTK_PASS_CTOC);

  if (track_boxes) {
    for (size_t i = 0; i < location[0].size(); ++i) {
      position_x[i] = location[0][i];
      position_y[i] = location[1][i];
      position_z[i] = location[2][i];
    }
  }

  // Write to pt_loc_foo and pt_vel_foo
  const std::array<std::vector<CCTK_REAL>, Loop::dim> &velocity =
      g_punctures->getVelocity();
  for (size_t i = 0; i < location[0].size(); ++i) {
    pt_loc_x[i] = location[0][i];
    pt_loc_y[i] = location[1][i];
    pt_loc_z[i] = location[2][i];
    pt_vel_x[i] = velocity[0][i];
    pt_vel_y[i] = velocity[1][i];
    pt_vel_z[i] = velocity[2][i];
  }
}

extern "C" void PunctureTracker_CheckShift(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_PunctureTracker_CheckShift;
  DECLARE_CCTK_PARAMETERS;

  const int level = ilogb(CCTK_REAL(cctk_levfac[0]));
  const int finest_lvl = num_levels[0] - 1;

  if (level == finest_lvl) {
    for (int n = 0; n < max_num_tracked; ++n) {
      if (track[n]) {
        const Arith::vect<CCTK_REAL, Loop::dim> loc_vec = {
            pt_loc_x[n], pt_loc_y[n], pt_loc_z[n]};

        grid.loop_all_device<0, 0, 0>(
            grid.nghostzones,
            [=] CCTK_DEVICE(
                const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
              if (maximum(abs(p.X - loc_vec)) <= (1 / pow(2, level))) {
                printf("Shift at level %d near puncture #%d is {%g, %g, %g} at "
                       "coords {%g, %g, %g}.",
                       level, n, betax(p.I), betay(p.I), betaz(p.I), p.x, p.y,
                       p.z);
              }
            });
      }
    }
  }
}

} // namespace PunctureTracker
