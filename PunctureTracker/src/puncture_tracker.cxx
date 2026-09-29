#include "puncture.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <ctype.h>

namespace PunctureTracker {

static PunctureContainer *g_punctures = nullptr;

const int max_num_tracked = 10;

// `BoxInBox::positions` is a vector grid scalar with one element per
// refinement region, and both `PunctureTracker_Setup` and
// `PunctureTracker_Track` declare `WRITES: BoxInBox::positions`. A WRITES
// clause names the whole vector group, so CarpetX poisons every region
// immediately before the routine and checks every region immediately after it.
// Writing only the tracked punctures leaves the remaining regions as nans
// and aborts in `valid.cxx`.  Write all of them: the tracked punctures where
// there is one, and BoxInBox's own parameters (i.e. what `BoxInBox_Init` put
// there) everywhere else. This also caps the write at the number of regions
// that actually exist, which the bare `n < nPunctures` loops did not do.

namespace {

const int max_num_boxes = 3; // BoxInBox::max_num_regions

int get_num_boxes() {
  const int gi = CCTK_GroupIndex("BoxInBox::positions");
  if (gi < 0)
    CCTK_VERROR("Could not find group BoxInBox::positions");
  cGroup gdata;
  const int ierr = CCTK_GroupData(gi, &gdata);
  if (ierr != 0)
    CCTK_VERROR("Could not query group BoxInBox::positions");
  if (gdata.vectorlength != max_num_boxes)
    CCTK_VERROR("BoxInBox::positions has %d regions, but PunctureTracker knows "
                "how to restore %d of them",
                gdata.vectorlength, max_num_boxes);
  return gdata.vectorlength;
}

// BoxInBox's `position_[xyz]_N` are private parameters, so they cannot be
// pulled in with `SHARES: BoxInBox` / `USES`. We read them at runtime instead.
CCTK_REAL get_box_position_param(const char *const component, const int box) {
  char name[32];
  snprintf(name, sizeof name, "position_%s_%d", component, box + 1);
  int type = -1;
  const void *const value = CCTK_ParameterGet(name, "BoxInBox", &type);
  if (!value)
    CCTK_VERROR("Could not read parameter BoxInBox::%s", name);
  if (type != PARAMETER_REAL)
    CCTK_VERROR("Parameter BoxInBox::%s is not a real", name);
  return *static_cast<const CCTK_REAL *>(value);
}

void set_box_positions(CCTK_REAL *restrict const position_x,
                       CCTK_REAL *restrict const position_y,
                       CCTK_REAL *restrict const position_z,
                       const int num_tracked_boxes,
                       const std::array<std::vector<CCTK_REAL>, Loop::dim>
                           &location) {
  const int num_boxes = get_num_boxes();
  for (int n = 0; n < num_boxes; ++n) {
    if (n < num_tracked_boxes) {
      position_x[n] = location[0][n];
      position_y[n] = location[1][n];
      position_z[n] = location[2][n];
    } else {
      position_x[n] = get_box_position_param("x", n);
      position_y[n] = get_box_position_param("y", n);
      position_z[n] = get_box_position_param("z", n);
    }
  }
}

} // namespace

extern "C" void PunctureTracker_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_PunctureTracker_Init;
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
}

extern "C" void PunctureTracker_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_PunctureTracker_Setup;
  DECLARE_CCTK_PARAMETERS;

  // Initialize PunctureContainer
  if (g_punctures == nullptr) {
    g_punctures = new PunctureContainer();

    for (int n = 0; n < max_num_tracked; ++n) {
      if (track[n]) {
        g_punctures->getTime().push_back(pt_loc_t[n]);
        g_punctures->getLocation()[0].push_back(pt_loc_x[n]);
        g_punctures->getLocation()[1].push_back(pt_loc_y[n]);
        g_punctures->getLocation()[2].push_back(pt_loc_z[n]);
        g_punctures->getVelocity()[0].push_back(pt_vel_x[n]);
        g_punctures->getVelocity()[1].push_back(pt_vel_y[n]);
        g_punctures->getVelocity()[2].push_back(pt_vel_z[n]);
      }
    }
  }

  const int nPunctures = g_punctures->getTime().size();
  g_punctures->getPreviousTime().resize(nPunctures);
  g_punctures->getBeta()[0].resize(nPunctures);
  g_punctures->getBeta()[1].resize(nPunctures);
  g_punctures->getBeta()[2].resize(nPunctures);
  g_punctures->setNumPunctures();
  assert(g_punctures->getNumPunctures() == nPunctures);

  // enabled if refinement regions should follow the punctures. The regions
  // that do not follow a puncture are still written, see `set_box_positions`
  const std::array<std::vector<CCTK_REAL>, Loop::dim> &location =
      g_punctures->getLocation();
  const int num_tracked_boxes =
      track_boxes ? std::min(nPunctures, max_num_boxes) : 0;
  if (verbose)
    for (int n = 0; n < num_tracked_boxes; ++n)
      CCTK_VINFO("Writing punc coords to box %d.", n);
  set_box_positions(position_x, position_y, position_z, num_tracked_boxes,
                    location);
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

  const int nPunctures = g_punctures->getNumPunctures();
  const std::array<std::vector<CCTK_REAL>, Loop::dim> &location =
      g_punctures->getLocation();
  const std::array<std::vector<CCTK_REAL>, Loop::dim> &velocity =
      g_punctures->getVelocity();
  const std::vector<CCTK_REAL> &time = g_punctures->getTime();

  if (verbose) {
    for (int n = 0; n < nPunctures; ++n) {
      CCTK_VINFO("Puncture #%d is at (%g,%g,%g)", n, double(location[0][n]),
                 double(location[1][n]), double(location[2][n]));
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
      for (int n = 0; n < nPunctures; ++n) {
        CCTK_VINFO("Shift at puncture #%d is at (%g,%g,%g)", n,
                   double(beta[0][n]), double(beta[1][n]), double(beta[2][n]));
      }
    }

    // Check for NaNs and large shift components
    for (int n = 0; n < nPunctures; ++n) {
      CCTK_REAL norm =
          sqrt(pow(beta[0][n], 2) + pow(beta[1][n], 2) + pow(beta[2][n], 2));

      if (!CCTK_isfinite(norm) || norm > shift_limit) {
        CCTK_VERROR("Shift at puncture #%d is (%g,%g,%g).  This likely "
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

  // Write to pt_loc_foo and pt_vel_foo. `pt_loc` and `pt_vel` are
  // max_num_tracked-element vector groups and this routine declares both as
  // WRITES, so the driver poisons all of them before the call and checks all
  // of them after it. The slots beyond the tracked punctures have to be given
  // the same zeros PunctureTracker_Init writes for an untracked puncture.
  for (int i = 0; i < max_num_tracked; ++i) {
    const bool tracked = i < nPunctures;
    pt_loc_t[i] = tracked ? time[i] : 0.0;
    pt_loc_x[i] = tracked ? location[0][i] : 0.0;
    pt_loc_y[i] = tracked ? location[1][i] : 0.0;
    pt_loc_z[i] = tracked ? location[2][i] : 0.0;
    pt_vel_t[i] = tracked ? time[i] : 0.0;
    pt_vel_x[i] = tracked ? velocity[0][i] : 0.0;
    pt_vel_y[i] = tracked ? velocity[1][i] : 0.0;
    pt_vel_z[i] = tracked ? velocity[2][i] : 0.0;
  }

  set_box_positions(position_x, position_y, position_z,
                    track_boxes ? std::min(nPunctures, max_num_boxes) : 0,
                    location);
}

} // namespace PunctureTracker
