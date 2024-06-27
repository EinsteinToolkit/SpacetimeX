#include "io.hxx"
#include "multipole.hxx"
#include "surface.hxx"

#include <assert.h>
#include <iomanip>
#include <sstream>
#include <stdio.h>
#include <string>

#include <loop_device.hxx>

#include <util_Table.h>

namespace Multipole {
using namespace std;

static Sphere *g_sphere = nullptr;

// Parsed variables
static vector<VariableParse> g_vars;
static vector<int> g_spin_weights;

// Function to output modes
static void outputModes(CCTK_ARGUMENTS, const VariableParse vars[],
                        const CCTK_REAL radii[], const ModeArray &modes) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (output_tsv && CCTK_MyProc(cctkGH) == 0) {
    for (int v = 0; v < modes.getNumVars(); ++v) {
      for (int i = 0; i < modes.getNumRadii(); ++i) {
        const CCTK_REAL rad = radii[i];
        for (int l = 0; l <= modes.getMaxL(); ++l) {
          for (int m = -l; m <= l; ++m) {
            std::ostringstream filename;
            filename << "mp_" << vars[v].name << "_l" << l << "_m" << m << "_r"
                     << std::fixed << std::setprecision(2) << rad << ".tsv";

            OutputComplexToFile(CCTK_PASS_CTOC, filename.str(),
                                modes(v, i, l, m, 0), modes(v, i, l, m, 1));
          }
        }
      }
    }
  }
  if (output_hdf5 && CCTK_MyProc(cctkGH) == 0) {
    OutputComplexToH5File(CCTK_PASS_CTOC, vars, radii, modes);
  }
}

static void fillVariable(int idx, const char *optString, void *callbackArg) {
  assert(idx >= 0);
  assert(callbackArg != nullptr);

  vector<VariableParse> *vars =
      static_cast<vector<VariableParse> *>(callbackArg);

  VariableParse var;
  var.realIndex = idx;

  // Initialize default values
  var.imagIndex = -1;
  var.spinWeight = 0;
  var.name = string(CCTK_VarName(var.realIndex));

  if (optString != nullptr) {
    int table = Util_TableCreateFromString(optString);

    if (table >= 0) {
      const int bufferLength = 256;
      char buffer[bufferLength];

      Util_TableGetInt(table, &var.spinWeight, "sw");
      if (Util_TableGetString(table, bufferLength, buffer, "cmplx") >= 0) {
        var.imagIndex = CCTK_VarIndex(buffer);
      }
      if (Util_TableGetString(table, bufferLength, buffer, "name") >= 0) {
        var.name = string(buffer);
      }

      const int ierr = Util_TableDestroy(table);
      if (ierr) {
        CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                    "Could not destroy table: %d", ierr);
      }
    }
  }

  vars->push_back(var);
}

extern "C" void Multipole_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  // Parse parameters and save them to g_vars
  int ierr = CCTK_TraverseString(string(variables).c_str(), fillVariable,
                                 &g_vars, CCTK_GROUP_OR_VAR);
  assert(ierr >= 0);

  // Get all different kinds of spin weights
  for (size_t i = 0; i < g_vars.size(); i++) {
    if (!int_in_array(g_vars[i].spinWeight, g_spin_weights)) {
      g_spin_weights.push_back(g_vars[i].spinWeight);
    }
  }

  // Initialize a sphere (whose radius can be changed later)
  if (g_sphere == nullptr) {
    g_sphere =
        new Sphere(ntheta, nphi, CCTK_Equals(integration_method, "midpoint"),
                   g_spin_weights, l_max);
  }

  CCTK_VINFO("initialized arrays");
}

extern "C" void Multipole_Finalize(CCTK_ARGUMENTS) {
  delete g_sphere;
  g_sphere = nullptr;
}

// This is the main scheduling file.  Because we are completely local here
// and do not use cactus arrays etc, we schedule only one function and then
// like program like one would in C, C++ with this function taking the
// place of int main(void).
// This function calls functions to accomplish 3 things:
//   1) Interpolate psi4 onto a sphere
//   2) Integrate psi4 with the ylm's over that sphere
//   3) Output the mode decomposed psi4
extern "C" void Multipole_Calc(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Multipole_Calc;
  DECLARE_CCTK_PARAMETERS;

  const int n_variables = g_vars.size();

  if (out_every == 0 || cctk_iteration % out_every != 0)
    return;

  ModeArray modes(n_variables, nradii, l_max);

  for (int v = 0; v < n_variables; v++) {
    int si = find_int_in_array(g_vars[v].spinWeight, g_spin_weights.data(),
                               g_spin_weights.size());
    assert(si != -1);

    for (int i = 0; i < nradii; i++) {

      // Set radius of the sphere
      g_sphere->setRadius(radius[i]);

      // Interpolate to the sphere
      g_sphere->interpolate(CCTK_PASS_CTOC, g_vars[v].realIndex,
                            g_vars[v].imagIndex);

      // Intergate of conj(sYlm)*F*sin(theta) over the sphere at radiusr[i]
      for (int l = 0; l <= l_max; l++) {
        for (int m = -l; m <= l; m++) {
          g_sphere->integrate(g_sphere->getRealY()[si][l][m + l],
                              g_sphere->getImagY()[si][l][m + l],
                              g_sphere->getRealF(), g_sphere->getImagF(),
                              &modes(v, i, l, m, 0), &modes(v, i, l, m, 1));

        } // loop over m
      } // loop over l

      g_sphere->output1D(CCTK_PASS_CTOC, g_vars[v], radius[i]);

    } // loop over radii
  } // loop over variables

  outputModes(CCTK_PASS_CTOC, g_vars.data(), radius, modes);
}

} // namespace Multipole
