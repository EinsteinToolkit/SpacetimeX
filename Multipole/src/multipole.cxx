#include "multipole.hxx"
#include "surface.hxx"
#include "utils.hxx"

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
static vector<VariableParse> vars;
static vector<int> spin_weights;

static void fillVariable(int idx, const char *optString, void *callbackArg) {
  assert(idx >= 0);
  assert(callbackArg != nullptr);

  vector<VariableParse> *vs = static_cast<vector<VariableParse> *>(callbackArg);
                                 //
  VariableParse v;
  v.index = idx;

  // Initialize default values
  v.imagIndex = -1;
  v.spinWeight = 0;
  v.name = string(CCTK_VarName(v.index));

  if (optString != nullptr) {
    int table = Util_TableCreateFromString(optString);

    if (table >= 0) {
      const int bufferLength = 256;
      char buffer[bufferLength];

      Util_TableGetInt(table, &v.spinWeight, "sw");
      if (Util_TableGetString(table, bufferLength, buffer, "cmplx") >= 0) {
        v.imagIndex = CCTK_VarIndex(buffer);
      }
      if (Util_TableGetString(table, bufferLength, buffer, "name") >= 0) {
        v.name = string(buffer);
      }

      const int ierr = Util_TableDestroy(table);
      if (ierr) {
        CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                    "Could not destroy table: %d", ierr);
      }
    }
  }
  vs->push_back(v);
}

static void parse_variables_string(const string &var_string,
                                   vector<VariableParse> &vs) {
  int ierr = CCTK_TraverseString(var_string.c_str(), fillVariable, &vs,
                                 CCTK_GROUP_OR_VAR);
  assert(ierr >= 0);
}

static void get_spin_weights(vector<VariableParse> vars,
                             vector<int> &spin_weights) {
  for (size_t i = 0; i < vars.size(); i++) {
    if (!int_in_array(vars[i].spinWeight, spin_weights)) {
      spin_weights.push_back(vars[i].spinWeight);
    }
  }
}

static void output_modes(CCTK_ARGUMENTS, const VariableParse vars[],
                         const CCTK_REAL radii[], const ModeArray &modes) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (output_tsv) {
    if (CCTK_MyProc(cctkGH) == 0) {
      for (int v = 0; v < modes.getNumVars(); v++) {
        for (int i = 0; i < modes.getNumRadii(); i++) {
          const CCTK_REAL rad = radii[i];
          for (int l = 0; l <= modes.getMaxL(); l++) {
            for (int m = -l; m <= l; m++) {
              ostringstream name;
              name << "mp_" << vars[v].name << "_l" << l << "_m" << m << "_r"
                   << setiosflags(ios::fixed) << setprecision(2) << rad
                   << ".tsv";
              OutputComplexToFile(CCTK_PASS_CTOC, name.str(),
                                  modes(v, i, l, m, 0), modes(v, i, l, m, 1));
            }
          }
        }
      }
    }
  }
  if (output_hdf5) {
    if (CCTK_MyProc(cctkGH) == 0) {
      OutputComplexToH5File(CCTK_PASS_CTOC, vars, radii, modes);
    }
  }
}

static void output_1d(CCTK_ARGUMENTS, const VariableParse *v, CCTK_REAL rad,
                      const CCTK_REAL *th, const CCTK_REAL *ph,
                      const CCTK_REAL *real, const CCTK_REAL *imag) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_MyProc(cctkGH) == 0 && output_tsv) {
    if (out_1d_every != 0 && cctk_iteration % out_1d_every == 0) {
      ostringstream real_base;
      real_base << "mp_" << string(CCTK_VarName(v->index)) << "_r"
                << setiosflags(ios::fixed) << setprecision(2) << rad;
      Output1D(CCTK_PASS_CTOC, real_base.str() + string(".th.tsv"),
               th, ph, mp_theta, real);
      Output1D(CCTK_PASS_CTOC, real_base.str() + string(".ph.tsv"),
               th, ph, mp_phi, real);

      if (v->imagIndex != -1) {
        ostringstream imag_base;
        imag_base << "mp_" << string(CCTK_VarName(v->imagIndex)) << "_r"
                  << setiosflags(ios::fixed) << setprecision(2) << rad;
        Output1D(CCTK_PASS_CTOC, imag_base.str() + string(".th.tsv"),
                 th, ph, mp_theta, imag);
        Output1D(CCTK_PASS_CTOC, imag_base.str() + string(".ph.tsv"),
                 th, ph, mp_phi, imag);
      }
    }
  }
}

extern "C" void Multipole_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  parse_variables_string(string(variables), vars);
  get_spin_weights(vars, spin_weights);

  if (g_sphere == nullptr) {
    g_sphere =
        new Sphere(ntheta, nphi, CCTK_Equals(integration_method, "midpoint"),
                   spin_weights, l_max);
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

  const int n_variables = vars.size();

  if (out_every == 0 || cctk_iteration % out_every != 0)
    return;

  ModeArray modes(n_variables, nradii, l_max);

  for (int v = 0; v < n_variables; v++) {
    // assert(vars[v].spinWeight == -2);
    int si = find_int_in_array(vars[v].spinWeight, spin_weights.data(),
                               spin_weights.size());
    assert(si != -1);

    for (int i = 0; i < nradii; i++) {

      // Set radius of the sphere
      g_sphere->setRadius(radius[i]);

      // Interpolate to the sphere
      g_sphere->interpolate(CCTK_PASS_CTOC, vars[v].index, vars[v].imagIndex);

      // Intergate of conj(sYlm)*F*sin(theta) over the sphere at radiusr[i]
      for (int l = 0; l <= l_max; l++) {
        for (int m = -l; m <= l; m++) {
          g_sphere->integrate(g_sphere->getRealY()[si][l][m + l],
                              g_sphere->getImagY()[si][l][m + l],
                              g_sphere->getRealF(), g_sphere->getImagF(),
                              &modes(v, i, l, m, 0), &modes(v, i, l, m, 1));

        } // loop over m
      } // loop over l

      output_1d(CCTK_PASS_CTOC, &vars[v], radius[i],
                g_sphere->getTheta().data(), g_sphere->getPhi().data(),
                g_sphere->getRealF().data(), g_sphere->getImagF().data());

    } // loop over radii
  } // loop over variables

  output_modes(CCTK_PASS_CTOC, vars.data(), radius, modes);
}

} // namespace Multipole
