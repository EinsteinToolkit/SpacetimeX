#include "integrate.hxx"
#include "interpolate.hxx"
#include "multipole.hxx"
#include "sphericalharmonic.hxx"
#include "utils.hxx"

#include <assert.h>
#include <iomanip>
#include <sstream>
#include <stdio.h>
#include <string>

#include <loop_device.hxx>

namespace Multipole {
using namespace std;

static vector<CCTK_REAL> xs, ys, zs, xhat, yhat, zhat;
static vector<CCTK_REAL> th, ph;
static vector<CCTK_REAL> realV, imagV;

// Y[sw][l][m][i]
static vector<vector<vector<vector<CCTK_REAL> > > > realY, imagY;

// Parsed variables
static vector<VariableParse> vars;
static vector<int> spin_weights;

static const int max_vars = 10;

// since each variable can have at most one spin weight, max_vars is a good
// upper limit for the number of spin weights to expect
static const int max_spin_weights = max_vars;

// modes are 0...max_l_modes-1
static const int max_l_modes = 10;
static const int max_m_modes = 2 * max_l_modes + 1;

static void fillVariable(int idx, const char *optString, void *callbackArg) {
  assert(idx >= 0);
  assert(callbackArg != nullptr);

  vector<VariableParse> *vs = static_cast<vector<VariableParse> *>(callbackArg);

  assert(vs->size() < max_vars); // Ensure we don't exceed max_vars
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
      assert(spin_weights.size() < max_spin_weights);
      spin_weights.push_back(vars[i].spinWeight);
    }
  }
}

static void
setup_harmonics(const vector<int> spin_weights, int lmax, vector<CCTK_REAL> th,
                vector<CCTK_REAL> ph, int array_size,
                vector<vector<vector<vector<CCTK_REAL> > > > &realY,
                vector<vector<vector<vector<CCTK_REAL> > > > &imagY) {
  for (size_t si = 0; si < spin_weights.size(); si++) {
    int sw = spin_weights[si];
    vector<vector<vector<CCTK_REAL> > > realY_s, imagY_s;
    for (int l = 0; l <= lmax; l++) {
      vector<vector<CCTK_REAL> > realY_s_l, imagY_s_l;
      for (int m = -l; m <= l; m++) {
        vector<CCTK_REAL> realY_s_l_m(array_size), imagY_s_l_m(array_size);
        HarmonicSetup(sw, l, m, array_size, th, ph, realY_s_l_m, imagY_s_l_m);
        realY_s_l.push_back(realY_s_l_m);
        imagY_s_l.push_back(imagY_s_l_m);
      }
      realY_s.push_back(realY_s_l);
      imagY_s.push_back(imagY_s_l);
    }
    realY.push_back(realY_s);
    imagY.push_back(imagY_s);
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
                      CCTK_REAL *th, CCTK_REAL *ph, CCTK_REAL *real,
                      CCTK_REAL *imag, int array_size) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_MyProc(cctkGH) == 0 && output_tsv) {
    if (out_1d_every != 0 && cctk_iteration % out_1d_every == 0) {
      ostringstream real_base;
      real_base << "mp_" << string(CCTK_VarName(v->index)) << "_r"
                << setiosflags(ios::fixed) << setprecision(2) << rad;
      Output1D(CCTK_PASS_CTOC, real_base.str() + string(".th.tsv"), array_size,
               th, ph, mp_theta, real);
      Output1D(CCTK_PASS_CTOC, real_base.str() + string(".ph.tsv"), array_size,
               th, ph, mp_phi, real);

      if (v->imagIndex != -1) {
        ostringstream imag_base;
        imag_base << "mp_" << string(CCTK_VarName(v->imagIndex)) << "_r"
                  << setiosflags(ios::fixed) << setprecision(2) << rad;
        Output1D(CCTK_PASS_CTOC, imag_base.str() + string(".th.tsv"),
                 array_size, th, ph, mp_theta, imag);
        Output1D(CCTK_PASS_CTOC, imag_base.str() + string(".ph.tsv"),
                 array_size, th, ph, mp_phi, imag);
      }
    }
  }
}

extern "C" void Multipole_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  const int array_size = (ntheta + 1) * (nphi + 1);

  realV.resize(array_size);
  imagV.resize(array_size);
  th.resize(array_size);
  ph.resize(array_size);
  xs.resize(array_size);
  ys.resize(array_size);
  zs.resize(array_size);
  xhat.resize(array_size);
  yhat.resize(array_size);
  zhat.resize(array_size);

  parse_variables_string(string(variables), vars);
  get_spin_weights(vars, spin_weights);
  CoordSetup(xhat.data(), yhat.data(), zhat.data(), th.data(), ph.data());
  setup_harmonics(spin_weights, l_max, th, ph, array_size, realY, imagY);
  CCTK_VINFO("initialized arrays");
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

  const int array_size = (ntheta + 1) * (nphi + 1);
  const int n_variables = vars.size();

  if (out_every == 0 || cctk_iteration % out_every != 0)
    return;

  int lmax = l_max;

  assert(lmax + 1 <= max_l_modes);

  ModeArray modes(n_variables, nradii, lmax);

  for (int v = 0; v < n_variables; v++) {
    // assert(vars[v].spinWeight == -2);
    int si = find_int_in_array(vars[v].spinWeight, spin_weights.data(),
                               spin_weights.size());
    assert(si != -1);

    for (int i = 0; i < nradii; i++) {
      // Compute x^i = r * \hat x^i
      ScaleCartesian(ntheta, nphi, radius[i], xhat.data(), yhat.data(),
                     zhat.data(), xs.data(), ys.data(), zs.data());

      // Interpolate Psi4r and Psi4i
      Interp(CCTK_PASS_CTOC, xs.data(), ys.data(), zs.data(), vars[v].index,
             vars[v].imagIndex, realV, imagV);

      for (int l = 0; l <= lmax; l++) {
        for (int m = -l; m <= l; m++) {
          // Integrate sYlm (realV + i imagV) over the sphere at radius r
          Integrate(array_size, ntheta, realY[si][l][m + l],
                    imagY[si][l][m + l], realV, imagV, th.data(), ph.data(),
                    &modes(v, i, l, m, 0), &modes(v, i, l, m, 1));

        } // loop over m
      } // loop over l
      output_1d(CCTK_PASS_CTOC, &vars[v], radius[i], th.data(), ph.data(),
                realV.data(), imagV.data(), array_size);
    } // loop over radii
  } // loop over variables

  output_modes(CCTK_PASS_CTOC, vars.data(), radius, modes);
}

} // namespace Multipole
