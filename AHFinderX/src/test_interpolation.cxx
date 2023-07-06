#include <cctk.h>
#include <cctk_Arguments.h>
#include <util_Table.h>

#include <array>
#include <cassert>
#include <vector>

namespace AHFinder {
using namespace std;

extern "C" void AHFinder_test_interpolation(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AHFinder_test_interpolation;

  const vector<CCTK_INT> all_operations{0, 1, 2, 3, 11, 12, 13, 22, 23, 33};

  const vector<CCTK_INT> all_varinds{
      CCTK_VarIndex("Coordinates::vcoordx"),
      CCTK_VarIndex("Coordinates::vcoordy"),
      CCTK_VarIndex("Coordinates::vcoordz"),
  };

  const int nvars = all_varinds.size() * all_operations.size();

  vector<CCTK_INT> varinds;
  vector<CCTK_INT> operations;
  varinds.reserve(nvars);
  operations.reserve(nvars);
  for (const auto op : all_operations) {
    for (const auto var : all_varinds) {
      varinds.push_back(var);
      operations.push_back(op);
    }
  }
  assert(int(varinds.size()) == nvars);
  assert(int(operations.size()) == nvars);

  constexpr int npoints = 4;
  const array<vector<CCTK_REAL>, 3> coords{
      vector<CCTK_REAL>{0.1, 0.2, 0.3, 0.1},
      vector<CCTK_REAL>{0.1, 0.2, 0.3, 0.2},
      vector<CCTK_REAL>{0.1, 0.2, 0.3, 0.3},
  };
  for (int d = 0; d < 3; ++d)
    assert(coords[d].size() == npoints);

  vector<vector<CCTK_REAL> > results(nvars);
  vector<vector<CCTK_REAL> > results2(nvars);
  vector<CCTK_REAL *> resultptrs(nvars);
  vector<CCTK_REAL *> resultptrs2(nvars);
  for (int n = 0; n < nvars; ++n) {
    results.at(n).resize(npoints);
    resultptrs.at(n) = results.at(n).data();
    results2.at(n).resize(npoints);
    resultptrs2.at(n) = results2.at(n).data();
  }

  Interpolate(cctkGH, npoints, coords[0].data(), coords[1].data(),
              coords[2].data(), nvars, varinds.data(), operations.data(),
              resultptrs.data());

  CCTK_INT const N_dims = 3;

  const void *interp_coords[N_dims];
  interp_coords[0] = coords[0].data();
  interp_coords[1] = coords[1].data();
  interp_coords[2] = coords[2].data();

  void *const *output_array = (void *const *)&resultptrs2[0];

  /* DriverInterpolate arguments that aren't currently used */
  CCTK_INT const coord_system_handle = 0;
  CCTK_INT const interp_coords_type_code = 0;
  CCTK_INT const output_array_type_codes[1] = {0};
  int interp_handle = 0;

  /* Table generation */

  int param_table_handle;
  int operands[nvars];
  for (int i = 0; i < nvars; i++) {
    operands[i] = i;
  }

  param_table_handle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);
  if (param_table_handle < 0)
    CCTK_ERROR("Can't create parameter table!");
  if (Util_TableSetInt(param_table_handle, 1, "order") < 0)
    CCTK_ERROR("Can't set order in parameter table!");
  if (Util_TableSetIntArray(param_table_handle, nvars, operands,
                            "operand_indices") < 0)
    CCTK_ERROR("Can't set operand_indices array in parameter table!");
  if (Util_TableSetIntArray(param_table_handle, nvars, operations.data(),
                            "operation_codes") < 0)
    CCTK_ERROR("Can't set operation_codes array in parameter table!");

  int ierr = DriverInterpolate(
      cctkGH, N_dims, interp_handle, param_table_handle, coord_system_handle,
      npoints, interp_coords_type_code, interp_coords, nvars, varinds.data(),
      nvars, output_array_type_codes, output_array);

  const auto chop{[](const auto x) { return fabs(x) <= 1.0e-12 ? 0 : x; }};
  for (int v = 0; v < nvars; ++v) {
    for (int n = 0; n < npoints; ++n) {
      const int d = v % all_varinds.size();
      assert(d >= 0 && d < 3);
      const int op = operations.at(v);
      switch (op) {
      case 0:
        assert(chop(results.at(v).at(n) - coords[d].at(n)) == 0);
        assert(chop(results2.at(v).at(n) - coords[d].at(n)) == 0);
        assert(chop(results.at(v).at(n) - results2.at(v).at(n)) == 0);
        break;
      case 1:
      case 2:
      case 3:
        assert(chop(results.at(v).at(n) - (op - 1 == d ? 1 : 0)) == 0);
        assert(chop(results2.at(v).at(n) - (op - 1 == d ? 1 : 0)) == 0);
        assert(chop(results2.at(v).at(n) - results.at(v).at(n)) == 0);
        break;
      case 11:
      case 12:
      case 13:
      case 22:
      case 23:
      case 33:
        assert(chop(results.at(v).at(n)) == 0);
        assert(chop(results2.at(v).at(n)) == 0);
        break;
      default:
        assert(0);
      }
    }
  }
}

} // namespace AHFinder
