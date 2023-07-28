#include <cmath>
#include <ctype.h>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>

#include <loop_device.hxx>
#include <array>

namespace TestPT {
using namespace std;
using namespace Loop;

extern "C" void TestPT_init_shift_circle(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestPT_init_shift_circle;
  DECLARE_CCTK_PARAMETERS;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout(cctkGH, indextype);

  const GF3D2<CCTK_REAL> betax_(layout, betax);
  const GF3D2<CCTK_REAL> betay_(layout, betay);
  const GF3D2<CCTK_REAL> betaz_(layout, betaz);

  const GridDescBaseDevice grid(cctkGH);
  grid.loop_all_device<0, 0, 0>(grid.nghostzones,
                                [=] CCTK_DEVICE(const PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      betax_(p.I) = p.y;
                                      betay_(p.I) = -p.x;
                                      betaz_(p.I) = p.z;
                                    });
}

extern "C" void TestPT_init_shift_ycircle(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestPT_init_shift_ycircle;
  DECLARE_CCTK_PARAMETERS;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout(cctkGH, indextype);

  const GF3D2<CCTK_REAL> betax_(layout, betax);
  const GF3D2<CCTK_REAL> betay_(layout, betay);
  const GF3D2<CCTK_REAL> betaz_(layout, betaz);

  const GridDescBaseDevice grid(cctkGH);
  grid.loop_all_device<0, 0, 0>(grid.nghostzones,
                                [=] CCTK_DEVICE(const PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      betax_(p.I) = p.z;
                                      betay_(p.I) = 0;
                                      betaz_(p.I) = -p.x;
                                    });
}

extern "C" void TestPT_init_shift_3D(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestPT_init_shift_3D;
  DECLARE_CCTK_PARAMETERS;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout(cctkGH, indextype);

  const GF3D2<CCTK_REAL> betax_(layout, betax);
  const GF3D2<CCTK_REAL> betay_(layout, betay);
  const GF3D2<CCTK_REAL> betaz_(layout, betaz);

  const GridDescBaseDevice grid(cctkGH);
  grid.loop_all_device<0, 0, 0>(grid.nghostzones,
                                [=] CCTK_DEVICE(const PointDesc &p)
                                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                      betax_(p.I) = (1.0 / 50.0) * p.x * p.x + p.y * p.y * p.y;
                                      betay_(p.I) = -(1.0 / 50.0) * p.x * p.x * p.x + p.y * p.y;
                                      betaz_(p.I) = p.z * (p.x * p.x + p.y * p.y);
                                    });
}

extern "C" void TestPT_interp_test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestPT_interp_test;
  DECLARE_CCTK_PARAMETERS;
	
  // Interpolate

  // Dimensions
  const int dim = 3;

	// Number of interpolation variables
	int const num_vars = 3;

  const int operator_handle = 0;

  // Interpolation parameter table
  CCTK_INT operations[1][dim];
	for (int var = 0 ; var < num_vars; var++) {
	  operations[0][var] = 0;
	}

  int operands[1][dim];
	for (int var = 0 ; var < num_vars; var++) {
		operands[0][var] = var;
	}

	int ierr;
  // const int param_table_handle = Util_TableCreateFromString("order=4");
  int param_table_handle = Util_TableCreate(UTIL_TABLE_FLAGS_DEFAULT);
  if (param_table_handle < 0)
    CCTK_VERROR("Can't create parameter table: %d", param_table_handle);
  if ((ierr = Util_TableSetInt(param_table_handle, 1, "order")) < 0)
    CCTK_VERROR("Can't set order in parameter table: %d", ierr);
  if ((ierr = Util_TableSetIntArray(param_table_handle, num_vars, (int const*const)operands,
                            "operand_indices")) < 0)
    CCTK_VERROR("Can't set operand_indices array in parameter table: %d", ierr);
  if ((ierr = Util_TableSetIntArray(param_table_handle, num_vars, (int const*const)operations,
                            "operation_codes")) < 0)
    CCTK_VERROR("Can't set operation_codes array in parameter table: %d", ierr);

  {

    // Interpolation coordinate system: Not used in CarpetX_DriverInterpolate
    const int coordsys_handle = 0;
		CCTK_INT const interp_coords_type_code = 0;

    // Only processor 0 interpolates
    const int num_points = CCTK_MyProc(cctkGH) == 0 ? 2 : 0;

    // Interpolation coordinates
    assert(dim == 3);
    const CCTK_REAL coords[dim][2] = {
			{4, 4}, {0, 0}, {-1, 1}
		};
		const void* interp_coords[dim] = {
			coords[0], coords[1], coords[2]
		};


    // Interpolated variables
    assert(num_vars == 3);
    int input_array_indices[3];
    input_array_indices[0] = CCTK_VarIndex("ADMBase::betax");
    input_array_indices[1] = CCTK_VarIndex("ADMBase::betay");
    input_array_indices[2] = CCTK_VarIndex("ADMBase::betaz");

    // Interpolation result types: Not used by CarpetX DriverInterp
		CCTK_INT const output_array_type_codes[1] = {0};

    // Interpolation result
    CCTK_REAL pt_betax[2];
    CCTK_REAL pt_betay[2];
    CCTK_REAL pt_betaz[2];

    assert(num_vars == 3);
    CCTK_POINTER output_arrays[3];
    output_arrays[0] = pt_betax;
    output_arrays[1] = pt_betay;
    output_arrays[2] = pt_betaz;

    // Interpolate
    int ierr;
    // Use CarpetX Funtion:
		ierr = DriverInterpolate(
		cctkGH, dim, operator_handle, param_table_handle, coordsys_handle,
		num_points, interp_coords_type_code, interp_coords, num_vars, (int const * const)input_array_indices,
		num_vars, output_array_type_codes, output_arrays);

    if (CCTK_MyProc(cctkGH) == 0) {

			CCTK_VINFO("Shift at z=-1 interpolated to be (%g,%g,%g)", 
								 double(pt_betax[0]), double(pt_betay[0]),
								 double(pt_betaz[0]));
			CCTK_VINFO("Shift at z=+1 interpolated to be (%g,%g,%g)", 
								 double(pt_betax[1]), double(pt_betay[1]),
								 double(pt_betaz[1]));
		}
	}
}
	
} //namespace TestPT
