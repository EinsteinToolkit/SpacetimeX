// io.cc -- I/O routines for this thorn
// $Header$
//
// input_gridfn - read an angular grid function from an input file
// input_gridfn__explicit_name - ... with the input file explicitly named
//
// setup_h_files - set up to write horizon-shape and similar output files
/// output_OpenDX_control_file - write an OpenDX control file for an output file
// output_gridfn - write an angular grid function to an output file
// output_Jacobians - write a Jacobian matrix or matrices to an output file
//
/// io_ASCII_file_name - compute file name for angular-gridfn I/O file
/// io_HDF5_file_name - compute file name for angular-gridfn I/O file
//

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "util_Table.h"
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "config.h"
#include "stdc.h"
#include "../jtutil/util.hh"
#include "../jtutil/array.hh"
#include "../jtutil/cpm_map.hh"
#include "../jtutil/linear_map.hh"
using jtutil::error_exit;

#include "../patch/coords.hh"
#include "../patch/grid.hh"
#include "../patch/fd_grid.hh"
#include "../patch/patch.hh"
#include "../patch/patch_edge.hh"
#include "../patch/patch_interp.hh"
#include "../patch/ghost_zone.hh"
#include "../patch/patch_system.hh"

#include "../elliptic/Jacobian.hh"

#include "../gr/gfns.hh"
#include "../gr/gr.hh"

#include "horizon_sequence.hh"
#include "BH_diagnostics.hh"
#include "driver.hh"

// all the code in this file is inside this namespace
namespace AHFinderDirect
	  {

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// prototypes for functions local to this file
//

namespace {
void output_OpenDX_control_file(const patch_system& ps,
				const struct IO_info& IO_info,
				int hn);
const char* io_ASCII_file_name(const struct IO_info& IO_info,
                               const char base_file_name[],
                               int min_digits,
                               int hn, int AHF_iteration = 0);
const char* io_HDF5_file_name(const struct IO_info& IO_info,
                              const char base_file_name[],
                              int min_digits,
                              int hn, int AHF_iteration = 0);
void create_h_directory(const struct IO_info& IO_info);
	  }

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function inputs a gridfn from a data file, with the file name
// chosen automagically.
//
// We assume that this gridfn is ghosted, but the ghost zones are *not*
// present in the data file.
//
void input_gridfn(patch_system& ps, int unknown_gfn,
		  const struct IO_info& IO_info, const char base_file_name[],
		  int min_digits,
		  int hn, bool print_msg_flag, int AHF_iteration /* = 0 */)
{
const char* file_name = io_ASCII_file_name(IO_info, base_file_name, min_digits,
                                           hn, AHF_iteration);

input_gridfn__explicit_name(ps, unknown_gfn,
			    IO_info, base_file_name, print_msg_flag);
}

//******************************************************************************

//
// This function inputs a gridfn from a data file, with the file name
// specified explicitly.
//
// We assume that this gridfn is ghosted, but the ghost zones are *not*
// present in the data file.
//
void input_gridfn__explicit_name(patch_system& ps, int unknown_gfn,
				 const struct IO_info& IO_info,
				 const char file_name[], bool print_msg_flag)
{
if (print_msg_flag)
   then {
	if (unknown_gfn == gfns::gfn__h)
	   then CCTK_VInfo(CCTK_THORNSTRING,
			   "   reading initial guess from \"%s\"", file_name);
	}

if (IO_info.output_ASCII_files)
   then	{
	ps.read_ghosted_gridfn(unknown_gfn,
			       file_name,
			       false);		// no ghost zones in data file
	}

else if (IO_info.output_HDF5_files)
   then {
	CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
"input_gridfn__explicit_name(): reading from HDF5 data files not implemented yet!");
	}

else
	{
	CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   input_gridfn__explicit_name():\n"
"        (this should never happen!)");				/*NOTREACHED*/
	}
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function sets up to output horizon-shape or similar data files:
// - it creates the output directory if this doesn't already exist
// - it writes the OpenDX control files if this is desired
//
void setup_h_files(patch_system& ps, const struct IO_info& IO_info, int hn)
{
// create the output directory (if it doesn't already exist)
create_h_directory(IO_info);
output_OpenDX_control_file(ps, IO_info, hn);
}

//******************************************************************************

//
// This function creates the h ooutput directory if it does not exist
//
namespace {
void create_h_directory(const struct IO_info& IO_info)
{
// create the output directory (if it doesn't already exist)
const int status = CCTK_CreateDirectory(IO_info.default_directory_permission,
					IO_info.h_directory);
if (status < 0)
   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
"\n"
"   setup_h_files():\n"
"        error %d trying to create output directory\n"
"        \"%s\"!"
		   ,
		   status,
		   IO_info.h_directory);			/*NOTREACHED*/
}
	  }

//******************************************************************************

//
// This function outputs an OpenDX control file to allow Thomas Radke's
// OpenDX macros to read in a "ASCII (gnuplot)" data file produced by
//  output_gridfn() .
//
namespace {
void output_OpenDX_control_file(const patch_system& ps,
				const struct IO_info& IO_info,
				int hn)
{
if (! IO_info.output_OpenDX_control_files) return;
static char file_name_buffer[IO_info::file_name_buffer_size];
snprintf(file_name_buffer, IO_info::file_name_buffer_size,
	 "%s/%s.ah%d.%s",
	 IO_info.h_directory, IO_info.h_base_file_name,
         hn, IO_info.OpenDX_control_file_name_extension);

FILE *fileptr = fopen(file_name_buffer, "w");
if (fileptr == NULL)
   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
		"output_OpenDX_control_file(): can't open output file \"%s\"!",
		   file_name_buffer);				/*NOTREACHED*/

fprintf(fileptr, "# list the size of each patch (N_rho x N_sigma)\n");
	for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
	{
	const patch& p = ps.ith_patch(pn);
	fprintf(fileptr, "object \"%s patch\" class array ", p.name());
	fprintf(fileptr, "type int rank 1 shape 2 items 1 data follows %d %d\n",
		p.effective_N_irho(IO_info.output_ghost_zones_for_h),
		p.effective_N_isigma(IO_info.output_ghost_zones_for_h));
	}
fprintf(fileptr, "\n");

fprintf(fileptr, "# collect all patch sizes into a single OpenDX group\n");
fprintf(fileptr, "# for the ImportAHFinderDirectGnuplot macro to read\n");
fprintf(fileptr, "object \"patchsizes\" class group\n");
	for (int pn = 0 ; pn < ps.N_patches() ; ++pn)
	{
	const patch& p = ps.ith_patch(pn);
	fprintf(fileptr, "member %d value \"%s patch\"\n", pn, p.name());
	}

fclose(fileptr);
}
	  }

//******************************************************************************

//
// This function outputs a gridfn from a data file.
//
// If the gridfn is h, then we also write out the xyz positions of the
// horizon surface points.
//
// FIXME: if the gridfn is not h, we assume that it's nominal-grid.
//
void output_gridfn(patch_system& ps, int unknown_gfn,
                   const char gfn_name[], const cGH *cctkGH,
		   const struct IO_info& IO_info, const char base_file_name[],
                   int min_digits,
		   int hn, bool print_msg_flag, int AHF_iteration /* = 0 */)
{

if (IO_info.output_ASCII_files)
   then	{
	const char* file_name
		= io_ASCII_file_name(IO_info, base_file_name, min_digits,
				     hn, AHF_iteration);
	// create the output directory (if it doesn't already exist)
	create_h_directory(IO_info);

	switch	(unknown_gfn)
		{
	case gfns::gfn__h:
		if (print_msg_flag)
		   then CCTK_VInfo(CCTK_THORNSTRING,
				   "writing h to \"%s\"", file_name);
		ps.print_ghosted_gridfn_with_xyz
		      (unknown_gfn,
		       true, gfns::gfn__h,
		       file_name,
		       IO_info.output_ghost_zones_for_h); // should we include
							  // ghost zones?
		break;
	case gfns::gfn__Theta:
		if (print_msg_flag)
		   then CCTK_VInfo(CCTK_THORNSTRING,
				   "writing Theta to \"%s\"", file_name);
		ps.print_gridfn(unknown_gfn, file_name);
		break;
	case gfns::gfn__mean_curvature:
		if (print_msg_flag)
		   then CCTK_VInfo(CCTK_THORNSTRING,
				   "writing mean curvature to \"%s\"", file_name);
		ps.print_gridfn(unknown_gfn, file_name);
		break;
	case gfns::gfn__Delta_h:
		if (print_msg_flag)
		   then CCTK_VInfo(CCTK_THORNSTRING,
				   "writing Delta_h to \"%s\"", file_name);
		ps.print_gridfn(unknown_gfn, file_name);
		break;
	default:
		if (print_msg_flag)
		   then CCTK_VInfo(CCTK_THORNSTRING,
				   "writing gfn=%d to \"%s\"",
				   unknown_gfn, file_name);
		ps.print_gridfn(unknown_gfn, file_name);
		break;
		}
	}

if (IO_info.output_HDF5_files)
   then	{
	const char* file_name
		= io_HDF5_file_name(IO_info, base_file_name, min_digits,
				     hn, AHF_iteration);
	// create the output directory (if it doesn't already exist)
	create_h_directory(IO_info);

	switch	(unknown_gfn)
		{
	case gfns::gfn__h:
		if (print_msg_flag)
		   then CCTK_VInfo(CCTK_THORNSTRING,
				   "writing h to \"%s\"", file_name);
		ps.output_ghosted_gridfn_with_xyz
		      (unknown_gfn, gfn_name, cctkGH,
		       true, gfns::gfn__h,
		       file_name,
		       IO_info.output_ghost_zones_for_h); // should we include
							  // ghost zones?
		break;
	case gfns::gfn__Theta:
		if (print_msg_flag)
		   then CCTK_VInfo(CCTK_THORNSTRING,
				   "writing Theta to \"%s\"", file_name);
		ps.output_gridfn(unknown_gfn, gfn_name, cctkGH,
                                 file_name);
		break;
	case gfns::gfn__Delta_h:
		if (print_msg_flag)
		   then CCTK_VInfo(CCTK_THORNSTRING,
				   "writing Delta_h to \"%s\"", file_name);
		ps.output_gridfn(unknown_gfn, gfn_name, cctkGH,
                                 file_name);
		break;
	default:
		if (print_msg_flag)
		   then CCTK_VInfo(CCTK_THORNSTRING,
				   "writing gfn=%d to \"%s\"",
				   unknown_gfn, file_name);
		ps.output_gridfn(unknown_gfn, gfn_name, cctkGH,
                                 file_name);
		break;
		}
	}
}

//******************************************************************************

//
// This function prints one or two Jacobian matrices (and their difference
// in the latter case) to a named output file.
//
// Bugs:
// - The file format is hardwired to ASCII.
//
// Arguments:
// A null Jacobian pointer means to skip that Jacobian.
//
void output_Jacobians(const patch_system& ps,
                      const Jacobian* Jac_NP_ptr,
                      const Jacobian* Jac_SD_FDdr_ptr,
		      const struct IO_info& IO_info, const char base_file_name[],
                      int min_digits,
		      int hn, bool print_msg_flag, int AHF_iteration /* = 0 */)
{
if (Jac_NP_ptr == NULL)
   then {
	CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
		   "output_Jacobians(): Jac_NP_ptr == NULL is not (yet) supported!");
	return;						// *** ERROR RETURN ***
	}

const char* file_name = io_ASCII_file_name(IO_info, base_file_name, min_digits,
                                           hn, AHF_iteration);
// create the output directory (if it doesn't already exist)
create_h_directory(IO_info);

if (print_msg_flag)
   then CCTK_VInfo(CCTK_THORNSTRING,
		   "   writing %s to \"%s\"",
		   ((Jac_SD_FDdr_ptr == NULL) ? "NP Jacobian"
					      : "NP and SD_FDdr Jacobians"),
		   file_name);

FILE *fileptr = fopen(file_name, "w");
if (fileptr == NULL)
   then CCTK_VWarn(FATAL_ERROR, __LINE__, __FILE__, CCTK_THORNSTRING,
		   "output_Jacobians(): can't open output file \"%s\"!",
		   file_name);					/*NOTREACHED*/

fprintf(fileptr, "# column 1 = x patch name\n");
fprintf(fileptr, "# column 2 = x patch number\n");
fprintf(fileptr, "# column 3 = x irho\n");
fprintf(fileptr, "# column 4 = x isigma\n");
fprintf(fileptr, "# column 5 = x II\n");
fprintf(fileptr, "# column 6 = y patch name\n");
fprintf(fileptr, "# column 7 = y patch number\n");
fprintf(fileptr, "# column 8 = y irho\n");
fprintf(fileptr, "# column 9 = y isigma\n");
fprintf(fileptr, "# column 10 = y JJ\n");
fprintf(fileptr, "# column 11 = Jac_NP.element(II,JJ)\n");
if (Jac_SD_FDdr_ptr != NULL)
   then {
	fprintf(fileptr, "# column 12 = Jac_SD_FDdr.element(II,JJ)\n");
	fprintf(fileptr, "# column 13 = abs error in Jac_SD_FDdr\n");
	fprintf(fileptr, "# column 14 = rel error in Jac_SD_FDdr\n");
	}

    for (int xpn = 0 ; xpn < ps.N_patches() ; ++xpn)
    {
    const patch& xp = ps.ith_patch(xpn);

    for (int x_irho = xp.min_irho() ; x_irho <= xp.max_irho() ; ++x_irho)
    {
    for (int x_isigma = xp.min_isigma() ;
	 x_isigma <= xp.max_isigma() ;
	 ++x_isigma)
    {
    const int II = ps.gpn_of_patch_irho_isigma(xp, x_irho,x_isigma);

	for (int ypn = 0 ; ypn < ps.N_patches() ; ++ypn)
	{
	const patch& yp = ps.ith_patch(ypn);

	for (int y_irho = yp.min_irho() ;
	     y_irho <= yp.max_irho() ;
	     ++y_irho)
	{
	for (int y_isigma = yp.min_isigma() ;
	     y_isigma <= yp.max_isigma() ;
	     ++y_isigma)
	{
	const int JJ = ps.gpn_of_patch_irho_isigma(yp, y_irho,y_isigma);

	if (! Jac_NP_ptr->is_explicitly_stored(II,JJ))
	   then continue;			// skip sparse points

	const fp NP = Jac_NP_ptr->element(II,JJ);
	const fp SD_FDdr = (Jac_SD_FDdr_ptr == NULL)
			   ? 0.0
			   : Jac_SD_FDdr_ptr->element(II,JJ);

	if ((NP == 0.0) && (SD_FDdr == 0.0))
	   then continue;		// skip zero values (if == )

	fprintf(fileptr,
		"%s %d %d %d %d\t%s %d %d %d %d\t%#.10g",
		xp.name(), xpn, x_irho, x_isigma, II,
		yp.name(), ypn, y_irho, y_isigma, JJ,
		double(NP));

	if (Jac_SD_FDdr_ptr != NULL)
	   then {
		const fp abs_NP      = jtutil::abs(NP     );
		const fp abs_SD_FDdr = jtutil::abs(SD_FDdr);
		const fp max_abs = jtutil::max(abs_SD_FDdr, abs_NP);
		const fp SD_FDdr_abs_error = SD_FDdr - NP;
		const fp SD_FDdr_rel_error = SD_FDdr_abs_error / max_abs;

		fprintf(fileptr,
			"\t%#.10g\t%e\t%e",
			double(SD_FDdr),
			double(SD_FDdr_abs_error), double(SD_FDdr_rel_error));
		}

	fprintf(fileptr, "\n");
	}
	}
	}
    }
    }
    }

fclose(fileptr);
}

//******************************************************************************
//******************************************************************************
//******************************************************************************

//
// This function encapsulates our file-naming conventions for angular-gridfn
// output files (those used for h, H, and other angular grid functions).
//
// Arguments:
// base_file_name[] = from the parameter file
// hn = the horizon number
// AHF_iteration = the apparent horizon finder's internal iteration
//		   number (>= 1) if this is an intermediate iterate,
//		   or the default (0) if this is a final computed
//		   horizon position
//
// Results:
// This function returns (a pointer to) the file name.  The returned
// result points into a private static buffer; the usual caveats apply.
//
namespace {
const char* io_ASCII_file_name(const struct IO_info& IO_info,
                               const char base_file_name[], int min_digits,
                               int hn, int AHF_iteration /* = 0 */)
{
static char file_name_buffer[IO_info::file_name_buffer_size];

const char* file_name_extension
	= IO_info.ASCII_gnuplot_file_name_extension;

if (AHF_iteration == 0)
   then snprintf(file_name_buffer, IO_info::file_name_buffer_size,
		 "%s/%s.t%0*d.ah%d.%s",
		 IO_info.h_directory, base_file_name,
		 min_digits, IO_info.time_iteration, hn,
		 file_name_extension);
   else snprintf(file_name_buffer, IO_info::file_name_buffer_size,
		 "%s/%s.t%0*d.ah%d.it%d.%s",
		 IO_info.h_directory, base_file_name,
		 min_digits, IO_info.time_iteration, hn, AHF_iteration,
		 file_name_extension);

return file_name_buffer;
}
	  }

//******************************************************************************

//
// This function encapsulates our file-naming conventions for angular-gridfn
// output files (those used for h, H, and other angular grid functions).
//
// Arguments:
// base_file_name[] = from the parameter file
// hn = the horizon number
// AHF_iteration = the apparent horizon finder's internal iteration
//		   number (>= 1) if this is an intermediate iterate,
//		   or the default (0) if this is a final computed
//		   horizon position
//
// Results:
// This function returns (a pointer to) the file name.  The returned
// result points into a private static buffer; the usual caveats apply.
//
namespace {
const char* io_HDF5_file_name(const struct IO_info& IO_info,
                              const char base_file_name[], int min_digits,
                              int hn, int AHF_iteration /* = 0 */)
{
static char file_name_buffer[IO_info::file_name_buffer_size];

const char* file_name_extension
	= IO_info.HDF5_file_name_extension;

if (AHF_iteration == 0)
   then snprintf(file_name_buffer, IO_info::file_name_buffer_size,
		 "%s/%s.ah%d.%s",
		 IO_info.h_directory, base_file_name,
		 hn,
		 file_name_extension);
   else snprintf(file_name_buffer, IO_info::file_name_buffer_size,
		 "%s/%s.ah%d.it%d.%s",
		 IO_info.h_directory, base_file_name,
		 hn, AHF_iteration,
		 file_name_extension);

return file_name_buffer;
}
	  }

//******************************************************************************

	  }	// namespace AHFinderDirect
