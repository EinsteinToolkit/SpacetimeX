// driver.hh -- header file for driver code
// $Header$

//
// prerequisites:
//	<stdio.h>
//	"horizon_sequence.hh"
//	"BH_diagnostics.hh"
//

// everything in this file is inside this namespace
namespace AHFinderDirect
	  {

//******************************************************************************

// Cactus likes 3-element arrays for xyz axes
enum	{ X_AXIS=0, Y_AXIS=1, Z_AXIS=2 };

//
// this enum holds the decoded  method  parameter, i.e. it specifies
// our top-level method
//
enum	method
	{
	method__evaluate_expansions,
	method__test_expansion_Jacobians,
	method__find_horizons // no comma
	};

//
// this enum holds the decoded  verbose_method  parameter, i.e. it
// specifies which (how many) informational messages we should print
//
enum	verbose_level
	{
	verbose_level__physics_highlights,
	verbose_level__physics_details,
	verbose_level__algorithm_highlights,
	verbose_level__algorithm_details,
	verbose_level__algorithm_debug // no comma
	};

//
// this enum holds the (a) decoded  initial_guess_method  parameter,
// i.e.  it specifies how we should set up the initial guess for a
// single apparent horizon
//
enum	initial_guess_method
	{
	initial_guess__read_from_named_file,
	initial_guess__read_from_h_file,
	initial_guess__Kerr_Kerr,
	initial_guess__Kerr_KerrSchild,
	initial_guess__coord_sphere,
	initial_guess__coord_ellipsoid // no comma
	};

//******************************************************************************

//
// This struct holds parameters for setting up the initial guess
// for a single apparent horizon.
//
struct	initial_guess_info
	{
	enum initial_guess_method method;
        bool reset_horizon_after_not_finding;

	// parameters for method == initial_guess__read_from_named_file
	struct	{
		const char* file_name;
		} read_from_named_file_info;

	// parameters for method == initial_guess__Kerr_Kerr
	struct	{
		fp x_posn, y_posn, z_posn;
		fp mass, spin;
		} Kerr_Kerr_info;

	// parameters for method == initial_guess__Kerr_KerrSchild
	struct	{
		fp x_posn, y_posn, z_posn;
		fp mass, spin;
		} Kerr_KerrSchild_info;

	// parameters for method == initial_guess__coord_sphere
	struct	{
		fp x_center, y_center, z_center;
		fp radius;
		} coord_sphere_info;

	// parameters for method == initial_guess__coord_ellipsoid
	struct	{
		fp x_center, y_center, z_center;
		fp x_radius, y_radius, z_radius;
		} coord_ellipsoid_info;
	};

//
// This struct holds parameters for solving the Theta(h) = 0 equations.
//
struct	solver_info
	{
	bool debugging_output_at_each_Newton_iteration;
	struct linear_solver_pars linear_solver_pars;
	int max_Newton_iterations__initial,
	    max_Newton_iterations__subsequent;
	fp max_allowable_Delta_h_over_h;
	fp Theta_norm_for_convergence;
	fp max_allowable_Theta;
        fp max_allowable_Theta_growth_iterations;
        fp max_allowable_Theta_nonshrink_iterations;
	fp *max_allowable_horizon_radius;	// --> new[]-allocated array
						//     of size  N_horizons+1 ,
						//     subscripted by hn
	bool want_expansion_gradients;
	};

//
// This struct holds info for I/O
//
struct	IO_info
	{
	// buffer size for file names
	// ... file names longer than this will be truncated, but
	//     they will NOT overflow the buffer (we always use snprintf(3))
	enum { file_name_buffer_size = 512 };

	// default permissions for newly created directories
	// = full access to everyone; the user's umask
	//   should cut this down to something reasonable
	enum { default_directory_permission = 0777 };

        bool output_ASCII_files;
        bool output_HDF5_files;
	bool output_initial_guess;
	int output_h_every, output_Theta_every;
	int output_mean_curvature_every;
	// based on the above, do we want to output things now (this time step)?
	bool output_h, output_Theta, output_mean_curvature;

	bool output_BH_diagnostics;
	const char* BH_diagnostics_directory;
	const char* BH_diagnostics_base_file_name;
	const char* BH_diagnostics_file_name_extension;

	bool output_ghost_zones_for_h;
	const char* ASCII_gnuplot_file_name_extension;
	const char* HDF5_file_name_extension;
	const char* h_directory;
	const char* h_base_file_name;
	const char* Theta_base_file_name;
	const char* mean_curvature_base_file_name;
	const char* Delta_h_base_file_name;
        int h_min_digits;

	const char* Jacobian_base_file_name;

	bool output_OpenDX_control_files;
	const char* OpenDX_control_file_name_extension;

	// this is used to choose file names
	int time_iteration;	// the Cactus time interation number
				// (cctk_iteration)
	fp time;		// the Cactus time coordinate (cctk_time)
	};

//
// This struct holds info describing how verbose we should be
//
struct	verbose_info
	{
	// decoded from verbose_level parameter
	enum verbose_level verbose_level;

	// derived Boolean flags saying whether or not
	// verbose_level is >= the appropriate level
	bool print_physics_highlights;
	bool print_physics_details;
	bool print_algorithm_highlights;
	bool print_algorithm_details;
	bool print_algorithm_debug;
	};

//
// This struct holds info for setting a mask grid function
// based on each apparent horizon shape.
//
struct	mask_info
	{
	bool  set_mask_for_any_horizon;
	bool* set_mask_for_this_horizon;
	CCTK_REAL radius_multiplier, radius_offset;
	CCTK_REAL buffer_thickness;
	bool mask_is_noshrink;
	CCTK_REAL min_horizon_radius_points_for_mask;
	bool set_old_style_mask, set_new_style_mask;
	struct	old_style_mask_info
		{
		const char* gridfn_name;
		int         gridfn_varindex;
		CCTK_REAL*  gridfn_dataptr;
		CCTK_REAL   inside_value, buffer_value, outside_value;
		} old_style_mask_info;
	struct	new_style_mask_info
		{
		const char* gridfn_name;
		int         gridfn_varindex;
		CCTK_INT*   gridfn_dataptr;
		const char* bitfield_name;
		CCTK_INT    bitfield_bitmask;
		const char* inside_value;
		const char* buffer_value;
		const char* outside_value;
		CCTK_INT    inside_bitvalue, buffer_bitvalue, outside_bitvalue;
		} new_style_mask_info;
	};

//
// This struct holds interprocessor-communication buffers for broadcasting
// status information from each active processor to all processors.
//
struct	iteration_status_buffers
	{
	// --> high-level buffers for broadcast-level semantics in Newton()
	// ... broadcast_status() [in Newton.c] sets these all to point to
	//     new-allocated arrays of size N_active_procs on its first call
	int* hn_buffer;
	int* iteration_buffer;
	enum expansion_status* expansion_status_buffer;
	fp* mean_horizon_radius_buffer;
	fp* Theta_infinity_norm_buffer;
	bool* found_horizon_buffer;

	// --> low-level buffers for CCTK_Reduce()
	// ... broadcast_status() [in Newton.c] sets these to point to
	//     new-allocated objects on its first call
	jtutil::array2d<CCTK_REAL> *send_buffer_ptr;
	jtutil::array2d<CCTK_REAL> *receive_buffer_ptr;

	iteration_status_buffers()
		: hn_buffer(NULL), iteration_buffer(NULL),
		  expansion_status_buffer(NULL),
		  mean_horizon_radius_buffer(NULL),
		  Theta_infinity_norm_buffer(NULL),
		  found_horizon_buffer(NULL),
		  send_buffer_ptr(NULL), receive_buffer_ptr(NULL)
		{ }
	};

//
// This struct holds interprocessor-communication buffers for broadcasting
// the BH diagnostics and horizon shape from the processor which finds a
// given horizon, to all processors.
//
struct	horizon_buffers
	{
	int N_buffer;		// number of CCTK_REAL values in each buffer
	CCTK_REAL* send_buffer;
	CCTK_REAL* receive_buffer;

	horizon_buffers()
		: N_buffer(0),
		  send_buffer(NULL),
		  receive_buffer(NULL)
		{ }
	};

//******************************************************************************

//
// (A single copy of) this struct holds all of our information about
// a single apparent horizon.
//
struct	AH_data
	{
	//
	// Any given horizon is allocated to a single processor.
	// On that processor (where we actually find the horizon)
	// we keep a "full-fledged" patch system and a Jacobian.
	// On other processors (where we only use the horizon shape
	// to (optionally) set the BH excision mask), we keep only a
	// "skeletal" patch system, and no Jacobian.
	//
	patch_system* ps_ptr;
	Jacobian* Jac_ptr;
        what_to_compute compute_info;

        bool move_origins;

        bool use_pretracking;
        int pretracking_max_iterations;

        fp pretracking_value;
        fp pretracking_minimum_value;
        fp pretracking_maximum_value;
        fp pretracking_delta;
        fp pretracking_minimum_delta;
        fp pretracking_maximum_delta;

        int depends_on;
        fp desired_value_factor;
        fp desired_value_offset;

        fp shiftout_factor;
        fp smoothing_factor;

	// are we finding this horizon "from scratch"?
	// ... true if this is the first time we've tried to find it,
	//          or if we've tried before but failed to find it the
	//          last time
	//     false if we've tried before and succeeded the last time
	//           (so we have that position as a very good initial guess)
	// ... also false if we're not finding this horizon on this processor
	bool initial_find_flag;
        // is this really the first time in this simulation that we are
        // trying to find a horizon?
        // ... true initially if this is a genuine horizon,
        //     false after the first time that a horizon has been found
	bool really_initial_find_flag;

	// used only if we're finding this horizon on this processor
	struct initial_guess_info initial_guess_info;

	bool search_flag;	// did we search for this horizon
	bool found_flag;	// did we find this horizon (successfully)
	bool h_files_written;	// have we written horizon-shape or similar
				// files for this horizon yet?

	struct BH_diagnostics BH_diagnostics;
	FILE *BH_diagnostics_fileptr;

	// interprocessor-communication buffers
	// for this horizon's BH diagnostics and (optionally) horizon shape
	struct horizon_buffers horizon_buffers;
	};

//
// (A single copy of) this struct holds all of our state on this processor
// that's persistent across Cactus scheduler calls.  This copy lives in
// "state.cc".
//
struct	state
	{
	enum method method;
	struct error_info error_info;
	struct verbose_info verbose_info;
	int timer_handle;
	int N_procs;			// total number of processors
	int my_proc;			// processor number of this processor
					// (0 to N_procs-1)

	int N_horizons;			// total number of genuine horizons
					// being searched for
	int N_active_procs;		// total number of active processors
					// (the active processors are processor
					//  numbers 0 to N_active_procs-1)

	struct cactus_grid_info cgi;
	struct geometry_info gi;
	struct Jacobian_info Jac_info;
	struct solver_info solver_info;
	struct IO_info IO_info;

	struct BH_diagnostics_info BH_diagnostics_info;
	struct mask_info mask_info;

	// interprocessor-communication buffers for broadcasting
	// Newton-iteration status from active processors to all processors
	struct iteration_status_buffers isb;

	horizon_sequence *my_hs;	// --> new-allocated object describing
					//     the sequence of genuine horizons
					//     assigned to this processor

	// horizon numbers ("hn") run from 1 to N_horizons inclusive
	struct AH_data** AH_data_array;	// --> new[]-allocated array of size
					//     N_horizons+1, subscripted by hn,
					//     of --> info
	};

//******************************************************************************

//
// prototypes for functions visible outside their source files
//

// setup.cc
// ... called from Cactus Scheduler
extern "C"
  void AHFinderDirect_setup(CCTK_ARGUMENTS);

// find_horizons.cc
// ... called from Cactus Scheduler
extern "C"
  void AHFinderDirect_find_horizons(CCTK_ARGUMENTS);

// announce.cc
// ... called from Cactus Scheduler
extern "C"
  void AHFinderDirect_announce(CCTK_ARGUMENTS);

// mask.cc
// ... called from Cactus Scheduler
extern "C"
  void AHFinderDirect_maybe_do_masks(CCTK_ARGUMENTS);

// aliased_functions.cc
// ... called from other thorns via the Cactus flesh function-aliasing mechanism
extern "C"
  CCTK_INT AHFinderDirect_local_coordinate_origin
    (CCTK_INT horizon_number,
     CCTK_REAL* origin_x_ptr, CCTK_REAL* origin_y_ptr, CCTK_REAL* origin_z_ptr);
extern "C"
  CCTK_INT AHFinderDirect_horizon_was_found(CCTK_INT horizon_number);
extern "C"
  CCTK_INT AHFinderDirect_radius_in_direction
    (CCTK_INT horizon_number,
     CCTK_INT N_points,
     const CCTK_REAL* const x, const CCTK_REAL* const y, const CCTK_REAL* const 
z,
     CCTK_REAL* const radius);

// initial_guess.cc
void setup_initial_guess(patch_system& ps,
			 const struct initial_guess_info& igi,
			 const struct IO_info& IO_info,
			 int hn, int N_horizons,
			 const struct verbose_info& verbose_info);
enum initial_guess_method
  decode_initial_guess_method(const char initial_guess_method_string[]);


void set_initial_guess_parameters(struct AH_data& AH_data, const int hn,
                                  const fp origin_x, const fp origin_y, const fp origin_z);

// Newton.cc
// returns true for success, false for failure to converge
void Newton(const cGH* GH,
	    int N_procs, int N_active_procs, int my_proc,
	    horizon_sequence& hs, struct AH_data* const AH_data_array[],
	    const struct cactus_grid_info& cgi,
	    const struct geometry_info& gi,
	    const struct Jacobian_info& Jacobian_info,
	    const struct solver_info& solver_info,
	    const struct IO_info& IO_info,
	    const struct BH_diagnostics_info& BH_diagnostics_info,
	    bool broadcast_horizon_shape,
	    const struct error_info& error_info,
	    const struct verbose_info& verbose_info,
	    struct iteration_status_buffers& isb);

// Tracks coordinate origin
void track_origin(const cGH* const cctkGH, patch_system& ps, 
                  struct AH_data* const AH_data_ptr, 
                  const int hn, const bool print_algorithm_details);

// io.cc
void input_gridfn(patch_system& ps, int unknown_gfn,
		  const struct IO_info& IO_info, const char base_file_name[],
		  int min_digits,
		  int hn, bool print_msg_flag, int AHF_iteration = 0);
void input_gridfn__explicit_name(patch_system& ps, int unknown_gfn,
				 const struct IO_info& IO_info,
				 const char file_name[], bool print_msg_flag);
void setup_h_files(patch_system& ps, const struct IO_info& IO_info, int hn);
void output_gridfn(patch_system& ps, int unknown_gfn,
                   const char gfn_name[], const cGH *cctkGH,
		   const struct IO_info& IO_info, const char base_file_name[],
                   int min_digits,
		   int hn, bool print_msg_flag, int AHF_iteration = 0);
void output_Jacobians(const patch_system& ps,
		      const Jacobian* Jac_NP_ptr,
		      const Jacobian* Jac_SD_FDdr_ptr,
		      const struct IO_info& IO_info, const char base_file_name[],
                      int min_digits,
		      int hn, bool print_msg_flag, int AHF_iteration = 0);

// misc-driver.cc
int Cactus_gridfn_varindex(const char gridfn_name[]);

//******************************************************************************

	  }	// namespace AHFinderDirect
