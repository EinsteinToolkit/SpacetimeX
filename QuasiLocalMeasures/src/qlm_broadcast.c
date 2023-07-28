#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#ifdef CCTK_MPI
#  include <mpi.h>
#endif



/* Broadcast a vector element of a grid group */
static void
bcast (cGH const * restrict const cctkGH,
       char const * restrict const group,
       int const vecind,
       int const root)
{
#ifdef CCTK_MPI
  
  int ierr;
  
  assert (group);
  
  MPI_Comm comm;
  if (CCTK_IsFunctionAliased ("GetMPICommWorld"))
  {
    comm = * (MPI_Comm const *) GetMPICommWorld (cctkGH);
  }
  else
  {
    comm = MPI_COMM_WORLD;
  }
  
  int num_procs;
  MPI_Comm_size (comm, & num_procs);
  
  assert (root >= 0 && root < num_procs);
  
  int const gi = CCTK_GroupIndex (group);
  assert (gi >= 0);
  
  cGroup data;
  ierr = CCTK_GroupData (gi, & data);
  assert (! ierr);
  
  assert (vecind >= 0 && vecind < data.vectorlength);
  
  if (data.numvars == 0) return;
  
  MPI_Datatype mpitype;
  int items;
  switch (data.vartype)
  {
  case CCTK_VARIABLE_INT:
    if (sizeof (CCTK_INT) == sizeof (int)) {
      mpitype = MPI_INT;
    } else if (sizeof (CCTK_INT) == sizeof (long)) {
      mpitype = MPI_LONG;
    } else if (sizeof (CCTK_INT) == sizeof (long long)) {
      mpitype = MPI_LONG_LONG;
    } else {
      CCTK_ERROR("Unsupported CCTK_INT type");
    }
    items = 1;
    break;
  case CCTK_VARIABLE_REAL:
    assert (sizeof (CCTK_REAL) == sizeof (double));
    mpitype = MPI_DOUBLE;
    items = 1;
    break;
  case CCTK_VARIABLE_COMPLEX:
    assert (sizeof (CCTK_COMPLEX) == 2 * sizeof (double));
    mpitype = MPI_DOUBLE;
    items = 2;
    break;
  default:
    assert (0);
  }
  assert (items > 0);
  
  cGroupDynamicData dyndata;
  ierr = CCTK_GroupDynamicData (cctkGH, gi, & dyndata);
  assert (! ierr);
  
  int elems = 1;
  {
    int d;
    for (d = 0; d < dyndata.dim; ++ d)
    {
      elems *= dyndata.lsh[d];
    }
  }
  assert (elems >= 0);
  
  int const v0 = CCTK_FirstVarIndexI (gi);
  assert (v0 >= 0);
  
  int size;
  MPI_Type_size (mpitype, & size);
  assert (items * size == CCTK_VarTypeSize (data.vartype));
  
  int const numvectors = data.numvars / data.vectorlength;
  
  {
    int vec;
    for (vec = 0; vec < numvectors; ++ vec)
    {
      int const vi = v0 + vecind + data.vectorlength * vec;
      void * restrict const ptr = CCTK_VarDataPtrI (cctkGH, 0, vi);
      assert (ptr);
      
      MPI_Bcast (ptr, elems * items, mpitype, root, comm);
    }
  }

#endif  /* ifdef CCTK_MPI */
}


void CCTK_FCALL
CCTK_FNAME(qlm_broadcast) (CCTK_POINTER_TO_CONST * restrict cctkGH_);

void CCTK_FCALL
CCTK_FNAME(qlm_broadcast) (CCTK_POINTER_TO_CONST * restrict const cctkGH_)
{
  cGH const * restrict const cctkGH = * (cGH const * const *) cctkGH_;
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (veryverbose)
  {
    CCTK_INFO ("Broadcasting quasi-local measures");
  }
  
  int const num_procs = CCTK_nProcs (cctkGH);
  int hn; 
  for (hn = 0; hn < num_surfaces; ++ hn)
  {
    int const root = hn % num_procs;
    
    bcast (cctkGH, "QuasiLocalMeasures::qlm_state"      , hn, root);
    bcast (cctkGH, "QuasiLocalMeasures::qlm_state_p"    , hn, root);
    bcast (cctkGH, "QuasiLocalMeasures::qlm_grid_int"   , hn, root);
    bcast (cctkGH, "QuasiLocalMeasures::qlm_grid_real"  , hn, root);
    bcast (cctkGH, "QuasiLocalMeasures::qlm_grid_real_p", hn, root);
    
    bcast (cctkGH, "QuasiLocalMeasures::qlm_shapes"               , hn, root);
    bcast (cctkGH, "QuasiLocalMeasures::qlm_shapes_p"             , hn, root);
    bcast (cctkGH, "QuasiLocalMeasures::qlm_coordinates"          , hn, root);
    bcast (cctkGH, "QuasiLocalMeasures::qlm_tetrad_l"             , hn, root);
    bcast (cctkGH, "QuasiLocalMeasures::qlm_tetrad_n"             , hn, root);
    bcast (cctkGH, "QuasiLocalMeasures::qlm_tetrad_m"             , hn, root);
    bcast (cctkGH, "QuasiLocalMeasures::qlm_newman_penrose"       , hn, root);
    bcast (cctkGH, "QuasiLocalMeasures::qlm_weyl_scalars"         , hn, root);
    bcast (cctkGH, "QuasiLocalMeasures::qlm_ricci_scalars"        , hn, root);
    bcast (cctkGH, "QuasiLocalMeasures::qlm_twometric"            , hn, root);
    bcast (cctkGH, "QuasiLocalMeasures::qlm_killing_vector"       , hn, root);
    bcast (cctkGH, "QuasiLocalMeasures::qlm_killed_twometric"     , hn, root);
    bcast (cctkGH, "QuasiLocalMeasures::qlm_invariant_coordinates", hn, root);
    
    bcast (cctkGH, "QuasiLocalMeasures::qlm_multipole_moments", hn, root);
    
    bcast (cctkGH, "QuasiLocalMeasures::qlm_3determinant", hn, root);
    
    bcast (cctkGH, "QuasiLocalMeasures::qlm_scalars"  , hn, root);
    bcast (cctkGH, "QuasiLocalMeasures::qlm_scalars_p", hn, root);
    
  } /* for hn */
}
