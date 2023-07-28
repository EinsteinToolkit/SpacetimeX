! $Header$

#include "cctk.h"

module cactus
  implicit none
  public
  
  external     cctk_pointerto
  CCTK_POINTER cctk_pointerto
  
  interface
     
     ! from Cactus:
     
     subroutine cctk_barrier (ierr, cctkgh)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
     end subroutine cctk_barrier
     
     subroutine cctk_coordrange (ierr, cctkgh, lower, upper, dir, name, &
          &                      systemname)
       implicit none
       integer ierr
       CCTK_POINTER cctkgh
       CCTK_REAL    lower
       CCTK_REAL    upper
       integer      dir
       character(*) name
       character(*) systemname
     end subroutine cctk_coordrange
     
     subroutine cctk_coordsystemhandle (handle, coordsystem)
       implicit none
       integer      handle
       character(*) coordsystem
     end subroutine cctk_coordsystemhandle
     
     function cctk_equals (arg1, arg2)
       implicit none
       integer      cctk_equals
       CCTK_POINTER arg1
       character(*) arg2
     end function cctk_equals
     
     subroutine cctk_fortranstring (clength, cstring, fortstring)
       implicit none
       CCTK_INT     clength     ! intent(out)
       CCTK_POINTER cstring
       character(*) fortstring
     end subroutine cctk_fortranstring
     
     subroutine cctk_groupbboxgn (ierr, cctkgh, nbbox, bbox, groupname)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       integer      nbbox
       integer      bbox(nbbox)
       character(*) groupname
     end subroutine cctk_groupbboxgn
     
     subroutine cctk_groupgshgn (ierr, cctkgh, ngsh, gsh, groupname)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       integer      ngsh
       integer      gsh(ngsh)
       character(*) groupname
     end subroutine cctk_groupgshgn
     
     subroutine cctk_grouplbndgn (ierr, cctkgh, nlbnd, lbnd, groupname)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       integer      nlbnd
       integer      lbnd(nlbnd)
       character(*) groupname
     end subroutine cctk_grouplbndgn
     
     subroutine cctk_grouplshgn (ierr, cctkgh, nlsh, lsh, groupname)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       integer      nlsh
       integer      lsh(nlsh)
       character(*) groupname
     end subroutine cctk_grouplshgn
     
     subroutine cctk_groupnghostzonesgn &
          (ierr, cctkgh, nnghostzones, nghostzones, groupname)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       integer      nnghostzones
       integer      nghostzones(nnghostzones)
       character(*) groupname
     end subroutine cctk_groupnghostzonesgn
     
     subroutine cctk_groupubndgn (ierr, cctkgh, nubnd, ubnd, groupname)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       integer      nubnd
       integer      ubnd(nubnd)
       character(*) groupname
     end subroutine cctk_groupubndgn
     
     subroutine cctk_info (thorn, message)
       implicit none
       character(*) thorn
       character(*) message
     end subroutine cctk_info
     
     subroutine cctk_interpgridarrays &
          (ierr, cctkGH, N_dims, &
          local_interp_handle, param_table_handle, coord_system_handle, &
          N_interp_points, interp_coords_type, interp_coords, &
          N_input_arrays, input_array_indices, &
          N_output_arrays, output_array_types, output_arrays)
       implicit none
       integer      ierr
       CCTK_POINTER cctkGH
       integer      N_dims
       integer      local_interp_handle
       integer      param_table_handle
       integer      coord_system_handle
       integer      N_interp_points
       integer      interp_coords_type
       CCTK_POINTER interp_coords(N_dims)
       integer      N_input_arrays
       CCTK_INT     input_array_indices(N_input_arrays)
       integer      N_output_arrays
       CCTK_INT     output_array_types(N_output_arrays)
       CCTK_POINTER output_arrays(N_output_arrays)
     end subroutine cctk_interpgridarrays
     
     subroutine cctk_interphandle (handle, interp)
       implicit none
       integer      handle
       character(*) interp
     end subroutine cctk_interphandle
     
     function cctk_isthornactive (name)
       implicit none
       integer      cctk_isthornactive
       character(*) name
     end function cctk_isthornactive
     
     function cctk_myproc (cctkgh)
       implicit none
       integer      cctk_myproc
       CCTK_POINTER cctkgh
     end function cctk_myproc
     
     function cctk_nprocs (cctkgh)
       implicit none
       integer      cctk_nprocs
       CCTK_POINTER cctkgh
     end function cctk_nprocs
     
     function cctk_nullpointer ()
       implicit none
       CCTK_POINTER cctk_nullpointer
     end function cctk_nullpointer
     
     subroutine cctk_paramwarn (thorn, message)
       implicit none
       character(*) thorn
       character(*) message
     end subroutine cctk_paramwarn
     
!!$     function cctk_pointerto (var)
!!$       implicit none
!!$       CCTK_POINTER cctk_pointerto
!!$       ???          var
!!$     end function cctk_pointerto
     
     subroutine cctk_reductionarrayhandle (handle, reduction)
       implicit none
       integer      handle
       character(*) reduction
     end subroutine cctk_reductionarrayhandle
     
     subroutine cctk_registerbanner (ierr, banner)
       implicit none
       integer      ierr
       character(*) banner
     end subroutine cctk_registerbanner
     
     subroutine cctk_syncgroup (ierr, cctkgh, groupname)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       character(*)  groupname
     end subroutine cctk_syncgroup
     
     subroutine cctk_syncgroupi (ierr, cctkgh, group)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       integer      group
     end subroutine cctk_syncgroupi
     
     subroutine cctk_syncgroupwithvar (ierr, cctkgh, varname)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       character(*) varname
     end subroutine cctk_syncgroupwithvar
     
     subroutine cctk_syncgroupwithvari (ierr, cctkgh, var)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       integer      var
     end subroutine cctk_syncgroupwithvari
     
     subroutine cctk_varindex (vi, vname)
       implicit none
       integer      vi
       character(*) vname
     end subroutine cctk_varindex
     
     subroutine cctk_vartypei (vartype, varindex)
       implicit none
       integer vartype
       integer varindex
     end subroutine cctk_vartypei
     
     subroutine cctk_warn (level, line, file, thorn, message)
       implicit none
       integer      level
       integer      line
       character(*) file
       character(*) thorn
       character(*) message
     end subroutine cctk_warn
     
     
     
     ! from thorn CactusBase/Boundary:
     subroutine bndcopygn (ierr, cctkgh, sw, togroupname, fromgroupname)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       integer      sw(3)
       character(*) togroupname
       character(*) fromgroupname
     end subroutine bndcopygn
     
     subroutine bndflatgn (ierr, cctkgh, sw, groupname)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       integer      sw(3)
       character(*) groupname
     end subroutine bndflatgn
     
     subroutine bndradiativegn (ierr, cctkgh, sw, var0, v0, &
          &                     togroupname, fromgroupname)
       implicit none
       integer       ierr
       CCTK_POINTER  cctkgh
       integer       sw(3)
       CCTK_REAL     var0
       CCTK_REAL     v0
       character(*)  togroupname
       character(*)  fromgroupname
     end subroutine bndradiativegn
     
     subroutine bndradiativevn (ierr, cctkgh, sw, var0, v0, &
          &                     tovarname, fromvarname)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       integer      sw(3)
       CCTK_REAL    var0
       CCTK_REAL    v0
       character(*) tovarname
       character(*) fromvarname
     end subroutine bndradiativevn
     
     subroutine bndrobinvn (ierr, cctkgh, sw, finf, npow, varname)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       integer      sw(3)
       CCTK_REAL    finf
       integer      npow
       character(*) varname
     end subroutine bndrobinvn
     
     subroutine bndscalargn (ierr, cctkgh, sw, var0, groupname)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       integer      sw(3)
       CCTK_REAL    var0
       character(*) groupname
     end subroutine bndscalargn
     
     subroutine bndscalarvn (ierr, cctkgh, sw, var0, varname)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       integer      sw(3)
       CCTK_REAL    var0
       character(*) varname
     end subroutine bndscalarvn
     
     
     
     ! from thorn CactusBase/CartGrid3D:
     
     subroutine cartsymgn (ierr, cctkgh, groupname)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       character(*) groupname
     end subroutine cartsymgn
     
     subroutine setcartsymvn (ierr, cctkgh, sw, varname)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       integer      sw(3)
       character(*) varname
     end subroutine setcartsymvn
     
     
     
     ! from thorn AlphaThorns/Cart3D:
     
     subroutine cart3dsymvi (ierr, cctkgh, vi)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       integer      vi
     end subroutine cart3dsymvi
     
     subroutine cart3dsymvn (ierr, cctkgh, vn)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       character(*) vn
     end subroutine cart3dsymvn
     
     subroutine cart3dsymgi (ierr, cctkgh, gi)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       integer      gi
     end subroutine cart3dsymgi
     
     subroutine cart3dsymgn (ierr, cctkgh, gn)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       character(*) gn
     end subroutine cart3dsymgn
     
     subroutine cart3dsettensortypevi (ierr, cctkgh, nvars, vi, typename)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       integer      nvars
       integer      vi(nvars)
       character(*) typename
     end subroutine cart3dsettensortypevi
     
     subroutine cart3dsettensortypevn (ierr, cctkgh, vn, typename)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       character(*) vn
       character(*) typename
     end subroutine cart3dsettensortypevn
     
     subroutine cart3dgetsymmetriesvi (ierr, cctkgh, sym, vi)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       integer      sym(6)
       integer      vi
     end subroutine cart3dgetsymmetriesvi
     
     subroutine cart3dgetsymmetriesvn (ierr, cctkgh, sym, vn)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       integer      sym(6)
       character(*) vn
     end subroutine cart3dgetsymmetriesvn
     
     subroutine cart3dgetsymmetryboundaries (ierr, cctkgh, symbnd)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       integer      symbnd(6)
     end subroutine cart3dgetsymmetryboundaries
     
     subroutine cart3dgetsymmetric (ierr, cctkgh, symmetric)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       integer      symmetric(6)
     end subroutine cart3dgetsymmetric
     
     subroutine cart3dgetstaggered (ierr, cctkgh, staggered)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       integer      staggered(6)
     end subroutine cart3dgetstaggered
     
     subroutine cart3dgetperiodic (ierr, cctkgh, periodic)
       implicit none
       integer      ierr
       CCTK_POINTER cctkgh
       integer      periodic(6)
     end subroutine cart3dgetperiodic
     
  end interface
  
end module cactus
