#include <fstream>
#include <iomanip>
#include <iostream>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

using namespace std;

extern "C"
void TwoPuncturesX_Metadata (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_TwoPuncturesX_Metadata;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_MyProc(cctkGH) == 0)
  {
    ofstream o;
    o.open(string(string(out_dir) + "/TwoPuncturesX.bbh").c_str());

    o << setprecision(17);

    o << "\
# ==================================\n\
# Numerical Relativity Metadata file\n\
# ==================================\n\
#\n\
# This file contains information about the simulation provided by the\n\
# TwoPuncturesX thorn.  The format is described in the NR Data Format Document\n\
# http://arxiv.org/abs/0709.0093 [draft SVN r707].\n\
" << "\n";

    o << "[metadata]" << "\n";
    o << "initial-ADM-energy            = " << *E << "\n";
    o << "initial-ADM-angular-momentumx = " << *J1 << "\n";
    o << "initial-ADM-angular-momentumy = " << *J2 << "\n";
    o << "initial-ADM-angular-momentumz = " << *J3 << "\n";
    o << "initial-separation            = " << par_b * 2 << "\n";
    o << "initial-data-type             = Bowen-York" << "\n";
    o << "initial-data-bibtex-keys      = Bowen:1980yu Brandt:1997tf Ansorg:2004ds" << "\n";
    o << "initial-bh-position1x         = " << par_b + center_offset[0] << "\n";
    o << "initial-bh-position1y         = " << center_offset[1] << "\n";
    o << "initial-bh-position1z         = " << center_offset[2] << "\n";
    o << "initial-bh-position2x         = " << -par_b + center_offset[0] << "\n";
    o << "initial-bh-position2y         = " << center_offset[1] << "\n";
    o << "initial-bh-position2z         = " << center_offset[2] << "\n";
    o << "initial-bh-momentum1x         = " << par_P_plus[0] << "\n";
    o << "initial-bh-momentum1y         = " << par_P_plus[1] << "\n";
    o << "initial-bh-momentum1z         = " << par_P_plus[2] << "\n";
    o << "initial-bh-momentum2x         = " << par_P_minus[0] << "\n";
    o << "initial-bh-momentum2y         = " << par_P_minus[1] << "\n";
    o << "initial-bh-momentum2z         = " << par_P_minus[2] << "\n";
    o << "initial-bh-spin1x             = " << par_S_plus[0] << "\n";
    o << "initial-bh-spin1y             = " << par_S_plus[1] << "\n";
    o << "initial-bh-spin1z             = " << par_S_plus[2] << "\n";
    o << "initial-bh-spin2x             = " << par_S_minus[0] << "\n";
    o << "initial-bh-spin2y             = " << par_S_minus[1] << "\n";
    o << "initial-bh-spin2z             = " << par_S_minus[2] << "\n";
    o << "initial-bh-puncture-adm-mass1 = " << *mp_adm << "\n";
    o << "initial-bh-puncture-adm-mass2 = " << *mm_adm << "\n";
    o << "initial-bh-puncture-bare-mass1 = " << *mp << "\n";
    o << "initial-bh-puncture-bare-mass2 = " << *mm << "\n";
  }
}
