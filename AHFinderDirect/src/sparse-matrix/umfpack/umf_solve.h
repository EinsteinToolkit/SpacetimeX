/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

GLOBAL Int UMF_solve
(
    Int sys,
    const Int Ap [ ],
    const Int Ai [ ],
    const double Ax [ ],
    double Xx [ ],
    const double Bx [ ],
#ifdef COMPLEX
    const double Az [ ],
    double Xz [ ],
    const double Bz [ ],
#endif
    NumericType *Numeric,
    Int irstep,
    double Info [UMFPACK_INFO],
    Int Pattern [ ],
    double SolveWork [ ]
) ;

