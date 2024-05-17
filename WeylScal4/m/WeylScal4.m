SetEnhancedTimes[False];
SetSourceLanguage["C"];

$RecursionLimit=1000;

(****************************************************************************
  Derivatives
 ****************************************************************************)

derivatives =
{
  PDstandard2nd[i_] -> StandardCenteredDifferenceOperator[1,1,i],
  PDstandard2nd[i_, i_] -> StandardCenteredDifferenceOperator[2,1,i],
  PDstandard2nd[i_, j_] -> StandardCenteredDifferenceOperator[1,1,i] *
    StandardCenteredDifferenceOperator[1,1,j],

  PDstandard4th[i_] -> StandardCenteredDifferenceOperator[1,2,i],
  PDstandard4th[i_, i_] -> StandardCenteredDifferenceOperator[2,2,i],
  PDstandard4th[i_, j_] -> StandardCenteredDifferenceOperator[1,2,i] *
    StandardCenteredDifferenceOperator[1,2,j],

  PDstandard[i_] -> StandardCenteredDifferenceOperator[1,fdOrder/2,i],
  PDstandard[i_, i_] -> StandardCenteredDifferenceOperator[2,fdOrder/2,i],
  PDstandard[i_, j_] -> StandardCenteredDifferenceOperator[1,fdOrder/2,i] *
    StandardCenteredDifferenceOperator[1,fdOrder/2,j]
};

(****************************************************************************
  Tensors 
 ****************************************************************************)

(* Register all the tensors that will be used with TensorTools *)
Map[DefineTensor, 
{
  R, gamma,g, gInv, k, ltet, n, rm, im, rmbar, imbar, tn, va, vb, vc,
  wa, wb, wc, ea, eb, ec,
  Ro, Rojo, R4p, Psi0r, Psi0i, Psi1r, Psi1i, Psi2r, Psi2i, Psi3r,
  Psi3i, Psi4r,Psi4i, curvIr, curvIi, curvJr, curvJi, curvJ1, curvJ2, curvJ3, curvJ4
}];

(* Psi0,2,4 behave as (pseudo)scalars *)
SetTensorAttribute[Psi0r, TensorParity, +1];
SetTensorAttribute[Psi2r, TensorParity, +1];
SetTensorAttribute[Psi4r, TensorParity, +1];
SetTensorAttribute[Psi0i, TensorParity, -1];
SetTensorAttribute[Psi2i, TensorParity, -1];
SetTensorAttribute[Psi4i, TensorParity, -1];
(* these are 'special' and require a patch to the symmetry thorns *)
SetTensorAttribute[Psi1r, TensorManualCartesianParities, {1,1,-1}];
SetTensorAttribute[Psi3r, TensorManualCartesianParities, {1,1,-1}];
SetTensorAttribute[Psi1i, TensorManualCartesianParities, {-1,-1,1}];
SetTensorAttribute[Psi3i, TensorManualCartesianParities, {-1,-1,1}];

(* The I and J invariants have parities which can be derived from the
   parities of the Psis *)
SetTensorAttribute[curvIr, TensorParity, +1];
SetTensorAttribute[curvIi, TensorParity, -1];
SetTensorAttribute[curvJr, TensorParity, +1];
SetTensorAttribute[curvJi, TensorParity, -1];

(* TODO: set the parities of the invariants J1, J2, J3 and J4 *)

(* Register the TensorTools symmetries (this is very simplistic) *)
Map[AssertSymmetricDecreasing, 
{
  k[la,lb], g[la,lb]
}];

AssertSymmetricDecreasing[gamma[ua,lb,lc], lb, lc];

Tensor[R, i_, j_, k_, l_] /; i > j := -R[j, i, k, l];
Tensor[R, i_, j_, k_, l_] /; i == j := 0;
Tensor[R, i_, j_, k_, l_] /; k > l := -R[i, j, l, k];
Tensor[R, i_, j_, k_, l_] /; k == l := 0;
Tensor[R, i_, j_, k_, l_] /; i > k := R[k, l, i, j];
Tensor[R, i_, j_, k_, l_] /; i == k && j > l := R[k, l, i, j];

Tensor[R4p, i_, j_, k_, l_] /; i > j := -R4p[j, i, k, l];
Tensor[R4p, i_, j_, k_, l_] /; i == j := 0;
Tensor[R4p, i_, j_, k_, l_] /; k > l := -R4p[i, j, l, k];
Tensor[R4p, i_, j_, k_, l_] /; k == l := 0;
Tensor[R4p, i_, j_, k_, l_] /; i > k := R4p[k, l, i, j];
Tensor[R4p, i_, j_, k_, l_] /; i == k && j > l := R4p[k, l, i, j];

(* Determinants of the metrics in terms of their components
   (Mathematica symbolic expressions) *)
gDet = Det[MatrixOfComponents[g[la,lb]]];

(****************************************************************************
  Groups
 ****************************************************************************)

(* Cactus group definitions *)

scalars = {Psi0r, Psi0i, Psi1r, Psi1i, Psi2r, Psi2i, Psi3r, Psi3i, Psi4r,Psi4i, curvIr, curvIi, curvJr, curvJi, curvJ1, curvJ2, curvJ3, curvJ4};

scalarGroups = Map[CreateGroupFromTensor, scalars];

(* We need extra timelevels so that interpolation onto the extraction sphere
   works properly with mesh refinement *)
scalarGroups = Map[AddGroupExtra[#, Timelevels -> 3] &, scalarGroups];

admGroups = 
  {{"admbase::metric", {gxx,gxy,gxz,gyy,gyz,gzz}},
   {"admbase::curv", {kxx,kxy,kxz,kyy,kyz,kzz}},
   {"admbase::lapse", {alp}},
   {"admbase::shift", {betax,betay,betaz}}};

declaredGroups = scalarGroups;
declaredGroupNames = Map[First, declaredGroups];

groups = Join[declaredGroups, admGroups];

(****************************************************************************
  Shorthands
 ****************************************************************************)

shorthands = 
{
  gamma[ua,lb,lc], R[la,lb,lc,ld], invdetg, detg, third, detgmthirdroot,
  gInv[ua,ub], Ro[la,lb,lc], Rojo[la,lb], R4p[li,lj,lk,ll],
  omega11, omega22, omega33, omega12, omega13, omega23, va[ua], vb[ua], vc[ua],
  wa[ua], wb[ua], wc[ua], ea[ua], eb[ua], ec[ua],
  tn[ua], nn, ltet[ua],n[ua],rm[ua],im[ua],rmbar[ua],imbar[ua], isqrt2, xmoved,
  ymoved, zmoved
};

k11=kxx; k21=kxy; k22=kyy; k31=kxz; k32=kyz; k33=kzz;
g11=gxx; g21=gxy; g22=gyy; g31=gxz; g32=gyz; g33=gzz;

realParameters = {{Name -> offset, Default -> 10^(-15)},xorig,yorig,zorig};

psi4Eqs[PD_] := Flatten@{
  detg -> gDet,
  invdetg -> 1 / detg,
  gInv[ua,ub] -> invdetg gDet MatrixInverse[g[ua,ub]],
  gamma[ua, lb, lc] -> 1/2 gInv[ua,ud] (PD[g[lb,ld], lc] + PD[g[lc,ld], lb] - PD[g[lb,lc],ld]),

  (****************************************************************************
    Offset the origin
   ****************************************************************************)

  xmoved -> x - xorig,
  ymoved -> y - yorig,
  zmoved -> z - zorig,

  (****************************************************************************
    Compute the local tetrad
   ****************************************************************************)
  (* azmuthal *)
  va1 -> -ymoved, va2 -> xmoved+offset, va3 -> 0,

  (* radial *)
  vb1 -> xmoved+offset, vb2 -> ymoved, vb3 -> zmoved,

  (* polar *)
  vc[ua] -> Sqrt[detg] gInv[ua,ud] Eps[ld,lb,lc] va[ub] vb[uc],

  (* Orthonormalize using Gram-Schmidt*)
  (* Orthonormalize in the order phi, r, theta *)
  wa[ua] -> va[ua],
  omega11 -> wa[ua] wa[ub] g[la,lb],
  ea[ua] -> wa[ua] / Sqrt[omega11],

  omega12 -> ea[ua] vb[ub] g[la,lb],
  wb[ua] -> vb[ua] - omega12 ea[ua],
  omega22 -> wb[ua] wb[ub] g[la,lb],
  eb[ua] -> wb[ua]/Sqrt[omega22],

  omega13 -> ea[ua] vc[ub] g[la,lb],
  omega23 -> eb[ua] vc[ub] g[la,lb],
  wc[ua] -> vc[ua] - omega13 ea[ua] - omega23 eb[ua],
  omega33 -> wc[ua] wc[ub] g[la,lb],
  ec[ua] -> wc[ua]/Sqrt[omega33],

  (* Create Spatial Portion of Null Tetrad *)
  isqrt2  ->  0.7071067811865475244,
  ltet[ua] -> isqrt2 eb[ua],
  n[ua] -> - isqrt2 eb[ua],
  rm[ua] -> isqrt2 ec[ua],
  im[ua] -> isqrt2 ea[ua],
  rmbar[ua] -> isqrt2 ec[ua],
  imbar[ua] -> -isqrt2 ea[ua],

  (* nn here is the projection of both l^a and n^a with u^a (the time-like unit
     vector normal to the hypersurface).  We do NOT save the t component of the
     tetrads in this code to avoid unnecessary factors of lapse and shift. *)
  nn -> isqrt2,


  (****************************************************************************
    Compute the NP pseudoscalars
   ****************************************************************************)

  (* Calculate the relevant Riemann Quantities *)
		
  (* The 3-Riemann *)
  R[la,lb,lc,ld] -> 1/2 ( PD[g[la,ld],lc,lb] + PD[g[lb,lc],ld,la] )
			      - 1/2 ( PD[g[la,lc],lb,ld] + PD[g[lb,ld],la,lc] )
			      + g[lj,le] gamma[uj,lb,lc] gamma[ue,la,ld]
			      - g[lj,le] gamma[uj,lb,ld] gamma[ue,la,lc],

  (* The 4-Riemann projected into the slice on all its indices. The Gauss equation. *)
  R4p[li,lj,lk,ll] -> R[li,lj,lk,ll] + 2 AntiSymmetrize[k[li,lk] k[ll,lj], lk, ll],

  (* The 4-Riemann projected in the unit normal direction on one
     index, then into the slice on the remaining indices. The Codazzi equation. *)
  Ro[lj,lk,ll] -> - 2 AntiSymmetrize[ PD[k[lj,lk],ll], lk,ll]
			      - 2 AntiSymmetrize[ gamma[up,lj,lk] k[ll,lp], lk,ll],

  (* The 4-Riemann projected in the unit normal direction on two
     indices, and into the slice on the remaining two. *)
  Rojo[lj,ll] ->  gInv[uc,ud] (R[lj,lc,ll,ld] ) - k[lj,lp] gInv[up,ud] k[ld,ll]
			      + k[lc,ld] gInv[uc,ud] k[lj,ll],

  (* Calculate End Quantities
     NOTE: In writing this, I assume m[0]=0!! to save lots of work *)

  Psi4r -> R4p[li,lj,lk,ll] n[ui] n[uk] ( rmbar[uj] rmbar[ul] - imbar[uj] imbar[ul] )
		 + 2 Ro[lj,lk,ll] n[uk] nn ( rmbar[uj] rmbar[ul] - imbar[uj] imbar[ul] )
		 + Rojo[lj,ll] nn nn ( rmbar[uj] rmbar[ul] - imbar[uj] imbar[ul] 
         (* + terms in mbar^0 == 0 *) ),

  Psi4i -> R4p[la,lb,lc,ld] n[ua] n[uc] ( - rm[ub] im[ud] - im[ub] rm[ud] )
		 + 2 Ro[la,lb,lc] n[ub] nn ( - rm[ua] im[uc] - im[ua] rm[uc] )
		 + Rojo[la,lb] nn nn ( - rm[ua] im[ub] - im[ua] rm[ub] )
};

otherPsiEqs = {
  Psi3r -> R4p[la,lb,lc,ld] ltet[ua] n[ub] rm[uc] n[ud]
		 + Ro[la,lb,lc] ( nn (n[ua]-ltet[ua]) rm[ub] n[uc] - nn rm[ua] ltet[ub] n[uc] )
		 - Rojo[la,lb] nn (n[ua]-ltet[ua]) nn rm[ub],

  Psi3i -> - R4p[la,lb,lc,ld] ltet[ua] n[ub] im[uc] n[ud]
		 - Ro[la,lb,lc] ( nn (n[ua]-ltet[ua]) im[ub] n[uc] - nn im[ua] ltet[ub] n[uc] )
		 + Rojo[la,lb] nn (n[ua]-ltet[ua]) nn im[ub],

  Psi2r -> R4p[la,lb,lc,ld] ltet[ua] n[ud] (rm[ub] rm[uc] + im[ub] im[uc])
		 + Ro[la,lb,lc] nn ( n[uc] (rm[ua] rm[ub] + im[ua] im[ub]) - ltet[ub] (rm[ua] rm[uc] + im[ua] im[uc]) )
		 - Rojo[la,lb] nn nn (rm[ua] rm[ub] + im[ua] im[ub]),

  Psi2i -> R4p[la,lb,lc,ld] ltet[ua] n[ud] (im[ub] rm[uc] - rm[ub] im[uc])
		 + Ro[la,lb,lc] nn ( n[uc] (im[ua] rm[ub] - rm[ua] im[ub]) - ltet[ub] (rm[ua] im[uc] - im[ua] rm[uc]) )
		 - Rojo[la,lb] nn nn (im[ua] rm[ub] - rm[ua] im[ub]),

  Psi1r -> R4p[la,lb,lc,ld] n[ua] ltet[ub] rm[uc] ltet[ud]
		 + Ro[la,lb,lc] ( nn ltet[ua] rm[ub] ltet[uc] - nn rm[ua] n[ub] ltet[uc] - nn n[ua] rm[ub] ltet[uc] )
		 + Rojo[la,lb] nn nn ( n[ua] rm[ub] - ltet[ua] rm[ub] ),

  Psi1i -> R4p[la,lb,lc,ld] n[ua] ltet[ub] im[uc] ltet[ud]
		 + Ro[la,lb,lc] ( nn ltet[ua] im[ub] ltet[uc] - nn im[ua] n[ub] ltet[uc] - nn n[ua] im[ub] ltet[uc] )
		 + Rojo[la,lb] nn nn ( n[ua] im[ub] - ltet[ua] im[ub] ),

  Psi0r -> R4p[la,lb,lc,ld] ltet[ua] ltet[uc] (rm[ub] rm[ud] - im[ub] im[ud])
		 + 2 Ro[la,lb,lc] nn ltet[ub] (rm[ua] rm[uc] - im[ua] im[uc])
		 + Rojo[la,lb] nn nn (rm[ua] rm[ub] - im[ua] im[ub]),

  Psi0i -> R4p[la,lb,lc,ld] ltet[ua] ltet[uc] (rm[ub] im[ud] + im[ub] rm[ud])
		 + 2 Ro[la,lb,lc] nn ltet[ub] (rm[ua] im[uc] + im[ua] rm[uc])
		 + Rojo[la,lb] nn nn (rm[ua] im[ub] + im[ua] rm[ub])
};

invariantEqs = {
  (* Scalar invariants I and J as defined in (2.2a) and (2.2b) of arXiv:gr-qc/0407013 *)
  curvIr -> ComplexExpand[Re[3 (Psi2r+I Psi2i)^2 - 4 (Psi1r+I Psi1i) (Psi3r + I Psi3i) + (Psi4r + I Psi4i) (Psi0r + I Psi0i)]],
  curvIi -> ComplexExpand[Im[3 (Psi2r+I Psi2i)^2 - 4 (Psi1r+I Psi1i) (Psi3r + I Psi3i) + (Psi4r + I Psi4i) (Psi0r + I Psi0i)]],
  curvJr -> ComplexExpand[Re[Det[{{Psi4r+I Psi4i,Psi3r+I Psi3i,Psi2r+I Psi2i},
                                  {Psi3r+I Psi3i,Psi2r+I Psi2i,Psi1r+I Psi1i},
                                  {Psi2r+I Psi2i,Psi1r+I Psi1i,Psi0r+I Psi0i}}]]],
  curvJi -> ComplexExpand[Im[Det[{{Psi4r+I Psi4i,Psi3r+I Psi3i,Psi2r+I Psi2i},
                                  {Psi3r+I Psi3i,Psi2r+I Psi2i,Psi1r+I Psi1i},
                                  {Psi2r+I Psi2i,Psi1r+I Psi1i,Psi0r+I Psi0i}}]]],

  (* Scalar invariants J1, J2, J3 and J4 of the Narlikar and Karmarkar basis as defined
     in B5-B8 of arXiv:0704.1756. Computed from Weyl tensor expressions using xTensor. *)
  curvJ1 -> -16(3 Psi2i^2-3 Psi2r^2-4 Psi1i Psi3i+4 Psi1r Psi3r+Psi0i Psi4i-Psi0r Psi4r),
  curvJ2 -> 96(-3 Psi2i^2 Psi2r+Psi2r^3+2 Psi1r Psi2i Psi3i+2 Psi1i Psi2r Psi3i-Psi0r Psi3i^2+2 Psi1i Psi2i Psi3r-2 Psi1r Psi2r Psi3r
    -2 Psi0i Psi3i Psi3r+Psi0r Psi3r^2-2 Psi1i Psi1r Psi4i+Psi0r Psi2i Psi4i+Psi0i Psi2r Psi4i-Psi1i^2 Psi4r+Psi1r^2 Psi4r
    +Psi0i Psi2i Psi4r-Psi0r Psi2r Psi4r),
  curvJ3 -> 64(9 Psi2i^4-54 Psi2i^2 Psi2r^2+9 Psi2r^4-24 Psi1i Psi2i^2 Psi3i+48 Psi1r Psi2i Psi2r Psi3i+24 Psi1i Psi2r^2 Psi3i
    +16 Psi1i^2 Psi3i^2-16 Psi1r^2 Psi3i^2+24 Psi1r Psi2i^2 Psi3r+48 Psi1i Psi2i Psi2r Psi3r-24 Psi1r Psi2r^2 Psi3r
    -64 Psi1i Psi1r Psi3i Psi3r-16 Psi1i^2 Psi3r^2+16 Psi1r^2 Psi3r^2+6 Psi0i Psi2i^2 Psi4i-12 Psi0r Psi2i Psi2r Psi4i
    -6 Psi0i Psi2r^2 Psi4i-8 Psi0i Psi1i Psi3i Psi4i+8 Psi0r Psi1r Psi3i Psi4i+8 Psi0r Psi1i Psi3r Psi4i
    +8 Psi0i Psi1r Psi3r Psi4i+Psi0i^2 Psi4i^2-Psi0r^2 Psi4i^2-6 Psi0r Psi2i^2 Psi4r-12 Psi0i Psi2i Psi2r Psi4r+6 Psi0r Psi2r^2 Psi4r
    +8 Psi0r Psi1i Psi3i Psi4r+8 Psi0i Psi1r Psi3i Psi4r+8 Psi0i Psi1i Psi3r Psi4r-8 Psi0r Psi1r Psi3r Psi4r-4 Psi0i Psi0r Psi4i Psi4r-Psi0i^2 Psi4r^2+Psi0r^2 Psi4r^2),
  curvJ4 -> -640(-15 Psi2i^4 Psi2r+30 Psi2i^2 Psi2r^3-3 Psi2r^5+10 Psi1r Psi2i^3 Psi3i+30 Psi1i Psi2i^2 Psi2r Psi3i-30 Psi1r Psi2i Psi2r^2 Psi3i
    -10 Psi1i Psi2r^3 Psi3i-16 Psi1i Psi1r Psi2i Psi3i^2-3 Psi0r Psi2i^2 Psi3i^2-8 Psi1i^2 Psi2r Psi3i^2+8 Psi1r^2 Psi2r Psi3i^2
    -6 Psi0i Psi2i Psi2r Psi3i^2+3 Psi0r Psi2r^2 Psi3i^2+4 Psi0r Psi1i Psi3i^3+4 Psi0i Psi1r Psi3i^3+10 Psi1i Psi2i^3 Psi3r
    -30 Psi1r Psi2i^2 Psi2r Psi3r-30 Psi1i Psi2i Psi2r^2 Psi3r+10 Psi1r Psi2r^3 Psi3r-16 Psi1i^2 Psi2i Psi3i Psi3r
    +16 Psi1r^2 Psi2i Psi3i Psi3r-6 Psi0i Psi2i^2 Psi3i Psi3r+32 Psi1i Psi1r Psi2r Psi3i Psi3r+12 Psi0r Psi2i Psi2r Psi3i Psi3r
    +6 Psi0i Psi2r^2 Psi3i Psi3r+12 Psi0i Psi1i Psi3i^2 Psi3r-12 Psi0r Psi1r Psi3i^2 Psi3r+16 Psi1i Psi1r Psi2i Psi3r^2
    +3 Psi0r Psi2i^2 Psi3r^2+8 Psi1i^2 Psi2r Psi3r^2-8 Psi1r^2 Psi2r Psi3r^2+6 Psi0i Psi2i Psi2r Psi3r^2-3 Psi0r Psi2r^2 Psi3r^2
    -12 Psi0r Psi1i Psi3i Psi3r^2-12 Psi0i Psi1r Psi3i Psi3r^2-4 Psi0i Psi1i Psi3r^3+4 Psi0r Psi1r Psi3r^3-6 Psi1i Psi1r Psi2i^2 Psi4i
    +2 Psi0r Psi2i^3 Psi4i-6 Psi1i^2 Psi2i Psi2r Psi4i+6 Psi1r^2 Psi2i Psi2r Psi4i+6 Psi0i Psi2i^2 Psi2r Psi4i
    +6 Psi1i Psi1r Psi2r^2 Psi4i-6 Psi0r Psi2i Psi2r^2 Psi4i-2 Psi0i Psi2r^3 Psi4i+12 Psi1i^2 Psi1r Psi3i Psi4i-4 Psi1r^3 Psi3i Psi4i
    -2 Psi0r Psi1i Psi2i Psi3i Psi4i-2 Psi0i Psi1r Psi2i Psi3i Psi4i-2 Psi0i Psi1i Psi2r Psi3i Psi4i
    +2 Psi0r Psi1r Psi2r Psi3i Psi4i-2 Psi0i Psi0r Psi3i^2 Psi4i+4 Psi1i^3 Psi3r Psi4i-12 Psi1i Psi1r^2 Psi3r Psi4i
    -2 Psi0i Psi1i Psi2i Psi3r Psi4i+2 Psi0r Psi1r Psi2i Psi3r Psi4i+2 Psi0r Psi1i Psi2r Psi3r Psi4i
    +2 Psi0i Psi1r Psi2r Psi3r Psi4i-2 Psi0i^2 Psi3i Psi3r Psi4i+2 Psi0r^2 Psi3i Psi3r Psi4i+2 Psi0i Psi0r Psi3r^2 Psi4i
    -Psi0r Psi1i^2 Psi4i^2-2 Psi0i Psi1i Psi1r Psi4i^2+Psi0r Psi1r^2 Psi4i^2+2 Psi0i Psi0r Psi2i Psi4i^2+Psi0i^2 Psi2r Psi4i^2
    -Psi0r^2 Psi2r Psi4i^2-3 Psi1i^2 Psi2i^2 Psi4r+3 Psi1r^2 Psi2i^2 Psi4r+2 Psi0i Psi2i^3 Psi4r+12 Psi1i Psi1r Psi2i Psi2r Psi4r
    -6 Psi0r Psi2i^2 Psi2r Psi4r+3 Psi1i^2 Psi2r^2 Psi4r-3 Psi1r^2 Psi2r^2 Psi4r-6 Psi0i Psi2i Psi2r^2 Psi4r+2 Psi0r Psi2r^3 Psi4r
    +4 Psi1i^3 Psi3i Psi4r-12 Psi1i Psi1r^2 Psi3i Psi4r-2 Psi0i Psi1i Psi2i Psi3i Psi4r+2 Psi0r Psi1r Psi2i Psi3i Psi4r
    +2 Psi0r Psi1i Psi2r Psi3i Psi4r+2 Psi0i Psi1r Psi2r Psi3i Psi4r-Psi0i^2 Psi3i^2 Psi4r+Psi0r^2 Psi3i^2 Psi4r
    -12 Psi1i^2 Psi1r Psi3r Psi4r+4 Psi1r^3 Psi3r Psi4r+2 Psi0r Psi1i Psi2i Psi3r Psi4r+2 Psi0i Psi1r Psi2i Psi3r Psi4r
    +2 Psi0i Psi1i Psi2r Psi3r Psi4r-2 Psi0r Psi1r Psi2r Psi3r Psi4r+4 Psi0i Psi0r Psi3i Psi3r Psi4r+Psi0i^2 Psi3r^2 Psi4r-Psi0r^2 Psi3r^2 Psi4r
    -2 Psi0i Psi1i^2 Psi4i Psi4r+4 Psi0r Psi1i Psi1r Psi4i Psi4r+2 Psi0i Psi1r^2 Psi4i Psi4r+2 Psi0i^2 Psi2i Psi4i Psi4r
    -2 Psi0r^2 Psi2i Psi4i Psi4r-4 Psi0i Psi0r Psi2r Psi4i Psi4r+Psi0r Psi1i^2 Psi4r^2+2 Psi0i Psi1i Psi1r Psi4r^2-Psi0r Psi1r^2 Psi4r^2
    -2 Psi0i Psi0r Psi2i Psi4r^2-Psi0i^2 Psi2r Psi4r^2+Psi0r^2 Psi2r Psi4r^2)
};

Psi4Calc[fdOrder_, PD_] :=
{
  Name -> "WeylScal4_psi4_calc_" <> fdOrder,
  Where -> Interior,
  After -> "ADMBase_SetADMVars",
  ConditionalOnKeywords -> {{"fd_order", fdOrder}, {"calc_scalars", "psi4"}},
  Shorthands -> shorthands,
  Equations -> psi4Eqs[PD]
};

PsisCalc[fdOrder_, PD_] :=
{
  Name -> "WeylScal4_psis_calc_" <> fdOrder,
  Where -> Interior,
  After -> "ADMBase_SetADMVars",
  ConditionalOnKeywords -> {{"fd_order", fdOrder}, {"calc_scalars", "psis"}},
  Shorthands -> shorthands,
  Equations -> Join[psi4Eqs[PD], otherPsiEqs]
};

InvariantsCalc[fdOrder_, PD_] :=
{
  Name -> "WeylScal4_invars_calc_" <> fdOrder,
  Where -> Everywhere,
  After -> "WeylScal4_psis_calc_" <> fdOrder <> "_group",
  ConditionalOnKeywords -> {{"fd_order", fdOrder}, {"calc_scalars", "psis"}, {"calc_invariants", "always"}},
  Equations -> invariantEqs
};


(****************************************************************************
  Construct the thorn
 ****************************************************************************)

fdOrderParam = 
{
  Name -> "fd_order",
  Default -> "Nth",
  AllowedValues -> {"Nth", "2nd", "4th"}
};

calcScalarsParam = {
  Name -> "calc_scalars",
  Description -> "Which scalars to calculate",
  AllowedValues -> {"psi4", "psis"},
  Default -> "psi4"
};

calcInvariantsParam = {
  Name -> "calc_invariants",
  Description -> "Compute invariants",
  AllowedValues -> {"always", "never"},
  Default -> "never"
};

keywordParameters = 
{
  fdOrderParam, calcScalarsParam, calcInvariantsParam
};

intParameters =
{
  {
    Name -> fdOrder,
    Default -> 2,
    AllowedValues -> {2,4,6,8}
  }
};


calculations = 
{
  Psi4Calc["Nth", PDstandard],
  Psi4Calc["2nd", PDstandard2nd],
  Psi4Calc["4th", PDstandard4th],
  PsisCalc["Nth", PDstandard],
  PsisCalc["2nd", PDstandard2nd],
  PsisCalc["4th", PDstandard4th],
  InvariantsCalc["Nth", PDstandard],
  InvariantsCalc["2nd", PDstandard2nd],
  InvariantsCalc["4th", PDstandard4th]
};

CreateKrancThornTT[groups, ".", "WeylScal4", 
  Calculations -> calculations,
  DeclaredGroups -> declaredGroupNames,
  PartialDerivatives -> derivatives,
  KeywordParameters -> keywordParameters,
  RealParameters -> realParameters,
  IntParameters -> intParameters,
  InheritedImplementations -> {"admbase", "methodoflines"},
  UseJacobian -> False,
  UseLoopControl -> False,
  UseVectors -> False];
