# ellipsoid.maple -- compute equations for offset ellipsoid setup
# $Header: /usr/local/svn/cvs-repositories/numrelcvs/AEIThorns/AHFinderDirect/archive/ellipsoid.mm,v 1.1 2002-07-18 17:37:11 jthorn Exp $

#
# ellipsoid has center (A,B,C), radius (a,b,c)
# angular coordinate system has center (U,V,W)
#
# direction cosines wrt angular coordinate center are (alpha,beta,gamma)
# but Maple predefines gamma = Euler's constant, so we use (xcos,ycos,zcos)
# instead, i.e. a point has coordinates (U+xcos*r, V+ycos*r, W+zcos*r)
#
# then the equation of the ellipsoid is
#	(U+xcos*r - A)^2     (V+ycos*r - B)^2     (W+zcos*r - C)^2
#	-----------------  +  ----------------  +  -----------------  =  1
#	        a^2                  b^2                   c^2
#
# to solve this, we introduce intermediate variables
#	AU = A - U
#	BV = B - V
#	CW = C - W
#
eqn := (xcos*r - AU)^2/a^2 + (ycos*r - BV)^2/b^2 + (zcos*r - CW)^2/c^2 = 1;

read "../maple/util.mm";
read "../maple/codegen2.mm";

[solve(eqn, r)];
map(simplify, %);
[r_plus = %[1], r_minus = %[2]];
solnlist := [codegen[optimize](%)];

ftruncate("ellipsoid.c");
print_name_list_dcl(temps_in_eqnlist(solnlist, [r_plus,r_minus]),
		    "fp", "ellipsoid.c");
codegen[C](solnlist, filename="ellipsoid.c");
