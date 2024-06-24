
# Multipole

|Authors          |Ian Hinder and Andrew Knapp (modernized by Liwei Ji, Samuel Cupp, Roland Haas, Allen Wen, Yosef Zlochower) |
|:----------------|:---------------------------|
|Maintainer       |Liwei Ji
|Licence          |GNU GPL version 2
|Documentation    |http://einsteintoolkit.org/thornguide/EinsteinAnalysis/Multipole/documentation.html


## Purpose of This Thorn

The Multipole thorn performs spherical harmonic mode decomposition of
Cactus grid functions on coordinate spheres.  It can decompose
multiple grid functions with any spin-weight on multiple spheres.
This thorn uses the interpolator interface to access grid functions,
so works with mesh-refinement and multi-patch.

## Required thorns

* An interpolator

## Related thorns

* WeylScal4 can be used to compute the Weyl scalars, which can then be
  decomposed into modes on coordinate spheres by Multipole.

## Publications

Multipole has been used in the following publications:

## Copyright

This thorn is copyright (C) 2007-2011 (C) by Ian Hinder and Andrew
Knapp.

This thorn is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This thorn is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this thorn (see the file COPYING in this directory);
if not, write to the Free Software Foundation, Inc., 59 Temple
Place, Suite 330, Boston, MA  02111-1307  USA
