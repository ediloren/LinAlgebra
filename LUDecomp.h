/***************************************************************************
*                                                                          *
*   Copyright (c) 2017                                                     *
*   FastFieldSolvers S.R.L.  http://www.fastfieldsolvers.com               *
*                                                                          *
*   This program is free software; you can redistribute it and/or modify   *
*   it under the terms of the GNU Lesser General Public License (LGPL)     *
*   as published by the Free Software Foundation; either version 2 of      *
*   the License, or (at your option) any later version.                    *
*   for detail see the LICENCE text file.                                  *
*                                                                          *
*   This program is distributed in the hope that it will be useful,        *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
*   GNU Library General Public License for more details.                   *
*                                                                          *
*   You should have received a copy of the GNU Library General Public      *
*   License along with this program; if not, write to the Free Software    *
*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307   *
*   USA                                                                    *
*                                                                          *
***************************************************************************/


// ludecomp.cpp header file

// matrix LU decomposition class
//

#ifndef LUDECOMP_H
#define LUDECOMP_H

#include "Linalg.h"
#include "Vect.h"
#include "Mtx.h"

// Solve system of linear equations Ax = b.
//
//  Typical usage:
//
//    CLin_Matrix A;
//    CLin_Vector ipiv;
//    CLin_Vector b;
//
//    1)  CLin_LU_Factor(A,ipiv);
//    2)  CLin_LU_Solve(A,ipiv,b);
//
//   Now b has the solution x.  Note that both A and b
//   are overwritten.  If these values need to be preserved, 
//   one can make temporary copies 


int CLin_LU_factor( CLin_Matrix &A, CLin_Vector &indx);
int CLin_LU_solve(const CLin_Matrix &A, const CLin_Vector &indx, CLin_Vector &b);

#endif // LUDECOMP_H
