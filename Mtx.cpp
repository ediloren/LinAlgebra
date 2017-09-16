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


// mtx.cpp

// basic numerical matrix class
//

#include "stdafx.h"

#include "Mtx.h"

using namespace std;

//
// constructors
//

std::ostream& operator<<(std::ostream &s, const CLin_Matrix &A)
{
    CLin_subscript M;
    CLin_subscript N;

    M = A.num_rows();
    N = A.num_cols();

	s << M << " " << N << endl;

	for (CLin_subscript i=0; i<M; i++)
	{
		for (CLin_subscript j=0; j<N; j++)
		{
			s << A[i][j] << " ";
		}
		s << endl;
	}

	return s;
}

std::istream& operator>>(std::istream &s, CLin_Matrix &A)
{
	CLin_subscript M, N;

	s >> M >> N;

	if ( !(M == A.num_rows() && N == A.num_cols() ))
	{
		A.newsize(M,N);
	}

	for (CLin_subscript i=0; i<M; i++)
		for (CLin_subscript j=0; j<N; j++)
		{
			s >>  A[i][j];
		}

	return s;
}



