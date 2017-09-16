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


// vect.cpp

// basic numerical vector class
//

#include "stdafx.h"

#include "Vect.h"

using namespace std;

// the items between 's' and 'e', included, are referenced by the range object
CLin_Range::CLin_Range(CLin_Vector &A, CLin_subscript s, CLin_subscript e)
{
	CLin_assert(A.m_dV != NULL);
	CLin_assert(s <= A.m_iN);
	CLin_assert(e <= A.m_iN);

	// point to the start of the range inside the vector
	m_dV = A.m_dV + s;
	// compute new length
	m_iN = e-s+1;
}

// the items between 's' and the end of the array, included, are referenced by the range object
CLin_Range::CLin_Range(CLin_Vector &A, CLin_subscript s)
{
	CLin_assert(A.m_dV != NULL);
	CLin_assert(s <= A.m_iN);

	// point to the start of the range inside the vector
	m_dV = A.m_dV + s;
	// compute new length
	m_iN = A.m_iN-s;
}

CLin_Range& CLin_Range::operator=( CLin_Vector &A)
{
	CLin_assert(A.m_dV != NULL);

	m_dV = A.m_dV;
	m_iN = A.m_iN;

	return *this;
}

//  I/O

std::ostream& operator<<(std::ostream &s, const CLin_Vector &A)
{
    CLin_subscript N=A.dim();

    s <<  N << endl;

    for (CLin_subscript i=0; i<N; i++)
        s   << A[i] << " " << endl;
    s << endl;

    return s;
}

std::istream & operator>>(std::istream &s, CLin_Vector &A)
{
    CLin_subscript N;

    s >> N;

    if ( !(N == A.size() ))
    {
        A.newsize(N);
    }

    for (CLin_subscript i=0; i<N; i++)
            s >>  A[i];

    return s;
}


