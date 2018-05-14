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


// range.cpp header file

// basic vector 1D range class
//

#ifndef RANGE_H
#define RANGE_H

#include "Linalg.h"
#include "Range.h"
#include <math.h>

class CLin_Range 
{

public:
    typedef         double*  iterator;
    typedef         double&  reference;
    typedef const   double&  const_reference;

	//
	// constructors
	//

	// default constructor
	inline CLin_Range() : m_dV(0), m_iN(0)  
	{
		// vectors default to dimension 3, to use in 3D representations 
		initialize(3);
		// do not zero by default
		//set(0.0);
	}

	// construct from another already-existing vector 
	inline CLin_Vector(const CLin_Vector &pVect) : m_dV(0), m_iN(0)
	{
		initialize(pVect.m_iN);
		copy(pVect.m_dV);
	}

	// construct vector of lenght 'N'  
	inline CLin_Vector(CLin_subscript N) :  m_dV(0), m_iN(0)
	{
		initialize(N);
	}

	// construct vector of lenght 'N' and fill it up with value 'value'  
	inline CLin_Vector(CLin_subscript N, double value = 0.0) :  m_dV(0), m_iN(0)
	{
		initialize(N);
		set(value);
	}

	// construct vector of lenght 'N' and initialize with values from array of doubles 'v' 
	inline CLin_Vector(CLin_subscript N, const double* v) :  m_dV(0), m_iN(0)
	{
		initialize(N);
		copy(v);
	}

	//
	// destructor
	//

	inline ~CLin_Vector() 
	{
		destroy();
	}

	//
	// access
	//

    inline iterator begin() { return m_dV;}
    inline iterator end()   { return m_dV + m_iN; }
    inline CLin_subscript dim() const { return  m_iN; }
    inline CLin_subscript size() const { return  m_iN; }

	//
	// methods
	//

	inline CLin_Vector& newsize(CLin_subscript N)
	{
		if (m_iN == N) return *this;
		
		destroy();
		initialize(N);
		
		return *this;
	}

	//
    // operators
    //

	inline CLin_Vector& operator=(const CLin_Vector &A)
	{
		if (m_dV == A.m_dV)
			return *this;
		
		// no need to re-alloc
		if (m_iN == A.m_iN)         
			copy(A.m_dV);
		else {
			destroy();
			initialize(A.m_iN);
			copy(A.m_dV);
		}
		
		return *this;
	}
	
	inline CLin_Vector& operator=(const double& scalar)
	{ 
		set(scalar);  
		return *this;
	}

    inline reference operator[](CLin_subscript i)
    { 
#ifdef LINALG_BOUNDS_CHECK
        CLin_assert(0<=i);
        CLin_assert(i < m_iN) ;
#endif
        return m_dV[i]; 
    }

    inline const_reference operator[](CLin_subscript i) const
    {
#ifdef LINALG_BOUNDS_CHECK
        CLin_assert(0<=i);
        CLin_assert(i < m_iN) ;
#endif
        return m_dV[i]; 
    }


protected:
	// initialize a new vector (that is, allocate memory and initialize members)
	inline void initialize(CLin_subscript N)
	{
		// allocate array
		CLin_assert(m_dV == NULL);
		m_dV = new double[N];
		CLin_assert(m_dV  != NULL);
		// init lenght member
		m_iN = N;
	}

	// copy from an array into the vector
	inline void copy(const double* v)
	{
		CLin_subscript i;
		
		for (i=0; i<m_iN; i++)
			m_dV[i] = v[i];
	}

	// fill array with value 'val'
	inline void set(const double val)
	{
		CLin_subscript i;
		
		for (i=0; i<m_iN; i++)
			m_dV[i] = val;    
	}
	
	inline void destroy()
	{     
		// do nothing, if no memory has been previously allocated 
		if (m_dV == NULL) return ;
		
		delete [] (m_dV);     
		
		m_dV = NULL;
	}

    double* m_dV;                  
    CLin_subscript m_iN;
};

//
//  I/O  
//

std::ostream& operator<<(std::ostream &s, const CLin_Vector &A);
std::istream& operator>>(std::istream &s, CLin_Vector &A);

//
// basic vector operations
//

inline CLin_Vector operator+(const CLin_Vector &A, const CLin_Vector &B)
{
    CLin_subscript N = A.dim();

    CLin_assert(N==B.dim());

    CLin_Vector tmp(N);
    CLin_subscript i;

    for (i=0; i<N; i++)
		tmp[i] = A[i] + B[i];

    return tmp;
}

inline CLin_Vector operator-(const CLin_Vector &A, const CLin_Vector &B)
{
    CLin_subscript N = A.dim();

    CLin_assert(N==B.dim());

    CLin_Vector tmp(N);
    CLin_subscript i;

    for (i=0; i<N; i++)
		tmp[i] = A[i] - B[i];

    return tmp;
}

inline CLin_Vector operator*(const CLin_Vector &A, const CLin_Vector &B)
{
    CLin_subscript N = A.dim();

    CLin_assert(N==B.dim());

    CLin_Vector tmp(N);
    CLin_subscript i;

    for (i=0; i<N; i++)
		tmp[i] = A[i] * B[i];

    return tmp;
}

inline CLin_Vector operator*(CLin_Vector A, const double b)
{
    CLin_subscript N = A.dim();
    CLin_subscript i;

    for (i=0; i<N; i++)
		A[i] *= b;

    return A;
}

inline CLin_Vector operator*(const double b, CLin_Vector A)
{
    CLin_subscript N = A.dim();
    CLin_subscript i;

    for (i=0; i<N; i++)
		A[i] *= b;

    return A;
}

inline CLin_Vector operator/(const CLin_Vector &A, const double &B)
{
    CLin_Vector tmp(A);
    CLin_subscript i;
	CLin_subscript N = A.dim();

    for (i=0; i<N; i++)
		tmp[i] = A[i] / B;

    return tmp;
}

inline double dot_prod(const CLin_Vector &A, const CLin_Vector &B)
{
    CLin_subscript N = A.dim();
    CLin_assert(N == B.dim());

    CLin_subscript i;
    double sum = 0;

    for (i=0; i<N; i++)
		sum += A[i] * B[i];

    return sum;
}

inline double mod(const CLin_Vector &A)
{
	return sqrt(dot_prod(A,A));
}

// Cross product of two 3D vectors - should be much faster
// than a general-purpose 'cross' function, because is specialized
// for the three dimensions
inline CLin_Vector cross3d(const CLin_Vector &A, const CLin_Vector &B) 
{
	CLin_subscript N = A.dim();
	CLin_subscript M = B.dim();

    CLin_assert(N == 3);
    CLin_assert(M == 3);

	CLin_Vector tmp(N);

	tmp[0] = A[1]*B[2]-A[2]*B[1];
	tmp[1] = -A[0]*B[2]+A[2]*B[0];
	tmp[2] = A[0]*B[1]-A[1]*B[0];

	return tmp;
}


#endif // RANGE_H
