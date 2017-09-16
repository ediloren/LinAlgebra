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


// vect.cpp header file

// basic numerical vector class
//

// Update 2010/07/15
// modified destroy() zeroing m_iN

#ifndef VECT_H
#define VECT_H

#ifdef MS_VS
// when using MS VisualC++
//#include "afx.h"
#endif

#include "Linalg.h"
#include <iostream>
#include <math.h>

class CLin_Vector;

class CLin_Range
{

public:
	friend class CLin_Vector;

	typedef         double&  reference;
    typedef const   double&  const_reference;

	//
	// constructors
	//

	// default constructor
	inline CLin_Range() : m_dV(NULL), m_iN(0) {}
	// other constructor;
	// cannot be defined inline because the function body is in a .cpp file;
	// and cannot be defined in the header because accesses members of the frienc class that are
	// not yet defined (see CLin_Vector definition)
	CLin_Range(CLin_Vector &v, CLin_subscript s, CLin_subscript e);
	CLin_Range(CLin_Vector &v, CLin_subscript s);

    inline reference operator[](CLin_subscript i)
    {
#ifdef LINALG_BOUNDS_CHECK
//        CLin_assert(0 <= i);
        CLin_assert(i < m_iN);
#endif
		CLin_assert(m_dV != NULL);

       return m_dV[i];
    }

    inline const_reference operator[](CLin_subscript i) const
    {
#ifdef LINALG_BOUNDS_CHECK
//        CLin_assert(0 <= i);
        CLin_assert(i < m_iN);
#endif
		CLin_assert(m_dV != NULL);

        return m_dV[i];
    }

	inline CLin_Range& operator=(const CLin_Range &A)
	{
		CLin_assert(A.m_dV != NULL);

		m_dV = A.m_dV;
		m_iN = A.m_iN;

		return *this;
	}

	// cannot be defined inline because the function body is in a .cpp file;
	// and cannot be defined in the header because accesses members of the friend class that are
	// not yet defined (see CLin_Vector definition)
	CLin_Range& operator=( CLin_Vector &A);

    inline CLin_subscript dim() const { return  m_iN; }
    inline CLin_subscript size() const { return  m_iN; }
    inline double* array() { return m_dV; }

protected:
    double* m_dV;
    CLin_subscript m_iN;
};


class CLin_Vector
{

public:
	friend class CLin_Range;

    typedef         double*  iterator;
    typedef         double&  reference;
    typedef const   double&  const_reference;

	//
	// constructors
	//

	// default constructor
	inline CLin_Vector() : m_dV(NULL), m_iN(0)
	{
		// vectors default to dimension 3, to use in 3D representations
		//initialize(3);
		// do not zero by default
		//set(0.0);
	}

	// construct from another already-existing vector
	inline CLin_Vector(const CLin_Vector &pVect) : m_dV(NULL), m_iN(0)
	{
		initialize(pVect.m_iN);
		copy(pVect.m_dV);
	}

	// construct vector of lenght 'N'
	inline CLin_Vector(CLin_subscript N) :  m_dV(NULL), m_iN(0)
	{
		initialize(N);
	}

	// construct vector of lenght 'N' and fill it up with value 'value'
	inline CLin_Vector(CLin_subscript N, double value) :  m_dV(NULL), m_iN(0)
	{
		initialize(N);
		set(value);
	}

	// construct vector of lenght 'N' and initialize with values from array of doubles 'v'
	inline CLin_Vector(CLin_subscript N, const double* v) :  m_dV(NULL), m_iN(0)
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
    inline double* array() { return m_dV; }

	//
	// methods
	//

	inline bool newsize(CLin_subscript N)
	{
		bool ret;

		ret = true;

		if (m_iN != N) {
			destroy();
			ret = initialize(N);
		}

		return ret;
	}

	//
    // operators
    //

	inline CLin_Vector& operator=(const CLin_Vector &A)
	{
		if (m_dV == A.m_dV)
			return *this;

		// need to re-alloc
		if (m_iN != A.m_iN) {
			destroy();
			initialize(A.m_iN);
		}

		copy(A.m_dV);

		return *this;
	}

	inline CLin_Vector& operator=(const double& scalar)
	{
		set(scalar);
		return *this;
	}


	inline CLin_Vector& operator=(const CLin_Range &R)
	{
		// need to re-alloc
		if (m_dV == NULL || m_iN != R.m_iN) {
			destroy();
			initialize(R.m_iN);
		}

		copy(R.m_dV);

		return *this;
	}


    inline reference operator[](CLin_subscript i)
    {
#ifdef LINALG_BOUNDS_CHECK
//        CLin_assert(0<=i);
        CLin_assert(i < m_iN) ;
#endif
        return m_dV[i];
    }

    inline const_reference operator[](CLin_subscript i) const
    {
#ifdef LINALG_BOUNDS_CHECK
//        CLin_assert(0<=i);
        CLin_assert(i < m_iN) ;
#endif
        return m_dV[i];
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
		m_iN = 0;
	}

protected:
	// initialize a new vector (that is, allocate memory and initialize members)
	inline bool initialize(CLin_subscript N)
	{
		// allocate array
		CLin_assert(m_dV == NULL);
		if(m_dV != NULL) {
			return false;
		}
		try {
			m_dV = new double[N];
		}
		catch( bad_alloc& ) {
			m_iN = 0;
			return false;
		}
		catch(...) {
			m_iN = 0;
			return false;
		}
		if(m_dV == NULL) {
			m_iN = 0;
			return false;
		}
		// init length member
		m_iN = N;

		return true;
	}

	// copy from an array into the vector
	inline void copy(const double* v)
	{
		CLin_subscript i;

		for (i=0; i<m_iN; i++)
			m_dV[i] = v[i];
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



#endif // VECT_H
