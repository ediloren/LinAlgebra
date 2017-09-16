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


// mtx.cpp header file

// basic numerical matrix class
//

// 2010/07/12 update
// added Frobenius norm
// 2010/08/18 update
// corrected memory leak bug in destroy(): wrong 'if' statement

#ifndef MTX_H
#define MTX_H

#include "Linalg.h"
#include "Vect.h"

#include <iostream>

class CLin_Matrix
{
public:

	//
	// constructors
	//

	// default constructor
	inline CLin_Matrix() : m_iM(0), m_iN(0), m_iMN(0), m_pV(0), m_pRow(0)
	{
	}

	// construct from another already-existing matrix
	inline CLin_Matrix(const CLin_Matrix &A)
	{
		initialize(A.m_iM, A.m_iN);
		copy(A.m_pV);
	}

	// construct matrix of dim 'M' x 'N' and fill it up with value 'value'
	inline CLin_Matrix(CLin_subscript M, CLin_subscript N, double value = 0.0)
	{
		initialize(M,N);
		set(value);
	}

	// construct matrix of dim 'M' x 'N' and initialize with values from matrix of doubles 'v'
	inline CLin_Matrix(CLin_subscript M, CLin_subscript N, const double* v)
	{
		initialize(M,N);
		copy(v);
	}

	//
	// destructor
	//
	inline ~CLin_Matrix()
	{
		destroy();
	}

	//
	// access
	//

	inline CLin_subscript size() const { return m_iMN; }

	// returns # of rows or cols according to parameter d (1=rows, 2=cols, other=0)
	inline CLin_subscript dim(CLin_subscript d)
	{
#ifdef LINALG_BOUNDS_CHECK
		CLin_assert( d >= 1);
		CLin_assert( d <= 2);
#endif
		if(d==1) {
			return m_iM;
		}
		else if(d==2) {
			return m_iN;
		}
		else {
			return 0;
		}
	}

	inline CLin_subscript num_rows() const { return m_iM; }

	inline CLin_subscript num_cols() const { return m_iN; }

	//
	// methods
	//

	inline CLin_Matrix& newsize(CLin_subscript M, CLin_subscript N)
	{
		if (num_rows() == M && num_cols() == N)
			return *this;

		destroy();
		initialize(M,N);

		return *this;
	}

	//
	// operators
	//

	inline CLin_Matrix& operator=(const CLin_Matrix &A)
	{
		if (m_pV == A.m_pV)
			return *this;

		// no need to re-alloc
		if (m_iM == A.m_iM  && m_iN == A.m_iN)
			copy(A.m_pV);
		else
		{
			destroy();
			initialize(A.m_iM, A.m_iN);
			copy(A.m_pV);
		}

		return *this;
	}

	inline CLin_Matrix& operator=(const double scalar)
	{
		set(scalar);
		return *this;
	}

	inline operator double**(){ return	m_pRow; }

	inline operator double**() const { return m_pRow; }

	inline double* operator[](CLin_subscript i)
	{
#ifdef LINALG_BOUNDS_CHECK
//		CLin_assert(0<=i);
		CLin_assert(i < m_iM) ;
#endif
		return m_pRow[i];
	}
/*
	inline const double* operator[](CLin_subscript i) const
	{
#ifdef LINALG_BOUNDS_CHECK
		CLin_assert(0<=i);
		CLin_assert(i < m_iM) ;
#endif
		return m_pRow[i];
	}
*/
	inline double operator()(CLin_subscript i, CLin_subscript j)
	{
#ifdef LINALG_BOUNDS_CHECK
//		CLin_assert(i>=0);
		CLin_assert(i < m_iM) ;
//		CLin_assert(j>=0);
		CLin_assert(j < m_iN);
#endif
		return	m_pRow[i][j];
	}


#ifdef LINALG_USE_REGIONS

	typedef Region2D<CLin_Matrix<T> > Region;


	Region operator()(const Index1D &I, const Index1D &J)
	{
		return Region(*this, I,J);
	}


	typedef const_Region2D< CLin_Matrix<T> > const_Region;
	const_Region operator()(const Index1D &I, const Index1D &J) const
	{
		return const_Region(*this, I,J);
	}

#endif

	inline void destroy()
	{
		// do nothing, if no memory has been previously allocated
		if (m_pV != NULL) delete [] (m_pV);
		if (m_pRow != NULL) delete [] (m_pRow);
		m_iM = 0;
		m_iN = 0;
		m_iMN = 0;
		m_pV = NULL;
		m_pRow = NULL;
	}

protected:
	// initialize a new matrix (that is, allocate memory and initialize members)
	inline void initialize(CLin_subscript M, CLin_subscript N)
	{
		m_iMN = M*N;
		m_iM = M;
		m_iN = N;

		// allocate memory
		m_pV = new double[m_iMN];
		m_pRow = new double*[M];

		CLin_assert(m_pV  != NULL);
		CLin_assert(m_pRow	!= NULL);

		// initialize vector of pointers to matrix rows
		double* p = m_pV;
		for (CLin_subscript i=0; i<M; i++)
		{
			m_pRow[i] = p;
			p += N ;
		}
	}

	// copy from an array into the matrix
	inline void copy(const double*  v)
	{
		CLin_subscript i;

		for (i=0; i < m_iM * m_iN; i++)
			m_pV[i] = v[i];
	}

	// fill matrix with value 'val'
	inline void set(const double val)
	{
		CLin_subscript i;

		for (i=0; i < m_iM * m_iN; i++)
			m_pV[i] = val;
	}

	// rows
	CLin_subscript m_iM;
	// cols
	CLin_subscript m_iN;
	// total size
	CLin_subscript m_iMN;
	// the actual matrix
	double* m_pV;
	// array of pointers to matrix rows
	double** m_pRow;
};


//	I/O

std::ostream& operator<<(std::ostream &s, const CLin_Matrix &A);
std::istream& operator>>(std::istream &s, CLin_Matrix &A);

// basic matrix operations

inline CLin_Matrix operator+(const CLin_Matrix &A, const CLin_Matrix &B)
{
	CLin_subscript M, N;

	M = A.num_rows();
	N = A.num_cols();

	CLin_assert(M==B.num_rows());
	CLin_assert(N==B.num_cols());

	CLin_Matrix tmp(M,N);
	CLin_subscript i,j;

	for (i=0; i<M; i++)
		for (j=0; j<N; j++)
			tmp[i][j] = A[i][j] + B[i][j];

	return tmp;
}

inline CLin_Matrix operator-(const CLin_Matrix &A, const CLin_Matrix &B)
{
	CLin_subscript M, N;

	M = A.num_rows();
	N = A.num_cols();

	CLin_assert(M==B.num_rows());
	CLin_assert(N==B.num_cols());

	CLin_Matrix tmp(M,N);
	CLin_subscript i,j;

	for (i=0; i<M; i++)
		for (j=0; j<N; j++)
			tmp[i][j] = A[i][j] - B[i][j];

	return tmp;
}

inline CLin_Matrix mult_element(const CLin_Matrix &A, const CLin_Matrix &B)
{
	CLin_subscript M, N;

	M = A.num_rows();
	N = A.num_cols();

	CLin_assert(M==B.num_rows());
	CLin_assert(N==B.num_cols());

	CLin_Matrix tmp(M,N);
	CLin_subscript i,j;

	for (i=0; i<M; i++)
		for (j=0; j<N; j++)
			tmp[i][j] = A[i][j] * B[i][j];

	return tmp;
}

inline CLin_Matrix transpose(const CLin_Matrix &A)
{
	CLin_subscript M, N;

	M = A.num_rows();
	N = A.num_cols();

	CLin_Matrix S(N,M);
	CLin_subscript i, j;

	for (i=0; i<M; i++)
		for (j=0; j<N; j++)
			S[j][i] = A[i][j];

	return S;
}

inline CLin_Vector matmult(const CLin_Matrix  &A, const CLin_Vector &x)
{

#ifdef LINALG_BOUNDS_CHECK
	CLin_assert(A.num_cols() == x.dim());
#endif

	CLin_subscript M, N;
	CLin_subscript i, j;

	M = A.num_rows();
	N = A.num_cols();

	CLin_Vector tmp(M);
	double sum;

	for (i=0; i<M; i++)
	{
		sum = 0;
		const double* rowi = A[i];
		for (j=0; j<N; j++)
			sum = sum +  rowi[j] * x[j];

		tmp[i] = sum;
	}

	return tmp;
}

inline CLin_Matrix matmult(const CLin_Matrix &A, const CLin_Matrix &B)
{

#ifdef LINALG_BOUNDS_CHECK
	CLin_assert(A.num_cols() == B.num_rows());
#endif

	CLin_subscript M, N, K;
	CLin_subscript i, j, k;
	double sum;

	M = A.num_rows();
	N = A.num_cols();
	K = B.num_cols();

	CLin_Matrix tmp(M,K);

	for (i=0; i<M; i++)
		for (k=0; k<K; k++)
		{
			sum = 0;
			for (j=0; j<N; j++)
				sum = sum +  A[i][j] * B[j][k];

			tmp[i][k] = sum;
		}

	return tmp;
}

inline int matmult(CLin_Matrix &C, const CLin_Matrix &A, const CLin_Matrix &B)
{
	CLin_subscript M, N, K;
	CLin_subscript i, j, k;
	double sum;
	const double* row_i;
	const double* col_k;

	CLin_assert(A.num_cols() == B.num_rows());

	M = A.num_rows();
	N = A.num_cols();
	K = B.num_cols();

	C.newsize(M,K);

	for (i=0; i<M; i++)
		for (k=0; k<K; k++)
		{
			row_i  = &(A[i][0]);
			col_k  = &(B[0][k]);
			sum = 0;
			for (j=0; j<N; j++)
			{
				sum  += *row_i * *col_k;
				row_i++;
				col_k += K;
			}
			C[i][k] = sum;
		}

	return 0;
}

inline CLin_Vector operator*(const CLin_Matrix	&A, const CLin_Vector &x)
{
	return matmult(A,x);
}

inline CLin_Matrix operator*(const CLin_Matrix &A, const CLin_Matrix &B)
{
	return matmult(A,B);
}

// Frobenius norm
inline double FNorm(const CLin_Matrix &A)
{
	CLin_subscript M, N;
	CLin_subscript i, j;
	const double* rowi;
	double norm;

	M = A.num_rows();
	N = A.num_cols();

	norm = 0;
	for (i=0; i<M; i++) {
		rowi = A[i];
		for (j=0; j<N; j++) {
			norm += rowi[j]*rowi[j];
		}
	}

	norm = sqrt(norm);

	return norm;
}

// Frobenius norm, complex matrix
inline double FNorm(const CLin_Matrix &ARe, const CLin_Matrix &AIm)
{
	CLin_subscript M, N;
	CLin_subscript i, j;
	const double *rowiRe, *rowiIm;
	double norm;

	M = ARe.num_rows();
	N = ARe.num_cols();

	CLin_assert(AIm.num_rows() == M);
	CLin_assert(AIm.num_cols() == N);

	if(AIm.num_rows() != M || AIm.num_cols() != N) {
		// impossible value to signal an issue (norm should be positive)
		return -1.0;
	}

	norm = 0;
	for (i=0; i<M; i++) {
		rowiRe = ARe[i];
		rowiIm = AIm[i];
		for (j=0; j<N; j++) {
			norm += rowiRe[j]*rowiRe[j] + rowiIm[j]*rowiIm[j];
		}
	}

	norm = sqrt(norm);

	return norm;
}


#endif // MTX_H
