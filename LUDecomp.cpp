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

#include "stdafx.h"

#include <cmath>

#include "LUDecomp.h"

// right-looking LU factorization algorithm (unblocked)
//
//   Factors matrix A into lower and upper  triangular matrices 
//   (L and U respectively) in solving the linear equation Ax=b.  
//
//
// Args:
//
// A        (input/output) Matrix(1:n, 1:n)  In input, matrix to be
//                  factored.  On output, overwritten with lower and 
//                  upper triangular factors.
//
// indx     (output) Vector(1:n)    Pivot vector. Describes how
//                  the rows of A were reordered to increase
//                  numerical stability.
//
// Return value:
//
// int      (0 if successful, 1 otherwise)
//
//

int CLin_LU_factor( CLin_Matrix &A, CLin_Vector &indx)
{
	CLin_subscript M, N, minMN;
	CLin_subscript i, j, k, jp;
	CLin_subscript ii,jj;
	double t, recp;

	M = A.num_rows();
	N = A.num_cols();

	if (M == 0 || N == 0) return 0;

	if (indx.dim() != M)
		indx.newsize(M);

	i = j = k = jp =0;

	 // min(M,N);
	minMN =  (M < N ? M : N) ;       

	for (j=0; j < minMN; j++)
	{
        // find pivot in column j and test for singularity.

        jp = j;
        t = fabs(A[j][j]);
        for (i=j+1; i<M; i++)
            if ( fabs(A[i][j]) > t)
            {
                jp = i;
                t = fabs(A[i][j]);
            }

        indx[j] = jp;

        // jp now has the index of maximum element 
        // of column j, below the diagonal

        if ( A[jp][j] == 0 )        
			// factorization failed because of zero pivot
            return 1;       


		// if pivot not already on the diagonal 
        if (jp != j)  
			// swap rows j and jp
            for (k=0; k<N; k++)
            {
                t = A[j][k];
                A[j][k] = A[jp][k];
                A[jp][k] = t;
            }
		
		// divide elements j+1:M of jth column by pivot element 
		// (for the first M-1 cols, last col's below-the-diagonal element
		// is already the pivot)
        if (j<M-1)                
        {
            // note that A(j,j) was previously A(jp,p), which was
            // guaranteed not to be zero
            recp =  1.0 / A[j][j];

            for (k=j+1; k<M; k++)
                A[k][j] *= recp;
        }


        if (j < minMN-1)
        {
            // rank-1 update to trailing submatrix:   E = E - x*y;
            //
            // E is the region A(j+1:M, j+1:N)
            // x is the column vector A(j+1:M,j)
            // y is row vector A(j,j+1:N)

            for (ii=j+1; ii<M; ii++)
                for (jj=j+1; jj<N; jj++)
                    A[ii][jj] -= A[ii][j]*A[j][jj];
        }
    }

    return 0;
}   

int CLin_LU_solve(const CLin_Matrix &A, const CLin_Vector &indx, CLin_Vector &b)
{

    CLin_subscript i, ii, ip, j;
    CLin_subscript n;
	double sum;
	bool flag;

	n = b.dim();
	flag = false;
    sum = 0.0;

    for (i=0; i<n; i++) 
    {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];

        if (flag)
            for (j=ii; j<=i-1; j++) 
                sum -= A[i][j]*b[j];
		else if (sum) {
			ii = i;
			flag = true;
		}

        b[i] = sum;
    }
    for (i=n-1;i>=0;i--) 
    {
        sum=b[i];
        for (j=i+1;j<n;j++) 
            sum -= A[i][j]*b[j];
        b[i] = sum/A[i][i];
    }

    return 0;
}
