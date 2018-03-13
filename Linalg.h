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


// LinAlg header file


#ifndef LINALG_H
#define LINALG_H

#include <cassert>
#include <iostream>

#ifdef MS_VS
#include <crtdbg.h>
#else
#  include <wx/debug.h>
#  define _ASSERT wxASSERT
#endif

#define LINALG_MAJOR_VERSION    '0'
#define LINALG_MINOR_VERSION    '1'
#define LINALG_SUBMINOR_VERSION '0'
#define LINALG_VERSION_STRING "0.1.0"


// Define the data type used for matrix and vector Subscripts.
#ifndef LINALG_SUBSCRIPT_TYPE
#define LINALG_SUBSCRIPT_TYPE unsigned long
#endif

typedef LINALG_SUBSCRIPT_TYPE CLin_subscript;

// Define this macro if you want LINALG to ensure all references
// are within the bounds of the array.  This encurs a run-time
// overhead, of course, but is recommended while developing
// code.  It can be turned off for production runs.
//
//       #define LINALG_BOUNDS_CHECK

#define LINALG_BOUNDS_CHECK
#ifdef LINALG_NO_BOUNDS_CHECK
#undef LINALG_BOUNDS_CHECK
#endif

#include <exception>
using namespace std;

// TBC warning: assertion routine, to be customized
//#define CLin_assert assert
#ifdef _DEBUG
inline void CLin_assert(bool a) { _ASSERT(a); }
#else
inline void CLin_assert(bool ) { }
#endif


#endif // LINALG_H

