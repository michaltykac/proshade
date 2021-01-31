/***************************************************************************
  **************************************************************************
  
  S2kit 1.0
  A lite version of Spherical Harmonic Transform Kit

  Copyright (c) 2004 Peter Kostelec, Dan Rockmore

  This file is part of S2kit.

  S2kit is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  S2kit is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
  
  See the accompanying LICENSE file for details.
  
  ************************************************************************
  ************************************************************************/

/*
  header file for the function that generates the
  weights for a bandwidth bw legendre transform.

  see makeweights.c for details.

*/


#ifndef _MAKEWEIGHTS_H
#define _MAKEWEIGHTS_H

#if defined ( _WIN64 ) || defined ( _WIN32 )
extern void __declspec(dllexport) makeweights( int ,
             double * ) ;
#else
extern void makeweights( int ,
             double * ) ;
#endif

extern void makeweights2( int ,
			  double * ) ;

#if defined ( _WIN64 ) || defined ( _WIN32 )
extern void __declspec(dllexport) releaseSOFTMemory ( double* relMe );
#else
extern void releaseSOFTMemory ( double* relMe );
#endif

#if defined ( _WIN64 ) || defined ( _WIN32 )
extern void __declspec(dllexport) releaseSOFTMemoryMulti ( double** relMe, int size );
#else
extern void releaseSOFTMemoryMulti ( double** relMe, int size );
#endif

#endif  /* _MAKEWEIGHTS_H */
