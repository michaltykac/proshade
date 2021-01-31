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

/* source code for generating cosine transforms 
   of Pml and Gml functions */

#include <math.h>
#include <string.h>   /* to declare memcpy */
#include <stdlib.h>
#include <stdio.h>

#include "fftw3.h"
#include "s2_primitive.h"
/* #include "s2_pmls.h" */


/************************************************************************/
/* utility functions for table management */
/************************************************************************/
/* Computes the number of non-zero entries in a table containing
   cosine series coefficients of all the P(m,l) or G(m,l) functions
   necessary for computing the seminaive transform for a given
   bandwidth bw and order m.  Works specifically for tables
   generated by CosPmlTableGen() 
*/

int TableSize(int m,
	      int bw)
{

  int k ;
  int fudge, fudge2 ;
  int a1, a2, a3 ;

  if ( bw % 2 )  /* if the bandwidth is odd */
    {
      k = bw/2 ;
      fudge = (m+1)%2 ;

      a1 = k*(k+1);
      a2 = fudge*(k+1);

      fudge2 = m/2;
      a3 = fudge2*(fudge2+1);

    }
  else /* bandwidth is even */
    {
      k = bw/2 ;
      fudge = m%2 ;

      a1 = (k-fudge)*(k-fudge+1);
      a2 = fudge*k;

      fudge2 = m/2;
      a3 = fudge2*(fudge2+1);

    }
    
  return(a1+a2-a3);
      
}
/************************************************************************/
/* Spharmonic_TableSize(bw) returns an integer value for
   the amount of space necessary to fill out an entire spharmonic
   table.  Note that in the above TableSize() formula, 
   you need to sum this formula over m as m ranges from 0 to
   (bw-1).  The critical closed form that you need is that

   \sum_{k=0}^n = \frac{(n(n+1)(2n+1)}{6}

   You also need to account for integer division.
   From this you should derive an upper bound on the
   amount of space. 

   Some notes - because of integer division, you need to account for
   a fudge factor - this is the additional " + bw" at the
   end.  This gaurantees that you will always have slightly more
   space than you need, which is clearly better than underestimating!
   Also, if bw > 512, the closed form
   fails because of the bw*bw*bw term (at least on Sun Sparcstations)
   so the loop computation is used instead.

   Also, the transpose is exactly the same size, obviously.

*/

int Spharmonic_TableSize(int bw)
{
  int m, sum;
  
  if (bw > 512)
    {
      sum = 0;
      
      for (m=0; m<bw; m++)
	sum += TableSize(m,bw);
      
      return sum;
    }
  else
    {
      return (
	      (((4*(bw*bw*bw)) + (6*(bw*bw)) - (8*bw))/24)
	      + bw
	      );
    }
}

/************************************************************************/
/* Reduced_Spharmonic_TableSize(bw,m) returns an integer value for
   the amount of space necessary to fill out a spharmonic table
   if interesting in using it only for orders up to (but NOT
   including) order m.
   This will be used in the hybrid algorithm's call of the
   semi-naive algorithm (which won't need the full table ... hopefully
   this'll cut down on the memory usage).

   Also, the transpose is exactly the same size, obviously.

   This is a "reduced" version of Spharmonic_TableSize(m).

*/

#if defined ( _WIN64 ) || defined ( _WIN32 )
int __declspec(dllexport) Reduced_SpharmonicTableSize(int bw,
				int m)
#else
int Reduced_SpharmonicTableSize(int bw,
                int m)
#endif
{
  
  int i, sum;

  sum = 0;
  
  for (i=0; i<m; i++)
    sum += TableSize(i,bw);

  return sum;
}


/************************************************************************/
/* For an array containing cosine series coefficients of Pml or Gml
   functions, computes the location of the first coefficient of Pml.
   This supersedes the TableOffset() function.
   Assumes table is generated by CosPmlTableGen()
*/

int NewTableOffset(int m,
		   int l)
{
  int offset;
  int tm, tl;
    
  if ( m % 2 )
    {
      tl = l-1;
      tm = m-1;
    }
  else
    {
      tl = l;
      tm = m;
    }

  offset = ((tl/2)*((tl/2)+1)) - ((tm/2)*((tm/2)+1));
  if (tl % 2)
    offset += (tl/2)+1;

  return offset;
}

/************************************************************************/
/* generate all of the Pmls for a specified value of m.  

   storeplm points to a double array of size 2 * bw * (bw - m);

   Workspace needs to be
   16 * bw 

   P(m,l,j) respresents the associated Legendre function P_l^m
   evaluated at the j-th Chebyshev point (for the bandwidth bw)
   Cos((2 * j + 1) * PI / (2 * bw)).

   The array is placed in storeplm as follows:

   P(m,m,0)    P(m,m,1)  ... P(m,m,2*bw-1)
   P(m,m+1,0)  P(m,m+1,1)... P(m,m+1,2*bw-1)
   P(m,m+2,0)  P(m,m+2,1)... P(m,m+2,2*bw-1)
   ...
   P(m,bw-1,0)   P(m,bw-1,1) ... P(m,bw-1,2*bw-1)

   This array will eventually be used by the naive transform algorithm.
   This function will precompute the arrays necessary for the algorithm.
*/

void PmlTableGen(int bw,
		 int m,
		 double *storeplm,
		 double *workspace)
{
  double *prev, *prevprev;
  double *temp1, *temp2, *temp3, *temp4;
  double *x_i, *eval_args;
  int i;
  
  prevprev = workspace;
  prev = prevprev + (2*bw);
  temp1 = prev + (2*bw);
  temp2 = temp1 + (2*bw);
  temp3 = temp2 + (2*bw);
  temp4 = temp3 + (2*bw);
  x_i = temp4 + (2*bw);
  eval_args = x_i + (2*bw);
  

  /* get the evaluation nodes */
  EvalPts(2*bw,x_i);
  ArcCosEvalPts(2*bw,eval_args);
  
  /* set initial values of first two Pmls */
  for (i=0; i<2*bw; i++) 
    prevprev[i] = 0.0;
  if (m == 0)
    for (i=0; i<2*bw; i++)
      prev[i] = 0.707106781186547 ;  /* 1/sqrt(2) */
  else 
    Pmm_L2(m, eval_args, 2*bw, prev);

  memcpy(storeplm, prev, sizeof(double) * 2 * bw);

  for(i = 0; i < bw - m - 1; i++)
    {
      vec_mul(L2_cn(m,m+i),prevprev,temp1,2*bw);
      vec_pt_mul(prev, x_i, temp2, 2*bw);
      vec_mul(L2_an(m,m+i), temp2, temp3, 2*bw);
      vec_add(temp3, temp1, temp4, 2*bw); /* temp4 now contains P(m,m+i+1) */
      
      storeplm += (2 * bw);
      memcpy(storeplm, temp4, sizeof(double) * 2 * bw);
      memcpy(prevprev, prev, sizeof(double) * 2 * bw);
      memcpy(prev, temp4, sizeof(double) * 2 * bw);
    }
}






/************************************************************************/
/*
  generate all of the cosine series for L2-normalized Pmls or Gmls for
  a specified value of m. Note especially that since series are
  zero-striped, all zeroes have been removed.  

  tablespace points to a double array of size TableSize(m,bw);
  
  Workspace needs to be
  9 * bw 

  Let P(m,l,j) represent the jth coefficient of the
  cosine series representation of Pml.  The array
  stuffed into tablespace is organized as follows:

  P(m,m,0)    P(m,m,2)  ...  P(m,m,m)
  P(m,m+1,1)  P(m,m+1,3)...  P(m,m+1,m+1)
  P(m,m+2,0)  P(m,m+2,2) ... P(m,m+2,m+2)
  
  etc.  Appropriate modifications are made for m odd (Gml functions).


  NOTE that the Pmls or Gmls are being sampled at bw-many points,
  and not 2*bw-many points. I can get away with this. HOWEVER, I
  need to multiply the coefficients by sqrt(2), because the expected
  input of the seminaive transform of bandwidth bw will be sampled
  at 2-bw many points. So the sqrt(2) is a scaling factor.


*/

void CosPmlTableGen(int bw,
		    int m,
		    double *tablespace,
		    double *workspace)
{
  double *prev, *prevprev, *temp1, *temp2, *temp3, *temp4;
  double *x_i, *eval_args;
  double *tableptr, *cosres ;
  int i, j, k;

  /* fftw stuff now */
  double fudge ;
  fftw_plan p ;

  prevprev = workspace;
  prev = prevprev + bw;
  temp1 = prev + bw;
  temp2 = temp1 + bw;
  temp3 = temp2 + bw;
  temp4 = temp3 + bw;
  x_i = temp4 + bw;
  eval_args = x_i + bw;
  cosres = eval_args + bw;

  tableptr = tablespace;

  /* make fftw plan */
  p = fftw_plan_r2r_1d( bw, temp4, cosres,
			FFTW_REDFT10, FFTW_ESTIMATE ) ;

  /* main loop */

  /* Set the initial number of evaluation points to appropriate
     amount */

  /* now get the evaluation nodes */
  EvalPts(bw,x_i);
  ArcCosEvalPts(bw,eval_args);
    
  /* set initial values of first two Pmls */
  for (i=0; i<bw; i++) 
    prevprev[i] = 0.0;

  if (m == 0)
    for (i=0; i<bw; i++)
      prev[i] = 0.707106781186547; /* sqrt(1/2) */
  else 
    Pmm_L2(m, eval_args, bw, prev);


  if ( m % 2 ) /* need to divide out sin x */
    for (i=0; i<bw; i++)
      prev[i] /= sin(eval_args[i]);
  

  /* set k to highest degree coefficient */
  if ((m % 2) == 0)
    k = m;
  else
    k = m-1;	

  /* now compute cosine transform */
  memcpy( temp4, prev, sizeof(double) * bw );
  fftw_execute( p );
  cosres[0] *= 0.707106781186547 ;
  fudge = 1. / sqrt(((double) bw ) );
  for ( i = 0 ; i < bw ; i ++ )
    cosres[i] *= fudge ;

  /* store what I've got so far */
  for (i=0; i<=k; i+=2)
    tableptr[i/2] = cosres[i];

  /* update tableptr */
  tableptr += k/2+1;

  /* now generate remaining pmls  */
  for (i=0; i<bw-m-1; i++)
    {
      vec_mul(L2_cn(m,m+i),prevprev,temp1,bw);
      vec_pt_mul(prev, x_i, temp2, bw);
      vec_mul(L2_an(m,m+i), temp2, temp3, bw);
      vec_add(temp3, temp1, temp4, bw); /* temp4 now contains P(m,m+i+1) */

      /* compute cosine transform */
      fftw_execute( p );
      cosres[0] *= 0.707106781186547 ;
      for ( j = 0 ; j < bw ; j ++ )
	cosres[j] *= fudge ;

      /* update degree counter */
      k++;

      /* now put decimated result into table */
      if ( i % 2 )
	for (j=0; j<=k; j+=2)
	  tableptr[j/2] = cosres[j];
      else
	for (j=1; j<=k; j+=2)
	  tableptr[j/2] = cosres[j];
      
      /* update tableptr */
      tableptr += k/2+1;

      /* now update Pi and P(i+1) */
      memcpy(prevprev, prev, sizeof(double) * bw);
      memcpy(prev, temp4, sizeof(double) * bw);
    }

  fftw_destroy_plan( p );

}


/************************************************************************/
/* RowSize returns the number of non-zero coefficients in a row of the
   cospmltable if were really in matrix form.  Helpful in transpose
   computations.  It is helpful to think of the parameter l as
   the row of the corresponding matrix.
*/

int RowSize(int m,
	    int l)
{
  if (l < m)
    return 0;
  else
    {
      if ((m % 2) == 0)
	return ((l/2)+1);
      else
	return (((l-1)/2)+1);
    }
}
/************************************************************************/
/* Transposed row size returns the number of non-zero coefficients
   in the transposition of the matrix representing a cospmltable.
   Used for generating arrays for inverse seminaive transform.
   Unlike RowSize, need to know the bandwidth bw.  Also, in
   the cospml array, the first m+1 rows are empty, but in
   the transpose, all rows have non-zero entries, and the first
   m+1 columns are empty.  So the input parameters are a bit different
   in the you need to specify the row you want.

*/

int Transpose_RowSize(int row,
		      int m,
		      int bw)
{
  /* my version might be longer, but at least I understand
     it better, and it's only minimally recursive */

  if ( bw % 2 )
    {
      if ( m % 2 )
	{
	  if ( m == 1 )
	    return( (bw-row)/2 );
	  else if ( row < m - 1 )
	    return ( (bw-m+1)/2 );
	  else
	    return ( Transpose_RowSize(row, 1, bw) ) ;
	}
      else
	{
	  if ( m == 0 )
	    return( (bw-row)/2 + ((row+1)%2) );
	  else if ( row < m )
	    return ( (bw-m)/2 + ((row+1)%2) );
	  else
	    return ( Transpose_RowSize(row, 0, bw) ) ;
	}
    }
  else
    {
      if ( m % 2 )
	{
	  if ( m == 1 )
	    return( (bw-row)/2 );
	  else if ( row < m - 1 )
	    return ( (bw-m+1)/2 - (row%2) );
	  else
	    return ( Transpose_RowSize(row, 1, bw) ) ;
	}
      else
	{
	  if ( m == 0 )
	    return( (bw-row)/2 + (row%2) );
	  else if ( row < m )
	    return ( (bw-m)/2 );
	  else
	    return ( Transpose_RowSize(row, 0, bw) ) ;
	}
    }
  

 
  /*** original version

  if (row >= bw)
  return 0;
  else if ((m % 2) == 0)
  {
  if (row <= m)
  return ( ((bw-m)/2) );
  else
  return ( ((bw-row-1)/2) + 1);
  }
  else
  {
  if (row == (bw-1))
  return 0;
  else if (row >= m)
  return (Transpose_RowSize(row+1,m-1,bw));
  else
  return (Transpose_RowSize(row+1,m-1,bw) - (row % 2));
  }

  ***/

}

/************************************************************************/
/* Inverse transform is transposition of forward transform.
   Thus, need to provide transposed version of table
   returned by CosPmlTableGen.  This function does that
   by taking as input a cos_pml_table for a particular value
   of bw and m, and loads the result as a
   transposed, decimated version of it for use by an inverse 
   seminaive transform computation.

   result needs to be of size TableSize(m,bw)

*/

void Transpose_CosPmlTableGen(int bw,
			      int m,
			      double *cos_pml_table,
			      double *result)
{
  /* recall that cospml_table has had all the zeroes
     stripped out, and that if m is odd, then it is
     really a Gml function, which affects indexing a bit.
  */
  
  double *trans_tableptr, *tableptr;
  int i, row, rowsize, stride, offset, costable_offset;

  /* note that the number of non-zero entries is the same
     as in the non-transposed case */

  trans_tableptr = result;
  
  /* now traverse the cos_pml_table , loading appropriate values
     into the rows of transposed array */

  if ( m == bw - 1 )
    memcpy( result, cos_pml_table, sizeof(double)*TableSize(m,bw));
  else
    {

      for (row = 0; row < bw; row++)
	{
	  /* if m odd, no need to do last row - all zeroes */
	  if (row == (bw-1))
	    {
	      if ( m % 2 )
		return;
	    }

	  /* get the rowsize for the transposed array */
	  rowsize = Transpose_RowSize(row, m, bw);

	  /* compute the starting point for values in cos_pml_table */
	  if (row <= m)
	    {
	      if ((row % 2) == 0)
		tableptr = cos_pml_table + (row/2);
	      else
		tableptr = cos_pml_table + (m/2) + 1 + (row/2);
	    }
	  else
	    {
	      /* if row > m, then the highest degree coefficient
		 of P(m,row) should be the first coefficient loaded
		 into the transposed array, so figure out where
		 this point is.
	      */
	      offset = 0;
	      if ( (m%2) == 0 )
		{
		  for (i=m; i<=row; i++)
		    offset += RowSize(m, i);
		}
	      else
		{
		  for (i=m;i<=row+1;i++)
		    offset += RowSize(m, i);
		}
	      /* now we are pointing one element too far, so decrement */
	      offset--;

	      tableptr = cos_pml_table + offset;
	    }

	  /* stride is how far we need to jump between
	     values in cos_pml_table, i.e., to traverse the columns of the
	     cos_pml_table.  Need to set initial value.  Stride always
	     increases by 2 after that 
	  */
	  if (row <= m)
	    stride = m + 2 - (m % 2) + (row % 2);
	  else
	    stride = row + 2;

	  /* now load up this row of the transposed table */
	  costable_offset = 0;
	  for (i=0; i < rowsize; i++)
	    {
	      trans_tableptr[i] = tableptr[costable_offset];
	      costable_offset += stride;
	      stride += 2;

	    } /* closes i loop */

	  trans_tableptr += rowsize;

	} /* closes row loop */
    }

}
/************************************************************************/
/* This is a function that returns all of the (cosine transforms of)
   Pmls and Gmls necessary
   to do a full spherical harmonic transform, i.e., it calls
   CosPmlTableGen for each value of m less than bw, returning a
   table of tables ( a pointer of type (double **), which points
   to an array of size m, each containing a (double *) pointer
   to a set of cospml or cosgml values, which are the (decimated)
   cosine series representations of Pml (even m) or Gml (odd m)
   functions.  See CosPmlTableGen for further clarification.

   Inputs - the bandwidth bw of the problem
   resultspace - need to allocate Spharmonic_TableSize(bw) for storing results
   workspace - needs to be (16 * bw)

   Note that resultspace is necessary and contains the results/values
   so one should be careful about when it is OK to re-use this space.
   workspace, though, does not have any meaning after this function is
   finished executing

*/

double **Spharmonic_Pml_Table(int bw,
			      double *resultspace,
			      double *workspace)
{

  int i;
  double **spharmonic_pml_table;

  /* allocate an array of double pointers */
  spharmonic_pml_table = (double **) malloc(sizeof(double *) * bw);

  /* traverse the array, assigning a location in the resultspace
     to each pointer */

  spharmonic_pml_table[0] = resultspace;

  for (i=1; i<bw; i++)
    {
      spharmonic_pml_table[i] = spharmonic_pml_table[i-1] +
	TableSize(i-1,bw);
    }
  
  /* now load up the array with CosPml and CosGml values */
  for (i=0; i<bw; i++)
    {
      CosPmlTableGen(bw, i, spharmonic_pml_table[i], workspace);
    }

  /* that's it */

  return spharmonic_pml_table;
}


/************************************************************************/
/* For the inverse semi-naive spharmonic transform, need the "transpose"
   of the spharmonic_pml_table.  Need to be careful because the
   entries in the spharmonic_pml_table have been decimated, i.e.,
   the zeroes have been stripped out.

   Inputs are a spharmonic_pml_table generated by Spharmonic_Pml_Table
   and the bandwidth bw

   Allocates memory for the (double **) result
   also allocates memory

   resultspace - need to allocate Spharmonic_TableSize(bw) for storing results
   workspace - not needed, but argument added to avoid
               confusion wth Spharmonic_Pml_Table

*/

double **Transpose_Spharmonic_Pml_Table(double **spharmonic_pml_table, 
					int bw,
					double *resultspace,
					double *workspace)
{
  
  int i;
  double **transpose_spharmonic_pml_table;

  /* allocate an array of double pointers */
  transpose_spharmonic_pml_table = (double **) malloc(sizeof(double *) * bw);

  /* now need to load up the transpose_spharmonic_pml_table by transposing
     the tables in the spharmonic_pml_table */

  transpose_spharmonic_pml_table[0] = resultspace;

  for (i=0; i<bw; i++)
    {
      Transpose_CosPmlTableGen(bw, 
			       i, 
			       spharmonic_pml_table[i],
			       transpose_spharmonic_pml_table[i]);

      if (i != (bw-1))
	{
	  transpose_spharmonic_pml_table[i+1] =
	    transpose_spharmonic_pml_table[i] + TableSize(i, bw);
	}
    }

  return transpose_spharmonic_pml_table;
}


/************************************************************************/
/* Reduced_Naive_TableSize(bw,m) returns an integer value for
   the amount of space necessary to fill out a reduced naive table
   of pmls if interested in using it only for orders m through bw - 1.

*/

#if defined ( _WIN64 ) || defined ( _WIN32 )
int __declspec(dllexport) Reduced_Naive_TableSize(int bw,
			    int m)
#else
int Reduced_Naive_TableSize(int bw,
                int m)
#endif
{
  
  int i, sum;

  sum = 0;
  
  for (i=m; i<bw; i++)
    sum += ( 2 * bw * (bw - i));

  return sum;

}

/*************************************************************

  just like Spharmonic_Pml_Table(), except generates a
  table for use with the semi-naive and naive algorithms.

  m is the cutoff order, where to switch from semi-naive to
  naive algorithms

  bw = bandwidth of problem
  m  = where to switch algorithms (order where naive is FIRST done)
  resultspace = where to store results, must be of
                size Reduced_Naive_TableSize(bw, m) +
		     Reduced_SpharmonicTableSize(bw, m);

***********************************************************/

#if defined ( _WIN64 ) || defined ( _WIN32 )
double __declspec(dllexport) **SemiNaive_Naive_Pml_Table(int bw,
				   int m,
				   double *resultspace,
				   double *workspace)
#else
double **SemiNaive_Naive_Pml_Table(int bw,
                   int m,
                   double *resultspace,
                   double *workspace)
#endif
{
  int i;
  double **seminaive_naive_table;
  int  lastspace;

  seminaive_naive_table = (double **) malloc(sizeof(double) * (bw+1));

  seminaive_naive_table[0] = resultspace;

  
  for (i=1; i<m; i++)
    {
      seminaive_naive_table[i] = seminaive_naive_table[i - 1] +
	TableSize(i-1,bw);
    }

  if( m == 0)
    {
      lastspace = 0;
      for (i=m+1; i<bw; i++)
	{ 
	  seminaive_naive_table[i] = seminaive_naive_table[i - 1] +
	    (2 * bw * (bw - (i - 1)));
	}
    }
  else
    {
      lastspace = TableSize(m-1,bw);
      seminaive_naive_table[m] = seminaive_naive_table[m-1] +
	lastspace;
      for (i=m+1; i<bw; i++)
	{ 
	  seminaive_naive_table[i] = seminaive_naive_table[i - 1] +
	    (2 * bw * (bw - (i - 1)));
	}
    }

  /* now load up the array with CosPml and CosGml values */
  for (i=0; i<m; i++)
    {
      CosPmlTableGen(bw, i, seminaive_naive_table[i], workspace);
    }

  /* now load up pml values */
  for(i=m; i<bw; i++)
    {
      PmlTableGen(bw, i, seminaive_naive_table[i], workspace);
    }

  /* that's it */

  return seminaive_naive_table;

}


/************************************************************************/
/* For the inverse seminaive_naive transform, need the "transpose"
   of the seminaive_naive_pml_table.  Need to be careful because the
   entries in the seminaive portion have been decimated, i.e.,
   the zeroes have been stripped out.

   Inputs are a seminaive_naive_pml_table generated by SemiNaive_Naive_Pml_Table
   and the bandwidth bw and cutoff order m

   Allocates memory for the (double **) result
   also allocates memory

   resultspace - need to allocate Reduced_Naive_TableSize(bw, m) +
		     Reduced_SpharmonicTableSize(bw, m) for storing results
   workspace - size 16 * bw 

*/

double **Transpose_SemiNaive_Naive_Pml_Table(double **seminaive_naive_pml_table, 
					     int bw,
					     int m,
					     double *resultspace,
					     double *workspace)
{
  
  int i, lastspace;
  double **trans_seminaive_naive_pml_table;

  /* allocate an array of double pointers */
  trans_seminaive_naive_pml_table = (double **) malloc(sizeof(double *) * (bw+1));

  /* now need to load up the transpose_seminaive_naive_pml_table by transposing
     the tables in the seminiave portion of seminaive_naive_pml_table */

  trans_seminaive_naive_pml_table[0] = resultspace;

  
  for (i=1; i<m; i++)
    {
      trans_seminaive_naive_pml_table[i] =
	trans_seminaive_naive_pml_table[i - 1] +
	TableSize(i-1,bw);
    }

  if( m == 0 )
    {
      lastspace = 0;
      for (i=m+1; i<bw; i++)
	{ 
	  trans_seminaive_naive_pml_table[i] =
	    trans_seminaive_naive_pml_table[i - 1] +
	    (2 * bw * (bw - (i - 1)));
	}
    }
  else
    {
      lastspace = TableSize(m-1,bw);
      trans_seminaive_naive_pml_table[m] =
	trans_seminaive_naive_pml_table[m-1] +
	lastspace;

      for (i=m+1; i<bw; i++)
	{ 
	  trans_seminaive_naive_pml_table[i] =
	    trans_seminaive_naive_pml_table[i - 1] +
	    (2 * bw * (bw - (i - 1)));
	}
    }

  for (i=0; i<m; i++)
    {
      Transpose_CosPmlTableGen(bw, 
			       i, 
			       seminaive_naive_pml_table[i],
			       trans_seminaive_naive_pml_table[i]);

      if (i != (bw-1))
	{
	  trans_seminaive_naive_pml_table[i+1] =
	    trans_seminaive_naive_pml_table[i] + TableSize(i, bw);
	}
    }

  /* now load up pml values */
  for(i=m; i<bw; i++)
    {
      PmlTableGen(bw, i, trans_seminaive_naive_pml_table[i], workspace);
    }

  return trans_seminaive_naive_pml_table;
}
