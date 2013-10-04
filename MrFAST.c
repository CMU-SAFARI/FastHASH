/*
 * Copyright (c) <2008 - 2020>, University of Washington, Simon Fraser University, 
 * Bilkent University and Carnegie Mellon University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this list
 * of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or other
 *   materials provided with the distribution.
 * - Neither the names of the University of Washington, Simon Fraser University, 
 *   Bilkent University, Carnegie Mellon University,
 *   nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without specific
 *   prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

  Authors: 
  Farhad Hormozdiari
	  farhadh AT uw DOT edu
  Faraz Hach
	  fhach AT cs DOT sfu DOT ca
  Can Alkan
	  calkan AT gmail DOT com
  Hongyi Xin
	  gohongyi AT gmail DOT com
  Donghyuk Lee
	  bleups AT gmail DOT com
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dirent.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <mmintrin.h>

#include "Common.h"
#include "Reads.h"
#include "HashTable.h"
#include "Output.h"
#include "MrFAST.h"
#include "RefGenome.h"


#define min(a,b) ((a)>(b)?(b):(a))
#define min3(a,b,c) ((a)>(b)?(b>c?c:b):(a>c?c:a))
#define CHARCODE(a) (a=='A' ? 0 : (a=='C' ? 1 : (a=='G' ? 2 : (a=='T' ? 3 : 4))))

double binomial_coefficient(int n, int k);
float calculateScore(int index, char *seq, char *qual, char *md);


/* ************************************************************************** */
/* ************************************************************************** */
/* 								GLOBAL VARIABLES                              */
/* ************************************************************************** */
/* ************************************************************************** */

/* _msf_ stands for mr fast */
unsigned char	 mrFAST			= 1;
char			*versionNumberF = "0.0";


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* This is only used in verifySingleEnd function, which is never called
 * anywhere. */
long long	 verificationCnt = 0;

/* Number of successful mappings of k-mers to the reference genome. */
long long	 mappingCnt		 = 0;

/* Number of mapped read sequences */
long long	 mappedSeqCnt	 = 0;

/* Only utilized in mapSingleEndSeq function. When we map a read for maxHits,
 * this is incremented. Actual purpose of this could not be discovered yet. */
long long	 completedSeqCnt = 0;

/* This seems not to be used in normal mode mapping. It is used in paired end
 * mapping. */
char		*mappingOutput;


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* reference genome related information. */
char		*_msf_refGen	   = NULL;
int			 _msf_refGenLength = 0;
int			 _msf_refGenOffset = 0;
char		*_msf_refGenName   = NULL;
int			 _msf_refGenBeg;
int			 _msf_refGenEnd;
IHashTable	*_msf_hashTable	   = NULL;


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Sampling locations. */
int *_msf_samplingLocs;			/* Starting position of sampling locations. */
int *_msf_samplingLocsEnds;		/* Ending position of sampling locations. */
int	 _msf_samplingLocsSize;


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Read sequence list info. */
Read	*_msf_seqList;
int		 _msf_seqListSize;


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Sorted reads wrt hash values of the complete read sequence. */
Pair	*_msf_sort_seqList = NULL;
int		*_msf_map_sort_seqList;


/* -------------------------------------------------------------------------- */
ReadIndexTable	*_msf_rIndex = NULL;
int				 _msf_rIndexSize;
int				 _msf_rIndexMax;


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Output related variables. These are used almost in every output function. */
SAM			 _msf_output;
OPT_FIELDS	*_msf_optionalFields; 
char		*_msf_op;


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Verified locations is for storing if a read's segment has been processed in
 * the reference genome location. It is of size |reference genome
 * length|. Entries are the readIds (related with read numbers - two entries for
 * each read, FORWARD and BACKWARD). Each read's k-mers operate on this
 * array. */
int			*_msf_verifiedLocs = NULL;


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Used in paired-end mode for storing mapping information for each read. Linked
 * list, as much as number of mapping locations. */
MappingInfo *_msf_mappingInfo;

char	_msf_numbers[200][3];
char	_msf_cigar[5];


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Used in paired-end mode. Of size #reads. Used to store the number of the
 * sequence mappings.*/
int *_msf_seqHits;

/* #open files, max no of left and right mapped sequences for reads. They are
 * summed up each time mapPairedEndSeq() function called. They are used in
 * outputPairedEnd function. The _msf_openFiles is incremented each time
 * mapPairedEndSeq() is called. */
int	 _msf_openFiles = 0;
int	 _msf_maxLSize	= 0;
int	 _msf_maxRSize	= 0;


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Best mapping information, used in best mapping mode. */
BestFullMappingInfo *bestHitMappingInfo;


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Utilized in outputPairedEnd. */
/* _msf_maxFile is incremented each time outputPairedEnd is called.  */
int		_msf_maxFile = 0;
char	_msf_fileName[4000][200][2][FILE_NAME_LENGTH];
int		_msf_fileCount[4000];


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Boolean array of size #reads/2 storing information about if a read has
 * concordant mapping or not. */
char	*_msf_readHasConcordantMapping;


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* These are used in paired-end mapping mode. Their size are equal to number of
 * reads. So an integer is stored for each read. */

/* Number of OEA mappings. A read can have two pairs, accessed thru 2*i and
 * 2*i+1. */
int *_msf_oeaMapping;

/* Number of discordant mappings. A read can have two pairs, accessed thru 2*i
 * and 2*i+1. */
int *_msf_discordantMapping;


/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Throughout the code, F stands for "Forward" and B stands for "Backward" */
/* These are heavily used in edit distance computation */
int scoreF[SEQ_MAX_LENGTH][SEQ_MAX_LENGTH];
int scoreB[SEQ_MAX_LENGTH][SEQ_MAX_LENGTH];
int score[SEQ_MAX_LENGTH][SEQ_MAX_LENGTH];
int direction1[SEQ_MAX_LENGTH][SEQ_MAX_LENGTH]; /* For left segment */
int direction2[SEQ_MAX_LENGTH][SEQ_MAX_LENGTH]; /* For right segment */


/* -------------------------------------------------------------------------- */
/* Need to clarify all the locations this is used. */
__m128i MASK;

/* ************************************************************************** */
/* ************************************************************************** */
/* 								GLOBAL VARIABLES                              */
/* ************************************************************************** */
/* ************************************************************************** */



/*------------------------------------------------------------------------------
 * Runs the edit distance for small mappings (forward). Called when length <= e.
 *----------------------------------------------------------------------------*/
int
smallEditDistanceF ( char *a, int lena, char *b, int lenb )
{
	/* Runs the full edit distance without SSE instructions. */
	int matrix[20][20];
	int i = 0;
	int j = 0;

	for(i = 0; i <= lena; i++)
    {
		matrix[0][i] = i;
    }

	for(i = 0; i <= lenb; i++)
	{
		matrix[i][0] = i;
	}

	for(i = 1; i <= lenb; i++)
	{
		for(j = 1; j <= lena; j++)
		{
			matrix[i][j] = min3(matrix[i-1][j-1]+ (a[j-1] != b[i-1]),matrix[i][j-1]+1 ,matrix[i-1][j]+1);
	  	}
	}
	return (matrix[lenb][lena] > errThreshold ? -1 : matrix[lenb][lena]);
}



/*------------------------------------------------------------------------------
 * Runs the edit distance for small mappings (backward). Called when length <= e
 *----------------------------------------------------------------------------*/
int
smallEditDistanceB ( char *a, int lena, char *b, int lenb )
{
	int matrix[20][20];
	int i = 0;
	int j = 0;

	for(i = 0; i <= lena; i++)
	{
		matrix[0][i] = i;
	}

	for(i = 0; i <= lenb; i++)
	{
		matrix[i][0] = i;
	}

	for(i = 1; i <= lenb; i++)
	{
		for(j = 1; j <= lena; j++)
	  	{
			matrix[i][j] = min3(matrix[i-1][j-1]+ (*(a-j+1) != *(b-i+1)),matrix[i][j-1]+1 ,matrix[i-1][j]+1);
	  	}
	}

	return (matrix[lenb][lena] > errThreshold ? -1 : matrix[lenb][lena]);
}



/*------------------------------------------------------------------------------
 * This is NOT called. XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 *----------------------------------------------------------------------------*/
char
fastEditDistance ( int per1, int per2 )
{
	int		i	  = 0;
	int		j	  = 0;
	char	str1[7];
	char	str2[7];
	int		val1  = per1;
	int		val2  = per2;
	int		index = 0;
	int		mod	  = 0;
	int		matrix[7][7];
	int		min	  = 20;

	while(index < 6)
	{
		mod			  = val1 % 5;
		str1[5-index] = (mod==0 ? 'A':(mod==1?'C':mod==2?'G':(mod==3)?'T':'N'));
		val1		  = val1 / 5;
		index++;
	}

	str1[6] = '\0';

	index = 0;
	while(index < 6)
	{
		mod			  = val2 % 5;
		str2[5-index] = (mod==0 ? 'A':(mod==1?'C':mod==2?'G':(mod==3)?'T':'N'));
		val2		  = val2 / 5;
		index++;
	}
	
	str2[6] = '\0';

	for(i = 0; i < 7; i++)
	{
		matrix[0][i] = i;
		matrix[i][0] = i;
	}

	for(i = 1; i < 7; i++)
	{
		for(j = 1; j < 7; j++)
		{
			matrix[i][j] = min3(matrix[i-1][j-1]+ (str1[i-1] != str2[j-1]),matrix[i][j-1]+1 ,matrix[i-1][j]+1);
	  	}
	}

	for(i = 0; i < 7; i++)
	{
		if (matrix[i][6] < min)
	  		min = matrix[i][6];
	}

	for(i = 0; i < 7; i++)
	{
		if (matrix[6][i] < min)
	  		min = matrix[6][i];
	}
	
	return min;
}



/*------------------------------------------------------------------------------
 * Initializes the global SSE2 variable MASK and scoreF and scoreB matrices.
 *----------------------------------------------------------------------------*/
void
initLookUpTable (   )
{
	int i = 0;

	MASK = _mm_insert_epi16(MASK,1,0);
	MASK = _mm_insert_epi16(MASK,1,1);
	MASK = _mm_insert_epi16(MASK,1,2);
	MASK = _mm_insert_epi16(MASK,1,3);
	MASK = _mm_insert_epi16(MASK,1,4);
	MASK = _mm_insert_epi16(MASK,0,5);
	MASK = _mm_insert_epi16(MASK,0,6);
	MASK = _mm_insert_epi16(MASK,0,7);

	for(i = 0; i < errThreshold + 1; i++)
	{
		scoreF[0][i] = i;
		scoreF[i][0] = i;
	}

	for(i = 0 ; i < errThreshold + 1; i++)
	{
		scoreB[0][i] = i;
		scoreB[i][0] = i;
	}
}



/*------------------------------------------------------------------------------
 * -FULLY COMMENTED-
 * 	a		: ref genome
 *	lena	: length of current portion of the ref genome that matching will be
 *			  performed
 *	b		: read sequence
 * 	lenb	: length the of the read sequence

 *	Called when e = 2 (verifySingleEndEditDistance2)
 *----------------------------------------------------------------------------*/
int
forwardEditDistance2SSE2 ( char *a, int lena, char *b, int lenb )
{
	if (lenb == 0 || lena == 0)
	  	return 0;

	int i0 = 0;
	int i1 = 0;

	int error;					//0: if the two character are equal 1: if not, good!
	int i		   = 0;			//loop index
	int e		   = 2;			//error bound
	int totalError = 0;

	/* 128-bit variables */
	/* R0 stores the result of even-numbered antidiagonals, R1 stores the result
	 * of odd-numbered antidiagonals. */
	__m128i R0;
	__m128i R1;

	/* These are helping variables, while computing main diagonal of the edit
	 * distance matrix. Diag stores the newly computed edit distance result. */
	__m128i Side1, Side;		//side matrix
	__m128i Down1, Down2, Down; //down matrix
	__m128i Diag;

	__m128i tmp;
	__m128i ERROR_REACH;
	
	/* R0 and R1 store two latest computed consecutive antidiagonals */
	R0			= _mm_setzero_si128();
	R1			= _mm_setzero_si128();

	/* These 6 variables are for the edge of the banded diagonal */
	Diag		= _mm_setzero_si128();
	Side1		= _mm_setzero_si128();
	Down1		= _mm_setzero_si128();
	Down2		= _mm_setzero_si128();
	Side		= _mm_setzero_si128();
	Down		= _mm_setzero_si128();

	/* Used for if the error threshold has been exceeded. */
	tmp			= _mm_setzero_si128();
	ERROR_REACH = _mm_setzero_si128();
	
	/* e = 2, if length of the read sequence is <= 2, just run small edit
	 * distance without bothering. */
	if (lenb <= e)
	{
	  	return smallEditDistanceF(a, lena, b, lenb);
	}

	/* _mm_set_epi16 sets the 8 16-bit signed integer values. */
	ERROR_REACH = _mm_set_epi16(0, 0, 0, 0, 0, e, e, e);

	/* Initialize R0 and R1 using a_0 and b_0  */
	R0 = _mm_insert_epi16(R0,0,0);
	R1 = _mm_insert_epi16(R1,1,0);
	R1 = _mm_insert_epi16(R1,1,1);

	Diag  = _mm_set_epi16(0, 0, 0, 0, 0, 2 * e, ((a[0]) != (b[0])), 2 *e);
	Side1 = _mm_set_epi16(0, 0, 0, 0, 0, 2 * e, 1, 1);
	Down1 = _mm_set_epi16(0, 0, 0, 0, 0, 2 * e, 1, 1);
	Down2 = _mm_set_epi16(0, 0, 0, 0, 0, 1, 1, 2 * e);

	tmp = _mm_slli_si128(R1,2);
	R0	= _mm_min_epi16(R1 + Side1, R0 + Diag);
	R0	= _mm_min_epi16(R0, tmp + Down2);

	for (i = 3; i < 2 * lena; i++)
	{		
		if (i % 2 == 1)			/* Odd numbered indexes - R1 */
		{
			/* Diag stores the currently computed two values */
			Diag  = _mm_xor_si128(Diag, Diag);
			
			error = ((a[(i + 1) / 2 - 1]) != (b[(i - 1) / 2 - 1]));
			Diag  = _mm_insert_epi16(Diag,error,0);
						
			error = ((a[(i - 1) / 2 - 1]) != (b[(i + 1) / 2 - 1]));
			Diag  = _mm_insert_epi16(Diag,error,1);
						
			/* Diag = _mm_set_epi16(0, 0, 0, 0, 0, 0, ((a[(i-1)/2-1]) !=
			 * (b[(i+1)/2-1])) ,((a[(i+1)/2-1]) != (b[(i-1)/2-1]))); */
			
			tmp = _mm_srli_si128(R0,2);

			/* Update two entries of R1 (stored in the 32 LSB bits of R1) */
			R1 = _mm_min_epi16(tmp + Side1, R1 + Diag);
			R1 = _mm_min_epi16(R1, R0 + Down1);

			if (i > 2 * lenb - 2)
			{
				i1 = _mm_extract_epi16(R1, 1);
		  		totalError = min(totalError, i1);
			}
		}
		else if (i % 2 == 0)	/* Even numbered indexes - R0 */
		{
			/* Diag stores the currently computed three values */
			error = ((a[i / 2]) != (b[i / 2 - 2]));
			Diag  = _mm_insert_epi16(Diag,error,0);
			
			error = ((a[i / 2 - 1]) != (b[i / 2 - 1]));
			Diag  = _mm_insert_epi16(Diag,error,1);
			
			error = ((a[i / 2 - 2]) != (b[i / 2]));
			Diag  = _mm_insert_epi16(Diag,error,2);			
						
			/*  Diag = _mm_set_epi16(0, 0, 0, 0, 0, ((a[i/2-2]) != (b[i/2])) ,
			 *  ((a[i/2-1]) != (b[i/2-1])) , ((a[i/2]) != (b[i/2-2])) ); */

			/* Update three entries of R0 (stored in the 48 LSB
			 * bits of R0) */
			tmp = _mm_slli_si128(R1,2);

			R0 = _mm_min_epi16(R1 + Side1, R0 + Diag);
			R0 = _mm_min_epi16(R0, tmp + Down2);

			tmp = _mm_sub_epi16(ERROR_REACH, R0);
			i0	= _mm_movemask_epi8(tmp); /* This is signed (63 = 111111) */

			/* Quit if we exceed the error threshold, it checks R1, and R0. */
			if (i0 == 63 &&
				_mm_extract_epi16(R1,0) > errThreshold &&
				_mm_extract_epi16(R1,1) > errThreshold &&
				i < 2 * lenb - 2) /* Did not reach the end yet */
			{
		  		return -1;
			}
			
			if (i == 2 * lenb - 2)
			{				
		  		totalError = _mm_extract_epi16(R0, 2);
			}
		}
	}
	
	Down1 = _mm_insert_epi16(Down1,2*e,0);


	/* ********************************************************************* */
	/* First part of the error */
	error = ((a[i / 2]) != (b[i / 2 - 2]));
	Diag  = _mm_insert_epi16(Diag,error,0);
	
	error = ((a[i / 2 - 1]) != (b[i / 2 - 1]));
	Diag  = _mm_insert_epi16(Diag,error,1);
	
	Diag  = _mm_insert_epi16(Diag,2*e,2);
	/*        Diag = _mm_set_epi16(0, 0, 0, 0, 0, 2*e , ((a[i/2-1]) !=
	 *        (b[i/2-1])) , ((a[i/2]) != (b[i/2-2])) ); */

	R0		   = _mm_min_epi16(R1 + Side1, R0 + Diag);
	R0		   = _mm_min_epi16(R0, _mm_slli_si128(R1,2) + Down1);
	i1		   = _mm_extract_epi16(R0, 1);
	totalError = min(totalError, i1);

	
	/* ********************************************************************* */
	/* Second part of the error */
	i++;

	Diag  = _mm_xor_si128(Diag, Diag);
	Diag  = _mm_insert_epi16(Diag,2*e,0);
	error = ((a[i / 2]) != (b[lenb - 1]));
	
	Diag  = _mm_insert_epi16(Diag,error,1);
	Diag  = _mm_insert_epi16(Diag,2*e,2);
	/*        Diag = _mm_set_epi16(0, 0, 0, 0, 0, 2*e , ((a[i/2]) !=
	 *        (b[lenb-1])) , 2*e ); */

	R1		   = _mm_min_epi16(R0 + Side1, _mm_slli_si128(R1,2) + Diag);
	R1		   = _mm_min_epi16(R1, _mm_slli_si128(R0,2) + Down1);
	i1		   = _mm_extract_epi16(R1, 1);
	totalError = min(totalError, i1);

	
	/* ********************************************************************* */
	/* Final element of the matrix */
	i++;
	
	Diag  = _mm_xor_si128(Diag, Diag);
	error = ((a[i / 2]) != (b[lenb - 1]));
	
	Diag  = _mm_insert_epi16(Diag,error,0);

	/*        Diag = _mm_set_epi16(0, 0, 0, 0, 0, 0 , 0 , ((a[i/2]) !=
	 *        (b[lenb-1])) ); */

	Down = _mm_insert_epi16(Down,1,0);
	Side = _mm_insert_epi16(Side,1,0);
	tmp	 = _mm_srli_si128(R1,2);

	R0		   = _mm_min_epi16(R1 + Down, _mm_srli_si128(R0,2) + Diag);
	R0		   = _mm_min_epi16(R0, tmp + Side);
	i0		   = _mm_extract_epi16(R0, 0);
	totalError = min(totalError, i0);

	/* If final error is greater than error threshold, return
	 * -1  */
	if (totalError > e)
	  	return -1;

	return totalError;
}



/*------------------------------------------------------------------------------
 * -FULLY COMMENTED-
 * 	a		: ref genome
 *	lena	: length of current portion of the ref genome that matching will be
 *			  performed
 *	b		: read sequence
 * 	lenb	: length the of the read sequence

 * Called when e = 2 (verifySingleEndEditDistance2)

 * This is actually exactly same with the forwardEditDistance2SSE2. The only
 * difference is that in each iteration we go backwards in a & b instead of
 * going forward. So it is like a - ... and b - ... instead of a[] & b[] which
 * translate to a + ... & b + ...
 *----------------------------------------------------------------------------*/
int
backwardEditDistance2SSE2 ( char *a, int lena, char *b, int lenb )
{
	if (lenb == 0 || lena == 0)
	  	return 0;

	int i0 = 0;
	int i1 = 0;

	int error; /* 0: if the two character are equal 1: if not, thanks for this
				* comment, it made my day. */
	int i		   = 0;			//loop index
	int e		   = 2;			//error bound
	int totalError = 0;


	/* 128-bit variables */
	/* R0 stores the result of even-numbered antidiagonals, R1 stores the result
	 * of odd-numbered antidiagonals. */
	__m128i R0;
	__m128i R1;

	/* These are helping variables, while computing main diagonal of the edit
	 * distance matrix. Diag stores the newly computed edit distance result. */	
	__m128i Side1, Side2, Side; //side matrix
	__m128i Down1, Down2, Down; //down matrix
	__m128i Diag; //diag matrix

	__m128i tmp;
	__m128i ERROR_REACH;

	
	/* R0 and R1 store two latest computed consecutive antidiagonals */
	R0			= _mm_setzero_si128();
	R1			= _mm_setzero_si128();

	/* Stores newly computed matching at each iteration */
	Diag		= _mm_setzero_si128();

	/* These 6 variables are for the edge of the banded diagonal */
	Side1		= _mm_setzero_si128();
	Side2		= _mm_setzero_si128();
	Side		= _mm_setzero_si128();
	Down1		= _mm_setzero_si128();
	Down2		= _mm_setzero_si128();
	Down		= _mm_setzero_si128();

	/* Used for if the error threshold has been exceeded. */
	ERROR_REACH = _mm_setzero_si128();
	tmp			= _mm_setzero_si128();
	
	/* e = 2, if length of the read sequence is <= 2, just run small edit
	 * distance without bothering. */
	if (lenb <= e)
	{
	  	return smallEditDistanceB(a, lena, b, lenb);
	}

	/* Initialization */
	ERROR_REACH = _mm_set_epi16(0, 0, 0, 0, 0, e, e, e);
	
	R0			= _mm_insert_epi16(R0,0,0);
	R1			= _mm_insert_epi16(R1,1,0);
	R1			= _mm_insert_epi16(R1,1,1);
	
	error		= ((a[0]) != (b[0]));

	/* Why didn't you use _mm_set_epi16 as you did in the forward version? */
	Diag		= _mm_insert_epi16(Diag,2*e,0);
	Diag		= _mm_insert_epi16(Diag,error,1);
	Diag		= _mm_insert_epi16(Diag,2*e,2);
	
	Side1		= _mm_insert_epi16(Side1,1,0);
	Side1		= _mm_insert_epi16(Side1,1,1);
	Side1		= _mm_insert_epi16(Side1,2*e,2);
	
	Side2		= _mm_insert_epi16(Side2,2*e,0);
	Side2		= _mm_insert_epi16(Side2,1,1);
	Side2		= _mm_insert_epi16(Side2,1,2);
	
	Down1		= _mm_insert_epi16(Down1,1,0);
	Down1		= _mm_insert_epi16(Down1,1,1);
	Down1		= _mm_insert_epi16(Down1,2*e,2);
	
	Down2		= _mm_insert_epi16(Down2,2*e,0);
	Down2		= _mm_insert_epi16(Down2,1,1);
	Down2		= _mm_insert_epi16(Down2,1,2);
	
	tmp	= _mm_slli_si128(R1,2);	
	R0	= _mm_min_epi16(R1 + Side1, R0 + Diag);
	R0	= _mm_min_epi16(R0, tmp + Down2);
	
	for (i = 3; i < 2 * lena; i++)
	{		
		if (i % 2 == 1)			/* Odd numbered indexes - R1 */
		{
			/* Diag stores the currently computed two values */
			/* We move backwards (subtraction from the given character.) */
			Diag  = _mm_sub_epi8(Diag, Diag);
			error = (*(a - ((i + 1) / 2 - 1)) != *(b - ((i - 1) / 2 - 1)));
			Diag  = _mm_insert_epi16(Diag,error,0);
			
			error = (*(a - ((i - 1) / 2 - 1)) != *(b - ((i + 1) / 2 - 1)));
			Diag  = _mm_insert_epi16(Diag,error,1);

			/* Update two entries of R1 (stored in the 32 LSB bits of R1) */
			tmp	= _mm_srli_si128(R0,2);
			R1	= _mm_min_epi16(tmp + Side1, R1 + Diag);
			R1	= _mm_min_epi16(R1, R0 + Down1);

			if (i > 2 * lenb - 2)
			{
		  		i1		   = _mm_extract_epi16(R1, 1);
		  		totalError = min(totalError, i1);
			}
		}
		else if (i % 2 == 0)	/* Even numbered indexes - R0 */
		{
			/* Diag stores the currently computed three values */
			/* Again we are going backwards. */
			error = (*(a - (i / 2)) != *(b - (i / 2 - 2)));
			Diag  = _mm_insert_epi16(Diag,error,0);
			
			error = (*(a - (i / 2 - 1)) != *(b - (i / 2 - 1)));
			Diag  = _mm_insert_epi16(Diag,error,1);
			
			error = (*(a - (i / 2 - 2)) != *(b - (i / 2)));
			Diag  = _mm_insert_epi16(Diag,error,2);

			tmp = _mm_slli_si128(R1,2);

			/* Update three entries of R0 (stored in the 48 LSB bits of R0) */
			R0	= _mm_min_epi16(R1 + Side1, R0 + Diag);
			R0	= _mm_min_epi16(R0, tmp + Down2);
			tmp = _mm_sub_epi16(ERROR_REACH, R0);
			i0	= _mm_movemask_epi8(tmp);

			/* Quit if we exceed the error threshold, it checks R1, and R0. */
			if (i0 == 63 &&
				_mm_extract_epi16(R1,0) > errThreshold &&
				_mm_extract_epi16(R1,1) > errThreshold &&
				i < 2 * lenb - 2)	/* Did not reach the end yet */
			{
		  		return -1;
			}

			if (i == 2 * lenb - 2)
			{
		  		totalError = _mm_extract_epi16(R0, 2);
			}
		}
	}
	
	Down1 = _mm_insert_epi16(Down1,2*e,0);
	
	/* ********************************************************************* */
	/* First part of the error */
	error = (*(a - (i / 2)) != *(b - (i / 2 - 2)));
	Diag  = _mm_insert_epi16(Diag,error,0);
	
	error = (*(a - (i / 2 - 1)) != *(b - (i / 2 - 1)));
	Diag = _mm_insert_epi16(Diag,error,1);
	
	Diag = _mm_insert_epi16(Diag,2*e,2);

	R0 = _mm_min_epi16(R1 + Side1, R0 + Diag);
	R0 = _mm_min_epi16(R0, _mm_slli_si128(R1,2) + Down1);

	i0 = _mm_extract_epi16(R0, 0); /* Not used here */
	i1 = _mm_extract_epi16(R0, 1);

	totalError = min(totalError, i1);

	
	/* ********************************************************************* */
	/* Second part of the error */
	i++;
	
	Diag = _mm_sub_epi8(Diag, Diag);
	Diag = _mm_insert_epi16(Diag,2*e,0);
	error = (*(a - (i / 2)) != *(b - (lenb - 1)));
	
	Diag = _mm_insert_epi16(Diag,error,1);
	Diag = _mm_insert_epi16(Diag,2*e,2);

	R1 = _mm_min_epi16(R0 + Side1, _mm_slli_si128(R1,2) + Diag);
	R1 = _mm_min_epi16(R1, _mm_slli_si128(R0,2) + Down1);
	i0 = _mm_extract_epi16(R1, 0); /* Not used here */
	i1 = _mm_extract_epi16(R1, 1);

	totalError = min(totalError, i1);

	
	/* ********************************************************************* */
	/* Final element of the matrix */
	i++;

	Diag = _mm_sub_epi8(Diag, Diag);
	error = (*(a - (i / 2)) != *(b - (lenb - 1)));
	
	Diag = _mm_insert_epi16(Diag,error,0);

	Down = _mm_insert_epi16(Down,1,0);
	Side = _mm_insert_epi16(Side,1,0);
	tmp = _mm_srli_si128(R1,2);

	R0 = _mm_min_epi16(R1 + Down, _mm_srli_si128(R0,2) + Diag);
	R0 = _mm_min_epi16(R0, tmp + Side);
	i0 = _mm_extract_epi16(R0, 0);
	totalError = min(totalError, i0);

	/* If final error is greater than error threshold, return
	 * -1  */
	if (totalError > e)
	  	return -1;
	
	return totalError;
}



/*------------------------------------------------------------------------------
 * -FULLY COMMENTED-
 * 	a		: ref genome
 *	lena	: length of current portion of the ref genome that matching will be
 *			  performed
 *	b		: read sequence
 * 	lenb	: length the of the read sequence
 
 * 	Called for e = 3 (verifySingleEndEditDistance)
 -----------------------------------------------------------------------------*/
int
forwardEditDistanceSSE2Odd ( char *a, int lena, char *b, int lenb )
{
	if (lenb == 0 || lena == 0)
	  	return 0;

	int		i		 = 0;
	int		j		 = 0;
	int		k		 = 0;
	int		e		 = errThreshold;
	int		minError = 2 * e;
	char	flag	 = 0;

	/* e = 3, if length of the read sequence is <= 3, just run small edit
	 * distance without bothering. */
	if (lenb <= e)
	{
	  	return smallEditDistanceF(a, lena, b, lenb);
	}

	/* 128-bit variables */
	/* R0 stores the result of even-numbered antidiagonals, R1 stores the result
	 * of odd-numbered antidiagonals. */
	__m128i R0, R1;

	/* These are helping variables, while computing main diagonal of the edit
	 * distance matrix. Diag stores the newly computed edit distance result. */
	__m128i Diag;
	__m128i Side1, Side2;
	__m128i Down1, Down2;
	__m128i Error;
	__m128i tmp;

	/* R0 and R1 store two latest computed consecutive antidiagonals */
	R0	  = _mm_setzero_si128();
	R1	  = _mm_setzero_si128();

	/* These 6 variables are for the edge of the banded diagonal */
	Diag  = _mm_setzero_si128();
	Side1 = _mm_setzero_si128();
	Side2 = _mm_setzero_si128();
	Down1 = _mm_setzero_si128();
	Down2 = _mm_setzero_si128();

	/* Used for if the error threshold has been exceeded. */
	Error = _mm_setzero_si128();
	tmp	  = _mm_setzero_si128();
	
	/* I thought these were already set to 0 */
	R1	  = _mm_xor_si128(R1, R1);
	R0	  = _mm_xor_si128(R0, R0);
	Diag  = _mm_xor_si128(Diag, Diag);
	Side1 = _mm_xor_si128(Side1, Side1);
	Down1 = _mm_xor_si128(Down1, Down1);

	/* Initialize helping variables */
	Diag  = _mm_insert_epi16(Diag,2*e,0);
	Side1 = _mm_insert_epi16(Side1,1,0);
	Side1 = _mm_insert_epi16(Side1,2*e,1);
	Down1 = _mm_insert_epi16(Down1,2*e,0);
	Down1 = _mm_insert_epi16(Down1,1,1);
	Down1 = _mm_insert_epi16(Down1,2*e,2);

	/* Initialize R0 and R1 */
	R0 = _mm_insert_epi16(R0,0,0);
	R1 = _mm_insert_epi16(R1,1,0);
	R1 = _mm_insert_epi16(R1,1,1);

	/* This is not the main loop yet, just goes up to 'e'. Computes R0 and R1
	 * till the i becomes large enough to handle a pattern. */
	for (i = 2; i <= e; i++)
	{		
		/* In each iteration, shift Side1 left by 16 bits and insert 1 to the
		 * least significant 16 bits:
		   	0	0	0	0	0	0	2e	1	before starting loop
			0	0	0	0	0	2e	1	1
			0	0	0	0	2e	1	1	1
			0	0	0	2e	1	1	1	1
			and so on
		*/
		Side1 = _mm_slli_si128(Side1,2);
		Side1 = _mm_insert_epi16(Side1,1,0);

		/* In each iteration, insert 1 to the least significant 16 bits of
		 * Down1, shift left by 16 bits, then insert 2e to the least significant
		 * 16 bits:
		 	0	0	0	0	0	2e	1	2e	before starting loop
			0	0	0	0	2e	1	1	2e
			0	0	0	2e	1	1	1	2e
			0	0	2e	1	1	1	1	2e
			and so on

		 */
		Down1 = _mm_insert_epi16(Down1,1,0);
		Down1 = _mm_slli_si128(Down1,2);
		Down1 = _mm_insert_epi16(Down1,2*e,0);

		/* Diag = 	0	0	0	0	0	0	0	0	0 */
		Diag = _mm_xor_si128(Diag, Diag);
		if (i % 2 == 0)			/* EVEN */
		{	
			Diag = _mm_insert_epi16(Diag,2*e,0);

			for (j = 1; j <= i - 1; j++)
			{	
		  		Diag = _mm_slli_si128(Diag, 2);
		  		Diag = _mm_insert_epi16(Diag, b[i/2-1+(i/2-j)] != a[i/2-1-(i/2-j)],0);
			}
			
			Diag = _mm_slli_si128(Diag, 2);
			Diag = _mm_insert_epi16(Diag, 2*e,0);
			
			R0 = _mm_min_epi16(R1 + Side1, _mm_slli_si128(R0,2) + Diag);
			R0 = _mm_min_epi16(R0, _mm_slli_si128(R1,2) + Down1);
		}
		else					/* ODD */
		{
			Diag = _mm_insert_epi16(Diag,2*e,0);
			
			for (j = i / 2 - 1; j >= -i / 2; j--)
			{
		  		Diag = _mm_slli_si128(Diag, 2);
		  		Diag = _mm_insert_epi16(Diag, b[(i+1)/2+j-1] != a[(i-1)/2-j-1],0);
			}
			
			Diag = _mm_slli_si128(Diag, 2);
			Diag = _mm_insert_epi16(Diag, 2*e,0);
			
			R1 = _mm_min_epi16(R0 + Side1, _mm_slli_si128(R1,2) + Diag);
			R1 = _mm_min_epi16(R1, _mm_slli_si128(R0,2) + Down1);
		}
	}

	/* Set all to 0 */
	Error = _mm_xor_si128(Error, Error);
	Side2 = _mm_xor_si128(Side2, Side2);
	Side1 = _mm_xor_si128(Side1, Side1);
	Down2 = _mm_xor_si128(Down2, Down2);
	Down1 = _mm_xor_si128(Down1, Down1);

	/* Initialize helping variables */
	Error = _mm_insert_epi16(Error,e,0);
	Side2 = _mm_insert_epi16(Side2,2*e,0);
	Side1 = _mm_insert_epi16(Side2,2*e,0);
	Down1 = _mm_insert_epi16(Down1,2*e,0);

	for (j = 0; j < e; j++)
	{
		Side2 = _mm_slli_si128(Side2, 2);
		Side2 = _mm_insert_epi16(Side2,1,0);

		Side1 = _mm_slli_si128(Side1, 2);
		Side1 = _mm_insert_epi16(Side1,1,0);

		Down1 = _mm_slli_si128(Down1, 2);
		Down1 = _mm_insert_epi16(Down1,1,0);

		Down2 = _mm_slli_si128(Down2, 2);
		Down2 = _mm_insert_epi16(Down2,1,0);

		Error = _mm_slli_si128(Error, 2);
		Error = _mm_insert_epi16(Error, e, 0);
	}

	Down2 = _mm_slli_si128(Down2, 2);
	Down2 = _mm_insert_epi16(Down2,2*e,0);
	
	/* We continue with the old i value = (e + 1) */
	for (; i <= 2 * lenb - (e - 1); i++)
	{
		flag = 0;
		Diag = _mm_xor_si128(Diag, Diag);
		if (i % 2 == 0)			/* EVEN */
		{			
			for (j = e / 2; j >= -e / 2; j--)
			{
		  		Diag = _mm_slli_si128(Diag, 2);
		  		Diag = _mm_insert_epi16(Diag, b[i/2-1+j] != a[i/2-1-j],0);
			}			

			R0 = _mm_min_epi16(_mm_srli_si128(R1,2) + Side1, R0 + Diag);
			R0 = _mm_min_epi16(R0, R1 + Down1);			

			/* Continue to process if we did not exceed the error threshold. */
			if (_mm_extract_epi16(R0,0) <= e)
		  		flag = 1;

			/* Check all values in R0 if we did exceed the error threshold or
			 * not. */
			tmp = _mm_srli_si128(R0,2);
			for (j = 0; j < e - 1; j++)
			{
				if (_mm_extract_epi16(tmp,0) <= e)
				  	flag = 1;
				tmp = _mm_srli_si128(tmp,2);
			}

			/* We seem to be doing the same thing for R1. */
			/* Farhad 28/02/2012 */
			if (_mm_extract_epi16(R1,0) <= e)
		  		flag = 1;

			tmp = _mm_srli_si128(R1,2);
			for (j = 0; j < e - 1; j++)
			{
				if (_mm_extract_epi16(tmp,0) <= e)
				  	flag = 1;
				tmp = _mm_srli_si128(tmp,2);
			}
			/* Farhad 28/02/2012 end */
			
			if (flag == 0)
		  		return -1;

			/* If we hit the end, update minError */
			if (i == 2 * lenb - (e - 1))
			{
				tmp = _mm_srli_si128(R0,2);
				for (k = 0; k < e - 2; k++)
				  	tmp = _mm_srli_si128(tmp,2);
				minError = _mm_extract_epi16(tmp,0);
			}
		}
		else					/* ODD */
		{		
			for (j = e / 2; j >= -e / 2 - 1; j--)
			{
		  		Diag = _mm_slli_si128(Diag, 2);
		  		Diag = _mm_insert_epi16(Diag, b[(i+1)/2+j-1] != a[(i)/2-j-1],0);
			}

			R1 = _mm_min_epi16(R0 + Side2, R1 + Diag);
			R1 = _mm_min_epi16(R1, _mm_slli_si128(R0,2) + Down2);

			/* If we reached the end, update the error value. */
			if (i >= 2 * lenb - e)
			{
		  		tmp = _mm_srli_si128(R1,2);
		  		for (k = 0; k < e - 1; k++)
					tmp = _mm_srli_si128(tmp,2);
		  		minError = min(minError, _mm_extract_epi16(tmp,0));
			}
		}
	}

	/* Compute the remaining cells manually */
	
	//first cell
	Diag	 = _mm_xor_si128(Diag, Diag);
	Diag	 = _mm_insert_epi16(Diag, b[lenb-3] != a[lena], 0);
	Diag	 = _mm_insert_epi16(Diag, b[lenb-2] != a[lena-1], 1);
	Diag	 = _mm_insert_epi16(Diag, b[lenb-1] != a[lena-2], 2);
	Diag	 = _mm_insert_epi16(Diag, 2*e, 3);
	R1		 = _mm_min_epi16(R0 + Side2, R1 + Diag);
	R1		 = _mm_min_epi16(R1, _mm_slli_si128(R0,2) + Down2);
	minError = min(minError, _mm_extract_epi16(R1,2));

	//second cell
	Diag	 = _mm_xor_si128(Diag, Diag);
	Diag	 = _mm_insert_epi16(Diag, b[lenb-2] != a[lena], 0);
	Diag	 = _mm_insert_epi16(Diag, b[lenb-1] != a[lena-1], 1);
	Diag	 = _mm_insert_epi16(Diag, 2*e, 2);
	R0		 = _mm_min_epi16(_mm_srli_si128(R1,2) + Side1, R0 + Diag);
	R0		 = _mm_min_epi16(R0, R1 + Down1);
	minError = min(minError, _mm_extract_epi16(R0,1));

	//third cell
	Diag	 = _mm_xor_si128(Diag, Diag);
	Diag	 = _mm_insert_epi16(Diag, b[lenb-2] != a[lena+1], 0);
	Diag	 = _mm_insert_epi16(Diag, b[lenb-1] != a[lena], 1);
	Diag	 = _mm_insert_epi16(Diag, 2*e, 2);
	R1		 = _mm_min_epi16(R0 + Side2, R1 + Diag);
	R1		 = _mm_min_epi16(R1, _mm_slli_si128(R0,2) + Down2);
	minError = min(minError, _mm_extract_epi16(R1,1));

	//forth
	Diag	 = _mm_xor_si128(Diag, Diag);
	Diag	 = _mm_insert_epi16(Diag, b[lenb-1] != a[lena+1], 0);
	Diag	 = _mm_insert_epi16(Diag, 2*e, 1);
	R0		 = _mm_min_epi16(_mm_srli_si128(R1,2) + Side1, R0 + Diag);
	R0		 = _mm_min_epi16(R0, R1 + Down1);
	minError = min(minError, _mm_extract_epi16(R0,0));

	//fifth
	Diag	 = _mm_xor_si128(Diag, Diag);
	Diag	 = _mm_insert_epi16(Diag, b[lenb-1] != a[lena+2], 0);
	Diag	 = _mm_insert_epi16(Diag, 2*e, 1);
	R1		 = _mm_min_epi16(R0 + Side2, R1 + Diag);
	R1		 = _mm_min_epi16(R1, _mm_slli_si128(R0,2) + Down2);
	minError = min(minError, _mm_extract_epi16(R1,0));

	/* Check if the error threshold is exceeded */
	if (minError > e)
	  	return -1;
	
	return minError;
}



/*------------------------------------------------------------------------------
 * -SKIPPED-
 * 	a		: ref genome
 *	lena	: length of current portion of the ref genome that matching will be
 *			  performed
 *	b		: read sequence
 * 	lenb	: length the of the read sequence
 
 * 	Called for e = 3 (verifySingleEndEditDistance)
 
 * 	I did not go into details of this as it should be very similar to
 * 	forwardEditDistanceSSE2Odd basically. There should be subtraction instead of
 * 	addition as in the backwardEditDistance2SSE2.
 -----------------------------------------------------------------------------*/
int
backwardEditDistanceSSE2Odd ( char *a, int lena, char *b, int lenb )
{
  if (lenb == 0 || lena == 0)
    return 0;

  int i = 0;
  int j = 0;
  int k = 0;

  int e = errThreshold;

  char flag = 0;

  int minError = 2 * e;

  __m128i R0, R1;
  __m128i Diag;
  __m128i Side1, Side2;
  __m128i Down1, Down2;
  __m128i Error;
  __m128i tmp;

  /* initialize */
  R0 = _mm_setzero_si128();
  R1 = _mm_setzero_si128();
  Diag = _mm_setzero_si128();
  Side1 = _mm_setzero_si128();
  Side2 = _mm_setzero_si128();
  Down1 = _mm_setzero_si128();
  Down2 = _mm_setzero_si128();
  Error = _mm_setzero_si128();
  tmp = _mm_setzero_si128();
  /* end initialize */

  if(lenb <= e) 
    {
      return smallEditDistanceB(a, lena, b, lenb);
    }

  R1 = _mm_xor_si128(R1, R1);
  R0 = _mm_xor_si128(R0, R0);

  Diag = _mm_xor_si128(Diag, Diag);
  Side1 = _mm_xor_si128(Side1, Side1);
  Down1 = _mm_xor_si128(Down1, Down1);

  Diag = _mm_insert_epi16(Diag,2*e,0);

  Side1 = _mm_insert_epi16(Side1,1,0);
  Side1 = _mm_insert_epi16(Side1,2*e,1);

  Down1 = _mm_insert_epi16(Down1,2*e,0);
  Down1 = _mm_insert_epi16(Down1,1,1);
  Down1 = _mm_insert_epi16(Down1,2*e,2);

  R0 = _mm_insert_epi16(R0,0,0);

  R1 = _mm_insert_epi16(R1,1,0);
  R1 = _mm_insert_epi16(R1,1,1);

  for(i = 2; i <= e; i++)
    {
      //set side
      Side1 = _mm_slli_si128(Side1,2);
      Side1 = _mm_insert_epi16(Side1,1,0);

      Down1 = _mm_insert_epi16(Down1,1,0);
      Down1 = _mm_slli_si128(Down1,2);
      Down1 = _mm_insert_epi16(Down1,2*e,0);

      Diag = _mm_xor_si128(Diag, Diag);
      if (i % 2 == 0) {
	Diag = _mm_insert_epi16(Diag,2*e,0);

	for (j = 1; j <= i - 1; j++) {
	  Diag = _mm_slli_si128(Diag, 2);
	  Diag =
	    _mm_insert_epi16(Diag, *(b-(i/2-1+(i/2-j))) != *(a-(i/2-1-(i/2-j))),0);
	}
	Diag = _mm_slli_si128(Diag, 2);
	Diag = _mm_insert_epi16(Diag, 2*e,0);

	R0 = _mm_min_epi16(R1 + Side1, _mm_slli_si128(R0,2) + Diag);
	R0 = _mm_min_epi16(R0, _mm_slli_si128(R1,2) + Down1);
      }

      else {
	Diag = _mm_insert_epi16(Diag,2*e,0);
	for (j = i / 2 - 1; j >= -i / 2; j--) {
	  Diag = _mm_slli_si128(Diag, 2);
	  Diag =
	    _mm_insert_epi16(Diag, *(b-((i+1)/2+j-1)) != *(a-((i-1)/2-j-1)),0);
	}
	Diag = _mm_slli_si128(Diag, 2);
	Diag = _mm_insert_epi16(Diag, 2*e,0);

	R1 = _mm_min_epi16(R0 + Side1, _mm_slli_si128(R1,2) + Diag);
	R1 = _mm_min_epi16(R1, _mm_slli_si128(R0,2) + Down1);
      }
    }
  Error = _mm_xor_si128(Error, Error);
  Side2 = _mm_xor_si128(Side2, Side2);
  Down2 = _mm_xor_si128(Down2, Down2);
  Down1 = _mm_xor_si128(Down1, Down1);

  Error = _mm_insert_epi16(Error,e,0);
  Side1 = _mm_insert_epi16(Side2,2*e,0);
  Side2 = _mm_insert_epi16(Side2,2*e,0);
  Down1 = _mm_insert_epi16(Down1,2*e,0);

  for (j = 0; j < e; j++) {
    Side2 = _mm_slli_si128(Side2, 2);
    Side2 = _mm_insert_epi16(Side2,1,0);

    Side1 = _mm_slli_si128(Side1, 2);
    Side1 = _mm_insert_epi16(Side1,1,0);

    Down1 = _mm_slli_si128(Down1, 2);
    Down1 = _mm_insert_epi16(Down1,1,0);

    Down2 = _mm_slli_si128(Down2, 2);
    Down2 = _mm_insert_epi16(Down2,1,0);

    Error = _mm_slli_si128(Error, 2);
    Error = _mm_insert_epi16(Error, e, 0);
  }

  Down2 = _mm_slli_si128(Down2, 2);
  Down2 = _mm_insert_epi16(Down2,2*e,0);

  for (; i <= 2 * lenb - (e - 1); i++) {
    flag = 0;
    Diag = _mm_xor_si128(Diag, Diag);
    if (i % 2 == 0) {
      for (j = e / 2; j >= -e / 2; j--) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag =
	  _mm_insert_epi16(Diag, *(b-(i/2-1+j)) != *(a-(i/2-1-j)),0);
      }

      R0 = _mm_min_epi16(_mm_srli_si128(R1,2) + Side1, R0 + Diag);
      R0 = _mm_min_epi16(R0, R1 + Down1);

      if (_mm_extract_epi16(R0,0) <= e)
	flag = 1;
      tmp = _mm_srli_si128(R0,2);
      for (j = 0; j <= e; j++) {
	if (_mm_extract_epi16(tmp,0) <= e)
	  flag = 1;
	tmp = _mm_srli_si128(tmp,2);
      }

      if (_mm_extract_epi16(R1,1) <= e)
	flag = 1;
      tmp = _mm_srli_si128(R1,2);
      for (j = 0; j < e-1; j++) {
	if (_mm_extract_epi16(tmp,0) <= e)
	  flag = 1;
	tmp = _mm_srli_si128(tmp,2);
      }

      if (flag == 0)
	return -1;
      /* farhad: lenb1 */
      if (i == 2 * lenb - (e - 1)) {
	tmp = _mm_srli_si128(R0,2);
	for (k = 0; k < e - 2; k++)
	  tmp = _mm_srli_si128(tmp,2);
	minError = _mm_extract_epi16(tmp,0);
      }

    }

    else {
      for (j = e / 2; j >= -e / 2 - 1; j--) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag =
	  _mm_insert_epi16(Diag, *(b-((i+1)/2+j-1)) != *(a-((i)/2-j-1)),0);
      }

      R1 = _mm_min_epi16(R0 + Side2, R1 + Diag);

      R1 = _mm_min_epi16(R1, _mm_slli_si128(R0,2) + Down2);
      /* farhad: lenb2 */
      if (i >= 2 * lenb - e) {
	tmp = _mm_srli_si128(R1,2);
	for (k = 0; k < e - 1; k++)
	  tmp = _mm_srli_si128(tmp,2);
	minError = min(minError, _mm_extract_epi16(tmp,0));
      }
    }
  }

  //first cell
  Diag = _mm_xor_si128(Diag, Diag);
  Diag = _mm_insert_epi16(Diag, *(b-(lenb-3)) != *(a-lena), 0);
  Diag = _mm_insert_epi16(Diag, *(b-(lenb-2)) != *(a-(lena-1)), 1);
  Diag = _mm_insert_epi16(Diag, *(b-(lenb-1)) != *(a-(lena-2)), 2);
  Diag = _mm_insert_epi16(Diag, 2*e, 3);
  R1 = _mm_min_epi16(R0 + Side2, R1 + Diag);
  R1 = _mm_min_epi16(R1, _mm_slli_si128(R0,2) + Down2);

  minError = min(minError, _mm_extract_epi16(R1,2));

  //second cell
  Diag = _mm_xor_si128(Diag, Diag);
  Diag = _mm_insert_epi16(Diag, *(b-(lenb-2)) != *(a-(lena)), 0);
  Diag = _mm_insert_epi16(Diag, *(b-(lenb-1)) != *(a-(lena-1)), 1);
  Diag = _mm_insert_epi16(Diag, 2*e, 2);

  R0 = _mm_min_epi16(_mm_srli_si128(R1,2) + Side1, R0 + Diag);
  R0 = _mm_min_epi16(R0, R1 + Down1);

  minError = min(minError, _mm_extract_epi16(R0,1));

  //third cell
  Diag = _mm_xor_si128(Diag, Diag);
  Diag = _mm_insert_epi16(Diag, *(b-(lenb-2)) != *(a-(lena+1)), 0);
  Diag = _mm_insert_epi16(Diag, *(b-(lenb-1)) != *(a-(lena)), 1);
  Diag = _mm_insert_epi16(Diag, 2*e, 2);

  R1 = _mm_min_epi16(R0 + Side2, R1 + Diag);
  R1 = _mm_min_epi16(R1, _mm_slli_si128(R0,2) + Down2);

  minError = min(minError, _mm_extract_epi16(R1,1));

  //forth
  Diag = _mm_xor_si128(Diag, Diag);
  Diag = _mm_insert_epi16(Diag, *(b-(lenb-1)) != *(a-(lena+1)), 0);
  Diag = _mm_insert_epi16(Diag, 2*e, 1);

  R0 = _mm_min_epi16(_mm_srli_si128(R1,2) + Side1, R0 + Diag);
  R0 = _mm_min_epi16(R0, R1 + Down1);

  minError = min(minError, _mm_extract_epi16(R0,0));

  //fifth
  Diag = _mm_xor_si128(Diag, Diag);
  Diag = _mm_insert_epi16(Diag, *(b-(lenb-1)) != *(a-(lena+2)), 0);
  Diag = _mm_insert_epi16(Diag, 2*e, 1);

  R1 = _mm_min_epi16(R0 + Side2, R1 + Diag);
  R1 = _mm_min_epi16(R1, _mm_slli_si128(R0,2) + Down2);

  minError = min(minError, _mm_extract_epi16(R0,0));

  if (minError > e)
    return -1;
  return minError;
}



/*------------------------------------------------------------------------------
 * -SKIPPED-
 * 	a		: ref genome
 *	lena	: length of current portion of the ref genome that matching will be
 *			  performed
 *	b		: read sequence
 * 	lenb	: length the of the read sequence
 
 * 	Called for e = 3 (verifySingleEndEditDistance) when the error threshold is
 * 	even. This should not be called in any case.
 *----------------------------------------------------------------------------*/
int
backwardEditDistanceSSE2G ( char *a, int lena, char *b, int lenb )
{
  if (lenb == 0 || lena == 0)
    return 0;

  int i = 0;
  int j = 0;
  int k = 0;

  int e = errThreshold;

  char flag = 0;

  int minError = 2 * e;

  __m128i R0, R1;
  __m128i Diag;
  __m128i Side1, Side2;
  __m128i Down1, Down2;
  __m128i Error;
  __m128i tmp;

  /* initialize */
  R0 = _mm_setzero_si128();
  R1 = _mm_setzero_si128();
  Diag = _mm_setzero_si128();
  Side1 = _mm_setzero_si128();
  Side2 = _mm_setzero_si128();
  Down1 = _mm_setzero_si128();
  Down2 = _mm_setzero_si128();
  Error = _mm_setzero_si128();
  tmp = _mm_setzero_si128();
  /* end initialize */

  if (lenb <= e) {
    return smallEditDistanceB(a, lena, b, lenb);
  }

  R1 = _mm_xor_si128(R1, R1);
  R0 = _mm_xor_si128(R0, R0);

  Diag = _mm_xor_si128(Diag, Diag);
  Side1 = _mm_xor_si128(Side1, Side1);
  Down1 = _mm_xor_si128(Down1, Down1);

  Diag = _mm_insert_epi16(Diag,2*e,0);

  Side1 = _mm_insert_epi16(Side1,1,0);
  Side1 = _mm_insert_epi16(Side1,2*e,1);

  Down1 = _mm_insert_epi16(Down1,2*e,0);
  Down1 = _mm_insert_epi16(Down1,1,1);
  Down1 = _mm_insert_epi16(Down1,2*e,2);

  R0 = _mm_insert_epi16(R0,0,0);

  R1 = _mm_insert_epi16(R1,1,0);
  R1 = _mm_insert_epi16(R1,1,1);

  for (i = 2; i <= e; i++) {
    //set side
    Side1 = _mm_slli_si128(Side1,2);
    Side1 = _mm_insert_epi16(Side1,1,0);

    Down1 = _mm_insert_epi16(Down1,1,0);
    Down1 = _mm_slli_si128(Down1,2);
    Down1 = _mm_insert_epi16(Down1,2*e,0);

    Diag = _mm_xor_si128(Diag, Diag);
    if (i % 2 == 0) {
      Diag = _mm_insert_epi16(Diag,2*e,0);

      for (j = 1; j <= i - 1; j++) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag =
	  _mm_insert_epi16(Diag, *(b-(i/2-1+(i/2-j))) != *(a-(i/2-1-(i/2-j))),0);
      }
      Diag = _mm_slli_si128(Diag, 2);
      Diag = _mm_insert_epi16(Diag, 2*e,0);

      R0 = _mm_min_epi16(R1 + Side1, _mm_slli_si128(R0,2) + Diag);
      R0 = _mm_min_epi16(R0, _mm_slli_si128(R1,2) + Down1);
    }

    else {
      Diag = _mm_insert_epi16(Diag,2*e,0);
      for (j = i / 2 - 1; j >= -i / 2; j--) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag =
	  _mm_insert_epi16(Diag, *(b-((i+1)/2+j-1)) != *(a-((i-1)/2-j-1)),0);
      }
      Diag = _mm_slli_si128(Diag, 2);
      Diag = _mm_insert_epi16(Diag, 2*e,0);

      R1 = _mm_min_epi16(R0 + Side1, _mm_slli_si128(R1,2) + Diag);
      R1 = _mm_min_epi16(R1, _mm_slli_si128(R0,2) + Down1);
    }
  }
  Error = _mm_xor_si128(Error, Error);
  Side2 = _mm_xor_si128(Side2, Side2);
  Down2 = _mm_xor_si128(Down2, Down2);
  Down1 = _mm_xor_si128(Down1, Down1);

  Error = _mm_insert_epi16(Error,e,0);
  Side2 = _mm_insert_epi16(Side2,2*e,0);
  Down1 = _mm_insert_epi16(Down1,2*e,0);

  for (j = 0; j < e; j++) {
    Side2 = _mm_slli_si128(Side2, 2);
    Side2 = _mm_insert_epi16(Side2,1,0);

    Down1 = _mm_slli_si128(Down1, 2);
    Down1 = _mm_insert_epi16(Down1,1,0);

    Down2 = _mm_slli_si128(Down2, 2);
    Down2 = _mm_insert_epi16(Down2,1,0);

    Error = _mm_slli_si128(Error, 2);
    Error = _mm_insert_epi16(Error, e, 0);
  }

  Down2 = _mm_slli_si128(Down2, 2);
  Down2 = _mm_insert_epi16(Down2,2*e,0);

  for (; i <= 2 * lenb - (e - 1); i++) {
    flag = 0;
    Diag = _mm_xor_si128(Diag, Diag);
    if (i % 2 == 0) {
      for (j = e / 2; j >= -e / 2; j--) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag =
	  _mm_insert_epi16(Diag, *(b-(i/2-1+j)) != *(a-(i/2-1-j)),0);
      }

      R0 = _mm_min_epi16(R1 + Side2, R0 + Diag);
      R0 = _mm_min_epi16(R0, _mm_slli_si128(R1,2) + Down2);

      if (_mm_extract_epi16(R0,0) <= e)
	flag = 1;
      tmp = _mm_srli_si128(R0,2);
      for (j = 0; j <= e; j++) {
	if (_mm_extract_epi16(tmp,0) <= e)
	  flag = 1;
	tmp = _mm_srli_si128(tmp,2);
      }

      if (flag == 0)
	return -1;

      if (i == 2 * lenb - e) {
	tmp = _mm_srli_si128(R0,2);
	for (k = 0; k < e - 1; k++)
	  tmp = _mm_srli_si128(tmp,2);
	minError = _mm_extract_epi16(tmp,0);
      }

    }

    else {
      for (j = -e / 2 + 1; j <= e / 2; j++) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag =
	  _mm_insert_epi16(Diag, *(b-((i+1)/2-j-1)) != *(a-((i-1)/2+j-1)),0);
      }

      R1 = _mm_min_epi16(_mm_srli_si128(R0,2) + Side1, R1 + Diag);
      R1 = _mm_min_epi16(R1, R0 + Down1);

      if (i >= 2 * lenb - e) {
	tmp = _mm_srli_si128(R1,2);
	for (k = 0; k < e - 2; k++)
	  tmp = _mm_srli_si128(tmp,2);
	minError = min(minError, _mm_extract_epi16(tmp,0));
      }
    }
  }

  j = 0;
  int tmpE = e;
  for (; j < 2 * (e - 2) + 1; j++) {

    Diag = _mm_xor_si128(Diag, Diag);
    //set the first element
    if (j == 0) {
      for (k = 0; k <= e - 1; k++) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag =
	  _mm_insert_epi16(Diag, *(b-(lenb-1-k)) != *(a-((i-lenb)-1+k)),0);
      }

      R0 = _mm_min_epi16(R1 + Side2, R0 + Diag);
      R0 = _mm_min_epi16(R0, _mm_slli_si128(R1,2) + Down2);

      tmpE--;
      tmp = _mm_srli_si128(R0,2);
      for (k = 0; k < e - 2; k++)
	tmp = _mm_srli_si128(tmp,2);
      minError = min(minError, _mm_extract_epi16(tmp,0));
    } else if (j % 2 == 0) {
      for (k = 0; k < tmpE; k++) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag =
	  _mm_insert_epi16(Diag, *(b-(lenb-1-k)) != *(a-((i-lenb)-1+k)),0);
      }

      R0 = _mm_min_epi16(R1 + Side2, R0 + Diag);
      R0 = _mm_min_epi16(R0, _mm_slli_si128(R1,2) + Down2);

      tmpE--;

      tmp = _mm_srli_si128(R0,2);
      for (k = 0; k < tmpE - 1; k++)
	tmp = _mm_srli_si128(tmp,2);
      minError = min(minError, _mm_extract_epi16(tmp,0));
    }

    else {
      for (k = 0; k < tmpE; k++) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag =
	  _mm_insert_epi16(Diag, *(b-(lenb-1-k)) != *(a-((i-lenb)-1+k)),0);
      }

      R1 = _mm_min_epi16(_mm_srli_si128(R0,2) + Side1, R1 + Diag);
      R1 = _mm_min_epi16(R1, R0 + Down1);

      tmp = _mm_srli_si128(R1,2);
      for (k = 0; k < tmpE - 2; k++)
	tmp = _mm_srli_si128(tmp,2);
      minError = min(minError, _mm_extract_epi16(tmp,0));
    }
    i++;
  }
  //Diag

  Diag = _mm_xor_si128(Diag, Diag);
  Diag = _mm_insert_epi16(Diag, 2*e, 0);
  Diag = _mm_insert_epi16(Diag, *(a-(lenb+e-2)) != *(b-(lenb-1)), 1);

  Side1 = _mm_insert_epi16(Side1,1,0);
  Side1 = _mm_insert_epi16(Side1,1,1);

  Down1 = _mm_insert_epi16(Down1, 2*e, 0);
  Down1 = _mm_insert_epi16(Down1, 1, 1);

  R1 = _mm_min_epi16(R0 + Side1, _mm_slli_si128(R1,2) + Diag);
  R1 = _mm_min_epi16(R1, _mm_slli_si128(R0,2) + Down1);

  minError = min(minError, _mm_extract_epi16(R1,1));

  Diag = _mm_insert_epi16(Diag, *(a-(lenb+e-1)) != *(b-(lenb-1)), 0);
  Down1 = _mm_insert_epi16(Down1, 1, 0);

  R0 = _mm_min_epi16(R1 + Down1, R0 + Diag);
  R0 = _mm_min_epi16(R0, _mm_srli_si128(R1,2) + Side1);

  minError = min(minError, _mm_extract_epi16(R0,0));

  if (minError > e)
    return -1;
  return minError;
}



/*------------------------------------------------------------------------------
 * -SKIPPED-
 * 	a		: ref genome
 *	lena	: length of current portion of the ref genome that matching will be
 *			  performed
 *	b		: read sequence
 * 	lenb	: length the of the read sequence
 
 * 	Called for e = 3 (verifySingleEndEditDistance) when the error threshold is
 * 	even. This should not be called in any case.
 *----------------------------------------------------------------------------*/
int
forwardEditDistanceSSE2G ( char *a, int lena, char *b, int lenb )
{
  if (lenb == 0 || lena == 0)
    return 0;

  int i = 0;
  int j = 0;
  int k = 0;

  int e = errThreshold;

  int minError = 2 * e;

  char flag = 0;

  if (lenb <= e) {
    return smallEditDistanceF(a, lena, b, lenb);
  }

  __m128i R0, R1;
  __m128i Diag;
  __m128i Side1, Side2;
  __m128i Down1, Down2;
  __m128i Error;
  __m128i tmp;

  /* initialize */
  R0 = _mm_setzero_si128();
  R1 = _mm_setzero_si128();
  Diag = _mm_setzero_si128();
  Side1 = _mm_setzero_si128();
  Side2 = _mm_setzero_si128();
  Down1 = _mm_setzero_si128();
  Down2 = _mm_setzero_si128();
  Error = _mm_setzero_si128();
  tmp = _mm_setzero_si128();
  /* end initialize */

  R1 = _mm_xor_si128(R1, R1);
  R0 = _mm_xor_si128(R0, R0);

  Diag = _mm_xor_si128(Diag, Diag);
  Side1 = _mm_xor_si128(Side1, Side1);
  Down1 = _mm_xor_si128(Down1, Down1);

  Diag = _mm_insert_epi16(Diag,2*e,0);

  Side1 = _mm_insert_epi16(Side1,1,0);
  Side1 = _mm_insert_epi16(Side1,2*e,1);

  Down1 = _mm_insert_epi16(Down1,2*e,0);
  Down1 = _mm_insert_epi16(Down1,1,1);
  Down1 = _mm_insert_epi16(Down1,2*e,2);

  R0 = _mm_insert_epi16(R0,0,0);

  R1 = _mm_insert_epi16(R1,1,0);
  R1 = _mm_insert_epi16(R1,1,1);

  for (i = 2; i <= e; i++) {
    //set side
    Side1 = _mm_slli_si128(Side1,2);
    Side1 = _mm_insert_epi16(Side1,1,0);

    Down1 = _mm_insert_epi16(Down1,1,0);
    Down1 = _mm_slli_si128(Down1,2);
    Down1 = _mm_insert_epi16(Down1,2*e,0);

    Diag = _mm_xor_si128(Diag, Diag);
    if (i % 2 == 0) {
      Diag = _mm_insert_epi16(Diag,2*e,0);

      for (j = 1; j <= i - 1; j++) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag =
	  _mm_insert_epi16(Diag, b[i/2-1+(i/2-j)] != a[i/2-1-(i/2-j)],0);
      }
      Diag = _mm_slli_si128(Diag, 2);
      Diag = _mm_insert_epi16(Diag, 2*e,0);

      R0 = _mm_min_epi16(R1 + Side1, _mm_slli_si128(R0,2) + Diag);
      R0 = _mm_min_epi16(R0, _mm_slli_si128(R1,2) + Down1);
    }

    else {
      Diag = _mm_insert_epi16(Diag,2*e,0);
      for (j = i / 2 - 1; j >= -i / 2; j--) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag =
	  _mm_insert_epi16(Diag, b[(i+1)/2+j-1] != a[(i-1)/2-j-1],0);
      }
      Diag = _mm_slli_si128(Diag, 2);
      Diag = _mm_insert_epi16(Diag, 2*e,0);

      R1 = _mm_min_epi16(R0 + Side1, _mm_slli_si128(R1,2) + Diag);
      R1 = _mm_min_epi16(R1, _mm_slli_si128(R0,2) + Down1);
    }
  }
  Error = _mm_xor_si128(Error, Error);
  Side2 = _mm_xor_si128(Side2, Side2);
  Down2 = _mm_xor_si128(Down2, Down2);
  Down1 = _mm_xor_si128(Down1, Down1);

  Error = _mm_insert_epi16(Error,e,0);
  Side2 = _mm_insert_epi16(Side2,2*e,0);
  Down1 = _mm_insert_epi16(Down1,2*e,0);

  for (j = 0; j < e; j++) {
    Side2 = _mm_slli_si128(Side2, 2);
    Side2 = _mm_insert_epi16(Side2,1,0);

    Down1 = _mm_slli_si128(Down1, 2);
    Down1 = _mm_insert_epi16(Down1,1,0);

    Down2 = _mm_slli_si128(Down2, 2);
    Down2 = _mm_insert_epi16(Down2,1,0);

    Error = _mm_slli_si128(Error, 2);
    Error = _mm_insert_epi16(Error, e, 0);
  }

  Down2 = _mm_slli_si128(Down2, 2);
  Down2 = _mm_insert_epi16(Down2,2*e,0);

  for (; i <= 2 * lenb - (e - 1); i++) {
    flag = 0;
    Diag = _mm_xor_si128(Diag, Diag);
    if (i % 2 == 0) {
      for (j = e / 2; j >= -e / 2; j--) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag = _mm_insert_epi16(Diag, b[i/2-1+j] != a[i/2-1-j],0);
      }

      R0 = _mm_min_epi16(R1 + Side2, R0 + Diag);
      R0 = _mm_min_epi16(R0, _mm_slli_si128(R1,2) + Down2);

      if (_mm_extract_epi16(R0,0) <= e)
	flag = 1;

      tmp = _mm_srli_si128(R0,2);
      for (j = 0; j < e - 1; j++) {
	if (_mm_extract_epi16(tmp,0) <= e)
	  flag = 1;
	tmp = _mm_srli_si128(tmp,2);
      }

      if (flag == 0)
	return -1;

      if (i == 2 * lenb - e) {
	tmp = _mm_srli_si128(R0,2);
	for (k = 0; k < e - 1; k++)
	  tmp = _mm_srli_si128(tmp,2);
	minError = _mm_extract_epi16(tmp,0);
      }

    }

    else {
      for (j = -e / 2 + 1; j <= e / 2; j++) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag =
	  _mm_insert_epi16(Diag, b[(i+1)/2-j-1] != a[(i-1)/2+j-1],0);
      }

      R1 = _mm_min_epi16(_mm_srli_si128(R0,2) + Side1, R1 + Diag);
      R1 = _mm_min_epi16(R1, R0 + Down1);

      if (i >= 2 * lenb - e) {
	tmp = _mm_srli_si128(R1,2);
	for (k = 0; k < e - 2; k++)
	  tmp = _mm_srli_si128(tmp,2);
	minError = min(minError, _mm_extract_epi16(tmp,0));
      }
    }
  }

  j = 0;
  int tmpE = e;
  for (; j < 2 * (e - 2) + 1; j++) {

    Diag = _mm_xor_si128(Diag, Diag);
    //set the first element
    if (j == 0) {
      for (k = 0; k <= e - 1; k++) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag = _mm_insert_epi16(Diag, b[lenb-1-k] != a[(i-lenb)-1+k],0);
      }

      R0 = _mm_min_epi16(R1 + Side2, R0 + Diag);
      R0 = _mm_min_epi16(R0, _mm_slli_si128(R1,2) + Down2);

      tmpE--;

      tmp = _mm_srli_si128(R0,2);
      for (k = 0; k < e - 2; k++)
	tmp = _mm_srli_si128(tmp,2);
      minError = min(minError, _mm_extract_epi16(tmp,0));
    } else if (j % 2 == 0) {
      for (k = 0; k < tmpE; k++) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag = _mm_insert_epi16(Diag, b[lenb-1-k] != a[(i-lenb)-1+k],0);
      }

      R0 = _mm_min_epi16(R1 + Side2, R0 + Diag);
      R0 = _mm_min_epi16(R0, _mm_slli_si128(R1,2) + Down2);

      tmpE--;

      tmp = _mm_srli_si128(R0,2);
      for (k = 0; k < tmpE - 1; k++)
	tmp = _mm_srli_si128(tmp,2);
      minError = min(minError, _mm_extract_epi16(tmp,0));
    }

    else {
      for (k = 0; k < tmpE; k++) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag = _mm_insert_epi16(Diag, b[lenb-1-k] != a[(i-lenb)-1+k],0);
      }

      R1 = _mm_min_epi16(_mm_srli_si128(R0,2) + Side1, R1 + Diag);
      R1 = _mm_min_epi16(R1, R0 + Down1);

      tmp = _mm_srli_si128(R1,2);
      for (k = 0; k < tmpE - 1; k++)
	tmp = _mm_srli_si128(tmp,2);
      minError = min(minError, _mm_extract_epi16(tmp,0));
    }
    i++;
  }
  //Diag

  Diag = _mm_xor_si128(Diag, Diag);
  Diag = _mm_insert_epi16(Diag, 2*e, 0);
  Diag = _mm_insert_epi16(Diag, a[lenb+e-2] != b[lenb-1], 1);

  Side1 = _mm_insert_epi16(Side1,1,0);
  Side1 = _mm_insert_epi16(Side1,1,1);

  Down1 = _mm_insert_epi16(Down1, 2*e, 0);
  Down1 = _mm_insert_epi16(Down1, 1, 1);

  R1 = _mm_min_epi16(R0 + Side1, _mm_slli_si128(R1,2) + Diag);
  R1 = _mm_min_epi16(R1, _mm_slli_si128(R0,2) + Down1);

  minError = min(minError, _mm_extract_epi16(R1,1));

  Diag = _mm_insert_epi16(Diag, a[lenb+e-1] != b[lenb-1], 1);
  Down1 = _mm_insert_epi16(Down1, 1, 0);

  R0 = _mm_min_epi16(R1 + Down1, R0 + Diag);
  R0 = _mm_min_epi16(R0, _mm_srli_si128(R1,2) + Side1);

  minError = min(minError, _mm_extract_epi16(R0,0));

  if (minError > e)
    return -1;
  return minError;
}



/*------------------------------------------------------------------------------
 * -SKIPPED-
 * 	a		: ref genome
 *	lena	: length of current portion of the ref genome that matching will be
 *			  performed
 *	b		: read sequence
 * 	lenb	: length the of the read sequence

 *	Called for e = 4 (verifySingleEndEditDistance4)

 * 	The logic is very similar, see forwardEditDistance2SSE2 for detailed
 * 	comments.
 *----------------------------------------------------------------------------*/
int
forwardEditDistance4SSE2 ( char *a, int lena, char *b, int lenb )
{
  if (lenb == 0 || lena == 0)
    return 0;

  int i = 0;
  int j = 0;
  int k = 0;

  int i0 = 0;
  int i1 = 0;
  int i2 = 0;
  int i4 = 0;
  int i5 = 0;

  int e = errThreshold;

  int minError = 2 * e;
  int index = 0;
  int tmpValue = 0;

  if (lenb <= e) {
    return smallEditDistanceF(a, lena, b, lenb);
  }

  register __m128i R0, R1;
  __m128i Diag;
  __m128i Side1, Side2;
  __m128i Down1, Down2;
  __m128i tmp;
  register __m128i SeqA, SeqB;
  __m128i Result;

  __m128i tmpSeqA;
  __m128i tmpSeqB;

  /* initialize */
  R0 = _mm_setzero_si128();
  R1 = _mm_setzero_si128();
  Diag = _mm_setzero_si128();
  Side1 = _mm_setzero_si128();
  Side2 = _mm_setzero_si128();
  Down1 = _mm_setzero_si128();
  Down2 = _mm_setzero_si128();
  SeqA = _mm_setzero_si128();
  SeqB = _mm_setzero_si128();
  Result = _mm_setzero_si128();
  /* end initialize */

  R1 = _mm_xor_si128(R1, R1);
  R0 = _mm_xor_si128(R0, R0);

  Diag = _mm_xor_si128(Diag, Diag);
  Diag = _mm_insert_epi16(Diag,2*e,0);

  i0 = (a[0] != b[0]);
  i1 = min(i0, (a[1]!=b[0])) + 1;
  i2 = min(i0,(a[0]!=b[1])) + 1;

  i0 = min3(i0+(a[1]!=b[1]),i1+1,i2+1);
  i4 = min(i1, (a[2]!=b[0])+1) + 1;
  i5 = min(i2, (a[0]!=b[2])+1) + 1;

  R1 = _mm_insert_epi16(R1, 3, 0);
  R1 = _mm_insert_epi16(R1, i1, 1);
  R1 = _mm_insert_epi16(R1, i2, 2);
  R1 = _mm_insert_epi16(R1, 3, 3);

  R0 = _mm_insert_epi16(R0, 4, 0);
  R0 = _mm_insert_epi16(R0, i4, 1);
  R0 = _mm_insert_epi16(R0, i0, 2);
  R0 = _mm_insert_epi16(R0, i5, 3);
  R0 = _mm_insert_epi16(R0, 4, 4);

  Side2 = _mm_xor_si128(Side2, Side2);
  Down2 = _mm_xor_si128(Down2, Down2);
  Down1 = _mm_xor_si128(Down1, Down1);
  Side1 = _mm_xor_si128(Side1, Side1);

  Side2 = _mm_insert_epi16(Side2,2*e,0);
  Down1 = _mm_insert_epi16(Down1,2*e,0);

  Side1 = _mm_insert_epi16(Side1,1,0);

  index = 0;
  for (j = 0; j < e; j++) {
    Side2 = _mm_slli_si128(Side2, 2);
    Side2 = _mm_insert_epi16(Side2,1,0);

    Down1 = _mm_slli_si128(Down1, 2);
    Down1 = _mm_insert_epi16(Down1,1,0);

    Down2 = _mm_slli_si128(Down2, 2);
    Down2 = _mm_insert_epi16(Down2,1,0);

    Side1 = _mm_slli_si128(Side1, 2);
    Side1 = _mm_insert_epi16(Side1,1,0);

    SeqA = _mm_slli_si128(SeqA, 2);
    SeqB = _mm_slli_si128(SeqB, 2);
    SeqA = _mm_insert_epi16(SeqA,a[index],0);
    SeqB = _mm_insert_epi16(SeqB,b[index],0);
    index++;
  }

  Down2 = _mm_slli_si128(Down2, 2);
  Down2 = _mm_insert_epi16(Down2,2*e,0);

  index = 4;
  i = 5;

  int loopEnd = 2 * lenb - (e - 1);
  for (; i <= loopEnd; i++) {
    //Diag = _mm_xor_si128(Diag, Diag);
    if (i % 2 == 0) {
      tmpSeqA = _mm_slli_si128(SeqA, 2);
      tmpSeqB = _mm_slli_si128(SeqB, 2);
      SeqA = _mm_insert_epi16(tmpSeqA,a[index],0);
      SeqB = _mm_insert_epi16(tmpSeqB,b[index],0);

      index++;

      tmp = _mm_shufflelo_epi16(SeqB,27);
      tmp = _mm_slli_si128(tmp, 2);
      tmpValue = _mm_extract_epi16(tmp, 5);
      tmp = _mm_insert_epi16(tmp, tmpValue, 0);

      Result = _mm_cmpeq_epi16(SeqA, tmp);
      Diag = _mm_andnot_si128(Result, MASK);

      R0 = _mm_min_epi16(R1 + Side2, R0 + Diag);
      R0 = _mm_min_epi16(R0, _mm_slli_si128(R1,2) + Down2);

      if (_mm_extract_epi16(R0, 0) > e && _mm_extract_epi16(R0, 1) > e
	  && _mm_extract_epi16(R0, 2) > e
	  && _mm_extract_epi16(R0, 3) > e
	  && _mm_extract_epi16(R0, 4) > e
	  && _mm_extract_epi16(R1, 0) > e
	  && _mm_extract_epi16(R1, 1) > e
	  && _mm_extract_epi16(R1, 2) > e
	  && _mm_extract_epi16(R1, 3) > e)
	return -1;

      if (i == 2 * lenb - e) {
	tmp = _mm_srli_si128(R0,2);
	for (k = 0; k < e - 1; k++)
	  tmp = _mm_srli_si128(tmp,2);
	minError = _mm_extract_epi16(tmp,0);
      }

    }

    else {
      Result = _mm_cmpeq_epi16(SeqA, _mm_shufflelo_epi16(SeqB,27));
      Diag = _mm_andnot_si128(Result, MASK);

      R1 = _mm_min_epi16(_mm_srli_si128(R0,2) + Side1, R1 + Diag);
      R1 = _mm_min_epi16(R1, R0 + Down1);

      if (i >= 2 * lenb - e) {
	tmp = _mm_srli_si128(R1,2);
	for (k = 0; k < e - 2; k++)
	  tmp = _mm_srli_si128(tmp,2);
	minError = min(minError, _mm_extract_epi16(tmp,0));
      }
    }

  }
  j = 0;
  int tmpE = e;
  for (; j < 2 * (e - 2) + 1; j++) {

    Diag = _mm_xor_si128(Diag, Diag);
    //set the first element
    if (j == 0) {
      for (k = 0; k <= e - 1; k++) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag = _mm_insert_epi16(Diag, b[lenb-1-k] != a[(i-lenb)-1+k],0);
      }

      R0 = _mm_min_epi16(R1 + Side2, R0 + Diag);
      R0 = _mm_min_epi16(R0, _mm_slli_si128(R1,2) + Down2);

      tmpE--;

      tmp = _mm_srli_si128(R0,2);
      for (k = 0; k < e - 2; k++)
	tmp = _mm_srli_si128(tmp,2);
      minError = min(minError, _mm_extract_epi16(tmp,0));
    } else if (j % 2 == 0) {
      for (k = 0; k < tmpE; k++) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag = _mm_insert_epi16(Diag, b[lenb-1-k] != a[(i-lenb)-1+k],0);
      }

      R0 = _mm_min_epi16(R1 + Side2, R0 + Diag);
      R0 = _mm_min_epi16(R0, _mm_slli_si128(R1,2) + Down2);

      tmpE--;

      tmp = _mm_srli_si128(R0,2);
      for (k = 0; k < tmpE - 1; k++)
	tmp = _mm_srli_si128(tmp,2);
      minError = min(minError, _mm_extract_epi16(tmp,0));
    }

    else {
      for (k = 0; k < tmpE; k++) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag = _mm_insert_epi16(Diag, b[lenb-1-k] != a[(i-lenb)-1+k],0);
      }

      R1 = _mm_min_epi16(_mm_srli_si128(R0,2) + Side1, R1 + Diag);
      R1 = _mm_min_epi16(R1, R0 + Down1);

      tmp = _mm_srli_si128(R1,2);
      for (k = 0; k < tmpE - 2; k++)
	tmp = _mm_srli_si128(tmp,2);
      minError = min(minError, _mm_extract_epi16(tmp,0));
    }
    i++;
  }
  //Diag

  Diag = _mm_xor_si128(Diag, Diag);
  Diag = _mm_insert_epi16(Diag, 2*e, 0);
  Diag = _mm_insert_epi16(Diag, a[lenb+e-2] != b[lenb-1], 1);

  Side1 = _mm_insert_epi16(Side1,1,0);
  Side1 = _mm_insert_epi16(Side1,1,1);

  Down1 = _mm_insert_epi16(Down1, 2*e, 0);
  Down1 = _mm_insert_epi16(Down1, 1, 1);

  R1 = _mm_min_epi16(R0 + Side1, _mm_slli_si128(R1,2) + Diag);
  R1 = _mm_min_epi16(R1, _mm_slli_si128(R0,2) + Down1);

  minError = min(minError, _mm_extract_epi16(R1,1));

  Diag = _mm_insert_epi16(Diag, a[lenb+e-1] != b[lenb-1], 0);
  Down1 = _mm_insert_epi16(Down1, 1, 0);

  R0 = _mm_min_epi16(R1 + Down1, R0 + Diag);
  R0 = _mm_min_epi16(R0, _mm_srli_si128(R1,2) + Side1);

  minError = min(minError, _mm_extract_epi16(R0,0));

  if (minError > e)
    return -1;
  return minError;
}



/*------------------------------------------------------------------------------
 * -SKIPPED-
 * 	a		: ref genome
 *	lena	: length of current portion of the ref genome that matching will be
 *			  performed
 *	b		: read sequence
 * 	lenb	: length the of the read sequence

 *	Called for e = 4 (verifySingleEndEditDistance4)

 * 	The logic is very similar, see backwardEditDistance2SSE2 for detailed
 * 	comments.
 *----------------------------------------------------------------------------*/
int
backwardEditDistance4SSE2 ( char *a, int lena, char *b, int lenb )
{
  if (lenb == 0 || lena == 0)
    return 0;

  int i = 0;
  int j = 0;
  int k = 0;

  int i0;
  int i1;
  int i2;
  int i4;
  int i5;

  int e = errThreshold;

  int minError = 2 * e;
  int index = 0;
  int tmpValue = 0;

  if (lenb <= e) {
    return smallEditDistanceB(a, lena, b, lenb);
  }

  __m128i R0, R1;
  __m128i Diag;
  __m128i Side1, Side2;
  __m128i Down1, Down2;
  __m128i tmp;
  __m128i SeqA, SeqB;
  __m128i Result;

  /* initialize */
  R0 = _mm_setzero_si128();
  R1 = _mm_setzero_si128();
  Diag = _mm_setzero_si128();
  Side1 = _mm_setzero_si128();
  Side2 = _mm_setzero_si128();
  Down1 = _mm_setzero_si128();
  Down2 = _mm_setzero_si128();
  SeqA = _mm_setzero_si128();
  SeqB = _mm_setzero_si128();
  Result = _mm_setzero_si128();
  /* end initialize */

  R1 = _mm_xor_si128(R1, R1);
  R0 = _mm_xor_si128(R0, R0);

  Diag = _mm_xor_si128(Diag, Diag);
  Diag = _mm_insert_epi16(Diag,2*e,0);

  i0 = (a[0] != b[0]);
  i1 = min(i0, ( *(a-1)!=*b) ) + 1;
  i2 = min(i0,( a[0] != *(b-1) ) ) + 1;

  i0 = min3( i0+ ( *(a-1)!=*(b-1) ),i1+1,i2+1);
  i4 = min(i1, ( *(a-2)!=b[0] )+1) + 1;
  i5 = min(i2, (a[0] != *(b-2))+1) + 1;

  R1 = _mm_insert_epi16(R1, 3, 0);
  R1 = _mm_insert_epi16(R1, i1, 1);
  R1 = _mm_insert_epi16(R1, i2, 2);
  R1 = _mm_insert_epi16(R1, 3, 3);

  R0 = _mm_insert_epi16(R0, 4, 0);
  R0 = _mm_insert_epi16(R0, i4, 1);
  R0 = _mm_insert_epi16(R0, i0, 2);
  R0 = _mm_insert_epi16(R0, i5, 3);
  R0 = _mm_insert_epi16(R0, 4, 4);

  Side2 = _mm_xor_si128(Side2, Side2);
  Down2 = _mm_xor_si128(Down2, Down2);
  Down1 = _mm_xor_si128(Down1, Down1);
  Side1 = _mm_xor_si128(Side1, Side1);

  Side2 = _mm_insert_epi16(Side2,2*e,0);
  Down1 = _mm_insert_epi16(Down1,2*e,0);

  Side1 = _mm_insert_epi16(Side1,1,0);

  index = 0;
  for (j = 0; j < e; j++) {
    Side2 = _mm_slli_si128(Side2, 2);
    Side2 = _mm_insert_epi16(Side2,1,0);

    Down1 = _mm_slli_si128(Down1, 2);
    Down1 = _mm_insert_epi16(Down1,1,0);

    Down2 = _mm_slli_si128(Down2, 2);
    Down2 = _mm_insert_epi16(Down2,1,0);

    Side1 = _mm_slli_si128(Side1, 2);
    Side1 = _mm_insert_epi16(Side1,1,0);

    SeqA = _mm_slli_si128(SeqA, 2);
    SeqB = _mm_slli_si128(SeqB, 2);
    SeqA = _mm_insert_epi16(SeqA,*(a-index),0);
    SeqB = _mm_insert_epi16(SeqB,*(b-index),0);
    index++;
  }

  Down2 = _mm_slli_si128(Down2, 2);
  Down2 = _mm_insert_epi16(Down2,2*e,0);

  index = 4;
  i = 5;
  int loopEnd = 2 * lenb - (e - 1);
  for (; i <= loopEnd; i++) {

    Diag = _mm_xor_si128(Diag, Diag);
    if (i % 2 == 0) {
      SeqA = _mm_slli_si128(SeqA, 2);
      SeqB = _mm_slli_si128(SeqB, 2);
      SeqA = _mm_insert_epi16(SeqA,*(a-(index)),0);
      SeqB = _mm_insert_epi16(SeqB,*(b-(index)),0);

      index++;

      tmp = _mm_shufflelo_epi16(SeqB,27);
      tmp = _mm_slli_si128(tmp, 2);
      tmpValue = _mm_extract_epi16(tmp, 5);
      tmp = _mm_insert_epi16(tmp, tmpValue, 0);

      Result = _mm_cmpeq_epi16(SeqA, tmp);
      Diag = _mm_andnot_si128(Result, MASK);

      R0 = _mm_min_epi16(R1 + Side2, R0 + Diag);
      R0 = _mm_min_epi16(R0, _mm_slli_si128(R1,2) + Down2);

      if (_mm_extract_epi16(R0, 0) > e && _mm_extract_epi16(R0, 1) > e
	  && _mm_extract_epi16(R0, 2) > e
	  && _mm_extract_epi16(R0, 3) > e
	  && _mm_extract_epi16(R0, 4) > e
	  && _mm_extract_epi16(R1, 0) > e
	  && _mm_extract_epi16(R1, 1) > e
	  && _mm_extract_epi16(R1, 2) > e
	  && _mm_extract_epi16(R1, 3) > e)
	return -1;

      if (i == 2 * lenb - e) {
	tmp = _mm_srli_si128(R0,2);
	for (k = 0; k < e - 1; k++)
	  tmp = _mm_srli_si128(tmp,2);
	minError = _mm_extract_epi16(tmp,0);
      }

    }

    else {
      Result = _mm_cmpeq_epi16(SeqA, _mm_shufflelo_epi16(SeqB,27));
      Diag = _mm_andnot_si128(Result, MASK);

      R1 = _mm_min_epi16(_mm_srli_si128(R0,2) + Side1, R1 + Diag);
      R1 = _mm_min_epi16(R1, R0 + Down1);

      if (i >= 2 * lenb - e) {
	tmp = _mm_srli_si128(R1,2);
	for (k = 0; k < e - 2; k++)
	  tmp = _mm_srli_si128(tmp,2);
	minError = min(minError, _mm_extract_epi16(tmp,0));
      }
    }

  }

  j = 0;

  int tmpE = e;

  for (; j < 2 * (e - 2) + 1; j++) {

    Diag = _mm_xor_si128(Diag, Diag);
    //set the first element
    if (j == 0) {
      for (k = 0; k <= e - 1; k++) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag =
	  _mm_insert_epi16(Diag, *(b-(lenb-1-k)) != *(a-((i-lenb)-1+k)),0);
      }

      R0 = _mm_min_epi16(R1 + Side2, R0 + Diag);
      R0 = _mm_min_epi16(R0, _mm_slli_si128(R1,2) + Down2);

      tmpE--;

      tmp = _mm_srli_si128(R0,2);
      for (k = 0; k < e - 2; k++)
	tmp = _mm_srli_si128(tmp,2);
      minError = min(minError, _mm_extract_epi16(tmp,0));
    } else if (j % 2 == 0) {
      for (k = 0; k < tmpE; k++) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag =
	  _mm_insert_epi16(Diag, *(b-(lenb-1-k)) != *(a-((i-lenb)-1+k)),0);
      }

      R0 = _mm_min_epi16(R1 + Side2, R0 + Diag);
      R0 = _mm_min_epi16(R0, _mm_slli_si128(R1,2) + Down2);

      tmpE--;

      tmp = _mm_srli_si128(R0,2);
      for (k = 0; k < tmpE - 1; k++)
	tmp = _mm_srli_si128(tmp,2);
      minError = min(minError, _mm_extract_epi16(tmp,0));
    }

    else {
      for (k = 0; k < tmpE; k++) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag =
	  _mm_insert_epi16(Diag, *(b-(lenb-1-k)) != *(a-((i-lenb)-1+k)),0);
      }

      R1 = _mm_min_epi16(_mm_srli_si128(R0,2) + Side1, R1 + Diag);
      R1 = _mm_min_epi16(R1, R0 + Down1);

      tmp = _mm_srli_si128(R1,2);
      for (k = 0; k < tmpE - 2; k++)
	tmp = _mm_srli_si128(tmp,2);
      minError = min(minError, _mm_extract_epi16(tmp,0));
    }
    i++;
  }
  //Diag

  Diag = _mm_xor_si128(Diag, Diag);
  Diag = _mm_insert_epi16(Diag, 2*e, 0);
  Diag = _mm_insert_epi16(Diag, *(a-(lenb+e-2)) != *(b-(lenb-1)), 1);

  Side1 = _mm_insert_epi16(Side1,1,0);
  Side1 = _mm_insert_epi16(Side1,1,1);

  Down1 = _mm_insert_epi16(Down1, 2*e, 0);
  Down1 = _mm_insert_epi16(Down1, 1, 1);

  R1 = _mm_min_epi16(R0 + Side1, _mm_slli_si128(R1,2) + Diag);
  R1 = _mm_min_epi16(R1, _mm_slli_si128(R0,2) + Down1);

  minError = min(minError, _mm_extract_epi16(R1,1));

  Diag = _mm_insert_epi16(Diag, *(a-(lenb+e-1)) != *(b-(lenb-1)), 0);
  Down1 = _mm_insert_epi16(Down1, 1, 0);

  R0 = _mm_min_epi16(R1 + Down1, R0 + Diag);
  R0 = _mm_min_epi16(R0, _mm_srli_si128(R1,2) + Side1);

  minError = min(minError, _mm_extract_epi16(R0,0));

  if (minError > e)
    return -1;
  return minError;
}



/*------------------------------------------------------------------------------
 * -FULLY COMMENTED-
 *	a		: ref genome
 *	lena	: length of current portion of the ref genome that matching will be
 *			  performed
 *	b		: read sequence
 * 	lenb	: length the of the read sequence
 
 * 	Called for e > 4 (verifySingleEndEditDistanceExtension)
 -----------------------------------------------------------------------------*/
inline int
forwardEditDistanceSSE2Extension ( char *a, int lena, char *b, int lenb )
{
	if (lenb == 0 || lena == 0)
	  	return 0;

	int i		 = 0;
	int j		 = 0;
	int k		 = 0;
	int i0		 = 0;
	int i1		 = 0;
	int i2		 = 0;
	int i4		 = 0;
	int i5		 = 0;
	int mismatch = errThreshold;
	int e		 = 4;
	int minError = 4 * mismatch + 1;
	int index	 = 0;
	int tmpValue = 0;

	/* e > 4, if length of the read sequence is < 4, just run small edit
	 * distance without bothering. */
	if (lenb <= e)
	{
	  	return smallEditDistanceF(a, lena, b, lenb);
	}


	/* 128-bit variables */
	/* R0 stores the result of even-numbered antidiagonals, R1 stores the result
	 * of odd-numbered antidiagonals. */
	register __m128i	R0, R1;

	/* These are helping variables, while computing main diagonal of the edit
	 * distance matrix. Diag stores the newly computed edit distance result. */
	__m128i				Diag;
	__m128i				Side1, Side2;
	__m128i				Down1, Down2;
	__m128i				tmp;

	/* Yet to come. */
	register __m128i	SeqA, SeqB;
	__m128i				Result;
	__m128i				tmpSeqA;
	__m128i				tmpSeqB;

	
	/* R0 and R1 store two latest computed consecutive antidiagonals */
	R0 = _mm_setzero_si128();
	R1 = _mm_setzero_si128();

	/* These 6 variables are for the edge of the banded diagonal */
	Diag  = _mm_setzero_si128();
	Side1 = _mm_setzero_si128();
	Side2 = _mm_setzero_si128();
	Down1 = _mm_setzero_si128();
	Down2 = _mm_setzero_si128();

	/* SeqA and SeqB stores the integer values of characters A G C T. In each
	 * iteration it is shifted left and new char is added to the least
	 * significant integer.*/
	SeqA   = _mm_setzero_si128();
	SeqB   = _mm_setzero_si128();
	Result = _mm_setzero_si128();
	
	/* Set R1 and R0 to 0. Already done so? */
	R1 = _mm_xor_si128(R1, R1);
	R0 = _mm_xor_si128(R0, R0);

	/* Initialize Diag */
	Diag = _mm_xor_si128(Diag, Diag);
	Diag = _mm_insert_epi16(Diag,minError,0);

	/* Compute edit distance values for the first two anti-diagonals */
	i0 = (a[0] != b[0]);
	i1 = min(i0, (a[1]!=b[0])) + 1;
	i2 = min(i0,(a[0]!=b[1])) + 1;
	i0 = min3(i0+(a[1]!=b[1]),i1+1,i2+1);
	i4 = min(i1, (a[2]!=b[0])+1) + 1;
	i5 = min(i2, (a[0]!=b[2])+1) + 1;

	/* Update R0 and R1 to store the results of first two anti-diagonals. */
	R1 = _mm_insert_epi16(R1, 3, 0);
	R1 = _mm_insert_epi16(R1, i1, 1);
	R1 = _mm_insert_epi16(R1, i2, 2);
	R1 = _mm_insert_epi16(R1, 3, 3);
	R0 = _mm_insert_epi16(R0, 4, 0);
	R0 = _mm_insert_epi16(R0, i4, 1);
	R0 = _mm_insert_epi16(R0, i0, 2);
	R0 = _mm_insert_epi16(R0, i5, 3);
	R0 = _mm_insert_epi16(R0, 4, 4);

	/* Initialize Side and Down helping variables. */
	Side2 = _mm_xor_si128(Side2, Side2);
	Down2 = _mm_xor_si128(Down2, Down2);
	Down1 = _mm_xor_si128(Down1, Down1);
	Side1 = _mm_xor_si128(Side1, Side1);
	Side2 = _mm_insert_epi16(Side2,minError,0);
	Down1 = _mm_insert_epi16(Down1,minError,0);
	Side1 = _mm_insert_epi16(Side1,1,0);

	/* Continue initializing variables wrt error value.  */
	index = 0;
	for (j = 0; j < e; j++)
	{
		Side2 = _mm_slli_si128(Side2, 2);
		Side2 = _mm_insert_epi16(Side2,1,0);

		Down1 = _mm_slli_si128(Down1, 2);
		Down1 = _mm_insert_epi16(Down1,1,0);

		Down2 = _mm_slli_si128(Down2, 2);
		Down2 = _mm_insert_epi16(Down2,1,0);

		Side1 = _mm_slli_si128(Side1, 2);
		Side1 = _mm_insert_epi16(Side1,1,0);

		SeqA = _mm_slli_si128(SeqA, 2);
		SeqB = _mm_slli_si128(SeqB, 2);
		SeqA = _mm_insert_epi16(SeqA,a[index],0);
		SeqB = _mm_insert_epi16(SeqB,b[index],0);
		index++;
	}

	Down2 = _mm_slli_si128(Down2, 2);
	Down2 = _mm_insert_epi16(Down2,minError,0);
	
	index = 4;
	i	  = 5;

	int loopEnd = 2 * lenb - (e - 1);
	for (; i <= loopEnd; i++)
	{	
		if (i % 2 == 0)
		{
			/* Shift SeqA and SeqB 16-bits left and insert a[index] and b[index]
			 * to least-significant integer. */
			tmpSeqA = _mm_slli_si128(SeqA, 2);
			tmpSeqB = _mm_slli_si128(SeqB, 2);
			SeqA	= _mm_insert_epi16(tmpSeqA,a[index],0);
			SeqB	= _mm_insert_epi16(tmpSeqB,b[index],0);

			index++;

			/* Take last four chars in SeqB, reverse them and store that into
			 * tmp. */
			tmp		 = _mm_shufflelo_epi16(SeqB,27);
			tmp		 = _mm_slli_si128(tmp, 2);
			tmpValue = _mm_extract_epi16(tmp, 5);
			tmp		 = _mm_insert_epi16(tmp, tmpValue, 0);

			/* Compare tmp and SeqA to find out Diag. Note that using MASK, we
			 * only compute lower 5 integer values in Diag. In EVEN iterations
			 * we compare 5 values from reference genome and read seq*/
			Result = _mm_cmpeq_epi16(SeqA, tmp);
			Diag   = _mm_andnot_si128(Result, MASK);

			/* Find new value of R0 */
			R0 = _mm_min_epi16(R1 + Side2, R0 + Diag);
			R0 = _mm_min_epi16(R0, _mm_slli_si128(R1,2) + Down2);

		    /* Extract error values and return if any of them exceeds error
			 * threshold. */
			if (_mm_extract_epi16(R0, 0) > errThreshold	&&
				_mm_extract_epi16(R0, 1) > errThreshold	&&
				_mm_extract_epi16(R0, 2) > errThreshold	&&
				_mm_extract_epi16(R0, 3) > errThreshold	&&
				_mm_extract_epi16(R0, 4) > errThreshold	&&
				_mm_extract_epi16(R1, 0) > errThreshold	&&
				_mm_extract_epi16(R1, 1) > errThreshold	&&
				_mm_extract_epi16(R1, 2) > errThreshold	&&
				_mm_extract_epi16(R1, 3) > errThreshold)
		  		return -1;

			/* Update minError if we are at the end of the loop */
			if (i == 2 * lenb - e)
			{
				tmp = _mm_srli_si128(R0,2);
				for (k = 0; k < e - 1; k++)
				  	tmp = _mm_srli_si128(tmp,2);
				minError = _mm_extract_epi16(tmp,0);
			}
		}
		else
		{			
			/* Compare SeqA and shuffled SeqB (last 4 chars reversed) to find
			 * out Diag. Note that using MASK, we only compute lower 5 integer
			 * values in Diag. In ODD iterations we compare 4 values from
			 * reference genome and read seq*/
			Result = _mm_cmpeq_epi16(SeqA, _mm_shufflelo_epi16(SeqB,27));
			Diag   = _mm_andnot_si128(Result, MASK);

			/* Find new value of R1 */
			R1 = _mm_min_epi16(_mm_srli_si128(R0,2) + Side1, R1 + Diag);
			R1 = _mm_min_epi16(R1, R0 + Down1);

			/* Update minError if we are at the end of the loop */	
			if (i >= 2 * lenb - e)
			{
				tmp = _mm_srli_si128(R1,2);
				for (k = 0; k < e - 2; k++)
				  	tmp = _mm_srli_si128(tmp,2);
				minError = min(minError, _mm_extract_epi16(tmp,0));
			}
		}
	}

	/* Fill the remaining elements */
	j		 = 0;
	int tmpE = e;
	for (; j < 2 * (e - 2) + 1; j++)
	{		
		Diag = _mm_xor_si128(Diag, Diag);
		
		//set the first element -- OK
		if (j == 0)
		{
			for (k = 0; k <= e - 1; k++)
			{				
		  		Diag = _mm_slli_si128(Diag, 2);
		  		Diag = _mm_insert_epi16(Diag, b[lenb-1-k] != a[(i-lenb)-1+k],0);
			}

			R0 = _mm_min_epi16(R1 + Side2, R0 + Diag);
			R0 = _mm_min_epi16(R0, _mm_slli_si128(R1,2) + Down2);

			tmpE--;

			tmp = _mm_srli_si128(R0,2);
			for (k = 0; k < e - 2; k++)
		  		tmp = _mm_srli_si128(tmp,2);
			minError = min(minError, _mm_extract_epi16(tmp,0));
		}		
		else if (j % 2 == 0)	/* EVEN -- R0 */
		{
			for (k = 0; k < tmpE; k++)
			{				
		  		Diag = _mm_slli_si128(Diag, 2);
		  		Diag = _mm_insert_epi16(Diag, b[lenb-1-k] != a[(i-lenb)-1+k],0);
			}

			R0 = _mm_min_epi16(R1 + Side2, R0 + Diag);
			R0 = _mm_min_epi16(R0, _mm_slli_si128(R1,2) + Down2);

			tmpE--;

			tmp = _mm_srli_si128(R0,2);
			for (k = 0; k < tmpE - 1; k++)
		  		tmp = _mm_srli_si128(tmp,2);
			minError = min(minError, _mm_extract_epi16(tmp,0));
		}
		else					/* ODD -- R1 */
		{
			for (k = 0; k < tmpE; k++)
			{				
		  		Diag = _mm_slli_si128(Diag, 2);
		  		Diag = _mm_insert_epi16(Diag, b[lenb-1-k] != a[(i-lenb)-1+k],0);
			}

			R1 = _mm_min_epi16(_mm_srli_si128(R0,2) + Side1, R1 + Diag);
			R1 = _mm_min_epi16(R1, R0 + Down1);

			tmp = _mm_srli_si128(R1,2);
			for (k = 0; k < tmpE - 2; k++)
		  		tmp = _mm_srli_si128(tmp,2);
			minError = min(minError, _mm_extract_epi16(tmp,0));
		}

		i++;
	}
	
	//Diag

	
	/* Almost there. */
	Diag = _mm_xor_si128(Diag, Diag);
	Diag = _mm_insert_epi16(Diag, minError, 0);
	Diag = _mm_insert_epi16(Diag, a[lenb+e-2] != b[lenb-1], 1);	

	Side1 = _mm_insert_epi16(Side1,1,0);
	Side1 = _mm_insert_epi16(Side1,1,1);

	Down1 = _mm_insert_epi16(Down1, minError, 0);
	Down1 = _mm_insert_epi16(Down1, 1, 1);

	R1 = _mm_min_epi16(R0 + Side1, _mm_slli_si128(R1,2) + Diag);
	R1 = _mm_min_epi16(R1, _mm_slli_si128(R0,2) + Down1);

	minError = min(minError, _mm_extract_epi16(R1,1));

	Diag = _mm_insert_epi16(Diag, a[lenb+e-1] != b[lenb-1], 0);
	Down1 = _mm_insert_epi16(Down1, 1, 0);

	R0 = _mm_min_epi16(R1 + Down1, R0 + Diag);
	R0 = _mm_min_epi16(R0, _mm_srli_si128(R1,2) + Side1);

	minError = min(minError, _mm_extract_epi16(R0,0));

	/* Check error and return  */
	if (minError > mismatch)
	  	return -1;
	
	return minError;
}



/*------------------------------------------------------------------------------
 * -SKIPPED-
 *	a		: ref genome
 *	lena	: length of current portion of the ref genome that matching will be
 *			  performed
 *	b		: read sequence
 * 	lenb	: length the of the read sequence
 
 * 	Called for e > 4 (verifySingleEndEditDistanceExtension)

 * Did not into the details of this as it should be very similar to
 * forwardEditDistanceSSE2Extension.
 -----------------------------------------------------------------------------*/
inline int
backwardEditDistanceSSE2Extension ( char *a, int lena, char *b, int lenb )
{
  if (lenb == 0 || lena == 0)
    return 0;

  int i = 0;
  int j = 0;
  int k = 0;

  int i0;
  int i1;
  int i2;
  int i4;
  int i5;

  int e = 4;
  int mismatch = errThreshold;

  int minError = 2 * errThreshold;
  int index = 0;
  int tmpValue = 0;

  if (lenb <= e) {
    return smallEditDistanceB(a, lena, b, lenb);
  }

  __m128i R0, R1;
  __m128i Diag;
  __m128i Side1, Side2;
  __m128i Down1, Down2;
  __m128i tmp;
  __m128i SeqA, SeqB;
  __m128i Result;

  /* initialize */
  R0 = _mm_setzero_si128();
  R1 = _mm_setzero_si128();
  Diag = _mm_setzero_si128();
  Side1 = _mm_setzero_si128();
  Side2 = _mm_setzero_si128();
  Down1 = _mm_setzero_si128();
  Down2 = _mm_setzero_si128();
  SeqA = _mm_setzero_si128();
  SeqB = _mm_setzero_si128();
  Result = _mm_setzero_si128();
  /* end initialize */

  R1 = _mm_xor_si128(R1, R1);
  R0 = _mm_xor_si128(R0, R0);

  Diag = _mm_xor_si128(Diag, Diag);
  Diag = _mm_insert_epi16(Diag,minError,0);

  i0 = (a[0] != b[0]);
  i1 = min(i0, ( *(a-1)!=*b) ) + 1;
  i2 = min(i0,( a[0] != *(b-1) ) ) + 1;

  i0 = min3( i0+ ( *(a-1)!=*(b-1) ),i1+1,i2+1);
  i4 = min(i1, ( *(a-2)!=b[0] )+1) + 1;
  i5 = min(i2, (a[0] != *(b-2))+1) + 1;

  R1 = _mm_insert_epi16(R1, 3, 0);
  R1 = _mm_insert_epi16(R1, i1, 1);
  R1 = _mm_insert_epi16(R1, i2, 2);
  R1 = _mm_insert_epi16(R1, 3, 3);

  R0 = _mm_insert_epi16(R0, 4, 0);
  R0 = _mm_insert_epi16(R0, i4, 1);
  R0 = _mm_insert_epi16(R0, i0, 2);
  R0 = _mm_insert_epi16(R0, i5, 3);
  R0 = _mm_insert_epi16(R0, 4, 4);

  Side2 = _mm_xor_si128(Side2, Side2);
  Down2 = _mm_xor_si128(Down2, Down2);
  Down1 = _mm_xor_si128(Down1, Down1);
  Side1 = _mm_xor_si128(Side1, Side1);

  Side2 = _mm_insert_epi16(Side2,minError,0);
  Down1 = _mm_insert_epi16(Down1,minError,0);

  Side1 = _mm_insert_epi16(Side1,1,0);

  index = 0;
  for (j = 0; j < e; j++) {
    Side2 = _mm_slli_si128(Side2, 2);
    Side2 = _mm_insert_epi16(Side2,1,0);

    Down1 = _mm_slli_si128(Down1, 2);
    Down1 = _mm_insert_epi16(Down1,1,0);

    Down2 = _mm_slli_si128(Down2, 2);
    Down2 = _mm_insert_epi16(Down2,1,0);

    Side1 = _mm_slli_si128(Side1, 2);
    Side1 = _mm_insert_epi16(Side1,1,0);

    SeqA = _mm_slli_si128(SeqA, 2);
    SeqB = _mm_slli_si128(SeqB, 2);
    SeqA = _mm_insert_epi16(SeqA,*(a-index),0);
    SeqB = _mm_insert_epi16(SeqB,*(b-index),0);
    index++;
  }

  Down2 = _mm_slli_si128(Down2, 2);
  Down2 = _mm_insert_epi16(Down2,minError,0);

  index = 4;
  i = 5;

  int loopEnd = 2 * lenb - (e - 1);
  for (; i <= loopEnd; i++) {

    Diag = _mm_xor_si128(Diag, Diag);
    if (i % 2 == 0) {
      SeqA = _mm_slli_si128(SeqA, 2);
      SeqB = _mm_slli_si128(SeqB, 2);
      SeqA = _mm_insert_epi16(SeqA,*(a-(index)),0);
      SeqB = _mm_insert_epi16(SeqB,*(b-(index)),0);

      index++;

      tmp = _mm_shufflelo_epi16(SeqB,27);
      tmp = _mm_slli_si128(tmp, 2);
      tmpValue = _mm_extract_epi16(tmp, 5);
      tmp = _mm_insert_epi16(tmp, tmpValue, 0);

      Result = _mm_cmpeq_epi16(SeqA, tmp);
      Diag = _mm_andnot_si128(Result, MASK);

      R0 = _mm_min_epi16(R1 + Side2, R0 + Diag);
      R0 = _mm_min_epi16(R0, _mm_slli_si128(R1,2) + Down2);

      if (_mm_extract_epi16(R0, 0) > errThreshold
	  && _mm_extract_epi16(R0, 1) > errThreshold
	  && _mm_extract_epi16(R0, 2) > errThreshold
	  && _mm_extract_epi16(R0, 3) > errThreshold
	  && _mm_extract_epi16(R0, 4) > errThreshold
	  && _mm_extract_epi16(R1, 0) > errThreshold
	  && _mm_extract_epi16(R1, 1) > errThreshold
	  && _mm_extract_epi16(R1, 2) > errThreshold
	  && _mm_extract_epi16(R1, 3) > errThreshold)
	return -1;

      if (i == 2 * lenb - e) {
	tmp = _mm_srli_si128(R0,2);
	for (k = 0; k < e - 1; k++)
	  tmp = _mm_srli_si128(tmp,2);
	minError = _mm_extract_epi16(tmp,0);
      }

    }

    else {
      Result = _mm_cmpeq_epi16(SeqA, _mm_shufflelo_epi16(SeqB,27));
      Diag = _mm_andnot_si128(Result, MASK);

      R1 = _mm_min_epi16(_mm_srli_si128(R0,2) + Side1, R1 + Diag);
      R1 = _mm_min_epi16(R1, R0 + Down1);

      if (i >= 2 * lenb - e) {
	tmp = _mm_srli_si128(R1,2);
	for (k = 0; k < e - 2; k++)
	  tmp = _mm_srli_si128(tmp,2);
	minError = min(minError, _mm_extract_epi16(tmp,0));
      }
    }

  }

  j = 0;
  int tmpE = e;
  for (; j < 2 * (e - 2) + 1; j++) {

    Diag = _mm_xor_si128(Diag, Diag);
    //set the first element
    if (j == 0) {
      for (k = 0; k <= e - 1; k++) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag =
	  _mm_insert_epi16(Diag, *(b-(lenb-1-k)) != *(a-((i-lenb)-1+k)),0);
      }

      R0 = _mm_min_epi16(R1 + Side2, R0 + Diag);
      R0 = _mm_min_epi16(R0, _mm_slli_si128(R1,2) + Down2);

      tmpE--;

      tmp = _mm_srli_si128(R0,2);
      for (k = 0; k < e - 2; k++)
	tmp = _mm_srli_si128(tmp,2);
      minError = min(minError, _mm_extract_epi16(tmp,0));
    } else if (j % 2 == 0) {
      for (k = 0; k < tmpE; k++) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag =
	  _mm_insert_epi16(Diag, *(b-(lenb-1-k)) != *(a-((i-lenb)-1+k)),0);
      }

      R0 = _mm_min_epi16(R1 + Side2, R0 + Diag);
      R0 = _mm_min_epi16(R0, _mm_slli_si128(R1,2) + Down2);

      tmpE--;

      tmp = _mm_srli_si128(R0,2);
      for (k = 0; k < tmpE - 1; k++)
	tmp = _mm_srli_si128(tmp,2);
      minError = min(minError, _mm_extract_epi16(tmp,0));
    }

    else {
      for (k = 0; k < tmpE; k++) {
	Diag = _mm_slli_si128(Diag, 2);
	Diag =
	  _mm_insert_epi16(Diag, *(b-(lenb-1-k)) != *(a-((i-lenb)-1+k)),0);
      }

      R1 = _mm_min_epi16(_mm_srli_si128(R0,2) + Side1, R1 + Diag);
      R1 = _mm_min_epi16(R1, R0 + Down1);

      tmp = _mm_srli_si128(R1,2);
      for (k = 0; k < tmpE - 2; k++)
	tmp = _mm_srli_si128(tmp,2);
      minError = min(minError, _mm_extract_epi16(tmp,0));
    }
    i++;
  }
  //Diag

  Diag = _mm_xor_si128(Diag, Diag);
  Diag = _mm_insert_epi16(Diag, 2*errThreshold, 0);
  Diag = _mm_insert_epi16(Diag, *(a-(lenb+e-2)) != *(b-(lenb-1)), 1);

  Side1 = _mm_insert_epi16(Side1,1,0);
  Side1 = _mm_insert_epi16(Side1,1,1);

  Down1 = _mm_insert_epi16(Down1, 2*errThreshold, 0);
  Down1 = _mm_insert_epi16(Down1, 1, 1);

  R1 = _mm_min_epi16(R0 + Side1, _mm_slli_si128(R1,2) + Diag);
  R1 = _mm_min_epi16(R1, _mm_slli_si128(R0,2) + Down1);

  minError = min(minError, _mm_extract_epi16(R1,1));

  Diag = _mm_insert_epi16(Diag, *(a-(lenb+e-1)) != *(b-(lenb-1)), 0);
  Down1 = _mm_insert_epi16(Down1, 1, 0);

  R0 = _mm_min_epi16(R1 + Down1, R0 + Diag);
  R0 = _mm_min_epi16(R0, _mm_srli_si128(R1,2) + Side1);

  minError = min(minError, _mm_extract_epi16(R0,0));

  if (minError > mismatch)
    return -1;
  return minError;
}



/*------------------------------------------------------------------------------
 * Initializes best mapping information, note that FullMappingInfo and
 * BestFullMappingInfo only differ in (double tprob) fields. This is also called
 * in the paired-end mode.
 *----------------------------------------------------------------------------*/
void
initBestMapping ( int totalReadNumber )
{
	int i			   = 0;
	bestHitMappingInfo = getMem(totalReadNumber * sizeof(BestFullMappingInfo));
	for (i = 0; i < totalReadNumber; i++)
	{
	  	bestHitMappingInfo[i].loc	= -1;
	  	bestHitMappingInfo[i].tprob = 0.0; 
	}
}



/*------------------------------------------------------------------------------
 * This is called from baseFast.c, after all mapping has been finished. The
 * best mapping information is not immediately dumped to the disk after a read
 * has been mapped as in the version without best mapping mode. It is performed
 * after all reads have been mapped.
 *----------------------------------------------------------------------------*/
void
finalizeBestSingleMapping (   ) 
{
	int		 i		  = 0;
	char	*_tmpQual, *_tmpSeq;
	char	 rqual[SEQ_LENGTH + 1];
	rqual[SEQ_LENGTH] = '\0';

	for(i = 0; i < _msf_seqListSize; i++)
	{
		if(_msf_seqList[i].hits[0] != 0)
	  	{		
			if (bestHitMappingInfo[i].dir)
			{
				reverse(_msf_seqList[i].qual, rqual, SEQ_LENGTH);
				_tmpQual = rqual;
				_tmpSeq	 = _msf_seqList[i].rseq;
			}
			else
			{
				_tmpQual = _msf_seqList[i].qual;
				_tmpSeq	 = _msf_seqList[i].seq;
			}

			/* Fill out mapping information for output. */
			_msf_output.QNAME = _msf_seqList[i].name;
			_msf_output.FLAG  = 16 * bestHitMappingInfo[i].dir;
			_msf_output.RNAME = bestHitMappingInfo[i].chr;
			_msf_output.POS	  = bestHitMappingInfo[i].loc;

			if (seqFastq)
			  	_msf_output.MAPQ = mapQ(i);
			else
			  	_msf_output.MAPQ = 255;

			_msf_output.CIGAR  = bestHitMappingInfo[i].cigar;
			_msf_output.MRNAME = "*";
			_msf_output.MPOS   = 0;
			_msf_output.ISIZE  = 0;

			_msf_output.SEQ	 = _tmpSeq;
			_msf_output.QUAL = _tmpQual;

			_msf_output.optSize	  = 2;
			_msf_output.optFields = _msf_optionalFields;

			_msf_optionalFields[0].tag	= "NM";
			_msf_optionalFields[0].type = 'i';
			_msf_optionalFields[0].iVal = bestHitMappingInfo[i].err;

			_msf_optionalFields[1].tag	= "MD";
			_msf_optionalFields[1].type = 'Z';
			_msf_optionalFields[1].sVal = bestHitMappingInfo[i].md;

			/* Output mapping information */
			output(_msf_output);
	  	}
	}

	/* Let go */
	freeMem(bestHitMappingInfo, _msf_seqListSize * sizeof(FullMappingInfo));
}



/*------------------------------------------------------------------------------
 * Compare two given reads with respect to their hash values.
 *----------------------------------------------------------------------------*/
int
compare ( const void *a, const void *b )
{
  	return ((Pair *) a)->hv - ((Pair *) b)->hv;
}



/*------------------------------------------------------------------------------
 * Sort the reads wrt their hash values of the read sequences.
 *----------------------------------------------------------------------------*/
void
preProcessReads (   )
{
	int i = 0;

	/* First pass: compute hash values of each read sequence. */
	_msf_sort_seqList = getMem(_msf_seqListSize * sizeof(Pair));
	for (i = 0; i < _msf_seqListSize; i++)
	{
		/* Note that we compute hash value using the function in HashTable. The
		 * *hashValue and *rhashValue fields of reads are not utilized. */
	  	_msf_sort_seqList[i].hv = hashVal(_msf_seqList[i].seq);
	  	_msf_sort_seqList[i].readNumber = i;
	}

	qsort(_msf_sort_seqList, _msf_seqListSize, sizeof(Pair), compare);

	_msf_map_sort_seqList = getMem(_msf_seqListSize * sizeof(int));

	for (i = 0; i < _msf_seqListSize; i++)
	  	_msf_map_sort_seqList[_msf_seqList[i].readNumber] = i;
}



/*------------------------------------------------------------------------------
 * This is never called.
 *----------------------------------------------------------------------------*/
int
verifySingleEnd ( int index, char* seq, int offset )
{
  int curOff = 0;
  int i;

  char *ref;

  int err;
  int errCnt = 0;
  int errCntOff = 0;
  int NCntOff = 0;

  ref = _msf_refGen + index - 1;

  verificationCnt++;

  for (i = 0; i < SEQ_LENGTH; i++) {
    err = *ref != *seq;
    errCnt += err;
    if (errCnt > errThreshold) {

      return -1;
    }

    if (i >= _msf_samplingLocs[curOff]
	&& i <= _msf_samplingLocsEnds[curOff]) {
      errCntOff += err;
      NCntOff += (*seq == 'N');
    } else if (curOff < _msf_samplingLocsSize
	       && i >= _msf_samplingLocs[curOff + 1]) {

      if (errCntOff == 0 && NCntOff == 0 && offset > curOff) {
	return -1;
      }

      errCntOff = 0;
      NCntOff = 0;
      curOff++;

      if (i >= _msf_samplingLocs[curOff]) {
	errCntOff += err;
	NCntOff += (*seq == 'N');
      }
    }

    ref++;
    seq++;
  }
  return errCnt;
}



/*------------------------------------------------------------------------------
 * Initializes various parameters such as optional fields, sampling locations,
 * reference genome related stuff etc. It also separates normal mode mapping
 * from paired-end mapping as paired-end mapping has different data structure
 * requirements.
 *----------------------------------------------------------------------------*/
void
initFAST ( Read *seqList, int seqListSize, int *samplingLocs,
		   int samplingLocsSize, char *genFileName )
{
	int i;

	/* Initialize and allocate optional fields for result output. */
	if (_msf_optionalFields == NULL)
	{
		_msf_op = getMem(SEQ_LENGTH);

		/* More space for paired-end mapping. */
		if (pairedEndMode)
		{
		  	_msf_optionalFields = getMem(8 * sizeof(OPT_FIELDS));
		}
		else
		{
		  	_msf_optionalFields = getMem(2 * sizeof(OPT_FIELDS));
		}

		for (i = 0; i < 200; i++)
		{
		  	sprintf(_msf_numbers[i], "%d%c", i, '\0');
		}
		sprintf(_msf_cigar, "%dM", SEQ_LENGTH);
	}


	/* Initialize read sequence related information */
	if (_msf_samplingLocsEnds == NULL)
	{
		_msf_samplingLocs	  = samplingLocs;
		_msf_samplingLocsSize = samplingLocsSize;

		/* The endpoints of the sampling locations were not determined. */
		_msf_samplingLocsEnds = getMem(sizeof(int) * _msf_samplingLocsSize);
		for (i = 0; i < _msf_samplingLocsSize; i++)
		{
		  	_msf_samplingLocsEnds[i] = _msf_samplingLocs[i] + WINDOW_SIZE - 1;
		}

		_msf_seqList	 = seqList;
		_msf_seqListSize = seqListSize;

		/* Sort reads wrt to their hash values. */
		preProcessReads();

		/* OEA means one end anchored. These are related with paired-end mapping
		 * mode. Allocate one integer for OEA and discordant mapping for each
		 * read. */
		_msf_oeaMapping = getMem(_msf_seqListSize * sizeof(int));
		for (i = 0; i < _msf_seqListSize; i++)
		{
		  	_msf_oeaMapping[i] = 0;
		}

		_msf_discordantMapping = getMem(_msf_seqListSize * sizeof(int));
		for (i = 0; i < _msf_seqListSize; i++)
		{
		  	_msf_discordantMapping[i] = 0;
		}
	}

	/* Initialize reference genome related information. */
	if (_msf_refGenName == NULL)
	{
	  	_msf_refGenName = getMem(4 * SEQ_LENGTH);
	}
	
	_msf_refGen		  = getRefGenome();
	_msf_refGenLength = strlen(_msf_refGen);

	_msf_refGenOffset = getRefGenomeOffset();
	snprintf(_msf_refGenName, 4 * SEQ_LENGTH, "%s%c", getRefGenomeName(), '\0');
	_msf_refGenName[strlen(getRefGenomeName())] = '\0';

	/* Doing something about the verified locations. There is an entry for each
	 * reference genome location. */
	if (_msf_verifiedLocs != NULL)
	  	freeMem(_msf_verifiedLocs, sizeof(int) * (_msf_refGenLength + 1));

	_msf_verifiedLocs = (int *) getMem(sizeof(int) * (_msf_refGenLength + 1));

	/* Inits a value for each position of the reference genome. The initialized
	 * value of it is still a mystery. Nonetheless, the verified location
	 * entries are assigned readIds, so it makes sense he may have wanted to
	 * make it something different than range of readId values it can get. */
	for (i = 0; i <= _msf_refGenLength; i++)
	  	_msf_verifiedLocs[i] = _msf_seqListSize * 10 + 1;

	/* Paired-end mode related auxiliary data. */
	if (pairedEndMode && _msf_seqHits == NULL)
	{
		/* For each read, allocate a MappingInfo structure which stores the
		 * output information. It is actually a linked list containing
		 * information about each mapped location for the read. */
		_msf_mappingInfo = getMem(seqListSize * sizeof(MappingInfo));

		for (i = 0; i < seqListSize; i++)
		{
			_msf_mappingInfo[i].next = NULL;
			_msf_mappingInfo[i].size = 0;
		}

		/* For each read, number of read sequence hits? */
		_msf_seqHits = getMem((_msf_seqListSize) * sizeof(int));

		for (i = 0; i < _msf_seqListSize; i++)
		{
			_msf_seqHits[i] = 0;
		}

		/* The read has concordant mapping or no. Dont get why its size is equal
		 * to sequence list size / 2. */
		_msf_readHasConcordantMapping = getMem(
						   _msf_seqListSize / 2 * sizeof(char));
		for (i = 0; i < _msf_seqListSize / 2; i++)
		{
			_msf_readHasConcordantMapping[i] = 0;
		}

		/* Allocates space for reference genome. */
		initLoadingRefGenome(genFileName);
	}

	/* Can read more of reference genome later. */
	if (_msf_refGenOffset == 0)
	  	_msf_refGenBeg = 1;
	else
	  	_msf_refGenBeg = CONTIG_OVERLAP - SEQ_LENGTH + 2 + errThreshold;
	_msf_refGenEnd = _msf_refGenLength - SEQ_LENGTH + 1;
}



/*------------------------------------------------------------------------------
 * 	Freed data:
 *		_msf_seqHits
 *		_msf_refGenName
 *		_msf_map_sort_seqList
 *		_msf_sort_seqList
 *----------------------------------------------------------------------------*/
void
finalizeFAST (   )
{
	freeMem(_msf_seqHits, (_msf_seqListSize) * sizeof(int));
	freeMem(_msf_refGenName, 4 * SEQ_LENGTH);

	freeMem(_msf_map_sort_seqList, sizeof(Pair) * _msf_seqListSize);
	freeMem(_msf_sort_seqList, sizeof(int) * _msf_seqListSize);
}



/*------------------------------------------------------------------------------
 * This does not seem to be in the anti-diagonal order, it seems to proceed in
 * CBC (column-by-column) order. This is not used anywhere.

 * OLD COMMENTS (not @OGUZ-COMMENT)
 * Will apply the Levenshtein Dynamic programming.  Different from
 * verifySingleEndEditDistance fucntion as in this function only one dynamic
 * table is made while in verifySingleEndEditDistance two dynamic table is made
 * for each right and left string
 *----------------------------------------------------------------------------*/
int
editDistance ( int refIndex, char *seq, int seqLength, char *matrix )
{
	int i			   = 0;
	int size		   = 0;
	int error		   = 0;
	int rIndex		   = 0;
	int directionIndex = 0;

	int min		 = 0;
	int minIndex = 0;

	int tempUp	 = 0;
	int tempDown = 0;

	char	*ref;

	int errorString = 0;
	/*
	  1: Up
	  2: Side
	  3: Diagonal Match
	  4: Diagonal Mismatch
	*/

	int upValue;
	int diagValue;
	int sideValue;

	/* ref stores the current entry in ref genome. */
	ref = _msf_refGen + refIndex - 1;

	rIndex = 1;

	for (i = 0; i <= errThreshold; i++)
	{
		score[0][i] = i;
		score[i][0] = i;
	}

	/* Runs the DP? */
	while (rIndex <= seqLength + errThreshold)
	{	
		tempUp = ((rIndex - errThreshold) > 0 ?
				  ((rIndex > seqLength) ?
				   seqLength - errThreshold :
				   rIndex - errThreshold) : 1);
		tempDown = ((rIndex >= seqLength - errThreshold) ?
					seqLength + 1 :
					rIndex + errThreshold + 1);
		
		for (i = tempUp; i < tempDown; i++)
		{
			
			errorString = (*(ref + rIndex - 1) == *(seq + i - 1));
			upValue		= score[i - 1][rIndex] + 1;
			diagValue	= score[i - 1][rIndex - 1] + !errorString;
			sideValue	= score[i][rIndex - 1] + 1;

			if (i != tempUp && i != tempDown - 1)
		  		score[i][rIndex] = min3(sideValue, diagValue , upValue);
			else if ((i == ((rIndex - errThreshold) > 0 ? rIndex - errThreshold : 1))
				 && rIndex <= seqLength)
		  		score[i][rIndex] = min(sideValue, diagValue);
			else if (rIndex > seqLength && (i == seqLength - errThreshold))
		  		score[i][rIndex] = sideValue;
			else
		  		score[i][rIndex] = min(diagValue , upValue);

			if (i == tempUp)
		  		error = score[i][rIndex];
			else if (error > score[i][rIndex])
		  		error = score[i][rIndex];
		}
		rIndex++;
	}

	min		 = score[seqLength][seqLength + errThreshold];
	minIndex = seqLength + errThreshold;

	// Find the Best error for all the possible ways.
	for (i = 1; i <= 2 * errThreshold; i++)
	{
		if (min >= score[seqLength][seqLength + errThreshold - i]
		&& seqLength + errThreshold - i > 0)
		{
			min = score[seqLength][seqLength + errThreshold - i];
			minIndex = seqLength + errThreshold - i;
		}
	}

	error = score[seqLength][minIndex];

	/* Find out how we came here (I = Insert, D = Delete, M = Match) */
	directionIndex = seqLength;
	rIndex = minIndex;
	while (directionIndex != 0 || rIndex != 0)
	{
		if (rIndex == 0)
		{
			if (score[directionIndex][rIndex] - score[directionIndex - 1][rIndex] == 1)
			{
				matrix[size] = *(seq + directionIndex - 1);
				size++;
				matrix[size] = 'I';
				directionIndex--;
			}
		}
		else if (directionIndex == 0)
		{
			if (score[directionIndex][rIndex] - score[directionIndex][rIndex - 1] == 1)
			{
		  		matrix[size] = *(ref + rIndex - 1);
		  		size++;
		  		matrix[size] = 'D';
		  		rIndex--;
			}
		}
		else if (directionIndex - rIndex == errThreshold)
		{
			if (score[directionIndex][rIndex] - score[directionIndex - 1][rIndex] == 1)
			{
				matrix[size] = *(seq + directionIndex - 1);
				size++;
				matrix[size] = 'I';
				directionIndex--;
			}
			else if (score[directionIndex][rIndex] - score[directionIndex - 1][rIndex - 1] == 1)
			{
		  		matrix[size] = *(ref + rIndex - 1);
		  		rIndex--;
		  		directionIndex--;
			}
			else
			{
		  		matrix[size] = 'M';
		  		rIndex--;
		  		directionIndex--;
			}
		}
		else if (rIndex - directionIndex == errThreshold)
		{
			if (score[directionIndex][rIndex] - score[directionIndex][rIndex - 1] == 1)
			{
		  		matrix[size] = *(ref + rIndex - 1);
		  		size++;
		  		matrix[size] = 'D';
		  		rIndex--;
			}
			else if (score[directionIndex][rIndex] - score[directionIndex - 1][rIndex - 1] == 1)
			{
		  		matrix[size] = *(ref + rIndex - 1);
		  		rIndex--;
		  		directionIndex--;
			}
			else
			{
		  		matrix[size] = 'M';
		  		rIndex--;
		  		directionIndex--;
			}
		}
		else
		{
			if (score[directionIndex][rIndex] - score[directionIndex - 1][rIndex] == 1
				&& directionIndex != 0)
			{
				matrix[size] = *(seq + directionIndex - 1);
				size++;
				matrix[size] = 'I';
				directionIndex--;
			}
			else if (score[directionIndex][rIndex] - score[directionIndex][rIndex - 1] == 1 && rIndex != 0)
			{
				matrix[size] = *(ref + rIndex - 1);
				size++;
				matrix[size] = 'D';
				rIndex--;
			}
			else if (score[directionIndex][rIndex] - score[directionIndex - 1][rIndex - 1] == 1)
			{
				matrix[size] = *(ref + rIndex - 1);
				rIndex--;
				directionIndex--;
			}
			else
			{
		  		matrix[size] = 'M';
		  		rIndex--;
		  		directionIndex--;
			}
		}
		size++;
	}

	matrix[size] = '\0';

	char returnString[2 * SEQ_MAX_LENGTH];

	returnString[0] = '\0';
	reverse(matrix, returnString, size);
	sprintf(matrix, "%s", returnString);

	return error;
}



/*------------------------------------------------------------------------------
 * PENDING
 *----------------------------------------------------------------------------*/
/*
  Will apply the Levenshtein Dynamic programming.
  in both right and left direction as long as the
  threshould error is reached or end of string length
*/

int
msfHashVal ( char *seq )
{
  int i = 0;
  int val = 0, numericVal = 0;

  while (i < 6) {
    switch (seq[i]) {
    case 'A':
      numericVal = 0;
      break;
    case 'C':
      numericVal = 1;
      break;
    case 'G':
      numericVal = 2;
      break;
    case 'T':
      numericVal = 3;
      break;
    default:
      return -1;
      break;
    }
    val = (val << 2) | numericVal;
    i++;
  }
  return val;
}



/*------------------------------------------------------------------------------
 * -FULLY COMMENTED-
 *	Called when error threshold = 2
 *
 *              |--> This is where matched against ref genome
 * |--> lSeq    |    (refIndex + seqLength should have the same content)        
 * |            |            |--> rSeq
 * |            |            |
 * TGCTCCGCCCGAAAGCCTGGATCCTCAGGGCCCTCCCGCCTGGCCCTCAG
 * -------------+++++++++++++------------------------
 *   lSeqLength    matched         rSeqLength
 *	              segLength
 *
 *	RETURN
 *	matrix			: 1D array (for storing I = Insert, D = Delete, M = Match)
 *	map_location	:
 *	seqHashValues	: Not used
 *----------------------------------------------------------------------------*/
int
verifySingleEndEditDistance2 ( int refIndex, char *lSeq, int lSeqLength,
							   char *rSeq, int rSeqLength, int segLength,
							   char *matrix, int *map_location,
							   short *seqHashValue )
{
	int		 i = 0;
	char	*ref;
	char	*tempref;

	int rIndex = 0;				//reference Index

	int e		   = errThreshold;
	int error	   = 0;
	int error1	   = 0;
	int totalError = 0;

	/* Using tons of global variables and not defining const static global
	 * variables for Up/Side/Diagonal Match/Diagonal Mismatch... */
	/*
	  1: Up
	  2: Side
	  3: Diagonal Match
	  4: Diagonal Mismatch
	*/
	int minIndex1	   = 0;
	int minIndex2	   = 0;
	int directionIndex = 0;
	int size		   = 0;
	int startIndex1	   = 0;

	rIndex = 1;

	char matrixR[SEQ_MAX_LENGTH];
	char matrixL[SEQ_MAX_LENGTH];

	ref		= _msf_refGen + refIndex - 1;
	tempref = _msf_refGen + refIndex - 1;

	/* This is not used */
	int jumpIndex = 0;

	/* A return value of -1 means it did not match. If neither of them does not
	 * match, we return immediately. */
	if (rSeqLength != 0)		/* right seq */
	{	
		error1 = forwardEditDistance2SSE2(ref + segLength + jumpIndex,
						  rSeqLength - jumpIndex, rSeq + jumpIndex,
						  rSeqLength - jumpIndex);
		if (error1 == -1)
		{
		  	return -1;
		}
	}

	if (lSeqLength != 0)		/* left seq */
	{		
		error = backwardEditDistance2SSE2(ref - 1, lSeqLength,
						  lSeq + lSeqLength - 1, lSeqLength);
		if (error == -1)
		{
		  	return -1;
		}
	}

	/* We matched the read sequence to the reference genome */

	matrixL[0] = '\0';
	matrixR[0] = '\0';

	ref = _msf_refGen + refIndex - 1;

	/* If the edit distance of left seq and right seq summed is greater than the
	 * threshold, could not match */
	if (error1 + error > errThreshold)
	{
	  	return -1;
	}

	ref	   = _msf_refGen + refIndex - 1;
	rIndex = startIndex1 + 1;
	
	int i0 = 0;
	int i1 = 0;
	int i2 = 0;

	/* These are actually the same variables used in the edit distance
	 * computation, R0 and R1 are responsible for storing the most recent two
	 * antidiagonals computed. */
	__m128i R0;
	__m128i R1;

	/* To handle edges of the banded matrix */
	__m128i Side1, Side2, Side; //side matrix
	__m128i Down1, Down2, Down; //down matrix

	/* Stores the newly matching information for the current antidiagonal. */
	__m128i Diag;

	__m128i tmp;


	/* You go your way my love */
	R0	  = _mm_setzero_si128();
	R1	  = _mm_setzero_si128();
	Diag  = _mm_setzero_si128();
	Side1 = _mm_setzero_si128();
	Side2 = _mm_setzero_si128();
	Down1 = _mm_setzero_si128();
	Down2 = _mm_setzero_si128();
	Down  = _mm_setzero_si128();
	Side  = _mm_setzero_si128();
	tmp	  = _mm_setzero_si128();
	
	/* This array stores currently compared ref and seq match/mismatch
	 * status. Note that if the current iteration is ODD, we compare two chars
	 * and if it is EVEN, we compare three chars.  */
	int mismatch[3] = { 0, 0, 0 };

	
	/* Trace back left sequence. Left segment length should be greater than 0,
	 * so we know some edit distance computation has already been
	 * performed. Computes score and direction1 matrices, these are used later
	 * in computing I-D-M? */
	if (lSeqLength != 0)
	{		
		/* a and b point to reference genome and left sequence */
		char	*a;
		char	*b;
		a = ref - 1;
		b = lSeq + lSeqLength - 1;

		/* Init R0 and R1 */
		R0 = _mm_insert_epi16(R0,0,0);
		R1 = _mm_insert_epi16(R1,1,0);
		R1 = _mm_insert_epi16(R1,1,1);

		/* score and direction matrices */
		score[0][0]		 = 0;
		score[1][0]		 = 1;
		direction1[1][0] = 1;
		score[0][1]		 = 1;
		direction1[0][1] = 2;

		mismatch[0] = ((a[0]) != (b[0]));

		/* Initialize Diag, Side1, Side2, Down1 and Down2 exactly in the same
		 * way done as in edit distance computation */
		Diag  = _mm_insert_epi16(Diag,2*e,0);
		Diag  = _mm_insert_epi16(Diag,mismatch[0],1);
		Diag  = _mm_insert_epi16(Diag,2*e,2);
		Side1 = _mm_insert_epi16(Side1,1,0);
		Side1 = _mm_insert_epi16(Side1,1,1);
		Side1 = _mm_insert_epi16(Side1,2*e,2);
		Side2 = _mm_insert_epi16(Side2,2*e,0);
		Side2 = _mm_insert_epi16(Side2,1,1);
		Side2 = _mm_insert_epi16(Side2,1,2);
		Down1 = _mm_insert_epi16(Down1,1,0);
		Down1 = _mm_insert_epi16(Down1,1,1);
		Down1 = _mm_insert_epi16(Down1,2*e,2);
		Down2 = _mm_insert_epi16(Down2,2*e,0);
		Down2 = _mm_insert_epi16(Down2,1,1);
		Down2 = _mm_insert_epi16(Down2,1,2);

		/* Compute first value of R0 manually, R0 stores the even indexed
		 * iterations so it has three score values in it. */
		tmp = _mm_slli_si128(R1,2);
		R0	= _mm_min_epi16(R1 + Side1, R0 + Diag);
		R0	= _mm_min_epi16(R0, tmp + Down2);

		/* Get three values from R0 and store them in score matrix and update
		 * direction matrix as well */
		i0				 = _mm_extract_epi16(R0, 0);
		i1				 = _mm_extract_epi16(R0, 1);
		i2				 = _mm_extract_epi16(R0, 2);
		score[0][2]		 = i0;
		score[1][1]		 = i1;
		score[2][0]		 = i2;
		direction1[0][2] = 2;
		direction1[1][1] = ((mismatch[0] == 0) ? 3 : 4);
		direction1[2][0] = 1;

		/* Start looping */
		for (i = 3; i < 2 * lSeqLength; i++)
		{
			if (i % 2 == 1)		/* ODD */
			{
				/* Compute match/mismatch information, note that we have only
				 * two chars to compare in ODD iterations */
				Diag = _mm_sub_epi8(Diag, Diag);
				mismatch[0] = (*(a - ((i + 1) / 2 - 1)) != *(b - ((i - 1) / 2 - 1)));
				Diag = _mm_insert_epi16(Diag,mismatch[0],0);
				mismatch[1] = (*(a - ((i - 1) / 2 - 1)) != *(b - ((i + 1) / 2 - 1)));
				Diag = _mm_insert_epi16(Diag,mismatch[1],1);

				/* Find antidiagonal min values and store in R1 */
				tmp = _mm_srli_si128(R0,2);
				R1 = _mm_min_epi16(tmp + Side1, R1 + Diag);
				R1 = _mm_min_epi16(R1, R0 + Down1);

				/* Update score table (two values) */
				i0 = _mm_extract_epi16(R1, 0);
				i1 = _mm_extract_epi16(R1, 1);
				score[i / 2][i / 2 + 1] = i0;
				score[i / 2 + 1][i / 2] = i1;

				/* Update direction table (again, two values) */
				direction1[i / 2][i / 2 + 1] =
				  (score[i / 2][i / 2 + 1] == score[i / 2 - 1][i / 2]
				   && mismatch[0] == 0) ? 3 :
				  (score[i / 2][i / 2 + 1] - score[i / 2 - 1][i / 2 + 1]
				   == 1) ? 1 :
				  (score[i / 2][i / 2 + 1] - score[i / 2][i / 2] == 1) ?
				  2 : 4;

				direction1[i / 2 + 1][i / 2] =
				  (score[i / 2 + 1][i / 2] == score[i / 2][i / 2 - 1]
				   && mismatch[1] == 0) ? 3 :
				  (score[i / 2 + 1][i / 2] - score[i / 2][i / 2] == 1) ?
				  1 :
				  (score[i / 2 + 1][i / 2] - score[i / 2 + 1][i / 2 - 1]
				   == 1) ? 2 : 4;

				/* Update error if we came up to the end */
				if (i > 2 * lSeqLength - 2)
				{
					error = min(error, i1);
					if (error == i1)
					  minIndex1 = i - lSeqLength;
				}
			}
			else if (i % 2 == 0) /* EVEN */
			{
				/* Compute match/mismatch information, note that we have three
				 * chars to compare in EVEN iterations */
				mismatch[0] = (*(a - (i / 2)) != *(b - (i / 2 - 2)));
				Diag = _mm_insert_epi16(Diag,mismatch[0],0);
				mismatch[1] = (*(a - (i / 2 - 1)) != *(b - (i / 2 - 1)));
				Diag = _mm_insert_epi16(Diag,mismatch[1],1);
				mismatch[2] = (*(a - (i / 2 - 2)) != *(b - (i / 2)));
				Diag = _mm_insert_epi16(Diag,mismatch[2],2);

				/* Find antidiagonal min values and store in R0 */
				tmp = _mm_slli_si128(R1,2);
				R0 = _mm_min_epi16(R1 + Side1, R0 + Diag);
				R0 = _mm_min_epi16(R0, tmp + Down2);

				/* Update score table (three values) */
				i0							= _mm_extract_epi16(R0, 0);
				i1							= _mm_extract_epi16(R0, 1);
				i2							= _mm_extract_epi16(R0, 2);
				score[i / 2 - 1][i / 2 + 1] = i0;
				score[i / 2][i / 2]			= i1;
				score[i / 2 + 1][i / 2 - 1] = i2;

				/* Update direction table (again, three values) */
				direction1[i / 2 - 1][i / 2 + 1] =
				  (score[i / 2 - 1][i / 2 + 1] == score[i / 2 - 2][i / 2]
				   && mismatch[0] == 0) ? 3 :
				  (score[i / 2 - 1][i / 2 + 1] - score[i / 2 - 1][i / 2]
				   == 1) ? 2 : 4;

				direction1[i / 2][i / 2] =
				  (score[i / 2][i / 2] == score[i / 2 - 1][i / 2 - 1]
				   && mismatch[1] == 0) ? 3 :
				  (score[i / 2][i / 2] - score[i / 2 - 1][i / 2] == 1) ?
				  1 :
				  (score[i / 2][i / 2] - score[i / 2][i / 2 - 1] == 1) ?
				  2 : 4;

				direction1[i / 2 + 1][i / 2 - 1] =
				  (score[i / 2 + 1][i / 2 - 1] == score[i / 2][i / 2 - 2]
				   && mismatch[2] == 0) ? 3 :
				  (score[i / 2 + 1][i / 2 - 1] - score[i / 2][i / 2 - 1]
				   == 1) ? 1 : 4;

				/* Did not get it */
				if ((i / 2) % segLength == 0 && i1 == 0) /* the segment has been
														  * processed no need to
														  * process it again */
				{
					return -1;
				}

				/* Update error if we came up to the end */
				if (i == 2 * lSeqLength - 2)
				{
					error = i2;
					minIndex1 = i - lSeqLength;
				}
			}
		}


		/* Now filling remaining elements of the matrices as in the edit
		 * distance computation */
		Down1 = _mm_insert_epi16(Down1,2*e,0);

		
		/* ****************************************************************** */
		/* Filling first part of the error */
		mismatch[0] = (*(a - (i / 2)) != *(b - (i / 2 - 2)));
		Diag = _mm_insert_epi16(Diag,mismatch[0],0);
		mismatch[1] = (*(a - (i / 2 - 1)) != *(b - (i / 2 - 1)));
		Diag = _mm_insert_epi16(Diag,mismatch[1],1);
		Diag = _mm_insert_epi16(Diag,2*e,2);

		/* Store the minimmum values of the AD in R0. */
		R0 = _mm_min_epi16(R1 + Side1, R0 + Diag);
		R0 = _mm_min_epi16(R0, _mm_slli_si128(R1,2) + Down1);

		i0 = _mm_extract_epi16(R0, 0);
		i1 = _mm_extract_epi16(R0, 1);

		/* Update error */
		error = min(error, i1);
		if (error == i1)
		  	minIndex1 = i - lSeqLength;

		/* Update score and direction matrices */
		score[i / 2 - 1][i / 2 + 1] = i0;
		score[i / 2][i / 2] = i1;

		direction1[i / 2 - 1][i / 2 + 1] =
		  (score[i / 2 - 1][i / 2 + 1] == score[i / 2 - 2][i / 2]
		   && mismatch[0] == 0) ? 3 :
		  (score[i / 2 - 1][i / 2 + 1] - score[i / 2 - 1][i / 2]) ? 2 : 4;

		direction1[i / 2][i / 2] =
		  (score[i / 2][i / 2] == score[i / 2 - 1][i / 2 - 1]
		   && mismatch[1] == 0) ? 3 :
		  (score[i / 2][i / 2] - score[i / 2 - 1][i / 2] == 1) ? 1 :
		  (score[i / 2][i / 2] - score[i / 2][i / 2 - 1] == 1) ? 2 : 4;

		
		/* ****************************************************************** */
		/* Filling second part of the error */
		i++;
		Diag = _mm_sub_epi8(Diag, Diag);
		Diag = _mm_insert_epi16(Diag,2*e,0);
		mismatch[0] = (*(a - (i / 2)) != *(b - (lSeqLength - 1)));
		Diag = _mm_insert_epi16(Diag,mismatch[0],1);
		Diag = _mm_insert_epi16(Diag,2*e,2);

		/* Store the minimmum values of the AD in R1. */
		R1 = _mm_min_epi16(R0 + Side1, _mm_slli_si128(R1,2) + Diag);
		R1 = _mm_min_epi16(R1, _mm_slli_si128(R0,2) + Down1);

		i0 = _mm_extract_epi16(R1, 0);
		i1 = _mm_extract_epi16(R1, 1);

		/* Update error */
		error = min(error, i1);
		if (error == i1)
		  	minIndex1 = i - lSeqLength;

		/* Update score and direction matrices */
		score[i / 2 - 1][i / 2 + 2] = i0;
		score[i / 2][i / 2 + 1] = i1;

		direction1[i / 2 - 1][i / 2 + 2] =
		  (score[i / 2 - 1][i / 2 + 2] == score[i / 2 - 2][i / 2 + 1]
		   && mismatch[0] == 0) ? 3 :
		  (score[i / 2 - 1][i / 2 + 2] - score[i / 2 - 1][i / 2 + 1] == 1) ?
		  2 : 4;

		direction1[i / 2][i / 2 + 1] =
		  (score[i / 2][i / 2 + 1] == score[i / 2 - 1][i / 2]) ? 3 :
		  (score[i / 2][i / 2 + 1] - score[i / 2 - 1][i / 2 + 1] == 1) ?
		  1 :
		  (score[i / 2][i / 2 + 1] - score[i / 2][i / 2] == 1) ? 2 : 4;


		/* ****************************************************************** */
		/* Filling final element of the matrix */
		i++;
		Diag = _mm_sub_epi8(Diag, Diag);
		mismatch[0] = (*(a - (i / 2)) != *(b - (lSeqLength - 1)));
		Diag = _mm_insert_epi16(Diag,mismatch[0],0);

		/* Store the minimmum values of the AD in R0. */
		Down = _mm_insert_epi16(Down,1,0);
		Side = _mm_insert_epi16(Side,1,0);
		tmp = _mm_srli_si128(R1,2);
		R0 = _mm_min_epi16(R1 + Down, R0 + Diag);
		R0 = _mm_min_epi16(R0, tmp + Side);

		i0 = _mm_extract_epi16(R0, 0);

		/* Update error */
		error = min(error, i0);
		if (error == 0)			/* Why? */
		  	return -1;
		
		if (error == i0)
		  	minIndex1 = i - lSeqLength;

		/* Update direction for the final element */
		if (mismatch[0] == 0)
		  	direction1[lSeqLength][lSeqLength + errThreshold] = 3;
		else
		{
			if (score[lSeqLength][lSeqLength + errThreshold]
			- score[lSeqLength][lSeqLength + errThreshold - 1] == 1)
		  		direction1[lSeqLength][lSeqLength + errThreshold] = 2;
			else if (score[lSeqLength][lSeqLength + errThreshold]
				 - score[lSeqLength - 1][lSeqLength + errThreshold] == 1)
		  		direction1[lSeqLength][lSeqLength + errThreshold] = 1;
			else
		  		direction1[lSeqLength][lSeqLength + errThreshold] = 4;
		}
	} /* End left sequence */


	/* Update mapping location, by also taking indels into account. */
	error1		   = error;
	error		   = 0;
	directionIndex = lSeqLength;
	rIndex		   = minIndex1;
	*map_location = ((lSeqLength == 0) ? refIndex : refIndex - rIndex);

	
	/* ********************************************************************** */
	/* ********************************************************************** */
	/* ********************************************************************** */


	/* Trace back right sequence */
	ref = ref + segLength;
	if (rSeqLength <= e)		/* Skip already */
	{
		char *a;
		char *b;

		int tmp_index = 0;

		a = ref;
		b = rSeq;

		for (tmp_index = 0; tmp_index < rSeqLength; tmp_index++)
		{
		  	matrixR[tmp_index] = (a[tmp_index] == b[tmp_index]) ? 'M' : a[tmp_index];
		}
		
		matrixR[tmp_index] = '\0';
	}
	else if (rSeqLength != 0 && rSeqLength >= e)
	{		
		/* a and b point to reference genome and left sequence */
		char	*a;
		char	*b;
		a = ref;
		b = rSeq;

		/* Init R0 and R1 */
		R0 = _mm_sub_epi8(R0, R0);
		R1 = _mm_sub_epi8(R1, R1);
		R0 = _mm_insert_epi16(R0,0,0);
		R1 = _mm_insert_epi16(R1,1,0);
		R1 = _mm_insert_epi16(R1,1,1);

		/* score and direction matrices */
		score[0][0]		 = 0;
		score[1][0]		 = 1;
		direction2[1][0] = 1;
		score[0][1]		 = 1;
		direction2[0][1] = 2;

		mismatch[0] = ((a[0]) != (b[0]));

		/* Initialize Diag, Side1, Side2, Down1 and Down2 exactly in the same
		 * way done as in edit distance computation */
		Diag  = _mm_insert_epi16(Diag,2*e,0);
		Diag  = _mm_insert_epi16(Diag,mismatch[0],1);
		Diag  = _mm_insert_epi16(Diag,2*e,2);
		Side1 = _mm_insert_epi16(Side1,1,0);
		Side1 = _mm_insert_epi16(Side1,1,1);
		Side1 = _mm_insert_epi16(Side1,2*e,2);
		Side2 = _mm_insert_epi16(Side2,2*e,0);
		Side2 = _mm_insert_epi16(Side2,1,1);
		Side2 = _mm_insert_epi16(Side2,1,2);
		Down1 = _mm_insert_epi16(Down1,1,0);
		Down1 = _mm_insert_epi16(Down1,1,1);
		Down1 = _mm_insert_epi16(Down1,2*e,2);
		Down2 = _mm_insert_epi16(Down2,2*e,0);
		Down2 = _mm_insert_epi16(Down2,1,1);
		Down2 = _mm_insert_epi16(Down2,1,2);

		/* Compute first value of R0 manually, R0 stores the even indexed
		 * iterations so it has three score values in it. */
		tmp = _mm_slli_si128(R1,2);
		R0	= _mm_min_epi16(R1 + Side1, R0 + Diag);
		R0	= _mm_min_epi16(R0, tmp + Down2);

		/* Get three values from R0 and store them in score matrix and update
		 * direction matrix as well */
		i0				 = _mm_extract_epi16(R0, 0);
		i1				 = _mm_extract_epi16(R0, 1);
		i2				 = _mm_extract_epi16(R0, 2);
		score[0][2]		 = i0;
		score[1][1]		 = i1;
		score[2][0]		 = i2;
		direction2[0][2] = 2;
		direction2[1][1] = ((mismatch[0] == 0) ? 3 : 4);
		direction2[2][0] = 1;

		/* Start looping */
		for (i = 3; i < 2 * rSeqLength; i++)
		{
			if (i % 2 == 1)		/* ODD */
			{
				/* Compute match/mismatch information, note that we have only
				 * two chars to compare in ODD iterations */
				Diag = _mm_sub_epi8(Diag, Diag);
				mismatch[0] = ((a[(i + 1) / 2 - 1]) != (b[(i - 1) / 2 - 1]));
				Diag = _mm_insert_epi16(Diag,mismatch[0],0);
				mismatch[1] = ((a[(i - 1) / 2 - 1]) != (b[(i + 1) / 2 - 1]));
				Diag = _mm_insert_epi16(Diag,mismatch[1],1);

				/* Find antidiagonal min values and store in R1 */
				tmp = _mm_srli_si128(R0,2);
				R1	= _mm_min_epi16(tmp + Side1, R1 + Diag);
				R1	= _mm_min_epi16(R1, R0 + Down1);

				/* Update score table (two values) */
				i0 = _mm_extract_epi16(R1, 0);
				i1 = _mm_extract_epi16(R1, 1);
				score[i / 2][i / 2 + 1] = i0;
				score[i / 2 + 1][i / 2] = i1;

				/* Update direction table (again, two values) */
				direction2[i / 2][i / 2 + 1] =
				  (score[i / 2][i / 2 + 1] == score[i / 2 - 1][i / 2]
				   && mismatch[0] == 0) ? 3 :
				  (score[i / 2][i / 2 + 1] - score[i / 2 - 1][i / 2 + 1]
				   == 1) ? 1 :
				  (score[i / 2][i / 2 + 1] - score[i / 2][i / 2] == 1) ?
				  2 : 4;

				direction2[i / 2 + 1][i / 2] =
				  (score[i / 2 + 1][i / 2] == score[i / 2][i / 2 - 1]
				   && mismatch[1] == 0) ? 3 :
				  (score[i / 2 + 1][i / 2] - score[i / 2][i / 2] == 1) ?
				  1 :
				  (score[i / 2 + 1][i / 2] - score[i / 2 + 1][i / 2 - 1]
				   == 1) ? 2 : 4;

				/* Update error if we came up to the end */
				if (i > 2 * rSeqLength - 2)
				{
					error = min(error, i1);
					if (error == i1)
						minIndex2 = i - rSeqLength;
				}
			}
			else if (i % 2 == 0) /* EVEN */
			{
				/* Compute match/mismatch information, note that we have three
				 * chars to compare in EVEN iterations */
				mismatch[0] = ((a[i / 2]) != (b[i / 2 - 2]));
				Diag = _mm_insert_epi16(Diag,mismatch[0],0);
				mismatch[1] = ((a[i / 2 - 1]) != (b[i / 2 - 1]));
				Diag = _mm_insert_epi16(Diag,mismatch[1],1);
				mismatch[2] = ((a[i / 2 - 2]) != (b[i / 2]));
				Diag = _mm_insert_epi16(Diag,mismatch[2],2);

				/* Find antidiagonal min values and store in R0 */
				tmp = _mm_slli_si128(R1,2);
				R0 = _mm_min_epi16(R1 + Side1, R0 + Diag);
				R0 = _mm_min_epi16(R0, tmp + Down2);

				/* Update score table (three values) */
				i0							= _mm_extract_epi16(R0, 0);
				i1							= _mm_extract_epi16(R0, 1);
				i2							= _mm_extract_epi16(R0, 2);
				score[i / 2 - 1][i / 2 + 1] = i0;
				score[i / 2][i / 2]			= i1;
				score[i / 2 + 1][i / 2 - 1] = i2;

				/* Update direction table (again, three values) */
				direction2[i / 2 - 1][i / 2 + 1] =
				  (score[i / 2 - 1][i / 2 + 1] == score[i / 2 - 2][i / 2]
				   && mismatch[0] == 0) ? 3 :
				  (score[i / 2 - 1][i / 2 + 1] - score[i / 2 - 1][i / 2]
				   == 1) ? 2 : 4;

				direction2[i / 2][i / 2] =
				  (score[i / 2][i / 2] == score[i / 2 - 1][i / 2 - 1]
				   && mismatch[1] == 0) ? 3 :
				  (score[i / 2][i / 2] - score[i / 2 - 1][i / 2] == 1) ?
				  1 :
				  (score[i / 2][i / 2] - score[i / 2][i / 2 - 1] == 1) ?
				  2 : 4;

				direction2[i / 2 + 1][i / 2 - 1] =
				  (score[i / 2 + 1][i / 2 - 1] == score[i / 2][i / 2 - 2]
				   && mismatch[2] == 0) ? 3 :
				  (score[i / 2 + 1][i / 2 - 1] - score[i / 2][i / 2 - 1]
				   == 1) ? 1 : 4;
				
				if (i == 2 * rSeqLength - 2)
				{
				  	error = i2;
				  	minIndex2 = i - rSeqLength;
				}
			}
		}

		/* Now filling remaining elements of the matrices as in the edit
		 * distance computation */
		Down1 = _mm_insert_epi16(Down1,2*e,0);

		/* ****************************************************************** */
		/* Filling first part of the error */
		mismatch[0] = ((a[i / 2]) != (b[i / 2 - 2]));
		Diag = _mm_insert_epi16(Diag,mismatch[0],0);
		mismatch[1] = ((a[i / 2 - 1]) != (b[i / 2 - 1]));
		Diag = _mm_insert_epi16(Diag,mismatch[1],1);
		Diag = _mm_insert_epi16(Diag,2*e,2);

		/* Store the minimmum values of the AD in R0. */
		R0 = _mm_min_epi16(R1 + Side1, R0 + Diag);
		R0 = _mm_min_epi16(R0, _mm_slli_si128(R1,2) + Down1);

		i0 = _mm_extract_epi16(R0, 0);
		i1 = _mm_extract_epi16(R0, 1);

		/* Update error */
		error = min(error, i1);
		if (error == i1)
		  	minIndex2 = i - rSeqLength;

		/* Update score and direction matrices */
		score[i / 2 - 1][i / 2 + 1] = i0;
		score[i / 2][i / 2] = i1;

		direction2[i / 2 - 1][i / 2 + 1] =
		  (score[i / 2 - 1][i / 2 + 1] == score[i / 2 - 2][i / 2]
		   && mismatch[0] == 0) ? 3 :
		  (score[i / 2 - 1][i / 2 + 1] - score[i / 2 - 1][i / 2] == 1) ?
		  2 : 4;

		direction2[i / 2][i / 2] =
		  (score[i / 2][i / 2] == score[i / 2 - 1][i / 2 - 1]
		   && mismatch[1] == 0) ? 3 :
		  (score[i / 2][i / 2] - score[i / 2 - 1][i / 2] == 1) ? 1 :
		  (score[i / 2][i / 2] - score[i / 2][i / 2 - 1] == 1) ? 2 : 4;


		/* ****************************************************************** */
		/* Filling second part of the error */
		i++;
		Diag = _mm_sub_epi8(Diag, Diag);
		Diag = _mm_insert_epi16(Diag,2*e,0);
		mismatch[0] = ((a[i / 2]) != (b[rSeqLength - 1]));
		Diag = _mm_insert_epi16(Diag,mismatch[0],1);
		Diag = _mm_insert_epi16(Diag,2*e,2);

		/* Store the minimmum values of the AD in R1. */
		R1 = _mm_min_epi16(R0 + Side1, _mm_slli_si128(R1,2) + Diag);
		R1 = _mm_min_epi16(R1, _mm_slli_si128(R0,2) + Down1);

		i0 = _mm_extract_epi16(R1, 0);
		i1 = _mm_extract_epi16(R1, 1);

		/* Update error */
		error = min(error, i1);
		if (error == i1)
		  	minIndex2 = i - rSeqLength;

		/* Update score and direction matrices */
		score[i / 2 - 1][i / 2 + 2] = i0;
		score[i / 2][i / 2 + 1] = i1;

		direction2[i / 2 - 1][i / 2 + 2] =
		  (score[i / 2 - 1][i / 2 + 2] == score[i / 2 - 2][i / 2 + 1]
		   && mismatch[0] == 0) ? 3 :
		  (score[i / 2 - 1][i / 2 + 2] - score[i / 2 - 1][i / 2 + 1] == 1) ?
		  2 : 3;

		direction2[i / 2][i / 2 + 1] =
		  (score[i / 2][i / 2 + 1] == score[i / 2 - 1][i / 2]
		   && mismatch[0] == 0) ? 3 :
		  (score[i / 2][i / 2 + 1] - score[i / 2 - 1][i / 2 + 1] == 1) ?
		  1 :
		  (score[i / 2][i / 2 + 1] - score[i / 2][i / 2] == 1) ? 2 : 4;


		/* ****************************************************************** */
		/* Filling final element of the matrix */
		i++;
		Diag = _mm_sub_epi8(Diag, Diag);
		mismatch[0] = ((a[i / 2]) != (b[rSeqLength - 1]));
		Diag = _mm_insert_epi16(Diag,mismatch[0],0);

		/* Store the minimmum values of the AD in R0. */
		Down = _mm_sub_epi8(Down, Down);
		Down = _mm_insert_epi16(Down,1,0);
		Side = _mm_sub_epi8(Side, Side);
		Side = _mm_insert_epi16(Side,1,0);
		tmp	 = _mm_srli_si128(R1,2);
		R0	 = _mm_min_epi16(R1 + Down, R0 + Diag);
		R0	 = _mm_min_epi16(R0, tmp + Side);

		i0 = _mm_extract_epi16(R0, 0);

		/* Update error */
		error = min(error, i0);
		if (error == i0)
		  	minIndex2 = i - rSeqLength;

		/* Update direction for the final element */
		if (mismatch[0] == 0)
		  	direction2[rSeqLength][rSeqLength + errThreshold] = 3;
		else
		{
			if (score[rSeqLength][rSeqLength + errThreshold]
			- score[rSeqLength][rSeqLength + errThreshold - 1] == 1)
		  		direction2[lSeqLength][lSeqLength + errThreshold] = 2;
			else if (score[rSeqLength][rSeqLength + errThreshold]
				 - score[rSeqLength - 1][rSeqLength + errThreshold] == 1)
		  		direction2[rSeqLength][rSeqLength + errThreshold] = 1;
			else
		  		direction2[rSeqLength][rSeqLength + errThreshold] = 4;
		}
	}

	
	totalError	   = error1 + error;
	size		   = 0;
	directionIndex = rSeqLength;
	rIndex		   = minIndex2;

	

	/* Update matrixR (right sequence), which contains mapping information. M =
	 * Match, D = Delete, I = Insert. */
	if (rSeqLength > e)
	{	
		while (directionIndex != 0 || rIndex != 0)
		{

			if (direction2[directionIndex][rIndex] == 3)
			{
		  		matrixR[size] = 'M';
		  		rIndex--;
		  		directionIndex--;
			}
			else if (direction2[directionIndex][rIndex] == 4)
			{
		  		matrixR[size] = *(ref + rIndex - 1);
		  		rIndex--;
		  		directionIndex--;
			}
			else if (direction2[directionIndex][rIndex] == 2)
			{
		  		matrixR[size] = *(ref + rIndex - 1);
		  		size++;
		  		matrixR[size] = 'D';
		  		rIndex--;
			} else
			{
		  		matrixR[size] = *(rSeq + directionIndex - 1);
		  		size++;
		  		matrixR[size] = 'I';
		  		directionIndex--;
			}
			size++;
		}
		matrixR[size] = '\0';
	}
	size = 0;
	directionIndex = lSeqLength;
	rIndex = minIndex1;

	/* Update matrixL (left sequence), which contains mapping information. M =
	 * Match, D = Delete, I = Insert. */
	while (directionIndex != 0 || rIndex != 0)
	{	
		if (direction1[directionIndex][rIndex] == 3)
		{
			matrixL[size] = 'M';
			rIndex--;
			directionIndex--;
		}
		else if (direction1[directionIndex][rIndex] == 4)
		{
			matrixL[size] = *(tempref - rIndex);
			rIndex--;
			directionIndex--;
		}
		else if (direction1[directionIndex][rIndex] == 2)
		{
			matrixL[size] = 'D';
			size++;
			matrixL[size] = *(tempref - rIndex);
			rIndex--;
		}
		else
		{
		  	matrixL[size] = 'I';
		  	size++;
		  	matrixL[size] = *(lSeq + lSeqLength - directionIndex);
		  	directionIndex--;
		}

		size++;
	}

	
	matrixL[size] = '\0';

	char middle[SEQ_MAX_LENGTH];
	middle[0] = '\0';

	/* Middle segment already matched to the reference genome, so make it all
	 * M */
	for (i = 0; i < segLength; i++)
	  	middle[i] = 'M';
	middle[segLength] = '\0';

	/* Reverse right sequence values. */
	char rmatrixR[SEQ_MAX_LENGTH];
	reverse(matrixR, rmatrixR, strlen(matrixR));
	
	/* Store all of them in matrix and return the totalError */
	sprintf(matrix, "%s%s%s", matrixL, middle, rmatrixR);

	return totalError;
}



/*------------------------------------------------------------------------------
 * -FULLY COMMENTED-
 *	Called when error threshold = 3
 *
 *              |--> This is where matched against ref genome
 * |--> lSeq    |    (refIndex + seqLength should have the same content)        
 * |            |            |--> rSeq
 * |            |            |
 * TGCTCCGCCCGAAAGCCTGGATCCTCAGGGCCCTCCCGCCTGGCCCTCAG
 * -------------+++++++++++++------------------------
 *   lSeqLength    matched         rSeqLength
 *	              segLength
 *
 *	RETURN
 *	matrix			: 1D array (for storing I = Insert, D = Delete, M = Match)
 *	map_location	:
 *	seqHashValues	: Not used

 * Note:
 *	As opposed to verifySingleEndEditDistance2, this does not use SSE
 *	instructions after calling forward and backward edit distance computation.

 * 
 *----------------------------------------------------------------------------*/
int
verifySingleEndEditDistance( int refIndex, char *lSeq, int lSeqLength,
							 char *rSeq, int rSeqLength, int segLength,
							 char *matrix, int *map_location,
							 short *seqHashValue )
{
	int		 i			 = 0;
	char	*ref;
	char	*tempref;
	int		 rIndex		 = 0;	//reference Index
	int		 error		 = 0;
	int		 error1		 = 0;
	int		 error2		 = 0;
	int		 error3		 = 0;
	int		 totalError	 = 0;
	int		 ERROR_BOUND = errThreshold;
	//	int errorSegment = 0;	

	/*
	  1: Up
	  2: Side
	  3: Diagonal Match
	  4: Diagonal Mismatch
	*/

	int min			   = 0;
	int minIndex1	   = 0;
	int minIndex2	   = 0;
	int directionIndex = 0;
	int size		   = 0;

	ref		= _msf_refGen + refIndex - 1;
	tempref = _msf_refGen + refIndex - 1;

	/* Match right segment */
	if (rSeqLength != 0)
	{
		if (errThreshold % 2 == 1)
			error2 = forwardEditDistanceSSE2Odd(ref + segLength, rSeqLength,
												rSeq, rSeqLength);
		else					/* This should not be called, error threshold is
								 * 3. */
			error2 = forwardEditDistanceSSE2G(ref + segLength, rSeqLength,
											  rSeq, rSeqLength);
		
		if (error2 == -1)
		{
		  	return -1;
		}
	}

	/* Match left segment */
	if (lSeqLength != 0)
	{
		if (errThreshold % 2 == 1)
			error3 = backwardEditDistanceSSE2Odd(ref - 1, lSeqLength,
												 lSeq + lSeqLength - 1,
												 lSeqLength);
		else					/* This should not be called, error threshold is
								 * 3. */
			error3 = backwardEditDistanceSSE2G(ref - 1, lSeqLength,
											   lSeq + lSeqLength - 1,
											   lSeqLength);
		if (error3 == -1)
		{
		  	return -1;
		}
	}

	/* If the edit distance of left seq and right seq summed is greater than the
	 * threshold, could not match */
	if (error3 + error2 > errThreshold)
	{
	  	return -1;
	}
	

	/* Initialize score matrix for the left segment (scoreB) */
	for (i = 0; i < errThreshold + 1; i++)
	{
	  	scoreB[0][i] = i;
	  	scoreB[i][0] = i;
	}

	
	//int prevError = 0;
	rIndex			= 1;
	int tempUp		= 0;
	int tempDown	= 0;
	int errorString = 0;
	int upValue;
	int diagValue;
	int sideValue;

	
	/* Process the left segment to find score and match, note that this does not
	 * use SSE instructions. Simple and plain DP. rIndex is the anti-diagonal
	 * index. */
	while (rIndex <= lSeqLength + errThreshold && lSeqLength != 0)
	{
		/* Find starting and finishing indexes of this anti-diagonal. */
		tempUp = ((rIndex - ERROR_BOUND) > 0 ?
				  ((rIndex > lSeqLength) ?
				   lSeqLength - ERROR_BOUND : rIndex - ERROR_BOUND) : 1);
		tempDown = ((rIndex >= lSeqLength - ERROR_BOUND) ?
					lSeqLength + 1 : rIndex + ERROR_BOUND + 1);

		/* Compute anti-diagonal entries. */
		for (i = tempUp; i < tempDown; i++)
		{
			/* Find the value of the newly computed entry by choosing min of
			 * three values. */
			errorString = (*(ref - rIndex) == *(lSeq + lSeqLength - i));
			upValue		= scoreB[i - 1][rIndex] + 1;
			diagValue	= scoreB[i - 1][rIndex - 1] + !errorString;
			sideValue	= scoreB[i][rIndex - 1] + 1;

			/* Update score matrix */
			if (i != tempUp && i != tempDown - 1)
		  		scoreB[i][rIndex] = min3(sideValue, diagValue , upValue);

			else if ((i == ((rIndex - ERROR_BOUND) > 0 ?
							rIndex - ERROR_BOUND : 1)) &&
					 rIndex <= lSeqLength)
		  		scoreB[i][rIndex] = min(sideValue, diagValue);
			else if (rIndex > lSeqLength && (i == lSeqLength - ERROR_BOUND))
		  		scoreB[i][rIndex] = sideValue;
			else
		  		scoreB[i][rIndex] = min(diagValue , upValue);

			if (i == tempUp)
		  		error = scoreB[i][rIndex];
			else if (error > scoreB[i][rIndex])
		  		error = scoreB[i][rIndex];
		}

		/* I wouldnt leave that here. */
		if (rIndex <= lSeqLength)
		{
		  	//errorSegment = error-prevError;
		}
		rIndex++;
	}
	
	/* Find the final error value. */
	if (lSeqLength != 0)
	{
		min = scoreB[lSeqLength][lSeqLength + errThreshold];
		minIndex1 = lSeqLength + errThreshold;

		// Find the Best error for all the possible ways.
		for (i = 1; i <= 2 * errThreshold; i++)
		{
			if (min >= scoreB[lSeqLength][lSeqLength + errThreshold - i] && lSeqLength + errThreshold - i > 0)
			{
		  		min = scoreB[lSeqLength][lSeqLength + errThreshold - i];
		  		minIndex1 = lSeqLength + errThreshold - i;
			}
		}
		
		error = scoreB[lSeqLength][minIndex1];
	}

	/* Update mapping location, by also taking indels into account. */
	error1		   = error;
	error		   = 0;
	//errorSegment = 0;
	directionIndex = lSeqLength;
	rIndex		   = minIndex1;
	*map_location  = ((lSeqLength == 0) ? refIndex : refIndex - rIndex);

	
	/* ********************************************************************** */
	/* ********************************************************************** */
	/* ********************************************************************** */

	
	/* Trace back right sequence (scoreF). */
	ref = ref + segLength;
	if (rSeqLength != 0)
	{	
		for (i = 0; i < errThreshold + 1; i++)
		{
		  	scoreF[0][i] = i;
		  	scoreF[i][0] = i;
		}

		/* rIndex is the anti-diagonal index. */
		rIndex = 1;
		while (rIndex <= rSeqLength + errThreshold - error1)
		{
			/* Find starting and finishing indexes of this anti-diagonal. */
			tempUp = (rIndex - ERROR_BOUND) > 0 ?
				((rIndex > rSeqLength) ?
				 rSeqLength - ERROR_BOUND : rIndex - ERROR_BOUND) : 1;
			tempDown = ((rIndex >= rSeqLength - ERROR_BOUND) ?
						rSeqLength + 1 : rIndex + ERROR_BOUND + 1);

			/* Compute anti-diagonal entries. */
			for (i = tempUp; i < tempDown; i++)
			{
				/* Find the value of the newly computed entry by choosing min of
				 * three values. */
		  		errorString = (*(ref + rIndex - 1) == *(rSeq + i - 1));
		  		upValue		= scoreF[i - 1][rIndex] + 1;
		  		diagValue	= scoreF[i - 1][rIndex - 1] + !errorString;
		  		sideValue	= scoreF[i][rIndex - 1] + 1;

				/* Update score matrix */
		  		if (i != tempUp && i != tempDown - 1)
					scoreF[i][rIndex] = min3(sideValue, diagValue , upValue);
		  		else if ((i	== ((rIndex - ERROR_BOUND) > 0 ? rIndex - ERROR_BOUND : 1)) && rIndex <= rSeqLength)
					scoreF[i][rIndex] = min(sideValue, diagValue);
		  		else if (rIndex > rSeqLength && (i == rSeqLength - ERROR_BOUND))
					scoreF[i][rIndex] = sideValue;
		  		else
					scoreF[i][rIndex] = min(diagValue , upValue);

		  		if (i == tempUp)
					error = scoreF[i][rIndex];
		  		if (error > scoreF[i][rIndex])
					error = scoreF[i][rIndex];
			}

			/* I wouldnt leave that here. */
			if (rIndex <= rSeqLength)
			{
		  		//errorSegment = error;
			}
			rIndex++;
		}

		/* Find the final error value. */
		min		  = scoreF[rSeqLength][rSeqLength + errThreshold - error1];
		minIndex2 = rSeqLength + errThreshold - error1;

		// Find the Best error for all the possible ways.
		for (i = 1; i <= 2 * (errThreshold - error1); i++)
		{
			if (min > scoreF[rSeqLength][rSeqLength + errThreshold - error1 - i]
				&& rSeqLength + errThreshold - error1 - i > 0)
			{
		  		min = scoreF[rSeqLength][rSeqLength + errThreshold - error1 - i];
		  		minIndex2 = rSeqLength + errThreshold - error1 - i;
			}
		}
		error = scoreF[rSeqLength][minIndex2];
	} /* End right segment. */


	/* Update total error value by adding error values of left and right
	 * segments. */
	totalError = error + error1;

	/* This is for debugging purposes. */
	if (debugMode && totalError != error2 + error3
		&& totalError > errThreshold)
	{
		for (i = 0; i < lSeqLength; i++)
		  	printf("%c", *(tempref - 1 - i));
		
		printf("\n");
		
		for (i = 0; i < lSeqLength; i++)
		  	printf("%c", *(lSeq + i));
		
		printf("\n");

		for (i = 0; i < rSeqLength; i++)
		  	printf("%c", *(tempref + segLength + i));
		
		printf("\n");

		for (i = 0; i < rSeqLength; i++)
		  	printf("%c", *(rSeq + i));
		
		printf("\n");

		printf("SSEF=%d SSEB%d\n", error2, error3);
		printf("F=%d B=%d\n", error, error1);
		//scanf("%d", &i);
	}


	/* Now compute matching information. */
	char matrixR[SEQ_MAX_LENGTH];
	char matrixL[SEQ_MAX_LENGTH];

	matrixR[0] = '\0';
	matrixL[0] = '\0';

	size		   = 0;
	directionIndex = rSeqLength;
	rIndex		   = minIndex2;
	
	
	/* Compute matching information of right segment. */
	while (directionIndex != 0 || rIndex != 0)
	{
		if (directionIndex - rIndex == errThreshold)
		{
			if (scoreF[directionIndex][rIndex]
			- scoreF[directionIndex - 1][rIndex] == 1)
			{
		  		matrixR[size] = *(rSeq + directionIndex - 1);
		  		size++;
		  		matrixR[size] = 'I';
		  		directionIndex--;
			}
			else if (scoreF[directionIndex][rIndex]
			   - scoreF[directionIndex - 1][rIndex - 1] == 1)
			{
		  		matrixR[size] = *(ref + rIndex - 1);
		  		rIndex--;
		  		directionIndex--;
			}
			else
			{
		  		matrixR[size] = 'M';
		  		rIndex--;
		  		directionIndex--;
			}

		}
		else if (rIndex - directionIndex == errThreshold)
		{
			if (scoreF[directionIndex][rIndex]
			- scoreF[directionIndex][rIndex - 1] == 1)
			{
		  		matrixR[size] = *(ref + rIndex - 1);
		  		size++;
		  		matrixR[size] = 'D';
		  		rIndex--;
			}
			else if (scoreF[directionIndex][rIndex]
			   - scoreF[directionIndex - 1][rIndex - 1] == 1)
			{
		  		matrixR[size] = *(ref + rIndex - 1);
		  		rIndex--;
		  		directionIndex--;
			}
			else
			{
		  		matrixR[size] = 'M';
		  		rIndex--;
		  		directionIndex--;
		  }
		}
		else
		{
			if (scoreF[directionIndex][rIndex] -
				scoreF[directionIndex - 1][rIndex] == 1
			&& directionIndex != 0)
			{
		  		matrixR[size] = *(rSeq + directionIndex - 1);
		  		size++;
		  		matrixR[size] = 'I';
		  		directionIndex--;
			}
			else if (scoreF[directionIndex][rIndex]
			   - scoreF[directionIndex][rIndex - 1] == 1 && rIndex != 0)
			{
		  		matrixR[size] = *(ref + rIndex - 1);
		  		size++;
		  		matrixR[size] = 'D';
		  		rIndex--;
			}
			else if (scoreF[directionIndex][rIndex]
			   - scoreF[directionIndex - 1][rIndex - 1] == 1)
			{
		  		matrixR[size] = *(ref + rIndex - 1);
		  		rIndex--;
		  		directionIndex--;
			}
			else
			{
		  		matrixR[size] = 'M';
		  		rIndex--;
		  		directionIndex--;
			}
		}
		size++;
	}
	
	matrixR[size] = '\0';

	
	/* Compute matching information of left segment. */
	size		   = 0;
	directionIndex = lSeqLength;
	rIndex		   = minIndex1;

	while (directionIndex != 0 || rIndex != 0)
	{
		if (directionIndex - rIndex == errThreshold)
		{
			if (scoreB[directionIndex][rIndex] -
				scoreB[directionIndex - 1][rIndex] == 1)
			{
		  		matrixL[size] = 'I';
		  		size++;
		  		matrixL[size] = *(lSeq + lSeqLength - directionIndex);
		  		directionIndex--;
			}
			else if (scoreB[directionIndex][rIndex] -
					 scoreB[directionIndex - 1][rIndex - 1] == 1)
			{
		  		matrixL[size] = *(tempref - rIndex);
		  		rIndex--;
		  		directionIndex--;
			}
			else
			{
		  		matrixL[size] = 'M';
		  		rIndex--;
		  		directionIndex--;
			}

		}
		else if (rIndex - directionIndex == errThreshold)
		{
			if (scoreB[directionIndex][rIndex] -
				scoreB[directionIndex][rIndex - 1] == 1)
			{
		  		matrixL[size] = 'D';
		  		size++;
		  		matrixL[size] = *(tempref - rIndex);
		  		rIndex--;
			}
			else if (scoreB[directionIndex][rIndex] -
					 scoreB[directionIndex - 1][rIndex - 1] == 1)
			{
		  		matrixL[size] = *(tempref - rIndex);
		  		rIndex--;
		  		directionIndex--;
			}
			else
			{
		  		matrixL[size] = 'M';
		  		rIndex--;
		  		directionIndex--;
			}
		}
		else
		{
			if (scoreB[directionIndex][rIndex] -
				scoreB[directionIndex - 1][rIndex] == 1 &&
				directionIndex != 0)
			{
		  		matrixL[size] = 'I';
		  		size++;
		  		matrixL[size] = *(lSeq + lSeqLength - directionIndex);
		  		directionIndex--;
			}
			else if (scoreB[directionIndex][rIndex]  -
					 scoreB[directionIndex][rIndex - 1] == 1 &&
					 rIndex != 0)
			{
		  		matrixL[size] = 'D';
		  		size++;
		  		matrixL[size] = *(tempref - rIndex);
		  		rIndex--;
			}
			else if (scoreB[directionIndex][rIndex]  -
					 scoreB[directionIndex - 1][rIndex - 1] == 1)
			{
		  		matrixL[size] = *(tempref - rIndex);
		  		rIndex--;
		  		directionIndex--;
			}
			else
			{
		  		matrixL[size] = 'M';
		  		rIndex--;
		  		directionIndex--;
			}
		}
		size++;
	}
	
	matrixL[size] = '\0';


	/* Middle segment already matched to the reference genome, so make it all
	 * M */
	char middle[SEQ_MAX_LENGTH];
	middle[0] = '\0';
	for (i = 0; i < segLength; i++)
	  	middle[i] = 'M';
	middle[segLength] = '\0';

	/* Reverse right sequence values. */
	char rmatrixR[SEQ_MAX_LENGTH];
	reverse(matrixR, rmatrixR, strlen(matrixR));
	
	sprintf(matrix, "%s%s%s", matrixL, middle, rmatrixR);

	return totalError;
}



/*------------------------------------------------------------------------------
 * -SKIPPED-
 *	e = 4
 *----------------------------------------------------------------------------*/
int
verifySingleEndEditDistance4 ( int refIndex, char *lSeq, int lSeqLength,
							   char *rSeq, int rSeqLength, int segLength,
							   char *matrix, int *map_location,
							   short *seqHashValue )
{
  int i = 0;

  char * ref;
  char * tempref;

  int rIndex = 0; //reference Index

  int error = 0;
  int error1 = 0;

  int error2 = 0;
  int error3 = 0;
  int totalError = 0;
  //int errorSegment = 0;

  int ERROR_BOUND = errThreshold;

  /*
    1: Up
    2: Side
    3: Diagonal Match
    4: Diagonal Mismatch
  */

  int min = 0;
  int minIndex1 = 0;
  int minIndex2 = 0;

  int directionIndex = 0;

  int size = 0;

  ref = _msf_refGen + refIndex - 1;
  tempref = _msf_refGen + refIndex - 1;

  if (lSeqLength != 0) {
    error3 = backwardEditDistance4SSE2(ref - 1, lSeqLength,
				       lSeq + lSeqLength - 1, lSeqLength);
    if (error3 == -1) {
      return -1;
    }
  }

  if (rSeqLength != 0) {
    error2 = forwardEditDistance4SSE2(ref + segLength, rSeqLength, rSeq,
				      rSeqLength);
    if (error2 == -1)
      return -1;
  }

  if (error2 + error3 > errThreshold)
    return -1;

  rIndex = 1;

  //int prevError = 0;

  int tempUp = 0;
  int tempDown = 0;

  int errorString = 0;

  int upValue;
  int diagValue;
  int sideValue;

  while (rIndex <= lSeqLength + errThreshold && lSeqLength != 0) {
    tempUp =
      ((rIndex - ERROR_BOUND) > 0 ?
       ((rIndex > lSeqLength) ?
	lSeqLength - ERROR_BOUND : rIndex - ERROR_BOUND) :
       1);
    tempDown = (
		(rIndex >= lSeqLength - ERROR_BOUND) ?
		lSeqLength + 1 : rIndex + ERROR_BOUND + 1);
    for (i = tempUp; i < tempDown; i++) {
      errorString = (*(ref - rIndex) == *(lSeq + lSeqLength - i));

      upValue = scoreB[i - 1][rIndex] + 1;
      diagValue = scoreB[i - 1][rIndex - 1] + !errorString;
      sideValue = scoreB[i][rIndex - 1] + 1;

      if (i != tempUp && i != tempDown - 1)
	scoreB[i][rIndex] = min3(sideValue, diagValue , upValue);

      else if ((i
		== ((rIndex - ERROR_BOUND) > 0 ? rIndex - ERROR_BOUND : 1))
	       && rIndex <= lSeqLength)
	scoreB[i][rIndex] = min(sideValue, diagValue);
      else if (rIndex > lSeqLength && (i == lSeqLength - ERROR_BOUND))
	scoreB[i][rIndex] = sideValue;
      else
	scoreB[i][rIndex] = min(diagValue , upValue);

      if (i == tempUp)
	error = scoreB[i][rIndex];
      else if (error > scoreB[i][rIndex])
	error = scoreB[i][rIndex];
    }
    if (rIndex <= lSeqLength) {
      //errorSegment = error-prevError;
    }
    rIndex++;
  }

  if (lSeqLength != 0) {
    min = scoreB[lSeqLength][lSeqLength + errThreshold];
    minIndex1 = lSeqLength + errThreshold;

    // Find the Best error for all the possible ways.
    for (i = 1; i <= 2 * errThreshold; i++) {
      if (min >= scoreB[lSeqLength][lSeqLength + errThreshold - i]
	  && lSeqLength + errThreshold - i > 0) {
	min = scoreB[lSeqLength][lSeqLength + errThreshold - i];
	minIndex1 = lSeqLength + errThreshold - i;
      }
    }
    error = scoreB[lSeqLength][minIndex1];
  }

  error1 = error;

  error = 0;
  //errorSegment = 0;

  directionIndex = lSeqLength;
  rIndex = minIndex1;

  *map_location = ((lSeqLength == 0) ? refIndex : refIndex - rIndex);

  ref = ref + segLength;

  if (rSeqLength != 0) {
    rIndex = 1;
    while (rIndex <= rSeqLength + errThreshold - error1) {
      tempUp =
	(rIndex - ERROR_BOUND) > 0 ?
	((rIndex > rSeqLength) ?
	 rSeqLength - ERROR_BOUND :
	 rIndex - ERROR_BOUND) :
	1;
      tempDown = (
		  (rIndex >= rSeqLength - ERROR_BOUND) ?
		  rSeqLength + 1 : rIndex + ERROR_BOUND + 1);
      for (i = tempUp; i < tempDown; i++) {
	errorString = (*(ref + rIndex - 1) == *(rSeq + i - 1));

	upValue = scoreF[i - 1][rIndex] + 1;
	diagValue = scoreF[i - 1][rIndex - 1] + !errorString;
	sideValue = scoreF[i][rIndex - 1] + 1;

	if (i != tempUp && i != tempDown - 1)
	  scoreF[i][rIndex] = min3(sideValue, diagValue , upValue);
	else if ((i
		  == ((rIndex - ERROR_BOUND) > 0 ?
		      rIndex - ERROR_BOUND : 1))
		 && rIndex <= rSeqLength)
	  scoreF[i][rIndex] = min(sideValue, diagValue);
	else if (rIndex > rSeqLength && (i == rSeqLength - ERROR_BOUND))
	  scoreF[i][rIndex] = sideValue;
	else
	  scoreF[i][rIndex] = min(diagValue , upValue);

	if (i == tempUp)
	  error = scoreF[i][rIndex];
	if (error > scoreF[i][rIndex])
	  error = scoreF[i][rIndex];
      }
      if (rIndex <= rSeqLength) {
	//errorSegment = error;
      }

      rIndex++;
    }

    min = scoreF[rSeqLength][rSeqLength + errThreshold - error1];
    minIndex2 = rSeqLength + errThreshold - error1;

    // Find the Best error for all the possible ways.
    for (i = 1; i <= 2 * (errThreshold - error1); i++) {
      if (min > scoreF[rSeqLength][rSeqLength + errThreshold - error1 - i]
	  && rSeqLength + errThreshold - error1 - i > 0) {
	min =
	  scoreF[rSeqLength][rSeqLength + errThreshold - error1
			     - i];
	minIndex2 = rSeqLength + errThreshold - error1 - i;
      }
    }
    error = scoreF[rSeqLength][minIndex2];
  }

  totalError = error + error1;

  if (errThreshold > 4)
    printf("ERROR in errorThreshold.\n");

  if (totalError != error2 + error3 && totalError > errThreshold) {
    printf("ErrorF=%d, ErrorB=%d Error=%d Error=%d\n", error2, error3,
	   error1, error);

    //scanf("%d", &i);
  }

  char matrixR[SEQ_MAX_LENGTH];
  char matrixL[SEQ_MAX_LENGTH];

  matrixR[0] = '\0';
  matrixL[0] = '\0';

  size = 0;
  directionIndex = rSeqLength;
  rIndex = minIndex2;

  while (directionIndex != 0 || rIndex != 0) {
    if (directionIndex - rIndex == errThreshold) {
      if (scoreF[directionIndex][rIndex]
	  - scoreF[directionIndex - 1][rIndex] == 1) {
	matrixR[size] = *(rSeq + directionIndex - 1);
	size++;
	matrixR[size] = 'I';
	directionIndex--;
      } else if (scoreF[directionIndex][rIndex]
		 - scoreF[directionIndex - 1][rIndex - 1] == 1) {
	matrixR[size] = *(ref + rIndex - 1);
	rIndex--;
	directionIndex--;
      } else {
	matrixR[size] = 'M';
	rIndex--;
	directionIndex--;
      }

    } else if (rIndex - directionIndex == errThreshold) {
      if (scoreF[directionIndex][rIndex]
	  - scoreF[directionIndex][rIndex - 1] == 1) {
	matrixR[size] = *(ref + rIndex - 1);
	size++;
	matrixR[size] = 'D';
	rIndex--;
      } else if (scoreF[directionIndex][rIndex]
		 - scoreF[directionIndex - 1][rIndex - 1] == 1) {
	matrixR[size] = *(ref + rIndex - 1);
	rIndex--;
	directionIndex--;
      } else {
	matrixR[size] = 'M';
	rIndex--;
	directionIndex--;
      }
    } else {
      if (scoreF[directionIndex][rIndex]
	  - scoreF[directionIndex - 1][rIndex] == 1
	  && directionIndex != 0) {
	matrixR[size] = *(rSeq + directionIndex - 1);
	size++;
	matrixR[size] = 'I';
	directionIndex--;
      } else if (scoreF[directionIndex][rIndex]
		 - scoreF[directionIndex][rIndex - 1] == 1 && rIndex != 0) {
	matrixR[size] = *(ref + rIndex - 1);
	size++;
	matrixR[size] = 'D';
	rIndex--;
      } else if (scoreF[directionIndex][rIndex]
		 - scoreF[directionIndex - 1][rIndex - 1] == 1) {
	matrixR[size] = *(ref + rIndex - 1);
	rIndex--;
	directionIndex--;
      } else {
	matrixR[size] = 'M';
	rIndex--;
	directionIndex--;
      }
    }
    size++;
  }
  matrixR[size] = '\0';

  size = 0;
  directionIndex = lSeqLength;
  rIndex = minIndex1;

  while (directionIndex != 0 || rIndex != 0) {
    if (directionIndex - rIndex == errThreshold) {
      if (scoreB[directionIndex][rIndex]
	  - scoreB[directionIndex - 1][rIndex] == 1) {
	matrixL[size] = 'I';
	size++;
	matrixL[size] = *(lSeq + lSeqLength - directionIndex);
	directionIndex--;
      } else if (scoreB[directionIndex][rIndex]
		 - scoreB[directionIndex - 1][rIndex - 1] == 1) {
	matrixL[size] = *(tempref - rIndex);
	rIndex--;
	directionIndex--;
      } else {
	matrixL[size] = 'M';
	rIndex--;
	directionIndex--;
      }

    } else if (rIndex - directionIndex == errThreshold) {
      if (scoreB[directionIndex][rIndex]
	  - scoreB[directionIndex][rIndex - 1] == 1) {
	matrixL[size] = 'D';
	size++;
	matrixL[size] = *(tempref - rIndex);
	rIndex--;
      } else if (scoreB[directionIndex][rIndex]
		 - scoreB[directionIndex - 1][rIndex - 1] == 1) {
	matrixL[size] = *(tempref - rIndex);
	rIndex--;
	directionIndex--;
      } else {
	matrixL[size] = 'M';
	rIndex--;
	directionIndex--;
      }
    } else {
      if (scoreB[directionIndex][rIndex]
	  - scoreB[directionIndex - 1][rIndex] == 1
	  && directionIndex != 0) {
	matrixL[size] = 'I';
	size++;
	matrixL[size] = *(lSeq + lSeqLength - directionIndex);
	directionIndex--;
      } else if (scoreB[directionIndex][rIndex]
		 - scoreB[directionIndex][rIndex - 1] == 1 && rIndex != 0) {
	matrixL[size] = 'D';
	size++;
	matrixL[size] = *(tempref - rIndex);
	rIndex--;
      } else if (scoreB[directionIndex][rIndex]
		 - scoreB[directionIndex - 1][rIndex - 1] == 1) {
	matrixL[size] = *(tempref - rIndex);
	rIndex--;
	directionIndex--;
      } else {
	matrixL[size] = 'M';
	rIndex--;
	directionIndex--;
      }
    }

    size++;
  }

  matrixL[size] = '\0';
  char middle[SEQ_MAX_LENGTH];
  middle[0] = '\0';

  for (i = 0; i < segLength; i++)
    middle[i] = 'M';
  middle[segLength] = '\0';

  char rmatrixR[SEQ_MAX_LENGTH];

  reverse(matrixR, rmatrixR, strlen(matrixR));

  sprintf(matrix, "%s%s%s", matrixL, middle, rmatrixR);

  return totalError;
}



/*------------------------------------------------------------------------------
 * -FULLY COMMENTED-
 *	e > 4
 *----------------------------------------------------------------------------*/
int
verifySingleEndEditDistanceExtension ( int refIndex, char *lSeq, int lSeqLength,
									   char *rSeq, int rSeqLength,
									   int segLength, char *matrix,
									   int *map_location, short *seqHashValue )
{
	int		 i			= 0;
	char	*ref;
	char	*tempref;
	int		 rIndex		= 0;	//reference Index
	int		 error		= 0;
	int		 error1		= 0;
	int		 error2		= 0;
	int		 error3		= 0;
	int		 totalError = 0;
	//int errorSegment = 0;

	int ERROR_BOUND = min(4, errThreshold);

	/*
	  1: Up
	  2: Side
	  3: Diagonal Match
	  4: Diagonal Mismatch
	*/

	int min			   = 0;
	int minIndex1	   = 0;
	int minIndex2	   = 0;
	int directionIndex = 0;
	int size		   = 0;

	ref		= _msf_refGen + refIndex - 1;
	tempref = _msf_refGen + refIndex - 1;

	/* Match left segment */
	if (lSeqLength != 0)
	{	
		error3 = backwardEditDistanceSSE2Extension(ref - 1, lSeqLength,
												   lSeq + lSeqLength - 1,
												   lSeqLength);
		if (error3 == -1)
		{
		  	return -1;
		}
	}

	/* Match right segment */
	if (rSeqLength != 0)
	{		
		error2 = forwardEditDistanceSSE2Extension(ref + segLength, rSeqLength,
												  rSeq, rSeqLength);
		if (error2 == -1)
		{
		  	return -1;
		}
	}
	
	/* If the edit distance of left seq and right seq summed is greater than the
	 * threshold, could not match */
	if (error2 + error3 > errThreshold)
	{	
	  	return -1;
	}

	//int prevError = 0;
	rIndex			= 1;
	int tempUp		= 0;
	int tempDown	= 0;
	int errorString = 0;
	int upValue;
	int diagValue;
	int sideValue;

	/* Process the left segment to find score and match, note that this does not
	 * use SSE instructions. Simple and plain DP. rIndex is the anti-diagonal
	 * index. */
	if (lSeqLength > ERROR_BOUND)
	{
		
		while (rIndex <= lSeqLength + ERROR_BOUND && lSeqLength != 0)
		{
			/* Find starting and finishing indexes of this anti-diagonal. */
			tempUp = ((rIndex - ERROR_BOUND) > 0 ?
					  ((rIndex > lSeqLength) ?
					   lSeqLength - ERROR_BOUND : rIndex - ERROR_BOUND) : 1);
			tempDown = ((rIndex >= lSeqLength - ERROR_BOUND) ?
						lSeqLength + 1 : rIndex + ERROR_BOUND + 1);

			/* Compute anti-diagonal entries. */
			for (i = tempUp; i < tempDown; i++)
			{
				/* Find the value of the newly computed entry by choosing min of
				 * three values. */
				errorString = (*(ref - rIndex) == *(lSeq + lSeqLength - i));
				upValue		= scoreB[i - 1][rIndex] + 1;
				diagValue	= scoreB[i - 1][rIndex - 1] + !errorString;
				sideValue	= scoreB[i][rIndex - 1] + 1;

				/* Update score matrix */
				if (i != tempUp && i != tempDown - 1)
				  	scoreB[i][rIndex] = min3(sideValue, diagValue , upValue);

				else if ((i == ((rIndex - ERROR_BOUND) > 0 ?
								rIndex - ERROR_BOUND : 1)) &&
						 rIndex <= lSeqLength)
				  	scoreB[i][rIndex] = min(sideValue, diagValue);
				else if (rIndex > lSeqLength && (i == lSeqLength - ERROR_BOUND))
				  	scoreB[i][rIndex] = sideValue;
				else
				  	scoreB[i][rIndex] = min(diagValue , upValue);

				if (i == tempUp)
				  	error = scoreB[i][rIndex];
				else if (error > scoreB[i][rIndex])
				  	error = scoreB[i][rIndex];
			}

			/* I wouldnt leave that here. */
			if (rIndex <= lSeqLength)
			{
		  		//errorSegment = error-prevError;
			}
			rIndex++;
		}

		/* Find the final error value. */
		if (lSeqLength != 0)
		{
			min = scoreB[lSeqLength][lSeqLength + ERROR_BOUND];
			minIndex1 = lSeqLength + ERROR_BOUND;

			// Find the Best error for all the possible ways.
			for (i = 1; i <= 2 * ERROR_BOUND; i++)
			{
				if (min >= scoreB[lSeqLength][lSeqLength + ERROR_BOUND - i]
					&& lSeqLength + ERROR_BOUND - i > 0)
				{
				  	min = scoreB[lSeqLength][lSeqLength + ERROR_BOUND - i];
				  	minIndex1 = lSeqLength + ERROR_BOUND - i;
				}
			}
			error = scoreB[lSeqLength][minIndex1];
		}
	}
	/* seq length is smaller than the ERROR_BOUND, run a simpler DP. */
	else
	{
		int j = 0;
		for (i = 1; i <= lSeqLength; i++)
		{
		  	for (j = 1; j <= lSeqLength; j++)
			{
				scoreB[i][j] = min3(scoreB[i-1][j-1]+ (*(ref-j) != *(lSeq+lSeqLength-i) ),
									scoreB[i][j-1]+1 ,
									scoreB[i-1][j]+1);
		  	}
		}
		
		error	  = scoreB[lSeqLength][lSeqLength];
		minIndex1 = lSeqLength;
	}

	/* Update mapping location, by also taking indels into account. */
	error1		   = error;
	error		   = 0;
	//errorSegment = 0;
	directionIndex = lSeqLength;
	rIndex		   = minIndex1;
	*map_location  = ((lSeqLength == 0) ? refIndex : refIndex - rIndex);

	
	/* ********************************************************************** */
	/* ********************************************************************** */
	/* ********************************************************************** */
	

	/* Trace back right sequence (scoreF). */
	ref = ref + segLength;
	if (rSeqLength != 0 && rSeqLength > ERROR_BOUND)
	{	
		ERROR_BOUND = min(ERROR_BOUND, rSeqLength);

		if (rSeqLength == ERROR_BOUND)
		{
		  	for (i = 0; i < 2 * ERROR_BOUND; i++)
				scoreF[0][i] = i;
		}

		/* rIndex is the anti-diagonal index. */
		rIndex = 1;
		while (rIndex <= rSeqLength + ERROR_BOUND)
		{
			/* Find starting and finishing indexes of this anti-diagonal. */
			tempUp = (rIndex - ERROR_BOUND) > 0 ?
				((rIndex > rSeqLength) ?
				 rSeqLength - ERROR_BOUND : rIndex - ERROR_BOUND) : 1;
			tempDown = ((rIndex >= rSeqLength - ERROR_BOUND) ?
						rSeqLength + 1 : rIndex + ERROR_BOUND + 1);

			/* Compute anti-diagonal entries. */
			for (i = tempUp; i < tempDown; i++)
			{
				/* Find the value of the newly computed entry by choosing min of
				 * three values. */
				errorString = (*(ref + rIndex - 1) == *(rSeq + i - 1));
				upValue		= scoreF[i - 1][rIndex] + 1;
				diagValue	= scoreF[i - 1][rIndex - 1] + !errorString;
				sideValue	= scoreF[i][rIndex - 1] + 1;

				/* Update score matrix */
				if (i != tempUp && i != tempDown - 1)
				  	scoreF[i][rIndex] = min3(sideValue, diagValue , upValue);
				else if ((i == ((rIndex - ERROR_BOUND) > 0 ?
								rIndex - ERROR_BOUND : 1)) &&
						 rIndex <= rSeqLength)
				  	scoreF[i][rIndex] = min(sideValue, diagValue);
				else if (rIndex > rSeqLength && (i == rSeqLength - ERROR_BOUND))
				  	scoreF[i][rIndex] = sideValue;
				else
				  	scoreF[i][rIndex] = min(diagValue , upValue);

				if (i == tempUp)
				  	error = scoreF[i][rIndex];
				if (error > scoreF[i][rIndex])
				  	error = scoreF[i][rIndex];
			}

			/* Not to mention how if branches are ruining me performance. */
			if (rIndex <= rSeqLength)
			{
		  		//errorSegment = error;
			}
			rIndex++;
		}

		/* Find the final error value. */
		min		  = scoreF[rSeqLength][rSeqLength + ERROR_BOUND];
		minIndex2 = rSeqLength + ERROR_BOUND;

		// Find the Best error for all the possible ways.
		for (i = 1; i <= 2 * ERROR_BOUND; i++)
		{
			if (min > scoreF[rSeqLength][rSeqLength + ERROR_BOUND - i] &&
				rSeqLength + ERROR_BOUND - i > 0)
			{
		  		min = scoreF[rSeqLength][rSeqLength + ERROR_BOUND - i];
		  		minIndex2 = rSeqLength + ERROR_BOUND - i;
			}
		}
		error = scoreF[rSeqLength][minIndex2];
	}
	else						/* segment length < ERROR_BOUND*/
	{
		int j = 0;
		for (i = 1; i <= rSeqLength; i++)
		{
			for (j = 1; j <= rSeqLength; j++)
			{
		  		scoreF[i][j] = min3(scoreF[i-1][j-1]+ (*(ref+j-1) != *(rSeq+i-1) ),
									scoreF[i][j-1]+1 ,
									scoreF[i-1][j]+1);
			}
		}
		
		error = scoreF[rSeqLength][rSeqLength];
		minIndex2 = rSeqLength;
	}

	/* Update total error value by adding error values of left and right
	 * segments. */
	totalError = error + error1;

	/* Farhad 08/07/2012 */
	if(totalError > errThreshold)
	  	return -1;
	/* Farhad 08/07/2012 */

	/* This is for debugging purposes. */
	if (debugMode && totalError != error2 + error3)
	{
		for (i = 0; i < lSeqLength; i++)
		  	printf("%c", *(tempref - 1 - i));
		
		printf("\n");
		
		for (i = 0; i < lSeqLength; i++)
		  	printf("%c", *(lSeq + i));
		
		printf("\n");

		for (i = 0; i < rSeqLength; i++)
		  	printf("%c", *(tempref + segLength + i));
		
		printf("\n");

		for (i = 0; i < rSeqLength; i++)
		  	printf("%c", *(rSeq + i));
		
		printf("\n");

		printf("ERROR=%d\n", totalError);
		printf("ERROR_SSE=%d\n", error3 + error2);

		printf("ERROR_SSE_back=%d E_SSE_forw=%d\n", error3, error2);
		printf("ERROR_back=%d E_forw=%d\n", error1, error);
	}

	/* Now compute matching information. */
	char matrixR[SEQ_MAX_LENGTH];
	char matrixL[SEQ_MAX_LENGTH];

	matrixR[0] = '\0';
	matrixL[0] = '\0';

	size		   = 0;
	directionIndex = rSeqLength;
	rIndex		   = minIndex2;

	
	/* Farhad 29/11/2012 */
	/* Compute matching information of right segment. */
	while (directionIndex > 0 || rIndex > 0)
	{
		if (directionIndex - rIndex == errThreshold)
		{
			if (scoreF[directionIndex][rIndex] -
				scoreF[directionIndex - 1][rIndex] == 1)
			{
		  		matrixR[size] = *(rSeq + directionIndex - 1);
		  		size++;
		  		matrixR[size] = 'I';
		  		directionIndex--;
			}
			else if (scoreF[directionIndex][rIndex] -
					 scoreF[directionIndex - 1][rIndex - 1] == 1)
			{
		  		matrixR[size] = *(ref + rIndex - 1);
		  		rIndex--;
		  		directionIndex--;
			}
			else
			{
		  		matrixR[size] = 'M';
		  		rIndex--;
		  		directionIndex--;
			}
		}
		else if (rIndex - directionIndex == errThreshold)
		{
			if (scoreF[directionIndex][rIndex] -
				scoreF[directionIndex][rIndex - 1] == 1)
			{
		  		matrixR[size] = *(ref + rIndex - 1);
		  		size++;
		  		matrixR[size] = 'D';
		  		rIndex--;
			}
			else if (scoreF[directionIndex][rIndex] -
					 scoreF[directionIndex - 1][rIndex - 1] == 1)
			{
		  		matrixR[size] = *(ref + rIndex - 1);
		  		rIndex--;
		  		directionIndex--;
			}
			else
			{
		  		matrixR[size] = 'M';
		  		rIndex--;
		  		directionIndex--;
			}
		}
		else
		{
			if (scoreF[directionIndex][rIndex] -
				scoreF[directionIndex - 1][rIndex] == 1 && directionIndex != 0)
			{
		  		matrixR[size] = *(rSeq + directionIndex - 1);
		  		size++;
		  		matrixR[size] = 'I';
		  		directionIndex--;
			}
			else if (scoreF[directionIndex][rIndex] -
					 scoreF[directionIndex][rIndex - 1] == 1 && rIndex != 0)
			{
		  		matrixR[size] = *(ref + rIndex - 1);
		  		size++;
		  		matrixR[size] = 'D';
		  		rIndex--;
			}
			else if (scoreF[directionIndex][rIndex] -
					 scoreF[directionIndex - 1][rIndex - 1] == 1)
			{
		  		matrixR[size] = *(ref + rIndex - 1);
		  		rIndex--;
		  		directionIndex--;
			}
			else
			{
		  		matrixR[size] = 'M';
		  		rIndex--;
		  		directionIndex--;
			}
		}
		
		size++;
	}
	
	matrixR[size] = '\0';

	/* Compute matching information of left segment. */
	size		   = 0;
	directionIndex = lSeqLength;
	rIndex		   = minIndex1;

	while (directionIndex != 0 || rIndex != 0)
	{
		if (directionIndex - rIndex == errThreshold)
		{
			if (scoreB[directionIndex][rIndex] -
				scoreB[directionIndex - 1][rIndex] == 1)
			{
		  		matrixL[size] = 'I';
		  		size++;
		  		matrixL[size] = *(lSeq + lSeqLength - directionIndex);
		  		directionIndex--;
			}
			else if (scoreB[directionIndex][rIndex] -
					 scoreB[directionIndex - 1][rIndex - 1] == 1)
			{
		  		matrixL[size] = *(tempref - rIndex);
		  		rIndex--;
		  		directionIndex--;
			}
			else
			{
		  		matrixL[size] = 'M';
		  		rIndex--;
		  		directionIndex--;
			}
		}
		else if (rIndex - directionIndex == errThreshold)
		{
			if (scoreB[directionIndex][rIndex] -
				scoreB[directionIndex][rIndex - 1] == 1)
			{
		  		matrixL[size] = 'D';
		  		size++;
		  		matrixL[size] = *(tempref - rIndex);
		  		rIndex--;
			}
			else if (scoreB[directionIndex][rIndex] -
					 scoreB[directionIndex - 1][rIndex - 1] == 1)
			{
		  		matrixL[size] = *(tempref - rIndex);
		  		rIndex--;
		  		directionIndex--;
			}
			else
			{
		  		matrixL[size] = 'M';
		  		rIndex--;
		  		directionIndex--;
			}
		}
		else
		{
			if (scoreB[directionIndex][rIndex] -
				scoreB[directionIndex - 1][rIndex] == 1	&&
				directionIndex != 0)
			{
		  		matrixL[size] = 'I';
		  		size++;
		  		matrixL[size] = *(lSeq + lSeqLength - directionIndex);
		  		directionIndex--;
			}
			else if (scoreB[directionIndex][rIndex] -
					 scoreB[directionIndex][rIndex - 1] == 1 &&
					 rIndex != 0)
			{
		  		matrixL[size] = 'D';
		  		size++;
		  		matrixL[size] = *(tempref - rIndex);
		  		rIndex--;
			}
			else if (scoreB[directionIndex][rIndex] -
					 scoreB[directionIndex - 1][rIndex - 1] == 1)
			{
		  		matrixL[size] = *(tempref - rIndex);
		  		rIndex--;
		  		directionIndex--;
			}
			else
			{
		  		matrixL[size] = 'M';
		  		rIndex--;
		  		directionIndex--;
			}
		}
		
		size++;
	}
	
	matrixL[size] = '\0';

	
	/* Middle segment already matched to the reference genome, so make it all
	 * M */
	char middle[SEQ_MAX_LENGTH];
	middle[0] = '\0';
	for (i = 0; i < segLength; i++)
	  	middle[i] = 'M';
	middle[segLength] = '\0';

	/* Reverse right sequence values. */
	char rmatrixR[SEQ_MAX_LENGTH];
	reverse(matrixR, rmatrixR, strlen(matrixR));

	sprintf(matrix, "%s%s%s", matrixL, middle, rmatrixR);

	return totalError;
}



/*------------------------------------------------------------------------------
 * Helping function to the generateCigar function.
 *----------------------------------------------------------------------------*/
int
addCigarSize ( int cnt )
{
	if (cnt < 10)
	  	return 1;
	else if (cnt < 100)
	  	return 2;
	return 3;
}



/*------------------------------------------------------------------------------
 * Generate CIGAR output using the matching information consisting of Match,
 * Delete and Insert values in the matrix variable. Write output to the cigar
 *----------------------------------------------------------------------------*/
void
generateCigar ( char *matrix, int matrixLength, char *cigar )
{
	int i		  = 0;
	int counterM  = 0;
	int counterI  = 0;
	int counterD  = 0;
	int cigarSize = 0;

	cigar[0] = '\0';

	while (i < matrixLength)
	{
		if (matrix[i] == 'M')
		{
			counterM++;
			if (counterI != 0)
			{
		  		sprintf(cigar, "%s%dI", cigar, counterI);
		  		cigarSize += addCigarSize(counterI) + 1;
		  		cigar[cigarSize] = '\0';
		  		counterI = 0;
			}
			else if (counterD != 0)
			{
		  		sprintf(cigar, "%s%dD", cigar, counterD);
		  		cigarSize += addCigarSize(counterD) + 1;
		  		cigar[cigarSize] = '\0';
		  		counterD = 0;
			}
		}
		else if (matrix[i] == 'I')
		{
			if (counterM != 0)
			{
		  		sprintf(cigar, "%s%dM", cigar, counterM);
		  		cigarSize += addCigarSize(counterM) + 1;
		  		cigar[cigarSize] = '\0';
		  		counterM = 0;
			}
			else if (counterD != 0)
			{
		  		sprintf(cigar, "%s%dD", cigar, counterD);
		  		cigarSize += addCigarSize(counterD) + 1;
		  		cigar[cigarSize] = '\0';
		  		counterD = 0;
			}
			counterI++;
			i++;
		}
		else if (matrix[i] == 'D')
		{
			if (counterM != 0)
			{
		  		sprintf(cigar, "%s%dM", cigar, counterM);
		  		cigarSize += addCigarSize(counterM) + 1;
		  		cigar[cigarSize] = '\0';
		  		counterM = 0;
			}
			else if (counterI != 0)
			{
		  		sprintf(cigar, "%s%dI", cigar, counterI);
		  		cigarSize += addCigarSize(counterI) + 1;
		  		cigar[cigarSize] = '\0';
		  		counterI = 0;
			}

			counterD++;
			i++;
		}
		else
		{
			counterM++;
			if (counterI != 0)
			{
		  		sprintf(cigar, "%s%dI", cigar, counterI);
		  		cigarSize += addCigarSize(counterI) + 1;
		  		cigar[cigarSize] = '\0';
		  		counterI = 0;
			}
			else if (counterD != 0)
			{
		  		sprintf(cigar, "%s%dD", cigar, counterD);
		  		cigarSize += addCigarSize(counterD) + 1;
		  		cigar[cigarSize] = '\0';
		  		counterD = 0;
			}
		}
		i++;
	}

	if (counterM != 0)
	{
		sprintf(cigar, "%s%dM", cigar, counterM);
		cigarSize += addCigarSize(counterM) + 1;
		cigar[cigarSize] = '\0';
		counterM = 0;
	}
	else if (counterI != 0)
	{
		sprintf(cigar, "%s%dI", cigar, counterI);
		cigarSize += addCigarSize(counterI) + 1;
		cigar[cigarSize] = '\0';
		counterI = 0;
	}
	else if (counterD != 0)
	{
		sprintf(cigar, "%s%dD", cigar, counterD);
		cigarSize += addCigarSize(counterD) + 1;
		cigar[cigarSize] = '\0';
		counterD = 0;
	}

	cigar[cigarSize] = '\0';
}



/*------------------------------------------------------------------------------
 * Generates CIGAR output from MD information
 *----------------------------------------------------------------------------*/
void
generateCigarFromMD ( char *mismatch, int mismatchLength, char *cigar )
{
	int i		  = 0;
	int j		  = 0;
	int start	  = 0;
	int cigarSize = 0;
	cigar[0]	  = '\0';

	while (i < mismatchLength)
	{
		if (mismatch[i] >= '0' && mismatch[i] <= '9')
		{
			start = i;
			while (mismatch[i] >= '0' && mismatch[i] <= '9' && i < mismatchLength)
		  		i++;

			int value = atoi(mismatch + start);
			for (j = 0; j < value - 1; j++)
			{
		  		cigar[cigarSize] = 'M';
		  		cigarSize++;
			}
			cigar[cigarSize] = 'M';
		}
		else if (mismatch[i] == '^')
		{
		  	cigar[cigarSize] = 'I';
		  	i++;
		}
		else if (mismatch[i] == '\'')
		{
		  	cigar[cigarSize] = 'D';
		  	i++;
		}
		else
		{
		  	cigar[cigarSize] = 'M';
		  	cigarSize++;
		}
		cigarSize++;
		i++;
	}
	cigar[cigarSize] = '\0';
}



/*------------------------------------------------------------------------------
 * Generate SAM output using the matching information consisting of Match,
 * Delete and Insert values in the matrix variable. Write output to the
 * outputSNP
 *----------------------------------------------------------------------------*/
void
generateSNPSAM ( char *matrix, int matrixLength, char *outputSNP )
{
	int		i		 = 0;
	int		counterM = 0;
	int		counterD = 0;
	char	delete[100];
	int		snpSize	 = 0;
	
	outputSNP[0] = '\0';
	delete[0]	 = '\0';

	/* Generate output information using the matching information M, D, I
	 * (Match, Delete, Insert) */
	while (i < matrixLength)
	{
		if (matrix[i] == 'M')
		{
			counterM++;
			if (counterD != 0)
			{
				delete[counterD] = '\0';
				counterD = 0;
				sprintf(outputSNP, "%s^%s", outputSNP, delete);
				snpSize += strlen(delete) + 1;
				outputSNP[snpSize] = '\0';
				delete[0] = '\0';
			}
		}
		else if (matrix[i] == 'D')
		{
			if (counterM != 0)
			{
				sprintf(outputSNP, "%s%d", outputSNP, counterM);
				snpSize += addCigarSize(counterM);
				outputSNP[snpSize] = '\0';
				counterM = 0;
				delete[counterD] = matrix[i + 1];
				i++;
				counterD++;
			}
			else if (counterD != 0)
			{
		  		delete[counterD] = matrix[i + 1];
		  		counterD++;
		  		i++;
			}
			else
			{
		  		delete[counterD] = matrix[i + 1];
		  		counterD++;
		  		i++;
			}
		}
		else if (matrix[i] == 'I')
		{
			if (counterM != 0)
			{
		  		// sprintf(outputSNP, "%s%d\0", outputSNP, counterM);
		  		//counterM++;
			}
			else if (counterD != 0)
			{
		  		delete[counterD] = '\0';
		  		sprintf(outputSNP, "%s^%s", outputSNP, delete);
		  		snpSize += strlen(delete) + 1;
		  		outputSNP[snpSize] = '\0';
		  		counterD = 0;
		  		delete[0] = '\0';
			}
			i++;
		}
		else
		{
			if (counterM != 0)
			{
		  		sprintf(outputSNP, "%s%d", outputSNP, counterM);
		  		snpSize += addCigarSize(counterM);
		  		outputSNP[snpSize] = '\0';
		  		counterM = 0;
			}
			if (counterD != 0)
			{
		  		delete[counterD] = '\0';
		  		counterD = 0;
		  		sprintf(outputSNP, "%s^%s", outputSNP, delete);
		  		snpSize += strlen(delete) + 1;
		  		outputSNP[snpSize] = '\0';
		  		delete[0] = '\0';
			}
			sprintf(outputSNP, "%s%c", outputSNP, matrix[i]);
			snpSize += 1;
			outputSNP[snpSize] = '\0';
		}
		i++;
	}

	if (counterM != 0)
	{
		sprintf(outputSNP, "%s%d", outputSNP, counterM);
		snpSize += addCigarSize(counterM);
		outputSNP[snpSize] = '\0';
		counterM = 0;
	}
	else if (counterD != 0)
	{
		delete[counterD] = '\0';
		sprintf(outputSNP, "%s^%s", outputSNP, delete);
		snpSize += strlen(delete) + 1;
		outputSNP[snpSize] = '\0';
		counterD = 0;
	}

	outputSNP[snpSize] = '\0';
}



/*------------------------------------------------------------------------------
 * Parameters:
 *	target_coor	: A position in the reference genome.
 *	entry_coor	: Locations of the current given read-segment.
 *	entry_size	: Length of the entry_size

 * OLD COMMENTS (not @OGUZ-COMMENT)
 *	MrFAST with fastHASH: searchKey()
 *----------------------------------------------------------------------------*/
int
searchKey ( int target_coor, unsigned int* entry_coor, int entry_size )
{
	if (entry_size <= 0)
	  	return 0;
	
	int lower_bound = 1;
	int upper_bound = entry_size;
	int mid			= lower_bound + entry_size / 2;

	while (lower_bound < upper_bound)
	{
		if (entry_coor[mid] <= target_coor + errThreshold
			&& entry_coor[mid] >= target_coor - errThreshold)
		  	break;
		else if (entry_coor[mid] < target_coor)
		  	lower_bound = mid + 1;
		else
		  	upper_bound = mid - 1;
		mid = lower_bound + (upper_bound - lower_bound) / 2;
	}

	if (entry_coor[mid] <= target_coor + errThreshold
		&& entry_coor[mid] >= target_coor - errThreshold)
	{
	  	return 1;
	}
	else
	  	return 0;
}



/*------------------------------------------------------------------------------
 * Parameters:
 *	l1					: The locations of this k-mer in the reference genome.
 *	s1					: Number of such locations.
 *	readNumber			: The id of this read in the global read seq array.
 *	readSegment			: The id of the current segment in this specific read
 						  (note that they were sorted so that this does not have
						  to be in order).
 *	direction			: 0 = FORWARD, 1 = BACKWARD
 *	index				: Which readSegment we are processing? (in order)
 *	keys_input			: All keys (readSegment or k-mer) of this read.
 *	potential_key_num	: Number of read-segments to be processed.

 * OLD COMMENTS (not @OGUZ-COMMENT)
 *	direction = 0 forward
 				1 backward
 * MrFAST with fastHASH: compareEntrySize()
 *----------------------------------------------------------------------------*/
void
mapSingleEndSeq ( unsigned int *l1, int s1, int readNumber, int readSegment,
				  int direction, int index, key_struct* keys_input,
				  int potential_key_number )
{
	int		 j	  = 0;
	int		 z	  = 0;
	int		*locs = (int *) l1;
	char	*_tmpSeq, *_tmpQual;
	char	 rqual[SEQ_LENGTH + 1];

	/* Left, middle and right segment lengths. */
	int genLoc			= 0;
	int leftSeqLength	= 0;
	int rightSeqLength	= 0;
	int middleSeqLength = 0;

	char	matrix[SEQ_MAX_LENGTH];
	char	editString[2 * SEQ_MAX_LENGTH];
	char	cigar[MAX_CIGAR_SIZE];

	short	*_tmpHashValue;
	int		 key_number = SEQ_LENGTH / WINDOW_SIZE;
	int		 realLoc, readId;

	rqual[SEQ_LENGTH] = '\0';

	/* With respect to chosen direction, either choose seq or rseq in the
	 * original read */
	if (direction)				/* 1 = BACKWARD */
	{
		reverse(_msf_seqList[readNumber].qual, rqual, SEQ_LENGTH);
		_tmpQual	  = rqual;
		_tmpSeq		  = _msf_seqList[readNumber].rseq;
		_tmpHashValue = _msf_seqList[readNumber].rhashValue;
	}
	else						/* 0 = FORWARD */
	{
		_tmpQual	  = _msf_seqList[readNumber].qual;
		_tmpSeq		  = _msf_seqList[readNumber].seq;
		_tmpHashValue = _msf_seqList[readNumber].hashValue;
	}

	/* readId is for _msf_verifiedLocs. We need to store two entries per read
	 * (FORWARD and BACKWARD), so x2. */
	readId = 2 * readNumber + direction;

	/* For each location this k-mer can be mapped to. */
	for ( z = 0; z < s1; z++)
	{
		int map_location = 0;
		int a			 = 0;
		int o			 = index;

		/* Get the current location of that read-segment will be mapped to in
		 * reference genome. */
		genLoc = locs[z];

		/* Assume you have mapped a k-mer to a location in ref genome. This
		 * location is genLoc. realLoc becomes the start of the whole read
		 * sequence at the ref genome with respect to genLoc. */
		
		//hxin: If the read is at the beginning of the contig, and there are insertions
		//hxin: then genLoc - _msf_sampleingLoc[0] might be smaller than _refGenBeg
		realLoc = genLoc - _msf_samplingLocs[o];

		/* Some corners cases about the beginning and end of the reference
		 * genome. */
		if (genLoc < _msf_samplingLocs[o]
			|| genLoc - _msf_samplingLocs[o] < _msf_refGenBeg
			|| genLoc - _msf_samplingLocs[o] > _msf_refGenEnd)
		{
		  	if (genLoc - _msf_samplingLocs[o] > _msf_refGenBeg - errThreshold)
				realLoc =  _msf_refGenBeg;
		  	else
				continue;
		}

		/* Only compare the start of the read sequence at the reference genome
		 * to ensure if this read was verified previously. */
		if (_msf_verifiedLocs[realLoc] == readId)
		{
		  	continue;
		}

		//Begin of long-K (long k-mer?)
		/* Tries to merge read-segments? It skips the current locs entry
		 * inspected if something is validated (Idx = index). */
		int mergeIdx = 2 * errThreshold + 1 - index;
		
		if (mergeIdx < potential_key_number)
		{
			if (!searchKey(
				   genLoc
				   + (keys_input[mergeIdx].key_number
				  - keys_input[o].key_number) * WINDOW_SIZE,
				   keys_input[mergeIdx].key_entry,
				   keys_input[mergeIdx].key_entry_size))
			{
			  	continue;
			}
		}
		//End of long-K

		// Adjacency Filtering Start ---------------------------------
		int skip_edit_distance = 0;
		int diff_num		   = 0;
		int ix				   = 0;
		for (ix = 0; ix < potential_key_number; ix++)
		{
			if (ix >= key_number - errThreshold)
			{
		  		break;
			}
			
			if (ix != o && ix != mergeIdx)
			{ // Changed with long-K
		  		if (!searchKey(
					 genLoc
					 + (keys_input[ix].key_number
					- keys_input[o].key_number)
					 * WINDOW_SIZE, keys_input[ix].key_entry,
					 keys_input[ix].key_entry_size))
		  		{
					diff_num++;
					if (diff_num > errThreshold)
					{
			  			skip_edit_distance = 1;
			  			break;
					}
		  		}
			}
		}
		// Adjacency Filtering End -----------------------------------

		
		int err		 = -1;
		map_location = 0;

		/* The variable "a" stores the starting point of right sequence. */
		leftSeqLength	= _msf_samplingLocs[o];
		middleSeqLength = WINDOW_SIZE;
		a				= leftSeqLength + middleSeqLength; 
		rightSeqLength	= SEQ_LENGTH - a;

		/* CALKAN: skip alignment if it is a perfect match 
		   this has to be re-addressed by Donghyuk

		   if (diff_num == 0) {
		   skip_edit_distance = 1;
		   err = 0;
		   sprintf(cigar, "%dM", SEQ_LENGTH);
		   sprintf(editString, "%d", SEQ_LENGTH);
		   }*/

		
		/* This is where the matching is done. */
		if (skip_edit_distance == 0)
		{
			if (errThreshold == 2) /* OK */
			{
		  		err = verifySingleEndEditDistance2(genLoc, _tmpSeq,
							 leftSeqLength, _tmpSeq + a, rightSeqLength,
							 middleSeqLength, matrix, &map_location, _tmpHashValue);
			}
			else if (errThreshold == 4) /* Skipped */
			{
		  		err = verifySingleEndEditDistance4(genLoc, _tmpSeq,
							 leftSeqLength, _tmpSeq + a, rightSeqLength,
							 middleSeqLength, matrix, &map_location, _tmpHashValue);
			}
			else if (errThreshold == 3) /* OK */
			{
		  		err = verifySingleEndEditDistance(genLoc, _tmpSeq,
							leftSeqLength, _tmpSeq + a, rightSeqLength,
							middleSeqLength, matrix, &map_location, _tmpHashValue);
			}
			else				/* OK */
			{
		  		err = verifySingleEndEditDistanceExtension(genLoc, _tmpSeq,
								 leftSeqLength, _tmpSeq + a, rightSeqLength,
								 middleSeqLength, matrix, &map_location, _tmpHashValue);
			}
		} 


		/* Update verified locations after an alignment for this k-mer. */
		for (j = -errThreshold+1; j < errThreshold; j++)
		{
			if(genLoc-(readSegment*WINDOW_SIZE)+j >= _msf_refGenBeg &&
		   		genLoc-(readSegment*WINDOW_SIZE)+j <= _msf_refGenEnd)
			{
		  		_msf_verifiedLocs[genLoc-(readSegment*WINDOW_SIZE)+j] = readId;
			}
		}
		

		/* Write mapping information if we could successfully map this read to
		 * the currently inspected location. */
		if (err != -1)
		{
			/* Generate SAM output information using the matrix containing
			 * matching information and store the output in editString. */
			generateSNPSAM(matrix, strlen(matrix), editString);

			/* Generate CIGAR output in a similar manner. */
			generateCigar(matrix, strlen(matrix), cigar);

			if (!bestMode)
			{
				/* Increment mapping count. This is for each k-mer, about how
				 * many times it was mapped in the reference genome. */
				mappingCnt++;

				/* Output information. */
				_msf_seqList[readNumber].hits[0]++;
				_msf_output.QNAME  = _msf_seqList[readNumber].name;
				_msf_output.FLAG   = 16 * direction;
				_msf_output.RNAME  = _msf_refGenName;
				_msf_output.POS	   = map_location + _msf_refGenOffset;
				_msf_output.MAPQ   = 255;
				_msf_output.CIGAR  = cigar;
				_msf_output.MRNAME = "*";
				_msf_output.MPOS   = 0;
				_msf_output.ISIZE  = 0;
				_msf_output.SEQ	   = _tmpSeq;
				_msf_output.QUAL   = _tmpQual;

				_msf_output.optSize	  = 2;
				_msf_output.optFields = _msf_optionalFields;

				_msf_optionalFields[0].tag	= "NM";
				_msf_optionalFields[0].type = 'i';
				_msf_optionalFields[0].iVal = err;

				_msf_optionalFields[1].tag	= "MD";
				_msf_optionalFields[1].type = 'Z';
				_msf_optionalFields[1].sVal = editString;

				/* Write the output to file. Note that _msf_output is of SAM
				 * type. */
				output(_msf_output);

				/* Increment global mappedSeqCnt if this read has been
				 * mappped. Note that hits is incremented every time we could
				 * map a k-mer to the ref genome but this is incremented only
				 * once by checking if it is equal to 1. */
				if ( _msf_seqList[readNumber].hits[0] == 1 )
				  	mappedSeqCnt++;

				/* purpose of maxHits - skip  */
				if ( maxHits == 0 )
				  	_msf_seqList[readNumber].hits[0] = 2;

				if ( maxHits != 0 &&
					 _msf_seqList[readNumber].hits[0] == maxHits )
				{
				  	completedSeqCnt++;
				  	break;
				}
			}
			else				/* best mode, this does not immediately dump the
								 * mapping information to the file. */
			{
				/* Increment mapping count. */
				mappingCnt++;

				/* Update hits. */
				_msf_seqList[readNumber].hits[0]++;

				/* Increment number of reads mapped. */
				if ( _msf_seqList[readNumber].hits[0] == 1 )
				  	mappedSeqCnt++;

				/* Need to resolve what maxHits is for */
				if ( maxHits == 0 )
				  	_msf_seqList[readNumber].hits[0] = 2;

				/* Update probability of this read, dont get what this
				 * probability is about. */
				if ( seqFastq )
				  	bestHitMappingInfo[readNumber].tprob +=
						mapProb(readNumber, editString, direction, err);

				/* Update full mapping information for the best mapping for the
				 * current read. Note that this _not_ performed only for
				 * once. Each time a best mapping is found, this information
				 * overwrites the previous one - for the current read if the
				 * editDistance of the current mapping is better than the
				 * previous ones OR it is the first mapping information written
				 * down. */
				if ( err  < bestHitMappingInfo[readNumber].err ||
					 bestHitMappingInfo[readNumber].loc == -1 )
					setFullMappingInfo(readNumber, map_location + _msf_refGenOffset,
									   direction, err, 0, editString, _msf_refGenName, cigar);
			}
		}
	}
}



/*------------------------------------------------------------------------------
 * This is for sorting the key_struct with respect to its key_entry_size field,
 * which basically is the number of locations this k-mer of the read is mapped
 * to.

 * OLD COMMENTS (not @OGUZ-COMMENT)
 *	MrFAST with fastHASH: compareEntrySize()
 *----------------------------------------------------------------------------*/
int
compareEntrySize ( const void *a, const void *b )
{
	return ((*(key_struct*) a).key_entry_size
		- (*(key_struct*) b).key_entry_size);
}



/*------------------------------------------------------------------------------
 * Maps all reads. For both forward and backward direction.

 * OLD COMMENTS (not @OGUZ-COMMENT)
 *	MrFAST with fastHASH: mapAllSingleEndSeq()
 *----------------------------------------------------------------------------*/
int
mapAllSingleEndSeq (   )
{
	int				 i			= 0;
	int				 j			= 0;
	int				 k			= 0;
	int				 it			= 0;
	unsigned int	*locs		= NULL;

	/* key_number is something related to sampling locations. For instance if
	 * the window size is 13 and read length is 50, then there are 3
	 * non-overlapping portions of the read we can start mapping to candidate
	 * locations in the ref genome. */
	int				 key_number = SEQ_LENGTH / WINDOW_SIZE;
	key_struct*		 sort_input = getMem(key_number * sizeof(key_struct));
	
	// Forward Mode
	/* Process reads in sorted order wrt their hash values. This ensures spatial
	 * coherency in the hash table. */
	for (i = 0; i < _msf_seqListSize; i++)
	{
		/* Get the read number to access the original read list. */
		k = _msf_sort_seqList[i].readNumber;
		int available_key_num = 0;

		/* Get locations in the ref genome for each sampling location. */
		for (it = 0; it < key_number; it++)
		{
			locs = getCandidates(
					 hashVal(_msf_seqList[k].seq + it * WINDOW_SIZE));
			if (locs != NULL)
			{
		  		sort_input[available_key_num].key_number	 = it;
		  		sort_input[available_key_num].key_entry		 = locs;
		  		sort_input[available_key_num].key_entry_size = locs[0];
		  		available_key_num++;
			}
		}

		/* available_key_num can ben smaller than the key_number if some locs =
		 * NULL. So we set the number of operation_key_num to this value found
		 * just above. */
		int operating_key_num = _msf_samplingLocsSize;
		if (available_key_num < operating_key_num)
		{
		  	operating_key_num = available_key_num;
		}

		/* Sort the available k-mers (whose locations are not NULL) with respect
		 * to their locs size. */
		qsort(sort_input, available_key_num, sizeof(key_struct),
		  compareEntrySize);
		
		for (j = 0; j < operating_key_num; j++)
		{			
			/* Determine the current sampling location of the read */
			_msf_samplingLocs[j] = sort_input[j].key_number * WINDOW_SIZE;

			/* Maps a single k-mer to all possible locations. */
			mapSingleEndSeq(sort_input[j].key_entry + 1,
					sort_input[j].key_entry_size, k, sort_input[j].key_number,
					0, j, sort_input, available_key_num);
		}
	}

	
	/* Do all of the above for the reverse sequence. */
	
	/* Process reads in sorted order wrt their hash values. This ensures spatial
	 * coherency in the hash table. */
	for (i = 0; i < _msf_seqListSize; i++)
	{
		/* Get the read number to access the original read list. */
		k					  = _msf_sort_seqList[i].readNumber;
		int available_key_num = 0;

		/* Get locations in the ref genome for each sampling location. */
		for (it = 0; it < key_number; it++)
		{
			locs = getCandidates(
					 hashVal(_msf_seqList[k].rseq + it * WINDOW_SIZE));
			if (locs != NULL)
			{
		  		sort_input[available_key_num].key_number	 = it;
		  		sort_input[available_key_num].key_entry		 = locs;
		  		sort_input[available_key_num].key_entry_size = locs[0];
		  		available_key_num++;
			}
		}		

		/* Sort the available k-mers (whose locations are not NULL) with respect
		 * to their locs size. */
		qsort(sort_input, available_key_num, sizeof(key_struct),
		  compareEntrySize);

		/* available_key_num can ben smaller than the key_number if some locs =
		 * NULL. So we set the number of operation_key_num to this value found
		 * just above. */
		int operating_key_num = _msf_samplingLocsSize;		
		
		if (available_key_num < operating_key_num)
		{
		  	operating_key_num = available_key_num;
		}

		for (j = 0; j < operating_key_num; j++)
		{
			/* Determine the current sampling location of the read */
			_msf_samplingLocs[j] = sort_input[j].key_number * WINDOW_SIZE;

			/* Maps a single k-mer to all possible locations. */
			mapSingleEndSeq(sort_input[j].key_entry + 1,
					sort_input[j].key_entry_size, k, sort_input[j].key_number,
					1, j, sort_input, available_key_num);
		}
	}

	/* Cleanup. */
	freeMem(sort_input, key_number * sizeof(key_struct));
	
	return 1;
}



/*------------------------------------------------------------------------------
 * Used for sorting FullMappingInfo structs with respect to their locations.
 *----------------------------------------------------------------------------*/
int
compareOut ( const void *a, const void *b )
{
	FullMappingInfo *aInfo = (FullMappingInfo *) a;
	FullMappingInfo *bInfo = (FullMappingInfo *) b;
	return aInfo->loc - bInfo->loc;
}




/*------------------------------------------------------------------------------
 * Parameters:
 *	l1					: The locations of this k-mer in the reference genome.
 *	s1					: Number of such locations.
 *	readNumber			: The id of this read in the global read seq array.
 *	readSegment			: The id of the current segment in this specific read
 						  (note that they were sorted so that this does not have
						  to be in order).
 *	direction			: 0 = FORWARD, 1 = BACKWARD
 *	index				: Which readSegment we are processing? (in order)
 *	keys_input			: All keys (readSegment or k-mer) of this read.
 *	potential_key_num	: Number of read-segments to be processed.

 * The function signature is the same with mapSingleEndSeq

 * OLD COMMENTS (not @OGUZ-COMMENT)
 *	direction = 0 forward
 				1 backward
 * MrFAST with fastHASH: mapPairEndSeqList()

 * Om - On the Mountain At Dawn
 *----------------------------------------------------------------------------*/
void
mapPairEndSeqList ( unsigned int *l1, int s1, int readNumber, int readSegment,
					int direction, int index, key_struct* keys_input,
					int potential_key_number)
{
	int		 z	  = 0;
	int		*locs = (int *) l1;
	char	*_tmpSeq;

	/* matrix variable stores the matching information.  */
	char	matrix[SEQ_MAX_LENGTH];
	char	editString[2 * SEQ_MAX_LENGTH];
	char	cigar[MAX_CIGAR_SIZE];

	short *_tmpHashValue;

	/* Left, middle and right segment lengths. */
	int leftSeqLength	= 0;
	int rightSeqLength	= 0;
	int middleSeqLength = 0;

	int genLoc = 0;
	int r	   = readNumber;

	/* if d = -1 it is BACKWARD, if it is then FORWARD. */
	char d = (direction == 1) ? -1 : 1;

	/* readId is for _msf_verifiedLocs. We need to store two entries per read
	 * (FORWARD and BACKWARD), so x2. */
	int readId	   = 2 * readNumber + direction;
	int key_number = SEQ_LENGTH / WINDOW_SIZE;

	/* With respect to chosen direction, either choose seq or rseq in the
	 * original read */
	if (d == -1)				/* -1 = BACKWARD */
	{
		_tmpSeq = _msf_seqList[readNumber].rseq;
		_tmpHashValue = _msf_seqList[readNumber].rhashValue;
	} else						/* 1 = FORWARD */
	{
		_tmpSeq = _msf_seqList[readNumber].seq;
		_tmpHashValue = _msf_seqList[readNumber].hashValue;
	}

	/* For each location this k-mer can be mapped to. */
	for (z = 0; z < s1; z++)
	{
		int map_location = 0;
		int a			 = 0;
		int o			 = index;

		/* Current location of that read-segment will be mapped to in reference
		 * genome. */
		genLoc			 = locs[z];

		/* We have three pieces in the read sequence - left, middle and right
		 * sequence. */
		leftSeqLength	= _msf_samplingLocs[o];
		middleSeqLength = WINDOW_SIZE;
		a				= leftSeqLength + middleSeqLength;
		rightSeqLength	= SEQ_LENGTH - a;

		/*    
		//hxinPE: If the read is at the beginning of the contig, and there are insertions
		//hxinPE: then genLoc - _msf_sampleingLoc[0] might be smaller than _refGenBeg
		int realLoc = genLoc - _msf_samplingLocs[o];
		if (genLoc < _msf_samplingLocs[o]
		|| genLoc - _msf_samplingLocs[o] < _msf_refGenBeg
		|| genLoc - _msf_samplingLocs[o] > _msf_refGenEnd) {
		  if (genLoc - _msf_samplingLocs[o] > _msf_refGenBeg - errThreshold)
		realLoc =  _msf_refGenBeg;
		  else
		continue;
		}

		if (_msf_verifiedLocs[realLoc] == readId || _msf_verifiedLocs[realLoc] == -readId )
		  continue;
		*/

		/* Some corners cases about the beginning and end of the reference
		 * genome. Different from normal mode mapping, the verified locations
		 * are checked here and skipped this read has been processed before. */
		if (genLoc - leftSeqLength < _msf_refGenBeg ||
			genLoc + rightSeqLength + middleSeqLength > _msf_refGenEnd ||
			_msf_verifiedLocs[genLoc - _msf_samplingLocs[o]] == readId ||
			_msf_verifiedLocs[genLoc - _msf_samplingLocs[o]] == -readId)
		{
		  	continue;
		}

		//Begin of long-K
		/* Tries to merge read-segments? It skips the current locs entry
		 * inspected if something is validated (Idx = index). */
		int mergeIdx = 2 * errThreshold + 1 - index;
		
		if (mergeIdx < potential_key_number)
		{
			if (!searchKey(
				   genLoc
				   + (keys_input[mergeIdx].key_number
				  - keys_input[o].key_number) * WINDOW_SIZE,
				   keys_input[mergeIdx].key_entry,
				   keys_input[mergeIdx].key_entry_size))
			{
		  		continue;
			}
		}
		//End of long-K
		
		// Adjacency Filtering Start ---------------------------------
		int skip_edit_distance = 0;
		int diff_num		   = 0;
		int ix				   = 0;
		for (ix = 0; ix < potential_key_number; ix++)
		{
			if (ix >= key_number - errThreshold)
			{
		  		break;
			}
			
			if (ix != o && ix != mergeIdx)
			{ // Changed with long-K
				if (!searchKey(
						   genLoc
						   + (keys_input[ix].key_number
						  - keys_input[o].key_number)
						   * WINDOW_SIZE, keys_input[ix].key_entry,
						   keys_input[ix].key_entry_size))
				{
				  	diff_num++;
				  	if (diff_num > errThreshold)
					{
						skip_edit_distance = 1;
						break;
				  	}
				}
			}
		}
		// Adjacency Filtering End -----------------------------------
		
		int err		 = -1;
		map_location = 0;

		/* This is where the matching is done. */
		if (skip_edit_distance == 0)
		{			
			if (errThreshold == 2) /* OK */
			{
				err = verifySingleEndEditDistance2(genLoc, _tmpSeq,
								   leftSeqLength, _tmpSeq + a, rightSeqLength,
								   middleSeqLength, matrix, &map_location, _tmpHashValue);
			}
			else if (errThreshold == 4) /* Skipped */
			{
				err = verifySingleEndEditDistance4(genLoc, _tmpSeq,
								   leftSeqLength, _tmpSeq + a, rightSeqLength,
								   middleSeqLength, matrix, &map_location, _tmpHashValue);
			}
			else if (errThreshold == 3) /* OK */
			{
				err = verifySingleEndEditDistance(genLoc, _tmpSeq,
								  leftSeqLength, _tmpSeq + a, rightSeqLength,
								  middleSeqLength, matrix, &map_location, _tmpHashValue);
			}
			else				/* OK */
			{
				err = verifySingleEndEditDistanceExtension(genLoc, _tmpSeq,
									   leftSeqLength, _tmpSeq + a, rightSeqLength,
									   middleSeqLength, matrix, &map_location, _tmpHashValue);
			}
		}
		else
		{
		  	err = -1;
		}


		/* Update verified locations after an alignment for this k-mer. */
		int j = 0;
		for (j = -errThreshold+1; j < errThreshold; j++)
		{
			if(genLoc-(readSegment*WINDOW_SIZE)+j >= _msf_refGenBeg &&
			   genLoc-(readSegment*WINDOW_SIZE)+j <= _msf_refGenEnd)
		  		_msf_verifiedLocs[genLoc-(readSegment*WINDOW_SIZE)+j] = readId;
		}


		/* Output information, this is different from normal mapping mode. */
		if (err != -1)
		{
			int i = 0;

			/* calkan counter */
			mappingCnt++;

			/* Generate SAM output information using the matrix containing
			 * matching information and store the output in editString. */
			generateSNPSAM(matrix, strlen(matrix), editString);

			/* Generate CIGAR output in a similar manner. */
			generateCigar(matrix, strlen(matrix), cigar);

			/* r is the current readId. */
			MappingLocations	*parent = NULL;
			MappingLocations	*child	= _msf_mappingInfo[r].next;
			genLoc						= map_location + _msf_refGenOffset;

			/* A single MappingLocations stores MAP_CHUNKS number of
			 * results. This finds the current MappingLocations that will be
			 * used to store latest mapping result. */
			for (i = 0; i < (_msf_mappingInfo[r].size / MAP_CHUNKS); i++)
			{
				parent = child;
				child  = child->next;
			}

			if (child == NULL)	/* A new MappingLocations struct required. */
			{
				MappingLocations		*tmp = getMem(sizeof(MappingLocations));
				tmp->next					 = NULL;

				/* This is the first entry in the newly created
				 * MappingLocations, so it is 0th element. Recall that d was
				 * either 1 or -1. Purpose not clear at the moment, but it is
				 * probably for the clarification of direction during output. */
				tmp->loc[0]					 = genLoc *	 d;	// d is required: DHL
				tmp->err[0]					 = err;
				tmp->cigarSize[0]			 = strlen(cigar);
				sprintf(tmp->cigar[0], "%s", cigar);
				tmp->mdSize[0] = strlen(editString);
				sprintf(tmp->md[0], "%s", editString);

				if (parent == NULL)
				  	_msf_mappingInfo[r].next = tmp;
				else
				  	parent->next = tmp;
			}
			else				/* Utilize an already allocated MappingLocations
								 * struct. */
			{
				/* Error checking. */
				if (strlen(cigar) > SEQ_LENGTH ||
					strlen(editString) > SEQ_LENGTH)
				{
					printf(
					   "ERROR in %d read size(After mapping) exceeds cigar=%d md =%d cigar=%s md =%s\n",
					   r, (int) strlen(cigar), (int) strlen(editString),
					   cigar, editString);
				}

				/* Write output information. Multiplication of genLoc with d
				 * probably related with clarification of direction during
				 * output. */
				child->loc[_msf_mappingInfo[r].size % MAP_CHUNKS] = genLoc * d;
				child->err[_msf_mappingInfo[r].size % MAP_CHUNKS] = err;
				child->cigarSize[_msf_mappingInfo[r].size % MAP_CHUNKS] = strlen(cigar);
				sprintf(child->cigar[_msf_mappingInfo[r].size % MAP_CHUNKS], "%s", cigar);
				child->mdSize[_msf_mappingInfo[r].size % MAP_CHUNKS] = strlen(editString);
				sprintf(child->md[_msf_mappingInfo[r].size % MAP_CHUNKS], "%s", editString);
			}

			/* Note that the size is not equal to MappingLocations structs but
			 * its multiplication with MAP_CHUNKS, i.e., the number of
			 * successful mappings. */
			_msf_mappingInfo[r].size++;
		}
	}
}



/*------------------------------------------------------------------------------
 * Maps all reads. For both forward and backward direction. It performs
 * additional processing related to paired-end mode compared to normal
 * mode. Actually, the initial code where mapping is done in forward and
 * backward is the same with mapAllSingleEndSeq. But the rest differs since
 * additional processing done later here.

 * OLD COMMENTS (not @OGUZ-COMMENT)
 *	MrFAST with fastHASH: mapPairedEndSeq()
 *----------------------------------------------------------------------------*/
void
mapPairedEndSeq (   )
{
	// DHL: Changed Start
	int				 i			= 0;
	int				 j			= 0;
	int				 k			= 0;
	unsigned int	*locs		= NULL;

	/* key_number is something related to sampling locations. For instance if
	 * the window size is 13 and read length is 50, then there are 3
	 * non-overlapping portions of the read we can start mapping to candidate
	 * locations in the ref genome. */
	int				 key_number = SEQ_LENGTH / WINDOW_SIZE;
	key_struct*		 sort_input = getMem(key_number * sizeof(key_struct));

	// Forward Mode
	/* Process reads in sorted order wrt their hash values. This ensures spatial
	 * coherency in the hash table. */
	for (i = 0; i < _msf_seqListSize; i++)
	{
		/* Get the read number to access the original read list. */
		k					  = _msf_sort_seqList[i].readNumber;
		int available_key_num = 0;

		/* Get locations in the ref genome for each sampling location. */
		int it				  = 0;
		for (it = 0; it < key_number; it++)
		{
			int key_hash = hashVal(_msf_seqList[k].seq + it * WINDOW_SIZE);
			locs		 = getCandidates(key_hash);
			if (locs != NULL)
			{
		  		sort_input[available_key_num].key_number	 = it;
		  		sort_input[available_key_num].key_entry		 = locs;
		  		sort_input[available_key_num].key_entry_size = locs[0];
		  		available_key_num++;
			}
		}

		/* available_key_num can ben smaller than the key_number if some locs =
		 * NULL. So we set the number of operation_key_num to this value found
		 * just above. */
		int operating_key_num = _msf_samplingLocsSize;
		if (available_key_num < operating_key_num)
		{
		  operating_key_num = available_key_num;
		}

		/* Sort the available k-mers (whose locations are not NULL) with respect
		 * to their locs size. */
		qsort(sort_input, available_key_num, sizeof(key_struct),
		  compareEntrySize);
		
		for (j = 0; j < operating_key_num; j++)
		{			
			/* Determine the current sampling location of the read */
			_msf_samplingLocs[j] = sort_input[j].key_number * WINDOW_SIZE;

			/* Maps a single k-mer to all possible locations. */
			mapPairEndSeqList(sort_input[j].key_entry + 1,
				  sort_input[j].key_entry_size, k, sort_input[j].key_number,
				  0, j, sort_input, available_key_num);
		}
	}

	/* Do all of the above for the reverse sequence. */	
	// Reverse Mode
	/* Process reads in sorted order wrt their hash values. This ensures spatial
	 * coherency in the hash table. */
	for (i = 0; i < _msf_seqListSize; i++)
	{
		/* Get the read number to access the original read list. */
		k					  = _msf_sort_seqList[i].readNumber;
		int key_number		  = SEQ_LENGTH / WINDOW_SIZE;
		int available_key_num = 0;

		/* Get locations in the ref genome for each sampling location. */
		int it				  = 0;
		for (it = 0; it < key_number; it++)
		{
			int key_hash = hashVal(_msf_seqList[k].rseq + it * WINDOW_SIZE);
			locs = getCandidates(key_hash);
			if (locs != NULL)
			{
		  		sort_input[available_key_num].key_number	 = it;
		  		sort_input[available_key_num].key_entry		 = locs;
		  		sort_input[available_key_num].key_entry_size = locs[0];
		  		available_key_num++;
			}
		}

		/* available_key_num can ben smaller than the key_number if some locs =
		 * NULL. So we set the number of operation_key_num to this value found
		 * just above. */
		int operating_key_num = _msf_samplingLocsSize;
		if (available_key_num < operating_key_num)
		{
		  	operating_key_num = available_key_num;
		}

		/* Sort the available k-mers (whose locations are not NULL) with respect
		 * to their locs size. */
		qsort(sort_input, available_key_num, sizeof(key_struct),
		  compareEntrySize);

		for (j = 0; j < operating_key_num; j++)
		{
			/* Determine the current sampling location of the read */
		  	_msf_samplingLocs[j] = sort_input[j].key_number * WINDOW_SIZE;

			/* Maps a single k-mer to all possible locations. */
		  	mapPairEndSeqList(sort_input[j].key_entry + 1,
					sort_input[j].key_entry_size, k, sort_input[j].key_number,
					1, j, sort_input, available_key_num);
		}
	}

	/* Cleanup. */
	freeMem(sort_input, key_number * sizeof(key_struct));
	// DHL: Changed End


	/* ********************************************************************** */
	/* ********************************************************************** */
	/* ********************************************************************** */

	
	/* This is where it differs from normal mode mapping. */
	char				 fname1[FILE_NAME_LENGTH];
	char				 fname2[FILE_NAME_LENGTH];
	MappingLocations	*cur;
	int					 tmpOut;
	int					 lmax = 0, rmax = 0;

	sprintf(fname1, "%s__%s__%s__%d__1.tmp", mappingOutputPath, _msf_refGenName,
		mappingOutput, _msf_openFiles);
	sprintf(fname2, "%s__%s__%s__%d__2.tmp", mappingOutputPath, _msf_refGenName,
		mappingOutput, _msf_openFiles);

	FILE* out;
	FILE* out1 = fileOpen(fname1, "w");
	FILE* out2 = fileOpen(fname2, "w");

	_msf_openFiles++;

	for (i = 0; i < _msf_seqListSize; i++)
	{
		/* Recall that in paired-end mode, read[i] and read[i+1] store a read's
		 * normal and complement sequence. Find the max number of left and right
		 * mappings among all reads. Also determine to which file to write. It
		 * can be either one of the ...__1 or ..._2 files. These are binary
		 * files.  */
		if (i % 2 == 0)
		{
			out = out1;
			if (lmax < _msf_mappingInfo[i].size)
		  		lmax = _msf_mappingInfo[i].size;
		}
		else
		{
			out = out2;
			if (rmax < _msf_mappingInfo[i].size)
		  		rmax = _msf_mappingInfo[i].size;
		}

		/* Output #mappings for the current read.  */
		tmpOut = fwrite(&(_msf_mappingInfo[i].size), sizeof(int), 1, out);
		if (_msf_mappingInfo[i].size > 0)
		{
			/* Find the MappingLocations struct to write. */
			cur = _msf_mappingInfo[i].next;
			for (j = 0; j < _msf_mappingInfo[i].size; j++)
			{
				
				if (j > 0 && j % MAP_CHUNKS == 0)
				{
				  	cur = cur->next;
				}

				if(debugMode &&
				   (cur->cigarSize[j % MAP_CHUNKS] > SEQ_LENGTH ||
					cur->mdSize[j % MAP_CHUNKS] > SEQ_LENGTH))
				{
					printf("ERROR in %d read size exceeds cigar=%d md =%d cigar=%s md =%s\n", i,  cur->cigarSize[j % MAP_CHUNKS], cur->mdSize[j % MAP_CHUNKS], cur->cigar[j % MAP_CHUNKS], cur->md[j % MAP_CHUNKS]);	
				}

				/* Dump the CIGAR and MD output along with other details. */
				tmpOut = fwrite(&(cur->loc[j % MAP_CHUNKS]),
								sizeof(int), 1, out);
				tmpOut = fwrite(&(cur->err[j % MAP_CHUNKS]),
								sizeof(int), 1, out);
				tmpOut = fwrite(&(cur->cigarSize[j % MAP_CHUNKS]),
								sizeof(int), 1, out);
				tmpOut = fwrite((cur->cigar[j % MAP_CHUNKS]),
								sizeof(char), (cur->cigarSize[j % MAP_CHUNKS]),
								out);
				tmpOut = fwrite(&(cur->mdSize[j % MAP_CHUNKS]),
								sizeof(int), 1, out);
				tmpOut = fwrite((cur->md[j % MAP_CHUNKS]),
								sizeof(char), (cur->mdSize[j % MAP_CHUNKS]),
								out);
			}

			/* It will futher be used. */
			_msf_mappingInfo[i].size = 0;
		}
	}

	/* Update global left and right max #mappings. */
	_msf_maxLSize += lmax;
	_msf_maxRSize += rmax;
	tmpOut++;					/* What is going on here? */

	fclose(out1);
	fclose(out2);
}



/*------------------------------------------------------------------------------
 * Outputs a read pair information, read/1 and read/2. This is only called from
 * finalizeBestConcordantDiscordant
 *----------------------------------------------------------------------------*/
void
outputPairFullMappingInfo ( FILE *fp, int readNumber )
{
	char *seq1, *seq2, *rseq1, *rseq2, *qual1, *qual2;
	char rqual1[SEQ_LENGTH + 1], rqual2[SEQ_LENGTH + 1];

	rqual1[SEQ_LENGTH] = rqual2[SEQ_LENGTH] = '\0';

	seq1  = _msf_seqList[readNumber * 2].seq;
	rseq1 = _msf_seqList[readNumber * 2].rseq;
	qual1 = _msf_seqList[readNumber * 2].qual;

	reverse(_msf_seqList[readNumber * 2].qual, rqual1, SEQ_LENGTH);

	seq2  = _msf_seqList[readNumber * 2 + 1].seq;
	rseq2 = _msf_seqList[readNumber * 2 + 1].rseq;
	qual2 = _msf_seqList[readNumber * 2 + 1].qual;

	reverse(_msf_seqList[readNumber * 2 + 1].qual, rqual2, SEQ_LENGTH);

	/* None of the /1 /2 of the read could not be mapped, so there should not be
	 * concordant/discordant mapping. Thus return. */
	if (bestHitMappingInfo[readNumber * 2].loc == -1 &&
		bestHitMappingInfo[readNumber * 2 + 1].loc == -1)
	  	return;
	else
	{
		char	*seq;
		char	*qual;
		char	 d1;
		char	 d2;
		int		 isize;
		int		 proper = 0;
		
		// ISIZE CALCULATION
		// The distance between outer edges
		isize = abs(
			bestHitMappingInfo[readNumber * 2].loc -
			bestHitMappingInfo[readNumber * 2 + 1].loc
			)
			+ SEQ_LENGTH - 2;

		if (bestHitMappingInfo[readNumber * 2].loc -
			bestHitMappingInfo[readNumber * 2 + 1].loc > 0)
		{
		  	isize *= -1;
		}

		/* Direction */
		d1 = (bestHitMappingInfo[readNumber * 2].dir == -1) ? 1 : 0;
		d2 = (bestHitMappingInfo[readNumber * 2 + 1].dir == -1) ? 1 : 0;
		if (d1)
		{
		  	seq = rseq1;
		  	qual = rqual1;
		}
		else
		{
		  	seq = seq1;
		  	qual = qual1;
		}

		if (
			(bestHitMappingInfo[readNumber * 2].loc <
			 bestHitMappingInfo[readNumber * 2 + 1].loc && !d1 && d2) ||
			(bestHitMappingInfo[readNumber * 2].loc >
			 bestHitMappingInfo[readNumber * 2 + 1].loc && d1
			&& !d2)
			)
		{
		  	proper = 2;
		}
		else
		{
		  	proper = 0;
		}

		/* proper variable is used in FLAG field. */
		_msf_output.POS	  = bestHitMappingInfo[readNumber * 2].loc;
		_msf_output.MPOS  = bestHitMappingInfo[readNumber * 2 + 1].loc;
		_msf_output.FLAG  = 1 + proper + 16 * d1 + 32 * d2 + 64;
		_msf_output.ISIZE = isize;
		_msf_output.SEQ	  = seq;
		_msf_output.QUAL  = qual;
		_msf_output.QNAME = _msf_seqList[readNumber * 2].name;
		_msf_output.RNAME = bestHitMappingInfo[readNumber * 2].chr;
		
		if (seqFastq)
		  	_msf_output.MAPQ = mapQ(readNumber * 2) + mapQ(readNumber * 2 + 1);
		else
		  	_msf_output.MAPQ = 255;
		
		_msf_output.CIGAR  = bestHitMappingInfo[readNumber * 2].cigar;
		_msf_output.MRNAME = "=";

		_msf_output.optSize	  = 2;
		_msf_output.optFields = _msf_optionalFields;

		_msf_optionalFields[0].tag	= "NM";
		_msf_optionalFields[0].type = 'i';
		_msf_optionalFields[0].iVal = bestHitMappingInfo[readNumber * 2].err;

		_msf_optionalFields[1].tag	= "MD";
		_msf_optionalFields[1].type = 'Z';
		_msf_optionalFields[1].sVal = bestHitMappingInfo[readNumber * 2].md;

		/* Output read/1 */
		output(_msf_output);


		/* Now do the same stuff for /2 */
		if (d2)
		{
		  	seq	 = rseq2;
		  	qual = rqual2;
		}
		else
		{
		  	seq	 = seq2;
		  	qual = qual2;
		}

		_msf_output.POS	  = bestHitMappingInfo[readNumber * 2 + 1].loc;
		_msf_output.MPOS  = bestHitMappingInfo[readNumber * 2].loc;
		_msf_output.FLAG  = 1 + proper + 16 * d2 + 32 * d1 + 128;
		_msf_output.ISIZE = -isize;
		_msf_output.SEQ	  = seq;
		_msf_output.QUAL  = qual;
		_msf_output.QNAME = _msf_seqList[readNumber * 2].name;
		_msf_output.RNAME = bestHitMappingInfo[readNumber * 2].chr;
		
		if (seqFastq)
			_msf_output.MAPQ = mapQ(readNumber * 2) + mapQ(readNumber * 2 + 1);
		else
		  	_msf_output.MAPQ = 255;

		_msf_output.CIGAR  = bestHitMappingInfo[readNumber * 2 + 1].cigar;
		_msf_output.MRNAME = "=";

		_msf_output.optSize	  = 2;
		_msf_output.optFields = _msf_optionalFields;

		_msf_optionalFields[0].tag	= "NM";
		_msf_optionalFields[0].type = 'i';
		_msf_optionalFields[0].iVal = bestHitMappingInfo[readNumber * 2 + 1].err;

		_msf_optionalFields[1].tag	= "MD";
		_msf_optionalFields[1].type = 'Z';
		_msf_optionalFields[1].sVal = bestHitMappingInfo[readNumber * 2 + 1].md;

		/* Output read/2 */
		output(_msf_output);
	}
}



/*------------------------------------------------------------------------------
 * -SKIPPED-
 *----------------------------------------------------------------------------*/
/*
  Find the closet one to the c
  @return 0: if the x1 is closer to c
  1: if the x2 is closer to c
  2: if both distance are equal
  -1: if error
*/
int
findNearest ( int x1, int x2, int c )
{

  	if (abs(x1 - c) > abs(x2 - c))
    	return 0;
  	else if (abs(x1 - c) < abs(x2 - c))
    	return 1;
  	else if (abs(x1 - c) == abs(x2 - c))
    	return 2;
  	else
    	return -1;
}



/*------------------------------------------------------------------------------
 * Writes the best concordant or best discordant information. This was stored
 * in bestHitMappingInfo.
 *----------------------------------------------------------------------------*/
void
finalizeBestConcordantDiscordant (   )
{
	int i = 0;
	for (i = 0; i < _msf_seqListSize / 2; i++)
	{
	  	outputPairFullMappingInfo(NULL, i);
	}

	/* bestHitMappingInfo is released here. This should actually be
	 * sizeof(BestFullMappingInfo), it is allocated in that
	 * way. sizeof(BestFullMappingInfo) > sizeof(FullMappingInfo)
	 * anyways. Possible leak but who cares. */
	freeMem(bestHitMappingInfo, _msf_seqListSize * sizeof(FullMappingInfo));
}



/*------------------------------------------------------------------------------
 * Not sure but computes a probability value related with the read using its
 * mapping information.
 *----------------------------------------------------------------------------*/
double
mapProb ( int readNumber, char *md, int dir, int err )
{
	int		i	  = 0;
	int		mdlen = strlen(md);
	char	buf[MAX_CIGAR_SIZE];
	int		j	  = 0;

	double	phred  = 0.0;
	int		errloc = 0;
	int		errcnt = 0;			//since I cannot calculate deletion base quality

	buf[0] = 0;

	if (err == 0) 
	  	return 1.0;

	while (i<mdlen)
	{
		if (isdigit(md[i]))
			buf[j++]=md[i++];
		
		else if (isalpha(md[i]))
		{
			/* mismatch */
			errcnt++;
			buf[j] = '\0'; 
			if (j != 0)
		  		errloc += atoi(buf);
			else if (i!=0)
		  		errloc++;

			j=0;
			buf[0]=0;

			if (dir)
		  		phred += (double) (_msf_seqList[readNumber].qual[SEQ_LENGTH-errloc-1] - 33);
			else
		  		phred += (double) (_msf_seqList[readNumber].qual[errloc] - 33);

			i++;
		}

		else if (md[i]=='^')
		{
			/* insertion to the read / deletion from reference  */
			if (j!=0)
			{
		  		buf[j]=0;
		  		errloc += atoi(buf);
		  		buf[0] = 0;
			}
			
			j=0; 
			i++; /* pass ^ */
			while (isalpha(md[i++]))
				j++;
			errloc += j;
			j = 0;
		}
	}

	double indel_prob = 1; 
	if (errcnt != err)
	  indel_prob = 0.0002 * (err - errcnt);

	return pow(10, -1 * (phred / 10)) * indel_prob;
}



/*------------------------------------------------------------------------------
 * -SKIPPED-
 *----------------------------------------------------------------------------*/
int
mapQ ( int readNumber )
{
  	int mapqual;
  	double mapprob;

  	mapprob = mapProb(readNumber, bestHitMappingInfo[readNumber].md, 
		    bestHitMappingInfo[readNumber].dir, bestHitMappingInfo[readNumber].err); 

  	if (mapprob == bestHitMappingInfo[readNumber].tprob)
    	mapqual = 40;
	else 
    	mapqual =  (int) (round(-10.0 * log10(1 - (mapprob / bestHitMappingInfo[readNumber].tprob))));
  
  	if (mapqual > 40) mapqual = 40;

  	return mapqual;
}



/*------------------------------------------------------------------------------
 * Sets attributes related to found best mapping.
 *----------------------------------------------------------------------------*/
void
setFullMappingInfo ( int readNumber, int loc, int dir, int err, int score,
					 char *md, char * refName, char *cigar )
{
	bestHitMappingInfo[readNumber].loc	 = loc;
	bestHitMappingInfo[readNumber].dir	 = dir;
	bestHitMappingInfo[readNumber].err	 = err;
	bestHitMappingInfo[readNumber].score = score;

	strncpy(bestHitMappingInfo[readNumber].md, md, strlen(md) + 1);

	/*
	if (bestHitMappingInfo[readNumber].chr == NULL)
	  bestHitMappingInfo[readNumber].chr = (char *) getMem(sizeof(char) * (strlen(refName)+1));
	else if (strlen(bestHitMappingInfo[readNumber].chr) < strlen(refName)){
	  freeMem(bestHitMappingInfo[readNumber].chr, (strlen(bestHitMappingInfo[readNumber].chr)+1));
	  bestHitMappingInfo[readNumber].chr = (char *) getMem(sizeof(char) * (strlen(refName)+1));
	}
	*/

	strncpy(bestHitMappingInfo[readNumber].chr, refName, strlen(refName) + 1);
	strncpy(bestHitMappingInfo[readNumber].cigar, cigar, strlen(cigar) + 1);
}



/*------------------------------------------------------------------------------
 * Sets attributes of the best paired-end mapping (seems to be for concordant
 * mappings + discordant mappings also). Performs it for both elements of the
 * pair. They are stored in 2*readNumber and 2*readNumber+1.
 *----------------------------------------------------------------------------*/
void
setPairFullMappingInfo ( int readNumber, FullMappingInfo mi1,
						 FullMappingInfo mi2 )
{
	bestHitMappingInfo[readNumber * 2].loc	 = mi1.loc;
	bestHitMappingInfo[readNumber * 2].dir	 = mi1.dir;
	bestHitMappingInfo[readNumber * 2].err	 = mi1.err;
	bestHitMappingInfo[readNumber * 2].score = mi1.score;

	/*
	if (bestHitMappingInfo[readNumber * 2].chr == NULL)
	  bestHitMappingInfo[readNumber * 2].chr = (char *) getMem(sizeof(char) * (strlen(_msf_refGenName)+1));
	else if (strlen(bestHitMappingInfo[readNumber * 2].chr) < strlen(_msf_refGenName)){
	  freeMem(bestHitMappingInfo[readNumber * 2].chr, (strlen(bestHitMappingInfo[readNumber * 2].chr)+1));
	  bestHitMappingInfo[readNumber * 2].chr = (char *) getMem(sizeof(char) * (strlen(_msf_refGenName)+1));
	  }
	*/

	snprintf(bestHitMappingInfo[readNumber * 2].chr,
			 strlen(_msf_refGenName)+1, "%s", _msf_refGenName);

	strncpy(bestHitMappingInfo[readNumber * 2].md, mi1.md, strlen(mi1.md) + 1);
	strncpy(bestHitMappingInfo[readNumber * 2].cigar, mi1.cigar,
		strlen(mi1.cigar) + 1);

	bestHitMappingInfo[readNumber * 2 + 1].loc	 = mi2.loc;
	bestHitMappingInfo[readNumber * 2 + 1].dir	 = mi2.dir;
	bestHitMappingInfo[readNumber * 2 + 1].err	 = mi2.err;
	bestHitMappingInfo[readNumber * 2 + 1].score = mi2.score;

	/*
	if (bestHitMappingInfo[readNumber * 2 + 1].chr == NULL)
	  bestHitMappingInfo[readNumber * 2 + 1].chr = (char *) getMem(sizeof(char) * (strlen(_msf_refGenName)+1));
	else if (strlen(bestHitMappingInfo[readNumber * 2 + 1].chr) < strlen(_msf_refGenName)){
	  freeMem(bestHitMappingInfo[readNumber * 2 + 1].chr, (strlen(bestHitMappingInfo[readNumber * 2 + 1].chr)+1));
	  bestHitMappingInfo[readNumber * 2 + 1].chr = (char *) getMem(sizeof(char) * (strlen(_msf_refGenName)+1));
	}
	*/

	snprintf(bestHitMappingInfo[readNumber * 2 + 1].chr,
			 strlen(_msf_refGenName)+1, "%s", _msf_refGenName);

	strncpy(bestHitMappingInfo[readNumber * 2 + 1].md, mi2.md,
		strlen(mi2.md) + 1);
	strncpy(bestHitMappingInfo[readNumber * 2 + 1].cigar, mi2.cigar,
		strlen(mi2.cigar) + 1);
}



/*------------------------------------------------------------------------------
 * For each read, it reads the previously written down FullMappingInfo data and
 * dumps OEA and DISCORDANT mapping information. It also outputs concordant
 * mapping data. The best concordant mapping information for a read (if not
 * concardant mapping can be found for a read, then best discordant mapping
 * information) is also stored in BestFullMappingInfo for later writing.

 * pairedEndMode
 * pairedEndDiscordantMode

 * If some location is < 0, it is reversed.
 *----------------------------------------------------------------------------*/
void
outputPairedEnd (   )
{
	int		i = 0;
	char	cigar[MAX_CIGAR_SIZE];
	int		tmpOut;

	FILE* in1[_msf_openFiles];
	FILE* in2[_msf_openFiles];

	char fname1[_msf_openFiles][FILE_NAME_LENGTH];
	char fname2[_msf_openFiles][FILE_NAME_LENGTH];

	// discordant
	FILE *out = NULL, *out1 = NULL;

	char fname3[FILE_NAME_LENGTH];
	char fname4[FILE_NAME_LENGTH];

	int meanDistanceMapping = 0;

	char	rqual1[SEQ_LENGTH + 1];
	char	rqual2[SEQ_LENGTH + 1];
	int		tmp = 0;

	/* Ref genome needed. */
	loadRefGenome(&_msf_refGen, &_msf_refGenName, &tmpOut);

	/* If in pairedEndDiscordantMode, need to write discordant and
	 * one-end-anchored mappings by appending them to the old files. */
	if (pairedEndDiscordantMode)
	{
		sprintf(fname3, "%s__%s__disc", mappingOutputPath, mappingOutput);
		sprintf(fname4, "%s__%s__oea", mappingOutputPath, mappingOutput);
		out	 = fileOpen(fname3, "a");
		out1 = fileOpen(fname4, "a");
	}

	/* Allocate FullMappingInfo for maximum number of mapping sizes. The max is
	 * here used obviously because a read can have at most that much sequence
	 * mappings. */
	FullMappingInfo *mi1 = getMem(sizeof(FullMappingInfo) * _msf_maxLSize);
	FullMappingInfo *mi2 = getMem(sizeof(FullMappingInfo) * _msf_maxRSize);

	/* Construct file names and open them to read. The convention is that the 1
	 * suffix is used for the first read of the paired-end read and the 2 suffix
	 * is used for the second read of the paired-end read. */
	_msf_fileCount[_msf_maxFile] = 0;
	for (i = 0; i < _msf_openFiles; i++)
	{
		sprintf(fname1[i], "%s__%s__%s__%d__1.tmp", mappingOutputPath,
			_msf_refGenName, mappingOutput, i);
		sprintf(_msf_fileName[_msf_maxFile][_msf_fileCount[_msf_maxFile]][0],
			"%s", fname1[i]);

		sprintf(fname2[i], "%s__%s__%s__%d__2.tmp", mappingOutputPath,
			_msf_refGenName, mappingOutput, i);
		sprintf(_msf_fileName[_msf_maxFile][_msf_fileCount[_msf_maxFile]][1],
			"%s", fname2[i]);

		in1[i] = fileOpen(fname1[i], "r");
		in2[i] = fileOpen(fname2[i], "r");
		_msf_fileCount[_msf_maxFile]++;
	}
	_msf_maxFile++;

	int size;
	int j, k;
	int size1, size2;

	/* meanDistanceMapping changes whether the mode is PE mode or PE discordant
	 * mode. */
	meanDistanceMapping = (pairedEndDiscordantMode == 1) ?
	  (minPairEndedDiscordantDistance + maxPairEndedDiscordantDistance) / 2
		+ SEQ_LENGTH :
	  (minPairEndedDistance + maxPairEndedDistance) / 2 + SEQ_LENGTH;

	/* Process each read two by two. Note that the paired-end reads are in
	 * consecutive locations in the array. */
	for (i = 0; i < _msf_seqListSize / 2; i++)
	{		
		/* Load information regarding 1st read of the current paired-end
		 * read. Load the mapping information into the FullMappingInfo mi1. */
		size1 = size2 = 0;
		for (j = 0; j < _msf_openFiles; j++)
		{
			/* size = number of locations this 1st pair mapped to. */
		  	tmpOut = fread(&size, sizeof(int), 1, in1[j]);
			
		  	if (size > 0)
		  	{
				for (k = 0; k < size; k++)
				{
					/* direction, location and edit distance. */
		  			mi1[size1 + k].dir = 1;
		  			tmpOut = fread(&(mi1[size1 + k].loc),
								   sizeof(int), 1, in1[j]);
		  			tmpOut = fread(&(mi1[size1 + k].err),
								   sizeof(int), 1, in1[j]);

					/* CIGAR. */
		  			tmpOut = fread(&(mi1[size1 + k].cigarSize),
								   sizeof(int), 1, in1[j]);
		  			tmpOut = fread((mi1[size1 + k].cigar),
								   sizeof(char), mi1[size1 + k].cigarSize,
								   in1[j]);
		  			mi1[size1 + k].cigar[mi1[size1 + k].cigarSize] = '\0';

					/* MD. */
		  			tmpOut = fread(&(mi1[size1 + k].mdSize),
								   sizeof(int), 1, in1[j]);
		  			tmpOut = fread((mi1[size1 + k].md),
								   sizeof(char), (mi1[size1 + k].mdSize),
								   in1[j]);
		  			mi1[size1 + k].md[mi1[size1 + k].mdSize] = '\0';

					/* If this is not mapped anywhere, make loc -loc and dir
					 * -1. ??? */
		  			if (mi1[size1 + k].loc < 1)
		  			{
						mi1[size1 + k].loc *= -1;
						mi1[size1 + k].dir = -1;
		  			}
				}

				/* Sort output with respect to their mapping locations. */
				qsort(mi1 + size1, size, sizeof(FullMappingInfo), compareOut);
				size1 += size;
		  	}
		}

		/* Load information regarding 2nd read of the current paired-end
		 * read. Load the mapping information into the FullMappingInfo mi2. */
		for (j = 0; j < _msf_openFiles; j++)
		{
			/* size = number of locations this 2nd pair mapped to. */
			tmpOut = fread(&size, sizeof(int), 1, in2[j]);
			if (size > 0)
			{
		  		for (k = 0; k < size; k++)
		  		{
					/* direction, location and edit distance. */
					mi2[size2 + k].dir = 1;
					tmpOut = fread(&(mi2[size2 + k].loc),
								   sizeof(int), 1, in2[j]);
					tmpOut = fread(&(mi2[size2 + k].err),
								   sizeof(int), 1, in2[j]);

					/* CIGAR. */
					tmpOut = fread(&(mi2[size2 + k].cigarSize),
								   sizeof(int), 1, in2[j]);
					tmpOut = fread((mi2[size2 + k].cigar),
								   sizeof(char), mi2[size2 + k].cigarSize,
								   in2[j]);
					mi2[size2 + k].cigar[mi2[size2 + k].cigarSize] = '\0';

					/* MD. */
					tmpOut = fread(&(mi2[size2 + k].mdSize),
								   sizeof(int), 1, in2[j]);
					tmpOut = fread((mi2[size2 + k].md),
								   sizeof(char), mi2[size2 + k].mdSize, in2[j]);
					mi2[size2 + k].md[mi2[size2 + k].mdSize] = '\0';

					/* If this is not mapped anywhere, make loc -loc and dir
					 * -1. ??? */
					if (mi2[size2 + k].loc < 1)
					{
			  			mi2[size2 + k].loc *= -1;
			  			mi2[size2 + k].dir = -1;
					}
		  		}

				/* Sort output with respect to their mapping locations. */
				qsort(mi2 + size2, size, sizeof(FullMappingInfo), compareOut);
		  		size2 += size;
			}
		}
		

		/* This if branch writes OEA information to file. */
		int lm, ll, rl, rm;
		int pos = 0;
		if (pairedEndDiscordantMode)
		{			
			/* Process the mapping locations of the 1st entry of the
			 * pair. Search in a Cartesian product way. Each of the pairs in the
			 * set (1st mapping locations) x (2nd mapping locations) are
			 * searched through.  This does not write it yet to the file, just
			 * marks the reads using _msf_seqList[i].hits[0] field. */
			for (j = 0; j < size1; j++)
			{
				lm = mi1[j].loc - maxPairEndedDiscordantDistance + 1;
				ll = mi1[j].loc - minPairEndedDiscordantDistance + 1;
				rl = mi1[j].loc + minPairEndedDiscordantDistance - 1;
				rm = mi1[j].loc + maxPairEndedDiscordantDistance - 1;

				/* Search for the 2nd enrty of the pair in the given
				 * interval. Stop when we hit to the right of lm. */
				while (pos < size2 && mi2[pos].loc < lm)
				{
				  	pos++;
				}				

				k = pos;
				while (k < size2 && mi2[k].loc <= rm)
				{							
					if (mi2[k].loc <= ll || mi2[k].loc >= rl)
					{
						/* Found it! This should be a concordant mapping? */
						if ((mi1[j].loc < mi2[k].loc && mi1[j].dir == 1 &&
							 mi2[k].dir == -1) ||
							(mi1[j].loc > mi2[k].loc && mi1[j].dir == -1 &&
							 mi2[k].dir == 1))
						{							
							/* Indicate the found information. */
						  	_msf_seqList[i * 2].hits[0]		= 1;
						  	_msf_seqList[i * 2 + 1].hits[0] = 1;

						  	if (nosamMode != 0)
						  	{
								size1 = 0;
								size2 = 0;
						  	}

						  	break;
						}
					}
					k++;
				}
			}

			/* Number of sequence mappings for the 1st and 2nd elements of the
			 * current read. */
			_msf_seqHits[i * 2]		+= size1;
			_msf_seqHits[i * 2 + 1] += size2;

			/* If the number of cartesian mapping locations is higher than some
			 * value and SAM mode is off, do something I do not understand.  */
			if (_msf_seqHits[i * 2 + 1] * _msf_seqHits[i * 2] >
				DISCORDANT_CUT_OFF && nosamMode != 0)
			{
				_msf_seqList[i * 2].hits[0]		= 1;
				_msf_seqList[i * 2 + 1].hits[0] = 1;
				size1							= 0;
				size2							= 0;
			}

			int		rNo = 0;
			int		loc = 0;
			int		err = 0;
			float	sc	= 0;
			char	l	= 0;

			/* If one of the pair elements could not be mapped, it means we have
			 * one-end anchored reads. 1st element mapped nowhere, output the
			 * mapping locations of the 2nd element of the pair. */
			if (_msf_seqHits[i * 2] == 0)
			{	
				/* 2nd element we are writing down. */
				for (k = 0;
					 k < size2 && _msf_oeaMapping[i * 2 + 1] < maxOEAOutput;
					 k++)
				{
					/* read number, location, error, score. */
					rNo = i * 2 + 1;
					loc = mi2[k].loc * mi2[k].dir;
					err = mi2[k].err;
					sc	= mi2[k].score;

					l = strlen(_msf_refGenName);

					tmp = fwrite(&rNo, sizeof(int), 1, out1);

					tmp = fwrite(&l, sizeof(char), 1, out1);
					tmp = fwrite(_msf_refGenName, sizeof(char), l, out1);

					tmp = fwrite(&loc, sizeof(int), 1, out1);
					tmp = fwrite(&err, sizeof(int), 1, out1);
					tmp = fwrite(&sc, sizeof(float), 1, out1);

					/* CIGAR */
					if (mi2[k].cigarSize > SEQ_LENGTH || mi2[k].cigarSize <= 0)
					  printf("ERROR  CIGAR size=%d %s\n", mi2[k].cigarSize,
						 _msf_seqList[i * 2 + 1].seq);

					tmp = fwrite(&(mi2[k].cigarSize), sizeof(int), 1, out1);
					tmp = fwrite((mi2[k].cigar), sizeof(char), mi2[k].cigarSize,
							 out1);

					/* MD */
					tmp = fwrite(&(mi2[k].mdSize), sizeof(int), 1, out1);
					tmp = fwrite((mi2[k].md), sizeof(char), mi2[k].mdSize,
							 out1);

					/* Increment the number of OEA mappings for this pair of the
					 * read. */
					_msf_oeaMapping[i * 2 + 1]++;
				}
			}

			/* If one of the pair elements could not be mapped, it means we have
			 * one-end anchored reads. 2nd element mapped nowhere, output the
			 * mapping locations of the 1sr element of the pair. */
			if (_msf_seqHits[i * 2 + 1] == 0)
			{
				/* 1st element we are writing down. */
				for (j = 0; j < size1 && _msf_oeaMapping[i * 2] < maxOEAOutput;
					 j++)
				{
					/* read number, location, error, score. */
					rNo = i * 2;
					loc = mi1[j].loc * mi1[j].dir;
					err = mi1[j].err;
					sc	= mi1[j].score;

					l = strlen(_msf_refGenName);

					tmp = fwrite(&rNo, sizeof(int), 1, out1);

					tmp = fwrite(&l, sizeof(char), 1, out1);
					tmp = fwrite(_msf_refGenName, sizeof(char), l, out1);

					tmp = fwrite(&loc, sizeof(int), 1, out1);
					tmp = fwrite(&err, sizeof(int), 1, out1);
					tmp = fwrite(&sc, sizeof(float), 1, out1);

					/* CIGAR */
					if (mi1[j].cigarSize > SEQ_LENGTH || mi1[j].cigarSize <= 0)
					  printf("ERROR %d %s\n", mi1[j].cigarSize,
						 _msf_seqList[i * 2 + 1].seq);

					tmp = fwrite(&(mi1[j].cigarSize), sizeof(int), 1, out1);
					tmp = fwrite((mi1[j].cigar), sizeof(char), mi1[j].cigarSize,
							 out1);

					/* MD */
					tmp = fwrite(&(mi1[j].mdSize), sizeof(int), 1, out1);
					tmp = fwrite((mi1[j].md), sizeof(char), mi1[j].mdSize,
							 out1);

					/* Increment the number of OEA mappings for this pair of the
					 * read. */
					_msf_oeaMapping[i * 2]++;
				}
			}
		} /* End paired-end discordant mode. */


		/* Information regarding 1st and 2nd pair. */
		char *seq1, *seq2, *rseq1, *rseq2, *qual1, *qual2;

		rqual1[SEQ_LENGTH] = '\0';
		rqual2[SEQ_LENGTH] = '\0';
		rqual1[0]		   = '\0';
		rqual2[0]		   = '\0';

		seq1  = _msf_seqList[i * 2].seq;
		rseq1 = _msf_seqList[i * 2].rseq;
		qual1 = _msf_seqList[i * 2].qual;
		strncpy(rqual1, _msf_seqList[i * 2].qual, SEQ_LENGTH);

		seq2  = _msf_seqList[i * 2 + 1].seq;
		rseq2 = _msf_seqList[i * 2 + 1].rseq;
		qual2 = _msf_seqList[i * 2 + 1].qual;
		strncpy(rqual2, _msf_seqList[i * 2 + 1].qual, SEQ_LENGTH);

		/* Calculate score of each mapping for both elements of the pair. Score
		 * is computed using quality of the reads. */
		if (pairedEndDiscordantMode)
		{
			for (k = 0; k < size1; k++)
			{
		  		mi1[k].score = calculateScore(
					mi1[k].loc, (mi1[k].dir == -1) ? rseq1 : seq1,
					(mi1[k].dir == -1) ? rqual1 : qual1, mi1[k].cigar);
			}

			for (k = 0; k < size2; k++)
			{
		  		mi2[k].score = calculateScore(
					mi2[k].loc, (mi2[k].dir == -1) ? rseq2 : seq2,
					(mi2[k].dir == -1) ? rqual2 : qual2, mi2[k].cigar);
			}
		}

		/* CALKAN MAPQ FOR PE */
		/* Compute tprob fields. */
		if (seqFastq)
		{
			for (j = 0; j < size1; j++)
			{
				if (mi1[j].err != 0)
				{
				  	bestHitMappingInfo[i*2].tprob +=
						mapProb(i*2, mi1[j].md, mi1[j].dir, mi1[j].err);
				}
			}
			for (k = 0; k < size2; k++)
			{
				if (mi2[k].err != 0)
				{
				  	bestHitMappingInfo[i*2+1].tprob +=
						mapProb((i*2+1), mi2[k].md, mi2[k].dir, mi2[k].err);
				}
			}
		}


		/* Check for the concordant and discordant mappings. Save them into
		 * bestHitMappingInfo. */
		if (pairedEndDiscordantMode)
		{
			/* For each cartesian product of the 1st element x 2nd element
			 * mapping locations, check for the discordant and concordant
			 * mappings. */
			for (j = 0; j < size1; j++)
			{
				for (k = 0; k < size2; k++)
				{
					/* Check for the pair if it is CONCORDANT. If it is, write
					 * to bestHitMappingInfo. If there are some discordant
					 * mappings found before finding the first concordant
					 * mapping, they had been written down into the
					 * BestFullMappingInfo struct for this read and here that
					 * information will be overwritten and the discordant
					 * mappings found hereafter will be disregarded.*/
					if (
						(pairedEndModePE && 
						 ((mi2[k].loc - mi1[j].loc >=
						   minPairEndedDiscordantDistance &&
						   mi2[k].loc - mi1[j].loc <=
						   maxPairEndedDiscordantDistance &&
						   mi1[j].dir > 0 && mi2[k].dir < 0) ||
						  (mi1[j].loc - mi2[k].loc >=
						   minPairEndedDiscordantDistance &&
						   mi1[j].loc - mi2[k].loc <=
						   maxPairEndedDiscordantDistance &&
						   mi1[j].dir < 0 && mi2[k].dir > 0)))
						||  // CALKAN MPPE
						(pairedEndModeMP &&
						 ((mi2[k].loc - mi1[j].loc >=
						   minPairEndedDiscordantDistance &&
						   mi2[k].loc - mi1[j].loc <=
						   maxPairEndedDiscordantDistance &&
						   mi1[j].dir < 0 && mi2[k].dir > 0) ||
						  (mi1[j].loc - mi2[k].loc >=
						   minPairEndedDiscordantDistance &&
						   mi1[j].loc - mi2[k].loc <=
						   maxPairEndedDiscordantDistance &&
						   mi1[j].dir > 0 && mi2[k].dir < 0)))
						)
					{
						//POSSIBLE CONCORDANT
						/* Store the concordant mapping information in
						 * bestHitMappingInfo. Both elements of the pair are
						 * stored (2*readNumber and 2*readNumber+1). */
						if(_msf_readHasConcordantMapping[i] == 0)
						{
							setPairFullMappingInfo(i, mi1[j], mi2[k]);
							_msf_readHasConcordantMapping[i] = 1;
							_msf_seqList[i * 2].hits[0]		 = 1;
							_msf_seqList[i * 2 + 1].hits[0]	 = 1;
						}
						else
						{
							/* If the found concordant mapping has a better edit
							 * distance than the processed ones, overwrite
							 * bestHitMappingInfo to store this better one. */
							if (bestHitMappingInfo[i * 2].err +
								bestHitMappingInfo[i * 2 + 1].err >=
								mi1[j].err + mi2[k].err)
							{
								/* If they have the same edit distance, it
								 * checks something else. */
								if (bestHitMappingInfo[i * 2].err +
									bestHitMappingInfo[i * 2 + 1].err ==
									mi1[j].err + mi2[k].err &&
									findNearest(
										abs(bestHitMappingInfo[i * 2 + 1].loc -
											bestHitMappingInfo[i * 2].loc),
										abs(mi2[k].loc - mi1[j].loc),
										meanDistanceMapping) == 0)
								{
									continue;
								}
								setPairFullMappingInfo(i, mi1[j], mi2[k]);
							}
						}
					}
					//DISCORDANT TO TEMP FILE FOR POST PROCESSING
					/* This writes down the discordant mappings to temporary
					 * file for post-processing until a concordant mapping has
					 * been found.  Also checks the elements of the pairs so
					 * both of them should be mapped somewhere in the
					 * genome. Recall that if single of them could be mapped, it
					 * becomes OEA mapping. */
					else if (_msf_readHasConcordantMapping[i] == 0 &&
							 _msf_seqHits[i * 2] != 0 &&
							 _msf_seqHits[i * 2 + 1] != 0)
					{
						int		rNo = i;
						int		loc = mi1[j].loc * mi1[j].dir;
						int		err = mi1[j].err;
						float	sc	= mi1[j].score;
						char	l	= strlen(_msf_refGenName);

						/* If we are still to write. */
						if (_msf_discordantMapping[i * 2] < maxDiscordantOutput)
						{
							/* 1st element of the pair. */
							tmp = fwrite(&rNo, sizeof(int), 1, out);
							tmp = fwrite(&l, sizeof(char), 1, out);
							tmp = fwrite(_msf_refGenName, sizeof(char), l, out);
							tmp = fwrite(&loc, sizeof(int), 1, out);
							tmp = fwrite(&err, sizeof(int), 1, out);
							tmp = fwrite(&sc, sizeof(float), 1, out);
							tmp = fwrite(&(mi1[j].cigarSize), sizeof(int), 1,
								 out);
							tmp = fwrite((mi1[j].cigar), sizeof(char),
								 mi1[j].cigarSize, out);
							tmp = fwrite(&(mi1[j].mdSize), sizeof(int), 1, out);
							tmp = fwrite((mi1[j].md), sizeof(char),
								 mi1[j].mdSize, out);

							/* 2nd element of the pair. */
							loc = mi2[k].loc * mi2[k].dir;
							err = mi2[k].err;
							sc	= mi2[k].score;
							tmp = fwrite(&loc, sizeof(int), 1, out);
							tmp = fwrite(&err, sizeof(int), 1, out);
							tmp = fwrite(&sc, sizeof(float), 1, out);
							tmp = fwrite(&(mi2[k].cigarSize), sizeof(int), 1,
								 out);
							tmp = fwrite((mi2[k].cigar), sizeof(char),
								 mi2[k].cigarSize, out);
							tmp = fwrite(&(mi2[k].mdSize), sizeof(int), 1, out);
							tmp = fwrite((mi2[k].md), sizeof(char),
								 mi2[k].mdSize, out);

							/* Increment discordant mapping count. */
							_msf_discordantMapping[i * 2]++;
						}
						
						//SET THE BEST DISCORDANT
						//BEGIN {Farhad Hormozdiari}
						/* This is the first discordant mapping pair, save it to
						 * the bestHitMappingInfo. Note that if a concordant
						 * mapping has been found beforehand, does nothing.*/
						if (bestHitMappingInfo[i * 2].loc == -1 &&
							bestHitMappingInfo[i * 2 + 1].loc == -1 &&
							_msf_readHasConcordantMapping[i] == 0)
						{
						  	setPairFullMappingInfo(i, mi1[j], mi2[k]);
						  	_msf_seqList[i * 2].hits[0] = 1;
						  	_msf_seqList[i * 2 + 1].hits[0] = 1;
						}
						/* Check if this is a better discordant mapping from the
						 * previously found discordant mappings. */
						else if (bestHitMappingInfo[i * 2].err +
								 bestHitMappingInfo[i * 2 + 1].err >=
								 mi1[j].err + mi2[k].err &&
								 _msf_readHasConcordantMapping[i] == 0)
						{
							/* If they have the same edit distance, check
							 * something else. */
							if (bestHitMappingInfo[i * 2].err +
								bestHitMappingInfo[i * 2 + 1].err ==
								mi1[j].err + mi2[k].err &&
								findNearest(
								   abs(
									   bestHitMappingInfo[i * 2 + 1].loc
									   - bestHitMappingInfo[i * 2].loc),
								   abs(mi1[j].loc - mi2[k].loc),
								   meanDistanceMapping) == 0)
							{
						  		continue;
							}
							setPairFullMappingInfo(i, mi1[j], mi2[k]);
						}
						//END {Farhad Hormozdiari}
					}	/* Discordant mapping. */
				} /* Map. locations of the 2nd element. */
			} /* Map. locations of the 1st element. */
		} /* End pairedEndDiscordantMode */
		else					/* Pure paired-end mode. */
		{
			/* For each cartesian product of the 1st element x 2nd element
			 * mapping locations. */
			for (j = 0; j < size1; j++)
			{
				for (k = 0; k < size2; k++)
				{
					/* Is it paired-end?  */
					if ((mi2[k].loc - mi1[j].loc >= minPairEndedDistance &&
						 mi2[k].loc - mi1[j].loc <= maxPairEndedDistance &&
						 mi1[j].dir > 0 && mi2[k].dir < 0)
						||
						(mi1[j].loc - mi2[k].loc >= minPairEndedDistance &&
						 mi1[j].loc - mi2[k].loc <= maxPairEndedDistance &&
						 mi1[j].dir < 0 && mi2[k].dir > 0))
					{
						char	*seq;
						char	*qual;
						char	 d1;
						char	 d2;
						int		 isize;
						int		 proper = 0;
						
						// ISIZE CALCULATION
						// The distance between outer edges
						isize = abs(mi1[j].loc - mi2[k].loc) + SEQ_LENGTH - 2;
						if (mi1[j].loc - mi2[k].loc > 0)
						  	isize *= -1;

						d1 = (mi1[j].dir == -1) ? 1 : 0;
						d2 = (mi2[k].dir == -1) ? 1 : 0;

						//SET THE READ HAS CONCORDANT MAPPING
						_msf_readHasConcordantMapping[i] = 1;

						/* Choisir normal or reverse sequence. */
						if (d1)
						{
						  	seq = rseq1;
						  	qual = rqual1;
						}
						else
						{
						  	seq = seq1;
						  	qual = qual1;
						}

						/* proper is used right below. */
						if ((mi1[j].loc < mi2[k].loc && !d1 && d2) ||
							(mi1[j].loc > mi2[k].loc && d1 && !d2))
						  	proper = 2;
						else
						  	proper = 0;

						/* Fill out the output related information (SAM and
						 * optional fields). The first elem of the pair is taken
						 * as the basis. */
						_msf_output.POS	   = mi1[j].loc;
						_msf_output.MPOS   = mi2[k].loc;
						_msf_output.FLAG   = 1 + proper + 16 * d1 + 32 * d2 + 64;
						_msf_output.ISIZE  = isize;
						_msf_output.SEQ	   = seq;
						_msf_output.QUAL   = qual;
						_msf_output.QNAME  = _msf_seqList[i * 2].name;
						_msf_output.RNAME  = _msf_refGenName;
						_msf_output.MAPQ   = 255;
						_msf_output.CIGAR  = cigar;
						_msf_output.MRNAME = "=";

						_msf_output.optSize	  = 2;
						_msf_output.optFields = _msf_optionalFields;

						_msf_optionalFields[0].tag	= "NM";
						_msf_optionalFields[0].type = 'i';
						_msf_optionalFields[0].iVal = mi1[j].err;

						_msf_optionalFields[1].tag	= "MD";
						_msf_optionalFields[1].type = 'Z';
						_msf_optionalFields[1].sVal = mi1[j].md;


						/* Output if not in bestMode. */
						if (!bestMode)
						  	output(_msf_output);

						if (d2)
						{
						  	seq	 = rseq2;
						  	qual = rqual2;
						}
						else
						{
						  	seq	 = seq2;
						  	qual = qual2;
						}

						/* Different from above info basically in POS and MPOS
						 * fields. The second elem of the pair is taken as the
						 * basis. */
						_msf_output.POS	   = mi2[k].loc;
						_msf_output.MPOS   = mi1[j].loc;
						_msf_output.FLAG   = 1 + proper + 16 * d2 + 32 * d1 + 128;
						_msf_output.ISIZE  = -isize;
						_msf_output.SEQ	   = seq;
						_msf_output.QUAL   = qual;
						_msf_output.QNAME  = _msf_seqList[i * 2].name;
						_msf_output.RNAME  = _msf_refGenName;
						_msf_output.MAPQ   = 255;
						_msf_output.CIGAR  = cigar;
						_msf_output.MRNAME = "=";

						_msf_output.optSize	  = 2;
						_msf_output.optFields = _msf_optionalFields;

						_msf_optionalFields[0].tag	= "NM";
						_msf_optionalFields[0].type = 'i';
						_msf_optionalFields[0].iVal = mi2[k].err;


						_msf_optionalFields[1].tag	= "MD";
						_msf_optionalFields[1].type = 'Z';
						_msf_optionalFields[1].sVal = mi2[k].md;

						/* Output if not in bestMode. */
						if (!bestMode)
						  	output(_msf_output);
						
						//SET THE BEST CONCORDANT
						//BEGIN {Farhad Hormozdiari}
						/* Store the concordant mapping information. First found
						 * concordant mapping for this read if both loc = -1. */
						if (bestHitMappingInfo[i * 2].loc == -1 &&
							bestHitMappingInfo[i * 2 + 1].loc == -1)
						{
						  	setPairFullMappingInfo(i, mi1[j], mi2[k]);
						}
						else
						{
							/* If the found concordant mapping has a better edit
							 * distance than the processed ones, overwrite
							 * bestHitMappingInfo to store this better one. */
							if (bestHitMappingInfo[i * 2].err +
								bestHitMappingInfo[i * 2 + 1].err >=
								mi1[j].err + mi2[k].err)
							{
								/* If they have the same edit distance, it
								 * checks something else. */
								if (bestHitMappingInfo[i * 2].err +
									bestHitMappingInfo[i * 2 + 1].err ==
									mi1[j].err + mi2[k].err &&
									findNearest(
										abs(
											bestHitMappingInfo[i * 2+ 1].loc
											- bestHitMappingInfo[i* 2].loc),
										   abs(mi2[k].loc - mi1[j].loc),
										   meanDistanceMapping) == 0)
								{
								  	continue;
								}
								setPairFullMappingInfo(i, mi1[j], mi2[k]);
							}
						}
						//END   {Farhad Hormozdiari}
					}	/* Is it paired end = concordant mapping? */
				}		/* 2nd elem */
			}			/* 1st elem */
		}				/* Paired end mode */
	}


	/* Close files */
	if (pairedEndDiscordantMode)
	{
	  	fclose(out);
	  	fclose(out1);
	}

	for (i = 0; i < _msf_openFiles; i++)
	{
		fclose(in1[i]);
		fclose(in2[i]);

		unlink(fname1[i]);
		unlink(fname2[i]);
	}

	tmp++;

	freeMem(mi1, sizeof(FullMappingInfo) * _msf_maxLSize);
	freeMem(mi2, sizeof(FullMappingInfo) * _msf_maxRSize);

	_msf_openFiles = 0;

	/* calkan counter */
	int unmappedCnt = 0;
	for (i = 0; i < _msf_seqListSize; i++)
	{
	  	if (_msf_seqHits[i] == 0) unmappedCnt++;
	}

	/* mappedSeqCnt is global. */
	mappedSeqCnt = _msf_seqListSize - unmappedCnt;
}




/*------------------------------------------------------------------------------
 * Simple util function.
 *----------------------------------------------------------------------------*/
float
str2int ( char *str, int index1, int index2 )
{
	char	tmp[SEQ_MAX_LENGTH];
	strncpy(tmp, &str[index1], index2 - index1);
	tmp[index2 - index1] = '\0';
	return atol(tmp);
}



/*------------------------------------------------------------------------------
 * Simple util function.
 *----------------------------------------------------------------------------*/
double
binomial_coefficient ( int n, int k )
{
	double	ret;
	int		i;
	ret = 1.0;

	for (i=0; i<k; i++)
	{
	  	ret *= (n - i);
	  	ret /= (k - i);
	}

	return ret; 
}



/*------------------------------------------------------------------------------
 * Calculates quality of the given read sequence using matching information.
 *----------------------------------------------------------------------------*/
float
calculateScore ( int index, char *seq, char *qual, char *md )
{
	int		 i;
	int		 j;
	char	*ref;
	char	*ver;

	ref			  = _msf_refGen + index - 1;
	ver			  = seq;
	float	score = 1;

	char	tmp[2 * SEQ_MAX_LENGTH];
	int		value  = 0;
	int		end	   = 0;
	int		index1 = 0;
	int		index2 = 0;

	/* Fill in tmp variable using matching information. */
	i = 0;
	while (1)
	{
		if (i >= strlen(md))
		  	break;

		index1 = i;

		while (md[i] >= '0' && md[i] <= '9')
		{
		  	i++;
		}

		index2 = i;

		value = str2int(md, index1, index2);

		if (md[i] == 'M')
		{
			for (j = 0; j < value; j++)
			{
				tmp[end] = 'M';
				end++;
			}
		}
		else if (md[i] == 'I')
		{
			for (j = 0; j < value; j++)
			{
		  		tmp[end] = 'I';
		  		end++;
			}
		}
		else if (md[i] == 'D')
		{
			for (j = 0; j < value; j++)
			{
		  		tmp[end] = 'D';
		  		end++;
			}
		}
		i++;
	}

	tmp[end] = '\0';

	j = 0;

	/* Calculate score. */
	for (i = 0; i < end; i++)
	{
		if (tmp[i] == 'M')
		{
			if (*ref != *ver)
			{
		  		score *= 0.001 + 1 / pow(10, ((qual[j] - 33) / 10.0));
			}

			ref++;
			ver++;
			j++;
		}
		else if (tmp[i] == 'I')
		{
			ver++;
			j++;
			score *= 0.0003;  // 0.0001 + 0.0002;  0.0001: indel rate in normal human, 0.0002: indel error rate in Illumina
		}
		else if (tmp[i] == 'D')
		{
			ref++;
			score *= 0.0003; // 0.0001 + 0.0002
		}
	}

	return score;
}



/*------------------------------------------------------------------------------
 * Pffff.
 *----------------------------------------------------------------------------*/
int
matoi ( char *str, int start, int end )
{
	int		i = 0;
	char	tmp[SEQ_MAX_LENGTH];

	for (i = 0; i < end - start; i++)
	  	tmp[i] = str[start + i];
	tmp[i] = '\0';

	return atoi(tmp);
}



/*------------------------------------------------------------------------------
 * Using CIGAR information, fills the matrix data (I, M, D).
 * Used for outputPairedEndDiscPP.
 *----------------------------------------------------------------------------*/
void
convertCigarToMatrix ( char *cigar, int cigar_size, char * matrix )
{
	int i	  = 0;
	int j	  = 0;
	int start = 0;
	int size  = 0;
	matrix[0] = '\0';

	while (i < cigar_size)
	{
		if (cigar[i] >= '0' && cigar[i] <= '9')
		{
			start = i;
			while (cigar[i] >= '0' && cigar[i] <= '9' && i < cigar_size)
		  		i++;

			int value = matoi(cigar, start, i);
			for (j = 0; j < value; j++)
			{
				if (cigar[i] == 'M')
				  matrix[size] = 'M';
				else if (cigar[i] == 'D')
				  matrix[size] = 'D';
				else if (cigar[i] == 'I')
				  matrix[size] = 'I';
				size++;
			}
		}
		i++;
	}
	matrix[size] = '\0';
}



/*------------------------------------------------------------------------------
 * Using MD information, fills the matrix data (I, M, D).
 * Used for outputPairedEndDiscPP.
 *----------------------------------------------------------------------------*/
void
convertMDToMatrix ( char *md, int md_size, char * matrix )
{
	int i	  = 0;
	int j	  = 0;
	int start = 0;
	int size  = 0;
	matrix[0] = '\0';

	while (i < md_size)
	{
		if (md[i] >= '0' && md[i] <= '9')
		{
			start = i;

			while (md[i] >= '0' && md[i] <= '9' && i < md_size)
		  		i++;

			int value = matoi(md, start, i);
			for (j = 0; j < value; j++)
			{
		  		matrix[size] = 'M';
		  		size++;
			}
			i--;
		}
		else if (md[i] == '^')
		{
		  	matrix[size] = 'D';
		  	size++;
		}
		else
		{
		  	matrix[size] = md[i];
		  	size++;
		}
		i++;
	}
	matrix[size] = '\0';
}



/*------------------------------------------------------------------------------
 * Uses both MD and CIGAR information to generate matching information writing
 * it into matrix
 * Used for outputPairedEndDiscPP. 
 *----------------------------------------------------------------------------*/
void
convertMDCigarToMatrix ( char *cigar, int cigar_size, char *md, int md_size,
						 char *matrix)
{
	int i = 0;
	int j = 0;

	int size = 0;

	char tmp1[SEQ_MAX_LENGTH];
	char tmp2[SEQ_MAX_LENGTH];
	
	convertMDToMatrix(md, md_size, tmp2);
	convertCigarToMatrix(cigar, cigar_size, tmp1);

	while (i < strlen(tmp1))
	{
		if (tmp1[i] == 'M')
		{
			if (j < strlen(tmp2))
			{
				if (tmp2[j] == 'M')
				{
				  	matrix[size] = 'M';
				  	size++;
				}
				if (tmp2[j] != 'M')
				{
				  	matrix[size] = tmp2[j];
				  	size++;
				}
			}
			else
			{
		  		matrix[size] = 'M';
		  		size++;
			}
		}
		else if (tmp1[i] == 'D')
		{
			matrix[size] = 'D';
			size++;
			j++;
			matrix[size] = tmp2[j];
			size++;
		}
		else if (tmp1[i] == 'I')
		{
		  	matrix[size] = 'I';
		  	size++;
		}

		i++;
		if (j < strlen(tmp2))
		  	j++;
	}

	if (strlen(tmp1))
	  	matrix[size] = '\0';
}



/*------------------------------------------------------------------------------
 *	PENDING
 *	Used for outputPairedEndDiscPP. Did not go into the details of it.
 *----------------------------------------------------------------------------*/
void
convertInsertion ( char * in_matrix, char * seq, char *out_matrix )
{
	int i	 = 0;
	int j	 = 0;
	int size = 0;

	while (i < strlen(in_matrix))
	{
		if (in_matrix[i] == 'M')
		{
		  	out_matrix[size] = 'M';
		  	size++;
		  	j++;
		}
		else if (in_matrix[i] == 'D')
		{
		  	out_matrix[size] = 'D';
		  	size++;
		  	i++;
		  	j++;
		  	out_matrix[size] = seq[j];
		  	j++;
		  	size++;
		}
		else if (in_matrix[i] == 'I')
		{
		  	out_matrix[size] = 'I';
		  	size++;
		  	out_matrix[size] = seq[j];
		  	size++;
		  	j++;
		}
		else
		{
		  	out_matrix[size] = in_matrix[i];
		  	size++;
		  	j++;
		}
		i++;
	}
	out_matrix[size] = '\0';
}



/*------------------------------------------------------------------------------
 *	Reads already dumped information from ...__disc file and writes ..._DIVET.vh
 *	file.
 *----------------------------------------------------------------------------*/
void
outputPairedEndDiscPP (   )
{
	char tmp_matrix1[SEQ_MAX_LENGTH];
	char tmp_matrix2[SEQ_MAX_LENGTH];

	char matrix1[SEQ_MAX_LENGTH];
	char matrix2[SEQ_MAX_LENGTH];

	char cigar1[MAX_CIGAR_SIZE];
	char editString1[2 * SEQ_MAX_LENGTH];

	char cigar2[MAX_CIGAR_SIZE];
	char editString2[2 * SEQ_MAX_LENGTH];
	char seq1[SEQ_LENGTH + 1];

	char seq2[SEQ_LENGTH + 1];

	char	 genName[SEQ_LENGTH];
	char	 fname1[FILE_NAME_LENGTH];
	char	 fname2[FILE_NAME_LENGTH];
	char	 l;
	int		 l_size;
	int		 loc1, loc2;
	int		 err1, err2;
	char	 dir1, dir2;
	float	 sc1, sc2, lsc = 0;
	int		 flag		   = 0;
	int		 rNo, lrNo	   = -1;
	int		 tmp;
	FILE	*in, *out;

	sprintf(fname1, "%s__%s__disc", mappingOutputPath, mappingOutput);
	sprintf(fname2, "%s%s_DIVET.vh", mappingOutputPath, mappingOutput);

	in	= fileOpen(fname1, "r");
	out = fileOpen(fname2, "w");

	if (in != NULL)
	{
	  	flag = fread(&rNo, sizeof(int), 1, in);
	}
	else
	{
	  	flag = 0;
	}

	seq1[SEQ_LENGTH] = '\0';
	seq2[SEQ_LENGTH] = '\0';

	/* Read the binary file __disc and generate vh file. */
	while (flag)
	{
		tmp				 = fread(&l, sizeof(char), 1, in);
		tmp				 = fread(genName, sizeof(char), l, in);
		genName[(int) l] = '\0';
		tmp				 = fread(&loc1, sizeof(int), 1, in);
		tmp				 = fread(&err1, sizeof(int), 1, in);
		tmp				 = fread(&sc1, sizeof(float), 1, in);

		tmp					 = fread(&l_size, sizeof(int), 1, in);
		tmp					 = fread(cigar1, sizeof(char), l_size, in);
		cigar1[(int) l_size] = '\0';

		tmp						  = fread(&l_size, sizeof(int), 1, in);
		tmp						  = fread(editString1, sizeof(char), l_size, in);
		editString1[(int) l_size] = '\0';

		tmp = fread(&loc2, sizeof(int), 1, in);
		tmp = fread(&err2, sizeof(int), 1, in);
		tmp = fread(&sc2, sizeof(float), 1, in);

		tmp					 = fread(&l_size, sizeof(int), 1, in);
		tmp					 = fread(cigar2, sizeof(char), l_size, in);
		cigar2[(int) l_size] = '\0';

		tmp						  = fread(&l_size, sizeof(int), 1, in);
		tmp						  = fread(editString2, sizeof(char), l_size, in);
		editString2[(int) l_size] = '\0';

		convertMDCigarToMatrix(cigar1, strlen(cigar1), editString1,
							   strlen(editString1), tmp_matrix1);
		convertMDCigarToMatrix(cigar2, strlen(cigar2), editString2,
							   strlen(editString2), tmp_matrix2);


		/* CALKAN: GO OVER THIS VERY CAREFULLY FOR PE vs MP */

		/* Read should not have concordant mapping. */
		if (_msf_readHasConcordantMapping[rNo] == 0 &&
			_msf_discordantMapping[rNo * 2] < maxDiscordantOutput )
		{
			dir1 = dir2 = 'F';
			strncpy(seq1, _msf_seqList[rNo * 2].seq, SEQ_LENGTH);
			strncpy(seq2, _msf_seqList[rNo * 2 + 1].seq, SEQ_LENGTH);

			if (loc1 < 0)
			{
		  		dir1 = 'R';
		  		loc1 = -loc1;
		  		strncpy(seq1, _msf_seqList[rNo * 2].rseq, SEQ_LENGTH);
			}

			if (loc2 < 0)
			{
		  		dir2 = 'R';
		  		loc2 = -loc2;
		  		strncpy(seq2, _msf_seqList[rNo * 2 + 1].rseq, SEQ_LENGTH);
			}

			convertInsertion(tmp_matrix1, seq1, matrix1);
			convertInsertion(tmp_matrix2, seq2, matrix2);

			if (rNo != lrNo)
			{
				int j;
				for (j = 0; j < SEQ_LENGTH; j++)
				{
				  	lsc += _msf_seqList[rNo * 2].qual[j] +
						_msf_seqList[rNo * 2 + 1].qual[j];
				}
				lsc /= 2 * SEQ_LENGTH;
				lsc -= 33;
				lrNo = rNo;
			}

			char event = '\0';

			/* Dont know what event is here. */
			if (dir1 == dir2)
			{
		  		event = 'V';
			} 
			else
			{
		  		if (pairedEndModePE && loc1 < loc2 && dir1 == 'R' && dir2 == 'F') 
					event = 'E';
		  		else if (pairedEndModeMP && loc1 < loc2 && dir1 == 'F' && dir2 == 'R') 
					event = 'E';
		  		else if (pairedEndModePE && loc2 < loc1 && dir1 == 'F' && dir2 == 'R') 
					event = 'E';
		  		else if (pairedEndModeMP && loc2 < loc1 && dir1 == 'R' && dir2 == 'F') 
					event = 'E';
		  		else if (abs(loc2 - loc1) >= maxPairEndedDiscordantDistance) 
					event = 'D';
		  		else 
					event = 'I';	    	 
			}

			/* Write out here. Why do you make hits[0] = 2? */
			_msf_seqList[rNo * 2].hits[0] = 2;
			fprintf(out,
					"%s\t%s\t%d\t%d\t%c\t=\t%d\t%d\t%c\t%c\t%d\t%0.0f\t%e\n",
					_msf_seqList[rNo * 2].name, genName, loc1,
					(loc1 + SEQ_LENGTH - 1), dir1, loc2,
					(loc2 + SEQ_LENGTH - 1), dir2, event, (err1 + err2),
					lsc, sc1 * sc2);

			//	      lsc, sc1 * sc2 * binomial_coefficient(2 * SEQ_LENGTH, (err1 + err2)));

		}
		flag = fread(&rNo, sizeof(int), 1, in);
	}

	/* Hallelujah! I saw the light! */
	tmp++;

	fclose(in);
	fclose(out);

	unlink(fname1);
}



/*------------------------------------------------------------------------------
 * Reads the already written ..._oea file in outputPairedEnd and infers
 * information to write the ...__OEA.sam file.
 *----------------------------------------------------------------------------*/
void
finalizeOEAReads ( char *fileName )
{
	FILE	*fp_out1;
	FILE	*in;

	char	genName[SEQ_LENGTH];
	char	fname1[FILE_NAME_LENGTH];
	char	fname2[FILE_NAME_LENGTH];

	char	 l		   = 0;
	int		 loc1	   = 0;
	int		 err1;
	char	 d;
	float	 sc1	   = 0;
	int		 flag	   = 0;
	int		 rNo	   = -1;
	int		 tmp	   = 0;
	int		 cigarSize = 0;
	int		 mdSize	   = 0;
	char	 cigar[SEQ_LENGTH + 1];
	char	 md[SEQ_LENGTH + 1];
	char	*seq1, *seq2, *qual1, *qual2;
	char	 rqual1[SEQ_LENGTH + 1];

	seq1	  = NULL;
	seq2	  = NULL;
	qual1	  = NULL;
	qual2	  = NULL;
	rqual1[0] = '\0';

	sprintf(fname1, "%s%s_OEA.sam", mappingOutputPath, mappingOutput);
	fp_out1 = fileOpen(fname1, "w");

	/* This information was actually saved to ...__oea file in
	 * pairedEndMode. Here we read it and write the corresponding SAM file? */
	SAMheaderTX(fp_out1, 0);
	in = NULL;
	if (pairedEndDiscordantMode)
	{
	  	sprintf(fname2, "%s__%s__oea", mappingOutputPath, mappingOutput);
	  	in = fileOpen(fname2, "r");
	}

	/* The first entry is the read id. */
	if (in != NULL)
	{
	  	flag = fread(&rNo, sizeof(int), 1, in);
	}
	else
	{
	  	flag = 0;
	}

	/* Got any more reads? */
	while (flag)
	{
		cigar[0] = '\0';
		md[0]	 = '\0';

		/* l is the length of the reference genome name. */
		tmp = fread(&l, sizeof(char), 1, in);
		tmp = fread(genName, sizeof(char), l, in);

		genName[(int) l] = '\0';

		/* location, error and score. */
		tmp = fread(&loc1, sizeof(int), 1, in);
		tmp = fread(&err1, sizeof(int), 1, in);
		tmp = fread(&sc1, sizeof(float), 1, in);

		/* CIGAR */
		tmp = fread(&cigarSize, sizeof(int), 1, in);
		tmp = fread(cigar, sizeof(char), cigarSize, in);
		cigar[cigarSize] = '\0';

		/* MD */
		tmp = fread(&mdSize, sizeof(int), 1, in);
		tmp = fread(md, sizeof(char), mdSize, in);
		md[mdSize] = '\0';

		/* If loc < 0, need to reverse. seq1 = /1 of the read.  */
		d = 1;
		if (loc1 < 0)
		{
			d	  = -1;
			loc1 *= -1;
			seq1  = _msf_seqList[rNo].rseq;
			reverse(_msf_seqList[rNo].qual, rqual1, SEQ_LENGTH);
			rqual1[SEQ_LENGTH] = '\0';
		}
		else
		{
		  	seq1  = _msf_seqList[rNo].seq;
		  	qual1 = _msf_seqList[rNo].qual;
		}

		/* seq2 = /2 of the read. */
		if (rNo % 2 == 0)
		{
		  	seq2  = _msf_seqList[rNo + 1].seq;
		  	qual2 = _msf_seqList[rNo + 1].qual;
		}
		else
		{
		  	seq2  = _msf_seqList[rNo - 1].seq;
		  	qual2 = _msf_seqList[rNo - 1].qual;
		}

		/* Output if this is a OEA read and it could be mapped. Uses _msf_output
		 * global variable. */
		if ( _msf_seqHits[rNo] != 0 && _msf_seqHits[rNo] < maxOEAOutput &&
			 _msf_seqHits[(rNo % 2 == 0) ? rNo + 1 : rNo - 1] == 0)
		{
			_msf_output.POS	   = loc1;
			_msf_output.MPOS   = 0;
			_msf_output.FLAG   = (rNo % 2 == 0) ?
				1 + 4 + 32 * d + 128 : 1 + 8 + 16 * d + 64;
			_msf_output.ISIZE  = 0;
			_msf_output.SEQ	   = seq1;
			_msf_output.QUAL   = qual1;
			_msf_output.QNAME  = _msf_seqList[rNo].name;
			_msf_output.RNAME  = genName;
			_msf_output.MAPQ   = 255;
			_msf_output.CIGAR  = cigar;
			_msf_output.MRNAME = "=";

			_msf_output.optSize = 4;
			_msf_output.optFields = _msf_optionalFields;

			_msf_optionalFields[0].tag = "NM";
			_msf_optionalFields[0].type = 'i';
			_msf_optionalFields[0].iVal = err1;

			_msf_optionalFields[1].tag = "MD";
			_msf_optionalFields[1].type = 'Z';
			_msf_optionalFields[1].sVal = md;

			//for the OEA reads
			_msf_optionalFields[2].tag = "NS";
			_msf_optionalFields[2].type = 'Z';
			_msf_optionalFields[2].sVal = seq2;

			_msf_optionalFields[3].tag = "NQ";
			_msf_optionalFields[3].type = 'Z';
			_msf_optionalFields[3].sVal = qual2;

			/* Output dem dataz. */
			outputSAM(fp_out1, _msf_output);

			/* Why does it still alter hits[0] fields? Hope this
			 * will come clear later. */
			_msf_seqList[rNo].hits[0] = -1;
			_msf_seqList[(rNo % 2 == 0) ? rNo + 1 : rNo - 1].hits[0] = -1;
		}
		else if(_msf_seqHits[rNo] != 0 && _msf_seqHits[(rNo % 2 == 0) ? rNo + 1 : rNo - 1] == 0)
		{						/* This is /2 */
			_msf_seqList[rNo].hits[0] = -1;
			_msf_seqList[(rNo % 2 == 0) ? rNo + 1 : rNo - 1].hits[0] = -1;
		}

		/* Get the other read sequence. */
		flag = fread(&rNo, sizeof(int), 1, in);
	}

	/* Yes this makes very much sense to me! */
	tmp++;

	fclose(in);
	unlink(fname2);

	fclose(fp_out1);
}



/* disabled until completed 


   void outputTransChromosomal(char *fileName1, char *fileName2, FILE * fp_out)
   {
   int i = 0;
   int j = 0;
   int k = 0;

   char *index;

   int size1 = 0;
   int size2 = 0;

   FILE *fp1 = NULL;
   FILE *fp2 = NULL;

   char geneFileName1[FILE_NAME_LENGTH];
   char geneFileName2[FILE_NAME_LENGTH];

   FullMappingInfoLink *miL = getMem(_msf_seqListSize * sizeof(FullMappingInfoLink));
   FullMappingInfoLink *miR = getMem(_msf_seqListSize * sizeof(FullMappingInfoLink));


   if(fileName1 != NULL && fileName2 != NULL)
   {

   fp1 = fileOpen(fileName1, "r");
   fp2 = fileOpen(fileName2, "r");

   index = strstr(fileName1, "__");
   strncpy(geneFileName1, index + 2 * sizeof(char), strstr(index + 2, "__") - index - 2);
   geneFileName1[strstr(index + 2, "__") - index - 2] = '\0';

   index = strstr(fileName2, "__");
   strncpy(geneFileName2, index + 2 * sizeof(char), strstr(index + 2, "__") - index - 2);
   geneFileName2[strstr(index + 2, "__") - index - 2] = '\0';


   for(i = 0; i < _msf_seqListSize / 2; i++)
   {
   fread(&size1, sizeof(int), 1, fp1);
   fread(&size2, sizeof(int), 1, fp2);

   miL[i].mi = getMem(size1 * sizeof(FullMappingInfo) );
   miR[i].mi = getMem(size2 * sizeof(FullMappingInfo) );

   miL[i].size = size1;
   miR[i].size = size2;

   for(j = 0; j < size1; j++)
   {
   fread(&(miL[i].mi[j].loc), sizeof(int), 1, fp1);

   fread (&(miL[i].mi[j].err), sizeof(int), 1, fp1);

   fread (&(miL[i].mi[j].cigarSize), sizeof(int), 1, fp1);
   fread ((miL[i].mi[j].cigar), sizeof(char), miL[i].mi[j].cigarSize+1, fp1);

   fread (&(miL[i].mi[j].mdSize), sizeof(int), 1, fp1);
   fread ((miL[i].mi[j].md), sizeof(char), miL[i].mi[j].mdSize+1, fp1);

   miL[i].mi[j].dir = 1;
   if(miL[i].mi[j].loc < 1)
   {
   miL[i].mi[j].loc *= -1;
   miL[i].mi[j].dir = -1;
   }
   }
   for(k = 0; k < size2; k++)
   {
   fread(&(miR[i].mi[k].loc), sizeof(int), 1, fp2);

   fread (&(miR[i].mi[k].err), sizeof(int), 1, fp2);

   fread (&(miR[i].mi[k].cigarSize), sizeof(int), 1, fp2);
   fread ((miR[i].mi[k].cigar), sizeof(char), miR[i].mi[k].cigarSize+1, fp2);

   fread (&(miR[i].mi[k].mdSize), sizeof(int), 1, fp2);
   fread ((miR[i].mi[k].md), sizeof(char), miR[i].mi[k].mdSize+1, fp2);

   miR[i].mi[k].dir = 1;
   if(miR[i].mi[k].loc < 1)
   {
   miR[i].mi[k].loc *= -1;
   miR[i].mi[k].dir = -1;
   }
   }
   if(_msf_readHasConcordantMapping[i] == 0 && size1 != 0 && size2 != 0 && (size1 * size2 < MAX_TRANS_CHROMOSAL_OUTPUT))
   {
   int d1 = 0;
   int d2 = 0;
   char *seq, *qual;
   char *seq1, *seq2, *rseq1, *rseq2, *qual1, *qual2;
   char rqual1[SEQ_LENGTH+1], rqual2[SEQ_LENGTH+1];
   rqual1[SEQ_LENGTH] = rqual2[SEQ_LENGTH] = '\0';
   seq1 = _msf_seqList[i*2].seq;
   rseq1 = _msf_seqList[i*2].rseq;
   qual1 = _msf_seqList[i*2].qual;
   reverse(_msf_seqList[i*2].qual, rqual1, SEQ_LENGTH);

   seq2 = _msf_seqList[i*2+1].seq;
   rseq2 = _msf_seqList[i*2+1].rseq;
   qual2 = _msf_seqList[i*2+1].qual;
   reverse(_msf_seqList[i*2+1].qual, rqual2, SEQ_LENGTH);

   for(j = 0; j < size1; j++)
   {
   d1 = (miL[i].mi[j].dir == -1)?1:0;

   if ( d1 )
   {
   seq = rseq1;
   qual = rqual1;
   }
   else
   {
   seq = seq1;
   qual = qual1;
   }

   for(k = 0; k < size2; k++)
   {

   d2 = (miR[i].mi[k].dir == -1)?1:0;

   _msf_output.POS                 = miL[i].mi[j].loc;
   _msf_output.MPOS                = miR[i].mi[k].loc;
   _msf_output.FLAG                = 0;
   _msf_output.ISIZE               = 0;
   _msf_output.SEQ                 = seq;
   _msf_output.QUAL                = qual;
   _msf_output.QNAME               = _msf_seqList[i*2].name;
   _msf_output.RNAME               = geneFileName1;
   _msf_output.MAPQ                = 255;
   _msf_output.CIGAR               = miL[i].mi[j].cigar;
   _msf_output.MRNAME              = "=";

   _msf_output.optSize     = 2;
   _msf_output.optFields   = _msf_optionalFields;

   _msf_optionalFields[0].tag = "NM";
   _msf_optionalFields[0].type = 'i';
   _msf_optionalFields[0].iVal = miL[i].mi[j].err;

   _msf_optionalFields[1].tag = "MD";
   _msf_optionalFields[1].type = 'Z';
   _msf_optionalFields[1].sVal = miL[i].mi[j].md;


   if ( d2 )
   {
   seq = rseq2;
   qual = rqual2;
   }
   else
   {
   seq = seq2;
   qual = qual2;
   }

   outputSAM(fp_out, _msf_output);


   _msf_output.POS                 = miR[i].mi[k].loc;
   _msf_output.MPOS                = miL[i].mi[j].loc;
   _msf_output.FLAG                = 0;
   _msf_output.ISIZE               = 0;
   _msf_output.SEQ                 = seq;
   _msf_output.QUAL                = qual;
   _msf_output.QNAME               = _msf_seqList[i*2+1].name;
   _msf_output.RNAME               = geneFileName2;
   _msf_output.MAPQ                = 255;
   _msf_output.CIGAR               = miR[i].mi[k].cigar;
   _msf_output.MRNAME              = "=";

   _msf_output.optSize     = 2;
   _msf_output.optFields   = _msf_optionalFields;

   _msf_optionalFields[0].tag = "NM";
   _msf_optionalFields[0].type = 'i';
   _msf_optionalFields[0].iVal = miR[i].mi[k].err;

   _msf_optionalFields[1].tag = "MD";
   _msf_optionalFields[1].type = 'Z';
   _msf_optionalFields[1].sVal = miR[i].mi[k].md;

   outputSAM(fp_out, _msf_output);

   }
   }
   }
   }

   }

   for(i = 0; i < _msf_seqListSize / 2; i++)
   {
   freeMem(miL[i].mi, miL[i].size * sizeof(FullMappingInfo));
   freeMem(miR[i].mi, miR[i].size * sizeof(FullMappingInfo));
   }

   freeMem(miL, _msf_seqListSize * sizeof(FullMappingInfoLink));
   freeMem(miR, _msf_seqListSize * sizeof(FullMappingInfoLink));

   fclose(fp1);
   fclose(fp2);
   }

*/

/*
  if flag is 1 it will output all the possible trans chromsal mapping
  otherwise only tmp file will be delete

*/


/*------------------------------------------------------------------------------
 * This is disabled, although it is called from baseFast, it just returns.
 *----------------------------------------------------------------------------*/
void
outputAllTransChromosomal ( int flag )
{
  return;
  /* disabled until completed

     int i = 0;
     int j = 0;
     int k = 0;
     int l = 0;

     FILE *fp_out = NULL;
     char fname1[FILE_NAME_LENGTH];

     if(flag)
     {
     fp_out = fileOpen(fname1, "w");

     sprintf(fname1, "%s%s_TRANSCHROMOSOMAL", mappingOutputPath, mappingOutput);

     i  = 0;
     for(j = i+1; j < _msf_maxFile; j++)
     {
     if(i != j)
     {
     for(k = 0; k < _msf_fileCount[i]; k++)
     {
     for(l = 0; l < _msf_fileCount[j]; l++)
     {
     outputTransChromosomal(_msf_fileName[i][k][0], _msf_fileName[j][l][1], fp_out);
     }// for l
     }// for k
     }// if
     }// for j
     }

     for(i = 0; i < _msf_maxFile; i++)
     {
     for(j = 0; j < _msf_fileCount[i]; j++)
     {
     unlink(_msf_fileName[i][j][0]);
     unlink(_msf_fileName[i][j][1]);
     }
     }
     if(flag)
     fclose(fp_out);
  */
}

