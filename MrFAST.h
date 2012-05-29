/*
 * Copyright (c) <2008 - 2020>, University of Washington, Simon Fraser University, Bilkent University
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
 */


/*
  Authors: 
  Farhad Hormozdiari
  Faraz Hach
  Can Alkan
  Emails: 
  farhadh AT uw DOT edu
  fhach AT cs DOT sfu DOT ca
  calkan AT cs DOT bilkent DOT edu DOT tr
*/



#ifndef __MR_FAST__
#define __MR_FAST__

#include "Reads.h"

#define MAP_CHUNKS 15
#define MAX_CIGAR_SIZE 100


// Pair is used to pre-processing and making the read index table
typedef struct
{
  int hv;
  //char hv[50];
  int readNumber;
} Pair;

typedef struct
{
  int hv;
  unsigned int *seqInfo;
} ReadIndexTable;


typedef struct 
{
  int loc;
  char dir;
  int err;
  float score;
  char md[MAX_CIGAR_SIZE];
  char chr[20];
  char cigar[MAX_CIGAR_SIZE];
  int cigarSize;
  int mdSize;
} FullMappingInfo;

typedef struct
{
  int loc;
  char dir;
  int err;
  float score;
  char md[MAX_CIGAR_SIZE];
  char chr[20];
  char cigar[MAX_CIGAR_SIZE];
  int cigarSize;
  int mdSize;
} BestFullMappingInfo;

typedef struct lc
{
  char md[MAP_CHUNKS][MAX_CIGAR_SIZE];
  int mdSize[MAP_CHUNKS];

  char cigar[MAP_CHUNKS][MAX_CIGAR_SIZE];
  int cigarSize[MAP_CHUNKS];

  int err[MAP_CHUNKS];
  int loc[MAP_CHUNKS];
  struct lc *next;
} MappingLocations;

typedef struct inf
{
  int size;
  MappingLocations *next;
} MappingInfo;


typedef struct 
{
  FullMappingInfo *mi;
  int size;
} FullMappingInfoLink;


extern long long			verificationCnt;
extern long long			mappingCnt;
extern long long			mappedSeqCnt;
extern long long			completedSeqCnt;

void initFAST(	Read *seqList,
		int seqListSize,
		int *samplingLocs,
		int samplingLocsSize, 
		char *fileName);

void initVerifiedLocs();
void initLookUpTable();
void initBestMapping();
void initBestConcordantDiscordant(int readNumber);

void finalizeFAST();
void finalizeBestSingleMapping();
void finalizeBestConcordantDiscordant();
void finalizeOEAReads(char *);


int mapAllSingleEndSeq();
//void mapSingleEndSeq(unsigned int *l1, int s1, int readNumber, int readSegment, int direction);
//void mapPairedEndSeqList(unsigned int *l1, int s1, int readNumber, int readSegment, int direction);

void mapPairedEndSeq();

void outputPairedEnd();
void outputPairedEndDiscPP();


void outputPairFullMappingInfo(FILE *fp, int readNumber);
void setPairFullMappingInfo(int readNumber, FullMappingInfo mi1, FullMappingInfo mi2);
void setFullMappingInfo(int readNumber, int loc, int dir, int err, int score, char *md, char * refName, char *cigar);

void outputAllTransChromosomal();
/*
  void outputTransChromosomal(char *fileName1, char *fileName2, FILE * fp_out);
*/

void generateSNPSAM(char *matrix, int matrixLength, char *outputSNP);
void generateCigar(char *matrix, int matrixLength, char *cigar);
void generateCigarFromMD(char *mistmatch, int mismatchLength, char *cigar);

int msfHashVal(char *seq);

int backwardEditDistance2SSE2(char *a, int lena, char *b,int lenb);
int forwardEditDistance2SSE2(char *a, int lena, char *b,int lenb);

int forwardEditDistanceSSE2G(char *a, int lena, char *b,int lenb);
int backwardEditDistanceSSE2G(char *a, int lena, char *b,int lenb);

int forwardEditDistance4SSE2(char *a, int lena, char *b,int lenb);
int backwardEditDistance4SSE2(char *a, int lena, char *b,int lenb);

int forwardEditDistanceSSE2Extension(char *a, int lena, char *b,int lenb);
int backwardEditDistanceSSE2Extension(char *a, int lena, char *b,int lenb);


/***********************************/

int editDistance(int refIndex, char *seq, int seqLength, char *matrix);

int verifySingleEndEditDistance(int refIndex, char *lSeq, int lSeqLength, char *rSeq, int rSeqLength, int segLength, 
				char *matrix, int *map_location, short *seqHashValue);

int verifySingleEndEditDistance2(int refIndex, char *lSeq, int lSeqLength, char *rSeq, int rSeqLength, int segLength, 
				 char *matrix, int *map_location, short *seqHashValue);

int verifySingleEndEditDistance4(int refIndex, char *lSeq, int lSeqLength, char *rSeq, int rSeqLength, int segLength, 
				 char *matrix, int *map_location, short *seqHashValue);

int verifySingleEndEditDistanceExtension(int refIndex, char *lSeq, int lSeqLength, char *rSeq, int rSeqLength, int segLength,
					 char *matrix, int *map_location, short *seqHashValue);

// for fastHASH 
int compareEntrySize (const void *a, const void *b);											// fastHASH()
void mapSingleEndSeq(unsigned int *l1, int s1, int readNumber, int readSegment, int direction,	// fastHASH()
                     int index, key_struct* keys_input, int potential_key_number); 				// fastHASH()
void mapPairEndSeqList(unsigned int *l1, int s1, int readNumber, int readSegment, int direction,// fastHASH()
                       int index, key_struct* keys_input, int potential_key_number); 			// fastHASH()

#endif
