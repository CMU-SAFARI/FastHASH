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
#include "Common.h"
#include "CommandLineParser.h"
#include "Reads.h"
#include "Output.h"
#include "HashTable.h"
#include "MrFAST.h"

char 				*versionNumber = "2.6";			// Current Version
unsigned char		seqFastq;

extern int
main ( int argc, char *argv[] )
{
	if (!parseCommandLine(argc, argv))
	  	return 1;

	/* Determined wrt WINDOW_SIZE */
	configHashTable();
	
	/****************************************************
	 * INDEXING
	 ***************************************************/
	if (indexingMode)
	{
		/********************************
		 * Regular Mode
		 ********************************/
		configHashTable();
		generateHashTable(fileName[0], fileName[1]);
	}
	/****************************************************
	 * SEARCHING
	 ***************************************************/
	else
	{
		Read			*seqList;
		unsigned int	 seqListSize;
		int				 fc;
		int				 samplingLocsSize;
		int				*samplingLocs;
		double			 totalLoadingTime = 0;
		double			 totalMappingTime = 0;
		double			 startTime;
		double			 loadingTime;
		double			 mappingTime;
		double			 lstartTime;
		double			 ppTime			  = 0.0;
		double			 tmpTime;;
		char			*prevGen		  = getMem(CONTIG_NAME_SIZE);
		prevGen[0]						  =	'\0';
		char			*curGen;
		int				 flag;
		double			 maxMem			  =	0;
		char			 fname1[FILE_NAME_LENGTH];
		char			 fname2[FILE_NAME_LENGTH];
		char			 fname3[FILE_NAME_LENGTH];
		char			 fname4[FILE_NAME_LENGTH];
		char			 fname5[FILE_NAME_LENGTH];

		char outputFileName[FILE_NAME_LENGTH];

		
		startTime = getTime();

		/* Loads all reads with all fields of them initialized into seqList, of
		 * size seqListSize. Also handles paired end mode. */
		if (!readAllReads(seqFile1, seqFile2, seqCompressed, &seqFastq, pairedEndMode, &seqList, &seqListSize))
			return 1;

		/* Loads sampling locations for reads. */
		loadSamplingLocations(&samplingLocs, &samplingLocsSize);
		totalLoadingTime += getTime()-startTime;

		/* Simple stuff related to paired end mode. */
		if (pairedEndMode)
	    {
			//Switching to Inferred Size 
			minPairEndedDistance = minPairEndedDistance - SEQ_LENGTH + 2;
			maxPairEndedDistance = maxPairEndedDistance - SEQ_LENGTH + 2;
			if (pairedEndDiscordantMode)
			{
				maxPairEndedDiscordantDistance = maxPairEndedDiscordantDistance - SEQ_LENGTH + 2;
				minPairEndedDiscordantDistance = minPairEndedDiscordantDistance - SEQ_LENGTH + 2;
			}

			/* These files are used for temporary output information. */
			sprintf(fname1, "__%s__1", mappingOutput);
			sprintf(fname2, "__%s__2", mappingOutput);
			sprintf(fname3, "__%s__disc", mappingOutput);
			sprintf(fname4, "__%s__oea1", mappingOutput);
			sprintf(fname5, "__%s__oea2", mappingOutput);
			unlink(fname1);
			unlink(fname2);
			unlink(fname3);
			unlink(fname4);
			unlink(fname5);
	    }

		sprintf(outputFileName, "%s%s",mappingOutputPath , mappingOutput);

		// Preparing output
		initOutput(outputFileName, outCompressed);

		fprintf(stdout, "-----------------------------------------------------------------------------------------------------------\n");
		fprintf(stdout, "| %15s | %15s | %15s | %15s | %15s %15s |\n","Seq. Name","Loading Time", "Mapping Time", "Memory Usage(M)","Total Mappings","Mapped reads");
		fprintf(stdout, "-----------------------------------------------------------------------------------------------------------\n");

		/********************************
		 * Regular Mode
		 ********************************/
		if (!pairedEndMode)
	    {
			/* SSE stuff */
			initLookUpTable();

			/* FullMappingInfo and BestFullMappingInfo only differ in (double
			 * tprob) field which exists in BestFullMappingInfo and does not
			 * exist in FullMappingInfo. */
			if(bestMode)
			  	initBestMapping(seqListSize);

			/* Initialize but do not load hash table. */
			if (!initLoadingHashTable(fileName[1]))
			{
				return 1;
			}

			mappingTime = 0;
			loadingTime = 0;
			prevGen[0]	= '\0';
			flag		= 1;

			do
			{
				/* Note that there may be multiple indexes for different
				 * reference genomes and they might be saved in different binary
				 * files. */
				flag = loadHashTable ( &tmpTime, errThreshold);  			// Reading a fragment
				curGen = getRefGenomeName();

				// First Time
				if (flag && prevGen[0]== '\0')
			    {
					sprintf(prevGen, "%s", curGen);
			    }

				/* Informing. */
				if ( !flag || strcmp(prevGen, curGen)!=0)
			    {
					fprintf(stdout, "| %15s | %15.2f | %15.2f | %15.2f | %15lld %15lld |\n",
						prevGen,loadingTime, mappingTime, maxMem, mappingCnt , mappedSeqCnt);
					fflush(stdout);

					totalMappingTime += mappingTime;
					totalLoadingTime += loadingTime;

					loadingTime = 0;
					mappingTime = 0;
					maxMem = 0;

					if (!flag)
					{
						break;
					}
			   	}
			   	else if (progressRep && mappingTime != 0)
			   	{
					fprintf(stdout, "| %15s | %15.2f | %15.2f | %15.2f | %15lld %15lld |\n",
						prevGen,loadingTime, mappingTime, maxMem, mappingCnt , mappedSeqCnt);
					fflush(stdout);
			  	}

				sprintf(prevGen, "%s", curGen);
				loadingTime += tmpTime;
				//						lstartTime = getTime();

				/* Init reference genome, sampling locations, etc. */
				initFAST(seqList, seqListSize, samplingLocs, samplingLocsSize, fileName[0]);
				
				lstartTime = getTime();

				/* This is where all stuff is done */
				mapAllSingleEndSeq();
				
				mappingTime += getTime() - lstartTime;
				if (maxMem < getMemUsage())
			  	{
					maxMem = getMemUsage();					 
			  	}
			} while (flag);

			if(bestMode)
			  	finalizeBestSingleMapping();
			
			finalizeFAST();
			finalizeLoadingHashTable();
	    }
		else					/* Pairedend Mapping Mode */
	  	{
			/* SSE stuff */
			initLookUpTable();

			/* This is actually best mapping mode output information which is
			 * also used for paired-end mode. */
			if(pairedEndMode)
			  	initBestMapping(seqListSize);

			/* Initialize but do not load hash table. */
			if (!initLoadingHashTable(fileName[1]))
				return 1;
			
			mappingTime = 0;
			loadingTime = 0;
			prevGen[0]	= '\0';
			flag		= 1;

			do
			{
				/* Note that there may be multiple indexes for different
				 * reference genomes and they might be saved in different binary
				 * files. */
				flag   = loadHashTable ( &tmpTime , errThreshold);
				curGen = getRefGenomeName();

				// First Time
				if (flag && prevGen[0]== '\0')
			  	{
					sprintf(prevGen, "%s", curGen);
			  	}

				/* Inform the user. */
				if ( !flag || strcmp(prevGen, curGen)!=0)
			  	{
					// DISCORDANT
					lstartTime = getTime();
					outputPairedEnd(); /* Be careful, dumping data here. */
					mappingTime += getTime() - lstartTime;
					//DISCORDANT			

					fprintf(stdout, "| %15s | %15.2f | %15.2f | %15.2f | %15lld %15lld |\n",
						prevGen,loadingTime, mappingTime, maxMem, mappingCnt , mappedSeqCnt);
					fflush(stdout);

					totalMappingTime += mappingTime;
					totalLoadingTime += loadingTime;

					loadingTime = 0;
					mappingTime = 0;
					maxMem = 0;

					if (!flag)
					{
						break;
					}
			  	}
				else if (progressRep && mappingTime != 0)
			  	{
					fprintf(stdout, "| %15s | %15.2f | %15.2f | %15.2f | %15lld %15lld |\n",
						prevGen,loadingTime, mappingTime, maxMem, mappingCnt , mappedSeqCnt);
					fflush(stdout);
			  	}

				sprintf(prevGen, "%s", curGen);

				loadingTime += tmpTime;
				lstartTime = getTime();

				/* Init reference genome, sampling locations, etc. */
				initFAST(seqList, seqListSize, samplingLocs, samplingLocsSize, fileName[0]);

				/* This is where all stuff is done */
				mapPairedEndSeq();

				mappingTime += getTime() - lstartTime;
				if (maxMem < getMemUsage())
			  	{
					maxMem = getMemUsage();
			  	}

		   	} while (flag);


			/* Need to write OEA and best concordant/discordant information if
			 * it is paired end mode. */
			if(pairedEndMode)
			{
				sprintf(outputFileName, "%s%s",mappingOutputPath , mappingOutput);

				/* Writes the SAM file regarding OEA reads. */
				finalizeOEAReads(outputFileName);

				/* This literally does nothing, it is disabled. */
				outputAllTransChromosomal(transChromosomal);

				/* Writes best concordant/discordant mapping information. A
				 * single entry in the file for each read/1 or read/2. */
				finalizeBestConcordantDiscordant();
			}

			finalizeLoadingHashTable();

			
			if (pairedEndDiscordantMode)
			{
				lstartTime = getTime();

				/* Writes out discordant mappings for VH. */
				outputPairedEndDiscPP();
				ppTime = getTime() - lstartTime;
			}
			finalizeFAST();
	  	} //else

		finalizeOutput();

		fprintf(stdout, "-----------------------------------------------------------------------------------------------------------\n");
		fprintf(stdout, "%19s%16.2f%18.2f\n\n", "Total:",totalLoadingTime, totalMappingTime);
		if (pairedEndDiscordantMode)
	  		fprintf(stdout, "Post Processing Time: %18.2f \n", ppTime);
		fprintf(stdout, "%-30s%10.2f\n","Total Time:", totalMappingTime+totalLoadingTime);
		fprintf(stdout, "%-30s%10d\n","Total No. of Reads:", seqListSize);
		fprintf(stdout, "%-30s%10lld\n","Total No. of Mappings:", mappingCnt);
		fprintf(stdout, "%-30s%10.0f\n\n","Avg No. of locations verified:", ceil((float)verificationCnt/seqListSize));
		
		int cof = (pairedEndMode)?2:1;

		/* maxHits used here! */
		if (progressRep && maxHits != 0)
	  	{
			int frequency[maxHits+1];
			int i;
			for ( i=0 ; i <= maxHits; i++)
			{
				frequency[i] = 0;
			}

			for (fc = 0; fc < seqListSize; fc++)
			{
				frequency[(int)(*(seqList[fc*cof].hits))]++;
			}
			frequency[maxHits] = completedSeqCnt;
			for ( i=0 ; i <= maxHits; i++)
			{
				fprintf(stdout, "%-30s%10d%10d%10.2f%%\n","Reads Mapped to ", i, frequency[i], 100*(float)frequency[i]/(float)seqListSize);
			}
	  	}

		finalizeReads(unmappedOutput);
		freeMem(prevGen, CONTIG_NAME_SIZE);
	}	/* Searching mode */

	return 0;
}
