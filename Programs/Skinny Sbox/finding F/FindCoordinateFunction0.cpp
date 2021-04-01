/*
* -----------------------------------------------------------------
* COMPANY : Ruhr University Bochum
* AUTHOR  : Amir Moradi amir.moradi@rub.de Aein Rezaei Shahmirzadi aein.rezaeishahmirzadi@rub.de
* DOCUMENT: "Second-Order SCA Security with almost no Fresh Randomness" TCHES 2021, Issue 3
* -----------------------------------------------------------------
*
* Copyright c 2021, Amir Moradi, Aein Rezaei Shahmirzadi
*
* All rights reserved.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTERS BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* INCLUDING NEGLIGENCE OR OTHERWISE ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
* Please see LICENSE and README for license and further instructions.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <omp.h>

#define	MaxNumberOfTerms 1000000000L

unsigned char   *HW;
unsigned int    *NumberOfIndexesInHW;
unsigned long  **IndexesInHW;

#define ONEL ((unsigned __int64)1)

void FillANFTables(unsigned int NumberOfCells)
{
	unsigned int	i, Bit;
	unsigned char	NumberOfBits = log((double)NumberOfCells) / log((double)2);

	HW = (unsigned char *)malloc(NumberOfCells * sizeof(unsigned char));

	NumberOfIndexesInHW = (unsigned int *)malloc((NumberOfBits + 1) * sizeof(unsigned int));
	IndexesInHW = (unsigned long **)malloc((NumberOfBits + 1) * sizeof(unsigned long *));

	for (Bit = 0; Bit < (NumberOfBits + 1); Bit++)
	{
		IndexesInHW[Bit] = (unsigned long *)malloc((NumberOfCells / 2) * sizeof(unsigned long));
		NumberOfIndexesInHW[Bit] = 0;
	}

	for (i = 0; i < NumberOfCells; i++)
	{
		HW[i] = 0;

		for (Bit = 0; Bit < NumberOfBits; Bit++)
			HW[i] += ((i & (1L << Bit)) >> Bit);

		IndexesInHW[HW[i]][NumberOfIndexesInHW[HW[i]]] = i;
		NumberOfIndexesInHW[HW[i]]++;
	}
}


struct FunctionStruct
{
	unsigned int		NumberOfTerms;
	unsigned __int64   *Term;
};

void DefineFunction(FunctionStruct *F, unsigned short NumberOfTerms)
{
	F->Term = (unsigned __int64 *)malloc(NumberOfTerms * sizeof(unsigned __int64));
}

void FreeFunction(FunctionStruct *F)
{
	free(F->Term);
}

unsigned char EvalFunction(FunctionStruct *Function, unsigned __int64 Input)
{
	unsigned char	Res = 0;
	unsigned int	i;

	for (i = 0; i < Function->NumberOfTerms; i++)
		Res ^= ((Function->Term[i] & Input) == Function->Term[i]);

	return (Res);
}

void AddTermtoFunction(FunctionStruct *Function, unsigned __int64 Symbol)
{
	Function->Term[Function->NumberOfTerms] = Symbol;
	Function->NumberOfTerms++;
}

void MakeANF(unsigned char* Table, unsigned int NumberOfElements, FunctionStruct* Functions, unsigned char NumberOfOutputBits, unsigned char Invert)
{
	unsigned char   NumberOfBits = log((double)NumberOfElements) / log((double)2);
	unsigned char	Bit, HWIndex;
	unsigned int	i;

	for (Bit = 0; Bit < NumberOfOutputBits; Bit++)
	{
		Functions[Bit].NumberOfTerms = 0;

		for (HWIndex = 0; HWIndex < (NumberOfBits + 1); HWIndex++)
			for (i = 0; i < NumberOfIndexesInHW[HWIndex]; i++)
				if (IndexesInHW[HWIndex][i] < NumberOfElements)
					if (EvalFunction(&Functions[Bit], IndexesInHW[HWIndex][i]) != (((Table[IndexesInHW[HWIndex][i]] ^ Invert) >> Bit) & 1))
						AddTermtoFunction(&Functions[Bit], IndexesInHW[HWIndex][i]);
	}
}

unsigned short* Term2Index;
unsigned short* Index2Term;

void FillANFTermTables()
{
	unsigned short i;
	unsigned char  c;

	Term2Index = (unsigned short*)malloc(512 * sizeof(unsigned short));
	Index2Term = (unsigned short*)malloc(512 * sizeof(unsigned short));

	// input order: (MSB) a1, a2, a3, b1, b2, b3, c1, c2, c3(LBS)

	c = 0;
	for (i = 1; i < 512; i <<= 1)  // linear terms
	{
		Term2Index[i] = c;
		Index2Term[c] = i;
		c++;
	}
}

unsigned char EvalFunctionSpecial(uint64_t Func, unsigned short Input)
{
	unsigned char	Res = 0;
	unsigned int	i;

	for (i = 0; Func; Func >>= 1)
	{
		if (Func & 1)
			Res ^= ((Index2Term[i] & Input) == Index2Term[i]);
		i++;
	}

	return (Res);
}

uint64_t MakeANFSpecial(unsigned char* Table, unsigned char NumberOfBits)
{
	unsigned int	NumberOfElements = (1 << NumberOfBits);
	unsigned char	HWIndex;
	unsigned char	i;
	uint64_t		Res = 0;
	uint64_t		ONE = 1;

	for (HWIndex = 1; HWIndex <= 1; HWIndex++)
		for (i = 0; i < NumberOfIndexesInHW[HWIndex]; i++)
			if (IndexesInHW[HWIndex][i] < NumberOfElements)
			{
				if (EvalFunctionSpecial(Res, IndexesInHW[HWIndex][i]) != Table[IndexesInHW[HWIndex][i]])
					Res |= (ONE << Term2Index[IndexesInHW[HWIndex][i]]);
			}
			else
				break;

	return Res;
}


unsigned char AlgDegree(FunctionStruct Function)
{
	unsigned int		i;
	unsigned char		Degree;

	Degree = 1;
	for (i = 0; i < Function.NumberOfTerms; i++)
		if (HW[Function.Term[i]] > Degree)
			Degree = HW[Function.Term[i]];

	return(Degree);
}

void SPrintFunction(FunctionStruct Function, char* ResStr, char InputStr[][5])
{
	unsigned int		i, j;
	unsigned __int64	ONE;

	ResStr[0] = 0;

	for (i = 0; i < Function.NumberOfTerms; i++)
	{
		if (i)
			strcat(ResStr, " + ");

		if (Function.Term[i])
		{
			ONE = 1;

			for (j = 0; j < 64; j++)
			{
				if (Function.Term[i] & ONE)
					sprintf(ResStr, "%s%s", ResStr, InputStr[j]);

				ONE <<= 1;
			}
		}
		else
			strcat(ResStr, "1");
	}
}

//---------------------------------------------------------------------------

int FillF3to1(unsigned char** F3to1[3], unsigned char** &F3to1_full)
{
	int				NumberOfF3to1;
	unsigned int	i, j;
	char			index;
	char			a, b, c1, c2, c3;

	for (index = 0; index < 3; index++)
		if (F3to1[index] == NULL)
			F3to1[index] = (unsigned char**)calloc(16, sizeof(unsigned char*)); // to fill every cell with NULL

	if (F3to1_full == NULL)
		F3to1_full = (unsigned char**)calloc(16, sizeof(unsigned char*)); // to fill every cell with NULL

	for (NumberOfF3to1 = 0; NumberOfF3to1 < 16; NumberOfF3to1++)
	{
		if (F3to1_full[NumberOfF3to1] == NULL)
			F3to1_full[NumberOfF3to1] = (unsigned char*)malloc(32 * sizeof(unsigned char));

		memset(F3to1_full[NumberOfF3to1], 0, 32);

		for (j = 0;j < 32;j++)
		{
			a = (j >> 0) & 1;
			b = (j >> 1) & 1;
			c1 = (j >> 2) & 1;
			c2 = (j >> 3) & 1;
			c3 = (j >> 4) & 1;

			F3to1_full[NumberOfF3to1][j] |= a;
			F3to1_full[NumberOfF3to1][j] |= b << 1;

			if (((NumberOfF3to1 >> 2) & 3) == 1)
				F3to1_full[NumberOfF3to1][j] |= c1 << 2;

			if (((NumberOfF3to1 >> 2) & 3) == 2)
				F3to1_full[NumberOfF3to1][j] |= c2 << 2;

			if (((NumberOfF3to1 >> 2) & 3) == 3)
				F3to1_full[NumberOfF3to1][j] |= c3 << 2;
		}

		//*************************************//

		for (index = 0; index < 3; index++)
		{
			if (F3to1[index][NumberOfF3to1] == NULL)
				F3to1[index][NumberOfF3to1] = (unsigned char*)malloc(32 * sizeof(unsigned char));

			memset(F3to1[index][NumberOfF3to1], 0, 32);

			for (j = 0;j < 32;j++)
			{
				a = (j >> 0) & 1;
				b = (j >> 1) & 1;
				c1 = (j >> 2) & 1;
				c2 = (j >> 3) & 1;
				c3 = (j >> 4) & 1;

				F3to1[index][NumberOfF3to1][j] ^= a & b;

				if ((NumberOfF3to1 >> 0) & 1)
					F3to1[index][NumberOfF3to1][j] ^= a;

				if ((NumberOfF3to1 >> 1) & 1)
					F3to1[index][NumberOfF3to1][j] ^= b;

				if (((NumberOfF3to1 >> 2) & 3) == 1)
					F3to1[index][NumberOfF3to1][j] ^= c1;

				if (((NumberOfF3to1 >> 2) & 3) == 2)
					F3to1[index][NumberOfF3to1][j] ^= c2;

				if (((NumberOfF3to1 >> 2) & 3) == 3)
					F3to1[index][NumberOfF3to1][j] ^= c3;

				F3to1[index][NumberOfF3to1][j] <<= index;
			}
		}
	}

	return(NumberOfF3to1);
}

struct ToupleStruct
{
	unsigned short	i[3];
	uint64_t		ANF;
};

void FillTableIndexes(unsigned char**& InputTableIndex, unsigned short*& OutputTableIndex, unsigned char*& Unmasking, unsigned char*& UnmaskedInputTable, unsigned char**& ExtraCheckTable)
{
	char			index;
	unsigned short  a, a1, a2, a3;
	unsigned short	b, b1, b2, b3;
	unsigned short	c, c1, c2, c3;
	unsigned short	abc;
	unsigned short	a1b1c1;
	unsigned short	a2b2c2;
	unsigned short	a3b3c3;
	unsigned short	Masked_InputIndex;
	unsigned short	OutIndex;
	unsigned short	j, k;

	InputTableIndex = (unsigned char**)malloc(9 * sizeof(unsigned char*));
	for (index = 0; index < 9; index++)
		InputTableIndex[index] = (unsigned char*)malloc(512 * sizeof(unsigned char));

	OutputTableIndex = (unsigned short*)malloc(512 * sizeof(unsigned short));
	UnmaskedInputTable = (unsigned char*)malloc(512 * sizeof(unsigned char));

	ExtraCheckTable = (unsigned char**)malloc(12 * sizeof(unsigned char*));
	for (index = 0; index < 12; index++)
		ExtraCheckTable[index] = (unsigned char*)malloc(512 * sizeof(unsigned char));

	Masked_InputIndex = 0;
	for (abc = 0; abc < 8; abc++)
	{
		a = (abc >> 0) & 1;
		b = (abc >> 1) & 1;
		c = (abc >> 2) & 1;

		for (a1b1c1 = 0; a1b1c1 < 8; a1b1c1++)
		{
			a1 = (a1b1c1 >> 0) & 1;
			b1 = (a1b1c1 >> 1) & 1;
			c1 = (a1b1c1 >> 2) & 1;

			for (a2b2c2 = 0; a2b2c2 < 8; a2b2c2++)
			{
				a3b3c3 = abc ^ a1b1c1 ^ a2b2c2;

				a2 = (a2b2c2 >> 0) & 1;
				b2 = (a2b2c2 >> 1) & 1;
				c2 = (a2b2c2 >> 2) & 1;

				a3 = (a3b3c3 >> 0) & 1;
				b3 = (a3b3c3 >> 1) & 1;
				c3 = (a3b3c3 >> 2) & 1;

				InputTableIndex[0][Masked_InputIndex] = (c3 << 4) | (c2 << 3) | (c1 << 2) | (b1 << 1) | (a1 << 0);
				InputTableIndex[1][Masked_InputIndex] = (c3 << 4) | (c2 << 3) | (c1 << 2) | (b2 << 1) | (a1 << 0);
				InputTableIndex[2][Masked_InputIndex] = (c3 << 4) | (c2 << 3) | (c1 << 2) | (b3 << 1) | (a1 << 0);

				InputTableIndex[3][Masked_InputIndex] = (c3 << 4) | (c2 << 3) | (c1 << 2) | (b1 << 1) | (a2 << 0);
				InputTableIndex[4][Masked_InputIndex] = (c3 << 4) | (c2 << 3) | (c1 << 2) | (b2 << 1) | (a2 << 0);
				InputTableIndex[5][Masked_InputIndex] = (c3 << 4) | (c2 << 3) | (c1 << 2) | (b3 << 1) | (a2 << 0);

				InputTableIndex[6][Masked_InputIndex] = (c3 << 4) | (c2 << 3) | (c1 << 2) | (b1 << 1) | (a3 << 0);
				InputTableIndex[7][Masked_InputIndex] = (c3 << 4) | (c2 << 3) | (c1 << 2) | (b2 << 1) | (a3 << 0);
				InputTableIndex[8][Masked_InputIndex] = (c3 << 4) | (c2 << 3) | (c1 << 2) | (b3 << 1) | (a3 << 0);

				OutIndex = (a1 << 8) | (a2 << 7) | (a3 << 6) | (b1 << 5) | (b2 << 4) | (b3 << 3) | (c1 << 2) | (c2 << 1) | (c3 << 0);
				OutputTableIndex[Masked_InputIndex] = OutIndex;

				UnmaskedInputTable[Masked_InputIndex] = (c << 2) | (b << 1) | (a << 0);

				ExtraCheckTable[0][Masked_InputIndex] = (c1 << 1) | (a1 << 0);
				ExtraCheckTable[1][Masked_InputIndex] = (c2 << 1) | (a1 << 0);
				ExtraCheckTable[2][Masked_InputIndex] = (c3 << 1) | (a1 << 0);
				ExtraCheckTable[3][Masked_InputIndex] = (c1 << 1) | (a2 << 0);
				ExtraCheckTable[4][Masked_InputIndex] = (c2 << 1) | (a2 << 0);
				ExtraCheckTable[5][Masked_InputIndex] = (c3 << 1) | (a2 << 0);
				ExtraCheckTable[6][Masked_InputIndex] = (c1 << 1) | (a3 << 0);
				ExtraCheckTable[7][Masked_InputIndex] = (c2 << 1) | (a3 << 0);
				ExtraCheckTable[8][Masked_InputIndex] = (c3 << 1) | (a3 << 0);

				ExtraCheckTable[9][Masked_InputIndex] =  (c1 << 1) | (b1 << 0);
				ExtraCheckTable[10][Masked_InputIndex] = (c2 << 1) | (b2 << 0);
				ExtraCheckTable[11][Masked_InputIndex] = (c3 << 1) | (b3 << 0);
				Masked_InputIndex++;
			}
		}
	}

	Unmasking = (unsigned char*)calloc(256, sizeof(unsigned char));

	for (j = 0; j < 256; j++)
		for (k = 0; k < 8; k++)
			Unmasking[j] ^= (j & (1 << k)) ? 1 : 0;
}


void MakeLargeTables(unsigned char** F3to1[3], unsigned char** F3to1_full, int NumberOfF3to1, char Type[9],
	unsigned char** InputTableIndex, unsigned short* OutputTableIndex, unsigned char* Unmasking, unsigned char** ExtraCheckTable,
	ToupleStruct*& FTuples1, ToupleStruct*& FTuples2, ToupleStruct*& FTuples3, unsigned int NumberOfFTuples[3],
	unsigned int**& FTouple3InANF, unsigned short*& NumberOfANF3)
{
	char			index;
	unsigned short	Masked_InputIndex;
	short 			i_last;
	unsigned short	i[4] = { 0 };
	unsigned char	Masked_Output;
	unsigned char	Unmasked_Output;
	unsigned char	UniformityCounter[2];
	unsigned char	Distribution2Old[4];
	unsigned char	Distribution3Old[8];
	unsigned char	Distribution2[4];
	unsigned char	Distribution3[8];
	unsigned char	DisOld[15][64];
	unsigned char	Dis[15][64];
	unsigned char	Uniform;
	unsigned char   full[3];
	unsigned char	MadeTable[512];
	uint64_t		MadeANF;
	unsigned short  j;
	char			k;
	unsigned int* TempBuff;
	unsigned short  Counter;

	if (FTuples1 == NULL)
		FTuples1 = (ToupleStruct*)malloc(4096 * sizeof(ToupleStruct)); // at most can be 4096

	Counter = 0;
	NumberOfFTuples[0] = 0;
	#pragma omp parallel for schedule(guided) private(i, index, MadeTable, MadeANF, UniformityCounter, Distribution2Old, Distribution3Old, Distribution2, Distribution3, full, DisOld, Dis, Masked_InputIndex, Masked_Output, Unmasked_Output, Uniform, j, k)
	for (i_last = 0; i_last < NumberOfF3to1; i_last++)
	{
		i[0] = 0;
		i[1] = 0;
		i[2] = i_last;

		while (i[2] == i_last)
		{
			UniformityCounter[0] = 0;
			UniformityCounter[1] = 0;
			Uniform = 1;
			memset(Distribution2, 0, 4);
			memset(Distribution3, 0, 8);
			for (k = 0; k < 15; k++)
				memset(Dis[k], 0, 64);

			for (Masked_InputIndex = 0; Masked_InputIndex < 512; Masked_InputIndex++)
			{
				Masked_Output = 0;
				for (index = 0; index < 3; index++)
				{
					Masked_Output |= F3to1[index][i[index]][InputTableIndex[Type[index]][Masked_InputIndex]];
					full[index] = F3to1_full[i[index]][InputTableIndex[Type[index]][Masked_InputIndex]];
				}

				Unmasked_Output = Unmasking[Masked_Output];
				MadeTable[OutputTableIndex[Masked_InputIndex]] = Unmasked_Output;
				UniformityCounter[Unmasked_Output]++;

				Distribution2[Masked_Output >> 1]++;
				Distribution3[Masked_Output >> 0]++;

				Dis[0][(full[0] << 3) | Masked_Output]++;
				Dis[1][(full[1] << 3) | Masked_Output]++;
				Dis[2][(full[2] << 3) | Masked_Output]++;

				Dis[3][(ExtraCheckTable[0][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[4][(ExtraCheckTable[1][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[5][(ExtraCheckTable[2][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[6][(ExtraCheckTable[3][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[7][(ExtraCheckTable[4][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[8][(ExtraCheckTable[5][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[9][(ExtraCheckTable[6][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[10][(ExtraCheckTable[7][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[11][(ExtraCheckTable[8][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[12][(ExtraCheckTable[9][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[13][(ExtraCheckTable[10][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[14][(ExtraCheckTable[11][Masked_InputIndex] << 3) | Masked_Output]++;

				if ((Masked_InputIndex & 0x3F) == 0x3F) // every 64
				{
					Uniform &= (UniformityCounter[0] == UniformityCounter[1]);

					UniformityCounter[0] = 0;
					UniformityCounter[1] = 0;

					//---------

					if (Masked_InputIndex == 0x3F) // the first one
					{
						memcpy(Distribution2Old, Distribution2, 4);
						memcpy(Distribution3Old, Distribution3, 8);

						for (k = 0; k < 15; k++)
							memcpy(DisOld[k], Dis[k], 64);
					}
					else
					{
						for (j = 0; j < 4; j++)
							if (Distribution2[j] != Distribution2Old[j])
							{
								i[0] = NumberOfF3to1;
								break;
							}

						if (j < 4)
							break;

						for (j = 0; j < 8; j++)
							if (Distribution3[j] != Distribution3Old[j])
								break;

						if (j < 8)
							break;

						//***********************//

						for (j = 0; j < 64; j++)
						{
							for (k = 0; k < 15; k++)
								if (Dis[k][j] != DisOld[k][j])
									break;

							if (k < 15)
								break;
						}

						if (j < 64)
							break;
					}

					memset(Distribution2, 0, 4);
					memset(Distribution3, 0, 8);

					for (k = 0; k < 15; k++)
						memset(Dis[k], 0, 64);
				}
			}

			if (Uniform & (Masked_InputIndex == 512))
			{
				#pragma omp critical (found1)
				{
					FTuples1[NumberOfFTuples[0]].i[0] = i[0];
					FTuples1[NumberOfFTuples[0]].i[1] = i[1];
					FTuples1[NumberOfFTuples[0]].i[2] = i[2];
					FTuples1[NumberOfFTuples[0]].ANF = MakeANFSpecial(MadeTable, 9);
					NumberOfFTuples[0]++;
				}
			}

			j = 0;
			do
			{
				if ((++i[j] >= NumberOfF3to1) & (j < 3))
					i[j] = 0;
			} while (i[j++] == 0);
		}

		#pragma omp critical (print)
		{
			Counter++;
			for (j = 0; j < 9; j++)
				printf("%d", Type[j]);
			printf(" Touple1 %d / %d %d\n", Counter, NumberOfF3to1, NumberOfFTuples[0]);
		}
	}

	//-------------------------------

	if (FTuples2 == NULL)
		FTuples2 = (ToupleStruct*)malloc(4096 * sizeof(ToupleStruct)); // at most can be 4096

	Counter = 0;
	NumberOfFTuples[1] = 0;
	#pragma omp parallel for schedule(guided) private(i, index, MadeTable, MadeANF, UniformityCounter, Distribution2Old, Distribution3Old, Distribution2, Distribution3, full, DisOld, Dis, Masked_InputIndex, Masked_Output, Unmasked_Output, Uniform, j, k)
	for (i_last = 0; i_last < NumberOfF3to1; i_last++)
	{
		i[0] = 0;
		i[1] = 0;
		i[2] = i_last;

		while (i[2] == i_last)
		{
			UniformityCounter[0] = 0;
			UniformityCounter[1] = 0;
			Uniform = 1;
			memset(Distribution2, 0, 4);
			memset(Distribution3, 0, 8);
			for (k = 0; k < 15; k++)
				memset(Dis[k], 0, 64);

			for (Masked_InputIndex = 0; Masked_InputIndex < 512; Masked_InputIndex++)
			{
				Masked_Output = 0;
				for (index = 0; index < 3; index++)
				{
					Masked_Output |= F3to1[index][i[index]][InputTableIndex[Type[index + 3]][Masked_InputIndex]];
					full[index] = F3to1_full[i[index]][InputTableIndex[Type[index + 3]][Masked_InputIndex]];
				}

				Unmasked_Output = Unmasking[Masked_Output];
				MadeTable[OutputTableIndex[Masked_InputIndex]] = Unmasked_Output;
				UniformityCounter[Unmasked_Output]++;

				Distribution2[Masked_Output >> 1]++;
				Distribution3[Masked_Output >> 0]++;

				Dis[0][(full[0] << 3) | Masked_Output]++;
				Dis[1][(full[1] << 3) | Masked_Output]++;
				Dis[2][(full[2] << 3) | Masked_Output]++;

				Dis[3][(ExtraCheckTable[0][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[4][(ExtraCheckTable[1][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[5][(ExtraCheckTable[2][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[6][(ExtraCheckTable[3][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[7][(ExtraCheckTable[4][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[8][(ExtraCheckTable[5][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[9][(ExtraCheckTable[6][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[10][(ExtraCheckTable[7][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[11][(ExtraCheckTable[8][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[12][(ExtraCheckTable[9][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[13][(ExtraCheckTable[10][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[14][(ExtraCheckTable[11][Masked_InputIndex] << 3) | Masked_Output]++;

				if ((Masked_InputIndex & 0x3F) == 0x3F) // every 64
				{
					Uniform &= (UniformityCounter[0] == UniformityCounter[1]);

					UniformityCounter[0] = 0;
					UniformityCounter[1] = 0;

					//---------

					if (Masked_InputIndex == 0x3F) // the first one
					{
						memcpy(Distribution2Old, Distribution2, 4);
						memcpy(Distribution3Old, Distribution3, 8);

						for (k = 0; k < 15; k++)
							memcpy(DisOld[k], Dis[k], 64);
					}
					else
					{
						for (j = 0; j < 4; j++)
							if (Distribution2[j] != Distribution2Old[j])
							{
								i[0] = NumberOfF3to1;
								break;
							}

						if (j < 4)
							break;

						for (j = 0; j < 8; j++)
							if (Distribution3[j] != Distribution3Old[j])
								break;

						if (j < 8)
							break;

						//***********************//

						for (j = 0; j < 64; j++)
						{
							for (k = 0; k < 15; k++)
								if (Dis[k][j] != DisOld[k][j])
									break;

							if (k < 15)
								break;
						}

						if (j < 64)
							break;
					}

					memset(Distribution2, 0, 4);
					memset(Distribution3, 0, 8);

					for (k = 0; k < 15; k++)
						memset(Dis[k], 0, 64);
				}
			}

			if (Uniform & (Masked_InputIndex == 512))
			{
				#pragma omp critical (found2)
				{
					FTuples2[NumberOfFTuples[1]].i[0] = i[0];
					FTuples2[NumberOfFTuples[1]].i[1] = i[1];
					FTuples2[NumberOfFTuples[1]].i[2] = i[2];
					MadeANF = MakeANFSpecial(MadeTable, 9);
					FTuples2[NumberOfFTuples[1]].ANF = MadeANF;
					NumberOfFTuples[1]++;
				}
			}

			j = 0;
			do
			{
				if ((++i[j] >= NumberOfF3to1) & (j < 3))
					i[j] = 0;
			} while (i[j++] == 0);
		}

		#pragma omp critical (print)
		{
			Counter++;
			for (j = 0; j < 9; j++)
				printf("%d", Type[j]);
			printf(" Touple2 %d / %d %d\n", Counter, NumberOfF3to1, NumberOfFTuples[1]);
		}
	}

	//-------------------------------

	if (FTuples3 == NULL)
		FTuples3 = (ToupleStruct*)malloc(4096 * sizeof(ToupleStruct)); // at most can be 4096

	if (NumberOfANF3 == NULL)
		NumberOfANF3 = (unsigned short*)calloc((1 << 9), sizeof(unsigned short));
	else
		for (unsigned int i = 0; i < (1 << 9); i++)
		{
			free(FTouple3InANF[i]);
			FTouple3InANF[i] = 0;
			NumberOfANF3[i] = 0;
		}

	if (FTouple3InANF == NULL)
		FTouple3InANF = (unsigned int**)calloc(4096, sizeof(unsigned int*));

	Counter = 0;
	NumberOfFTuples[2] = 0;
	#pragma omp parallel for schedule(guided) private(i, index, MadeTable, MadeANF, UniformityCounter, Distribution2Old, Distribution3Old, Distribution2, Distribution3, full, DisOld, Dis, Masked_InputIndex, Masked_Output, Unmasked_Output, Uniform, j, k)
	for (i_last = 0; i_last < NumberOfF3to1; i_last++)
	{
		i[0] = 0;
		i[1] = 0;
		i[2] = i_last;

		while (i[2] == i_last)
		{
			UniformityCounter[0] = 0;
			UniformityCounter[1] = 0;
			Uniform = 1;
			memset(Distribution2, 0, 4);
			memset(Distribution3, 0, 8);
			for (k = 0; k < 15; k++)
				memset(Dis[k], 0, 64);

			for (Masked_InputIndex = 0; Masked_InputIndex < 512; Masked_InputIndex++)
			{
				Masked_Output = 0;
				for (index = 0; index < 3; index++)
				{
					Masked_Output |= F3to1[index][i[index]][InputTableIndex[Type[index + 6]][Masked_InputIndex]];
					full[index] = F3to1_full[i[index]][InputTableIndex[Type[index + 6]][Masked_InputIndex]];
				}

				Unmasked_Output = Unmasking[Masked_Output];
				MadeTable[OutputTableIndex[Masked_InputIndex]] = Unmasked_Output;
				UniformityCounter[Unmasked_Output]++;

				Distribution2[Masked_Output >> 1]++;
				Distribution3[Masked_Output >> 0]++;

				Dis[0][(full[0] << 3) | Masked_Output]++;
				Dis[1][(full[1] << 3) | Masked_Output]++;
				Dis[2][(full[2] << 3) | Masked_Output]++;

				Dis[3][(ExtraCheckTable[0][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[4][(ExtraCheckTable[1][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[5][(ExtraCheckTable[2][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[6][(ExtraCheckTable[3][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[7][(ExtraCheckTable[4][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[8][(ExtraCheckTable[5][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[9][(ExtraCheckTable[6][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[10][(ExtraCheckTable[7][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[11][(ExtraCheckTable[8][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[12][(ExtraCheckTable[9][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[13][(ExtraCheckTable[10][Masked_InputIndex] << 3) | Masked_Output]++;
				Dis[14][(ExtraCheckTable[11][Masked_InputIndex] << 3) | Masked_Output]++;

				if ((Masked_InputIndex & 0x3F) == 0x3F) // every 64
				{
					Uniform &= (UniformityCounter[0] == UniformityCounter[1]);

					UniformityCounter[0] = 0;
					UniformityCounter[1] = 0;

					//---------

					if (Masked_InputIndex == 0x3F) // the first one
					{
						memcpy(Distribution2Old, Distribution2, 4);
						memcpy(Distribution3Old, Distribution3, 8);

						for (k = 0; k < 15; k++)
							memcpy(DisOld[k], Dis[k], 64);
					}
					else
					{
						for (j = 0; j < 4; j++)
							if (Distribution2[j] != Distribution2Old[j])
							{
								i[0] = NumberOfF3to1;
								break;
							}

						if (j < 4)
							break;

						for (j = 0; j < 8; j++)
							if (Distribution3[j] != Distribution3Old[j])
								break;

						if (j < 8)
							break;

						//***********************//

						for (j = 0; j < 64; j++)
						{
							for (k = 0; k < 15; k++)
								if (Dis[k][j] != DisOld[k][j])
									break;

							if (k < 15)
								break;
						}

						if (j < 64)
							break;
					}

					memset(Distribution2, 0, 4);
					memset(Distribution3, 0, 8);

					for (k = 0; k < 15; k++)
						memset(Dis[k], 0, 64);
				}
			}

			if (Uniform & (Masked_InputIndex == 512))
			{
				#pragma omp critical (found3)
				{
					FTuples3[NumberOfFTuples[2]].i[0] = i[0];
					FTuples3[NumberOfFTuples[2]].i[1] = i[1];
					FTuples3[NumberOfFTuples[2]].i[2] = i[2];

					MadeANF = MakeANFSpecial(MadeTable, 9);
					FTuples3[NumberOfFTuples[2]].ANF = MadeANF;
					TempBuff = (unsigned int*)malloc((NumberOfANF3[MadeANF] + 1) * sizeof(unsigned int));
					memcpy(TempBuff, FTouple3InANF[MadeANF], NumberOfANF3[MadeANF] * sizeof(unsigned int));
					free(FTouple3InANF[MadeANF]);
					FTouple3InANF[MadeANF] = TempBuff;

					FTouple3InANF[MadeANF][NumberOfANF3[MadeANF]] = NumberOfFTuples[2];
					NumberOfANF3[MadeANF]++;
					NumberOfFTuples[2]++;
				}
			}

			j = 0;
			do
			{
				if ((++i[j] >= NumberOfF3to1) & (j < 3))
					i[j] = 0;
			} while (i[j++] == 0);
		}

		#pragma omp critical (print)
		{
			Counter++;
			for (j = 0; j < 9; j++)
				printf("%d", Type[j]);
			printf(" Touple3 %d / %d %d\n", Counter, NumberOfF3to1, NumberOfFTuples[2]);
		}
	}
}

void FindCombination(char* FileName, unsigned short* OutputTableIndex,
	char SelectedVarIndexes[3], char SelectedVars[3][5], char invert)
{
	unsigned short	Masked_InputIndex;
	int				i;

	unsigned char	Masked_Output1;
	unsigned char	Masked_Output2;
	unsigned char	Masked_Output3;
	unsigned char	MadeTable1[512];
	unsigned char	MadeTable2[512];
	unsigned char	MadeTable3[512];
	unsigned short*	TableIndex;
	unsigned short  j;
	short			abcd;
	short			a1b1c1d1;
	short			a2b2c2d2;
	short			v[4];
	short			v1[4];
	short			v2[4];
	short			InputIndex;
	char			index;
	FILE*			F;
	char			InputStr[9][5] = {
		{ "x1"},
		{ "x1"},
		{ "x1"},
		{ "x2"},
		{ "x2"},
		{ "x2"},
		{ "x3"},
		{ "x3"},
		{ "x3"} };

	for (i = 0; i < 9; i++)
	{
		InputStr[i][0] = SelectedVars[0][0];
	}

	TableIndex = (unsigned short*)malloc((1 << 12) * sizeof(short));

	InputIndex = 0;
	for (abcd=0;abcd<16;abcd++)
		for (a1b1c1d1 = 0;a1b1c1d1 < 16;a1b1c1d1++)
			for (a2b2c2d2 = 0;a2b2c2d2 < 16;a2b2c2d2++)
			{
				for (i = 0;i < 4;i++)
				{
					v[i]  = (abcd >> i) & 1;
					v1[i] = (a1b1c1d1 >> i) & 1;
					v2[i] = (a2b2c2d2 >> i) & 1;
				}

				TableIndex[InputIndex] = 0;

				for (i = 0;i < 3;i++)
				{
					TableIndex[InputIndex] <<= 1;
					TableIndex[InputIndex] |= (i == 2 ? v[SelectedVarIndexes[2-i]] : 0);
				}

				for (i = 0;i < 3;i++)
				{
					TableIndex[InputIndex] <<= 1;
					TableIndex[InputIndex] |= (i == 2 ? v1[SelectedVarIndexes[2-i]] : 0);
				}

				for (i = 0;i < 3;i++)
				{
					TableIndex[InputIndex] <<= 1;
					TableIndex[InputIndex] |= (i == 2 ? v2[SelectedVarIndexes[2-i]] : 0);
				}

				InputIndex++;
			}

	/******************************************/


	for (Masked_InputIndex = 0; Masked_InputIndex < 512; Masked_InputIndex++)
	{
		Masked_Output1 = (OutputTableIndex[Masked_InputIndex] >> 8) & 1;
		Masked_Output2 = (OutputTableIndex[Masked_InputIndex] >> 7) & 1;
		Masked_Output3 = (OutputTableIndex[Masked_InputIndex] >> 6) & 1;

		MadeTable1[Masked_InputIndex] = Masked_Output1;
		MadeTable2[Masked_InputIndex] = Masked_Output2;
		MadeTable3[Masked_InputIndex] = Masked_Output3;
	}

	F = fopen(FileName, "at");
	fprintf(F, "012345678");

	fprintf(F, ", 0, ");
	for (index = 0; index < 3; index++)
		fprintf(F, "%d, %s%s, ", 0, index ? "0" : InputStr[index], (!index) & invert ? " + 1" : "");

	for (j = 0; j < (1 << 12); j++)
		fprintf(F, "%x", MadeTable1[TableIndex[j]]);

	//----------------

	fprintf(F, ", 0, ");
	for (index = 0; index < 3; index++)
		fprintf(F, "%d, %s, ", 0, index ? "0" : InputStr[index + 3]);

	for (j = 0; j < (1 << 12); j++)
		fprintf(F, "%x", MadeTable2[TableIndex[j]]);

	//----------------

	fprintf(F, ", 0, ");
	for (index = 0; index < 3; index++)
		fprintf(F, "%d, %s, ", 0, index ? "0" : InputStr[index + 6]);

	for (j = 0; j < (1 << 12); j++)
		fprintf(F, "%x", MadeTable3[TableIndex[j]]);
	fprintf(F, ",\n");
	fclose(F);

	free(TableIndex);
}

int main()
{
	unsigned char**		InputTableIndex = NULL;
	unsigned short*		OutputTableIndex = NULL;
	unsigned char*		Unmasking = NULL;
	unsigned char*		UnmaskedInputTable = NULL;
	unsigned char**		ExtraCheckTable = NULL;

	unsigned char**		F3to1[3] = { NULL };
	unsigned char**		F3to1_full = NULL;
	int					NumberOfF3to1;
	ToupleStruct*		FTouple1 = NULL;
	ToupleStruct*		FTouple2 = NULL;
	ToupleStruct*		FTouple3 = NULL;
	unsigned int		NumberOfFTouples[3];
	unsigned int**		FTouple3InANF = NULL;
	unsigned short*		NumberOfANF3 = NULL;
	FILE*				F;

	unsigned char		x, y, z;
	unsigned char		xyz;
	unsigned char		TargetFunc[8];
	FunctionStruct		ANFTargetFunc;
	unsigned char		Invert;
	unsigned char		i, j;
	char				FileName[100];
	unsigned char		OutputBit;
	char			    OrigVars[4][5] = { "a", "b", "c", "d" };
	char				SelectedVars[3][5];
	char				SelectedVarIndexes[3];

	FillANFTables(512);
	FillANFTermTables();
	FillTableIndexes(InputTableIndex, OutputTableIndex, Unmasking, UnmaskedInputTable, ExtraCheckTable);

	for (OutputBit = 0; OutputBit <= 0; OutputBit++)
	{
		sprintf(FileName, "Res_%d.csv", OutputBit);
		F = fopen(FileName, "wt");
		fclose(F);

		//---------------------------------------

		SelectedVarIndexes[0] = 2;  //c + 1
		strcpy(SelectedVars[0], OrigVars[SelectedVarIndexes[0]]);

		FindCombination(FileName, OutputTableIndex,	SelectedVarIndexes, SelectedVars, 1);
	}

	printf("done");

	return 0;
}
