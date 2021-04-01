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

const char	path[500] = "Res_";

const unsigned char	NumberOfOutputBits = 4;
const unsigned char Bits[NumberOfOutputBits] = { 1,0,2,3 };

const unsigned char ReduceTables = 1;
const unsigned char	ReorderTabes = 0;
const unsigned char CheckCouples = 0;
const unsigned char	SameType = 0;
const unsigned char NumberOfCPUs = 22;

const long long		MaxCouplesPossible = 100000000L;

//-----------------------------------------------------

unsigned char*	Unmasking;

void FillTables()
{
	unsigned short	j, k;

	Unmasking = (unsigned char*)calloc(512, sizeof(unsigned char));

	for (j = 0; j < 512; j++)
		for (k = 0; k < 9; k++)
			Unmasking[j] ^= (j & (1 << k)) ? 1 : 0;
}

unsigned char CheckUniformity(unsigned char* TablesAddress[NumberOfOutputBits], unsigned char NumberOfTables,
	unsigned short** UniformityCounter)
{
	unsigned int	Masked_InputIndex;
	unsigned short	Masked_Output;
	unsigned short	i;
	unsigned short	j;
	unsigned short	k;
	unsigned short	ExpectedCounterValue[NumberOfOutputBits - 1];
	unsigned short	UniformityTable[NumberOfOutputBits - 1][1024];
	unsigned short	NumberOfUniformityTable;

	for (j = 0; j < NumberOfTables - 1; j++)
		ExpectedCounterValue[j] = 1024 / (1 << (2 * (j + 2)));

	NumberOfUniformityTable = 0;
	for (Masked_InputIndex = 0;Masked_InputIndex < (1 << 12);Masked_InputIndex++)
	{
		Masked_Output = 0;
		for (j = 0; j < NumberOfTables; j++)
		{
			Masked_Output <<= 3;
			Masked_Output |= TablesAddress[j][Masked_InputIndex];

			if (j)
			{
				UniformityTable[j - 1][NumberOfUniformityTable] = Masked_Output;

				if (++UniformityCounter[j - 1][Masked_Output] > ExpectedCounterValue[j - 1])
					break;
			}
		}

		if ((j < NumberOfTables) || ((Masked_InputIndex & 0xFF) == 0xFF)) // every 256
		{
			for (k = 0;k < NumberOfTables - 1;k++)
				for (i = 0;i < NumberOfUniformityTable;i++)
					UniformityCounter[k][UniformityTable[k][i]] = 0;

			NumberOfUniformityTable = 0;

			if (j < NumberOfTables)
				return(j);
		}
		else
			NumberOfUniformityTable++;
	}

	return(0); // is uniform
}

unsigned char CheckDistributions(unsigned char** FullTablesAddress[NumberOfOutputBits], unsigned char** SmallTablesAddress[NumberOfOutputBits], unsigned char NumberOfTables, unsigned short** Dis, unsigned short** DisOld)
{
	unsigned int	Masked_InputIndex;
	unsigned short	i, ii;
	unsigned short	j, jj;
	unsigned short	c, k;
	unsigned char   RejectIndex[630];

	c = 0;
	for (i = 1;i < NumberOfTables; i++)
	{
		for (ii = 0;ii < 3;ii++)
			for (j = 0; j < i; j++)
			{
				for (jj = 0;jj < 3;jj++)
					RejectIndex[c++] = i;

				for (jj = 0;jj < 9;jj++)
					RejectIndex[c++] = i;
			}

		for (ii = 0;ii < 9;ii++)
			for (j = 0; j < i; j++)
				for (jj = 0;jj < 3;jj++)
					RejectIndex[c++] = i;
	}

	//******************************//

	for (Masked_InputIndex = 0;Masked_InputIndex < (1 << 12);Masked_InputIndex++)
	{
		if (((Masked_InputIndex & 0xFF) == 0)) // every 256
			for (k = 0;k < c;k++)
				memset(Dis[k], 0, 64 * sizeof(unsigned short));

		c = 0;
		for (i = 1;i < NumberOfTables; i++)
		{
			for (ii = 0;ii < 3;ii++)
				for (j = 0; j < i; j++)
				{
					for (jj = 0;jj < 3;jj++)
						Dis[c++][(FullTablesAddress[j][jj][Masked_InputIndex] << 3) | FullTablesAddress[i][ii][Masked_InputIndex]]++;

					for (jj = 0;jj < 9;jj++)
						Dis[c++][(SmallTablesAddress[j][jj][Masked_InputIndex] << 3) | FullTablesAddress[i][ii][Masked_InputIndex]]++;
				}

			for (ii = 0;ii < 9;ii++)
				for (j = 0; j < i; j++)
					for (jj = 0;jj < 3;jj++)
						Dis[c++][(FullTablesAddress[j][jj][Masked_InputIndex] << 3) | SmallTablesAddress[i][ii][Masked_InputIndex]]++;
		}

		if (((Masked_InputIndex & 0xFF) == 0xFF)) // every 256
		{
			if (Masked_InputIndex == 0xFF) // the first one
			{
				for (k = 0;k < c;k++)
					memcpy(DisOld[k], Dis[k], 64 * sizeof(unsigned short));
			}
			else
			{
				for (k = 0;k < c;k++)
					if (memcmp(DisOld[k], Dis[k], 64 * sizeof(unsigned short)))
						return(RejectIndex[k]);
			}
		}
	}

	return(0); // distributions are identical
}



struct UniqueTablestruct
{
	unsigned char*	FullTable[3];
	unsigned char*	Table;
	unsigned char*	SmallTable[9];
	unsigned int	Type;
	unsigned int	StartIndex;
	unsigned int	Count;
};


int main()
{
	FILE*				F;
	char				FilePath[500];
	char*				TempStr;
	UniqueTablestruct  *UniqueTables[NumberOfOutputBits];
	unsigned int		NumberOfUniqueTables[NumberOfOutputBits];
	unsigned int		Type;
	uint64_t			ANF_old;
	uint64_t			ANF;
	char**				Strings[18][NumberOfOutputBits];
	unsigned int		NumberOfTables[NumberOfOutputBits];
	unsigned char*		Flags[NumberOfOutputBits];
	unsigned char		OnPattern[NumberOfOutputBits];
	unsigned char		OffPattern[NumberOfOutputBits];
	unsigned int		NumberOfReducedTables[NumberOfOutputBits];
	unsigned int*		ReducedTablesIndex[NumberOfOutputBits];

	unsigned int		TempSize[NumberOfOutputBits];
	unsigned int		min;
	unsigned int		max;
	unsigned char		TablesOrder[NumberOfOutputBits];
	unsigned char		TablesOrderInv[NumberOfOutputBits];
	char				TempChar;
	int					j;
	unsigned int		k, l;
	char				index;
	char				index2;
	int					i_last;
	unsigned int		i[NumberOfOutputBits + 1];
	int					i0;
	int					i1;
	unsigned char*		TablesAddress[NumberOfOutputBits];
	unsigned char**		FullTablesAddress[NumberOfOutputBits];
	unsigned char**		SmallTablesAddress[NumberOfOutputBits];
	long long			NumberOfFound;
	int					Counter;
	char				Last;
	char				CheckBit;
	char				Res;
	char**				CouplesPossible[NumberOfOutputBits][NumberOfOutputBits] = { NULL };
	int					ThreadNum;
	unsigned short***	UniformityCounter;
	unsigned short***	Dis;
	unsigned short***	DisOld;
	char                VarStr[12][4] = { "a2","b2","c2","d2", "a1","b1","c1","e1", "a3","b3","c3","d3" };
	char				Dependency[3];
	char				NumberOfDependencies;
	short				abcd;
	short				a1b1c1d1;
	short				a2b2c2d2;
	short				a3b3c3d3;
	short				Value[12];
	short				InputIndex;

	FillTables();
	TempStr = (char*)malloc(100000 * sizeof(char));

	for (index = 0; index < NumberOfOutputBits; index++)
	{
		printf("getting size of table %d", index);

		sprintf(FilePath, "%s%d.csv", path, Bits[index]);
		F = fopen(FilePath, "rt");
		NumberOfTables[index] = 0;
		Last = '\n';

		while (!feof(F))
		{
			j = fread(TempStr, 1, 100000, F);

			for (k = 0; k < j; k++)
			{
				TempChar = TempStr[k];
				if (((TempChar == '\n') | (TempChar == '\r')) & (Last != '\n') & (Last != '\r'))
					NumberOfTables[index]++;
				Last = TempChar;
			}
		}

		printf(": %d, reading", NumberOfTables[index]);
		fseek(F, 0, SEEK_SET);

		UniqueTables[index] = (UniqueTablestruct*)malloc(NumberOfTables[index] * sizeof(UniqueTablestruct));
		for (k = 0; k < 18; k++)
			Strings[k][index] = (char**)malloc(NumberOfTables[index] * sizeof(char*));

		ANF_old = 0;
		NumberOfUniqueTables[index] = 0;
		for (j = 0; j < NumberOfTables[index]; j++)
		{
			fscanf(F, "%s", TempStr);
			TempStr[strlen(TempStr) - 1] = 0;
			//Type = atoi(TempStr);
			Type = 0;

			fscanf(F, "%s", TempStr);
			TempStr[strlen(TempStr) - 1] = 0;
			sscanf(TempStr, "%" SCNu64 , &ANF);

			for (k = 0; k < 6; k++)
			{
				l = 0;
				do
				{
					TempChar = fgetc(F);
					if (TempChar != '\n')
						TempStr[l++] = TempChar;
				} while (TempChar != ',');

				TempStr[l] = 0;
				Strings[k][index][j] = (char*)malloc((l + 1) * sizeof(char));
				strcpy(Strings[k][index][j], TempStr);
			}

			fscanf(F, "%s", TempStr);
			if (1) //(ANF != ANF_old)
			{
				UniqueTables[index][NumberOfUniqueTables[index]].StartIndex = j;
				UniqueTables[index][NumberOfUniqueTables[index]].Count = 1;
				UniqueTables[index][NumberOfUniqueTables[index]].Type = Type;
				UniqueTables[index][NumberOfUniqueTables[index]].Table = (unsigned char*)malloc((1 << 12) * sizeof(unsigned char));
				UniqueTables[index][NumberOfUniqueTables[index]].FullTable[0] = (unsigned char*)malloc((1 << 12) * sizeof(unsigned char));

				for (k = 0; k < (1 << 12); k++)
				{
					UniqueTables[index][NumberOfUniqueTables[index]].FullTable[0][k] = (TempStr[k] <= '9') ? (TempStr[k] - '0') : (TempStr[k] - 'a' + 10);
					UniqueTables[index][NumberOfUniqueTables[index]].Table[k] = Unmasking[UniqueTables[index][NumberOfUniqueTables[index]].FullTable[0][k]];
				}

				//****************************//

				for (k = 0;k < 3;k++)
				{
					UniqueTables[index][NumberOfUniqueTables[index]].SmallTable[k] = (unsigned char*)malloc((1 << 12) * sizeof(unsigned char));

					NumberOfDependencies = 0;
					for (l = 0;l < 12;l++)
						if (strstr(Strings[k * 2 + 1][index][j], VarStr[l]) != NULL)
							Dependency[NumberOfDependencies++] = l;

					InputIndex = 0;
					for (abcd = 0; abcd < 16; abcd++)
						for (a1b1c1d1 = 0; a1b1c1d1 < 16; a1b1c1d1++)
							for (a2b2c2d2 = 0; a2b2c2d2 < 16; a2b2c2d2++)
							{
								a3b3c3d3 = abcd ^ a1b1c1d1 ^ a2b2c2d2;
								for (l = 0; l < 4; l++)
								{
									Value[l + 0] = (a2b2c2d2 >> l) & 1;
									Value[l + 4] = (a1b1c1d1 >> l) & 1;
									Value[l + 8] = (a3b3c3d3 >> l) & 1;
								}

								UniqueTables[index][NumberOfUniqueTables[index]].SmallTable[k][InputIndex] = 0;
								for (l = 0;l < NumberOfDependencies;l++)
									UniqueTables[index][NumberOfUniqueTables[index]].SmallTable[k][InputIndex] |= Value[Dependency[l]] << l;

								InputIndex++;
							}
				}
			}
			else
				UniqueTables[index][NumberOfUniqueTables[index]].Count++;

			fscanf(F, "%s", TempStr); // skip the second ANF

			for (k = 6; k < 12; k++)
			{
				l = 0;
				do
				{
					TempChar = fgetc(F);
					if (TempChar != '\n')
						TempStr[l++] = TempChar;
				} while (TempChar != ',');

				TempStr[l] = 0;
				Strings[k][index][j] = (char*)malloc((l + 1) * sizeof(char));
				strcpy(Strings[k][index][j], TempStr);
			}

			fscanf(F, "%s", TempStr);
			if (1) //ANF != ANF_old)
			{
				UniqueTables[index][NumberOfUniqueTables[index]].FullTable[1] = (unsigned char*)malloc((1 << 12) * sizeof(unsigned char));

				for (k = 0; k < (1 << 12); k++)
				{
					UniqueTables[index][NumberOfUniqueTables[index]].FullTable[1][k] = (TempStr[k] <= '9') ? (TempStr[k] - '0') : (TempStr[k] - 'a' + 10);
					UniqueTables[index][NumberOfUniqueTables[index]].Table[k] |= Unmasking[UniqueTables[index][NumberOfUniqueTables[index]].FullTable[1][k]] << 1;
				}

				//****************************//

				for (k = 0;k < 3;k++)
				{
					UniqueTables[index][NumberOfUniqueTables[index]].SmallTable[3 + k] = (unsigned char*)malloc((1 << 12) * sizeof(unsigned char));

					NumberOfDependencies = 0;
					for (l = 0;l < 12;l++)
						if (strstr(Strings[k * 2 + 7][index][j], VarStr[l]) != NULL)
							Dependency[NumberOfDependencies++] = l;

					InputIndex = 0;
					for (abcd = 0; abcd < 16; abcd++)
						for (a1b1c1d1 = 0; a1b1c1d1 < 16; a1b1c1d1++)
							for (a2b2c2d2 = 0; a2b2c2d2 < 16; a2b2c2d2++)
							{
								a3b3c3d3 = abcd ^ a1b1c1d1 ^ a2b2c2d2;
								for (l = 0;l < 4;l++)
								{
									Value[l + 0] = (a2b2c2d2 >> l) & 1;
									Value[l + 4] = (a1b1c1d1 >> l) & 1;
									Value[l + 8] = (a3b3c3d3 >> l) & 1;
								}

								UniqueTables[index][NumberOfUniqueTables[index]].SmallTable[3 + k][InputIndex] = 0;
								for (l = 0;l < NumberOfDependencies;l++)
									UniqueTables[index][NumberOfUniqueTables[index]].SmallTable[3 + k][InputIndex] |= Value[Dependency[l]] << l;

								InputIndex++;
							}
				}
			}

			fscanf(F, "%s", TempStr); // skip the third ANF

			for (k = 12; k < 18; k++)
			{
				l = 0;
				do
				{
					TempChar = fgetc(F);
					if (TempChar != '\n')
						TempStr[l++] = TempChar;
				} while (TempChar != ',');

				TempStr[l] = 0;
				Strings[k][index][j] = (char*)malloc((l + 1) * sizeof(char));
				strcpy(Strings[k][index][j], TempStr);
			}

			fscanf(F, "%s", TempStr);
			if (1) //ANF != ANF_old)
			{
				UniqueTables[index][NumberOfUniqueTables[index]].FullTable[2] = (unsigned char*)malloc((1 << 12) * sizeof(unsigned char));

				for (k = 0; k < (1 << 12); k++)
				{
					UniqueTables[index][NumberOfUniqueTables[index]].FullTable[2][k] = (TempStr[k] <= '9') ? (TempStr[k] - '0') : (TempStr[k] - 'a' + 10);
					UniqueTables[index][NumberOfUniqueTables[index]].Table[k] |= Unmasking[UniqueTables[index][NumberOfUniqueTables[index]].FullTable[2][k]] << 2;
				}

				//****************************//

				for (k = 0;k < 3;k++)
				{
					UniqueTables[index][NumberOfUniqueTables[index]].SmallTable[6 + k] = (unsigned char*)malloc((1 << 12) * sizeof(unsigned char));

					NumberOfDependencies = 0;
					for (l = 0;l < 12;l++)
						if (strstr(Strings[k * 2 + 13][index][j], VarStr[l]) != NULL)
							Dependency[NumberOfDependencies++] = l;

					InputIndex = 0;
					for (abcd = 0; abcd < 16; abcd++)
						for (a1b1c1d1 = 0; a1b1c1d1 < 16; a1b1c1d1++)
							for (a2b2c2d2 = 0; a2b2c2d2 < 16; a2b2c2d2++)
							{
								a3b3c3d3 = abcd ^ a1b1c1d1 ^ a2b2c2d2;
								for (l = 0; l < 4; l++)
								{
									Value[l + 0] = (a2b2c2d2 >> l) & 1;
									Value[l + 4] = (a1b1c1d1 >> l) & 1;
									Value[l + 8] = (a3b3c3d3 >> l) & 1;
								}

								UniqueTables[index][NumberOfUniqueTables[index]].SmallTable[6 + k][InputIndex] = 0;
								for (l = 0;l < NumberOfDependencies;l++)
									UniqueTables[index][NumberOfUniqueTables[index]].SmallTable[6 + k][InputIndex] |= Value[Dependency[l]] << l;

								InputIndex++;
							}
				}

				//****************************//

				NumberOfUniqueTables[index]++;
			}

			ANF_old = ANF;
		}

		fclose(F);

		printf(", unique: %d\n", NumberOfUniqueTables[index]);
	}

	//====================================================

	omp_set_num_threads(NumberOfCPUs);
	UniformityCounter = (unsigned short***)malloc(omp_get_max_threads() * sizeof(unsigned short**));
	Dis = (unsigned short***)malloc(omp_get_max_threads() * sizeof(unsigned short**));
	DisOld = (unsigned short***)malloc(omp_get_max_threads() * sizeof(unsigned short**));

	for (j = 0;j < omp_get_max_threads();j++)
	{
		UniformityCounter[j] = (unsigned short**)malloc((NumberOfOutputBits - 1) * sizeof(unsigned short*));
		Dis[j] = (unsigned short**)malloc(630 * sizeof(unsigned short*));
		DisOld[j] = (unsigned short**)malloc(630 * sizeof(unsigned short*));

		for (index = 0; index < NumberOfOutputBits - 1; index++)
			UniformityCounter[j][index] = (unsigned short*)calloc((1 << (3 * (index + 2))), sizeof(unsigned short));

		for (i0 = 0; i0 < 630; i0++)
		{
			Dis[j][i0] = (unsigned short*)malloc(64 * sizeof(unsigned short));
			DisOld[j][i0] = (unsigned short*)malloc(64 * sizeof(unsigned short));
		}
	}

	//====================================================

	memcpy(TempSize, NumberOfUniqueTables, NumberOfOutputBits * sizeof(unsigned int));
	for (index = 0; index < NumberOfOutputBits; index++)
	{
		if (ReorderTabes)
		{
			min = 0;
			for (j = 1; j < NumberOfOutputBits; j++)
				if (TempSize[min] > TempSize[j])
					min = j;

			TablesOrder[index] = min;
			TempSize[min] = -1;
		}
		else
			TablesOrder[index] = index;
	}

	for (index = 0; index < NumberOfOutputBits; index++)
		TablesOrderInv[TablesOrder[index]] = index;

	//====================================================

	for (index = 0; index < NumberOfOutputBits; index++)
	{
		OnPattern[index] = 0;
		for (index2 = 0; index2 <= index; index2++)
			OnPattern[index] |= (1 << index2);

		OffPattern[index] = ~(1 << index);
	}

	for (index = 0; index < NumberOfOutputBits; index++)
		Flags[index] = (unsigned char*)calloc(NumberOfUniqueTables[index], sizeof(unsigned char));

	printf("\n");

	for (index = 0; index < NumberOfOutputBits; index++)
		NumberOfReducedTables[TablesOrder[index]] = NumberOfUniqueTables[TablesOrder[index]];

	for (index = 0; index < NumberOfOutputBits; index++)
	{
		if (ReduceTables)
		{
			printf("Checking table %d: ", TablesOrder[index]);

			NumberOfReducedTables[TablesOrder[index]] = 0;
			for (i[0] = 0; i[0] < NumberOfUniqueTables[TablesOrder[index]]; i[0]++)
				if (Flags[TablesOrder[index]][i[0]] == (OnPattern[index] & OffPattern[index]))
					NumberOfReducedTables[TablesOrder[index]]++;
			printf("%d (", NumberOfReducedTables[TablesOrder[index]]);

			for (index2 = index + 1; index2 < NumberOfOutputBits; index2++)
			{
				if (index2 != index + 1)
					printf(", ");
				printf("with %d ", TablesOrder[index2]);

				if (((long long)NumberOfReducedTables[TablesOrder[index]]) * ((long long)NumberOfReducedTables[TablesOrder[index2]]) < MaxCouplesPossible)
				{
					#pragma omp parallel for schedule(guided) private(ThreadNum, i1, TablesAddress, FullTablesAddress, SmallTablesAddress, CheckBit)
					for (i0 = 0; i0 < NumberOfUniqueTables[TablesOrder[index]]; i0++)
					{
						ThreadNum = omp_get_thread_num();

						if (Flags[TablesOrder[index]][i0] == (OnPattern[index2 - 1] & OffPattern[index]))
						{
							TablesAddress[0] = UniqueTables[TablesOrder[index]][i0].Table;
							FullTablesAddress[0]= UniqueTables[TablesOrder[index]][i0].FullTable;
							SmallTablesAddress[0]= UniqueTables[TablesOrder[index]][i0].SmallTable;
							CheckBit = 0;

							for (i1 = 0; i1 < NumberOfUniqueTables[TablesOrder[index2]]; i1++)
								if (((!SameType) ||
									(UniqueTables[TablesOrder[index]][i0].Type == UniqueTables[TablesOrder[index2]][i1].Type)) &&
										((((index == 0) & (Flags[TablesOrder[index2]][i1] == 0)) ||
									((index != 0) & (Flags[TablesOrder[index2]][i1] == (OnPattern[index - 1] & OffPattern[index2]))) ||
											((Flags[TablesOrder[index2]][i1] == (OnPattern[index] & OffPattern[index2])) & (!CheckBit)))))
								{
									TablesAddress[1] = UniqueTables[TablesOrder[index2]][i1].Table;
									FullTablesAddress[1] = UniqueTables[TablesOrder[index2]][i1].FullTable;
									SmallTablesAddress[1] = UniqueTables[TablesOrder[index2]][i1].SmallTable;

									if (!CheckDistributions(FullTablesAddress, SmallTablesAddress, 2, Dis[ThreadNum], DisOld[ThreadNum]))
										if (!CheckUniformity(TablesAddress, 2, UniformityCounter[ThreadNum]))
										{
											CheckBit = 1;
											#pragma omp atomic
											Flags[TablesOrder[index2]][i1] |= (1 << index);
										}
								}

							if (CheckBit)
							{
								#pragma omp atomic
								Flags[TablesOrder[index]][i0] |= (1 << index2);
							}
						}
					}
				}
				else
				{
					for (i0 = 0; i0 < NumberOfUniqueTables[TablesOrder[index]]; i0++)
						Flags[TablesOrder[index]][i0] |= (1 << index2);

					for (i1 = 0; i1 < NumberOfUniqueTables[TablesOrder[index2]]; i1++)
						Flags[TablesOrder[index2]][i1] |= (1 << index);
				}

				if (index2 != NumberOfOutputBits - 1)
				{
					NumberOfReducedTables[TablesOrder[index]] = 0;
					for (i[0] = 0; i[0] < NumberOfUniqueTables[TablesOrder[index]]; i[0]++)
						if (Flags[TablesOrder[index]][i[0]] == (OnPattern[index2] & OffPattern[index]))
							NumberOfReducedTables[TablesOrder[index]]++;
					printf("-> %d", NumberOfReducedTables[TablesOrder[index]]);
				}
			}

			NumberOfReducedTables[TablesOrder[index]] = 0;
			for (i[0] = 0; i[0] < NumberOfUniqueTables[TablesOrder[index]]; i[0]++)
				if (Flags[TablesOrder[index]][i[0]] == (OnPattern[NumberOfOutputBits - 1] & OffPattern[index]))
					NumberOfReducedTables[TablesOrder[index]]++;

			printf(") -> %d\n", NumberOfReducedTables[TablesOrder[index]]);
		}
		else
			NumberOfReducedTables[TablesOrder[index]] = NumberOfUniqueTables[TablesOrder[index]];

		ReducedTablesIndex[TablesOrder[index]] = (unsigned int*)malloc(NumberOfReducedTables[TablesOrder[index]] * sizeof(unsigned int));
		NumberOfReducedTables[TablesOrder[index]] = 0;
		for (i[0] = 0; i[0] < NumberOfUniqueTables[TablesOrder[index]]; i[0]++)
			if ((ReduceTables == 0) ||
				(Flags[TablesOrder[index]][i[0]] == (OnPattern[NumberOfOutputBits - 1] & OffPattern[index])))
				ReducedTablesIndex[TablesOrder[index]][NumberOfReducedTables[TablesOrder[index]]++] = i[0];
	}

	printf("\n");

	//====================================================

	memcpy(TempSize, NumberOfReducedTables, NumberOfOutputBits * sizeof(unsigned int));

	for (index = 0; index < NumberOfOutputBits; index++)
	{
		if (ReorderTabes)
		{
			max = 0;
			for (j = 1; j < NumberOfOutputBits; j++)
				if (TempSize[max] < TempSize[j])
					max = j;

			TablesOrder[index] = max;
			TempSize[max] = 0;
		}
		else
			TablesOrder[index] = index;
	}

	for (index = 0; index < NumberOfOutputBits; index++)
		TablesOrderInv[TablesOrder[index]] = index;

	//====================================================

	for (index2= NumberOfOutputBits-1;index2>=1;index2--)
		for (index = index2 - 1; index >=0; index--)
			if (CheckCouples)
			{
				printf("Checking possible couples (%d,%d) ", TablesOrder[index2], TablesOrder[index]);

				if (((long long)NumberOfReducedTables[TablesOrder[index2]]) * ((long long)NumberOfReducedTables[TablesOrder[index]]) < MaxCouplesPossible)
				{
					if (NumberOfReducedTables[TablesOrder[index2]])
						CouplesPossible[TablesOrder[index2]][TablesOrder[index]] = (char**)malloc(NumberOfReducedTables[TablesOrder[index2]] * sizeof(char*));
					else
						CouplesPossible[TablesOrder[index2]][TablesOrder[index]] = NULL;

					#pragma omp parallel for schedule(guided) private(ThreadNum, i0, TablesAddress, FullTablesAddress, SmallTablesAddress)
					for (i1 = 0; i1 < NumberOfReducedTables[TablesOrder[index2]]; i1++)
					{
						ThreadNum = omp_get_thread_num();

						TablesAddress[0] = UniqueTables[TablesOrder[index2]][ReducedTablesIndex[TablesOrder[index2]][i1]].Table;
						FullTablesAddress[0] = UniqueTables[TablesOrder[index2]][ReducedTablesIndex[TablesOrder[index2]][i1]].FullTable;
						SmallTablesAddress[0] = UniqueTables[TablesOrder[index2]][ReducedTablesIndex[TablesOrder[index2]][i1]].SmallTable;
						CouplesPossible[TablesOrder[index2]][TablesOrder[index]][i1] = (char*)malloc(NumberOfReducedTables[TablesOrder[index]] * sizeof(char));

						for (i0 = 0; i0 < NumberOfReducedTables[TablesOrder[index]]; i0++)
						{
							TablesAddress[1] = UniqueTables[TablesOrder[index]][ReducedTablesIndex[TablesOrder[index]][i0]].Table;
							FullTablesAddress[1] = UniqueTables[TablesOrder[index]][ReducedTablesIndex[TablesOrder[index]][i0]].FullTable;
							SmallTablesAddress[1] = UniqueTables[TablesOrder[index]][ReducedTablesIndex[TablesOrder[index]][i0]].SmallTable;

							if (!CheckDistributions(FullTablesAddress, SmallTablesAddress, 2, Dis[ThreadNum], DisOld[ThreadNum]))
								CouplesPossible[TablesOrder[index2]][TablesOrder[index]][i1][i0] = CheckUniformity(TablesAddress, 2, UniformityCounter[ThreadNum]);
							else
								CouplesPossible[TablesOrder[index2]][TablesOrder[index]][i1][i0] = 1;
						}
					}

					printf("done\n");
				}
				else
				{
					CouplesPossible[TablesOrder[index2]][TablesOrder[index]] = NULL;
					printf("too large\n");
				}
			}
			else
				CouplesPossible[TablesOrder[index2]][TablesOrder[index]] = NULL;

	//====================================================

	strcpy(FilePath, path);
	for (index = 0; index < NumberOfOutputBits; index++)
		sprintf(FilePath, "%s%d", FilePath, Bits[index]);

	strcat(FilePath, ".csv");
	F = fopen(FilePath, "wt");
	fclose(F);

	NumberOfFound = 0;
	Counter = 0;
	#pragma omp parallel for schedule(guided) private(ThreadNum, j, i, index, index2, TablesAddress, FullTablesAddress, SmallTablesAddress, Res)
	for (i_last = 0; i_last < NumberOfReducedTables[TablesOrder[NumberOfOutputBits - 1]]; i_last++)
	{
		ThreadNum = omp_get_thread_num();

		#pragma omp atomic
		Counter++;

		#pragma omp critical (print)
		{
			printf("%d / %d\n", Counter, NumberOfReducedTables[TablesOrder[NumberOfOutputBits - 1]]);
		}

		for (j = 0; j < NumberOfOutputBits - 1; j++)
			i[j] = 0;
		i[NumberOfOutputBits - 1] = i_last;

		while (i[NumberOfOutputBits - 1] == i_last)
		{
			for (index = NumberOfOutputBits - 1; index >= 0; index--)
			{
				if (SameType && (index < (NumberOfOutputBits - 1)) &&
					(UniqueTables[TablesOrder[index]][ReducedTablesIndex[TablesOrder[index]][i[index]]].Type !=
						UniqueTables[TablesOrder[index + 1]][ReducedTablesIndex[TablesOrder[index + 1]][i[index + 1]]].Type))
					break;

				for (index2 = NumberOfOutputBits - 1; index2 > index; index2--)
					if (CouplesPossible[TablesOrder[index2]][TablesOrder[index]] &&
						CouplesPossible[TablesOrder[index2]][TablesOrder[index]][i[index2]][i[index]])
						break;

				if (index2 > index)
					break;

				FullTablesAddress[NumberOfOutputBits - 1 - index] = UniqueTables[TablesOrder[index]][ReducedTablesIndex[TablesOrder[index]][i[index]]].FullTable;
				SmallTablesAddress[NumberOfOutputBits - 1 - index] = UniqueTables[TablesOrder[index]][ReducedTablesIndex[TablesOrder[index]][i[index]]].SmallTable;
				TablesAddress[NumberOfOutputBits - 1 - index] = UniqueTables[TablesOrder[index]][ReducedTablesIndex[TablesOrder[index]][i[index]]].Table;
			}

			if (index >= 0)
			{
				for (j = index - 1; j >= 0; j--)
					i[j] = NumberOfReducedTables[TablesOrder[j]];
			}
			else
			{
				Res = CheckDistributions(FullTablesAddress, SmallTablesAddress, NumberOfOutputBits, Dis[ThreadNum], DisOld[ThreadNum]);
				if (!Res)
					Res = CheckUniformity(TablesAddress, NumberOfOutputBits, UniformityCounter[ThreadNum]);

				if (Res)
				{
					for (j = NumberOfOutputBits - 2 - Res; j >= 0; j--)
						i[j] = NumberOfReducedTables[TablesOrder[j]];
				}
				else
				{
					#pragma omp critical (found)
					{
						unsigned int small_i[NumberOfOutputBits + 1];

						for (index = 0; index < NumberOfOutputBits + 1; index++)
							small_i[index] = 0;

						F = fopen(FilePath, "at");

						NumberOfFound++;
						if ((NumberOfFound & 0xff) == 0xff)
						{
							printf("%d / %d found %I64d ", Counter, NumberOfReducedTables[TablesOrder[NumberOfOutputBits - 1]], NumberOfFound);
							for (index = 0; index < NumberOfOutputBits; index++)
								printf("%d ", i[index]);
							printf("\n");
						}

						while (!small_i[NumberOfOutputBits])
						{
							NumberOfFound++;
							if ((NumberOfFound & 0xff) == 0xff)
							{
								printf("%d / %d found %I64d ", Counter, NumberOfReducedTables[TablesOrder[NumberOfOutputBits - 1]], NumberOfFound);
								for (index = 0; index < NumberOfOutputBits; index++)
									printf("%d ", i[index]);
								printf("\n");
							}

							for (index = 0; index < NumberOfOutputBits; index++)
							{
								//fprintf(F, "%08d, ", UniqueTables[index][ReducedTablesIndex[index][i[TablesOrderInv[index]]]].Type);
								for (k = 0; k < 18; k++)
									fprintf(F, "%s ", Strings[k][index][UniqueTables[index][ReducedTablesIndex[index][i[TablesOrderInv[index]]]].StartIndex + small_i[index]]);
								fprintf(F, ", ");
							}

							fprintf(F, "\n");


							j = 0;
							do
							{
								small_i[j]++;
								if (j < NumberOfOutputBits)
									if (small_i[j] >= UniqueTables[j][ReducedTablesIndex[j][i[TablesOrderInv[j]]]].Count)
										small_i[j] = 0;
							} while (small_i[j++] == 0);
						}
						fclose(F);
					}
				}
			}

			j = 0;
			do
			{
				i[j]++;
				if (j < NumberOfOutputBits)
					if (i[j] >= NumberOfReducedTables[TablesOrder[j]])
						i[j] = 0;
			} while (i[j++] == 0);
		}
	}

	printf("Total found %I64d\n\n", NumberOfFound);

	//--------------------------------------------------

	for (index = 0; index < NumberOfOutputBits; index++)
	{
		for (index2 = 0; index2 < NumberOfOutputBits; index2++)
			if (CouplesPossible[TablesOrder[index]][TablesOrder[index2]])
			{
				for (i0 = 0; i0 < NumberOfReducedTables[TablesOrder[index]]; i0++)
					free(CouplesPossible[TablesOrder[index]][TablesOrder[index2]][i0]);

				free(CouplesPossible[TablesOrder[index]][TablesOrder[index2]]);
			}

		if (ReducedTablesIndex[index])
			free(ReducedTablesIndex[index]);
		if (Flags[index])
			free(Flags[index]);

		for (j = 0; j < NumberOfUniqueTables[index]; j++)
			free(UniqueTables[index][j].Table);
		free(UniqueTables[index]);

		for (j = 0; j < NumberOfTables[index]; j++)
		{
			for (k = 0; k < 18; k++)
				free(Strings[k][index][j]);
		}

		for (k = 0; k < 18; k++)
			free(Strings[k][index]);
	}

	printf("done");

	return 0;
}
