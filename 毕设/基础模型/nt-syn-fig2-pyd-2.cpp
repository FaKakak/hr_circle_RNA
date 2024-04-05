// nt-syn-fig2-pyd-2.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
// nt-syn-fig2-pyd-1.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
// nt-syn-fig2a-pyd.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
//"Nucleotide synthetase ribozymes may have emerged first in the RNA World"
//by Wentao Ma et al.
//-----C source codes for the simulation program
//----Parameters are for the case in figure 2a of the article

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <conio.h>
#include <ctime>

/******* RANDOM NUMBER GENERATOR BY ZIFF **********/
#define A1 471
#define B1 1586
#define C1 6988
#define D1 9689
#define M 16383
#define RIMAX 2147483648.0        // = 2^31 
#define RandomInteger (++nd, ra[nd & M] = ra[(nd-A1) & M] ^ ra[(nd-B1) & M] ^ ra[(nd-C1) & M] ^ ra[(nd-D1) & M])
void seed(long seed);  // random number initialization 
static long ra[M + 1], nd;
/**************************************************/

#define LEN sizeof(struct rna)
#define C 2
#define G 3
#define A 1
#define U 4
#define STEPNUM 3000000    // Total time steps of Monte Carlo simulation
#define STAREC 0          // The step to start record
#define RECINT 10000       // The interval steps of recording
#define MAX_RNA_LENGTH 100    // Defining maximum RNA length allowed in the simulation
#define LONG_CHAIN_LEN 30     // Defining long chains for recording

#define SD 5   //555
#define N 30                // The N length of the two-dimensional grid
#define TOTAL_MATERIAL 40000  // Total materials in the system
#define NRSEQ A,G,C,A,U,G,C,U   // The presumed specific sequence with which a polynucleotide could act as a nt-synthetase
#define CTR1SEQ C,U,C,U,A,G,A,G  // 
#define CTR2SEQ C,G,U,U,A,A,C,G  // 
#define CTR3SEQ A,U,C,G,C,G,A,U  // 

#define INOCUSEQ NRSEQ
#define INOCUNUM 10
#define INOCUSEQ1 CTR1SEQ
#define INOCUSEQ2 CTR2SEQ   
#define INOCUSEQ3 CTR3SEQ   

#define INOCUSTEP 10000
#define PSBP 0.9             // Probability of separation of a base-pair
#define PBB 0.000001  //0.000005         // Probability of breaking of a phosphodiester bond
#define PLMC 0.000002        // Probability of ligation of two unit in a cell with mineral catalysis 
#define PAT 0.1              // Probability of attracting a substrate by a template when the substrate could base-pair with the template
#define PLT 0.9              // Probability of a template-directed ligation 
#define PMOVR 0.01           // Probability of movement of raw
#define PMOV (PMOVR/2)       // Probability of movement of a mononucleotide
#define PMF 0.0002           // probability of mononucleotide formation from raw materials
#define PMFS 0.9             // probability of mononucleotide formation from raw materials under the catalysis of nt-synthetase
#define TNSS 1               // Turn of nt-synthesis by nt-synthetase each step
#define PMD 0.01             // probability of mononucleotide decay to raw materials

#define PNDE  0.00002  // Ma --- Probability of nucleotide decaying into its precursor at RNA's chain end


#define PFP 0.001   //0.01             // probability of false base-pairing
#define FDMOV (pow(p->length1+p->length2,1/3.0))  
							 // The factor defining the relationship between probability of moving and molecular weight
#define CELLNUM (N*N)  // Cell numbers in the grid

long randl(long);         // random number between 0 and parameter 
double randd(void);       // random double between 0 and 1         
void avail_xy_init(void); // Initialization for xy_choose
void xy_choose(void);     // Picks a cell at random 
void fresh_unit(void);    // Updating a unit for the next time step
int findseq(char seq[], int seqlength, struct rna* p); //find a specific subsequence in a sequence 
void inits(void);         // initialization of the system
void inoculate(void);
void unit_case(void);     // Action of units (molecules) in the system
int record(void);         // Data recording
void freepool(void);      // Memory releasing

struct rna                // A mononucleotide or polynucleotide
{
	char information[2][MAX_RNA_LENGTH];
	int length1;
	int length2;
	int nick;
	struct rna* next;
	struct rna* prior;
};
struct rna* cell_head[2][N][N];
struct rna* p, * p1, * p2, * p3;

static char nrseq[50] = { NRSEQ };        // Presumed nt_synthetase sequence
static char ctr1seq[50] = { CTR1SEQ };
static char ctr2seq[50] = { CTR2SEQ };
static char ctr3seq[50] = { CTR3SEQ };
static char inocuseq[50] = { INOCUSEQ };
static char inocuseq1[50] = { INOCUSEQ1 };
static char inocuseq2[50] = { INOCUSEQ2 };
static char inocuseq3[50] = { INOCUSEQ3 };
static int raw_arr[N][N];

int over_max_len = 0;
int x, y;                 // The coordinate of cells in the grid 
int nrlength, ctr1length, ctr2length, ctr3length;
int inoculength, inoculength1, inoculength2, inoculength3;
int randcase, randcase1, randcaser, randcaser1, g = 0, h = 0;
int flag, flag1, flag2, flag3;
long i;                  // Cycle variable for Monte Carlo steps
long available;
long availabl[CELLNUM];
long recstep[(STEPNUM - STAREC) / RECINT + 1];    // Record steps
float nr_num[(STEPNUM - STAREC) / RECINT + 1];   // Record number of nt-synthetases in steps
float ctr1_num[(STEPNUM - STAREC) / RECINT + 1];  // 
float ctr2_num[(STEPNUM - STAREC) / RECINT + 1];  // 
float ctr3_num[(STEPNUM - STAREC) / RECINT + 1];  // 
float total_mat_num[(STEPNUM - STAREC) / RECINT + 1];  // Record number of total materials in steps
float RNA_num[(STEPNUM - STAREC) / RECINT + 1];  // Record number of units in steps
float raw_num[(STEPNUM - STAREC) / RECINT + 1];  // Record number of raw in steps
time_t time_start, time_end, time1, time2, test;
time_t* ptimer = &test;
// Number of replicases and polynucleotides including the complementary sequence of presumed replicase sequence 


/***********************************************************
 * Random generator initialization                         *
 *                      by a simple Congruential generator *
 ***********************************************************/
void seed(long seed)
{
	int a;

	if (seed < 0) { puts("SEED error."); exit(1); }
	ra[0] = (long)fmod(16807.0 * (double)seed, 2147483647.0);
	for (a = 1; a <= M; a++)
	{
		ra[a] = (long)fmod(16807.0 * (double)ra[a - 1], 2147483647.0);
	}
}

//------------------------------------------------------------------------------
long randl(long num)      /* random integer number between 0 and num-1 */
{
	return(RandomInteger % num);
}

//------------------------------------------------------------------------------
double randd(void)        /* random real number between 0 and 1 */
{
	return((double)RandomInteger / RIMAX);
}

//------------------------------------------------------------------------------
void avail_xy_init(void)   // Initialization for xy_choose
{
	int j;
	for (j = 0; j < CELLNUM; j++)
	{
		availabl[j] = j + 1;
	}
	available = CELLNUM;
}

//------------------------------------------------------------------------------
void xy_choose(void)       // Picks a cell at random
{
	long rl, s;
	rl = randl(available);
	s = availabl[rl];
	x = (s - 1) % N;
	y = (s - 1) / N;
	availabl[rl] = availabl[available - 1];
	available--;
}

//------------------------------------------------------------------------------
void fresh_unit(void)     // Updating a unit for the next time step
{
	p1 = p->prior;
	p2 = p->next;
	p3 = cell_head[!h][y][x]->next;
	cell_head[!h][y][x]->next = p;
	p->next = p3;
	p->prior = cell_head[!h][y][x];
	if (p3 != cell_head[!h][y][x])p3->prior = p;
	p1->next = p2;
	if (p2 != cell_head[h][y][x])p2->prior = p1;
	p = p1;
}

int findseq(char subseq[], int subseqlength, struct rna* p)  // Find a specific subsequence in an RNA
{
	int a, b, flag1;

	flag1 = 1;   // Assuming the RNA does not contain the subsequence
	if (p->length1 >= subseqlength)
	{
		for (b = 0; p->length1 - subseqlength - b >= 0; b++)
		{
			for (a = 0; a < subseqlength; a++)
			{
				if (p->information[0][b + a] != subseq[a]) break;
			}
			if (a == subseqlength) { flag1 = 0; break; }  // The subsequence has been found in the RNA
		}
	}

	if (flag1 == 0)return(0);   // Yes, the RNA contains the subsequence
	else return(1);   // No, the RNA does not contain the subsequence
}

//------------------------------------------------------------------------------
void inits(void)         // initialization of the system
{
	int j, m, k;
	seed(SD);

	nrlength = 0;
	for (j = 0; nrseq[j] != 0; j++)
		nrlength++;

	ctr1length = 0;
	for (j = 0; ctr1seq[j] != 0; j++)
		ctr1length++;

	ctr2length = 0;
	for (j = 0; ctr2seq[j] != 0; j++)
		ctr2length++;

	ctr3length = 0;
	for (j = 0; ctr3seq[j] != 0; j++)
		ctr3length++;
	///////////////////////////////

	inoculength = 0;
	for (j = 0; inocuseq[j] != 0; j++)
		inoculength++;

	inoculength1 = 0;
	for (j = 0; inocuseq1[j] != 0; j++)
		inoculength1++;

	inoculength2 = 0;
	for (j = 0; inocuseq2[j] != 0; j++)
		inoculength2++;

	inoculength3 = 0;
	for (j = 0; inocuseq3[j] != 0; j++)
		inoculength3++;

	for (m = 0; m < 2; m++)
	{
		for (y = 0; y < N; y++)
		{
			for (x = 0; x < N; x++)
			{
				p1 = (struct rna*)malloc(LEN);
				if (!p1) { printf("\tinit1--memeout\n"); exit(0); }
				cell_head[m][y][x] = p1;
				p1->next = cell_head[m][y][x];
			}
		}
	}

	for (k = 0; k < TOTAL_MATERIAL; k++)  //initial distribution of raw material
	{
		x = randl(N);
		y = randl(N);
		raw_arr[y][x]++;
	}
}

//------------------------------------------------------------------------------
void inoculate(void)
{
	int k, k1;

	for (k = 0; k < INOCUNUM; k++) {
		x = randl(N);
		y = randl(N);
		p2 = (struct rna*)malloc(LEN);
		if (!p2) { printf("\t%dform_monomer--memeout\n", k + 1); exit(0); }
		for (k1 = 0; k1 < inoculength; k1++) p2->information[0][k1] = inocuseq[k1];
		p2->information[0][k1] = 0;
		p2->information[1][0] = 0;

		p2->length1 = inoculength;
		p2->length2 = 0;
		p2->nick = 0;
		p2->next = cell_head[h][y][x]->next;
		if (p2->next != cell_head[h][y][x])(p2->next)->prior = p2;
		cell_head[h][y][x]->next = p2;
		p2->prior = cell_head[h][y][x];
	}

	for (k = 0; k < INOCUNUM; k++) {
		x = randl(N);
		y = randl(N);
		p2 = (struct rna*)malloc(LEN);
		if (!p2) { printf("\t%dform_monomer--memeout\n", k + 1); exit(0); }
		for (k1 = 0; k1 < inoculength1; k1++) p2->information[0][k1] = inocuseq1[k1];
		p2->information[0][k1] = 0;
		p2->information[1][0] = 0;

		p2->length1 = inoculength1;
		p2->length2 = 0;
		p2->nick = 0;
		p2->next = cell_head[h][y][x]->next;
		if (p2->next != cell_head[h][y][x])(p2->next)->prior = p2;
		cell_head[h][y][x]->next = p2;
		p2->prior = cell_head[h][y][x];
	}

	for (k = 0; k < INOCUNUM; k++) {
		x = randl(N);
		y = randl(N);
		p2 = (struct rna*)malloc(LEN);
		if (!p2) { printf("\t%dform_monomer--memeout\n", k + 1); exit(0); }
		for (k1 = 0; k1 < inoculength2; k1++) p2->information[0][k1] = inocuseq2[k1];
		p2->information[0][k1] = 0;
		p2->information[1][0] = 0;

		p2->length1 = inoculength2;
		p2->length2 = 0;
		p2->nick = 0;
		p2->next = cell_head[h][y][x]->next;
		if (p2->next != cell_head[h][y][x])(p2->next)->prior = p2;
		cell_head[h][y][x]->next = p2;
		p2->prior = cell_head[h][y][x];
	}

	for (k = 0; k < INOCUNUM; k++) {
		x = randl(N);
		y = randl(N);
		p2 = (struct rna*)malloc(LEN);
		if (!p2) { printf("\t%dform_monomer--memeout\n", k + 1); exit(0); }
		for (k1 = 0; k1 < inoculength3; k1++) p2->information[0][k1] = inocuseq3[k1];
		p2->information[0][k1] = 0;
		p2->information[1][0] = 0;

		p2->length1 = inoculength3;
		p2->length2 = 0;
		p2->nick = 0;
		p2->next = cell_head[h][y][x]->next;
		if (p2->next != cell_head[h][y][x])(p2->next)->prior = p2;
		cell_head[h][y][x]->next = p2;
		p2->prior = cell_head[h][y][x];
	}
}

//------------------------------------------------------------------------------
void unit_case(void)      // Action of units (molecules) in the system
{
	int a, b, d, j, k, m, n, randnt, raw_bef, nt_turn;
	double f, rtdaddlig, rtdaddphili;

	avail_xy_init();      // Initialization for xy_choose
	for (d = 0; d < CELLNUM; d++)
	{
		xy_choose();    // Picks a cell at random 
		raw_bef = raw_arr[y][x];
		for (k = 0; k < raw_bef; k++)
		{
			randcaser = randl(2);
			switch (randcaser)
			{
			case 0:  //forming nt
				if (randd() < PMF)
				{
					raw_arr[y][x]--;
					p3 = (struct rna*)malloc(LEN);
					if (!p3) { printf("\t%d form_monomer--memeout\n", k + 1); exit(0); }
					randnt = randl(4) + 1;
					switch (randnt)
					{
					case 1:  p3->information[0][0] = A; break;
					case 2:  p3->information[0][0] = C; break;
					case 3:  p3->information[0][0] = G; break;
					case 4:  p3->information[0][0] = U; break;
					default: printf("form randnt error");
					}
					p3->information[0][1] = 0;
					p3->information[1][0] = 0;

					p3->length1 = 1;
					p3->length2 = 0;
					p3->nick = 0;

					p3->prior = cell_head[!h][y][x];
					p3->next = cell_head[!h][y][x]->next;
					if (p3->next != cell_head[!h][y][x])(p3->next)->prior = p3;
					cell_head[!h][y][x]->next = p3;
				}
				break;

			case 1:   // raw moving
				if (randd() < PMOVR)
				{
					randcaser1 = randl(4);   // Four possible directions
					switch (randcaser1)
					{
					case 0:
						if (x > 0)
						{
							raw_arr[y][x]--;
							raw_arr[y][x - 1]++;
						}
						break;

					case 1:
						if (x < N - 1)
						{
							raw_arr[y][x]--;
							raw_arr[y][x + 1]++;
						}
						break;

					case 2:
						if (y > 0)
						{
							raw_arr[y][x]--;
							raw_arr[y - 1][x]++;
						}
						break;

					case 3:
						if (y < N - 1)
						{
							raw_arr[y][x]--;
							raw_arr[y + 1][x]++;
						}
						break;

					default:printf("raw moving error");
					}
				}
				break;

			default:printf("raw case error");
			}
		}


		for (p = cell_head[h][y][x]->next; p != cell_head[h][y][x]; p = p->next)
		{
			randcase = randl(6);
			switch (randcase)
			{
			case 0:                        // Chain ligation with mineral catalysis
				for (p3 = p->next; p3 != p; p3 = p3->next)
				{
					if (p3 == cell_head[h][y][x]) { p3 = cell_head[h][y][x]->next; if (p3 == p)break; }
					if (p3->length2 == 0)
					{
						if (randd() < PLMC / p3->length1)
						{
							if (p->length1 + p3->length1 > MAX_RNA_LENGTH - 1)
							{
								over_max_len++; continue;
							}

							for (a = 0; a < p3->length1; a++)
								p->information[0][a + p->length1] = p3->information[0][a];
							p->information[0][p->length1 + p3->length1] = 0;
							p->length1 = p->length1 + p3->length1;

							(p3->prior)->next = p3->next;
							if (p3->next != cell_head[h][y][x])(p3->next)->prior = p3->prior;
							free(p3);

							break;
						}
					}
				}
				fresh_unit();
				break;

			case 1:            // Decay and degradation 
				if (p->length1 == 1)  // Decay of mononucleotide
				{
					if (p->length2 == 0 && randd() < PMD)
					{
						raw_arr[y][x]++;
						(p->prior)->next = p->next;
						if (p->next != cell_head[h][y][x])(p->next)->prior = p->prior;
						p3 = p;
						p = p->prior;
						free(p3); break;
					}
				}
				else                  //Degradation of chain
				{
					if (p->length1 > p->length2 && randd() < PNDE) // Ma--start    Nucleotide residue at RNA end decaying
					{
						raw_arr[y][x]++;
						p->information[0][p->length1 - 1] = 0;
						p->length1--;
						if (p->length1 == 1) {fresh_unit(); break; }
					} // Ma-end

					f = PBB;
					for (j = p->length1; j > 1; j--)
					{
						if (j <= p->length2)   // Falling into double chain region
						{
							if (p->nick == 0)
							{
								m = j - 1;
								n = p->length2 - j + 1;
								k = (m < n) ? m : n;
								f = PBB * k * PBB * k;
							}
							else
							{
								if (j == p->nick + 1)
								{
									m = j - 1;
									n = p->length2 - j + 1;
									k = (m < n) ? m : n;
									f = PBB * k;
								}
								else if (j > p->nick + 1)
								{
									m = j - p->nick - 1;
									n = p->length2 - j + 1;
									k = (m < n) ? m : n;
									f = PBB * k * PBB * k;
								}
								else
								{
									m = j - 1;
									n = p->nick - j + 1;
									k = (m < n) ? m : n;
									f = PBB * k * PBB * k;
								}
							}
						}

						if (randd() < f)
						{
							p3 = (struct rna*)malloc(LEN);
							if (!p3) { printf("\t%ddeg--memeout\n", i); exit(0); }

							for (b = 0; b < p->length1 - j + 1; b++)
								p3->information[0][b] = p->information[0][b + j - 1];
							p3->information[0][p->length1 - j + 1] = 0;
							p->information[0][j - 1] = 0;
							p3->length1 = p->length1 - j + 1;
							p->length1 = j - 1;

							if (p->length2 > j - 1)
							{
								for (b = 0; b < p->length2 - j + 1; b++)
									p3->information[1][b] = p->information[1][b + j - 1];
								p3->information[1][p->length2 - j + 1] = 0;
								p->information[1][j - 1] = 0;
								p3->length2 = p->length2 - j + 1;
								p->length2 = j - 1;
							}
							else
							{
								p3->information[1][0] = 0;
								p3->length2 = 0;
							}

							if (p->nick > j - 1) { p3->nick = p->nick - j + 1; p->nick = 0; }
							else if (p->nick == j - 1) { p3->nick = 0; p->nick = 0; }
							else p3->nick = 0;

							p3->prior = cell_head[!h][y][x];
							p3->next = cell_head[!h][y][x]->next;
							if (p3->next != cell_head[!h][y][x])(p3->next)->prior = p3;
							cell_head[!h][y][x]->next = p3;
							break;
						}
					}
				}
				fresh_unit();
				break;

			case 2:                         //Template-directed addition
				if (p->nick == 0)                      //Template-directed attraction of substrates
				{
					for (p3 = p->next; p3 != p; p3 = p3->next)
					{
						if (p3 == cell_head[h][y][x]) { p3 = cell_head[h][y][x]->next; if (p3 == p)break; }
						if (p3->length2 == 0)
						{
							if (p3->length1 <= p->length1 - p->length2)
							{
								for (flag = 0, b = 0; b < p3->length1; b++)
								{
									if ((p3->information[0][p3->length1 - 1 - b] + p->information[0][p->length2 + b]) == 5)continue;
									else if (randd() < PFP)continue;
									else { flag = 1; break; }
								}
								if (flag == 0)
								{
									rtdaddphili = randd();
									if (rtdaddphili < PAT)
									{
										for (a = 0; a < p3->length1; a++)
											p->information[1][p->length2 + a] = p3->information[0][p3->length1 - 1 - a];
										p->information[1][p->length2 + p3->length1] = 0;
										if (p->length2 != 0)p->nick = p->length2;
										p->length2 = p->length2 + p3->length1;

										(p3->prior)->next = p3->next;
										if (p3->next != cell_head[h][y][x])(p3->next)->prior = p3->prior;
										free(p3);
										break;
									}
								}
							}
						}
					}
				}
				else                     //Template-directed ligation
				{
					rtdaddlig = randd();
					if (rtdaddlig < PLT)
					{
						p->nick = 0;
					}
				}
				fresh_unit();
				break;

			case 3:                           // Separation
				if (p->length2 != 0)    // Separation of double chain  
				{
					if (randd() < pow(PSBP, p->length2 - p->nick))
					{
						p3 = (struct rna*)malloc(LEN);
						if (!p3) { printf("\t%dsep--memeout\n", i); exit(0); }
						for (b = 0; b < p->length2 - p->nick; b++)
							p3->information[0][b] = p->information[1][p->length2 - 1 - b];
						p->information[1][p->nick] = 0;

						p3->information[0][p->length2 - p->nick] = 0;
						p3->information[1][0] = 0;

						p3->length1 = p->length2 - p->nick;
						p->length2 = p->nick;
						p->nick = 0;
						p3->length2 = 0;
						p3->nick = 0;

						p3->prior = cell_head[!h][y][x];
						p3->next = cell_head[!h][y][x]->next;
						if (p3->next != cell_head[!h][y][x])(p3->next)->prior = p3;
						cell_head[!h][y][x]->next = p3;
					}
				}
				fresh_unit();
				break;

				//-------------------------------------------------------------------------------------
			case 4:
				if (p->length2 == 0 && p->length1<1.5* nrlength)  //Ma-2
				{
					flag = findseq(nrseq, nrlength, p);
					if (flag == 0)    // nt-synthetase catalyses the synthesis of nt. 
					{
						nt_turn = TNSS;
						raw_bef = raw_arr[y][x];
						for (k = 0; k < raw_bef; k++)
						{
							if (nt_turn <= 0)break;
							nt_turn--;
							if (randd() < PMFS)
							{
								raw_arr[y][x]--;

								p3 = (struct rna*)malloc(LEN);
								if (!p3) { printf("\t%dsyn form_monomer--memeout\n", k + 1); exit(0); }
								randnt = randl(4) + 1;
								switch (randnt)
								{
								case 1:  p3->information[0][0] = A; break;
								case 2:  p3->information[0][0] = C; break;
								case 3:  p3->information[0][0] = G; break;
								case 4:  p3->information[0][0] = U; break;
								default: printf("syn randnt error");
								}
								p3->information[0][1] = 0;
								p3->information[1][0] = 0;

								p3->length1 = 1;
								p3->length2 = 0;
								p3->nick = 0;

								p3->prior = cell_head[!h][y][x];
								p3->next = cell_head[!h][y][x]->next;
								if (p3->next != cell_head[!h][y][x])(p3->next)->prior = p3;
								cell_head[!h][y][x]->next = p3;
							}
						}
					}
				}
				fresh_unit();
				break;

			case 5:            //moving to another adjacent cell
				if (randd() * FDMOV < PMOV)
				{
					randcase1 = randl(4);   // Four possible directions
					switch (randcase1)
					{
					case 0:
						if (x > 0)
						{
							p1 = p->prior;
							p2 = p->next;

							p3 = cell_head[!h][y][x - 1]->next;
							cell_head[!h][y][x - 1]->next = p;
							p->next = p3;
							p->prior = cell_head[!h][y][x - 1];
							if (p3 != cell_head[!h][y][x - 1])p3->prior = p;

							p1->next = p2;
							if (p2 != cell_head[h][y][x])p2->prior = p1;
							p = p1;
						}
						else fresh_unit();
						break;

					case 1:
						if (x < N - 1)
						{
							p1 = p->prior;
							p2 = p->next;

							p3 = cell_head[!h][y][x + 1]->next;
							cell_head[!h][y][x + 1]->next = p;
							p->next = p3;
							p->prior = cell_head[!h][y][x + 1];
							if (p3 != cell_head[!h][y][x + 1])p3->prior = p;

							p1->next = p2;
							if (p2 != cell_head[h][y][x])p2->prior = p1;
							p = p1;
						}
						else fresh_unit();
						break;

					case 2:
						if (y > 0)
						{
							p1 = p->prior;
							p2 = p->next;

							p3 = cell_head[!h][y - 1][x]->next;
							cell_head[!h][y - 1][x]->next = p;
							p->next = p3;
							p->prior = cell_head[!h][y - 1][x];
							if (p3 != cell_head[!h][y - 1][x])p3->prior = p;

							p1->next = p2;
							if (p2 != cell_head[h][y][x])p2->prior = p1;
							p = p1;
						}
						else fresh_unit();
						break;


					case 3:
						if (y < N - 1)
						{
							p1 = p->prior;
							p2 = p->next;

							p3 = cell_head[!h][y + 1][x]->next;
							cell_head[!h][y + 1][x]->next = p;
							p->next = p3;
							p->prior = cell_head[!h][y + 1][x];
							if (p3 != cell_head[!h][y + 1][x])p3->prior = p;

							p1->next = p2;
							if (p2 != cell_head[h][y][x])p2->prior = p1;
							p = p1;
						}
						else fresh_unit();
						break;

					default:printf("rna moving error");
					}
				}
				else fresh_unit();
				break;

			default: printf("rna case error");
			}
		}
	}
}

//------------------------------------------------------------------------------
int record(void)               // Data recording
{
	FILE* fptxt, * fptxt1;  //Ma-start
	errno_t err, err1;
	err = fopen_s(&fptxt, "monitor.txt", "at");
	if (err != 0) { printf("cannot open file");  exit(-1); }
	err1 = fopen_s(&fptxt1, "picture.txt", "at");
	if (err1 != 0) { printf("cannot open file1");  exit(-1); }

	int ch_num[MAX_RNA_LENGTH], dch_num[MAX_RNA_LENGTH], long_chain_num, si; //Ma-end

	nr_num[g] = 0;
	ctr1_num[g] = 0;
	ctr2_num[g] = 0;
	ctr3_num[g] = 0;

	total_mat_num[g] = 0;
	RNA_num[g] = 0;
	raw_num[g] = 0;
	for (y = 0; y < N; y++)
	{
		for (x = 0; x < N; x++)
		{
			raw_num[g] += raw_arr[y][x];
			for (p = cell_head[h][y][x]->next; p != cell_head[h][y][x]; p = p->next)
			{
				RNA_num[g]++;
				total_mat_num[g] += p->length1 + p->length2;
				flag1 = findseq(nrseq, nrlength, p);
				if (flag1 == 0) nr_num[g]++; //Nr sequence in chain1
				flag1 = findseq(ctr1seq, ctr1length, p);
				if (flag1 == 0) ctr1_num[g]++; //
				flag1 = findseq(ctr2seq, ctr2length, p);
				if (flag1 == 0) ctr2_num[g]++; //
				flag1 = findseq(ctr3seq, ctr3length, p);
				if (flag1 == 0) ctr3_num[g]++; //
			}
		}
	}
	total_mat_num[g] += raw_num[g];

	//Ma
	fprintf(fptxt1, "step=%d: nr=%d, ctr1=%d, ctr2=%d, ctr3=%d, RNA=%d  (tn=%d, r=%d)\n", i,
		(int)nr_num[g], (int)ctr1_num[g], (int)ctr2_num[g], (int)ctr3_num[g], (int)RNA_num[g], (int)total_mat_num[g], (int)raw_num[g]);

	for (si = 0; si < LONG_CHAIN_LEN; si++)  //Ma---start
	{
		ch_num[si] = 0;
		dch_num[si] = 0;   // double chains
	}

	long_chain_num = 0;
	for (y = 0; y < N; y++)
	{
		for (x = 0; x < N; x++)
		{
			for (p = cell_head[h][y][x]->next; p != cell_head[h][y][x]; p = p->next)
			{
				if (p->length1 > LONG_CHAIN_LEN) p - long_chain_num++;
				else
				{
					ch_num[p->length1 - 1]++;
					if (p->length2 != 0) dch_num[p->length1 - 1]++;
				}
			}
		}
	}

	printf("\nstep=%d: nr=%d, ctr1=%d, ctr2=%d, ctr3=%d, RNA=%d  (tn=%d, r=%d)\n", i,
		(int)nr_num[g], (int)ctr1_num[g], (int)ctr2_num[g], (int)ctr3_num[g], (int)RNA_num[g], (int)total_mat_num[g], (int)raw_num[g]);
	fprintf(fptxt, "\nstep=%d: nr=%d, ctr1=%d, ctr2=%d, ctr3=%d, RNA=%d  (tn=%d, r=%d)\n", i,
		(int)nr_num[g], (int)ctr1_num[g], (int)ctr2_num[g], (int)ctr3_num[g], (int)RNA_num[g], (int)total_mat_num[g], (int)raw_num[g]);
	printf(" chains longer than %d = %d\n", LONG_CHAIN_LEN, long_chain_num);
	fprintf(fptxt, " chains longer than %d = %d\n", LONG_CHAIN_LEN, long_chain_num);

	for (si = 0; si < LONG_CHAIN_LEN; si++)
	{
		printf("%dnt-%d(%d), ", si + 1, ch_num[si], dch_num[si]);
		fprintf(fptxt, "%dnt-%d(%d), ", si + 1, ch_num[si], dch_num[si]);
	}

	printf("\n");
	fprintf(fptxt, "\n");
	//Ma---end

	g++;
	fclose(fptxt);
	fclose(fptxt1);
	return(0);
}

//------------------------------------------------------------------------------
void freepool(void)        // Memory releasing  
{
	int m;
	for (m = 0; m < 2; m++)
	{
		for (y = 0; y < N; y++)
		{
			for (x = 0; x < N; x++)
			{
				while (1)
				{
					if (cell_head[m][y][x]->next != cell_head[m][y][x])
					{
						p = cell_head[m][y][x]->next;
						cell_head[m][y][x]->next = p->next;
						free(p);
					}
					else break;
				}
				free(cell_head[m][y][x]);
			}
		}
	}
}

//------------------------------------------------------------------------------ 
int main()
{
	inits();        // initialization of the system

	time1 = time(ptimer);

	for (i = 0; i <= STEPNUM; i++)      // Monte-Carlo cycle
	{
		if (i == INOCUSTEP)inoculate();
		if (i >= STAREC && i % RECINT == 0)
		{
			record();
		}
		unit_case();
		h = !h;
	}
	
	time2 = time(ptimer);

	FILE* fptxt;  //Ma-start
	errno_t err;
	err = fopen_s(&fptxt, "picture.txt", "at");
	fprintf(fptxt,"execution time:%i, seed:%i\n\n",time2-time1,SD);

	freepool();
	return (0);
}
//========================================================  End of the program

