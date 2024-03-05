// hr-circle-4-9-2-8-5-8-5.cpp: 主项目文件。
// hr-circle-4-9-2-8-5-8-4.cpp: 主项目文件。
// hr-circle-4-9-2-8-5-8.cpp: 主项目文件。
// hr-circle-4-9-2-8-5-7.cpp: 主项目文件。
// hr-circle-4-9-2-8-5-6.cpp: 主项目文件。
// hr-circle-4-9-2-8-5-5.cpp: 主项目文件。
// hr-circle-4-9-2-8-5-4.cpp: 主项目文件。
// hr-circle-4-9-2-8-5-3.cpp: 主项目文件。
// hr-circle-4-9-2-8-5-2.cpp: 主项目文件。
// hr-circle-4-9-2-8-5-1.cpp: 主项目文件。
// hr-circle-4-9-2-8-5.cpp: 主项目文件。
// hr-circle-4-9-2-8-4.cpp: 主项目文件。
// hr-circle-4-9-2-8-3.cpp: 主项目文件。
// hr-circle-4-9-2-8-2.cpp: 主项目文件。
// hr-circle-4-9-2-8-1.cpp: 主项目文件。
// hr-circle-4-9-2-8.cpp: 主项目文件。
// hr-circle-4-9-2-7-2.cpp: 主项目文件。
// hr-circle-4-9-2-7.cpp: 主项目文件。
// hr-circle-4-9-2-5.cpp: 主项目文件。
// hr-circle-4-9-2-4.cpp: 主项目文件。
// hr-circle-4-9-2-3.cpp: 主项目文件。
// hr-circle-4-9-2-2.cpp: 主项目文件。
// hr-circle-4-9-2-1.cpp: 主项目文件。

// hhr-circle-4-9-2-1.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
// hhr-circle-4-9-2.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
// hhr-circle-4-9-1.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
// hhr-circle-4-9.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
// hhr-circle-4-8.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
// hhr-circle-4-7.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
// hhr-circle-4-6.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
// hhr-circle-4-5.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
// hhr-circle-4-4.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
// hhr-circle-4-3.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
// hhr-circle-4-2.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
// hhr-circle-4.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
// hhr-circle-3.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
// hhr-circle-2.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
// hhr-circle-1.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
// circle-1.cpp: 主项目文件。

//#include <stdafx.h>
// 
//"Nucleotide synthetase ribozymes may have emerged first in the RNA World"
//by Wentao Ma
//-----C source codes for the simulation program
//----Parameters are for the case in figure 2a of the article

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <conio.h>
#include <string.h>

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
#define STEPNUM 10000000     // Total time steps of Monte Carlo simulation
#define STAREC 0          // The step to start record
#define RECINT 10000       // The interval steps of recording
#define MAX_RNA_LENGTH 100    // Defining maximum RNA length allowed in the simulation
#define LONG_CHAIN_LEN 30     // Defining long chains for recording

#define SD 11  //13   //5
#define N 30  //20    //Ma2-1            // The side length of the two-dimensional grid
#define TOTAL_MATERIAL 80000 //40000 // Total materials in the system
#define HRSEQ A,G,U,C
/*
#define NRSEQ U,C,G,C,G,A    // The presumed specific sequence with which a polynucleotide could act as a nt-synthetase
#define CONTRSEQ1 U,C,U,A,G,A  // 
#define CONTRSEQ2 G,U,U,A,A,C  // 
#define CONTRSEQ3 A,C,A,U,G,U  // 
#define HR_NRSEQ U,C,G,C,G,A,A,G,U,C
#define HR_CONTRSEQ1 U,C,U,A,G,A,A,G,U,C
#define HR_CONTRSEQ2 G,U,U,A,A,C,A,G,U,C
#define NHR_CONTRSEQ3 A,C,A,U,G,U,C,G,C,A
*/
//Ma2-1 start
#define MAX_CHAR_LENGTH 20   //RRRRR-ma3

#define NRSEQ U,G,A,U,G,C,A,G   // The presumed specific sequence with which a polynucleotide could act as a nt-synthetase
#define CONTRSEQ1 A,C,U,G,A,C,G,U  //  
#define CONTRSEQ2 A,A,C,G,C,U,C,G  // 
#define CONTRSEQ3 C,G,A,U,C,A,A,U
//Ma2-1 end

#define INOCUSEQ NRSEQ                                         //new
#define INOCUNUM 20   //Ma2-1
#define INOCUSEQ1 CONTRSEQ1 // Ma2-1
#define INOCUNUM1 INOCUNUM
#define INOCUSEQ2 CONTRSEQ2 // Ma2-1 
#define INOCUNUM2 INOCUNUM
#define INOCUSEQ3 CONTRSEQ3 // Ma2-1
#define INOCUNUM3 INOCUNUM
#define INOCUSTEP 10000
#define PSBP 0.5  //0.9   // Ma2-1           // Probability of separation of a base-pair
#define PBB 0.00001              // Probability of breaking of a phosphodiester bond
#define PLMC 0.0000001  //0.000002 //Ma2-1        // Probability of ligation of two unit in a cell with mineral catalysis 
#define PEL  0.0000001  //0.00001  //Ma2-1  

#define MINI_CIRCLE_LENGTH 8 //Ma2
#define FDA 5  //Ma2: Factor concerning de novo attraction of substrate by RNA template
#define FHR 1   //1000  // Ma2-1: Factor concerning self-cleavage or self-ligation of HR
#define FLT 0.5   //1  //0.1  //0.01  //Ma2-1 Factor concerning linear RNA acting as template

#define PAT 0.5               // Probability of attracting a substrate by a template when the substrate could base-pair with the template
#define PLT 0.2   //0.5  //0.9    //Ma2-1         // Probability of a template-directed ligation 
#define PMR 0.01  // Ma2-1
#define PMV 0.002  // Ma2-1           // Probability of movement of raws or nucleotides
#define PMF 0.001 //0.0002 // Ma2-1          // probability of mononucleotide formation from raw materials
#define PMFS 0.9            // probability of mononucleotide formation from raw materials under the catalysis of nt-synthetase
#define TNSS 1              // Turn of nt-synthesis by nt-synthetase each step
#define PMD 0.05 //0.01   // Ma2-1          // probability of mononucleotide decay to raw materials
#define PNDE 0.001  //0.0001 // Ma2-1   // Probability of a nucleotide residue decaying at RNA?s chain end
#define PFP 0.01             // probability of false base-pairing
#define RMRW (pow(p->length1+p->length2,1/2.0))  //Ma2: Zimm model  
							 // The relationship between the movement of RNAs and their molecular weight
#define CELLNUM (N*N)   // Cell numbers in the grid

long randl(long);      // random number between 0 and parameter 
double randd(void);    // random double between 0 and 1         
void avail_xy_init(void); // Initialization for xy_choose
void xy_choose(void);     // Picks a cell at random 
void fresh_unit(void);    // Updating a unit for the next time step
int findseq(char seq[], int seqlength, struct rna* p); //find a specific subsequence in a sequence 
void inits(void);         // initialization of the system
void inoculate(char);
void unit_case(void);     // Action of units (molecules) in the system
void record(void);         // Data recording
void freepool(void);      // Memory releasing

struct rna                // A unit of mononucleotide or polynucleotide
{
	char information[2][MAX_RNA_LENGTH];
	int length1;
	int length2;
	struct c2_frag* chain2;
	char type1;
	char type2;
	struct rna* next;
	struct rna* prior;
};
struct rna* room_head[2][N][N];
struct rna* p, * p1, * p2, * p3, * p4, * ps, * ps1, * ps2;

struct c2_frag
{
	int start;
	int length;
	struct c2_frag* next;
	struct c2_frag* prior;
};
struct c2_frag* chain2;
struct c2_frag* c2f, * c2f1, * c2f2, * c2f3, * c2f4;

static char nrseq[50] = { NRSEQ };        // Presumed nt_synthetase sequence
static char hrseq[50] = { HRSEQ };		// Presumed hammerhead ribozyme sequence
static char contrseq1[50] = { CONTRSEQ1 };
static char contrseq2[50] = { CONTRSEQ2 };
static char contrseq3[50] = { CONTRSEQ3 };
//static char hr_nrseq[50] = { HR_NRSEQ };
//static char hr_contrseq1[50] = { HR_CONTRSEQ1 };
//static char hr_contrseq2[50] = { HR_CONTRSEQ2 };
//static char nhr_contrseq3[50] = { NHR_CONTRSEQ3 };
static char inocuseq[50] = { INOCUSEQ };
static char inocuseq1[50] = { INOCUSEQ1 };
static char inocuseq2[50] = { INOCUSEQ2 };
static char inocuseq3[50] = { INOCUSEQ3 };
static int raw_arr[N][N];
char temp_information[2][MAX_RNA_LENGTH];

int over_max_len = 0;
int x, y;                 // The coordinate of cells in the grid 
int nrlength, hrlength, contrlength1, contrlength2, contrlength3;
int inoculength, inoculength1, inoculength2, inoculength3;
int randcase, randcase1, randcaser, randcaser1, length3, g = 0, h = 0, g_end = 0, gi, g_stop = 0;
int flag, flag1, flag2, flag3, flag4, flag5;
long i;                  // Cycle variable for Monte Carlo steps
long available;
long availabl[CELLNUM];
long recstep[(STEPNUM - STAREC) / RECINT + 1];    // Record steps
float nr[(STEPNUM - STAREC) / RECINT + 1];  // Record number of nt-synthetases in steps
float cir_nr[(STEPNUM - STAREC) / RECINT + 1];
float hr[(STEPNUM - STAREC) / RECINT + 1];  // Record number of hr in steps
float cir_hr[(STEPNUM - STAREC) / RECINT + 1];
float hr_nr[(STEPNUM - STAREC) / RECINT + 1];  // Record number of hr_nr in steps
float cir_hr_nr[(STEPNUM - STAREC) / RECINT + 1];

float ctr1[(STEPNUM - STAREC) / RECINT + 1];  // 
float cir_ctr1[(STEPNUM - STAREC) / RECINT + 1];  //
float ctr2[(STEPNUM - STAREC) / RECINT + 1];  // 
float cir_ctr2[(STEPNUM - STAREC) / RECINT + 1];  // 
float ctr3[(STEPNUM - STAREC) / RECINT + 1];  // 
float cir_ctr3[(STEPNUM - STAREC) / RECINT + 1];  // 

float total_mat_num[(STEPNUM - STAREC) / RECINT + 1];  // Record number of total materials in steps
float unit[(STEPNUM - STAREC) / RECINT + 1];  // Record number of units in steps
float cir_unit[(STEPNUM - STAREC) / RECINT + 1];  // 

float raw_num[(STEPNUM - STAREC) / RECINT + 1];  // Record number of raw in steps
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
	p3 = room_head[!h][y][x]->next;
	room_head[!h][y][x]->next = p;
	p->next = p3;
	p->prior = room_head[!h][y][x];
	if (p3 != room_head[!h][y][x])p3->prior = p;
	p1->next = p2;
	if (p2 != room_head[h][y][x])p2->prior = p1;
	p = p1;
}

//------------------------------------------------------------------------------
int findseq(char seq[], int seqlength, struct rna* p)  // Find a specific subsequence in a sequence
{
	int flag2, a, b;
//	char inf[MAX_RNA_LENGTH];
//	for (a = 0; a < MAX_RNA_LENGTH; a++)inf[a] = 0;

	char inf[MAX_RNA_LENGTH + MAX_CHAR_LENGTH];    //RRRRR-ma3,  #define MAX_CHAR_LENGTH 20, --- to avoid using array out of bounds in the "dangerous block"
	for (a = 0; a < MAX_RNA_LENGTH + MAX_CHAR_LENGTH; a++)inf[a] = 0;  //RRRRR-ma3, 

	for (a = 0; a < p->length1 + seqlength; a++)  // dangerous block
	{
		if (a < p->length1) inf[a] = p->information[0][a];
		else inf[a] = p->information[0][a - p->length1];
	}

	flag2 = 0;
	// search for the subsequence
	if (p->length1 >= seqlength)
	{
		if (p->type1 == 0)
		{
			for (b = 0; p->length1 - seqlength - b >= 0; b++)
			{
				flag2 = 0;
				for (a = 0; a < seqlength; a++)
				{
					if (inf[b + a] == seq[a])continue;
					else { flag2 = 1; break; }  //this location has not this subsequence
				}
				if (flag2 == 0)break; // this location has this subsequence
			}
		}
		else if (p->type1 == 1)
		{

			for (b = 0; b <= p->length1 - 1; b++)
			{
				flag2 = 0;
				for (a = 0; a < seqlength; a++)
				{
					if (inf[b + a] == seq[a])continue;
					else { flag2 = 1; break; }
				}
				if (flag2 == 0)break;
			}
		}
	}
	else flag2 = 1;

	if (flag2 == 0)return(0);   //Yes, the sequence contains the subsequence
	else return(1);   //no, the sequence does not contain the subsequence
}

//------------------------------------------------------------------------------
void inits(void)         // initialization of the system
{
	int j, m, k;
	seed(SD);

	nrlength = 0;
	for (j = 0; nrseq[j] != 0; j++)
		nrlength++;

	hrlength = 0;                                               //new
	for (j = 0; hrseq[j] != 0; j++)
		hrlength++;

	contrlength1 = 0;
	for (j = 0; contrseq1[j] != 0; j++)
		contrlength1++;

	contrlength2 = 0;
	for (j = 0; contrseq2[j] != 0; j++)
		contrlength2++;

	contrlength3 = 0;
	for (j = 0; contrseq3[j] != 0; j++)
		contrlength3++;

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
				room_head[m][y][x] = p1;
				p1->next = room_head[m][y][x];
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
void inoculate(char ch_typ)   //Ma2-1   inoculating circular or linear RNAs
{
	int k, k1;

	for (k = 0; k < INOCUNUM; k++) {
		x = randl(N);
		y = randl(N);
		p2 = (struct rna*)malloc(LEN);
		if (!p2) { printf("\t%dform_monomer--memeout\n", k + 1); exit(0); }
		memset(p2->information, 0, sizeof(p2->information));
		for (k1 = 0; k1 < inoculength; k1++) p2->information[0][k1] = inocuseq[k1];
		p2->information[0][k1] = 0;
		p2->information[1][0] = 0;

		p2->length1 = inoculength;
		p2->length2 = 0;
		p2->type1 = ch_typ;
		p2->type2 = 0;
		p2->next = room_head[h][y][x]->next;
		if (p2->next != room_head[h][y][x])(p2->next)->prior = p2;
		room_head[h][y][x]->next = p2;
		p2->prior = room_head[h][y][x];

		c2f = (struct c2_frag*)malloc(sizeof(struct c2_frag));
		if (!c2f) { printf("\tinit1--memeout\n"); exit(0); }
		p2->chain2 = c2f;
		c2f->next = p2->chain2;
		c2f->prior = p2->chain2;
	}

	for (k = 0; k < INOCUNUM1; k++) {
		x = randl(N);
		y = randl(N);
		p2 = (struct rna*)malloc(LEN);
		if (!p2) { printf("\t%dform_monomer--memeout\n", k + 1); exit(0); }
		memset(p2->information, 0, sizeof(p2->information));
		for (k1 = 0; k1 < inoculength1; k1++) p2->information[0][k1] = inocuseq1[k1];
		p2->information[0][k1] = 0;
		p2->information[1][0] = 0;

		p2->length1 = inoculength1;
		p2->length2 = 0;
		p2->type1 = ch_typ;
		p2->type2 = 0;
		p2->next = room_head[h][y][x]->next;
		if (p2->next != room_head[h][y][x])(p2->next)->prior = p2;
		room_head[h][y][x]->next = p2;
		p2->prior = room_head[h][y][x];

		c2f = (struct c2_frag*)malloc(sizeof(struct c2_frag));
		if (!c2f) { printf("\tinit1--memeout\n"); exit(0); }
		p2->chain2 = c2f;
		c2f->next = p2->chain2;
		c2f->prior = p2->chain2;
	}

	for (k = 0; k < INOCUNUM2; k++) {
		x = randl(N);
		y = randl(N);
		p2 = (struct rna*)malloc(LEN);
		if (!p2) { printf("\t%dform_monomer--memeout\n", k + 1); exit(0); }
		memset(p2->information, 0, sizeof(p2->information));
		for (k1 = 0; k1 < inoculength2; k1++) p2->information[0][k1] = inocuseq2[k1];
		p2->information[0][k1] = 0;
		p2->information[1][0] = 0;

		p2->length1 = inoculength2;
		p2->length2 = 0;
		p2->type1 = ch_typ;
		p2->type2 = 0;
		p2->next = room_head[h][y][x]->next;
		if (p2->next != room_head[h][y][x])(p2->next)->prior = p2;
		room_head[h][y][x]->next = p2;
		p2->prior = room_head[h][y][x];

		c2f = (struct c2_frag*)malloc(sizeof(struct c2_frag));
		if (!c2f) { printf("\tinit1--memeout\n"); exit(0); }
		p2->chain2 = c2f;
		c2f->next = p2->chain2;
		c2f->prior = p2->chain2;
	}

	for (k = 0; k < INOCUNUM3; k++) {
		x = randl(N);
		y = randl(N);
		p2 = (struct rna*)malloc(LEN);
		if (!p2) { printf("\t%dform_monomer--memeout\n", k + 1); exit(0); }
		memset(p2->information, 0, sizeof(p2->information));
		for (k1 = 0; k1 < inoculength3; k1++) p2->information[0][k1] = inocuseq3[k1];
		p2->information[0][k1] = 0;
		p2->information[1][0] = 0;

		p2->length1 = inoculength3;
		p2->length2 = 0;
		p2->type1 = ch_typ;
		p2->type2 = 0;
		p2->next = room_head[h][y][x]->next;
		if (p2->next != room_head[h][y][x])(p2->next)->prior = p2;
		room_head[h][y][x]->next = p2;
		p2->prior = room_head[h][y][x];

		c2f = (struct c2_frag*)malloc(sizeof(struct c2_frag));
		if (!c2f) { printf("\tinit1--memeout\n"); exit(0); }
		p2->chain2 = c2f;
		c2f->next = p2->chain2;
		c2f->prior = p2->chain2;
	}
}

//------------------------------------------------------------------------------
void unit_case(void)      // Action of units (molecules) in the system
{
	int a, b, c, d, j, k, randnt, raw_bef, nt_turn, length, randseq;
	double f, f1, rtdaddlig, rtdaddphili;

	avail_xy_init();      // Initialization for xy_choose
	for (d = 0; d < CELLNUM; d++)
	{
		xy_choose();    // Picks a cell at random 
		raw_bef = raw_arr[y][x];     //Events of raw materials
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
					memset(p3->information, 0, sizeof(p3->information));
					randnt = randl(4) + 1;
					switch (randnt)
					{
					case 1:  p3->information[0][0] = A; break;
					case 2:  p3->information[0][0] = C; break;
					case 3:  p3->information[0][0] = G; break;
					case 4:  p3->information[0][0] = U; break;
					default: printf("form randnt error");
					}
					//p3->information[0][1] = 0;  //Ma2
					//p3->information[1][0] = 0;  //Ma2

					p3->length1 = 1;
					p3->length2 = 0;
					p3->type1 = 0;
					p3->type2 = 0;

					p3->prior = room_head[!h][y][x];
					p3->next = room_head[!h][y][x]->next;
					if (p3->next != room_head[!h][y][x])(p3->next)->prior = p3;
					room_head[!h][y][x]->next = p3;

					c2f = (struct c2_frag*)malloc(sizeof(struct c2_frag));
					if (!c2f) { printf("\tinit1--memeout\n"); exit(0); }
					p3->chain2 = c2f;
					c2f->next = p3->chain2;
					c2f->prior = p3->chain2;
				}
				break;

			case 1:   // raw moving
				if (randd() < PMR)   //Ma2: with toroidal topology to avoid edge effects
				{
					raw_arr[y][x]--;   //Ma2
					randcaser1 = randl(4);   // Four possible directions
					switch (randcaser1)
					{
					case 0:
						raw_arr[y][(N + x - 1) % N]++;  //Na2: toroidal topology
						break;
					case 1:
						raw_arr[y][(x + 1) % N]++;
						break;
					case 2:
						raw_arr[(N + y - 1) % N][x]++;
						break;
					case 3:
						raw_arr[(y + 1) % N][x]++;
						break;
					default: printf("raw moving error");
					}
				}
				break;

			default:printf("raw case error");
			}
		}

		for (p = room_head[h][y][x]->next; p != room_head[h][y][x]; p = p->next)
		{
			randcase = randl(6);
			switch (randcase)
			{
			case 0:                        // Chain ligation with mineral catalysis
				if (p->type1 == 0)
				{
					if (p->length1 >= MINI_CIRCLE_LENGTH)      //Ma2: considering self-ligation before ligation between RNAs //new  
					{
						f = PEL;
						flag = 0;  //Ma2
						for (a = 0; a < hrlength / 2; a++)  //Ma2
						{
							if (p->information[0][p->length1 - hrlength / 2 + a] == hrseq[a] && p->information[1][p->length1 - hrlength / 2 + a] == 0)continue;
							else { flag = 1; break; }
						}
						if (flag == 0)
						{
							for (a = 0; a < hrlength / 2; a++)
							{
								if (p->information[0][a] == hrseq[hrlength / 2 + a] && p->information[1][a] == 0)continue;
								else { flag = 1; break; }
							}
						}
						if (flag == 0)f = f * FHR;  //Ma2
						if (randd() < f)     //Ma2
						{
							p->type1 = 1;
							fresh_unit();
							break;
						}
					}

					for (p3 = p->next; p3 != p; p3 = p3->next)
					{
						if (p3 == room_head[h][y][x]) { p3 = room_head[h][y][x]->next; if (p3 == p)break; }
						if (p3->type1 == 0)   //Ma2---start: the p3 needs not to be single chain
						{
							if (randd() < PLMC / (p->length1 * p3->length1))   //Ma2: the random ligation should be influenced by the lengths of both strands
							{
								if (p->length1 + p3->length1 > MAX_RNA_LENGTH - 1)
								{
									over_max_len++; continue;
								}

								for (a = 0; a < p3->length1; a++)
								{
									p->information[0][a + p->length1] = p3->information[0][a];
									p->information[1][a + p->length1] = p3->information[1][a];  //Ma2
								}

								for (c2f3 = p3->chain2->next; c2f3 != p3->chain2; c2f3 = c2f3->next)  //Ma2---start
								{
									c2f3->start += p->length1;
								}
								p->chain2->prior->next = p3->chain2->next;
								p3->chain2->next->prior = p->chain2->prior;
								p->chain2->prior = p3->chain2->prior;
								p3->chain2->prior->next = p->chain2;   //Ma2---end

								p->length1 = p->length1 + p3->length1;
								p->length2 = p->length2 + p3->length2;  //Ma2

								(p3->prior)->next = p3->next;
								if (p3->next != room_head[h][y][x])(p3->next)->prior = p3->prior;
								free(p3->chain2);
								free(p3);
								break;
							}
						}
					}
				}
				fresh_unit();
				break;

			case 1:            // Decay and degradation 
				if (p->length1 == 1)  // Decay of mononucleotide
				{
					//if (p->information[0][0] == 0) { printf("unexpected error on p-nucleotide"); exit(0); }
					if (p->length2 == 0)
					{
						if (randd() < PMD)
						{
							raw_arr[y][x]++;
							(p->prior)->next = p->next;
							if (p->next != room_head[h][y][x])(p->next)->prior = p->prior;
							p3 = p;
							p = p->prior;
							free(p3->chain2);
							free(p3); break;
						}
					}
					else if (p->length2 == 1)      //Ma2: The decay of the paired nucleotides at the same time
					{
						if (randd() < PMD * sqrt(PMD))
						{
							raw_arr[y][x] += 2;
							(p->prior)->next = p->next;
							if (p->next != room_head[h][y][x])(p->next)->prior = p->prior;
							p3 = p;
							p = p->prior;
							free(p3->chain2);
							free(p3); break;
						}
					}
					else { printf("unexpected error on chain-length"); exit(0); }
				}
				else                  //Degradation of chain
				{
					if (p->type1 == 0)
					{
						c2f1 = p->chain2->prior;
						
						// Nucleotide residue decaying at the end of RNA
						if (p->information[1][p->length1 - 1] == 0)         // Single chain at the end 
						{
							if (randd() < PNDE)
							{
								p->information[0][p->length1 - 1] = 0;
								p->length1--;
								raw_arr[y][x]++;
							}
						}
						else if (p->information[1][p->length1 - 1] != 0)  //Double chain at the end
						{
							if (randd() < PNDE * sqrt(PNDE))                 // The decay of the paired residues at the same time
							{
								p->information[0][p->length1 - 1] = 0;
								raw_arr[y][x]++;

								p->information[1][p->length1 - 1] = 0;
								raw_arr[y][x]++;

								p->length1--;
								p->length2--;
								c2f1->length--;

								if (c2f1->length == 0)
								{
									(c2f1->prior)->next = c2f1->next;
									(c2f1->next)->prior = c2f1->prior;
									free(c2f1);
								}
							}
						}
						//else { printf("unexpected error on chain-length"); exit(0); }

						if (p->length1 == 1)
						{
							fresh_unit();
							break;
						}

						// Nucleotide residue decaying at the start of RNA
						c2f2 = p->chain2->next;
						if (p->information[1][0] == 0)      //Single chain at the start
						{
							if (randd() < PNDE)             // The decay of the start residue on this single chain
							{
								for (b = 1; b < p->length1; b++)
								{
									p->information[0][b - 1] = p->information[0][b];
									p->information[1][b - 1] = p->information[1][b];
								}
								p->information[0][p->length1 - 1] = 0;
								p->information[1][p->length1 - 1] = 0;
								p->length1--;
								raw_arr[y][x]++;

								for (c2f3 = c2f2; c2f3 != p->chain2; c2f3 = c2f3->next)    //Ma2
								{
									c2f3->start--;
								}
							}
						}
						else if (p->information[1][0] != 0) //Double chain at the start
						{
							if (randd() < PNDE * sqrt(PNDE)) { // The decay of the paired residues at the same time

								for (b = 1; b < p->length1; b++)
								{
									p->information[0][b - 1] = p->information[0][b];
									p->information[1][b - 1] = p->information[1][b];
								}
								p->information[0][p->length1 - 1] = 0;
								p->information[1][p->length1 - 1] = 0;
								raw_arr[y][x]++;
								raw_arr[y][x]++;

								p->length1--;
								p->length2--;
								c2f2->length--;

								for (c2f3 = c2f2->next; c2f3 != p->chain2; c2f3 = c2f3->next)
								{
									c2f3->start--;
								}

								if (c2f2->length == 0)
								{
									(c2f2->prior)->next = c2f2->next;
									(c2f2->next)->prior = c2f2->prior;
									free(c2f2);
								}
							}
						}
						//else { printf("unexpected error on chain-length"); exit(0); }

						if (p->length1 == 1)
						{
							fresh_unit();
							break;
						}

						while (1)
						{
							for (j = p->length1; j > 1; j--)
							{
								f = PBB;
								for (c2f1 = p->chain2->prior; c2f1 != p->chain2; c2f1 = c2f1->prior)
								{
									if (j > c2f1->start + 1)
									{
										if (j <= c2f1->length + c2f1->start)         // Falling into double chain region
										{
											f = f * sqrt(f);
										}
										break;
									}
								}

								if (j - 1 >= hrlength / 2 && p->length1 - j + 1 >= hrlength / 2)       //new
								{
									flag = 0;
									for (a = 0; a < hrlength; a++)
									{
										if (p->information[0][j - 1 - hrlength / 2 + a] == hrseq[a] && p->information[1][j - 1 - hrlength / 2 + a] == 0)continue;
										else { flag = 1; break; }
									}
								}
								else flag = 1;
								if (flag == 0)f = f * FHR;   //Ma2

								c2f1 = p->chain2->prior;
								c2f2 = p->chain2->next;
								if (randd() < f)
								{
									p3 = (struct rna*)malloc(LEN);
									if (!p3) { printf("\t%ddeg--memeout\n", i); exit(0); }
									memset(p3->information, 0, sizeof(p3->information));
									c2f = (struct c2_frag*)malloc(sizeof(struct c2_frag));
									if (!c2f) { printf("\tinit1--memeout\n"); exit(0); }
									p3->chain2 = c2f;
									c2f->next = p3->chain2;
									c2f->prior = p3->chain2;

									for (b = 0; b < p->length1 - j + 1; b++)
									{
										p3->information[0][b] = p->information[0][b + j - 1];
										p->information[0][b + j - 1] = 0;
									}
									p3->length1 = p->length1 - j + 1;
									p->length1 = j - 1;
									p3->length2 = 0;
									p3->type1 = 0;
									p3->type2 = 0;  //Ma2

									//if (p->length2 == 0 || c2f1->start + c2f1->length <= j - 1)    //Ma2
									//{
									//	p3->information[1][0] = 0;
									//}
									if (c2f1 != p->chain2 && c2f1->start + c2f1->length > j - 1)   //Ma2
									{
										//p3->type2 = 0;  //Ma2

										for (b = 0; b < c2f1->start + c2f1->length - j + 1; b++)
										{
											p3->information[1][b] = p->information[1][b + j - 1];
											p->information[1][b + j - 1] = 0;
										}

										while (c2f1 != p->chain2)
										{
											c2f3 = (struct c2_frag*)malloc(sizeof(struct c2_frag));
											if (!c2f3) { printf("\tinit1--memeout\n"); exit(0); }
											c2f3->prior = p3->chain2;
											c2f3->next = p3->chain2->next;
											p3->chain2->next->prior = c2f3;
											p3->chain2->next = c2f3;

											if (j <= c2f1->start + 1)
											{
												c2f3->length = c2f1->length;
												c2f3->start = c2f1->start - j + 1;
												p3->length2 = p3->length2 + c2f3->length;
												p->length2 = p->length2 - c2f3->length;

												(c2f1->prior)->next = c2f1->next;
												(c2f1->next)->prior = c2f1->prior;
												c2f = c2f1;
												c2f1 = c2f1->prior;
												free(c2f);

												if (p3->length1 < p3->length2) { printf("unexpected error on p3-chain-length"); exit(0); }
												if (c2f3->start + c2f3->length > p3->length1) { printf("unexpected error on p3-chain-length"); exit(0); }

												if (j >= c2f1->start + c2f1->length + 1)
												{
													break;
												}
											}
											else
											{
												c2f3->length = c2f1->start + c2f1->length - j + 1;
												c2f3->start = 0;
												c2f1->length = c2f1->length - c2f3->length;

												p3->length2 = p3->length2 + c2f3->length;
												p->length2 = p->length2 - c2f3->length;

												if (p->length1 < p->length2) { printf("unexpected error on p-chain-length"); exit(0); }
												if (p3->length1 < p3->length2) { printf("unexpected error on p3-chain-length"); exit(0); }
												if (c2f3->start + c2f3->length > p3->length1) { printf("unexpected error on p-chain-length"); exit(0); }
												break;
											}
										}
									}
									p3->prior = room_head[!h][y][x];
									p3->next = room_head[!h][y][x]->next;
									if (p3->next != room_head[!h][y][x])(p3->next)->prior = p3;
									room_head[!h][y][x]->next = p3;
									break;     //Bond break occurs
								}
							}
							if (j == 1) break;
						}
					}
					else if (p->type1 == 1)
					{
						flag5 = 0;  //Ma2: a flag labeling whether the circle chain breaking has happened
						for (j = p->length1; j > 0; j--)
						{
							f = PBB; flag1 = 0; //Ma2
							for (c2f1 = p->chain2->prior; c2f1 != p->chain2; c2f1 = c2f1->prior)
							{
								if (j > c2f1->start + 1)
								{
									if (j <= c2f1->length + c2f1->start)         // Falling into double chain region
									{
										f = f * sqrt(f); flag1 = 1;  //Ma2
									}
									break;
								}
								//else if (c2f1 == p->chain2->next)break;  //Ma2
							}
							if (j == 1 && p->type2 == 1)f = f * sqrt(f);

							for (flag = 0, a = 0; a < hrlength; a++)                    //new
							{
								if (j - 1 - hrlength / 2 + a >= p->length1)
								{
									if (p->information[0][j - 1 - hrlength / 2 + a - p->length1] == hrseq[a] && p->information[1][j - 1 - hrlength / 2 + a - p->length1] == 0)continue;
									else { flag = 1; break; }
								}
								else if (j - 1 - hrlength / 2 + a < 0)
								{
									if (p->information[0][p->length1 + j - 1 - hrlength / 2 + a] == hrseq[a] && p->information[1][p->length1 + j - 1 - hrlength / 2 + a] == 0)continue;
									else { flag = 1; break; }
								}
								else
								{
									if (p->information[0][j - 1 - hrlength / 2 + a] == hrseq[a] && p->information[1][j - 1 - hrlength / 2 + a] == 0)continue;
									else { flag = 1; break; }
								}
							}
							if (flag == 0)f = f * FHR;   //Ma2

							if (randd() < f)
							{
								flag5 = j;  // Ma2: the location of the circle chain breaking 
								p->type1 = 0;
								p->type2 = 0;
								if (j == 1)break;

								for (b = 0; b < j - 1; b++)
								{
									temp_information[0][b] = p->information[0][b];
									temp_information[1][b] = p->information[1][b];
								}
								for (b = 0; b < p->length1; b++)
								{
									if (b < p->length1 - j + 1)
									{
										p->information[0][b] = p->information[0][b + j - 1];
										p->information[1][b] = p->information[1][b + j - 1];
									}
									else
									{
										p->information[0][b] = temp_information[0][b - p->length1 + j - 1];
										p->information[1][b] = temp_information[1][b - p->length1 + j - 1];
									}
								}
								if (p->chain2->next == p->chain2)break;  //Ma2

								if (flag1 == 1) //f == PBB * sqrt(PBB)) //Ma2
								{
									c2f3 = (struct c2_frag*)malloc(sizeof(struct c2_frag));
									if (!c2f3) { printf("\tinit1--memeout\n"); exit(0); }
									c2f3->prior = c2f1;
									c2f3->next = c2f1->next;
									c2f1->next->prior = c2f3;
									c2f1->next = c2f3;

									c2f3->start = j - 1;
									c2f3->length = c2f1->start + c2f1->length - j + 1;
									c2f1->length = c2f1->length - c2f3->length;
								}
								else
								{
									c2f3 = c2f1->next;   //Ma2
								}

								//Ma2: two "if" sentences were deleted

								for (c2f2 = c2f3; c2f2 != p->chain2; c2f2 = c2f2->next) c2f2->start = c2f2->start - j + 1;
								for (c2f2 = p->chain2->next; c2f2 != c2f1->next; c2f2 = c2f2->next) c2f2->start = c2f2->start + p->length1 - j + 1;

								//if (p->chain2->next->next == p->chain2)break;  //Ma2

								if (c2f3 != p->chain2 && c2f1 != p->chain2)   //Ma2
								{
									p->chain2->prior->next = p->chain2->next;
									p->chain2->next->prior = p->chain2->prior;
									p->chain2->next = c2f3;
									p->chain2->prior = c2f1;
									c2f1->next = p->chain2;
									c2f3->prior = p->chain2;
								}

								if (p->length1 < p->length2) { printf("unexpected error on p-chain-length5"); exit(0); }
								//if (c2f1->start + c2f1->length > p->length1) { printf("unexpected error on p-chain-length6"); exit(0); }  //Ma2
								break;
							}
						}
						if (flag5 != 0) //Ma2--start: Circlar chain breaking has happened
						{
							flag4 = p->length1 - flag5 + 1;
							while (1)
							{
								for (j = p->length1; j > flag4; j--)  //Ma2: Only consider those bonds that have not been considered before the circular chain breaking
								{
									f = PBB;
									for (c2f1 = p->chain2->prior; c2f1 != p->chain2; c2f1 = c2f1->prior)
									{
										if (j > c2f1->start + 1)
										{
											if (j <= c2f1->length + c2f1->start)         // Falling into double chain region
											{
												f = f * sqrt(f);
											}
											break;
										}
									}

									if (j - 1 >= hrlength / 2 && p->length1 - j + 1 >= hrlength / 2)       //new
									{
										flag = 0;
										for (a = 0; a < hrlength; a++)
										{
											if (p->information[0][j - 1 - hrlength / 2 + a] == hrseq[a] && p->information[1][j - 1 - hrlength / 2 + a] == 0)continue;
											else { flag = 1; break; }
										}
									}
									else flag = 1;
									if (flag == 0)f = f * FHR;   //Ma2

									c2f1 = p->chain2->prior;
									c2f2 = p->chain2->next;
									if (randd() < f)
									{
										p3 = (struct rna*)malloc(LEN);
										if (!p3) { printf("\t%ddeg--memeout\n", i); exit(0); }
										memset(p3->information, 0, sizeof(p3->information));
										c2f = (struct c2_frag*)malloc(sizeof(struct c2_frag));
										if (!c2f) { printf("\tinit1--memeout\n"); exit(0); }
										p3->chain2 = c2f;
										c2f->next = p3->chain2;
										c2f->prior = p3->chain2;

										for (b = 0; b < p->length1 - j + 1; b++)
										{
											p3->information[0][b] = p->information[0][b + j - 1];
											p->information[0][b + j - 1] = 0;
										}
										p3->length1 = p->length1 - j + 1;
										p->length1 = j - 1;
										p3->length2 = 0;
										p3->type1 = 0;
										p3->type2 = 0;  //Ma2

										//if (p->length2 == 0 || c2f1->start + c2f1->length <= j - 1)    //Ma2
										//{
										//	p3->information[1][0] = 0;
										//}
										if (c2f1 != p->chain2 && c2f1->start + c2f1->length > j - 1)   //Ma2
										{
											//p3->type2 = 0;  //Ma2

											for (b = 0; b < c2f1->start + c2f1->length - j + 1; b++)
											{
												p3->information[1][b] = p->information[1][b + j - 1];
												p->information[1][b + j - 1] = 0;
											}

											while (c2f1 != p->chain2)
											{
												c2f3 = (struct c2_frag*)malloc(sizeof(struct c2_frag));
												if (!c2f3) { printf("\tinit1--memeout\n"); exit(0); }
												c2f3->prior = p3->chain2;
												c2f3->next = p3->chain2->next;
												p3->chain2->next->prior = c2f3;
												p3->chain2->next = c2f3;

												if (j <= c2f1->start + 1)
												{
													c2f3->length = c2f1->length;
													c2f3->start = c2f1->start - j + 1;
													p3->length2 = p3->length2 + c2f3->length;
													p->length2 = p->length2 - c2f3->length;

													(c2f1->prior)->next = c2f1->next;
													(c2f1->next)->prior = c2f1->prior;
													c2f = c2f1;
													c2f1 = c2f1->prior;
													free(c2f);

													if (p3->length1 < p3->length2) { printf("unexpected error on p3-chain-length"); exit(0); }
													if (c2f3->start + c2f3->length > p3->length1) { printf("unexpected error on p3-chain-length"); exit(0); }

													if (j >= c2f1->start + c2f1->length + 1)
													{
														break;
													}
												}
												else
												{
													c2f3->length = c2f1->start + c2f1->length - j + 1;
													c2f3->start = 0;
													c2f1->length = c2f1->length - c2f3->length;

													p3->length2 = p3->length2 + c2f3->length;
													p->length2 = p->length2 - c2f3->length;

													if (p->length1 < p->length2) { printf("unexpected error on p-chain-length"); exit(0); }
													if (p3->length1 < p3->length2) { printf("unexpected error on p3-chain-length"); exit(0); }
													if (c2f3->start + c2f3->length > p3->length1) { printf("unexpected error on p-chain-length"); exit(0); }
													break;
												}
											}
										}
										p3->prior = room_head[!h][y][x];
										p3->next = room_head[!h][y][x]->next;
										if (p3->next != room_head[!h][y][x])(p3->next)->prior = p3;
										room_head[!h][y][x]->next = p3;
										break;    //Bond break occurs
									}
								}
								if (j == flag4) break;
							}
						}//Ma2--end
					}
				}
				fresh_unit();
				break;

			case 2:                         //Template-directed addition
				for (c2f2 = p->chain2->next; c2f2 != p->chain2; c2f2 = c2f2->next)         //Template-directed ligation
				{
					rtdaddlig = randd();
					c2f3 = c2f2->next;
					if (c2f3 != p->chain2)
					{
						if (c2f2->length + c2f2->start == c2f3->start && rtdaddlig < PLT)
						{
							//if (p->information[1][c2f3->start] == 0 || p->information[1][c2f3->start - 1] == 0) { printf("unexpected error on add-nucleotide"); exit(0); }
							c2f2->length = c2f2->length + c2f3->length;
							c2f2->next = c2f3->next;
							(c2f3->next)->prior = c2f2;
							free(c2f3);
							break;
						}
					}
					else
					{
						if (p->type1 == 1 && c2f2->length + c2f2->start == p->length1 && p->chain2->next->start == 0 && rtdaddlig < PLT)
						{
							//if (p->information[1][0] == 0 && p->information[1][p->length1 - 1] == 0) { printf("unexpected error on add-nucleotide2"); exit(0); }
							c2f1 = p->chain2->prior;
							c2f4 = p->chain2->next;   //Ma2
							if (c2f1 == c2f4)
							{
								p->type2 = 1;
								break;
							}

							for (b = 0; b < c2f1->length; b++)
							{
								temp_information[0][b] = p->information[0][b + p->length1 - c2f1->length];
								temp_information[1][b] = p->information[1][b + p->length1 - c2f1->length];
							}
							for (b = p->length1 - 1; b >= 0; b--)
							{
								if (b >= c2f1->length)
								{
									p->information[0][b] = p->information[0][b - c2f1->length];
									p->information[1][b] = p->information[1][b - c2f1->length];
								}
								else
								{
									p->information[0][b] = temp_information[0][b];
									p->information[1][b] = temp_information[1][b];
								}
							}

							c2f4->length = c2f4->length + c2f1->length;   //Ma2
							for (c2f4 = c2f4->next; c2f4 != p->chain2; c2f4 = c2f4->next)  //Ma2
							{
								c2f4->start = c2f4->start + c2f1->length;
							}
							c2f1->prior->next = p->chain2;
							p->chain2->prior = c2f1->prior;
							free(c2f1);
							break;
						}
					}
				}

				//if (p->type1 == 0) { fresh_unit(); break; }  //Ma2-1  Linear RNA cannot attract substrates
				//Template-directed attraction of substrates	
				if (p->type1 == 1) f1=PAT;  //Ma2-1
				else f1 = PAT * FLT;   //Ma2-1
				for (p3 = p->next; p3 != p; p3 = p3->next)
				{
					if (p3 == room_head[h][y][x])
					{
						p3 = room_head[h][y][x]->next;
						if (p3 == p)break;
					}
					if (p3->length2 == 0 && p3->length1 <= p->length1)
					{
						c2f2 = p->chain2->next;
						if (p->type1 == 0 && (p->length2 == 0 || c2f2->start != 0))
						{
							if (p->length2 == 0)
							{
								length = p->length1;
							}
							else if (c2f2->start != 0)
							{
								length = c2f2->start;
							}

							if (p3->length1 <= length)
							{
								for (c = 0; c <= length - p3->length1; c++)
								{
									for (flag = 0, b = 0; b < p3->length1; b++)
									{
										if ((p3->information[0][p3->length1 - 1 - b] + p->information[0][c + b]) == 5)continue;
										else if (randd() < PFP)continue;
										else { flag = 1; break; }
									}
									if (flag == 0)break;
								}
								if (flag == 0)
								{
									rtdaddphili = randd() * FDA;  //Ma2
									if (rtdaddphili < f1)
									{
										for (a = 0; a < p3->length1; a++)
										{
											p->information[1][c + a] = p3->information[0][p3->length1 - 1 - a];
											if (p->information[1][c + a] == 0) { printf("unexpected error on add-nucleotide3"); exit(0); };
										}

										c2f3 = (struct c2_frag*)malloc(sizeof(struct c2_frag));
										if (!c2f3) { printf("\tinit1--memeout\n"); exit(0); }
										c2f3->prior = p->chain2;
										c2f3->next = p->chain2->next;
										p->chain2->next->prior = c2f3;
										p->chain2->next = c2f3;

										c2f3->start = c;
										c2f3->length = p3->length1;

										p->length2 = p->length2 + p3->length1;
										//p->type2 = 0;    //Ma2

										(p3->prior)->next = p3->next;
										if (p3->next != room_head[h][y][x])(p3->next)->prior = p3->prior;
										free(p3->chain2);
										free(p3);
										break;
									}
								}
							}
						}
						else if (p->type1 == 1 && p->length2 == 0)
						{
							if (p->length1 >= p3->length1)
							{
								for (c = 0; c <= p->length1 - 1; c++)
								{
									for (flag = 0, b = 0; b < p3->length1; b++)
									{
										if (c + b < p->length1)
										{
											if ((p3->information[0][p3->length1 - 1 - b] + p->information[0][c + b]) == 5)continue;
											else if (randd() < PFP)continue;
											else { flag = 1; break; }
										}
										else
										{
											if ((p3->information[0][p3->length1 - 1 - b] + p->information[0][c + b - p->length1]) == 5)continue;
											else if (randd() < PFP)continue;
											else { flag = 1; break; }
										}
									}
									if (flag == 0)break;
								}
								if (flag == 0)
								{
									rtdaddphili = randd() * FDA;   //Ma2
									if (rtdaddphili < f1)
									{
										c2f3 = (struct c2_frag*)malloc(sizeof(struct c2_frag));
										if (!c2f3) { printf("\tinit1--memeout\n"); exit(0); }
										c2f3->prior = p->chain2;
										c2f3->next = p->chain2->next;
										p->chain2->next->prior = c2f3;
										p->chain2->next = c2f3;
										c2f3->length = p3->length1;

										if (p3->length1 + c <= p->length1)
										{
											for (a = 0; a < p3->length1; a++) {
												p->information[1][c + a] = p3->information[0][p3->length1 - 1 - a];
												if (p->information[1][c + a] == 0) { printf("unexpected error on add-nucleotide41"); exit(0); };
											}

											c2f3->start = c;
										}
										else
										{
											for (b = 0; b < p->length1 - c; b++)
												temp_information[0][b] = p->information[0][b + c];

											for (b = p->length1 - 1; b >= 0; b--)
											{
												if (b >= p->length1 - c) {
													p->information[0][b] = p->information[0][b - p->length1 + c];
													if (p->information[0][b] == 0) { printf("unexpected error on add-nucleotide42"); exit(0); };
												}
												else {
													p->information[0][b] = temp_information[0][b];
													if (p->information[0][b] == 0) { printf("unexpected error on add-nucleotide43"); exit(0); };
												}
											}
											for (a = 0; a < p3->length1; a++) { p->information[1][a] = p3->information[0][p3->length1 - 1 - a]; }

											c2f3->start = 0;
										}
										p->length2 = p->length2 + p3->length1;
										//p->type2 = 0;

										(p3->prior)->next = p3->next;
										if (p3->next != room_head[h][y][x])(p3->next)->prior = p3->prior;
										free(p3->chain2);
										free(p3);
										break;
									}
								}
							}
						}
						else if (p->type1 == 1 && c2f2->start != 0)
						{
							c2f1 = p->chain2->prior;
							length = c2f2->start + p->length1 - c2f1->start - c2f1->length;

							if (length >= p3->length1)
							{
								for (c = 0; c <= length - p3->length1; c++)
								{
									for (flag = 0, b = 0; b < p3->length1; b++)
									{
										if (c + c2f1->start + c2f1->length + b < p->length1)   //Ma2
										{
											if ((p3->information[0][p3->length1 - 1 - b] + p->information[0][c + c2f1->start + c2f1->length + b]) == 5)continue;  //Ma2
											else if (randd() < PFP)continue;
											else { flag = 1; break; }
										}
										else
										{
											if ((p3->information[0][p3->length1 - 1 - b] + p->information[0][c + c2f1->start + c2f1->length + b - p->length1]) == 5)continue;  //Ma2
											else if (randd() < PFP)continue;
											else { flag = 1; break; }
										}
									}
									if (flag == 0)break;
								}
								if (flag == 0)
								{
									rtdaddphili = randd();
									if (c != 0)rtdaddphili = rtdaddphili * FDA;  //Ma2
									if (rtdaddphili < f1)
									{
										c2f3 = (struct c2_frag*)malloc(sizeof(struct c2_frag));
										if (!c2f3) { printf("\tinit1--memeout\n"); exit(0); }
										c2f3->length = p3->length1;

										if (p3->length1 + c <= p->length1 - c2f1->start - c2f1->length)
										{
											for (a = 0; a < p3->length1; a++)
												p->information[1][c2f1->start + c2f1->length + c + a] = p3->information[0][p3->length1 - 1 - a];

											c2f3->next = p->chain2;
											c2f3->prior = p->chain2->prior;
											p->chain2->prior->next = c2f3;
											p->chain2->prior = c2f3;
											c2f3->start = c + c2f1->start + c2f1->length;
										}
										else if (c >= p->length1 - c2f1->start - c2f1->length)
										{
											for (a = 0; a < p3->length1; a++)
												p->information[1][c - p->length1 + c2f1->start + c2f1->length + a] = p3->information[0][p3->length1 - 1 - a];

											c2f3->prior = p->chain2;
											c2f3->next = p->chain2->next;
											p->chain2->next->prior = c2f3;
											p->chain2->next = c2f3;
											c2f3->start = c - p->length1 + c2f1->start + c2f1->length;
										}
										else
										{
											for (b = 0; b < p->length1 - c2f1->start - c2f1->length - c; b++)
											{
												temp_information[0][b] = p->information[0][b + c2f1->start + c2f1->length + c];
												temp_information[1][b] = p->information[1][b + c2f1->start + c2f1->length + c];
											}
											for (b = p->length1 - 1; b >= 0; b--)
											{
												if (b >= p->length1 - c2f1->start - c2f1->length - c)
												{
													p->information[0][b] = p->information[0][b - p->length1 + c2f1->start + c2f1->length + c];
													p->information[1][b] = p->information[1][b - p->length1 + c2f1->start + c2f1->length + c];
													if (p->information[0][b] == 0) { printf("unexpected error on add-nucleotide51"); exit(0); }
												}
												else
												{
													p->information[0][b] = temp_information[0][b];
													p->information[1][b] = temp_information[1][b];
													if (p->information[0][b] == 0) { printf("unexpected error on add-nucleotide52"); exit(0); };
												}
											}
											for (a = 0; a < p3->length1; a++) { p->information[1][a] = p3->information[0][p3->length1 - 1 - a]; }

											for (c2f2 = p->chain2->next; c2f2 != p->chain2; c2f2 = c2f2->next)
												c2f2->start = c2f2->start + p->length1 - c2f1->start - c2f1->length - c;

											c2f3->prior = p->chain2;
											c2f3->next = p->chain2->next;
											p->chain2->next->prior = c2f3;
											p->chain2->next = c2f3;
											c2f3->start = 0;
										}
										p->length2 = p->length2 + p3->length1;

										(p3->prior)->next = p3->next;
										if (p3->next != room_head[h][y][x])(p3->next)->prior = p3->prior;
										free(p3->chain2);
										free(p3);
										break;
									}
								}
							}
						}
						flag3 = 0;  //Ma2:  a flag labeling whether the attraction has happened
						for (c2f2 = p->chain2->next; c2f2 != p->chain2; c2f2 = c2f2->next)
						{
							if (c2f2->next == p->chain2)
							{
								if (p->type1 == 1 && p->chain2->next->start != 0) break;  //Ma2:  "p->length2 == 0 ||" is deleted
								length = p->length1 - c2f2->start - c2f2->length;
							}
							else
							{
								length = c2f2->next->start - c2f2->start - c2f2->length;
							}

							if (p3->length1 <= length)
							{
								for (c = 0; c <= length - p3->length1; c++)
								{
									for (flag = 0, b = 0; b < p3->length1; b++)
									{
										if ((p3->information[0][p3->length1 - 1 - b] + p->information[0][c2f2->start + c2f2->length + c + b]) == 5)continue;
										else if (randd() < PFP)continue;
										else { flag = 1; break; }
									}
									if (flag == 0)break;
								}
								if (flag == 0)
								{
									rtdaddphili = randd();
									if (c != 0)rtdaddphili = rtdaddphili * FDA; //Ma2
									if (rtdaddphili < f1)
									{
										for (a = 0; a < p3->length1; a++) {
											p->information[1][c2f2->start + c2f2->length + c + a] = p3->information[0][p3->length1 - 1 - a];
											if (p->information[1][c2f2->start + c2f2->length + c + a] == 0) { printf("unexpected error on add-nucleotide6"); exit(0); };
										}

										c2f3 = (struct c2_frag*)malloc(sizeof(struct c2_frag));
										if (!c2f3) { printf("\tinit1--memeout\n"); exit(0); }
										c2f3->prior = c2f2;
										c2f3->next = c2f2->next;
										c2f2->next->prior = c2f3;
										c2f2->next = c2f3;

										c2f3->start = c2f2->start + c2f2->length + c;
										c2f3->length = p3->length1;

										p->length2 = p->length2 + p3->length1;

										(p3->prior)->next = p3->next;
										if (p3->next != room_head[h][y][x])(p3->next)->prior = p3->prior;
										free(p3->chain2);
										free(p3);
										flag3 = 1;   //Ma2:  the attraction has happened.  -- "flag=2" is changed.
										break;
									}
								}
							}
						}
						if (flag3 == 1) break;  //Ma2 "flag==2" is changed
						//{
						//	flag = 0;
						//	break;
						//}
					}
				}
				fresh_unit();
				break;

			case 3:                           // Separation
				if (p->length2 != 0)    // Separation of double chain  
				{
					j = 0;
					for (c2f2 = p->chain2->next; c2f2 != p->chain2; c2f2 = c2f2->next) j++;
					randseq = randl(j);
					c2f2 = p->chain2->next;
					for (j = 1; j <= randseq; j++)
					{
						c2f2 = c2f2->next;
					}

					if (randd() < pow(PSBP, sqrt(c2f2->length *1.0)))  //Ma2-1
					{
						p3 = (struct rna*)malloc(LEN);
						if (!p3) { printf("\t%dsep--memeout\n", i); exit(0); }
						memset(p3->information, 0, sizeof(p3->information));
						c2f = (struct c2_frag*)malloc(sizeof(struct c2_frag));
						if (!c2f) { printf("\tinit1--memeout\n"); exit(0); }
						p3->chain2 = c2f;
						c2f->next = p3->chain2;
						c2f->prior = p3->chain2;

						for (b = 0; b < c2f2->length; b++)
						{
							p3->information[0][b] = p->information[1][c2f2->start + c2f2->length - 1 - b];
							if (p3->information[0][b] == 0) { printf("unexpected error on seperation"); exit(0); }
							p->information[1][c2f2->start + c2f2->length - 1 - b] = 0;
						}
						p3->length1 = c2f2->length;
						p->length2 = p->length2 - c2f2->length;
						p3->length2 = 0;
						if (p->type2 == 1)
						{
							p3->type1 = 1;
							p->type2 = 0;    //Ma2
						}
						else p3->type1 = 0;
						p3->type2 = 0;      //Ma2

						(c2f2->prior)->next = c2f2->next;
						(c2f2->next)->prior = c2f2->prior;
						free(c2f2);

						p3->prior = room_head[!h][y][x];
						p3->next = room_head[!h][y][x]->next;
						if (p3->next != room_head[!h][y][x])(p3->next)->prior = p3;
						room_head[!h][y][x]->next = p3;
					}
				}
				fresh_unit();
				break;

				//-------------------------------------------------------------------------------------
			case 4:    // nr catalyses the synthesis of nt
				if (p->type1 == 0 && p->length2 == 0 && p->length1 >= nrlength && p->length1 < 2 * nrlength)  //Ma2-1, p->length1 < 2 * nrlength:: Nr cannot be much longer than its characteristic domain
				{
					flag = findseq(nrseq, nrlength, p);
					if (flag == 0)     //Ma2-1
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
								if (!p3) { printf("\t%dnr form_monomer--memeout\n", k + 1); exit(0); }
								memset(p3->information, 0, sizeof(p3->information));
								randnt = randl(4) + 1;
								switch (randnt)
								{
								case 1:  p3->information[0][0] = A; break;
								case 2:  p3->information[0][0] = C; break;
								case 3:  p3->information[0][0] = G; break;
								case 4:  p3->information[0][0] = U; break;
								default: printf("nr randnt error");
								}
								//p3->information[0][1] = 0;  //Ma2
								//p3->information[1][0] = 0;  //Ma2

								p3->length1 = 1;
								p3->length2 = 0;
								p3->type1 = 0;
								p3->type2 = 0;

								p3->prior = room_head[!h][y][x];
								p3->next = room_head[!h][y][x]->next;
								if (p3->next != room_head[!h][y][x])(p3->next)->prior = p3;
								room_head[!h][y][x]->next = p3;

								c2f = (struct c2_frag*)malloc(sizeof(struct c2_frag));
								if (!c2f) { printf("\tinit1--memeout\n"); exit(0); }
								p3->chain2 = c2f;
								c2f->next = p3->chain2;
								c2f->prior = p3->chain2;
							}
						}
					}
				}
				fresh_unit();
				break;

			case 5:            //moving to another adjacent cell
				if (randd() * RMRW < PMV)    //Ma2: with toroidal topology to avoid edge effects
				{
					randcase1 = randl(4);   // Four possible directions
					switch (randcase1)
					{
					case 0:
						p1 = p->prior;
						p2 = p->next;

						p3 = room_head[!h][y][(N + x - 1) % N]->next;
						room_head[!h][y][(N + x - 1) % N]->next = p;
						p->next = p3;
						p->prior = room_head[!h][y][(N + x - 1) % N];
						if (p3 != room_head[!h][y][(N + x - 1) % N])p3->prior = p;

						p1->next = p2;
						if (p2 != room_head[h][y][x])p2->prior = p1;
						p = p1;
						break;

					case 1:
						p1 = p->prior;
						p2 = p->next;

						p3 = room_head[!h][y][(x + 1) % N]->next;
						room_head[!h][y][(x + 1) % N]->next = p;
						p->next = p3;
						p->prior = room_head[!h][y][(x + 1) % N];
						if (p3 != room_head[!h][y][(x + 1) % N])p3->prior = p;

						p1->next = p2;
						if (p2 != room_head[h][y][x])p2->prior = p1;
						p = p1;
						break;

					case 2:
						p1 = p->prior;
						p2 = p->next;

						p3 = room_head[!h][(N + y - 1) % N][x]->next;
						room_head[!h][(N + y - 1) % N][x]->next = p;
						p->next = p3;
						p->prior = room_head[!h][(N + y - 1) % N][x];
						if (p3 != room_head[!h][(N + y - 1) % N][x])p3->prior = p;

						p1->next = p2;
						if (p2 != room_head[h][y][x])p2->prior = p1;
						p = p1;
						break;

					case 3:
						p1 = p->prior;
						p2 = p->next;

						p3 = room_head[!h][(y + 1) % N][x]->next;
						room_head[!h][(y + 1) % N][x]->next = p;
						p->next = p3;
						p->prior = room_head[!h][(y + 1) % N][x];
						if (p3 != room_head[!h][(y + 1) % N][x])p3->prior = p;

						p1->next = p2;
						if (p2 != room_head[h][y][x])p2->prior = p1;
						p = p1;
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
void record(void)               // Data recording
{
	int ch_num[MAX_RNA_LENGTH], ch_nr_num[MAX_RNA_LENGTH], ch_hr_num[MAX_RNA_LENGTH], ch_hrnr_num[MAX_RNA_LENGTH], ch_0_num[MAX_RNA_LENGTH], long_chain_num, si;

	FILE* fptxt, * fptxt1;
	errno_t err, err1;
	err = fopen_s(&fptxt, "monitor.txt", "at");
	if (err != 0) { printf("cannot open file");  exit(-1); }
	err1 = fopen_s(&fptxt1, "picture.txt", "at");
	if (err1 != 0) { printf("cannot open file1");  exit(-1); }

	nr[g] = 0;
	cir_nr[g] = 0;

	hr[g] = 0;
	cir_hr[g] = 0;

	hr_nr[g] = 0;
	cir_hr_nr[g] = 0;

	ctr1[g] = 0;
	ctr2[g] = 0;
	ctr3[g] = 0;

	cir_ctr1[g] = 0;
	cir_ctr2[g] = 0;
	cir_ctr3[g] = 0;

	cir_unit[g] = 0;

	total_mat_num[g] = 0;
	unit[g] = 0;
	raw_num[g] = 0;
	for (y = 0; y < N; y++)
	{
		for (x = 0; x < N; x++)
		{
			raw_num[g] += raw_arr[y][x];
			for (p = room_head[h][y][x]->next; p != room_head[h][y][x]; p = p->next)
			{
				unit[g]++;
				total_mat_num[g] += p->length1 + p->length2;

				flag = findseq(nrseq, nrlength, p);
				flag4 = findseq(hrseq, hrlength, p);
				if (flag == 0) nr[g]++;
				if (flag4 == 0) hr[g]++;
				if (flag == 0 && flag4 == 0) hr_nr[g]++;

				flag1 = findseq(contrseq1, contrlength1, p);
				if (flag1 == 0) ctr1[g]++; //
				flag2 = findseq(contrseq2, contrlength2, p);
				if (flag2 == 0) ctr2[g]++; //
				flag3 = findseq(contrseq3, contrlength3, p);
				if (flag3 == 0) ctr3[g]++; //

				if (p->type1 == 1)
				{
					cir_unit[g]++;
					if (flag == 0) cir_nr[g]++;
					if (flag4 == 0) cir_hr[g]++;
					if (flag == 0 && flag4 == 0) cir_hr_nr[g]++;
					if (flag1 == 0) cir_ctr1[g]++;
					if (flag2 == 0) cir_ctr2[g]++;
					if (flag3 == 0) cir_ctr3[g]++;
				}
			}
		}
	}
	total_mat_num[g] += raw_num[g];

	printf("----- step=%d: nr=%d (%d), hr=%d (%d), hr_nr=%d (%d), ctr1=%d (%d), ctr2=%d (%d), ctr3=%d (%d), unit=%d (%d)  (tn=%d, r=%d)\n", i,
		(int)nr[g], (int)cir_nr[g], (int)hr[g], (int)cir_hr[g], (int)hr_nr[g], (int)cir_hr_nr[g],
		(int)ctr1[g], (int)cir_ctr1[g], (int)ctr2[g], (int)cir_ctr2[g], (int)ctr3[g], (int)cir_ctr3[g],
		(int)unit[g], (int)cir_unit[g], (int)total_mat_num[g], (int)raw_num[g]);
	fprintf(fptxt, "----- step=%d: nr=%d (%d), hr=%d (%d), hr_nr=%d (%d), ctr1=%d (%d), ctr2=%d (%d), ctr3=%d (%d), unit=%d (%d)  (tn=%d, r=%d)\n", i,
		(int)nr[g], (int)cir_nr[g], (int)hr[g], (int)cir_hr[g], (int)hr_nr[g], (int)cir_hr_nr[g],
		(int)ctr1[g], (int)cir_ctr1[g], (int)ctr2[g], (int)cir_ctr2[g], (int)ctr3[g], (int)cir_ctr3[g],
		(int)unit[g], (int)cir_unit[g], (int)total_mat_num[g], (int)raw_num[g]);
	fprintf(fptxt1, "step=%d: nr=%d (%d), hr=%d (%d), hr_nr=%d (%d), ctr1=%d (%d), ctr2=%d (%d), ctr3=%d (%d), unit=%d (%d)  (tn=%d, r=%d)\n", i,
		(int)nr[g], (int)cir_nr[g], (int)hr[g], (int)cir_hr[g], (int)hr_nr[g], (int)cir_hr_nr[g],
		(int)ctr1[g], (int)cir_ctr1[g], (int)ctr2[g], (int)cir_ctr2[g], (int)ctr3[g], (int)cir_ctr3[g],
		(int)unit[g], (int)cir_unit[g], (int)total_mat_num[g], (int)raw_num[g]);

	for (si = 0; si < LONG_CHAIN_LEN; si++)  //Ma2---start
	{
		ch_num[si] = 0;
		ch_nr_num[si] = 0;
		ch_hr_num[si] = 0;
		ch_hrnr_num[si] = 0;
		ch_0_num[si] = 0;
	}
	long_chain_num = 0;
	for (y = 0; y < N; y++)
	{
		for (x = 0; x < N; x++)
		{
			for (p = room_head[0][y][x]->next; p != room_head[0][y][x]; p = p->next)
			{
				ch_num[p->length1 - 1]++;
				if (p->length1 > LONG_CHAIN_LEN)
				{
					long_chain_num++;
					for (int t = 0; t < p->length1; t++)
					{
						switch (p->information[0][t])
						{
						case 1: printf("A"); fprintf(fptxt, "A"); break;
						case 2: printf("C"); fprintf(fptxt, "C"); break;
						case 3: printf("G"); fprintf(fptxt, "G"); break;
						case 4: printf("U"); fprintf(fptxt, "U"); break;
						default: printf("error");
						}
					}
					printf("\t%d\t%d\n", p->length1, p->type1);
					fprintf(fptxt, "\t%d\t%d\n", p->length1, p->type1);
					for (int t = 0; t < p->length1; t++)
					{
						switch (p->information[1][t])
						{
						case 0: printf("-"); fprintf(fptxt, "-"); break;
						case 1: printf("A"); fprintf(fptxt, "A"); break;
						case 2: printf("C"); fprintf(fptxt, "C"); break;
						case 3: printf("G"); fprintf(fptxt, "G"); break;
						case 4: printf("U"); fprintf(fptxt, "U"); break;
						default: printf("error");
						}
					}
					printf("\t\t%d\n", p->type2);
					fprintf(fptxt, "\t\t%d\n", p->type2);
				}
				flag1 = findseq(nrseq, nrlength, p);
				//flag3 = findseq(hrseq, hrlength, p);
				flag3 = findseq(contrseq1, contrlength1, p);
				if (flag1 == 0 && flag3 != 0)ch_nr_num[p->length1 - 1]++;
				else if (flag1 != 0 && flag3 == 0)ch_hr_num[p->length1 - 1]++;
				else if (flag1 == 0 && flag3 == 0)ch_hrnr_num[p->length1 - 1]++;
				else ch_0_num[p->length1 - 1]++;
			}
		}
	}

	for (si = 0; si < LONG_CHAIN_LEN; si++)
	{
		printf("%dnt-%d(%d|%d/%d^%d), ", si+1, ch_num[si], ch_hr_num[si], ch_nr_num[si], ch_hrnr_num[si], ch_0_num[si]);
		fprintf(fptxt, "%dnt-%d(%d|%d/%d^%d), ", si+1, ch_num[si], ch_hr_num[si], ch_nr_num[si], ch_hrnr_num[si], ch_0_num[si]);
	}
	printf("\nChains over %d = %d   step=%d\n\n", LONG_CHAIN_LEN, long_chain_num, i);
	fprintf(fptxt, "\nChains over %d = %d   step=%d\n\n", LONG_CHAIN_LEN, long_chain_num, i);  //Ma2---end

	g++;

	fclose(fptxt);
	fclose(fptxt1);
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
					if (room_head[m][y][x]->next != room_head[m][y][x])
					{
						p = room_head[m][y][x]->next;
						room_head[m][y][x]->next = p->next;
						while (1)
						{
							if (p->chain2->next != p->chain2)
							{
								c2f = p->chain2->next;
								p->chain2->next = c2f->next;
								free(c2f);
							}
							else break;
						}
						free(p->chain2);
						free(p);
					}
					else break;
				}
				free(room_head[m][y][x]);
			}
		}
	}
}

//------------------------------------------------------------------------------ 
int main()
{
	inits();        // initialization of the system

	for (i = 0; i <= STEPNUM; i++)      // Monte-Carlo cycle
	{
		if (i == INOCUSTEP)inoculate(0);  // Ma2-1
		if (i >= STAREC && i % RECINT == 0)
		{
			record();
		}
		unit_case();
		h = !h;
	}

	freepool();

	return (0);
}
//========================================================  End of the program






