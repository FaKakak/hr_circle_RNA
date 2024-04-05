//Monte Carlo simulation of emergence and subsequent evolution of RNA replicases in a nucleotide pool
//-----C source codes for the simulation program (replicase.c)
//record computing time
//erase bug of "raw--"
//no xy choose
//write fresh_unit into unit_case
//write find_seq into unit_case

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <iostream>
using namespace std;

/******* RANDOM NUMBER GENERATOR BY ZIFF **********/
#define A1 471
#define B1 1586
#define C1 6988
#define D1 9689
#define M 16383
#define RIMAX 2147483648.0        // = 2^31 
#define RandomInteger (++nd, ra[nd & M] = ra[(nd-A1) & M] ^ ra[(nd-B1) & M] ^ ra[(nd-C1) & M] ^ ra[(nd-D1) & M])
void seed(long seed);  // random number initialization 
static long ra[M+1], nd;
/**************************************************/

#define RECDATA1 "rep-1.dat"
#define RECDATA2 "rep-2.dat"

#define LEN sizeof(struct rna)
#define C 2
#define G 3
#define A 1
#define U 4
#define STEPNUM 300000     // Total time steps of Monte Carlo simulation
#define STAREC 0          // The step to start record
#define RECINT 30000       // The interval steps of recording
#define MAX_RNA_LENGTH 100    // Defining maximum RNA length allowed in the simulation

#define SIDE 3                // The side length of the two-dimensional grid
#define TOTAL_MATERIAL 3000  // Total materials in the system
#define REPSEQ G,A,C,G,U,C   // The presumed specific sequence with which a polynucleotide could act as a replicase
#define REPCOMSEQ G,A,C,G,U,C // The complementary sequence of the presumed replicase sequence
const double PSBP=0.9;              // Probability of separation of a base-pair
const double PBB=0.0001    ;         // Probability of breaking of a phosphodiester bond
const double PLMC= 0.0002 ;          // Probability of ligation of two unit in a cell with mineral catalysis 
const double PAT= 0.01   ;            // Probability of attracting a substrate by a template when the substrate could base-pair with the template
const double PATR= 0.9  ;            // Probability of attracting a substrate by a template in the existing of a binding replicase
const double PLT= 0.005;             // Probability of a template-directed ligation 
const double PLTR= 0.9      ;        // probability of a template-directed ligation in the existing of a binding replicase
const double PRB= 0.95     ;         // probability of binding of a replicase onto a template
const double PRD= 0.05    ;          // probability of dropping of replicase from template while the copying of complementary chain has not completed
const double PMOV= 0.9   ;          // Probability of movement of a mononucleotide
const double PMF= 0.0001;            // probability of mononucleotide formation from raw materials
const double PMD= 0.001;             // probability of mononucleotide decay to raw materials
const double PFP= 0.01;             // probability of false base-pairing
#define REPCOVERSEQ (sqrt(replength)+2)  // The cover sequence length of a replicase binding on the template
#define FDMOV (pow(p->length1+p->length2+length3,1/3.0))  
                             // The factor defining the relationship between probability of moving and molecular weight
//#define CELLNUM (SIDE*SIDE)   // Cell numbers in the grid

long randl(long);      // random number between 0 and parameter 
double randd(void);    // random double between 0 and 1         
void inits(void);         // initialization of the system
void form_monomer(void);  // Forming of mononucleotide from raw materials
void unit_case(void);     // Action of units (molecules) in the system
int record(void);         // Data recording
void save_result(void);   // Data saving 
void freepool(void);      // Memory releasing

struct rna                // A unit of mononucleotide or polynucleotide
{char information[2][MAX_RNA_LENGTH];
 int length1;
 int length2;
 int nick;
 struct rna *chain3;
 struct rna *next;
 struct rna *prior;
 };
struct rna *cell_head[2][SIDE][SIDE];
struct rna *p,*p1,*p2,*p3,*p4;

static char repseq[50]={REPSEQ};        // Presumed replicase sequence
static char repcomseq[50]={REPCOMSEQ};
int raw=TOTAL_MATERIAL;  // Initializating the raw materials
int over_max_len=0;
int x,y;                 // The coordinate of cells in the grid 
int replength,randcase,length3,g=0,h=0,g_end=0,gi,g_stop=0;
int flag,flag1,flag2,flag3,flagn,flagn1,flagn2,flagn3,flag4=0,flag5;
long i;                  // Cycle variable for Monte Carlo steps
//long available;
//long availabl[CELLNUM];
long recstep[(STEPNUM-STAREC)/RECINT+1];    // Record steps
float rep_num[(STEPNUM-STAREC)/RECINT+1];  // Record number of replicases in steps
  // Number of replicases and polynucleotides including the complementary sequence of presumed replicase sequence 
time_t time_start, time_end, time1, time2, time3, test;
time_t *ptimer=&test;

/***********************************************************
 * Random generator initialization                         *
 *                      by a simple Congruential generator *
 ***********************************************************/
void seed(long seed)      
{
  int a;
 
  if(seed<0) { puts("SEED error."); exit(1); }
  ra[0]= (long) fmod(16807.0*(double)seed, 2147483647.0);
  for(a=1; a<=M; a++)
  {
    ra[a] = (long)fmod( 16807.0 * (double) ra[a-1], 2147483647.0);
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
   return((double) RandomInteger / RIMAX);
  }

//------------------------------------------------------------------------------
void inits(void)         // initialization of the system
{
	int j,m;
	seed(555);       

	replength=0;
	for(j=0;repseq[j]!=0;j++)
	replength++;

	for(m=0;m<2;m++)
	{
		for(y=0;y<SIDE;y++)
		{
			for(x=0;x<SIDE;x++)
			{
				p1=(struct rna *)malloc(LEN);          
				if(!p1){printf("\tinit1--memeout\n"); exit(0);}
				cell_head[m][y][x]=p1;
				p1->next=cell_head[m][y][x];
			}
		}
	}
}

//------------------------------------------------------------------------------
void form_monomer(void)   // Forming of mononucleotide from raw materials
{
 int k,randnt,raw_bef=raw;

 for(k=0;k<raw_bef;k++)    
 {
	 if(randd()<PMF)
	 {
		 raw--;
		 x=randl(SIDE);
		 y=randl(SIDE);
 
//		 p1=cell_head[h][y][x];
//		 for(p=cell_head[h][y][x]->next; p!=cell_head[h][y][x]; p=p->next)   
//			 p1=p;

		 p2=(struct rna *)malloc(LEN);    
		  if(!p2){printf("\t%dform_monomer--memeout\n",k+1); exit(0);}
         randnt=randl(4)+1;
		   switch(randnt)
		   {
			   case 1:  p2->information[0][0]=A; break;
			   case 2:  p2->information[0][0]=C; break;
			   case 3:  p2->information[0][0]=G; break;
			   case 4:  p2->information[0][0]=U; break;
			   default: printf("randnt error");
		   }
          p2->information[0][1]=0;
		  p2->information[1][0]=0;

		  p2->length1=1;
		  p2->length2=0;
		  p2->nick=0;
		  p2->chain3=NULL;
		  p2->next=cell_head[!h][y][x]->next;
		  if(p2->next!=cell_head[!h][y][x])(p2->next)->prior=p2;
		  cell_head[!h][y][x]->next=p2;
		  p2->prior=cell_head[!h][y][x];
//		  p1->next=p2;
	 }
 }
}

//------------------------------------------------------------------------------
void unit_case(void)      // Action of units (molecules) in the system
{ 
	int a,b,c,j,k,m,n;
	double f,rtdaddlig, rtdaddphili;
	
	for(y=0;y<SIDE;y++)
		{
			for(x=0;x<SIDE;x++)
			{

			for(p=cell_head[h][y][x]->next; p!=cell_head[h][y][x]; p=p->next)
		  {
			  randcase=randl(6);    
		   switch(randcase)
		   {case 0:                        // Chain ligation with mineral catalysis
		    for(p3=p->next;p3!=cell_head[h][y][x];p3=p3->next)   
			{
              if(p3->length2==0&&p3->chain3==NULL)
				  {if(randd()<PLMC/p3->length1)
				   {
						if(p->length1+p3->length1>MAX_RNA_LENGTH-1)          
						{   
							over_max_len++;continue;
						}
						
						for(a=0;a<p3->length1;a++)
						  p->information[0][a+p->length1]=p3->information[0][a];
						p->information[0][p->length1+p3->length1]=0;
						p->length1=p->length1+p3->length1;

						(p3->prior)->next=p3->next;
						if(p3->next!=cell_head[h][y][x])(p3->next)->prior=p3->prior;
						free(p3);
			
						break;
				   }
				  }
			}
			p1=p->prior;                    //fresh_unit      
			p2=p->next;
			p3=cell_head[!h][y][x]->next;           
			cell_head[!h][y][x]->next=p;
			p->next=p3;
			p->prior=cell_head[!h][y][x];
			if(p3!=cell_head[!h][y][x])p3->prior=p;
			p1->next=p2;                              
			if(p2!=cell_head[h][y][x])p2->prior=p1;
			p=p1;                                       
			break;

			case 1:            // Decay and degradation 
			if(p->length1==1)  // Decay of mononucleotide
			{
				if(p->length2==0&&p->chain3==NULL&&randd()<PMD)
				{
					raw++;
					(p->prior)->next=p->next;
					if(p->next!=cell_head[h][y][x])(p->next)->prior=p->prior;
					p3=p;
					p=p->prior;
					free(p3); break;
				}
			}
			else{                 //Degradation of chain
				if(p->chain3!=NULL){
					if(p->length2==0) c=1;             
					else if(p->nick==0) c=p->length2;
					else c=p->nick;
				}

				f=PBB;
				for(j=p->length1;j>1;j--){
					if(j<=p->length2){   // Falling into double chain region
						if(p->nick==0){
							m=j-1;
							n=p->length2-j+1;
							k=(m<n)?m:n;
							f=PBB*k*PBB*k;
						}
						else{
							if(j==p->nick+1){
								m=j-1;
								n=p->length2-j+1;
								k=(m<n)?m:n;
								f=PBB*k;
							}
							else if(j>p->nick+1){
								m=j-p->nick-1;
								n=p->length2-j+1;
								k=(m<n)?m:n;
								f=PBB*k*PBB*k;
							}
							else{
								m=j-1;
								n=p->nick-j+1;
								k=(m<n)?m:n;
								f=PBB*k*PBB*k;
							}
						}
					}


					if((p->chain3!=NULL&&(j<=c||j>=c+REPCOVERSEQ)||p->chain3==NULL)&&randd()<f){									
       				p3=(struct rna *)malloc(LEN);
						if(!p3){printf("\t%ddeg--memeout\n",i); exit(0);}

						for(b=0;b<p->length1-j+1;b++)
							p3->information[0][b]=p->information[0][b+j-1];
						p3->information[0][p->length1-j+1]=0;
						p->information[0][j-1]=0;
						p3->length1=p->length1-j+1;
						p->length1=j-1;

						if(p->length2>j-1){
							for(b=0;b<p->length2-j+1;b++)
								p3->information[1][b]=p->information[1][b+j-1];
 							p3->information[1][p->length2-j+1]=0;
							p->information[1][j-1]=0;
							p3->length2=p->length2-j+1; 
							p->length2=j-1;
						}
						else{
							p3->information[1][0]=0;
							p3->length2=0;
						}
 						  
						if(p->nick>j-1) {p3->nick=p->nick-j+1; p->nick=0;}
						else if(p->nick==j-1) {p3->nick=0; p->nick=0;}
  						else p3->nick=0;

						if(p->chain3!=NULL&&j<=c){
							p3->chain3=p->chain3;
							p->chain3=NULL;
						}
						else p3->chain3=NULL;
						
						p3->prior=cell_head[!h][y][x];    
						p3->next=cell_head[!h][y][x]->next;
						if(p3->next!=cell_head[!h][y][x])(p3->next)->prior=p3;
						cell_head[!h][y][x]->next=p3;
						break;
					}
				}
			}
	p1=p->prior;                          
	p2=p->next;
	p3=cell_head[!h][y][x]->next;           
	cell_head[!h][y][x]->next=p;
	p->next=p3;
	p->prior=cell_head[!h][y][x];
	if(p3!=cell_head[!h][y][x])p3->prior=p;
	p1->next=p2;                              
	if(p2!=cell_head[h][y][x])p2->prior=p1;
	p=p1;                                       
   			break;

			case 2:                         //Template-directed addition
			if(p->nick==0){                      //Template-directed attraction of substrates
				for(p3=p->next;p3!=cell_head[h][y][x];p3=p3->next){
					if(p3->length2==0&&p3->chain3==NULL){
						if(p3->length1<=p->length1-p->length2){
							for(flag=0,b=0;b<p3->length1;b++){
								if((p3->information[0][p3->length1-1-b]+p->information[0][p->length2+b])!=5){
									if(randd()>PFP){flag=1;	break;}
								}
							}
							if(flag==0){
								rtdaddphili=randd();
								if(p->chain3==NULL&&rtdaddphili<PAT||p->chain3!=NULL&&rtdaddphili<PATR){
									for(a=0;a<p3->length1;a++)
										p->information[1][p->length2+a]=p3->information[0][p3->length1-1-a];
									p->information[1][p->length2+p3->length1]=0;
								   if(p->length2!=0)p->nick=p->length2;
								   p->length2=p->length2+p3->length1;
								   
								   (p3->prior)->next=p3->next;
								   if(p3->next!=cell_head[h][y][x])(p3->next)->prior=p3->prior;
								   free(p3);
								   break;
								}
							}
						}
					}
				}
			}
			else{                    //Template-directed ligation
				rtdaddlig=randd();
			   if(p->chain3==NULL&&rtdaddlig<PLT||p->chain3!=NULL&&rtdaddlig<PLTR)
					p->nick=0;
			}
			p1=p->prior;                          
			p2=p->next;
			p3=cell_head[!h][y][x]->next;           
			cell_head[!h][y][x]->next=p;
			p->next=p3;
			p->prior=cell_head[!h][y][x];
			if(p3!=cell_head[!h][y][x])p3->prior=p;
			p1->next=p2;                              
			if(p2!=cell_head[h][y][x])p2->prior=p1;
			p=p1;                                       
			break;

			case 3:                           // Separation
			if(p->chain3!=NULL)                // Separation of replicase and template
			{if(randd()<PRD||(p->length2==p->length1&&p->nick==0))  
			 {  
				p3=p->chain3;           
				p->chain3=NULL;
   
				p3->prior=cell_head[!h][y][x];    
				p3->next=cell_head[!h][y][x]->next;
				if(p3->next!=cell_head[!h][y][x])(p3->next)->prior=p3;
				cell_head[!h][y][x]->next=p3;
			 }
			}
			else                  // Separation of double chain   
			{                 
				if(p->length2!=0)
				{
					if(randd()<pow(PSBP,p->length2-p->nick))            
					{
						p3=(struct rna *)malloc(LEN);
						if(!p3){printf("\t%dsep--memeout\n",i); exit(0);}
						for(b=0;b<p->length2-p->nick;b++)
						   p3->information[0][b]=p->information[1][p->length2-1-b];
						p->information[1][p->nick]=0;
						  
					   p3->information[0][p->length2-p->nick]=0;
						p3->information[1][0]=0;

						p3->length1=p->length2-p->nick;
						p->length2=p->nick;
						p->nick=0;
						p3->length2=0;
						p3->nick=0;
						p3->chain3=NULL;

						p3->prior=cell_head[!h][y][x];    
						p3->next=cell_head[!h][y][x]->next;
						if(p3->next!=cell_head[!h][y][x])(p3->next)->prior=p3;
						cell_head[!h][y][x]->next=p3;
					}
				}
			}
	p1=p->prior;                          
	p2=p->next;
	p3=cell_head[!h][y][x]->next;           
	cell_head[!h][y][x]->next=p;
	p->next=p3;
	p->prior=cell_head[!h][y][x];
	if(p3!=cell_head[!h][y][x])p3->prior=p;
	p1->next=p2;                              
	if(p2!=cell_head[h][y][x])p2->prior=p1;
	p=p1;                                       
			break;

			case 4:                            // Binding of replicase to template
			if(p->length2==0&&p->chain3==NULL){
				flag=1;                      // search for the repsequence in p
				if(p->length1>=replength){
					for(int b=0;b<=p->length1-replength;b++){  //site b
						flag2=0;
						for(int a=0;a<replength;a++){
							if(p->information[0][b+a]!=repseq[a]){
								flag2=1;break;         //not in this site, to next site
							}  
						}
						if(flag2==0){flag=0;break;} //yes in this site, finish finding
					}
				}

				if(flag==0){
					for(p3=p->next;p3!=cell_head[h][y][x];p3=p3->next){
						flag1=1;                      // search for the repsequence in p3
						if(p3->length1>=replength){
							for(int b=0;b<=p3->length1-replength;b++){  //site b
								flag2=0;
								for(int a=0;a<replength;a++){
									if(p3->information[0][b+a]!=repseq[a]){
										flag2=1;break;         //not in this site, to next site
									}  
								}
								if(flag2==0){flag1=0;break;} //yes in this site, finish finding
							}
						}
						if(flag1==0){ // Template should include replicase sequence or its complementary sequence
							if(randd()<PRB&&p3->length1>=REPCOVERSEQ&&(p3->length2<p3->length1||p3->length2==p3->length1&&p3->nick!=0)&&p3->chain3==NULL){
								(p3->prior)->next=p3->next;
								if(p3->next!=cell_head[h][y][x])(p3->next)->prior=p3->prior;
								(p->prior)->next=p3;
								if(p->next!=cell_head[h][y][x])(p->next)->prior=p3;
								p3->prior=p->prior;
								p3->next=p->next;
								p3->chain3=p;
								p=p3;
								break;
							}
						}
						else{
							flagn1=1;                      // search for the repcomsequence in p3
							if(p3->length1>=replength){
								for(int b=0;b<=p3->length1-replength;b++){  //site b
									flag2=0;
									for(int a=0;a<replength;a++){
										if(p3->information[0][b+a]!=repcomseq[a]){
											flag2=1;break;         //not in this site, to next site
										}  
									}
									if(flag2==0){flagn1=0;break;} //yes in this site, finish finding
								}
							}
							if(flagn1==0){ // Template should include replicase sequence or its complementary sequence
								if(randd()<PRB&&p3->length1>=REPCOVERSEQ&&(p3->length2<p3->length1||p3->length2==p3->length1&&p3->nick!=0)&&p3->chain3==NULL){
									(p3->prior)->next=p3->next;
									if(p3->next!=cell_head[h][y][x])(p3->next)->prior=p3->prior;
									(p->prior)->next=p3;
									if(p->next!=cell_head[h][y][x])(p->next)->prior=p3;
									p3->prior=p->prior;
									p3->next=p->next;
									p3->chain3=p;
									p=p3;
									break;
								}
							}
						}
					}
				}
			}
			p1=p->prior;                          
			p2=p->next;
			p3=cell_head[!h][y][x]->next;           
			cell_head[!h][y][x]->next=p;
			p->next=p3;
			p->prior=cell_head[!h][y][x];
			if(p3!=cell_head[!h][y][x])p3->prior=p;
			p1->next=p2;                              
			if(p2!=cell_head[h][y][x])p2->prior=p1;
			p=p1;                                       
			break;

			case 5:            //moving to another adjacent cell
			if(p->chain3==NULL)length3=0;
			else length3=(p->chain3)->length1;
    		if(randd()*FDMOV<PMOV)              
			{ 
			   randcase=randl(4);   // Four possible directions
			   switch(randcase)
			   {
				case 0:
					if(x>0)
					{
						p1=p->prior;
						p2=p->next;

						p3=cell_head[!h][y][x-1]->next;           
						cell_head[!h][y][x-1]->next=p;
						p->next=p3;
						p->prior=cell_head[!h][y][x-1];
						if(p3!=cell_head[!h][y][x-1])p3->prior=p;

						p1->next=p2;                              
						if(p2!=cell_head[h][y][x])p2->prior=p1;
						p=p1;                                      
					}
					else {	p1=p->prior;                          
	p2=p->next;
	p3=cell_head[!h][y][x]->next;           
	cell_head[!h][y][x]->next=p;
	p->next=p3;
	p->prior=cell_head[!h][y][x];
	if(p3!=cell_head[!h][y][x])p3->prior=p;
	p1->next=p2;                              
	if(p2!=cell_head[h][y][x])p2->prior=p1;
	p=p1;                                       
}
					break;
				
				case 1:
					if(x<SIDE-1)
					{
						p1=p->prior;
						p2=p->next;

						p3=cell_head[!h][y][x+1]->next;           
						cell_head[!h][y][x+1]->next=p;
						p->next=p3;
						p->prior=cell_head[!h][y][x+1];
						if(p3!=cell_head[!h][y][x+1])p3->prior=p;

						p1->next=p2;                              
						if(p2!=cell_head[h][y][x])p2->prior=p1;
						p=p1;                                 
					}
					else {	p1=p->prior;                          
	p2=p->next;
	p3=cell_head[!h][y][x]->next;           
	cell_head[!h][y][x]->next=p;
	p->next=p3;
	p->prior=cell_head[!h][y][x];
	if(p3!=cell_head[!h][y][x])p3->prior=p;
	p1->next=p2;                              
	if(p2!=cell_head[h][y][x])p2->prior=p1;
	p=p1;                                       
}
					break;
				
				case 2:
					if(y>0)
					{
						p1=p->prior;
						p2=p->next;

						p3=cell_head[!h][y-1][x]->next;           
						cell_head[!h][y-1][x]->next=p;
						p->next=p3;
						p->prior=cell_head[!h][y-1][x];
						if(p3!=cell_head[!h][y-1][x])p3->prior=p;

						p1->next=p2;                              
						if(p2!=cell_head[h][y][x])p2->prior=p1;
						p=p1;                                        
					}
					else {	p1=p->prior;                          
	p2=p->next;
	p3=cell_head[!h][y][x]->next;           
	cell_head[!h][y][x]->next=p;
	p->next=p3;
	p->prior=cell_head[!h][y][x];
	if(p3!=cell_head[!h][y][x])p3->prior=p;
	p1->next=p2;                              
	if(p2!=cell_head[h][y][x])p2->prior=p1;
	p=p1;                                       
}
					break;


				case 3:
					if(y<SIDE-1)
					{
						p1=p->prior;
						p2=p->next;

						p3=cell_head[!h][y+1][x]->next;           
						cell_head[!h][y+1][x]->next=p;
						p->next=p3;
						p->prior=cell_head[!h][y+1][x];
						if(p3!=cell_head[!h][y+1][x])p3->prior=p;

						p1->next=p2;                             
						if(p2!=cell_head[h][y][x])p2->prior=p1;
						p=p1;                                     
					}
					else {	p1=p->prior;                          
	p2=p->next;
	p3=cell_head[!h][y][x]->next;           
	cell_head[!h][y][x]->next=p;
	p->next=p3;
	p->prior=cell_head[!h][y][x];
	if(p3!=cell_head[!h][y][x])p3->prior=p;
	p1->next=p2;                              
	if(p2!=cell_head[h][y][x])p2->prior=p1;
	p=p1;                                       
}
					break;

				default :printf("randmove error");
			   } 
			}  
			else {	p1=p->prior;                          
	p2=p->next;
	p3=cell_head[!h][y][x]->next;           
	cell_head[!h][y][x]->next=p;
	p->next=p3;
	p->prior=cell_head[!h][y][x];
	if(p3!=cell_head[!h][y][x])p3->prior=p;
	p1->next=p2;                              
	if(p2!=cell_head[h][y][x])p2->prior=p1;
	p=p1;                                       
}
			break;

			default: printf("rna case error");
		   }    
       }  
		}
	}
}

//------------------------------------------------------------------------------
int record(void)               // Data recording
{
	recstep[g]=i;

	rep_num[g]=0; 
	for(y=0;y<SIDE;y++){
		for(x=0;x<SIDE;x++){
			for(p=cell_head[h][y][x]->next; p!=cell_head[h][y][x]; p=p->next){
				flag1=1;                      // search for the repsequence in p
				if(p->length1>=replength){
					for(int b=0;b<=p->length1-replength;b++){  //site b
						flag2=0;
						for(int a=0;a<replength;a++){
							if(p->information[0][b+a]!=repseq[a]){
								flag2=1;break;         //not in this site, to next site
							}  
						}
						if(flag2==0){flag1=0;break;} //yes in this site, finish finding
					}
				}
				if(flag1==0) rep_num[g]++; //replicase sequence as template
				if(p->chain3!=NULL) rep_num[g]++; //replicase sequence as ribozyme
			}
		} 
	} 
	g++;
	return(0);
}

//------------------------------------------------------------------------------
void save_result(void)    // Data saving
{
	FILE *fp1,*fp2; 

	printf("over_max_length = %d times\n",over_max_len);
	printf("ready to write into file:");
	printf(RECDATA1);
	printf("\n");

	if((fp1=fopen(RECDATA1,"wb"))==NULL){printf("cannot open file1");exit(0);}
	if((fp2=fopen(RECDATA2,"wb"))==NULL){printf("cannot open file2");exit(0);}

	for(g=0;g<(STEPNUM-STAREC)/RECINT+1;g++)
	{
	  gi=g*RECINT+STAREC;
	  cout<<"step"<<gi<<":rep_num="<<rep_num[g]<<endl;
	  if(fwrite(&gi,sizeof(int),1,fp1)!=1)printf("file write error1\n");
	}
	for(g=0;g<(STEPNUM-STAREC)/RECINT+1;g++)
	{
	  if(fwrite(&rep_num[g],sizeof(float),1,fp2)!=1)printf("file write error2\n");
	}

	printf("\nend\n");

	fclose(fp1);
	fclose(fp2);
}

//------------------------------------------------------------------------------
void freepool(void)        // Memory releasing  
{
	int m;
	for(m=0;m<2;m++)
	{
		for(y=0;y<SIDE;y++)
		{
			for(x=0;x<SIDE;x++)
			{
				while(1) 
				{
					if(cell_head[m][y][x]->next!=cell_head[m][y][x])
					 {
						 p=cell_head[m][y][x]->next;
						 cell_head[m][y][x]->next=p->next;
						 if(p->chain3!=NULL)free(p->chain3);
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
	time_start=time(ptimer);

	inits();        // initialization of the system

	time1=time(ptimer);

	for(i=0; i<=STEPNUM; i++)      // Monte-Carlo cycle
	{
		if(i>=STAREC&&i%RECINT==0) 
		{
		 record();
		}
		printf("go%ld\n",i);
		form_monomer();
		unit_case();
		h=!h;
	}   
	time2=time(ptimer);

	save_result();
	time3=time(ptimer);

	freepool();
	time_end=time(ptimer);

	cout<<"init_time="<<time1-time_start<<endl;
	cout<<"monte_cycle_time="<<time2-time1<<endl;
	cout<<"save_result_time="<<time3-time2<<endl;
	cout<<"freepool_time="<<time_end-time3<<endl;
	cin.get();

	return (0);
}
//========================================================  End of the program
