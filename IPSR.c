/* A code to get the IPSR (instantaneous population spike rate) for the constituent cells */

#include <stdio.h>   /* Standard header file for the Input and Output */
#include <stdlib.h>  /* Standard library header file */
#include <math.h>    /* Math header file */
#include <conio.h>   /* Header file for the IO on the console */

#define N 1000       /* No. of cells:
                        The index of the array begins from i=1.
                        --> The size of array is N+1. */
#define NOD 10000    /* No. of spikes in the input raster plot */

unsigned long int RKCNT; /* No. of calling the Runge-Kutta intergration routine 
                            FLOW() to avoid the calculation error in time */

int SN,TN,PN;
int NOSi[N+1];
float Time,SPTime[N+1][20000];
float IISR[N+1],IPSR;
float h,BW;
char FINN[20],FOUTN[20];
FILE *FIN,*FOUT;

/* General Routines */

void KDE();                     /* A routine to obtain the individual and 
                                   population spiking rate using the kernel 
                                   density estimation */
float GaussKernel(float t);     /* A routine to return value of the Gaussian kernel */

void main()
{
    int i;

    /* Band width for the weighting function */
    BW=10.;

    /* Calculation Conditions to obtain the raster plot:
       h: RK time interval,
       TN: Transient time to eliminate the transient behavior,
       PN: Plotting time. */
    SN=10; h=1./SN; 
    TN=0; PN=2000;

    /* Input names of input and output files, and open them */
    printf("\n INPUT THE NAME OF THE INPUT FILE\n");scanf("%s",FINN);
    if((FIN = fopen(FINN,"r"))==NULL){printf("FILE OPEN ERROR...\n");exit(-1);}
    printf("\n INPUT THE NAME OF THE OUTPUT FILE\n");scanf("%s",FOUTN);
    if((FOUT = fopen(FOUTN,"w"))==NULL){printf("FILE OPEN ERROR...\n");exit(-1);}

    KDE();

    fclose(FIN); fclose(FOUT);
    printf("<***** End *****>\n");
    getche();
}

/* A routine to obtain the individual and population spiking rate using the 
   kernel density estimation */
void KDE()
{
    int i,j,k,ni,tmpi;
    float ftmp,ftmp2;

    /* Initialize No. of spike per each neuron */
    for(i=1;i<=N;i++) NOSi[i]=0;

    /* Read the spike time from the input raster plot */
    for(i=1;i<=NOD;i++) {
        fscanf(FIN,"%f %d",&Time,&ni);
        NOSi[ni]=NOSi[ni]+1;
        SPTime[ni][NOSi[ni]]=Time;
    }

    /* Calculate the individual and population spiking rate kernel estimation */
    for(i=1;i<=PN*SN;i++) {
        Time=((float)i+(float)TN*SN)*h;

        /* Calculate the individual spiking rate kernel estimation */
        for(j=1;j<=N;j++) {
            IISR[j]=0.;    /* Initialize */
            for(k=1;k<=NOSi[j];k++) {
                if(SPTime[j][k]>(Time-TN)) {
                    IISR[j]=IISR[j]+GaussKernel((Time-SPTime[j][k])/BW)/BW;
                }
                if(SPTime[j][k]>(Time+TN)) break;
            }
            IISR[j]=1000.*IISR[j];  /* 1000x: msec --> sec */
        }

        /* Calculate the population spiking rate kernel estimation */
        IPSR=0.;   /* Initialize */
        for(j=1;j<=N;j++) {
            IPSR=IPSR+IISR[j];
        }
        IPSR=IPSR/N;

        fprintf(FOUT,"%12.8f %12.8f ",Time,IPSR);
    }
}

/* A routine to return value of the rectangular kernel */
float RectKernel(float t)
{
    float KERNEL;
    
    if (fabs(t)<1.) KERNEL=0.5;
    else KERNEL=0.;
    
    return KERNEL;
}

/* A routine to return value of the triangular kernel */
float TriangKernel(float t)
{
    float KERNEL;
    
    if (fabs(t)<1.) KERNEL=1.-fabs(t);
    else KERNEL=0.;
    
    return KERNEL;
}

/* A routine to return value of the Gaussian kernel */
float GaussKernel(float t)
{
    float KERNEL;
    
    KERNEL=(1./sqrt(2.*M_PI))*exp(-t*t/2.);

    return KERNEL;
}

