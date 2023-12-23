/* A code to get the raster plots of spikes for the constituent cells */

#include <stdio.h>      /* Standard header file for the Input and Output */
#include <stdlib.h>     /* Standard library header file */
#include <math.h>       /* Math header file */
#include <conio.h>      /* Header file for the IO on the console */
#include <time.h>       /* Header file for initializing the random number generator */
                   
unsigned long int RKCNT;  /* No. of calling the Runge-Kutta intergration routine FLOW() to avoid the calculation error in time */

#define NC      1024    /* No. of GR cluster */
#define NGR     50      /* No. of GR cells in each GR cluster */
#define NPC     16      /* NPC: No. of PCs */
#define NMFVN   100

#define MaxSPIKECount  1000	    /* Maximum No. of spikes used to calculate the synaptic conductance for each neuron to save the computational memory 
                                   --> The size of array for the spike time is N+1 x MaxSPIKECount+1 */

int SN,TN,PN,NCycle;
int NOSPIKEGR[NC+1][NGR+1],NOSPIKEMF[NC+1][NGR+1],NOSPIKEGO[NC+1],NOSPIKEBC[NPC+1],NOSPIKEPC[NPC+1];
int NOSPIKEVN,NOSPIKEMFVN[NMFVN+1],NOSPIKEIO,NOSPIKEES;
int NoPreGOGR[NC+1],IPreGOGR[NC+1][10*NC],iPreGOGR[NC+1][10*NC];
int NoPreGLGO[NC+1],IPreGLGO[NC+1][10*NC];
float t,h,vGR[NC+1][NGR+1],vGO[NC+1],vBC[NPC+1],vPC[NPC+1],vVN,vIO;
float fvGR[NC+1][NGR+1],fvGO[NC+1],fvBC[NPC+1],fvPC[NPC+1],fvVN,fvIO;
float CGR,gLGR,VLGR,gAHPGR,tauAHPGR,VAHPGR,VthGR;
float CGO,gLGO,VLGO,gAHPGO,tauAHPGO,VAHPGO,VthGO;
float CBC,gLBC,VLBC,gAHPBC,tauAHPBC,VAHPBC,VthBC;
float CPC,gLPC,VLPC,gAHPPC,tauAHPPC,VAHPPC,VthPC;
float CVN,gLVN,VLVN,gAHPVN,tauAHPVN,VAHPVN,VthVN;
float CIO,gLIO,VLIO,gAHPIO,tauAHPIO,VAHPIO,VthIO;
float SPTIMEGR[NC+1][NGR+1][MaxSPIKECount+1],SPTIMEMF[NC+1][NGR+1][MaxSPIKECount+1];
float SPTIMEGO[NC+1][MaxSPIKECount+1],SPTIMEBC[NPC+1][MaxSPIKECount+1];
float SPTIMEPC[NPC+1][MaxSPIKECount+1],SPTIMEVN[MaxSPIKECount+1];
float SPTIMEMFVN[NMFVN+1][MaxSPIKECount+1],SPTIMEIO[MaxSPIKECount+1],SPTIMEES[MaxSPIKECount+1];
float fAVGMF,fT,fAVGES;
float IdcPC,IdcVN;
float ALPHAbound,CutOffTimebound;
    /* DEbound: Bound for the double exponential function,
       CutOffTimebound: Error bound for the cut-off time (ms). */
float CutOffTimeGRMFAMPA,CutOffTimeGRMFNMDA,VGRMFSyn;
float CutOffTimeGOGRAMPA,CutOffTimeGOGRNMDA,VGOGRSyn;
float CutOffTimeGRGOGABA,VGRGOSyn;
float CutOffTimeBCGRAMPA,VBCGRSyn;
float CutOffTimePCGRAMPA,CutOffTimePCIOAMPA,CutOffTimePCBCGABA,VPCGRSyn,VPCIOSyn,VPCBCSyn;
float CutOffTimeVNMFAMPA,CutOffTimeVNMFNMDA,VVNMFSyn;
float CutOffTimeVNPCGABA,VVNPCSyn;
float CutOffTimeIOESAMPA,VIOESSyn,CutOffTimeIOVNGABA,VIOVNSyn;
float tauGRMFAMPA,tauGRMFNMDA,gGRMFAMPA,gGRMFNMDA,JGRMF;
float tauGOGRAMPA,tauGOGRNMDA1,tauGOGRNMDA2,AGOGRNMDA1,AGOGRNMDA2,gGOGRAMPA,gGOGRNMDA,JGOGR;
float tauGRGOGABA1,tauGRGOGABA2,AGRGOGABA1,AGRGOGABA2,gGRGOGABA,JGRGO;
float tauBCGRAMPA,gBCGRAMPA,JBCGR;
float tauPCGRAMPA,gPCGRAMPA,JPCGR[NPC+1][NC+1][NGR+1],JPCGR0,tauPCBCGABA,gPCBCGABA,JPCBC;
float tauPCIOAMPA,gPCIOAMPA,JPCIO;
float tauVNMFAMPA,tauVNMFNMDA,gVNMFAMPA,gVNMFNMDA,JVNMF,tauVNPCGABA,gVNPCGABA,JVNPC;
float tauIOESAMPA,gIOESAMPA,JIOES,tauIOVNGABA,gIOVNGABA,JIOVN;
float PGOGR,PGLGO;
float deltaLTD,deltaLTP,CutOffTimeLTD,dTl,dTr,A,B,t0,sigma;

char FOUTN1[20],FOUTN2[20],FOUTN3[20],FOUTN4[20],FOUTN5[20],FOUTN6[20],FOUTN7[20];
FILE *FOUT1,*FOUT2,*FOUT3,*FOUT4,*FOUT5,*FOUT6,*FOUT7;

/* General Routines */
void TS();              /* A routine to obtain the global and individual outputs 
                           and the raster plot of spikings */
void FLOW();		    /* A routine to integrate the differential equations */
float CutOff1(float tau);     /* A Routine to find the cut-off time of the alpha
                                function for the local synapse by using the bisection method */
float CutOff2(float A1, float A2, float tau1, float tau2);     /* A Routine to find the cut-off time of the alpha
                       function for the local synapse by using the bisection method */
float ALPHA1(float t, float tau);  /* A Routine to return the value of the alpha
                       function for the local synaptic current */
float ALPHA2(float t, float A1, float A2, float tau1, float tau2);  /* A Routine to return the value of the alpha
                       function for the local synaptic current */
float UniRandom(); 	    /* A routine to generate a uniform random number on (0,1) 
                           using the C-supplied rand function */

/* Routines dependent on the system */
void F();               /* A routine to return the function values in the 
                           differential equations */

/* Cerebellar network in Yamazaki & Nagao (2012) */

void main()
{
    int i,j,k,l,m;
    int Idx;
    float eta;
    time_t ttt;

    /* Initialize the random number generator */
    srand((unsigned) time(&ttt));

    /* Parameters for the each cell */
    CGR=3.1; gLGR=0.43; VLGR=-58.0; gAHPGR=1.0; tauAHPGR=5.0; VAHPGR=-82.0; VthGR=-35.0;        /* GR */
    CGO=28.0; gLGO=2.3; VLGO=-55.0; gAHPGO=20.0; tauAHPGO=5.0; VAHPGO=-72.7; VthGO=-52.0;       /* GO */
    CBC=107.0; gLBC=2.32; VLBC=-68.0; gAHPBC=100.0; tauAHPBC=2.5; VAHPBC=-70.0; VthBC=-55.0;    /* BC */
    CPC=107.0; gLPC=2.32; VLPC=-68.0; gAHPPC=100.0; tauAHPPC=5.0; VAHPPC=-70.0; VthPC=-55.0;    /* PC */
    CVN=122.3; gLVN=1.63; VLVN=-56.0; gAHPVN=50.0; tauAHPVN=2.5; VAHPVN=-70.0; VthVN=-38.8;     /* VN */
    CIO=10.; gLIO=0.67; VLIO=-60.0; gAHPIO=1.0; tauAHPIO=10.; VAHPIO=-75.0; VthIO=-50.;         /* IO */

    /* Parameters for the synaptic currents */
    tauGRMFAMPA=1.2; gGRMFAMPA=0.18; VGRMFSyn=0.; tauGRMFNMDA=52.0; gGRMFNMDA=0.025; JGRMF=8.;  /* MF to GR */
    tauGOGRAMPA=1.5; gGOGRAMPA=45.5; VGOGRSyn=0.; tauGOGRNMDA1=31.0; tauGOGRNMDA2=170.0; 
    AGOGRNMDA1=0.33; AGOGRNMDA2=(1.-AGOGRNMDA1); gGOGRNMDA=30.0; JGOGR=0.00004;                 /* GR to GO */
    tauGRGOGABA1=7.0; tauGRGOGABA2=59.0; gGRGOGABA=0.028; VGRGOSyn=-82.0; 
    AGRGOGABA1=0.43; AGRGOGABA2=(1.-AGRGOGABA1); JGRGO=10.0;                                    /* GO to GR */
    tauBCGRAMPA=8.3; gBCGRAMPA=0.7; VBCGRSyn=0.; JBCGR=0.006;                                   /* GR to BC (PF) */
    tauPCGRAMPA=8.3; gPCGRAMPA=0.7; VPCGRSyn=0.; JPCGR0=0.006;                                  /* GR to PC (PF) */
    for(i=1;i<=NPC;i++) {
        for(j=1;j<=NC;j++) {
            for(k=1;k<=NGR;k++) JPCGR[i][j][k]=JPCGR0;
        }
    }
    tauPCIOAMPA=8.3; gPCIOAMPA=0.7; VPCIOSyn=0.; JPCIO=1.0;                                     /* IO to PC */
    tauPCBCGABA=10.0; gPCBCGABA=1.0; VPCBCSyn=-75.0; JPCBC=5.3;                                 /* BC to PC */
    tauVNMFAMPA=9.9; gVNMFAMPA=50.0; VVNMFSyn=0.; tauVNMFNMDA=30.6; gVNMFNMDA=25.8; JVNMF=0.002;    /* MF to VN */
    tauVNPCGABA=42.3; gVNPCGABA=30.0; VVNPCSyn=-88.0; JVNPC=0.008;                              /* PC to VN */
    tauIOESAMPA=10.0; gIOESAMPA=1.0; VIOESSyn=0.; tauIOVNGABA=10.0; gIOVNGABA=0.18; JIOVN=5.0;  /* VN to IO */
    JIOES=1.0;

    /* Parameters for the synaptic plasticity rule */
    deltaLTD=0.005; deltaLTP=0.0005; CutOffTimeLTD=50.;
    dTl=-117.5; dTr=277.5;
    A=-0.09; B=0.31; t0=80.0; sigma=180.;

    /* Constant stimuli for PC and VN */
    IdcPC=250.;
    IdcVN=700.;

    /* ALPHAbound: Bound for the alpha function,
       CutOffTimebound: Error bound for the cut-off time (ms). */
    ALPHAbound=1.e-5; CutOffTimebound=1.e-5;
    CutOffTimeGRMFAMPA=CutOff1(tauGRMFAMPA);    /* Find the cut-off time for the local synaptic current */
    CutOffTimeGRMFNMDA=CutOff1(tauGRMFNMDA);    /* Find the cut-off time for the local synaptic current */
    CutOffTimeGOGRAMPA=CutOff1(tauGOGRAMPA);    /* Find the cut-off time for the local synaptic current */
    CutOffTimeGOGRNMDA=CutOff2(AGOGRNMDA1,AGOGRNMDA2,tauGOGRNMDA1,tauGOGRNMDA2);    /* Find the cut-off time for the local synaptic current */
    CutOffTimeGRGOGABA=CutOff2(AGRGOGABA1,AGRGOGABA2,tauGRGOGABA1,tauGRGOGABA2);    /* Find the cut-off time for the local synaptic current */
    CutOffTimeBCGRAMPA=CutOff1(tauBCGRAMPA);    /* Find the cut-off time for the local synaptic current */
    CutOffTimePCGRAMPA=CutOff1(tauPCGRAMPA);    /* Find the cut-off time for the local synaptic current */
    CutOffTimePCIOAMPA=CutOff1(tauPCIOAMPA);    /* Find the cut-off time for the local synaptic current */
    CutOffTimePCBCGABA=CutOff1(tauPCBCGABA);    /* Find the cut-off time for the local synaptic current */
    CutOffTimeVNMFAMPA=CutOff1(tauVNMFAMPA);    /* Find the cut-off time for the local synaptic current */
    CutOffTimeVNMFNMDA=CutOff1(tauVNMFNMDA);    /* Find the cut-off time for the local synaptic current */
    CutOffTimeVNPCGABA=CutOff1(tauVNPCGABA);    /* Find the cut-off time for the local synaptic current */
    CutOffTimeIOESAMPA=CutOff1(tauIOESAMPA);    /* Find the cut-off time for the local synaptic current */
    CutOffTimeIOVNGABA=CutOff1(tauIOVNGABA);    /* Find the cut-off time for the local synaptic current */

    /* Parameters of the Poisson spike train generator tor for the post-movement
       sensory signal via MF
       fT (kHz): Frequency of the post-movement sensory signal,
       fAVGMF (kHz): Mean firing rate of sensory signal of MF */
    fAVGMF=15./1000.;fT=0.5/1000.; fAVGES=1.5/1000.;

    /* Connection probability */
    PGOGR=0.1;          /* GR to GO */
    PGLGO=0.06;         /* GO to GR */

    /* Calculation Conditions:
       SN: No. of Runge-Kutta (RK) steps per unit time (1ms),
       h: RK time interval,
       TN: Transient time to eliminate the transient behavior,
       PN: Plotting time. */
    SN=1; h=1./SN; TN=0; PN=2000;

    /* Input names of output files, and open them */
    printf("\n INPUT THE NAME OF THE OUTPUT FILE 1 (Raster plot for GR) \n");scanf("%s",FOUTN1);
    if((FOUT1 = fopen(FOUTN1,"w"))==NULL){printf("FILE OPEN ERROR...\n");exit(-1);}
    printf("\n INPUT THE NAME OF THE OUTPUT FILE 2 (Raster plot for GO) \n");scanf("%s",FOUTN2);
    if((FOUT2 = fopen(FOUTN2,"w"))==NULL){printf("FILE OPEN ERROR...\n");exit(-1);}
    printf("\n INPUT THE NAME OF THE OUTPUT FILE 3 (Raster plot for PC) \n");scanf("%s",FOUTN3);
    if((FOUT3 = fopen(FOUTN3,"w"))==NULL){printf("FILE OPEN ERROR...\n");exit(-1);}
    printf("\n INPUT THE NAME OF THE OUTPUT FILE 4 (Raster plot for BC) \n");scanf("%s",FOUTN4);
    if((FOUT4 = fopen(FOUTN4,"w"))==NULL){printf("FILE OPEN ERROR...\n");exit(-1);}
    printf("\n INPUT THE NAME OF THE OUTPUT FILE 5 (Raster plot for VN) \n");scanf("%s",FOUTN5);
    if((FOUT5 = fopen(FOUTN5,"w"))==NULL){printf("FILE OPEN ERROR...\n");exit(-1);}
    printf("\n INPUT THE NAME OF THE OUTPUT FILE 6 (Raster plot for IO) \n");scanf("%s",FOUTN6);
    if((FOUT6 = fopen(FOUTN6,"w"))==NULL){printf("FILE OPEN ERROR...\n");exit(-1);}
    printf("\n INPUT THE NAME OF THE OUTPUT FILE 7 (Active PF-PC synaptic weight) \n");scanf("%s",FOUTN7);
    if((FOUT7 = fopen(FOUTN7,"w"))==NULL){printf("FILE OPEN ERROR...\n");exit(-1);}

    /* Gerenation of connection from GR to GO */
    for(i=1;i<=NC;i++) {
        NoPreGOGR[i]=0;
        for(j=i-24;j<=i+24;j++) {
            if(j<1) Idx=j+NC;
            if(j>NC) Idx=j-NC;
            for(k=1;k<=NGR;k++) {
                eta=UniRandom();
                if(eta<PGOGR) {
                    NoPreGOGR[i]=NoPreGOGR[i]+1;
                    IPreGOGR[i][NoPreGOGR[i]]=Idx;
                    iPreGOGR[i][NoPreGOGR[i]]=k;
                }
            }
        }
    }

    /* Gerenation of connection from GO to GR for GL */
    for(i=1;i<=NC;i++) {
        NoPreGLGO[i]=0;
        for(j=i-40;j<=i+40;j++) {
            if(j<1) Idx=j+NC;
            if(j>NC) Idx=j-NC;
            eta=UniRandom();
            if(eta<PGLGO) {
                NoPreGLGO[i]=NoPreGLGO[i]+1;
                IPreGLGO[i][NoPreGLGO[i]]=Idx;
            }
        }
    }

    /* Initialize variables */
    for(i=1;i<=NC;i++) {
        for(j=1;j<=NGR;j++) vGR[i][j]=VLGR+10.*(UniRandom()-0.5);
    }
    for(i=1;i<=NC;i++) vGO[i]=VLGO+10.*(UniRandom()-0.5);
    for(i=1;i<=NPC;i++) vBC[i]=VLBC+10.*(UniRandom()-0.5);
    for(i=1;i<=NPC;i++) vPC[i]=VLPC+10.*(UniRandom()-0.5);
    vVN=VLVN+10.*(UniRandom()-0.5);
    vIO=VLIO+10.*(UniRandom()-0.5);

    /* Begin the Main Calculation */
    for(NCycle=1;NCycle<=400;NCycle++) TS();

    fclose(FOUT1);fclose(FOUT2);fclose(FOUT3);fclose(FOUT4);fclose(FOUT5);fclose(FOUT6);fclose(FOUT7);
    printf("<***** End *****>\n");
    getche();
}

/* A routine to obtain the global and individual outputs and the raster plot 
   of spikings */
void TS()
{
    int i,j,k,l,m,Idx;
    int CF[NPC+1],PF[NC+1][NGR+1],DF;
    float fPoissonMF,pr,eta,fPoissonES,prES;
    float tCFeff,tPFeff,dT,dWJ;

    /* Initialize variables */
    for(i=1;i<=NC;i++) {
        for(j=1;j<=NGR;j++) NOSPIKEMF[i][j]=0;
    }
    for(i=1;i<=NC;i++) {
        for(j=1;j<=NGR;j++) {
            NOSPIKEGR[i][j]=0; SPTIMEGR[i][j][1]=-10.;
        }
    }
    for(i=1;i<=NC;i++) {
        NOSPIKEGO[i]=0; SPTIMEGO[i][1]=-10.;
    }
    for(i=1;i<=NPC;i++) {
        NOSPIKEBC[i]=0; SPTIMEBC[i][1]=-10.;
    }
    for(i=1;i<=NPC;i++) {
        NOSPIKEPC[i]=0; SPTIMEPC[i][1]=-10.;
    }
    for(i=1;i<=NMFVN;i++) {NOSPIKEMFVN[i]=0;}
    NOSPIKEVN=0; SPTIMEVN[1]=-10.;
    NOSPIKEES=0;
    NOSPIKEIO=0; SPTIMEIO[1]=-10.;

    /* Initialize the No. of call the Runge-Kutta intergration routine FLOW() */
    RKCNT=0; t=0.;

    for(i=1;i<=PN*SN;i++) {
        FLOW();

        /* Initialize */
        for(j=1;j<=NPC;j++) CF[j]=0;
        for(j=1;j<=NC;j++) {
            for(k=1;k<=NGR;k++) PF[j][k]=0;
        }

        /* Generation of MF sensory signal */
        fPoissonMF=-fAVGMF*cos(2.*M_PI*fT*t)+fAVGMF;
        pr=fPoissonMF*h;
        for(j=1;j<=NC;j++) {
            for(k=1;k<=NGR;k++) {
                eta=UniRandom();
                if(eta<pr) {
                    /* Increase in the number of spikes of MF */
                    NOSPIKEMF[j][k]=NOSPIKEMF[j][k]+1;
                    if(NOSPIKEMF[j][k]>MaxSPIKECount) NOSPIKEMF[j][k]=MaxSPIKECount;
                    for(l=NOSPIKEMF[j][k];l>1;l--) {SPTIMEMF[j][k][l]=SPTIMEMF[j][k][l-1];}
                    SPTIMEMF[j][k][1]=t;
                }
            }
        }
                
        /* Find a spike in GR cell */
        for(j=1;j<=NC;j++) {
            for(k=1;k<=NGR;k++) {
                if(vGR[j][k]>=VthGR) {
                    fprintf(FOUT1,"%4d %12.8f %4d %4d\n",NCycle,t,j,k); 
                    for(l=1;l<=NPC;l++) fprintf(FOUT7,"%4d %12.8f %12.8f %4d %4d %4d\n",NCycle,t,JPCGR[j][k][l],j,k,l);
                    NOSPIKEGR[j][k]=NOSPIKEGR[j][k]+1;
                    if(NOSPIKEGR[j][k]>MaxSPIKECount) NOSPIKEGR[j][k]=MaxSPIKECount;
                    for(l=NOSPIKEGR[j][k];l>1;l--) {SPTIMEGR[j][k][l]=SPTIMEGR[j][k][l-1];}
                    SPTIMEGR[j][k][1]=t;
                    PF[j][k]=1;
                }
            }
        }

        /* Find a spike in GO cell */
        for(j=1;j<=NC;j++) {
            if(vGO[j]>=VthGO) {
                fprintf(FOUT2,"%4d %12.8f %4d\n",NCycle,t,j); 
                NOSPIKEGO[j]=NOSPIKEGO[j]+1;
                if(NOSPIKEGO[j]>MaxSPIKECount) NOSPIKEGO[j]=MaxSPIKECount;
                for(k=NOSPIKEGO[j];k>1;k--) {SPTIMEGO[j][k]=SPTIMEGO[j][k-1];}
                SPTIMEGO[j][1]=t;
            }
        }

        /* Find a spike in BC */
        for(j=1;j<=NPC;j++) {
            if(vBC[j]>=VthBC) {
                fprintf(FOUT4,"%4d %12.8f %4d\n",NCycle,t,j); 
                NOSPIKEBC[j]=NOSPIKEBC[j]+1;
                if(NOSPIKEBC[j]>MaxSPIKECount) NOSPIKEBC[j]=MaxSPIKECount;
                for(k=NOSPIKEBC[j];k>1;k--) {SPTIMEBC[j][k]=SPTIMEBC[j][k-1];}
                SPTIMEBC[j][1]=t;
            }
        }
        
        /* Find a spike in PC */
        for(j=1;j<=NPC;j++) {
            if(vPC[j]>=VthPC) {
                fprintf(FOUT3,"%4d %12.8f %4d\n",NCycle,t,j); 
                NOSPIKEPC[j]=NOSPIKEPC[j]+1;
                if(NOSPIKEPC[j]>MaxSPIKECount) NOSPIKEPC[j]=MaxSPIKECount;
                for(k=NOSPIKEPC[j];k>1;k--) {SPTIMEPC[j][k]=SPTIMEPC[j][k-1];}
                SPTIMEPC[j][1]=t;
                CF[j]=1;
            }
        }

        /* Generation of MF sensory signal into the VN */
        for(j=1;j<=NMFVN;j++) {
            eta=UniRandom();
            if(eta<pr) {
                /* Increase in the number of spikes of MF */
                NOSPIKEMFVN[j]=NOSPIKEMFVN[j]+1;
                if(NOSPIKEMFVN[j]>MaxSPIKECount) NOSPIKEMFVN[j]=MaxSPIKECount;
                for(k=NOSPIKEMFVN[j];k>1;k--) {SPTIMEMFVN[j][k]=SPTIMEMFVN[j][k-1];}
                SPTIMEMFVN[j][1]=t;
            }
        }

        /* Find a spike in VN */
        if(vVN>=VthVN) {
            fprintf(FOUT5,"%4d %12.8f %4d\n",NCycle,t,1);
            NOSPIKEVN=NOSPIKEVN+1;
            if(NOSPIKEVN>MaxSPIKECount) NOSPIKEVN=MaxSPIKECount;
            for(j=NOSPIKEVN;j>1;j--) {SPTIMEVN[j]=SPTIMEVN[j-1];}
            SPTIMEVN[1]=t;
        }

        /* Generation of MF sensory signal into the VN */
        fPoissonES=-fAVGES*cos(2.*M_PI*fT*t)+fAVGES;
        prES=fPoissonES*h;
        eta=UniRandom();
        if(eta<prES) {
            /* Increase in the number of spikes of MF */
            NOSPIKEES=NOSPIKEES+1;
            if(NOSPIKEES>MaxSPIKECount) NOSPIKEES=MaxSPIKECount;
            for(j=NOSPIKEES;j>1;j--) {SPTIMEES[j]=SPTIMEES[j-1];}
            SPTIMEES[1]=t;
        }

        /* Find a spike in IO */
        if(vIO>=VthIO) {
            fprintf(FOUT6,"%4d %12.8f %4d\n",NCycle,t,1);
            NOSPIKEIO=NOSPIKEIO+1;
            if(NOSPIKEIO>MaxSPIKECount) NOSPIKEIO=MaxSPIKECount;
            for(j=NOSPIKEIO;j>1;j--) {SPTIMEIO[j]=SPTIMEIO[j-1];}
            SPTIMEIO[1]=t; 
            for(j=1;j<=NPC;j++) CF[j]=1;
        }
        
        /* Plasticity rule */
        for(j=1;j<=NPC;j++) {
            if(CF[j]==1) {      /* CF firing --> LTD */
                /* Calculate the effective range */
                tCFeff=t-dTr;
                if(tCFeff<0.) tCFeff=0.;
                for(k=(NC-143)+(j-1)*64;k<=(NC-143)+(j-1)*64+287;k++) { /* Find the PF firing of pre-synaptic GR cells in the effective range */
                    if(k<1) Idx=k+NC;
                    if(k>NC) Idx=k-NC;
                    for(l=1;l<=NGR;l++) {
                        for(m=1;m<=NOSPIKEGR[k][l];m++) {
                            if(SPTIMEGR[k][l][m] < tCFeff) break;
                            dT=t-SPTIMEGR[k][l][m];
                            dWJ=A+B*exp(-pow((dT-t0)/sigma,2.));
                            JPCGR[j][k][l]=JPCGR[j][k][l]-deltaLTD*JPCGR[j][k][l]*dWJ;
                        }
                    }
                }
            }
            else {
                for(k=(NC-143)+(j-1)*64;k<=(NC-143)+(j-1)*64+287;k++) { /* Find the PF firing of pre-synaptic GR cells */
                    if(k<1) Idx=k+NC;
                    if(k>NC) Idx=k-NC;
                    for(l=1;l<=NGR;l++) {
                        if(PF[k][l]==1) {   /* PF firing */
                            DF=0;
                            /* Calculate the effective range */
                            tPFeff=t+dTl;
                            if(tPFeff<0.) tPFeff=0.;
                            for(m=1;m<=NOSPIKEPC[j];m++) {
                                if(SPTIMEPC[j][m] < tPFeff) break;
                                dT=t-SPTIMEPC[j][m];
                                dWJ=A+B*exp(-pow((dT-t0)/sigma,2.));
                                JPCGR[j][k][l]=JPCGR[j][k][l]-deltaLTD*JPCGR[j][k][l]*dWJ;
                                DF=1;
                            }
                            if(DF==0) {
                                JPCGR[j][k][l]=JPCGR[j][k][l]+deltaLTP*(JPCGR0-JPCGR[j][k][l]);
                            }
                        }
                        /* No PF firing --> No Plasticity */
                    }
                }
            }
        }
    }
}

/* A routine to integrate the differential equations */
void FLOW()
{
    int i,j;
    float t0,vGR0[NC+1][NGR+1],vGO0[NC+1],vBC0[NPC+1],vPC0[NPC+1],vVN0,vIO0;
    float fvGR1[NC+1][NGR+1],fvGR2[NC+1][NGR+1],fvGO1[NC+1],fvGO2[NC+1];
    float fvBC1[NPC+1],fvBC2[NPC+1],fvPC1[NPC+1],fvPC2[NPC+1];
    float fvVN1,fvVN2,fvIO1,fvIO2;

    t0=(float)RKCNT/SN;
    for(i=1;i<=NC;i++) {
        for(j=1;j<=NGR;j++) vGR0[i][j]=vGR[i][j];
    }
    for(i=1;i<=NC;i++) vGO0[i]=vGO[i];
    for(i=1;i<=NPC;i++) vBC0[i]=vBC[i];
    for(i=1;i<=NPC;i++) vPC0[i]=vPC[i];
    vVN0=vVN; vIO0=vIO;

    t=t0;
    F();
    for(i=1;i<=NC;i++) {
        for(j=1;j<=NGR;j++) {
            fvGR1[i][j]=fvGR[i][j];
            vGR[i][j]=vGR0[i][j]+h*fvGR1[i][j];
        }
    }
    for(i=1;i<=NC;i++) {
        fvGO1[i]=fvGO[i];
        vGO[i]=vGO0[i]+h*fvGO1[i];
    }
    for(i=1;i<=NPC;i++) {
        fvBC1[i]=fvBC[i];
        vBC[i]=vBC0[i]+h*fvBC1[i];
    }
    for(i=1;i<=NPC;i++) {
        fvPC1[i]=fvPC[i];
        vPC[i]=vPC0[i]+h*fvPC1[i];
    }
    fvVN1=fvVN; vVN=vVN0+h*fvVN1;
    fvIO1=fvIO; vIO=vIO0+h*fvIO1;

    t=t0+h; RKCNT++;
    F();
    for(i=1;i<=NC;i++) {
        for(j=1;j<=NGR;j++) {
            fvGR2[i][j]=fvGR[i][j];
            vGR[i][j]=vGR0[i][j]+h*(fvGR1[i][j]+fvGR2[i][j])/2.;
        }
    }
    for(i=1;i<=NC;i++) {
        fvGO2[i]=fvGO[i];
        vGO[i]=vGO0[i]+h*(fvGO1[i]+fvGO2[i])/2.;
    }
    for(i=1;i<=NPC;i++) {
        fvBC2[i]=fvBC[i];
        vBC[i]=vBC0[i]+h*(fvBC1[i]+fvBC2[i])/2.;
    }
    for(i=1;i<=NPC;i++) {
        fvPC2[i]=fvPC[i];
        vPC[i]=vPC0[i]+h*(fvPC1[i]+fvPC2[i])/2.;
    }
    fvVN2=fvVN; vVN=vVN0+h*(fvVN1+fvVN2)/2.;
    fvIO2=fvIO; vIO=vIO0+h*(fvIO1+fvIO2)/2.;
}

/* A routine to return the function values in the differential equations */
void F()
{
    int i,j,k,l,m,Idx;
    float ILGR,IAHPGR,ILGO,IAHPGO,ILBC,IAHPBC,ILPC,IAHPPC,ILVN,IAHPVN,ILIO,IAHPIO;
    float SGRMFAMPA[NC+1][NGR+1],SGRMFNMDA[NC+1][NGR+1],gGRMFAMPASyn,gGRMFNMDASyn,IGRMFAMPASyn,IGRMFNMDASyn;
    float SGOGRAMPA[NC+1][NGR+1],SGOGRNMDA[NC+1][NGR+1],gGOGRAMPASyn,gGOGRNMDASyn,IGOGRAMPASyn,IGOGRNMDASyn;
    float SGRGOGABA[NC+1],gGRGOGABASyn,IGRGOGABASyn;
    float SBCGRAMPA[NC+1][NGR+1],gBCGRAMPASyn,IBCGRAMPASyn;
    float SPCGRAMPA[NC+1][NGR+1],gPCGRAMPASyn,IPCGRAMPASyn;
    float SPCIOAMPA,gPCIOAMPASyn,IPCIOAMPASyn;
    float SPCBCGABA[NPC+1],gPCBCGABASyn,IPCBCGABASyn;
    float SVNMFAMPA[NMFVN+1],SVNMFNMDA[NMFVN+1],gVNMFAMPASyn,gVNMFNMDASyn,IVNMFAMPASyn,IVNMFNMDASyn;
    float SVNPCGABA[NPC+1],gVNPCGABASyn,IVNPCGABASyn;
    float SIOESAMPA,gIOESAMPASyn,IIOESAMPASyn,SIOVNGABA,gIOVNGABASyn,IIOVNGABASyn;

    /* Calculate the summation of s[i][j] for MF to GR */
	for(i=1;i<=NC;i++) {
 	    for(j=1;j<=NGR;j++) {
            SGRMFAMPA[i][j]=0.;
            for(k=1;k<=NOSPIKEMF[i][j];k++) {
                if((t-SPTIMEMF[i][j][k]) > CutOffTimeGRMFAMPA) break;
                /* Spikes in (t-tau_l-CutOffTime)<SPTIME[i][j]<t-tau_l are 
                   considered for the synaptic current. */
                SGRMFAMPA[i][j]=SGRMFAMPA[i][j]+ALPHA1(t-SPTIMEMF[i][j][k],tauGRMFAMPA);
            }
        }
    }
	for(i=1;i<=NC;i++) {
 	    for(j=1;j<=NGR;j++) {
            SGRMFNMDA[i][j]=0.;
            for(k=1;k<=NOSPIKEMF[i][j];k++) {
                if((t-SPTIMEMF[i][j][k]) > CutOffTimeGRMFNMDA) break;
                /* Spikes in (t-tau_l-CutOffTime)<SPTIME[i][j]<t-tau_l are 
                   considered for the synaptic current. */
                SGRMFNMDA[i][j]=SGRMFNMDA[i][j]+ALPHA1(t-SPTIMEMF[i][j][k],tauGRMFNMDA);
            }
        }
    }

    /* Calculate the summation of s[i][j] for GR to GO */
	for(i=1;i<=NC;i++) {
 	    for(j=1;j<=NGR;j++) {
            SGOGRAMPA[i][j]=0.;
            for(k=1;k<=NOSPIKEGR[i][j];k++) {
                if((t-SPTIMEGR[i][j][k]) > CutOffTimeGOGRAMPA) break;
                /* Spikes in (t-tau_l-CutOffTime)<SPTIME[i][j]<t-tau_l are 
                   considered for the synaptic current. */
                SGOGRAMPA[i][j]=SGOGRAMPA[i][j]+ALPHA1(t-SPTIMEGR[i][j][k],tauGOGRAMPA);
            }
        }
    }
	for(i=1;i<=NC;i++) {
	    for(j=1;j<=NGR;j++) {
            SGOGRNMDA[i][j]=0.;
            for(k=1;k<=NOSPIKEGR[i][j];k++) {
                if((t-SPTIMEGR[i][j][k]) > CutOffTimeGOGRNMDA) break;
                /* Spikes in (t-tau_l-CutOffTime)<SPTIME[i][j]<t-tau_l are 
                   considered for the synaptic current. */
                SGOGRNMDA[i][j]=SGOGRNMDA[i][j]+ALPHA2(t-SPTIMEGR[i][j][k],AGOGRNMDA1,AGOGRNMDA2,tauGOGRNMDA1,tauGOGRNMDA2);
            }
        }
    }

    /* Calculate the summation of s[i][j] for GO to GR */
	for(i=1;i<=NC;i++) {
        SGRGOGABA[i]=0.;
        for(j=1;j<=NOSPIKEGO[i];j++) {
            if((t-SPTIMEGO[i][j]) > CutOffTimeGRGOGABA) break;
            /* Spikes in (t-tau_l-CutOffTime)<SPTIME[i][j]<t-tau_l are 
               considered for the synaptic current. */
            SGRGOGABA[i]=SGRGOGABA[i]+ALPHA2(t-SPTIMEGO[i][j],AGRGOGABA1,AGRGOGABA2,tauGRGOGABA1,tauGRGOGABA2);
        }
    }

    /* Calculate the summation of s[i][j] for GR to BC */
	for(i=1;i<=NC;i++) {
 	    for(j=1;j<=NGR;j++) {
            SBCGRAMPA[i][j]=0.;
            for(k=1;k<=NOSPIKEGR[i][j];k++) {
                if((t-SPTIMEGR[i][j][k]) > CutOffTimeBCGRAMPA) break;
                /* Spikes in (t-tau_l-CutOffTime)<SPTIME[i][j]<t-tau_l are 
                   considered for the synaptic current. */
                SBCGRAMPA[i][j]=SBCGRAMPA[i][j]+ALPHA1(t-SPTIMEGR[i][j][k],tauBCGRAMPA);
            }
        }
    }

    /* Calculate the summation of s[i][j] for GR to PC */
	for(i=1;i<=NC;i++) {
 	    for(j=1;j<=NGR;j++) {
            SPCGRAMPA[i][j]=0.;
            for(k=1;k<=NOSPIKEGR[i][j];k++) {
                if((t-SPTIMEGR[i][j][k]) > CutOffTimePCGRAMPA) break;
                /* Spikes in (t-tau_l-CutOffTime)<SPTIME[i][j]<t-tau_l are 
                   considered for the synaptic current. */
                SPCGRAMPA[i][j]=SPCGRAMPA[i][j]+ALPHA1(t-SPTIMEGR[i][j][l],tauPCGRAMPA);
            }
        }
    }

    /* Calculate the summation of s[i][j] for CF to PC */
    SPCIOAMPA=0.;
    for(i=1;i<=NOSPIKEIO;i++) {
        if((t-SPTIMEIO[i]) > CutOffTimePCIOAMPA) break;
        /* Spikes in (t-tau_l-CutOffTime)<SPTIME[i][j]<t-tau_l are 
           considered for the synaptic current. */
        SPCIOAMPA=SPCIOAMPA+ALPHA1(t-SPTIMEIO[i],tauPCIOAMPA);
    }

    /* Calculate the summation of s[i][j] for BC to PC */
	for(i=1;i<=NPC;i++) {
        SPCBCGABA[i]=0.;
        for(j=1;j<=NOSPIKEBC[i];j++) {
            if((t-SPTIMEBC[i][j]) > CutOffTimePCBCGABA) break;
            /* Spikes in (t-tau_l-CutOffTime)<SPTIME[i][j]<t-tau_l are 
               considered for the synaptic current. */
            SPCBCGABA[i]=SPCBCGABA[i]+ALPHA1(t-SPTIMEBC[i][j],tauPCBCGABA);
        }
    }

    /* Calculate the summation of s[i][j] for MF to VN */
	for(i=1;i<=NMFVN;i++) {
        SVNMFAMPA[i]=0.;
        for(j=1;j<=NOSPIKEMFVN[i];j++) {
            if((t-SPTIMEMFVN[i][j]) > CutOffTimeVNMFAMPA) break;
            /* Spikes in (t-tau_l-CutOffTime)<SPTIME[i][j]<t-tau_l are 
               considered for the synaptic current. */
            SVNMFAMPA[i]=SVNMFAMPA[i]+ALPHA1(t-SPTIMEMFVN[i][j],tauVNMFAMPA);
        }
    }
	for(i=1;i<=NMFVN;i++) {
        SVNMFNMDA[i]=0.;
        for(j=1;j<=NOSPIKEMFVN[i];j++) {
            if((t-SPTIMEMFVN[i][j]) > CutOffTimeVNMFNMDA) break;
            /* Spikes in (t-tau_l-CutOffTime)<SPTIME[i][j]<t-tau_l are 
               considered for the synaptic current. */
            SVNMFNMDA[i]=SVNMFNMDA[i]+ALPHA1(t-SPTIMEMFVN[i][j],tauVNMFNMDA);
        }
    }

    /* Calculate the summation of s[i][j] for PC to VN */
	for(i=1;i<=NPC;i++) {
        SVNPCGABA[i]=0.;
        for(j=1;j<=NOSPIKEPC[i];j++) {
            if((t-SPTIMEPC[i][j]) > CutOffTimeVNPCGABA) break;
            /* Spikes in (t-tau_l-CutOffTime)<SPTIME[i][j]<t-tau_l are 
               considered for the synaptic current. */
            SVNPCGABA[i]=SVNPCGABA[i]+ALPHA1(t-SPTIMEPC[i][j],tauVNPCGABA);
        }
    }

    /* Calculate the summation of s[i][j] for ES to IO */
    SIOESAMPA=0.;
    for(i=1;i<=NOSPIKEES;i++) {
        if((t-SPTIMEES[i]) > CutOffTimeIOESAMPA) break;
        /* Spikes in (t-tau_l-CutOffTime)<SPTIME[i][j]<t-tau_l are 
           considered for the synaptic current. */
        SIOESAMPA=SIOESAMPA+ALPHA1(t-SPTIMEES[i],tauIOESAMPA);
    }

    /* Calculate the summation of s[i][j] for VN to IO */
    SIOVNGABA=0.;
    for(i=1;i<=NOSPIKEVN;i++) {
        if((t-SPTIMEVN[i]) > CutOffTimeIOVNGABA) break;
        /* Spikes in (t-tau_l-CutOffTime)<SPTIME[i][j]<t-tau_l are 
           considered for the synaptic current. */
        SIOVNGABA=SIOVNGABA+ALPHA1(t-SPTIMEVN[i],tauIOVNGABA);
    }

    /* Governing Eqs. for GR cells */
    for(i=1;i<=NC;i++) {
        for(j=1;j<=NGR;j++) {
            ILGR=gLGR*(vGR[i][j]-VLGR);
            if(SPTIMEGR[i][j][1]>0.) IAHPGR=gAHPGR*exp(-(t-SPTIMEGR[i][j][1])/tauAHPGR)*(vGR[i][j]-VAHPGR);
            else IAHPGR=0.;

            gGRMFAMPASyn=0.; 
            for(k=i;k<=(i+1);k++) gGRMFAMPASyn=gGRMFAMPASyn+SGRMFAMPA[k][j];
            gGRMFAMPASyn=gGRMFAMPA*JGRMF*gGRMFAMPASyn;
            IGRMFAMPASyn=gGRMFAMPASyn*(vGR[i][j]-VGRMFSyn);

            gGRMFNMDASyn=0.; 
            for(k=i;k<=(i+1);k++) gGRMFNMDASyn=gGRMFNMDASyn+SGRMFNMDA[k][j];
            gGRMFNMDASyn=gGRMFNMDA*JGRMF*gGRMFNMDASyn;
            IGRMFNMDASyn=gGRMFNMDASyn*(vGR[i][j]-VGRMFSyn);

            gGRGOGABASyn=0.; 
            for(k=i;k<=(i+1);k++) {
                for(l=1;l<=NoPreGLGO[k];m++) {
                    gGRGOGABASyn=gGRGOGABASyn+SGRGOGABA[IPreGLGO[k][m]];
                }
            }
            gGRGOGABASyn=gGRGOGABA*JGRGO*gGRGOGABASyn;
            IGRGOGABASyn=gGRGOGABASyn*(vGR[i][j]-VGRGOSyn);

            fvGR[i][j]=(-ILGR-IAHPGR-IGRMFAMPASyn-IGRMFNMDASyn-IGRGOGABASyn)/CGR;  /* Note the (1/CGR) coefficient */
        }
    }

    /* Governing Eqs. for GO cells */
    for(i=1;i<=NC;i++) {
        ILGO=gLGO*(vGO[i]-VLGO);
        if(SPTIMEGO[i][1]>0.) IAHPGO=gAHPGO*exp(-(t-SPTIMEGO[i][1])/tauAHPGO)*(vGO[i]-VAHPGO);
        else IAHPGO=0.;

        gGOGRAMPASyn=0.; 
        for(j=1;j<=NoPreGOGR[i];j++) {
            gGOGRAMPASyn=gGOGRAMPASyn+SGOGRAMPA[IPreGOGR[i][j]][iPreGOGR[i][j]];
        }
        gGOGRAMPASyn=gGOGRAMPA*JGOGR*gGOGRAMPASyn;
        IGOGRAMPASyn=gGOGRAMPASyn*(vGO[i]-VGOGRSyn);

        gGOGRNMDASyn=0.; 
        for(j=1;j<=NoPreGOGR[i];j++) {
            gGOGRNMDASyn=gGOGRNMDASyn+SGOGRNMDA[IPreGOGR[i][j]][iPreGOGR[i][j]];
        }
        gGOGRNMDASyn=gGOGRNMDA*JGOGR*gGOGRNMDASyn;
        IGOGRNMDASyn=gGOGRNMDASyn*(vGO[i]-VGOGRSyn);

        fvGO[i]=(-ILGO-IAHPGO-IGOGRAMPASyn-IGOGRNMDASyn)/CGO;  /* Note the (1/CGO) coefficient */
    }
    
    /* Governing Eqs. for BCs */
    for(i=1;i<=NPC;i++) {
        ILBC=gLBC*(vBC[i]-VLBC);
        if(SPTIMEBC[i][1]>0.) IAHPBC=gAHPBC*exp(-(t-SPTIMEBC[i][1])/tauAHPBC)*(vBC[i]-VAHPBC);
        else IAHPBC=0.;

        gBCGRAMPASyn=0.; 
        for(j=(NC-143)+(i-1)*64;j<=(NC-143)+(i-1)*64+287;j++) {
            if(j<1) Idx=j+NC;
            if(j>NC) Idx=j-NC;
            for(k=1;k<=NGR;k++) {
                gBCGRAMPASyn=gBCGRAMPASyn+SBCGRAMPA[Idx][k];
            }
        }
        gBCGRAMPASyn=gBCGRAMPA*JBCGR*gBCGRAMPASyn;
        IBCGRAMPASyn=gBCGRAMPASyn*(vBC[i]-VBCGRSyn);

        fvBC[i]=(-ILBC-IAHPBC-IBCGRAMPASyn)/CBC;  /* Note the (1/CBC) coefficient */
    }
    
    /* Governing Eqs. for PCs */
    for(i=1;i<=NPC;i++) {
        ILPC=gLPC*(vPC[i]-VLPC);
        if(SPTIMEPC[i][1]>0.) IAHPPC=gAHPPC*exp(-(t-SPTIMEPC[i][1])/tauAHPPC)*(vPC[i]-VAHPPC);
        else IAHPPC=0.;

        gPCGRAMPASyn=0.; 
        for(j=(NC-143)+(i-1)*64;j<=(NC-143)+(i-1)*64+287;j++) {
            if(j<1) Idx=j+NC;
            if(j>NC) Idx=j-NC;
            for(k=1;k<=NGR;k++) {
                gPCGRAMPASyn=gPCGRAMPASyn+SPCGRAMPA[Idx][k];
            }
        }
        gPCGRAMPASyn=gPCGRAMPA*gPCGRAMPASyn;
        IPCGRAMPASyn=gPCGRAMPASyn*(vPC[i]-VPCGRSyn);

        gPCIOAMPASyn=gPCIOAMPA*JPCIO*SPCIOAMPA;
        IPCIOAMPASyn=gPCIOAMPASyn*(vPC[i]-VPCIOSyn);
    
        gPCBCGABASyn=0.;
        for(j=(i-1);j<=(i+1);j++) {
            if(j<1) Idx=j+NC;
            if(j>NC) Idx=j-NC;
            gPCBCGABASyn=gPCBCGABASyn+SPCBCGABA[Idx];
        }
        gPCBCGABASyn=gPCBCGABA*JPCBC*gPCBCGABASyn;
        IPCBCGABASyn=gPCBCGABASyn*(vPC[i]-VPCBCSyn);

        fvPC[i]=(-ILPC-IAHPPC+IdcPC-IPCGRAMPASyn-IPCIOAMPASyn-IPCBCGABASyn)/CPC;  /* Note the (1/CPC) coefficient */
    }

    /* Governing Eqs. for VN */
    ILVN=gLVN*(vVN-VLVN);
    if(SPTIMEVN[1]>0.) IAHPVN=gAHPVN*exp(-(t-SPTIMEVN[1])/tauAHPVN)*(vVN-VAHPVN);
    else IAHPVN=0.;

    gVNMFAMPASyn=0.; 
    for(i=1;i<=NMFVN;i++) gVNMFAMPASyn=gVNMFAMPASyn+SVNMFAMPA[i];
    gVNMFAMPASyn=gVNMFAMPA*JVNMF*gVNMFAMPASyn;
    IVNMFAMPASyn=gVNMFAMPASyn*(vVN-VVNMFSyn);

    gVNPCGABASyn=0.; 
    for(i=1;i<=NPC;i++) gVNPCGABASyn=gVNPCGABASyn+SVNPCGABA[i];
    gVNPCGABASyn=gVNPCGABA*JVNPC*gVNPCGABASyn;
    IVNPCGABASyn=gVNPCGABASyn*(vVN-VVNPCSyn);

    fvVN=(-ILVN-IAHPVN+IdcVN-IVNMFAMPASyn-IVNMFNMDASyn-IVNPCGABASyn)/CVN;  /* Note the (1/CVN) coefficient */

    /* Governing Eqs. for IO */
    ILIO=gLIO*(vIO-VLIO);
    if(SPTIMEIO[1]>0.) IAHPIO=gAHPIO*exp(-(t-SPTIMEIO[1])/tauAHPIO)*(vIO-VAHPIO);
    else IAHPIO=0.;

    gIOESAMPASyn=gIOESAMPA*JIOES*SIOESAMPA;
    IIOESAMPASyn=gIOESAMPASyn*(vIO-VIOESSyn);
    gIOVNGABASyn=gIOVNGABA*JIOVN*SIOVNGABA;
    IIOVNGABASyn=gIOVNGABASyn*(vIO-VIOVNSyn);

    fvIO=(-ILIO-IAHPIO-IIOESAMPASyn-IIOVNGABASyn)/CIO;  /* Note the (1/CIO) coefficient */
}

/* A Routine to find the cut-off time of the alpha function 
   for the local synapse by using the bisection method */
float CutOff1(float tau)
{
    float XT,XL,XC,XR,FTN,FTNC;

    XT=0.;
    do{
        XT=XT+0.01;
        FTN=ALPHA1(XT,tau);
    }while(FTN>ALPHAbound);

    /* For t=XL (or t=XR), the double exponential function > (or <) DEbound. */
    XL=XT-0.1; XR=XT;

    /* Find the cut-off time by using bisection method */
    do{
        XC=(XL+XR)/2.;
        FTNC=ALPHA1(XC,tau)-ALPHAbound;
        if(FTNC>0.) XL=XC;
        else XR=XC;
    }while((XR-XL)>CutOffTimebound);
    XC=XT;

    return XC;
}

/* A Routine to find the cut-off time of the alpha function 
   for the local synapse by using the bisection method */
float CutOff2(float A1, float A2, float tau1, float tau2)
{
    float XT,XL,XC,XR,FTN,FTNC;

    XT=0.;
    do{
        XT=XT+0.01;
        FTN=ALPHA2(XT,A1,A2,tau1,tau2);
    }while(FTN>ALPHAbound);

    /* For t=XL (or t=XR), the double exponential function > (or <) DEbound. */
    XL=XT-0.1; XR=XT;

    /* Find the cut-off time by using bisection method */
    do{
        XC=(XL+XR)/2.;
        FTNC=ALPHA2(XC,A1,A2,tau1,tau2)-ALPHAbound;
        if(FTNC>0.) XL=XC;
        else XR=XC;
    }while((XR-XL)>CutOffTimebound);
    XC=XT; 

    return XC;
}

/* A Routine to return the value of the alpha function for the local synaptic current */
float ALPHA1(float t, float tau)
{
    float alp1;

    alp1=exp(-t/tau);
    return alp1;
}

/* A Routine to return the value of the alpha function for the local synaptic current */
float ALPHA2(float t, float A1, float A2, float tau1, float tau2)
{
    float alp1;

    alp1=A1*exp(-t/tau1)+A2*exp(-t/tau2);
    return alp1;
}

/* A routine to generate a uniform random number on (0,1) 
   using the C-supplied rand function */
float UniRandom()
{   
    float R;

    R=((float) rand()+1.)/((float) RAND_MAX+2.);
    return R;
}

