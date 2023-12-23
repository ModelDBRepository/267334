/* A code to get the conjunction index */

#include <stdio.h>      /* Standard header file for the Input and Output */
#include <stdlib.h>     /* Standard library header file */
#include <math.h>       /* Math header file */
#include <conio.h>      /* Header file for the IO on the console */
#include <time.h>       /* Header file for initializing the random number generator */
                   
#define NOD      2000   /* No. of data in the input file */

char FINN1[20],FINN2[20],FOUTN[20];
FILE *FIN1,*FIN2,*FOUT;

/* Cerebellar network in Yamazaki & Nagao (2012) */

void main()
{
    int i;
    float RG[NOD+1],RGI[NOD+1],dRG[NOD+1],dRGI[NOD+1];
    float AVGRG,AVGRGI,COV,SDRG,SDRGI,CI,tmpt,tmpf;

    /* Input names of output files, and open them */
    printf("\n INPUT THE NAME OF THE INPUT FILE 1 (R_G) \n");scanf("%s",FINN1);
    if((FIN1 = fopen(FINN1,"r"))==NULL){printf("FILE OPEN ERROR...\n");exit(-1);}
    printf("\n INPUT THE NAME OF THE INPUT FILE 2 (R_G^I) \n");scanf("%s",FINN2);
    if((FIN2 = fopen(FINN2,"r"))==NULL){printf("FILE OPEN ERROR...\n");exit(-1);}
    printf("\n INPUT THE NAME OF THE OUTPUT FILE \n");scanf("%s",FOUTN);
    if((FOUT = fopen(FOUTN,"w"))==NULL){printf("FILE OPEN ERROR...\n");exit(-1);}

    for(i=1;i<=NOD;i++) {
        fscanf(FIN1,"%f %f\n",&tmpt,&tmpf);
        RG[i]=tmpf; 
        fscanf(FIN2,"%f %f\n",&tmpt,&tmpf);
        RGI[i]=tmpf; 
    }

    AVGRG=0.; AVGRGI=0.;
    for(i=1;i<=NOD;i++) {AVGRG=AVGRG+RG[i]; AVGRGI=AVGRGI+RGI[i];}
    AVGRG=AVGRG/NOD; AVGRGI=AVGRGI/NOD;
    
    for(i=1;i<=NOD;i++) {dRG[i]=RG[i]-AVGRG; dRGI[i]=RGI[i]-AVGRGI;}
    
    COV=0.; SDRG=0.; SDRGI=0.;
    for(i=1;i<=NOD;i++) {
        COV=COV+dRG[i]*dRGI[i];
        SDRG=SDRG+dRG[i]*dRG[i];
        SDRGI=SDRGI+dRGI[i]*dRGI[i];
    }
    COV=COV/NOD; SDRG=sqrt(SDRG/NOD); SDRGI=sqrt(SDRGI/NOD);
    CI=COV/(SDRG*SDRGI);
    
    printf("Conjunction index = %12.8f",CI);
    fprintf(FOUT,"%12.8f\n",CI);

    fclose(FIN1);fclose(FIN2);fclose(FOUT);
    printf("<***** End *****>\n");
    getche();
}
