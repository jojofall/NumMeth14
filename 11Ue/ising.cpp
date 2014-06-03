#define PROGNAME "mc_Ising_2D"
#define VERSION "0.2a"
#define DATE  "today"
#define AUTHOR "DSJH"

/* Metropolis Monte Carlo for 2D-Ising Model with periodic boundary conditions */
/* Code adapted from http://gold.cchem.berkeley.edu/~acpan/220A/2DIsing.c and N. Bluemer*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

/*************************
 * Constant Declarations *
 ************************/
#define NMAX 5           /* Lattice size (i.e. NMAX x NMAX) */

#define J 1.0             /* Ising lattice coupling constant
			    J > 0: Ferromagnetic 
			    J < 0: Antiferromagnetic */
#define H 0.0              /* Ising lattice field strength */

void error(char error_text[])
/* standard error handler */
{
    fprintf(stderr,"\n%s run-time error\n",PROGNAME);
    fprintf(stderr,"--%s--\n",error_text);
    fprintf(stderr,"for general help use option -h\n");
    fprintf(stderr,"...now exiting to system...\n");
    exit(1);
}

void random_init(int mag[NMAX][NMAX], int L)
{
  int i,j;
  
  for(i = 0; i < L; i++){
    for(j = 0; j < L; j++)
      {
          if(drand48() > 0.5)
              mag[i][j] = 1;
          else
              mag[i][j] = -1;
       }
    }
  }



void outputmag(int mag[NMAX][NMAX], int L)
{
  int i,j;

  printf("\n");

  for(i = 0; i < L; i++)
    {
      for(j = 0; j < L; j++)
	{
	  if(mag[i][j] == -1)
	    printf("O ");
	  else
	    printf("X ");
	}
      printf("\n");
    }
}

double mag_av(int mag[NMAX][NMAX], int L)
/* computes average magnetization */
{
  int i,j;
  double x;

  x=0.0;

	for(i = 0; i < L; i++){
		for(j = 0; j < L; j++){
			x+=1.0*mag[i][j];
  		}
	}
	x=sqrt(x*x)/(L*L);		

	return(x);
}

double calcEnergy(int mag[NMAX][NMAX], int L){ //Function for calculating the Energy
    double Energy=0.;
    int Inn[4] = {1,-1,0,0};       /* Nearest neighbor array I */
    int Jnn[4] = {0,0,1,-1};       /* Nearest neighbor array J */
    int Inew, Jnew;                /* Nearest neighbot indices */
    for(int i=0;i<L;i++){
        for(int j=0;j<L;j++){
            for(int k=0;k<4;k++){
                Inew = i + Inn[k];
                Jnew = j + Jnn[k];
                if(Inew < 0){ // Periodic boundary
                    Inew = L-1;
                }
                else if(Inew >= L){
                    Inew = 0;
                }
                if(Jnew < 0){
                    Jnew = L-1;
                }
                else if(Jnew >= L){
                    Jnew = 0;
                }
                Energy += -J * mag[i][j] * mag[Inew][Jnew];
            }
            Energy += 2*H*mag[i][j]; //Field contribution
        }
    }
    Energy /= 2.0;

    return Energy;
}

void flipAndCheck(int mag[NMAX][NMAX], int L, double Energy, double Td){
    int xi, yi;
    int Inn[4] = {1,-1,0,0};       /* Nearest neighbor array I */
    int Jnn[4] = {0,0,1,-1};       /* Nearest neighbor array J */
    int Inew, Jnew;                /* Nearest neighbot indices */
    double DeltaE=0.,r,z;
    for (int b=0;b<L*L;b++){
        xi=int(drand48()*L);
        yi=int(drand48()*L);
        mag[xi][yi]=-mag[xi][yi];
        DeltaE=0.;
        /* Loop over nearest neighbors */
        for(int k=0;k<4;k++){
            Inew = xi + Inn[k];
            Jnew = yi + Jnn[k];
            /* Check periodic boundary conditions */
            if(Inew < 0){
                Inew = L-1;}
            else if(Inew >= L){
                Inew = 0;}
            if(Jnew < 0){
                Jnew = L-1;}
            else if(Jnew >= L){
                Jnew = 0;}
            DeltaE += -2*J * mag[xi][yi] * mag[Inew][Jnew]; // only nearest neighbors important
        }
        /*Calculate the contribution from the field H */
        DeltaE += 2*H*mag[xi][yi];
        r=exp(-DeltaE/Td);
        z=drand48();
        
        if(z<r || r>1){ //check flipping criteria
            Energy-=DeltaE;
        }
        else{ //flip back
            mag[xi][yi]=-mag[xi][yi];
        }
	}
    
}


void printhelp ()
{
 printf("**********************************************************\n");
 printf("%s: Metropolis Monte Carlo for 2D Ising model\n",PROGNAME);
 printf("Version: %s, %s by %s\n",VERSION,DATE,AUTHOR);
 printf("based on http://gold.cchem.berkeley.edu/~acpan/220A/2DIsing.c\n");
 printf("options: -T temperature (really: k_B T/|J|)\n");
 printf("         -L linear dimension of lattice (L<=%d)\n",NMAX);
 printf("         -w# number of initialization sweeps\n");
 printf("         -n# number of measurement sweeps\n");
 printf("         -c print configurations\n");
 printf("         -h this help\n");
}

/***********************************************************/
int main(int argc, char *argv[])
{
  /*************************
   * Variable Declarations *
   ************************/
  int mag[NMAX][NMAX],magS[NMAX][NMAX];           /* 2D Ising Lattice */
  int i, j, k;                   /* Loop counters */
  int s,d;                       /* Lattice spin variables */  
  double Energy;                 /* Total lattice energy */
  int Inn[4] = {1,-1,0,0};       /* Nearest neighbor array I */
  int Jnn[4] = {0,0,1,-1};       /* Nearest neighbor array J */
  int Inew, Jnew;                /* Nearest neighbot indices */ 
  double Etemp, deltaE;          /* Temp energy variables for MC moves */ 
  int accept = 0;                /* Number of accepted moves */
  int move = 0;                  /* Number of moves total */ 

  char c;
  double T;                      /* temperature (in units of J/k_B) */
  int sweeps;                    /* number of measurement sweeps */
  int warm;                      /* number of warm-up sweeps */
  int print_conf;                /* flag for printing configurations */
  int L;                         /* lattice dimension */
    double magAv=0., eAv=0.,magAvSq,magAvCu;
    vector<double> magV,eV,magVSq,magVCu;

  /***************************
   * Initialization          *
   ***************************/
  //T=10;
    sweeps=10000;
    warm=1000;
    print_conf=0;
    L=NMAX;


  while (--argc > 0 && (*++argv)[0] == '-'){
    while (c== *++argv[0]){
	   switch (c) {
	   case 'n':
 	     sscanf(++argv[0],"%d\n",&sweeps); 
             break;
	   case 'w':
 	     sscanf(++argv[0],"%d\n",&warm); 
             break;
	   case 'L':
 	     sscanf(++argv[0],"%d\n",&L); 
             //if (L>NMAX) error ("L too large");
             break;
           case 'T':
 	     sscanf(++argv[0],"%lf\n",&T); 
             break;
           case 'c':
             print_conf=1;
             break;
	   case 'h':
	     printhelp();
	     exit(0);
             break;
/*  	   default:  */
/*  	     error("No valid choice");  */
	   }
}}

    fstream out;
    out.open("Wert.dat",std::fstream::out);

    /* Seed the random number generator */
    srand48((unsigned int) time(NULL) - (time(NULL)/100000)*100000);
    out<< "T \t  magAv \t stdmagAv \t m^2 \t stdm^2 \t m^4 \t stdm^4 \t eAv\t stdeAv\t Binder \t std Binder"<< endl;

    double Td,r,DeltaE,z;
    double stdE,stdM,stdM2,stdM4,binder,stdBinder;
//    int xi, yi;
    
    for(T=1;T<50;T++){ //Loop over temperatures
        Td=T/10.;

        random_init(mag, L);
        Energy=calcEnergy(mag, L);
 
        for(int a=0;a<warm;a++){ //warm-up runs
            flipAndCheck(mag,L,Energy,Td);
        }
        
        magAv=magAvSq=magAvCu=eAv=0.; // empty averages and vectors
        magV.clear();
        magVSq.clear();
        magVCu.clear();
        eV.clear();

        for(int a=0;a<sweeps;a++){ // run!
            flipAndCheck(mag,L,Energy,Td);
            magAv+=mag_av(mag,L);
            magAvSq+=pow(mag_av(mag,L),2);
            magAvCu+=pow(mag_av(mag,L),4);
            Energy=calcEnergy(mag,L);
      
            eAv+=Energy/(L*L);
            magV.push_back(mag_av(mag,L));
            magVSq.push_back(pow(mag_av(mag,L),2));
            magVCu.push_back(pow(mag_av(mag,L),4));
            eV.push_back(Energy/(L*L));
        }
        

        eAv=eAv/sweeps;
        magAv=magAv/sweeps;
        magAvSq=magAvSq/sweeps;
        magAvCu=magAvCu/sweeps;
        binder=1.-(magAvCu/(3*(pow(magAvSq,2))));
        
        //Standarderror-calculation:
        stdE=stdM=stdM2=stdM4=0.;
        for(int i=0;i<sweeps;i++){
            stdE+=pow(eV[i]-eAv,2);
            stdM+=pow(magV[i]-magAv,2);
            stdM2+=pow(magVSq[i]-magAvSq,2);
            stdM4+=pow(magVCu[i]-magAvCu,2);
        }
      
        stdE=sqrt(stdE/(sweeps-1));
        stdM=sqrt(stdM/(sweeps-1));
        stdM2=sqrt(stdM2/(sweeps-1));
        stdM4=sqrt(stdM4/(sweeps-1));
        stdBinder=sqrt(pow(-1./(3*pow(magAvSq,2)),2))*stdM2+sqrt(pow(2.*magAvSq/(3*pow(magAvCu,3)),2))*stdM4; // propagation of uncertainty for M^2 and M^4
      
        cout<<scientific;
        cout<<"T= "<<Td<<"\tMagAv= "<<magAv<<"\teAv= "<<eAv<<endl;
        out<< Td << "\t" << magAv << "\t" << stdM << "\t" << magAvSq << "\t" << stdM2 << "\t" << magAvCu << "\t" << stdM4 << "\t" << eAv << "\t" << stdE << "\t" << binder << "\t" << stdBinder << endl;
        
  
    }

    return 0;
}

