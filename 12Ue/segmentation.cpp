// MR image processing tool
// Daniel Spitzbart and Johannes Herbst
// June 2014
//
//

#include <sstream>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#define RangeX 254
#define RangeY 333


using namespace std;

int readfile(string filename, int mat[RangeX][RangeY]){
  int x, y, temp;
  fstream file;
  string completeline;
  file.open(filename, std::ios::in);
  if(file==1)cout<<"Reading file " << filename << endl;
  else {
    cout << "File " << filename << " missing. Abort.";
    return 0;
  }
  while(!file.eof()){
     stringstream sstr;
     getline(file,completeline);
     if(completeline.length()>0){
        sstr<<completeline;
        sstr>>x>>y>>temp;
        mat[x-1][y-1]=temp;
     }
  }
  file.close();
  return 1;
}

double calcError(int a[RangeX][RangeY], int b[RangeX][RangeY], double error[5]){
  int count[5], right[5],e=0;
  for (int i=0; i<5;i++){
    count[i]=right[i]=0;
    error[i]=0.;
  }
  for(int i=0; i<RangeX;i++){
     for(int j=0; j<RangeY;j++){
       right[b[i][j]-1]++;
       if(a[i][j]!=b[i][j])count[b[i][j]-1]++;
     }
   }
   for (int i=0; i<5;i++){
     error[i]=double(count[i])/right[i];
     e+=count[i];
   }
   return double(e)/(RangeX*RangeY);
}

void seg(int tis[RangeX][RangeY], int grey[RangeX][RangeY], int type[5], int std[5]){
  int minInt;
  fstream out;
  out.open("classification.dat",std::ios::out);
   for(int i=0; i<RangeX;i++){
      for(int j=0; j<RangeY;j++){
         double tempD,min=0;
         for(int k=0; k<5; k++){
            tempD=exp(-pow((grey[i][j]-type[k]),2)/(2*pow(std[k],2))+log(std[k]));
            if(tempD>min){
               min=tempD;
               minInt=k;
            }
         }
         tis[i][j]=minInt+1;
         out<<i+1<<"\t"<<j+1<< "\t" << minInt+1 << endl;
         minInt=0;
      }
      out << endl;
   }
   out.close();
}

void SA(double Tbegin, double Tfinal, double lambda, int tissue[RangeX][RangeY], int grey[RangeX][RangeY], int type[5], int std[5]){
  int sweeps=1000;
  int Inew, Jnew,rand1,rand2,rand3,temp;
  int Inn[4] = {1,-1,0,0};       /* Nearest neighbor array I */
  int Jnn[4] = {0,0,1,-1};       /* Nearest neighbor array J */

  double DeltaE, coupling=1.1,r,z;

  cout << endl << "Simulated annealing with " << sweeps << " sweeps per temp" << endl;
  cout << "Starting temp: " << Tbegin << endl;
  cout << "End temp: " << Tfinal << endl;
  cout << "Stepsize: " << lambda << endl << endl;

  for (double T=Tbegin;T>=Tfinal;T/=lambda){
    cout << "Processing temperature " << T << endl;
    for(int a=0;a<sweeps;a++){
      for(int i=0;i<RangeX*RangeY;i++){
        rand1=int(drand48()*RangeX); //find random position to flip
        rand2=int(drand48()*RangeY);
        rand3=tissue[rand1][rand2];
        while(tissue[rand1][rand2]==rand3){ //create a random "spin" that isn't equal to the actual spin
          rand3=int(drand48()*5)+1;
        }
        temp=tissue[rand1][rand2];
        tissue[rand1][rand2]=rand3;

        DeltaE=0.;
        /* Loop over nearest neighbors */
        for(int k=0;k<4;k++){
          Inew = rand1 + Inn[k];
          Jnew = rand2 + Jnn[k];
          /* Check periodic boundary conditions */
          if(Inew < 0){
            Inew = RangeX-1;
          }
          else if(Inew >= RangeX){
            Inew = 0;
          }
          if(Jnew < 0){
            Jnew = RangeY-1;
          }
          else if(Jnew >= RangeY){
            Jnew = 0;
          }
          if(tissue[rand1][rand2]==tissue[Inew][Jnew]) DeltaE += -coupling;// * tis[rand1][rand2] * tis[Inew][Jnew]; // only nearest neighbors important
          else DeltaE+=coupling;
        }
        /*Calculate the contribution from the field H */

        DeltaE += pow((grey[rand1][rand2]-type[tissue[rand1][rand2]-1]),2)/(2*pow(std[tissue[rand1][rand2]-1],2))+log(std[tissue[rand1][rand2]-1]);

        z=drand48();
        if ((DeltaE/T)>20.) r=0;
        else{
           if ((DeltaE/T)<0.05) r=1.01;
           else r=exp(-DeltaE/T);
        }
        if(z<r || r>=1){ //check flipping criteria
          DeltaE=0.;
        }
        else{ //flip back
          tissue[rand1][rand2]=temp;
        }
      }
    }
  }
}

int main(){

  int type[5]={30,426,602,1223,167};
  int std[5]={30,59,102,307,69};
  int grey[RangeX][RangeY], tis[RangeX][RangeY], cor[RangeX][RangeY];

  double error[5];

  fstream image,correct,out2;

  if(readfile("SimMRimage.dat", grey)==0)return 0;
  if(readfile("CorrectSegImage.dat", cor)==0)return 0;

  cout << "Starting segmentation" << endl;

  seg(tis, grey, type, std); //segmentation

  cout << endl << "Error after segmentation: " << calcError(tis,cor,error)*100 << "%" << endl;
  for (int i=0; i<5;i++)cout << "Tissue #" << i+1 << ": "<< error[i]*100 << "%" << endl;

  SA(1.,.5,1.1,tis,grey, type, std); // first SA

  cout << "Error after SA: " << calcError(tis,cor,error)*100 << "%" << endl;
  for (int i=0; i<5;i++)cout << "Tissue #" << i+1 << ": "<< error[i]*100 << "%" << endl;

  SA(0.8,0.5,1.1,tis,grey,type,std); // re-heating and cooling

  cout << "Error after SA: " << calcError(tis,cor,error)*100 << "%" << endl;
  for (int i=0; i<5;i++)cout << "Tissue #" << i+1 << ": "<< error[i]*100 << "%" << endl;

  out2.open("saSweep.dat",std::ios::out);
  for(int i=0;i<RangeX;i++){
    for(int j=0;j<RangeY;j++){
      out2 << i+1<< "\t" << j+1 << "\t" << tis[i][j] << endl;
    }
    out2<<endl;
  }

  cout << "Segmentation good :)" << endl;
  return 1;
}
