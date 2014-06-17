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
#include <ctime>


using namespace std;

double fun(double x){
  return sin(x)*cos(x)/pow(x,2);
}

double dist(double x){
  return (10./9.)/pow(x,2);
}

double funReduced(double x){
  return fun(x)/dist(x);
}

int main(){
  double sum=0.,sum2=0.,ran1,ran2,value;
  int counter=0,max=200000,sweeps=200;
  srand48(time(NULL));
  for(int j=0;j<sweeps;j++){
    for(int i=0;i<max;i++){
      ran1=drand48()*9.+1.;
      ran2=drand48();
      if(ran2<dist(ran1)){
        sum+=funReduced(ran1);
        counter++;
      }
      else if(dist(ran1)>1){
        sum+=funReduced(ran1);
        counter++;
      }
    }
    value+=sum/counter;
    sum=0.;
    counter=0;
  }
  cout << endl << "Integrated value: " << value/sweeps << endl;

  return 1;
}
