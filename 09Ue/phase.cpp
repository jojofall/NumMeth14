#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

double vPot(double x,double vZero,double rZero,double aZero,double vOne,double rOne,double aOne) {
  double vPott;
  vPott=-vZero/(1.+exp((x-rZero)/aZero));
    vPott=vPott+vOne*exp(-pow((x-rOne),2)/aOne);
  return vPott;
}

double wr(double r,double energy,int l, double Vpotential,double mhb2){
  double wrr;
  wrr=mhb2*energy-l*(l+1)/pow(r,2)-mhb2*Vpotential;
  return wrr;
}

double getU(double h, double w, double z){
  return (1+pow(h,2)*w/12.)*z;
}

double getZ(double h, double w, double u){
  return u/(1+pow(h,2)*w/12.);
}

double getNextZ(double z1, double z0, double h, double w){
  return 2*z1-z0-pow(h,2)*w*z1;
}

double Neumann(double x, int l){
  if(l==0) return cos(x);
  else if (l==1) return cos(x)/x+sin(x);
  else return -1.;
}

double Bessel(double x, int l){
  if(l==0) return sin(x);
  else if (l==1) return sin(x)/x-cos(x);
  else return -1.;
}

double getDelta(double u1, double u2, double k, double r1, double r2, int l){
  return atan((u1*Bessel(k*r2,l)-u2*Bessel(k*r1,l))/(u1*Neumann(k*r2,l)-u2*Neumann(k*r1,l)));
}

int main(){

  double v0,r0,A,a,v1,r1,a1,pi,mRed,mn=939.56563,mA,amu=931.49432,mhb2,hbarc=197.327053,dr,Rmax,rr, Energy, dE=0.1;
  int l,Nrmax=2000,NEmax=2000;

  fstream in,outPot;
  in.open("Werte.dat",std::fstream::in);
  in>>v0>>rr>>A>>a>>v1>>r1>>a1>>l>>dr>>Rmax>>mA;
  in.close();

  mRed=mn*mA*amu/(mA*amu+mn);
   mhb2=2.*mRed/pow(hbarc,2);
   r0=rr*pow(A,(1./3.));

 outPot.open("Potential.dat",std::fstream::out);
 int mxf,n=0;
 mxf=int(Rmax/dr);

  for(int f=0;f<=mxf ;f++){
   outPot<<f*dr<< " " <<vPot(f*dr,v0,r0,a,v1,r1,a1)<<endl;
 }

   outPot.close();

  vector<double> z,u,delta;
  double vnew,k, singleDelta;

  pi=2.*asin(1.);
  Energy=200.;
  ofstream fout;
  fout.open("test.dat");
  //cout << getZ(dr,0.,0.1) << endl;
  for(int i=0; i<NEmax; i++){
    z.push_back(getZ(dr,0.,0.));
    z.push_back(getZ(dr,0.,0.1));
    u.push_back(0.);
    u.push_back(0.1);
    k=sqrt(mhb2*Energy);
    for(int j=2; j<Nrmax; j++){
      vnew=vPot(j*dr,v0,r0,a,v1,r1,a1);
      z.push_back(getNextZ(z[j-1],z[j-2],dr,wr(j*dr,Energy,l,vnew,mhb2)));
      u.push_back(getU(dr,wr(j*dr,Energy,l,vnew,mhb2),getNextZ(z[j-1],z[j-2],dr,wr(j*dr,Energy,l,vnew,mhb2))));
    }

    singleDelta=getDelta(u[1980],u[1981],k,1980*dr,1981*dr,l);
    delta.push_back(-singleDelta);
    cout << i-1 << " " << delta[i-1] << " " << i << " " << delta[i] << endl;
    if((delta[i-1]-delta[i])>0.7*pi){
      n++;
    }
    fout << Energy << " " << -singleDelta+n*pi << endl;// << " " << u[1980] << " " << u[1981] << " " << 1980*dr << " " << 1981*dr << endl;

    Energy-=dE;
    z.clear();
    u.clear();
  }

  fout.close();


}
