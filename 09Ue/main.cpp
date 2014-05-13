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



int main(){

	double v0,r0,A,a,v1,r1,a1,pi,mRed,mn=939.56563,mA,amu=931.49432,mhb2,hbarc=197.327053,dr,Rmax,Nr,rr;
	int l,Nrmax=2000,NEmax=2000;

	fstream in,outPot,outU;	
	in.open("Werte.dat",std::fstream::in);	
	in>>v0>>rr>>A>>a>>v1>>r1>>a1>>l>>dr>>Rmax>>mA;
	in.close();
	
	pi=2.*asin(1.);

 	mRed=mn*mA*amu/(mA*amu+mn);
  	mhb2=2.*mRed/pow(hbarc,2);
  	r0=rr*pow(A,(1./3.));

	outPot.open("Potential.dat",std::fstream::out);
	int mxf,len,i;
	mxf=int(Rmax/dr);

	for(int f=0;f<=mxf ;f++){
		outPot<<f*dr<< " " <<vPot(f*dr,v0,r0,a,v1,r1,a1)<<endl;
	}

	outPot.close();
	cout << "Shooting..." << endl;
	vector<double> z,u;
	double vnew,znew,stepsize,upast=1.;
	double E=0.,ud,norm=0.;
	int f=0,hitter=1,totalsteps=0;
	outU.open("U.dat",std::fstream::out);
  vector<double>::iterator it;
  stepsize=1.;
	while (sqrt(upast*upast)>0.000001){
		norm=0.;
		z.push_back((2.*exp(-sqrt(mhb2*-E))*(Rmax))/((1+(pow(dr,2)/12.)*wr(Rmax,E,l,vnew,mhb2))));
		u.push_back(2.*exp(-sqrt(mhb2*-E))*(Rmax));
		z.push_back((2.*exp(-sqrt(mhb2*(-E))*(Rmax-dr)))/(1+(pow(dr,2)/12.)*wr(Rmax-dr,E,l,vnew,mhb2)));
		u.push_back(2.*exp(-sqrt(mhb2*-E))*(Rmax-dr));
		norm=(u[0]+u[1])*dr/2.;
		
	
		for(int i=1;i<=mxf;i++){
			vnew=vPot(Rmax-i*dr,v0,r0,a,v1,r1,a1);
			znew=2*z[i]-z[i-1]-pow(dr,2)*wr(Rmax-i*dr,E,l,vnew,mhb2)*z[i];
			ud=(1+(pow(dr,2)/12.)*wr(Rmax-i*dr,E,l,vnew,mhb2))*znew;
			z.push_back(znew);		
			u.push_back(ud);
			norm=norm+(u[i]+u[i-1])*dr/2.;

		}
		norm=sqrt(norm*norm);

		for(int i=0;i<=mxf;i++){
			u[i]=-u[i]/norm;
		}
    
    
    //Alternativer Algorithmus zur Nullstellensuche (Sekantenverfahren), konvergiert aber bei l=1 nicht daher ->
    //x=-u[mxf]*dr/(u[mxf-1]-u[mxf]);
    //if(x*x>xpast*xpast&&f>0){
    //  stepsize=-stepsize/3.;
    //}
    //xpast=x;
		
    //Dieser Algorithmus
    if(u[mxf]*u[mxf]>upast*upast&&f>hitter){
      stepsize=-stepsize/3.;
    }
    upast=u[mxf];
    
		cout << ".";
    //cout << f << " " << E << " " << u[mxf] << endl;
		f++;
    totalsteps++;
    E=E-stepsize;
    if(E>0){ //Da der Algorithmus nicht immer gleich konvergiert, muss diese Vorkehrung getroffen werden (verkleinerung der Schritte und Umkehr erst nach einer Mindest-Schrittanzahl)
      hitter++;
      stepsize=1.;
      upast=1.;
      E=0.;
      f=0;
    }
		z.clear();
		u.clear();
	}
  E=E+stepsize;
  cout << endl << endl << "Energieeigenwert: " << E << endl << "Konvergiert nach " << totalsteps << " Schritten" << endl;
  for(int k=0;k<=mxf;k++){
    u[i]=u[i]/norm;
    outU<<Rmax-k*dr<<" "<<u[k]<<endl;
    }

  outU.close();
	cout<<"><><><><"<<endl;
	return 0;
}
