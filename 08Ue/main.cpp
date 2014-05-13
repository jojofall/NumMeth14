#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "gnuplot_i.hpp"

using namespace std;

double z2(double x0, double x1, double x2, double eta, double w, double f, double omega){
       return -2*eta*x2-pow(w,2)*sin(x1)+f*sin(omega*x0);
}
       
int main(){
  const double gravity=9.80665;
  double length, mass, thetaZero,etha,f,Force,Omega,omegaZero,theta,thetaDot,t;
  double k1,k2,k3,k4,Phi,h,runtime,Pi;
	vector<double> thetaDotVec,thetaVec,tVec; 
	cout<<"laufzeit variable t[s]: ";
	cin>>runtime;
        
/*     Daten einlesen       */ 
  std::fstream in,out;
  in.open("Werte.dat", std::fstream::in);
  in>>length>>mass>>thetaZero>>etha>>Force>>Omega;
  in.close();

	cout<<"Eingelesene Werte: "<<endl<<"Länge: "<<length<<endl<<"Masse: "<<mass<<endl<<"Theta(t=0): "<<thetaZero<<endl<<"Etha: "<<etha<<endl<<"Kraft: "<<Force<<endl<<"Omega(Erreger): "<<Omega<<endl<<endl;

/*     werte                */
  h=0.01;
  Pi=asin(1.)*2.;
  omegaZero=sqrt(gravity/length);
	cout << "OmegaNull: " << omegaZero << endl;
  f=Force/(length*mass);
  theta=thetaZero;
  thetaDot=0;
  t=0;

	out.open("Ausgaben.dat",std::fstream::out);	

	for(int i=0;i<runtime/h;i++){
    k1=z2(t,theta,thetaDot,etha,omegaZero,f,Omega);
	  k2=z2(t+h/2,theta+h/2,thetaDot+h*k1/2,etha,omegaZero,f,Omega);
    k3=z2(t+h/2,theta+h/2,thetaDot+h*k2/2,etha,omegaZero,f,Omega);
    k4=z2(t+h,theta+h,thetaDot+h*k3,etha,omegaZero,f,Omega);
    Phi=1./6*(k1+2*k2+2*k3+k4);		
	  thetaDot=thetaDot+h*Phi;
		thetaDotVec.push_back(thetaDot);

		k1=thetaDot;
		k2=thetaDot+h*k1/2;
		k3=thetaDot+h*k2/2;
		k4=thetaDot+h*k3;	
		
	  Phi=1./6*(k1+2*k2+2*k3+k4);	
		theta=theta+h*Phi;
		t=t+h;
		thetaVec.push_back(theta);
		tVec.push_back(t);
		
		out<<t<<" "<<theta<<" "<<thetaDot<<endl;
	}
	
	out.close();
  cout << "Lösung der Differentialgleichung abgeschlossen. Fouriertransformation..." << endl << endl;
  
	//Gnuplot g1 = Gnuplot(tVec, thetaVec, "Schwingung", "lines", "x-Achse", "y-Achse");
	//Gnuplot g2 = Gnuplot(thetaVec, thetaDotVec, "Phasespace", "lines", "Theta", "ThetaDot");
	fstream outFourier;
	outFourier.open("Fourier.dat",std::fstream::out);
	vector<double>::iterator thetaIt=thetaVec.begin();
	vector<double>::iterator tIt=tVec.begin();
	vector<double> omegaVec, fourierVec;
	double FT=0.,om;
	for (int i=1;i<10000;i++){
		om=i*0.01;
		omegaVec.push_back(om);
		thetaIt=thetaVec.begin();
		tIt=tVec.begin();
		FT=-cos(om*(*tIt))*(*thetaIt);
		while(thetaIt!=thetaVec.end()){
			FT=FT+2*cos(om*(*tIt))*(*thetaIt);
			thetaIt++;
			tIt++;
		}
		thetaIt--;
		tIt--;
		FT=FT-cos(om*(*tIt))*(*thetaIt);
		FT=h*FT/(2*3.1415);
		fourierVec.push_back(FT*FT);
		outFourier<<om<<" "<<FT<<" "<<FT*FT<<endl;
	}
	//Gnuplot g3 = Gnuplot("lines");
	//g3.set_xrange(0.,25.);
	//g3.plot_xy(omegaVec,fourierVec,"F(omega)**2");
	
  cout << "Fouriertransformation abgeschlossen" << endl << endl;

  //Berechnungen für Weyl
  //

  vector<double>::iterator thetaDotIt=thetaDotVec.begin();
  thetaIt=thetaVec.begin();
  double thetaMax=0.,thetaDotMax=0.;
  while(thetaIt!=thetaVec.end()){
    if((*thetaIt)>thetaMax)thetaMax=(*thetaIt);
    if((*thetaDotIt)>thetaDotMax)thetaDotMax=(*thetaDotIt);
    thetaIt++;
    thetaDotIt++;
  }
  
  double hP=0.0002;//Wert des Wirkungsquantums frei gewählt
  theta=acos(1.-hP*omegaZero/(4*Pi*mass*gravity*length));
  
  cout << "Theta(t=0) für Grundzustand: " << theta << endl;
  cout << "qMax: " << thetaMax*length << endl;
  cout << "pMax: " << thetaDotMax*length*mass << endl;
  cout << "Phasenraumvolumen (pMax*qMax*Pi): " << thetaMax*thetaDotMax*Pi*mass*length*length << endl << "Zum Vergleich, h/2 in unserem Fall: " << hP/2. << endl;
  
  cout<<"Beendet. Daten liegen in 'Ausgaben.dat' & 'Fourier.dat'"<<endl;
       	int i;
	cin>>i; //damit gnuplots nicht gleich wieder geschlossen werden
	return 0;
}
