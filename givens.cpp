////////////////////////////////////////////////
///////// GIVENS Rotation berechnen ////////////
////////////////////////////////////////////////
//////////////////i = 4 quits///////////////////


#include <iostream>
#include <cmath>

using namespace std;

int main(){
    double M[3][3],MatrixBegin[3][3],MatrixBeginTemp[3][3],cosine,sine;
    int i,j,ii,jj;
    
    //  Anfangscredits
    
    cout << endl << "////////////////////////////////////////////////" << endl << "///////// GIVENS Rotation berechnen ////////////" << endl << "////////////////////////////////////////////////"<< endl<<"//////////////////i = 4 quits///////////////////" << endl;

    
    //  Anfangsmatrix befüllen
    
    for(int i=1;i<=3;i++){
        for(int j=1;j<=3;j++){
            MatrixBegin[i-1][j-1]=j*i;
            MatrixBeginTemp[i-1][j-1]=j*i;
        }
    }
    
    
    //HIER MUSS SCHLEIFE BEGINNEN
    
    while(i!=4){
    
        //  i und j einlesen
    
        cout<<"i: ";
        cin>>ii;
        if (ii==4){return 88;}
        else{
        cout<<endl<<"j: ";
        cin>>jj;
        
        cosine=MatrixBegin[jj-1][jj-1]/sqrt(pow(MatrixBegin[jj-1][jj-1],2)+pow(MatrixBegin[ii-1][jj-1],2));
        sine=-MatrixBegin[ii-1][jj-1]/sqrt(pow(MatrixBegin[jj-1][jj-1],2)+pow(MatrixBegin[ii-1][jj-1],2));

        //  GivensMatrix befüllen
        
        for(int i=0;i<=2;i++){
            for(int j=0;j<=2;j++){
                if (i==j && (i==ii-1 || j==jj-1)){
                    M[i][j]=cosine;
                }
                else if (i==j && (i!=ii-1 || j!=jj-1)){
                    M[i][j]=1;
                }
                else if (i==ii-1 && j==jj-1){
                    M[i][j]=sine;
                }
                else if (j==ii-1 && i==jj-1){
                    M[i][j]=-sine;
                }
                else{
                    M[i][j]=0.;
                }
            }
        
        }
   
        cout << "i= " << ii << ",    " << "j= " << jj << endl;
    
        //  Anfangsmatrix ausgeben
    
        cout << "Anfangsmatrix: " << endl;
    
        for(i=0;i<3;i++){
            for(j=0;j<3;j++){
                cout << "  " << MatrixBegin[i][j]  <<  " ";
            }
            cout  << endl;
        }
    
        //Givens matrix ausgeben
        
        cout << "Givensmatrix: " << endl;
        
        for(i=0;i<3;i++){
            for(j=0;j<3;j++){
                cout << "  " << M[i][j]  <<  " ";
            }
        
            cout  << endl;
        
            if (i==2){
                cout << endl;}
        }
    
        // Schleife 1 Test:
    
        double matrixSchleife1[i][j];
    
        for(i=0; i<3; i++){
            for(j=0; j<3; j++){
                matrixSchleife1[i][j] = 0;
                for(int k=0; k<3; k++){
                    matrixSchleife1[i][j] += M[i][k] * MatrixBegin[k][j];
                }
                MatrixBegin[i][j]=matrixSchleife1[i][j];
            }
        }
    
        //1. Multiplizierte matrix ausgeben
    
        cout << "1. Multiplizierte matrix: " << endl;
    
        for(i=0;i<3;i++){
            for(j=0;j<3;j++){
                cout << "  " << matrixSchleife1[i][j]  <<  " ";
            }
        
            cout  << endl;
        
            if (i==2){
                cout << endl;
            }
        }
        }
    }


}