//GIVENS

#include <iostream>
#include <cmath>

using namespace std;

int main(){
    double M[3][3],MatrixBegin[3][3],cosine,sine;
    int i,j,ii,jj;
    
    cout<<"geben sie i ein"<<endl;
    cin>>ii;
    cout<<"geben sie j ein"<<endl;
    cin>>jj;
    
    for(int i=1;i<=3;i++){
        for(int j=1;j<=3;j++){
            MatrixBegin[i-1][j-1]=j*i;
        }
    }
 
    for(int i=0;i<=2;i++){
        for(int j=0;j<=2;j++){
            if(i==j){
                M[i][j]=1.;
            }
            else{
                M[i][j]=0.;
            }
        }
        
    }

    cosine=M[jj][jj]/sqrt(pow(M[jj][jj],2)+pow(M[ii][jj],2));
    sine=-M[ii][jj]/sqrt(pow(M[jj][jj],2)+pow(M[ii][jj],2));
    
    
    
    return 0;
}

void dot(double matrix1[3][3],double matrix2[3][3],double matrix3[3][3]){
    for(int g=0;g<=2;g++){
        for(int h=0;h<=2;g++){
            matrix3[g][h]=matrix1[g][h]*matrix2[h][g]+matrix3[g][h];
        }
    }
}
