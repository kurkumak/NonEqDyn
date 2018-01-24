#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <cassert>
using namespace std;
//#include "comandi.h"


#define p1 3
#define p2 4

FILE * f11;
FILE * f12;
FILE * f13;
FILE * f14;


vector<double> mu;
vector<vector<double> >C, R, Q, F;

char type;
double eps=0.;
double ciccio=1;
double h;
//double Td=.805166;
double Tp=.805;
double T;
double q;
int n1,ip;
double Bp;
//double B=1/T;
//double Tp=1/Bp;

double R0=1.;



double mypow(double x,int p){
  double ff;
  ff=1;
  for(int i=0;i<p;i++)ff*=x;
  return ff;
}



double f(double q){
  double ff;
  //  ff=.5*q*q*(ciccio+eps*q*q);
  ff=.5*(ciccio*mypow(q,p1)+eps*mypow(q,p2));
  return ff;
}

double f1(double q){
  double ff;
  //    ff=q*(ciccio+eps*2.*q*q);
  ff=.5*(ciccio*p1*mypow(q,p1-1)+eps*p2*mypow(q,p2-1));
  return ff;
}

double f2(double q){
  double ff;
  //   ff=(ciccio+eps*6.*q*q);
  ff=.5*(ciccio*p1*(p1-1)*mypow(q,p1-2)+eps*p2*(p2-1)*mypow(q,p2-2));
  return ff;
}

double If2RR(int i, int j){
  double I=0;
  int k;
  for ( k=j+1;k<i;k++){
    I+=f2( C[i][k] )*R[i][k]*R[j][k];
////      I+=1./T*(f1(C[i][k])-f1( C[i][k]-T*R[i][k]*h ))*R[j][k];
  }
  k=j;
  //  I+=.5*f2( C[i][k] )*R[i][k]*R[j][k];
  k=i;
  //  I+=.5*f2( C[i][k] )*R[i][k]*R[j][k];
  I*=h;
  return I;
}


double If2RC(int i, int j){
    double I=0;
    int k;
  for ( k=0;k<i;k++){
      I+=f2( C[i][k] )*R[i][k]*C[j][k];
////      I+=1./T*(f1(C[i][k])-f1( C[i][k]-T*R[i][k]*h ))*C[j][k];
  }
  k=0;
  //  I-=.5*f2( C[i][k] )*R[i][k]*C[j][k];
  k=i;
  //  I+=.5*f2( C[i][k] )*R[i][k]*C[j][k];
  I*=h;
  return I;
}

double If1R(int i, int j){
  double I=0;
  int k;
  for (k=0;k<j;k++){
    I+=f1( C[i][k] )*R[j][k];
  }
  k=0;
  //   I-=.5*f1( C[i][k] )*R[j][k];
  k=j;
  //   I+=.5*f1( C[i][k] )*R[j][k];
  I*=h;
  return I;
}

void prog1(){
  int i,j,k;
  //  double B,Bp;
  //  B=1/T;
  //  Bp=1/Tp;Bp=0.; //Achtung
  for(i=0;i<n1;i++){for(j=i;j<n1;j++){
      C[i][j]=0; R[i][j]=0;
    } }
  C[0][0]=1;  R[0][0]=0;
  mu[0]=T+Bp*f1(1.);
  for (i=0;i<n1;i++){
    if(i % 100 ==0) { printf("\r%d",i); fflush(stdout); }
    R[i][i]=0; C[i][i]=1;
    for(j=0;j<i+1;j++){
      R[i+1][j]=R[i][j]+h*(-mu[i]*R[i][j]+If2RR(i,j));
      R[j][i+1]=R[i+1][j];
      C[i+1][j]=C[i][j]+h*(-mu[i]*C[i][j]+If2RC(i,j)+If1R(i,j)+Bp*f1(C[i][0])*C[j][0] );
      C[j][i+1]=C[i+1][j];
        if( j>0 && C[i][j]<q && C[i][j-1]>=q){
	/////fprintf(f13,"%d %d %f %f %f\n",i,j,/*T*R[i][j], (C[i][j]-C[i][j-1])/h,*/C[i][j],q,C[i][j-1]);
        }
     }
    //    R[i+1][i]=1+h*(-mu[i]+If2RR(i,i));
    R[i+1][i]=1;                              //THIS PASSAGE IS FUNDAMENTAL FOR THE STABILITY
    R[i][i+1]=R[i+1][i];                      //
    C[i+1][i+1]=1;
    mu[i+1]=If1R(i+1,i+1)+If2RC(i+1,i+1)+T+Bp*f1(C[i+1][0])*C[i+1][0]; // Questa e' la prescrizione migliore per mu
  

    // mu[i+1]=(If1R(i,i+1)+If2RC(i,i+1)+T+Bp*f1(C[i+1][0])*C[i+1][0]);
    // mu[i+1]+=.5*(If1R(i,i+1)+If2RC(i,i+1)+T+Bp*f1(C[i][0])*C[i][0]);

    double ene=-Bp*f(C[i][0])-If1R(i,i);
    fprintf(f11,"%f %f %f %f %f\n",i*h,C[i][0],R[i][0],ene,mu[i]);
 }
  fprintf(f11,"\n\n");
  //  fprintf(f12,"\n\n");
  //  fprintf(f13,"\n\n");
  /*for(i=n1;i>n1/5-1;i-=n1/5){
    for(j=i-1;j>=0;j--){
      if( j>0 && C[i][j]>q && C[i][j-1]<=q){
	 ////fprintf(f13,"%d %d %f %f %f %f %f\n",i,j,T*R[i][j], (C[i][j]-C[i][j-1])/h,C[i][j],q,C[i][j-1]);
      }
      double chi=0;
	for(k=j;k<i+1;k++)chi+=R[i][k];
	fprintf(f14,"%f %f %f %f %f %f\n",h*i,h*j,C[i][j],R[i][j],(C[i][j+1]-C[i][j])/h,chi*h);
    }
    fprintf(f14,"\n\n");
  }*/
}

double If2QC(int i, int j){
  double I=0;
  int k;
  for ( k=0;k<i;k++){
    I+=1./T*f2( C[i][k] )*(Q[i][k+1]-Q[i][k])*C[j][k];
  }
  return I;
}

double If2QCb(int i, int j){
  double I=0;
  int k;
  for ( k=0;k<j;k++){
    I+=1./T*f2( C[i][k] )*(Q[i][k+1]-Q[i][k])*C[j][k];
  }
  return I;
}

double If1Q(int i, int j){
    double I=0;
    int k;
  for ( k=0;k<j;k++){
      I+=1./T*f1( C[i][k] )*(Q[j][k+1]-Q[j][k]);
  }
  return I;
}

double If1Qb(int i, int j){
    double I=0;
    int k;
  for ( k=j;k<i;k++){
      I+=1./T*f1( C[i][k] )*(Q[k+1][j]-Q[k][j]);
  }
  I*=h;
  return I;
}

double If1C(int i, int j){
    double I=0;
    int k;
  for ( k=j;k<i;k++){
      I+=1./T*f1( C[i][k] )*(C[j][k+1]-C[j][k]);
  }
  return I;
}

double If2QQ(int i, int j){
  double I=0;
  int k;
  for ( k=j;k<i;k++){
      I+=1./T*f2( C[i][k] )*(Q[i][k+1]-Q[i][k])*(1.-Q[k][j]);
  }
  return I;
}

double I1C(int i, int j){
	double I=0;
    int k;
    for(k=j;k<i;k++) {
        I+=0.5*(f1(C[i][k+1])+f1(C[i][k]))*(C[j][k+1]-C[j][k]);
    }
    return I;
}

double I2C(int i, int j){
	double I=0;
    int k;
    for(k=j;k<i;k++) {
        I+=0.25*(f2(C[i][k+1])+f2(C[i][k]))*(Q[i][k+1]-Q[i][k])*(C[j][k+1]+C[j][k]);
    }
    return I;
}

double I3C(int i, int j){
	double I=0;
    int k;
    for(k=0;k<j;k++) {
        I+=0.5*(f1(C[i][k+1])+f1(C[i][k]))*(Q[j][k+1]-Q[j][k]);
    }
    return I;
}

double I4C(int i, int j){
	double I=0;
    int k;
    for(k=0;k<j;k++) {
        I+=0.25*(f2(C[i][k+1])+f2(C[i][k]))*(Q[i][k+1]-Q[i][k])*(C[j][k+1]+C[j][k]);
    }
    return I;
}

double I1Q(int i, int j){
	double I=0;
    int k;
    for(k=j;k<i;k++) {
        I+=0.5*(f1(C[i][k+1])+f1(C[i][k]))*(Q[j][k+1]-Q[j][k]);
    }
    return I;
}

double I2Q(int i, int j){
	double I=0;
    int k;
    for(k=j;k<i;k++) {
        I+=0.25*(f2(C[i][k+1])+f2(C[i][k]))*(Q[i][k+1]-Q[i][k])*(2.-Q[j][k+1]-Q[j][k]);
    }
    return I;
}

void prog2(){
    double ene;
    double i1C,i2C,i3C,i4C,i1Q,i2Q,i3Q,i4Q;
    int i,j,k;

    for(i=0;i<n1;i++){
        for(j=i;j<n1;j++){
            C[i][j]=0; Q[i][j]=0;
        } 
    }

    C[0][0]=1; Q[0][0]=0;
    mu[0]=T+(Bp-1./T)*f1(1.);

    for(i=0;i<n1;i++){
        if(i % 100 ==0) { printf("\r%d",i); fflush(stdout); }
        C[i][i]=1; Q[i][i]=0;

        for(j=0;j<i+1;j++){
            i1C = I1C(i,j);
            i2C = I2C(i,j);
            i3C = I3C(i,j);
            i4C = I4C(i,j);
            i1Q = I1Q(i,j);
            i2Q = I2Q(i,j);
            i3Q = i3C;
            i4Q = i4C;
      
            Q[i+1][j]=Q[i][j]+h*(mu[i]*(1.-Q[i][j])-T*R0+(-i1Q-i2Q-i3Q-i4Q)/T-(Bp-1./T)*f1(C[i][0])*C[j][0]);
            Q[j][i+1]=Q[i+1][j];
            C[i+1][j]=C[i][j]+h*(-mu[i]*C[i][j]+(-i1C+i2C+i3C+i4C)/T+(Bp-1./T)*f1(C[i][0])*C[j][0]);
            C[j][i+1]=C[i+1][j];
        }
        C[i+1][i+1]=1; Q[i+1][i+1]=0;
    
        mu[i+1]=(I3C(i+1,i+1)+I4C(i+1,i+1))/T+T+(Bp-1./T)*f1(C[i+1][0])*C[i+1][0]; // Questa e' la prescrizione migliore per mu
        ene = -(Bp-1./T)*f(C[i][0])-f(1.)/T-I3C(i+1,i+1)/T; //-Bp*f(C[i][0])-If1R(i,i);
        fprintf(f12,"%f %f %f %f %f\n",i*h,C[i][0],(Q[i][1]+C[i][1]-Q[i][0]-C[i][0])/h/T,ene,mu[i]);     
    }
    fprintf(f12,"\n\n");
}

double If2F1F(int i, int j){
  double I=0;
  int k;
  for ( k=j+1;k<i;k++){
    I+=f2( C[i][k] )*(F[i][k+1]-F[i][k])*F[j][k];
////      I+=1./T*(f1(C[i][k])-f1( C[i][k]-T*R[i][k]*h ))*R[j][k];
  }
  k=j;
  //  I+=.5*f2( C[i][k] )*R[i][k]*R[j][k];
  k=i;
  //  I+=.5*f2( C[i][k] )*R[i][k]*R[j][k];
  //I*=h;
  return I;
}


double If2F1C(int i, int j){
    double I=0;
    int k;
  for ( k=0;k<i;k++){
      I+=f2( C[i][k] )*(F[i][k+1]-F[i][k])*C[j][k];
////      I+=1./T*(f1(C[i][k])-f1( C[i][k]-T*R[i][k]*h ))*C[j][k];
  }
  k=0;
  //  I-=.5*f2( C[i][k] )*R[i][k]*C[j][k];
  k=i;
  //  I+=.5*f2( C[i][k] )*R[i][k]*C[j][k];
  //I*=h;
  return I;
}

double If1F1(int i, int j){
  double I=0;
  int k;
  for (k=0;k<j;k++){
    I+=f1( C[i][k] )*(F[j][k+1]-F[j][k]);
  }
  k=0;
  //   I-=.5*f1( C[i][k] )*R[j][k];
  k=j;
  //   I+=.5*f1( C[i][k] )*R[j][k];
  //I*=h;
  return I;
}

void prog3(){
  int i,j,k;
  //  double B,Bp;
  //  B=1/T;
  //  Bp=1/Tp;Bp=0.; //Achtung
  for(i=0;i<n1;i++){for(j=i;j<n1;j++){
      C[i][j]=0; F[i][j]=0;
    } }
  C[0][0]=1;  F[0][0]=0;
  mu[0]=T+Bp*f1(1.);
  for (i=0;i<n1;i++){
    if(i % 100 ==0){ printf("\r%d",i); fflush(stdout); }
    F[i][i]=0; C[i][i]=1;
    for(j=0;j<i+1;j++){
      F[i+1][j]=F[i][j]+h*(-1.-mu[i]*F[i][j]+If2F1F(i,j));
      F[j][i+1]=F[i+1][j];
      C[i+1][j]=C[i][j]+h*(-mu[i]*C[i][j]+If2F1C(i,j)+If1F1(i,j)+Bp*f1(C[i][0])*C[j][0]);
      C[j][i+1]=C[i+1][j];
        if( j>0 && C[i][j]<q && C[i][j-1]>=q){
  /////fprintf(f13,"%d %d %f %f %f\n",i,j,/*T*R[i][j], (C[i][j]-C[i][j-1])/h,*/C[i][j],q,C[i][j-1]);
        }
     }
    //    R[i+1][i]=1+h*(-mu[i]+If2RR(i,i));
    //F[i+1][i]=-h;                              //THIS PASSAGE IS FUNDAMENTAL FOR THE STABILITY
    //F[i][i+1]=F[i+1][i];                      //
    C[i+1][i+1]=1;
    mu[i+1]=If1F1(i+1,i+1)+If2F1C(i+1,i+1)+T+Bp*f1(C[i+1][0])*C[i+1][0]; // Questa e' la prescrizione migliore per mu
  

    // mu[i+1]=(If1R(i,i+1)+If2RC(i,i+1)+T+Bp*f1(C[i+1][0])*C[i+1][0]);
    // mu[i+1]+=.5*(If1R(i,i+1)+If2RC(i,i+1)+T+Bp*f1(C[i][0])*C[i][0]);

    double ene=-Bp*f(C[i][0])-If1F1(i,i);
    fprintf(f13,"%f %f %f %f %f\n",i*h,C[i][0],(F[i][1]-F[i][0])/h,ene,mu[i]);
 }
  fprintf(f13,"\n\n");
}


int main(int argc, char *argv[]){
    //Comandi cmd( argc, argv);
    printf("type(R,F or Q) eps, B, Bp, n1, dt ?\n");
    type=*argv[1];
    eps=atof(argv[2]);
    T=1./atof(argv[3]);
    Bp=atof(argv[4]);
    n1=atoi(argv[5])+1;
    h=atof(argv[6]);
    string nome=string(argv[1])+"_"+string(argv[2])+"_"+string(argv[3])+"_"+string(argv[4])+"_"+string(argv[5])+"_"+string(argv[6]);
    f11=fopen((nome+".dat").c_str(),"w");
    //f14=fopen((nome+"cr.dat").c_str(),"w");
    cout<<"#####INPUTVALUES "<<eps<<" "<<n1<<" "<<1./T<<" "<<Bp<<" "<<h<<endl;
    cout<<"eps = "<<eps<<endl;
    cout<<"dt = "<<h<<endl;

    mu.resize(n1+2);
    C.resize(n1+2,mu);

    if(type=='R') {
      R.resize(n1+2,mu);
      cout << endl << "R" << endl;
      prog1();
      h*=2.;
      cout << endl << "dt x2" << endl;
      prog1();    
    }    
    else if(type=='Q') {
      Q.resize(n1+2,mu);
      cout << endl << "Q" << endl;
      prog2();
      h*=2.;
      cout << endl << "dt x2" << endl;
      prog2(); 
    } 
    else if(type=='F') {
      F.resize(n1+2,mu);
      cout << endl << "F" << endl;
      prog3();
      h*=2.;
      cout << endl << "dt x2" << endl;
      prog3(); 
    }

    cout << endl <<"#####OUTPUT_FORMAT " << "time | correlation | response | energy | diagonal_term -depend on R,F,Q" << endl;
}
