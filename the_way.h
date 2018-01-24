#ifndef _THE_WAY_LIBRARY_
#define _THE_WAY_LIBRARY_

#define N_PRINT
/**************************PrInT****************************/
#ifndef N_PRINT
#include "lib/share_memory.h"

void print_vector(double *sh, double **R, int l) {

	int j,k;
	double temp_sh;

	temp_sh = *(sh-2);
	*(sh-2) = 0.;
	for(k=0;k<l;k++) {
		for(j=0;j<l;j++) { *(sh+k*l+j) = R[k][j]; }
	}
	*(sh-2) = temp_sh+1.;
}
#endif
/**********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <math.h>
#include <sys/stat.h>

#define Pi  3.141592653589793
#define Ntmax 1 << 20

#define D1 3./2.
#define D2 -2.
#define D3 1./2.
#define I1 5./12.
#define I2 8./12.
#define I3 -1./12.

struct pmct{
	int Nt,Ntexp,Nt2,Nc,rpt,den,itr;	//Nt: #(grid elements), Ntexp: log2(Nt), rpt, Nt2: Nt/2, Nsh: #(short-grid elements), den: Nt/Nsh , itr: #(cycles)
	double t,t0,eps,alpha;
	double p,s_eps,s,T,beta;			//p: (p+s)-spin, s_eps: weight of s-interaction, s: (p+s)-spin, T: temperature, beta: inverse temperature of the initial equilibrium
};

struct parr{							//see appendix C
	double dt,dmu;
	double *mu,*E;
	double **C,**X,**f1,**f2;			//X can stend for R,F or Q
	double **dCh,**dXh,**df1h,**df2h;
	double **dCv,**dXv,**df1v,**df2v;
};

struct psys{
	char file[100],dir[100];			//file and directory name
	double *sh,*sh2;
};

/**************************GrIdMaNiPuLaTiOn****************************/
void contract(struct pmct z,struct parr *px,double *dt,double *dmu);
void Ih(struct parr *px, int i, int j);
void Iv(struct parr *px, int i, int j);
void Mirroring(struct parr *px, int i, int j);
void Extrapolation(struct parr *px, struct pmct z, int i);
void Banal_Extrapolation(struct parr *px, struct pmct z, int i);
void Advanced_Extrapolation(struct parr *px, struct pmct z, int i);
/**************************GeNeRaLfUnCtIoNs****************************/
double power(double x,int p);
double f(double x,struct pmct z);
double f_1(double x,struct pmct z);
double f_2(double x,struct pmct z);
/**************************ArRay****************************/
void array_initialization(struct pmct z,struct parr *px);
/**************************PaRaMeTeRs****************************/
void parameters_initialization(struct pmct *pz,struct parr *px,struct psys *pw,char *argv[]);
/**************************ScReEn****************************/
void save_config(struct pmct z,struct parr x,struct psys w);
void open_config(struct pmct *pz,struct parr *px,struct psys *pw,char *dir);
/**************************OuTpUt****************************/
void write_parameters(struct pmct z,struct parr x,struct psys w);
void write_C_0(struct parr x,struct psys w,int ini,int ifi);
void write_C(struct pmct z,struct parr x,struct psys w,int j);

/**************************GrIdMaNiPuLaTiOn****************************/

void contract(struct pmct z,struct parr *px,double *dt,double *dmu){
	int i,j,k;
	double Dl;
	i = z.Nt;
	for(k=z.Nt-z.Nc*2+1;k<=z.Nt-z.Nc;k++){
		Dl=(px->X[i][k]-px->X[i][k-1])*(I3*(px->f2[i][k+1]*px->C[i][k+1]+px->f1[i][k+1])+
						I2*(px->f2[i][k]*px->C[i][k]+px->f1[i][k])+
						I1*(px->f2[i][k-1]*px->C[i][k-1]+px->f1[i][k-1]));
	(*dmu) += Dl;
	}
	for(i=0;i<=z.Nt2;i++){
		for(j=0;j<z.Nt2;j++){
			px->dCv[i][j]= 0.5*(px->dCv[2*i][2*j+1]+px->dCv[2*i][2*j]);
			px->dXv[i][j]= 0.5*(px->dXv[2*i][2*j+1]+px->dXv[2*i][2*j]);
			px->df1v[i][j]= 0.5*(px->df1v[2*i][2*j+1]+px->df1v[2*i][2*j]);
			px->df2v[i][j]= 0.5*(px->df2v[2*i][2*j+1]+px->df2v[2*i][2*j]);
		}
	}
	for(i=1;i<=z.Nt2;i++){
		for(j=0;j<=z.Nt2;j++){
			px->dCh[i][j]= 0.5*(px->dCh[2*i][2*j]+px->dCh[2*i-1][2*j]);
			px->dXh[i][j]= 0.5*(px->dXh[2*i][2*j]+px->dXh[2*i-1][2*j]);
			px->df1h[i][j]= 0.5*(px->df1h[2*i][2*j]+px->df1h[2*i-1][2*j]);
			px->df2h[i][j]= 0.5*(px->df2h[2*i][2*j]+px->df2h[2*i-1][2*j]);
		}
	}
	for(i=0;i<=z.Nt2;i++){
		px->mu[i]= px->mu[2*i];
		for(j=0;j<=i;j++){
			px->C[i][j]= px->C[2*i][2*j];
			px->X[i][j]= px->X[2*i][2*j];
 			px->f1[i][j]= px->f1[2*i][2*j];
			px->f2[i][j]= px->f2[2*i][2*j];
		}
	}
	(*dt) *= 2.0;
}

void Ih(struct parr *px, int i, int j) {
	if(i-2>=j) {
		px->dCh[i-1][j]=I1*px->C[i][j]+I2*px->C[i-1][j]+I3*px->C[i-2][j];
		px->dXh[i-1][j]=I1*px->X[i][j]+I2*px->X[i-1][j]+I3*px->X[i-2][j];
		px->df1h[i-1][j]=I1*px->f1[i][j]+I2*px->f1[i-1][j]+I3*px->f1[i-2][j];
		px->df2h[i-1][j]=I1*px->f2[i][j]+I2*px->f2[i-1][j]+I3*px->f2[i-2][j];
	} else if(i-1==j) {
		px->dCh[i-1][j]=0.5*px->C[i][j]+0.5*px->C[i-1][j];
		px->dXh[i-1][j]=0.5*px->X[i][j]+0.5*px->X[i-1][j];
		px->df1h[i-1][j]=0.5*px->f1[i][j]+0.5*px->f1[i-1][j];
		px->df2h[i-1][j]=0.5*px->f2[i][j]+0.5*px->f2[i-1][j];
	}
}

void Iv(struct parr *px, int i, int j) {
	if(j+2<=i) {
		px->dCv[i][j]= I1*px->C[i][j]+I2*px->C[i][j+1]+I3*px->C[i][j+2];
		px->dXv[i][j]= I1*px->X[i][j]+I2*px->X[i][j+1]+I3*px->X[i][j+2];
		px->df1v[i][j]= I1*px->f1[i][j]+I2*px->f1[i][j+1]+I3*px->f1[i][j+2];
		px->df2v[i][j]= I1*px->f2[i][j]+I2*px->f2[i][j+1]+I3*px->f2[i][j+2];
	} else if (j+1==i) {
		px->dCv[i][j]= 0.5*px->C[i][j]+0.5*px->C[i][j+1];
		px->dXv[i][j]= 0.5*px->X[i][j]+0.5*px->X[i][j+1];
		px->df1v[i][j]= 0.5*px->f1[i][j]+0.5*px->f1[i][j+1];
		px->df2v[i][j]= 0.5*px->f2[i][j]+0.5*px->f2[i][j+1];
	}
}

void Mirroring(struct parr *px, int i, int j) {
	px->C[j][i] = px->C[i][j];
	px->X[j][i] = px->X[i][j];
	px->f1[j][i] = px->f1[i][j];
	px->f2[j][i] = px->f2[i][j];

	px->dCv[j][i-1] = px->dCh[i-1][j];
	px->dXv[j][i-1] = px->dXh[i-1][j];
	px->df1v[j][i-1] = px->df1h[i-1][j];
	px->df2v[j][i-1] = px->df2h[i-1][j];

	px->dCh[j][i] = px->dCv[i][j];
	px->dXh[j][i] = px->dXv[i][j];
	px->df1h[j][i] = px->df1v[i][j];
	px->df2h[j][i] = px->df2v[i][j];
}

void Banal_Extrapolation(struct parr *px, struct pmct z, int i) {

	int j;	
	for(j=i-z.Nc;j<=i;j++){
		px->C[i][j]= px->C[i-1][j-1];
		px->X[i][j]= px->X[i-1][j-1];
		px->f1[i][j]= px->f1[i-1][j-1];
		px->f2[i][j]= px->f2[i-1][j-1];
	}
	for(j=i-z.Nc;j<i;j++){
		px->dCh[i][j]= px->dCh[i-1][j-1];
		px->dXh[i][j]= px->dXh[i-1][j-1];
		px->df1h[i][j]= px->df1h[i-1][j-1];
		px->df2h[i][j]= px->df2h[i-1][j-1];
		px->dCv[i][j]= px->dCv[i-1][j-1];
		px->dXv[i][j]= px->dXv[i-1][j-1];
		px->df1v[i][j]= px->df1v[i-1][j-1];
		px->df2v[i][j]= px->df2v[i-1][j-1];
	}

  // (2) Prepare test values
	for(j=0;j<=i-z.Nc-1;j++){
		px->C[i][j]= px->C[i-1][j];
		px->X[i][j]= px->X[i-1][j];
		px->f1[i][j]= f_1(px->C[i][j],z);
		px->f2[i][j]= f_2(px->C[i][j],z);
	}
	for(j=0;j<=i-z.Nc-1;j++){
		Ih(px,i,j);
		Iv(px,i,j);
	}
}

void Advanced_Extrapolation(struct parr *px, struct pmct z, int i) {

	int j;	
	for(j=i-z.Nc;j<=i;j++){
		px->C[i][j]= px->C[i-1][j-1];
		px->X[i][j]= px->X[i-1][j-1];
		px->f1[i][j]= px->f1[i-1][j-1];
		px->f2[i][j]= px->f2[i-1][j-1];
	}
	for(j=i-z.Nc;j<i;j++){
		px->dCh[i][j]= px->dCh[i-1][j-1];
		px->dXh[i][j]= px->dXh[i-1][j-1];
		px->df1h[i][j]= px->df1h[i-1][j-1];
		px->df2h[i][j]= px->df2h[i-1][j-1];
		px->dCv[i][j]= px->dCv[i-1][j-1];
		px->dXv[i][j]= px->dXv[i-1][j-1];
		px->df1v[i][j]= px->df1v[i-1][j-1];
		px->df2v[i][j]= px->df2v[i-1][j-1];
	}

  // (2) Prepare test values
	for(j=0;j<=i-z.Nc-1;j++){
		px->C[i][j]= (D1+1.)*px->C[i-1][j-1]+D2*px->C[i-2][j-2]+D3*px->C[i-3][j-3];
		px->X[i][j]= (D1+1.)*px->X[i-1][j-1]+D2*px->X[i-2][j-2]+D3*px->X[i-3][j-3];
		px->f1[i][j]= (D1+1.)*px->f1[i-1][j-1]+D2*px->f1[i-2][j-2]+D3*px->f1[i-3][j-3];
		px->f2[i][j]= (D1+1.)*px->f2[i-1][j-1]+D2*px->f2[i-2][j-2]+D3*px->f2[i-3][j-3];
	}
	for(j=0;j<=i-z.Nc-1;j++){
		Ih(px,i,j);
		Iv(px,i,j);
	}
}

void Extrapolation(struct parr *px, struct pmct z, int i) {
	int j;
	for(j=i;j>0;j--) {
	if(j-3>=0) {
		px->C[i][j]= (D1+1.)*px->C[i-1][j-1]+D2*px->C[i-2][j-2]+D3*px->C[i-3][j-3];
		px->X[i][j]= (D1+1.)*px->X[i-1][j-1]+D2*px->X[i-2][j-2]+D3*px->X[i-3][j-3];
		px->f1[i][j]= (D1+1.)*px->f1[i-1][j-1]+D2*px->f1[i-2][j-2]+D3*px->f1[i-3][j-3];
		px->f2[i][j]= (D1+1.)*px->f2[i-1][j-1]+D2*px->f2[i-2][j-2]+D3*px->f2[i-3][j-3];
		px->dCh[i][j]= (D1+1.)*px->dCh[i-1][j-1]+D2*px->dCh[i-2][j-2]+D3*px->dCh[i-3][j-3];
		px->dXh[i][j]= (D1+1.)*px->dXh[i-1][j-1]+D2*px->dXh[i-2][j-2]+D3*px->dXh[i-3][j-3];
		px->df1h[i][j]= (D1+1.)*px->df1h[i-1][j-1]+D2*px->df1h[i-2][j-2]+D3*px->df1h[i-3][j-3];
		px->df2h[i][j]= (D1+1.)*px->df2h[i-1][j-1]+D2*px->df2h[i-2][j-2]+D3*px->df2h[i-3][j-3];
		px->dCv[i][j]= (D1+1.)*px->dCv[i-1][j-1]+D2*px->dCv[i-2][j-2]+D3*px->dCv[i-3][j-3];
		px->dXv[i][j]= (D1+1.)*px->dXv[i-1][j-1]+D2*px->dXv[i-2][j-2]+D3*px->dXv[i-3][j-3];
		px->df1v[i][j]= (D1+1.)*px->df1v[i-1][j-1]+D2*px->df1v[i-2][j-2]+D3*px->df1v[i-3][j-3];
		px->df2v[i][j]= (D1+1.)*px->df2v[i-1][j-1]+D2*px->df2v[i-2][j-2]+D3*px->df2v[i-3][j-3];
	} else if (j-2==0) {
		px->C[i][j]= (D1+1.)*px->C[i-1][2]+D2*px->C[i-2][2]+D3*px->C[i-3][2];
		px->X[i][j]= (D1+1.)*px->X[i-1][2]+D2*px->X[i-2][2]+D3*px->X[i-3][2];
		px->f1[i][j]= (D1+1.)*px->f1[i-1][2]+D2*px->f1[i-2][2]+D3*px->f1[i-3][2];
		px->f2[i][j]= (D1+1.)*px->f2[i-1][2]+D2*px->f2[i-2][2]+D3*px->f2[i-3][2];
		px->dCh[i][j]= (D1+1.)*px->dCh[i-1][2]+D2*px->dCh[i-2][2]+D3*px->dCh[i-3][2];
		px->dXh[i][j]= (D1+1.)*px->dXh[i-1][2]+D2*px->dXh[i-2][2]+D3*px->dXh[i-3][2];
		px->df1h[i][j]= (D1+1.)*px->df1h[i-1][2]+D2*px->df1h[i-2][2]+D3*px->df1h[i-3][2];
		px->df2h[i][j]= (D1+1.)*px->df2h[i-1][2]+D2*px->df2h[i-2][2]+D3*px->df2h[i-3][2];
		px->dCv[i][j]= (D1+1.)*px->dCv[i-1][2]+D2*px->dCv[i-2][2]+D3*px->dCv[i-3][2];
		px->dXv[i][j]= (D1+1.)*px->dXv[i-1][2]+D2*px->dXv[i-2][2]+D3*px->dXv[i-3][2];
		px->df1v[i][j]= (D1+1.)*px->df1v[i-1][2]+D2*px->df1v[i-2][2]+D3*px->df1v[i-3][2];
		px->df2v[i][j]= (D1+1.)*px->df2v[i-1][2]+D2*px->df2v[i-2][2]+D3*px->df2v[i-3][2];
	} else if (j-1==0) {
		px->C[i][j]= (D1+1.)*px->C[i-1][1]+D2*px->C[i-2][1]+D3*px->C[i-3][1];
		px->X[i][j]= (D1+1.)*px->X[i-1][1]+D2*px->X[i-2][1]+D3*px->X[i-3][1];
		px->f1[i][j]= (D1+1.)*px->f1[i-1][1]+D2*px->f1[i-2][1]+D3*px->f1[i-3][1];
		px->f2[i][j]= (D1+1.)*px->f2[i-1][1]+D2*px->f2[i-2][1]+D3*px->f2[i-3][1];
		px->dCh[i][j]= (D1+1.)*px->dCh[i-1][1]+D2*px->dCh[i-2][1]+D3*px->dCh[i-3][1];
		px->dXh[i][j]= (D1+1.)*px->dXh[i-1][1]+D2*px->dXh[i-2][1]+D3*px->dXh[i-3][1];
		px->df1h[i][j]= (D1+1.)*px->df1h[i-1][1]+D2*px->df1h[i-2][1]+D3*px->df1h[i-3][1];
		px->df2h[i][j]= (D1+1.)*px->df2h[i-1][1]+D2*px->df2h[i-2][1]+D3*px->df2h[i-3][1];
		px->dCv[i][j]= (D1+1.)*px->dCv[i-1][1]+D2*px->dCv[i-2][1]+D3*px->dCv[i-3][1];
		px->dXv[i][j]= (D1+1.)*px->dXv[i-1][1]+D2*px->dXv[i-2][1]+D3*px->dXv[i-3][1];
		px->df1v[i][j]= (D1+1.)*px->df1v[i-1][1]+D2*px->df1v[i-2][1]+D3*px->df1v[i-3][1];
		px->df2v[i][j]= (D1+1.)*px->df2v[i-1][1]+D2*px->df2v[i-2][1]+D3*px->df2v[i-3][1];
	}  else if (j==0) {
		px->C[i][j]= (D1+1.)*px->C[i-1][0]+D2*px->C[i-2][0]+D3*px->C[i-3][0];
		px->X[i][j]= (D1+1.)*px->X[i-1][0]+D2*px->X[i-2][0]+D3*px->X[i-3][0];
		px->f1[i][j]= (D1+1.)*px->f1[i-1][0]+D2*px->f1[i-2][0]+D3*px->f1[i-3][0];
		px->f2[i][j]= (D1+1.)*px->f2[i-1][0]+D2*px->f2[i-2][0]+D3*px->f2[i-3][0];
		px->dCh[i][j]= (D1+1.)*px->dCh[i-1][0]+D2*px->dCh[i-2][0]+D3*px->dCh[i-3][0];
		px->dXh[i][j]= (D1+1.)*px->dXh[i-1][0]+D2*px->dXh[i-2][0]+D3*px->dXh[i-3][0];
		px->df1h[i][j]= (D1+1.)*px->df1h[i-1][0]+D2*px->df1h[i-2][0]+D3*px->df1h[i-3][0];
		px->df2h[i][j]= (D1+1.)*px->df2h[i-1][0]+D2*px->df2h[i-2][0]+D3*px->df2h[i-3][0];
		px->dCv[i][j]= (D1+1.)*px->dCv[i-1][0]+D2*px->dCv[i-2][0]+D3*px->dCv[i-3][0];
		px->dXv[i][j]= (D1+1.)*px->dXv[i-1][0]+D2*px->dXv[i-2][0]+D3*px->dXv[i-3][0];
		px->df1v[i][j]= (D1+1.)*px->df1v[i-1][0]+D2*px->df1v[i-2][0]+D3*px->df1v[i-3][0];
		px->df2v[i][j]= (D1+1.)*px->df2v[i-1][0]+D2*px->df2v[i-2][0]+D3*px->df2v[i-3][0];
	}
	}
}

/**************************GeNeRaLfUnCtIoNs****************************/

double power(double x,int p){
	double ff=1.;
	int i;
	for(i=0;i<p;i++) { ff *= x; }
	return ff;
}

double f(double x, struct pmct z){
	return 0.5*power(x,z.p) + z.s_eps*0.5*power(x,z.s);
}

double f_1(double x, struct pmct z){
	return 0.5*(z.p)*power(x,z.p-1.0) + z.s_eps*0.5*(z.s)*power(x,z.s-1.0);
}

double f_2(double x, struct pmct z){
	return 0.5*(z.p)*(z.p-1.0)*power(x,z.p-2.0) + z.s_eps*0.5*(z.s)*(z.s-1.0)*power(x,z.s-2.0);
}

/**************************ArRay****************************/

void array_initialization(struct pmct z,struct parr *px) {
	int i;

	px->mu = (double *) calloc ((z.Nt+1),sizeof(double));
	px->E = (double *) calloc ((z.Nt+1),sizeof(double));
	px->C  = (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->X  = (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->f1  = (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->f2  = (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->dCh= (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->dXh= (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->df1h= (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->df2h= (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->dCv= (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->dXv= (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->df1v= (double **) malloc ((z.Nt+1) * sizeof(double*));
	px->df2v= (double **) malloc ((z.Nt+1) * sizeof(double*));


	for(i=0;i<(z.Nt+1);i++) {
		px->C[i]  = (double *) calloc ((z.Nt+1),sizeof(double));
		px->X[i]  = (double *) calloc ((z.Nt+1),sizeof(double));
		px->f1[i]  = (double *) calloc ((z.Nt+1),sizeof(double));
		px->f2[i]  = (double *) calloc ((z.Nt+1),sizeof(double));
		px->dCh[i]= (double *) calloc ((z.Nt+1),sizeof(double));
		px->dXh[i]= (double *) calloc ((z.Nt+1),sizeof(double));
		px->df1h[i]= (double *) calloc ((z.Nt+1),sizeof(double));
		px->df2h[i]= (double *) calloc ((z.Nt+1),sizeof(double));
		px->dCv[i]= (double *) calloc ((z.Nt+1),sizeof(double));
		px->dXv[i]= (double *) calloc ((z.Nt+1),sizeof(double));
		px->df1v[i]= (double *) calloc ((z.Nt+1),sizeof(double));
		px->df2v[i]= (double *) calloc ((z.Nt+1),sizeof(double));
	}
}

/**************************PaRaMeTeRs****************************/

void parameters_initialization(struct pmct *pz,struct parr *px, struct psys *pw, char *argv[]) {
	/* system parameters */
	pz->T = 1./atof(argv[2]);		//  T
	pz->beta = atof(argv[3]);		//  beta'
	pz->p     = 3.0;			//  f1(x) = x^p/2
	pz->s_eps = atof(argv[1]);		//	weight of s-spin
	pz->s     = 4.0;			//  f2(x) = x^s/2

	/* mct parameters */
	pz->Ntexp = atoi(argv[4]);         	//  Nt=2^{Ntexp}  [5..10]
	pz->Nt = 1 << pz->Ntexp;
	pz->Nt2= (int)pz->Nt/2;

	pz->den  = 1 << 5;			// Nt/Nc
	pz->Nc = 1 << 2;
	pz->itr = 34;				// number of cycle.

	pz->t0 = atof(argv[5]);			// Initial time window
	pz->eps= 1.0E-12; 			// Accepted error distance between the solution and the discrete equation ->
	pz->rpt= 1000;				// 		in the iteration procedure with z.rpt the maximum number of iterations

	pz->alpha = 1;				// Coefficient of self-consistence iteration

	/* array parameters */
	px->dt = pz->t0/pz->Nt2;		// Grid time set to Initial Grid time
	px->dmu = 0.;				//IMPORTANT: local upgrade of mu

	/* output files */
	sprintf(pw->file,"0.dat");
	sprintf(pw->dir,"eps%.2f_1suT%.5f_beta%.5fNt%dt0%.2e",pz->s_eps,1./pz->T,pz->beta,pz->Nt,pz->t0);
    mkdir(pw->dir, 0700);
}

/**************************ScReEn****************************/

void save_config(struct pmct z,struct parr x,struct psys w) {
	int i,j;
	FILE *f;
	char fn[100];
	sprintf(fn,"%s/config",w.dir);
	if((f=fopen(fn, "w"))==NULL){
		fprintf(stderr,"save_config: Cannot open a outfile\n");
		exit(1);
	}
	//pmct z
	fprintf(f,"%d %d %d %d %d %d %d\n",z.Nt,z.Ntexp,z.Nt2,z.Nc,z.rpt,z.den,z.itr);
	fprintf(f,"%.5e %.5e %.5e %.5e\n",z.t,z.t0,z.eps,z.alpha);
	fprintf(f,"%.5e %.5e %.5e %.10e %.10e\n",z.p,z.s_eps,z.s,z.T,z.beta);
	//psys w
	fprintf(f,"%s %s\n",w.file,w.dir);
	//parr x
	fprintf(f,"%.5e %.5e\n",x.dt,x.dmu);
	
	for (i=0; i<z.Nt+1; i++) {
		fprintf(f,"%.5e ",x.mu[i]);
	}
	fprintf(f,"\n");

	for (i=0; i<z.Nt+1; i++) {
		for (j=0; j<z.Nt+1; j++) {
			fprintf(f,"%.5e %.5e %.5e %.5e ",x.C[i][j],x.X[i][j],x.f1[i][j],x.f2[i][j]);
			fprintf(f,"%.5e %.5e %.5e %.5e ",x.dCh[i][j],x.dXh[i][j],x.df1h[i][j],x.df2h[i][j]);
			fprintf(f,"%.5e %.5e %.5e %.5e ",x.dCv[i][j],x.dXv[i][j],x.df1v[i][j],x.df2v[i][j]);
		}	
	}
	fprintf(f,"\n");
}

void open_config(struct pmct *pz,struct parr *px,struct psys *pw, char *dir) {
	int i,j;
	FILE *f;
	char fn[100];
	sprintf(fn,"%s/config",dir);
	if((f=fopen(fn, "r"))==NULL){
		fprintf(stderr,"open_config: Cannot open a outfile\n");
		exit(1);
	}
	//pmct z
	fscanf(f,"%d %d %d %d %d %d %d\n",&pz->Nt,&pz->Ntexp,&pz->Nt2,&pz->Nc,&pz->rpt,&pz->den,&pz->itr);
	fscanf(f,"%lf %lf %lf %lf\n",&pz->t,&pz->t0,&pz->eps,&pz->alpha);
	fscanf(f,"%lf %lf %lf %lf %lf\n",&pz->p,&pz->s_eps,&pz->s,&pz->T,&pz->beta);
	//psys w
	fscanf(f,"%s %s\n",pw->file,pw->dir);
	//parr x
	fscanf(f,"%lf %lf\n",&px->dt,&px->dmu);

	array_initialization(*pz,px);
	
	for (i=0; i<pz->Nt+1; i++) {
		fscanf(f,"%lf ",&px->mu[i]);
	}
	fscanf(f,"\n");

	for (i=0; i<pz->Nt+1; i++) {
		for (j=0; j<pz->Nt+1; j++) {
			fscanf(f,"%lf %lf %lf %lf ",&px->C[i][j],&px->X[i][j],&px->f1[i][j],&px->f2[i][j]);
			fscanf(f,"%lf %lf %lf %lf ",&px->dCh[i][j],&px->dXh[i][j],&px->df1h[i][j],&px->df2h[i][j]);
			fscanf(f,"%lf %lf %lf %lf ",&px->dCv[i][j],&px->dXv[i][j],&px->df1v[i][j],&px->df2v[i][j]);
		}	
	}
	fscanf(f,"\n");
}

/**************************OuTpUt****************************/

void write_parameters(struct pmct z,struct parr x,struct psys w){
  FILE *fout;
  char fn[100];
  sprintf(fn,"%s/par.txt",w.dir);
  if((fout=fopen(fn, "w"))==NULL){
    fprintf(stderr,"write_parameters: Cannot open a outfile\n");
    exit(1);
  }

  fprintf(fout,"####[Parameters]##########################");
  fprintf(fout,"##########################################\n");
  fprintf(fout,"# Program name:   '%s'\n",thisfile);
  fprintf(fout,"# This directory name: '%s'\n",w.dir);
  fprintf(fout,"####[System-related parameters]############");
  fprintf(fout,"##########################################\n");
  fprintf(fout,"# z.T      =%.5f\n", z.T);
  fprintf(fout,"# z.beta   =%.5f\n", z.beta);
  fprintf(fout,"# z.Nt     =2^%d=%d\n",   z.Ntexp,z.Nt);
  fprintf(fout,"# z.Nc     =%d   \n",   z.Nc);
  fprintf(fout,"# z.p      =%.0f  \n",   z.p);
  fprintf(fout,"# z.s_eps  =%.2e  \n",   z.s_eps);
  fprintf(fout,"# z.s      =%.0f  \n",   z.s);
  fprintf(fout,"# z.itr    =%d    \n",   z.itr);
  fprintf(fout,"# z.t0     =%.2e \n",   z.t0);
  fprintf(fout,"# z.rpt    =%d   \n",   z.rpt);
  fprintf(fout,"# z.eps    =%.1e \n",   z.eps);
  fprintf(fout,"####[Time stamp]##########################");
  fprintf(fout,"##########################################\n");
  fclose(fout);
}

void write_C_0(struct parr x, struct psys w, int ini, int ifi){
	FILE *fout;
	int i;
	char fn[100];
	sprintf(fn,"%s/%s",w.dir,w.file);
	if((fout=fopen(fn, "at"))==NULL){
		fprintf(stderr,"write_C_0: Cannot open a outfile\n");
		exit(1);
	}

	for(i=ini;i<ifi;i++) {
		fprintf(fout,"%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.10e\t%1.10e\n",
			(double)i*x.dt,x.C[i][0],x.dCv[i][0],x.dCh[i][0],x.X[i][0],x.dXv[i][0],x.dXh[i][0],x.mu[i],x.E[i]);
  	}
	fclose(fout);
}

void write_C(struct pmct z, struct parr x, struct psys w, int j){
	FILE *fout;
	int i;
	char fn[100];
	sprintf(fn,"%s/%d_%.2e.dat",w.dir,j,j*x.dt);
	if((fout=fopen(fn, "w"))==NULL){
		fprintf(stderr," write_C: Cannot open a outfile\n");
		exit(1);
	}

	for(i=1;i<=z.Nt;i++) {
		fprintf(fout,"%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\n",(double)(i-j)*x.dt,x.C[i][j],x.dCv[i][j],x.dCh[i][j],x.X[i][j],x.dXv[i][j],x.dXh[i][j]);
  	}
	fclose(fout);
}



#endif