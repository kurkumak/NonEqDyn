#define thisfile "TQW.c"
#include "../the_way.h"

//INSTRUCTIONS gcc TQW.c lib/share_memory.c -o TQW
/*****************************************************
p-spin model. Algorithm for chi_4. Similar to the aging (Kim&Latz, 2000)
written:  08/13/05
upgraded: 03/30/06
******************************************************
Reference article (appendix C): Spontaneous and induced dynamic correlations in glass formers. II. - L.Berthier, G.Biroli, G.-P.Bouchaud, W.Kob, K.Miyazaki, D.R.Reichman 
modified: 22/03/17
******************************************************/

/*void contract(struct pmct z,struct parr *x,double *dt,double *dmu){
  int i,j;
  double Dl;
  i=z.Nt;
  for(j=z.Nt-z.Nc*2+1;j<=z.Nt-z.Nc;j++){
    Dl=(x->X[i][j]-x->X[i][j-1])*(I3*(x->f1[i][j+1]+x->f2[i][j+1]*x->C[i][j+1])
				+I2*(x->f1[i][j  ]+x->f2[i][j  ]*x->C[i][j  ])
				+I1*(x->f1[i][j-1]+x->f2[i][j-1]*x->C[i][j-1]));
    (*dmu) += Dl;
  }
  for(i=1;i<=z.Nt2;i++){
    for(j=0;j<=i-1;j++){
      x->dCh[i][j]= 0.5*(x->dCh[2*i][2*j]+x->dCh[2*i-1][2*j]);
      x->dXh[i][j]= 0.5*(x->dXh[2*i][2*j]+x->dXh[2*i-1][2*j]);
      x->df1h[i][j]	= 0.5*(x->df1h[2*i][2*j]+x->df1h[2*i-1][2*j]);
      x->df2h[i][j]= 0.5*(x->df2h[2*i][2*j]+x->df2h[2*i-1][2*j]);
      x->dCv[i][j]= 0.5*(x->dCv[2*i][2*j+1]+x->dCv[2*i][2*j]);
      x->dXv[i][j]= 0.5*(x->dXv[2*i][2*j+1]+x->dXv[2*i][2*j]);
      x->df1v[i][j]= 0.5*(x->df1v[2*i][2*j+1]+x->df1v[2*i][2*j]);
      x->df2v[i][j]= 0.5*(x->df2v[2*i][2*j+1]+x->df2v[2*i][2*j]);
    }
  }
  for(i=0;i<=z.Nt2;i++){
    for(j=0;j<=i;j++){
      x->C[i][j]= x->C[2*i][2*j];
      x->X[i][j]= x->X[2*i][2*j];
      x->f1[i][j]= x->f1[2*i][2*j];
      x->f2[i][j]= x->f2[2*i][2*j];
    }
  }
  (*dt) *= 2.0;
}*/

/**************************GeNeRaL****************************/
void mct(struct pmct z,struct parr *px,struct psys w);
void initialarray(struct pmct z,struct parr *px);
int step(int i,struct pmct z,struct parr *px,struct psys w);
/**************************SeLf-CoNsIsTeNcE****************************/
double SC2(double *gC,double *gF,double D,double mu,struct pmct z,struct parr x,int i,int j);
double I1C(struct parr x,int i,int j,int m);
double I2C(struct parr x,int i,int j,int m);
double I3C(struct parr x,int i,int j);
double I4C(struct parr x,int i,int j);
double I1Q(struct parr x,int i,int j,int m);
double I2Q(struct parr x,int i,int j,int m);

double initial_dmu(struct pmct z,struct parr x,int i);
double mu_t(struct pmct z,struct parr x,int i);
double E_t(struct pmct z,struct parr x,int i);

/*******************************MaIn*********************************/

int main(int argc, char *argv[]){

	#ifdef N_PRINT
		printf("\nWithout printing!\nTo change see #define N_PRINT\n");
	#endif

	if(argc!=6) { printf("Input 5 parameters! $./FP <eps_s> <B> <Be> <log2(Nt)> <dt0>\n"); exit(0); }

	struct pmct z;
	struct parr x;
	struct psys w;

	parameters_initialization(&z,&x,&w,argv);

	//open_config(&z,&x,&w,w.dir);

	write_parameters(z,x,w);

	/* shared memory */
	#ifndef N_PRINT
		int KEY = z.Ntexp*12;							//  for the printing procedure
		w.sh = start_double_shared_memory(KEY,z.Nt*z.Nt);			//PRINT_INITIALIZE
		w.sh2 = start_double_shared_memory(KEY+1,z.Nt*z.Nt);		//PRINT_INITIALIZE
	#endif

	printf("\nDIRECTORY_OUTPUT: %s\n", w.dir);

	/*run*/
	printf("\n-----------------------------------------------------START-------------------------------------------------------\n");
	printf("------SYSTEM: (%d+eps*%d)-spin glass, eps = %2.3f------------------------------------------------------------------\n",(int)z.p,(int)z.s,z.s_eps);
	printf("------PARAMETERS: T = %1.3f, T' = %1.3f, grid dimension = %d (Nsh = %d), initial time window = %.2e-----\n",z.T,1./z.beta,z.Nt,z.Nc,z.t0);
	printf("-----------------------------------------------------------------------------------------------------------------\n\n");
	mct(z,&x,w);
	printf("\n-------------------------------------------------------END-------------------------------------------------------\n");

	return(0);
}


/**************************GeNeRaL****************************/

void mct(struct pmct z,struct parr *px,struct psys w){

/*------------------------------array initialization------------------------------*/
	int i,j,itr;
	int *scmaxx;	//variables #iteration (itr) reached in the self consitence procedure
	int scMAX,i_scMAX;
	
	array_initialization(z,px);

	scmaxx = (int *) calloc ((z.Nt+1),sizeof(int));
	//for(j=0;j<=z.Nt;j++) scmaxx[j] = 0;
/*------------------------------end array initialization------------------------------*/

	itr=0;

	initialarray(z,px); 	// prepare the array btwn 0 <= i,j <= Nt/2
	//px->dmu = initial_dmu(z,*px,i);

	write_C_0(*px,w,0,z.Nt2);
	write_C(z,*px,w,1);

	#ifndef N_PRINT
		print_vector(w.sh,px->C,z.Nt);		//PRINT
		print_vector(w.sh2,px->X,z.Nt);		//PRINT
	#endif

	while(itr < z.itr){

		//if(scMAX>10) { save_config(z,*px,w); }

		scMAX = 0;
		for(i=z.Nt2+1;i<=z.Nt;i++){
/*------------------------------------------------------------*/
			scmaxx[i]=step(i,z,px,w);	// propagate the solution the array btwn Nt/2 <= i,j <= Nt
			if(scmaxx[i]>scMAX) { scMAX = scmaxx[i]; i_scMAX = i; }
/*------------------------------------------------------------*/
		}

		write_C_0(*px,w,z.Nt2+1,z.Nt);
		write_C(z,*px,w,1);

	#ifndef N_PRINT
		print_vector(w.sh,px->C,z.Nt);		//PRINT
		print_vector(w.sh2,px->X,z.Nt);		//PRINT
	#endif

/*------------------------------------------------------------*/
		contract(z,px,&(px->dt),&(px->dmu));	// double the size of the system
/*------------------------------------------------------------*/
		printf(" %dth/%d cycle - self_iter %d/%d: %d - window_time: %.2e \n",itr,z.itr,i_scMAX,z.Nt,scMAX,px->dt*z.Nt);
		itr++;
	}

	for(j=2;j<z.Nt;j++) {
		write_C(z,*px,w,j);
	}
}

void initialarray(struct pmct z,struct parr *px){

	int i,j;

	for(i=0;i<=z.Nt2;i++){
		px->mu[i] = 0.0;
		for(j=0;j<=i;j++){
			px->C[i][j]= 1.0 - (double)(i-j)*px->dt*z.T;	// very short time expansion --> time homogeneous initial C
			px->X[i][j]= 0.0;							// very short time expansion --> FDT initially respected (Q=0)
			px->f1[i][j]= f_1(px->C[i][j],z);
			px->f2[i][j]= f_2(px->C[i][j],z);
		}
	}
	for(i=1;i<=z.Nt2;i++){
		for(j=0;j<i;j++){
			px->dCh[i][j]= 0.5*(px->C[i-1][j]+px->C[i][j]);
			px->dXh[i][j]= 0.5*(px->X[i-1][j]+px->X[i][j]);
			px->df1h[i][j]= 0.5*(px->f1[i-1][j]+px->f1[i][j]);
			px->df2h[i][j]= 0.5*(px->f2[i-1][j]+px->f2[i][j]);
			px->dCv[i][j]= 0.5*(px->C[i][j+1]+px->C[i][j]);
			px->dXv[i][j]= 0.5*(px->X[i][j+1]+px->X[i][j]);
			px->df1v[i][j]= 0.5*(px->f1[i][j+1]+px->f1[i][j]);
			px->df2v[i][j]= 0.5*(px->f2[i][j+1]+px->f2[i][j]);
		}
	}

}

int step(int i,struct pmct z,struct parr *px,struct psys w){
	int j,scmax,err2_j=99999;
	double D,err2=1.0,err_temp2;
	double *gC,*gQ;
	gC = (double *) calloc ((z.Nt),sizeof(double));
	gQ = (double *) calloc ((z.Nt),sizeof(double));

	// (1) copy the value for the top Nsh from the previous column
	px->mu[i] = (D1+1.)*px->mu[i-1]+D2*px->mu[i-2]+D3*px->mu[i-3];

	for(j=3;j<=i;j++){
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
	}

	j=2;
		px->C[i][j]= 2.*px->C[i-1][j-1]-px->C[i-2][j-2];
		px->X[i][j]= 2.*px->X[i-1][j-1]-px->X[i-2][j-2];
		px->f1[i][j]= 2.*px->f1[i-1][j-1]-px->f1[i-2][j-2];
		px->f2[i][j]= 2.*px->f2[i-1][j-1]-px->f2[i-2][j-2];
		px->dCh[i][j]= 2.*px->dCh[i-1][j-1]-px->dCh[i-2][j-2];
		px->dXh[i][j]= 2.*px->dXh[i-1][j-1]-px->dXh[i-2][j-2];
		px->df1h[i][j]= 2.*px->df1h[i-1][j-1]-px->df1h[i-2][j-2];
		px->df2h[i][j]= 2.*px->df2h[i-1][j-1]-px->df2h[i-2][j-2];
		px->dCv[i][j]= 2.*px->dCv[i-1][j-1]-px->dCv[i-2][j-2];
		px->dXv[i][j]= 2.*px->dXv[i-1][j-1]-px->dXv[i-2][j-2];
		px->df1v[i][j]= 2.*px->df1v[i-1][j-1]-px->df1v[i-2][j-2];
		px->df2v[i][j]= 2.*px->df2v[i-1][j-1]-px->df2v[i-2][j-2];
	j=1;
		px->C[i][j]= px->C[i-1][j-1];
		px->X[i][j]= px->X[i-1][j-1];
		px->f1[i][j]= px->f1[i-1][j-1];
		px->f2[i][j]= px->f2[i-1][j-1];
		px->dCh[i][j]= px->dCh[i-1][j-1];
		px->dXh[i][j]= px->dXh[i-1][j-1];
		px->df1h[i][j]= px->df1h[i-1][j-1];
		px->df2h[i][j]= px->df2h[i-1][j-1];
		px->dCv[i][j]= px->dCv[i-1][j-1];
		px->dXv[i][j]= px->dXv[i-1][j-1];
		px->df1v[i][j]= px->df1v[i-1][j-1];
		px->df2v[i][j]= px->df2v[i-1][j-1];
	j=0;
		px->C[i][j]= 2.*px->C[i-1][0]-px->C[i-2][0];
		px->X[i][j]= 2.*px->X[i-1][0]-px->X[i-2][0];
		px->f1[i][j]= 2.*px->f1[i-1][0]-px->f1[i-2][0];
		px->f2[i][j]= 2.*px->f2[i-1][0]-px->f2[i-2][0];
		px->dCh[i][j]= 2.*px->dCh[i-1][0]-px->dCh[i-2][0];
		px->dXh[i][j]= 2.*px->dXh[i-1][0]-px->dXh[i-2][0];
		px->df1h[i][j]= 2.*px->df1h[i-1][0]-px->df1h[i-2][0];
		px->df2h[i][j]= 2.*px->df2h[i-1][0]-px->df2h[i-2][0];
		px->dCv[i][j]= 2.*px->dCv[i-1][0]-px->dCv[i-2][0];
		px->dXv[i][j]= 2.*px->dXv[i-1][0]-px->dXv[i-2][0];
		px->df1v[i][j]= 2.*px->df1v[i-1][0]-px->df1v[i-2][0];
		px->df2v[i][j]= 2.*px->df2v[i-1][0]-px->df2v[i-2][0];


  // (2) Prepare test values
	/*for(j=0;j<=i-z.Nc-1;j++){
		px->C[i][j]= px->C[i-1][j];
		px->X[i][j]= px->X[i-1][j];
		px->f1[i][j]= f_1(px->C[i][j],w);
		px->f2[i][j]= f_2(px->C[i][j],w);
	}
	for(j=0;j<=i-z.Nc-1;j++){
		px->dCh[i][j] = (-1.*px->C[i-2][j]+8.*px->C[i-1][j]+5.*px->C[i][j])/12.;
		px->dXh[i][j] = (-1.*px->X[i-2][j]+8.*px->X[i-1][j]+5.*px->X[i][j])/12.;
		px->df1h[i][j]=(-1.*px->f1[i-2][j]+8.*px->f1[i-1][j]+5.*px->f1[i][j])/12.;
		px->df2h[i][j]=(-1.*px->f2[i-2][j]+8.*px->f2[i-1][j]+5.*px->f2[i][j])/12.;
		px->dCv[i][j] = (-1.*px->C[i][j+2]+8.*px->C[i][j+1]+5.*px->C[i][j])/12.;
		px->dXv[i][j] = (-1.*px->X[i][j+2]+8.*px->X[i][j+1]+5.*px->X[i][j])/12.;
		px->df1v[i][j]=(-1.*px->f1[i][j+2]+8.*px->f1[i][j+1]+5.*px->f1[i][j])/12.;
		px->df2v[i][j]=(-1.*px->f2[i][j+2]+8.*px->f2[i][j+1]+5.*px->f2[i][j])/12.;
	}*/

	// (3) Go to the SC (self-consistence) loop
	scmax = 0;
	err2 =1.0; 

	j=1;

	while( err2 >= z.eps*z.eps && scmax < z.rpt){
	err2 =0.0;

	//****** PART ----> j<i-1

	while(j>=0){

		D = 1.5/px->dt  + px->mu[i] + px->df1v[i][i-1]/z.T;

		err_temp2 = SC2(gC,gQ,D,px->mu[i],z,*px,i,j);
		if (err_temp2>err2) { err2=err_temp2; err2_j=j; }
		// renew all variable
		px->C[i][j]+= gC[j]*z.alpha;
		px->X[i][j]+= gQ[j]*z.alpha;
		//if(x->C[i][j] >= z.Cmin) x->X[i][j]+= gQ[j];
		//else { x->C[i][j] = 0.0; }
		px->f1[i][j] = f_1(px->C[i][j],z);
		px->f2[i][j] = f_2(px->C[i][j],z);

		px->dCh[i][j]=(-1.*px->C[i-2][j]+8.*px->C[i-1][j]+5.*px->C[i][j])/12.;
		px->dXh[i][j]=(-1.*px->X[i-2][j]+8.*px->X[i-1][j]+5.*px->X[i][j])/12.;
		px->df1h[i][j]=(-1.*px->f1[i-2][j]+8.*px->f1[i-1][j]+5.*px->f1[i][j])/12.;
		px->df2h[i][j]=(-1.*px->f2[i-2][j]+8.*px->f2[i-1][j]+5.*px->f2[i][j])/12.;
		px->dCv[i][j]= (-1.*px->C[i][j+2]+8.*px->C[i][j+1]+5.*px->C[i][j])/12.;
		px->dXv[i][j]= (-1.*px->X[i][j+2]+8.*px->X[i][j+1]+5.*px->X[i][j])/12.;
		px->df1v[i][j]= (-1.*px->f1[i][j+2]+8.*px->f1[i][j+1]+5.*px->f1[i][j])/12.;
		px->df2v[i][j]= (-1.*px->f2[i][j+2]+8.*px->f2[i][j+1]+5.*px->f2[i][j])/12.;

		px->mu[i] = mu_t(z,*px,i);

		j--;
	}

	scmax++;
	if(scmax%10==0) { printf("\r%d/%d %d\r",j,i,scmax); fflush(stdout); }
	j=i-z.Nc;
		
	}

	px->E[i] = E_t(z,*px,i);

	free(gC);
	free(gQ);
	return scmax;
}

/**************************SeLf-CoNsIsTeNcE****************************/

double SC2(double *gC, double *gQ, double D, double mu, struct pmct z, struct parr x, int i, int j){
  int m;
  double i1C,i2C,i3C,i4C,i1Q,i2Q,i3Q,i4Q;
  m = (int)(0.5*(i+j));
  i1C = I1C(x,i,j,m);
  i2C = I2C(x,i,j,m);
  i3C = I3C(x,i,j);
  i4C = I4C(x,i,j);
  i1Q = I1Q(x,i,j,m);
  i2Q = I2Q(x,i,j,m);
  i3Q = i3C;
  i4Q = i4C;

  gC[j] = -0.5/x.dt*x.C[i-2][j]+2.0/x.dt*x.C[i-1][j]+(-i1C+i2C+i3C+i4C)/z.T;
  gC[j]-= (1./z.T-z.beta)*f_1(x.C[i][0],z)*x.C[j][0];
  gC[j]/= D;
  gC[j]-= x.C[i][j];
  gQ[j] = -z.T+mu-0.5/x.dt*x.X[i-2][j]+2.0/x.dt*x.X[i-1][j]+(-i1Q-i2Q-i3Q-i4Q)/z.T;
  gQ[j]+= (1./z.T-z.beta)*f_1(x.C[i][0],z)*x.C[j][0];
  gQ[j]/= D;
  gQ[j]-= x.X[i][j];
	
  return gC[j]*gC[j]+gQ[j]*gQ[j];
}

double I1C(struct parr x,int i,int j,int m){
  int l;
  double sum;
  sum = 0.0;
  sum += x.f1[i][m]*x.C[m][j]-x.f1[i][j]*x.C[j][j];
  for(l=m+1;l<=i-1;l++){
    sum += x.df1v[i][l-1]*(x.C[l][j]-x.C[l-1][j]);
  }
  sum += x.df1v[i][i-1]*(-x.C[i-1][j]);
  for(l=j+1;l<=m;l++){
    sum -= (x.f1[i][l]-x.f1[i][l-1])*x.dCh[l][j];
  }
  return sum;
}

double I2C(struct parr x,int i,int j,int m){
  int l;
  double sum;
  sum = 0.0;
  for(l=m+1;l<=i;l++) sum += x.df2v[i][l-1]*(x.X[i][l]-x.X[i][l-1])*(x.C[l][j]+x.C[l-1][j]);
  for(l=j+1;l<=m;l++) sum += (x.f2[i][l]+x.f2[i][l-1])*(x.X[i][l]-x.X[i][l-1])*x.dCh[l][j];
  return 0.5*sum;
}

double I3C(struct parr x,int i,int j){
  int l;
  double sum;
  sum = 0.0;
  sum += x.f1[i][j]*x.X[j][j]-x.f1[i][0]*x.X[j][0];
  for(l=1;l<=j;l++) sum -= (x.f1[i][l]-x.f1[i][l-1])*x.dXv[j][l-1];
  return sum;
}

double I4C(struct parr x,int i,int j){
  int l;
  double sum;
  sum = 0.0;
  for(l=1;l<=j;l++) sum += (x.f2[i][l]+x.f2[i][l-1])*(x.X[i][l]-x.X[i][l-1])*x.dCv[j][l-1];
  return 0.5*sum;
}

double I1Q(struct parr x,int i,int j,int m){
  int l;
  double sum;
  sum = 0.0;
  sum += x.f1[i][m]*x.X[m][j]-x.f1[i][j]*x.X[j][j];
  for(l=m+1;l<=i-1;l++){
    sum += x.df1v[i][l-1]*(x.X[l][j]-x.X[l-1][j]);
  }
  sum += x.df1v[i][i-1]*(-x.X[i-1][j]);
  for(l=j+1;l<=m;l++){
    sum -= (x.f1[i][l]-x.f1[i][l-1])*x.dXh[l][j];
  }
  return sum;
}

double I2Q(struct parr x,int i,int j,int m){
  int l;
  double sum;
  sum = 0.0;
  for(l=m+1;l<=i;l++) sum += x.df2v[i][l-1]*(x.X[i][l]-x.X[i][l-1])*(2.0-x.X[l][j]-x.X[l-1][j]);
  for(l=j+1;l<=m;l++) sum += (x.f2[i][l]+x.f2[i][l-1])*(x.X[i][l]-x.X[i][l-1])*(1.0-x.dXh[l][j]);
  return 0.5*sum;
}

double mu_t(struct pmct z,struct parr x,int i){
	int l;
	double mu;
	mu = z.T + x.dmu;
	for(l=1;l<=i-z.Nc;l++){
		mu+=(x.X[i][l]-x.X[i][l-1])*(I3*(x.f1[i][l+1]+x.f2[i][l+1]*x.C[i][l+1])
			+I2*(x.f1[i][l  ]+x.f2[i][l  ]*x.C[i][l  ])
			+I1*(x.f1[i][l-1]+x.f2[i][l-1]*x.C[i][l-1])
		)/z.T;
	}
	mu -= (1./z.T-z.beta)*f_1(x.C[i][0],z)*x.C[i][0];
	return mu;
}

double E_t(struct pmct z,struct parr x,int i){
	int k;
	double E = 0.;
	for(k=0;k<i;k++){
		E -= x.df1v[i][k]*(x.C[i][k+1]-x.C[i][k]+x.X[i][k+1]-x.X[i][k])/z.T;
	}
	E -= z.beta*f(x.C[i][0],z);
	return E;
}
