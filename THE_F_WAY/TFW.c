#define thisfile "TFW.c"
#include "../the_way.h"

//INSTRUCTIONS gcc TFW.c lib/share_memory.c -o TFW
/*****************************************************
p-spin model. Algorithm for chi_4. Similar to the aging (Kim&Latz, 2000)
written:  08/13/05
upgraded: 03/30/06
******************************************************
Reference article (appendix C): Spontaneous and induced dynamic correlations in glass formers. II. - L.Berthier, G.Biroli, G.-P.Bouchaud, W.Kob, K.Miyazaki, D.R.Reichman 
modified: 22/03/17
******************************************************/

/**************************GeNeRaL****************************/
void mct(struct pmct z,struct parr *px,struct psys w);
void initialarray(struct pmct z,struct parr *px);
int step(int i,struct pmct z,struct parr *px,struct psys w);
/**************************SeLf-CoNsIsTeNcE****************************/
double SC2(double *gC,double *gF,double D,double mu,struct pmct z,struct parr x,int i,int j);
double grid_If2F1F(struct parr x,int i,int j);
double grid_If2F1C(struct parr x,int i,int j);
double grid_If1F1(struct parr x,int i,int j);
double grid_reduced_If2F1F(struct parr x,int i,int j);
double grid_reduced_If2F1C(struct parr x,int i,int j);

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
	px->dmu = initial_dmu(z,*px,i);

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

		write_C_0(*px,w,z.Nt2,z.Nt);
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
	double if2F1C,if2F1F,if1F1,bf1C,bfC;

	px->C[0][0]=1;
	px->X[0][0]=0;
	px->f1[0][0] = f_1(px->C[0][0],z);
	px->f2[0][0] = f_2(px->C[0][0],z);

	px->mu[0]=z.T+z.beta*f_1(1.,z);
	px->E[0]=-z.beta*f(px->C[0][0],z);

	for(i=0;i<z.Nt2;i++){

		px->C[i+1][i+1]=1;
		px->X[i+1][i+1]=0;
		px->f1[i+1][i+1] = f_1(px->C[i+1][i+1],z);
		px->f2[i+1][i+1] = f_2(px->C[i+1][i+1],z);

		for(j=i;j>=0;j--){

			if2F1C = grid_If2F1C(*px,i,j); //If2F1C(z,*px,i,j); //
			if2F1F = grid_If2F1F(*px,i,j); //If2F1F(z,*px,i,j); //
			if1F1 = grid_If1F1(*px,i,j); //If1F1(z,*px,i,j); //
			bf1C = z.beta*px->f1[i][0]*px->C[j][0]; //z.beta*f1(px->C[i][0],z)*px->C[j][0]; 

			px->C[i+1][j] = px->C[i][j]+px->dt*(-px->mu[i]*px->C[i][j]+if2F1C+if1F1+bf1C);
			px->X[i+1][j] = px->X[i][j]+px->dt*(-1.-px->mu[i]*px->X[i][j]+if2F1F);
			px->f1[i+1][j] = f_1(px->C[i+1][j],z);
			px->f2[i+1][j] = f_2(px->C[i+1][j],z);

			Ih(px,i+1,j);
			Iv(px,i+1,j);

			Mirroring(px,i+1,j);
		}

		if2F1C = grid_If2F1C(*px,i+1,i+1); //If2C(z,*px,i+1,i+1); //
		if1F1 = grid_If1F1(*px,i+1,i+1); //If1F1(z,*px,i+1,i+1); //
		bf1C = z.beta*px->f1[i+1][0]*px->C[i+1][0]; //z.beta*f1(px->C[i+1][0],z)*px->C[i+1][0]; 
		bfC = z.beta*f(px->C[i+1][0],z);

		px->mu[i+1]=if2F1C+if1F1+z.T+bf1C; // Questa e' la prescrizione migliore per mu
		px->E[i+1]=-bfC-if1F1;

		//px->C[i+1][i+1]=1;
		//px->X[i][i]=0;
	}
}

int step(int i,struct pmct z,struct parr *px,struct psys w){
	int j,k,scmax=0,dj=-1,err2_j=99999;
	double D,err2=1.0,err_temp2;

	double *gC,*gF;
	gC = (double *) calloc ((z.Nt),sizeof(double));
	gF = (double *) calloc ((z.Nt),sizeof(double));

	// (1) First extrapolation to begin the self-consistence loop
	Banal_Extrapolation(px,z,i);
	for(j=i;j>0;j--) { Mirroring(px,i,j); }
	 //px->mu[i] = px->mu[i-1]; //(D1+1.)*px->mu[i-1]+D2*px->mu[i-2]+D3*px->mu[i-3];

	// (2) Go to the SC (self-consistence) loop
	j=i-z.Nc-1;
	//k=j-1;

	while(err2 >= z.eps*z.eps && scmax < z.rpt){

	err2=0.; 
	//****** PART ----> j<i-1
		while(j>=0 && j<i-z.Nc){
			
			px->mu[i] = mu_t(z,*px,i);
			D = D1/px->dt + px->mu[i] - 0.5*px->df2v[i][i-1]*(px->X[i][i]-px->X[i][i-1]);

			err_temp2 = SC2(gC,gF,D,px->mu[i],z,*px,i,j);
			if (err_temp2>err2) { err2=err_temp2; err2_j=j; }

			// renew all variable
			px->C[i][j]+= gC[j]*z.alpha;
			px->X[i][j]+= gF[j]*z.alpha;
			px->f1[i][j] = f_1(px->C[i][j],z);
			px->f2[i][j] = f_2(px->C[i][j],z);

			Ih(px,i,j);
			Iv(px,i,j);
			Mirroring(px,i,j);
			
			/*Iv(p-x,i,j+1);
			Mirroring(px,i,j+1);
			if(j>0) {
				Iv(px,i,j-1);
				Mirroring(px,i,j-1);
			}*/	

			j+=dj;
			//if(j==k) { j=i-3; k--; }
		} 

		//printf("%d %e\n",err2_j,err2); fflush(stdout); getchar();

		scmax++;
		if(scmax%10==0) { printf("\r%d %d\r",i,scmax); fflush(stdout); }

		j=i-z.Nc-1;

		//j-=dj; //other option loop ...
		//dj*=-1;
	}

	px->E[i] = E_t(z,*px,i);

	free(gC);
	free(gF);

	return scmax;
}


/**************************SeLf-CoNsIsTeNcE****************************/

double SC2(double *gC, double *gF, double D, double mu, struct pmct z, struct parr x, int i, int j){
	//int m = (int)(0.5*(i+j));
	double if2F1C,if2F1F,if1F1;

	if2F1C = grid_reduced_If2F1C(x,i,j);
	if2F1F = grid_reduced_If2F1F(x,i,j);
	if1F1 = grid_If1F1(x,i,j);

	gC[j] = -D3/x.dt*x.C[i-2][j]-D2/x.dt*x.C[i-1][j]+if2F1C+if1F1;
	gC[j]+= z.beta*x.f1[i][0]*x.C[j][0];
	gC[j]+= 0.5*x.df2v[i][i-1]*(x.X[i][i]-x.X[i][i-1])*x.C[j][i-1];
	gC[j]/= D;
	gC[j]-= x.C[i][j];
	gF[j] = -D3/x.dt*x.X[i-2][j]-D2/x.dt*x.X[i-1][j]-1.+if2F1F;
	gF[j]+= 0.5*x.df2v[i][i-1]*(x.X[i][i]-x.X[i][i-1])*x.X[j][i-1];
	gF[j]/= D;
	gF[j]-= x.X[i][j];

 	//printf("%d %d %f %e %e\r",j,i,x.C[i][j],gC[j],gF[j]); fflush(stdout); getchar();

	return gC[j]*gC[j]+gF[j]*gF[j];
}

double grid_If2F1F(struct parr x,int i,int j){
	double I=0;
	int m = (i+j)/2;
	int k;
	for(k=j;k<=m;k++){
		I+=0.5*(x.f2[i][k+1]+x.f2[i][k])*(x.X[i][k+1]-x.X[i][k])*x.dXv[j][k];
	}
	for(k=m+1;k<i;k++){
		I+=0.5*(x.df2v[i][k])*(x.X[i][k+1]-x.X[i][k])*(x.X[j][k+1]+x.X[j][k]);
	}
	return I;
}

double grid_reduced_If2F1F(struct parr x,int i,int j){
	double I=0;
	int m = (i+j)/2;
	int k;
	for(k=j;k<=m;k++){
		I+=0.5*(x.f2[i][k+1]+x.f2[i][k])*(x.X[i][k+1]-x.X[i][k])*x.dXv[j][k];
	}
	for(k=m+1;k<i-1;k++){
		I+=0.5*(x.df2v[i][k])*(x.X[i][k+1]-x.X[i][k])*(x.X[j][k+1]+x.X[j][k]);
	}
	return I;
}

double grid_If2F1C(struct parr x,int i,int j){
	double I=0;
	int m = (i+j)/2;
	int k;
	for(k=0;k<=m;k++){
		I+=0.5*(x.f2[i][k+1]+x.f2[i][k])*(x.X[i][k+1]-x.X[i][k])*x.dCv[j][k];
	}
	for(k=m+1;k<i;k++){
		I+=0.5*(x.df2v[i][k])*(x.X[i][k+1]-x.X[i][k])*(x.C[j][k+1]+x.C[j][k]);
	}
	return I;
}

double grid_reduced_If2F1C(struct parr x,int i,int j){
	double I=0;
	int m = (i+j)/2;
	int k;
	for(k=0;k<=m;k++){
		I+=0.5*(x.f2[i][k+1]+x.f2[i][k])*(x.X[i][k+1]-x.X[i][k])*x.dCv[j][k];
	}
	for(k=m+1;k<i-1;k++){
		I+=0.5*(x.df2v[i][k])*(x.X[i][k+1]-x.X[i][k])*(x.C[j][k+1]+x.C[j][k]);
	}
	return I;
}

double grid_If1F1(struct parr x,int i,int j){
	double I=0;
	int k;
	for(k=0;k<j;k++){
		I+=x.df1v[i][k]*(x.X[j][k+1]-x.X[j][k]);
	}
	return I;
}

double initial_dmu(struct pmct z,struct parr x,int i){
	int k;
	double i_dmu=0.;
	double Dl;
	i = z.Nt;
	for(k=z.Nt-z.Nc+1;k<=z.Nt;k++){
		Dl=(x.X[i][k]-x.X[i][k-1])*(I3*(x.f2[i][k+1]*x.C[i][k+1]+x.f1[i][k+1])+
									I2*(x.f2[i][k]*x.C[i][k]+x.f1[i][k])+
									I1*(x.f2[i][k-1]*x.C[i][k-1]+x.f1[i][k-1]));
		i_dmu += Dl;
	}
	return i_dmu;
}

double mu_t(struct pmct z,struct parr x,int i){
	int k;
	double mu;
	mu = z.T + z.beta*x.f1[i][0]*x.C[i][0]+x.dmu;
	for(k=1;k<=i-z.Nc;k++){
		mu+=(x.X[i][k]-x.X[i][k-1])*(I3*(x.f2[i][k+1]*x.C[i][k+1]+x.f1[i][k+1])+
									I2*(x.f2[i][k]*x.C[i][k]+x.f1[i][k])+
									I1*(x.f2[i][k-1]*x.C[i][k-1]+x.f1[i][k-1]));
	}
	return mu;
}

double E_t(struct pmct z,struct parr x,int i){
	int k;
	double E = 0.;
	for(k=0;k<i;k++){
		E -= x.df1v[i][k]*(x.X[i][k+1]-x.X[i][k]);
	}
	E -= z.beta*f(x.C[i][0],z);
	return E;
}