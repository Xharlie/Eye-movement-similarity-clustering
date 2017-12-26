#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "random_gen.h"
#include "read_lib_ORIGINAL.h"


//DATA ACQUISITION
double **get_vx(char *FILENAME);
double **get_vy(char *FILENAME);
void normalize_velocity(double **vx, double **vy);
vector get_mean_velocity(int viewer, double **vx, double **vy);
double **get_Cij(double **vx, double **vy);



#define PI 3.14159265359
#define C 3.0						//tunable constant
#define Num_of_MC_Sweeps 200000
#define Num_Meas 200000
#define eps 0.01



typedef struct structAngle{
	double theta;
}angle;

void init_config(angle *S, int N) {
	int i;
	for(i = 1; i <= N; i++)
		S[i].theta = PI*(2.*FRANDOM - 1.);
}

void Init_Jij(double **J, int N) {
	int i, j;
	for(i = 0; i <= N; i++) {
		for(j = 0; j <= N; j++) {
			if(j != i) {
				J[i][j] = 0.1;
			}
		}
	}
}

void Init_H(vector *H, int N) {
	int i;
	for(i = 0; i <= N; i++) {
		H[i].x = 0.1;
		H[i].y = 0.1;
	}
}

double get_energy_diff(angle *S, int N, int i, double **J, vector *H, double new_angle) {
	int j;
	double de = 0.0;
	for(j = 1; j <= N; j++){
		if(j != i) {
			de += 0.5 * J[i][j] * ( cos(S[i].theta - S[j].theta) - cos(new_angle - S[j].theta) );
		}
	}
	de += H[i].x*(cos(S[i].theta) - cos(new_angle)) + H[i].y*(sin(S[i].theta) - sin(new_angle));
	
	return de;
}

int OneMCstep(angle *S, int N, double **J, vector *H, double beta) {
	int toBeUpdated = (int)(N*FRANDOM)+1;
	double de;
	double new_theta;
	int accept = 0;
	
	new_theta = fmod(S[toBeUpdated].theta + C*(FRANDOM - 0.5), PI);
	
	new_theta = S[toBeUpdated].theta + C*(FRANDOM - 0.5);
	if(new_theta > 0) {
		new_theta = (new_theta > PI) ? (new_theta - 2*PI): new_theta;
	} else {
		new_theta = (new_theta < -PI) ? (new_theta + 2*PI) : new_theta;
	}
	
	if(fabs(new_theta) > PI)
		printf("wrong update %f\n", new_theta);
	
	de = get_energy_diff(S, N, toBeUpdated, J, H, new_theta);
	
	if(de < 0.0) {
		S[toBeUpdated].theta = new_theta;
		accept = 1;
	} else {
		if(FRANDOM < exp(-beta*de)) { 
			S[toBeUpdated].theta = new_theta;
			accept = 1;
		}
	}
	return accept;
}

double OneMCsweep(angle *S, int N, double **J, vector *H, double beta) {
	int n=1;
	double accept_ratio = 0;
	while(n<=N) {
		accept_ratio += (double)OneMCstep(S, N, J, H, beta);
		n++;
	}
	return (accept_ratio)/N;
}


void Init_mag(vector *mag, int N) {
	int i;
	for(i = 0; i <= N; i++) {
		mag[i].x = 0.0;
		mag[i].y = 0.0;
	}
}

void Init_G(double **G, int N) {
	int i, j;
	for(i = 0; i <= N; i++) {
		for(j = 0; j <= N; j++) {
			G[i][j] = 0.0;
		}
	}
}

double get_magnetization(angle *S, int N) {
	double mx = 0.0;
	double my = 0.0;
	int i;
	for(i = 1; i <= N; i++) {
		mx += cos(S[i].theta);
		my += sin(S[i].theta);
	}
	mx /= N;
	my /= N;
	
	return sqrt(mx*mx +my*my);
}

double get_energy(angle *S, int N, double **J, vector *H, double beta) {
	int i, j;
	double E_link = 0.0;
	double E_site = 0.0;
	
	for(i = 1; i <= N; i++)
		E_site += (H[i].x * cos(S[i].theta)) + (H[i].y * sin(S[i].theta));

	for(i = 1; i <= N; i++) {
		for(j = i + 1; j <= N; j++) {
			E_link += J[i][j] * (cos(S[i].theta - S[j].theta));
		}
	}	
	return (-(E_link + E_site));
}






int main(int argc, char *argv[]) {
	/* SEED */
	FILE *devran = fopen("/dev/random","r");
	fread(&myrand, 4, 1, devran);
	fclose(devran);
	initRandom();
	/***********************/	
	
	int i;
	double **vx, **vy, **Cij;
	vector *mi;
	int num_viewers, N, video_num;
	char fname_vx[MAXSTRING];
	char fname_vy[MAXSTRING];
	
	video_num = atoi(argv[1]);
	
	sprintf(fname_vx, "vxo%d.txt", video_num);
	sprintf(fname_vy, "vyo%d.txt", video_num);
	
	fprintf(stdout, "Num_viewers = %d\n", get_viewers(fname_vx));
	
	vx = get_vx(fname_vx);
	vy = get_vy(fname_vy);
	normalize_velocity(vx, vy);
	
	Cij = get_Cij(vx, vy);     //Experimental Cij
	mi  = get_mag(vx, vy);     //Experimental mi
	
	num_viewers = vx[0][0];
	N = num_viewers;
	
	
	angle *S;
	double **J, **DeltaJ;
	vector *H, *DeltaH;
	
	S = (angle *)calloc(N+1, sizeof(angle));
	J = (double **)calloc(N+1,sizeof(double *));
	DeltaJ = (double **)calloc(N+1,sizeof(double *));
	for(i = 0; i<=N; i++) {
		J[i] = (double *)calloc(N+1, sizeof(double));
		DeltaJ[i] = (double *)calloc(N+1, sizeof(double));
	}
	H = (vector *)calloc(N + 1, sizeof(vector));
	DeltaH = (vector *)calloc(N + 1, sizeof(vector));
	
	
	int j, t, iter, time_meas, tau, tau_max, Samples  ;
	double beta, rate, eta, alpha, maxDiff_C;
	
	double Jtemp, mag_ave;
	vector Htemp;
	
	vector *mag, *m;
	double **G;
	
	mag = (vector *)calloc(N + 1, sizeof(vector));
	m = (vector *)calloc(N + 1, sizeof(vector));
	
	G =(double **)calloc(N + 1, sizeof(double **));
	for(i = 0; i <= N; i++) {
		G[i] = (double *)calloc(N+1, sizeof(double));
	}
	
	
	//MONTECARLO
	Init_Jij(J, N);
	Init_H(H, N);
	init_config(S, N);
	
	beta = 1.0;
	eta = 0.01;
	tau_max = 1000;
	alpha = 0.7;
	time_meas = 5;
	maxDiff_C = 0.0;
	
	
	
	for(tau = 1; tau <= tau_max ; ++tau) {
		
		//THERMALIZATION
		rate = 0.0;
		for(iter = 1; iter <= Num_of_MC_Sweeps; ++iter) 
			rate += OneMCsweep(S, N, J, H, beta);
		
		
		//MEASUREMENTS
		Init_G(G, N);
		Init_mag(mag, N);
		Samples = 0;
		
		for(iter = 1; iter <= Num_Meas; ++iter){
			OneMCsweep(S, N, J, H, beta);
			
			if( (iter % time_meas) == 0) {
				Samples++;
				for(i = 1; i <= N; i++) {
					mag[i].x += cos(S[i].theta);
					mag[i].y += sin(S[i].theta);
				}
				
				for(i = 1; i <= N; i++) {
					for(j = 1; j <= N; j++) {
						G[i][j] += cos(S[i].theta - S[j].theta);
					}
				}
			}
		}
		//END MEASUREMENTS
		
		for(i = 1; i <= N; i++) {
			mag[i].x /= Samples;
			mag[i].y /= Samples;
		}
		for(i = 1; i <= N; i++) {
			for(j = 1; j <= N; j++) {
				G[i][j] /= Samples;
				G[i][j] = G[i][j] - (mag[i].x * mag[j].x + mag[i].y*mag[j].y);
			}
		}
		
		//GAUGE Jij and Hi
		for(i = 1; i <= N; i++) {
			for(j = 1; j <= N; j++) {
				Jtemp = J[i][j];
				
				J[i][j] = J[i][j] - eta*( G[i][j] - Cij[i][j]) + alpha*DeltaJ[i][j];
				DeltaJ[i][j] = J[i][j] - Jtemp;
			}
			Htemp.x = H[i].x;
			Htemp.y = H[i].y;
			
			H[i].x = H[i].x - eta*( mag[i].x - mi[i].x) + alpha*DeltaH[i].x;
			H[i].y = H[i].y - eta*( mag[i].y - mi[i].y) + alpha*DeltaH[i].y;
			
			DeltaH[i].x = H[i].x - Htemp.x;
			DeltaH[i].y = H[i].y - Htemp.y;
			
		}
		
		eta -= 0.01/((double)(tau_max+1));
		
	    maxDiff_C = 0.0;
		
		for(i = 1; i <= N; i++) {
			for(j = i + 1; j <= N; j++) {
				if(fabs(G[i][j] - Cij[i][j]) > maxDiff_C)
					maxDiff_C = fabs(G[i][j] - Cij[i][j]);
			}
		}
		
		char fname_corr[MAXSTRING];
		char fname_mag[MAXSTRING];
		char fname_J[MAXSTRING];
		char fname_H[MAXSTRING];
		
		sprintf(fname_corr, "corr_corr_v0_video%d.txt", video_num);
		sprintf(fname_mag, "mag_mag_v0_video%d.txt", video_num);
		sprintf(fname_J, "Jinf_v0_video%d.txt", video_num);
		sprintf(fname_H, "Hinf_v0_video%d.txt", video_num);

		FILE *fcorr = fopen(fname_corr, "w");
		FILE *fmag  = fopen(fname_mag, "w");
		FILE *fJ    = fopen(fname_J, "w");
		FILE *fH    = fopen(fname_H, "w");		
		
		//Cij_MC VS Cij_REAL
		for(i = 1; i <= N; i++) {
			for(j = i + 1; j <= N; j++) {
				fprintf(fcorr, "%f %f\n", G[i][j], Cij[i][j]);
				fflush(fcorr);
			}
		}
		fclose(fcorr);
		
		//MC_mag VS REAL_mag
		for(i = 1; i <= N; i++) {
			fprintf(fmag, "%f %f %f %f\n", mag[i].x, mi[i].x, mag[i].y, mi[i].y);
			fflush(fmag);
		}
		fclose(fmag);
		
		double meanJ=0.0;
		//J inferred
		int pair = 0;
		for(i = 1; i <= N; i++) {
			for(j = i + 1; j <= N; j++) {
				pair++;
				meanJ += J[i][j];
				fprintf(fJ, "%d %f\n", pair, J[i][j]);
				fflush(fJ);
			}
		}
		fprintf(fJ, "%f\n", meanJ/pair);
		fclose(fJ);
		
		//H inferred
		for(i = 1; i <= N; i++) {
			fprintf(fH, "%d %f %f\n", i, H[i].x, H[i].y);
			fflush(fH);
		}
		fclose(fH);
		
		//Measuring CV
		if( (maxDiff_C < eps) || (tau == tau_max)) {
			
			char fname_Cv[MAXSTRING];
			sprintf(fname_Cv, "Cv_v0_video%d.txt", video_num);
			FILE *fCv    = fopen(fname_Cv, "w");
			
			double Mag, M_square, E_mean, E_square, magnetization, energy, Heat, chi;
			
			init_config(S, N);
			for(beta = 0.7; beta <= 6.0; beta += 0.05) {
				rate = 0.0;
				//THERMALIZATION
				for(iter = 1; iter <= Num_of_MC_Sweeps; ++iter) 
					rate +=  OneMCsweep(S, N, J, H, beta);
				
				// MEASUREMENTS
				Mag = 0.0;
                M_square = 0.0;
				E_mean = 0.0;
				E_square = 0.0;
				
				Samples = 0;
				for(iter = 1; iter <= Num_Meas; ++iter){
					OneMCsweep(S, N, J, H, beta);
	
					if( (iter % time_meas) == 0) { 
						Samples++;
						magnetization = get_magnetization(S, N);    //modulus of magnetization
						Mag += magnetization;
                        M_square += (magnetization * magnetization);
						energy = get_energy(S, N, J, H, beta);
						E_mean += energy;
						E_square += (energy * energy);
					}
				}
	
				M_square /= Samples;
                Mag /= Samples;
				E_mean /= Samples;
				E_square /= Samples;
				Heat = (E_square - E_mean*E_mean)*pow(beta,2.0);
                chi = (M_square - Mag*Mag)*pow(beta,2.0);
				
				fprintf(fCv, "%f %f %f %f %f\n", beta, Mag, rate/Num_of_MC_Sweeps, Heat, chi);
				fflush(fCv);
			}
			fclose(fCv);
			fprintf(stdout, "tau = %d, maxDiff_C = %f\n", tau, maxDiff_C);
			break;
		}
	}
}
	
	
	
	
	











