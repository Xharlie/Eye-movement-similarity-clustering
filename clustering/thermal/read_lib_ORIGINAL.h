#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define MAXSTRING 100000

typedef struct structVector{
	double x;
	double y;
}vector;


int get_viewers(char *FILENAME){
	int i, n, viewer, t, T;
	double time;
	int num_viewers;
	char *fname = FILENAME;
	char line[MAXSTRING], *start;
	FILE *list;
	
	
	//count number of times
	list = fopen(fname, "r");
	T = 0;
	while( fgets(line, MAXSTRING, list) != NULL) {
		start = line;
		num_viewers = 0;
		while( sscanf(start, "%lf%n", &time, &n) == 1) {
			start += n;
			num_viewers++;
		}
		T++;
	}
	fclose(list);
	
	return num_viewers;
}
	





double **get_vx(char *FILENAME) {
	
	int i, n, viewer, t, T;
	double time;
	int num_viewers = get_viewers(FILENAME);
	char *fname = FILENAME;
	char line[MAXSTRING], *start;
	FILE *list;
	
	
	//count number of times
	list = fopen(fname, "r");
	T = 0;
	while( fgets(line, MAXSTRING, list) != NULL) {
		start = line;
		while( sscanf(start, "%lf%n", &time, &n) == 1) {
			start += n;
		}
		T++;
	}
	fclose(list);
	
	double **vx;
	
	vx = (double **)calloc(num_viewers + 1, sizeof(double *));
	for(i = 0; i <= num_viewers; i++)
		vx[i] = (double *)calloc(T + 1, sizeof(double));
	
	list = fopen(fname, "r");
	t = 1;
	while( fgets(line, MAXSTRING, list) != NULL) {
		start = line;
		viewer = 1;
		while( sscanf(start, "%lf%n", &time, &n) == 1) {
			start += n;
			vx[viewer][t] = time;
			viewer++;
		}
		t++;
	}
	fclose(list);
	
	vx[0][0] = num_viewers;
	vx[0][1] = T;
	
	return vx;
}


double **get_vy(char *FILENAME) {
	
	int i, n, viewer, t, T;
	double time;
	int num_viewers = get_viewers(FILENAME);
	char *fname = FILENAME;
	char line[MAXSTRING], *start;
	FILE *list;
	
	
	//count number of times
	list = fopen(fname, "r");
	T = 0;
	while( fgets(line, MAXSTRING, list) != NULL) {
		start = line;
		while( sscanf(start, "%lf%n", &time, &n) == 1) {
			start += n;
		}
		T++;
	}
	fclose(list);
	
	double **vy;
	
	vy = (double **)calloc(num_viewers + 1, sizeof(double *));
	for(i = 0; i <= num_viewers; i++)
		vy[i] = (double *)calloc(T + 1, sizeof(double));
	
	list = fopen(fname, "r");
	t = 1;
	while( fgets(line, MAXSTRING, list) != NULL) {
		start = line;
		viewer = 1;
		while( sscanf(start, "%lf%n", &time, &n) == 1) {
			start += n;
			vy[viewer][t] = time;
			viewer++;
		}
		t++;
	}
	fclose(list);
	
	vy[0][0] = num_viewers;
	vy[0][1] = T;
	
	
	return vy;
}

void normalize_velocity(double **vx, double **vy) {
	int i, t, num_viewers, T;
	double mod;
	num_viewers = vx[0][0];
	T = vx[0][1];
	
	for(i = 1; i <= num_viewers; i++) {
		for(t = 1; t <= T; t++) {
			
			mod = sqrt( pow(vx[i][t],2.0) + pow(vy[i][t],2.0) );
			vx[i][t] /= mod;
			vy[i][t] /= mod;
		}
	}
}

vector get_mean_velocity(int viewer, double **vx, double **vy) {
	int t, T;
	double mean_vx, mean_vy;
	vector mean_velocity;
	T = vx[0][1];
	
	mean_vx = 0.0;
	mean_vy = 0.0;
	for(t = 1; t <= T; t++) {
		mean_vx += vx[viewer][t];
		mean_vy += vy[viewer][t];
	}
	
	mean_velocity.x = (mean_vx/T);
	mean_velocity.y = (mean_vy/T);
	
	return mean_velocity;
}


double **get_Cij(double **vx, double **vy) {
	int i, j, t, T;
	int num_viewers;
	
	num_viewers = vx[0][0];
	T = vx[0][1];
	
	double **corr;
	corr = (double **)calloc(num_viewers + 1, sizeof(double *));
	for(i = 0; i <= num_viewers; i++) {
		corr[i] = (double *)calloc(num_viewers + 1, sizeof(double));
	}
	
	double Cij;
	vector mi, mj;
	
	for(i = 1; i <= num_viewers; i++) {
		for(j = 1; j <= num_viewers; j++) {
			Cij = 0.0;
			mi = get_mean_velocity(i, vx, vy);
			mj = get_mean_velocity(j, vx, vy);
			
			for(t = 1; t <= T; t++) {
				Cij += (vx[i][t] * vx[j][t]) + (vy[i][t] * vy[j][t]);
			}
			Cij /= T;
			Cij = Cij - ( (mi.x * mj.x) + (mi.y * mj.y) );
			corr[i][j] = Cij;
		}
	}
	return corr;
}



vector *get_mag(double **vx, double **vy) {
	int i, j, t, T;
	int num_viewers;

	num_viewers = vx[0][0];
	T = vx[0][1];
	
	vector *mag;
	mag = (vector *)calloc(num_viewers + 1, sizeof(vector));
	
	for(i = 1; i <= num_viewers; i++) 
		mag[i] = get_mean_velocity(i, vx, vy);
			
	return mag;
}




















