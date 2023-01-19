#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double factor;

void calc_next_step(int, double, double *,double(*)(int, double, double *), double);

double func(int i, double t, double x[]){
	switch(i){
		case 0:
			return (factor*x[2] / pow(sqrt(x[2]*x[2] + x[3]*x[3]),3));
		case 1:
			return (factor * x[3] / pow(sqrt(x[2]*x[2] + x[3]*x[3]),3));
		case 2:
			return (x[0]);
		case 3:
			return (x[1]);
			break;
	}
	return 0;
}

int main(int argc, char *argv[]){
	double E_alpha = 8.954;
	double Z_alpha = 2;
	double Z_Au =79;
	double m_alpha = 3755.7;
	double c = 2.9979E+23;
	double e2hc = 1/137.036;
	double hc = 197.327;
	double dt =500.0;
	double t, x[4];

	factor = Z_alpha * Z_Au *hc *e2hc/m_alpha;
	t = 0.0;
	x[0] = sqrt(2.0*E_alpha/m_alpha);
	x[1] = 0.0;
	x[2] = -2000.0;
	x[3] = atof(argv[1]);
	while( fabs(x[2]) < 3000.0 && fabs(x[3]) < 3000.0){
		printf("%8.0f%12.4f%12.4f\n",t,x[2],x[3]);
		calc_next_step(4,t,x,func,dt);
		t += dt;
	}
}


void calc_next_step(int n, double x, double *y, double (*dydx)(int, double, double *), double h){
	double *yt, *k1, *k2, *k3, *k4;
	yt = malloc(n*sizeof(double));
	k1 = malloc(n*sizeof(double));
	k2 = malloc(n*sizeof(double));
	k3 = malloc(n*sizeof(double));
	k4 = malloc(n*sizeof(double));

	for(int i=0; i<n; i++) yt[i] = y[i];
	for(int i=0; i<n; i++) k1[i] = h * dydx(i, x, yt);
	for(int i=0; i<n; i++) yt[i] = y[i] + 0.5*k1[i];
	for(int i=0; i<n; i++) k2[i] = h * dydx(i, x+0.5*h, yt);
	for(int i=0; i<n; i++) yt[i] = y[i] + 0.5*k2[i];
	for(int i=0; i<n; i++) k3[i] = h * dydx(i, x + 0.5*h, yt);
	for(int i=0; i<n; i++) k4[i] = h * dydx(i, x, yt);
	for(int i=0; i<n; i++) y[i] += ((k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])/6.0);
	free(yt);free(k1);free(k2);free(k3);free(k4);
}