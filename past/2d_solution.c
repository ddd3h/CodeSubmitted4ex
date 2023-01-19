#include <stdio.h>
#include <math.h>

double shiki_1(double, double);
double shiki_2(double, double);
double shiki_1a(double, double);

int main(void){
	double a = 1000.0;
	double b = 7.0;

	printf("%.6e %.6e\n",shiki_1(a, b),shiki_1(a, b)*shiki_2(a, b) - b);
	printf("%.6e %.6e\n",shiki_1a(a, b),shiki_1a(a, b)*shiki_2(a, b) - b);

	return 0;
}


double shiki_1(double a, double b){
	return ((-a+sqrt(pow(a, 2) - 4*b))/2);
}

double shiki_2(double a, double b){
	return ((-a-sqrt(pow(a, 2) - 4*b))/2);
}

double shiki_1a(double a, double b){
	return (2*b/(-a-sqrt(pow(a,2) - 4*b)));
}