#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

double dimen(int);

int main(void){
	printf("%.6e\n",dimen(3));
	return 0;
}


double dimen(int dim){
	double x;
	int N = 1000000;
	int m = 0;
	double volume;
	double square_sum;
	double average = 0.0;

	int seed = time(NULL);
	srand(seed);

	for(int i = 1; i <= N; i ++){
		square_sum = 0;
		for(int j = 0; j < dim; j++){
			x = drand48();
			square_sum += x*x;
		}

		if(square_sum < 1){
			m++;
		}

		if(i%100==0){
			volume = pow(2,dim)*(double)m/i;
			if(N*0.9<i){
				average += volume;
			}
		}
	}

	double ans = average/(N*0.1/100);
	return ans;
}


