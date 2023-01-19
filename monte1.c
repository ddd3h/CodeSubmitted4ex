#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

int main(void){
	double x, y, z;
	int N = 100000;
	int m = 0;
	double volume;

	int seed = time(NULL);
	srand(seed);

	for(int i = 1; i <= N; i ++){
		x = drand48();
		y = drand48();
		z = drand48();

		if (x*x + y*y + z*z < 1){
			m++;
		}

		if(i%100==0){
			volume = pow(2,3)*(double)m/(double)i;
			printf("%6d %.5e\n",i,volume);
		}
	}


	return 0;
}