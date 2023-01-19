#include <stdio.h>
#include <math.h>

int main(void){
	float s_float = 0.0;
	double s_double = 0.0;
	double n = pow(10,10);

	for(double i =1; i<=n; i++){
		s_float += 1/(i*(i+1));
		s_double += 1/(i*(i+1));
	}

	printf("%.6e %.6e",s_double, s_float);
}