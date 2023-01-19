#include <stdio.h>
#include <math.h>
#include <float.h>

int main(void){
	double a, _e;
	double e = 1.0;
	do{
		_e = e;
		e = e/2.0;
		a = 1.0+e;
	}while(a > 1.0);
	printf("%.5e\n",_e);
}