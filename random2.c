#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int irand(int num, int A, int C, int M){
	num = (num*A+C)%M;
	return (num);
}

void kadai1(void){
	int num = 10;
	for(int i=1; i<=256; i++){
		printf("%4d ",num);
		num =irand(num, 57, 1, 256);
		printf("%4d\n",num);
	}
}

void kadai2(void){
	int num = 1;
	for(int i=1; i<=256; i++){
		printf("%5d ",num);
		num =irand(num, 12869, 6925, 32768);
		printf("%5d\n",num);
	}
}

void kadai3(void){
	int seed=1533627;
	double r = drand48();
	for(int i=1;i<=256;i++){
		printf("%f ",r);
		r = drand48();
		printf("%f\n",r);
	}
}


int main(void){
	kadai1();

	return 0;
}

