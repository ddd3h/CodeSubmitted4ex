#include <stdio.h>
#include <math.h>

int irand(int,int,int,int);

int main(void){
	int num = 0;
	int M = 16;
	int Maxstep = 100000;
	printf("[No1] ");
	for(int i = 0; i < Maxstep; i++){
		num = irand(num, 5, 1, M);
		if(i%M == 0){
			printf("%d ",num);
		}	
	}


	num = 1;
	M = 32768;
	printf("\n\n[No2] ");
	for(int i = 0; i < Maxstep; i++){
		num = irand(num, 12869, 6925, M);
		if(i%M == 0){
			printf("%d ",num);
		}
	}


	num = 10;
	M = 256;
	printf("\n\n[No3] ");
	for(int i = 0; i < Maxstep; i++){
		num = irand(num, 57, 1, M);
		if(i%M == 0){
			printf("%d ",num);
		}
	}

	return 0;
}



int irand(int num, int A, int C, int M){
	num = (num*A+C)%M;
	return (num);
}