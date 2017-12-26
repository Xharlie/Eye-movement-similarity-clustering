#include<stdio.h>
#include<stdlib.h>
#include<math.h>

/*********************************************************************/
/* RANDOM NUMBER GENERATOR */
#define FNORM   (2.3283064365e-10)
#define RANDOM  ((ira[ip++] = ira[ip1++] + ira[ip2++]) ^ ira[ip3++])
#define FRANDOM (FNORM * RANDOM)
#define pm (FRANDOM > 0.5 ? 1:-1)
unsigned int myrand, ira[256];
unsigned char ip, ip1, ip2, ip3;
unsigned int randForInit(void) {
	unsigned long long int y;
	y = myrand * 16807LL;
	myrand = (y & 0x7fffffff) + (y >> 31);
	if (myrand & 0x80000000) {
		myrand = (myrand & 0x7fffffff) + 1;
	}
	return myrand;
}
void initRandom(void) {
	int i;
	ip = 128;    
	ip1 = ip - 24;    
	ip2 = ip - 55;    
	ip3 = ip - 61;
	for (i = ip3; i < ip; i++) {
		ira[i] = randForInit();
	}
}
/*********************************************************************/

/* Generator of gaussian random variables */

float gaussRan(void) {
	static int iset = 0;
	static float gset;
	float fac, rsq, v1, v2;
	
	if (iset == 0) {
		do {
			v1 = 2.0 * FRANDOM - 1.0;
			v2 = 2.0 * FRANDOM - 1.0;
			rsq = v1 * v1 + v2 * v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac = sqrt(-2.0 * log(rsq) / rsq);
		gset = v1 * fac;
		iset = 1;
		return v2 * fac;
	} else {
		iset = 0;
		return gset;
	}
}

/*********************************************************************/
