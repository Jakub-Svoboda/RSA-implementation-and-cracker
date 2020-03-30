/*
	RSA implementation and cracker
	Author: Jakub Svoboda
	Date 27.3.2020
	Login: xsvobo0z
	Email: xsvobo0z@stud.fit.vutbr.cz
*/

/*-------------------- Imports --------------------*/
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <gmp.h>
#include <stdlib.h>
#include <time.h>
#include <random>
#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <ctime>

using namespace std;


/*-------------------- Macros --------------------*/
#define BITLENGTH 1024 				//a length in bits of generated p and q random numbers
#define SOLSTRASITERS 1000			//the number of iterations in the Solovay–Strassen primality test


/*-------------------- Function Declarations --------------------*/
void getRandPrime(mpz_t* rand);
bool solStras(mpz_t*, uint32_t, gmp_randstate_t);
unsigned long long int getSeed();


/*-------------------- Function Definitions --------------------*/
/**
 * @brief Generate a random prime number of a large size
 * 
 * @param rand Pointer to a mpt_z variable holding the random number object
 */
void getRandPrime(mpz_t* rand, gmp_randstate_t* rstate){	
	mpz_urandomb(*rand, *rstate, BITLENGTH);	//Generate random 1024 bit length number
	int i = 0;
	while(!solStras(rand, SOLSTRASITERS, *rstate)){
		i++;
		mpz_urandomb(*rand, *rstate, BITLENGTH);		 
	}	
}


/**
 * @brief Calculates the Jacobi symbol for n|k
 * 
 * @param nn N
 * @param kk K
 * @return N|K 
 */
int jacobi(mpz_t nn, mpz_t kk){
	mpz_t tmp, tmp2, r, n, k;
	mpz_inits(tmp, tmp2, r, n, k, NULL);
	mpz_set(n,nn);
	mpz_set(k,kk);

	assert(mpz_cmp_ui(k, 0) >= 1);			//assert k > 0
	mpz_mod_ui(tmp,k,2);					//assert k % 2 == 1
	assert(mpz_cmp_ui(tmp,1) == 0);			

	mpz_mod(n,n,k);					//n = n % k
	int t = 1;						//t = 1

	while(mpz_cmp_ui(n,0) != 0){			//while n != 0:	
		mpz_mod_ui(tmp,n,2);
		while(mpz_cmp_ui(tmp,0) == 0){			//while n%2 == 0
			mpz_fdiv_q_ui(n,n,2);					//n = n/2
			mpz_mod_ui(r, k, 8);						//r = k % 8
			if(mpz_cmp_ui(r,3) == 0 || mpz_cmp_ui(r,5) == 0){		//if r == 3 || r == 5
				t = -t;												//t=-t					
			}
			mpz_mod_ui(tmp,n,2);				//loop condition reevaluation
		}
		mpz_set(tmp,n);			// n = k 
		mpz_set(n,k);
		mpz_set(k, tmp);			// k = n

		mpz_mod_ui(tmp, n, 4);
		mpz_mod_ui(tmp2, k, 4);

		if(mpz_cmp_ui(tmp,3) == 0 && mpz_cmp_ui(tmp2,3) == 0){	//if n % 4 == 3 and k % 4 == 3
			t = -t;												//t=-t		
		}
		mpz_mod(n,n,k);										//n = n%k
	}

	if(mpz_cmp_ui(k,1)==0){		//if k == 0
		mpz_clears(tmp, tmp2, r, n, k, NULL);
		return t;
	}else{
		mpz_clears(tmp, tmp2, r, n, k, NULL);
		return 0;
	}
}

/**
 * @brief Solovay-Strassen primality test
 * 
 * @param n Number of which the primality is to be tested
 * @param k Number of iteration of the solovay-strassen algorithm
 * @param rstate Random state generator
 * @return true When the number is probably prime
 * @return false When the number is not prime
 */
bool solStras(mpz_t* nn, uint32_t k, gmp_randstate_t rstate){
	mpz_t a,n,tmp,left,right;
	mpz_inits(a,n,tmp,left,right, NULL);
	mpz_set(n, *nn);
	if(mpz_cmp_ui(n,0) == 0) return false;		// 0 and 1 are not prime
	if(mpz_cmp_ui(n,1) == 0) return false;
	if(mpz_cmp_ui(n,2) == 0) return true;		// 2 is prime
	mpz_mod_ui(tmp,n,2);
	if(mpz_cmp_ui(tmp,0) == 0) return false;		//even numbers are not prime 

	for (uint32_t i = 0U; i < k; i++){				//repeat K times
		mpz_set(a,n);								//init a to n, to make sure the cond below matches
		while(mpz_cmp(a,n) >= 0 or (mpz_cmp_ui(a,2) < 0)){			//choose a randomly in the range [2,n − 1]
			mpz_urandomb(a, rstate, BITLENGTH);						//generate until smaller	
		}
		int x = jacobi(a, n);									// x = a/n
		mpz_sub_ui(left,n,1U);										//left = n-1
		mpz_fdiv_q_ui(left,left,2U);								//left = (n-1)/2
		mpz_powm(left, a, left, n);									//left = a^((n-1)/2) MOD n
		mpz_set_ui(tmp, (uint32_t)x);	//tmp = x

		mpz_mod(right,tmp,n);											//right = x MOD n

		if(x == -1) mpz_sub_ui(tmp,n,1);
		if(x == 1) mpz_add_ui(tmp,n,1);
		
		mpz_mod(tmp,tmp,n);;

		if(x == 0 || (mpz_cmp(left,tmp) != 0)){
			mpz_clears(a,n,tmp,left,right,NULL);
			return false;
		}
	}
	mpz_clears(a,n,tmp,left,right,NULL);
	return true;
}

/**
 * @brief Get a random number representing seed based on system time
 * 
 * @return unsigned long long representing seed.
 */
unsigned long long getSeed(){
	std::srand(std::time(nullptr)); // use current time as seed for random generator
	return std::rand();				//return seed
}


/**
 * @brief Euclid algorithm for finding the greatest common divisor
 * 
 * @param res Pointer to an object storing result.
 * @param ee Numerator 
 * @param phii Denominator
 */
void euclidGDC(mpz_t res, mpz_t ee, mpz_t phii){
	mpz_t e,phi,phiNew;
	mpz_inits(e,phi,phiNew,NULL);	
	
	mpz_set(e,ee);
	mpz_set(phi, phii);

	if(mpz_cmp_ui(e,0) == 0){					//if e == 0
		mpz_set(res,phi);								//return e
	}else if(mpz_cmp_ui(phi,0) == 0){
		mpz_set(res,e);
	}else{					
		//cout << phi << " " << e <<" \n";
		assert(mpz_cmp_ui(phi,0) != 0);					
		mpz_mod(phiNew,e,phi);
		euclidGDC(res, phi, phiNew);					//else gdc(b, a mod b)
	}
	mpz_clears(e, phi, phiNew, NULL);	
	return;	
}


void getE(mpz_t e, mpz_t phi, gmp_randstate_t* rstate){
	mpz_t gdcRes;
	mpz_init(gdcRes);
	while(true){												//while e not between 1 and phi or gdc(e,phi) is not 1
		mpz_urandomb(e, *rstate, BITLENGTH);						// e = random
		if(mpz_cmp(e,phi) > 0 ) continue;							// e is too large, go again
		if(mpz_cmp_ui(e,1) <= 0) continue;							// e is too small, go again
		euclidGDC(gdcRes,e,phi);									//get GDC
		if(mpz_cmp_ui(gdcRes, 1) != 0){								// gdc(e,phi) != 1, go again
			continue;
		}else{
			break;	
		}																	
	}		
	mpz_clears(gdcRes,NULL);
	return;
}


int main() {
	mpz_t p,q,n,phi,e;								
	mpz_inits(p,q,n,phi,e,NULL);

	gmp_randstate_t rstate;	 				// initilze the state object for the random generator functions
	gmp_randinit_default (rstate); 			//Initialize state with a default algorithm.
	gmp_randseed_ui(rstate, getSeed()); 	// create the generator seed for the random engine to reference

	getRandPrime(&p, &rstate);				//Randomly generate 2 large prime numbers p and q
	getRandPrime(&q, &rstate);

	mpz_mul(n, p, q);						//n = p * q
	mpz_sub_ui(p,p,1);						//p = p-1
	mpz_sub_ui(q,q,1);						//q = q-1
	mpz_mul(phi,p,q);						//phi(n) = (p - 1) * (q - 1)
	getE(e, phi, &rstate);






	gmp_randclear(rstate);					// empty the memory location for the random generator state
	mpz_clears(p,q,n,NULL);					//clear variables

	return 0;
}