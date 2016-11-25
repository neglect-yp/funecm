#include <stdio.h>
#include "gmp.h"
#include "point.h"

void double_add(PROJECTIVE_POINT R, PROJECTIVE_POINT P, const mpz_t N)
{
	mpz_t B,C,D,E,F,H,J;

	mpz_inits(B,C,D,E,F,H,J,NULL);

	mpz_add(B,P->X,P->Y); //B = X1+Y1
	mpz_pow_ui(B,B,2); //B = (X1+Y1)^2
	mpz_mod(B,B,N);
	mpz_pow_ui(C,P->X,2); //C = X1^2
	mpz_mod(C,C,N);
	mpz_pow_ui(D,P->Y,2); //D = Y1^2 
	mpz_mod(D,D,N);
	mpz_mul_si(E,C,-1); //E=-C
	mpz_add(F,E,D); //F = E+D
	mpz_pow_ui(H,P->Z,2); //H = Z1^2
	mpz_mod(H,H,N);
	mpz_sub(J,F,H); //J = F-H
	mpz_sub(J,J,H); //J = F-2H
	
	/* calcurate X_3 */
	mpz_sub(R->X,B,C); //X3=B-C
	mpz_sub(R->X,R->X,D); //X3=B-C-D
	mpz_mul_mod(R->X,R->X,J,N); //X3
	
	/* calcurate Y_3 */
	mpz_sub(R->Y,E,D); //Y3 = E-D
	mpz_mul_mod(R->Y,R->Y,F,N); //Y3

	/* calcurate Z_3 */
	mpz_mul_mod(R->Z,F,J,N); //Z3

	mpz_clears(B,C,D,E,F,H,J,NULL);
}

void dedicated_doubling(EXTENDED_POINT R, const EXTENDED_POINT P, const mpz_t N)
{
	mpz_t A,B,C,D,E,G,F,H;

	mpz_inits(A,B,C,D,E,G,F,H,NULL);

	mpz_pow_ui(A, P->X, 2);
	mpz_mod(A, A, N);

	mpz_pow_ui(B, P->Y, 2);
	mpz_mod(B, B, N);

	mpz_pow_ui(C, P->Z, 2);
	mpz_mod(C, C, N);
	mpz_mul_ui(C, C, 2);
	mpz_mod(C, C, N);
	
	mpz_mul_si(D, A, -1);

	mpz_add(E, P->X, P->Y);
	mpz_pow_ui(E, E, 2);
	mpz_mod(E, E, N);
	mpz_sub(E, E, A);
	mpz_sub(E, E, B);

	mpz_add(G, D, B);
	
	mpz_sub(F, G, C);
	
	mpz_sub(H, D, B);

	mpz_mul_mod(R->X, E, F, N);

	mpz_mul_mod(R->Y, G, H, N);
	
	mpz_mul_mod(R->T, E, H, N);

	mpz_mul_mod(R->Z, F, G, N);

	mpz_clears(A,B,C,D,E,G,F,H,NULL);
}
