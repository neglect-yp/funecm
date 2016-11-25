#include <stdio.h>
#include "gmp.h"
#include "point.h"

void extended_dedicated_add(EXTENDED_POINT R, EXTENDED_POINT P, EXTENDED_POINT Q, const mpz_t N)
{
	mpz_t A,B,C,D,E,F,G,H,tmp;
	
	mpz_inits(A,B,C,D,E,F,G,H,tmp,NULL);
	
	/* A<-(Y1-X1)*(Y2+X2) */
	mpz_sub(A, P->Y, P->X);
	mpz_add(tmp, Q->Y, Q->X);
	mpz_mul_mod(A, A, tmp, N);

	/* B<-(Y1+X1)*(Y2-X2) */
	mpz_add(B, P->Y, P->X);
	mpz_sub(tmp, Q->Y, Q->X);
	mpz_mul_mod(B, B, tmp, N);

	/* C<-2*Z1*T2 */
	mpz_mul_ui(C, P->Z, 2);
	mpz_mul_mod(C, C, Q->T, N);

	/* D<-2*T1*Z2 */
	mpz_mul_ui(D, P->T, 2);
	mpz_mod(D, D, N);
	
	/* E<-D+C */
	mpz_add(E, D, C);
	mpz_mod(E, E, N);
	
	/* F<-B-A */
	mpz_sub(F, B, A);
	mpz_mod(F, F, N);

	/* G<-B+A */
	mpz_add(G, B, A);
	mpz_mod(G, G, N);
	
	/* H<-D-C */
	mpz_sub(H, D, C);
	mpz_mod(H, H, N);

	/* X3<-E*F */
	mpz_mul_mod(R->X, E, F, N);
	
	/* Y3<-G*H */
	mpz_mul_mod(R->Y, G, H, N);

	/* Z3<-F*G */
	mpz_mul_mod(R->Z, F, G, N);

	mpz_clears(A,B,C,D,E,F,G,H,tmp,NULL);
}
