#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "gmp.h"
#include "point.h"

void ecm(mpz_t f, const mpz_t N, const mpz_t X, const mpz_t Y, mpz_t d, const unsigned long int B1, const unsigned long int B2, FILE *fp, const int window_size)
{
	PROJECTIVE_POINT P;
	int e;
	int i;

	mpz_t tmp;
	mpz_t tmp2;
	mpz_t inv;
	mpz_init(tmp);
	mpz_init(tmp2);
	mpz_init(inv);

	projective_point_init(P);

	/* set P */
	if (X == NULL)
		mpz_set_ui(P->X, 2);
	else
		mpz_set(P->X, X);
	mpz_set(P->Y, Y);
	mpz_set_ui(P->Z, 1);

	/* calcurate d if atkin_flag isn't set*/
	if (X == NULL) {
		mpz_init(d);
		mpz_pow_ui(tmp,P->X,2); //tmp = x^2
		mpz_mod(tmp,tmp,N);
		mpz_pow_ui(tmp2,P->Y,2); //tmp2 = y^2
		mpz_mod(tmp2,tmp2,N);
		mpz_sub(d,tmp2,tmp); //d = y^2-x^2
		mpz_sub_ui(d,d,1); //d = y^2-x^2-1
		mpz_mul_mod(tmp,tmp,tmp2,N); //tmp = x^2y^2
		mpz_invert(tmp,tmp,N);
		mpz_mul_mod(d,d,tmp,N); //dの値
	}

	/* set prime number */
	unsigned long int p = 2;
	mpz_t prime;
	mpz_init(prime);
	mpz_set_ui(prime, p);

	double stage1_time, stage2_time = -1;
	double start, end;
	start = omp_get_wtime();
	/* stage1 */
	while (p <= B1) {
		/* e = log p k */
		e = (int)(log(B1) / log(p));
		for (i = 1; i <= e; i++) {
			/* P_Z <- 1 */
			mpz_invert(inv, P->Z, N);
			mpz_mul_mod(P->X, P->X, inv, N);
			mpz_mul_mod(P->Y, P->Y, inv, N);
			mpz_set_ui(P->Z, 1);

			scalar(P, P, p, d, window_size, N);
			mpz_gcd(f, P->X, N);
			if (mpz_cmp_ui(f,1) != 0) {
				end = omp_get_wtime();
				stage1_time = end - start;
				goto FACTOR_FOUND;
			}
		}

		mpz_nextprime(prime, prime);
		p = mpz_get_ui(prime);
	}
	end = omp_get_wtime();
	stage1_time = end - start;

	start = omp_get_wtime();
	/* stage2 */
	mpz_t product;
	mpz_init(product);
	mpz_set_ui(product, 1);
	while (p <= B2) {
		/* P_Z <- 1 */
		mpz_invert(inv, P->Z, N);
		mpz_mul_mod(P->X, P->X, inv, N);
		mpz_mul_mod(P->Y, P->Y, inv, N);
		mpz_set_ui(P->Z, 1);

		scalar(P, P, p, d, window_size, N);
		mpz_mul_mod(product, product, P->Z, N);

		mpz_nextprime(prime, prime);
		p = mpz_get_ui(prime);
	}
	mpz_gcd(f, product, N);
	mpz_clear(product);
	end = omp_get_wtime();
	stage2_time = end - start;

FACTOR_FOUND:
	gmp_fprintf(fp,"Stage1: d = %Zd\n", d);
	fprintf(fp, "Stage1 time: %f seconds\n", stage1_time);
	if (stage2_time != -1)
		fprintf(fp, "Stage2 time: %f seconds\n", stage2_time);
	else
		fprintf(fp, "Stage2 time: ----\n");

	projective_point_clear(P);
	mpz_clear(tmp);
	mpz_clear(tmp2);
	mpz_clear(inv);
	mpz_clear(prime);
}
