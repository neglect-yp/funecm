#include <stdio.h>
#include "gmp.h"
#include <time.h>
#include "point.h"

/* convert affine point to projective point */
void afftopro(PROJECTIVE_POINT R, const AFFINE_POINT P, const mpz_t N)
{
	mpz_set(R->X, P->x);
	mpz_set(R->Y, P->y);
	mpz_set_ui(R->Z, 1);
}

/* convert projective point to affine point */
void protoaff(AFFINE_POINT R, const PROJECTIVE_POINT P, const mpz_t N)
{
	mpz_t inv;
	mpz_init(inv);

	mpz_invert(inv, P->Z, N);

	mpz_mul(R->x, P->X, inv);
	mpz_mul(R->y, P->Y, inv);

	mpz_mod(R->x, R->x, N);
	mpz_mod(R->y, R->y, N);

	mpz_clear(inv);
}

/* convert projective point to extended point */
void protoext(EXTENDED_POINT R, const PROJECTIVE_POINT P, const mpz_t N)
{
	/* (X:Y:Z) -> (XZ:YZ:XY:Z^2) */
	mpz_mul(R->X, P->X, P->Z);
	mpz_mod(R->X, R->X, N);
	mpz_mul(R->Y, P->Y, P->Z);
	mpz_mod(R->Y, R->Y, N);
	mpz_mul(R->T, P->X, P->Y);
	mpz_mod(R->T, R->T, N);
	mpz_pow_ui(R->Z, P->Z, 2);
	mpz_mod(R->Z, R->Z, N);
}

/* convert extended point to projective point */
void exttopro(PROJECTIVE_POINT R, const EXTENDED_POINT P, const mpz_t N)
{
	mpz_set(R->X, P->X);
	mpz_set(R->Y, P->Y);
	mpz_set(R->Z, P->Z);
}
