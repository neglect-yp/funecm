#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "point.h"

/* for debug */
void print_bit(unsigned long int n)
{
	while (n != 0) {
		printf("%d", n & 1);
		n >>= 1;
	}
	printf("\n");
}

static long int count_bit(unsigned long int n)
{
	long int count = 0;
	while (n != 0) {
		n >>= 1;
		count++;
	}
	return count;
}

/* scalar multiplication */
void scalar(PROJECTIVE_POINT R, PROJECTIVE_POINT P, const unsigned long int k, const mpz_t D, const int window_size, const mpz_t N)
{
	EXTENDED_POINT tP;
	extended_point_init(tP);
	protoext(tP, P, N);

	EXTENDED_POINT eP;
	extended_point_init(eP);
	extended_point_set(eP, tP);

	long int m = count_bit(k);
	char *bit = (char *)malloc(m);

	long int i = 0;
	for (i = 0; i < m; i++) {
		bit[i] = (k >> i) & 1;
	}

	/* pre-calcurate for sliding-window method */
	EXTENDED_POINT *Parray = (EXTENDED_POINT *)malloc((1<<window_size)*sizeof(EXTENDED_POINT)); // allocate 2^window_size
	extended_point_init(Parray[1]);
	extended_point_set(Parray[1], eP);
	extended_point_init(Parray[2]);
	dedicated_doubling(Parray[2], eP, N);
	mpz_t inv;
	mpz_init(inv);
	mpz_invert(inv, Parray[2]->Z, N);
	mpz_mul_mod(Parray[2]->X, Parray[2]->X, inv, N);
	mpz_mul_mod(Parray[2]->Y, Parray[2]->Y, inv, N);
	mpz_mul_mod(Parray[2]->T, Parray[2]->X, Parray[2]->Y, N);
	mpz_set_ui(Parray[2]->Z, 1);
	for (i = 1; i <= (1<<(window_size-1))-1; i++) { // from 1 to 2^(window_size-1) - 1
		extended_point_init(Parray[2*i+1]);
		extended_dedicated_add(Parray[2*i+1],Parray[2*i-1],Parray[2],N);
		mpz_invert(inv, Parray[2*i+1]->Z, N);
		mpz_mul_mod(Parray[2*i+1]->X, Parray[2*i+1]->X, inv, N);
		mpz_mul_mod(Parray[2*i+1]->Y, Parray[2*i+1]->Y, inv, N);
		mpz_mul_mod(Parray[2*i+1]->T, Parray[2*i+1]->X, Parray[2*i+1]->Y, N);
		mpz_set_ui(Parray[2*i+1]->Z, 1);
	}
	mpz_clear(inv);

	/* sliding window method */
	i = m - 1;
	long int t = ((i-(window_size-1) > 0) ? i-(window_size-1) : 0);
	while (!bit[t])
		t++;
	int h = 0;
	while (i >= t) {
		h *= 2;
		if (bit[i])
			h++;
		i--;
	}
	extended_point_set(tP, Parray[h]);

	while (i >= 0) {
		if (!bit[i]) {
			dedicated_doubling(tP, tP, N);
			i--;
		} else {
			t = ((i-(window_size-1) > 0) ? i-(window_size-1) : 0);
			while (!bit[t])
				t++;
			h = 0;
			while (i >= t) {
				h *= 2;
				if (bit[i])
					h++;
				dedicated_doubling(tP, tP, N);
				i--;
			}
			extended_dedicated_add(tP, tP, Parray[h], N);
		}
	}

	exttopro(R, tP, N);

	free(bit);
	extended_point_clear(tP);
	extended_point_clear(eP);
	extended_point_clear(Parray[1]);
	extended_point_clear(Parray[2]);
	for (i = 1; i <= (1<<(window_size-1))-1; i++) {
		extended_point_clear(Parray[2*i+1]);
	}
	free(Parray);
}
