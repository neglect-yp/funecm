/*
 * point.c
 * Affine座標系とProjective座標系の点の変換を行う
 *
 * 更新履歴
 * 2014/10/30 新規作成
 * 2014/11/02 引数の一部をconstに変更
 *
 */

#include <stdio.h>
#include "gmp.h"
#include <time.h>
#include "point.h"

/*
 * Affine座標系からProjective座標系に変換する
 * PROJECTIVE_POINT R  :Projective座標系の点を格納する
 * const AFFINE_POINT P:変換元のAffine座標の点
 * const mpz_t N       :Mod N
 */
void afftopro(PROJECTIVE_POINT R, const AFFINE_POINT P, const mpz_t N)
{
	gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state, (unsigned long int)time(NULL));

	mpz_urandomm(R->Z, state, N);
	if (mpz_cmp_ui(R->Z, 0) == 0)
		mpz_add_ui(R->Z, R->Z, 1);

	mpz_mul(R->X, P->x, R->Z);
	mpz_mul(R->Y, P->y, R->Z);

	mpz_mod(R->X, R->X, N);
	mpz_mod(R->Y, R->Y, N);

	gmp_randclear(state);
}

/*
 * Projective座標系からAffine座標系に変換する
 * AFFINE_POINT R           :格納するAffine座標系の点
 * const PROJECTIVE_POINT P :変換元のProjective座標の点
 * const mpz_t N            :Mod N
 */
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

