/*
 * point.c
 * Affine座標系とProjective座標系の点の変換を行う
 *
 * 更新履歴
 * 2014/10/30 新規作成
 * 2014/11/02 引数の一部をconstに変更
 * 2016/05/25 Projective座標系からExtended Twisted Edwards Coodinate座標系への変換を追加
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
	mpz_set(R->X, P->x);
	mpz_set(R->Y, P->y);
	mpz_set_ui(R->Z, 1);
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

void protoext(EXTENDED_POINT R, const PROJECTIVE_POINT P, const mpz_t N)
{
	mpz_set(R_X, P_X);
	mpz_set(R_Y, P_Y);
	mpz_set(R_Z, P_Z);

	mpz_t inv;
	mpz_init(inv);

	mpz_invert(inv, P->Z, N);
	
	mpz_mul(R->T, P->X, P->Y);
	mpz_mul(R->T, R->T, inv);
	mpz_mod(R->T, R->T, N);

	mpz_clear(inv);
}

void exttopro(PROJECTIVE_POINT R, const EXTENDED_POINT P, const mpz_t N)
{
	mpz_set(R->X, P->X);
	mpz_set(R->Y, P->Y);
	mpz_set(R->Z, P->Z);
}
