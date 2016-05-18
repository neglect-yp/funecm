/*
 * point.h
 *
 * 更新履歴
 * 2014/10/21 新規作成
 * 2014/10/30 マクロ追加
 *            projective_point_cmpのバグを修正 //そもそも役に立たないかもしれない
 * 2014/11/02 関数の宣言の一部の引数をconstに変更
 */

#ifndef GMP_H
#define GMP_H
#include "gmp.h"
#endif

#ifndef POINT_H
#define POINT_H

typedef struct
{
	mpz_t x;
	mpz_t y;
} AFFINE_POINT[1];

typedef struct
{
	mpz_t X;
	mpz_t Y;
	mpz_t Z;
} PROJECTIVE_POINT[1];

void afftopro(PROJECTIVE_POINT R, const AFFINE_POINT P, const mpz_t N);
void protoaff(AFFINE_POINT R, const PROJECTIVE_POINT P, const mpz_t N);

#define affine_point_init(P) do { mpz_init(P->x); mpz_init(P->y); } while (0)
#define affine_point_clear(P) do { mpz_clear(P->x); mpz_clear(P->y); } while (0)
#define affine_point_set(R, P) do { mpz_set(R->x, P->x); mpz_set(R->y, P->y); } while (0)
#define affine_point_cmp(P, Q) (((mpz_cmp(P->x, Q->x) == 0) && (mpz_cmp(P->y, Q->y) == 0) ) ? 0 : 1)

#define projective_point_init(P) do{mpz_init(P->X); mpz_init(P->Y); mpz_init(P->Z);}while(0)
#define projective_point_clear(P) do { mpz_clear(P->X); mpz_clear(P->Y); mpz_clear(P->Z); } while (0)
#define projective_point_set(R, P) do { mpz_set(R->X, P->X); mpz_set(R->Y, P->Y); mpz_set(R->Z, P->Z); } while (0)
#define projective_point_cmp(P, Q) (((mpz_cmp(P->X, Q->X) == 0) && (mpz_cmp(P->Y, Q->Y) == 0) && (mpz_cmp(P->Z, Q->Z) == 0)) ? 0 : 1)

#endif
