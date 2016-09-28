#include <stdio.h>
#include "gmp.h"
#include "point.h"

void atkin_double_add(PROJECTIVE_POINT R, PROJECTIVE_POINT P, const mpz_t N);

void atkin_moraine (mpz_t X, mpz_t Y, mpz_t d, mpz_t s, mpz_t t, const mpz_t N) {
	mpz_t alpha, beta, tmp, tmp2;
	mpz_inits(alpha, beta, tmp, tmp2, NULL);

	/* calculate alpha  */
	mpz_sub_ui(alpha, s, 8);
	mpz_add_ui(tmp, t, 25);
	mpz_invert(tmp, tmp, N);
	mpz_mul(alpha, alpha, tmp);
	mpz_mod(alpha, alpha, N);

	/* calculate beta */
	/* beta <- 2*alpha * (4*alpha + 1) */
	mpz_mul_ui(beta, alpha, 4);
	mpz_add_ui(beta, beta, 1);
	mpz_mul(beta, beta, alpha);
	mpz_mod(beta, beta, N);
	mpz_mul_ui(beta, beta, 2);
	/* tmp <- 8 * alpha^2 -1 */
	mpz_pow_ui(tmp, alpha, 2);
	mpz_mod(tmp, tmp, N);
	mpz_mul_ui(tmp, tmp, 8);
	mpz_mod(tmp, tmp, N);
	mpz_sub_ui(tmp, tmp, 1);
	/* beta <- beta / tmp */
	mpz_invert(tmp, tmp, N);
	mpz_mul(beta, beta, tmp);
	mpz_mod(beta, beta, N);

	/* calculate d */
	/* tmp2 <- 2 * beta - 1 */
	mpz_mul_ui(tmp2, beta, 2);
	mpz_sub_ui(tmp2, tmp2, 1);
	/* d <- 2 * tmp2^2 - 1 */
	mpz_pow_ui(d, tmp2, 2);
	mpz_mod(d, d, N);
	mpz_mul_ui(d, d, 2);
	mpz_sub_ui(d, d, 1);
	/* tmp <- (2*beta - 1)^4 */
	mpz_pow_ui(tmp, tmp2, 4);
	mpz_mod(tmp, tmp, N);
	/* d <- d / tmp */
	mpz_invert(tmp, tmp, N);
	mpz_mul(d, d, tmp);
	mpz_mod(d, d, N);

	/* calculate X */
	/* X <- tmp2 * (4*beta - 3) */
	mpz_mul_ui(X, beta, 4);
	mpz_sub_ui(X, X, 3);
	mpz_mul(X, X, tmp2);
	mpz_mod(X, X, N);
	/* tmp <- 6*beta - 4 */
	mpz_mul_ui(tmp, beta, 6);
	mpz_sub_ui(tmp, tmp, 4);
	/* X <- X / tmp */
	mpz_invert(tmp, tmp, N);
	mpz_mul(X, X, tmp);
	mpz_mod(X, X, N);

	/* calculate Y */
	/* Y <- tmp2 * (t^2 + 50*t - 2*s^3 + 27*s^2 - 104) */
	mpz_pow_ui(Y, t, 2);
	mpz_mod(Y, Y, N);
	mpz_mul_ui(tmp, t, 50);
	mpz_add(Y, Y, tmp);
	mpz_pow_ui(tmp, s, 3);
	mpz_mod(tmp, tmp, N);
	mpz_mul_ui(tmp, tmp, 2);
	mpz_sub(Y, Y, tmp);
	mpz_pow_ui(tmp, s, 2);
	mpz_mod(tmp, tmp, N);
	mpz_mul_ui(tmp, tmp, 27);
	mpz_add(Y, Y, tmp);
	mpz_sub_ui(Y, Y, 104);
	mpz_mul(Y, Y, tmp2);
	mpz_mod(Y, Y, N);
	/* tmp <- t + 3*s - 2 */
	mpz_mul_ui(tmp, s, 3);
	mpz_add(tmp, tmp, t);
	mpz_sub_ui(tmp, tmp, 2);
	/* tmp2 <- t + s + 16 */
	mpz_add(tmp2, t, s);
	mpz_add_ui(tmp2, tmp2, 16);
	/* tmp <- tmp * tmp2 */
	mpz_mul(tmp, tmp, tmp2);
	mpz_mod(tmp, tmp, N);
	/* Y <- Y / tmp */
	mpz_invert(tmp, tmp, N);
	mpz_mul(Y, Y, tmp);
	mpz_mod(Y, Y, N);

	/* (s, t) <- 2(s, t) */
	AFFINE_POINT aP;
	affine_point_init(aP);
	PROJECTIVE_POINT P, R;
	projective_point_init(P);
	projective_point_init(R);
	mpz_set(aP->x, s);
	mpz_set(aP->y, t);
	afftopro(P, aP, N);
	atkin_double_add(R, P, N);
	mpz_set(s, R->X);
	mpz_set(t, R->Y);
	affine_point_clear(aP);
	projective_point_clear(P);
	projective_point_clear(R);

	mpz_clears(alpha, beta, tmp, tmp2, NULL);
}

/*
 * Projective座標系の2倍算を行う
 * PROJECTIVE_POINT R :Pの2倍算の結果を格納する点
 * PROJECTIVE_POINT P :2倍算を行う点
 * const mpz_t N :mod N
 */
void atkin_double_add(PROJECTIVE_POINT R, PROJECTIVE_POINT P, const mpz_t N)
{
	mpz_t w, s, B, h;

	PROJECTIVE_POINT tP;

	/* 一時変数 */
	mpz_t tmp;
	mpz_t tmp2;

	/* 初期化 */
	mpz_inits(w, s, B, h, tmp, tmp2, NULL);
	projective_point_init(tP);

	projective_point_set(tP, P);

	/* wの計算 */
	mpz_pow_ui(w, tP->Z, 2);
	mpz_mod(w, w, N);
	mpz_mul_si(w, w, -8);
	mpz_pow_ui(tmp, tP->X, 2);
	mpz_mod(tmp, tmp, N);
	mpz_mul_ui(tmp, tmp, 3);
	mpz_add(w, w, tmp);
	mpz_mod(w, w, N);

	/* sの計算 */
	mpz_mul(s, tP->Y, tP->Z);
	mpz_mod(s, s, N);

	/* Bの計算 */
	mpz_mul(B, tP->X, tP->Y);
	mpz_mod(B, B, N);
	mpz_mul(B, B, s);
	mpz_mod(B, B, N);

	/* hの計算 */
	mpz_pow_ui(h, w, 2);
	mpz_mod(h, h, N);
	mpz_mul_ui(tmp, B, 8);
	mpz_sub(h, h, tmp);
	mpz_mod(h, h, N);

	/* Xrの計算 */
	mpz_mul_ui(R->X, h, 2);
	mpz_mul(R->X, R->X, s);
	mpz_mod(R->X, R->X, N);

	/* Yrの計算 */
	mpz_mul_ui(R->Y, B, 4);
	mpz_sub(R->Y, R->Y, h);
	mpz_mul(R->Y, R->Y, w);
	mpz_mod(R->Y, R->Y, N);
	mpz_pow_ui(tmp, tP->Y, 2);
	mpz_mod(tmp, tmp, N);
	mpz_pow_ui(tmp2, s, 2);
	mpz_mod(tmp2, tmp2, N);
	mpz_mul(tmp, tmp, tmp2);
	mpz_mul_ui(tmp, tmp, 8);
	mpz_sub(R->Y, R->Y, tmp);
	mpz_mod(R->Y, R->Y, N);

	/* Zrの計算 */
	mpz_pow_ui(R->Z, s, 3);
	mpz_mod(R->Z, R->Z, N);
	mpz_mul_ui(R->Z, R->Z, 8);
	mpz_mod(R->Z, R->Z, N);

	/* 使用した変数の開放 */
	mpz_clears(w, s, B, h, tmp, tmp2, NULL);

	projective_point_clear(tP);
}
