/*
 * ecm.c
 * 楕円曲線法によって因数分解を行う
 *
 * 更新履歴
 * 2014/10/29 新規作成
 * 2014/10/30 mpz_t f追加
 *            gcd処理追加
 * 2014/11/02 引数の一部をconstに変更
 * 2015/10/** a値を廃止しd値を追加
 * 2015/11/17 ファイル処理
 */

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "gmp.h"
#include "point.h"

/* logメモ */
/* logab= logb/loga */

/*
 * 楕円曲線法にて因数分解を行う関数
 * mpz_t f              :発見した素因数
 * const mpz_t N              :素因数分解したい合成数
 * const unsigned long int Y  :ベースポイントPのY座標
 * const unsigned long int B1  :スカラー倍
 */
void ecm(mpz_t f, const mpz_t N, const mpz_t X, const mpz_t Y, mpz_t d, const unsigned long int B1, const unsigned long int B2, FILE *fp, const int window_size)
{
	/* 使用変数・構造体の宣言 */
	PROJECTIVE_POINT P;
	int e;
	int i;

	/* 一時変数 */
	mpz_t tmp;
	mpz_t tmp2;
	mpz_t inv;
	mpz_init(tmp);
	mpz_init(tmp2);
	mpz_init(inv);

	/* Pの初期化 */
	projective_point_init(P);

	/* Pの点の座標を指定 */
	if (X == NULL)
		mpz_set_ui(P->X, 2);
	else
		mpz_set(P->X, X);
	mpz_set(P->Y, Y);
	mpz_set_ui(P->Z, 1);

	/* dの決定 */
	if (X == NULL) {
		mpz_init(d);
		mpz_pow_ui(tmp,P->X,2); //tmp = x^2
		mpz_mod(tmp,tmp,N);
		mpz_pow_ui(tmp2,P->Y,2); //tmp2 = y^2
		mpz_mod(tmp2,tmp2,N);
		mpz_sub(d,tmp2,tmp); //d = y^2-x^2
		mpz_sub_ui(d,d,1); //d = y^2-x^2-1
		mpz_mul(tmp,tmp,tmp2); //tmp = x^2y^2
		mpz_mod(tmp,tmp,N);
		mpz_invert(tmp,tmp,N);
		mpz_mul(d,d,tmp); //dの値
		mpz_mod(d,d,N);
	}

	/* 素数の決定 */
	unsigned long int p = 2;
	mpz_t prime;
	mpz_init(prime);
	mpz_set_ui(prime, p);

	double stage1_time, stage2_time = -1;
	double start, end;
	start = omp_get_wtime();
	/* stage1 */
	while (p <= B1) {
		/* e = log p kを決める */
		e = (int)(log(B1) / log(p));
		for (i = 1; i <= e; i++) {
			/* Zを1にするための処理 */
			mpz_invert(inv, P->Z, N);
			mpz_mul(P->X, P->X, inv);
			mpz_mod(P->X, P->X, N);
			mpz_mul(P->Y, P->Y, inv);
			mpz_mod(P->Y, P->Y, N);
			mpz_set_ui(P->Z, 1);

			scalar(P, P, p, d, window_size, N);
			mpz_gcd(f, P->X, N);
			if (mpz_cmp_ui(f,1) != 0) {
				end = omp_get_wtime();
				stage1_time = end - start;
				goto FACTOR_FOUND;
			}
		}

		/* pを次の素数に */
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
		/* Zを1にするための処理 */
		mpz_invert(inv, P->Z, N);
		mpz_mul(P->X, P->X, inv);
		mpz_mod(P->X, P->X, N);
		mpz_mul(P->Y, P->Y, inv);
		mpz_mod(P->Y, P->Y, N);
		mpz_set_ui(P->Z, 1);

		scalar(P, P, p, d, window_size, N);
		mpz_mul(product, product, P->Z);
		mpz_mod(product, product, N);

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

	/* 使用変数・関数の開放 */
	projective_point_clear(P);
	mpz_clear(tmp);
	mpz_clear(tmp2);
	mpz_clear(inv);
	mpz_clear(prime);
}
