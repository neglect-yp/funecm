/*
 * scalar.c
 * スカラー倍を行う
 *
 * !CAUTION!
 * kが１, 2の時未定義
 * 更新履歴
 * 2014/ 6/16 新規作成					// !?
 * 2014/10/29 ｋの値による分岐を削除
 * 2014/11/01 ローカルに点を取るように変更
 * 2014/11/02 引数の一部をconstに変更
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "point.h"

/* DEBUG*/
void print_bit(unsigned long int n)
{
	while (n != 0) {
		printf("%d", n & 1);
		n >>= 1;
	}
	printf("\n");
}

/*
 * ビット数をカウントする関数
 * unsigned long int n :ビット数をカウントする
 */
static long int count_bit(unsigned long int n)
{
	long int count = 0;
	while (n != 0) {
		n >>= 1;
		count++;
	}
	return count;
}

/* PROJECTIVE_POINT R  :スカラー倍の計算を格納する点
 * PROJECTIVE_POINT P  :スカラー倍を行う点
 * const unsigned long int k :スカラー倍k
 * const mpz_t D             :楕円曲線の係数
 * const mpz_t N             :mod N
 * */
void scalar(PROJECTIVE_POINT R, PROJECTIVE_POINT P, const unsigned long int k, const mpz_t D, const mpz_t N)
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

	/* 移動窓法のための事前計算 */	
	EXTENDED_POINT *Parray = (EXTENDED_POINT *)malloc(16*sizeof(EXTENDED_POINT));
	extended_point_init(Parray[1]);
	extended_point_set(Parray[1], eP);
	extended_point_init(Parray[2]);
	dedicated_doubling(Parray[2], eP, N);
	mpz_t inv;
	mpz_init(inv);
	mpz_invert(inv, Parray[2]->Z, N);
	mpz_mul(Parray[2]->X, Parray[2]->X, inv);
	mpz_mod(Parray[2]->X, Parray[2]->X, N);
	mpz_mul(Parray[2]->Y, Parray[2]->Y, inv);
	mpz_mod(Parray[2]->Y, Parray[2]->Y, N);
	mpz_mul(Parray[2]->T, Parray[2]->X, Parray[2]->Y);
	mpz_mod(Parray[2]->T, Parray[2]->T, N);
	mpz_set_ui(Parray[2]->Z, 1);
	for (i = 1; i <= 7; i++) {
		extended_point_init(Parray[2*i+1]);
		extended_dedicated_add(Parray[2*i+1],Parray[2*i-1],Parray[2],N);
		mpz_invert(inv, Parray[2*i+1]->Z, N);
		mpz_mul(Parray[2*i+1]->X, Parray[2*i+1]->X, inv);
		mpz_mod(Parray[2*i+1]->X, Parray[2*i+1]->X, N);
		mpz_mul(Parray[2*i+1]->Y, Parray[2*i+1]->Y, inv);
		mpz_mod(Parray[2*i+1]->Y, Parray[2*i+1]->Y, N);
		mpz_mul(Parray[2*i+1]->T, Parray[2*i+1]->X, Parray[2*i+1]->Y);
		mpz_mod(Parray[2*i+1]->T, Parray[2*i+1]->T, N);
		mpz_set_ui(Parray[2*i+1]->Z, 1);
	}
	mpz_clear(inv);

	i = m - 1;
	/* バイナリー法で計算を行う */
	/*
	while (i > 0) {
		i--;
		dedicated_doubling(tP, tP, N);
		if (bit[i] == 1) {
			extended_dedicated_add(tP, tP, eP, N);
		}
	}
	*/
	while (i > 0) {
		if (!bit[i]) {
			dedicated_doubling(tP, tP, N);
			i--;
		} else {
			long int t = ((i-3 > 0) ? i-3 : 0);
			while (!bit[t])
				t++;
			int h = 0;
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

	free(bit);

	exttopro(R, tP, N);

	/* メモリーの解放 */
	extended_point_clear(tP);
	extended_point_clear(eP);
	extended_point_clear(Parray[1]);
	extended_point_clear(Parray[2]);
	for (i = 1; i <= 7; i++)
		extended_point_clear(Parray[2*i+1]);
	free(Parray);
}
