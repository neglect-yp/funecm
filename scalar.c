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
static unsigned long int count_bit(unsigned long int n)
{
	unsigned long int count = 0;
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
	PROJECTIVE_POINT tP;
	PROJECTIVE_POINT tmp;
	EXTENDED_POINT eP;
	EXTENDED_POINT doubled;

	projective_point_init(tP);
	projective_point_init(tmp);
	extended_point_init(eP);
	extended_point_init(doubled);

	projective_point_set(tP, P);
	projective_point_set(tmp, P);
	protoext(eP, P, N);

	unsigned long int m = count_bit(k);
	char *bit = (char *)malloc(m);

	unsigned long int i = 0;
	for (i = 0; i < m; i++) {
		bit[i] = (k >> i) & 1;
	}

	i = m - 1;
	/* バイナリー法で計算を行う */
	while (i > 0) {
		i--;
		double_add(tP, tP, N);
		if (bit[i] == 1) {
			protoext(doubled, tP, N);
			extended_normal_add(doubled, doubled, eP, D, N);
			exttopro(tP, doubled, N);
		}
	}

	free(bit);

	projective_point_set(R, tP);

	/* メモリーの解放 */
	projective_point_clear(tP);
	projective_point_clear(tmp);
	extended_point_clear(eP);
	extended_point_clear(doubled);
}
