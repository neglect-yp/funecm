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
	EXTENDED_POINT tP;
	extended_point_init(tP);
	protoext(tP, P, N);

	EXTENDED_POINT eP;
	extended_point_init(eP);
	extended_point_set(eP, tP);

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
		dedicated_doubling(tP, tP, N);
		if (bit[i] == 1) {
			extended_dedicated_add(tP, tP, eP, N);
		}
	}

	free(bit);

	exttopro(R, tP, N);

	/* メモリーの解放 */
	extended_point_clear(tP);
	extended_point_clear(eP);
}
