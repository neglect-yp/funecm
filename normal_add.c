 /*
 * 更新履歴
 * 2014/10/17 新規作成
 * 2014/10/19 バグ修正
 * 2014/11/01 各宣言をマクロに修正
 * 2014/11/02 引数の一部をconstに変更
 * 2015/10/** 前年度より計算式を変更
 * 2016/05/27 a=-1に変更 Extended Twisted Edwards Coodinatesに対応
 */
#include <stdio.h>
#include "gmp.h"
#include "point.h"

/* 
 * バイナリ法における加算公式を用いる部分では、Pを足すことになるが、P->Zが1のときは
 * mixed addition アルゴリズムを用いることができる。
 * R <- P + Q
 */

void extended_dedicated_add(EXTENDED_POINT R, EXTENDED_POINT P, EXTENDED_POINT Q, const mpz_t N) {
	mpz_t A,B,C,D,E,F,G,H,tmp;
	
	/*初期化*/
	mpz_inits(A,B,C,D,E,F,G,H,tmp,NULL);
	
	/* A<-(Y1-X1)*(Y2+X2) */
	mpz_sub(A, P->Y, P->X);
	mpz_add(tmp, Q->Y, Q->X);
	mpz_mul(A, A, tmp);
	mpz_mod(A, A, N);

	/* B<-(Y1+X1)*(Y2-X2) */
	mpz_add(B, P->Y, P->X);
	mpz_sub(tmp, Q->Y, Q->X);
	mpz_mul(B, B, tmp);
	mpz_mod(B, B, N);

	/* C<-2*Z1*T2 */
	mpz_mul_ui(C, P->Z, 2);
	mpz_mul(C, C, Q->T);
	mpz_mod(C, C, N);

	/* D<-2*T1*Z2 */
	mpz_mul_ui(D, P->T, 2);
	//mpz_mul(D, D, Q->Z); //Z2 = 1
	mpz_mod(D, D, N);
	
	/* E<-D+C */
	mpz_add(E, D, C);
	mpz_mod(E, E, N);
	
	/* F<-B-A */
	mpz_sub(F, B, A);
	mpz_mod(F, F, N);

	/* G<-B+A */
	mpz_add(G, B, A);
	mpz_mod(G, G, N);
	
	/* H<-D-C */
	mpz_sub(H, D, C);
	mpz_mod(H, H, N);

	/* X3<-E*F */
	mpz_mul(R->X, E, F);
	mpz_mod(R->X, R->X, N);
	
	/* Y3<-G*H */
	mpz_mul(R->Y, G, H);
	mpz_mod(R->Y, R->Y, N);

	/* T3<-E*H */
	//mpz_mul(R->T, E, H);
	//mpz_mod(R->T, R->T, N);

	/* Z3<-F*G */
	mpz_mul(R->Z, F, G);
	mpz_mod(R->Z, R->Z, N);

	mpz_clears(A,B,C,D,E,F,G,H,tmp,NULL);
}
