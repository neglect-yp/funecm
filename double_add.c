/*
 * double_add.c
 * Projective座標系の2倍算を行う
 *
 * R=P+Pをするメソッド
 * P:とある点
 * R:点Pを二倍算した点
 * a:楕円曲線の係数
 *
 * 更新履歴
 * 2014/10/17 新規作成
 * 2014/10/19 一部バグ修正
 *            mod処理を追加
 * 2014/10/30 コメント追加
 * 2014/11/01 宣言をマクロに変更
 * 2014/11/02 引数の一部をconstに変更
 *            メモリリーク修正
 * 2015/10/** 前年度より計算式を変更
 * 2016/03/23 mul処理の下にmod処理を追加
 * 2016/05/27 a=-1に変更
 */
#include <stdio.h>
#include "gmp.h"
#include "point.h"

/*
 * Projective座標系の2倍算を行う
 * PROJECTIVE_POINT R :Pの2倍算の結果を格納する点
 * PROJECTIVE_POINT P :2倍算を行う点
 * const mpz_t N            :mod N
 */
void double_add(PROJECTIVE_POINT R, PROJECTIVE_POINT P, const mpz_t N)
{
	mpz_t B,C,D,E,F,H,J;

	/* 初期化 */
	mpz_inits(B,C,D,E,F,H,J,NULL);

	mpz_add(B,P->X,P->Y); //B = X1+Y1
	mpz_pow_ui(B,B,2); //B = (X1+Y1)^2
	mpz_mod(B,B,N);
	mpz_pow_ui(C,P->X,2); //C = X1^2
	mpz_mod(C,C,N);
	mpz_pow_ui(D,P->Y,2); //D = Y1^2 
	mpz_mod(D,D,N);
	mpz_mul_si(E,C,-1); //E=-C
	mpz_add(F,E,D); //F = E+D
	mpz_pow_ui(H,P->Z,2); //H = Z1^2
	mpz_mod(H,H,N);
	mpz_sub(J,F,H); //J = F-H
	mpz_sub(J,J,H); //J = F-2H
	
	/*X3の計算*/
	mpz_sub(R->X,B,C); //X3=B-C
	mpz_sub(R->X,R->X,D); //X3=B-C-D
	mpz_mul(R->X,R->X,J); //X3
	mpz_mod(R->X,R->X,N);
	
	/*Y3の計算*/
	mpz_sub(R->Y,E,D); //Y3 = E-D
	mpz_mul(R->Y,R->Y,F); //Y3
	mpz_mod(R->Y,R->Y,N);

	/*Z3の計算*/
	mpz_mul(R->Z,F,J); //Z3
	mpz_mod(R->Z,R->Z,N);

	mpz_clears(B,C,D,E,F,H,J,NULL);
}

void dedicated_doubling(EXTENDED_POINT R, const EXTENDED_POINT P, const mpz_t N){
	
	mpz_t A,B,C,D,E,G,F,H;

	mpz_inits(A,B,C,D,E,G,F,H,NULL);

	mpz_pow_ui(A, P->X, 2);
	mpz_mod(A, A, N);

	mpz_pow_ui(B, P->Y, 2);
	mpz_mod(B, B, N);

	mpz_pow_ui(C, P->Z, 2);
	mpz_mod(C, C, N);
	mpz_mul_ui(C, C, 2);
	mpz_mod(C, C, N);
	
	mpz_mul_si(D, A, -1);

	mpz_add(E, P->X, P->Y);
	mpz_pow_ui(E, E, 2);
	mpz_mod(E, E, N);
	mpz_sub(E, E, A);
	mpz_sub(E, E, B);

	mpz_add(G, D, B);
	
	mpz_sub(F, G, C);
	
	mpz_sub(H, D, B);

	mpz_mul(R->X, E, F);
	mpz_mod(R->X, R->X, N);

	mpz_mul(R->Y, G, H);
	mpz_mod(R->Y, R->Y, N);
	
	mpz_mul(R->T, E, H);
	mpz_mod(R->T, R->T, N);

	mpz_mul(R->Z, F, G);
	mpz_mod(R->Z, R->Z, N);

	mpz_clears(A,B,C,D,E,G,F,H,NULL);
}
