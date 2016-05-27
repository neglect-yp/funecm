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
#include "gmp.h"
#include "point.h"

/*
 * Projective座標系の2倍算を行う
 * PROJECTIVE_POINT R :Pの2倍算の結果を格納する点
 * PROJECTIVE_POINT P :2倍算を行う点
 * const mpz_t a            :楕円曲線の係数
 * const mpz_t N            :mod N
 */
void double_add(PROJECTIVE_POINT R, PROJECTIVE_POINT P, /*const mpz_t a,*/ const mpz_t N)
{
	mpz_t B;
	mpz_t C;
	mpz_t D;
	mpz_t E;
	mpz_t F;
	mpz_t H;
	mpz_t J;

	PROJECTIVE_POINT tP; //Pの座標を格納する

	/* 初期化 */
	mpz_init(B);
	mpz_init(C);
	mpz_init(D);
	mpz_init(E);
	mpz_init(F);
	mpz_init(H);
	mpz_init(J);

	projective_point_init(tP);

	projective_point_set(tP, P);

	mpz_add(B,tP->X,tP->Y); //B = X1+Y1
	mpz_pow_ui(B,B,2); //B = (X1+Y1)^2
	mpz_mod(B,B,N);
	mpz_pow_ui(C,tP->X,2); //C = X1^2
	mpz_mod(C,C,N);
	mpz_pow_ui(D,tP->Y,2); //D = Y1^2 
	mpz_mod(D,D,N);
	mpz_mul_ui(E,C,-1); //E=-C
	mpz_add(F,C,D); //F = E+D
	mpz_pow_ui(H,tP->Z,2); //H = Z1^2
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

	mpz_clear(B);
	mpz_clear(C);
	mpz_clear(D);
	mpz_clear(E);
	mpz_clear(F);
	mpz_clear(H);
	mpz_clear(J);

	projective_point_clear(tP);
}

