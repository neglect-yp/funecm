 /*
 * 更新履歴
 * 2014/10/17 新規作成
 * 2014/10/19 バグ修正
 * 2014/11/01 各宣言をマクロに修正
 * 2014/11/02 引数の一部をconstに変更
 * 2015/10/** 前年度より計算式を変更
 */

#include "gmp.h"
#include "point.h"

/*
 * 異なる点の加算を行う関数
 * PROJECTIVE_POINT R:P + Qの結果を格納する点
 * PROJECTIVE_POINT P:演算される点
 * PROJECTIVE_POINT Q:演算される点
 *const mpz_t d:楕円曲線の係数
 * const mpz_t N           :mod N
 */
void normal_add(PROJECTIVE_POINT R, PROJECTIVE_POINT P, PROJECTIVE_POINT Q, const mpz_t d,const mpz_t N)
{
	PROJECTIVE_POINT tP;
	PROJECTIVE_POINT tQ;
	mpz_t A;
	mpz_t B;
	mpz_t C;
	mpz_t D;
	mpz_t F;
	mpz_t E;
	mpz_t G;
	mpz_t H;
	mpz_t I;

	/*初期化*/
	mpz_init(A);
	mpz_init(B);
	mpz_init(C);
	mpz_init(D);
	mpz_init(E);
	mpz_init(F);
	mpz_init(G);
	mpz_init(H);
	mpz_init(I);

	projective_point_init(tP);
	projective_point_init(tQ);

	projective_point_set(tP, P);
	projective_point_set(tQ, Q);

	mpz_mul(A,tP->Z,tQ->Z); //A = Z1Z2
	mpz_mod(A,A,N);
	mpz_pow_ui(B,A,2); //B = A^2
	mpz_mod(B,B,N);
	/*mpz_mul(B,B,d); B = dA^2*/
	mpz_mul(C,tP->X,tQ->X); //C = X1X2
	mpz_mod(C,C,N);
	mpz_mul(D,tP->Y,tQ->Y); //D = Y1Y2
	mpz_mod(D,D,N);
	mpz_mul(E,C,D); //E = CD
	mpz_mod(E,E,N);
	mpz_mul(E,E,d); //E = dCD
	mpz_mod(E,E,N);
	mpz_sub(F,B,E); //F = B-E
	mpz_add(G,B,E); //F = B+E
	mpz_add(H,tP->X,tP->Y); //H =X1+Y1
	mpz_add(I,tQ->X,tQ->Y); //I =X2+Y2

	/*X3の計算*/
	mpz_mul(R->X,H,I); //X3 = HI
	mpz_mod(R->X,R->X,N);
	mpz_sub(R->X,R->X,C); //X3 = HI-C
	mpz_sub(R->X,R->X,D); //X3 = HI-C-D
	mpz_mul(R->X,R->X,A); //X3 = A(HI-C-D)
	mpz_mod(R->X,R->X,N);
	mpz_mul(R->X,R->X,F); //X3
	mpz_mod(R->X,R->X,N);

	/*Y3の計算*/
	mpz_sub(R->Y,D,C); //Y3 = D-C
	mpz_mul(R->Y,R->Y,A); //Y3 = A(D-C)
	mpz_mod(R->Y,R->Y,N);
	mpz_mul(R->Y,R->Y,G); //Y3
	mpz_mod(R->Y,R->Y,N);

	/*Z3の計算*/
	mpz_mul(R->Z,F,G); //Z3
	mpz_mod(R->Z,R->Z,N);

	/*使用した変数の開放*/
	mpz_clear(A);
	mpz_clear(B);
	mpz_clear(C);
	mpz_clear(D);
	mpz_clear(E);
	mpz_clear(F);
	mpz_clear(H);
	mpz_clear(G);
	mpz_clear(I);

	projective_point_clear(tP);
	projective_point_clear(tQ);
}
