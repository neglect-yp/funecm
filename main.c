/*
 * main.c
 *
 * 更新履歴
 * 2014/10/30 新規作成
 * 2014/11/01 誤字修正(null, unsinedなど)
 *            define修正
 * 2014/11/07 出力, loop処理追加
 * 2014/11/08 並列化実装
 * 2015/11/17 結果をファイルに出力 loop処理追加
 * 2015/12/01 メモリリーク修正
 */

#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include "point.h"

#define Y_LOOP 12000

/* !	適当です	!
 * 素因数が見つからなかった    : 0
 * エラー終了                  : 1
 * 因数が素数で余因数が素数    : 2
 * 因数が素数で余因数が合成数  : 4
 * 因数が合成数で余因数が素数  : 8
 * 因数が合成数で余因数が合成数: 16
 */
int main (int argc, char *argv[])
{
	if (argc <= 3) {
		fprintf (stderr, "Error: Need three Argument\n");
		fprintf (stderr, "Usage: funecm [options] <composite number> <k> <filename>\n");
		return 1;
	}

	/* オプション処理 */
	int opt;
	int loop = 0;
	while ((opt = getopt (argc, argv, "hl")) != -1) {
		switch (opt) {
			case 'h':
				fprintf(stdout, "Usage: funecm [options] <composite number> <k> <filename>\n");
				fprintf(stdout, "-h: help\n");
				return 0;
				break;
			case 'l':
				loop = 1;
				break;
			default:
				fprintf(stderr, "No such option\n");
				fprintf(stdout, "Usage: funecm [options] <composite number> <k> <filename>\n");
				return 1;
		}
	}

	/* ARGUMENT CONVERSION */
	mpz_t N;
	mpz_init_set_str(N, argv[optind++], 10);
	unsigned long int k;
	k = (unsigned long int)strtol(argv[optind++], NULL, 10);

	/* 修正予定 */
	if (k <= 2)
		return 0;

	switch (mpz_probab_prime_p (N, 25)) {
		case 2:
			gmp_printf("%Zd is definitely prime\n", N);
			return 0;
			break;
		case 1:
			gmp_printf("%Zd is probably prime\n", N);
			return 0;
			break;
		case 0:
			gmp_printf("%Zd is definitely composite\n", N);
			break;
		default:
			break;
	}

	/*
	AFFINE_POINT P;
	affine_point_init(P);
	*/
	
	/*コマンドライン引数からファイル名取得、拡張子が無ければ付ける*/
	FILE *fp;
	char f_name[64];
	strcpy(f_name,argv[optind++]);
	
	char ext[5]=".txt";
	if (strstr(f_name, ext) == NULL ) strcat(f_name, ext);
	
	if ((fp = fopen(f_name, "w")) == NULL) {
		fprintf(stderr,"file open error!!\n");
		return 1;
	}

	char digits[1000];
	
	double total_start;
	double total_end;
	unsigned long int Y;
	int found;
	
	do{
		gmp_fprintf(fp,"Input number: %Zd  ", N);
		mpz_get_str(digits, 10, N);
		fprintf(fp,"digits: %d\n", strlen(digits));
		gmp_fprintf(fp,"k: %ld\n", k);
		
		found = 0;
		total_start = omp_get_wtime();
		//int n = 1;
		int n = omp_get_max_threads();
		fprintf(fp,"threads = %d\n", n);
		//#pragma omp parallel num_threads(n) shared(found)
		#pragma omp parallel num_threads(n) shared(found)
		{
			//#pragma omp for private(factor)
			#pragma omp for
			for (Y = 2; Y < Y_LOOP; Y++) {
				printf("Y=%d\n",Y);
				mpz_t factor;
				mpz_init(factor);
				mpz_t cofactor;
				mpz_init(cofactor);
				double A_start;
				double A_end;
				if (found == 0) {
					A_start = omp_get_wtime();

					ecm(factor, N, Y, k, fp);
					mpz_divexact(cofactor, N, factor);
					/* 因数が1又はNだった場合係数を変えてやり直す */
					if (mpz_cmp_ui(factor, 1) == 0 || mpz_cmp(factor, N) == 0) {
						A_end = omp_get_wtime();
						fprintf(fp,"stage1 time: %.3lf seconds\nY=%d\nfactor not found\n--------------------------------------------------\n", (A_end - A_start),Y);
						
						/*
						printf("stage1 time: %.3lf seconds\n", (A_end - A_start));
						printf("factor not found\n");
						printf("--------------------------------------------------\n");
						*/
						mpz_clear(factor);
						mpz_clear(cofactor);
						continue;
					} else {
						found = 1;
					}
					mpz_get_str(digits, 10, factor);
					/* 終了ステータス */
					switch (mpz_probab_prime_p(factor, 25)) {
						case 2:
							gmp_fprintf(fp,"@ definite prime factor found: %Zd  digits: %d cofactor: %Zd\n", factor, strlen(digits), cofactor);
							break;
						case 1:
							gmp_fprintf(fp,"@ probable prime factor found: %Zd  digits: %d cofactor: %Zd\n", factor, strlen(digits), cofactor);
							break;
						case 0:
							gmp_fprintf(fp,"@ composite prime factor found: %Zd  digits: %d cofactor: %Zd\n", factor, strlen(digits), cofactor);
							break;
						default:
							break;
					}
					
					/*ループ用の初期化*/
					if(loop==1){
						mpz_set(N,cofactor);
						switch (mpz_probab_prime_p(cofactor, 25)) {
							case 2:
								gmp_fprintf(fp,"@ definite prime factor found: %Zd \n", cofactor);
								loop=0;
								break;
							case 1:
								loop=0;
								gmp_fprintf(fp,"@ probable prime factor found: %Zd \n", cofactor);
								break;
							case 0:
								gmp_fprintf(fp,"composite prime factor found: %Zd  \n", cofactor);
								break;
							default:
								break;
						}
						
					}
				}
				mpz_clear(factor);
				mpz_clear(cofactor);
			}
		   }
		total_end = omp_get_wtime();
		fprintf(fp,"total time: %.3lf seconds\n\n\n", (total_end - total_start));
	}while(loop==1&&found==1);
	
	/* メモリの解放*/
	/*affine_point_clear(P);*/
	fclose(fp);

	return 0;
}
