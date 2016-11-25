#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include "point.h"

int main (int argc, char *argv[])
{
	int opt;
	int optc = 0;
	int loop = 0;
	int number_of_elliptic_curves = 4000;
	int window_size = 4;
	int atkin_flag = 0;
	while ((opt = getopt (argc, argv, "hlc:w:a")) != -1) {
		switch (opt) {
			case 'h':
				fprintf(stdout, "Usage: funecm [options] <composite number> <B1> <B2> <filename>\n");
				fprintf(stdout, "-h: help\n");
				fprintf(stdout, "-c <number>: number of elliptic curves\n");
				fprintf(stdout, "-w <bit>: window size\n");
				fprintf(stdout, "-a: use Atkin-Moraine ECPP\n");
				return 0;
				break;
			case 'l':
				loop = 1;
				break;
			case 'c':
				number_of_elliptic_curves = atoi(optarg);
				if (number_of_elliptic_curves == 0) {
					fprintf(stderr, "Error: the argument of 'c' option must be number\n");
					return 1;
				}
				break;
			case 'w':
				window_size = atoi(optarg);
				if (window_size <= 1 || 31 <= window_size) {
					fprintf(stderr, "Error: the argument of 'w' option must be [2-30]\n");
					return 1;
				}
				break;
			case 'a':
				atkin_flag = 1;
				break;
			default:
				fprintf(stderr, "No such option\n");
				fprintf(stdout, "Usage: funecm [options] <composite number> <B1> <B2> <filename>\n");
				return 1;
		}
		optc++;
	}

	if ((argc - optc) <= 3) {
		fprintf (stderr, "Error: Need three Argument\n");
		fprintf (stderr, "Usage: funecm [options] <composite number> <B1> <B2> <filename>\n");
		return 1;
	}

	mpz_t N;
	mpz_init_set_str(N, argv[optind++], 10);
	unsigned long int B1;
	unsigned long int B2;
	B1 = (unsigned long int)strtol(argv[optind++], NULL, 10);
	B2 = (unsigned long int)atol(argv[optind++]);
	if (B2 == 0)
		B2 = B1 * 100;

	if (B1 <= 2)
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
	int i;
	int found;

	gmp_fprintf(fp,"Input number: %Zd  ", N);
	mpz_get_str(digits, 10, N);
	fprintf(fp,"digits: %d\n", strlen(digits));
	gmp_fprintf(fp,"B1=%ld\n", B1);
	gmp_fprintf(fp,"B2=%ld\n", B2);
	fprintf(fp,"elliptic curves: %d\n", number_of_elliptic_curves);
	fprintf(fp,"window size: %d\n", window_size);

	found = 0;
	total_start = omp_get_wtime();
	int n = omp_get_max_threads();
	fprintf(fp,"threads = %d\n", n);

	/* variables for atkin-moraine ECPP */
	mpz_t s, t;
	mpz_inits(s, t, NULL);
	mpz_set_ui(s, 12);
	mpz_set_ui(t, 40);

	#pragma omp parallel num_threads(n) shared(found,s,t)
	{
		#pragma omp for
		for (i = 0; i < number_of_elliptic_curves; i++) {
			mpz_t X, Y, d;
			mpz_inits(X, Y, d, NULL);

			if (atkin_flag) {
				/* set X, Y, d by atkin-moraine ECPP */
				atkin_moraine(X, Y, d, s, t, N);
				gmp_printf("X=%Zd, ",X);
			} else {
				/* use random Y */
				gmp_randstate_t state;
				gmp_randinit_default(state);
				gmp_randseed_ui(state, (unsigned long int)time(NULL)+i);
				mpz_urandomm(Y, state, N);
				while (mpz_cmp_ui(Y, 2) < 0)
					mpz_add_ui(Y, Y, 1);
			} 

			gmp_printf("Y=%Zd\n",Y);

			mpz_t factor;
			mpz_init(factor);
			mpz_t cofactor;
			mpz_init(cofactor);
			double A_start;
			double A_end;
			if (found == 0) {
				A_start = omp_get_wtime();

				if (atkin_flag)
					ecm(factor, N, X, Y, d, B1, B2, fp, window_size);
				else
					ecm(factor, N, NULL, Y, d, B1, B2, fp, window_size);
				A_end = omp_get_wtime();

				/* output log to file */
				if (atkin_flag)
					gmp_fprintf(fp,"total time: %.3lf seconds\nX=%Zd\nY=%Zd\n--------------------------------------------------\n", (A_end - A_start),X,Y);
				else
					gmp_fprintf(fp,"total time: %.3lf seconds\nY=%Zd\n--------------------------------------------------\n", (A_end - A_start),Y);

				mpz_divexact(cofactor, N, factor);
				/* retry if factor is 1 or N */
				if (mpz_cmp_ui(factor, 1) == 0 || mpz_cmp(factor, N) == 0) {
					fprintf(fp, "factor not found\n\n");
					mpz_clears(X, Y, d, NULL);
					mpz_clear(factor);
					mpz_clear(cofactor);
					continue;
				} else {
					found = 1;
				}
				mpz_get_str(digits, 10, factor);

				switch (mpz_probab_prime_p(factor, 25)) {
					case 2:
						gmp_fprintf(fp,"@ definite prime factor found: %Zd  digits: %d cofactor: %Zd\n\n", factor, strlen(digits), cofactor);
						break;
					case 1:
						gmp_fprintf(fp,"@ probable prime factor found: %Zd  digits: %d cofactor: %Zd\n\n", factor, strlen(digits), cofactor);
						break;
					case 0:
						gmp_fprintf(fp,"@ composite prime factor found: %Zd  digits: %d cofactor: %Zd\n\n", factor, strlen(digits), cofactor);
						break;
					default:
						break;
				}
			}
			mpz_clears(X, Y, d, NULL);
			mpz_clear(factor);
			mpz_clear(cofactor);
		}
	}
	total_end = omp_get_wtime();
	fprintf(fp,"total time: %.3lf seconds\n\n\n", (total_end - total_start));
	mpz_clears(s, t, NULL);

	fclose(fp);

	return 0;
}
