#define USE_MATH_DEFINES

#define BAND (100)
#define BASIS (2)
#define PRM (2)
#define H (BASIS*BASIS)

#define COMPLEX2(complex) (pow(creal(complex), 2) + pow(cimag(complex), 2))
#define DOT(v, k, g, r) ((creal(v) * cos((k+g) * r) + cimag(v) * sin((k+g) * r)) - (creal(v) * sin((k+g) * r) - cimag(v) * cos((k+g) * r)) * I)

#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lapacke.h>

void CalcK(double *k) {
	int i;
	
	for(i=0; i<BAND; i++) {
		k[i] = M_PI * i / BAND;
	}
}

void TightBinding1D(double *k, double *w, lapack_complex_double *v) {
	char jobz = 'V', uplo = 'L';
	double rwork[3*BASIS-2], w_tmp[BASIS];
	lapack_int ln = BASIS, lda = BASIS, lwork = 2*BASIS-1, info;
	lapack_complex_double work[lwork], v_tmp[H];
	
	int i, j;
	lapack_complex_double tbb[BAND*H];

	memset(tbb, 0, sizeof(tbb));

	for(i=0; i<BAND; i++) {
		tbb[H*i + 1] = -(1 + cos(2*k[i]) + sin(2*k[i])*I);
	}

	for(i=0; i<BAND; i++) {
		for(j=0; j<H; j++) {
			v_tmp[j] = tbb[H*i + j];
		}

		LAPACK_zheev(&jobz, &uplo, &ln, v_tmp, &lda, w_tmp, work, &lwork, rwork, &info);
		if(info != 0) {
			printf("LAPACK_zheev FAIL\n");
			exit(1);
		}

		for(j=0; j<BASIS; j++) {
			w[BASIS*i + j] = w_tmp[j];
		}
		for(j=0; j<H; j++) {
			v[H*i + j] = v_tmp[j];
		}	
	}
}

void MakeFold(double *w) {
	FILE *f;
	char f_name[16];

	sprintf(f_name, "band.txt");

	if((f = fopen(f_name, "w")) == NULL) {
		printf("%s fopen FAIL\n", f_name);
		exit(1);
	}

	int i, j;

	for(i=0; i<BAND; i++) {
		fprintf(f, "%4d", i);
		for(j=0; j<BASIS; j++) {
			fprintf(f, "%12f", w[BASIS*i + j]);
		}
		fprintf(f, "\n");
	}

	fclose(f);
}

void MakeUnfold(double *k, double *w, lapack_complex_double *v) {
	FILE *f;
	char f_name[16];

	sprintf(f_name, "band_unfold.txt");

	if((f = fopen(f_name, "w")) == NULL) {
		printf("%s fopen FAIL\n", f_name);
		exit(1);
	}

	int i, j, l;
	double r[2] = {0, 1}, g = 0, p2;
	lapack_complex_double p;

	//double r[2] = {0, 1}, g[2] = {0, M_PI}, p2[2];
	//lapack_complex_double p[2];

	for(i=0; i<BAND; i++) {
		fprintf(f, "%4d", i);
		
		for(j=0; j<BASIS; j++) {
			p = 0;

			for(l=0; l<PRM; l++) {
				p += DOT(v[BASIS*(BASIS*i + j) + l], k[i], g, r[l]);
			}
			p2 = COMPLEX2(p) / PRM;

			if(p2 > 0.5) {
				fprintf(f, "%12f", w[BASIS*i + j] * p2);
			}
		}
		fprintf(f, "\n");
	}

	fclose(f);
}

int main() {
	double k[BAND], w[BAND*BASIS];
	lapack_complex_double v[BAND*H];

	CalcK(k);

	TightBinding1D(k, w, v);

	//MakeFold(w);
	MakeUnfold(k, w, v);

	return 0;
}
