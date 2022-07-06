// hf3.c : universal functions for calculating tight-binding Hamiltonian of 3 band models

#include "hf3.h" 

void CalcK(Vector *v) {
	int i;

	for(i=0; i<K3; i++) {
		v[i].x = -M_PI + 2*M_PI * (i / (K*K)) / K;
		v[i].y = -M_PI + 2*M_PI * ((i / K) % K) / K;
		v[i].z = -M_PI + 2*M_PI * (i % K) / K;
	}
}	

void CalcTB(char *f_name, const int basis, int k_num, int tb_size, Vector *v, Vector q, lapack_complex_double *tb) {
	FILE *f;

	if((f = fopen(f_name, "wb")) == NULL) {
		printf("%s fopen FAIL\n", f_name);
		exit(1);
	}

	time_t t0 = time(NULL);

	int i;
	void (*Fourier)(const int, int, Vector, Vector, lapack_complex_double*);

	memset(tb, 0, tb_size); 

	if(strstr(f_name, "_f")) {
		Fourier = FourierF; 
	}
	else {
		Fourier = FourierA; 
	}

	for(i=0; i<k_num; i++) {
		Fourier(basis, i, v[i], q, tb);
	}
	
	fwrite(tb, tb_size, sizeof(lapack_complex_double), f);
	fclose(f);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, f_name, t1 - t0);
}

void CalcEigenTB(const int basis, lapack_complex_double *tbb, double *w, lapack_complex_double *v) {
	char jobz = 'V', uplo = 'L';
	double rwork[3*basis-2], w_tmp[basis];
	lapack_int ln = basis, lda = basis, lwork = 2*basis-1, info;
	lapack_complex_double work[lwork], v_tmp[H(basis)];
	
	int i, j;

	for(i=0; i<BAND; i++) {
		for(j=0; j<H(basis); j++) {
			v_tmp[j] = tbb[H(basis)*i + j];
		}

		LAPACK_zheev(&jobz, &uplo, &ln, v_tmp, &lda, w_tmp, work, &lwork, rwork, &info);
		if(info != 0) {
			printf("LAPACK_zheev FAIL\n");
			exit(1);
		}

		for(j=0; j<basis; j++) {
			w[basis*i + j] = w_tmp[j];
		}
		for(j=0; j<H(basis); j++) {
			v[H(basis)*i + j] = v_tmp[j];
		}	
	}
}

void MakeTB(char *type, const int basis, lapack_complex_double *tbb) {
	FILE *f;
	char f_name[16];

	sprintf(f_name, "band_%s.txt", type);

	if((f = fopen(f_name, "w")) == NULL) {
		printf("%s fopen FAIL\n", f_name);
		exit(1);
	}

	char jobz = 'V', uplo = 'L';
	double rwork[3*basis-2], w[basis];
	lapack_int ln = basis, lda = basis, lwork = 2*basis-1, info;
	lapack_complex_double work[lwork], v[H(basis)];
	
	int i, j;

	for(i=0; i<BAND; i++) {
		for(j=0; j<H(basis); j++) {
			v[j] = tbb[H(basis)*i + j];
		}

		LAPACK_zheev(&jobz, &uplo, &ln, v, &lda, w, work, &lwork, rwork, &info);
		if(info != 0) {
			printf("LAPACK_zheev FAIL\n");
			exit(1);
		}

		fprintf(f, "%4d", i);
		for(j=0; j<basis; j++) {
			fprintf(f, "%12f", w[j]);
		}
		fprintf(f, "\n");
	}

	fclose(f);
}
