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

void CalcTB(char *fo_name, const int basis, const int k_num, const int tb_size, Vector *v, lapack_complex_double *tb) {
	FILE *fi, *fo;
	char fi_name[32];

	sprintf(fi_name, "lattice.txt");

	if((fi = fopen(fi_name, "r")) == NULL) {
		printf("%s fopen FAIL\n", fi_name);
		exit(1);
	}
	if((fo = fopen(fo_name, "wb")) == NULL) {
		printf("%s fopen FAIL\n", fo_name);
		exit(1);
	}

	time_t t0 = time(NULL);

	int i;
	void (*Fourier)(FILE*, const int, const int, Vector, Vector, lapack_complex_double*);
	Vector q;

	memset(tb, 0, tb_size); 

	if(strstr(fo_name, "_f")) {
		Fourier = FourierF; 
		q.x = 0;
		q.y = 0;
		q.z = 0;
	}
	else if(strstr(fo_name, "_a")) {
		Fourier = FourierA; 
		q.x = M_PI;
		q.y = 0;
		q.z = 0;
	}
	else if(strstr(fo_name, "_c")) {
		Fourier = FourierA; 
		q.x = M_PI;
		q.y = M_PI;
		q.z = 0;
	}
	else if(strstr(fo_name, "_g")) {
		Fourier = FourierA; 
		q.x = M_PI;
		q.y = M_PI;
		q.z = M_PI;
	}
	else {
		printf("\"%s\" is not available\n", fo_name);
		exit(1);
	}

	for(i=0; i<k_num; i++) {
		Fourier(fi, i, basis, v[i], q, tb);
		rewind(fi);
	}
	
	fwrite(tb, tb_size, sizeof(lapack_complex_double), fo);

	fclose(fi);
	fclose(fo);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fo_name, t1 - t0);
}

void MakeTB(char *type, int basis, lapack_complex_double *tbb) {
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
