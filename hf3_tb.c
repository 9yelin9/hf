// hf3.c : universal functions for calculating tight-binding Hamiltonian of 3 band models

#include "hf3.h" 

void ReadPath(Vector *path_k, Vector *path_b) {
	FILE *fk, *fb;
	char fks[32], fbs[32];

	sprintf(fks, "path_k%d.bin", K);
	sprintf(fbs, "path_b.bin");

	if((fk = fopen(fks, "rb")) == NULL) {
		printf("%s fopen FAIL\n", fks);
		exit(1);
	}
	if((fb = fopen(fbs, "rb")) == NULL) {
		printf("%s fopen FAIL\n", fbs);
		exit(1);
	}

	fread(path_k, K3, sizeof(Vector), fk);
	fread(path_b, BAND, sizeof(Vector), fb);

	fclose(fk);
	fclose(fb);
}

void CalcPathK() {
	FILE *f;
	char fs[32];

	sprintf(fs, "path_k%d.bin", K);

	if((f = fopen(fs, "wb")) == NULL) {
		printf("%s fopen FAIL\n", fs);
		exit(1);
	}

	int i;
	Vector path[K3];

	for(i=0; i<K3; i++) {
		path[i].x = -M_PI + 2*M_PI * (i / (K*K)) / K;
		path[i].y = -M_PI + 2*M_PI * ((i / K) % K) / K;
		path[i].z = -M_PI + 2*M_PI * (i % K) / K;
	}

	fwrite(path, sizeof(path), 1, f);
	fclose(f);
}	

void CalcEigenTB(char *type, double *w, lapack_complex_double *v) {
	FILE *f;
	char fs[32];

	sprintf(fs, "tbb_%s.bin", type);

	if((f = fopen(fs, "rb")) == NULL) {
		printf("%s fopen FAIL\n", fs);
		exit(1);
	}

	const int basis = strstr(type, "f") ? BASIS1 : BASIS2;

	lapack_complex_double tb[HB(basis)];
	fread(tb, sizeof(tb), 1, f);
	fclose(f);

	char jobz = 'V', uplo = 'L';
	double rwork[3*basis-2], w_tmp[basis];
	lapack_int ln = basis, lda = basis, lwork = 2*basis-1, info;
	lapack_complex_double work[lwork], v_tmp[H(basis)];
	
	int i, j;

	for(i=0; i<BAND; i++) {
		for(j=0; j<H(basis); j++) v_tmp[j] = tb[H(basis)*i + j];

		LAPACK_zheev(&jobz, &uplo, &ln, v_tmp, &lda, w_tmp, work, &lwork, rwork, &info);
		if(info != 0) {
			printf("LAPACK_zheev FAIL\n");
			exit(1);
		}

		for(j=0; j<basis; j++) w[basis*i + j] = w_tmp[j];
		for(j=0; j<H(basis); j++) v[H(basis)*i + j] = v_tmp[j];
	}
}

void MakeTB(char *type, char *fs, Vector *path, Vector q) {
	FILE *f;

	if((f = fopen(fs, "wb")) == NULL) {
		printf("%s fopen FAIL\n", fs);
		exit(1);
	}

	time_t t0 = time(NULL);

	const int basis = strstr(type, "f") ? BASIS1 : BASIS2;
	const int path_num = strstr(fs, "tbk") ? K3 : BAND;

	int i, j;
	lapack_complex_double tb[path_num * H(basis)];
	void (*Fourier)(char*, int, Vector, Vector, lapack_complex_double*);

	if(strstr(type, "f")) Fourier = FourierF; 
	else if(strstr(type, "s")) Fourier = FourierSubA; 
	else Fourier = FourierA;

	memset(tb, 0, sizeof(tb));

	for(i=0; i<path_num; i++) Fourier(type, i, path[i], q, tb);
	
	fwrite(tb, sizeof(tb), 1, f);
	fclose(f);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fs, t1 - t0);

	if(strstr(fs, "tbk")) return;
	else {
		FILE *fb;
		char fbs[32];

		sprintf(fbs, "band_%s.txt", type);

		if((fb = fopen(fbs, "w")) == NULL) {
			printf("%s fopen FAIL\n", fbs);
			exit(1);
		}

		double w[path_num * basis];
		lapack_complex_double v[path_num * H(basis)];

		CalcEigenTB(type, w, v);

		for(i=0; i<BAND; i++) {
			fprintf(fb, "%4d", i);
			for(j=0; j<basis; j++) fprintf(fb, "%12f", w[basis*i + j]);
			fprintf(fb, "\n");
		}

		fclose(fb);
	}
}

