// lib/com.c : common functions

#include "hf.h" 

void ReadConfig(Config *c) {
	char fn[256];
	sprintf(fn, "input/%s/config_%c.txt", c->strain, c->type[0]);

	FILE *f0 = fopen("input/config.txt", "r"), *f = fopen(fn, "r");
	int i;
	char buf[1024];

	while(!feof(f0)) {
		fgets(buf, sizeof(buf), f0);
		if     (strstr(buf, "Nkg1")) sscanf(buf, "Nkg1 %d", &c->Nkg1);
		else if(strstr(buf, "Nkb"))  sscanf(buf, "Nkb %d", &c->Nkb);
	}
	fclose(f0);

	c->Nkg = pow(c->Nkg1, 3);

	while(!feof(f)) {
		fgets(buf, sizeof(buf), f);
		if     (strstr(buf, "Ni"))  sscanf(buf, "Ni %d", &c->Ni);
		else if(strstr(buf, "Q"))   sscanf(buf, "Q %lf%lf%lf", &c->Q[0], &c->Q[1], &c->Q[2]);
		else if(strstr(buf, "Lat")) sscanf(buf, "Lat %s", c->lat);
		else if(strstr(buf, "Sym")) sscanf(buf, "Sym %s", c->sym);
	}
	fclose(f);

	c->Ns = c->Ni * Nc;
	c->Nb = c->Ni * Nc * 2;
	for(i=0; i<DIM; i++) c->Q[i] *= M_PI;
}

void CalcEigen(Config c, double *ev, lapack_complex_double *es) {
	char jobz = 'V', uplo = 'L';
	double rwork[3*c.Nb-2];
	lapack_int ln = c.Nb, lda = c.Nb, lwork = 2*c.Nb-1, info;
	lapack_complex_double work[lwork];

	LAPACK_zheev(&jobz, &uplo, &ln, es, &lda, ev, work, &lwork, rwork, &info);
	if(info != 0) {
		printf("LAPACK_zheev FAIL\n");
		exit(1);
	}
}

void ShowProgress(int i, int i_max) {
	int x, ticks=50;

	printf("\r[");
	for(x=0; x<ticks*(i/i_max); x++) printf("=");
	for(   ; x<ticks;           x++) printf(" ");
	printf("] %10d/%d", i+1, i_max);
	fflush(stdout);
}

void ReplaceStr(char *in, char *org, char *rep, char *out) {
	int n;
	char *tmp;

	tmp = strstr(in, org);
	n = tmp - in;

	strncpy(out, in, n);  
	sprintf(out+n, "%s%s", rep, tmp+strlen(org));
}

void GetDimsH5(char *fn, char *dn, hsize_t *dims) {
	hid_t file_id, dataset_id, dataspace_id;

	file_id      = H5Fopen(fn, H5F_ACC_RDONLY, H5P_DEFAULT); 
	dataset_id   = H5Dopen(file_id, dn, H5P_DEFAULT);
	dataspace_id = H5Dget_space(dataset_id);

	H5Sget_simple_extent_dims(dataspace_id, dims, NULL);

	H5Fclose(file_id);
	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);
}

void ReadH5(char *fn, char *dn, double *val) {
	hid_t file_id, dataset_id;

	file_id    = H5Fopen(fn, H5F_ACC_RDONLY, H5P_DEFAULT); 
	dataset_id = H5Dopen(file_id, dn, H5P_DEFAULT);

	H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, val);

	H5Fclose(file_id);
	H5Dclose(dataset_id);
}

void WriteH5(char *fn, char *dn, int dim, hsize_t *dims, double *val) {
	hid_t file_id, dataset_id, dataspace_id;

	file_id      = H5Fcreate(fn, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	dataspace_id = H5Screate_simple(dim, dims, NULL);
	dataset_id   = H5Dcreate(file_id, dn, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, val);

	H5Fclose(file_id);
	H5Sclose(dataspace_id);
	H5Dclose(dataset_id);
}
