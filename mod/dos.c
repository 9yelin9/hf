// mod/dos.c : generate dos dataset

#define Nbs      12   // num of bases
#define Nhsp     4    // num of high symmetry points
#define BINS_MAX 8193 // num of energy bins

#define E_MIN -8.1
#define E_MAX  8.1

#include "hf3.h" 

typedef struct BandFeatures {
	double e[Nbs * Nhsp];
	double w[Nbs * Nhsp];
} BFeats;

void GetDimsH5(char *fn, char *dn, hsize_t *dims) {
	hid_t f, f_dset, f_size;

	f      = H5Fopen(fn, H5F_ACC_RDONLY, H5P_DEFAULT); 
	f_dset = H5Dopen(f, dn, H5P_DEFAULT);
	f_size = H5Dget_space(f_dset);

	H5Sget_simple_extent_dims(f_size, dims, NULL);

	H5Fclose(f);
	H5Dclose(f_dset);
	H5Sclose(f_size);
}

void ReadH5(char *fn, char *dn, double *val) {
	hid_t f, f_dset;

	f      = H5Fopen(fn, H5F_ACC_RDONLY, H5P_DEFAULT); 
	f_dset = H5Dopen(f, dn, H5P_DEFAULT);

	H5Dread(f_dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, val);

	H5Fclose(f);
	H5Dclose(f_dset);
}

void WriteH5(char *fn, char *dn, int dim, hsize_t *dims, double *val) {
	hid_t f, f_size, f_dset;

	f      = H5Fcreate(fn, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	f_size = H5Screate_simple(dim, dims, NULL);
	f_dset = H5Dcreate(f, dn, H5T_NATIVE_DOUBLE, f_size, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	H5Dwrite(f_dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, val);

	H5Fclose(f);
	H5Sclose(f_size);
	H5Dclose(f_dset);
}

void GenEnergy(char *dir) {
	char fn[128];
	sprintf(fn, "%s/energy.h5", dir);

	time_t t0 = time(NULL);

	int i;
	double e[BINS_MAX];

	for(i=0; i<BINS_MAX; i++) e[i] = E_MIN + (E_MAX - E_MIN) * i / BINS_MAX;

	hsize_t dims[1] = {BINS_MAX};
	WriteH5(fn, "/energy", 1, dims, e);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fn, t1 - t0);
}

void DOS0(BFeats b, char *dtype, double eta, double *e, double *dos) { // Local DOS
}

#define D_D_IDX (BINS_MAX*hsp + i)
#define D_B_IDX (Nbs*hsp + j)

void DOS1(BFeats b, char *dtype, double eta, double *e, double *dos) { // k-projected DOS at high symmetry points
	int i, j, hsp;
	double weight[BINS_MAX], broad[BINS_MAX], sum;

	// options
	if(strstr(dtype, "f")) { // ignore energy over fermi level
		for(i=0; i<BINS_MAX; i++) {
			if(e[i] > 0) weight[i] = 0;
			else         weight[i] = 1;
		}
	}
	else for(i=0; i<BINS_MAX; i++) weight[i] = 1;

	if(strstr(dtype, "b")) { // add linear broadening regard to distance from fermi level
		for(i=0; i<BINS_MAX; i++) {
			broad[i] = eta + fabs(e[i]) / 10;
		}
	}
	else for(i=0; i<BINS_MAX; i++) broad[i] = eta;

	for(hsp=0; hsp<Nhsp; hsp++) {
		for(i=0; i<BINS_MAX; i++) {
			sum = 0;
			for(j=0; j<Nbs; j++) sum += (broad[i] / (pow(e[i] - b.e[D_B_IDX], 2) + pow(broad[i], 2))) * b.w[D_B_IDX] * weight[i];
			dos[D_D_IDX] = sum / M_PI;
		}
	}
}

void GenDOS(char *dir, char *dtype, double eta, double *e) {
	FILE *fb;
	char fbn[128], fdn[128];

	sprintf(fbn, "%s/band.txt", dir);
	sprintf(fdn, "%s/dos_dt%s_eta%.2f.h5", dir, dtype, eta);
	fb = OpenFile(fbn, "r");

	time_t t0 = time(NULL);

	int i, j, fb_len;
	char fb_len_c[16], tmp[2*(Nbs * Nhsp)][16];
	BFeats b;

	fgets(fb_len_c, sizeof(fb_len_c), fb);
	fb_len = atoi(fb_len_c);
	double dos[fb_len][BINS_MAX * Nhsp];

	void (*DOS)(BFeats, char*, double, double*, double*);
	if     (strstr(dtype, "0")) DOS = DOS0; // Local DOS
	else if(strstr(dtype, "1")) DOS = DOS1; // k-projected DOS at high symmetry points
	else {
		printf("'%s' is wrong dtype\n", dtype);
		exit(1);
	}

	for(i=0; i<2*(Nbs * Nhsp); i++) fscanf(fb, "%s", tmp[i]); // ignore first line

	for(i=0; i<fb_len; i++) {
		for(j=0; j<Nbs * Nhsp; j++) fscanf(fb, "%lf", &b.e[j]);
		for(j=0; j<Nbs * Nhsp; j++) fscanf(fb, "%lf", &b.w[j]);

		DOS(b, dtype, eta, e, dos[i]);
	}	
	fclose(fb);

	hsize_t dims[2] = {fb_len, BINS_MAX * Nhsp};
	WriteH5(fdn, "/dos", 2, dims, (double*)dos);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fdn, t1 - t0);
}

#define R_RD_IDX (bins*hsp + j)
#define R_D_IDX  (BINS_MAX*hsp + Nclus*j + k)
#define R_E_IDX  (Nclus*j + k)

void ReduceDOS(char *dir, char *dtype, double eta, double *e, int bins) {
	char fdn[128], frn[128];

	sprintf(fdn, "%s/dos_dt%s_eta%.2f.h5", dir, dtype, eta);
	sprintf(frn, "%s/rdos_dt%s_eta%.2f_bins%d.h5", dir, dtype, eta, bins);

	hsize_t fd_dims[2];
	GetDimsH5(fdn, "/dos", fd_dims);

	double dos[fd_dims[0]][fd_dims[1]];
	ReadH5(fdn, "/dos", (double*)dos);

	time_t t0 = time(NULL);

	int i, j, k, hsp, Nclus = BINS_MAX / bins;
	double rdos[fd_dims[0]][bins * Nhsp];

	memset(rdos, 0, sizeof(rdos));

	for(i=0; i<fd_dims[0]; i++) {
		for(hsp=0; hsp<Nhsp; hsp++) {
			for(j=0; j<bins; j++) {
				for(k=0; k<Nclus; k++) {
					rdos[i][R_RD_IDX] += dos[i][R_D_IDX] * (e[R_E_IDX+1] - e[R_E_IDX]);
				}
			}
		}
	}

	hsize_t fr_dims[2] = {fd_dims[0], bins * Nhsp};
	WriteH5(frn, "/rdos", 2, fr_dims, (double*)rdos);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, frn, t1 - t0);
}

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("%s <energy/dos/rdos> <(d/r)dtype> <(d/r)eta> <(r)bins> : generate dos dataset\n", argv[0]);
		exit(1);
	}
	omp_set_num_threads(1); 
	char *dir = "output/baoso3_ms/magstr";

	if(strstr(argv[1], "e")) GenEnergy(dir);
	else {
		char *dtype = argv[2];
		double eta  = atof(argv[3]);

		char fen[128];
		sprintf(fen, "%s/energy.h5", dir);

		hsize_t dims[1];
		GetDimsH5(fen, "/energy", dims);

		double e[dims[0]];
		ReadH5(fen, "/energy", e);

		if(strstr(argv[1], "r")) ReduceDOS(dir, dtype, eta, e, atoi(argv[4]));
		else                     GenDOS(dir, dtype, eta, e);
	}

	return 0;
}
