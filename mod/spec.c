// mod/spec.c : generate spec dataset

#define BINS_MAX 1024 // maximum bins
#define Nkb      1024 // num of kpoints
#define Nbs      12   // num of bases
#define Nhsp     4    // num of high symmetry points

#include "hf3.h" 

void ShowProgress(int i, int i_max) {
	int x, ticks=50;

	printf("\r[");
	for(x=0; x<i/(i_max/ticks); x++) printf("=");
	for(   ; x<ticks;           x++) printf(" ");
	printf("] %10d/%d", i+1, i_max);
	fflush(stdout);
}

void GenDatabase(char *save, char *dtype) {
	FILE *fb, *fu;
	char fbn[256], fun[256], fen[256], fwn[256];
	
	sprintf(fbn, "%s/list_b_%c.txt", save, dtype[0]);
	sprintf(fun, "%s/list_u_%c.txt", save, dtype[0]);
	sprintf(fen, "%s/band_e_%c.h5",  save, dtype[0]);
	sprintf(fwn, "%s/band_w_%c.h5",  save, dtype[0]);
	fb = OpenFile(fbn, "r");
	fu = OpenFile(fun, "r");

	time_t t0 = time(NULL);

	FILE *b, *u;
	int i, j, k, f_len;
	char buf[1024], bn[256], un[256];
	double fermi;

	fgets(buf, sizeof(buf), fb);
	fgets(buf, sizeof(buf), fu);
	f_len = atoi(buf);
	double band_e[f_len][Nkb * Nbs], band_w[f_len][Nkb * Nbs];

	for(i=0; i<f_len; i++) {
		ShowProgress(i, f_len);

		fscanf(fb, "%s%lf", bn, &fermi);
		fscanf(fu, "%s",    un);

		b = OpenFile(bn, "r");
		u = OpenFile(un, "r");

		fgets(buf, sizeof(buf), b);
		fgets(buf, sizeof(buf), u);

		for(j=0; j<Nkb; j++) {
			for(k=0; k<Nbs; k++) {
				fscanf(b, "%lf", &band_e[i][Nbs*j + k]);
				fscanf(u, "%lf", &band_w[i][Nbs*j + k]);
				band_e[i][Nbs*j + k] -= fermi;
			}
		}

		fclose(b);
		fclose(u);
	}

	fclose(fb);
	fclose(fu);

	hsize_t dims[2] = {f_len, Nkb * Nbs};
	WriteH5(fen, "/energy", 2, dims, (double*)band_e);
	WriteH5(fwn, "/weight", 2, dims, (double*)band_w);

	time_t t1 = time(NULL);
	printf("\n%s(%s, %s) : %lds\n", __func__, fen, fwn, t1 - t0);
}

void Spec(double *band_e, double *band_w, double *e, double *broad, double *spec) {
	int i, j, k;
	double sum;

	for(i=0; i<BINS_MAX; i++) {
		for(j=0; j<Nkb; j++) {
			sum = 0;
			for(k=0; k<Nbs; k++) {
				sum += (broad[i] / (pow(e[i] - band_e[Nbs*j + k], 2) + pow(broad[i], 2))) * band_w[Nbs*j + k];
			}
			spec[Nkb*i + j] = sum / M_PI;
		}
	}
}

void GenSpec(char *save, char *dtype, double eta, double *e) {
	char fen[256], fwn[256], fsn[256];
	
	sprintf(fen, "%s/band_e_%c.h5", save, dtype[0]);
	sprintf(fwn, "%s/band_w_%c.h5", save, dtype[0]);
	sprintf(fsn, "%s/spec_%s_eta%.2f.h5", save, dtype, eta);

	hsize_t f_dims[2];
	GetDimsH5(fen, "/energy", f_dims);
	double band_e[f_dims[0]][f_dims[1]], band_w[f_dims[0]][f_dims[1]];
	ReadH5(fen, "/energy", (double*)band_e);
	ReadH5(fwn, "/weight", (double*)band_w);

	time_t t0 = time(NULL);

	int i, f_len = f_dims[0];
	double broad[BINS_MAX];
	double spec[f_len][BINS_MAX * Nkb];

	// options
	/*
	if(strstr(dtype, "f")) { // ignore energy over fermi level
		for(i=0; i<BINS_MAX; i++) {
			if(e[i] > 0) weight[i] = -1;
			else         weight[i] =  1;
		}
	}
	else for(i=0; i<BINS_MAX; i++) weight[i] = 1;
	*/

	if(strstr(dtype, "b")) { // add linear broadening regard to distance from fermi level
		for(i=0; i<BINS_MAX; i++) {
			broad[i] = 0.1 + fabs(e[i] / e[0]) * eta;
		}
	}
	else for(i=0; i<BINS_MAX; i++) broad[i] = eta;

#pragma omp parallel for ordered
	for(i=0; i<f_len; i++) {
		Spec(band_e[i], band_w[i], e, broad, spec[i]);
	}

	hsize_t dims[2] = {f_len, BINS_MAX * Nkb};
	WriteH5(fsn, "/spec", 2, dims, (double*)spec);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fsn, t1 - t0);
}

/*
void ReduceSpec(char *save, char *dtype, double eta, double *e, int bins) {
	char fsn[256], frn[256];

	sprintf(fsn, "%s/spec_%s_eta%.2f_bins%d.h5", save, dtype, eta, BINS_MAX);
	sprintf(frn, "%s/spec_%s_eta%.2f_bins%d.h5", save, dtype, eta, bins);

	// read spec.h5
	hsize_t fs_dims[2];
	GetDimsH5(fsn, "/spec", fs_dims);
	double spec[fs_dims[0]][fs_dims[1]];
	ReadH5(fsn, "/spec", (double*)spec);
	
	time_t t0 = time(NULL);

	int i, j, m, n, Nclus = BINS_MAX / bins;
	double rspec[fs_dims[0]][bins * Nkb];
	memset(rspec, 0, sizeof(rspec));

	for(i=0; i<fs_dims[0]; i++) {
		for(j=0; j<bins; j++) {
			for(m=0; m<Nclus; m++) {
				for(n=0; n<Nkb; n++) {
					rspec[i][Nkb*j + n] += spec[i][(Nkb*Nclus)*j + Nkb*m + n];
				}
			}
		}
	}	

	hsize_t fr_dims[2] = {fs_dims[0], bins * Nkb};
	WriteH5(frn, "/spec", 2, fr_dims, (double*)rspec);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, frn, t1 - t0);
}
*/

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("%s <save> <base/spec> <dtype> <(s)eta> : generate spec dataset\n", argv[0]);
		exit(1);
	}

	char save[256], *dtype = argv[3];
	sprintf(save, "output/%s/magspec", argv[1]);

	if(strstr(argv[2], "b")) {
		omp_set_num_threads(1); 
		GenDatabase(save, dtype);
		return 0;
	}
	else omp_set_num_threads(16);

	// read energy.h5
	char en[256];
	sprintf(en, "%s/energy.h5", save);
	hsize_t dims[1];
	GetDimsH5(en, "/energy", dims);
	double e[dims[0]];
	ReadH5(en, "/energy", e);

	GenSpec(save, dtype, atof(argv[4]), e);

	return 0;
}
