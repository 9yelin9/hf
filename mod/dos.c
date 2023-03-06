// mod/dos.c : generate dos dataset

#define Nbs      12   // num of bases
#define Nhsp     4    // num of high symmetry points
#define BINS_MAX 1024 // num of energy bins
#define PEAK_MAX 2    // num of peaks

#include "hf3.h" 

typedef struct BandFeatures {
	double e[Nbs * Nhsp];
	double w[Nbs * Nhsp];
} BFeats;

/*
void GenEnergy(char *save) {
	char fn[128];
	sprintf(fn, "%s/energy.h5", save);

	time_t t0 = time(NULL);

	int i;
	double e[BINS_MAX];

	for(i=0; i<BINS_MAX; i++) e[i] = E_MIN + (E_MAX - E_MIN) * i / BINS_MAX;

	hsize_t dims[1] = {BINS_MAX};
	WriteH5(fn, "/energy", 1, dims, e);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fn, t1 - t0);
}
*/

#define D_D_IDX (BINS_MAX*hsp + i)
#define D_B_IDX (Nbs*hsp + j)

void DOS(BFeats b, double *e, double *weight, double *broad, double *dos) {
	int i, j, hsp;
	double sum;

	for(hsp=0; hsp<Nhsp; hsp++) {
		for(i=0; i<BINS_MAX; i++) {
			sum = 0;
			for(j=0; j<Nbs; j++) {
				sum += (broad[i] / (pow(e[i] - b.e[D_B_IDX], 2) + pow(broad[i], 2))) * b.w[D_B_IDX] * weight[i];
			}
			dos[D_D_IDX] = sum / M_PI;
		}
	}
}

void GenDOSK(char *save, char *dtype, double eta, double *e) { // k-projected DOS at high symmetry points
	FILE *fb;
	char fbn[128], fdn[128];

	sprintf(fbn, "%s/band_%c.txt", save, dtype[0]);
	sprintf(fdn, "%s/dos_%s_eta%.2f.h5", save, dtype, eta);
	fb = OpenFile(fbn, "r");

	time_t t0 = time(NULL);

	int i, j, fb_len = 0;
	char buf[2*(Nbs * Nhsp)][16], tmp[1048576];
	double weight[BINS_MAX], broad[BINS_MAX];
	BFeats b;

	fgets(tmp, sizeof(tmp), fb);
	while(1) {
		fgets(tmp, sizeof(tmp), fb);
		if(feof(fb)) break;
		fb_len++;
	}
	rewind(fb);
	double dos[fb_len][BINS_MAX * Nhsp];

	// options
	if(strstr(dtype, "f")) { // ignore energy over fermi level
		for(i=0; i<BINS_MAX; i++) {
			if(e[i] > 0) weight[i] = -1;
			else         weight[i] =  1;
		}
	}
	else for(i=0; i<BINS_MAX; i++) weight[i] = 1;

	if(strstr(dtype, "b")) { // add linear broadening regard to distance from fermi level
		for(i=0; i<BINS_MAX; i++) {
			broad[i] = 0.1 + fabs(e[i] / e[0]) * eta;
		}
	}
	else for(i=0; i<BINS_MAX; i++) broad[i] = eta;

	for(i=0; i<2*(Nbs * Nhsp); i++) fscanf(fb, "%s", buf[i]); // ignore first line

	for(i=0; i<fb_len; i++) {
		for(j=0; j<Nbs * Nhsp; j++) fscanf(fb, "%lf", &b.e[j]);
		for(j=0; j<Nbs * Nhsp; j++) fscanf(fb, "%lf", &b.w[j]);

		DOS(b, e, weight, broad, dos[i]);
	}	
	fclose(fb);

	hsize_t dims[2] = {fb_len, BINS_MAX * Nhsp};
	WriteH5(fdn, "/dos", 2, dims, (double*)dos);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fdn, t1 - t0);
}

#define R_RD_K_IDX (bins*hsp + j)
#define R_D_K_IDX  (BINS_MAX*hsp + Nclus*j + k)
#define R_D_L_IDX  (Nclus*j + k)

void ReduceDOS(char *save, char *dtype, double eta, int bins) {
	char fdn[128], frn[128];

	sprintf(fdn, "%s/dos_%s_eta%.2f.h5", save, dtype, eta);
	sprintf(frn, "%s/dos_%s_eta%.2f_bins%d.h5", save, dtype, eta, bins);

	// read dos.h5
	hsize_t fd_dims[2];
	GetDimsH5(fdn, "/dos", fd_dims);
	double dos[fd_dims[0]][fd_dims[1]];
	ReadH5(fdn, "/dos", (double*)dos);

	time_t t0 = time(NULL);

	int i, j, k, hsp, Nclus = BINS_MAX / bins;

	if(strstr(dtype, "K")) {
		double rdos[fd_dims[0]][bins * Nhsp];
		memset(rdos, 0, sizeof(rdos));

		for(i=0; i<fd_dims[0]; i++) {
			for(hsp=0; hsp<Nhsp; hsp++) {
				for(j=0; j<bins; j++) {
					for(k=0; k<Nclus; k++) rdos[i][R_RD_K_IDX] += dos[i][R_D_K_IDX];
				}
			}
		}

		hsize_t fr_dims[2] = {fd_dims[0], bins * Nhsp};
		WriteH5(frn, "/dos", 2, fr_dims, (double*)rdos);
	}
	else if(strstr(dtype, "L")) {
		double rdos[fd_dims[0]][bins];
		memset(rdos, 0, sizeof(rdos));

		for(i=0; i<fd_dims[0]; i++) {
			for(j=0; j<bins; j++) {
				for(k=0; k<Nclus; k++) {
					rdos[i][j] += dos[i][R_D_L_IDX];
				}
			}
		}

		hsize_t fr_dims[2] = {fd_dims[0], bins};
		WriteH5(frn, "/dos", 2, fr_dims, (double*)rdos);
	}
	else printf("'%s' is wrong dtype\n", dtype);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, frn, t1 - t0);
}

double AvgEnergy(int Nclus, int idx, double *e) {
	int i;
	double energy = 0;
	for(i=0; i<Nclus; i++) energy += e[Nclus*idx + i];
	return energy / Nclus;
}

#define P_K_IDX(j) (bins*hsp + j)
#define GRADIENT_K_E(j) ((dos[i][P_K_IDX(j)-1] - dos[i][P_K_IDX(j)]) / (e[j-1] - e[j]))
#define GRADIENT_K_H(j) ((dos[i][P_K_IDX(j)+1] - dos[i][P_K_IDX(j)]) / (e[j+1] - e[j]))
#define GRADIENT_L_E(j) ((dos[i][j-1] - dos[i][j]) / (e[j-1] - e[j]))
#define GRADIENT_L_H(j) ((dos[i][j+1] - dos[i][j]) / (e[j+1] - e[j]))

void GenPeak(char *save, char *dtype, double eta, double *e, int bins) {
	char fdn[128], fpn[128];

	if(bins == BINS_MAX) {
		sprintf(fdn, "%s/dos_%s_eta%.2f.h5",  save, dtype, eta);
		sprintf(fpn, "%s/peak_%s_eta%.2f.h5", save, dtype, eta);
	}
	else {
		sprintf(fdn, "%s/dos_%s_eta%.2f_bins%d.h5",  save, dtype, eta, bins);
		sprintf(fpn, "%s/peak_%s_eta%.2f_bins%d.h5", save, dtype, eta, bins);
	}

	hsize_t fd_dims[2];
	GetDimsH5(fdn, "/dos", fd_dims);

	double dos[fd_dims[0]][fd_dims[1]];
	ReadH5(fdn, "/dos", (double*)dos);

	time_t t0 = time(NULL);

	int i, j, grd, grd_old, Nclus = BINS_MAX / bins;
	double tol_g = -1e-3, tol_w=0.1; // tol_w = 5.5 일때 제일 좋았음
	
	if(strstr(dtype, "D")) tol_w = 1.0;

	if(strstr(dtype, "K")) {
		int hsp;
		double peak[fd_dims[0]][PEAK_MAX * Nhsp];

		for(i=0; i<fd_dims[0]; i++) {
			for(j=0; j<PEAK_MAX * Nhsp; j++) peak[i][j] = 99;
		}

		for(i=0; i<fd_dims[0]; i++) {
			for(hsp=0; hsp<Nhsp; hsp++) {
				// electron
				grd_old = GRADIENT_K_E(bins/2) > -tol_g;
				for(j=bins/2; j>0; j--) {
					grd = GRADIENT_K_E(j) > -tol_g;
					if(grd - grd_old > 0  && dos[i][P_K_IDX(j)] > tol_w) {
						peak[i][PEAK_MAX*hsp + 0] = AvgEnergy(Nclus, j, e);
						break;
					}
					grd_old = grd;
				}

				// hole
				grd_old = GRADIENT_K_H(bins/2) < tol_g;
				for(j=bins/2; j<bins; j++) {
					grd = GRADIENT_K_H(j) < tol_g;
					if(grd - grd_old > 0 && dos[i][P_K_IDX(j)] > tol_w) {
						peak[i][PEAK_MAX*hsp + 1] = AvgEnergy(Nclus, j, e);
						break;
					}
					grd_old = grd;
				}
			}
		}

		hsize_t fp_dims[2] = {fd_dims[0], PEAK_MAX * Nhsp};
		WriteH5(fpn, "/peak", 2, fp_dims, (double*)peak);
	}
	else if(strstr(dtype, "L")) {
		double peak[fd_dims[0]][PEAK_MAX];

		for(i=0; i<fd_dims[0]; i++) {
			for(j=0; j<PEAK_MAX; j++) peak[i][j] = 99;
		}

		for(i=0; i<fd_dims[0]; i++) {
			// electron
			grd_old = GRADIENT_L_E(0) > -tol_g;
			for(j=bins/2; j>0; j--) {
				grd = GRADIENT_L_E(j) > -tol_g;
				if(grd - grd_old > 0  && dos[i][j] > tol_w) {
					peak[i][PEAK_MAX + 0] = AvgEnergy(Nclus, j, e);	
					break;
				}
				grd_old = grd;
			}

			// hole
			grd_old = GRADIENT_L_H(bins/2) < tol_g;
			for(j=bins/2; j<bins; j++) {
				grd = GRADIENT_L_H(j) < tol_g;
				if(grd - grd_old > 0 && dos[i][j] > tol_w) {
					peak[i][PEAK_MAX + 1] = AvgEnergy(Nclus, j, e);	
					break;
				}
				grd_old = grd;
			}
		}

		hsize_t fp_dims[2] = {fd_dims[0], PEAK_MAX};
		WriteH5(fpn, "/peak", 2, fp_dims, (double*)peak);
	}
	else printf("'%s' is wrong dtype\n", dtype);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fpn, t1 - t0);
}

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("%s <save> <dos/rdos/peak/spec> <dtype> <eta> <(r)bins> : generate dos dataset\n", argv[0]);
		exit(1);
	}
	omp_set_num_threads(1); 

	char *dtype = argv[3];
	double eta  = atof(argv[4]);

	char save[64], fen[128];
	sprintf(save, "output/%s/magstr", argv[1]);
	sprintf(fen, "%s/energy.h5", save);

	hsize_t dims[1];
	GetDimsH5(fen, "/energy", dims);

	double e[dims[0]];
	ReadH5(fen, "/energy", e);

	if     (strstr(argv[2], "r")) ReduceDOS(save, dtype, eta, atoi(argv[5]));
	else if(strstr(argv[2], "p")) GenPeak(save, dtype, eta, e, atoi(argv[5]));
	else                          GenDOSK(save, dtype, eta, e);

	return 0;
}
