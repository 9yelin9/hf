// lib/tb.c : functions for tight-binding Hamiltonian

#include "hf.h"

void FourierN(Config c, int Nl, Lattice lat, double *k, lapack_complex_double *tb) {
	int i, j, v = 0;
	double dot;
	lapack_complex_double tb0[c.Nb][c.Nb], t, e;

	memset(tb0, 0, sizeof(tb0));

	for(i=0; i<Nl; i++) {
		dot = 0;
		for(j=0; j<DIM; j++) {
			dot += lat.site[i][j] * k[j];
		}

		t = lat.t[i][0] + lat.t[i][1] * I;
		e = cos(dot)    + sin(dot)    * I;

		tb0[lat.obt[i][0]-1][lat.obt[i][1]-1] += t * e;
	}

	for(i=0; i<c.Ns; i++) {
		for(j=0; j<c.Ns; j++) {
			tb0[i+c.Ns][j+c.Ns] = tb0[i][j];
		}
	}

	for(i=0; i<c.Nb; i++) {
		for(j=0; j<c.Nb; j++) {
			tb[v] = tb0[i][j];
			v++;
		}
	}
}	

void FourierQ(Config c, int Nl, Lattice lat, double *k, lapack_complex_double *tb) {
	int i, j, v = 0;
	double dot, dotQ;
	lapack_complex_double tb0[c.Nb][c.Nb], t, e, eQ;

	memset(tb0, 0, sizeof(tb0));

	for(i=0; i<Nl; i++) {
		dot = dotQ = 0;
		for(j=0; j<DIM; j++) {
			dot  += lat.site[i][j] * (k[j]);
			dotQ += lat.site[i][j] * (k[j] + c.Q[j]);
		}

		t  = lat.t[i][0] + lat.t[i][1] * I;
		e  = cos(dot)    + sin(dot)    * I;
		eQ = cos(dotQ)   + sin(dotQ)   * I;

		tb0[lat.obt[i][0]-1       ][lat.obt[i][1]-1       ] += t * e;
		tb0[lat.obt[i][0]-1 + c.Nc][lat.obt[i][1]-1 + c.Nc] += t * eQ;
	}

	for(i=0; i<c.Ns; i++) {
		for(j=0; j<c.Ns; j++) {
			tb0[i+c.Ns][j+c.Ns] = tb0[i][j];
		}
	}

	for(i=0; i<c.Nb; i++) {
		for(j=0; j<c.Nb; j++) {
			tb[v] = tb0[i][j];
			v++;
		}
	}
}

void GenTB(Config c, char *ktype, void (*Fourier)()) {
	int Nk = strstr(ktype, "g") ? Nkg : Nkb;
	char fln[256], fkn[256], ftn[256];
	sprintf(fln, "input/%s.txt",           c.lat);
	sprintf(fkn, "input/k%s_Nk%d.txt",     ktype, Nk);
	sprintf(ftn, "input/tb%s_Nk%d_%s.bin", ktype, Nk, c.type);

	time_t t0 = time(NULL);

	FILE *fl = fopen(fln, "r"), *fk = fopen(fkn, "r"), *ft = fopen(ftn, "wb");
	int i, Nl;
	char buf[1024];
	double k[Nk][DIM];
	lapack_complex_double tb[Nk][c.Nb*c.Nb];

	// lat
	fgets(buf, sizeof(buf), fl); // skip header
	Nl = atoi(buf);
	Lattice lat = {
		.site = calloc(Nl, sizeof(*lat.site)),
		.obt  = calloc(Nl, sizeof(*lat.obt)),
		.t    = calloc(Nl, sizeof(*lat.t))
	};

	for(i=0; i<Nl; i++) fscanf(fl, "%d%d%d%d%d%lf%lf", &lat.site[i][0], &lat.site[i][1], &lat.site[i][2], &lat.obt[i][0], &lat.obt[i][1], &lat.t[i][0], &lat.t[i][1]);
	fclose(fl);
	
	// k
	fgets(buf, sizeof(buf), fk); // skip header
	for(i=0; i<Nk; i++) fscanf(fk, "%lf%lf%lf", &k[i][0], &k[i][1], &k[i][2]);
	fclose(fk);

	// tb
#pragma omp parallel for ordered
	for(i=0; i<Nk; i++) Fourier(c, Nl, lat, k[i], tb[i]);
	fwrite(tb, sizeof(tb), 1, ft);
	fclose(ft);

	free(lat.site);
	free(lat.obt);
	free(lat.t);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, ftn, t1 - t0);
}

void GenTBBand(Config c) {
	char ftn[256], fbn[256];
	sprintf(ftn, "input/tbb_Nk%d_%s.bin",  Nkb, c.type);
	sprintf(fbn, "input/band_Nk%d_%s.txt", Nkb, c.type);

	time_t t0 = time(NULL);

	FILE *ft = fopen(ftn, "rb"), *fb = fopen(fbn, "w");
	int i, j;
	double ev[Nkb][c.Nb];
	lapack_complex_double tb[Nkb][c.Nb*c.Nb];

	// tb
	fread(tb, sizeof(tb), 1, ft);
	fclose(ft);

	// band
	for(i=0; i<c.Nb; i++) fprintf(fb, "%20s%02d", "e", i);
	fprintf(fb, "\n");

	for(i=0; i<Nkb; i++) {
		CalcEigen(c, ev[i], tb[i]);

		for(j=0; j<c.Nb; j++) fprintf(fb, "%22.16f", ev[i][j]);
		fprintf(fb, "\n");
	}
	fclose(fb);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fbn, t1 - t0);
}
