// lib/tb.c : functions for tight-binding

void FourierN(Config c, int Nl, Lattice *lat, double *k, lapack_complex_double *tb) {
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

void FourierQ(Config c, int Nl, Lattice *lat, double *k, lapack_complex_double *tb) {
	int i, j, v = 0;
	double dot, dotQ;
	lapack_complex_double tb0[c.Nb][c.Nb], t, e, eQ;

	memset(tb0, 0, sizeof(tb0));

	for(i=0; i<Nl; i++) {
		dot = dotQ = 0;
		for(j=0; j<DIM; j++) {
			dot  += lat.c[i][j] * (k[j]);
			dotQ += lat.c[i][j] * (k[j] + c.Q[j]);
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

void CalcTB(Config c, char *ktype, void (*Fourier)()) {
	int Nk = strstr(ktype, "g") ? Nkg : Nkb;
	char fln[256], fkn[256], ftn[256];
	sprintf(fln, "input/%s.txt",           c.lat);
	sprintf(fkn, "input/k%s_Nk%d.txt",     ktype, Nk);
	sprintf(ftn, "input/tb%s_Nk%d_%s.bin", ktype, Nk, c.type);

	time_t t0 = time(NULL);

	FILE *fl = fopen(fln, "r"), *fk = fopen(fkn, "r"), *ft = fopen(ftn, "wb");
	int i, j, Nl;
	char buf[1024];
	double k[Nk][DIM];
	lapack_complex_double tb[Nk][c.Nb*c.Nb];

	// lat
	fgets(buf, sizeof(buf), fl);
	Nl = atoi(buf);
	Lattice lat = {
		.site = malloc(Nl * sizeof(lat.site)),
		.obt  = malloc(Nl * sizeof(lat.obt)),
		.t    = malloc(Nl * sizeof(lat.t))
	};

	for(i=0; i<Nl; i++) fscanf(f, "%lf%lf%lf%d%d%lf%lf", &lat.site[i][0], &lat.site[i][1], &lat.site[i][2], &lat.obt[i][0], &lat.obt[i][1], &lat.t[i][0], &lat.t[i][1]);
	fclose(fl);
	
	// k
	fgets(buf, sizeof(buf), fk); // skip header
	for(i=0; i<Nk; i++) {
		fgets(buf, sizeof(buf), fk);
		for(j=0; j<DIM; j++) sscanf(buf, "%lf", &k[i][j]);
	}
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

	/*
	// make band
	FILE *f;
	char fn[64];
	
	sprintf(fn, "input/%s/band_%c.txt", c.name, c.type[0]);
	f  = OpenFile(fn,  "w");

	t0 = time(NULL);

	LAPACK lp = {
		.jobz = 'V',
		.uplo = 'L',
		.rwork = (double*)malloc(sizeof(double) * (3*c.Nb-2)),
		.ln = c.Nb,
		.lda = c.Nb,
		.lwork = 2*c.Nb-1,
		.work = (lapack_complex_double*)malloc(sizeof(lapack_complex_double) * lp.lwork)
	};

	double ev0[c.Nb];
	lapack_complex_double es0[c.Nbb];

	for(i=0; i<Nkb; i++) {
		for(j=0; j<c.Nbb; j++) es0[j] = tbb[c.Nbb*i + j];

		LAPACK_zheev(&lp.jobz, &lp.uplo, &lp.ln, es0, &lp.lda, ev0, lp.work, &lp.lwork, lp.rwork, &lp.info);
		if(lp.info != 0) {
			printf("LAPACK_zheev FAIL\n");
			exit(1);
		}

		for(j=0; j<c.Nb; j++) fprintf(f, "%12f", ev0[j]);
		fprintf(f, "\n");
	}
	fclose(f);

	t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fn, t1 - t0);
	*/
}

