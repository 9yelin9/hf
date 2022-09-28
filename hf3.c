// hf3.c : universal functions for calculating hf3 model

#include "hf3.h" 

FILE* OpenFile(char *fs, char *ftype) {
	FILE *f;

	if((f = fopen(fs, ftype)) == NULL) {
		printf("%s fopen FAIL", fs);
		exit(1);
	}

	return f;
}

int GetLen(char *fs) {
	FILE *f = OpenFile(fs, "r");

	int len = 0;

	fscanf(f, "%d", &len);
	fclose(f);

	return len;
}

void ShowBlock(int basis, lapack_complex_double *block) {
	int i;

	for(i=0; i<basis*basis; i++) {
		if(i % basis == 0) printf("\n");
		printf("%.3f\t", creal(block[i]));
	}
	printf("\n");
}

void GenName(Solution *s, char *data_type, char *fs) {
	sprintf(fs,\
			"output/JU%.2f_SOC%.2f/%s_%s_N%.1f_U%.1f_n%f_m%f_fermi%f_e%f_%s.txt",\
			s->JU, s->SOC, data_type, s->type, s->N, s->U, s->ntot, s->mtot, s->fermi, s->e, s->runtime);
}

void DotProd(int g_len, Vector *g, lapack_complex_double *coef) {
	int i, j, k;
	double dot;
	Vector r[SUPER] = {{{0, 0, 0}}, {{1, 0, 0}}};

	for(i=0; i<g_len; i++) {
		for(j=0; j<SUPER; j++) {
			dot = 0;
			for(k=0; k<3; k++) dot += r[j].c[k] * g[i].c[k];
			coef[2*i + j] = cos(dot) - sin(dot) * I;
		}
	}
}

void TransBasis(lapack_complex_double *coef, lapack_complex_double *es) {
	int i, j, k;
	lapack_complex_double es1, es2;

	for(i=0; i<BAND; i++) {
		for(j=0; j<DOUBLE; j++) {
			for(k=0; k<SINGLE; k++) {
				es1 = sqrt(0.5) * coef[2*i + 1] * (es[STATE_IDX + k] + es[STATE_IDX + k+OBT]);
				es2 = sqrt(0.5) * coef[2*i + 0] * (es[STATE_IDX + k] - es[STATE_IDX + k+OBT]);

				es[STATE_IDX + k] = es1;
				es[STATE_IDX + k+OBT] = es2;
			}
		}
	}
}

void CalcGauss() {
	FILE *fg = OpenFile("input/gg.bin", "wb"), *fw = OpenFile("input/wg.bin", "wb");

	int i, i0, i1, i2;
	double g[GAUSS], w[GAUSS], wg[GAUSS3];
	Vector gg[GAUSS3];
	gsl_integration_glfixed_table *t = gsl_integration_glfixed_table_alloc(GAUSS);

	for(i=0; i<GAUSS; i++) gsl_integration_glfixed_point(-M_PI, M_PI, i, &g[i], &w[i], t);

	for(i=0; i<GAUSS3; i++) {
		i0 = i / (GAUSS * GAUSS);
		i1 =(i / GAUSS) % GAUSS;
		i2 = i % GAUSS;

		gg[i].c[0] = g[i0];
		gg[i].c[1] = g[i1];
		gg[i].c[2] = g[i2];
		wg[i] = w[i0] * w[i1] * w[i2];
	}

	fwrite(gg, sizeof(gg), 1, fg);
	fwrite(wg, sizeof(wg), 1, fw);
	fclose(fg);
	fclose(fw);
}

void CalcCoef() {
	FILE *fg = OpenFile("input/coefg.bin", "wb"), *fb = OpenFile("input/coefb.bin", "wb");
	Vector gg[GAUSS3], gb[BAND];
	lapack_complex_double coefg[2 * GAUSS3], coefb[2 * BAND];

	ReadBin("input/gg.bin", sizeof(Vector) * GAUSS3, gg);
	ReadBin("input/gb.bin", sizeof(Vector) * BAND, gb);

	DotProd(GAUSS3, gg, coefg);
	DotProd(BAND, gb, coefb);

	fwrite(coefg, sizeof(coefg), 1, fg);
	fwrite(coefb, sizeof(coefb), 1, fb);
	fclose(fg);
	fclose(fb);
}

void CalcTB(char *type) {
	FILE *fg, *fb, *f;
	char fgs[64], fbs[64], fs[64];

	sprintf(fgs, "input/tbg_%s.bin", type);
	sprintf(fbs, "input/tbb_%s.bin", type);
	sprintf(fs, "input/band_%s.txt", type);

	fg = OpenFile(fgs, "wb");
	fb = OpenFile(fbs, "wb");
	f  = OpenFile(fs,  "w");

	time_t t0 = time(NULL);

	int i, j;
	char ctype[8];
	Vector q, gg[GAUSS3], gb[BAND];
	void (*Fourier)(int, int, Lattice*, Vector, Vector, lapack_complex_double*);

	ReadInfo(type, ctype, &q);
	ReadBin("input/gg.bin", sizeof(Vector) * GAUSS3, gg);
	ReadBin("input/gb.bin", sizeof(Vector) * BAND, gb);

	const int l_len = GetLen("input/lattice.txt");
	Lattice l[l_len];
	ReadLattice(l_len, l);

	const int basis = strstr(ctype, "s") ? SINGLE : DOUBLE;
	lapack_complex_double tbg[basis*basis * GAUSS3], tbb[basis*basis * BAND];
	Fourier = strstr(ctype, "s") ? FourierS : FourierD;

#pragma omp parallel for
	for(i=0; i<GAUSS3; i++) Fourier(l_len, i, l, gg[i], q, tbg);
#pragma omp parallel for
	for(i=0; i<BAND; i++) Fourier(l_len, i, l, gb[i], q, tbb);

	fwrite(tbg, sizeof(tbg), 1, fg);
	fwrite(tbb, sizeof(tbb), 1, fb);
	fclose(fg);
	fclose(fb);

	time_t t1 = time(NULL);
	printf("%s(%s, %s) : %lds\n", __func__, fgs, fbs, t1 - t0);

	// make band
	LParameter lp = {
		.jobz = 'V',
		.uplo = 'L',
		.rwork = (double*)malloc((3*basis-2) * sizeof(double)),
		.ln = basis,
		.lda = basis,
		.lwork = 2*basis-1,
		.work = (lapack_complex_double*)malloc(lp.lwork * sizeof(lapack_complex_double))
	};
	double ev[basis * BAND];
	lapack_complex_double es[basis*basis * BAND];

	for(i=0; i<BAND; i++) {
		for(j=0; j<basis*basis; j++) es[j] = tbb[basis*basis*i + j];

		LAPACK_zheev(&lp.jobz, &lp.uplo, &lp.ln, es, &lp.lda, ev, lp.work, &lp.lwork, lp.rwork, &lp.info);
		if(lp.info != 0) {
			printf("LAPACK_zheev FAIL\n");
			exit(1);
		}

		fprintf(f, "%4d", i);
		for(j=0; j<basis; j++) fprintf(f, "%12f", ev[j]);
		fprintf(f, "\n");
	}
	fclose(f);
}

void CalcEigen(Solution *s, LParameter *lp, int tb_len, lapack_complex_double *tb, double *ev, lapack_complex_double *es, Energy *e) {
	e->min =  100;
	e->max = -100;

	int i, j;
	double ev_block[s->basis];
	lapack_complex_double tb_block[s->basis*s->basis];
	void (*Interaction)(Solution*, lapack_complex_double*);

	Interaction = strstr(s->ctype, "s") ? InteractionS : InteractionD;

	for(i=0; i<tb_len; i++) {
		for(j=0; j<s->basis*s->basis; j++) tb_block[j] = tb[s->basis*s->basis*i + j];

		Interaction(s, tb_block);

		LAPACK_zheev(&lp->jobz, &lp->uplo, &lp->ln, tb_block, &lp->lda, ev_block, lp->work, &lp->lwork, lp->rwork, &lp->info);
		if(lp->info != 0) {
			printf("LAPACK_zheev FAIL\n");
			exit(1);
		}

		for(j=0; j<s->basis; j++) ev[s->basis*i + j] = ev_block[j];
		for(j=0; j<s->basis*s->basis; j++) es[s->basis*s->basis*i + j] = tb_block[j];

		if(ev_block[0] < e->min) e->min = ev_block[0];
		if(ev_block[s->basis-1] > e->max) e->max = ev_block[s->basis-1];
	}
}

void CalcSolution(Solution *s, LParameter *lp, lapack_complex_double *tbg) {
	FILE *f;
	char fs[256];

	GenName(s, "sol", fs);
	f = OpenFile(fs, "w");

	time_t t0 = time(NULL);

	int itr, i, j, is_cvg;
	double wg[GAUSS3], ev[s->basis * GAUSS3], fermi0, fermi1, e, n[OBT], m[OBT], ntot = 0, mtot = 0, cvg[OBT][3];
	Energy energy;
	lapack_complex_double es[s->basis*s->basis * GAUSS3];
	void (*GaussQuad)(double, double*, double*, lapack_complex_double*, double*, double*, double*);

	GaussQuad = strstr(s->ctype, "s") ? GaussQuadS : GaussQuadD;
	ReadBin("input/wg.bin", sizeof(double) * GAUSS3, wg);

	for(i=0; i<OBT; i++) {
		for(j=0; j<3; j++) cvg[i][j] = 100;
	}

	fprintf(f, "%5s%12s%12s", "#itr", "fermi", "e");
	for(i=0; i<OBT; i++) fprintf(f, "%10s%02d%10s%02d", "n", i+1, "m", i+1);
	fprintf(f, "%12s%12s\n", "ntotal", "mtotal");

	fprintf(f, "%5d%12f%12f", 0, s->fermi, s->e);
	for(i=0; i<OBT; i++) fprintf(f, "%12f%12f", s->n[i], s->m[i]);
	fprintf(f, "%12f%12f\n", s->ntot, s->mtot);

	for(itr=0; itr<50; itr++) {
		CalcEigen(s, lp, GAUSS3, tbg, ev, es, &energy);
		fermi0 = energy.min;
		fermi1 = energy.max;
		s->fermi = 0.5 * (fermi0 + fermi1);

		while(fabs(s->fermi - fermi1) > 1e-8) {
			ntot = 0;
			mtot = 0;

			GaussQuad(s->fermi, wg, ev, es, n, m, &e);
		
			for(i=0; i<OBT; i++) {	
				n[i] /= pow(2*M_PI , 3);
				m[i] /= pow(2*M_PI , 3);
				ntot += n[i];
				mtot += m[i];
			}
			e /= pow(2*M_PI, 3);

			if(fabs(ntot - s->N) < 1e-6) break;
			else {
				if(ntot < s->N) fermi0 = s->fermi;
				else fermi1 = s->fermi;
				s->fermi = 0.5 * (fermi0 + fermi1);
			}
		}

		for(i=0; i<OBT; i++) {
			s->n[i] = n[i];
			s->m[i] = m[i];
		}
		s->ntot = ntot;
		s->mtot = mtot;

		fprintf(f, "%5d%12f%12f", itr+1, s->fermi, s->e);
		for(i=0; i<OBT; i++) fprintf(f, "%12f%12f", s->n[i], s->m[i]);
		fprintf(f, "%12f%12f\n", s->ntot, s->mtot);
			
		is_cvg = 0;
		for(i=0; i<OBT; i++) {
			cvg[i][itr % 3] = m[i];
			if(fabs((cvg[i][0] + cvg[i][1] + cvg[i][2])/3 - cvg[i][itr % 3]) < 1e-3) is_cvg++;
		}
		if(is_cvg == 3) break;
	}

	fclose(f);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fs, t1 - t0);
}

void ReadPathInfo(int sym_len, int *sym) {
	FILE *f = OpenFile("input/info.txt", "r");

	int i, path, tmp = 0;
	char buf[1024];

	fgets(buf, sizeof(buf), f);

	for(i=0; i<sym_len; i++) {
		fscanf(f, "%d", &path);	
		sym[i] = path - tmp;
		tmp = path;
	}
	fclose(f);
}

void ReadLattice(int l_len, Lattice *l) {
	FILE *f = OpenFile("input/lattice.txt", "r");

	int i;
	char buf[1024];

	fgets(buf, sizeof(buf), f);

	for(i=0; i<l_len; i++) {
		fscanf(f,\
				"%d%d%d%d%d%lf%lf",\
				&l[i].c[0], &l[i].c[1], &l[i].c[2], &l[i].obt1, &l[i].obt2, &l[i].tre, &l[i].tim);
	}
	fclose(f);
}

void ReadInfo(char *type, char *ctype, Vector *q) {
	FILE *f = OpenFile("input/info.txt", "r");

	int i;
	char buf[1024], tmp[8];

	while(fgets(buf, sizeof(buf), f)) {
		if(strstr(buf, type)) {
			sscanf(buf, "%s%s%lf%lf%lf", tmp, ctype, &q->c[0], &q->c[1], &q->c[2]);
			break;
		}
	}
	fclose(f);

	for(i=0; i<3; i++) q->c[i] *= M_PI;
}

void ReadBin(char *fs, int bin_size, void *bin) {
	FILE *f = OpenFile(fs, "rb");

	fread(bin, bin_size, 1, f);
	fclose(f);
}

void MakeBand(Solution *s, LParameter *lp, lapack_complex_double *tbb) {
	FILE *f;
	char fs[256];

	GenName(s, "band", fs);
	f = OpenFile(fs, "w");

	time_t t0 = time(NULL);

	int i, j;
	double ev[s->basis * BAND];
	Energy energy;
	lapack_complex_double es[s->basis*s->basis * BAND];

	CalcEigen(s, lp, BAND, tbb, ev, es, &energy);

	fprintf(f, "%4s", "#");
	for(i=0; i<s->basis; i++) fprintf(f, "%10s%02d", "e", i+1);
	fprintf(f, "\n");

	for(i=0; i<BAND; i++) {
		fprintf(f, "%4d", i);
		for(j=0; j<s->basis; j++) fprintf(f, "%12f", ev[s->basis*i + j]);
		fprintf(f, "\n");
	}

	fclose(f);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fs, t1 - t0);

	if(!strstr(s->ctype, "s")) MakeUFW(s, ev, es);
}

void MakeUFW(Solution *s, double *ev, lapack_complex_double *es) {
	FILE *f;
	char fs[256];

	GenName(s, "ufw", fs);
	f = OpenFile(fs, "w");

	time_t t0 = time(NULL);

	int i, j, k;
	double p2;
	lapack_complex_double coefb[2 * BAND], p; 

	ReadBin("input/coefb.bin", sizeof(lapack_complex_double) * 2 * BAND, coefb);
	TransBasis(coefb, es);

	fprintf(f, "%4s", "#");
	for(i=0; i<s->basis; i++) fprintf(f, "%10s%02d", "w", i+1);
	fprintf(f, "\n");

	for(i=0; i<BAND; i++) {
		fprintf(f, "%4d", i);

		for(j=0; j<DOUBLE; j++) {
			p2 = 0;

			for(k=0; k<SINGLE; k++) {
				p = es[STATE_IDX + k] * coefb[2*i + 0] + es[STATE_IDX + k+OBT] * coefb[2*i + 1];
				p2 += CSQR(p) / SUPER;
			}
			fprintf(f, "%12f", p2); 
		}
		fprintf(f, "\n");
	}

	fclose(f);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fs, t1 - t0);
}

void MakeDOS(Solution *s, LParameter *lp, lapack_complex_double *tbg) {
	FILE *f;
	char fs[256];

	GenName(s, "dos", fs);
	f = OpenFile(fs, "w");

	time_t t0 = time(NULL);

	int itv, i, j;
	double ev[s->basis * GAUSS3], e, dos[s->basis];
	Energy energy;
	lapack_complex_double es[s->basis*s->basis * GAUSS3];
	
	CalcEigen(s, lp, GAUSS3, tbg, ev, es, &energy);

	if(!strstr(s->ctype, "s")) {
		lapack_complex_double coefg[2* GAUSS3];
		ReadBin("input/coefg.bin", sizeof(lapack_complex_double) * 2 * GAUSS3, coefg);
		TransBasis(coefg, es);
	}

	energy.min -= 1.0;
	energy.max += 1.0;
	
	fprintf(f, "%12s", "#e");
	for(i=0; i<s->basis; i++) fprintf(f, "%10s%02d", "d", i+1);
	fprintf(f, "\n");

	for(itv=0; itv<256; itv++) {
		memset(dos, 0, sizeof(dos));
		e = energy.min + (energy.max - energy.min) * itv / 256;

		for(i=0; i<GAUSS3*s->basis; i++) {
			for(j=0; j<s->basis; j++) {
				dos[j] += CSQR(es[s->basis*i + j]) * GREEN(i);
			}
		}

		fprintf(f, "%12f", e);
		for(i=0; i<s->basis; i++) {
			fprintf(f, "%12f", dos[i] / pow(2*M_PI, 3));
		}
		fprintf(f, "\n");
	}

	fclose(f);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fs, t1 - t0);
}

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("%s in : make reciprocal lattice vectors and basis transform coefficients\n%s tb <type> : make tight-binding matrix\n%s <type> <J/U> <SOC> <N> <U> : make Hartree-Fock approximated 3-band Hubbard model\n", argv[0], argv[0], argv[0]);
		exit(1);
	}

	omp_set_num_threads(OMP_THREAD);

	if(strstr(argv[1], "in")) {
		CalcGauss();
		CalcBandPath();
		CalcCoef();
		return 0;
	}
	else if(strstr(argv[1], "tb")) {
		CalcTB(argv[2]);
		return 0;
	}
	else {
		time_t t = time(NULL);
		struct tm *tm = localtime(&t);

		Solution s = {
			.type = argv[1],
			.JU = atof(argv[2]),
			.SOC = atof(argv[3]),
			.N = atof(argv[4]),
			.U = atof(argv[5]),
			.J = s.JU * s.U,

			.n = {s.N/3, s.N/3, s.N/3},
			.m = {0.1, 0.1, 0.1},
			.ntot = 100,
			.mtot = 100,
			.fermi = 100,
			.e = 100
		};

		sprintf(s.runtime, "%d%d%d%d", tm->tm_mon+1, tm->tm_mday, tm->tm_hour, tm->tm_min); 

		ReadInfo(s.type, s.ctype, &s.q);
		s.basis = strstr(s.ctype, "s") ? SINGLE : DOUBLE;

		LParameter lp = {
			.jobz = 'V',
			.uplo = 'L',
			.rwork = (double*)malloc((3*s.basis-2) * sizeof(double)),
			.ln = s.basis,
			.lda = s.basis,
			.lwork = 2*s.basis-1,
			.work = (lapack_complex_double*)malloc(lp.lwork * sizeof(lapack_complex_double))
		};
		lapack_complex_double tbg[s.basis*s.basis * GAUSS3], tbb[s.basis*s.basis * BAND];
		char tbgs[64], tbbs[64];

		sprintf(tbgs, "input/tbg_%s.bin", s.type);
		sprintf(tbbs, "input/tbb_%s.bin", s.type);

		ReadBin(tbgs, sizeof(lapack_complex_double) * s.basis*s.basis * GAUSS3, tbg);
		ReadBin(tbbs, sizeof(lapack_complex_double) * s.basis*s.basis * BAND, tbb);

		CalcSolution(&s, &lp, tbg);

		MakeBand(&s, &lp, tbb);
		MakeDOS(&s, &lp, tbg);

		free(lp.rwork);
		free(lp.work);

		return 0;
	}

	printf("nothing happened\n");
	exit(1);
}
