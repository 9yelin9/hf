// lib/func.c : functions for calculating Hartree-Fock approximated 3-band model

#define OBT_IDX ((i / s->basis) % OBT)
#define STATE_IDX (DOUBLE*DOUBLE*i + DOUBLE*j + OBT*(k/OBT))
#define INTER_N (0.5 * ((s->U) * n[OBT_IDX] + (s->U - 2*s->J) * n_[OBT_IDX] + (s->U - 3*s->J) * n_[OBT_IDX]))
#define INTER_M (0.5 * ((s->U) * m[OBT_IDX] + (s->U - 2*s->J) * m_[OBT_IDX] - (s->U - 3*s->J) * m_[OBT_IDX])) 

#define CSQR(c) (pow(creal(c), 2) + pow(cimag(c), 2))
#define GREEN(i) (0.05 / (pow(e - ev[i], 2) + pow(0.05, 2))) 

#include "hf3.h" 

FILE* OpenFile(char *fs, char *ftype) {
	FILE *f;

	if((f = fopen(fs, ftype)) == NULL) {
		printf("%s fopen FAIL\n", fs);
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

int GetB(char *name) {
	FILE *f;
	char fs[64];

	sprintf(fs, "input/%s/info.txt", name);
	f = OpenFile(fs, "r");

	int i, p, p_len = GetLen(fs);
	char buf[1024];

	fgets(buf, sizeof(buf), f);

	for(i=0; i<p_len; i++) fscanf(f, "%d", &p);
	fclose(f);

	return p;
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
	if(!strstr(data_type, "sol")) {
		sprintf(fs,\
				"output/%s/JU%.2f_SOC%.2f/%s_%s_N%.1f_U%.1f_n%f_m%f_e%f_gap%f_fermi%f_dntop%f_%s.txt",\
				s->name, s->JU, s->SOC, data_type, s->type, s->N, s->U, s->ntot, s->mtot, s->e, s->gap, s->fermi, s->dntop, s->runtime);
	}
	else {
		sprintf(fs,\
				"output/%s/JU%.2f_SOC%.2f/%s_%s_N%.1f_U%.1f_%s.txt",\
				s->name, s->JU, s->SOC, data_type, s->type, s->N, s->U, s->runtime);
	}
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

void TransBasis(int coef_len, lapack_complex_double *coef, lapack_complex_double *es) {
	int i, j, k;
	lapack_complex_double es1, es2;

	for(i=0; i<coef_len; i++) {
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

void CalcGauss(char *name) {
	FILE *fg, *fw;
	char fgs[64], fws[64];

	sprintf(fgs, "input/%s/gg.bin", name);
	sprintf(fws, "input/%s/wg.bin", name);

	fg = OpenFile(fgs, "wb");
	fw = OpenFile(fws, "wb");

	int i, i0, i1, i2;
	double g[G], w[G], wg[G3];
	Vector gg[G3];
	gsl_integration_glfixed_table *t = gsl_integration_glfixed_table_alloc(G);

	for(i=0; i<G; i++) gsl_integration_glfixed_point(-M_PI, M_PI, i, &g[i], &w[i], t);

	for(i=0; i<G3; i++) {
		i0 = i / (G * G);
		i1 =(i / G) % G;
		i2 = i % G;

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

void CalcCoef(int B, char *name) {
	FILE *fg, *fb;
	char fgs[64], fbs[64];

	sprintf(fgs, "input/%s/coefg.bin", name);
	sprintf(fbs, "input/%s/coefb.bin", name);

	fg = OpenFile(fgs, "wb");
	fb = OpenFile(fbs, "wb");

	char ggs[64], gbs[64];
	Vector gg[G3], gb[B];
	lapack_complex_double coefg[2 * G3], coefb[2 * B];

	sprintf(ggs, "input/%s/gg.bin", name);
	sprintf(gbs, "input/%s/gb.bin", name);

	ReadBin(ggs, sizeof(Vector) * G3, gg);
	ReadBin(gbs, sizeof(Vector) * B, gb);

	DotProd(G3, gg, coefg);
	DotProd(B, gb, coefb);

	fwrite(coefg, sizeof(coefg), 1, fg);
	fwrite(coefb, sizeof(coefb), 1, fb);
	fclose(fg);
	fclose(fb);
}

void CalcTB(int B, char *name, char *type) {
	FILE *fg, *fb, *f;
	char fgs[64], fbs[64], fs[64];

	sprintf(fgs, "input/%s/tbg_%s.bin", name, type);
	sprintf(fbs, "input/%s/tbb_%s.bin", name, type);
	sprintf(fs, "input/%s/band_%s.txt", name, type);

	fg = OpenFile(fgs, "wb");
	fb = OpenFile(fbs, "wb");
	f  = OpenFile(fs,  "w");

	time_t t0 = time(NULL);

	int i, j;
	char ctype[8], ggs[64], gbs[64];
	Vector q, gg[G3], gb[B];
	void (*Fourier)(int, int, Lattice*, Vector, Vector, lapack_complex_double*);

	sprintf(ggs, "input/%s/gg.bin", name);
	sprintf(gbs, "input/%s/gb.bin", name);

	ReadInfo(name, type, ctype, &q);
	ReadBin(ggs, sizeof(Vector) * G3, gg);
	ReadBin(gbs, sizeof(Vector) * B, gb);

	const int l_len = GetLen("input/lattice.txt");
	Lattice l[l_len];
	ReadLattice(name, l_len, l);

	const int basis = strstr(ctype, "s") ? SINGLE : DOUBLE;
	lapack_complex_double tbg[basis*basis * G3], tbb[basis*basis * B];
	Fourier = strstr(ctype, "s") ? FourierS : FourierD;

	omp_set_num_threads(OMP_THREAD);

#pragma omp parallel for ordered
	for(i=0; i<G3; i++) Fourier(l_len, i, l, gg[i], q, tbg);
#pragma omp parallel for ordered
	for(i=0; i<B; i++) Fourier(l_len, i, l, gb[i], q, tbb);

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
	double ev[basis * B];
	lapack_complex_double es[basis*basis * B];

	for(i=0; i<B; i++) {
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

void CalcGap(Solution *s, double *wg, double *ev, lapack_complex_double *es) {
	int i, j;
	double uplow = 100, dntop = -100, e = 0;

	// calc gap
	for(i=0; i<s->basis*G3; i++) {
		if(ev[i] > s->fermi) {
			if(ev[i] < uplow) uplow = ev[i];
		}
		else {
			if(ev[i] > dntop) dntop = ev[i];
		}
	}
	s->dntop = dntop;
	s->gap = uplow - dntop;

	// calc energy
	for(i=0; i<s->basis*G3; i++) {
		if(ev[i] < dntop) {
			for(j=0; j<s->basis; j++) e += CSQR(es[s->basis*i + j]) * wg[i/s->basis] * (ev[i] - dntop);
		}
	}
	s->e = e / pow(2*M_PI , 3);
}

inline void SymmetryF(double *n, double *m) { }
inline void SymmetryA(double *n, double *m) { n[2] = n[0];        m[2] = m[0]; }
inline void SymmetryC(double *n, double *m) { n[2] = n[1];        m[2] = m[1]; }
inline void SymmetryG(double *n, double *m) { n[2] = n[1] = n[0]; m[2] = m[1] = m[0]; }

void CalcSolution(Solution *s, LParameter *lp, lapack_complex_double *tbg) {
	FILE *f;
	char fs[256];

	GenName(s, "sol", fs);
	f = OpenFile(fs, "w");

	time_t t0 = time(NULL);

	int itr, i, j, is_cvg;
	char wgs[64], coefgs[64];
	double wg[G3], ev[s->basis * G3], fermi0, fermi1, n[OBT], m[OBT], ntot = 0, mtot = 0, cvg[OBT][3];
	Energy energy;
	lapack_complex_double es[s->basis*s->basis * G3], coefg[2 * G3];
	void (*GaussQuad)(double, double*, double*, lapack_complex_double*, double*, double*);
	void (*Symmetry)(double*, double*);

	GaussQuad = strstr(s->ctype, "s") ? GaussQuadS : GaussQuadD;

	if(strstr(s->type, "a")) Symmetry = SymmetryA;
	else if(strstr(s->type, "c")) Symmetry = SymmetryC;
	else if(strstr(s->type, "g")) Symmetry = SymmetryG;
	else Symmetry = SymmetryF;

	sprintf(wgs, "input/%s/wg.bin", s->name);
	sprintf(coefgs, "input/%s/coefg.bin", s->name);

	ReadBin(wgs, sizeof(double) * G3, wg);
	ReadBin(coefgs, sizeof(lapack_complex_double) * 2 * G3, coefg);

	for(i=0; i<OBT; i++) {
		for(j=0; j<3; j++) cvg[i][j] = 100;
	}

	fprintf(f, "%5s%12s", "#itr", "fermi");
	for(i=0; i<OBT; i++) fprintf(f, "%10s%02d%10s%02d", "n", i+1, "m", i+1);
	fprintf(f, "%12s%12s\n", "ntotal", "mtotal");

	fprintf(f, "%5d%12f", 0, s->fermi);
	for(i=0; i<OBT; i++) fprintf(f, "%12f%12f", s->n[i], s->m[i]);
	fprintf(f, "%12f%12f\n", s->ntot, s->mtot);

	for(itr=0; itr<50; itr++) {
		CalcEigen(s, lp, G3, tbg, ev, es, &energy);
		if(!strstr(s->ctype, "s")) TransBasis(G3, coefg, es);
		
		fermi0 = energy.min;
		fermi1 = energy.max;
		s->fermi = 0.5 * (fermi0 + fermi1);

		while(fabs(s->fermi - fermi1) > 1e-8) {
			ntot = 0;
			mtot = 0;

			GaussQuad(s->fermi, wg, ev, es, n, m);
		
			for(i=0; i<OBT; i++) {	
				n[i] /= pow(2*M_PI , 3);
				m[i] /= pow(2*M_PI , 3);
			}

			Symmetry(n, m);

			for(i=0; i<OBT; i++) {
				ntot += n[i];
				mtot += m[i];
			}

			if(fabs(ntot - s->N) < 1e-3) break;
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

		fprintf(f, "%5d%12f", itr+1, s->fermi);
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

	CalcGap(s, wg, ev, es);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fs, t1 - t0);
}

void ReadPathInfo(char *name, int p_len, int *p) {
	FILE *f;
	char fs[64];

	sprintf(fs, "input/%s/info.txt", name);
	f = OpenFile(fs, "r");

	int i, path, tmp = 0;
	char buf[1024];

	fgets(buf, sizeof(buf), f);

	for(i=0; i<p_len; i++) {
		fscanf(f, "%d", &path);	
		p[i] = path - tmp;
		tmp = path;
	}
	fclose(f);
}

void ReadLattice(char *name, int l_len, Lattice *l) {
	FILE *f;
	char fs[64];

	sprintf(fs, "input/%s/lattice.txt", name);
	f = OpenFile(fs, "r");

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

void ReadInfo(char *name, char *type, char *ctype, Vector *q) {
	FILE *f;
	char fs[64];

	sprintf(fs, "input/%s/info.txt", name);
	f = OpenFile(fs, "r");

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
	double ev[s->basis * s->B];
	Energy energy;
	lapack_complex_double es[s->basis*s->basis * s->B];

	CalcEigen(s, lp, s->B, tbb, ev, es, &energy);

	fprintf(f, "%4s", "#");
	for(i=0; i<s->basis; i++) fprintf(f, "%10s%02d", "e", i+1);
	fprintf(f, "\n");

	for(i=0; i<s->B; i++) {
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
	char coefbs[64];
	double p2;
	lapack_complex_double coefb[2 * s->B], p; 

	sprintf(coefbs, "input/%s/coefb.bin", s->name);
	ReadBin(coefbs, sizeof(lapack_complex_double) * 2 * s->B, coefb);
	TransBasis(s->B, coefb, es);

	fprintf(f, "%4s", "#");
	for(i=0; i<s->basis; i++) fprintf(f, "%10s%02d", "w", i+1);
	fprintf(f, "\n");

	for(i=0; i<s->B; i++) {
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
	double ev[s->basis * G3], e, dos[s->basis];
	Energy energy;
	lapack_complex_double es[s->basis*s->basis * G3];
	
	CalcEigen(s, lp, G3, tbg, ev, es, &energy);

	if(!strstr(s->ctype, "s")) {
		char coefgs[64];
		lapack_complex_double coefg[2 * G3];
		sprintf(coefgs, "input/%s/coefg.bin", s->name);
		ReadBin(coefgs, sizeof(lapack_complex_double) * 2 * G3, coefg);
		TransBasis(G3, coefg, es);
	}

	energy.min -= 1.0;
	energy.max += 1.0;
	
	fprintf(f, "%12s", "#e");
	for(i=0; i<s->basis; i++) fprintf(f, "%10s%02d", "d", i+1);
	fprintf(f, "\n");

	for(itv=0; itv<128; itv++) {
		memset(dos, 0, sizeof(dos));
		e = energy.min + (energy.max - energy.min) * itv / 128;

		for(i=0; i<G3*s->basis; i++) {
			for(j=0; j<s->basis; j++) dos[j] += CSQR(es[s->basis*i + j]) * GREEN(i);
		}

		fprintf(f, "%12f", e);
		for(i=0; i<s->basis; i++) fprintf(f, "%12f", dos[i] / pow(2*M_PI, 3));
		fprintf(f, "\n");
	}

	fclose(f);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fs, t1 - t0);
}

void FourierS(int l_len, int tb_num, Lattice *l, Vector g, Vector q, lapack_complex_double *tb) {
	int i, j, k = 0;
	double dot;
	lapack_complex_double tb_block[SINGLE][SINGLE] = {{0,}}, t, e;

	for(i=0; i<l_len; i++) {
		dot = 0;
		for(j=0; j<3; j++) dot += l[i].c[j] * g.c[j];

		t = l[i].tre + l[i].tim * I;
		e = cos(dot) + sin(dot) * I;

		tb_block[l[i].obt1-1][l[i].obt2-1] += t * e;
	}

	for(i=0; i<OBT; i++) {
		for(j=0; j<OBT; j++) {
			tb_block[i + OBT][j + OBT] = tb_block[i][j];
		}
	}

	for(i=0; i<SINGLE; i++) {
		for(j=0; j<SINGLE; j++) {
			tb[SINGLE*SINGLE*tb_num + k] = tb_block[i][j];
			k++;
		}
	}
}

void FourierD(int l_len, int tb_num, Lattice *l, Vector g, Vector q, lapack_complex_double *tb) {
	int i, j, k = 0;
	double dot_k, dot_q;
	lapack_complex_double tb_block[DOUBLE][DOUBLE] = {{0,}}, t, e_k, e_q;

	for(i=0; i<l_len; i++) {
		dot_k = dot_q = 0;
		for(j=0; j<3; j++) {
			dot_k += l[i].c[j] * (g.c[j]);
			dot_q += l[i].c[j] * (g.c[j] + q.c[j]);
		}

		t = l[i].tre + l[i].tim * I;
		e_k = cos(dot_k) + sin(dot_k) * I;
		e_q = cos(dot_q) + sin(dot_q) * I;

		tb_block[l[i].obt1-1][l[i].obt2-1] += t * e_k;
		tb_block[l[i].obt1-1 + OBT][l[i].obt2-1 + OBT] += t * e_q;
	}

	for(i=0; i<2*OBT; i++) {
		for(j=0; j<2*OBT; j++) {
			tb_block[i + 2*OBT][j + 2*OBT] = tb_block[i][j];
		}
	}

	for(i=0; i<DOUBLE; i++) {
		for(j=0; j<DOUBLE; j++) {
			tb[DOUBLE*DOUBLE*tb_num + k] = tb_block[i][j];
			k++;
		}
	}
}

void InteractionS(Solution *s, lapack_complex_double *tb_block) {
	int i, j;
	double n[OBT], n_[OBT], m[OBT], m_[OBT];
	
	memset(n_, 0, sizeof(n_));
	memset(m_, 0, sizeof(m_));

	for(i=0; i<OBT; i++) {
		n[i] = s->n[i];
		m[i] = s->m[i];

		for(j=0; j<OBT; j++) {
			if(j != i) { n_[i] += s->n[j]; m_[i] += s->m[j];
			}
		}		
	}

	for(i=0;              i<SINGLE*SINGLE;  i+=(SINGLE+1)) tb_block[i] += INTER_N;
	for(i=0;              i<SINGLE*OBT;     i+=(SINGLE+1)) tb_block[i] -= INTER_M;
	for(i=(SINGLE+1)*OBT; i<SINGLE*(2*OBT); i+=(SINGLE+1)) tb_block[i] += INTER_M;
}

void InteractionD(Solution *s, lapack_complex_double *tb_block) {
	int i, j;
	double n[OBT], n_[OBT], m[OBT], m_[OBT];
	
	memset(n_, 0, sizeof(n_));
	memset(m_, 0, sizeof(m_));

	for(i=0; i<OBT; i++) {
		n[i] = s->n[i] / 2;
		m[i] = s->m[i] / 2;

		for(j=0; j<OBT; j++) {
			if(j != i) {
				n_[i] += s->n[j] / 2;
				m_[i] += s->m[j] / 2;
			}
		}		
	}

	for(i=0;                        i<DOUBLE*DOUBLE;  i+=(DOUBLE+1)) tb_block[i] += INTER_N;
	for(i=OBT;                      i<DOUBLE*OBT;     i+=(DOUBLE+1)) tb_block[i] -= INTER_M;
	for(i=(DOUBLE+1)*OBT - OBT;     i<DOUBLE*(2*OBT); i+=(DOUBLE+1)) tb_block[i] -= INTER_M;
	for(i=(DOUBLE+1)*(2*OBT) + OBT; i<DOUBLE*(3*OBT); i+=(DOUBLE+1)) tb_block[i] += INTER_M;
	for(i=(DOUBLE+1)*(3*OBT) - OBT; i<DOUBLE*(4*OBT); i+=(DOUBLE+1)) tb_block[i] += INTER_M;

	/*
	for(i=0;                  i<DOUBLE*DOUBLE;  i+=(DOUBLE+1)) tb_block[i] += INTER_N;
	for(i=0;                  i<DOUBLE*OBT;     i+=(DOUBLE+1)) tb_block[i] -= INTER_M;
	for(i=(DOUBLE+1)*OBT;     i<DOUBLE*(2*OBT); i+=(DOUBLE+1)) tb_block[i] += INTER_M;
	for(i=(DOUBLE+1)*(2*OBT); i<DOUBLE*(3*OBT); i+=(DOUBLE+1)) tb_block[i] += INTER_M;
	for(i=(DOUBLE+1)*(3*OBT); i<DOUBLE*(4*OBT); i+=(DOUBLE+1)) tb_block[i] -= INTER_M;
	*/
}

void GaussQuadS(double fermi, double *wg, double *ev, lapack_complex_double *es, double *n, double *m) {
	int i, j;

	memset(n, 0, OBT * sizeof(double));
	memset(m, 0, OBT * sizeof(double));

	for(i=0; i<SINGLE*G3; i++) {
		if(ev[i] < fermi) {
			for(j=0; j<OBT; j++) {
				n[j] += (CSQR(es[SINGLE*i + j]) + CSQR(es[SINGLE*i + j+OBT])) * wg[i/SINGLE];
				m[j] += (CSQR(es[SINGLE*i + j]) - CSQR(es[SINGLE*i + j+OBT])) * wg[i/SINGLE];
			}
		}
	}
}

void GaussQuadD(double fermi, double *wg, double *ev, lapack_complex_double *es, double *n, double *m) {
	int i, j;

	memset(n, 0, OBT * sizeof(double));
	memset(m, 0, OBT * sizeof(double));

	for(i=0; i<DOUBLE*G3; i++) {
		if(ev[i] < fermi) {
			for(j=0; j<OBT; j++) {
				n[j] += (CSQR(es[DOUBLE*i + j]) + CSQR(es[DOUBLE*i + j+OBT]) + CSQR(es[DOUBLE*i + j+2*OBT]) + CSQR(es[DOUBLE*i + j+3*OBT])) * wg[i/DOUBLE];
				m[j] += (CSQR(es[DOUBLE*i + j]) - CSQR(es[DOUBLE*i + j+OBT])) * 2 * wg[i/DOUBLE];
			}
		}
	}
}
