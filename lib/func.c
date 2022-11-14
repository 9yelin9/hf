// lib/func.c : functions for calculating Hartree-Fock approxinameed 3-band model

#include "hf3.h" 

FILE* OpenFile(char *fn, char *mode) {
	FILE *f;

	if((f = fopen(fn, mode)) == NULL) {
		printf("%s fopen FAIL\n", fn);
		exit(1);
	}

	return f;
}

void CalcQuadPoints(Cell c) {
	FILE *fk, *fw;
	char fkn[64], fwn[64];

	sprintf(fkn, "input/%s/kg.bin", c.name);
	sprintf(fwn, "input/%s/wg.bin", c.name);
	fk = OpenFile(fkn, "wb");
	fw = OpenFile(fwn, "wb");

	time_t t0 = time(NULL);

	int i, i0, i1, i2;
	double k0[Nkg1], w0[Nkg1], w[Nkg];
	Coord k[Nkg];
	gsl_integration_glfixed_table *t = gsl_integration_glfixed_table_alloc(Nkg1);

	for(i=0; i<Nkg1; i++) gsl_integration_glfixed_point(-M_PI, M_PI, i, &k0[i], &w0[i], t);

	for(i=0; i<Nkg; i++) {
		i0 = i / (Nkg1 * Nkg1);
		i1 =(i / Nkg1) % Nkg1;
		i2 = i % Nkg1;

		k[i].c[0] = k0[i0];
		k[i].c[1] = k0[i1];
		k[i].c[2] = k0[i2];
		w[i] = w0[i0] * w0[i1] * w0[i2];
	}

	fwrite(k, sizeof(k), 1, fk);
	fwrite(w, sizeof(w), 1, fw);

	fclose(fk);
	fclose(fw);

	time_t t1 = time(NULL);
	printf("%s(%s, %s) : %lds\n", __func__, fkn, fwn, t1 - t0);
}

void ReadBin(char *fn, int size, void *v) {
	FILE *f;

	f = OpenFile(fn, "rb");
	fread(v, size, 1, f);
	fclose(f);
}

void ReadCell(Cell *c) {
	FILE *f;
	char fn[64];

	sprintf(fn, "input/%s/cell.txt", c->name);
	f = OpenFile(fn, "r");

	int i;
	char buf[1024], type[16];

	while(!feof(f)) {
		fgets(buf, sizeof(buf), f);
		sscanf(buf, "%s%s%s%d%d%lf%lf%lf",\
				type, c->sys, c->bas, &c->Ni, &c->Nc, &c->q.c[0], &c->q.c[1], &c->q.c[2]);
		if(strstr(type, c->type)) break;
	}
	fclose(f);

	c->Ns  = c->Ni * c->Nc;
	c->Nb  = c->Ni * c->Nc * 2;
	c->Nbb = c->Nb * c->Nb;
	for(i=0; i<3; i++) c->q.c[i] *= M_PI;
}

void ReadLat(char *name, int *Nl, Lattice *l) {
	FILE *f;
	char fn[64];

	sprintf(fn, "input/%s/lat.txt", name);
	f = OpenFile(fn, "r");

	int i;
	char buf[1024];

	fgets(buf, sizeof(buf), f);
	sscanf(buf, "%d", Nl);

	for(i=0; i<*Nl; i++) {
		fscanf(f, "%lf%lf%lf%d%d%lf%lf",\
				&l[i].d.c[0], &l[i].d.c[1], &l[i].d.c[2], &l[i].obi, &l[i].obf, &l[i].tre, &l[i].tim);
	}
	fclose(f);
}

void Fourier0(Cell c, int ith, Coord k, int Nl, Lattice *l, lapack_complex_double *tb) {
	int i, j, v = 0;
	double dot;
	lapack_complex_double tb0[c.Nb][c.Nb], t, e;

	memset(tb0, 0, sizeof(tb0));

	for(i=0; i<Nl; i++) {
		dot = 0;
		for(j=0; j<3; j++) {
			dot += l[i].d.c[j] * k.c[j];
		}

		t = l[i].tre + l[i].tim * I;
		e = cos(dot) + sin(dot) * I;

		tb0[l[i].obi-1][l[i].obf-1] += t * e;
	}

	for(i=0; i<c.Ns; i++) {
		for(j=0; j<c.Ns; j++) {
			tb0[i + c.Ns][j + c.Ns] = tb0[i][j];
		}
	}

	for(i=0; i<c.Nb; i++) {
		for(j=0; j<c.Nb; j++) {
			tb[c.Nbb*ith + v] = tb0[i][j];
			v++;
		}
	}
}	

void FourierQ(Cell c, int ith, Coord k, int Nl, Lattice *l, lapack_complex_double *tb) {
	int i, j, v = 0;
	double dot, dotQ;
	lapack_complex_double tb0[c.Nb][c.Nb], t, e, eQ;

	memset(tb0, 0, sizeof(tb0));

	for(i=0; i<Nl; i++) {
		dot = dotQ = 0;
		for(j=0; j<3; j++) {
			dot  += l[i].d.c[j] * (k.c[j]);
			dotQ += l[i].d.c[j] * (k.c[j] + c.q.c[j]);
		}

		t = l[i].tre + l[i].tim * I;
		e  = cos(dot)  + sin(dot)  * I;
		eQ = cos(dotQ) + sin(dotQ) * I;

		tb0[l[i].obi-1][l[i].obf-1] += t * e;
		tb0[l[i].obi-1 + c.Nc][l[i].obf-1 + c.Nc] += t * eQ;
	}

	for(i=0; i<c.Ns; i++) {
		for(j=0; j<c.Ns; j++) {
			tb0[i + c.Ns][j + c.Ns] = tb0[i][j];
		}
	}

	for(i=0; i<c.Nb; i++) {
		for(j=0; j<c.Nb; j++) {
			tb[c.Nbb*ith + v] = tb0[i][j];
			v++;
		}
	}
}

void DotProd(Cell c, int Nk, Coord *k, Coord *r, lapack_complex_double *cf) {
	int i, j, l;
	double dot;

	for(i=0; i<Nk; i++) {
		for(j=0; j<c.Ni; j++) {
			dot = 0;
			for(l=0; l<DIM; l++) dot += r[j].c[l] * k[i].c[l];
			cf[c.Ni*i + j] = cos(dot) - sin(dot) * I;
		}
	}
}

void CalcCoef(Cell c, Coord *r) {
	FILE *fg, *fb;
	char fgn[64], fbn[64];

	sprintf(fgn, "input/%s/cfg.bin", c.name);
	sprintf(fbn, "input/%s/cfb.bin", c.name);
	fg = OpenFile(fgn, "wb");
	fb = OpenFile(fbn, "wb");

	time_t t0 = time(NULL);

	char kgn[64], kbn[64];
	Coord kg[Nkg], kb[Nkb];
	lapack_complex_double cfg[c.Ni * Nkg], cfb[c.Ni * Nkb];

	sprintf(kgn, "input/%s/kg.bin", c.name);
	sprintf(kbn, "input/%s/kb.bin", c.name);

	ReadBin(kgn, sizeof(Coord) * Nkg, kg);
	ReadBin(kbn, sizeof(Coord) * Nkb, kb);

	DotProd(c, Nkg, kg, r, cfg);
	DotProd(c, Nkb, kb, r, cfb);

	fwrite(cfg, sizeof(cfg), 1, fg);
	fwrite(cfb, sizeof(cfb), 1, fb);

	fclose(fg);
	fclose(fb);

	time_t t1 = time(NULL);
	printf("%s(%s, %s) : %lds\n", __func__, fgn, fbn, t1 - t0);
}

void CalcTB(Cell c, void (*Fourier)()) {
	FILE *fg, *fb;
	char fgn[64], fbn[64];

	sprintf(fgn, "input/%s/tbg_%c.bin", c.name, c.type[0]);
	sprintf(fbn, "input/%s/tbb_%c.bin", c.name, c.type[0]);
	fg = OpenFile(fgn, "wb");
	fb = OpenFile(fbn, "wb");

	time_t t0 = time(NULL);

	int i, j;
	char kgn[64], kbn[64];
	Coord kg[Nkg], kb[Nkb];
	lapack_complex_double tbg[c.Nbb * Nkg], tbb[c.Nbb * Nkb];

	sprintf(kgn, "input/%s/kg.bin", c.name);
	sprintf(kbn, "input/%s/kb.bin", c.name);

	ReadBin(kgn, sizeof(Coord) * Nkg, kg);
	ReadBin(kbn, sizeof(Coord) * Nkb, kb);

	int Nl;
	Lattice l[Nl_MAX];
	ReadLat(c.name, &Nl, l);

	omp_set_num_threads(OMP_THREAD);
#pragma omp parallel for ordered
	for(i=0; i<Nkg; i++) Fourier(c, i, kg[i], Nl, l, tbg);
#pragma omp parallel for ordered              
	for(i=0; i<Nkb; i++) Fourier(c, i, kb[i], Nl, l, tbb);

	fwrite(tbg, sizeof(tbg), 1, fg);
	fwrite(tbb, sizeof(tbb), 1, fb);

	fclose(fg);
	fclose(fb);

	time_t t1 = time(NULL);
	printf("%s(%s, %s) : %lds\n", __func__, fgn, fbn, t1 - t0);

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
}

void DataName(Cell c, Solution *s, char *dtype, char *dn) {
	if(strstr(dtype, "sol")) {
		sprintf(dn, "output/%s/JU%.2f_SOC%.2f/%s_%s_N%.1f_U%.1f_%s.txt",\
				c.name, s->JU, s->SOC, dtype, c.type, s->N, s->U, s->runtime);
	}
	else {
		sprintf(dn, "output/%s/JU%.2f_SOC%.2f/%s_%s_N%.1f_U%.1f_n%f_m%f_e%f_gap%f_fermi%f_dntop%f_%s.txt",\
				c.name, s->JU, s->SOC, dtype, c.type, s->N, s->U, s->ns, s->ms, s->e, s->gap, s->fermi, s->dntop, s->runtime);
	}
}

void CalcEigen(Cell c, Solution *s, LAPACK *lp, int Nk, lapack_complex_double *tb, double *ev, lapack_complex_double *es, Energy *e, void (*Interaction)()) {
	e->min =  100;
	e->max = -100;

	int i, j;
	double ev0[c.Nb];
	lapack_complex_double es0[c.Nbb];

	for(i=0; i<Nk; i++) {
		for(j=0; j<c.Nbb; j++) es0[j] = tb[c.Nbb*i + j];

		Interaction(c, s, es0);

		LAPACK_zheev(&lp->jobz, &lp->uplo, &lp->ln, es0, &lp->lda, ev0, lp->work, &lp->lwork, lp->rwork, &lp->info);
		if(lp->info != 0) {
			printf("LAPACK_zheev FAIL\n");
			exit(1);
		}

		for(j=0; j<c.Nb;  j++) ev[c.Nb*i  + j] = ev0[j];
		for(j=0; j<c.Nbb; j++) es[c.Nbb*i + j] = es0[j];

		if(ev0[0]      < e->min) e->min = ev0[0];
		if(ev0[c.Nb-1] > e->max) e->max = ev0[c.Nb-1];
	}
}

void CalcGap(Cell c, Solution *s, double *ev, lapack_complex_double *es) {
	int i, j;
	char wgn[64];
	double wg[Nkg], uplow = 100, dntop = -100, e = 0, tol = 1e-5;

	sprintf(wgn, "input/%s/wg.bin", c.name);
	ReadBin(wgn, sizeof(double) * Nkg, wg);

	// calc gap
	for(i=0; i<c.Nb*Nkg; i++) {
		if(ev[i] - s->fermi > tol) {
			if(ev[i] < uplow) uplow = ev[i];
		}
		else if(ev[i] - s->fermi < -tol){
			if(ev[i] > dntop) dntop = ev[i];
		}
	}
	s->dntop = dntop;
	s->gap = uplow - dntop;

	// calc energy
	for(i=0; i<c.Nb*Nkg; i++) {
		if(ev[i] < dntop) {
			for(j=0; j<c.Nb; j++) e += CSQR(es[c.Nb*i + j]) * wg[i/c.Nb] * ev[i];
		}
	}
	s->e = e / pow(2*M_PI , 3);
}

#define OBT_IDX ((i / c.Nb) % c.Nc)
#define INTER_N (0.5 * ((s->U) * n[OBT_IDX] + (s->U - 2*s->J) * n_[OBT_IDX] + (s->U - 3*s->J) * n_[OBT_IDX]))
#define INTER_M (0.5 * ((s->U) * m[OBT_IDX] + (s->U - 2*s->J) * m_[OBT_IDX] - (s->U - 3*s->J) * m_[OBT_IDX])) 

void Interaction0(Cell c, Solution *s, lapack_complex_double *tb0) {
	int i, j;
	double n[c.Nc], n_[c.Nc], m[c.Nc], m_[c.Nc];
	
	memset(n_, 0, sizeof(n_));
	memset(m_, 0, sizeof(m_));

	for(i=0; i<c.Nc; i++) {
		n[i] = s->n[i] / c.Ni;
		m[i] = s->m[i] / c.Ni;

		for(j=0; j<c.Nc; j++) {
			if(j != i) {
				n_[i] += s->n[j] / c.Ni;
				m_[i] += s->m[j] / c.Ni;
			}
		}		
	}

	for(i=0; i<c.Nbb; i+=c.Nb+1) tb0[i] += INTER_N;

	for(i=0;             i<c.Nb*c.Ns; i+=c.Nb+1) tb0[i] -= INTER_M;
	for(i=(c.Nb+1)*c.Ns; i<c.Nbb;     i+=c.Nb+1) tb0[i] += INTER_M;
}

void InteractionS(Cell c, Solution *s, lapack_complex_double *tb0) {
	int i, j;
	double n[c.Nc], n_[c.Nc], m[c.Nc], m_[c.Nc];
	
	memset(n_, 0, sizeof(n_));
	memset(m_, 0, sizeof(m_));

	for(i=0; i<c.Nc; i++) {
		n[i] = s->n[i] / c.Ni;
		m[i] = s->m[i] / c.Ni;

		for(j=0; j<c.Nc; j++) {
			if(j != i) {
				n_[i] += s->n[j] / c.Ni;
				m_[i] += s->m[j] / c.Ni;
			}
		}		
	}

	for(i=0; i<c.Nbb; i+=c.Nb+1) tb0[i] += INTER_N;

	for(i=0;             i<c.Nb*c.Ns; i+=c.Nb+1) tb0[i] -= INTER_M * pow(-1, i / (c.Nb*c.Nc));
	for(i=(c.Nb+1)*c.Ns; i<c.Nbb;     i+=c.Nb+1) tb0[i] += INTER_M * pow(-1, (i + c.Nb*c.Ns)/ (c.Nb*c.Nc)) ;
}

void InteractionQ(Cell c, Solution *s, lapack_complex_double *tb0) {
	int i, j;
	double n[c.Nc], n_[c.Nc], m[c.Nc], m_[c.Nc];
	
	memset(n_, 0, sizeof(n_));
	memset(m_, 0, sizeof(m_));

	for(i=0; i<c.Nc; i++) {
		n[i] = s->n[i] / c.Ni;
		m[i] = s->m[i] / c.Ni;

		for(j=0; j<c.Nc; j++) {
			if(j != i) {
				n_[i] += s->n[j] / c.Ni;
				m_[i] += s->m[j] / c.Ni;
			}
		}		
	}

	for(i=0; i<c.Nbb; i+=c.Nb+1) tb0[i] += INTER_N;

	for(i=c.Nc;                 i<c.Nb*c.Nc; i+=c.Nb+1) tb0[i] -= INTER_M;
	for(i=(c.Nb+1)*c.Ns + c.Nc; i<c.Nbb;     i+=c.Nb+1) tb0[i] += INTER_M;
}

#define STATE_IDX (c.Nbb*i + c.Nb*j + c.Nc*(l/c.Nc) + l)

void Basis0(Cell c, int Nk, lapack_complex_double *cf, lapack_complex_double *es) {}
void BasisQ(Cell c, int Nk, lapack_complex_double *cf, lapack_complex_double *es) {
	int i, j, l;
	lapack_complex_double es0[c.Ni];

	for(i=0; i<Nk; i++) {
		for(j=0; j<c.Nb; j++) {
			for(l=0; l<c.Ns; l++) {
				es0[0] = sqrt(0.5) * cf[c.Ni*i + 1] * (es[STATE_IDX] + es[STATE_IDX + c.Nc]);
				es0[1] = sqrt(0.5) * cf[c.Ni*i + 0] * (es[STATE_IDX] - es[STATE_IDX + c.Nc]);

				es[STATE_IDX]        = es0[0];
				es[STATE_IDX + c.Nc] = es0[1];
			}
		}
	}
}

#define OC_IDX (c.Nc*i + j)

void Quadrature(Cell c, Solution *s, double *wg, double *ev, lapack_complex_double *es, double *n, double *m) {
	int i, j;
	double oc[c.Nb], mt[c.Nc], ms[c.Nc];

	memset(n,  0, sizeof(double) * c.Nc);
	memset(mt, 0, sizeof(double) * c.Nc);
	memset(ms, 0, sizeof(double) * c.Nc);
	memset(oc, 0, sizeof(double) * c.Nb);

	for(i=0; i<c.Nb*Nkg; i++) {
		if(ev[i] < s->fermi) {
			for(j=0; j<c.Nb; j++) {
				oc[j] += CSQR(es[c.Nb*i + j]) * wg[i / c.Nb];
			}
		}	
	}

	for(i=0; i<c.Nb; i++) oc[i] /= pow(2*M_PI, 3);

	for(i=0; i<c.Ni; i++) {
		for(j=0; j<c.Nc; j++) {
			n[j]  +=  oc[OC_IDX] + oc[OC_IDX + c.Ns];
			mt[j] +=  oc[OC_IDX] - oc[OC_IDX + c.Ns];
			ms[j] += (oc[OC_IDX] - oc[OC_IDX + c.Ns]) * pow(-1, i);
		}
	}

	if(strstr(c.type, "f")) for(i=0; i<c.Nc; i++) m[i] = mt[i];
	else                    for(i=0; i<c.Nc; i++) m[i] = ms[i];
}

void System0(double *n, double *m) {}
void SystemScA(double *n, double *m) {n[2] = n[0];        m[2] = m[0];}
void SystemScC(double *n, double *m) {n[2] = n[1];        m[2] = m[1];}
void SystemScG(double *n, double *m) {n[2] = n[1] = n[0]; m[2] = m[1] = m[0];}

void CalcSolution(Cell c, Solution *s, LAPACK *lp, void (*System)(), void (*Interaction)(), void (*Basis)()) {
	FILE *f;
	char fn[256];

	DataName(c, s, "sol", fn);
	f = OpenFile(fn, "w");

	time_t t0 = time(NULL);

	int itr, i, j, is_cvg;
	char wgn[64], cfgn[64], tbgn[64];
	double wg[Nkg], ev[c.Nb * Nkg], fermi0, fermi1, n[c.Nc], m[c.Nc], ns = 0, ms = 0, cvg[c.Nc][CVG_MAX], avg;
	Energy energy;
	lapack_complex_double cfg[c.Ni * Nkg], tbg[c.Nbb * Nkg], es[c.Nbb * Nkg];

	sprintf(wgn,  "input/%s/wg.bin", c.name);
	sprintf(cfgn, "input/%s/cfg.bin", c.name);
	sprintf(tbgn, "input/%s/tbg_%c.bin", c.name, c.type[0]);
	ReadBin(wgn,  sizeof(double) * Nkg, wg);
	ReadBin(cfgn, sizeof(lapack_complex_double) * c.Ni  * Nkg, cfg);
	ReadBin(tbgn, sizeof(lapack_complex_double) * c.Nbb * Nkg, tbg);

	for(i=0; i<c.Nc; i++) {
		for(j=0; j<CVG_MAX; j++) {
			cvg[i][j] = 100;
		}
	}

	fprintf(f, "%5s%12s", "#itr", "fermi");
	for(i=0; i<c.Nc; i++) fprintf(f, "%10s%02d%10s%02d", "n", i+1, "m", i+1);
	fprintf(f, "%12s%12s\n", "ns", "ms");

	fprintf(f, "%5d%12f", 0, s->fermi);
	for(i=0; i<c.Nc; i++) fprintf(f, "%12f%12f", s->n[i], s->m[i]);
	fprintf(f, "%12f%12f\n", s->ns, s->ms);

	for(itr=0; itr<25; itr++) {
		CalcEigen(c, s, lp, Nkg, tbg, ev, es, &energy, Interaction);
		Basis(c, Nkg, cfg, es);
		
		fermi0 = energy.min;
		fermi1 = energy.max;
		s->fermi = 0.5 * (fermi0 + fermi1);

		while(fabs(s->fermi - fermi1) > 1e-8) {
			ns = 0;
			ms = 0;

			Quadrature(c, s, wg, ev, es, n, m);
			System(n, m);

			for(i=0; i<c.Nc; i++) {
				ns += n[i];
				ms += m[i];
			}

			if(fabs(ns - s->N) < 1e-4) break;
			else {
				if(ns < s->N) fermi0 = s->fermi;
				else fermi1 = s->fermi;
				s->fermi = 0.5 * (fermi0 + fermi1);
			}
		}

		for(i=0; i<c.Nc; i++) {
			s->n[i] = n[i];
			s->m[i] = m[i];
		}
		s->ns = ns;
		s->ms = ms;

		fprintf(f, "%5d%12f", itr+1, s->fermi);
		for(i=0; i<c.Nc; i++) fprintf(f, "%12f%12f", s->n[i], s->m[i]);
		fprintf(f, "%12f%12f\n", s->ns, s->ms);
			
		is_cvg = 0;
		for(i=0; i<c.Nc; i++) {
			cvg[i][itr % CVG_MAX] = m[i];

			avg = 0;
			for(j=0; j<CVG_MAX; j++) avg += cvg[i][j];
			avg /= CVG_MAX;

			if(fabs(avg - m[i]) < 1e-6) is_cvg++;
		}
		if(is_cvg == c.Nc) break;
	}
	fclose(f);

	CalcGap(c, s, ev, es);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fn, t1 - t0);
}

void MakeBand(Cell c, Solution *s, LAPACK *lp, void (*Interaction)(), void (*Basis)()) {
	FILE *f;
	char fn[256];

	DataName(c, s, "band", fn);
	f = OpenFile(fn, "w");

	time_t t0 = time(NULL);

	int i, j;
	char tbbn[64];
	double ev[c.Nb * Nkb];
	Energy energy;
	lapack_complex_double tbb[c.Nbb * Nkb], es[c.Nbb * Nkb];

	sprintf(tbbn, "input/%s/tbb_%c.bin", c.name, c.type[0]);
	ReadBin(tbbn, sizeof(lapack_complex_double) * c.Nbb * Nkb, tbb);

	CalcEigen(c, s, lp, Nkb, tbb, ev, es, &energy, Interaction);

	fprintf(f, "%10s%02d", "#e", 0);
	for(i=1; i<c.Nb; i++) fprintf(f, "%10s%02d", "e", i+1);
	fprintf(f, "\n");

	for(i=0; i<Nkb; i++) {
		for(j=0; j<c.Nb; j++) fprintf(f, "%12f", ev[c.Nb*i + j]);
		fprintf(f, "\n");
	}

	fclose(f);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fn, t1 - t0);

	if(strstr(c.bas, "q")) MakeUFW(c, s, ev, es, Basis);
}

void MakeUFW(Cell c, Solution *s, double *ev, lapack_complex_double *es, void (*Basis)()) {
	FILE *f;
	char fn[256];

	DataName(c, s, "ufw", fn);
	f = OpenFile(fn, "w");

	time_t t0 = time(NULL);

	int i, j, l;
	char cfbn[64];
	double p2;
	lapack_complex_double cfb[c.Ni * Nkb], p; 

	sprintf(cfbn, "input/%s/cfb.bin", c.name);
	ReadBin(cfbn, sizeof(lapack_complex_double) * c.Ni * Nkb, cfb);
	Basis(c, Nkb, cfb, es);

	fprintf(f, "%10s%02d", "#w", 0);
	for(i=1; i<c.Nb; i++) fprintf(f, "%10s%02d", "w", i+1);
	fprintf(f, "\n");

	for(i=0; i<Nkb; i++) {
		for(j=0; j<c.Nb; j++) {
			p2 = 0;
			for(l=0; l<c.Ns; l++) {
				p = es[STATE_IDX] * cfb[c.Ni*i + 0] + es[STATE_IDX + c.Nc] * cfb[c.Ni*i + 1];
				p2 += CSQR(p) / c.Ni;
			}
			fprintf(f, "%12f", p2); 
		}
		fprintf(f, "\n");
	}

	fclose(f);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fn, t1 - t0);
}

#define WIDTH 0.05
#define GREEN(i) (WIDTH / (pow(e - ev[i], 2) + pow(WIDTH, 2))) 

void MakeDOS(Cell c, Solution *s, LAPACK *lp, void (*Interaction)(), void (*Basis)()) {
	FILE *f;
	char fn[256];

	DataName(c, s, "dos", fn);
	f = OpenFile(fn, "w");

	time_t t0 = time(NULL);

	int itv, i, j;
	char cfgn[64], tbgn[64];
	double ev[c.Nb * Nkg], e, dos[c.Nb];
	Energy energy;
	lapack_complex_double cfg[c.Ni * Nkg], tbg[c.Nbb * Nkg], es[c.Nbb * Nkg];
	
	sprintf(cfgn, "input/%s/cfg.bin", c.name);
	sprintf(tbgn, "input/%s/tbg_%c.bin", c.name, c.type[0]);
	ReadBin(cfgn, sizeof(lapack_complex_double) * c.Ni  * Nkg, cfg);
	ReadBin(tbgn, sizeof(lapack_complex_double) * c.Nbb * Nkg, tbg);

	CalcEigen(c, s, lp, Nkg, tbg, ev, es, &energy, Interaction);
	Basis(c, Nkg, cfg, es);

	energy.min -= 1.0;
	energy.max += 1.0;
	
	fprintf(f, "%12s", "#e");
	for(i=0; i<c.Nb; i++) fprintf(f, "%10s%02d", "d", i+1);
	fprintf(f, "\n");

	for(itv=0; itv<128; itv++) {
		memset(dos, 0, sizeof(dos));
		e = energy.min + (energy.max - energy.min) * itv / 128;

		for(i=0; i<c.Nb*Nkg; i++) {
			for(j=0; j<c.Nb; j++) dos[j] += CSQR(es[c.Nb*i + j]) * GREEN(i);
		}

		fprintf(f, "%12f", e);
		for(i=0; i<c.Nb; i++) fprintf(f, "%12f", dos[i] / pow(2*M_PI, 3));
		fprintf(f, "\n");
	}

	fclose(f);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fn, t1 - t0);
}
