// lib/hf.c : functions for Hartree-Fock approximation

#include "hf.h" 

void ReadInfo(Config *c) {
	FILE *f;
	char fn[64];

	sprintf(fn, "input/%s/config_%s.txt", c->name, c->type);
	f = OpenFile(fn, "r");

	int i;
	char buf[1024], type[16];

	while(!feof(f)) {
		fgets(buf, sizeof(buf), f);
		sscanf(buf, "%s%s%s%d%d%lf%lf%lf",\
				type, c->sys, c->bas, &c->Ni, &c->Nc, &c->q.c[0], &c->q.c[1], &c->q.c[2]);
		if(strstr(c->type, type)) break;
	}
	fclose(f);

	c->Ns  = c->Ni * c->Nc;
	c->Nb  = c->Ni * c->Nc * 2;
	c->Nbb = c->Nb * c->Nb;
	for(i=0; i<3; i++) c->q.c[i] *= M_PI;
}

void ReadLattice(char *name, int Nl, Lattice *l, char *ltype) {
	FILE *f;
	char fn[64];

	sprintf(fn, "input/%s/%slat.txt", name, ltype);
	f = OpenFile(fn, "r");

	int i;
	char buf[1024];

	fgets(buf, sizeof(buf), f);
	for(i=0; i<Nl; i++) {
		fscanf(f, "%lf%lf%lf%d%d%lf%lf",\
				&l[i].d.c[0], &l[i].d.c[1], &l[i].d.c[2], &l[i].obi, &l[i].obf, &l[i].tre, &l[i].tim);
	}
	fclose(f);
}

void FourierN(Config c, int ith, Coord k, int Nl, Lattice *l, lapack_complex_double *tb) {
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

void FourierQ(Config c, int ith, Coord k, int Nl, Lattice *l, lapack_complex_double *tb) {
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

void CalcTB(Config c, char *ltype, void (*Fourier)()) {
	FILE *fl, *fg, *fb;
	char fln[64], fgn[64], fbn[64];

	sprintf(fln, "input/%s/%slat.txt", c.name, ltype);
	sprintf(fgn, "input/%s/tbg_%c.bin", c.name, c.type[0]);
	sprintf(fbn, "input/%s/tbb_%c.bin", c.name, c.type[0]);
	fl = OpenFile(fln, "r");
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
	char buf[16];
	Lattice *l;

	fgets(buf, sizeof(buf), fl);
	sscanf(buf, "%d", &Nl);
	fclose(fl);

	l = (Lattice*)malloc(sizeof(Lattice) * Nl);
	ReadLat(c.name, Nl, l, ltype);
	
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

void DataName(Config c, Solution *s, char *dtype, char *dn) {
	if(strstr(dtype, "sol")) {
		sprintf(dn, "output/%s/JU%.2f_SOC%.2f/%s_%s_N%.1f_U%.1f_%s.txt",\
				c.save, s->JU, s->SOC, dtype, c.type, s->N, s->U, s->runtime);
	}
	else if(strstr(dtype, "dos")) {
		sprintf(dn, "output/%s/JU%.2f_SOC%.2f/%s_%s_eta%.2f_N%.1f_U%.1f_n%.16f_m%.16f_e%.16f_gap%.16f_fermi%.16f_dntop%.16f_%s.txt",\
				c.save, s->JU, s->SOC, dtype, c.type, c.eta, s->N, s->U, s->ns, s->ms, s->e, s->gap, s->fermi, s->dntop, s->runtime);
	}
	else {
		sprintf(dn, "output/%s/JU%.2f_SOC%.2f/%s_%s_N%.1f_U%.1f_n%.16f_m%.16f_e%.16f_gap%.16f_fermi%.16f_dntop%.16f_%s.txt",\
				c.save, s->JU, s->SOC, dtype, c.type, s->N, s->U, s->ns, s->ms, s->e, s->gap, s->fermi, s->dntop, s->runtime);
	}
}

void CalcEigen(Config c, Solution *s, LAPACK *lp, int Nk, lapack_complex_double *tb, double *ev, lapack_complex_double *es, Energy *e, void (*Interaction)()) {
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

void CalcGap(Config c, Solution *s, double *ev, lapack_complex_double *es) {
	int i, j;
	char wgn[64];
	double wg[Nkg], uplow = 100, dntop = -100, e = 0, tol = 1e-6;

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

void InteractionN(Config c, Solution *s, lapack_complex_double *tb0) {
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

void InteractionS(Config c, Solution *s, lapack_complex_double *tb0) {
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
	for(i=(c.Nb+1)*c.Ns; i<c.Nbb;     i+=c.Nb+1) tb0[i] += INTER_M * pow(-1, (i + c.Nb*c.Ns) / (c.Nb*c.Nc)) ;
}

void InteractionQ(Config c, Solution *s, lapack_complex_double *tb0) {
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

void BasisN(Config c, int Nk, lapack_complex_double *cf, lapack_complex_double *es) {}
void BasisQ(Config c, int Nk, lapack_complex_double *cf, lapack_complex_double *es) {
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

void Quadrature(Config c, Solution *s, double *wg, double *ev, lapack_complex_double *es, double *oc) {
	int i, j;

	memset(oc, 0, sizeof(double) * c.Nb);

	for(i=0; i<c.Nb*Nkg; i++) {
		if(ev[i] < s->fermi) {
			for(j=0; j<c.Nb; j++) {
				oc[j] += CSQR(es[c.Nb*i + j]) * wg[i / c.Nb];
			}
		}	
	}

	for(i=0; i<c.Nb; i++) oc[i] /= pow(2*M_PI, 3);
}

void SystemN(double *n, double *m) {}
void SystemScA(double *n, double *m) {n[2] = n[0];        m[2] = m[0];}
void SystemScC(double *n, double *m) {n[2] = n[1];        m[2] = m[1];}
void SystemScG(double *n, double *m) {n[2] = n[1] = n[0]; m[2] = m[1] = m[0];}

#define OC_IDX (c.Nc*i + j)

void CalcSolution(Config c, Solution *s, LAPACK *lp, void (*System)(), void (*Interaction)(), void (*Basis)()) {
	FILE *f;
	char fn[256];

	DataName(c, s, "sol", fn);
	f = OpenFile(fn, "w");

	time_t t0 = time(NULL);

	int itr, i, j, one, is_cvg;
	char wgn[64], cfgn[64], tbgn[64];
	double wg[Nkg], ev[c.Nb * Nkg], fermi0, fermi1, oc[c.Nb], n[c.Nc], m[c.Nc], ns = 0, ms = 0, cvg[c.Nc][CVG_MAX], avg;
	Energy energy;
	lapack_complex_double cfg[c.Ni * Nkg], tbg[c.Nbb * Nkg], es[c.Nbb * Nkg];

	if(strstr(c.type, "f")) one = 1;
	else one = -1;

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

	fprintf(f, "%22s", "#fermi");
	for(i=0; i<c.Nb; i++) fprintf(f, "%20s%02d", "oc", i+1);
	fprintf(f, "\n");

	for(itr=0; itr<25; itr++) {
		CalcEigen(c, s, lp, Nkg, tbg, ev, es, &energy, Interaction);
		Basis(c, Nkg, cfg, es);
		
		fermi0 = energy.min;
		fermi1 = energy.max;
		s->fermi = 0.5 * (fermi0 + fermi1);

		while(fabs(s->fermi - fermi1) > 1e-8) {
			memset(n, 0, sizeof(double) * c.Nc);
			memset(m, 0, sizeof(double) * c.Nc);
			ns = 0;
			ms = 0;

			Quadrature(c, s, wg, ev, es, oc);
			for(i=0; i<c.Ni; i++) {
				for(j=0; j<c.Nc; j++) {
					n[j] +=  oc[OC_IDX] + oc[OC_IDX + c.Ns];
					m[j] += (oc[OC_IDX] - oc[OC_IDX + c.Ns]) * pow(one, i);
				}
			}
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

		fprintf(f, "%22.16f", s->fermi);
		for(i=0; i<c.Nb; i++) fprintf(f, "%22.16f", oc[i]);
		fprintf(f, "\n");
			
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

void MakeBand(Config c, Solution *s, LAPACK *lp, void (*Interaction)(), void (*Basis)()) {
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

	fprintf(f, "%20s%02d", "#e", 0);
	for(i=0; i<c.Nb; i++) fprintf(f, "%20s%02d", "e", i+1);
	fprintf(f, "\n");

	for(i=0; i<Nkb; i++) {
		for(j=0; j<c.Nb; j++) fprintf(f, "%22.16f", ev[c.Nb*i + j]);
		fprintf(f, "\n");
	}

	fclose(f);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fn, t1 - t0);

	if(strstr(c.bas, "q")) MakeUFW(c, s, ev, es, Basis);
}

void MakeUFW(Config c, Solution *s, double *ev, lapack_complex_double *es, void (*Basis)()) {
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

	fprintf(f, "%20s%02d", "#w", 0);
	for(i=0; i<c.Nb; i++) fprintf(f, "%20s%02d", "w", i+1);
	fprintf(f, "\n");

	for(i=0; i<Nkb; i++) {
		for(j=0; j<c.Nb; j++) {
			p2 = 0;
			for(l=0; l<c.Ns; l++) {
				p = es[STATE_IDX] * cfb[c.Ni*i + 0] + es[STATE_IDX + c.Nc] * cfb[c.Ni*i + 1];
				p2 += CSQR(p) / c.Ni;
			}
			fprintf(f, "%22.16f", p2); 
		}
		fprintf(f, "\n");
	}

	fclose(f);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fn, t1 - t0);
}

#define GREEN(i) (c.eta / (pow(e - ev[i] + s->fermi, 2) + pow(c.eta, 2))) 

void MakeDOS(Config c, Solution *s, LAPACK *lp, void (*Interaction)(), void (*Basis)()) {
	FILE *f;
	char fn[256];

	DataName(c, s, "dos", fn);
	f = OpenFile(fn, "w");

	time_t t0 = time(NULL);

	int itv, v, i, j;
	char cfgn[64], tbgn[64];
	double ev[c.Nb * Nkg], dos[c.Nb], e, emin, emax, vol = pow(2*M_PI, 3) * pow(2*M_PI, 3) * M_PI; 
	Energy energy;
	lapack_complex_double cfg[c.Ni * Nkg], tbg[c.Nbb * Nkg], es[c.Nbb * Nkg];
	
	sprintf(cfgn, "input/%s/cfg.bin", c.name);
	sprintf(tbgn, "input/%s/tbg_%c.bin", c.name, c.type[0]);
	ReadBin(cfgn, sizeof(lapack_complex_double) * c.Ni  * Nkg, cfg);
	ReadBin(tbgn, sizeof(lapack_complex_double) * c.Nbb * Nkg, tbg);

	CalcEigen(c, s, lp, Nkg, tbg, ev, es, &energy, Interaction);
	Basis(c, Nkg, cfg, es);

	itv = 255;
	emin = -16;
	emax =  16;

	fprintf(f, "%22s", "#e");
	for(i=0; i<c.Nb; i++) fprintf(f, "%20s%02d", "d", i+1);
	fprintf(f, "\n");

	for(v=0; v<=itv; v++) {
		memset(dos, 0, sizeof(dos));
		e = emin + (emax - emin) * v / itv;

		for(i=0; i<c.Nb*Nkg; i++) {
			for(j=0; j<c.Nb; j++) dos[j] += GREEN(i) * CSQR(es[c.Nb*i + j]);
		}

		fprintf(f, "%22.16f", e);
		for(i=0; i<c.Nb; i++) fprintf(f, "%22.16f", dos[i] / vol);
		fprintf(f, "\n");
	}

	fclose(f);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fn, t1 - t0);
}
