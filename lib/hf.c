// lib/hf.c : functions for Hartree-Fock approximation

#include "hf.h" 

void FileName_(Solution *s, char *ftype, char *fn) {
	char save[1024];
	sprintf(save, "%s/%s", s->save, ftype);
	if(-access(save, 0)) mkdir(save, 0755);

	if(strstr(ftype, "sol")) {
		sprintf(fn, "%s/%s_N%.1f_U%.1f_v%s.txt",\
			   	ftype, s->N, s->U, s->runtime);
	}
	else {
		sprintf(fn, "%s/%s_N%.1f_U%.1f_n%.16f_m%.16f_e%.16f_gap%.16f_fermi%.16f_dntop%.16f_v%s.txt",\
			   	ftype, s->N, s->U, s->ns, s->ms, s->e, s->gap, s->fermi, s->dntop, s->runtime);
	}
}

void CalcEigen_(Config c, double *ev, lapack_complex_double *es) {
	char jobz = 'V', uplo = 'L';
	double rwork[3*c.Nb-2];
	lapack_int ln = lda = c.Nb, lwork = 2*c.Nb-1, info;
	lapack_complex_double work[lwork];

	//Interaction(c, s, tb);

	LAPACK_zheev(&jobz, &uplo, &ln, es, &lda, ev, work, &lwork, rwork, &info);
	if(info != 0) {
		printf("LAPACK_zheev FAIL\n");
		exit(1);
	}

	//if(ev0[0]      < e->min) e->min = ev0[0];
	//if(ev0[c.Nb-1] > e->max) e->max = ev0[c.Nb-1];
}

void CalcGap_(Config c, Solution *s, double *ev, lapack_complex_double *es) {
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

void ReadConfig(Config *c) {
	char fn[256];
	sprintf(fn, "input/config_%s.txt", c->type);

	FILE *f = fopen(fn, "r");
	int i;
	char buf[1024];

	while(!feof(f)) {
		fgets(buf, sizeof(buf), f);
		if      (strstr(buf, "Ni"))      sscanf(buf, "Ni %d", &c->Ni);
		else if (strstr(buf, "Nc"))      sscanf(buf, "Nc %d", &c->Nc);
		else if (strstr(buf, "Q"))       sscanf(buf, "Q %lf%lf%lf", &c->Q[0], &c->Q[1], &c->Q[2]);
		else if (strstr(buf, "Lattice")) sscanf(buf, "Lattice %s", c->lat);
	}
	fclose(f);

	c->Ns  = c->Ni * c->Nc;
	c->Nb  = c->Ni * c->Nc * 2;
	for(i=0; i<DIM; i++) c->Q[i] *= M_PI;
}

void InitSolution(Config c, Solution *s) {
	int i;

	for(i=0; i<c.Nc; i++) {
		s->n[i] = s->N / c.Nc;
		s->m[i] = M_INIT;
	}
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
				es0[0] = sqrt(0.5) * (cos(cf[c.Ni*i + 1]) - sin(cf[c.Ni*i + 1])*I) * (es[STATE_IDX] + es[STATE_IDX + c.Nc]);
				es0[1] = sqrt(0.5) * (cos(cf[c.Ni*i + 0]) - sin(cf[c.Ni*i + 0])*I) * (es[STATE_IDX] - es[STATE_IDX + c.Nc]);

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

#define OC_IDX (c.Nc*i + j)

void CalcSolution(Config c, Solution *s, LAPACK *lp, void (*System)(), void (*Interaction)(), void (*Basis)()) {
	FILE *f;
	char fn[256];

	FileName_(c, s, "sol", fn);
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

	FileName_(c, s, "band", fn);
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

	FileName_(c, s, "ufw", fn);
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

	FileName_(c, s, "dos", fn);
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
