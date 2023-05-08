// lib/hf.c : functions for Hartree-Fock approximation

#define OBT_IDX   ((i / c.Nb) % Nc)
#define ES_IDX    (c.Nb*i + Nc*(j/Nc) + j)
#define OC_IDX    (Nc*i + j)
#define INTER_N   (0.5 * ((c.U) * n[OBT_IDX] + (c.U - 2*c.J) * n_[OBT_IDX] + (c.U - 3*c.J) * n_[OBT_IDX]))
#define INTER_M   (0.5 * ((c.U) * m[OBT_IDX] + (c.U - 2*c.J) * m_[OBT_IDX] - (c.U - 3*c.J) * m_[OBT_IDX])) 
#define GREEN_IMG (ep / (pow(e - ev[i] + s->fermi, 2) + pow(ep, 2))) 

#include "hf.h" 

void FileName(Config c, Solution *s, char *ftype, char *fn) {
	char path_save[512];
	sprintf(path_save, "%s/%s", c.path_save, ftype);
	if(-access(path_save, 0)) mkdir(path_save, 0755);

	if(strstr(ftype, "sol")) {
		sprintf(fn, "%s/N%.1f_U%.1f.bin",\
			   	path_save, c.N, c.U);
	}
	else if(strstr(ftype, "oc")) {
		sprintf(fn, "%s/N%.1f_U%.1f_%s.txt",\
			   	path_save, c.N, c.U, c.runtime);
	}
	else {
		sprintf(fn, "%s/N%.1f_U%.1f_n%f_m%f_e%f_gap%f_fermi%f_dntop%f_%s.txt",\
			   	path_save, c.N, c.U, s->ns, s->ms, s->e, s->gap, s->fermi, s->dntop, c.runtime);
	}
}

void ReadSolution(Config c, Solution *s) {
	char fsn[256];
	sprintf(fsn, "%s/sol/%s.bin", c.path_save, c.sol);

	FILE *fs = fopen(fsn, "rb");
	fread(s, sizeof(Solution), 1, fs);
	fclose(fs);
}

void CalcGap(Config c, Solution *s, double *ev, lapack_complex_double *es, double *uplow, double *dntop) {
	int i;

	for(i=0; i<c.Nb; i++) {
		if     (ev[i] > s->fermi + FERMI_WIDTH / 2 && ev[i] < *uplow) *uplow = ev[i];
		else if(ev[i] < s->fermi - FERMI_WIDTH / 2 && ev[i] > *dntop) *dntop = ev[i];
	}
}

void CalcE(Config c, Solution *s, double w, double *ev, lapack_complex_double *es, double *e) {
	int i, j;

	for(i=0; i<c.Nb; i++) {
		if(ev[i] < s->fermi) {
			for(j=0; j<c.Nb; j++) *e += ev[i] * CSQR(es[c.Nb*i + j]) * w;
		}
	}
}

void CalcUFW(Config c, double *uf, lapack_complex_double *es, double *ufw) {
	int i, j, k;
	lapack_complex_double p;
	
	for(i=0; i<c.Nb; i++) {
		ufw[i] = 0;
		for(j=0; j<c.Ns; j++) {
			p = 0;
			for(k=0; k<c.Ni; k++) p += sqrt(1.0/c.Ni) * (cos(uf[k]) - sin(uf[k]) * I) * es[ES_IDX + Nc*k];
			ufw[i] += CSQR(p);
		}
	}
}

void InteractionN(Config c, Solution *s, lapack_complex_double *tb) {
	int i, j;
	double n[Nc], m[Nc], n_[Nc], m_[Nc];
	
	memset(n_, 0, sizeof(n_));
	memset(m_, 0, sizeof(m_));

	for(i=0; i<Nc; i++) {
		n[i] = s->n[i] / c.Ni;
		m[i] = s->m[i] / c.Ni;

		for(j=0; j<Nc; j++) {
			if(j != i) {
				n_[i] += s->n[j] / c.Ni;
				m_[i] += s->m[j] / c.Ni;
			}
		}		
	}

	for(i=0; i<c.Nb*c.Nb; i+=c.Nb+1) tb[i] += INTER_N;

	for(i=0;             i<c.Nb*c.Ns; i+=c.Nb+1) tb[i] -= INTER_M;
	for(i=(c.Nb+1)*c.Ns; i<c.Nb*c.Nb; i+=c.Nb+1) tb[i] += INTER_M;
}

void InteractionQ(Config c, Solution *s, lapack_complex_double *tb) {
	int i, j;
	double n[Nc], m[Nc], n_[Nc], m_[Nc];
	
	memset(n_, 0, sizeof(n_));
	memset(m_, 0, sizeof(m_));

	for(i=0; i<Nc; i++) {
		n[i] = s->n[i] / c.Ni;
		m[i] = s->m[i] / c.Ni;

		for(j=0; j<Nc; j++) {
			if(j != i) {
				n_[i] += s->n[j] / c.Ni;
				m_[i] += s->m[j] / c.Ni;
			}
		}		
	}

	for(i=0; i<c.Nb*c.Nb; i+=c.Nb+1) tb[i] += INTER_N;

	for(i=Nc;                 i<c.Nb*Nc;   i+=c.Nb+1) tb[i] -= INTER_M;
	for(i=Nc + (c.Nb+1)*c.Ns; i<c.Nb*c.Nb; i+=c.Nb+1) tb[i] += INTER_M;
}

/*
void InteractionS(Config c, Solution *s, lapack_complex_double *tb) {
	int i, j;
	double n[Nc], m[Nc], n_[Nc], m_[Nc];
	
	memset(n_, 0, sizeof(n_));
	memset(m_, 0, sizeof(m_));

	for(i=0; i<Nc; i++) {
		n[i] = s->n[i] / c.Ni;
		m[i] = s->m[i] / c.Ni;

		for(j=0; j<Nc; j++) {
			if(j != i) {
				n_[i] += s->n[j] / c.Ni;
				m_[i] += s->m[j] / c.Ni;
			}
		}		
	}

	for(i=0; i<c.Nb*c.Nb; i+=c.Nb+1) tb[i] += INTER_N;

	for(i=0;             i<c.Nb*c.Ns; i+=c.Nb+1) tb[i] -= INTER_M * pow(-1, i / (c.Nb*Nc));
	for(i=(c.Nb+1)*c.Ns; i<c.Nb*c.Nb; i+=c.Nb+1) tb[i] += INTER_M * pow(-1, (i + c.Nb*c.Ns) / (c.Nb*Nc)) ;
}
*/

void BasisN(Config c, double *uf, lapack_complex_double *es) {}
void BasisQ(Config c, double *uf, lapack_complex_double *es) {
	int i, j;
	lapack_complex_double es0[c.Ni];

	for(i=0; i<c.Nb; i++) {
		for(j=0; j<c.Ns; j++) {
			es0[0] = sqrt(1.0/c.Ni) * (cos(uf[1]) - sin(uf[1]) * I) * (es[ES_IDX] + es[ES_IDX + Nc]);
			es0[1] = sqrt(1.0/c.Ni) * (cos(uf[0]) - sin(uf[0]) * I) * (es[ES_IDX] - es[ES_IDX + Nc]);

			es[ES_IDX]      = es0[0];
			es[ES_IDX + Nc] = es0[1];
		}
	}
}

void Quadrature(Config c, Solution *s, double w, double *ev, lapack_complex_double *es, double *oc) {
	int i, j;

	for(i=0; i<c.Nb; i++) {
		if(ev[i] < s->fermi) {
			for(j=0; j<c.Nb; j++) oc[j] += CSQR(es[c.Nb*i + j]) * w;
		}
	}
}

void GenSolution(Config c, Solution *s, void (*Symmetry)(), void (*Interaction)(), void (*Basis)()) {
	char fkn[256], fun[256], ftn[256], fon[256], fsn[256];
	sprintf(fkn, "input/%s/tb/kg_Nk%d.txt",     c.strain, c.Nkg);
	sprintf(fun, "input/%s/tb/ufg_Nk%d_%c.txt", c.strain, c.Nkg, c.type[0]);
	sprintf(ftn, "input/%s/tb/tbg_Nk%d_%c.bin", c.strain, c.Nkg, c.type[0]);
	FileName(c, s, "oc",  fon);
	FileName(c, s, "sol", fsn);

	time_t t0 = time(NULL);

	FILE *fk = fopen(fkn, "r"), *fu = fopen(fun, "r"), *ft = fopen(ftn, "rb"), *fo = fopen(fon, "w"), *fs = fopen(fsn, "wb");
	int itr, i, j, is_cvg;
	char buf[1024];
	double tmp, w[c.Nkg], uf[c.Nkg][c.Ni], ev[c.Nkg][c.Nb], e_min, e_max, e = 0;
	double V = pow(2*M_PI, DIM), oc[c.Nb], oc_sum, cvg[Nc][CVG_MAX], avg;
	double uplow = 100, dntop = -100;
	lapack_complex_double tb[c.Nkg][c.Nb*c.Nb], es[c.Nkg][c.Nb*c.Nb];

	// w
	fgets(buf, sizeof(buf), fk); // skip header
	for(i=0; i<c.Nkg; i++) fscanf(fk, "%lf%lf%lf%lf", &tmp, &tmp, &tmp, &w[i]);
	fclose(fk);

	// uf
	fgets(buf, sizeof(buf), fu); // skip header
	for(i=0; i<c.Nkg; i++) {
		for(j=0; j<c.Ni; j++) fscanf(fu, "%lf", &uf[i][j]);
	}
	fclose(fu);

	// tb
	fread(tb, sizeof(tb), 1, ft);
	fclose(ft);

	for(i=0; i<Nc; i++) {
		for(j=0; j<CVG_MAX; j++) cvg[i][j] = 100;
	}

	fprintf(fo, "%22s", "fermi");
	for(i=0; i<c.Nb; i++) fprintf(fo, "%20s%02d", "oc", i+1);
	fprintf(fo, "\n");

	for(itr=0; itr<ITR_MAX; itr++) {
		e_min =  100;
		e_max = -100;

		for(i=0; i<c.Nkg; i++) {
			for(j=0; j<c.Nb*c.Nb; j++) es[i][j] = tb[i][j];

			Interaction(c, s, es[i]);
			CalcEigen(c, ev[i], es[i]);
			Basis(c, uf[i], es[i]);

			if(ev[i][0]      < e_min) e_min = ev[i][0];
			if(ev[i][c.Nb-1] > e_max) e_max = ev[i][c.Nb-1];
		}
		
		s->fermi = 0.5 * (e_min + e_max);

		while(fabs(s->fermi - e_max) > 1e-8) {
			memset(oc, 0, sizeof(oc));
			oc_sum = 0;

			for(i=0; i<c.Nkg; i++) Quadrature(c, s, w[i], ev[i], es[i], oc);

			for(i=0; i<c.Nb; i++) {
				oc[i] /= V;
				oc_sum += oc[i];
			}

			if(fabs(oc_sum - c.N) < 1e-4) break;
			else {
				if(oc_sum < c.N) e_min = s->fermi;
				else             e_max = s->fermi;
				s->fermi = 0.5 * (e_min + e_max);
			}
		}

		for(i=0; i<Nc; i++) {
			s->n[i] = 0;
			s->m[i] = 0;
		}
		s->ns = 0;
		s->ms = 0;

		for(i=0; i<c.Ni; i++) {
			for(j=0; j<Nc; j++) {
				s->n[j] +=  oc[OC_IDX] + oc[OC_IDX + c.Ns];
				s->m[j] += (oc[OC_IDX] - oc[OC_IDX + c.Ns]) * pow(-1, i);
			}
		}
		for(i=0; i<Nc; i++) {
			s->ns += s->n[i];
			s->ms += s->m[i];
		}

		fprintf(fo, "%22.16f", s->fermi);
		for(i=0; i<c.Nb; i++) fprintf(fo, "%22.16f", oc[i]);
		fprintf(fo, "\n");
			
		is_cvg = 0;
		for(i=0; i<Nc; i++) {
			cvg[i][itr % CVG_MAX] = s->m[i];

			avg = 0;
			for(j=0; j<CVG_MAX; j++) avg += cvg[i][j];
			avg /= CVG_MAX;

			if(fabs(avg - s->m[i]) < 1e-6) is_cvg++;
		}
		if(is_cvg == Nc) break;
		Symmetry(s->n, s->m);
	}
	fclose(fo);

	for(i=0; i<c.Nkg; i++) CalcGap(c, s, ev[i], es[i], &uplow, &dntop);
	s->gap   = uplow - dntop;
	s->dntop = dntop;

	for(i=0; i<c.Nkg; i++) CalcE(c, s, w[i], ev[i], es[i], &e);
	s->e = e / V;

	fwrite(s, sizeof(Solution), 1, fs);
	fclose(fs);

	time_t t1 = time(NULL);
	printf("%s(%s, %s) : %lds\n", __func__, fon, fsn, t1 - t0);

	GenSolDOS(c, s, EP, (double*)ev, (lapack_complex_double*)es);
	GenBand(c, s, Interaction, Basis);
}

void GenSolDOS(Config c, Solution *s, double ep, double *ev, lapack_complex_double *es) {
	char ftype[32], fdn[256];
	sprintf(ftype, "dos_ep%.2f", ep);
	FileName(c, s, ftype, fdn);

	time_t t0 = time(NULL);

	FILE *fd = fopen(fdn, "w");
	int itv, itv_max, i, j;
	double e, e_min, e_max, dos[c.Nb], V = pow(2*M_PI, 3)*pow(2*M_PI, 3)*M_PI; 
	
	itv_max = 256;
	e_min = -8;
	e_max =  8;

	fprintf(fd, "%22s", "e");
	for(i=0; i<c.Nb; i++) fprintf(fd, "%20s%02d", "dos", i+1);
	fprintf(fd, "\n");

	for(itv=0; itv<=itv_max; itv++) {
		memset(dos, 0, sizeof(dos));
		e = e_min + (e_max - e_min) * itv / itv_max;

		for(i=0; i<c.Nkg*c.Nb; i++) {
			for(j=0; j<c.Nb; j++) dos[j] += GREEN_IMG * CSQR(es[c.Nb*i + j]);
		}

		fprintf(fd, "%22.16f", e);
		for(i=0; i<c.Nb; i++) fprintf(fd, "%22.16f", dos[i] / V);
		fprintf(fd, "\n");
	}

	fclose(fd);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fdn, t1 - t0);
}

void GenBand(Config c, Solution *s, void (*Interaction)(), void (*Basis)()) {
	char fun[256], ftn[256], ftype[32], fbn[256];
	sprintf(fun, "input/%s/tb/ufb_Nk%d_%c.txt", c.strain, c.Nkb, c.type[0]);
	sprintf(ftn, "input/%s/tb/tbb_Nk%d_%c.bin", c.strain, c.Nkb, c.type[0]);
	sprintf(ftype, "band_Nk%d", c.Nkb);
	FileName(c, s, ftype, fbn);

	time_t t0 = time(NULL);

	FILE *fu = fopen(fun, "r"), *ft = fopen(ftn, "rb"), *fb = fopen(fbn, "w");
	int i, j;
	char buf[1024];
	double uf[c.Nkb][c.Ni], ev[c.Nkb][c.Nb], ufw[c.Nb];
	lapack_complex_double tb[c.Nkb][c.Nb*c.Nb], es[c.Nkb][c.Nb*c.Nb];

	// uf
	fgets(buf, sizeof(buf), fu); // skip header	
	for(i=0; i<c.Nkb; i++) {
		for(j=0; j<c.Ni; j++) fscanf(fu, "%lf", &uf[i][j]);
	}
	fclose(fu);

	// tb
	fread(tb, sizeof(tb), 1, ft);
	fclose(ft);

	for(i=0; i<c.Nb; i++) fprintf(fb, "%20s%02d", "e", i+1);
	for(i=0; i<c.Nb; i++) fprintf(fb, "%20s%02d", "w", i+1);
	fprintf(fb, "\n");

	for(i=0; i<c.Nkb; i++) {
		for(j=0; j<c.Nb*c.Nb; j++) es[i][j] = tb[i][j];

		Interaction(c, s, es[i]);
		CalcEigen(c, ev[i], es[i]);
		Basis(c, uf[i], es[i]);
		CalcUFW(c, uf[i], es[i], ufw);

		for(j=0; j<c.Nb; j++) fprintf(fb, "%22.16f", ev[i][j] - s->fermi);
		for(j=0; j<c.Nb; j++) fprintf(fb, "%22.16f", ufw[j]);
		fprintf(fb, "\n");
	}
	fclose(fb);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fbn, t1 - t0);
}

void GenDOS(Config c, Solution *s, double ep, void (*Interaction)(), void (*Basis)()) {
	char fkn[256], fun[256], ftn[256];
	sprintf(fkn, "input/%s/tb/kg_Nk%d.txt",     c.strain, c.Nkg);
	sprintf(fun, "input/%s/tb/ufg_Nk%d_%c.txt", c.strain, c.Nkg, c.type[0]);
	sprintf(ftn, "input/%s/tb/tbg_Nk%d_%c.bin", c.strain, c.Nkg, c.type[0]);

	FILE *fk = fopen(fkn, "r"), *fu = fopen(fun, "r"), *ft = fopen(ftn, "rb");
	int i, j;
	char buf[1024];
	double tmp, w[c.Nkg], uf[c.Nkg][c.Ni], ev[c.Nkg][c.Nb];
	lapack_complex_double tb[c.Nkg][c.Nb*c.Nb], es[c.Nkg][c.Nb*c.Nb];

	// w
	fgets(buf, sizeof(buf), fk); // skip header
	for(i=0; i<c.Nkg; i++) fscanf(fk, "%lf%lf%lf%lf", &tmp, &tmp, &tmp, &w[i]);
	fclose(fk);

	// uf
	fgets(buf, sizeof(buf), fu); // skip header
	for(i=0; i<c.Nkg; i++) {
		for(j=0; j<c.Ni; j++) fscanf(fu, "%lf", &uf[i][j]);
	}
	fclose(fu);

	// tb
	fread(tb, sizeof(tb), 1, ft);
	fclose(ft);

	for(i=0; i<c.Nkg; i++) {
		for(j=0; j<c.Nb*c.Nb; j++) es[i][j] = tb[i][j];

		Interaction(c, s, es[i]);
		CalcEigen(c, ev[i], es[i]);
		Basis(c, uf[i], es[i]);
	}

	GenSolDOS(c, s, ep, (double*)ev, (lapack_complex_double*)es);
}
