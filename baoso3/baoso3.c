// baoso3/baoso3.c : calculate BaOsO3 model

#include "../hf3.h"

void CalcBandPath() {
	FILE *f = OpenFile("input/gb.bin", "wb");

	int i = 0, j;
	Vector gb[BAND];

	const int sym_len = GetLen("input/info.txt");
	int sym[sym_len];
	ReadPathInfo(sym_len, sym);

	for(j=0; j<sym[1]; j++) { // G-X
		gb[i].c[0] = M_PI * j / sym[1];
		gb[i].c[1] = 0;
		gb[i].c[2] = 0;
		i++;
	}

	for(j=0; j<sym[2]; j++) { // X-M
		gb[i].c[0] = M_PI;
		gb[i].c[1] = M_PI * j / sym[2];
		gb[i].c[2] = 0;
		i++;
	}

	for(j=0; j<sym[3]; j++) { // M-G
		gb[i].c[0] = M_PI - M_PI * j / sym[3];
		gb[i].c[1] = M_PI - M_PI * j / sym[3];
		gb[i].c[2] = 0;
		i++;
	}

	for(j=0; j<sym[4]; j++) { // G-R
		gb[i].c[0] = M_PI * j / sym[4];
		gb[i].c[1] = M_PI * j / sym[4];
		gb[i].c[2] = M_PI * j / sym[4];
		i++;
	}

	fwrite(gb, sizeof(gb), 1, f);
	fclose(f);
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

void GaussQuadS(double fermi, double *wg, double *ev, lapack_complex_double *es, double *n, double *m, double *e) {
	int i, j;

	memset(n, 0, OBT * sizeof(double));
	memset(m, 0, OBT * sizeof(double));
	*e = 0;

	for(i=0; i<GAUSS3*SINGLE; i++) {
		if(ev[i] < fermi) {
			for(j=0; j<OBT; j++) {
				n[j] += (CSQR(es[SINGLE*i + j]) + CSQR(es[SINGLE*i + j+OBT])) * wg[i/SINGLE];
				m[j] += (CSQR(es[SINGLE*i + j]) - CSQR(es[SINGLE*i + j+OBT])) * wg[i/SINGLE];
				*e   += (CSQR(es[SINGLE*i + j]) + CSQR(es[SINGLE*i + j+OBT])) * wg[i/SINGLE] * ev[i];
			}
		}
	}
}

void GaussQuadD(double fermi, double *wg, double *ev, lapack_complex_double *es, double *n, double *m, double *e) {
	int i, j;

	memset(n, 0, OBT * sizeof(double));
	memset(m, 0, OBT * sizeof(double));
	*e = 0;

	for(i=0; i<GAUSS3*DOUBLE; i++) {
		if(ev[i] < fermi) {
			for(j=0; j<OBT; j++) {
				n[j] += (CSQR(es[DOUBLE*i + j]) + CSQR(es[DOUBLE*i + j+OBT]) + CSQR(es[DOUBLE*i + j+2*OBT]) + CSQR(es[DOUBLE*i + j+3*OBT])) * wg[i/DOUBLE];
				m[j] += (creal(es[DOUBLE*i + j]) * creal(es[DOUBLE*i + j+OBT]) + cimag(es[DOUBLE*i + j]) * cimag(es[DOUBLE*i + j+OBT])) * 4 * wg[i/DOUBLE];
				*e   += (CSQR(es[DOUBLE*i + j]) + CSQR(es[DOUBLE*i + j+OBT]) + CSQR(es[DOUBLE*i + j+2*OBT]) + CSQR(es[DOUBLE*i + j+3*OBT])) * wg[i/DOUBLE] * ev[i];
			}
		}
	}
}
