// baoso3/baoso3.c : calculate BaOsO3 model

#define OBT_IDX ((i / s->basis) % OBT)
#define INTER_N (0.5*s->U * n[OBT_IDX] + (s->U - 2.5*s->J) * n_[OBT_IDX])
#define INTER_M (0.5*s->U * m[OBT_IDX] + (s->U - 2.5*s->J) * m_[OBT_IDX])

#include "../hf3.h"

void CalcVB() {
	FILE *f;
	char fs[64];

	sprintf(fs, "input/vb.bin");

	if((f = fopen(fs, "wb")) == NULL) {
		printf("%s fopen FAIL\n", fs);
		exit(1);
	}

	int i = 0, j, *sym = ReadPath();
	Vector vb[BAND];

	for(j=0; j<sym[1]; j++) { // G-X
		vb[i].k[0] = M_PI * j / sym[1];
		vb[i].k[1] = 0;
		vb[i].k[2] = 0;
		i++;
	}

	for(j=0; j<sym[2]; j++) { // X-M
		vb[i].k[0] = M_PI;
		vb[i].k[1] = M_PI * j / sym[2];
		vb[i].k[2] = 0;
		i++;
	}

	for(j=0; j<sym[3]; j++) { // M-G
		vb[i].k[0] = M_PI - M_PI * j / sym[3];
		vb[i].k[1] = M_PI - M_PI * j / sym[3];
		vb[i].k[2] = 0;
		i++;
	}

	for(j=0; j<sym[4]; j++) { // G-R
		vb[i].k[0] = M_PI * j / sym[4];
		vb[i].k[1] = M_PI * j / sym[4];
		vb[i].k[2] = M_PI * j / sym[4];
		i++;
	}

	free(sym);

	fwrite(vb, sizeof(vb), 1, f); 
	fclose(f);
}

void FourierS(int l_len, int h_num, Vector v, Vector vq, Lattice *l, lapack_complex_double *h) {
	int i, j, k = 0;
	double dot;
	lapack_complex_double tmp[SINGLE][SINGLE], t, e;

	memset(tmp, 0, sizeof(tmp));

	for(i=0; i<l_len; i++) {
		dot = 0;
		for(j=0; j<3; j++) dot += l[i].r[j] * v.k[j];

		t = l[i].tre + l[i].tim * I;
		e = cos(dot) + sin(dot) * I;

		tmp[l[i].obt1-1][l[i].obt2-1] += t * e;
	}

	for(i=0; i<OBT; i++) {
		for(j=0; j<OBT; j++) {
			tmp[i + OBT][j + OBT] = tmp[i][j];
		}
	}

	for(i=0; i<SINGLE; i++) {
		for(j=0; j<SINGLE; j++) {
			h[SINGLE*SINGLE*h_num + k] = tmp[i][j];
			k++;
		}
	}
}

void FourierD(int l_len, int h_num, Vector v, Vector vq, Lattice *l, lapack_complex_double *h) {
	int i, j, k = 0;
	double dot, dot_sign, sign;
	lapack_complex_double tmp[DOUBLE][DOUBLE], t, e;

	memset(tmp, 0, sizeof(tmp));

	for(i=0; i<l_len; i++) {
		dot = dot_sign = 0;
		for(j=0; j<3; j++) {
			dot += l[i].r[j] * v.k[j];
			dot_sign += l[i].r[j] * vq.k[j];
		}

		t = l[i].tre + l[i].tim * I;
		e = cos(dot) + sin(dot) * I;
		sign = cos(dot_sign);

 		// change basis
		tmp[l[i].obt1-1      ][l[i].obt2-1      ] += 0.5 * t * e * (1 + sign);
		tmp[l[i].obt1-1 + OBT][l[i].obt2-1 + OBT] += 0.5 * t * e * (1 + sign);

		tmp[l[i].obt1-1      ][l[i].obt2-1 + OBT] += 0.5 * t * e * (1 - sign);
		tmp[l[i].obt1-1 + OBT][l[i].obt2-1      ] += 0.5 * t * e * (1 - sign);
	}

	for(i=0; i<2*OBT; i++) {
		for(j=0; j<2*OBT; j++) {
			tmp[i + 2*OBT][j + 2*OBT] = tmp[i][j];
		}
	}

	for(i=0; i<DOUBLE; i++) {
		for(j=0; j<DOUBLE; j++) {
			h[DOUBLE*DOUBLE*h_num + k] = tmp[i][j];
			k++;
		}
	}
}

void InteractionS(Solution *s, lapack_complex_double *v_tmp) {
	int i, j;
	double n[OBT], n_[OBT], m[OBT], m_[OBT];
	
	memset(n_, 0, sizeof(n_));
	memset(m_, 0, sizeof(m_));

	for(i=0; i<OBT; i++) {
		n[i] = s->n[i];
		m[i] = s->m[i];

		for(j=0; j<OBT; j++) {
			if(j != i) {
				n_[i] += s->n[j];
				m_[i] += s->m[j];
			}
		}		
	}

	for(i=0;              i<SINGLE*SINGLE;  i+=(SINGLE+1)) v_tmp[i] += INTER_N;
	for(i=0;              i<SINGLE*OBT;     i+=(SINGLE+1)) v_tmp[i] -= INTER_M;
	for(i=(SINGLE+1)*OBT; i<SINGLE*(2*OBT); i+=(SINGLE+1)) v_tmp[i] += INTER_M;
}
/*
void InteractionA(Solution *s, lapack_complex_double *v_tmp) {
	int i, j;
	double n[OBT], n_[OBT], m[OBT], m_[OBT];
	
	memset(n_, 0, sizeof(n_));
	memset(m_, 0, sizeof(m_));

	for(i=0; i<OBT; i++) {
		n[i] = s->n[i]/2;
		m[i] = s->m[i]/2;

		for(j=0; j<OBT; j++) {
			if(j != i) {
				n_[i] += s->n[j]/2;
				m_[i] += s->m[j]/2;
			}
		}		
	}

	for(i=0;                        i<DOUBLE*DOUBLE;  i+=(DOUBLE+1)) v_tmp[i] += INTER_N;
	for(i=OBT;                      i<DOUBLE*OBT;     i+=(DOUBLE+1)) v_tmp[i] -= INTER_M;
	for(i=(DOUBLE+1)*OBT - OBT;     i<DOUBLE*(2*OBT); i+=(DOUBLE+1)) v_tmp[i] -= INTER_M;
	for(i=(DOUBLE+1)*(2*OBT) + OBT; i<DOUBLE*(3*OBT); i+=(DOUBLE+1)) v_tmp[i] += INTER_M;
	for(i=(DOUBLE+1)*(3*OBT) - OBT; i<DOUBLE*(4*OBT); i+=(DOUBLE+1)) v_tmp[i] += INTER_M;
}
*/
void InteractionD(Solution *s, lapack_complex_double *v_tmp) {
	int i, j;
	double n[OBT], n_[OBT], m[OBT], m_[OBT];
	
	memset(n_, 0, sizeof(n_));
	memset(m_, 0, sizeof(m_));

	for(i=0; i<OBT; i++) {
		n[i] = s->n[i]/2;
		m[i] = s->m[i]/2;

		for(j=0; j<OBT; j++) {
			if(j != i) {
				n_[i] += s->n[j]/2;
				m_[i] += s->m[j]/2;
			}
		}		
	}

	for(i=0;                  i<DOUBLE*DOUBLE;  i+=(DOUBLE+1)) v_tmp[i] += INTER_N;
	for(i=0;                  i<DOUBLE*OBT;     i+=(DOUBLE+1)) v_tmp[i] -= INTER_M;
	for(i=(DOUBLE+1)*OBT;     i<DOUBLE*(2*OBT); i+=(DOUBLE+1)) v_tmp[i] += INTER_M;
	for(i=(DOUBLE+1)*(2*OBT); i<DOUBLE*(3*OBT); i+=(DOUBLE+1)) v_tmp[i] += INTER_M;
	for(i=(DOUBLE+1)*(3*OBT); i<DOUBLE*(4*OBT); i+=(DOUBLE+1)) v_tmp[i] -= INTER_M;
}

void OccupationS(double fermi, double *w, lapack_complex_double *v, double *n, double *m, double *e) {
	int i, j;

	memset(n, 0, OBT * sizeof(double));
	memset(m, 0, OBT * sizeof(double));
	*e = 0;

	for(i=0; i<K3*SINGLE; i++) {
		if(w[i] < fermi) {
			for(j=0; j<OBT; j++) {
				n[j] += CSQR(v[SINGLE*i + j]) + CSQR(v[SINGLE*i + j+OBT]);
				m[j] += CSQR(v[SINGLE*i + j]) - CSQR(v[SINGLE*i + j+OBT]);
				*e   +=(CSQR(v[SINGLE*i + j]) + CSQR(v[SINGLE*i + j+OBT])) * w[i];
			}
		}
	}
}

void OccupationD(double fermi, double *w, lapack_complex_double *v, double *n, double *m, double *e) {
	int i, j;

	memset(n, 0, OBT * sizeof(double));
	memset(m, 0, OBT * sizeof(double));
	*e = 0;

	for(i=0; i<K3*DOUBLE; i++) {
		if(w[i] < fermi) {
			for(j=0; j<OBT; j++) {
				n[j] += CSQR(v[DOUBLE*i + j]) + CSQR(v[DOUBLE*i + j+OBT]) + CSQR(v[DOUBLE*i + j+2*OBT]) + CSQR(v[DOUBLE*i + j+3*OBT]);
				m[j] +=(CSQR(v[DOUBLE*i + j]) - CSQR(v[DOUBLE*i + j+2*OBT])) * 2;
				*e   +=(CSQR(v[DOUBLE*i + j]) + CSQR(v[DOUBLE*i + j+OBT]) + CSQR(v[DOUBLE*i + j+2*OBT]) + CSQR(v[DOUBLE*i + j+3*OBT])) * w[i];
			}
		}
	}
}
