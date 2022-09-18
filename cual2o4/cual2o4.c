// cual2o4/cual2o4.c : calculate CuAl2O4 model

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

	for(j=0; j<sym[0]; j++) { // L-G
		vb[i].k[0] = M_PI - M_PI * j / sym[0];
		vb[i].k[1] = M_PI - M_PI * j / sym[0];
		vb[i].k[2] = M_PI - M_PI * j / sym[0];
		i++;
	}

	for(j=0; j<sym[1]; j++) { // G-X
		vb[i].k[0] = M_PI * j / sym[1];
		vb[i].k[1] = 0;
		vb[i].k[2] = M_PI * j / sym[1];
		i++;
	}

	for(j=0; j<sym[2]; j++) { // X-W
		vb[i].k[0] = M_PI;
		vb[i].k[1] = (M_PI/2) * j / sym[2];
		vb[i].k[2] = M_PI + (M_PI/2) * j / sym[2];
		i++;
	}

	for(j=0; j<sym[3]; j++) { // W-L
		vb[i].k[0] = M_PI;
		vb[i].k[1] = (M_PI/2) + (M_PI/2) * j / sym[3];
		vb[i].k[2] = (3*M_PI/2) - (M_PI/2) * j / sym[3];
		i++;
	}

	for(j=0; j<sym[4]; j++) { // L-K
		vb[i].k[0] = M_PI - (M_PI/4) * j / sym[4];
		vb[i].k[1] = M_PI + (M_PI/2) * j / sym[4];
		vb[i].k[2] = M_PI - (M_PI/4) * j / sym[4];
		i++;
	}

	for(j=0; j<sym[5]; j++) { // K-G
		vb[i].k[0] = (3*M_PI/4) - (3*M_PI/4) * j / sym[5];
		vb[i].k[1] = (3*M_PI/2) - (3*M_PI/2) * j / sym[5];
		vb[i].k[2] = (3*M_PI/4) - (3*M_PI/4) * j / sym[5];
		i++;
	}

	free(sym);

	fwrite(vb, sizeof(vb), 1, f); 
	fclose(f);
}

void FourierS(int l_len, int h_num, Vector v, Vector vq, Lattice *l, lapack_complex_double *h) {
	int i, j, k = 0;
	double dot;
	lapack_complex_double block[SINGLE][SINGLE], t, e;

	memset(block, 0, sizeof(block));

	for(i=0; i<l_len; i++) {
		dot = 0;
		for(j=0; j<3; j++) dot += l[i].r[j] * v.k[j];

		t = l[i].tre + l[i].tim * I;
		e = cos(dot) + sin(dot) * I;

		block[l[i].obt1-1][l[i].obt2-1] += t * e;
	}

	for(i=0; i<OBT; i++) {
		for(j=0; j<OBT; j++) {
			block[i + OBT][j + OBT] = block[i][j];
		}
	}

	for(i=0; i<SINGLE; i++) {
		for(j=0; j<SINGLE; j++) {
			h[SINGLE*SINGLE*h_num + k] = block[i][j];
			k++;
		}
	}
}

void FourierD(int l_len, int h_num, Vector v, Vector vq, Lattice *l, lapack_complex_double *h) {
	int i, j, k = 0;
	double dot_k, dot_q;
	lapack_complex_double block[DOUBLE][DOUBLE], t, e_k, e_q;

	memset(block, 0, sizeof(block));

	for(i=0; i<l_len; i++) {
		dot_k = dot_q = 0;
		for(j=0; j<3; j++) {
			dot_k += l[i].r[j] * (v.k[j]);
			dot_q += l[i].r[j] * (v.k[j] + vq.k[j]);
		}

		t = l[i].tre + l[i].tim * I;
		e_k = cos(dot_k) + sin(dot_k) * I;
		e_q = cos(dot_q) + sin(dot_q) * I;

		block[l[i].obt1-1][l[i].obt2-1] += t * e_k;
		block[l[i].obt1-1 + OBT][l[i].obt2-1 + OBT] += t * e_q;
	}

	for(i=0; i<2*OBT; i++) {
		for(j=0; j<2*OBT; j++) {
			block[i + 2*OBT][j + 2*OBT] = block[i][j];
		}
	}

	for(i=0; i<DOUBLE; i++) {
		for(j=0; j<DOUBLE; j++) {
			h[DOUBLE*DOUBLE*h_num + k] = block[i][j];
			k++;
		}
	}
}

void InteractionS(Solution *s, lapack_complex_double *v_block) {
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

	for(i=0;              i<SINGLE*SINGLE;  i+=(SINGLE+1)) v_block[i] += INTER_N;
	for(i=0;              i<SINGLE*OBT;     i+=(SINGLE+1)) v_block[i] -= INTER_M;
	for(i=(SINGLE+1)*OBT; i<SINGLE*(2*OBT); i+=(SINGLE+1)) v_block[i] += INTER_M;
}

void InteractionD(Solution *s, lapack_complex_double *v_block) {
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

	for(i=0;                        i<DOUBLE*DOUBLE;  i+=(DOUBLE+1)) v_block[i] += INTER_N;
	for(i=OBT;                      i<DOUBLE*OBT;     i+=(DOUBLE+1)) v_block[i] -= INTER_M;
	for(i=(DOUBLE+1)*OBT - OBT;     i<DOUBLE*(2*OBT); i+=(DOUBLE+1)) v_block[i] -= INTER_M;
	for(i=(DOUBLE+1)*(2*OBT) + OBT; i<DOUBLE*(3*OBT); i+=(DOUBLE+1)) v_block[i] += INTER_M;
	for(i=(DOUBLE+1)*(3*OBT) - OBT; i<DOUBLE*(4*OBT); i+=(DOUBLE+1)) v_block[i] += INTER_M;

	/*
	for(i=0;                  i<DOUBLE*DOUBLE;  i+=(DOUBLE+1)) v_block[i] += INTER_N;
	for(i=0;                  i<DOUBLE*OBT;     i+=(DOUBLE+1)) v_block[i] -= INTER_M;
	for(i=(DOUBLE+1)*OBT;     i<DOUBLE*(2*OBT); i+=(DOUBLE+1)) v_block[i] += INTER_M;
	for(i=(DOUBLE+1)*(2*OBT); i<DOUBLE*(3*OBT); i+=(DOUBLE+1)) v_block[i] += INTER_M;
	for(i=(DOUBLE+1)*(3*OBT); i<DOUBLE*(4*OBT); i+=(DOUBLE+1)) v_block[i] -= INTER_M;
	*/
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
				m[j] +=((creal(v[DOUBLE*i + j]) - cimag(v[DOUBLE*i + j]) * I) * v[DOUBLE*i + j+OBT]) * 4;
				//m[j] +=(CSQR(v[DOUBLE*i + j]) - CSQR(v[DOUBLE*i + j+2*OBT])) * 2;
				*e   +=(CSQR(v[DOUBLE*i + j]) + CSQR(v[DOUBLE*i + j+OBT]) + CSQR(v[DOUBLE*i + j+2*OBT]) + CSQR(v[DOUBLE*i + j+3*OBT])) * w[i];
			}
		}
	}
}
