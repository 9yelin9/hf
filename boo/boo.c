// boo/boo.c : calculate BaOsO3 model

#define OBT_IDX ((i / s->basis) % OBT)
#define INTER_N (0.5*s->U * n[OBT_IDX] + (s->U - 2.5*s->J) * n_[OBT_IDX])
#define INTER_M (0.5*s->U * m[OBT_IDX] + (s->U - 2.5*s->J) * m_[OBT_IDX])

#include "../hf3.h"
#include "boo.h"

void InteractionF(Solution *s, lapack_complex_double *v_tmp) {
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

	for(i=0; i<H(BASIS1); i+=(BASIS1+1)) v_tmp[i] += INTER_N;
	for(i=0; i<BASIS1*OBT; i+=(BASIS1+1)) v_tmp[i] -= INTER_M;
	for(i=(BASIS1+1)*OBT; i<BASIS1*(2*OBT); i+=(BASIS1+1)) v_tmp[i] += INTER_M;
}

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

	for(i=0; i<H(BASIS2); i+=(BASIS2+1)) v_tmp[i] += INTER_N;
	for(i=OBT; i<BASIS2*OBT; i+=(BASIS2+1)) v_tmp[i] -= INTER_M;
	for(i=(BASIS2+1)*OBT - OBT; i<BASIS2*(2*OBT); i+=(BASIS2+1)) v_tmp[i] -= INTER_M;
	for(i=(BASIS2+1)*(2*OBT) + OBT; i<BASIS2*(3*OBT); i+=(BASIS2+1)) v_tmp[i] += INTER_M;
	for(i=(BASIS2+1)*(3*OBT) - OBT; i<BASIS2*(4*OBT); i+=(BASIS2+1)) v_tmp[i] += INTER_M;
}

void InteractionSubA(Solution *s, lapack_complex_double *v_tmp) {
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

	for(i=0; i<H(BASIS2); i+=(BASIS2+1)) v_tmp[i] += INTER_N;
	for(i=0; i<BASIS2*OBT; i+=(BASIS2+1)) v_tmp[i] -= INTER_M;
	for(i=(BASIS2+1)*OBT; i<BASIS2*(2*OBT); i+=(BASIS2+1)) v_tmp[i] += INTER_M;
	for(i=(BASIS2+1)*(2*OBT); i<BASIS2*(3*OBT); i+=(BASIS2+1)) v_tmp[i] += INTER_M;
	for(i=(BASIS2+1)*(3*OBT); i<BASIS2*(4*OBT); i+=(BASIS2+1)) v_tmp[i] -= INTER_M;
}

void OccupationF(double fermi, double *w, lapack_complex_double *v, double *n, double *m, double *e) {
	int i, j;

	memset(n, 0, OBT*sizeof(double));
	memset(m, 0, OBT*sizeof(double));
	*e = 0;

	for(i=0; i<K3*BASIS1; i++) {
		if(w[i] < fermi) {
			for(j=0; j<OBT; j++) {
				n[j] += COMPLEX2(v[BASIS1*i + j]) + COMPLEX2(v[BASIS1*i + j+OBT]);
				m[j] += COMPLEX2(v[BASIS1*i + j]) - COMPLEX2(v[BASIS1*i + j+OBT]);
				*e += w[i] * (COMPLEX2(v[BASIS1*i + j]) + COMPLEX2(v[BASIS1*i + j+OBT]));
			}
		}
	}
}

void OccupationA(double fermi, double *w, lapack_complex_double *v, double *n, double *m, double *e) {
	int i, j;
	lapack_complex_double k1, k2;

	memset(n, 0, OBT*sizeof(double));
	memset(m, 0, OBT*sizeof(double));
	*e = 0;

	for(i=0; i<K3*BASIS2; i++) {
		if(w[i] < fermi) {
			for(j=0; j<OBT; j++) {
				k1 = creal(v[BASIS2*i + j]) - cimag(v[BASIS2*i + j]) * I;
				k2 = creal(v[BASIS2*i + j+2*OBT]) - cimag(v[BASIS2*i + j+2*OBT]) * I;

				n[j] += COMPLEX2(v[BASIS2*i + j]) + COMPLEX2(v[BASIS2*i + j+OBT]) + COMPLEX2(v[BASIS2*i + j+2*OBT]) + COMPLEX2(v[BASIS2*i + j+3*OBT]);
				m[j] += 2 * (k1 * v[BASIS2*i + j+OBT] - k2 * v[BASIS2*i + j+3*OBT]);
				*e += w[i] * (COMPLEX2(v[BASIS2*i + j]) + COMPLEX2(v[BASIS2*i + j+OBT]) + COMPLEX2(v[BASIS2*i + j+2*OBT]) + COMPLEX2(v[BASIS2*i + j+3*OBT]));
			}
		}
	}
}

void OccupationSubA(double fermi, double *w, lapack_complex_double *v, double *n, double *m, double *e) {
	int i, j;

	memset(n, 0, OBT*sizeof(double));
	memset(m, 0, OBT*sizeof(double));
	*e = 0;

	for(i=0; i<K3*BASIS2; i++) {
		if(w[i] < fermi) {
			for(j=0; j<OBT; j++) {
				n[j] += COMPLEX2(v[BASIS2*i + j]) + COMPLEX2(v[BASIS2*i + j+OBT]) + COMPLEX2(v[BASIS2*i + j+2*OBT]) + COMPLEX2(v[BASIS2*i + j+3*OBT]);
				m[j] += 2 * (COMPLEX2(v[BASIS2*i + j]) - COMPLEX2(v[BASIS2*i + j+2*OBT]));
				*e += w[i] * (COMPLEX2(v[BASIS2*i + j]) + COMPLEX2(v[BASIS2*i + j+OBT]) + COMPLEX2(v[BASIS2*i + j+2*OBT]) + COMPLEX2(v[BASIS2*i + j+3*OBT]));
			}
		}
	}
}

int main(int argc, char *argv[]) {
	if(argc != 7) {
		printf("Usage : %s <type> <J/U> <SOC> <N> <U> <is_unfold>\n", argv[0]);
		exit(1);
	}

	int is_unfold = atoi(argv[6]);

	Solution s = {
		.type = argv[1],
		.JU = atof(argv[2]),
		.SOC = atof(argv[3]),
		.N = atof(argv[4]),
		.U = atof(argv[5]),
		.J = s.JU * s.U,

		.n = {s.N/3, s.N/3, s.N/3},
		.m = {0.1, 0.1, 0.1},
		.n_total = 100,
		.m_total = 100,
		.fermi = 100,
		.e = 100
	};

	time_t t = time(NULL);
	struct tm *tm = localtime(&t);
	sprintf(s.runtime, "%d%d%d%d", tm->tm_mon+1, tm->tm_mday, tm->tm_hour, tm->tm_min); 

	s.basis = strstr(s.type, "f") ? 6 : 12,
	s.tbk = (lapack_complex_double*)malloc(HK(s.basis) * sizeof(lapack_complex_double));
	s.tbb = (lapack_complex_double*)malloc(HB(s.basis) * sizeof(lapack_complex_double));

	ReadTB(&s);
	CalcSolution(&s, is_unfold);

	MakeBand(&s, is_unfold);
	MakeDos(&s, is_unfold);

	free(s.tbk);
	free(s.tbb);

	return 0;
}
