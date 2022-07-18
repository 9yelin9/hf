// boo/boo.c : calculate BaOsO3 model

#include "../hf3.h"
#include "boo.h"

void InteractionF(Solution *s, lapack_complex_double *v_tmp) {
	int i, j;
	double n_diff[OBT], m_diff[OBT];
	
	memset(n_diff, 0, sizeof(n_diff));
	memset(m_diff, 0, sizeof(m_diff));

	for(i=0; i<OBT; i++) {
		for(j=0; j<OBT; j++) {
			if(j != i) {
				n_diff[i] += s->n[j];
				m_diff[i] += s->m[j];
			}
		}		
	}

	for(i=0; i<H(s->basis); i+=(s->basis+1)) {
		v_tmp[i] += N_TERM(i); 
	}
	for(i=0; i<OBT*(s->basis+1); i+=(s->basis+1)) {
		v_tmp[i] -= M_TERM(i); 
	}
	for(i=OBT*(s->basis+1); i<2*OBT*(s->basis+1); i+=(s->basis+1)) {
		v_tmp[i] += M_TERM(i); 
	}	
}

void InteractionA(Solution *s, lapack_complex_double *v_tmp) {
	int i, j;
	double n_diff[OBT], m_diff[OBT];
	
	memset(n_diff, 0, sizeof(n_diff));
	memset(m_diff, 0, sizeof(m_diff));

	for(i=0; i<OBT; i++) {
		for(j=0; j<OBT; j++) {
			if(j != i) {
				n_diff[i] += s->n[j];
				m_diff[i] += s->m[j];
			}
		}		
	}

	for(i=0; i<H(s->basis); i+=(s->basis+1)) {
		v_tmp[i] += (s->U * s->n[OBT_IDX(i)] + (s->U - 2.5*s->J) * n_diff[OBT_IDX(i)]); 
	}
	for(i=OBT; i<OBT*(s->basis+1); i+=(s->basis+1)) {
		v_tmp[i] -= M_TERM(i); 
	}
	for(i=OBT + 2*OBT*(s->basis+1); i<4*OBT*(s->basis+1); i+=(s->basis+1)) {
		v_tmp[i] += M_TERM(i); 
	}	
}

void OccupationF(int basis, double fermi, double *w, lapack_complex_double *v, double *n, double *m, double *e) {
	int i, j;

	memset(n, 0, OBT);
	memset(m, 0, OBT);
	*e = 0;

	for(i=0; i<K3*basis; i++) {
		if(w[i] < fermi) {
			for(j=0; j<OBT; j++) {
				n[j] += COMPLEX2(v[basis*i + j]) + COMPLEX2(v[basis*i + j+OBT]);
				m[j] += COMPLEX2(v[basis*i + j]) - COMPLEX2(v[basis*i + j+OBT]);
				*e += w[i] * (COMPLEX2(v[basis*i + j]) + COMPLEX2(v[basis*i + j+OBT]));
			}
		}
	}
}

void OccupationA(int basis, double fermi, double *w, lapack_complex_double *v, double *n, double *m, double *e) {
	int i, j;

	memset(n, 0, OBT);
	memset(m, 0, OBT);
	*e = 0;

	for(i=0; i<K3*basis; i++) {
		if(w[i] < fermi) {
			for(j=0; j<OBT; j++) {
				n[j] += COMPLEX2(v[basis*i + j]) + COMPLEX2(v[basis*i + j+OBT]) + COMPLEX2(v[basis*i + j+2*OBT]) + COMPLEX2(v[basis*i + j+3*OBT]);
				m[j] += 2 * (v[basis*i + j] * v[basis*i + j+OBT] - v[basis*i + j+2*OBT] * v[basis*i + j+3*OBT]);
				*e += w[i] * (COMPLEX2(v[basis*i + j]) + COMPLEX2(v[basis*i + j+OBT]) + COMPLEX2(v[basis*i + j+2*OBT]) + COMPLEX2(v[basis*i + j+3*OBT]));
			}
		}
	}
}

int main(int argc, char *argv[]) {
	if(argc != 6) {
		printf("Usage : %s <type> <J/U> <SOC> <N> <U>\n", argv[0]);
		exit(1);
	}

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

	ReadTB(s.type, s.basis, s.tbk, s.tbb);
	CalcSolution(&s, 0);

	MakeBand(&s);
	MakeDos(&s);

	free(s.tbk);
	free(s.tbb);

	return 0;
}
