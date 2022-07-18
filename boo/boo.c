// boo/boo.c : calculate BaOsO3 model

#define OBT_IDX ((i / s->basis) % OBT)

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

	for(i=0; i<H(BASIS1); i+=(BASIS1+1)) {
		v_tmp[i] += 0.5*s->U * s->n[OBT_IDX] + (s->U - 2.5*s->J) * n_diff[OBT_IDX]; 
	}
	for(i=0; i<OBT*(BASIS1+1); i+=(BASIS1+1)) {
		v_tmp[i] -= 0.5*s->U * s->m[OBT_IDX] - 0.5*s->J * m_diff[OBT_IDX]; 
	}
	for(i=OBT*(BASIS1+1); i<2*OBT*(BASIS1+1); i+=(BASIS1+1)) {
		v_tmp[i] += 0.5*s->U * s->m[OBT_IDX] - 0.5*s->J * m_diff[OBT_IDX]; 
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
	//for(i=0; i<H(BASIS2); i++) v_tmp[i] = 0;

	for(i=0; i<H(BASIS2); i+=(BASIS2+1)) {
		v_tmp[i] += 0.5*s->U * s->n[OBT_IDX] + (s->U - 2.5*s->J) * n_diff[OBT_IDX]; 
	}
	for(i=3*OBT; i<OBT*(BASIS2+1); i+=(BASIS2+1)) {
		v_tmp[i] -= 0.5*s->U * s->m[OBT_IDX] - 0.5*s->J * m_diff[OBT_IDX]; 
	}
	for(i=OBT + OBT*(BASIS2+1); i<2*OBT*(BASIS2+1); i+=(BASIS2+1)) {
		v_tmp[i] += 0.5*s->U * s->m[OBT_IDX] - 0.5*s->J * m_diff[OBT_IDX]; 
	}	
	/*for(j=0; j<H(BASIS2); j++) {
		if(j%BASIS2 == 0) printf("\n");
		printf("%.3f\t", creal(v_tmp[j]));
	}
	printf("\n");*/
}

void InteractionSubA(Solution *s, lapack_complex_double *v_tmp) {
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

	for(i=0; i<H(BASIS2); i+=(BASIS2+1)) {
		v_tmp[i] += 0.5*s->U * s->n[OBT_IDX]/2 + (s->U - 2.5*s->J) * n_diff[OBT_IDX]/2; 
	}
	for(i=0; i<2*OBT*(BASIS2+1); i+=(BASIS2+1)) {
		v_tmp[i] -= 0.5*s->U * s->m[OBT_IDX]/2 - 0.5*s->J * m_diff[OBT_IDX]/2; 
	}
	for(i=2*OBT*(BASIS2+1); i<4*OBT*(BASIS2+1); i+=(BASIS2+1)) {
		v_tmp[i] += 0.5*s->U * s->m[OBT_IDX]/2 - 0.5*s->J * m_diff[OBT_IDX]/2; 
	}	
}

void OccupationF(double fermi, double *w, lapack_complex_double *v, double *n, double *m, double *e) {
	int i, j;

	memset(n, 0, OBT);
	memset(m, 0, OBT);
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

	memset(n, 0, OBT);
	memset(m, 0, OBT);
	*e = 0;

	for(i=0; i<K3*BASIS2; i++) {
		if(w[i] < fermi) {
			for(j=0; j<OBT; j++) {
				k1 = creal(v[BASIS2*i + j]) - cimag(v[BASIS2*i + j]) * I;
				k2 = creal(v[BASIS2*i + j+2*OBT]) - cimag(v[BASIS2*i + j+2*OBT]) * I;

				n[j] += COMPLEX2(v[BASIS2*i + j]) + COMPLEX2(v[BASIS2*i + j+OBT]) + COMPLEX2(v[BASIS2*i + j+2*OBT]) + COMPLEX2(v[BASIS2*i + j+3*OBT]);
				m[j] += 2 * (k1 * v[BASIS2*i + j+3*OBT] - k2 * v[BASIS2*i + j+OBT]);
				*e += w[i] * (COMPLEX2(v[BASIS2*i + j]) + COMPLEX2(v[BASIS2*i + j+OBT]) + COMPLEX2(v[BASIS2*i + j+2*OBT]) + COMPLEX2(v[BASIS2*i + j+3*OBT]));
				//printf("%f\t%f\t%f\t%f\n", creal(v[BASIS2*i + j]), creal(v[BASIS2*i + j+OBT]), creal(v[BASIS2*i + j+2*OBT]) , creal(v[BASIS2*i + j+3*OBT]));
				//printf("%f\t%f\t%f\t%f\n", cimag(v[BASIS2*i + j]), cimag(v[BASIS2*i + j+OBT]), cimag(v[BASIS2*i + j+2*OBT]) , cimag(v[BASIS2*i + j+3*OBT]));
				//if(cimag(k2 * v[BASIS2*i + j+3*OBT]) > 1e-6) printf("%f\t%f\t%f\t%f\n", creal(k1 * v[BASIS2*i + j+OBT]), cimag(k1 * v[BASIS2*i + j+OBT]), creal(k2 * v[BASIS2*i + j+3*OBT]), cimag(k2 * v[BASIS2*i + j+3*OBT]));
			}
			//printf("\n");
		}
	}
}

void OccupationSubA(double fermi, double *w, lapack_complex_double *v, double *n, double *m, double *e) {
	int i, j;

	memset(n, 0, OBT);
	memset(m, 0, OBT);
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
