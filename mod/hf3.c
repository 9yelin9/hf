// mod/hf3.c : make Hartree-Fock approximated 3-band Hubbard model

#include "hf3.h" 

void Init0(double *m)  {m[0] = 1.0; m[1] = 0.1; m[2] = 0.1;}
void Init1(double *m)  {m[0] = 0.1; m[1] = 1.0; m[2] = 0.1;}
void Init2(double *m)  {m[0] = 0.1; m[1] = 0.1; m[2] = 1.0;}

void Init01(double *m) {m[0] = 1.0; m[1] = 1.0; m[2] = 0.1;}
void Init02(double *m) {m[0] = 1.0; m[1] = 0.1; m[2] = 1.0;}
void Init12(double *m) {m[0] = 0.1; m[1] = 1.0; m[2] = 1.0;}

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("%s <name> <type> <J/U> <SOC> <N> <U> : make Hartree-Fock approximated 3-band Hubbard model\n", argv[0]);
		exit(1);
	}

	time_t t = time(NULL);
	struct tm *tm = localtime(&t);

	Cell c = {
		.name = argv[1],
		.type = argv[2]
	};
	ReadCell(&c);

	Solution s = {
		.JU = atof(argv[3]),
		.SOC = atof(argv[4]),
		.N = atof(argv[5]),
		.U = atof(argv[6]),
		.J = s.JU * s.U,

		.n = (double*)malloc(sizeof(double) * c.Nc),
		.m = (double*)malloc(sizeof(double) * c.Nc),
		.ns = 100,
		.ms = 100,
		.fermi = 100,
		.dntop = 100,
		.gap = 100,
		.e = 100
	};
	sprintf(s.runtime, "%d%d%d%d", tm->tm_mon+1, tm->tm_mday, tm->tm_hour, tm->tm_min); 

	int i;
	for(i=0; i<c.Nc; i++) {
		s.n[i] = s.N / c.Nc;
		s.m[i] = 0.1;
	}

	if(strstr(c.sys, "sc")) {
		if     (strstr(c.type, "a1")) Init1(s.m);
		else if(strstr(c.type, "a2")) Init02(s.m);

		else if(strstr(c.type, "c1")) Init0(s.m);
		else if(strstr(c.type, "c2")) Init12(s.m);
	}
	if(strstr(c.sys, "fcc")) {
		if     (strstr(c.type, "f1")) Init0(s.m);
		else if(strstr(c.type, "f2")) Init1(s.m);
		else if(strstr(c.type, "f3")) Init2(s.m);
		else if(strstr(c.type, "f4")) Init01(s.m);
		else if(strstr(c.type, "f5")) Init02(s.m);
		else if(strstr(c.type, "f6")) Init12(s.m);

		else if(strstr(c.type, "a1")) Init0(s.m);
		else if(strstr(c.type, "a2")) Init1(s.m);
		else if(strstr(c.type, "a3")) Init2(s.m);
		else if(strstr(c.type, "a4")) Init01(s.m);
		else if(strstr(c.type, "a5")) Init02(s.m);
		else if(strstr(c.type, "a6")) Init12(s.m);
	}

	char dir[1024];
	sprintf(dir, "output/%s/JU%.2f_SOC%.2f", c.name, s.JU, s.SOC);
	if(-access(dir, 0)) mkdir(dir, 0755);

	LAPACK lp = {
		.jobz = 'V',
		.uplo = 'L',
		.rwork = (double*)malloc(sizeof(double) * (3*c.Nb-2)),
		.ln = c.Nb,
		.lda = c.Nb,
		.lwork = 2*c.Nb-1,
		.work = (lapack_complex_double*)malloc(sizeof(lapack_complex_double) * lp.lwork)
	};

	void (*System)(double*, double*);
	if(strstr(c.sys, "sc")) {
		if     (strstr(c.type, "a")) System = SystemScA;
		else if(strstr(c.type, "c")) System = SystemScC;
		else if(strstr(c.type, "g")) System = SystemScG;
		else System = System0;
	}
	else System = System0;

	void (*Interaction)(Cell, Solution*, lapack_complex_double*);
	if     (strstr(c.bas, "q")) Interaction = InteractionQ;
	else if(strstr(c.bas, "s")) Interaction = InteractionS;
	else Interaction = Interaction0;

	void (*Basis)(Cell, int, lapack_complex_double*, lapack_complex_double*);
	if (strstr(c.bas, "q")) Basis = BasisQ;
	else Basis = Basis0;

	// run
	CalcSolution(c, &s, &lp, System, Interaction, Basis);
	MakeBand(c, &s, &lp, Interaction, Basis);
	MakeDOS(c, &s, &lp, Interaction, Basis);

	free(s.n);
	free(s.m);
	free(lp.rwork);
	free(lp.work);

	return 0;
}
