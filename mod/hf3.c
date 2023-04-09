// mod/hf3.c : make Hartree-Fock approximated 3-band Hubbard model

#include "hf3.h" 

void Init0(Cell c, Solution *s) {
	int i;
	for(i=0; i<c.Nc; i++) {
		s->n[i] = s->N / c.Nc;
		s->m[i] = 0.1;
	}
}

// one magnetized orbital
void Init1(double *m) {m[0] = 1.0; m[1] = 0.1; m[2] = 0.1;}
void Init2(double *m) {m[0] = 0.1; m[1] = 1.0; m[2] = 0.1;}
void Init3(double *m) {m[0] = 0.1; m[1] = 0.1; m[2] = 1.0;}

// two magnetized orbitals
void Init4(double *m) {m[0] = 1.0; m[1] = 1.0; m[2] = 0.1;}
void Init5(double *m) {m[0] = 1.0; m[1] = 0.1; m[2] = 1.0;}
void Init6(double *m) {m[0] = 0.1; m[1] = 1.0; m[2] = 1.0;}

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("%s <name> <save> <J/U> <SOC> <type> <N> <U> : make Hartree-Fock approximated 3-band Hubbard model\n", argv[0]);
		exit(1);
	}
	omp_set_num_threads(1);

	time_t t = time(NULL);
	struct tm *tm = localtime(&t);

	Cell c = {
		.name = argv[1],
		.save = argv[2],
		.type = argv[5],
		.eta  = 0.1
	};
	ReadCell(&c);

	Solution s = {
		.JU  = atof(argv[3]),
		.SOC = atof(argv[4]),
		.N   = atof(argv[6]),
		.U   = atof(argv[7]),
		.J   = s.JU * s.U,

		.n     = (double*)malloc(sizeof(double) * c.Nc),
		.m     = (double*)malloc(sizeof(double) * c.Nc),
		.ns    = 100,
		.ms    = 100,
		.fermi = 100,
		.dntop = 100,
		.gap   = 100,
		.e     = 100
	};
	sprintf(s.runtime, "%d%d%d%d", tm->tm_mon+1, tm->tm_mday, tm->tm_hour, tm->tm_min); 

	Init0(c, &s);
	if     (strstr(c.type, "1")) Init1(s.m);
	else if(strstr(c.type, "2")) Init2(s.m);
	else if(strstr(c.type, "3")) Init3(s.m);
	else if(strstr(c.type, "4")) Init4(s.m);
	else if(strstr(c.type, "5")) Init5(s.m);
	else if(strstr(c.type, "6")) Init6(s.m);

	char dir[1024];
	sprintf(dir, "output/%s/JU%.2f_SOC%.2f_%s", c.save, s.JU, s.SOC, c.type);
	if(-access(dir, 0)) mkdir(dir, 0755);

	LAPACK lp = {
		.jobz  = 'V',
		.uplo  = 'L',
		.rwork = (double*)malloc(sizeof(double) * (3*c.Nb-2)),
		.ln    = c.Nb,
		.lda   = c.Nb,
		.lwork = 2*c.Nb-1,
		.work  = (lapack_complex_double*)malloc(sizeof(lapack_complex_double) * lp.lwork)
	};

	void (*System)(double*, double*);
	if(strstr(c.sys, "sc")) {
		if     (strstr(c.type, "a")) System = SystemScA;
		else if(strstr(c.type, "c")) System = SystemScC;
		else if(strstr(c.type, "g")) System = SystemScG;
		else System = SystemN;
	}
	else System = SystemN;

	void (*Interaction)(Cell, Solution*, lapack_complex_double*);
	if     (strstr(c.bas, "q")) Interaction = InteractionQ;
	else if(strstr(c.bas, "s")) Interaction = InteractionS;
	else Interaction = InteractionN;

	void (*Basis)(Cell, int, lapack_complex_double*, lapack_complex_double*);
	if (strstr(c.bas, "q")) Basis = BasisQ;
	else Basis = BasisN;

	// run
	CalcSolution(c, &s, &lp, System, Interaction, Basis);
	//MakeBand(c, &s, &lp, Interaction, Basis);
	//MakeDOS(c, &s, &lp, Interaction, Basis);

	free(s.n);
	free(s.m);
	free(lp.rwork);
	free(lp.work);

	return 0;
}
