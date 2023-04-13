// main/hf.c : Hartree-Fock approximation

#include "hf.h" 

// one magnetized orbital
void Init1(double *m) {m[0] = 1;      m[1] = M_INIT; m[2] = M_INIT;}
void Init2(double *m) {m[0] = M_INIT; m[1] = 1;      m[2] = M_INIT;}
void Init3(double *m) {m[0] = M_INIT; m[1] = M_INIT; m[2] = 1;}

// two magnetized orbitals
void Init4(double *m) {m[0] = 1;      m[1] = 1;      m[2] = M_INIT;}
void Init5(double *m) {m[0] = 1;      m[1] = M_INIT; m[2] = 1;}
void Init6(double *m) {m[0] = M_INIT; m[1] = 1;      m[2] = 1;}

// force symmetry
void SymmetryN(double *n, double *m) {}
void SymmetryA(double *n, double *m) {n[2] = n[0];        m[2] = m[0];}
void SymmetryC(double *n, double *m) {n[2] = n[1];        m[2] = m[1];}
void SymmetryG(double *n, double *m) {n[2] = n[1] = n[0]; m[2] = m[1] = m[0];}

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("%s <save> <type> <J/U> <SOC> <N> <U> <(dos)ep>\n", argv[0]);
		exit(1);
	}
	omp_set_num_threads(1);

	time_t t = time(NULL);
	struct tm *tm = localtime(&t);

	Config c = {.type = argv[2]};
	ReadConfig(&c);

	Solution s = {
		.save = argv[1],
		.type = argv[2],

		.JU   = atof(argv[3]),
		.SOC  = atof(argv[4]),
		.N    = atof(argv[5]),
		.U    = atof(argv[6]),
		.J    = s.JU * s.U,

		.n     = calloc(sizeof(*s.n), c.Nc),
		.m     = calloc(sizeof(*s.m), c.Nc),
		.ns    = 100,
		.ms    = 100,
		.fermi = 100,
		.dntop = 100,
		.gap   = 100,
		.e     = 100
	};
	sprintf(s.runtime, "%d%d%d%d", tm->tm_mon+1, tm->tm_mday, tm->tm_hour, tm->tm_min); 

	char save[1024];
	sprintf(save, "output/%s/%s_JU%.2f_SOC%.2f", s.save, s.type, s.JU, s.SOC);
	if(-access(save, 0)) mkdir(save, 0755);
	s.save = save;

	InitSolution(c, &s);
	if     (strstr(c.type, "1")) Init1(s.m);
	else if(strstr(c.type, "2")) Init2(s.m);
	else if(strstr(c.type, "3")) Init3(s.m);
	else if(strstr(c.type, "4")) Init4(s.m);
	else if(strstr(c.type, "5")) Init5(s.m);
	else if(strstr(c.type, "6")) Init6(s.m);

	void (*Symmetry)(double*, double*);
	if     (strstr(c.type, "A")) Symmetry = SymmetryA;
	else if(strstr(c.type, "C")) Symmetry = SymmetryC;
	else if(strstr(c.type, "G")) Symmetry = SymmetryG;
	else                         Symmetry = SymmetryN;

	int i, is_Q = 0;
	for(i=0; i<DIM; i++) if(c.Q[i]) is_Q = 1;

	void (*Interaction)(Config, Solution*, lapack_complex_double*);
	void (*Basis)(Config, double*, lapack_complex_double*);
	if(is_Q) {
		Interaction = InteractionQ;
		Basis       = BasisQ;
	}
	else {
		Interaction = InteractionN;
		Basis       = BasisN;
	}

	if(argv[7]) {
		double ep = atof(argv[7]);

		char dsn[256], ftype[32], fsn[256];
		sprintf(dsn, "%s/sol", s.save);
		sprintf(ftype, "N%.1f_U%.1f", s.N, s.U);

		DIR *d = opendir(dsn);
		struct dirent *f;
		while(!(f = readdir(d))) {
			if(strstr(f->d_name, ftype)) printf("%s\n", f->d_name);
		}

		//GenDOS(c, &s, fsn, ep, Interaction, Basis);
	}
	else GenSolution(c, &s, Symmetry, Interaction, Basis);

	free(s.n);
	free(s.m);

	return 0;
}
