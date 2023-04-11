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
		printf("%s <save> <type> <J/U> <SOC> <N> <U>\n", argv[0]);
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

	char save[1024];
	sprintf(save, "output/%s/%s_JU%.2f_SOC%.2f", s.save, s.type, s.JU, s.SOC);
	if(-access(save, 0)) mkdir(save, 0755);
	s.save = save;

	InitSolution(c, &s);
	if     (strstr(s.type, "1")) Init1(s.m);
	else if(strstr(s.type, "2")) Init2(s.m);
	else if(strstr(s.type, "3")) Init3(s.m);
	else if(strstr(s.type, "4")) Init4(s.m);
	else if(strstr(s.type, "5")) Init5(s.m);
	else if(strstr(s.type, "6")) Init6(s.m);

	void (*Symmetry)(double*, double*);
	if     (strstr(s.type, "A")) Symmetry = SymmetryA;
	else if(strstr(s.type, "C")) Symmetry = SymmetryC;
	else if(strstr(s.type, "G")) Symmetry = SymmetryG;
	else                         Symmetry = SymmetryN;

	void (*Interaction)(Config, Solution*, lapack_complex_double*);
	if     (strstr(c.bas, "q")) Interaction = InteractionQ;
	else if(strstr(c.bas, "s")) Interaction = InteractionS;
	else                        Interaction = InteractionN;

	void (*Basis)(Config, int, lapack_complex_double*, lapack_complex_double*);
	if (strstr(c.bas, "q")) Basis = BasisQ;
	else Basis = BasisN;

	// run
	CalcSolution(c, &s, &lp, System, Interaction, Basis);
	//MakeBand(c, &s, &lp, Interaction, Basis);
	//MakeDOS(c, &s, &lp, Interaction, Basis);

	free(s.n);
	free(s.m);
	return 0;
}
