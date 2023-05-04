// main/hf.c : Hartree-Fock approximation

#include "hf.h" 

// force symmetry
void SymmetryN(double *n, double *m) {}
void SymmetryA(double *n, double *m) {n[2] = n[0];        m[2] = m[0];}
void SymmetryC(double *n, double *m) {n[2] = n[1];        m[2] = m[1];}
void SymmetryG(double *n, double *m) {n[2] = n[1] = n[0]; m[2] = m[1] = m[0];}

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("%s <sol> <save> <strain> <type> <JU> <N> <U> [(band)Nk/(dos)ep]\n", argv[0]);
		exit(1);
	}
	omp_set_num_threads(1);

	time_t t = time(NULL);
	struct tm *tm = localtime(&t);
	char runtime[16];
	sprintf(runtime, "v%d%d%d%d", tm->tm_mon+1, tm->tm_mday, tm->tm_hour, tm->tm_min); 

	Solution s;
	Config c = {
		.runtime = runtime,
		.sol    = argv[1],
		.save   = argv[2],
		.strain = argv[3],
		.type   = argv[4],

		.JU  = atof(argv[5]),
		.N   = atof(argv[6]),
		.U   = atof(argv[7]),
		.J   = c.JU * c.U
	};

	double option = argv[8] ? atof(argv[8]) : -1;

	sprintf(c.path_save, "output/%s/%s_%s_JU%.2f", c.save, c.strain, c.type, c.JU);
	if(-access(c.path_save, 0)) mkdir(c.path_save, 0755);

	ReadConfig(&c);
	ReadSolution(c, &s);

	// Symmetry
	void (*Symmetry)(double*, double*);
	if     (strstr(c.sym, "A")) Symmetry = SymmetryA;
	else if(strstr(c.sym, "C")) Symmetry = SymmetryC;
	else if(strstr(c.sym, "G")) Symmetry = SymmetryG;
	else if(strstr(c.sym, "N")) Symmetry = SymmetryN;
	else                        Symmetry = SymmetryN;

	// Interaction, Basis
	void (*Interaction)(Config, Solution*, lapack_complex_double*);
	void (*Basis)(Config, double*, lapack_complex_double*);
	if(c.Q[0] + c.Q[1] + c.Q[2] > 1e-6) {
		Interaction = InteractionQ;
		Basis       = BasisQ;
	}
	else {
		Interaction = InteractionN;
		Basis       = BasisN;
	}

	if(option > 1) {
		c.Nkb = (int)option;
		GenBand(c, &s, Interaction, Basis);
	}
	else if(option > 0) GenDOS(c, &s, option, Interaction, Basis);
	else GenSolution(c, &s, Symmetry, Interaction, Basis);

	return 0;
}
