// main/hf.c : Hartree-Fock approximation

#include "hf.h" 

// one magnetized orbital
void Init1(double *m) {m[0] = M_MAX; m[1] = M_MIN; m[2] = M_MIN;}
void Init2(double *m) {m[0] = M_MIN; m[1] = M_MAX; m[2] = M_MIN;}
void Init3(double *m) {m[0] = M_MIN; m[1] = M_MIN; m[2] = M_MAX;}

// two magnetized orbitals
void Init4(double *m) {m[0] = M_MAX; m[1] = M_MAX; m[2] = M_MIN;}
void Init5(double *m) {m[0] = M_MAX; m[1] = M_MIN; m[2] = M_MAX;}
void Init6(double *m) {m[0] = M_MIN; m[1] = M_MAX; m[2] = M_MAX;}

// force symmetry
void SymmetryN(double *n, double *m) {}
void SymmetryA(double *n, double *m) {n[2] = n[0];        m[2] = m[0];}
void SymmetryC(double *n, double *m) {n[2] = n[1];        m[2] = m[1];}
void SymmetryG(double *n, double *m) {n[2] = n[1] = n[0]; m[2] = m[1] = m[0];}

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("%s <save> <strain> <type> <JU> <SOC> <N> <U> [(band)Nk/(dos)ep]\n", argv[0]);
		exit(1);
	}
	omp_set_num_threads(1);

	time_t t = time(NULL);
	struct tm *tm = localtime(&t);

	Config c = {.strain = argv[2], .type = argv[3]};
	ReadConfig(&c);

	Solution s = {
		.JU   = atof(argv[4]),
		.SOC  = atof(argv[5]),
		.N    = atof(argv[6]),
		.U    = atof(argv[7]),
		.J    = s.JU * s.U,

		.ns    = 100,
		.ms    = 100,
		.fermi = 100,
		.dntop = 100,
		.gap   = 100,
		.e     = 100
	};

	sprintf(s.runtime, "v%d%d%d%d", tm->tm_mon+1, tm->tm_mday, tm->tm_hour, tm->tm_min); 
	sprintf(s.type, "%s", argv[3]);
	sprintf(s.save, "output/%s/%s_%s_JU%.2f_SOC%.2f", argv[1], c.strain, c.type, s.JU, s.SOC);
	if(-access(s.save, 0)) mkdir(s.save, 0755);

	InitSolution(c, &s);
	if     (strstr(c.type, "1")) Init1(s.m);
	else if(strstr(c.type, "2")) Init2(s.m);
	else if(strstr(c.type, "3")) Init3(s.m);
	else if(strstr(c.type, "4")) Init4(s.m);
	else if(strstr(c.type, "5")) Init5(s.m);
	else if(strstr(c.type, "6")) Init6(s.m);

	void (*Symmetry)(double*, double*);
	if     (strstr(c.sym, "A")) Symmetry = SymmetryA;
	else if(strstr(c.sym, "C")) Symmetry = SymmetryC;
	else if(strstr(c.sym, "G")) Symmetry = SymmetryG;
	else if(strstr(c.sym, "N")) Symmetry = SymmetryN;
	else                        Symmetry = SymmetryN;

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

	if(argv[8]) {
		double opt = atof(argv[8]);

		char dsn[256], stype[32], fsn[256];
		sprintf(dsn, "%s/sol", s.save);
		sprintf(stype, "N%.1f_U%.1f", s.N, s.U);

		DIR *d = opendir(dsn);
		struct dirent *f;
		while((f = readdir(d)) != NULL) {
			if(strstr(f->d_name, stype)) break;
		}
		sprintf(fsn, "%s/sol/%s", s.save, f->d_name);

		if(opt > 1) {
			c.Nkb = (int)opt;
			GenBand(c, &s, fsn, Interaction, Basis);
		}
		else GenDOS(c, &s, fsn, opt, Interaction, Basis);
		closedir(d);
	}
	else GenSolution(c, &s, Symmetry, Interaction, Basis);

	return 0;
}
