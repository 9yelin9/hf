// hf3/hf3.c : make Hartree-Fock approximated 3-band Hubbard model

#include "hf3.h" 

inline void InitA1(double *m) { m[0] = 0.1; m[1] = 1.0; m[2] = 0.1; }
inline void InitA2(double *m) { m[0] = 1.0; m[1] = 0.1; m[2] = 1.0; }
inline void InitC1(double *m) { m[0] = 1.0; m[1] = 0.1; m[2] = 0.1; }
inline void InitC2(double *m) { m[0] = 0.1; m[1] = 1.0; m[2] = 1.0; }

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("%s <name> <type> <J/U> <SOC> <N> <U> : make Hartree-Fock approximated 3-band Hubbard model\n", argv[0]);
		exit(1);
	}

	time_t t = time(NULL);
	struct tm *tm = localtime(&t);

	Solution s = {
		.name = argv[1],
		.type = argv[2],
		.JU = atof(argv[3]),
		.SOC = atof(argv[4]),
		.N = atof(argv[5]),
		.U = atof(argv[6]),
		.J = s.JU * s.U,

		.n = {s.N/3, s.N/3, s.N/3},
		.m = {0.1, 0.1, 0.1},
		.ntot = 100,
		.mtot = 100,
		.fermi = 100,
		.dntop = 100,
		.gap = 100,
		.e = 100
	};

	// set initial m
	if(strstr(s.type, "a1")) InitA1(s.m);
	else if(strstr(s.type, "a2")) InitA2(s.m);
	else if(strstr(s.type, "c1")) InitC1(s.m);
	else if(strstr(s.type, "c2")) InitC2(s.m);

	// make directory
	char dir[128];
	argv[0] = &argv[0][2];
	sprintf(dir, "output/%s/JU%.2f_SOC%.2f", s.name, s.JU, s.SOC);
	if(-access(dir, 0)) mkdir(dir, 0755);

	// set runtime
	sprintf(s.runtime, "%d%d%d%d", tm->tm_mon+1, tm->tm_mday, tm->tm_hour, tm->tm_min); 

	// read info, set B and basis
	ReadInfo(s.name, s.type, s.ctype, &s.q);
	s.B = GetB(s.name);
	s.basis = strstr(s.ctype, "s") ? SINGLE : DOUBLE;

	// set LAPACK parameters
	LParameter lp = {
		.jobz = 'V',
		.uplo = 'L',
		.rwork = (double*)malloc((3*s.basis-2) * sizeof(double)),
		.ln = s.basis,
		.lda = s.basis,
		.lwork = 2*s.basis-1,
		.work = (lapack_complex_double*)malloc(lp.lwork * sizeof(lapack_complex_double))
	};

	// set tight-binding
	lapack_complex_double tbg[s.basis*s.basis * G3], tbb[s.basis*s.basis * s.B];
	char tbgs[64], tbbs[64];
	sprintf(tbgs, "input/%s/tbg_%c.bin", s.name, s.type[0]);
	sprintf(tbbs, "input/%s/tbb_%c.bin", s.name, s.type[0]);
	ReadBin(tbgs, sizeof(lapack_complex_double) * s.basis*s.basis * G3, tbg);
	ReadBin(tbbs, sizeof(lapack_complex_double) * s.basis*s.basis * s.B, tbb);

	// run
	CalcSolution(&s, &lp, tbg);
	MakeBand(&s, &lp, tbb);
	MakeDOS(&s, &lp, tbg);

	free(lp.rwork);
	free(lp.work);

	return 0;
}

