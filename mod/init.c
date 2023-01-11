// mod/init.c : initialize input files - Gauss-Legendre quadrature points and tight-binding Hamitonian

#include "hf3.h" 

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("%s <name> <type> <tb_only> : initialize input files\n", argv[0]);
		exit(1);
	}
	omp_set_num_threads(16);

	Cell c = {
		.name = argv[1],
		.type = argv[2]
	};
	ReadCell(&c);

	int tb_only = atoi(argv[3]);
	char ltype[8] = {};
	if(strstr(c.type, "0")) sprintf(ltype, "s");

	if(!tb_only) {
		CalcQuadPoints(c);

		const int Ni = c.Ni;
		Coord r[Ni];
		memset(r, 0, sizeof(r));

		if(strstr(c.sys, "sc")) {
			r[1].c[0] = 1;
		}
		CalcCoef(c, r);
	}

	void (*Fourier)(Cell, int, Coord, int, Lattice*, lapack_complex_double*);
	if(strstr(c.bas, "q")) Fourier = FourierQ;
	else                   Fourier = Fourier0;

	CalcTB(c, ltype, Fourier);

	return 0;
}
