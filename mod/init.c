// mod/init.c : initialize input files - Gauss-Legendre quadrature points and tight-binding Hamitonian

#include "hf3.h" 

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("%s <name> <type> : initialize input files\n", argv[0]);
		exit(1);
	}

	Cell c = {
		.name = argv[1],
		.type = argv[2]
	};
	ReadCell(&c);

	CalcQuadPoints(c);

	const int Ni = c.Ni;
	Coord r[Ni];
	memset(r, 0, sizeof(r));

	if(strstr(c.sys, "sc")) {
		r[1].c[0] = 1;
	}
	CalcCoef(c, r);

	void (*Fourier)(Cell, int, Coord, int, Lattice*, lapack_complex_double*);
	if(strstr(c.bas, "q")) Fourier = FourierQ;
	else Fourier = Fourier0;

	CalcTB(c, Fourier);

	return 0;
}
