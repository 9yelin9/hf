// main/tb.c : tight-binding Hamiltonian

#include "hf.h" 

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("%s <type>\n", argv[0]);
		exit(1);
	}
	omp_set_num_threads(72);

	Config c = {.type = argv[1]};
	ReadConfig(&c);

	int i, is_Q = 0;
	for(i=0; i<DIM; i++) if(c.Q[i]) is_Q = 1;

	void(*Fourier)(Config, int, Lattice, double*, lapack_complex_double*);
	if(is_Q) Fourier = FourierQ;
	else     Fourier = FourierN;

	GenTB(c, "g", Fourier);
	GenTB(c, "b", Fourier);
	GenTBBand(c);

	return 0;
}
