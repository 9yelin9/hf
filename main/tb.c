// main/tb.c : tight-binding

#include "hf.h" 

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("%s <type>\n", argv[0]);
		exit(1);
	}
	omp_set_num_threads(16);

	time_t t = time(NULL);
	struct tm *tm = localtime(&t);

	Config c = {.type = argv[1]};
	ReadConfig(&c);

	void(*Fourier)(Config);
	if(c.Ni > 2) Fourier = FourierN;
	else         Fourier = FourierQ;

	CalcTB(c, "g", Fourier);
	CalcTB(c, "b", Fourier);

	return 0;
}
