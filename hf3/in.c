// hf3/in.c : make band path points and basis transform coefficients

#include "hf3.h" 

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("%s <name> : make band path points and basis transform coefficients\n", argv[0]);
		exit(1);
	}

	char *name = argv[1];
	int B = GetB(name);

	if(strstr(name, "baoso3")) CalcPathBaOsO3(B);
	else if(strstr(name, "cual2o4")) CalcPathCuAl2O4(B);
	else {
		printf("%s is not available\n");
		exit(1);
	}

	CalcGauss(name);
	CalcCoef(B, name);

	return 0;
}
