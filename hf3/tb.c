// hf3/in.c : make tight-binding matrix

#include "hf3.h" 

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("%s tb <name> <type> : make tight-binding matrix\n", argv[0]);
		exit(1);
	}

	char *name = argv[1], *type = argv[2];
	int B = GetB(name);

	CalcTB(B, name, type);

	return 0;
}
