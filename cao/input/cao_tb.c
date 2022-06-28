// cao/input/input.c : calculate tight-binding Hamiltonian of CuAl2O4 model 

#include "../../hf3.h"
#include "../cao.h"

void CalcBAND(Vector *v) {
	int i = 0, j, high_sym[6] = {150, 200, 50, 100, 75, 225};

	for(j=0; j<high_sym[0]; j++) { // L-G
		v[i].x = M_PI - M_PI * j / high_sym[0];
		v[i].y = M_PI - M_PI * j / high_sym[0];
		v[i].z = M_PI - M_PI * j / high_sym[0];
		i++;
	}

	for(j=0; j<high_sym[1]; j++) { // G-X
		v[i].x = M_PI * j / high_sym[1];
		v[i].y = 0;
		v[i].z = M_PI * j / high_sym[1];
		i++;
	}

	for(j=0; j<high_sym[2]; j++) { // X-W
		v[i].x = M_PI;
		v[i].y = (M_PI/2) * j / high_sym[2];
		v[i].z = M_PI + (M_PI/2) * j / high_sym[2];
		i++;
	}

	for(j=0; j<high_sym[3]; j++) { // W-L
		v[i].x = M_PI;
		v[i].y = (M_PI/2) + (M_PI/2) * j / high_sym[3];
		v[i].z = (3*M_PI/2) - (M_PI/2) * j / high_sym[3];
		i++;
	}

	for(j=0; j<high_sym[4]; j++) { // L-K
		v[i].x = M_PI - (M_PI/4) * j / high_sym[4];
		v[i].y = M_PI + (M_PI/2) * j / high_sym[4];
		v[i].z = M_PI - (M_PI/4) * j / high_sym[4];
		i++;
	}

	for(j=0; j<high_sym[5]; j++) { // K-G
		v[i].x = (3*M_PI/4) - (3*M_PI/4) * j / high_sym[5];
		v[i].y = (3*M_PI/2) - (3*M_PI/2) * j / high_sym[5];
		v[i].z = (3*M_PI/4) - (3*M_PI/4) * j / high_sym[5];
		i++;
	}
}

void FourierF(FILE *f, const int num, const int basis, Vector v, Vector q, lapack_complex_double *tb) {
	int i, j, k, obt1, obt2;
	double tre, tim;
	lapack_complex_double tb_tmp[basis][basis];

	memset(tb_tmp, 0, sizeof(tb_tmp));

	while(!feof(f)) {
		fscanf(f, "%d%d%d%d%d%lf%lf\n", &i, &j, &k, &obt1, &obt2, &tre, &tim);
		tb_tmp[obt1-1][obt2-1] += (tre * cos(v.x*i + v.y*j + v.z*k) - tim * sin(v.x*i + v.y*j + v.z*k))
								+ (tre * sin(v.x*i + v.y*j + v.z*k) + tim * cos(v.x*i + v.y*j + v.z*k)) * I; 
	}

	k = 0;
	for(i=0; i<basis; i++) {
		for(j=0; j<basis; j++) {
			tb[H(basis)*num + k] = tb_tmp[i][j];
			k++;
		}
	}
}

void FourierA(FILE *f, const int num, const int basis, Vector v, Vector q, lapack_complex_double *tb) {
}

int main(int argc, char *argv[]) {
	if(argc != 2) {
		printf("Usage : %s <type>\n", argv[0]);
		exit(1);
	}

	char *type = argv[1];

	int basis;
	char fk_name[32], fb_name[32];
	Vector vk[K3], vb[BAND];

	sprintf(fk_name, "tb_%s_K%d.bin", type, K);
	sprintf(fb_name, "tb_%s_BAND.bin", type);

	CalcK(vk);
	CalcBAND(vb);

	if(strstr(type, "f")) {
		basis = 6;
	}
	else {
		basis = 12;
	}

	lapack_complex_double tbk[HK(K3, basis)], tbb[HB(BAND, basis)];

	CalcTB(fk_name, basis, K3, HK(K3, basis), vk, tbk);
	CalcTB(fb_name, basis, BAND, HB(BAND, basis), vb, tbb);
	MakeTB(type, basis, tbb);

	return 0;
}
