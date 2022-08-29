// cao/input/cao_tb.c : calculate tight-binding Hamiltonian of CuAl2O4 model 

#define DOT_F(lat) (path.x*lat.i + path.y*lat.j + path.z*lat.k)
#define DOT_A(lat) ((path.x + q.x)*lat.i + (path.y + q.y)*lat.j + (path.z + q.z)*lat.k)
#define DOT_SA(lat) (path.x*(lat.i-1) + path.y*lat.j + path.z*lat.k)
#define SUB_VEC(a, b) ((a.x - b.x) + (a.y - b.y) + (a.z - b.z))

#include "../../hf3.h"
#include "../cao.h"

void CalcPathB() {
	FILE *f;

	if((f = fopen("path_b.bin", "wb")) == NULL) {
		printf("path_b.bin fopen FAIL\n");
		exit(1);
	}

	int i = 0, j, high_sym[6] = {150, 200, 50, 100, 75, 225};
	Vector path[BAND];

	for(j=0; j<high_sym[0]; j++) { // L-G
		path[i].x = M_PI - M_PI * j / high_sym[0];
		path[i].y = M_PI - M_PI * j / high_sym[0];
		path[i].z = M_PI - M_PI * j / high_sym[0];
		i++;
	}

	for(j=0; j<high_sym[1]; j++) { // G-X
		path[i].x = M_PI * j / high_sym[1];
		path[i].y = 0;
		path[i].z = M_PI * j / high_sym[1];
		i++;
	}

	for(j=0; j<high_sym[2]; j++) { // X-W
		path[i].x = M_PI;
		path[i].y = (M_PI/2) * j / high_sym[2];
		path[i].z = M_PI + (M_PI/2) * j / high_sym[2];
		i++;
	}

	for(j=0; j<high_sym[3]; j++) { // W-L
		path[i].x = M_PI;
		path[i].y = (M_PI/2) + (M_PI/2) * j / high_sym[3];
		path[i].z = (3*M_PI/2) - (M_PI/2) * j / high_sym[3];
		i++;
	}

	for(j=0; j<high_sym[4]; j++) { // L-K
		path[i].x = M_PI - (M_PI/4) * j / high_sym[4];
		path[i].y = M_PI + (M_PI/2) * j / high_sym[4];
		path[i].z = M_PI - (M_PI/4) * j / high_sym[4];
		i++;
	}

	for(j=0; j<high_sym[5]; j++) { // K-G
		path[i].x = (3*M_PI/4) - (3*M_PI/4) * j / high_sym[5];
		path[i].y = (3*M_PI/2) - (3*M_PI/2) * j / high_sym[5];
		path[i].z = (3*M_PI/4) - (3*M_PI/4) * j / high_sym[5];
		i++;
	}

	fwrite(path, sizeof(path), 1, f);
	fclose(f);
}

void FourierF(char *type, int path_num, int lat_len, Vector path, Vector q, Lattice *lat, lapack_complex_double *tb) {
	int i, j, k;
	lapack_complex_double tb_tmp[BASIS1][BASIS1], t, exp;

	memset(tb_tmp, 0, sizeof(tb_tmp));

	for(i=0; i<lat_len; i++) {
		t = lat[i].tre + lat[i].tim * I;
		exp = cos(DOT_F(lat[i])) + sin(DOT_F(lat[i])) * I;

		tb_tmp[lat[i].obt1-1][lat[i].obt2-1] += t * exp;
	}

	k = 0;
	for(i=0; i<BASIS1; i++) {
		for(j=0; j<BASIS1; j++) {
			tb[H(BASIS1)*path_num + k] = tb_tmp[i][j];
			k++;
		}
	}
}

void FourierA(char *type, int path_num, int lat_len, Vector path, Vector q, Lattice *lat, lapack_complex_double *tb) {
	int i, j, k;
	lapack_complex_double tb_tmp[BASIS2][BASIS2], t, exp1, exp2;

	memset(tb_tmp, 0, sizeof(tb_tmp));

	for(i=0; i<lat_len; i++) {
		t = lat[i].tre + lat[i].tim * I;
		exp1 = cos(DOT_F(lat[i])) + sin(DOT_F(lat[i])) * I;
		exp2 = cos(DOT_A(lat[i])) + sin(DOT_A(lat[i])) * I;

		tb_tmp[lat[i].obt1-1][lat[i].obt2-1] += t * exp1;
		tb_tmp[lat[i].obt1-1 + OBT][lat[i].obt2-1 + OBT] += t * exp2;
	}

	for(i=0; i<2*OBT; i++) {
		for(j=0; j<2*OBT; j++) {
			tb_tmp[i + 2*OBT][j + 2*OBT] = tb_tmp[i][j];
		}
	}

	k = 0;
	for(i=0; i<BASIS2; i++) {
		for(j=0; j<BASIS2; j++) {
			tb[H(BASIS2)*path_num + k] = tb_tmp[i][j];
			k++;
		}
	}
}

int ConditionA(int i, int j, int k) { return i; }
int ConditionC(int i, int j, int k) { return i + j; }
int ConditionG(int i, int j, int k) { return i + j - k; }

void FourierSubA(char *type, int path_num, int lat_len, Vector path, Vector q, Lattice *lat, lapack_complex_double *tb) {
	int i, j, k;
	lapack_complex_double tb_tmp[BASIS2][BASIS2], t, exp;

	memset(tb_tmp, 0, sizeof(tb_tmp));

	int (*Condition)(int, int, int);

	if(strstr(type, "a")) Condition = ConditionA;
	else if(strstr(type, "c")) Condition = ConditionC;
	else if(strstr(type, "g")) Condition = ConditionG;
	else {
		printf("'%s' is wrong type\n", type);
		exit(1);
	}

	for(i=0; i<lat_len; i++) {
		t = lat[i].tre + lat[i].tim * I;
		exp = cos(DOT_F(lat[i])) + sin(DOT_F(lat[i])) * I;

		if(Condition(lat[i].i, lat[i].j, lat[i].k) % 2 == 0) {
			tb_tmp[lat[i].obt1-1][lat[i].obt2-1] += t * exp;
			tb_tmp[lat[i].obt1-1 + OBT][lat[i].obt2-1 + OBT] += t * exp;
		}
		else {
			exp = cos(DOT_SA(lat[i])) + sin(DOT_SA(lat[i])) * I;
			tb_tmp[lat[i].obt1-1][lat[i].obt2-1 + OBT] += t * exp;
		}
	}

	for(i=0; i<2*OBT; i++) {
		for(j=0; j<2*OBT; j++) {
			tb_tmp[i + 2*OBT][j + 2*OBT] = tb_tmp[i][j];
		}
	}

	k = 0;
	for(i=0; i<BASIS2; i++) {
		for(j=0; j<BASIS2; j++) {
			tb[H(BASIS2)*path_num + k] = tb_tmp[i][j];
			k++;
		}
	}
}

