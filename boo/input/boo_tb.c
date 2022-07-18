// boo/input/boo_tb.c : calculate tight-binding Hamiltonian of BaOsO3 model

#define DOT_F (path.x*i + path.y*j + path.z*k)
#define DOT_A ((path.x + q.x)*i + (path.y + q.y)*j + (path.z + q.z)*k)

#include "../../hf3.h"
#include "../boo.h"

void CalcPathB() {
	FILE *f;

	if((f = fopen("path_b.bin", "wb")) == NULL) {
		printf("path_b.bin fopen FAIL\n");
		exit(1);
	}

	int i = 0, j, high_sym[4] = {100, 100, 140, 170};
	Vector path[BAND];

	for(j=0; j<high_sym[0]; j++) { // G-X
		path[i].x = M_PI * j / high_sym[0];
		path[i].y = 0;
		path[i].z = 0;
		i++;
	}

	for(j=0; j<high_sym[1]; j++) { // X-M
		path[i].x = M_PI;
		path[i].y = M_PI * j / high_sym[1];
		path[i].z = 0;
		i++;
	}

	for(j=0; j<high_sym[2]; j++) { // M-G
		path[i].x = M_PI - M_PI * j / high_sym[2];
		path[i].y = M_PI - M_PI * j / high_sym[2];
		path[i].z = 0;
		i++;
	}

	for(j=0; j<high_sym[3]; j++) { // G-R
		path[i].x = M_PI * j / high_sym[3];
		path[i].y = M_PI * j / high_sym[3];
		path[i].z = M_PI * j / high_sym[3];
		i++;
	}

	fwrite(path, sizeof(path), 1, f); 
	fclose(f);
}

void FourierF(char *type, int path_num, Vector path, Vector q, lapack_complex_double *tb) {
	FILE *f; 

	if((f = fopen("lattice.txt", "r")) == NULL) {
		printf("lattice.txt fopen FAIL\n");
		exit(1);
	}

	int i, j, k, obt1, obt2;
	double tre, tim;
	lapack_complex_double tb_tmp[BASIS1][BASIS1], t, exp;

	memset(tb_tmp, 0, sizeof(tb_tmp));

	while(!feof(f)) {
		fscanf(f, "%d%d%d%d%d%lf%lf\n", &i, &j, &k, &obt1, &obt2, &tre, &tim);

		t = tre + tim * I;
		exp = cos(DOT_F) + sin(DOT_F) * I;

		tb_tmp[obt1-1][obt2-1] += t * exp;
	}
	fclose(f);

	for(i=0; i<OBT; i++) {
		for(j=0; j<OBT; j++) {
			tb_tmp[i + OBT][j + OBT] = tb_tmp[i][j];
		}
	}

	k = 0;
	for(i=0; i<BASIS1; i++) {
		for(j=0; j<BASIS1; j++) {
			tb[H(BASIS1)*path_num + k] = tb_tmp[i][j];
			k++;
		}
	}
}

void FourierA(char *type, int path_num, Vector path, Vector q, lapack_complex_double *tb) {
	FILE *f; 

	if((f = fopen("lattice.txt", "r")) == NULL) {
		printf("lattice.txt fopen FAIL\n");
		exit(1);
	}

	int i, j, k, obt1, obt2;
	double tre, tim;
	lapack_complex_double tb_tmp[BASIS2][BASIS2], t, exp1, exp2;

	memset(tb_tmp, 0, sizeof(tb_tmp));

	while(!feof(f)) {
		fscanf(f, "%d%d%d%d%d%lf%lf\n", &i, &j, &k, &obt1, &obt2, &tre, &tim);

		t = tre + tim * I;
		exp1 = cos(DOT_F) + sin(DOT_F) * I;
		exp2 = cos(DOT_A) + sin(DOT_A) * I;

		tb_tmp[obt1-1][obt2-1] += t * exp1;
		tb_tmp[obt1-1 + OBT][obt2-1 + OBT] += t * exp2;
	}
	fclose(f);

	for(i=0; i<OBT; i++) {
		for(j=0; j<OBT; j++) {
			tb_tmp[i + 2*OBT][j + 2*OBT] = tb_tmp[i][j];
			tb_tmp[i + 3*OBT][j + 3*OBT] = tb_tmp[i + OBT][j + OBT];
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

void FourierSubA(char *type, int path_num, Vector path, Vector q, lapack_complex_double *tb) {
	FILE *fon, *fsup, *frem; 

	if((fon = fopen("lattice_onsite.txt", "r")) == NULL) {
		printf("lattice_onsite.txt fopen FAIL\n");
		exit(1);
	}
	if((fsup = fopen("lattice_super.txt", "r")) == NULL) {
		printf("lattice_SUPER.txt fopen FAIL\n");
		exit(1);
	}
	if((frem = fopen("lattice_remain.txt", "r")) == NULL) {
		printf("lattice_remain.txt fopen FAIL\n");
		exit(1);
	}

	int i, j, k, obt1, obt2;
	double tre, tim;
	lapack_complex_double tb_tmp[BASIS2][BASIS2], t, exp;

	memset(tb_tmp, 0, sizeof(tb_tmp));

	while(!feof(fon)) {
		fscanf(fon, "%d%d%d%d%d%lf%lf\n", &i, &j, &k, &obt1, &obt2, &tre, &tim);
		tb_tmp[obt1-1][obt2-1] += tre + tim * I;
		tb_tmp[obt1-1 + OBT][obt2-1 + OBT] += tre + tim * I;
	}
	fclose(fon);

	while(!feof(fsup)) {
		fscanf(fsup, "%d%d%d%d%d%lf%lf\n", &i, &j, &k, &obt1, &obt2, &tre, &tim);
		tb_tmp[obt1-1][obt2-1 + OBT] += tre + tim * I;
	}
	fclose(fsup);

	int (*Condition)(int, int, int);

	if(strstr(type, "a")) Condition = ConditionA;
	else if(strstr(type, "c")) Condition = ConditionC;
	else if(strstr(type, "g")) Condition = ConditionG;
	else {
		printf("'%s' is wrong type\n", type);
		exit(1);
	}

	while(!feof(frem)) {
		fscanf(frem, "%d%d%d%d%d%lf%lf\n", &i, &j, &k, &obt1, &obt2, &tre, &tim);

		t = tre + tim * I;
		exp = cos(DOT_F) + sin(DOT_F) * I;

		if(Condition(i, j, k) % 2 == 0) {
			tb_tmp[obt1-1][obt2-1] += t * exp;
			tb_tmp[obt1-1 + OBT][obt2-1 + OBT] += t * exp;
		}
		else {
			exp = cos(path.x*(i-1) + path.y*j + path.z*k) + sin(path.x*(i-1) + path.y*j + path.z*k) * I;
			tb_tmp[obt1-1][obt2-1 + OBT] += t * exp;
		}
	}
	fclose(frem);

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

void Unfold(char *type, Vector *path, Vector q) {
	FILE  *f;
	char fs[32];

	sprintf(fs, "band_%s_unfold.txt", type);

	if((f = fopen(fs, "w")) == NULL) {
		printf("%s fopen FAIL\n", fs);
		exit(1);
	}

	int i, j, k, m, n;
	double p2;
	Vector r[SUPER] = {{1, 0, 0}, {0, 0, 0}};
	lapack_complex_double p, exp;

	double w[BAND*BASIS2];
	lapack_complex_double v[HB(BASIS2)];

	CalcEigenTB(type, w, v);

	for(i=0; i<BAND; i++) {
		fprintf(f, "%4d", i);

		for(j=0; j<BASIS2; j++) {
			p2 = 0;
			for(k=0; k<OBT; k++) {
				for(m=0; m<2; m++) {
					p = 0;
					for(n=0; n<SUPER; n++) {
						exp = cos(DOT_UNFOLD) - sin(DOT_UNFOLD) * I;	
						p += v[BASIS2*(BASIS2*i + j) + OBT*(2*m + n) + k] * exp;
					}
					p2 += COMPLEX2(p) / SUPER;
				}
			}
			
			if(p2 > 0.5) fprintf(f, "%12f", w[BASIS2*i + j]);
		}
		fprintf(f, "\n");
	}

	fclose(f);
}

int main(int argc, char *argv[]) {
	if(argc != 3) {
		printf("Usage : %s <type> <is_unfold>\n", argv[0]);
		exit(1);
	}

	char *type = argv[1];

	if(strstr(type, "path")) {
		CalcPathK();
		CalcPathB();
	}
	else {
		int is_unfold = atoi(argv[2]);
		Vector path_k[K3], path_b[BAND], q;

		ReadPath(path_k, path_b);

		if(strstr(type, "f")) {
			q.x = 0;
			q.y = 0;
			q.z = 0;
		}
		else if(strstr(type, "a")) {
			q.x = M_PI;
			q.y = 0;
			q.z = 0;
		}
		else if(strstr(type, "c")) {
			q.x = M_PI;
			q.y = M_PI;
			q.z = 0;
		}
		else if(strstr(type, "g")) {
			q.x = M_PI;
			q.y = M_PI;
			q.z = M_PI;
		}

		if(is_unfold) {
			Unfold(type, path_b, q);
		}
		else {
			char fks[32], fbs[32];

			sprintf(fks, "tbk%d_%s.bin", K, type);
			sprintf(fbs, "tbb_%s.bin", type);

			MakeTB(type, fks, path_k, q);
			MakeTB(type, fbs, path_b, q);
		}
	}

	return 0;
}
