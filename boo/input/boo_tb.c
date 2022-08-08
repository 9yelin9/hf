// boo/input/boo_tb.c : calculate tight-binding Hamiltonian of BaOsO3 model

#define DOT_F(lat) (path.x*lat.i + path.y*lat.j + path.z*lat.k)
#define DOT_A(lat) ((path.x + q.x)*lat.i + (path.y + q.y)*lat.j + (path.z + q.z)*lat.k)
#define DOT_SA(lat) (path.x*(lat.i-1) + path.y*lat.j + path.z*lat.k)
#define SUB_VEC(a, b) ((a.x - b.x) + (a.y - b.y) + (a.z - b.z))

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

void FourierF(char *type, int path_num, int lat_len, Vector path, Vector q, Lattice *lat, lapack_complex_double *tb) {
	int i, j, k;
	lapack_complex_double tb_tmp[BASIS1][BASIS1], t, exp;

	memset(tb_tmp, 0, sizeof(tb_tmp));

	for(i=0; i<lat_len; i++) {
		t = lat[i].tre + lat[i].tim * I;
		exp = cos(DOT_F(lat[i])) + sin(DOT_F(lat[i])) * I;

		tb_tmp[lat[i].obt1-1][lat[i].obt2-1] += t * exp;
	}

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

/*void Unfold(char *type, Vector *path, Vector q) {
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

	CalcEigenTB(type, "b", w, v);

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
}*/

/*void TransBasis(char *type, char *path_type) {
	FILE *f1, *f2, *f3;

	f1 = fopen("tbb_a.bin", "rb");
	f2 = fopen("tbb_sa.bin", "rb");
	f3 = fopen("path_b.bin", "rb");

	int i, j;
	Vector path[BAND], vk, vq, va, vb;
	lapack_complex_double a[HB(BASIS2)], sa[HB(BASIS2)];

	fread(a, sizeof(a), 1, f1);
	fread(sa, sizeof(sa), 1, f2);
	fread(path, sizeof(path), 1, f3);

	fclose(f1);
	fclose(f2);
	fclose(f3);

	for(i=0; i<3; i++) {
		vk = vq = va = vb = path[i];
		va.x = path[i].x * 2;

		[>printf("%f\t%f\t%f\n", vk.x, vk.y, vk.z);
		printf("%f\t%f\t%f\n", va.x, va.y, va.z);<]
		
		for(j=0; j<H(BASIS2); j++) {
			if(j % BASIS2 == 0) printf("\n");
			printf("%9.6f", creal(a[H(BASIS2)*i + j]));
		}	
		printf("\n");

		for(j=0; j<H(BASIS2); j++) {
			if(j % BASIS2 == 0) printf("\n");
			printf("%9.6f", creal(sa[H(BASIS2)*i + j]));
		}	
		printf("\n");


	}
}*/

/*void TransBasis(char *type, char *path_type) {
	FILE *f;
	char fs[16];

	if(strstr(path_type, "k")) sprintf(fs, "tk%d_%s.bin", K, type);
	else sprintf(fs, "tb_%s.bin", type);

	if((f = fopen(fs, "wb")) == NULL) {
		printf("%s fopen FAIL\n", fs);
		exit(1);
	}

	const int path_num = strstr(path_type, "k") ? K3 : BAND;

	int i, j, k;
	char typek[8], types[8];
	double wk[path_num * BASIS2], ws[path_num * BASIS2];
	lapack_complex_double vk[path_num * H(BASIS2)], vs[path_num * H(BASIS2)], vk_[H(BASIS2)], vs_[H(BASIS2)];
	lapack_complex_double t[path_num * 2*H(BASIS2)];

	lapack_int ln = BASIS2, lda = ln, lwork = ln, ipiv[ln], info;
	lapack_complex_double work[ln];

	sprintf(typek, "%s", type);
	sprintf(types, "s%s", type);

	CalcEigenTB(typek, path_type, wk, vk);
	CalcEigenTB(types, path_type, ws, vs);

	memset(t, 0, sizeof(t));

	[>FILE *fk = fopen("tbk16_a.bin", "rb"), *fsub = fopen("tbk16_sa.bin", "rb");
	lapack_complex_double hk[HK(BASIS2)], hs[HK(BASIS2)], v1[H(BASIS2)], v2[H(BASIS2)];
	fread(hk, sizeof(hk), 1, fk);
	fread(hs, sizeof(hs), 1, fsub);
	fclose(fk);
	fclose(fsub);<]

	for(i=0; i<path_num; i++) {
		for(j=0; j<H(BASIS2); j++) {
			vk_[j] = vk[H(BASIS2)*i + j];
			vs_[j] = vs[H(BASIS2)*i + j];
		}

		LAPACK_zgetrf(&ln, &ln, vk_, &lda, ipiv, &info);	
		LAPACK_zgetri(&ln, vk_, &lda, ipiv, work, &lwork, &info);

		LAPACK_zgetrf(&ln, &ln, vs_, &lda, ipiv, &info);	
		LAPACK_zgetri(&ln, vs_, &lda, ipiv, work, &lwork, &info);

		for(j=0; j<H(BASIS2); j++) {
			for(k=0; k<BASIS2; k++) {
				t[H(BASIS2)*(2*i)   + j] += vs_[(j/BASIS2)*BASIS2 + k] * vk[H(BASIS2)*i + (j%BASIS2) + BASIS2*k];
				t[H(BASIS2)*(2*i+1) + j] += vk_[(j/BASIS2)*BASIS2 + k] * vs[H(BASIS2)*i + (j%BASIS2) + BASIS2*k];
			}
		}

		[>memset(v1, 0, sizeof(v1));
		memset(v2, 0, sizeof(v2));

		for(j=0; j<H(BASIS2); j++) {
			for(k=0; k<BASIS2; k++) {
				v1[j] += t[H(BASIS2)*(2*i) + (j/BASIS2)*BASIS2 + k] * hk[H(BASIS2)*i + (j%BASIS2) + BASIS2*k];
			}
		}

		for(j=0; j<H(BASIS2); j++) {
			for(k=0; k<BASIS2; k++) {
				v2[j] += v1[(j/BASIS2)*BASIS2 + k] * t[H(BASIS2)*(2*i+1) + (j%BASIS2) + BASIS2*k];
			}
		}

		for(j=0; j<H(BASIS2); j++) {
			if(j % BASIS2 == 0) printf("\n");
			printf("%.2f\t", creal(v2[j]));
		}
		printf("\n");
		for(j=0; j<H(BASIS2); j++) {
			if(j % BASIS2 == 0) printf("\n");
			printf("%.2f\t", creal(hs[H(BASIS2)*i + j]));
		}
		printf("==\n");<]
	}

	fwrite(t, sizeof(t), 1, f);
	fclose(f);
}*/

int main(int argc, char *argv[]) {
	if(argc != 2) {
		printf("Usage : %s <type>\n", argv[0]);
		exit(1);
	}

	char *type = argv[1];
	Vector q;

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
	else {
		printf("'%s' is wrong type\n", type);
		exit(1);
	}

	if(strstr(type, "path")) {
		CalcPathK();
		CalcPathB();
	}
	/*else if(strstr(type, "basis")) {
		TransBasis("a", "k");	
		TransBasis("a", "b");	
		TransBasis("c", "k");	
		TransBasis("c", "b");	
		TransBasis("g", "k");	
		TransBasis("g", "b");
	}*/
	else {
		char fks[32], fbs[32];
		Vector path_k[K3], path_b[BAND];

		sprintf(fks, "tbk%d_%s.bin", K, type);
		sprintf(fbs, "tbb_%s.bin", type);

		ReadPath(path_k, path_b);

		MakeBandTB(type, fks, path_k, q);
		MakeBandTB(type, fbs, path_b, q);
	}

	return 0;
}
