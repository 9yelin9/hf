// boo/input/boo_tb.c : calculate tight-binding Hamiltonian of BaOsO3 model

#define DOT0(k, g, r) ((k.x + g.x)*r.x + (k.y + g.y)*r.y + (k.z + g.z)*r.z)
#define DOT(v, k, g, r) ((creal(v) * cos(DOT0(k, g, r)) + cimag(v) * sin(DOT0(k, g, r))) - (creal(v) * sin(DOT0(k, g, r)) - cimag(v) * cos(DOT0(k, g, r))) * I)

#include "../../hf3.h"
#include "../boo.h"

void CalcBAND(Vector *v) {
	int i = 0, j, high_sym[4] = {100, 100, 140, 170};

	for(j=0; j<high_sym[0]; j++) { // G-X
		v[i].x = M_PI * j / high_sym[0];
		v[i].y = 0;
		v[i].z = 0;
		i++;
	}

	for(j=0; j<high_sym[1]; j++) { // X-M
		v[i].x = M_PI;
		v[i].y = M_PI * j / high_sym[1];
		v[i].z = 0;
		i++;
	}

	for(j=0; j<high_sym[2]; j++) { // M-G
		v[i].x = M_PI - M_PI * j / high_sym[2];
		v[i].y = M_PI - M_PI * j / high_sym[2];
		v[i].z = 0;
		i++;
	}

	for(j=0; j<high_sym[3]; j++) { // G-R
		v[i].x = M_PI * j / high_sym[3];
		v[i].y = M_PI * j / high_sym[3];
		v[i].z = M_PI * j / high_sym[3];
		i++;
	}
}

void FourierF(const int basis, int num, Vector v, Vector q, lapack_complex_double *tb) {
	FILE *f; 

	if((f = fopen("lattice.txt", "r")) == NULL) {
		printf("lattice.txt fopen FAIL\n");
		exit(1);
	}

	int i, j, k, obt1, obt2;
	double tre, tim;
	lapack_complex_double tb_tmp[basis][basis];

	memset(tb_tmp, 0, sizeof(tb_tmp));

	while(!feof(f)) {
		fscanf(f, "%d%d%d%d%d%lf%lf\n", &i, &j, &k, &obt1, &obt2, &tre, &tim);
		tb_tmp[obt1-1][obt2-1] += (tre * cos(v.x*i + v.y*j + v.z*k) - tim * sin(v.x*i + v.y*j + v.z*k))
								+ (tre * sin(v.x*i + v.y*j + v.z*k) + tim * cos(v.x*i + v.y*j + v.z*k)) * I; 
	}
	fclose(f);

	for(i=0; i<OBT; i++) {
		for(j=0; j<OBT; j++) {
			tb_tmp[i + OBT][j + OBT] = tb_tmp[i][j];
		}
	}

	k = 0;
	for(i=0; i<basis; i++) {
		for(j=0; j<basis; j++) {
			tb[H(basis)*num + k] = tb_tmp[i][j];
			k++;
		}
	}
}

void FourierA(const int basis, int num, Vector v, Vector q, lapack_complex_double *tb) {
	FILE *f; 

	if((f = fopen("lattice.txt", "r")) == NULL) {
		printf("lattice.txt fopen FAIL\n");
		exit(1);
	}

	int i, j, k, obt1, obt2;
	double tre, tim;
	lapack_complex_double tb_tmp[basis][basis];

	memset(tb_tmp, 0, sizeof(tb_tmp));

	while(!feof(f)) {
		fscanf(f, "%d%d%d%d%d%lf%lf\n", &i, &j, &k, &obt1, &obt2, &tre, &tim);
		tb_tmp[obt1-1][obt2-1] += (tre * cos(v.x*i + v.y*j + v.z*k) - tim * sin(v.x*i + v.y*j + v.z*k))
								+ (tre * sin(v.x*i + v.y*j + v.z*k) + tim * cos(v.x*i + v.y*j + v.z*k)) * I; 
		tb_tmp[obt1-1 + OBT][obt2-1 + OBT] += (tre * cos((v.x + q.x)*i + (v.y + q.y)*j + (v.z + q.z)*k) - tim * sin((v.x + q.x)*i + (v.y + q.y)*j + (v.z + q.z)))
											+ (tre * sin((v.x + q.x)*i + (v.y + q.y)*j + (v.z + q.z)*k) + tim * cos((v.x + q.x)*i + (v.y + q.y)*j + (v.z + q.z))) * I; 
	}
	fclose(f);

	for(i=0; i<OBT; i++) {
		for(j=0; j<OBT; j++) {
			tb_tmp[i + 2*OBT][j + 2*OBT] = tb_tmp[i][j];
			tb_tmp[i + 3*OBT][j + 3*OBT] = tb_tmp[i + OBT][j + OBT];
		}
	}

	k = 0;
	for(i=0; i<basis; i++) {
		for(j=0; j<basis; j++) {
			tb[H(basis)*num + k] = tb_tmp[i][j];
			k++;
		}
	}

}

void Unfold(char *type, const int super, Vector q) {
	FILE *fi, *fo;
	char fi_name[32], fo_name[32];

	sprintf(fi_name, "tb_%s_BAND.bin", type);
	sprintf(fo_name, "band_%s_unfold.txt", type);

	if((fi = fopen(fi_name, "rb")) == NULL) {
		printf("%s fopen FAIL\n", fi_name);
		exit(1);
	}
	if((fo = fopen(fo_name, "w")) == NULL) {
		printf("%s fopen FAIL\n", fo_name);
		exit(1);
	}

	const int basis = BASIS2;
	int i, j, l, m;
	double p2up, p2dn;
	Vector vb[BAND], g, r[super];
	lapack_complex_double pup, pdn;

	double w[BAND*basis];
	lapack_complex_double tbb[HB(basis)], v[HB(basis)];

	g.x = 0;
	g.y = 0;
	g.z = 0;

	r[0].x = 0;
	r[0].y = 0;
	r[0].z = 0;

	r[1].x = 1;
	r[1].y = 0;
	r[1].z = 0;

	fread(tbb, HB(basis), sizeof(lapack_complex_double), fi); 
	CalcBAND(vb);
	CalcEigenTB(basis, tbb, w, v);

	for(i=0; i<BAND; i++) {
		fprintf(fo, "%4d", i);

		for(j=0; j<basis; j++) {
			p2up = 0;
			p2dn = 0;

			for(l=0; l<OBT; l++) {
				pup = 0;
				pdn = 0;

				for(m=0; m<super; m++) {
					pup += DOT(v[basis*(basis*i + j) + l + OBT*m], vb[i], g, r[m]);
					pdn += DOT(v[basis*(basis*i + j) + l + OBT*(m+2)], vb[i], g, r[m]);
				}
				p2up += COMPLEX2(pup) / super;
				p2dn += COMPLEX2(pdn) / super;
			}
			printf("%f\t%f\t%f\n", w[basis*i + j], p2up, p2dn);
			
			if(p2up > 0.1) {
				fprintf(fo, "%12f", w[basis*i + j] * p2up);
			}
		}
		fprintf(fo, "\n");
		printf("\n");
	}

	fclose(fi);
	fclose(fo);
}

int main(int argc, char *argv[]) {
	if(argc != 3) {
		printf("Usage : %s <type> <is_unfold>\n", argv[0]);
		exit(1);
	}

	char *type = argv[1];
	int is_unfold = atoi(argv[2]), super;
	Vector q;

	if(strstr(type, "f")) {
		q.x = 0;
		q.y = 0;
		q.z = 0;
	}
	else if(strstr(type, "a")) {
		super = 2;
		q.x = M_PI;
		q.y = 0;
		q.z = 0;
	}
	else if(strstr(type, "c")) {
		super = 2;
		q.x = M_PI;
		q.y = M_PI;
		q.z = 0;
	}
	else if(strstr(type, "g")) {
		super = 2;
		q.x = M_PI;
		q.y = M_PI;
		q.z = M_PI;
	}

	if(is_unfold) {
		Unfold(type, super, q);
	}
	else {
		int basis = strstr(type, "f") ? BASIS1 : BASIS2;
		char fk_name[32], fb_name[32];
		Vector vk[K3], vb[BAND];
		lapack_complex_double tbk[HK(basis)], tbb[HB(basis)];


		sprintf(fk_name, "tb_%s_K%d.bin", type, K);
		sprintf(fb_name, "tb_%s_BAND.bin", type);

		CalcK(vk);
		CalcBAND(vb);

		CalcTB(fk_name, basis, K3, HK(basis), vk, q, tbk);
		CalcTB(fb_name, basis, BAND, HB(basis), vb, q, tbb);
		MakeTB(type, basis, tbb);
	}

	return 0;
}
