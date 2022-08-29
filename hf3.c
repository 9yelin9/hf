// hf3.c : universal functions for calculating hf3 model

#define GREEN(i) (0.05 / (pow(energy - w[i], 2) + pow(0.05, 2))) 

#include "hf3.h" 

void CalcVK() {
	FILE *f;
	char fs[64];

	sprintf(fs, "input/vk%d.bin", K);

	if((f = fopen(fs, "wb")) == NULL) {
		printf("%s fopen FAIL\n", fs);
		exit(1);
	}

	int i;
	Vector vk[K3];

	for(i=0; i<K3; i++) {
		vk[i].k[0] = -M_PI + 2*M_PI * (i / (K*K)) / K;
		vk[i].k[1] = -M_PI + 2*M_PI * ((i / K) % K) / K;
		vk[i].k[2] = -M_PI + 2*M_PI * (i % K) / K;
	}

	fwrite(vk, sizeof(vk), 1, f);
	fclose(f);
}

int* ReadPath() {
	FILE *f;
	char fs[64];

	sprintf(fs, "input/info.txt");

	if((f = fopen(fs, "r")) == NULL) {
		printf("%s fopen FAIL\n", fs);
		exit(1);
	}

	int i, tmp = 0, sym_len, *sym;
	char buf[1024];

	fscanf(f, "%s", buf);
	sym_len = atoi(buf);
	sym = (int*)malloc(sym_len * sizeof(int));

	for(i=0; i<sym_len; i++) {
		fscanf(f, "%s", buf);	
		sym[i] = atoi(buf) - tmp;
		tmp = atoi(buf);
	}
	fclose(f);

	return sym;
}

Lattice* ReadLattice(int *l_len) {
	FILE *f;
	char fs[64];

	sprintf(fs, "input/lattice.txt");

	if((f = fopen(fs, "r")) == NULL) {
		printf("%s fopen FAIL\n", fs);
		exit(1);
	}

	int i;
	char buf[1024];

	// get length of lattice.txt
	while(!feof(f)) {
		fgets(buf, sizeof(buf), f);
		(*l_len)++;
	}
	rewind(f);

	*l_len -= 1;

	// read lattice.txt
	Lattice *l = (Lattice*)malloc(*l_len * sizeof(Lattice));
	for(i=0; i<*l_len; i++) {
		fscanf(f, "%d%d%d%d%d%lf%lf", &l[i].r[0], &l[i].r[1], &l[i].r[2], &l[i].obt1, &l[i].obt2, &l[i].tre, &l[i].tim);
	}
	fclose(f);

	return l;
}

void ReadInfo(char *type, char *cell_type, Vector *vq) {
	FILE *f;
	char fs[64];

	sprintf(fs, "input/info.txt");

	if((f = fopen(fs, "r")) == NULL) {
		printf("%s fopen FAIL\n", fs);
		exit(1);
	}

	int i;
	char buf[1024], tmp[8];

	while(!feof(f)) {
		fgets(buf, sizeof(buf), f);
		if(strstr(buf, type)) {
			sscanf(buf, "%s%s%lf%lf%lf", tmp, cell_type, &vq->k[0], &vq->k[1], &vq->k[2]);
			break;
		}
	}
	fclose(f);

	for(i=0; i<3; i++) vq->k[i] *= M_PI;
}

void ReadVector(Vector *vk, Vector *vb) {
	FILE *fk, *fb;
	char fks[64], fbs[64];

	sprintf(fks, "input/vk%d.bin", K);
	sprintf(fbs, "input/vb.bin");

	if((fk = fopen(fks, "rb")) == NULL) {
		printf("%s fopen FAIL\n", fks);
		exit(1);
	}
	if((fb = fopen(fbs, "rb")) == NULL) {
		printf("%s fopen FAIL\n", fbs);
		exit(1);
	}

	fread(vk, K3,   sizeof(Vector), fk);
	fread(vb, BAND, sizeof(Vector), fb);

	fclose(fk);
	fclose(fb);
}

void CalcTB(char *type) {
	FILE *fk, *fb, *f;
	char fks[64], fbs[64], fs[64];

	sprintf(fks, "input/hk%d_%s.bin", K, type);
	sprintf(fbs, "input/hb_%s.bin", type);
	sprintf(fs, "input/band_%s.txt", type);

	if((fk = fopen(fks, "wb")) == NULL) {
		printf("%s fopen FAIL\n", fks);
		exit(1);
	}
	if((fb = fopen(fbs, "wb")) == NULL) {
		printf("%s fopen FAIL\n", fbs);
		exit(1);
	}
	if((f = fopen(fs, "w")) == NULL) {
		printf("%s fopen FAIL\n", fs);
		exit(1);
	}

	time_t t0 = time(NULL);

	int i, j, l_len = 0;
	char cell_type[8];
	Vector vq, vk[K3], vb[BAND];
	void (*Fourier)(int, int, Vector, Vector, Lattice*, lapack_complex_double*);

	Lattice *l = ReadLattice(&l_len);
	ReadInfo(type, cell_type, &vq);
	ReadVector(vk, vb);

	const int basis = strstr(cell_type, "s") ? SINGLE : DOUBLE;
	lapack_complex_double hk[K3 * basis*basis], hb[BAND * basis*basis];

	if(strstr(cell_type, "s")) Fourier = FourierS;
	else Fourier = FourierD;

	for(i=0; i<K3;   i++) Fourier(l_len, i, vk[i], vq, l, hk);
	for(i=0; i<BAND; i++) Fourier(l_len, i, vb[i], vq, l, hb);

	fwrite(hk, sizeof(hk), 1, fk);
	fwrite(hb, sizeof(hb), 1, fb);

	fclose(fk);
	fclose(fb);

	time_t t1 = time(NULL);
	printf("%s(%s, %s) : %lds\n", __func__, fks, fbs, t1 - t0);

	// make band
	char jobz = 'V', uplo = 'L';
	double rwork[3*basis-2], w[basis];
	lapack_int ln = basis, lda = basis, lwork = 2*basis-1, info;
	lapack_complex_double work[lwork], v[basis*basis];

	for(i=0; i<BAND; i++) {
		for(j=0; j<basis*basis; j++) v[j] = hb[basis*basis*i + j];

		LAPACK_zheev(&jobz, &uplo, &ln, v, &lda, w, work, &lwork, rwork, &info);
		if(info != 0) {
			printf("LAPACK_zheev FAIL\n");
			exit(1);
		}

		fprintf(f, "%4d", i);
		for(j=0; j<basis; j++) fprintf(f, "%12f", w[j]);
		fprintf(f, "\n");
	}
	fclose(f);
}

void ReadTB(Solution *s) {
	FILE *fk, *fb;
	char fks[32], fbs[32];

	sprintf(fks, "input/hk%d_%s.bin", K, s->type);
	sprintf(fbs, "input/hb_%s.bin", s->type);

	if((fk = fopen(fks, "rb")) == NULL) {
		printf("%s fopen FAIL\n", fks);
		exit(1);
	}
	if((fb = fopen(fbs, "rb")) == NULL) {
		printf("%s fopen FAIL\n", fbs);
		exit(1);
	}

	fread(s->hk, K3   * s->basis*s->basis, sizeof(lapack_complex_double), fk);
	fread(s->hb, BAND * s->basis*s->basis, sizeof(lapack_complex_double), fb);

	fclose(fk);
	fclose(fb);
}

void GetName(Solution *s, char *data_type, char *fs) {
	sprintf(fs,\
		   	"output/K%d_JU%.2f_SOC%.2f/%s_%s_N%.1f_U%.1f_n%f_m%f_fermi%f_e%f_%s.txt",\
		   	K, s->JU, s->SOC, data_type, s->type, s->N, s->U, s->ntot, s->mtot, s->fermi, s->e, s->runtime);
}

Energy CalcEigen(Solution *s, int h_len, lapack_complex_double *h, double *w, lapack_complex_double *v) {
	Energy e = {
		.min =  100,
		.max = -100
	};

	char jobz = 'V', uplo = 'L';
	double rwork[3*s->basis-2], w_tmp[s->basis];
	lapack_int ln = s->basis, lda = s->basis, lwork = 2*s->basis-1, info;
	lapack_complex_double work[lwork], v_tmp[s->basis*s->basis];
	
	int i, j;
	void (*Interaction)(Solution*, lapack_complex_double*);

	if(strstr(s->cell_type, "s")) Interaction = InteractionS;
	else Interaction = InteractionD;

	for(i=0; i<h_len; i++) {
		//memset(v_tmp, 0, sizeof(v_tmp));
		for(j=0; j<s->basis*s->basis; j++) v_tmp[j] = h[s->basis*s->basis*i + j];

		Interaction(s, v_tmp);

		LAPACK_zheev(&jobz, &uplo, &ln, v_tmp, &lda, w_tmp, work, &lwork, rwork, &info);
		if(info != 0) {
			printf("LAPACK_zheev FAIL\n");
			exit(1);
		}

		for(j=0; j<s->basis; j++) {
			w[s->basis*i + j] = w_tmp[j];
		}
		for(j=0; j<s->basis*s->basis; j++) {
			v[s->basis*s->basis*i + j] = v_tmp[j];
		}

		if(w_tmp[0] < e.min) {
			e.min = w_tmp[0];
		}
		if(w_tmp[s->basis-1] > e.max) {
			e.max = w_tmp[s->basis-1];
		}
	}

	return e;
}

void CalcSolution(Solution *s) {
	FILE *f;
	char fs[256];
	GetName(s, "cvg", fs);

	if((f = fopen(fs, "w")) == NULL) {
		printf("%s fopen FAIL\n", fs);
		exit(1);
	}

	time_t t0 = time(NULL);

	double w[K3 * s->basis];
	lapack_complex_double v[K3 * s->basis*s->basis];

	int itr, i, is_cvg;
	double fermi0, fermi1, e, n[OBT], m[OBT], ntot = 0, mtot = 0, cvg[OBT*3];
	void (*Occupation)(double, double*, lapack_complex_double*, double*, double*, double*);
	Energy energy;

	if(strstr(s->cell_type, "s")) Occupation = OccupationS;
	else Occupation = OccupationD;

	for(i=0; i<OBT*3; i++) cvg[i] = 100;

	fprintf(f, "%5s%12s%12s", "#itr", "fermi", "e");
	for(i=0; i<OBT; i++) fprintf(f, "%10s%02d%10s%02d", "n", i+1, "m", i+1);
	fprintf(f, "%12s%12s\n", "ntotal", "mtotal");

	fprintf(f, "%5d%12f%12f", 0, s->fermi, s->e);
	for(i=0; i<OBT; i++) fprintf(f, "%12f%12f", s->n[i], s->m[i]);
	fprintf(f, "%12f%12f\n", s->ntot, s->mtot);

	for(itr=0; itr<50; itr++) {
		energy = CalcEigen(s, K3, s->hk, w, v);
		fermi0 = energy.min;
		fermi1 = energy.max;
		s->fermi = 0.5 * (fermi0 + fermi1);

		while(fabs(s->fermi - fermi1) > 1e-8) {
			ntot = 0;
			mtot = 0;

			Occupation(s->fermi, w, v, n, m, &e);
		
			for(i=0; i<OBT; i++) {	
				n[i] /= K3;
				m[i] /= K3;
				ntot += n[i];
				mtot += m[i];
			}
			s->e = e / K3;

			if(fabs(ntot - s->N) < 1e-6) break;
			else {
				if(ntot < s->N) fermi0 = s->fermi;
				else fermi1 = s->fermi;
				s->fermi = 0.5 * (fermi0 + fermi1);
			}
		}

		for(i=0; i<OBT; i++) {
			s->n[i] = n[i];
			s->m[i] = m[i];
		}
		s->ntot = ntot;
		s->mtot = mtot;

		fprintf(f, "%5d%12f%12f", itr+1, s->fermi, s->e);
		for(i=0; i<OBT; i++) fprintf(f, "%12f%12f", s->n[i], s->m[i]);
		fprintf(f, "%12f%12f\n", s->ntot, s->mtot);
			
		is_cvg = 0;
		for(i=0; i<OBT; i++) {
			cvg[3*(itr%3) + i] = m[i];
			if(fabs((cvg[3*0 + i] + cvg[3*1 + i] + cvg[3*2 + i])/3 - m[i]) < 1e-3) is_cvg++;
		}

		if(is_cvg == 3) break;
	}

	fclose(f);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fs, t1 - t0);
}

void MakeBand(Solution *s) {
	FILE *f;
	char fs[256];
	GetName(s, "band", fs);

	if((f = fopen(fs, "w")) == NULL) {
		printf("%s fopen FAIL\n", fs);
		exit(1);
	}

	time_t t0 = time(NULL);

	int i, j;
	double w[BAND * s->basis];
	lapack_complex_double v[BAND * s->basis*s->basis];

	CalcEigen(s, BAND, s->hb, w, v);

	fprintf(f, "%4s", "#");
	for(i=0; i<s->basis; i++) fprintf(f, "%10s%02d", "e", i+1);
	fprintf(f, "\n");

	for(i=0; i<BAND; i++) {
		fprintf(f, "%4d", i);
		for(j=0; j<s->basis; j++) fprintf(f, "%12f", w[s->basis*i + j]);
		fprintf(f, "\n");
	}

	fclose(f);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fs, t1 - t0);

	if(!strstr(s->cell_type, "s")) {
		GetName(s, "ufw", fs);

		if((f = fopen(fs, "w")) == NULL) {
			printf("%s fopen FAIL\n", fs);
			exit(1);
		}

		t0 = time(NULL);

		int k;
		double dot, p2;
		Vector r[SUPER] = {{{1, 0, 0}}, {{0, 0, 0}}};	
		lapack_complex_double exp[SUPER], p; 

		fprintf(f, "%4s", "#");
		for(i=0; i<s->basis; i++) fprintf(f, "%10s%02d", "w", i+1);
		fprintf(f, "\n");
		
		for(i=0; i<BAND; i++) {
			fprintf(f, "%4d", i);

			for(j=0; j<SUPER; j++) {
				dot = 0;
				for(k=0; k<3; k++) dot += s->vb[i].k[k] * r[j].k[k];
				exp[j] = cos(dot) - sin(dot) * I;
			}

			for(j=0; j<DOUBLE; j++) {
				p2 = 0;

				for(k=0; k<SINGLE; k++) {
					p = v[DOUBLE*DOUBLE*i + DOUBLE*j + OBT*(k/OBT) + k] * exp[0] + v[DOUBLE*DOUBLE*i + DOUBLE*j + OBT*(k/OBT) + k+OBT] * exp[1];	
					p2 += CSQR(p) / SUPER;
				}
				
				fprintf(f, "%12f", p2); 
			}
			fprintf(f, "\n");
		}

		fclose(f);

		t1 = time(NULL);
		printf("%s(%s) : %lds\n", __func__, fs, t1 - t0);
	}
}

void MakeDos(Solution *s) {
	FILE *f;
	char fs[256];
	GetName(s, "dos", fs);

	if((f = fopen(fs, "w")) == NULL) {
		printf("%s fopen FAIL\n", fs);
		exit(1);
	}

	time_t t0 = time(NULL);

	double w[K3 * s->basis];
	lapack_complex_double v[K3 * s->basis*s->basis];

	int itv, i, j;
	double energy, dos[s->basis];
	Energy e = CalcEigen(s, K3, s->hk, w, v);

	e.min -= 0.1;
	e.max += 0.1;
	
	fprintf(f, "%12s", "#e");
	for(i=0; i<s->basis; i++) fprintf(f, "%10s%02d", "d", i+1);
	fprintf(f, "\n");

	for(itv=0; itv<512; itv++) {
		for(i=0; i<s->basis; i++) {
			dos[i] = 0;
		}
		energy = e.min + (e.max - e.min) * itv / 512;

		for(i=0; i<K3*s->basis; i++) {
			for(j=0; j<s->basis; j++) {
				dos[j] += M_PI * CSQR(v[s->basis*i + j]) * GREEN(i);
			}
		}

		fprintf(f, "%12f", energy);
		for(i=0; i<s->basis; i++) {
			fprintf(f, "%12f", dos[i] / K3);
		}
		fprintf(f, "\n");
	}

	fclose(f);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, fs, t1 - t0);
}

int main(int argc, char *argv[]) {
	if(argc < 2) {
		printf("%s v : make reciprocal lattice vectors\n%s tb <type> : make tight-binding matrix\n%s <type> <J/U> <SOC> <N> <U> : make Hartree-Fock approximated 3-band Hubbard model\n", argv[0], argv[0], argv[0]);
		exit(1);
	}

	if(strstr(argv[1], "v")) {
		CalcVK();
		CalcVB();
		return 0;
	}
	else if(strstr(argv[1], "tb")) {
		CalcTB(argv[2]);
		return 0;
	}
	else {
		Solution s = {
			.type = argv[1],
			.JU = atof(argv[2]),
			.SOC = atof(argv[3]),
			.N = atof(argv[4]),
			.U = atof(argv[5]),
			.J = s.JU * s.U,

			.n = {s.N/3, s.N/3, s.N/3},
			.m = {0.1, 0.1, 0.1},
			.ntot = 100,
			.mtot = 100,
			.fermi = 100,
			.e = 100
		};

		ReadInfo(s.type, s.cell_type, &s.vq);
		ReadVector(s.vk, s.vb);
		s.basis = strstr(s.cell_type, "s") ? SINGLE : DOUBLE;

		time_t t = time(NULL);
		struct tm *tm = localtime(&t);
		sprintf(s.runtime, "%d%d%d%d", tm->tm_mon+1, tm->tm_mday, tm->tm_hour, tm->tm_min); 

		s.hk = (lapack_complex_double*)calloc(K3   * s.basis*s.basis, sizeof(lapack_complex_double));
		s.hb = (lapack_complex_double*)calloc(BAND * s.basis*s.basis, sizeof(lapack_complex_double));
		ReadTB(&s);

		CalcSolution(&s);

		MakeBand(&s);
		MakeDos(&s);

		free(s.hk);
		free(s.hb);

		return 0;
	}

	printf("nothing happened\n");
	exit(1);
}
