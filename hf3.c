// hf3.c : universal functions for calculating 3 band models

#include "hf3.h" 

void ReadTB(char *type, const int basis, lapack_complex_double *tbk, lapack_complex_double *tbb) {
	FILE *fk, *fb;
	char fk_name[32], fb_name[32];

	sprintf(fk_name, "input/tb_%s_K%d.bin", type, K);
	sprintf(fb_name, "input/tb_%s_BAND.bin", type);

	if((fk = fopen(fk_name, "rb")) == NULL) {
		printf("%s fopen FAIL\n", fk_name);
		exit(1);
	}
	if((fb = fopen(fb_name, "rb")) == NULL) {
		printf("%s fopen FAIL\n", fb_name);
		exit(1);
	}

	fread(tbk, HK(basis), sizeof(lapack_complex_double), fk);
	fread(tbb, HB(basis), sizeof(lapack_complex_double), fb);

	fclose(fk);
	fclose(fb);
}

Energy CalcEigen(Solution *s, int k_num, lapack_complex_double *tb, double *w, lapack_complex_double *v) {
	Energy e = {
		.min =  100,
		.max = -100
	};

	char jobz = 'V', uplo = 'L';
	double rwork[3*s->basis-2], w_tmp[s->basis];
	lapack_int ln = s->basis, lda = s->basis, lwork = 2*s->basis-1, info;
	lapack_complex_double work[lwork], v_tmp[H(s->basis)];
	
	int i, j;
	void (*Interaction)(Solution*, lapack_complex_double*);

	Interaction = strstr(s->type, "f") ? InteractionF : InteractionA;

	for(i=0; i<k_num; i++) {
		memset(v_tmp, 0, sizeof(v_tmp));

		for(j=0; j<H(s->basis); j++) {
			v_tmp[j] = tb[H(s->basis)*i + j];
		}

		Interaction(s, v_tmp);

		LAPACK_zheev(&jobz, &uplo, &ln, v_tmp, &lda, w_tmp, work, &lwork, rwork, &info);
		if(info != 0) {
			printf("LAPACK_zheev FAIL\n");
			exit(1);
		}

		for(j=0; j<s->basis; j++) {
			w[s->basis*i + j] = w_tmp[j];
		}
		for(j=0; j<H(s->basis); j++) {
			v[H(s->basis)*i + j] = v_tmp[j];
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

void CalcSolution(Solution *s, int is_unfold) {
	FILE *f;
	char f_name[256];

	if(!is_unfold) {
		sprintf(f_name, "output/K%d_JU%.2f_SOC%.2f/cvg_%s_N%.1f_U%.1f_%s.txt", K, s->JU, s->SOC, s->type, s->N, s->U, s->runtime);
	}
	else {
		sprintf(f_name, "output/K%d_JU%.2f_SOC%.2f_unfold/cvg_%s_N%.1f_U%.1f_%s.txt", K, s->JU, s->SOC, s->type, s->N, s->U, s->runtime);
	}

	if((f = fopen(f_name, "w")) == NULL) {
		printf("%s fopen FAIL\n", f_name);
		exit(1);
	}

	time_t t0 = time(NULL);

	double w[K3*s->basis];
	lapack_complex_double v[HK(s->basis)];

	int itr, i, is_cvg;
	double fermi0, fermi1, e, n[OBT], m[OBT], n_total = 0, m_total = 0, cvg[OBT*3];
	void (*Occupation)(int, double, double*, lapack_complex_double*, double*, double*, double*);
	Energy energy;

	Occupation = strstr(s->type, "f") ? OccupationF : OccupationA;

	for(i=0; i<OBT*3; i++) {
		cvg[i] = 100;
	}

	fprintf(f, "%3d%12f%12f", 0, s->fermi, s->e);
	for(i=0; i<OBT; i++) {
		fprintf(f, "%12f%12f", s->n[i], s->m[i]);
	}
	fprintf(f, "%12f%12f\n", s->n_total, s->m_total);

	for(itr=0; itr<50; itr++) {
		energy = CalcEigen(s, K3, s->tbk, w, v);
		fermi0 = energy.min;
		fermi1 = energy.max;
		s->fermi = 0.5 * (fermi0 + fermi1);

		while(fabs(s->fermi - fermi1) > 1e-8) {
			n_total = 0;
			m_total = 0;

			Occupation(s->basis, s->fermi, w, v, n, m, &e);
		
			for(i=0; i<OBT; i++) {	
				n[i] /= K3;
				m[i] /= K3;
				n_total += n[i];
				m_total += m[i];
			}
			s->e = e / K3;

			if(fabs(n_total - s->N) < 1e-6) {
				break;
			}
			else {
				if(n_total < s->N) {
					fermi0 = s->fermi;
				}
				else {
					fermi1 = s->fermi;
				}
				s->fermi = 0.5 * (fermi0 + fermi1);
			}
		}

		for(i=0; i<OBT; i++) {
			s->n[i] = n[i];
			s->m[i] = m[i];
		}
		s->n_total = n_total;
		s->m_total = m_total;

		fprintf(f, "%3d%12f%12f", itr+1, s->fermi, s->e);
		for(i=0; i<OBT; i++) {
			fprintf(f, "%12f%12f", s->n[i], s->m[i]);
		}
		fprintf(f, "%12f%12f\n", s->n_total, s->m_total);
			
		is_cvg = 0;
		for(i=0; i<OBT; i++) {
			cvg[3*(itr%3) + i] = m[i];
			if(fabs((cvg[3*0 + i] + cvg[3*1 + i] + cvg[3*2 + i])/3 - m[i]) < 1e-3) {
				is_cvg++;
			}
		}

		if(is_cvg == 3) {
			break;
		}
	}

	fclose(f);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, f_name, t1 - t0);
}

void MakeBand(Solution *s) {
	FILE *f;
	char f_name[256];

	sprintf(f_name, "output/K%d_JU%.2f_SOC%.2f/band_%s_N%.1f_U%.1f_n%f_m%f_fermi%f_e%f_%s.txt", K, s->JU, s->SOC, s->type, s->N, s->U, s->n_total, s->m_total, s->fermi, s->e, s->runtime);

	if((f = fopen(f_name, "w")) == NULL) {
		printf("%s fopen FAIL\n", f_name);
		exit(1);
	}

	time_t t0 = time(NULL);

	double w[BAND*s->basis];
	lapack_complex_double v[HB(s->basis)];

	int i, j;

	CalcEigen(s, BAND, s->tbb, w, v);
	
	for(i=0; i<BAND; i++) {
		fprintf(f, "%4d", i);
		for(j=0; j<s->basis; j++) {
			fprintf(f, "%12f", w[s->basis*i + j]);
		}
		fprintf(f, "\n");
	}

	fclose(f);

	time_t t1 = time(NULL);
	printf("%s(%s) : %lds\n", __func__, f_name, t1 - t0);
}

void MakeDos(Solution *s) {
	FILE *f;
	char f_name[256];

	sprintf(f_name, "output/K%d_JU%.2f_SOC%.2f/dos_%s_N%.1f_U%.1f_n%f_m%f_fermi%f_e%f_%s.txt", K, s->JU, s->SOC, s->type, s->N, s->U, s->n_total, s->m_total, s->fermi, s->e, s->runtime);

	if((f = fopen(f_name, "w")) == NULL) {
		printf("%s fopen FAIL\n", f_name);
		exit(1);
	}

	time_t t0 = time(NULL);

	double w[K3*s->basis];
	lapack_complex_double v[HK(s->basis)];

	int itv, i, j;
	double energy, dos[s->basis];
	Energy e = CalcEigen(s, K3, s->tbk, w, v);

	e.min -= 0.1;
	e.max += 0.1;
	
	for(itv=0; itv<=100; itv++) {
		for(i=0; i<s->basis; i++) {
			dos[i] = 0;
		}
		energy = e.min + (e.max - e.min) * itv / 100;

		for(i=0; i<K3*s->basis; i++) {
			for(j=0; j<s->basis; j++) {
				dos[j] += GREEN(energy, w[i]) * COMPLEX2(v[s->basis*i + j]);
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
	printf("%s(%s) : %lds\n", __func__, f_name, t1 - t0);
}
