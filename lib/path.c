// lib/path.c : functions for calculating band path points

#include "hf3.h" 

void CalcPathBaOsO3(int BAND) {
	FILE *f = OpenFile("input/baoso3/gb.bin", "wb");

	int i = 0, j;
	Vector gb[BAND];

	const int p_len = GetLen("input/baoso3/info.txt");
	int p[p_len];
	ReadPathInfo("baoso3", p_len, p);

	for(j=0; j<p[1]; j++) { // G-X
		gb[i].c[0] = M_PI * j / p[1];
		gb[i].c[1] = 0;
		gb[i].c[2] = 0;
		i++;
	}

	for(j=0; j<p[2]; j++) { // X-M
		gb[i].c[0] = M_PI;
		gb[i].c[1] = M_PI * j / p[2];
		gb[i].c[2] = 0;
		i++;
	}

	for(j=0; j<p[3]; j++) { // M-G
		gb[i].c[0] = M_PI - M_PI * j / p[3];
		gb[i].c[1] = M_PI - M_PI * j / p[3];
		gb[i].c[2] = 0;
		i++;
	}

	for(j=0; j<p[4]; j++) { // G-R
		gb[i].c[0] = M_PI * j / p[4];
		gb[i].c[1] = M_PI * j / p[4];
		gb[i].c[2] = M_PI * j / p[4];
		i++;
	}

	fwrite(gb, sizeof(gb), 1, f);
	fclose(f);
}

void CalcPathCuAl2O4(int BAND) {
	FILE *f = OpenFile("input/cual2o4/gb.bin", "wb");

	int i = 0, j;
	Vector gb[BAND];

	const int p_len = GetLen("input/cual2o4/info.txt");
	int p[p_len];
	ReadPathInfo("cual2o4", p_len, p);

	for(j=0; j<p[0]; j++) { // L-G
		gb[i].c[0] = M_PI - M_PI * j / p[0];
		gb[i].c[1] = M_PI - M_PI * j / p[0];
		gb[i].c[2] = M_PI - M_PI * j / p[0];
		i++;
	}

	for(j=0; j<p[1]; j++) { // G-X
		gb[i].c[0] = M_PI * j / p[1];
		gb[i].c[1] = 0;
		gb[i].c[2] = M_PI * j / p[1];
		i++;
	}

	for(j=0; j<p[2]; j++) { // X-W
		gb[i].c[0] = M_PI;
		gb[i].c[1] = (M_PI/2) * j / p[2];
		gb[i].c[2] = M_PI + (M_PI/2) * j / p[2];
		i++;
	}

	for(j=0; j<p[3]; j++) { // W-L
		gb[i].c[0] = M_PI;
		gb[i].c[1] = (M_PI/2) + (M_PI/2) * j / p[3];
		gb[i].c[2] = (3*M_PI/2) - (M_PI/2) * j / p[3];
		i++;
	}

	for(j=0; j<p[4]; j++) { // L-K
		gb[i].c[0] = M_PI - (M_PI/4) * j / p[4];
		gb[i].c[1] = M_PI + (M_PI/2) * j / p[4];
		gb[i].c[2] = M_PI - (M_PI/4) * j / p[4];
		i++;
	}

	for(j=0; j<p[5]; j++) { // K-G
		gb[i].c[0] = (3*M_PI/4) - (3*M_PI/4) * j / p[5];
		gb[i].c[1] = (3*M_PI/2) - (3*M_PI/2) * j / p[5];
		gb[i].c[2] = (3*M_PI/4) - (3*M_PI/4) * j / p[5];
		i++;
	}

	fwrite(gb, sizeof(gb), 1, f); 
	fclose(f);
}
