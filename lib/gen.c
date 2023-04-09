// lib/gen.c : functions for generating I/O files

void GenGaussQuad(Config c) {
	FILE *fk, *fw;
	char fkn[64], fwn[64];

	sprintf(fkn, "input/%s/gaussquad_k.bin", c.name);
	sprintf(fwn, "input/%s/gaussquad_w.bin", c.name);
	fk = OpenFile(fkn, "wb");
	fw = OpenFile(fwn, "wb");

	time_t t0 = time(NULL);

	int i, i0, i1, i2;
	double k0[Nkg1], w0[Nkg1], w[Nkg];
	Coord k[Nkg];
	gsl_integration_glfixed_table *t = gsl_integration_glfixed_table_alloc(Nkg1);

	for(i=0; i<Nkg1; i++) gsl_integration_glfixed_point(-M_PI, M_PI, i, &k0[i], &w0[i], t);

	for(i=0; i<Nkg; i++) {
		i0 = i / (Nkg1 * Nkg1);
		i1 =(i / Nkg1) % Nkg1;
		i2 = i % Nkg1;

		k[i].c[0] = k0[i0];
		k[i].c[1] = k0[i1];
		k[i].c[2] = k0[i2];
		w[i] = w0[i0] * w0[i1] * w0[i2];
	}

	fwrite(k, sizeof(k), 1, fk);
	fwrite(w, sizeof(w), 1, fw);

	fclose(fk);
	fclose(fw);

	time_t t1 = time(NULL);
	printf("%s(%s, %s) : %lds\n", __func__, fkn, fwn, t1 - t0);
}

void DotProd(Config c, int Nk, Coord *k, Coord *r, lapack_complex_double *cf) {
	int i, j, l;
	double dot;

	for(i=0; i<Nk; i++) {
		for(j=0; j<c.Ni; j++) {
			dot = 0;
			for(l=0; l<DIM; l++) dot += r[j].c[l] * k[i].c[l];
			cf[c.Ni*i + j] = cos(dot) - sin(dot) * I;
		}
	}
}

void CalcCoef(Config c, Coord *r) {
	FILE *fg, *fb;
	char fgn[64], fbn[64];

	sprintf(fgn, "input/%s/cfg.bin", c.name);
	sprintf(fbn, "input/%s/cfb.bin", c.name);
	fg = OpenFile(fgn, "wb");
	fb = OpenFile(fbn, "wb");

	time_t t0 = time(NULL);

	char kgn[64], kbn[64];
	Coord kg[Nkg], kb[Nkb];
	lapack_complex_double cfg[c.Ni * Nkg], cfb[c.Ni * Nkb];

	sprintf(kgn, "input/%s/kg.bin", c.name);
	sprintf(kbn, "input/%s/kb.bin", c.name);

	ReadBin(kgn, sizeof(Coord) * Nkg, kg);
	ReadBin(kbn, sizeof(Coord) * Nkb, kb);

	DotProd(c, Nkg, kg, r, cfg);
	DotProd(c, Nkb, kb, r, cfb);

	fwrite(cfg, sizeof(cfg), 1, fg);
	fwrite(cfb, sizeof(cfb), 1, fb);

	fclose(fg);
	fclose(fb);

	time_t t1 = time(NULL);
	printf("%s(%s, %s) : %lds\n", __func__, fgn, fbn, t1 - t0);
}

