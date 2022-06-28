// cao/cao.h

#ifndef CAO_H
#define CAO_H

const int K = 16;
const int K3 = 16*16*16;
const int BAND = 800;

void Interaction(Solution *s, lapack_complex_double *v_tmp); 
void FourierF(FILE *f, const int num, const int basis, Vector v, Vector q, lapack_complex_double *tb);
void FourierA(FILE *f, const int num, const int basis, Vector v, Vector q, lapack_complex_double *tb);
void OccupationF(int basis, double fermi, double *w, lapack_complex_double *v, double *n, double *m, double *e);
void OccupationA(int basis, double fermi, double *w, lapack_complex_double *v, double *n, double *m, double *e);

#endif
