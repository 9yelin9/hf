// boo/boo.h

#ifndef BOO_H
#define BOO_H

const int K = 16;
const int K3 = 16*16*16;
const int BAND = 510;

void Interaction(Solution *s, lapack_complex_double *v_tmp); 
void FourierF(char *type, int num, Vector v, Vector q, lapack_complex_double *tb);
void FourierA(char *type, int num, Vector v, Vector q, lapack_complex_double *tb);
void OccupationF(double fermi, double *w, lapack_complex_double *v, double *n, double *m, double *e);
void OccupationA(double fermi, double *w, lapack_complex_double *v, double *n, double *m, double *e);

#endif
