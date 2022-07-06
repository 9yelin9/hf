// boo/boo.h

#ifndef BOO_H
#define BOO_H

#define OBT_IDX(i) ((i / s->basis) % OBT)
#define N_TERM(i) (s->U * s->n[OBT_IDX(i)] + (s->U - 2.5*s->J) * n_diff[OBT_IDX(i)])
#define M_TERM(i) (s->U * s->m[OBT_IDX(i)] + (s->U - 2.5*s->J) * m_diff[OBT_IDX(i)])

const int K = 16;
const int K3 = 16*16*16;
const int BAND = 510;

void Interaction(Solution *s, lapack_complex_double *v_tmp); 
void FourierF(const int basis, int num, Vector v, Vector q, lapack_complex_double *tb);
void FourierA(const int basis, int num, Vector v, Vector q, lapack_complex_double *tb);
void OccupationF(const int basis, double fermi, double *w, lapack_complex_double *v, double *n, double *m, double *e);
void OccupationA(const int basis, double fermi, double *w, lapack_complex_double *v, double *n, double *m, double *e);

#endif
