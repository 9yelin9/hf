// hf3.h : headers for Hartree-Fock approximated 3-band Hubbard model

#ifndef HF3_H
#define HF3_H

#define USE_MATH_DEFINES
#define OMP_THREAD (16)

#define GAUSS3 (GAUSS*GAUSS*GAUSS)
#define OBT (3) // num of orbital
#define SUPER (2) // num of supercell in primitive cell (in reciprocal space)
#define SINGLE (6) // num of single cell basis
#define DOUBLE (2*SINGLE) // num of double cell basis

#define GAUSS_START ((GAUSS * (GAUSS + 1)) / 2)
#define OBT_IDX ((i / s->basis) % OBT)
#define STATE_IDX (DOUBLE*DOUBLE*i + DOUBLE*j + OBT*(k/OBT))

#define INTER_N (0.5 * ((s->U) * n[OBT_IDX] + (s->U - 2*s->J) * n_[OBT_IDX] + (s->U - 3*s->J) * n_[OBT_IDX]))
#define INTER_M (0.5 * ((s->U) * m[OBT_IDX] + (s->U - 2*s->J) * m_[OBT_IDX] - (s->U - 3*s->J) * m_[OBT_IDX])) 
#define CSQR(c) (pow(creal(c), 2) + pow(cimag(c), 2))
#define GREEN(i) (0.05 / (pow(e - ev[i], 2) + pow(0.05, 2))) 

#include <omp.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lapacke.h>
#include <gsl/gsl_integration.h>

typedef struct Vector {
	double c[3];
} Vector;

typedef struct Lattice {
	int c[3];
	int obt1;
	int obt2;
	double tre;
	double tim;
} Lattice;

typedef struct LAPACKParameter {
	char jobz;
	char uplo;
	double *rwork;
	lapack_int ln;
	lapack_int lda;
	lapack_int lwork;
	lapack_int info;
	lapack_complex_double *work;
} LParameter;

typedef struct SelfConsistentSolution {
	// input
	int basis; // num of bases
	char runtime[16]; // runtime
	char ctype[8]; // type of unit cell

	char *type; // type of magnetization
	double JU; // Hund coupling / Coulomb interaction
	double SOC; // spin-orbit coupling
	double N; // occupation
	double U; // Coulomb interaction
	double J; // Hund coupling
	Vector q; // order of magnetization

	// output
	double n[OBT]; // occupation per orbital
	double m[OBT]; // magnetization per orbital
	double ntot; // total occupation
	double mtot; // total magnetization
	double fermi; // Fermi level	
	double e; // energy
} Solution;

typedef struct Energy {
	double min;
	double max;
} Energy;

FILE* OpenFile(char *fs, char *ftype); // open file
int GetLen(char *fs); // get length of file
void ShowBlock(int basis, lapack_complex_double *block); // show block
void GenName(Solution *s, char *data_type, char *fs); // genenrate file name
void DotProd(int g_len, Vector *g, lapack_complex_double *coef); // dot product
void TransBasis(lapack_complex_double *coef, lapack_complex_double *es); // Transform basis (Q to sublattice) of matrices

void CalcGauss(); // calculate points and weights for Gauss-Legendre quadrature
void CalcBandPath(); // calculate reciprocal lattice vectors at band path
void CalcCoef(); // calculate basis transform coefficients
void CalcTB(char *type); // calculate tight-binding Hamiltonian matrices
void CalcEigen(Solution *s, LParameter *lp, int tb_len, lapack_complex_double *tb, double *ev, lapack_complex_double *es, Energy *e); // calculate eigenproblems 
void CalcSolution(Solution *s, LParameter *lp, lapack_complex_double *tbg); // calculate self-consistent solution

void ReadPathInfo(int sym_len, int *sym); // read path in info.txt
void ReadLattice(int l_len, Lattice *l); // read lattice in lattice.txt
void ReadInfo(char *type, char *ctype, Vector *q); // read info in info.txt
void ReadBin(char *fs, int bin_size, void *bin); // read binary files

void MakeBand(Solution *s, LParameter *lp, lapack_complex_double *tbb); // make band structure data
void MakeUFW(Solution *s, double *ev, lapack_complex_double *es); // make unfolding weight data
void MakeDOS(Solution *s, LParameter *lp, lapack_complex_double *tbg); // make density of states data

void FourierS(int l_len, int tb_num, Lattice *l, Vector g, Vector q, lapack_complex_double *tb); // Fourier transform of single cell
void FourierD(int l_len, int tb_num, Lattice *l, Vector g, Vector q, lapack_complex_double *tb); // Fourier transform of double cell

void InteractionS(Solution *s, lapack_complex_double *tb_block); // add interaction term to Hamiltonian of single cell
void InteractionD(Solution *s, lapack_complex_double *tb_block); // add interaction term to Hamiltonian of double cell

void GaussQuadS(double fermi, double *wg, double *ev, lapack_complex_double *es, double *n, double *m, double *e); // calculate occupation of single cell
void GaussQuadD(double fermi, double *wg, double *ev, lapack_complex_double *es, double *n, double *m, double *e); // calcultae occupation of double cell

#endif
