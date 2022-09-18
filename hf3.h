// hf3.h : headers for Hartree-Fock approximated 3-band Hubbard model

#ifndef HF3_H
#define HF3_H

#define USE_MATH_DEFINES
#define THREAD_NUM (16)

#define K (16) // num of k-points
#define K3 (K*K*K) // K^3
#define OBT (3) // num of orbital
#define SUPER (2) // num of supercell in primitive cell (in k-space)
#define SINGLE (6) // num of single cell basis
#define DOUBLE (2*SINGLE) // num of double cell basis

#define OBT_IDX ((i / s->basis) % OBT)
#define STATE_IDX (DOUBLE*DOUBLE*i + DOUBLE*j + OBT*(k/OBT))

#define INTER_N (0.5 * ((s->U) * n[OBT_IDX] + (s->U - 2*s->J) * n_[OBT_IDX] + (s->U - 3*s->J) * n_[OBT_IDX]))
#define INTER_M (0.5 * ((s->U) * m[OBT_IDX] + (s->U - 2*s->J) * m_[OBT_IDX] - (s->U - 3*s->J) * m_[OBT_IDX]))

#define CSQR(x) (pow(creal(x), 2) + pow(cimag(x), 2))

#include <omp.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lapacke.h>

typedef struct Lattice {
	int r[3];
	int obt1;
	int obt2;
	double tre;
	double tim;
} Lattice;

typedef struct Vector {
	double k[3];
} Vector;

typedef struct Energy {
	double min;
	double max;
} Energy;

typedef struct SelfConsistentSolution {
	// input
	int basis; // num of bases
	char runtime[16]; // runtime
	char cell_type[8]; // type of unit cell

	char *type; // type of magnetization
	double JU; // Hund coupling / Coulomb interaction
	double SOC; // spin-orbit coupling
	double N; // occupation
	double U; // Coulomb interaction
	double J; // Hund coupling

	Vector vq; // order of magnetization
	Vector vk[K3]; // vectors in k-space
	Vector vb[BAND]; // vectors in band path
	lapack_complex_double ek[2 * K3]; // basis transform coefficients in k-space
	lapack_complex_double eb[2 * BAND]; // basis transform coefficients in band path
	lapack_complex_double *hk; // tight-binding Hamiltonian matrices in k-space
	lapack_complex_double *hb; // tight-binding Hamiltonian matrices in band path

	// output
	double n[OBT]; // occupation per orbital
	double m[OBT]; // magnetization per orbital
	double ntot; // total occupation
	double mtot; // total magnetization
	double fermi; // Fermi level	
	double e; // energy
} Solution;

void CalcVK(); // calculate vectors in k-space
void CalcVB(); // calculate vectors in band path
void CalcE(); // calculate basis transform coefficients in k-space and band path
void CalcTB(char *type); // calculate tight-binding Hamiltonian matrices
Energy CalcEigen(Solution *s, int h_len, lapack_complex_double *h, double *w, lapack_complex_double *v); // calculate eigenproblems 
void CalcSolution(Solution *s); // calculate self-consistent solution

int* ReadPath(); // read path in info.txt
Lattice* ReadLattice(int *l_len); // read lattice in lattice.txt
void ReadInfo(char *type, char *cell, Vector *vq); // read info in info.txt
void ReadVector(Vector *vk, Vector *vb); // read vectors
void ReadE(Solution *s); // read basis transform coefficients in k-space and band path
void ReadTB(Solution *s); // read tight-binding Hamiltonian in k-space and band path

void GetName(Solution *s, char *data_type, char *fs); // get file name

void TransBasis(int e_len, lapack_complex_double *e, lapack_complex_double *v); // Transform basis (Q to sublattice) of matrices

void MakeBand(Solution *s); // make band structure data
void MakeUfw(Solution *s, double *w, lapack_complex_double *v); // make unfolding weight data
void MakeDos(Solution *s); // make density of states data

void FourierS(int l_len, int h_num, Vector v, Vector vq, Lattice *l, lapack_complex_double *h); // Fourier transform of single cell
void FourierD(int l_len, int h_num, Vector v, Vector vq, Lattice *l, lapack_complex_double *h); // Fourier transform of double cell

void InteractionS(Solution *s, lapack_complex_double *v_block); // add interaction term to Hamiltonian of single cell
void InteractionD(Solution *s, lapack_complex_double *v_block); // add interaction term to Hamiltonian of double cell

void OccupationS(double fermi, double *w, lapack_complex_double *v, lapack_complex_double *ek, double *n, double *m, double *e); // calculate occupation of single cell
void OccupationD(double fermi, double *w, lapack_complex_double *v, lapack_complex_double *ek, double *n, double *m, double *e); // calcultae occupation of double cell

#endif
