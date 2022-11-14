// hf3.h : headers for calculating Hartree-Fock approximated 3-band Hubbard model

#ifndef HF3_H
#define HF3_H

#define USE_MATH_DEFINES
#define OMP_THREAD 16

#define DIM 3 // dimension
#define Nkg1 32 // num of Gauss-Legendre quadrature points (in 1D)
#define Nkg 32768 // num of Gauss-Legendre quadrature points (in 3D) = (Nkg1)^3
#define Nkb 1024 // num of band path points

#define Nl_MAX 200000 // max len of lat.txt
#define CVG_MAX 3 // max count of checking convergence

#define CSQR(c) (pow(creal(c), 2) + pow(cimag(c), 2))

#include <omp.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <lapacke.h>
#include <sys/stat.h>
#include <gsl/gsl_integration.h>

typedef struct Coordinate {
	double c[DIM];
} Coord;

typedef struct Cell {
	char *name; // name of material
	char *type; // type of magnetic structure
	char sys[16]; // system
	char bas[16]; // basis
	int Ni; // num of atoms per unit cell
	int Nc; // num of orbitals per atom
	int Ns; // num of sites = Ni * Nc
	int Nb; // num of bases = Ni * Nc * 2(spin up and dn)
	int Nbb; // num of block = Nb * Nb
	Coord q; // ordering vector
} Cell;

typedef struct Lattice {
	Coord d; // distance
	int obi; // initial orbital
	int obf; // final orbital
	double tre; // real part of t
	double tim; // imag part of t
} Lattice;

typedef struct LAPACK {
	char jobz;
	char uplo;
	double *rwork;
	lapack_int ln;
	lapack_int lda;
	lapack_int lwork;
	lapack_int info;
	lapack_complex_double *work;
} LAPACK;

typedef struct SelfConsistentSolution {
	char runtime[16];

	// input
	double JU; // Hund coupling per Coulomb interaction
	double SOC; // spin-orbit coupling
	double N; // target occupation
	double U; // Coulomb interaction
	double J; // Hund coupling

	// output
	double *n; // occupation per orbital
	double *m; // magnetization per orbital
	double ns; // sum of occupation
	double ms; // sum of total magnetization
	double fermi; // Fermi level	
	double dntop; // dntop energy
	double gap; // band gap
	double e; // energy
} Solution;

typedef struct Energy {
	double min;
	double max;
} Energy;

FILE* OpenFile(char *fn, char *mode); // open file
void ReadBin(char *fn, int size, void *v); // read binary files
void ReadCell(Cell *c); // read cell info in cell.txt

// mod/init.c
void CalcQuadPoints(Cell c); // calculate points and weights for Gauss-Legendre quadrature
void DotProd(Cell c, int Nk, Coord *k, Coord *r, lapack_complex_double *cf); // dot product
void CalcCoef(Cell c, Coord *r); // calculate coefficients for basis transform
void ReadLat(char *name, int *Nl, Lattice *l); // read lattice info in lat.txt
void Fourier0(Cell c, int ith, Coord k, int Nl, Lattice *l, lapack_complex_double *tb); // Fourier transform
void FourierQ(Cell c, int ith, Coord k, int Nl, Lattice *l, lapack_complex_double *tb); // Fourier transform (Q basis)
void CalcTB(Cell c, void (*Fourier)()); // calculate tight-binding Hamiltonian matrices

// mod/hf3.c
void DataName(Cell c, Solution *s, char *dtype, char *dn); // data name
void CalcEigen(Cell c, Solution *s, LAPACK *lp, int Nk, lapack_complex_double *tb, double *ev, lapack_complex_double *es, Energy *e, void (*Interaction)()); // calculate eigenproblems 
void CalcGap(Cell c, Solution *s, double *ev, lapack_complex_double *es); // calculate band gap
void Interaction0(Cell c, Solution *s, lapack_complex_double *tb0); // add interaction term
void InteractionS(Cell c, Solution *s, lapack_complex_double *tb0); // add interaction term (Staggered m)
void InteractionQ(Cell c, Solution *s, lapack_complex_double *tb0); // add interaction term (Q basis)
void Basis0(Cell c, int Nk, lapack_complex_double *cf, lapack_complex_double *es); // basis transform
void BasisQ(Cell c, int Nk, lapack_complex_double *cf, lapack_complex_double *es); // basis transform (Q basis)
void Quadrature(Cell c, Solution *s, double *wg, double *ev, lapack_complex_double *es, double *n, double *m); // calculate occupation
void System0(double *n, double *m); // system
void SystemScA(double *n, double *m); // system (simple cubic - a)
void SystemScC(double *n, double *m); // system (simple cubic - c)
void SystemScG(double *n, double *m); // system (simple cubic - g)
void CalcSolution(Cell c, Solution *s, LAPACK *lp, void (*System)(), void (*Interaction)(), void (*Basis)()); // calculate self-consistent solution
void MakeBand(Cell c, Solution *s, LAPACK *lp, void (*Interaction)(), void (*Basis)()); // make band structure data
void MakeUFW(Cell c, Solution *s, double *ev, lapack_complex_double *es, void (*Basis)()); // make unfolding weight data
void MakeDOS(Cell c, Solution *s, LAPACK *lp, void (*Interaction)(), void (*Basis)()); // make density of states data

#endif
