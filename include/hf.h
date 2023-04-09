// hf.h : headers for Hartree-Fock approximation

#ifndef HF_H
#define HF_H

#define USE_MATH_DEFINES

#define DIM  3 // dimension
#define Nkg1 32 // num of Gauss-Legendre quadrature points (in 1D)
#define Nkg  32768 // num of Gauss-Legendre quadrature points (in 3D) = (Nkg1)^3
#define Nkb  1024 // num of band path points
//#define ETA  0.05 // broadening rate of DOS

#define CVG_MAX 3 // max count of checking convergence

#define CSQR(c) (pow(creal(c), 2) + pow(cimag(c), 2))

#include <omp.h>
#include <math.h>
#include <hdf5.h>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <lapack.h>
#include <sys/stat.h>
#include <gsl/gsl_integration.h>

typedef struct Coordinate {
	double c[DIM];
} Coord;

typedef struct Configure {
	char *name; // name of material
	char *save; // name of directory to save output
	char *type; // type of magnetic structure
	char sys[16]; // system
	char bas[16]; // basis
	int Ni; // num of atoms per unit cell
	int Nc; // num of orbitals per atom
	int Ns; // num of sites = Ni * Nc
	int Nb; // num of bases = Ni * Nc * 2(spin up and dn)
	int Nbb; // num of block = Nb * Nb
	double eta; // broadening rate of DOS
	Coord q; // ordering vector
} Config;

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

// lib/hf.c


// lib/gen.c
void GenGaussQuad(Config c); // generate points and weights for Gauss-Legendre quadrature


// lib/mod.c
void ReplaceStr(char *in, char *org, char *rep, char *out); // replace org->rep (out must be array, not pointer)
void ShowProgress(int i, int i_max); // show progress of iteration
void GetDimsH5(char *fn, char *dn, hsize_t *dims); // get dimenstion of hdf5 format
void ReadH5(char *fn, char *dn, double *val); // read hdf5 format
void WriteH5(char *fn, char *dn, int dim, hsize_t *dims, double *val); // write hdf5 format

void ReadConfig(Config *c); // read cell info in cell.txt

void DotProd(Config c, int Nk, Coord *k, Coord *r, lapack_complex_double *cf); // dot product
void CalcCoef(Config c, Coord *r); // calculate coefficients for basis transform
void ReadLat(char *name, int Nl, Lattice *l, char *ltype); // read lattice info in lat.txt
void Fourier0(Config c, int ith, Coord k, int Nl, Lattice *l, lapack_complex_double *tb); // Fourier transform
void FourierQ(Config c, int ith, Coord k, int Nl, Lattice *l, lapack_complex_double *tb); // Fourier transform (Q basis)
void CalcTB(Config c, char *ltype, void (*Fourier)()); // calculate tight-binding Hamiltonian matrices

void DataName(Config c, Solution *s, char *dtype, char *dn); // data name
void CalcEigen(Config c, Solution *s, LAPACK *lp, int Nk, lapack_complex_double *tb, double *ev, lapack_complex_double *es, Energy *e, void (*Interaction)()); // calculate eigenproblems 
void CalcGap(Config c, Solution *s, double *ev, lapack_complex_double *es); // calculate band gap
void Interaction0(Config c, Solution *s, lapack_complex_double *tb0); // add interaction term
void InteractionS(Config c, Solution *s, lapack_complex_double *tb0); // add interaction term (Staggered m)
void InteractionQ(Config c, Solution *s, lapack_complex_double *tb0); // add interaction term (Q basis)
void Basis0(Config c, int Nk, lapack_complex_double *cf, lapack_complex_double *es); // basis transform
void BasisQ(Config c, int Nk, lapack_complex_double *cf, lapack_complex_double *es); // basis transform (Q basis)
void Quadrature(Config c, Solution *s, double *wg, double *ev, lapack_complex_double *es, double *oc); // calculate occupation
void System0(double *n, double *m); // system
void SystemScA(double *n, double *m); // system (simple cubic - a)
void SystemScC(double *n, double *m); // system (simple cubic - c)
void SystemScG(double *n, double *m); // system (simple cubic - g)
void CalcSolution(Config c, Solution *s, LAPACK *lp, void (*System)(), void (*Interaction)(), void (*Basis)()); // calculate self-consistent solution
void MakeBand(Config c, Solution *s, LAPACK *lp, void (*Interaction)(), void (*Basis)()); // make band structure data
void MakeUFW(Config c, Solution *s, double *ev, lapack_complex_double *es, void (*Basis)()); // make unfolding weight data
void MakeDOS(Config c, Solution *s, LAPACK *lp, void (*Interaction)(), void (*Basis)()); // make density of states data

#endif
