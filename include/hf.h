// hf.h : headers for Hartree-Fock approximation

#ifndef HF_H
#define HF_H

#define USE_MATH_DEFINES

#define DIM 3      // dimension
#define M_INIT 0.1 // initial magnetization
#define EP  0.02   // broadening rate of DOS

#define Nkg1 32    // num of Gauss-Legendre quadrature points in 1D
#define Nkg  32768 // num of Gauss-Legendre quadrature points in 3D = (Nkg1)^3
#define Nkb  1024  // num of band path points

#define CVG_MAX 3 // maximum iteration of checking convergence

#define CSQR(c) (creal(c)*creal(c) + cimag(c)*cimag(c))

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

typedef struct Lattice {
	int (*site)[DIM]; // site
	int (*obt)[2];    // orbital
	double (*t)[2];   // t
} Lattice;

typedef struct Configure {
	char *type;    // type of configure
	char *lat;     // type of lattice
	int Ni;        // num of atoms per unit cell
	int Nc;        // num of orbitals per atom
	int Ns;        // num of sites = Ni * Nc
	int Nb;        // num of bases = Ni * Nc * 2(spin up & dn)
	double Q[DIM]; // ordering vector
} Config;

typedef struct SelfConsistentSolution {
	char *save;       // name of directory to save output
	char *type;       // type of configure
	char runtime[16]; // runtime

	// input
	double JU;  // Hund coupling over Coulomb interaction
	double SOC; // spin-orbit coupling
	double N;   // target occupation
	double U;   // Coulomb interaction
	double J;   // Hund coupling

	// output
	double *n;    // occupation per orbital
	double *m;    // magnetization per orbital
	double ns;    // sum of occupation
	double ms;    // sum of total magnetization
	double fermi; // Fermi level	
	double dntop; // dntop energy
	double gap;   // band gap
	double e;     // energy
} Solution;

// lib/hf.c
void ReadConfig(Config *c); // read cell info in cell.txt
void InteractionS(Config c, Solution *s, lapack_complex_double *tb0); // add interaction term (Staggered m)
void InteractionQ(Config c, Solution *s, lapack_complex_double *tb0); // add interaction term (Q basis)
void BasisN(Config c, int Nk, lapack_complex_double *cf, lapack_complex_double *es); // basis transform
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

// lib/tb.c
void FourierN(Config c, int Nl, Lattice *lat, double *k, lapack_complex_double *tb); // Fourier transform
void FourierQ(Config c, int Nl, Lattice *lat, double *k, lapack_complex_double *tb); // Fourier transform (Q basis)
void CalcTB(Config c, char *ktype, void (*Fourier)()); // calculate tight-binding Hamiltonian matrices void DataName(Config c, Solution *s, char *dtype, char *dn); // data name void CalcEigen(Config c, Solution *s, LAPACK *lp, int Nk, lapack_complex_double *tb, double *ev, lapack_complex_double *es, Energy *e, void (*Interaction)()); // calculate eigenproblems void CalcGap(Config c, Solution *s, double *ev, lapack_complex_double *es); // calculate band gap void Interaction0(Config c, Solution *s, lapack_complex_double *tb0); // add interaction term

// lib/mod.c
void ReplaceStr(char *in, char *org, char *rep, char *out); // replace org->rep (out must be array, not pointer)
void ShowProgress(int i, int i_max); // show progress of iteration
void GetDimsH5(char *fn, char *dn, hsize_t *dims); // get dimenstion of hdf5 format
void ReadH5(char *fn, char *dn, double *val); // read hdf5 format
void WriteH5(char *fn, char *dn, int dim, hsize_t *dims, double *val); // write hdf5 format

#endif
