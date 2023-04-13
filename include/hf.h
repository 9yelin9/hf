// hf.h : headers for Hartree-Fock approximation

#ifndef HF_H
#define HF_H

#define USE_MATH_DEFINES

#define DIM  3     // dimension
#define Nkg1 32    // num of Gauss-Legendre quadrature points in 1D
#define Nkg  32768 // num of Gauss-Legendre quadrature points in 3D = (Nkg1)^3
#define Nkb  1024  // num of band path points

#define M_INIT 0.1  // initial magnetization

#define FERMI_WIDTH 1e-4 // width of Fermi level
#define ITR_MAX     32   // maximum iteration of bisection
#define CVG_MAX     3    // maximum iteration of checking convergence

#define EP 0.02 // broadening rate of DOS

#define CSQR(c) (creal(c)*creal(c) + cimag(c)*cimag(c))

#include <omp.h>
#include <math.h>
#include <hdf5.h>
#include <time.h>
#include <stdio.h>
#include <dirent.h>
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
	char lat[16];     // type of lattice
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
void FileName(Solution *s, char *ftype, char *fn); // file name
void ReadConfig(Config *c); // read config_*.txt
void InitSolution(Config c, Solution *s); // initialize occupation and magnetization
void CalcEigen(Config c, double *ev, lapack_complex_double *es); // calculate eigenproblems
void CalcGap(Config c, Solution *s, double *ev, lapack_complex_double *es, double *uplow, double *dntop); // calculate band gap 
void CalcE(Config c, Solution *s, double w, double *ev, lapack_complex_double *es, double *e); // calculate sum of energy under dntop
void InteractionN(Config c, Solution *s, lapack_complex_double *tb); // add interaction term
void InteractionQ(Config c, Solution *s, lapack_complex_double *tb); // add interaction termi (Q basis)
void BasisN(Config c, double *uf, lapack_complex_double *es); // basis transform
void BasisQ(Config c, double *uf, lapack_complex_double *es); // basis transform (Q basis)
void Quadrature(Config c, Solution *s, double w, double *ev, lapack_complex_double *es, double *oc); // Gauss-Legendre quadrature
void GenSolution(Config c, Solution *s, void (*Symmetry)(), void (*Interaction)(), void (*Basis)()); // generate self-consistent solution by bisection
void GenSolBand(Config c, Solution *s, void (*Interaction)(), void (*Basis)()); // generate band structure from solution
void GenSolDOS(Config c, Solution *s, double ep, double *ev, lapack_complex_double *es); // generate density of states from solution
void GenDOS(Config c, Solution *s, char *fsn, double ep, void (*Interaction)(), void (*Basis)()); // generate density of states from solution

// lib/tb.c
void FourierN(Config c, int Nl, Lattice lat, double *k, lapack_complex_double *tb); // Fourier transform
void FourierQ(Config c, int Nl, Lattice lat, double *k, lapack_complex_double *tb); // Fourier transform (Q basis)
void GenTB(Config c, char *ktype, void (*Fourier)()); // generate tight-binding Hamiltonian
void GenTBBand(Config c); // generate band structure from tight-binding Hamiltonian

// lib/mod.c
void ReplaceStr(char *in, char *org, char *rep, char *out); // replace org->rep (out must be array, not pointer)
void ShowProgress(int i, int i_max); // show progress of iteration
void GetDimsH5(char *fn, char *dn, hsize_t *dims); // get dimenstion of hdf5 format
void ReadH5(char *fn, char *dn, double *val); // read hdf5 format
void WriteH5(char *fn, char *dn, int dim, hsize_t *dims, double *val); // write hdf5 format

#endif
