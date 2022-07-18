// hf3.h : universal headers for calculating 3 band models

#ifndef HF3_H
#define HF3_H

#define USE_MATH_DEFINES
#define OMP_THREAD (1)

#define OBT (3) // num of orbital
#define SUPER (2) // num of supercell
#define BASIS1 (6) // num of basis1
#define BASIS2 (12) // num of basis2
#define GREEN (0.05 / (pow(energy - w[i], 2) + pow(0.05, 2))) 
#define DOT_UNFOLD (path[i].x*r[n].x + path[i].y*r[n].y + path[i].z*r[n].z)

#define H(basis) (basis * basis) // size of Hamiltonian matrix at single k point
#define HK(basis) (K3 * H(basis)) // size of Hamiltonian matrix in k space
#define HB(basis) (BAND * H(basis)) // size of Hamiltonian matrix in band path
#define COMPLEX2(x) (pow(creal(x), 2) + pow(cimag(x), 2))

#include <omp.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lapacke.h>

typedef struct Vector {
	double x;
	double y;
	double z;
} Vector;

typedef struct SelfConsistentSolution {
	// input
	int basis; // num of bases
	char *type; // type of magnetic structure
	char runtime[16]; // runtime
	double JU; // Hund coupling / Coulomb interaction
	double SOC; // spin-orbit coupling
	double N; // occupation
	double U; // Coulomb interaction
	double J; // Hund coupling
	lapack_complex_double *tbk; // tight-binding Hamiltonian in k-space
	lapack_complex_double *tbb; // tight-binding Hamiltonian in band path

	// output
	double n[OBT]; // occupation per orbital
	double m[OBT]; // magnetization per orbital
	double n_total; // total occupation
	double m_total; // total magnetization
	double fermi; // Fermi level	
	double e; // energy
} Solution;

typedef struct Energy {
	double min;
	double max;
} Energy;

extern const int K; // num of k points in k space path
extern const int K3; // K*K*K
extern const int BAND; //num of kpoints in band path

extern void InteractionF(Solution *s, lapack_complex_double *v_tmp); // add interaction term to Hamiltonian of FM
extern void InteractionA(Solution *s, lapack_complex_double *v_tmp); // add interaction term to Hamiltonian of AF
extern void InteractionSubA(Solution *s, lapack_complex_double *v_tmp); // add interaction term to Hamiltonian of AF in sublattice basis
extern void FourierF(char *type, int num, Vector path, Vector q, lapack_complex_double *tb); // Fourier transform of FM
extern void FourierA(char *type, int num, Vector path, Vector q, lapack_complex_double *tb); // Fourier transform of AF
extern void FourierSubA(char *type, int num, Vector path, Vector q, lapack_complex_double *tb); // Fourier transform of AF in sublattice basis
extern void OccupationF(double fermi, double *w, lapack_complex_double *v, double *n, double *m, double *e); // calculate occupation of FM
extern void OccupationA(double fermi, double *w, lapack_complex_double *v, double *n, double *m, double *e); // calcultae occupation of AF
extern void OccupationSubA(double fermi, double *w, lapack_complex_double *v, double *n, double *m, double *e); // calcultae occupation of AF in sublattice basis

void CalcPathK(); // calculate k space path
void CalcPathB(); // calculate band path
void ReadPath(Vector *path_k, Vector *path_b); // read path binary data
void CalcEigenTB(char *type, double *w, lapack_complex_double *v); // calculate eigenproblems of tight-binding Hamiltonian
void MakeTB(char *type, char *fs, Vector *path, Vector q); // make band structure data (tight-binding only)

void ReadTB(Solution *s); // read tight-binding K, BAND Hamiltonian binary data
Energy CalcEigen(Solution *s, int k_num, lapack_complex_double *tb, double *w, lapack_complex_double *v); // calculate eigenproblems 
void CalcSolution(Solution *s, int is_unfold); // calculate self-consistent solution
void MakeBand(Solution *s, int is_unfold); // make band structure data
void MakeDos(Solution *s, int is_unfold); // make density of states data

#endif
