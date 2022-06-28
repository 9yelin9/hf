// hf3.h : universal headers for calculating 3 band models

#ifndef HF3_H
#define HF3_H

#define USE_MATH_DEFINES
#define OMP_THREAD (1)

#define OBT (3) // num of orbitals
#define H(basis) (basis * basis) // size of Hamiltonian matrix at single k point
#define HK(K3, basis) (K3 * H(basis)) // size of Hamiltonian matrix in k space
#define HB(BAND, basis) (BAND * H(basis)) // size of Hamiltonian matrix in band path
#define COMPLEX2(complex) (pow(creal(complex), 2) + pow(cimag(complex), 2))
#define GREEN(energy, eigenvalue) (0.05 / (pow(energy - eigenvalue, 2) + pow(0.05, 2))) 

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

extern void Interaction(Solution *s, lapack_complex_double *v_tmp); // add interaction term to Hamiltonian 
extern void FourierF(FILE *f, const int num, const int basis, Vector v, Vector q, lapack_complex_double *tb); // Fourier transform FM
extern void FourierA(FILE *f, const int num, const int basis, Vector v, Vector q, lapack_complex_double *tb); // Fourier transform AF
extern void OccupationF(int basis, double fermi, double *w, lapack_complex_double *v, double *n, double *m, double *e); // calculate occupation of FM
extern void OccupationA(int basis, double fermi, double *w, lapack_complex_double *v, double *n, double *m, double *e); // calcultae occupation of AF

void CalcK(Vector *v); // calculate k space path
void CalcBAND(Vector *v); // calculate band path
void CalcTB(char *fo_name, const int basis, const int k_num, const int tb_size, Vector *v, lapack_complex_double *tb); // calculate tight-binding Hamiltonian
void MakeTB(char *type, int basis, lapack_complex_double *tbb); // make band structure data (tight-binding only)

void ReadTB(char *type, int basis, lapack_complex_double *tbk, lapack_complex_double *tbb); // read tight-binding Hamiltonian binary data
Energy CalcEigen(Solution *s, const int k_num, lapack_complex_double *tb, double *w, lapack_complex_double *v); // calculate eigenproblems 
void CalcSolution(Solution *s); // calculate self-consistent solution
void MakeBand(Solution *s); // make band structure data
void MakeDos(Solution *s); // make density of states data

#endif
