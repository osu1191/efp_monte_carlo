/*-
 * Copyright (c) 2012-2015 Ilya Kaliman
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

#ifndef EFPMD_COMMON_H
#define EFPMD_COMMON_H

#include <assert.h>
#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef EFP_USE_MPI
#include <mpi.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
 
#include <efp.h>
#include <ff.h>
#include <mathutil.h>

#include "cfg.h"
#include "msg.h"
#include "phys.h"

#define NORETURN __attribute__((noreturn))
#define efpmd_ARRAY_SIZE(x) (sizeof(x)/sizeof(x[0]))

enum run_type {
	RUN_TYPE_SP,
	RUN_TYPE_GRAD,
	RUN_TYPE_HESS,
	RUN_TYPE_OPT,
	RUN_TYPE_MD,
	RUN_TYPE_MC, // SKP
	RUN_TYPE_EFIELD,
	RUN_TYPE_ELPOT,
	RUN_TYPE_FRAG_ELPOT,
	RUN_TYPE_CLIB_TEST, // SKP 7th June
	RUN_TYPE_GTEST,
	RUN_TYPE_ETEST
};

enum ensemble_type {
	ENSEMBLE_TYPE_NVE,
	ENSEMBLE_TYPE_NVT,
	ENSEMBLE_TYPE_NPT
};

// Added by SKP

typedef struct {
    int type;       // Particle type
    double x, y, z; // Position coordinates
} Particle;

double calculateEnergy(Particle particles[6]);
double calculatePairwiseEnergy(Particle particle1, Particle particle2);
void monteCarloMove(Particle* particles, int num_particles, double temperature);
void monteCarloSimulation(Particle* particles, int num_particles, double temperature);

// =============//

struct efpmd_frag {
	char *name;
	double coord[12];
	size_t n_atoms;
//	size_t chosen_frag;  // SKP for frag_elpot
	struct efp_atom *atoms;
	double vel[6];
	bool constraint_enable;
	vec_t constraint_xyz;
	double constraint_k;
};

// SKP 6th June
/*
struct Tensor{
    void* data;
    int64_t* sizes;
    int ndim;
    int type_id;
    int is_variable;
};


struct TenGrad {
    float *data;
    int64_t *sizes;
    int num_dims;
};


//struct Net Net;
struct Net *createNet();
void destroyNet(struct Net *model);
void forward(struct Net* model, const float *inputs, float *output, int input_size, int output_size);
void trainModelWrapper(struct Net *model, const float **input_data, const float *target_data, int num_samples, int num_epochs, float learning_rate);
void generateEnergyWrapper(struct Net *model, const float **input_data, int batch_size, int input_size);

*/

// SKP
/*
struct efpmd_frag_nn {
	char *nn_name;	
        int chosen_frag;
        double nn_coord[12];
        size_t nn_natoms;
};
*/
struct charge {
	double q;
	vec_t pos;
};

struct sys {
	size_t n_frags;
	struct efpmd_frag *frags;
	size_t n_charges;
	struct charge *charges;
};

struct state {
	struct efp *efp;
	struct ff *ff;
	struct cfg *cfg;
	struct sys *sys;
	double energy;
	double *grad;
};

void NORETURN die(const char *, ...);
void NORETURN error(const char *, ...);

void *xmalloc(size_t);
void *xcalloc(size_t, size_t);
void *xrealloc(void *, size_t);

void print_vec(const double *);
void print_geometry(struct efp *);
void print_geometry_pbc(struct efp *, int);
void print_energy(struct state *);
void print_gradient(struct state *);
void print_fragment(const char *, const double *, const double *);
void print_charge(double, double, double, double);
void print_vector(size_t, const double *);
void print_matrix(size_t, size_t, const double *);
void print_pair_energy(struct state *);

void check_fail(enum efp_result);
void compute_energy(struct state *, bool);
struct sys *parse_input(struct cfg *, const char *);
six_t box_from_str(const char *);
int efp_strcasecmp(const char *, const char *);
int efp_strncasecmp(const char *, const char *, size_t);


#endif /* EFPMD_COMMON_H */