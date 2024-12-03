/*-
 * Copyright (c) 2012-2017 Ilya Kaliman
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

#ifndef LIBEFP_UTIL_H
#define LIBEFP_UTIL_H

#include "mathutil.h"

struct efp;
struct frag;

int efp_skip_frag_pair(const struct efp *, size_t, size_t);
struct swf efp_make_swf(const struct efp *, const struct frag *,
    const struct frag *, int);
int efp_check_rotation_matrix(const mat_t *);
void efp_points_to_matrix(const double *, mat_t *);
const struct frag *efp_find_lib(struct efp *, const char *);
void efp_add_stress(const vec_t *, const vec_t *, mat_t *);
void efp_add_force(six_t *, const vec_t *, const vec_t *,
    const vec_t *, const vec_t *);
void efp_sub_force(six_t *, const vec_t *, const vec_t *,
    const vec_t *, const vec_t *);
void efp_move_pt(const vec_t *, const mat_t *, const vec_t *, vec_t *);
void efp_rotate_t2(const mat_t *, const double *, double *);
void efp_rotate_t3(const mat_t *, const double *, double *);
mat_t rotmat_2frags(const mat_t *rotmat1, const mat_t *rotmat2);
int efp_strcasecmp(const char *, const char *);
int efp_strncasecmp(const char *, const char *, size_t);
// find plane by 3 points
void find_plane(const vec_t, const vec_t, const vec_t, vec_t *, double);
// computes maximum cut_off distance for an arbitrary periodic cell
double max_cutoff(const six_t);
// computes rmsd between two fragment structures
double calc_rmsd(const struct frag *frag1, const struct frag *frag2);

#endif /* LIBEFP_UTIL_H */
