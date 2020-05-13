/*
 *	linux/drivers/cpufreq/cpufreq_adaptive.h
 *
 *	Header file for adaptive cpufreq governor template.
 *
 *	Copyright (C)  2020 Miłosz Malczak
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation.
 */


#include "cpufreq_governor.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define d_Ao 		(2)
#define d_Am		(3)
#define d_Rd		(1)
#define d_Sd		(1)
#define d_A		(2)
#define d_B		(1)

#define d_Bplus		(0)
#define d_Bmd		(0)
#define d_Bminus	(d_B)
#define d_R		(d_Sd+d_Bminus+d_Rd+d_Bplus)
#define d_S		(d_A + d_Rd-1 + d_Sd)
#define d_T		(d_Bmd + d_Ao)
#define d_D		MAX(d_R, d_Ao)
#define deg		(d_A+d_B+1)

#define d_ARd		(d_A+d_Rd)
#define d_BminusSd	(d_Sd+d_Bminus)
#define d_Sp		(d_ARd-1)
#define d_Rp		(d_BminusSd)
#define d_Acl		(d_ARd+d_BminusSd)
#define N		((d_Rp+1)+(d_Sp+1))

#if d_Acl != d_Ao+d_Am
#error
#endif


/* Controller */
#include <linux/string.h>

#define buf_length (1<<3)
#define BUF_IDX_ADD(BUF_IDX, VALUE) ((BUF_IDX+(VALUE))&(buf_length-1))

#define PARAMS_FILTER_LENGTH 1

struct adaptive_estimation_params
{
	int64_t P[deg*deg];
	int64_t theta[deg];		//value used for estimation
	int64_t theta_out[deg]; //filtered params or params set by the user
};

struct adaptive_controller_polynomials
{
	int64_t R[d_R+1];
	int64_t S[d_S+1];
	int64_t T[d_T+1];
	int64_t D[d_D+1];
};

struct adaptive_controller_buffers
{
	int idx;
	int64_t y[buf_length];
	int64_t u[buf_length];
	int64_t uc[buf_length];
	int64_t v[buf_length];
};

struct ParamsFilter
{
	int idx;
	int64_t buf[(deg)*PARAMS_FILTER_LENGTH];
};

/* Controller end*/

/* Fixed point arythmetics */
#include <linux/types.h>

#define POINT_POS (32)
#define FP(X) ((X)*((int64_t)1<<(POINT_POS)))
#define FP_MAX ((int64_t)0x7fffffffffffffff)
#define FP_MIN ((int64_t)0x8000000000000000)

static inline int64_t mult(int64_t a, int64_t b)
{
	int64_t a_h = a>>(POINT_POS/2);
	int64_t a_l = a-(a_h<<(POINT_POS/2));
	int64_t b_h = b>>(POINT_POS/2);
	int64_t b_l = b-(b_h<<(POINT_POS/2));

	int64_t c = (a>>(POINT_POS/2))*(b>>(POINT_POS/2));
	c += (a_h*b_l)>>(POINT_POS/2);
	c += (a_l*b_h)>>(POINT_POS/2);
	return c;
}

static inline int64_t division(int64_t a, int64_t b)
{
	int64_t c = 0;
	int counter = 0;
	while(abs(a)<((int64_t)1<<(62))) {
		a = a<<1;
		counter+=1;
		if(counter==POINT_POS)
			break;
	}
	c = a/(b)*((int64_t)1<<(POINT_POS-counter));
	return c;
}

/* Fixed point arythmetics end */

/* Maths */
int pow_10(int l)
{
	int i;
	int ret=1;
	for(i=0; i<l; i++)
		ret*=10;
	return ret;
}

/* Maths end*/

/* Matrix operations */
static void conv(int64_t *a, int64_t *b,int64_t *c, int a_l, int b_l)
{
	int i = 0;
	int j = 0;
	int l = 0;
	int up = 0;
	for(i=0; i<a_l+b_l-1; i++) {
		c[i] = 0;
		up = i;
		if(up>(a_l-1))
			up=(a_l-1);
		l = i-(b_l-1);
		if(l<0)
			l=0;
		for(j=up; j>=l; j--) {
			c[i] += mult(a[j], b[i-j]);
		}
	}
}

static inline int64_t multiply_vectors(int64_t *a, int64_t *b, int length)
{
	int i=0;
	int64_t y = 0;
	for(i=0; i<length; i++) {
		y += mult(a[i], b[i]);
	}
	return y;
}

static void multiply_matrices(int64_t *A, int64_t *B, int64_t *C,
				int A_row, int A_col, int B_row, int B_col)
{
	int i=0, row=0, col=0, A_idx=0, B_idx=0, C_idx=0;
	int64_t y = 0;
	for(row=0; row<A_row; row++) {
		for(col=0; col<B_col; col++) {
			C_idx = col + row * B_col;
			y = 0;
			for(i=0; i<A_col; i++) {
				A_idx = i + row * A_col;
				B_idx = i * B_col + col;;
				y += mult(A[A_idx], B[B_idx]);
			}
			C[C_idx] = y;
		}
	}
}

static inline void swap_vector_elements(int64_t *v, int pos_1, int pos_2)
{
	int64_t temp;
	temp = v[pos_1];
	v[pos_1] = v[pos_2];
	v[pos_2] = temp;
}

static void swap_rows_to_right_side(int64_t *A, int i_1, int i_2, int j)
{
	int64_t temp[N];
	int idx;
	int idx_2;
	int k;
	idx = i_1*N+j;
	for(k=j; k<N; k++) {
		temp[k] = A[idx];
		idx++;
	}
	idx = i_1*N+j;
	idx_2 = i_2*N+j;
	for(k=0; k<N-j; k++) {
		A[idx+k] = A[idx_2+k];
	}
	idx = i_2*N+j;
	for(k=j; k<N; k++) {
		A[idx] = temp[k];
		idx++;
	}
}
/* Matrix operations end*/


/* Linear solver */

static int partial_pivoting(int64_t *A, int64_t *b, int j)
{
	int i = 0;
	int64_t max = 0;
	int max_i = 0;
	int64_t temp = 0;
	for(i=j+1; i<N; i++) {
		temp = abs(A[i*N+j]);
		if(temp>max) {
			max = temp;
			max_i = i;
		}
	}
	if(max>abs(A[j*(N+1)])) {
		swap_rows_to_right_side(A, j, max_i, j);
		swap_vector_elements(b, j, max_i);
	}
	if(A[j*(N+1)]==0) {
		printk("ERROR, pivot == 0 \n");
		return -1;
	}
	return 0;
}

static void elimination_step(int64_t *A, int64_t *b, int j)
{
	int64_t pivot;
	int64_t multiplier = 0;
	int idx_top = 0;
	int idx = 0;
	int i = 0;
	int k = 0;
	idx = j*(N+1);
	pivot = A[idx];
	for(k=j; k<N; k++) {
		A[idx] = division(A[idx], pivot);
		idx++;
	}
	b[j] = division(b[j], pivot);

	for(i=j+1; i<N; i++) {
		idx = j+i*N;
		idx_top = j+j*N;
		multiplier = A[idx];
		b[i] = b[i] - mult(multiplier, b[j]);
		A[idx] = 0;
		for(k=j+1; k<N; k++)
		{
			idx++;
			idx_top++;
			A[idx] = A[idx] - mult(multiplier, A[idx_top]);
		}
	}
}

static int solve_linear_equation(int64_t *A_orig, int64_t *b_orig, int64_t *x)
{
	int err = 0;
	int i = 0;
	int j = 0;
	int idx = 0;
	int64_t A[N*N];
	int64_t b[N];

	/* Copy A and b matrices, otherwise original data will be lost*/
	for(idx=0; idx<N*N; idx++) {
		A[idx] = A_orig[idx];
	}
	for(idx=0; idx<N; idx++) {
		b[idx] = b_orig[idx];
	}
	/* elimination */
	for(j=0; j<N; j++) {
		err = partial_pivoting(A, b, j);
		if(err==-1)
			return err;
		elimination_step(A, b, j);
	}
	/* back-substitution */
	for(j=N-1; j>=0; j--) {
		for(i=j-1; i>=0; i--) {
			idx = j+i*N;
			b[i] = b[i] - mult(A[idx], b[j]);
			A[idx] = 0;
		}
	}
	for(idx=0; idx<N; idx++) {
		x[idx] = b[idx];
	}
	return err;
}

/* Linear solver end */

/* sysfs read/write */
ssize_t sscanf_fp(const char *buf, int64_t *value, int *buf_idx)
{
	int ret = 0;
	int negative = 0;
	int decimal = 0;
	int fractional = 0;
	int len = 0;
	*value = 0;
	if(buf[*buf_idx]=='-') {
		negative = 1;
		(*buf_idx)++;
	}
	ret = sscanf(buf+*buf_idx, "%d%n", &decimal, &len);
	if(ret!=1)
		return -EINVAL;
	*value += FP(decimal);
	*buf_idx += len;
	if(buf[*buf_idx] == '.' || buf[*buf_idx] == ',') {
		(*buf_idx)++;
		ret = sscanf(buf+*buf_idx, "%d%n", &fractional, &len);
		*buf_idx += len;
		if(ret!=1)
			return -EINVAL;
		*value += FP(fractional)/pow_10(len);
	}

	if(negative)
		*value = -1*(*value);
	return ret;
}

ssize_t sprintf_fp(char *buf, int64_t value, char end_char)
{
	int decimal;
	int fractional;
	int precision = 3;
	int64_t round = FP(0.999);
	/* Since we are working on fixed point arythmetics, the rounding does
	 * not work as expected. Division by the number a little less, but
	 * close to 1 gives expected results in rounding
	 */

	decimal = value/FP(1);
	value -= decimal*FP(1);
	value *= pow_10(precision);
	fractional = value/round;
	return sprintf(buf, "%d.%0*d%c", decimal, precision, abs(fractional),
			end_char);
}

/* sysfs read/write end */

struct adaptive_policy_dbs_info {
	struct policy_dbs_info policy_dbs;
	int err;
	int64_t phi_P_phi;
	struct adaptive_estimation_params est_params;
	struct adaptive_controller_polynomials poly;
	struct adaptive_controller_buffers buf;
	struct ParamsFilter p_filter;
};

static inline struct adaptive_policy_dbs_info *to_dbs_info(struct policy_dbs_info *policy_dbs)
{
	return container_of(policy_dbs, struct adaptive_policy_dbs_info, policy_dbs);
}

struct adaptive_dbs_tuners {
	int64_t lambda;
	int64_t theta_limit_up[deg];
	int64_t theta_limit_down[deg];
	int64_t Ao[d_Ao+1];
	int64_t Am[d_Am+1];
	int64_t Rd[d_Rd+1];
	int64_t Sd[d_Sd+1];
};
