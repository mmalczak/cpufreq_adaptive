// SPDX-License-Identifier: GPL-2.0-only
/*
 *	linux/drivers/cpufreq/cpufreq_adaptive.c
 *
 *	Copyright (C) 2020, NASK National Research Institute for Cybersecurity & AI
 *	Authors: 	Miłosz Malczak <milosz.malczak@nask.pl>
 *			Michał Karpowicz <michal.karpowicz@nask.pl>
 *			Michał Getka <michal.getka@nask.pl>
 */

#define TELEMETRY 1

#define pr_fmt(fmt) KBUILD_MODNAME ": " fmt

#include "cpufreq_adaptive.h"

#include <linux/kernel.h>

#if TELEMETRY
/* Header to support the kernel Driver Model */
#include <linux/device.h>
/* Required for the copy to user function */
#include <linux/uaccess.h>
/* Timesatmp handling */
#include <linux/time.h>

#include <linux/fs.h>

// Class name should be unique for kernel module module
#define  CLASS_NAME  ( KBUILD_MODNAME "tlm" )
#include "tlmt.h"

static int	  device_major_number;
static struct class*  tlmsrv_class	= NULL;
static struct device* tlmsrv_device = NULL;

/* Only one user-space application can open device */
static DEFINE_MUTEX(tlmsrv_mutex);
/* Data buffer access mutex */
static DEFINE_MUTEX(tlmbuff_mutex);

static void tlm_add_sample(struct tlm_sample *);
static int tlmsrv_init(void);
static void tlmsrv_exit(void);
#endif


/* Fixed point arythmetics */

static inline int64_t mult(int64_t a, int64_t b)
{
	int64_t a_h = a>>(POINT_POS / 2);
	int64_t a_l = a-(a_h<<(POINT_POS / 2));
	int64_t b_h = b>>(POINT_POS / 2);
	int64_t b_l = b-(b_h<<(POINT_POS / 2));

	int64_t c = (a>>(POINT_POS / 2)) * (b>>(POINT_POS / 2));
	c += (a_h * b_l)>>(POINT_POS / 2);
	c += (a_l * b_h)>>(POINT_POS / 2);
	return c;
}

static inline int64_t division(int64_t a, int64_t b)
{
	int64_t c = 0;
	int counter = 0;
	while (abs(a) < ((int64_t)1<<(62))) {
		a = a<<1;
		counter += 1;
		if(counter == POINT_POS)
			break;
	}
	c = a / (b) * ((int64_t)1<<(POINT_POS - counter));
	return c;
}

/* Fixed point arythmetics end */

/* Maths */
static inline int pow_10(int l)
{
	int i;
	int ret = 1;
	for (i = 0; i < l; i++)
		ret *= 10;
	return ret;
}

static void conv(int64_t *a, int64_t *b,int64_t *c, int a_l, int b_l)
{
	int i = 0;
	int j = 0;
	int l = 0;
	int up = 0;
	for (i = 0; i < a_l + b_l - 1; i++) {
		c[i] = 0;
		up = i;
		if (up > (a_l - 1))
			up = (a_l - 1);
		l = i - (b_l - 1);
		if (l < 0)
			l = 0;
		for (j = up; j >= l; j--) {
			c[i] += mult(a[j], b[i - j]);
		}
	}
}
/* Maths end*/

/* Matrix operations */
static inline int64_t dot_product(int N, int64_t *a, int64_t *b)
{
	int i = 0;
	int64_t y = 0;
	for (i = 0; i < N; i++) {
		y += mult(a[i], b[i]);
	}
	return y;
}

/*
 * M - number of rows of the matrix A
 * N - number of columns of the matrix B
 * K - number of columns of the matrix A equal
 * 	to the number of rows of the matrix B
 * A, B - input matrices
 * C - output matrix
 */
static void multiply_matrices(int M, int N, int K,
				int64_t *A, int64_t *B, int64_t *C)
{
	int i = 0;
	int row = 0;
	int col = 0;
	int A_idx = 0;
	int B_idx = 0;
	int C_idx = 0;
	int64_t y = 0;
	for (row = 0; row < M; row++) {
		for (col = 0; col < N; col++) {
			C_idx = col + row * N;
			y = 0;
			for (i = 0; i < K; i++) {
				A_idx = i + row * K;
				B_idx = i * N + col;;
				y += mult(A[A_idx], B[B_idx]);
			}
			C[C_idx] = y;
		}
	}
}
/* Matrix operations end*/

/* Linear solver */

static inline void swap_vector_elements(int64_t *v, int pos_1, int pos_2)
{
	int64_t temp;
	temp = v[pos_1];
	v[pos_1] = v[pos_2];
	v[pos_2] = temp;
}

static void swap_right_side_of_rows(int N, int64_t *A, int i_1, int i_2, int j)
{
	int64_t temp;
	int idx;
	int idx_2;
	int k;
	idx = i_1 * N + j;
	idx_2 = i_2 * N + j;
	for (k = j; k < N; k++) {
		temp = A[idx];
		A[idx] = A[idx_2];
		A[idx_2] = temp;
		idx++;
		idx_2++;
	}
}

static int partial_pivoting(int N, int64_t *A, int64_t *b, int j)
{
	int i = 0;
	int64_t max = 0;
	int max_i = 0;
	int64_t temp = 0;
	for (i = j + 1; i < N; i++) {
		temp = abs(A[i * N + j]);
		if (temp > max) {
			max = temp;
			max_i = i;
		}
	}
	if (max > abs(A[j * (N + 1)])) {
		swap_right_side_of_rows(N, A, j, max_i, j);
		swap_vector_elements(b, j, max_i);
	}
	if (A[j * (N + 1)] == 0) {
		printk("ERROR, pivot == 0 \n");
		return -1;
	}
	return 0;
}

static void elimination_step(int N, int64_t *A, int64_t *b, int j)
{
	int64_t pivot;
	int64_t multiplier = 0;
	int idx_top = 0;
	int idx = 0;
	int i = 0;
	int k = 0;
	idx = j * (N + 1);
	pivot = A[idx];
	for (k = j; k < N; k++) {
		A[idx] = division(A[idx], pivot);
		idx++;
	}
	b[j] = division(b[j], pivot);

	for (i = j + 1; i < N; i++) {
		idx = j + i * N;
		idx_top = j + j * N;
		multiplier = A[idx];
		b[i] = b[i] - mult(multiplier, b[j]);
		A[idx] = 0;
		for (k = j + 1; k < N; k++)
		{
			idx++;
			idx_top++;
			A[idx] = A[idx] - mult(multiplier, A[idx_top]);
		}
	}
}

/*
 * Solves a set of linear equations described by the matrix equation
 * Ab = x, where A is NxN matrix, b and x are vectors of length N.
 * The solver modifies the values of A and b. Therefore, if their values are
 * to be used afterwards, the copy of them has to be made before calling this
 * function
 */
static int solve_linear_equation_inplace(int N, int64_t *A, int64_t *b,
					int64_t *x)
{
	int err = 0;
	int i = 0;
	int j = 0;
	int idx = 0;

	/* elimination */
	for (j = 0; j < N; j++) {
		err = partial_pivoting(N, A, b, j);
		if (err == -1)
			return err;
		elimination_step(N, A, b, j);
	}
	/* back-substitution */
	for (j = N - 1; j >= 0; j--) {
		for (i = j - 1; i >= 0; i--) {
			idx = j + i * N;
			b[i] = b[i] - mult(A[idx], b[j]);
			A[idx] = 0;
		}
	}
	for (idx = 0; idx < N; idx++) {
		x[idx] = b[idx];
	}
	return err;
}

/* Linear solver end */

/* Controller */

void update_estimation(int64_t *theta, int64_t *P, int64_t lambda,
			int64_t *y_buf, int64_t *u_buf, int buf_idx,
			int64_t *phi_P_phi_);
int controller_synthesis(int64_t *A, int64_t *Bplus, int64_t *Bminus, int64_t *Bmd,
			struct adaptive_dbs_tuners *tuners,
			struct adaptive_controller_polynomials *poly);


void restart_controller(struct adaptive_policy_dbs_info *dbs_info)
{
	int i;
	int j;

	dbs_info->err = 0;
	dbs_info->phi_P_phi = 0;
	memset(&(dbs_info->est_params), 0, sizeof(struct adaptive_estimation_params));
	memset(&(dbs_info->poly), 0, sizeof(struct adaptive_controller_polynomials));
	memset(&(dbs_info->buf), 0, sizeof(struct adaptive_controller_buffers));
	memset(&(dbs_info->filt_params), 0, sizeof(struct filter_params));

	for (i = 0; i < deg; i++) {
		for (j = 0; j < deg; j++) {
			if (j == i)
				dbs_info->est_params.P[i * deg + j] = FP(100);
			else
				dbs_info->est_params.P[i * deg + j] = FP(0);
		}
	}
	for (i = 0; i < deg; i++) {
		if (i < d_A)
			dbs_info->est_params.theta[i] = FP(0.1);
		else
			dbs_info->est_params.theta[i] = FP(0.1);
	}
}

void filter_params(struct adaptive_policy_dbs_info *dbs_info, int64_t *params,
					int64_t *params_filtered)
{
	int i, j;
	for (i = 0; i < deg; i++)
		params_filtered[i] = 0;
	dbs_info->filt_params.idx += 1;
	if (dbs_info->filt_params.idx == PARAMS_FILTER_LENGTH)
		dbs_info->filt_params.idx = 0;
	for (i = 0; i < (deg); i++)
		dbs_info->filt_params.buf[dbs_info->filt_params.idx + i * PARAMS_FILTER_LENGTH] = params[i];

	for (i = 0; i < PARAMS_FILTER_LENGTH; i++) {
		for (j = 0; j < (deg); j++) {
			params_filtered[j] += dbs_info->filt_params.buf[i + j * PARAMS_FILTER_LENGTH];
		}
	}
	for (i = 0; i < (deg); i++) {
		params_filtered[i] = params_filtered[i] / PARAMS_FILTER_LENGTH;
	}

}

int64_t regulate(struct adaptive_dbs_tuners *tuners,
		struct adaptive_policy_dbs_info *dbs_info, int64_t y, int64_t u,
		int64_t uc)
{
	int i;
	int64_t v;

	int64_t Bplus[d_Bplus + 1] = {FP(1.0)};
	int64_t Bmd[d_Bmd + 1] = {FP(1.0)};
	int64_t A[d_A + 1];
	int64_t Bminus[d_Bminus + 1];

	// Update buffers
	dbs_info->buf.idx = BUF_IDX_ADD(dbs_info->buf.idx, -1);
	dbs_info->buf.y[dbs_info->buf.idx] = y;
	dbs_info->buf.u[dbs_info->buf.idx] = u;
	dbs_info->buf.uc[dbs_info->buf.idx] = uc;

	update_estimation(dbs_info->est_params.theta,
		dbs_info->est_params.P, tuners->lambda, dbs_info->buf.y,
		dbs_info->buf.u, dbs_info->buf.idx, &(dbs_info->phi_P_phi));
	for (i = 0; i < deg; i++) {
		if (dbs_info->est_params.theta[i] > tuners->theta_limit_up[i])
			dbs_info->est_params.theta[i] = tuners->theta_limit_up[i];
		else if (dbs_info->est_params.theta[i] < tuners->theta_limit_down[i])
			dbs_info->est_params.theta[i] = tuners->theta_limit_down[i];
	}
	filter_params(dbs_info, dbs_info->est_params.theta,
			dbs_info->est_params.theta_out);

	A[0] = FP(1.0);
	for (i = 0; i < d_A; i++)
		A[1 + i] = dbs_info->est_params.theta_out[i];
	for (i = 0; i < d_B + 1; i++)
		Bminus[i] = dbs_info->est_params.theta_out[d_A + i];
	dbs_info->err = controller_synthesis(A, Bplus, Bminus, Bmd,
						tuners, &(dbs_info->poly));
	if (dbs_info->err == -1)
		return dbs_info->buf.u[dbs_info->buf.idx];
	/*
	 * Rv=T*uc−S*y+D(v-u)
	 * with anti_windup(with assumption that R and Ao are monic)
	 */
	v = 0;
	for (i = 1; i < d_R+1; i++)
		v = v - mult(dbs_info->poly.R[i],
				dbs_info->buf.v[BUF_IDX_ADD(dbs_info->buf.idx, i)]);
	for (i = 0; i < d_S + 1; i++)
		v = v - mult(dbs_info->poly.S[i],
				dbs_info->buf.y[BUF_IDX_ADD(dbs_info->buf.idx, i)]);
	for (i = 0; i < d_T + 1; i++) {
		v = v + mult(dbs_info->poly.T[i],
				dbs_info->buf.uc[BUF_IDX_ADD(dbs_info->buf.idx, i)]);
	}
	for (i = 0; i < d_D; i++) {
		v = v + mult(dbs_info->poly.D[i + 1],
			(dbs_info->buf.v[BUF_IDX_ADD(dbs_info->buf.idx, i + 1)]-
			dbs_info->buf.u[BUF_IDX_ADD(dbs_info->buf.idx, i)]));
	}
	dbs_info->buf.v[dbs_info->buf.idx] = v;
	return v;
}

static void conditionally_update_theta(int64_t phi_P_phi, int64_t *theta,
						int64_t *K, int64_t epsilon)
{
	int i = 0;
	int64_t kappa = 0;

	if(phi_P_phi > FP(1000))
		kappa = 1;
	if(kappa)
		for (i = 0; i < deg; i++)
			theta[i] += mult(K[i], epsilon);
}

static void constant_trace_covariance_matrix(int64_t *P)
{
	int i = 0;
	int64_t trace = 0;
	for (i = 0; i < deg; i++) {
		trace += P[i * (deg + 1)];
	}
	for (i = 0; i < deg * deg; i++) {
		P[i] = division(P[i], trace);
		P[i] /= 1000;
	}
	for (i = 0; i < deg; i++) {
		P[i * (deg + 1)] = P[i * (deg + 1)] + FP(1)/(d_A + d_B + 1);
	}
}

void update_estimation(int64_t *theta, int64_t *P, int64_t lambda,
			int64_t *y_buf, int64_t *u_buf, int buf_idx,
			int64_t *phi_P_phi_)
{
	int i;
	int64_t phi[deg];
	int64_t P_phi[deg];
	int64_t phi_P_phi[1];
	int64_t K[deg];
	int64_t K_phi[deg*deg];
	int64_t K_phi_P[deg*deg];
	int64_t K_denominator;
	int64_t epsilon;

	for (i = 0; i < d_A; i++) {
		phi[i] = -y_buf[(buf_idx + 1 + i)&(buf_length - 1)];
	}
	for (i = 0; i < d_B + 1; i++) {
		phi[d_A + i] = u_buf[(buf_idx + i)&(buf_length - 1)];
	}
	multiply_matrices(deg, 1, deg, P, phi, P_phi);
	multiply_matrices(1, 1, deg, phi, P_phi, phi_P_phi);
	*phi_P_phi_ = phi_P_phi[0];
	K_denominator = lambda + phi_P_phi[0];

	for(i=0; i<deg; i++) {
		K[i] = division(P_phi[i], K_denominator);
	}

	epsilon = y_buf[buf_idx] - dot_product(deg, phi, theta);

	conditionally_update_theta(phi_P_phi[0], theta, K, epsilon);

	multiply_matrices(deg, deg, 1, K, phi, K_phi);
	multiply_matrices(deg, deg, deg, K_phi, P, K_phi_P);
	for (i = 0; i < deg * deg; i++) {
		P[i] = P[i] - K_phi_P[i];
		P[i] = division(P[i], lambda);
	}
	constant_trace_covariance_matrix(P);
}

static void combine_model_and_custom_filters(int64_t *A, int64_t *Rd, int64_t *ARd,
			int64_t *Bminus, int64_t *Sd, int64_t *BminusSd)
{
	conv(A, Rd, ARd, d_A + 1, d_Rd + 1);
	conv(Bminus, Sd, BminusSd, d_Bminus + 1, d_Sd + 1);
}

static inline void combine_observer_polynomial_and_modeled_dynamics(
					int64_t *Ao, int64_t *Am, int64_t *c)
{
	conv(Ao, Am, c, d_Ao + 1, d_Am + 1);
}

static void fill_autoregresive_coeff(int64_t *M, int N, int64_t *ARd, int i, int j)
{
	int idx = - j + i;
	if((idx < 0) || (idx > d_ARd))
		M[i * N + j] = 0;
	else
		M[i * N + j] = ARd[idx];
}

static void fill_moving_average_coeff(int64_t *M, int N, int64_t *BminusSd, int i, int j)
{
	int idx = - (j - (d_Rp + 1)) + i - (d_Acl - (d_BminusSd + d_Sp));
	if ((idx < 0) || (idx > d_BminusSd))
		M[i * N + j] = 0;
	else
		M[i * N + j] = BminusSd[idx];
}

static int normalise_controller_gain(int64_t *Am, int64_t *Ao, int64_t *Bminus,
					int64_t *Bmd, int64_t *T)
{
	int64_t BmdAo[d_Bmd + d_Ao + 1];
	int64_t Bm[d_Bminus+d_Bmd + 1];
	int64_t Am_sum = 0;
	int64_t Bm_sum = 0;
	int64_t beta = 0;
	int i;

	conv(Bminus, Bmd, Bm, d_Bminus + 1, d_Bmd + 1);
	conv(Bmd, Ao, BmdAo, d_Bmd + 1, d_Ao + 1);
	for (i = 0; i <= d_Am; i++)
		Am_sum += Am[i];
	for (i = 0; i <= d_Bminus + d_Bmd; i++)
		Bm_sum += Bm[i];
	if(Bm_sum == 0)
		return -1;
	beta = division(Am_sum, Bm_sum);

	for (i = 0; i <= d_Bmd + d_Ao; i++)
		T[i] = mult(BmdAo[i], beta);
	return 0;

}

static void calculate_anti_windup_polynomial(int64_t *D, int64_t *R, int64_t *Ao)
{
	int i;
	for (i = 0; i <= d_D; i++) {
		if (i <= d_R)
			D[i] = R[i];
		else
			D[i] = 0;
	}
	for (i = 0; i <= d_Ao; i++) {
		D[i] = D[i] - Ao[i];
	}
}

int controller_synthesis(int64_t *A, int64_t *Bplus, int64_t *Bminus, int64_t *Bmd,
			struct adaptive_dbs_tuners *tuners,
			struct adaptive_controller_polynomials *poly)
{
	int err = 0;
	int i = 0;
	int j = 0;

	int64_t M[d_M * d_M];
	int64_t c[d_M];
	int64_t temp[d_M];
	int64_t ARd[d_ARd + 1];
	int64_t BminusSd[d_BminusSd + 1];
	int64_t Rtemp[d_Rp + d_Rd + 1];

	combine_model_and_custom_filters(A, tuners->Rd, ARd,
					Bminus, tuners->Sd, BminusSd);
	for (i = 0; i < d_M; i++) {
		for (j = 0; j < d_M; j++) {
			if (j < d_Rp + 1)
				fill_autoregresive_coeff(M, d_M, ARd, i, j);
			else
				fill_moving_average_coeff(M, d_M, BminusSd, i, j);
		}
	}
	combine_observer_polynomial_and_modeled_dynamics(tuners->Ao,
							tuners->Am, c);
	err = solve_linear_equation_inplace(d_M, M, c, temp);
	if (err != 0)
		return err;;
	conv(temp, tuners->Rd, Rtemp, d_Rp + 1, d_Rd + 1);
	conv(Rtemp, Bplus, poly->R, d_Rp + d_Rd + 1, d_Bplus + 1);
	conv(&temp[d_Rp + 1], tuners->Sd, poly->S, d_Sp + 1, d_Sd + 1);
	err = normalise_controller_gain(tuners->Am, tuners->Ao, Bminus, Bmd,
					poly->T);
	if (err == -1)
		return err;
	calculate_anti_windup_polynomial(poly->D, poly->R, tuners->Ao);
	return err;
}

/* Controller end*/

static unsigned int adaptive_dbs_update(struct cpufreq_policy *policy)
{
	struct policy_dbs_info *policy_dbs = policy->governor_data;
	struct dbs_data *dbs_data = policy_dbs->dbs_data;
	struct adaptive_policy_dbs_info *dbs_info = to_dbs_info(policy_dbs);
	struct adaptive_dbs_tuners *tuners = dbs_data->tuners;

	unsigned int load = dbs_update(policy);

	int64_t uc;
	static uint64_t counter = 0;
	int u;
	static int buf_idx = 0;
	unsigned int freq_next;
	static int load_counter = 0;
	int64_t v;
	int i;
#if TELEMETRY
	struct tlm_sample sample;
#endif

	counter++;
	if ((counter / 200)%2)
		uc = 60;
	else
		uc = 90;

	buf_idx = BUF_IDX_ADD(buf_idx, 1);

#define scale 20000
	u = cpufreq_quick_get(0);
	v = regulate(tuners, dbs_info, FP(load), FP(u)/scale, FP(uc))/FP(1);
	v = v*scale;
	if (load == 100)
		load_counter++;
	else
		load_counter = 0;
	if (load_counter == 5)
		v = 2710000;

#if TELEMETRY
	sample.u_prev = u;
	sample.v = v;
	sample.load = load;
	sample.uc = uc;
	for (i = 0; i < deg; i++)
		sample.theta[i] = dbs_info->est_params.theta_out[i];
	for (i = 0; i < d_R + 1; i++)
		sample.R[i] = dbs_info->poly.R[i];
	for (i = 0; i < d_S + 1; i++)
		sample.S[i] = dbs_info->poly.S[i];
	for (i = 0; i < d_T + 1; i++)
		sample.T[i] = dbs_info->poly.T[i];
	for (i = 0; i < d_D + 1; i++)
		sample.D[i] = dbs_info->poly.D[i];
	for (i = 0; i < deg * deg; i++)
		sample.P[i] = dbs_info->est_params.P[i];
	sample.phi_P_phi = dbs_info->phi_P_phi;

	tlm_add_sample(&sample);
#endif

	if (dbs_info->err == -1)
			restart_controller(dbs_info);
	if (v < 0)
		freq_next = 0;
	else
		freq_next = (unsigned int)v;


	__cpufreq_driver_target(policy, freq_next, CPUFREQ_RELATION_C);

	return dbs_data->sampling_rate * policy_dbs->rate_mult;
}


/************************** sysfs interface ************************/

/* sysfs read/write */
static ssize_t sscanf_fp(const char *buf, int64_t *value, int *buf_idx)
{
	int ret = 0;
	int negative = 0;
	int decimal = 0;
	int fractional = 0;
	int len = 0;
	*value = 0;
	if (buf[*buf_idx] == '-') {
		negative = 1;
		(*buf_idx)++;
	}
	ret = sscanf(buf+*buf_idx, "%d%n", &decimal, &len);
	if (ret != 1)
		return -EINVAL;
	*value += FP(decimal);
	*buf_idx += len;
	if (buf[*buf_idx] == '.' || buf[*buf_idx] == ',') {
		(*buf_idx)++;
		ret = sscanf(buf+*buf_idx, "%d%n", &fractional, &len);
		*buf_idx += len;
		if (ret != 1)
			return -EINVAL;
		*value += FP(fractional)/pow_10(len);
	}

	if (negative)
		*value = -1 * (*value);
	return ret;
}

static size_t sprintf_fp(char *buf, int64_t value, char end_char)
{
	int decimal;
	int fractional;
	int precision = 3;
	int64_t round = FP(0.999);
	/*
	 * Since we are working on fixed point arythmetics, the rounding does
	 * not work as expected. Division by the number a little less, but
	 * close to 1 gives expected results in rounding
	 */

	decimal = value / FP(1);
	value -= decimal * FP(1);
	value *= pow_10(precision);
	fractional = value / round;
	return sprintf(buf, "%d.%0*d%c", decimal, precision, abs(fractional),
			end_char);
}

/* sysfs read/write end */

gov_show_one_common(sampling_rate);

gov_attr_rw(sampling_rate);

/* Specific parameters for cpufreq_adaptive */


static ssize_t store_lambda(struct gov_attr_set *attr_set, const char *buf,
				size_t count)
{
	struct dbs_data *dbs_data = to_dbs_data(attr_set);
	struct adaptive_dbs_tuners *tuners = dbs_data->tuners;
	int64_t input;
	int ret;
	int buf_idx=0;
	ret = sscanf_fp(buf, &input, &buf_idx);

	if (ret != 1) {
		return -EINVAL;
	}
	if (input > FP(1) || input <= FP(0)) {
		return -EINVAL;
	}
	tuners->lambda = input;

	return count;
}

gov_store_vector_adaptive(theta_limit_up, deg);
gov_store_vector_adaptive(theta_limit_down, deg);
gov_store_vector_adaptive(Ao, (d_Ao + 1));
gov_store_vector_adaptive(Am, (d_Am + 1));
gov_store_vector_adaptive(Rd, (d_Rd + 1));
gov_store_vector_adaptive(Sd, (d_Sd + 1));

gov_show_one_adaptive(lambda);
gov_show_vector_adaptive(theta_limit_up, deg);
gov_show_vector_adaptive(theta_limit_down, deg);
gov_show_vector_adaptive(Ao, d_Ao + 1);
gov_show_vector_adaptive(Am, d_Am + 1);
gov_show_vector_adaptive(Rd, d_Rd + 1);
gov_show_vector_adaptive(Sd, d_Sd + 1);

gov_attr_rw(lambda);
gov_attr_rw(theta_limit_up);
gov_attr_rw(theta_limit_down);
gov_attr_rw(Ao);
gov_attr_rw(Am);
gov_attr_rw(Rd);
gov_attr_rw(Sd);

static struct attribute *adaptive_attributes[] = {
	&sampling_rate.attr,
/* Specific parameters for cpufreq_adaptive */
	&lambda.attr,
	&theta_limit_up.attr,
	&theta_limit_down.attr,
	&Ao.attr,
	&Am.attr,
	&Rd.attr,
	&Sd.attr,
	NULL
};

/************************** sysfs end ************************/

static struct policy_dbs_info *adaptive_alloc(void)
{
	struct adaptive_policy_dbs_info *dbs_info;

	dbs_info = kzalloc(sizeof(*dbs_info), GFP_KERNEL);
	return dbs_info ? &dbs_info->policy_dbs : NULL;
}

static void adaptive_free(struct policy_dbs_info *policy_dbs)
{
	kfree(to_dbs_info(policy_dbs));
}

static void adaptive_dbs_tuners_init(struct adaptive_dbs_tuners *tuners)
{
	tuners->lambda = FP(0.99);
	tuners->theta_limit_up[0] = FP(-1);
	tuners->theta_limit_up[1] = FP(1000);
	tuners->theta_limit_up[2] = FP(-1);
	tuners->theta_limit_up[3] = FP(1000);
	tuners->theta_limit_down[0] = FP(-1000);
	tuners->theta_limit_down[1] = FP(-1000);
	tuners->theta_limit_down[2] = FP(-1000);
	tuners->theta_limit_down[3] = FP(-1000);
	tuners->Ao[0] = FP(1);
	tuners->Ao[1] = FP(0);
	tuners->Ao[2] = FP(0);
	tuners->Am[0] = FP(1);
	tuners->Am[1] = FP(-0.5);
	tuners->Am[2] = FP(0);
	tuners->Am[3] = FP(0);
	tuners->Rd[0] = FP(1);
	tuners->Rd[1] = FP(-1);
	tuners->Sd[0] = FP(1);
	tuners->Sd[1] = FP(0);
}

static int adaptive_init(struct dbs_data *dbs_data)
{
	struct adaptive_dbs_tuners *tuners;

	tuners = kzalloc(sizeof(*tuners), GFP_KERNEL);
	if (!tuners)
		return -ENOMEM;

	dbs_data->tuners = tuners;
	adaptive_dbs_tuners_init(dbs_data->tuners);
	return 0;
}

static void adaptive_exit(struct dbs_data *dbs_data)
{
	kfree(dbs_data->tuners);
}

static void adaptive_start(struct cpufreq_policy *policy)
{
	struct adaptive_policy_dbs_info *dbs_info = to_dbs_info(policy->governor_data);
	restart_controller(dbs_info);
}

static struct dbs_governor adaptive_dbs_gov = {
	.gov = CPUFREQ_DBS_GOVERNOR_INITIALIZER("adaptive"),
	.kobj_type = { .default_attrs = adaptive_attributes },
	.gov_dbs_update = adaptive_dbs_update,
	.alloc = adaptive_alloc,
	.free = adaptive_free,
	.init = adaptive_init,
	.exit = adaptive_exit,
	.start = adaptive_start,
};

#define CPU_FREQ_GOV_ADAPTIVE	(&adaptive_dbs_gov.gov)

static int __init cpufreq_gov_adaptive_init(void)
{
#if TELEMETRY
	tlmsrv_init();
#endif
	return cpufreq_register_governor(CPU_FREQ_GOV_ADAPTIVE);
}

static void __exit cpufreq_gov_adaptive_exit(void)
{
#if TELEMETRY
	tlmsrv_exit();
#endif
	cpufreq_unregister_governor(CPU_FREQ_GOV_ADAPTIVE);
}


MODULE_AUTHOR("Miłosz Malczak <milosz.malczak@nask.pl>");
MODULE_DESCRIPTION("CPUfreq policy governor 'adaptive'");
MODULE_LICENSE("GPL");

#ifdef CONFIG_CPU_FREQ_DEFAULT_GOV_ADAPTIVE
struct cpufreq_governor *cpufreq_default_governor(void)
{
	return CPU_FREQ_GOV_ADAPTIVE;
}

fs_initcall(cpufreq_gov_adaptive_init);
#else
module_init(cpufreq_gov_adaptive_init);
#endif
module_exit(cpufreq_gov_adaptive_exit);

#if TELEMETRY
static void tlm_add_sample(struct tlm_sample *sample)
{

  struct timespec64 ts;

  // update sample time
  ktime_get_real_ts64(&ts);
  sample->ts_sec = ts.tv_sec;
  sample->ts_nsec = ts.tv_nsec;

  mutex_lock(&tlmbuff_mutex);

  if (tlm.data_count >= TLM_BUFFER_SIZE) {
	tlm.data_count = 0;
	tlm.data_lost++;
  }

  memcpy(tlm.data + tlm.data_count, sample, sizeof(struct tlm_sample));

  tlm.data_count++;

  mutex_unlock(&tlmbuff_mutex);

}

static void tlm_buffer_reset(void)
{

  tlm.data_count = 0;
  tlm.data_lost = 0;

}

/*
 * Char device handling stuff
 */

static int dev_open(struct inode *inodep, struct file *filep)
{

  if (!mutex_trylock(&tlmsrv_mutex)) {
	pr_info("Device in use by another process");
	return -EBUSY;
  }

  return 0;

}

static ssize_t dev_read(struct file *filep, char *buffer, size_t len, loff_t *offset)
{

  int error_count = 0;

  /*
   * Specific tlmsrv should be read ONLY by appropriate tlm compiled along
   * with it. This check is not bulletproof, but it is the most that can be
   * done.
   */
  if (len != sizeof(struct tlm_private)) {
	pr_alert("Invalid data length requested. Use appropriate tlm for this tlmsrv.");
	return -EBADR;
  }

  mutex_lock(&tlmbuff_mutex);

  error_count = copy_to_user(buffer, &tlm, sizeof(struct tlm_private));

  tlm_buffer_reset();

  mutex_unlock(&tlmbuff_mutex);

  if (error_count == 0) {
	return 0;
  } else {
	pr_alert("Failed to provde ts data to user-space.");
	return -EFAULT;
  }

}

static ssize_t dev_write(struct file *filep, const char *buffer, size_t len, loff_t *offset)
{

  return -EACCES;

}

static int dev_release(struct inode *inodep, struct file *filep)
{
  mutex_unlock(&tlmsrv_mutex);
  return 0;
}

static struct file_operations fops =
{
  .open = dev_open,
  .read = dev_read,
  .write = dev_write,
  .release = dev_release,
};

/*
 * tlm life-cycle handling
 */

static int __init tlmsrv_init(void){

  int errno = 0;

  /* Try to dynamically allocate a major number for the device */
  device_major_number = register_chrdev(0, DEVICE_NAME, &fops);
  if (device_major_number < 0){
	pr_alert("Failed to register a major number\n");
	return device_major_number;
  }

  /* Register the device class */
  tlmsrv_class = class_create(THIS_MODULE, CLASS_NAME);
  if (IS_ERR(tlmsrv_class)) {
	errno = PTR_ERR(tlmsrv_class);
	pr_alert("Failed to register device class");
	goto out_clean_class;
  }

  /* Register the device driver */
  tlmsrv_device = device_create(tlmsrv_class, NULL, MKDEV(device_major_number, 0), NULL, DEVICE_NAME);
  if (IS_ERR(tlmsrv_device)){
	errno = PTR_ERR(tlmsrv_device);
	pr_alert("Failed to create the device");
	goto out_clean_device;
  }

  mutex_init(&tlmsrv_mutex);
  mutex_init(&tlmbuff_mutex);

  tlm_buffer_reset();

  /* Initialization successful */
  return 0;

out_clean_device:
  class_destroy(tlmsrv_class);
out_clean_class:
  unregister_chrdev(device_major_number, DEVICE_NAME);
  pr_alert("Failed to create the device\n");
  return errno;

}

static void __exit tlmsrv_exit(void){

  mutex_destroy(&tlmsrv_mutex);
  mutex_destroy(&tlmbuff_mutex);

  device_destroy(tlmsrv_class, MKDEV(device_major_number, 0));
  class_unregister(tlmsrv_class);
  class_destroy(tlmsrv_class);
  unregister_chrdev(device_major_number, DEVICE_NAME);

}
#endif
