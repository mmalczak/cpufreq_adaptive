// SPDX-License-Identifier: GPL-2.0-only
/*
 *	linux/drivers/cpufreq/cpufreq_adaptive.h
 *
 *	Header file for adaptive cpufreq governor template.
 *
 *	Copyright (C) 2020, NASK National Research Institute for Cybersecurity & AI
 *	Authors: 	Miłosz Malczak <milosz.malczak@nask.pl>
 *			Michał Karpowicz <michal.karpowicz@nask.pl>
 *			Michał Getka <michal.getka@nask.pl>

 *
 */


#include "cpufreq_governor.h"
#include <linux/string.h>
#include <linux/types.h>
/* In-kernel memory allocation */
#include <linux/slab.h>

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

#define d_M		((d_Rp+1)+(d_Sp+1))


#if d_Acl != d_Ao+d_Am
#error
#endif

#define buf_length (1<<3)
#define BUF_IDX_ADD(BUF_IDX, VALUE) ((BUF_IDX+(VALUE))&(buf_length-1))
#define PARAMS_FILTER_LENGTH 1

#define POINT_POS (32)
#define FP(X) ((X)*((int64_t)1<<(POINT_POS)))
#define FP_MAX ((int64_t)0x7fffffffffffffff)
#define FP_MIN ((int64_t)0x8000000000000000)

#define gov_show_one_adaptive(file_name)				\
static ssize_t show_##file_name						\
(struct gov_attr_set *attr_set, char *buf)				\
{									\
	struct dbs_data *dbs_data = to_dbs_data(attr_set);		\
	struct adaptive_dbs_tuners *tuners = dbs_data->tuners;		\
	return sprintf_fp(buf, tuners->file_name, '\n');		\
}

#define gov_show_vector_adaptive(file_name, length)			\
static ssize_t show_##file_name						\
(struct gov_attr_set *attr_set, char *buf)				\
{									\
	struct dbs_data *dbs_data = to_dbs_data(attr_set);		\
	struct adaptive_dbs_tuners *tuners = dbs_data->tuners;		\
	int i;								\
	int ret=0;							\
	for(i=0; i<length-1; i++)					\
		ret+=sprintf_fp(buf+ret, tuners->file_name[i], ' ');	\
	ret+=sprintf_fp(buf+ret, tuners->file_name[i], '\n');		\
	return ret;							\
}

#define gov_store_vector_adaptive(file_name, size)			\
static ssize_t store_##file_name					\
(struct gov_attr_set *attr_set, const char *buf, size_t count)		\
{									\
	struct dbs_data *dbs_data = to_dbs_data(attr_set);		\
	struct adaptive_dbs_tuners *tuners = dbs_data->tuners;		\
	int64_t input[size];						\
	int ret=0;							\
	int buf_idx=0;							\
	int i;								\
									\
	for(i=0; i<size; i++) {						\
		ret = sscanf_fp(buf, &input[i], &buf_idx);		\
		if (ret != 1) {						\
			return -EINVAL;					\
		}							\
	}								\
	for(i=0; i<size; i++)						\
		tuners->file_name[i] = input[i];			\
	return count;							\
}

struct adaptive_estimation_params
{
	int64_t P[deg*deg];
	int64_t theta[deg];		//value used for estimation
	int64_t theta_out[deg]; 	//filtered params
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

struct params_filter
{
	int idx;
	int64_t buf[(deg)*PARAMS_FILTER_LENGTH];
};

struct adaptive_policy_dbs_info {
	struct policy_dbs_info policy_dbs;
	int err;
	int64_t phi_P_phi;
	struct adaptive_estimation_params est_params;
	struct adaptive_controller_polynomials poly;
	struct adaptive_controller_buffers buf;
	struct params_filter p_filter;
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
