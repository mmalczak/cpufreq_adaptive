/*
 *	tlmt.h
 *
 * Datatypes definition. File common for LKM and user-space application.
 *
 *	Copyright (C)  2017 Michał Getka
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation.
 *
 */

#define DEVICE_NAME "adaptsrv"	   // Device path /dev/DEVICE_NAME

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#ifndef d_Ao
#define d_Ao 2
#define d_Am 3
#define d_Rd 1
#define d_Sd 1
#define d_A 2
#define d_B 1

#define d_Bplus (0)
#define d_Bmd (0)
#define d_Bminus (d_B)
#define d_R (d_Sd+d_Bminus+d_Rd+d_Bplus)
#define d_S (d_A + d_Rd-1 + d_Sd)
#define d_T (d_Bmd + d_Ao)
#define d_D MAX(d_R, d_Ao)
#define deg (d_A+d_B+1)
#endif

#ifdef __user_space__
#include <stdint.h>
#else
#include <linux/types.h>
#endif

#define __tlm_var_norm(var) {							\
	printf("%s%f", delimiter, (float)(sample->var)/((int64_t)1<<32));	\
}

#define __tlm_vector_norm(var,length) {							\
	printf("%s[", delimiter);							\
	for(i=0; i<length; i++)								\
	{										\
		if(i<length-1)								\
			printf("%.3f,", (float)(sample->var[i])/((int64_t)1<<32));	\
		else									\
			printf("%.3f", (float)(sample->var[i])/((int64_t)1<<32));	\
	}										\
	printf("]");									\
}


#define __tlm_var(var,fmt) {					\
	if (first_var) {					\
		first_var = 0;					\
		printf("%" fmt, sample->var);			\
	} else							\
		printf("%s%" fmt, delimiter, sample->var);	\
}

#define TLM_BUFFER_SIZE 4096

struct tlm_sample {

	long long ts_sec;
	long int ts_nsec;

	// Insert your timeservies signals here. Pointer data types are prohibited.
	int err;
	unsigned int load;
	int load_est;
	int64_t uc;
	int64_t u_prev;
	int64_t v;
	// Model vector
	int64_t theta[deg];
	// Controller polynomials
	int64_t R[d_R+1];
	int64_t S[d_S+1];
	int64_t T[d_T+1];
	int64_t D[d_D+1];
	int64_t P[deg*deg];
	int64_t phi_P_phi;
};

// Labels for timeseries signals
const char *labelsDef[] = {
	"err",
	"load",
	"load_est",
	"uc",
	"u_prev",
	"v",
	"theta",
	"R",
	"S",
	"T",
	"D",
	"P",
	"phi_P_phi",
};

// Only for user-space use
#ifndef MODULE

// Displaying timeseries data row to stdout
static inline void print_sample(struct tlm_sample *sample, char * delimiter, int long_ver)
{
	/*
	 * This function needs to print out appropriately formatted sample data
	 * arranged in the same way as they are specified in labels definition
	 * labelsDef.
	 */

	// For the use of __tlm_var macro
	unsigned int first_var = 1;
	int i;
	/*
	 * Macro __tlm_var(var, fmt) prints out single sample component.
	 *  var - sample variable name
	 *  fmt - printf format for the variable
	 */
	__tlm_var(err, "d");
	__tlm_var(load, "u");
	__tlm_var(load_est, "d");
	if(long_ver) printf("\t");
	__tlm_var_norm(uc);
	__tlm_var(u_prev, "ld");
	__tlm_var(v, "ld");
	if(long_ver && (sample->v<10000000)) printf("\t");
	__tlm_vector_norm(theta, deg);
	__tlm_vector_norm(R, d_R+1);
	__tlm_vector_norm(S, d_S+1);
	__tlm_vector_norm(T, d_T+1);
	__tlm_vector_norm(D, d_D+1);
	__tlm_vector_norm(P, deg*deg);
	__tlm_var_norm(phi_P_phi);
	printf(";");
}

#endif

struct tlm_private {
	// Number of available samples in data buffer
	unsigned int data_count;
	// Number of buffer reloads since last read action
	unsigned int data_lost;
	struct tlm_sample data[TLM_BUFFER_SIZE];
};

struct tlm_private tlm;
