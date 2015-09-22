/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: sqrt.h 638 2011-09-11 15:31:19Z jasonp_sf $
--------------------------------------------------------------------*/

#ifndef _GNFS_SQRT_SQRT_H_
#define _GNFS_SQRT_SQRT_H_

#include "gnfs.h"
#include <thread.h>

#ifdef __cplusplus
extern "C" {
#endif

/*------------------------------------------------------------------

Modification to msieve version 1.52

Method: sqrt_fb_light_init()

This is performs a light initialization of the mpz_poly structs
in the factor_base_t.

Used only in sqrt_data_init()

-------------------------------------------------------------------*/

void sqrt_fb_light_init (factor_base_t *fb);

/*------------------------------------------------------------------

Modification to msieve version 1.52

Method: sqrt_fb_light_copy_a_r_poly()

This is performs a light copy of the mpz_poly structs
in the factor_base_t.

Used only in sqrt_data_deep_copy()

-------------------------------------------------------------------*/

void sqrt_fb_light_copy_a_r_poly (factor_base_t *dst, factor_base_t *src);

/*------------------------------------------------------------------

Modification to msieve version 1.52

Struct: sqrt_data

Contains all of the data necessary for the sqrt stage, to be passed
to each thread.

-------------------------------------------------------------------*/

typedef struct {
	uint32 check_q;
	factor_base_t fb;
	mpz_poly_t monic_alg_poly;
	mpz_poly_t *rpoly;
	mpz_poly_t *apoly;
	mpz_t exponent, sqrt_r, sqrt_a;
	mpz_t c, tmp1, tmp2;
} sqrt_data;

/*------------------------------------------------------------------

Modification to msieve version 1.52

Method: sqrt_data_init()

This is performs a initialization of the sqrt_data struct.

Used in sqrt_data_deep_copy() and nfs_find_factors_threaded().

-------------------------------------------------------------------*/

void sqrt_data_init (sqrt_data *dat);

/*------------------------------------------------------------------

Modification to msieve version 1.52

Method: sqrt_data_deep_copy()

This is performs a deep copy of the sqrt_data struct. Only pass 
in an uninitialized, but memory allocated, sqrt_data pointer as
dst. The current method does not do the correct checking to
ensure memory-safety otherwise.

Used in nfs_find_factors_threaded() in order to make a copy of
sqrt_data for each of the threads.

-------------------------------------------------------------------*/

void sqrt_data_deep_copy (sqrt_data *dst, sqrt_data *src);

/*------------------------------------------------------------------

Modification to msieve version 1.52

Method: free_sqrt_data()

Frees sqrt_data.

Used in nfs_find_factors_threaded() and sqrt_thread_shutdown()

-------------------------------------------------------------------*/

void free_sqrt_data (sqrt_data *dat);

/*------------------------------------------------------------------

Modification to msieve version 1.52

Struct: sqrt_thread_data

Contains all of the data necessary for each thread to compute its
dependency.

-------------------------------------------------------------------*/

typedef struct {
	pthread_mutex_t *factor_found_mutex;
	pthread_mutex_t *status_mutex;
	pthread_cond_t *status_cond;
	pthread_mutex_t *count_mutex;
	int *status;
	int *count;
	uint32 num_relations;
	uint32 num_free_relations;
	relation_t *rlist;
	abpair_t *abpairs;
	sqrt_data *dat;
	uint32 dep_no;
	msieve_obj *obj;
	mpz_t n;
	factor_list_t *factor_list;

} sqrt_thread_data;

/*------------------------------------------------------------------

Modification to msieve version 1.52

Method: sqrt_thread_data_init()

Initializes the data object required by each thread.

Used in nfs_find_factors_threaded().

-------------------------------------------------------------------*/

void sqrt_thread_data_init(msieve_obj *obj, 
				sqrt_thread_data *thread_dat, 
				sqrt_data *dat, relation_lists_t *rl, 
				pthread_mutex_t *factor_found_mutex,
				pthread_mutex_t *status_mutex, 
				pthread_cond_t *status_cond, int *status, pthread_mutex_t *count_mutex, 
				int *count, mpz_t n, factor_list_t *factor_list);

/*------------------------------------------------------------------

Modification to msieve version 1.52

Method: free_sqrt_thread_data()

Frees the thread data object upon completion of the thread.

Used in nfs_find_factors_threaded() and sqrt_thread_shutdown().

-------------------------------------------------------------------*/

void free_sqrt_thread_data(sqrt_thread_data *data);

/*------------------------------------------------------------------

Modification to msieve version 1.52

Method: sqrt_thread_shutdown()

Performs necessary clean up of thread data after the task is complete.

Used in nfs_find_factors_threaded().

-------------------------------------------------------------------*/

void sqrt_thread_shutdown(void *arg, int thread_num);

uint32 get_prime_for_sqrt(mpz_poly_t *alg_poly,
			  uint32 min_value,
			  uint32 *q_out); 

void alg_square_root(msieve_obj *obj, mpz_poly_t *monic_alg_poly, 
			mpz_t n, mpz_t c, mpz_t m1, mpz_t m0,
			abpair_t *rlist, uint32 num_relations, 
			uint32 check_q, mpz_t sqrt_a);

#ifdef __cplusplus
}
#endif

#endif /* _GNFS_SQRT_SQRT_H_ */
