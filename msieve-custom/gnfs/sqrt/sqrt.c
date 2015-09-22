/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id: sqrt.c 904 2013-07-04 12:00:48Z jasonp_sf $
--------------------------------------------------------------------*/

#include "sqrt.h"

/* we will need to find primes q for which f(x) mod q
   is irreducible. This is the maximum number of q to try */

#define NUM_PRIME_RETRIES 100

/*--------------------------------------------------------------------*/
uint32 get_prime_for_sqrt(mpz_poly_t *alg_poly,
			  uint32 min_value,
			  uint32 *q_out) {

	uint32 i;
	uint32 status = 0;
	uint32 q = 0;
	prime_sieve_t prime_sieve;

	init_prime_sieve(&prime_sieve, min_value, (uint32)(-1));
	for (i = 0; i < NUM_PRIME_RETRIES; i++) {
		uint32 tmp_q = get_next_prime(&prime_sieve);
		if (is_irreducible(alg_poly, tmp_q)) {
			q = tmp_q;
			break;
		}
		else if (q == 0) {
			/* in rare cases, there is no q for which alg_poly
			   mod q is irreducible. Life becomes much more
			   difficult in this case, since square roots mod
			   q will not be unique. The only alternative when
			   this happens is to pick a q such that alg_poly 
			   mod q has no linear roots (so that all of the
			   relations mod q are relatively prime to alg_poly), 
			   then keep trying dependencies until by luck we 
			   start with the correct initial square root for 
			   the Newton iteration to succeed.

			   Buhler et. al. show that for polynomial degree d,
			   on average one in 2^(d/2) dependencies will lead to
			   a congruence of squares (and about half of those
			   will lead to a factor). Technically we also need to
			   check that alg_poly mod q is squarefree, but that
			   would require a full polynomial factoring routine;
			   I'm gambling that being squarefree is not rare. */

			uint32 roots[MAX_POLY_DEGREE];
			uint32 high_coeff;
			uint32 num_roots = poly_get_zeros(roots, alg_poly,
						tmp_q, &high_coeff, 1);
			if (high_coeff != 0 && num_roots == 0)
				q = tmp_q;
		}
	}
	free_prime_sieve(&prime_sieve);

	if (i == NUM_PRIME_RETRIES)
		status = 1;

	*q_out = q;
	return status;
}

/*--------------------------------------------------------------------*/
static void eval_poly_derivative(mpz_poly_t *poly, 
				mpz_t m1, mpz_t m0, 
				mpz_t n, mpz_t res) {
	uint32 i;
	mpz_t m1_pow;
	mpz_t tmp;

	mpz_init_set(m1_pow, m1);
	mpz_init(tmp);

	i = poly->degree;
	mpz_mul_ui(res, poly->coeff[i], i);

	mpz_neg(m0, m0);
	while (--i > 1) {
		mpz_mul(res, res, m0);
		mpz_mul_ui(tmp, poly->coeff[i], i);
		mpz_addmul(res, tmp, m1_pow);
		mpz_mul(m1_pow, m1_pow, m1);
	}
	mpz_mul(res, res, m0);
	mpz_addmul(res, poly->coeff[i], m1_pow);
	mpz_neg(m0, m0);

	mpz_clear(m1_pow);
	mpz_clear(tmp);
	mpz_mod(res, res, n);
}

/*--------------------------------------------------------------------*/
typedef struct {
	uint64 p;
	uint64 count;
} rat_prime_t;

static uint32 rat_square_root(relation_t *rlist, uint32 num_relations,
				mpz_t n, mpz_t sqrt_r) {
	uint32 i, j, num_primes;
	hashtable_t h;
	uint32 already_seen;
	uint32 array_size;
	mpz_t base, exponent, tmp;
	uint32 status = 0;
	rat_prime_t *curr;

	mpz_init(base);
	mpz_init(exponent);
	mpz_init(tmp);

	/* count up the number of times each prime factor in
	   rlist occurs */

	hashtable_init(&h, (uint32)WORDS_IN(rat_prime_t), 
				(uint32)WORDS_IN(uint64));

	for (i = 0; i < num_relations; i++) {
		relation_t *r = rlist + i;

		for (j = array_size = 0; j < r->num_factors_r; j++) {
			uint64 p = decompress_p(r->factors, &array_size);
			curr = (rat_prime_t *)hashtable_find(&h, &p, NULL,
							    &already_seen);
			if (!already_seen)
				curr->count = 1;
			else
				curr->count++;
		}
	}

	/* verify all such counts are even, and form the 
	   rational square root */

	mpz_set_ui(sqrt_r, 1);
	num_primes = hashtable_get_num(&h);
	curr = hashtable_get_first(&h);

	for (i = 0; i < num_primes; i++) {
		uint64 p = curr->p;
		uint64 count = curr->count;

		if (count % 2) {
			status = 1;
			break;
		}
		if (p > 0 && count > 0) {
			uint64_2gmp(p, base);
			uint64_2gmp(count / 2, exponent);
			mpz_powm(tmp, base, exponent, n);
			mpz_mul(sqrt_r, sqrt_r, tmp);
			mpz_tdiv_r(sqrt_r, sqrt_r, n);
		}
		curr = hashtable_get_next(&h, curr);
	}

	hashtable_free(&h);
	mpz_clear(base);
	mpz_clear(exponent);
	mpz_clear(tmp);
	return status;
}

/*--------------------------------------------------------------------*/
/* we will not do any computations involving the count,
   only verifying that it is even. Thus we can get away
   with storing only the low word of the count */

typedef struct {
	ideal_t ideal;
	uint32 count;
} alg_prime_t;

static uint32 verify_alg_ideal_powers(relation_t *rlist, 
					uint32 num_relations,
					uint32 *num_free_relations) {

	uint32 i, j, num_ideals;
	hashtable_t h;
	uint32 already_seen;
	alg_prime_t *curr;
	uint32 status = 0;

	/* count the multiplicity of each algebraic ideal (not
	   just the prime to which the ideal corresponds) in rlist */

	*num_free_relations = 0;

	hashtable_init(&h, (uint32)WORDS_IN(alg_prime_t),
			(uint32)WORDS_IN(ideal_t));

	for (i = 0; i < num_relations; i++) {
		relation_t *r = rlist + i;
		relation_lp_t rlp;

		find_large_ideals(r, &rlp, 0, 0);

		for (j = 0; j < rlp.ideal_count; j++) {
			ideal_t *curr_ideal = rlp.ideal_list + j;

			if (curr_ideal->rat_or_alg == RATIONAL_IDEAL)
				continue;

			curr = (alg_prime_t *)hashtable_find(&h, curr_ideal, 
						NULL, &already_seen);

			if (!already_seen)
				curr->count = 1;
			else
				curr->count++;
		}

		if (r->b == 0)
			(*num_free_relations)++;
	}

	/* verify each ideal occurs an even number of times */

	num_ideals = hashtable_get_num(&h);
	curr = hashtable_get_first(&h);

	for (i = 0; i < num_ideals; i++) {
		if (curr->count % 2) {
			status = 1;
			break;
		}
		curr = hashtable_get_next(&h, curr);
	}

	hashtable_free(&h);
	return status;
}

/*--------------------------------------------------------------------*/
uint32 nfs_find_factors(msieve_obj *obj, mpz_t n, 
			factor_list_t *factor_list) {

	/* external interface for the NFS square root */

	uint32 i, j;
	uint32 check_q;
	factor_base_t fb;
	mpz_poly_t monic_alg_poly;
	mpz_poly_t *rpoly;
	mpz_poly_t *apoly;
	mpz_t exponent, sqrt_r, sqrt_a;
	mpz_t c, tmp1, tmp2;
	uint32 dep_lower = 1;
	uint32 dep_upper = 64;
	uint32 factor_found = 0;
	time_t cpu_time;

	logprintf(obj, "\n");
	logprintf(obj, "commencing square root phase\n");

	memset(&fb, 0, sizeof(fb));
	apoly = &fb.afb.poly;
	rpoly = &fb.rfb.poly;
	mpz_poly_init(rpoly);
	mpz_poly_init(apoly);
	mpz_poly_init(&monic_alg_poly);
	mpz_init(exponent);
	mpz_init(sqrt_r);
	mpz_init(sqrt_a);
	mpz_init(c);
	mpz_init(tmp1);
	mpz_init(tmp2);

	/* read in the NFS polynomials */

	cpu_time = time(NULL);
	if (read_poly(obj, n, rpoly, apoly, NULL)) {
		logprintf(obj, "polynomials not found\n");
		goto finished;
	}

	/* find the values needed to convert the algebraic 
	   square root back to an integer */

	if (rpoly->degree != 1) {
		logprintf(obj, "cannot handle non-linear polynomials\n");
		goto finished;
	}

	/* construct a monic version of the algebraic poly,
	   saving off the leading coefficient separately */

	j = apoly->degree;
	if (mpz_cmp_ui(apoly->coeff[j], 0) < 0) {
		logprintf(obj, "cannot handle negative leading "
				"algebraic polynomial coefficient\n");
		goto finished;
	}

	mpz_set(c, apoly->coeff[j]);
	mpz_set(tmp1, c);
	mpz_set(monic_alg_poly.coeff[j-1], apoly->coeff[j-1]);
	monic_alg_poly.degree = j;
	mpz_set_ui(monic_alg_poly.coeff[j], 1);

	for (i = j - 2; (int32)i >= 0; i--) {
		mpz_mul(monic_alg_poly.coeff[i], apoly->coeff[i], tmp1);
		if (i > 0)
			mpz_mul(tmp1, tmp1, c);
	}
	get_prime_for_sqrt(&monic_alg_poly, (uint32)0x80000000, &check_q);

	/* determine the list of dependencies to compute */

	if (obj->nfs_args != NULL) {

		const char *tmp;
		const char *lower_limit;
		const char *upper_limit;

		tmp = strstr(obj->nfs_args, "dep_first=");
		if (tmp != NULL)
			dep_lower = strtoul(tmp + 10, NULL, 10);

		tmp = strstr(obj->nfs_args, "dep_last=");
		if (tmp != NULL)
			dep_upper = strtoul(tmp + 9, NULL, 10);

		/* old-style 'X,Y' format */

		upper_limit = strchr(obj->nfs_args, ',');
		if (upper_limit != NULL) {
			lower_limit = upper_limit - 1;
			while (lower_limit > obj->nfs_args &&
				isdigit(lower_limit[-1])) {
				lower_limit--;
			}
			upper_limit++;
			dep_lower = strtoul(lower_limit, NULL, 10);
			dep_upper = strtoul(upper_limit, NULL, 10);
		}

		dep_lower = MAX(dep_lower, 1);
		dep_upper = MAX(dep_upper, 1);
		dep_lower = MIN(dep_lower, 64);
		dep_upper = MIN(dep_upper, 64);
		logprintf(obj, "handling dependencies %u to %i\n",
				dep_lower, dep_upper);
	}

	/* for each dependency */

	for (i = dep_lower; i <= dep_upper; i++) {

		uint32 num_relations;
		uint32 num_free_relations;
		relation_t *rlist;
		abpair_t *abpairs;

		logprintf(obj, "reading relations for dependency %u\n", i);

		/* read in only the relations for dependency i */

		nfs_read_cycles(obj, &fb, NULL, NULL,
				&num_relations, &rlist, 0, i);

		if (num_relations == 0)
			continue;

		/* do some sanity checking, performing increasing
		   amounts of work as the dependency proves itself
		   to be valid */

		if (num_relations % 2) {
			/* the LA is supposed to force the number of 
			   relations in the dependency to be even. 
			   This isn't necessary if the leading coeff of
			   both NFS polynomials are squares, or if both
			   NFS polynomials are monic, since the 
			   corrections below that need the number of 
			   relations are avoided. But only a small 
			   minority of NFS jobs would satisfy this condition */

			logprintf(obj, "number of relations is not even\n");
			nfs_free_relation_list(rlist, num_relations);
			continue;
		}
		if (verify_alg_ideal_powers(rlist, 
				num_relations, &num_free_relations) != 0) {
			logprintf(obj, "algebraic side is not a square!\n");
			nfs_free_relation_list(rlist, num_relations);
			continue;
		}
		if (num_free_relations % 2) {
			logprintf(obj, "number of free relations (%u) is "
					"not even\n", num_free_relations);
			nfs_free_relation_list(rlist, num_relations);
			continue;
		}
		if (rat_square_root(rlist, num_relations, n, sqrt_r) != 0) {
			logprintf(obj, "rational side is not a square!\n");
			nfs_free_relation_list(rlist, num_relations);
			continue;
		}

		/* flatten the list of relations; each occurrence of
		   a relation gets its own abpair_t */

		abpairs = (abpair_t *)xmalloc(num_relations *
						sizeof(abpair_t));
		for (j = 0; j < num_relations; j++) {
			abpairs[j].a = rlist[j].a;
			abpairs[j].b = rlist[j].b;
		}
		nfs_free_relation_list(rlist, num_relations);

		/* perform the major work: the algebraic square root.
		   Note that to conserve memory, abpairs is freed in
		   the following call */

		mpz_set_ui(sqrt_a, 0);
		alg_square_root(obj, &monic_alg_poly, n, c, 
				rpoly->coeff[1], rpoly->coeff[0], 
				abpairs, num_relations, check_q, sqrt_a);
		if (mpz_sgn(sqrt_a) == 0) {
			logprintf(obj, "algebraic square root failed\n");
			continue;
		}

		/* an algebraic square root is available; move on
		   to the final congruence of squares. The arithmetic
		   is as given in Buhler et. al. with one exception:
		   when the rational poly is nonmonic there is a 
		   correction to the final square root value but the 
		   free relations *do not* figure into it. This latter
		   point is completely ignored in the literature! */

		eval_poly_derivative(apoly, rpoly->coeff[1], 
					rpoly->coeff[0], n, tmp1);
		mpz_mul(sqrt_r, sqrt_r, tmp1);
		mpz_mod(sqrt_r, sqrt_r, n);

		mpz_set_ui(exponent, 0);
		if (mpz_cmp_ui(c, 1) != 0) {
			mpz_set_ui(exponent, num_relations / 2 + 
						apoly->degree - 2);
			mpz_powm(tmp1, c, exponent, n);
			mpz_mul(sqrt_r, sqrt_r, tmp1);
			mpz_mod(sqrt_r, sqrt_r, n);
		}

		if (mpz_cmp_ui(rpoly->coeff[1], 1) != 0) {
			mpz_set_ui(exponent, (num_relations -
				       		num_free_relations) / 2);
			mpz_set(tmp1, rpoly->coeff[1]);
			if (mpz_sgn(tmp1) < 0)
				mpz_add(tmp1, tmp1, n);

			mpz_powm(tmp2, tmp1, exponent, n);
			mpz_mul(sqrt_a, sqrt_a, tmp2);
			mpz_mod(sqrt_a, sqrt_a, n);
		}

		/* a final sanity check: square the rational and algebraic 
		   square roots, expecting the same value modulo n */

		mpz_mul(tmp1, sqrt_r, sqrt_r);
		mpz_mul(tmp2, sqrt_a, sqrt_a);
		mpz_mod(tmp1, tmp1, n);
		mpz_mod(tmp2, tmp2, n);
		if (mpz_cmp(tmp1, tmp2) != 0) {
			logprintf(obj, "dependency does not form a "
					"congruence of squares!\n");
			continue;
		}

		/* look for a nontrivial factor of n */

		mpz_add(tmp1, sqrt_r, sqrt_a);
		mpz_gcd(tmp1, tmp1, n);
		if (mpz_cmp_ui(tmp1, 1) == 0) {
			logprintf(obj, "GCD is 1, no factor found\n");
		}
		else if (mpz_cmp(tmp1, n) == 0) {
			logprintf(obj, "GCD is N, no factor found\n");
		}
		else {
			/* factor found; add it to the list of factors. 
			   Stop trying dependencies if the remaining
			   composite is small enough that another method
			   will factor it faster.

			   Actually, we should be stopping when the remaining
			   composite is much larger (70-80 digits), but 
			   avoid doing this because the MPQS code will run
			   and wipe out all the NFS relations we've collected */

			uint32 composite_bits;
			mp_t junk;

			gmp2mp(tmp1, &junk);
			composite_bits = factor_list_add(obj, 
						factor_list, &junk);

			factor_found = 1;
			if (composite_bits < SMALL_COMPOSITE_CUTOFF_BITS) {
				break;
			}
			else {
				/* a single dependency could take hours,
				   and if N has more than two factors then
				   we'll need several dependencies to find
				   them all. So at least report the smallest
				   cofactor that we just found */

				mpz_divexact(tmp2, n, tmp1);
				gmp_sprintf(obj->mp_sprintf_buf, "%Zd",
						(mpz_cmp(tmp1, tmp2) < 0) ? 
						tmp1 : tmp2);
				logprintf(obj, "found factor: %s\n",
						obj->mp_sprintf_buf);
			}
		}
	}

finished:
	cpu_time = time(NULL) - cpu_time;
	logprintf(obj, "sqrtTime: %u\n", (uint32)cpu_time);

	mpz_poly_free(&fb.rfb.poly);
	mpz_poly_free(&fb.afb.poly);
	mpz_poly_free(&monic_alg_poly);
	mpz_clear(exponent);
	mpz_clear(sqrt_r);
	mpz_clear(sqrt_a);
	mpz_clear(c);
	mpz_clear(tmp1);
	mpz_clear(tmp2);

	return factor_found;
}

/*------------------------------------------------------------------

Modification to msieve version 1.52

Method: sqrt_fb_light_init()

This is performs a light initialization of the mpz_poly structs
in the factor_base_t.

Used only in sqrt_data_init()

-------------------------------------------------------------------*/

void sqrt_fb_light_init (factor_base_t *fb) {
	memset(fb, 0, sizeof(factor_base_t));
	mpz_poly_init(&(fb->rfb.poly));
	mpz_poly_init(&(fb->afb.poly));
}

/*------------------------------------------------------------------

Modification to msieve version 1.52

Method: sqrt_fb_light_copy_a_r_poly()

This is performs a light copy of the mpz_poly structs
in the factor_base_t.

Used only in sqrt_data_deep_copy()

-------------------------------------------------------------------*/

void sqrt_fb_light_copy_a_r_poly (factor_base_t *dst, factor_base_t *src) {
	mpz_poly_copy(&(dst->rfb.poly), &(src->rfb.poly));
	mpz_poly_copy(&(dst->afb.poly), &(src->afb.poly));
}

/*------------------------------------------------------------------

Modification to msieve version 1.52

Method: sqrt_data_init()

This is performs a initialization of the sqrt_data struct.

Used in sqrt_data_deep_copy() and nfs_find_factors_threaded().

-------------------------------------------------------------------*/

void sqrt_data_init (sqrt_data *dat) {
	sqrt_fb_light_init(&(dat->fb));
	dat->apoly = &(dat->fb.afb.poly);
	dat->rpoly = &(dat->fb.rfb.poly);
	mpz_poly_init(&(dat->monic_alg_poly));
	mpz_init(dat->exponent);
	mpz_init(dat->sqrt_r);
	mpz_init(dat->sqrt_a);
	mpz_init(dat->c);
	mpz_init(dat->tmp1);
	mpz_init(dat->tmp2);
}

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

void sqrt_data_deep_copy (sqrt_data *dst, sqrt_data *src) {
	sqrt_data_init(dst);
	dst->check_q = src->check_q;
	sqrt_fb_light_copy_a_r_poly(&(dst->fb), &(src->fb));
	mpz_poly_copy(&(dst->monic_alg_poly), &(src->monic_alg_poly));
	mpz_set(dst->exponent, src->exponent);
	mpz_set(dst->sqrt_r, src->sqrt_r);
	mpz_set(dst->sqrt_a, src->sqrt_a);
	mpz_set(dst->c, src->c);
	mpz_set(dst->tmp1, src->tmp1);
	mpz_set(dst->tmp2, src->tmp2);
}

/*------------------------------------------------------------------

Modification to msieve version 1.52

Method: free_sqrt_data()

Frees sqrt_data.

Used in nfs_find_factors_threaded() and sqrt_thread_shutdown()

-------------------------------------------------------------------*/

void free_sqrt_data (sqrt_data *dat) {
	mpz_poly_free(&(dat->fb.rfb.poly));
	mpz_poly_free(&(dat->fb.afb.poly));
	mpz_poly_free(&(dat->monic_alg_poly));
	mpz_clear(dat->exponent);
	mpz_clear(dat->sqrt_r);
	mpz_clear(dat->sqrt_a);
	mpz_clear(dat->c);
	mpz_clear(dat->tmp1);
	mpz_clear(dat->tmp2);

	free(dat);
}

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
				pthread_cond_t *status_cond, int *status, 
				pthread_mutex_t *count_mutex, int *count,
				mpz_t n, factor_list_t *factor_list) {

	thread_dat->dat = (sqrt_data *)xcalloc((size_t) 1, sizeof(sqrt_data));
	sqrt_data_deep_copy(thread_dat->dat, dat);
	thread_dat->dep_no = rl->dep_no;
	thread_dat->rlist = rl->rlist;
	thread_dat->num_relations = rl->num_relations;
	thread_dat->factor_found_mutex = factor_found_mutex;
	thread_dat->status_mutex = status_mutex;
	thread_dat->status_cond = status_cond;
	thread_dat->count_mutex = count_mutex;
	thread_dat->count = count;
	thread_dat->status = status;
	thread_dat->obj = obj;
	mpz_init(thread_dat->n);
	mpz_set(thread_dat->n, n);
	thread_dat->factor_list = factor_list;
}

/*------------------------------------------------------------------

Modification to msieve version 1.52

Method: free_sqrt_thread_data()

Frees the thread data object upon completion of the thread.

Used in nfs_find_factors_threaded() and sqrt_thread_shutdown().

-------------------------------------------------------------------*/

void free_sqrt_thread_data(sqrt_thread_data *data) {
	mpz_clear(data->n);
	free_sqrt_data(data->dat);
}

/*------------------------------------------------------------------

Modification to msieve version 1.52

Method: sqrt_thread_shutdown()

Performs necessary clean up of thread data after the task is complete.

Used in nfs_find_factors_threaded().

-------------------------------------------------------------------*/

void sqrt_thread_shutdown(void *arg, int thread_num) {
	sqrt_thread_data *data = (sqrt_thread_data *) arg;

	logprintf(data->obj, "Sqrt: Thread %u shutting down dep_task %u\n", 
					thread_num, data->dep_no);

	/* Update the number of completed dependencies */
	if (pthread_mutex_lock(data->count_mutex)) {
		perror("pthread_mutex_lock...");
	}
	*(data->count) = *(data->count) + 1;
	if (pthread_mutex_unlock(data->count_mutex)) {
		perror("pthread_mutex_unlock...");
	}

	/* Let the master check if we are done */
	if (pthread_mutex_lock(data->status_mutex)) {
		perror("pthread_mutex_lock...");
	}
	pthread_cond_broadcast(data->status_cond);
	if (pthread_mutex_unlock(data->status_mutex)) {
		perror("pthread_mutex_unlock...");
	}

	free_sqrt_thread_data(data);
}

/*------------------------------------------------------------------

Modification to msieve version 1.52

Method: sqrt_thread_find_factor()

Method that each thread runs to compute each dependency.

-------------------------------------------------------------------*/

void sqrt_thread_find_factor(void *data, int thread_num) {
	
	int i;
	sqrt_thread_data *init = (sqrt_thread_data *) data;
	msieve_obj *obj = init->obj;
	sqrt_data *dat = init->dat;
	pthread_mutex_t *factor_found_mutex = init->factor_found_mutex;
	pthread_mutex_t *status_mutex = init->status_mutex;
	int *status = init->status;
	uint32 num_relations = init->num_relations;
	uint32 num_free_relations;
	uint32 dep_no = init->dep_no;
	relation_t *rlist = init->rlist;
	abpair_t *abpairs = init->abpairs;
	factor_list_t *factor_list = init->factor_list;

	logprintf(obj, "Sqrt: Thread %u working on dependency %u\n", 
			  thread_num, dep_no);

	if (num_relations == 0)
		return;

	/* do some sanity checking, performing increasing
	   amounts of work as the dependency proves itself
	   to be valid */

	if (num_relations % 2) {
		/* the LA is supposed to force the number of 
		   relations in the dependency to be even. 
		   This isn't necessary if the leading coeff of
		   both NFS polynomials are squares, or if both
		   NFS polynomials are monic, since the 
		   corrections below that need the number of 
		   relations are avoided. But only a small 
		   minority of NFS jobs would satisfy this condition */

		logprintf(obj, "Sqrt: Dependency %u failed: number of relations" 
			      "is not even\n", dep_no);
		return;
	}


	if (verify_alg_ideal_powers(rlist, 
			num_relations, &num_free_relations) != 0) {
		logprintf(obj, "Sqrt: Dependency %u failed: algebraic side is not "
				  "a square!\n", dep_no);
		return;
	}
	if (num_free_relations % 2) {
		logprintf(obj, "Sqrt: Dependency %u failed: number of free relations"
		   		  " (%u) is not even\n", dep_no, num_free_relations);
		return;
	}
	if (rat_square_root(rlist, num_relations, init->n, dat->sqrt_r) != 0) {
		logprintf(obj, "Sqrt: Dependency %u failed: rational side is not "
			      "a square!\n", dep_no);
		return;
	}


	/* flatten the list of relations; each occurrence of
	   a relation gets its own abpair_t */

	abpairs = (abpair_t *)xmalloc(num_relations *
					sizeof(abpair_t));

	for (i = 0; i < num_relations; i++) {
		abpairs[i].a = rlist[i].a;
		abpairs[i].b = rlist[i].b;
	}

	nfs_free_relation_list(rlist, num_relations);

	/* perform the major work: the algebraic square root.
	   Note that to conserve memory, abpairs is freed in
	   the following call */
	logprintf(obj, "Sqrt: Thread %u beginning algebraic square root "
			  "for dependency %u\n", thread_num, dep_no);
	mpz_set_ui(dat->sqrt_a, 0);
	alg_square_root(obj, &(dat->monic_alg_poly), init->n, dat->c, 
							 dat->rpoly->coeff[1], dat->rpoly->coeff[0], 
							 abpairs, num_relations, dat->check_q, 
							 dat->sqrt_a);
	if (mpz_sgn(dat->sqrt_a) == 0) {
		logprintf(obj, "Sqrt: Thread %u: dependency %u, algebraic square "
			"root failed\n", thread_num, dep_no);
		return;
	}

	/* an algebraic square root is available; move on
	   to the final congruence of squares. The arithmetic
	   is as given in Buhler et. al. with one exception:
	   when the rational poly is nonmonic there is a 
	   correction to the final square root value but the 
	   free relations *do not* figure into it. This latter
	   point is completely ignored in the literature! */
	logprintf(obj, "Sqrt: Thread %u beginning congruence of squares "
			  "for dependency %u\n", thread_num, dep_no);

	eval_poly_derivative(dat->apoly, dat->rpoly->coeff[1], 
				dat->rpoly->coeff[0], init->n, dat->tmp1);
	mpz_mul(dat->sqrt_r, dat->sqrt_r, dat->tmp1);
	mpz_mod(dat->sqrt_r, dat->sqrt_r, init->n);

	mpz_set_ui(dat->exponent, 0);
	if (mpz_cmp_ui(dat->c, 1) != 0) {
		mpz_set_ui(dat->exponent, num_relations / 2 + 
					dat->apoly->degree - 2);
		mpz_powm(dat->tmp1, dat->c, dat->exponent, init->n);
		mpz_mul(dat->sqrt_r, dat->sqrt_r, dat->tmp1);
		mpz_mod(dat->sqrt_r, dat->sqrt_r, init->n);
	}

	if (mpz_cmp_ui(dat->rpoly->coeff[1], 1) != 0) {
		mpz_set_ui(dat->exponent, (num_relations -
			       		num_free_relations) / 2);
		mpz_set(dat->tmp1, dat->rpoly->coeff[1]);

		if (mpz_sgn(dat->tmp1) < 0)
			mpz_add(dat->tmp1, dat->tmp1, init->n);

		mpz_powm(dat->tmp2, dat->tmp1, dat->exponent, init->n);
		mpz_mul(dat->sqrt_a, dat->sqrt_a, dat->tmp2);
		mpz_mod(dat->sqrt_a, dat->sqrt_a, init->n);
	}

	/* a final sanity check: square the rational and algebraic 
	   square roots, expecting the same value modulo n */

	mpz_mul(dat->tmp1, dat->sqrt_r, dat->sqrt_r);
	mpz_mul(dat->tmp2, dat->sqrt_a, dat->sqrt_a);
	mpz_mod(dat->tmp1, dat->tmp1, init->n);
	mpz_mod(dat->tmp2, dat->tmp2, init->n);
	if (mpz_cmp(dat->tmp1, dat->tmp2) != 0) {
		logprintf(obj, "Sqrt: Thread %u: dependency %u does not form a "
				"congruence of squares!\n", thread_num, dep_no);
		return;
	}

	/* look for a nontrivial factor of n */

	mpz_add(dat->tmp1, dat->sqrt_r, dat->sqrt_a);
	mpz_gcd(dat->tmp1, dat->tmp1, init->n);

	if (mpz_cmp_ui(dat->tmp1, 1) == 0) {
		logprintf(obj, "Sqrt: Thread %u: dependency %u, "
			      "GCD is 1, no factor found\n",thread_num, dep_no);
	}
	else if (mpz_cmp(dat->tmp1, init->n) == 0) {
		logprintf(obj, "Sqrt: Thread %u: dependency %u, "
		          "GCD is N, no factor found\n", thread_num, dep_no);
	}
	else {
		/* factor found; add it to the list of factors. 
		   Stop trying dependencies if the remaining
		   composite is small enough that another method
		   will factor it faster.

		   Actually, we should be stopping when the remaining
		   composite is much larger (70-80 digits), but 
		   avoid doing this because the MPQS code will run
		   and wipe out all the NFS relations we've collected */

		/* If the factor has already been found from another thread,
		   stop calling factor_found to avoid confusing the rest of 
		   msieve code */

		if (pthread_mutex_lock(status_mutex)) {
			perror("pthread_mutex_lock for factor_found");
			return;
		}

		if (*status) {
			if (pthread_mutex_unlock(status_mutex)) {
				perror("pthread_mutex_lock for factor_found");
				return;
			}
			return;
		}

		if (pthread_mutex_unlock(status_mutex)) {
			perror("pthread_mutex_lock for factor_found");
			return;
		}

		/*  Lock out the factor_found portion to prevent 
			unspecified behavior from factor_list_add() 
			which is not known to be thread-safe */
		if (pthread_mutex_lock(factor_found_mutex)) {
			perror("pthread_mutex_lock for factor_found");
			return;
		}

		uint32 composite_bits;
		mp_t junk;

		gmp2mp(dat->tmp1, &junk);
		composite_bits = factor_list_add(obj, 
					factor_list, &junk);

		/* We are done if composite_bits condition is satisfied 
		   set status to 1 to prevent other threads from interfering */
		if (composite_bits < SMALL_COMPOSITE_CUTOFF_BITS) {
			if (pthread_mutex_lock(status_mutex)) {
				perror("pthread_mutex_lock for status_mutex");
				return;
			}

			*status = 1;
			
			if (pthread_mutex_unlock(status_mutex)) {
				perror("pthread_mutex_lock for status_mutex");
				return;
			}
		}
		else {
			/* a single dependency could take hours,
			   and if N has more than two factors then
			   we'll need several dependencies to find
			   them all. So at least report the smallest
			   cofactor that we just found */

			mpz_divexact(dat->tmp2, init->n, dat->tmp1);
			gmp_sprintf(obj->mp_sprintf_buf, "%Zd",
					(mpz_cmp(dat->tmp1, dat->tmp2) < 0) ? 
					dat->tmp1 : dat->tmp2);
			logprintf(obj, "found factor: %s\n",
					obj->mp_sprintf_buf);
		}

		if (pthread_mutex_unlock(factor_found_mutex)) {
			perror("pthread_mutex_unlock for factor found");
			return;
		}
	}
}


/*------------------------------------------------------------------

Modification to msieve version 1.52

Method: nfs_find_factors_threaded()

Implements a threaded version of nfs_find_factors(), where multiple
dependencies are worked on at the same time. This eliminates the 
non-deterministic runtime of nfs_find_factors(), where the time taken
would scale linearly with the number of "bad" dependencies it encounters.

Currently, all threads are given new tasks at the same time. This is
to allow all of the threads to finish first, and allow the master thread
to check if the factor has been found. 

A modification that could be made in the future to modify the threading
code to allow each thread to be killed at any time. This is a major
nightmare due to the fact that mutexes are involved and memory is freed
during the execution of thread, so I have not implemented it here.

As we must wait for all threads to complete, there is a ~10% overhead
here compared to the time it takes to work on one dependency only.

-------------------------------------------------------------------*/

uint32 nfs_find_factors_threaded(msieve_obj *obj, mpz_t n, 
			factor_list_t *factor_list) {

	/* external interface for the NFS square root */

	uint32 i, j;
	sqrt_data *dat; 
	uint32 *dep_lower;
	uint32 *dep_upper;
	time_t cpu_time;
	uint32 num_threads;
	thread_control_t control;
	relation_lists_t *rlists;
	struct threadpool *pool;
	pthread_mutex_t *status_mutex;
	pthread_mutex_t *factor_found_mutex;
	pthread_cond_t *status_cond;
	pthread_mutex_t *count_mutex;
	int *status;
	int *count;
	int result;
	sqrt_thread_data *thread_dat;

	logprintf(obj, "\n");
	logprintf(obj, "commencing square root phase\n");

	/* Allocate memory for sqrt_data, mutexes etc. */
	dat = (sqrt_data *)xcalloc((size_t) 1, sizeof(sqrt_data));
	sqrt_data_init(dat);

	status = (int *)xmalloc(sizeof(int));
	*status = 0;
	count = (int *)xmalloc(sizeof(int));
	*count = 0;
	dep_lower = (uint32 *)xmalloc(sizeof(uint32));
	*dep_lower = 1;
	dep_upper = (uint32 *)xmalloc(sizeof(uint32));
	*dep_upper = 64;
	factor_found_mutex = (pthread_mutex_t *)xmalloc(sizeof(pthread_mutex_t));
	status_mutex = (pthread_mutex_t *)xmalloc(sizeof(pthread_mutex_t));
	status_cond = (pthread_cond_t *)xmalloc(sizeof(pthread_cond_t));
	count_mutex = (pthread_mutex_t *)xmalloc(sizeof(pthread_mutex_t));

	/* read in the NFS polynomials */

	cpu_time = time(NULL);
	if (read_poly(obj, n, dat->rpoly, dat->apoly, NULL)) {
		logprintf(obj, "polynomials not found\n");
		goto finished;
	}

	/* find the values needed to convert the algebraic 
	   square root back to an integer */

	if (dat->rpoly->degree != 1) {
		logprintf(obj, "cannot handle non-linear polynomials\n");
		goto finished;
	}

	/* construct a monic version of the algebraic poly,
	   saving off the leading coefficient separately */

	j = dat->apoly->degree;
	if (mpz_cmp_ui(dat->apoly->coeff[j], 0) < 0) {
		logprintf(obj, "cannot handle negative leading "
				"algebraic polynomial coefficient\n");
		goto finished;
	}

	mpz_set(dat->c, dat->apoly->coeff[j]);
	mpz_set(dat->tmp1, dat->c);
	mpz_set(dat->monic_alg_poly.coeff[j-1], dat->apoly->coeff[j-1]);
	dat->monic_alg_poly.degree = j;
	mpz_set_ui(dat->monic_alg_poly.coeff[j], 1);

	for (i = j - 2; (int32)i >= 0; i--) {
		mpz_mul(dat->monic_alg_poly.coeff[i], dat->apoly->coeff[i], dat->tmp1);
		if (i > 0)
			mpz_mul(dat->tmp1, dat->tmp1, dat->c);
	}
	get_prime_for_sqrt(&(dat->monic_alg_poly), (uint32)0x80000000, 
		               &(dat->check_q));

	/* determine the list of dependencies to compute */

	if (obj->nfs_args != NULL) {

		const char *tmp;
		const char *lower_limit;
		const char *upper_limit;

		tmp = strstr(obj->nfs_args, "dep_first=");
		if (tmp != NULL)
			*dep_lower = strtoul(tmp + 10, NULL, 10);

		tmp = strstr(obj->nfs_args, "dep_last=");
		if (tmp != NULL)
			*dep_upper = strtoul(tmp + 9, NULL, 10);

		/* old-style 'X,Y' format */

		upper_limit = strchr(obj->nfs_args, ',');
		if (upper_limit != NULL) {
			lower_limit = upper_limit - 1;
			while (lower_limit > obj->nfs_args &&
				isdigit(lower_limit[-1])) {
				lower_limit--;
			}
			upper_limit++;
			*dep_lower = strtoul(lower_limit, NULL, 10);
			*dep_upper = strtoul(upper_limit, NULL, 10);
		}

		*dep_lower = MAX(*dep_lower, 1);
		*dep_upper = MAX(*dep_upper, 1);
		*dep_lower = MIN(*dep_lower, 64);
		*dep_upper = MIN(*dep_upper, 64);
		logprintf(obj, "Sqrt: Handling dependencies %u to %i\n",
				*dep_lower, *dep_upper);
	}

	logprintf(obj,"Sqrt: getting relations for each dependency\n");
	/* Grab all of the relations for each dependency */
	nfs_read_cycles_threaded(obj, &(dat->fb), &rlists, dep_lower, dep_upper);

	if (rlists == NULL) {
		logprintf(obj, "Sqrt: Could not process dependency cycle relations.\n");
		goto finished;
	}

	/* Time to distribute the tasks to the threads! */

	/* Initialize threadpool and mutexes */
	logprintf(obj,"Sqrt: initializing threads\n");
	num_threads = MIN(obj->num_threads, *dep_upper - *dep_lower + 1);
	logprintf(obj, "Sqrt using %u threads\n", num_threads);
	control.init = NULL;
    control.shutdown = NULL;
    control.data = NULL;
	pool = threadpool_init(num_threads, (*dep_upper - *dep_lower + 1) * 2,
						   &control);

	if (pthread_mutex_init(factor_found_mutex, NULL)) {
		perror("pthread mutex init...");
		goto finished;
	}
	if (pthread_mutex_init(status_mutex, NULL)) {
		perror("pthread mutex init...");
		goto finished;
	}
	if (pthread_cond_init(status_cond, NULL)) {
		perror("pthread cond init...");
		goto finished;
	}
	if (pthread_mutex_init(count_mutex, NULL)) {
		perror("pthread mutex init...");
		goto finished;
	}

	/* Prepare the thread data needed for each dependency task */
	thread_dat = (sqrt_thread_data *)xcalloc((size_t) (*dep_upper - 
											 		   *dep_lower + 1), 
									         sizeof(sqrt_thread_data));

	for (i = *dep_lower; i <= *dep_upper; i++) {
		relation_lists_t *rl = rlists + i - *dep_lower;
		sqrt_thread_data *data = thread_dat + i - *dep_lower;

		sqrt_thread_data_init(obj, data, dat, rl, 
				factor_found_mutex, status_mutex, 
				status_cond, status, count_mutex, 
				count, n, factor_list);
	}

	if (pthread_mutex_lock(status_mutex)) {
		perror("pthread mutex lock...");
		goto finished;
	}

	while (*status == 0) {
		if (pthread_mutex_lock(count_mutex)) {
			perror("pthread_mutex_lock");
			goto finished;
		}

		/* Add tasks in a batch fashion to avoid unusual thread conditions 

		   Tthis could be overhauled in the future if threadpool.c 
		   and sqrt_thread_find_factor supports task cancellation allowing
		   for all of the tasks to be added at the same time */

		if (*count == *dep_upper - *dep_lower + 1) {
			logprintf(obj, "Sqrt: All dependencies have been tried and have "
					  "failed. Exiting...\n");
			*status = -1;
		}

		/* If the batch is completed, add more tasks */
		if (*count % num_threads == 0) {
			uint32 start_dep = MIN(*dep_lower + *count, *dep_upper);
			uint32 end_dep = MIN(start_dep + num_threads - 1, *dep_upper);
			logprintf(obj, "Sqrt: Adding dependencies %u to %u "
			 		  "to the task pool\n", start_dep, end_dep);

			for (i = start_dep; i <= end_dep; i++) {
				sqrt_thread_data *data = thread_dat + i - *dep_lower;
				task_control_t t;
				t.init = NULL;
				t.run = sqrt_thread_find_factor;
				t.data = data;
				t.shutdown = sqrt_thread_shutdown;
				threadpool_add_task(pool, &t, 1);
			}
		}

		if (pthread_mutex_unlock(count_mutex)) {
			perror("pthread_mutex_lock");
			goto finished;
		}
		pthread_cond_wait(status_cond, status_mutex);
	}

	if (pthread_mutex_unlock(status_mutex)) {
		perror("pthread mutex unlock...");
		goto finished;
	}

	/* freeing up the threadpool and remaining unused dep data */
	threadpool_free(pool);
	for (i = *dep_lower + *count; i<=*dep_upper; i++) {
		free_sqrt_thread_data(thread_dat + i - *dep_lower);
	}
	free(thread_dat);

finished:
	result = *status == -1 ? 0 : 1;
	cpu_time = time(NULL) - cpu_time;
	
	free_sqrt_data(dat);
	pthread_mutex_destroy(status_mutex);
	pthread_mutex_destroy(factor_found_mutex);
	pthread_mutex_destroy(count_mutex);
	pthread_cond_destroy(status_cond);
	free(status_mutex);
	free(factor_found_mutex);
	free(count_mutex);
	free(status_cond);
	free(status);
	free(count);
	free(dep_lower);
	free(dep_upper);
	logprintf(obj, "sqrtTime: %u\n", (uint32)cpu_time);
	return result;
}
