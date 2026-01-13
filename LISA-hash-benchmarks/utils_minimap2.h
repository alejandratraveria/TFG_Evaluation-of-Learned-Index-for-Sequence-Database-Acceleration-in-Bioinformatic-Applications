#ifndef UTILS_MINIMAP2_H
#define UTILS_MINIMAP2_H

#include <stdint.h>
#include <stdio.h>
#include <fcntl.h>
#include <zlib.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <inttypes.h>

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;

#include "kseq.h"
#include "kvec.h"
#include "ksort.h"
#include "kdq.h"
#include "khash.h"
#include "metrics.h"
#include "ketopt.h"
#include "easyperf.h"
#include "perf_profiler.h"

#ifdef LISA_HASH
#include "lisa_hash.h"
extern lisa_hash<uint64_t, uint64_t> *lh;
#endif

#define idx_hash(a) ((a)>>1)
#define idx_eq(a, b) ((a)>>1 == (b)>>1)

KHASH_INIT(idx, uint64_t, uint64_t, 1, idx_hash, idx_eq)
KHASH_MAP_INIT_STR(str, uint32_t)
KDQ_INIT(int)
KSEQ_INIT2(, gzFile, gzread)
typedef khash_t(idx) idxhash_t;

#define kroundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))

#define mm_seq4_set(s, i, c) ((s)[(i)>>3] |= (uint32_t)(c) << (((i)&7)<<2))

// Options
typedef struct {
	short k, w, flag, bucket_bits;
	int64_t mini_batch_size;
	uint64_t batch_size;
} mm_idxopt_t;

typedef struct {
	int64_t flag;    // see MM_F_* macros
	int seed;
	int sdust_thres; // score threshold for SDUST; 0 to disable

	int max_qlen;    // max query length

	int bw, bw_long; // bandwidth
	int max_gap, max_gap_ref; // break a chain if there are no minimizers in a max_gap window
	int max_frag_len;
	int max_chain_skip, max_chain_iter;
	int min_cnt;         // min number of minimizers on each chain
	int min_chain_score; // min chaining score
	float chain_gap_scale;
	float chain_skip_scale;
	int rmq_size_cap, rmq_inner_dist;
	int rmq_rescue_size;
	float rmq_rescue_ratio;

	float mask_level;
	int mask_len;
	float pri_ratio;
	int best_n;      // top best_n chains are subjected to DP alignment

	float alt_drop;

	int a, b, q, e, q2, e2; // matching score, mismatch, gap-open and gap-ext penalties
	int sc_ambi; // score when one or both bases are "N"
	int noncan;      // cost of non-canonical splicing sites
	int junc_bonus;
	int zdrop, zdrop_inv;   // break alignment if alignment score drops too fast along the diagonal
	int end_bonus;
	int min_dp_max;  // drop an alignment if the score of the max scoring segment is below this threshold
	int min_ksw_len;
	int anchor_ext_len, anchor_ext_shift;
	float max_clip_ratio; // drop an alignment if BOTH ends are clipped above this ratio

	int rank_min_len;
	float rank_frac;

	int pe_ori, pe_bonus;

	float mid_occ_frac;  // only used by mm_mapopt_update(); see below
	float q_occ_frac;
	int32_t min_mid_occ, max_mid_occ;
	int32_t mid_occ;     // ignore seeds with occurrences above this threshold
	int32_t max_occ, max_max_occ, occ_dist;
	int64_t mini_batch_size; // size of a batch of query bases to process in parallel
	int64_t max_sw_mat;
	int64_t cap_kalloc;

	const char *split_prefix;
} mm_mapopt_t;

// minimap.h
#define MM_MAX_SEG       255

#define MM_F_ALL_CHAINS    (0x800000LL)
#define MM_F_NO_DIAG       (0x001LL) // no exact diagonal hit
#define MM_F_NO_DUAL       (0x002LL) // skip pairs where query name is lexicographically larger than target name
#define MM_F_NO_LJOIN      (0x400LL)
#define MM_I_HPC          0x1
#define MM_F_RMQ           (0x80000000LL)
#define MM_F_SR            (0x1000LL)
#define MM_F_FRAG_MODE     (0x2000LL)
#define MM_F_NO_PRINT_2ND  (0x4000LL)
#define MM_F_2_IO_THREADS  (0x8000LL)
#define MM_F_HEAP_SORT     (0x400000LL)
#define MM_F_SPLICE        (0x080LL) // splice mode
#define MM_F_SPLICE_FOR    (0x100LL) // match GT-AG
#define MM_F_SPLICE_REV    (0x200LL) // match CT-AC, the reverse complement of GT-AG
#define MM_F_SPLICE_FLANK  (0x40000LL)
#define MM_IDX_MAGIC   "MMI\2"
#define MM_I_NO_SEQ       0x2
#define MM_I_NO_NAME      0x4
#define MM_F_OUT_SAM       (0x008LL)
#define MM_F_NO_QUAL       (0x010LL)
#define MM_F_COPY_COMMENT  (0x2000000LL)


#define validate_int(x) (assert(INT_MIN <= x && x <= INT_MAX));
#define validate_uint(x) (assert(0 <= x && x <= UINT_MAX));
#define validate_ulong(x) (assert(0 <= x && x <= ULONG_MAX));

#define CHECK_PAIR_THRES 1000000

#define SD_WLEN 3
#define SD_WTOT (1<<(SD_WLEN<<1))
#define SD_WMSK (SD_WTOT - 1)

// Global variables
unsigned char seq_comp_table[256] = {
	  0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
	 16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
	 32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
	 48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
	 64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
	'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	 92,  93,  94,	95,
	 96, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
	'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127,
	128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
	144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
	160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
	176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
	192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
	208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
	224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
	240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255
};

unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

int mm_verbose = 3;
int mm_dbg_flag = 0;
double mm_realtime0;

typedef struct {
	int l_seq, rid;
	char *name, *seq, *qual, *comment;
} mm_bseq1_t;

struct mm_bseq_file_s {
	gzFile fp;
	kseq_t *ks;
	mm_bseq1_t s;
};

typedef struct mm_bseq_file_s mm_bseq_file_t;

// index reader
typedef struct {
	int is_idx, n_parts;
	int64_t idx_size;
	mm_idxopt_t opt;
	FILE *fp_out;
	union {
		struct mm_bseq_file_s *seq;
		FILE *idx;
	} fp;
} mm_idx_reader_t;

// minimap2 index
typedef struct {
	char *name;      // name of the db sequence
	uint64_t offset; // offset in mm_idx_t::S
	uint32_t len;    // length
	uint32_t is_alt;
} mm_idx_seq_t;

typedef struct {
	int32_t st, en, max; // max is not used for now
	int32_t score:30, strand:2;
} mm_idx_intv1_t;

typedef struct mm_idx_intv_s {
	int32_t n, m;
	mm_idx_intv1_t *a;
} mm_idx_intv_t;

typedef struct {
	int32_t b, w, k, flag;
	uint32_t n_seq;            // number of reference sequences
	int32_t index;
	int32_t n_alt;
	mm_idx_seq_t *seq;         // sequence name, length and offset
	uint32_t *S;               // 4-bit packed sequence
	struct mm_idx_bucket_s *B; // index (hidden)
	struct mm_idx_intv_s *I;   // intervals (hidden)
	void *km, *h;
} mm_idx_t;

class hash_entry {
	public: 
	uint64_t key;
	uint64_t n;
	uint64_t *p;
	hash_entry(uint64_t k, uint64_t n_, uint64_t *p_){
		key = k;
		n = n_;
		p = p_;
	}

};

bool key_sort( hash_entry i1, hash_entry i2)
{
    return (i1.key < i2.key);
}

// emulate 128-bit integers and arrays
typedef struct { uint64_t x, y; char seq[40]; } mm128_t;
//typedef struct { uint64_t x, y; } mm128_t;
typedef struct { size_t n, m; mm128_t *a; } mm128_v;

typedef struct mm_idx_bucket_s {
	mm128_v a;   // (minimizer, position) array
	int32_t n;   // size of the _p_ array
	uint64_t *p; // position array for minimizers appearing >1 times
	void *h;     // hash table indexing _p_ and minimizers appearing once
} mm_idx_bucket_t;

typedef struct header_t {
	size_t size;
	struct header_t *ptr;
} header_t;

typedef struct {
	void *par;
	size_t min_core_size;
	header_t base, *loop_head, *core_head; /* base is a zero-sized block always kept in the loop */
} kmem_t;

typedef struct {
	int mini_batch_size;
	uint64_t batch_size, sum_len;
	mm_bseq_file_t *fp;
	mm_idx_t *mi;
} pipeline_t;

typedef struct {
    int n_seq;
	mm_bseq1_t *seq;
	mm128_v a;
} step_t;

typedef struct { // a simplified version of kdq
	int front, count;
	int a[32];
} tiny_queue_t;

typedef struct {
	struct kt_for_t *t;
	long i;
} ktf_worker_t;

typedef struct kt_for_t {
	int n_threads;
	long n;
	ktf_worker_t *w;
	void (*func)(void*,long,int,long);
	void *data;
} kt_for_t;

#if (defined(WIN32) || defined(_WIN32)) && defined(_MSC_VER)
#define __sync_fetch_and_add(ptr, addend)     _InterlockedExchangeAdd((void*)ptr, addend)
#endif

typedef struct {
	struct ktp_t *pl;
	int64_t index;
	int step;
	void *data;
} ktp_worker_t;

typedef struct ktp_t {
	void *shared;
	void *(*func)(void*, int, void*);
	int64_t index;
	int n_workers, n_steps;
	ktp_worker_t *workers;
	pthread_mutex_t mutex;
	pthread_cond_t cv;
} ktp_t;

typedef struct {
	int start, finish;
	int r, l;
} perf_intv_t;

typedef kvec_t(perf_intv_t) perf_intv_v;
typedef kvec_t(uint64_t) uint64_v;

struct sdust_buf_s {
	kdq_t(int) *w;
	perf_intv_v P; // the list of perfect intervals for the current window, sorted by descending start and then by ascending finish
	uint64_v res;  // the result
	void *km;      // memory pool
};

typedef struct sdust_buf_s sdust_buf_t;

#define sort_key_128x(a) ((a).x)
KRADIX_SORT_INIT(128x, mm128_t, sort_key_128x, 8) 

#define sort_key_64(x) (x)
KRADIX_SORT_INIT(64, uint64_t, sort_key_64, 8)

KSORT_INIT_GENERIC(uint32_t)
KSORT_INIT_GENERIC(uint64_t)

// Set options ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void mm_idxopt_init(mm_idxopt_t *opt)
{
	memset(opt, 0, sizeof(mm_idxopt_t));
	opt->k = 20, opt->w = 10, opt->flag = 0;
	opt->bucket_bits = 14;
	opt->mini_batch_size = 50000000;
	opt->batch_size = 4000000000ULL;
}

void mm_mapopt_init(mm_mapopt_t *opt)
{
	memset(opt, 0, sizeof(mm_mapopt_t));
	opt->seed = 11;
	opt->mid_occ_frac = 2e-4f;
	opt->min_mid_occ = 10;
	opt->max_mid_occ = 1000000;
	opt->sdust_thres = 0; // no SDUST masking
	opt->q_occ_frac = 0.01f;

	opt->min_cnt = 3;
	opt->min_chain_score = 40;
	opt->bw = 500, opt->bw_long = 20000;
	opt->max_gap = 5000;
	opt->max_gap_ref = -1;
	opt->max_chain_skip = 25;
	opt->max_chain_iter = 5000;
	opt->rmq_inner_dist = 1000;
	opt->rmq_size_cap = 100000;
	opt->rmq_rescue_size = 1000;
	opt->rmq_rescue_ratio = 0.1f;
	opt->chain_gap_scale = 0.8f;
	opt->chain_skip_scale = 0.0f;
	opt->max_max_occ = 4095;
	opt->occ_dist = 500;

	opt->mask_level = 0.5f;
	opt->mask_len = INT_MAX;
	opt->pri_ratio = 0.8f;
	opt->best_n = 5;

	opt->alt_drop = 0.15f;

	opt->a = 2, opt->b = 4, opt->q = 4, opt->e = 2, opt->q2 = 24, opt->e2 = 1;
	opt->sc_ambi = 1;
	opt->zdrop = 400, opt->zdrop_inv = 200;
	opt->end_bonus = -1;
	opt->min_dp_max = opt->min_chain_score * opt->a;
	opt->min_ksw_len = 200;
	opt->anchor_ext_len = 20, opt->anchor_ext_shift = 6;
	opt->max_clip_ratio = 1.0f;
	opt->mini_batch_size = 500000000;
	opt->max_sw_mat = 100000000;
	opt->cap_kalloc = 1000000000;

	opt->rank_min_len = 500;
	opt->rank_frac = 0.9f;

	opt->pe_ori = 0; // FF
	opt->pe_bonus = 33;
}

int mm_set_opt(const char *preset, mm_idxopt_t *io, mm_mapopt_t *mo)
{
	if (preset == 0) {
		mm_idxopt_init(io);
		mm_mapopt_init(mo);
	} else if (strcmp(preset, "map-ont") == 0) { // this is the same as the default
	} else if (strcmp(preset, "ava-ont") == 0) {
		io->flag = 0, io->k = 19, io->w = 5;
		mo->flag |= MM_F_ALL_CHAINS | MM_F_NO_DIAG | MM_F_NO_DUAL | MM_F_NO_LJOIN;
		mo->min_chain_score = 100, mo->pri_ratio = 0.0f, mo->max_chain_skip = 25;
		mo->bw = mo->bw_long = 2000;
		mo->occ_dist = 0;
	} else if (strcmp(preset, "map10k") == 0 || strcmp(preset, "map-pb") == 0) {
		#if defined (PARALLEL_CHAINING) && (defined(__AVX2__)) && (!defined(__AVX512BW__))
		enable_vect_dp_chaining = false;
		#endif
		io->flag |= MM_I_HPC, io->k = 19;
	} else if (strcmp(preset, "ava-pb") == 0) {
		io->flag |= MM_I_HPC, io->k = 19, io->w = 5;
		mo->flag |= MM_F_ALL_CHAINS | MM_F_NO_DIAG | MM_F_NO_DUAL | MM_F_NO_LJOIN;
		mo->min_chain_score = 100, mo->pri_ratio = 0.0f, mo->max_chain_skip = 25;
		mo->bw_long = mo->bw;
		mo->occ_dist = 0;
	} else if (strcmp(preset, "map-hifi") == 0 || strcmp(preset, "map-ccs") == 0) {
		#if defined (PARALLEL_CHAINING) && (defined(__AVX2__)) && (!defined(__AVX512BW__))
		enable_vect_dp_chaining = false;
		#endif
		io->flag = 0, io->k = 19, io->w = 19;
		mo->max_gap = 10000;
		mo->a = 1, mo->b = 4, mo->q = 6, mo->q2 = 26, mo->e = 2, mo->e2 = 1;
		mo->occ_dist = 500;
		mo->min_mid_occ = 50, mo->max_mid_occ = 500;
		mo->min_dp_max = 200;
	} else if (strncmp(preset, "asm", 3) == 0) {
		io->flag = 0, io->k = 19, io->w = 19;
		mo->bw = 1000, mo->bw_long = 100000;
		mo->max_gap = 10000;
		mo->flag |= MM_F_RMQ;
		mo->min_mid_occ = 50, mo->max_mid_occ = 500;
		mo->min_dp_max = 200;
		mo->best_n = 50;
		if (strcmp(preset, "asm5") == 0) {
			mo->a = 1, mo->b = 19, mo->q = 39, mo->q2 = 81, mo->e = 3, mo->e2 = 1, mo->zdrop = mo->zdrop_inv = 200;
		} else if (strcmp(preset, "asm10") == 0) {
			mo->a = 1, mo->b = 9, mo->q = 16, mo->q2 = 41, mo->e = 2, mo->e2 = 1, mo->zdrop = mo->zdrop_inv = 200;
		} else if (strcmp(preset, "asm20") == 0) {
			mo->a = 1, mo->b = 4, mo->q = 6, mo->q2 = 26, mo->e = 2, mo->e2 = 1, mo->zdrop = mo->zdrop_inv = 200;
			io->w = 10;
		} else return -1;
	} else if (strcmp(preset, "short") == 0 || strcmp(preset, "sr") == 0) {
		io->flag = 0, io->k = 21, io->w = 11;
		mo->flag |= MM_F_SR | MM_F_FRAG_MODE | MM_F_NO_PRINT_2ND | MM_F_2_IO_THREADS | MM_F_HEAP_SORT;
		mo->pe_ori = 0<<1|1; // FR
		mo->a = 2, mo->b = 8, mo->q = 12, mo->e = 2, mo->q2 = 24, mo->e2 = 1;
		mo->zdrop = mo->zdrop_inv = 100;
		mo->end_bonus = 10;
		mo->max_frag_len = 800;
		mo->max_gap = 100;
		mo->bw = mo->bw_long = 100;
		mo->pri_ratio = 0.5f;
		mo->min_cnt = 2;
		mo->min_chain_score = 25;
		mo->min_dp_max = 40;
		mo->best_n = 20;
		mo->mid_occ = 1000;
		mo->max_occ = 5000;
		mo->mini_batch_size = 50000000;
	} else if (strncmp(preset, "splice", 6) == 0 || strcmp(preset, "cdna") == 0) {
		io->flag = 0, io->k = 20, io->w = 5;
		mo->flag |= MM_F_SPLICE | MM_F_SPLICE_FOR | MM_F_SPLICE_REV | MM_F_SPLICE_FLANK;
		mo->max_sw_mat = 0;
		mo->max_gap = 2000, mo->max_gap_ref = mo->bw = mo->bw_long = 200000;
		mo->a = 1, mo->b = 2, mo->q = 2, mo->e = 1, mo->q2 = 32, mo->e2 = 0;
		mo->noncan = 9;
		mo->junc_bonus = 9;
		mo->zdrop = 200, mo->zdrop_inv = 100; // because mo->a is halved
		if (strcmp(preset, "splice:hq") == 0)
			mo->junc_bonus = 5, mo->b = 4, mo->q = 6, mo->q2 = 24;
	} else return -1;
	return 0;
}

int32_t mm_idx_cal_max_occ(const mm_idx_t *mi, float f)
{
	int i;
	size_t n = 0;
	uint32_t thres;
	khint_t *a, k;
	if (f <= 0.) return INT32_MAX;
	for (i = 0; i < 1<<mi->b; ++i)
		if (mi->B[i].h) n += kh_size((idxhash_t*)mi->B[i].h);
	a = (uint32_t*)malloc(n * 4);
	for (i = n = 0; i < 1<<mi->b; ++i) {
		idxhash_t *h = (idxhash_t*)mi->B[i].h;
		if (h == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			if (!kh_exist(h, k)) continue;
			a[n++] = kh_key(h, k)&1? 1 : (uint32_t)kh_val(h, k);
		}
	}
	thres = ks_ksmall_uint32_t(n, a, (uint32_t)((1. - f) * n)) + 1;
	free(a);
	return thres;
}

void mm_mapopt_update(mm_mapopt_t *opt, const mm_idx_t *mi)
{
	if ((opt->flag & MM_F_SPLICE_FOR) || (opt->flag & MM_F_SPLICE_REV))
		opt->flag |= MM_F_SPLICE;
	if (opt->mid_occ <= 0) {
		opt->mid_occ = mm_idx_cal_max_occ(mi, opt->mid_occ_frac);
		if (opt->mid_occ < opt->min_mid_occ)
			opt->mid_occ = opt->min_mid_occ;
		if (opt->max_mid_occ > opt->min_mid_occ && opt->mid_occ > opt->max_mid_occ)
			opt->mid_occ = opt->max_mid_occ;
	}
	if (opt->bw_long < opt->bw) opt->bw_long = opt->bw;
}

// Open index file ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int64_t mm_idx_is_idx(const char *fn)
{
	int fd, is_idx = 0;
	int64_t ret, off_end;
	char magic[4];

	if (strcmp(fn, "-") == 0) return 0; // read from pipe; not an index
	fd = open(fn, O_RDONLY);
	if (fd < 0) return -1; // error
#ifdef WIN32
	if ((off_end = _lseeki64(fd, 0, SEEK_END)) >= 4) {
		_lseeki64(fd, 0, SEEK_SET);
#else
	if ((off_end = lseek(fd, 0, SEEK_END)) >= 4) {
		lseek(fd, 0, SEEK_SET);
#endif // WIN32
		ret = read(fd, magic, 4);
		if (ret == 4 && strncmp(magic, MM_IDX_MAGIC, 4) == 0)
			is_idx = 1;
	}
	close(fd);
	return is_idx? off_end : 0;
}

mm_bseq_file_t *mm_bseq_open(const char *fn)
{
	mm_bseq_file_t *fp;
	gzFile f;
	f = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(0, "r");
	if (f == 0) return 0;
	fp = (mm_bseq_file_t*)calloc(1, sizeof(mm_bseq_file_t));
	fp->fp = f;
	fp->ks = kseq_init(fp->fp);
	return fp;
}

mm_idx_reader_t *mm_idx_reader_open(const char *fn, const mm_idxopt_t *opt, const char *fn_out)
{
	int64_t is_idx;
	mm_idx_reader_t *r;
	// fn -> file name of the index file or ref_seq file 
	is_idx = mm_idx_is_idx(fn);	// is_idx = 0 if fn is a reference sequence file and not an index file
	//fprintf(stderr, "is_idx: %ld\n", is_idx);
	if (is_idx < 0) return 0; // failed to open the index
	// r: reader object
	r = (mm_idx_reader_t*)calloc(1, sizeof(mm_idx_reader_t));
	r->is_idx = is_idx;
	if (opt) r->opt = *opt;
	else mm_idxopt_init(&r->opt);
	if (r->is_idx) {	// open the index file
		r->fp.idx = fopen(fn, "rb");
		r->idx_size = is_idx;
	} else r->fp.seq = mm_bseq_open(fn);	// open the reference sequence file
	if (fn_out) r->fp_out = fopen(fn_out, "wb");
	return r;
}

// Index reader ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void mm_idx_dump(FILE *fp, const mm_idx_t *mi)
{
	uint64_t sum_len = 0;
	uint32_t x[5], i;

	x[0] = mi->w, x[1] = mi->k, x[2] = mi->b, x[3] = mi->n_seq, x[4] = mi->flag;
	fwrite(MM_IDX_MAGIC, 1, 4, fp);
	fwrite(x, 4, 5, fp);
	for (i = 0; i < mi->n_seq; ++i) {
		if (mi->seq[i].name) {
			uint8_t l = strlen(mi->seq[i].name);
			fwrite(&l, 1, 1, fp);
			fwrite(mi->seq[i].name, 1, l, fp);
		} else {
			uint8_t l = 0;
			fwrite(&l, 1, 1, fp);
		}
		fwrite(&mi->seq[i].len, 4, 1, fp);
		sum_len += mi->seq[i].len;
	}
	for (i = 0; i < 1<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		khint_t k;
		idxhash_t *h = (idxhash_t*)b->h;
		uint32_t size = h? h->size : 0;
		fwrite(&b->n, 4, 1, fp);
		fwrite(b->p, 8, b->n, fp);
		fwrite(&size, 4, 1, fp);
		if (size == 0) continue;
		for (k = 0; k < kh_end(h); ++k) {
			uint64_t x[2];
			if (!kh_exist(h, k)) continue;
			x[0] = kh_key(h, k), x[1] = kh_val(h, k);
			fwrite(x, 8, 2, fp);
		}
	}
	if (!(mi->flag & MM_I_NO_SEQ))
		fwrite(mi->S, 4, (sum_len + 7) / 8, fp);
	fflush(fp);
}

mm_idx_t *mm_idx_init(int w, int k, int b, int flag)
{
	mm_idx_t *mi;
	if (k*2 < b) b = k * 2;
	if (w < 1) w = 1;
	mi = (mm_idx_t*)calloc(1, sizeof(mm_idx_t));
	mi->w = w, mi->k = k, mi->b = b, mi->flag = flag;
	mi->B = (mm_idx_bucket_t*)calloc(1<<b, sizeof(mm_idx_bucket_t));
	if (!(mm_dbg_flag & 1)) mi->km = km_init();
	return mi;
}

mm_idx_t *mm_idx_load(FILE *fp)
{
	char magic[4];
	uint32_t x[5], i;
	uint64_t sum_len = 0;
	mm_idx_t *mi;

	if (fread(magic, 1, 4, fp) != 4) return 0;
	if (strncmp(magic, MM_IDX_MAGIC, 4) != 0) return 0;
	if (fread(x, 4, 5, fp) != 5) return 0;
	mi = mm_idx_init(x[0], x[1], x[2], x[4]);
	mi->n_seq = x[3];
	mi->seq = (mm_idx_seq_t*)kcalloc(mi->km, mi->n_seq, sizeof(mm_idx_seq_t));
	for (i = 0; i < mi->n_seq; ++i) {
		uint8_t l;
		mm_idx_seq_t *s = &mi->seq[i];
		fread(&l, 1, 1, fp);
		//klocwork fix
		assert(l >=0 && l <= 255);
		if (l) {
			s->name = (char*)kmalloc(mi->km, l + 1);
			fread(s->name, 1, l, fp);
			s->name[l] = 0;
		}
		fread(&s->len, 4, 1, fp);
		//klocwork fix
		validate_uint(s->len)
		s->offset = sum_len;
		s->is_alt = 0;
		sum_len += s->len;
	}
	for (i = 0; i < 1<<mi->b; ++i) {
		mm_idx_bucket_t *b = &mi->B[i];
		uint32_t j, size;
		khint_t k;
		idxhash_t *h;
		fread(&b->n, 4, 1, fp);
		//klocwork fix
		validate_ulong(b->n)	
		b->p = (uint64_t*)malloc(b->n * 8);
		fread(b->p, 8, b->n, fp);
		fread(&size, 4, 1, fp);
		//klocwork fix
		//assert(0 <= size && size <= UINT_MAX);
		validate_uint(size)
		if (size == 0) continue;
		b->h = h = kh_init(idx);
		kh_resize(idx, h, size);
		for (j = 0; j < size; ++j) {
			uint64_t x[2];
			int absent;
			fread(x, 8, 2, fp);
			k = kh_put(idx, h, x[0], &absent);
			assert(absent);
			kh_val(h, k) = x[1];
		}
	}
	if (!(mi->flag & MM_I_NO_SEQ)) {
		mi->S = (uint32_t*)malloc((sum_len + 7) / 8 * 4);
		fread(mi->S, 4, (sum_len + 7) / 8, fp);
	}
	return mi;
}

void *km_init2(void *km_par, size_t min_core_size)
{
	kmem_t *km;
	km = (kmem_t*)kcalloc(km_par, 1, sizeof(kmem_t));
	km->par = km_par;
	km->min_core_size = min_core_size > 0? min_core_size : 0x80000;
	return (void*)km;
}

void *km_init(void) { return km_init2(0, 0); }

void *kcalloc(void *_km, size_t count, size_t size)
{
	kmem_t *km = (kmem_t*)_km;
	void *p;
	if (size == 0 || count == 0) return 0;
	if (km == NULL) return calloc(count, size);
	p = kmalloc(km, count * size);
	memset(p, 0, count * size);
	return p;
}

static void panic(const char *s)
{
	fprintf(stderr, "%s\n", s);
	abort();
}

static header_t *morecore(kmem_t *km, size_t nu)
{
	header_t *q;
	size_t bytes, *p;
	nu = (nu + 1 + (km->min_core_size - 1)) / km->min_core_size * km->min_core_size; /* the first +1 for core header */
	bytes = nu * sizeof(header_t);
	q = (header_t*)kmalloc(km->par, bytes);
	if (!q) panic("[morecore] insufficient memory");
	q->ptr = km->core_head, q->size = nu, km->core_head = q;
	p = (size_t*)(q + 1);
	*p = nu - 1; /* the size of the free block; -1 because the first unit is used for the core header */
	kfree(km, p + 1); /* initialize the new "core"; NB: the core header is not looped. */
	return km->loop_head;
}

void *kmalloc(void *_km, size_t n_bytes)
{
	kmem_t *km = (kmem_t*)_km;
	size_t n_units;
	header_t *p, *q;

	if (n_bytes == 0) return 0;
	if (km == NULL) return malloc(n_bytes);
	n_units = (n_bytes + sizeof(size_t) + sizeof(header_t) - 1) / sizeof(header_t); /* header+n_bytes requires at least this number of units */

	if (!(q = km->loop_head)) /* the first time when kmalloc() is called, intialize it */
		km->base.ptr = &km->base;
		km->loop_head = km->base.ptr;
		q = km->loop_head; 
	for (p = q->ptr;; q = p, p = p->ptr) { /* search for a suitable block */
		if (p->size >= n_units) { /* p->size if the size of current block. This line means the current block is large enough. */
			if (p->size == n_units) q->ptr = p->ptr; /* no need to split the block */
			else { /* split the block. NB: memory is allocated at the end of the block! */
				p->size -= n_units; /* reduce the size of the free block */
				p += p->size; /* p points to the allocated block */
				*(size_t*)p = n_units; /* set the size */
			}
			km->loop_head = q; /* set the end of chain */
			return (size_t*)p + 1;
		}
		if (p == km->loop_head) { /* then ask for more "cores" */
			if ((p = morecore(km, n_units)) == 0) return 0;
		}
	}
}

void *krealloc(void *_km, void *ap, size_t n_bytes) // TODO: this can be made more efficient in principle
{
	kmem_t *km = (kmem_t*)_km;
	size_t cap, *p, *q;

	if (n_bytes == 0) {
		kfree(km, ap); return 0;
	}
	if (km == NULL) return realloc(ap, n_bytes);
	if (ap == NULL) return kmalloc(km, n_bytes);
	p = (size_t*)ap - 1;
	cap = (*p) * sizeof(header_t) - sizeof(size_t);
	if (cap >= n_bytes) return ap; /* TODO: this prevents shrinking */
	q = (size_t*)kmalloc(km, n_bytes);
	memcpy(q, ap, cap);
	kfree(km, ap);
	return q;
}

void kfree(void *_km, void *ap) /* kfree() also adds a new core to the circular list */
{
	header_t *p, *q;
	kmem_t *km = (kmem_t*)_km;
	
	if (!ap) return;
	if (km == NULL) {
		free(ap);
		return;
	}
	p = (header_t*)((size_t*)ap - 1);
	p->size = *((size_t*)ap - 1);
	/* Find the pointer that points to the block to be freed. The following loop can stop on two conditions:
	 *
	 * a) "p>q && p<q->ptr": @------#++++++++#+++++++@-------    @---------------#+++++++@-------
	 *    (can also be in    |      |                |        -> |                       |
	 *     two cores)        q      p           q->ptr           q                  q->ptr
	 *
	 *                       @--------    #+++++++++@--------    @--------    @------------------
	 *                       |            |         |         -> |            |
	 *                       q            p    q->ptr            q       q->ptr
	 *
	 * b) "q>=q->ptr && (p>q || p<q->ptr)":  @-------#+++++   @--------#+++++++     @-------#+++++   @----------------
	 *                                       |                |        |         -> |                |
	 *                                  q->ptr                q        p       q->ptr                q
	 *
	 *                                       #+++++++@-----   #++++++++@-------     @-------------   #++++++++@-------
	 *                                       |       |                 |         -> |                         |
	 *                                       p  q->ptr                 q       q->ptr                         q
	 */
	for (q = km->loop_head; !(p > q && p < q->ptr); q = q->ptr)
		if (q >= q->ptr && (p > q || p < q->ptr)) break;
	if (p + p->size == q->ptr) { /* two adjacent blocks, merge p and q->ptr (the 2nd and 4th cases) */
		p->size += q->ptr->size;
		p->ptr = q->ptr->ptr;
	} else if (p + p->size > q->ptr && q->ptr >= p) {
		panic("[kfree] The end of the allocated block enters a free block.");
	} else p->ptr = q->ptr; /* backup q->ptr */

	if (q + q->size == p) { /* two adjacent blocks, merge q and p (the other two cases) */
		q->size += p->size;
		q->ptr = p->ptr;
		km->loop_head = q;
	} else if (q + q->size > p && p >= q) {
		panic("[kfree] The end of a free block enters the allocated block.");
	} else km->loop_head = p, q->ptr = p; /* in two cores, cannot be merged; create a new block in the list */
}

int mm_bseq_eof(mm_bseq_file_t *fp)
{
	return (ks_eof(fp->ks->f) && fp->s.seq == 0);
}

static void mm_idx_add(mm_idx_t *mi, int n, const mm128_t *a)
{
	// Add the n minimizers (stored in a) to the index mi
	int i, mask = (1<<mi->b) - 1;
	// mask: a bitmask used to select the bucket for each minimizer
	for (i = 0; i < n; ++i) {
		// a[i].x: hash value of the minimizer
		mm128_v *p = &mi->B[a[i].x>>8&mask].a;
		kv_push(mm128_t, 0, *p, a[i]);	// Add the minimizer to the bucket
	}
}

static inline char *kstrdup(const kstring_t *s)
{
	char *t;
	t = (char*)malloc(s->l + 1);
	memcpy(t, s->s, s->l + 1);
	return t;
}

static inline void kseq2bseq(kseq_t *ks, mm_bseq1_t *s, int with_qual, int with_comment)
{
	int i;
	//if (ks->name.l == 0)
	//	fprintf(stderr, "[WARNING]\033[1;31m empty sequence name in the input.\033[0m\n");
	s->name = kstrdup(&ks->name);
	s->seq = kstrdup(&ks->seq);
	for (i = 0; i < (int)ks->seq.l; ++i) // convert U to T
		if (s->seq[i] == 'u' || s->seq[i] == 'U')
			--s->seq[i];
	s->qual = with_qual && ks->qual.l? kstrdup(&ks->qual) : 0;
	s->comment = with_comment && ks->comment.l? kstrdup(&ks->comment) : 0;
	s->l_seq = ks->seq.l;
}

static inline int mm_qname_len(const char *s)
{
	int l;
	l = strlen(s);
	return l >= 3 && s[l-1] >= '0' && s[l-1] <= '9' && s[l-2] == '/'? l - 2 : l;
}

static inline int mm_qname_same(const char *s1, const char *s2)
{
	int l1, l2;
	l1 = mm_qname_len(s1);
	l2 = mm_qname_len(s2);
	return (l1 == l2 && strncmp(s1, s2, l1) == 0);
}

mm_bseq1_t *mm_bseq_read3(mm_bseq_file_t *fp, int64_t chunk_size, int with_qual, int with_comment, int frag_mode, int *n_)
{
	int64_t size = 0;
	int ret;
	kvec_t(mm_bseq1_t) a = {0,0,0};
	kseq_t *ks = fp->ks;
	*n_ = 0;
	if (fp->s.seq) {
		kv_resize(mm_bseq1_t, 0, a, 256);
		kv_push(mm_bseq1_t, 0, a, fp->s);
		size = fp->s.l_seq;
		memset(&fp->s, 0, sizeof(mm_bseq1_t));
	}
	while ((ret = kseq_read(ks)) >= 0) {
		mm_bseq1_t *s;
		assert(ks->seq.l <= INT32_MAX);
		if (a.m == 0) kv_resize(mm_bseq1_t, 0, a, 256);
		kv_pushp(mm_bseq1_t, 0, a, &s);
		kseq2bseq(ks, s, with_qual, with_comment);
		size += s->l_seq;
		if (size >= chunk_size) {
			if (frag_mode && a.a[a.n-1].l_seq < CHECK_PAIR_THRES) {
				while ((ret = kseq_read(ks)) >= 0) {
					kseq2bseq(ks, &fp->s, with_qual, with_comment);
					if (mm_qname_same(fp->s.name, a.a[a.n-1].name)) {
						kv_push(mm_bseq1_t, 0, a, fp->s);
						memset(&fp->s, 0, sizeof(mm_bseq1_t));
					} else break;
				}
			}
			break;
		}
	}
	if (ret < -1) {
		if (a.n) fprintf(stderr, "[WARNING]\033[1;31m failed to parse the FASTA/FASTQ record next to '%s'. Continue anyway.\033[0m\n", a.a[a.n-1].name);
		else fprintf(stderr, "[WARNING]\033[1;31m failed to parse the first FASTA/FASTQ record. Continue anyway.\033[0m\n");
	}
	*n_ = a.n;
	return a.a;
}

mm_bseq1_t *mm_bseq_read2(mm_bseq_file_t *fp, int64_t chunk_size, int with_qual, int frag_mode, int *n_)
{
	return mm_bseq_read3(fp, chunk_size, with_qual, 0, frag_mode, n_);
}

mm_bseq1_t *mm_bseq_read(mm_bseq_file_t *fp, int64_t chunk_size, int with_qual, int *n_)
{
	return mm_bseq_read2(fp, chunk_size, with_qual, 0, n_);
}

static inline void tq_push(tiny_queue_t *q, int x)
{
	q->a[((q->count++) + q->front) & 0x1f] = x;
}

static inline int tq_shift(tiny_queue_t *q)
{
	int x;
	if (q->count == 0) return -1;
	x = q->a[q->front++];
	q->front &= 0x1f;
	--q->count;
	return x;
}

static inline uint64_t hash64(uint64_t key, uint64_t mask)
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

/**
 * Find symmetric (w,k)-minimizers on a DNA sequence
 *
 * @param km     thread-local memory pool; using NULL falls back to malloc()
 * @param str    DNA sequence
 * @param len    length of $str
 * @param w      find a minimizer for every $w consecutive k-mers
 * @param k      k-mer size
 * @param rid    reference ID; will be copied to the output $p array
 * @param is_hpc homopolymer-compressed or not
 * @param p      minimizers
 *               p->a[i].x = kMer<<8 | kmerSpan
 *               p->a[i].y = rid<<32 | lastPos<<1 | strand
 *               where lastPos is the position of the last base of the i-th minimizer,
 *               and strand indicates whether the minimizer comes from the top or the bottom strand.
 *               Callers may want to set "p->n = 0"; otherwise results are appended to p
 */
void mm_sketch(void *km, const char *str, int len, int w, int k, uint32_t rid, int is_hpc, mm128_v *p)
{
	// DEBUG
	// Al principio de mm_sketch:
	/*fprintf(stderr, "[mm_sketch][input] rid=%u len=%d w=%d k=%d is_hpc=%d seq[0:20]=%.20s ... seq[%d:]=%.20s\n",
        rid, len, w, k, is_hpc, str, len > 20 ? len - 20 : 0, len > 20 ? str + len - 20 : "");
	*/
	// DEBUG
	/*if (rid == 0){
		fprintf(stderr, "seq_nt4_table[0..127]:\n");
    	for (int i = 0; i < 128; ++i) {
        	fprintf(stderr, "%d ", seq_nt4_table[i]);
        	if ((i+1) % 32 == 0) fprintf(stderr, "\n");
    	}
    	fprintf(stderr, "\n");
	}
	*/

	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};
	int i, j, l, buf_pos, min_pos, kmer_span = 0;
	//mm128_t buf[256], min = { UINT64_MAX, UINT64_MAX };
	mm128_t buf[256];
	buf[0].x = UINT64_MAX;
	buf[0].y = UINT64_MAX;
	buf[0].seq[0] = '\0';
	mm128_t min;
	min.x = UINT64_MAX;
	min.y = UINT64_MAX;
	min.seq[0] = '\0';
	tiny_queue_t tq;

	assert(len > 0 && (w > 0 && w < 256) && (k > 0 && k <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
	//memset(buf, 0xff, w * 16);
	memset(buf, 0xff, w * sizeof(mm128_t));
	memset(&tq, 0, sizeof(tiny_queue_t));
	kv_resize(mm128_t, km, *p, p->n + len/w);

	for (i = l = buf_pos = min_pos = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)str[i]];
		mm128_t info;
		info.x = UINT64_MAX;
		info.y = UINT64_MAX;
		info.seq[0] = '\0';
		if (c < 4) { // not an ambiguous base
			int z;
			if (is_hpc) {
				int skip_len = 1;
				if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
					for (skip_len = 2; i + skip_len < len; ++skip_len)
						if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
							break;
					i += skip_len - 1; // put $i at the end of the current homopolymer run
				}
				tq_push(&tq, skip_len);
				kmer_span += skip_len;
				if (tq.count > k) kmer_span -= tq_shift(&tq);
			} else kmer_span = l + 1 < k? l + 1 : k;
			kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
			if (kmer[0] == kmer[1]) continue; // skip "symmetric k-mers" as we don't know it strand
			z = kmer[0] < kmer[1]? 0 : 1; // strand
			++l;
			if (l >= k && kmer_span < 256) {
				info.x = hash64(kmer[z], mask) << 8 | kmer_span;
				info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | z;

#ifdef COLLECT_MMINFO
				// Metrics: store the sequence of the minimizer
				strncpy(info.seq, &str[i - k + 1], k);	// initial position of the kmer: i - k + 1 
				info.seq[k] = '\0';
#endif
			}
		} else l = 0, tq.count = tq.front = 0, kmer_span = 0;
		buf[buf_pos] = info; // need to do this here as appropriate buf_pos and buf[buf_pos] are needed below
		if (l == w + k - 1 && min.x != UINT64_MAX) { // special case for the first window - because identical k-mers are not stored yet
			for (j = buf_pos + 1; j < w; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, km, *p, buf[j]);
			for (j = 0; j < buf_pos; ++j)
				if (min.x == buf[j].x && buf[j].y != min.y) kv_push(mm128_t, km, *p, buf[j]);
		}
		if (info.x <= min.x) { // a new minimum; then write the old min
			if (l >= w + k && min.x != UINT64_MAX) kv_push(mm128_t, km, *p, min);
			min = info, min_pos = buf_pos;
		} else if (buf_pos == min_pos) { // old min has moved outside the window
			if (l >= w + k - 1 && min.x != UINT64_MAX) kv_push(mm128_t, km, *p, min);
			for (j = buf_pos + 1, min.x = UINT64_MAX; j < w; ++j) // the two loops are necessary when there are identical k-mers
				if (min.x >= buf[j].x) min = buf[j], min_pos = j; // >= is important s.t. min is always the closest k-mer
			for (j = 0; j <= buf_pos; ++j)
				if (min.x >= buf[j].x) min = buf[j], min_pos = j;
			if (l >= w + k - 1 && min.x != UINT64_MAX) { // write identical k-mers
				for (j = buf_pos + 1; j < w; ++j) // these two loops make sure the output is sorted
					if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, km, *p, buf[j]);
				for (j = 0; j <= buf_pos; ++j)
					if (min.x == buf[j].x && min.y != buf[j].y) kv_push(mm128_t, km, *p, buf[j]);
			}
		}
		if (++buf_pos == w) buf_pos = 0;
	}
	if (min.x != UINT64_MAX)
		kv_push(mm128_t, km, *p, min);

	// DEBUG
    /*size_t n_added = p->n;
    uint64_t checksum = 0;
    for (size_t i = 0; i < n_added; ++i)
        checksum ^= p->a[i].x ^ p->a[i].y;

    // Imprime resumen de la llamada
    fprintf(stderr, "[mm_sketch] rid=%u len=%d w=%d k=%d is_hpc=%d minimizers=%zu checksum=0x%" PRIx64 "\n",
            rid, len, w, k, is_hpc, n_added, checksum);
	*/
}

static void *worker_pipeline(void *shared, int step, void *in)
{
	int i;
    pipeline_t *p = (pipeline_t*)shared;
    if (step == 0) { // step 0: read sequences
        step_t *s;
		if (p->sum_len > p->batch_size) return 0;	// Not to read more sequences than the batch size
        s = (step_t*)calloc(1, sizeof(step_t));
		// Read a mini-batch of n_seq sequences
		s->seq = mm_bseq_read(p->fp, p->mini_batch_size, 0, &s->n_seq); // read a mini-batch
		if (s->seq) {
			// DEBUG
			//fprintf(stderr, "[worker_pipeline][step 0] n_seq=%d\n", s->n_seq);
            //for (i = 0; i < s->n_seq && i < 5; ++i) // muestra solo las primeras 5 secuencias
            //    fprintf(stderr, "  seq[%d]: name=%s len=%d\n", i, s->seq[i].name, s->seq[i].l_seq);
			uint32_t old_m, m;
			assert((uint64_t)p->mi->n_seq + s->n_seq <= UINT32_MAX); // to prevent integer overflow
			// make room for p->mi->seq
			old_m = p->mi->n_seq, m = p->mi->n_seq + s->n_seq;
			kroundup32(m); kroundup32(old_m);
			if (old_m != m)
				p->mi->seq = (mm_idx_seq_t*)krealloc(p->mi->km, p->mi->seq, m * sizeof(mm_idx_seq_t));
			// make room for p->mi->S
			if (!(p->mi->flag & MM_I_NO_SEQ)) {
				uint64_t sum_len, old_max_len, max_len;
				for (i = 0, sum_len = 0; i < s->n_seq; ++i) sum_len += s->seq[i].l_seq;
				old_max_len = (p->sum_len + 7) / 8;
				max_len = (p->sum_len + sum_len + 7) / 8;
				kroundup64(old_max_len); kroundup64(max_len);
				if (old_max_len != max_len) {
					p->mi->S = (uint32_t*)realloc(p->mi->S, max_len * 4);
					//klocwork fix
					assert(p->mi->S != NULL);
					memset(&p->mi->S[old_max_len], 0, 4 * (max_len - old_max_len));
				}
			}
			// populate p->mi->seq
			for (i = 0; i < s->n_seq; ++i) {
				mm_idx_seq_t *seq = &p->mi->seq[p->mi->n_seq];
				uint32_t j;
				if (!(p->mi->flag & MM_I_NO_NAME)) {
					seq->name = (char*)kmalloc(p->mi->km, strlen(s->seq[i].name) + 1);
					strcpy(seq->name, s->seq[i].name);
				} else seq->name = 0;
				seq->len = s->seq[i].l_seq;
				seq->offset = p->sum_len;
				seq->is_alt = 0;
				// copy the sequence
				if (!(p->mi->flag & MM_I_NO_SEQ)) {
					for (j = 0; j < seq->len; ++j) { // TODO: this is not the fastest way, but let's first see if speed matters here
						uint64_t o = p->sum_len + j;
						int c = seq_nt4_table[(uint8_t)s->seq[i].seq[j]];
						  (p->mi->S, o, c);
					}
				}
				// update p->sum_len and p->mi->n_seq
				p->sum_len += seq->len;
				s->seq[i].rid = p->mi->n_seq++;
			}
			return s;
		} else free(s);	// s->seq is NULL (the sequences have not been read), free s
    } else if (step == 1) { // step 1: compute sketch (set of minimizers) for each read
        step_t *s = (step_t*)in;	// contains the sequences to be processed
		// DEBUG
		//fprintf(stderr, "[worker_pipeline][step 1] n_seq=%d\n", s->n_seq);
        //for (i = 0; i < s->n_seq && i < 5; ++i)
        //    fprintf(stderr, "  seq[%d]: name=%s minimizers=%zu\n", i, s->seq[i].name, s->a.n);		
		for (i = 0; i < s->n_seq; ++i) {	// iterate over the sequences
			mm_bseq1_t *t = &s->seq[i];
			if (t->l_seq > 0)
				// Compute the MINIMIZERS sketch for the sequence
				// The sketch is stored in s->a
				mm_sketch(0, t->seq, t->l_seq, p->mi->w, p->mi->k, t->rid, p->mi->flag&MM_I_HPC, &s->a);
			else if (mm_verbose >= 2)
				fprintf(stderr, "[WARNING] the length database sequence '%s' is 0\n", t->name);
			free(t->seq); free(t->name);
		}
		free(s->seq); s->seq = 0;
		return s;
    } else if (step == 2) { // step 2: dispatch sketch to buckets
        step_t *s = (step_t*)in;	// contains the sketches to be processed
		// DEBUG
		//fprintf(stderr, "[worker_pipeline][step 2] dispatching %zu minimizers\n", s->a.n);
		// s->a contains the sketches
		mm_idx_add(p->mi, s->a.n, s->a.a);	// Add the sketches to the index
		kfree(0, s->a.a); free(s);
	}
    return 0;
}

static void worker_post(void *g, long i, int tid, long thread_index)
{
	int n, n_keys;
	size_t j, start_a, start_p;
	idxhash_t *h;
	mm_idx_t *mi = (mm_idx_t*)g;
	mm_idx_bucket_t *b = &mi->B[i];	// b is the i-th bucket
	if (b->a.n == 0) return;	// If the bucket is empty, return (n: number of minimizers in the bucket)

	// sort by minimizer
	// a.a: array of minimizers (minimizer, position)
	radix_sort_128x(b->a.a, b->a.a + b->a.n);	// Sort the minimizers by their hash values

	// count and preallocate
	for (j = 1, n = 1, n_keys = 0, b->n = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x>>8 != b->a.a[j-1].x>>8) {
			++n_keys;	// Count the number of distinct minimizers (unique hash values)
			if (n > 1) b->n += n;	// Count for the total number of positions
			n = 1;
		} else ++n;
	}

	// DEBUG: print the number of minimizers in the bucket
	//fprintf(stderr, "[M::%s] bucket %zu has %d minimizers and %d distinct minimizers\n", __func__, i, b->a.n, n_keys);

	h = kh_init(idx);
	kh_resize(idx, h, n_keys);	// Resize the hash table to the number of distinct minimizers
	b->p = (uint64_t*)calloc(b->n, 8);

	// create the hash table
	for (j = 1, n = 1, start_a = start_p = 0; j <= b->a.n; ++j) {
		if (j == b->a.n || b->a.a[j].x>>8 != b->a.a[j-1].x>>8) {
			khint_t itr;
			int absent;
			mm128_t *p = &b->a.a[j-1];
			// DEBUG: imprime el rango del grupo actual
        	//fprintf(stderr, "DEBUG: j=%zu, start_a=%zu, n=%d, absent=%d, key=%" PRIu64 "\n", j, start_a, n, absent, p->x>>8>>mi->b<<1);
        	//for (size_t dbg = start_a; dbg < j; ++dbg) {
            //	fprintf(stderr, "  minimizer[%zu]: x=%" PRIu64 " y=%" PRIu64 "\n", dbg, b->a.a[dbg].x, b->a.a[dbg].y);
        	//}
        	itr = kh_put(idx, h, p->x>>8>>mi->b<<1, &absent);
        	// Depuración: imprime la condición del assert
        	//if (!(absent && j == start_a + n)) {
            //	fprintf(stderr, "ASSERT FAIL: absent=%d, j=%zu, start_a=%zu, n=%d, j==start_a+n=%d\n", absent, j, start_a, n, (j==start_a+n));
        	//}
        	assert(absent && j == start_a + n);

			if (n == 1) {
				kh_key(h, itr) |= 1;	// If the minimizer appears only once, set the least significant bit of the key
				kh_val(h, itr) = p->y; // Store directly the position in the value (y: mm information containing the position)
			} else {
				int k;
				for (k = 0; k < n; ++k)
					b->p[start_p + k] = b->a.a[start_a + k].y;	// Store the positions of the minimizers in the bucket
				radix_sort_64(&b->p[start_p], &b->p[start_p + n]); // sort by position; needed as in-place radix_sort_128x() is not stable
				kh_val(h, itr) = (uint64_t)start_p<<32 | n;	// Value: start position in the array and number of minimizers
				start_p += n;
			}
			start_a = j, n = 1;
		} else ++n;
	}
	b->h = h;
	assert(b->n == (int32_t)start_p);

	// deallocate and clear b->a
	kfree(0, b->a.a);
	b->a.n = b->a.m = 0, b->a.a = 0;
}

static inline long steal_work(kt_for_t *t)
{
	int i, min_i = -1;
	long k, min = LONG_MAX;
	for (i = 0; i < t->n_threads; ++i)
		if (min > t->w[i].i) min = t->w[i].i, min_i = i;
	k = __sync_fetch_and_add(&t->w[min_i].i, t->n_threads);
	return k >= t->n? -1 : k;
}

static void *ktf_worker(void *data)
{
	ktf_worker_t *w = (ktf_worker_t*)data;
	long i;
	long original_i = w->i;
	for (;;) {
		i = __sync_fetch_and_add(&w->i, w->t->n_threads);
		if (i >= w->t->n) break;
		w->t->func(w->t->data, i, w - w->t->w, original_i);
	}
	while ((i = steal_work(w->t)) >= 0)
		w->t->func(w->t->data, i, w - w->t->w, original_i);
	pthread_exit(0);
}

void kt_for(int n_threads, void (*func)(void*,long,int,long), void *data, long n)
{
	if (n_threads > 1) {
		int i;
		kt_for_t t;
		pthread_t *tid;
		t.func = func, t.data = data, t.n_threads = n_threads, t.n = n;
		t.w = (ktf_worker_t*)calloc(n_threads, sizeof(ktf_worker_t));
		tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));

		for (i = 0; i < n_threads; ++i){
			t.w[i].t = &t, t.w[i].i = i;
			//fprintf(stderr, "[M::%s] num_thread: %d\n", __func__, i);
		}
		for (i = 0; i < n_threads; ++i){
			pthread_create(&tid[i], 0, ktf_worker, &t.w[i]);
		} 
		for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);
		free(tid); free(t.w);
	} else {
		long j;
		for (j = 0; j < n; ++j) func(data, j, 0, 0);
	}
}

static void mm_idx_post(mm_idx_t *mi, int n_threads)
{
	kt_for(n_threads, worker_post, mi, 1<<mi->b);
}

void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps)
{
	ktp_t aux;
	//pthread_t *tid;
	int i;
	//int res;
	if (n_threads < 1) n_threads = 1;
	aux.n_workers = n_threads;
	aux.n_steps = n_steps;
	aux.func = func;
	aux.shared = shared_data;
	aux.index = 0;
	// klocwork fix;
	//res = pthread_mutex_init(&aux.mutex, 0);
	//assert((res == 0));
	//res = pthread_cond_init(&aux.cv, 0);
	//assert((res == 0));

	aux.workers = (ktp_worker_t*)calloc(n_threads, sizeof(ktp_worker_t));
	for (i = 0; i < n_threads; ++i) {
		ktp_worker_t *w = &aux.workers[i];
		w->step = 0; w->pl = &aux; w->data = 0;
		w->index = aux.index++;
	}

	//tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
	//for (i = 0; i < n_threads; ++i) pthread_create(&tid[i], 0, ktp_worker, &aux.workers[i]);
	//for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);
	
	// Single-threaded execution
    for (i = 0; i < n_threads; ++i) {
        ktp_worker_t *w = &aux.workers[i];
        while (w->step < aux.n_steps) {
            // Call the function directly for each step
            w->data = aux.func(aux.shared, w->step, w->step ? w->data : 0);

            // Update the step
            w->step = w->step == aux.n_steps - 1 || w->data ? (w->step + 1) % aux.n_steps : aux.n_steps;
            if (w->step == 0) w->index = aux.index++;
        }
    }
	
	//free(tid); 
	free(aux.workers);
	// klocwork fix;
	//res = pthread_mutex_destroy(&aux.mutex);
	//assert((res == 0));
	//res = pthread_cond_destroy(&aux.cv);
	//assert((res == 0));
}

double cputime()
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double realtime()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

mm_idx_t *mm_idx_gen(mm_bseq_file_t *fp, int w, int k, int b, int flag, int mini_batch_size, int n_threads, uint64_t batch_size)
{
	// fp: ref seq file, w: window size, k: k-mer size, b: bucket bits, flag: flags, mini_batch_size: mini batch size, n_threads: number of threads, batch_size: batch size
	pipeline_t pl;
	if (fp == 0 || mm_bseq_eof(fp)) return 0;
	memset(&pl, 0, sizeof(pipeline_t));	// initialize pl with 0
	
	pl.mini_batch_size = (uint64_t)mini_batch_size < batch_size? mini_batch_size : batch_size;
	pl.batch_size = batch_size;
	pl.fp = fp;
	pl.mi = mm_idx_init(w, k, b, flag);

	// Read the sequence, compute minimizers and add them to the buckets of the index
	kt_pipeline(n_threads < 3? n_threads : 3, worker_pipeline, &pl, 3);
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] collected minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));

	// Sort the minimizers and generate hash tables
	mm_idx_post(pl.mi, n_threads);
	if (mm_verbose >= 3)
		fprintf(stderr, "[M::%s::%.3f*%.2f] sorted minimizers\n", __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0));

	return pl.mi;
}

mm_idx_t *mm_idx_reader_read(mm_idx_reader_t *r, int n_threads)
{
	mm_idx_t *mi;	// index to be returned
	// r->is_idx: 1 if the input is an index file, 0 if the input is a reference sequence file
	if (r->is_idx) {	// <-- prebuilt index case
		fprintf(stderr, "[M::%s] loading prebuilt index...\n", __func__);
		// load prebuilt index
		mi = mm_idx_load(r->fp.idx);
		if (mi && mm_verbose >= 2 && (mi->k != r->opt.k || mi->w != r->opt.w || (mi->flag&MM_I_HPC) != (r->opt.flag&MM_I_HPC)))
			fprintf(stderr, "[WARNING]\033[1;31m Indexing parameters (-k, -w or -H) overridden by parameters used in the prebuilt index.\033[0m\n");
	} else {
	    fprintf(stderr, "[M::%s] building index from sequences...\n", __func__);
		// build index from sequences
		mi = mm_idx_gen(r->fp.seq, r->opt.w, r->opt.k, r->opt.bucket_bits, r->opt.flag, r->opt.mini_batch_size, n_threads, r->opt.batch_size);
	}
	if (mi) {
		// Write the index to a file if requested
		if (r->fp_out) mm_idx_dump(r->fp_out, mi);
		mi->index = r->n_parts++;
	}
	return mi;
}

void mm_bseq_close(mm_bseq_file_t *fp)
{
	kseq_destroy(fp->ks);
	gzclose(fp->fp);
	free(fp);
}

static inline void mm_revcomp_bseq(mm_bseq1_t *s)
{
	int i, t, l = s->l_seq;
	for (i = 0; i < l>>1; ++i) {
		t = s->seq[l - i - 1];
		s->seq[l - i - 1] = seq_comp_table[(uint8_t)s->seq[i]];
		s->seq[i] = seq_comp_table[t];
	}
	if (l&1) s->seq[l>>1] = seq_comp_table[(uint8_t)s->seq[l>>1]];
	if (s->qual)
		for (i = 0; i < l>>1; ++i)
			t = s->qual[l - i - 1], s->qual[l - i - 1] = s->qual[i], s->qual[i] = t;
}

sdust_buf_t *sdust_buf_init(void *km)
{
	sdust_buf_t *buf;
	buf = (sdust_buf_t*)kcalloc(km, 1, sizeof(sdust_buf_t));
	buf->km = km;
	buf->w = kdq_init(int, buf->km);
	kdq_resize(int, buf->w, 8);
	return buf;
}

void sdust_buf_destroy(sdust_buf_t *buf)
{
	if (buf == 0) return;
	kdq_destroy(int, buf->w);
	kfree(buf->km, buf->P.a); kfree(buf->km, buf->res.a); kfree(buf->km, buf);
}

static inline void save_masked_regions(void *km, uint64_v *res, perf_intv_v *P, int start)
{
	int i, saved = 0;
	perf_intv_t *p;
	if (P->n == 0 || P->a[P->n - 1].start >= start) return;
	p = &P->a[P->n - 1];
	if (res->n) {
		int s = res->a[res->n - 1]>>32, f = (uint32_t)res->a[res->n - 1];
		if (p->start <= f) // if overlapping with or adjacent to the previous interval
			saved = 1, res->a[res->n - 1] = (uint64_t)s<<32 | (f > p->finish? f : p->finish);
	}
	if (!saved) kv_push(uint64_t, km, *res, (uint64_t)p->start<<32|p->finish);
	for (i = P->n - 1; i >= 0 && P->a[i].start < start; --i); // remove perfect intervals that have falled out of the window
	P->n = i + 1;
}

static inline void shift_window(int t, kdq_t(int) *w, int T, int W, int *L, int *rw, int *rv, int *cw, int *cv)
{
	int s;
	if ((int)kdq_size(w) >= W - SD_WLEN + 1) { // TODO: is this right for SD_WLEN!=3?
		//klocwork fix
		//s = *kdq_shift(int, w);
		int *s_ptr = kdq_shift(int, w);
		assert(s_ptr != NULL);
		s = *s_ptr;
		// ------
		
		*rw -= --cw[s];
		if (*L > (int)kdq_size(w))
			--*L, *rv -= --cv[s];
	}
	kdq_push(int, w, t);
	++*L;
	*rw += cw[t]++;
	*rv += cv[t]++;
	if (cv[t] * 10 > T<<1) {
		do {
			// klocwork fix -- fence
			//_mm_lfence();
			s = kdq_at(w, kdq_size(w) - *L);
			*rv -= --cv[s];
			--*L;
		} while (s != t);
	}
}

static void find_perfect(void *km, perf_intv_v *P, const kdq_t(int) *w, int T, int start, int L, int rv, const int *cv)
{
	int c[SD_WTOT], r = rv, i, max_r = 0, max_l = 0;
	memcpy(c, cv, SD_WTOT * sizeof(int));
	for (i = (long)kdq_size(w) - L - 1; i >= 0; --i) {
		//klocwork fix -- fence
		//_mm_lfence();
		int j, t = kdq_at(w, i), new_r, new_l;
		r += c[t]++;
		new_r = r, new_l = kdq_size(w) - i - 1;
		if (new_r * 10 > T * new_l) {
			for (j = 0; j < (int)P->n && P->a[j].start >= i + start; ++j) { // find insertion position
				perf_intv_t *p = &P->a[j];
				if (max_r == 0 || p->r * max_l > max_r * p->l)
					max_r = p->r, max_l = p->l;
			}
			if (max_r == 0 || new_r * max_l >= max_r * new_l) { // then insert
				max_r = new_r, max_l = new_l;
				if (P->n == P->m) kv_resize(perf_intv_t, km, *P, P->n + 1);
				memmove(&P->a[j+1], &P->a[j], (P->n - j) * sizeof(perf_intv_t)); // make room
				++P->n;
				P->a[j].start = i + start, P->a[j].finish = kdq_size(w) + (SD_WLEN - 1) + start;
				P->a[j].r = new_r, P->a[j].l = new_l;
			}
		}
	}
}

const uint64_t *sdust_core(const uint8_t *seq, int l_seq, int T, int W, int *n, sdust_buf_t *buf)
{
	int rv = 0, rw = 0, L = 0, cv[SD_WTOT], cw[SD_WTOT];
	int i, start, l; // _start_: start of the current window; _l_: length of a contiguous A/C/G/T (sub)sequence
	unsigned t; // current word

	buf->P.n = buf->res.n = 0;
	buf->w->front = buf->w->count = 0;
	memset(cv, 0, SD_WTOT * sizeof(int));
	memset(cw, 0, SD_WTOT * sizeof(int));
	if (l_seq < 0) l_seq = strlen((const char*)seq);
	for (i = l = t = 0; i <= l_seq; ++i) {
		int b = i < l_seq? seq_nt4_table[seq[i]] : 4;
		if (b < 4) { // an A/C/G/T base
			++l, t = (t<<2 | b) & SD_WMSK;
			if (l >= SD_WLEN) { // we have seen a word
				start = (l - W > 0? l - W : 0) + (i + 1 - l); // set the start of the current window
				save_masked_regions(buf->km, &buf->res, &buf->P, start); // save intervals falling out of the current window?
				shift_window(t, buf->w, T, W, &L, &rw, &rv, cw, cv);
				if (rw * 10 > L * T)
					find_perfect(buf->km, &buf->P, buf->w, T, start, L, rv, cv);
			}
		} else { // N or the end of sequence; N effectively breaks input into pieces of independent sequences
			start = (l - W + 1 > 0? l - W + 1 : 0) + (i + 1 - l);
			while (buf->P.n) save_masked_regions(buf->km, &buf->res, &buf->P, start++); // clear up unsaved perfect intervals
			l = t = 0;
		}
	}
	*n = buf->res.n;
	return buf->res.a;
}

static int mm_dust_minier(void *km, int n, mm128_t *a, int l_seq, const char *seq, int sdust_thres)
{
	int n_dreg, j, k, u = 0;
	const uint64_t *dreg;
	sdust_buf_t *sdb;
	if (sdust_thres <= 0) return n;
	sdb = sdust_buf_init(km);
	dreg = sdust_core((const uint8_t*)seq, l_seq, sdust_thres, 64, &n_dreg, sdb);
	for (j = k = 0; j < n; ++j) { // squeeze out minimizers that significantly overlap with LCRs
		int32_t qpos = (uint32_t)a[j].y>>1, span = a[j].x&0xff;
		int32_t s = qpos - (span - 1), e = s + span;
		while (u < n_dreg && (int32_t)dreg[u] <= s) ++u;
		if (u < n_dreg && (int32_t)(dreg[u]>>32) < e) {
			int v, l = 0;
			for (v = u; v < n_dreg && (int32_t)(dreg[v]>>32) < e; ++v) { // iterate over LCRs overlapping this minimizer
				int ss = s > (int32_t)(dreg[v]>>32)? s : dreg[v]>>32;
				int ee = e < (int32_t)dreg[v]? e : (uint32_t)dreg[v];
				l += ee - ss;
			}
			if (l <= span>>1) a[k++] = a[j]; // keep the minimizer if less than half of it falls in masked region
		} else a[k++] = a[j];
	}
	sdust_buf_destroy(sdb);
	return k; // the new size
}

static void collect_minimizers(void *km, const mm_mapopt_t *opt, const mm_idx_t *mi, int n_segs, const int *qlens, const char **seqs, mm128_v *mv)
{
	int i, n, sum = 0;
	mv->n = 0;
	for (i = n = 0; i < n_segs; ++i) {
		size_t j;
		mm_sketch(km, seqs[i], qlens[i], mi->w, mi->k, i, mi->flag&MM_I_HPC, mv);
		for (j = n; j < mv->n; ++j)
			mv->a[j].y += sum << 1;
		if (opt->sdust_thres > 0) // mask low-complexity minimizers
			mv->n = n + mm_dust_minier(km, mv->n - n, mv->a + n, qlens[i], seqs[i], opt->sdust_thres);
		sum += qlens[i], n = mv->n;
	}
}

void mm_seed_mz_flt(void *km, mm128_v *mv, int32_t q_occ_max, float q_occ_frac)
{
	// -- DEBUG --
	//fprintf(stderr, "[%s] q_occ_max: %d, q_occ_frac: %.3f, mv->n: %zu\n", __func__, q_occ_max, q_occ_frac, mv->n);
	mm128_t *a;
	size_t i, j, st;
	if (mv->n <= q_occ_max || q_occ_frac <= 0.0f || q_occ_max <= 0) return;
	KMALLOC(km, a, mv->n);
	for (i = 0; i < mv->n; ++i)
		a[i].x = mv->a[i].x, a[i].y = i;
	radix_sort_128x(a, a + mv->n);
	for (st = 0, i = 1; i <= mv->n; ++i) {
		if (i == mv->n || a[i].x != a[st].x) {
			int32_t cnt = i - st;
			if (cnt > q_occ_max && cnt > mv->n * q_occ_frac)
				for (j = st; j < i; ++j)
					mv->a[a[j].y].x = 0;
			st = i;
		}
	}
	kfree(km, a);
	for (i = j = 0; i < mv->n; ++i)
		if (mv->a[i].x != 0)
			mv->a[j++] = mv->a[i];
	mv->n = j;
}

void mm_idx_stat(const mm_idx_t *mi)
{
	int n = 0, n1 = 0;
	uint32_t i;
	uint64_t sum = 0, len = 0;
	fprintf(stderr, "[M::%s] kmer size: %d; skip: %d; is_hpc: %d; #seq: %d\n", __func__, mi->k, mi->w, mi->flag&MM_I_HPC, mi->n_seq);
	for (i = 0; i < mi->n_seq; ++i)
		len += mi->seq[i].len;
	for (i = 0; i < 1U<<mi->b; ++i)
		if (mi->B[i].h) n += kh_size((idxhash_t*)mi->B[i].h);
	for (i = 0; i < 1U<<mi->b; ++i) {
		idxhash_t *h = (idxhash_t*)mi->B[i].h;
		khint_t k;
		if (h == 0) continue;
		for (k = 0; k < kh_end(h); ++k)
			if (kh_exist(h, k)) {
				sum += kh_key(h, k)&1? 1 : (uint32_t)kh_val(h, k);
				if (kh_key(h, k)&1) ++n1;
			}
	}
	//klocwork fix
	assert(n != 0);
	assert(sum != 0);
	// ---- 
	fprintf(stderr, "[M::%s::%.3f*%.2f] distinct minimizers: %d (%.2f%% are singletons); average occurrences: %.3lf; average spacing: %.3lf; total length: %ld\n",
			__func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), n, 100.0*n1/n, (double)sum / n, (double)len / sum, (long)len);
}

void mm_idx_dump_hash(const char* f_name, const mm_idx_t *mi)  
{
	std::vector<hash_entry> v_hash;
	fprintf(stderr, "Building sorted key-val map\n");

	uint32_t i;
	uint64_t num_values = 0;
	for (i = 0; i < 1U<<mi->b; ++i) {
		idxhash_t *h = (idxhash_t*)mi->B[i].h;
		khint_t k;
		if (h == 0) continue;
		for (k = 0; k < kh_end(h); ++k){
			if (kh_exist(h, k)) {
				uint64_t key = kh_key(h, k), bucket_id = i;
				key = key>>1;
				
				key = key<<mi->b | bucket_id;
				
				if(kh_key(h, k)&1)
				{
					v_hash.push_back(hash_entry(key, kh_val(h, k), NULL));
				}
				else
				{
					uint32_t n = (uint32_t)kh_val(h, k);
					v_hash.push_back(hash_entry(key, n, &mi->B[i].p[(kh_val(h, k)>>32) + 0]));
				}
			}
		}
	}
	sort(v_hash.begin(), v_hash.end(), key_sort);
	fprintf(stderr, "Storing hash to %s \n", f_name);

	// =====================================
	//  Exportar idx,key
	// =====================================
	/*string csv_name = (string)f_name + "_keys_ordered.csv";
	ofstream csv(csv_name);
	csv << "idx,key\n";

	for (size_t i = 0; i < v_hash.size(); i++) {
		csv << i << "," << v_hash[i].key << "\n";
	}
	csv.close();*/
	// =====================================

	vector<uint64_t> key_list;
	vector<uint64_t> val_list;
	vector<uint64_t> p_list;

	key_list.push_back(v_hash.size());
	int64_t itr_p = 0;
	uint64_t sum_pos = 0;
	string f1_name = (string)f_name + "_pos_bin";
	string f2_name = (string)f_name + "_val_bin";

	ofstream f1(f1_name, ios::out | ios::binary);
	if (!f1.is_open()) fprintf(stderr, "ERROR: Could not open %s\n", f1_name.c_str());
	ofstream f2(f2_name, ios::out | ios::binary);
	if (!f2.is_open()) fprintf(stderr, "ERROR: Could not open %s\n", f2_name.c_str());
	
	for( int i = 0; i < v_hash.size(); i++){
		key_list.push_back(v_hash[i].key);
		if(v_hash[i].p == NULL){	
			val_list.push_back(sum_pos<<32|(uint64_t)1);
			sum_pos+=1;
			p_list.push_back(v_hash[i].n);
			num_values++;
			continue;
		}
		
		val_list.push_back(sum_pos<<32|(uint64_t)v_hash[i].n);
		sum_pos+=v_hash[i].n;	
			
		num_values+=v_hash[i].n;


		for(int j = 0; j < v_hash[i].n; j++){
			p_list.push_back(v_hash[i].p[j]);
		}	
	}
	f1.write((char*)&val_list[0], (val_list.size())*sizeof(uint64_t));		
	f2.write((char*)&p_list[0], (p_list.size())*sizeof(uint64_t));		
	f1.close();
	f2.close();

	string size_file_name = (string) f_name + "_size";
	ofstream size_f(size_file_name);
	size_f<<v_hash.size()<<" "<<num_values;
	size_f.close();

	string prefix = (string)f_name + "_keys";	
	string keys_bin_file_name = prefix + ".uint64";		
	ofstream wf(keys_bin_file_name, ios::out | ios::binary);
	wf.write((char*)&key_list[0], (key_list.size())*sizeof(uint64_t));		
	wf.close();	

	key_list.clear();
	v_hash.clear();
}

void mm_idx_destroy_mm_hash(mm_idx_t *mi)
{
	//fprintf(stderr, "mm_destroy_hash\n");
	uint32_t i;
	if (mi == 0) return;
	if (mi->h) kh_destroy(str, (khash_t(str)*)mi->h);
	if (mi->B) {
		for (i = 0; i < 1U<<mi->b; ++i) {
			free(mi->B[i].p);
			free(mi->B[i].a.a);
			kh_destroy(idx, (idxhash_t*)mi->B[i].h);
		}
	}
}

void mm_idx_destroy_seq(mm_idx_t *mi)
{
	//fprintf(stderr, "mm_destroy_seq\n");

	uint32_t i;
	if (mi == 0) return;
	if (mi->I) {
		for (i = 0; i < mi->n_seq; ++i)
			free(mi->I[i].a);
		free(mi->I);
	}
	if (!mi->km) {
		for (i = 0; i < mi->n_seq; ++i)
			free(mi->seq[i].name);
		free(mi->seq);
	} else km_destroy(mi->km);
	free(mi->B); free(mi->S); free(mi);
}

// Extract minimizers from a query file ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int mm_extract_minimizers_file_simple(const mm_idx_t *idx, const char *fn, const mm_mapopt_t *opt, 
                                      mm128_v *all_minimizers, uint32_t **seq_ids, int **frag_info)
{
    mm_bseq_file_t *fp;
    int total_seqs = 0, total_frags = 0;
    size_t k;
    void *km;   // memory pool 
    std::vector<int> seg_off_vec, n_seg_vec;

    // Initialize the output minimizers array
    all_minimizers->n = 0;      // number of minimizers collected
    all_minimizers->m = 0;      // initial capacity
    all_minimizers->a = NULL;   // pointer to minimizers array 
    *seq_ids = NULL; *frag_info = NULL;

    // Open input file
    fp = mm_bseq_open(fn);
    if (!fp) {
        fprintf(stderr, "Failed to open query file: %s\n", fn);
        return -1;
    }
    
    // Initialize memory pool
    km = km_init();

    // Inicialization options for file read
    int with_qual = (!!(opt->flag & MM_F_OUT_SAM) && !(opt->flag & MM_F_NO_QUAL));
    int with_comment = !!(opt->flag & MM_F_COPY_COMMENT);
    int frag_mode = (1 > 1 || !!(opt->flag & MM_F_FRAG_MODE));  // n_fp=1, so first part is false
    
    int n_seqs = 0, i, j, n_frag = 0, qlen_sum = 0;
    uint32_t n_processed = 0;
    mm_bseq1_t *seqs;

    // Read all sequences -> not by batches
    // Output: seqs - array of sequences, n_seqs - number of sequences read

    do {
        seqs = mm_bseq_read3(fp, opt->mini_batch_size, with_qual, with_comment, frag_mode, &n_seqs);
        if (!seqs || n_seqs == 0) break;

        // Assign sequence IDs
        for (i = 0; i < n_seqs; ++i)
            seqs[i].rid = n_processed++;

        // Fragment segmentation
        int *seg_off = (int*)calloc(n_seqs, sizeof(int));
        int *n_seg = (int*)calloc(n_seqs, sizeof(int));
        n_frag = 0;
        for (i = 1, j = 0; i <= n_seqs; ++i) {
            if (i == n_seqs || !frag_mode || !mm_qname_same(seqs[i-1].name, seqs[i].name)) {
                n_seg[n_frag] = i - j;		 	// Number of reads/seqs in this fragment
                seg_off[n_frag++] = j;			// Offset of this fragment, n_frag is incremented
                j = i;
            }
        }

        // Process fragments
		// -- DEBUG --
		//fprintf(stderr, "n_frag: %d\n", n_frag);
        for (i = 0; i < n_frag; ++i) {
            int off = seg_off[i], pe_ori = opt->pe_ori;
            int qlens[MM_MAX_SEG];
            const char *qseqs[MM_MAX_SEG];
            mm128_v mv = {0,0,0};
			
			// Prepare query sequences with reverse complement
            for (j = 0; j < n_seg[i]; ++j) {
                if (n_seg[i] == 2 && ((j == 0 && (pe_ori>>1&1)) || (j == 1 && (pe_ori&1))))
                    mm_revcomp_bseq(&seqs[off + j]);
                qlens[j] = seqs[off + j].l_seq;
                qseqs[j] = seqs[off + j].seq;
            }

			// Compute the total length of query sequences
            for (j = 0, qlen_sum = 0; j < n_seg[i]; ++j)
                qlen_sum += qlens[j];
			// Validate fragment -> skip to the next fragment if invalid (empty or too long)
            if (qlen_sum == 0 || n_seg[i] <= 0 || n_seg[i] > MM_MAX_SEG) continue;
            if (opt->max_qlen > 0 && qlen_sum > opt->max_qlen) continue;
			
			// Extract minimizers for this fragment (from his query sequences)
			// -- DEBUG --
			//fprintf(stderr, "[Pre collect_minimizers] Processing fragment %d with %d sequences, total length %d\n", i, n_seg[i], qlen_sum);
            collect_minimizers(km, opt, idx, n_seg[i], qlens, qseqs, &mv);
			// -- DEBUG --
			//fprintf(stderr, "[Post collect_minimizers] Extracted %d minimizers for fragment with %d sequences, total length: %d\n", mv.n, n_seg[i], qlen_sum);
            if (opt->q_occ_frac > 0.0f) {	// Filter minimizers based on occurrence
                mm_seed_mz_flt(km, &mv, opt->mid_occ, opt->q_occ_frac);
			}
			// -- DEBUG --
			//fprintf(stderr, "[Post mm_seed_mz_flt] Actually %d minimizers for fragment with %d sequences, total length: %d\n", mv.n, n_seg[i], qlen_sum);

			// Store minimizers with sequence information in the output array
            if (mv.n > 0) {
                if (all_minimizers->n + mv.n > all_minimizers->m) {  // Resize if needed (capacity exceeded)
                    all_minimizers->m = all_minimizers->n + mv.n + 1000;  // Add extra new space
                    all_minimizers->a = (mm128_t*)realloc(all_minimizers->a, all_minimizers->m * sizeof(mm128_t));  // Resize minimizers array mantaining existing data
                    *seq_ids = (uint32_t*)realloc(*seq_ids, all_minimizers->m * sizeof(uint32_t));
                }
				// Copy minimizers and track which fragment they came from
                memcpy(&all_minimizers->a[all_minimizers->n], mv.a, mv.n * sizeof(mm128_t));
				// Assign sequence IDs to the minimizers
                for (k = 0; k < mv.n; ++k)
                    (*seq_ids)[all_minimizers->n + k] = seqs[off].rid;  // Track sequence ID
                all_minimizers->n += mv.n;  // Update total count of minimizers
            }

			// Restore original sequences (undo reverse complement)
            for (j = 0; j < n_seg[i]; ++j) {
                if (n_seg[i] == 2 && ((j == 0 && (pe_ori>>1&1)) || (j == 1 && (pe_ori&1))))
                    mm_revcomp_bseq(&seqs[off + j]);
            }
			// Free temporary minimizers from this sequence
            kfree(km, mv.a);
        }

        // Store fragment info for this batch
        for (i = 0; i < n_frag; ++i) {
            seg_off_vec.push_back(seg_off[i]);  // Fragment offset
            n_seg_vec.push_back(n_seg[i]);      // Fragment size
        }
        total_frags += n_frag;

        // Cleanup batch
        for (i = 0; i < n_seqs; ++i) {
            free(seqs[i].seq); free(seqs[i].name);
            if (seqs[i].qual) free(seqs[i].qual);
            if (seqs[i].comment) free(seqs[i].comment);
        }
        free(seqs); free(seg_off); free(n_seg);

        total_seqs += n_seqs;
    } while (n_seqs > 0);

    km_destroy(km); mm_bseq_close(fp);

    // Store fragment information for later use
    *frag_info = (int*)malloc(seg_off_vec.size() * 2 * sizeof(int));
	if (*frag_info == NULL) {
		fprintf(stderr, "Failed to allocate memory for fragment information.\n");
		return -1;  // Memory allocation failed
	}
    for (size_t i = 0; i < seg_off_vec.size(); ++i) {
        (*frag_info)[i*2] = seg_off_vec[i];
        (*frag_info)[i*2+1] = n_seg_vec[i];
    }

    //fprintf(stderr, "Extracted %d minimizers from %d fragments\n", (int)all_minimizers->n, n_frag);  // Ho llegueix en batches de 500000 sequencies/fragments
	fprintf(stderr, "Extracted %d minimizers from %d sequences\n", (int)all_minimizers->n, total_seqs);

    // Save minimizers keys to a file
    vector<uint64_t> key_list;
    key_list.push_back(all_minimizers->n);  // First element: number of keys/minimizers
    for (k = 0; k < all_minimizers->n; ++k){
        uint64_t hash = all_minimizers->a[k].x >> 8;  // Extract hash from minimizer
        key_list.push_back(hash);  // Store the key
    }
    // Write keys to a file
    string prefix = string(fn) + "_keys";  // Use input filename as prefix
    string keys_bin_file_name = prefix + ".uint64";  // Output file name
    ofstream wf(keys_bin_file_name, ios::out | ios::binary);
    wf.write((char*)&key_list[0], key_list.size() * sizeof(uint64_t));  // Write keys to file
    wf.close();  // Close the file
    key_list.clear();  // Clear the key list
    fprintf(stderr, "Minimizers keys stored successfully.\n");

    return total_seqs;  // Number of sequences processed
}

void km_destroy(void *_km)
{
	kmem_t *km = (kmem_t*)_km;
	void *km_par;
	header_t *p, *q;
	if (km == NULL) return;
	km_par = km->par;
	for (p = km->core_head; p != NULL;) {
		q = p->ptr;
		kfree(km_par, p);
		p = q;
	}
	kfree(km_par, km);
}

// Helper function to free minimizers
void free_minimizers_simple(mm128_v *minimizers)
{
    if (minimizers && minimizers->a) {
        free(minimizers->a);
        minimizers->a = NULL;
        minimizers->n = 0;
        minimizers->m = 0;
    }
}

void mm_idx_reader_close(mm_idx_reader_t *r)
{
	if (r->is_idx) fclose(r->fp.idx);
	else mm_bseq_close(r->fp.seq);
	if (r->fp_out) fclose(r->fp_out);
	free(r);
}

void mm_idx_destroy(mm_idx_t *mi)
{

	uint32_t i;
	if (mi == 0) return;
	if (mi->h) kh_destroy(str, (khash_t(str)*)mi->h);
	if (mi->B) {
		for (i = 0; i < 1U<<mi->b; ++i) {
			free(mi->B[i].p);
			free(mi->B[i].a.a);
			kh_destroy(idx, (idxhash_t*)mi->B[i].h);
		}
	}
	if (mi->I) {
		for (i = 0; i < mi->n_seq; ++i)
			free(mi->I[i].a);
		free(mi->I);
	}
	if (!mi->km) {
		for (i = 0; i < mi->n_seq; ++i)
			free(mi->seq[i].name);
		free(mi->seq);
	} else km_destroy(mi->km);
	free(mi->B); free(mi->S); free(mi);
}

const uint64_t *mm_idx_get(const mm_idx_t *mi, uint64_t minier, int *n, long thread_index, mm_metrics *mm_info)
{
	int mask = (1<<mi->b) - 1;	// mask to retrieve the bucket
	khint_t k;
	mm_idx_bucket_t *b = &mi->B[minier&mask];	// b is the bucket
	idxhash_t *h = (idxhash_t*)b->h;	// h: hash table
	
	*n = 0; // n: number of positions of the minimizer
	if (h == 0) return 0;

#if defined(COLLECT_METRICS) || defined(TOPDOWN)
    // Store the times we access the reference index
	if (!accessed[thread_index] && num_index_accesses[thread_index] == 0) {
		accessed[thread_index] = true;
		//fprintf(stderr, "[M::%s] INDEX ACCESSED (%d), num_accesses: %d\n", __func__, thread_index, num_index_accesses[thread_index]);
	}
	if (accessed[thread_index]){
		num_index_accesses[thread_index]++;
	}
#endif 
	// look up minimizer in the hash table -> return k: index of the minimizer
#ifdef COLLECT_MMINFO
	k = kh_get_wmetrics(idx, h, minier>>mi->b<<1, thread_index, mm_info);
#else
	k = kh_get(idx, h, minier>>mi->b<<1);
#endif

	// if the minimizer is not found
	if (k == kh_end(h)) return 0;

	// if the minimizer is found
	if (kh_key(h, k)&1) { // special casing when there is only one k-mer
		*n = 1;
		return &kh_val(h, k);
	} else {	// when there are multiple k-mers
		*n = (uint32_t)kh_val(h, k);
		return &b->p[kh_val(h, k)>>32];
	}
}

void update_perf_metrics(int group_counters, int thread_index, uint64_t* start, uint64_t* end) {
    switch (group_counters) {
        case 0:
            cycles[thread_index] += end[0] - start[0];
			instructions[thread_index] += end[1] - start[1];
            all_loads[thread_index] += end[2] - start[2];
            all_stores[thread_index] += end[3] - start[3];
            break;
        case 1:
            llc_references[thread_index] += end[0] - start[0];
            llc_misses[thread_index] += end[1] - start[1];
            cache_references[thread_index] += end[2] - start[2];
            cache_misses[thread_index] += end[3] - start[3];
            break;
        case 2:
            branch_instructions[thread_index] += end[0] - start[0];
            branch_misses[thread_index] += end[1] - start[1];
            retired_branch_instructions[thread_index] += end[2] - start[2];
            mispredicted_branch_instructions[thread_index] += end[3] - start[3];
            break;
        case 3:
            stlb_miss_loads[thread_index] += end[0] - start[0];
            stlb_miss_stores[thread_index] += end[1] - start[1];
            break;
        default:
            break;
    }
}

void update_topdown_metrics(int thread_index, uint64_t* results) {
    td_counters[thread_index].cycles += results[0];
    td_counters[thread_index].instructions += results[1];
    td_counters[thread_index].cache_misses += results[2];
    td_counters[thread_index].slots += results[3];
    td_counters[thread_index].uops_retired += results[4];
    td_counters[thread_index].backend_slots += results[5];
    td_counters[thread_index].bad_spec_slots += results[6];
    td_counters[thread_index].memory_slots += results[7];
}

uint64_t mm_collect_matches(const mm_idx_t *mi, const mm128_v *mv)
{
    // mi: reference index, mv: minimizers of the query
	// [metrics] - minimizers information
	mm_metrics *mm_info = (mm_metrics*)malloc(mv->n * sizeof(mm_metrics));

#ifdef COLLECT_METRICS
	// Perf profiling to obtain microarchitectural metrics [metrics]--------------------------------------------------------------------------------------------------------------------------------------------
    uint64_t start[num_counters], end[num_counters];
	if (group_counters != -1){
		perf_read_all(start);
	}
#endif
#ifdef TOPDOWN
	// Perf profiling to obtain the topdown plot [metrics]--------------------------------------------------------------------------------------------------------------------------------------------
	PerfCounters counters = start_profiling();
#endif

#ifdef COLLECT_MMINFO
	// Total number of minimizers [metrics]----------------------------------------------------------------
	total_num_mm[0] += mv->n;	// total number of minimizers [metrics]

	char seq[40];
	int64_t found, num_matches;
	
	// Open the file to write the information of the minimizers
    FILE *mm_file = fopen("minimizers_output.csv", "a");
    if (mm_file == NULL) {
        fprintf(stderr, "Error: The minimizers file could not be opened for writing.\n");
        return 1; 
    }
#ifdef LISA_HASH
	// Headers of the file
	fprintf(mm_file, "%ld distinct minimizers\nSequence,Key,Accesses,Matches,Range adjustments,Found\n", mv->n);
#else
	// Headers of the file
	fprintf(mm_file, "%ld distinct minimizers\nSequence,Key,Accesses,Matches,Collisions,Found\n", mv->n);
#endif
#endif

#ifdef LISA_HASH
//-----------------------------------
	uint64_t** cr_batch = (uint64_t**) malloc((mv->n)*sizeof(uint64_t*));
	int* t_batch = (int*)malloc((mv->n)*sizeof(int));
	uint64_t* minimizers = (uint64_t*) malloc((mv->n)*sizeof(uint64_t));
	int64_t* lisa_pos = (int64_t*) malloc((max(32, (int)mv->n))* sizeof(int64_t));

#ifdef F64
	double* key = (double*) malloc(mv->n * sizeof(double));
#endif
#ifdef UINT64
	uint64_t* key = (uint64_t*) malloc(mv->n * sizeof(uint64_t));
#endif
	if (!key) { 
		fprintf(stderr, "ERROR: Could not allocate memory for key array\n"); exit(1); 
		return 1;
	}

	// iterate over the minimizers of the query (mv)
	for (size_t i = 0; i < mv->n; i++) {
		mm128_t *p = &mv->a[i];
		minimizers[i] = p->x>>8;
		key[i] = minimizers[i];
	}

#ifdef ENABLE_BATCHING
	// minimizers: array of mm keys, mv->n: number of minimizers, 
	// lisa_pos: array to store the positions of mm, 
	// cr_batch: array to store the pointers to the mm in the reference index
	// t_batch: array to store the number of mm in the reference index (hits)

	// ----------------------------------------------------------------------
	// What I have to measure: seeding matches !!!!
	// ----------------------------------------------------------------------

	// això ho paralelitzem dins de la funcio
	lh->mm_idx_get_batched(minimizers, mv->n, lisa_pos, cr_batch, t_batch, 0, mm_info); 
#endif
#endif

	size_t i;
	int32_t k;

	// aquest for el paralelitzem amb la directiva pragma omp parallel sense el for 
	// perque dividirem nosaltres a ma els rangs perque found i num_matches es modifiquen en cada iteracio 
	// i es compartida (caldrà fer-la de manera que cada thread moidifiqui el valor del seu index/posicio) 
	
	// millor no utilitzar reduction i deixar-ho fet a ma com ho tenim en el main perque reduction implica molt overhead que no volem si fem un profiling 

	//iterperthr=mv->n/omp_get_max_threads();
	// pragma omp parallel shared(mv, mi) private(i, k){
	// int tid = omp_get_thread_num();
	// init = (tid*itpthr);
	// finish = (((id+1)itpthr)-1);
	// for (;init<finish;init++)
	// …
	// num_mm_matched[id]
	// found[id]
	// }

	// iterate over the minimizers of the query (mv)
	for (i = k = 0; i < mv->n; ++i) {
		const uint64_t *cr;
		mm128_t *p = &mv->a[i];
		int t;
		// ACCES THE REFERENCE INDEX
#ifdef LISA_HASH
#ifdef ENABLE_BATCHING
		// uses the batched lookup results
		t = t_batch[i];		
		cr = cr_batch[i];
#else
		cr = lh->get_hash_value(p->x>>8, &t, 0, &mm_info[i]);
#endif
#else	
		// look up the minimizer in the reference index (mi)
		// cr: pointer to the minimizer in the reference index
		cr = mm_idx_get(mi, p->x>>8, &t, 0, &mm_info[i]);
#endif
#ifdef COLLECT_MMINFO
		// Information for the seed's file
		if (t == 0){
			found = 0;
			num_matches = 0;
		}else{
			found = 1;	// if the minimizer is found in the reference index
			num_matches = t;	// the number of matches of the minimizer
			num_mm_matched[0]++;
		}
        num_exact_matches[0] += t;

		// Sequence of the minimizer
		for(int j = 0; j < strlen(p->seq); j++){
			seq[j] = p->seq[j];
		}

		// Write the information of the seed in the file
#ifdef LISA_HASH
#ifdef F64
		fprintf(mm_file, "%.*s,%0.2f,%lu,%d,%lu,%ld\n", mi->k, seq, key[i], mm_info[i].lookup_accesses, num_matches, mm_info[i].range_adjustments, found);
#endif
#ifdef UINT64
		fprintf(mm_file, "%.*s,%lu,%lu,%ld,%lu,%ld\n", mi->k, seq, key[i], mm_info[i].lookup_accesses, num_matches, mm_info[i].range_adjustments, found);
#endif
#else
#ifdef F64
		fprintf(mm_file, "%.*s,%0.2f,%lu,%d,%lu,%ld\n", mi->k, seq, mm_info[i].key, mm_info[i].lookup_accesses, num_matches, mm_info[i].num_collisions, found);
#endif
#ifdef UINT64
		fprintf(mm_file, "%.*s,%lu,%lu,%ld,%lu,%ld\n", mi->k, seq, mm_info[i].key, mm_info[i].lookup_accesses, num_matches, mm_info[i].num_collisions, found);
#endif
#endif
#endif
	}

#ifdef COLLECT_METRICS
    if (group_counters != -1) {
        perf_read_all(end);
        update_perf_metrics(group_counters, 0, start, end);
    }
#endif

#ifdef TOPDOWN
    uint64_t results[8] = {0};
    stop_profiling("mm_collect_matches", counters, results);
    update_topdown_metrics(0, results);
#endif

#ifdef LISA_HASH
	free(key);
	free(cr_batch);
	free(t_batch);
	free(minimizers);
	free(lisa_pos);
#endif

#ifdef COLLECT_MMINFO
	fclose(mm_file);
#endif

	return num_mm_matched[0];
}

#endif // UTILS_MINIMAP2_H