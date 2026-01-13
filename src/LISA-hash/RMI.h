/*************************************************************************************
MIT License
Copyright (c) 2020 Intel Labs

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Authors: Sanchit Misra <sanchit.misra@intel.com>; Saurabh Kalikar <saurabh.kalikar@intel.com>;
*****************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <stdint.h>
#include <limits.h>
#include <immintrin.h>
#include <fstream>
#include <cmath>
#include <cfloat> 
#include <sys/mman.h>
#include "metrics.h"

#define BATCH_SIZE 32

//typedef double rmi_key_t;
//typedef uint64_t rmi_key_t;

enum query_state{
    GUESS_RMI_ROOT,
    GUESS_RMI_LEAF,
    LAST_MILE
};

template<typename rmi_key_t>
class RMI{
    public:
    RMI(char *prefix);
    ~RMI();

    typedef struct batchMetadata
    {
        int64_t qid;
        query_state state;
        rmi_key_t key;
        int64_t modelIndex;
        int64_t first;
        int64_t m;
        int range_adjustments;
    }BatchMetadata;

    void get_random_keys(int64_t nq, rmi_key_t *key_array, int64_t *orig_pos_array);
    rmi_key_t get_element(int64_t index);
    uint64_t get_guess(rmi_key_t key, int64_t *err);
    uint64_t last_mile_search(rmi_key_t key, int64_t guess, size_t err, mm_metrics *mm_info);
    uint64_t last_mile_search_vectorized(rmi_key_t key, int64_t guess, size_t err, mm_metrics *mm_info);
    int64_t lookup(rmi_key_t key, long thread_index, mm_metrics *mm_info);
    void lookup_batched(rmi_key_t *key_array, int64_t num_queries, int64_t *pos_array, long thread_index, mm_metrics *mm_info);
    void print_stats();

    private:
    // klocwork fix -- valid
    RMI& operator=(const RMI&){ return *this;}
    RMI(const RMI& src){ /* do not create copies */ }
    int64_t n;
    rmi_key_t *sorted_array;

    // --- MODELO 16 SEGMENTOS VECTORIZADO ---
    alignas(64) double L0_SPLITS[16];
    // 16 Slopes
    alignas(64) double L0_SLOPES[16];
    // 16 Intercepts
    alignas(64) double L0_INTERCEPTS[16];

    int64_t L1_SIZE = 0;
    double *L1_PARAMETERS;

    bool load_sorted_array(char *prefix);
    bool load_rmi(char *prefix);
    inline int64_t FCLAMP(double inp, double bound);
    inline int64_t get_guess_root_step(rmi_key_t key);
    inline int64_t get_guess_leaf_step(rmi_key_t key, int64_t modelIndex, int64_t *err);
    inline void last_mile_search_one_step(rmi_key_t key, int64_t &first, int64_t &m);
    inline void last_mile_search_vectorized_step(rmi_key_t key, int64_t &first, int64_t &m);
    int process_query_one_step(BatchMetadata &meta, int64_t &pos, long thread_index);

#ifdef STATS
    double total_log_err = 0;
    double total_log_gap = 0;
    int64_t max_err = 0;
    int64_t max_gap = 0;
    int64_t err_hist[20];
    int64_t nq = 0;
#endif
};

template<typename rmi_key_t>
RMI<rmi_key_t>::RMI(char *prefix)
{
    // Inicializar arrays
    memset(L0_SPLITS, 0, 16 * sizeof(double));
    memset(L0_SLOPES, 0, 16 * sizeof(double));
    memset(L0_INTERCEPTS, 0, 16 * sizeof(double));

    load_sorted_array(prefix);
    load_rmi(prefix);
#ifdef STATS
    memset(err_hist, 0, 20*8);
#endif
}

template<typename rmi_key_t>
RMI<rmi_key_t>::~RMI()
{
    free(sorted_array);
    _mm_free(L1_PARAMETERS);
}

template<typename rmi_key_t>
void RMI<rmi_key_t>::get_random_keys(int64_t nq, rmi_key_t *key_array, int64_t *orig_pos_array)
{
    srand(0);
    for(int64_t i = 0; i < nq; i++)
    {
        double rand_num = rand() * 1.0 /(1L + RAND_MAX);
        int64_t query_id = rand_num * nq;
        key_array[i] = sorted_array[query_id];
        orig_pos_array[i] = query_id;
    }
}

template<typename rmi_key_t>
rmi_key_t RMI<rmi_key_t>::get_element(int64_t index)
{
    if(index < 0 || index >= n)
    {
        printf("index is out of bounds, index = %ld\n", index);
        exit(0);
    }
    return sorted_array[index];
}

template<typename rmi_key_t>
bool RMI<rmi_key_t>::load_sorted_array(char *prefix)
{
    std::string filename = prefix;
#ifdef UINT64
    filename = filename + ".uint64";
#endif
#ifdef F64
    filename = filename + ".f64";
#endif
    std::ifstream infile(filename, std::ios::in | std::ios::binary);
    if (!infile.good())
    {
        std::cout<<filename<<" file not found\n";
        exit(0);
    }
    infile.read((char *)&(this->n), sizeof(uint64_t));
    fprintf(stderr, "n = %ld, ", n);
#ifdef MMAP
    sorted_array = (rmi_key_t*) mmap(NULL, n * sizeof(rmi_key_t), PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB, -1, 0);
#else
    sorted_array = (rmi_key_t*) malloc(n * sizeof(rmi_key_t));
#endif
    if (sorted_array == NULL) return false;
    infile.read((char*)sorted_array, n * sizeof(rmi_key_t));
    if (!infile.good()) return false;
    return true;
}

template<typename rmi_key_t>
bool RMI<rmi_key_t>::load_rmi(char *prefix)
{
    std::string filename = prefix;
    filename = filename + ".rmi_PARAMETERS";
    std::ifstream infile(filename, std::ios::in | std::ios::binary);
    if (!infile.good())
    {
        std::cout<<filename<<" file not found\n";
        exit(0);
    }
    
    // 1. Splits (15)
    for (int i = 0; i < 15; i++) {
        infile.read((char *)&(this->L0_SPLITS[i]), sizeof(double));
    }
    this->L0_SPLITS[15] = DBL_MAX;

    // 2. Slopes (16)
    for (int i = 0; i < 16; i++) {
        infile.read((char *)&(this->L0_SLOPES[i]), sizeof(double));
    }

    // 3. Intercepts (16)
    for (int i = 0; i < 16; i++) {
        infile.read((char *)&(this->L0_INTERCEPTS[i]), sizeof(double));
    }

    // Leer tamaño L1
    infile.read((char *)&(this->L1_SIZE), sizeof(int64_t));
    fprintf(stderr,
            "L0 16-Seg SIMD-Ready loaded. Mid Split=%E\n", L0_SPLITS[7]);

    L1_PARAMETERS = (double*) _mm_malloc(L1_SIZE * 3 * sizeof(double), 64);
    if (L1_PARAMETERS == NULL) return false;
    infile.read((char*)L1_PARAMETERS, L1_SIZE * 3 * sizeof(double));
    if (!infile.good()) return false;

#ifdef COLLECT_MMINFO
    // Open the file to write the ERRORS of the models
    FILE *errs_file = fopen("L1_PARAMETERS_errors.csv", "w");
    if (errs_file == NULL) {
        fprintf(stderr, "Error: The file could not be opened for writing.\n");
        return NULL;
    }
    fprintf(errs_file, "model,err\n");
    for (int64_t i = 0; i < L1_SIZE; ++i) {
        uint64_t err = *((uint64_t*)(&L1_PARAMETERS[i * 3 + 2])); // Leer el error como uint64_t
        fprintf(errs_file, "%ld,%ld\n", i, err);
    }
    fclose(errs_file);
#endif

    return true;
}


template<typename rmi_key_t>
inline int64_t RMI<rmi_key_t>::get_guess_root_step(rmi_key_t key)
{
    double key_d = (double)key;
    int index = 0;

#if defined(__AVX512F__)
    __m512d vkey = _mm512_set1_pd(key_d);

    __m512d splits_lo = _mm512_load_pd(&L0_SPLITS[0]); // 0..7
    __m512d splits_hi = _mm512_load_pd(&L0_SPLITS[8]); // 8..15

    // 3. Comparación paralela: key > split?
    //Máscara de bits (1 si true, 0 si false)
    __mmask8 mask_lo = _mm512_cmp_pd_mask(vkey, splits_lo, _CMP_GT_OQ);
    __mmask8 mask_hi = _mm512_cmp_pd_mask(vkey, splits_hi, _CMP_GT_OQ);

    // 4. Contar los bits en 1. 
    index = _mm_popcnt_u32((uint32_t)mask_lo) + _mm_popcnt_u32((uint32_t)mask_hi);

#else
    for(int i = 0; i < 15; i++) {
        index += (key_d > L0_SPLITS[i]);
    }
#endif

    double m = L0_SLOPES[index];
    double b = L0_INTERCEPTS[index];
    
    double fpred = std::fma(m, key_d, b);

    return FCLAMP(fpred, L1_SIZE - 1.0);
}

template<typename rmi_key_t>
uint64_t RMI<rmi_key_t>::get_guess(rmi_key_t key, int64_t *err)
{
    int64_t modelIndex = get_guess_root_step(key);

    // leaf step
    double fpred = std::fma(L1_PARAMETERS[modelIndex * 3 + 1], key, L1_PARAMETERS[modelIndex * 3]);
    *err = *((uint64_t*) (L1_PARAMETERS + (modelIndex * 3 + 2)));
    uint64_t pos = FCLAMP(fpred, n - 1.0);
    return pos;
}

template<typename rmi_key_t>
uint64_t RMI<rmi_key_t>::last_mile_search_vectorized(rmi_key_t key, int64_t guess, size_t err, mm_metrics *mm_info)
{
    int64_t first = guess - err;
    if(first < 0) first = 0;
    int64_t last = guess + err + 1;
    if(last > n) last = n;
    int64_t m = last - first;

#ifdef COLLECT_MMINFO
    mm_info->range_adjustments=0;
#endif

    while(m > 8)
    {
   #ifdef COLLECT_MMINFO
        mm_info->range_adjustments++;
#endif
        int64_t half = m / 2;
        int64_t middle = first + half;
        if(key >= sorted_array[middle])
        {
            first = middle;
            m -= half;
        }
        else
        {
            m = half;
        }
    }

#ifdef F64
    __m512d key_vec = _mm512_set1_pd(key);
    __m512d element_vec = _mm512_loadu_pd(sorted_array + first);
    __mmask8 mask_lt = _mm512_cmplt_pd_mask(element_vec, key_vec);
#endif
#ifdef UINT64
    __m512i key_vec = _mm512_set1_epi64(key);
    __m512i element_vec = _mm512_loadu_si512(sorted_array + first);
    __mmask8 mask_lt = _mm512_cmplt_epi64_mask(element_vec, key_vec);
#endif

    int32_t numlt = _mm_popcnt_u32(mask_lt & ((1 << m) - 1));
    first = first + numlt;
    m = 0;
    return first;
}

template<typename rmi_key_t>
uint64_t RMI<rmi_key_t>::last_mile_search(rmi_key_t key, int64_t guess, size_t err, mm_metrics *mm_info)
{
    int64_t first = guess - err;
    if(first < 0) first = 0;
    int64_t last = guess + err + 1;
    if(last > n) last = n;
    int64_t m = last - first;

#ifdef COLLECT_MMINFO
    mm_info->range_adjustments=0;
#endif

    while(m > 1)
    {
   #ifdef COLLECT_MMINFO
        mm_info->range_adjustments++;
#endif
        int64_t half = m / 2;
        int64_t middle = first + half;
        if(key >= sorted_array[middle])
        {
            first = middle;
            m -= half;
        }
        else
        {
            m = half;
        }
    }
    return first;
}

template<typename rmi_key_t>
int64_t RMI<rmi_key_t>::lookup(rmi_key_t key, long thread_index, mm_metrics *mm_info)
{
#if defined(COLLECT_METRICS) || defined(TOPDOWN)
    if (!accessed[thread_index])
    {
        accessed[thread_index] = true;
    }
#endif

    int64_t err;
    int64_t guess = get_guess(key, &err);
#ifdef VECTORIZE
    int64_t pos = last_mile_search_vectorized(key, guess, err, mm_info);
#else
    int64_t pos = last_mile_search(key, guess, err, mm_info);
#endif

#ifdef COLLECT_MMINFO
    mm_info->lookup_accesses = mm_info->range_adjustments + 1;
    // index_accesses [metrics]--------------------------------------------------------------------------------
    if(mm_info->range_adjustments == 1){
        accesses_w0_range_adjustments[thread_index]++;
        num_total_index_accesses[thread_index]++;
    }else if (mm_info->range_adjustments == 2){
        accesses_w1_range_adjustment[thread_index]++;
        num_total_index_accesses[thread_index] += 2;   // 1 for the final access and 1 for the range adjustment
    }else if (mm_info->range_adjustments == 3){
        accesses_w2_range_adjustments[thread_index]++;
        num_total_index_accesses[thread_index] += 3;   // 1 for the final access and 2 for the range adjustments
    }else if (mm_info->range_adjustments == 4){
        accesses_w3_range_adjustments[thread_index]++;
        num_total_index_accesses[thread_index] += 4;   // 1 for the final access and 3 for the range adjustments
    }else if (mm_info->range_adjustments == 5){
        accesses_w4_range_adjustments[thread_index]++;
        num_total_index_accesses[thread_index] += 5;   // 1 for the final access and 4 for the range adjustments
    }else if (mm_info->range_adjustments == 6){
        accesses_w6_range_adjustments[thread_index]++;
        num_total_index_accesses[thread_index] += 6;   // 1 for the final access and 6 for the range adjustments
    }else if (mm_info->range_adjustments == 7){
        accesses_w6_range_adjustments[thread_index]++;
        num_total_index_accesses[thread_index] += 7;   // 1 for the final access and 6 for the range adjustments
    }else if (mm_info->range_adjustments == 8){
        accesses_w7_range_adjustments[thread_index]++;
        num_total_index_accesses[thread_index] += 8;   // 1 for the final access and 7 for the range adjustments
    }else if (mm_info->range_adjustments == 9){
        accesses_w8_range_adjustments[thread_index]++;
        num_total_index_accesses[thread_index] += 9;   // 1 for the final access and 8 for the range adjustments
    }else if (mm_info->range_adjustments > 9){
        accesses_wmore8_range_adjustments[thread_index]++;
        num_total_index_accesses[thread_index] += 10; // More than 1 for the final access and 9 for the range adjustments
    }
    num_index_accesses[thread_index]++; // increment the number of accesses to the index (final access: verification of the key)
#endif
    if(key != sorted_array[pos])
    {
        return -1;
    }
#ifdef STATS
    nq++;
    int64_t gap = labs(guess-pos);
    int64_t log2_err = log2(1.0 * err);
    err_hist[log2_err]++;
    if(max_err < err) max_err = err;
    if(max_gap < gap) max_gap = gap;
    total_log_err += log2(1.0 + err);
    total_log_gap += log2(1.0 + gap);
#endif
    return pos;
}

template<typename rmi_key_t>
inline int64_t RMI<rmi_key_t>::get_guess_leaf_step(rmi_key_t key, int64_t modelIndex, int64_t *err)
{
    double fpred = std::fma(L1_PARAMETERS[modelIndex * 3 + 1], key, L1_PARAMETERS[modelIndex * 3]);
    *err = *((uint64_t*) (L1_PARAMETERS + (modelIndex * 3 + 2)));
    int64_t guess = FCLAMP(fpred, n - 1.0);
    return guess;
}

template<typename rmi_key_t>
inline void RMI<rmi_key_t>::last_mile_search_one_step(rmi_key_t key, int64_t &first, int64_t &m)
{
    int64_t half = m >> 1;
    int64_t middle = first + half;
    int64_t cond = (key >= sorted_array[middle]);
    first = middle * cond + first * (1 - cond);
    m = (m - half) * cond + half * (1 - cond);
}

#if VECTORIZE && __AVX512BW__
template<typename rmi_key_t>
inline void RMI<rmi_key_t>::last_mile_search_vectorized_step(rmi_key_t key, int64_t &first, int64_t &m)
{
    #ifdef F64
    __m512d key_vec = _mm512_set1_pd(key);
    __m512d element_vec = _mm512_loadu_pd(sorted_array + first);
    __mmask8 mask_lt = _mm512_cmplt_pd_mask(element_vec, key_vec);
#endif
#ifdef UINT64
    __m512i key_vec = _mm512_set1_epi64(key);
    __m512i element_vec = _mm512_loadu_si512(sorted_array + first);
    __mmask8 mask_lt = _mm512_cmplt_epi64_mask(element_vec, key_vec);
#endif

    int32_t numlt = _mm_popcnt_u32(mask_lt & ((1 << m) - 1));
    first = first + numlt;
    m = 0;
}

template<typename rmi_key_t>
int RMI<rmi_key_t>::process_query_one_step(BatchMetadata &bm, int64_t &pos, long thread_index)
{
    if(bm.state == GUESS_RMI_ROOT)
    {
        bm.modelIndex = get_guess_root_step(bm.key);
        bm.state = GUESS_RMI_LEAF;
#ifdef ENABLE_PREFETCH
        _mm_prefetch((const char *)(&L1_PARAMETERS[bm.modelIndex * 3]), _MM_HINT_T0);
        _mm_prefetch((const char *)(&L1_PARAMETERS[bm.modelIndex * 3 + 2]), _MM_HINT_T0);
#endif
    }
    else if(bm.state == GUESS_RMI_LEAF)
    {
        int64_t err;
        int64_t guess = get_guess_leaf_step(bm.key, bm.modelIndex, &err);
        bm.first = guess - err;
        if(bm.first < 0) bm.first = 0;
        int64_t last = guess + err + 1;
        if(last > n) last = n;
        bm.m = last - bm.first;
        bm.state = LAST_MILE;

        int64_t middle = bm.m >> 1;
        int64_t cond = (bm.m > 8);
        int64_t pf0 = bm.first + cond * middle;
        int64_t pf1 = bm.first + (1 - cond) * 7;
#ifdef ENABLE_PREFETCH
        _mm_prefetch((const char *)(&sorted_array[pf0]), _MM_HINT_T0);
        _mm_prefetch((const char *)(&sorted_array[pf1]), _MM_HINT_T0);
#endif
    }
    else
    {
#if defined(COLLECT_METRICS) || defined(TOPDOWN)
        if (!accessed[thread_index])
        {
            accessed[thread_index] = true;
        }
#endif
        if(bm.m > 8)
        {
            last_mile_search_one_step(bm.key, bm.first, bm.m);
            int64_t middle = bm.m >> 1;
            int64_t cond = (bm.m > 8);
            int64_t pf0 = bm.first + cond * middle;
            int64_t pf1 = bm.first + (1 - cond) * 7;
#ifdef ENABLE_PREFETCH
            _mm_prefetch((const char *)(&sorted_array[pf0]), _MM_HINT_T0);
            _mm_prefetch((const char *)(&sorted_array[pf1]), _MM_HINT_T0);
#endif
#ifdef COLLECT_MMINFO
            bm.range_adjustments++;
#endif
        }
#if 1
        else
        {
            last_mile_search_vectorized_step(bm.key, bm.first, bm.m);
#ifdef COLLECT_MMINFO
            bm.range_adjustments++;
            // index_accesses [metrics]--------------------------------------------------------------------------------
            if(bm.range_adjustments == 1){
                accesses_w0_range_adjustments[thread_index]++;
                num_total_index_accesses[thread_index]++;
            }else if (bm.range_adjustments == 2){
                accesses_w1_range_adjustment[thread_index]++;
                num_total_index_accesses[thread_index] += 2;   // 1 for the final access and 1 for the range adjustment
            }else if (bm.range_adjustments == 3){
                accesses_w2_range_adjustments[thread_index]++;
                num_total_index_accesses[thread_index] += 3;   // 1 for the final access and 2 for the range adjustments
            }else if (bm.range_adjustments == 4){
                accesses_w3_range_adjustments[thread_index]++;
                num_total_index_accesses[thread_index] += 4;   // 1 for the final access and 3 for the range adjustments
            }else if (bm.range_adjustments == 5){
                accesses_w4_range_adjustments[thread_index]++;
                num_total_index_accesses[thread_index] += 5;   // 1 for the final access and 4 for the range adjustments
            }else if (bm.range_adjustments == 6){
                accesses_w6_range_adjustments[thread_index]++;
                num_total_index_accesses[thread_index] += 6;   // 1 for the final access and 6 for the range adjustments
            }else if (bm.range_adjustments == 7){
                accesses_w6_range_adjustments[thread_index]++;
                num_total_index_accesses[thread_index] += 7;   // 1 for the final access and 6 for the range adjustments
            }else if (bm.range_adjustments == 8){
                accesses_w7_range_adjustments[thread_index]++;
                num_total_index_accesses[thread_index] += 8;   // 1 for the final access and 7 for the range adjustments
            }else if (bm.range_adjustments == 9){
                accesses_w8_range_adjustments[thread_index]++;
                num_total_index_accesses[thread_index] += 9;   // 1 for the final access and 8 for the range adjustments
            }else if (bm.range_adjustments > 9){
                accesses_wmore8_range_adjustments[thread_index]++;
                num_total_index_accesses[thread_index] += 10; // More than 1 for the final access and 9 for the range adjustments
            }
            num_index_accesses[thread_index]++; // increment the number of accesses to the index (final access: verification of the key)
#endif
            pos = bm.first;
        if(sorted_array[pos] != bm.key)
            pos = -1;
            return 0;
        }
#else
        if(bm.m == 1)
        {
            pos = bm.first;
            return 0;
        }
#endif
    }
    return 1;
}

#else
template<typename rmi_key_t>
int RMI<rmi_key_t>::process_query_one_step(BatchMetadata &bm, int64_t &pos, long thread_index)
{
    if(bm.state == GUESS_RMI_ROOT)
    {
        bm.modelIndex = get_guess_root_step(bm.key);
        bm.state = GUESS_RMI_LEAF;
#ifdef ENABLE_PREFETCH
        _mm_prefetch((const char *)(&L1_PARAMETERS[bm.modelIndex * 3]), _MM_HINT_T0);
        _mm_prefetch((const char *)(&L1_PARAMETERS[bm.modelIndex * 3 + 2]), _MM_HINT_T0);
#endif
    }
    else if(bm.state == GUESS_RMI_LEAF)
    {
        int64_t err;
        int64_t guess = get_guess_leaf_step(bm.key, bm.modelIndex, &err);
        
        bm.first = guess - err;
        if(bm.first < 0) bm.first = 0;
        int64_t last = guess + err + 1;
        if(last > n) last = n;
        bm.m = last - bm.first;
        bm.state = LAST_MILE;

        int64_t middle = bm.m >> 1;
#ifdef ENABLE_PREFETCH
        _mm_prefetch((const char *)(&sorted_array[bm.first + middle]), _MM_HINT_T0);
#endif
    }
    else
    {
#if defined(COLLECT_METRICS) || defined(TOPDOWN)
        // the thread accessed to the index [metrics]----------------------------------------
        if (!accessed[thread_index])
        {
            accessed[thread_index] = true;
        }
#endif
        if(bm.m > 1)
        {
            last_mile_search_one_step(bm.key, bm.first, bm.m);
            int64_t middle = bm.m >> 1;
#ifdef ENABLE_PREFETCH
            _mm_prefetch((const char *)(&sorted_array[bm.first + middle]), _MM_HINT_T0);
#endif
#ifdef COLLECT_MMINFO
            bm.range_adjustments++;
#endif
        }
        if(bm.m == 1)   // Query is done with a range of 1 position
        {
#ifdef COLLECT_MMINFO
            // index_accesses [metrics]--------------------------------------------------------------------------------
            if(bm.range_adjustments == 1){
                accesses_w0_range_adjustments[thread_index]++;
                num_total_index_accesses[thread_index]++;
            }else if (bm.range_adjustments == 2){
                accesses_w1_range_adjustment[thread_index]++;
                num_total_index_accesses[thread_index] += 2;   // 1 for the final access and 1 for the range adjustment
            }else if (bm.range_adjustments == 3){
                accesses_w2_range_adjustments[thread_index]++;
                num_total_index_accesses[thread_index] += 3;   // 1 for the final access and 2 for the range adjustments
            }else if (bm.range_adjustments == 4){
                accesses_w3_range_adjustments[thread_index]++;
                num_total_index_accesses[thread_index] += 4;   // 1 for the final access and 3 for the range adjustments
            }else if (bm.range_adjustments == 5){
                accesses_w4_range_adjustments[thread_index]++;
                num_total_index_accesses[thread_index] += 5;   // 1 for the final access and 4 for the range adjustments
            }else if (bm.range_adjustments == 6){
                accesses_w6_range_adjustments[thread_index]++;
                num_total_index_accesses[thread_index] += 6;   // 1 for the final access and 6 for the range adjustments
            }else if (bm.range_adjustments == 7){
                accesses_w6_range_adjustments[thread_index]++;
                num_total_index_accesses[thread_index] += 7;   // 1 for the final access and 6 for the range adjustments
            }else if (bm.range_adjustments == 8){
                accesses_w7_range_adjustments[thread_index]++;
                num_total_index_accesses[thread_index] += 8;   // 1 for the final access and 7 for the range adjustments
            }else if (bm.range_adjustments == 9){
                accesses_w8_range_adjustments[thread_index]++;
                num_total_index_accesses[thread_index] += 9;   // 1 for the final access and 8 for the range adjustments
            }else if (bm.range_adjustments > 9){
                accesses_wmore8_range_adjustments[thread_index]++;
                num_total_index_accesses[thread_index] += 10; // More than 1 for the final access and 9 for the range adjustments
            }
            num_index_accesses[thread_index]++; // increment the number of accesses to the index (final access: verification of the key)
#endif
            pos = bm.first;
        if(sorted_array[pos] != bm.key)
            pos = -1;
            return 0;
        }
    }
    return 1;
}
#endif

#if 0
int RMI<rmi_key_t>::process_query_one_step(BatchMetadata &bm, int64_t &pos){
    if(bm.state == GUESS_RMI_ROOT)
    {
        bm.modelIndex = get_guess_root_step(bm.key);
        bm.state = GUESS_RMI_LEAF;
#ifdef ENABLE_PREFETCH
        _mm_prefetch((const char *)(&L1_PARAMETERS[bm.modelIndex * 3]), _MM_HINT_T0);
        _mm_prefetch((const char *)(&L1_PARAMETERS[bm.modelIndex * 3 + 2]), _MM_HINT_T0);
#endif
    }
    else if(bm.state == GUESS_RMI_LEAF)
    {
        int64_t err;
        int64_t guess = get_guess_leaf_step(bm.key, bm.modelIndex, &err);
        bm.first = guess - err;
        if(bm.first < 0) bm.first = 0;
        int64_t last = guess + err + 1;
        if(last > n) last = n;
        bm.m = last - bm.first;
        bm.state = LAST_MILE;

#ifdef ENABLE_PREFETCH
        if(bm.m > 8)
        {
            _mm_prefetch((const char *)(&sorted_array[bm.first + (bm.m >> 1)]), _MM_HINT_T0);
        }
        else
        {
            _mm_prefetch((const char *)(&sorted_array[bm.first]), _MM_HINT_T0);
            _mm_prefetch((const char *)(&sorted_array[bm.first + 7]), _MM_HINT_T0);
        }
#endif
    }
    else
    {
        if(bm.m > 8)
        {
            last_mile_search_one_step(bm.key, bm.first, bm.m);
#ifdef ENABLE_PREFETCH
            if(bm.m > 8)
            {
                _mm_prefetch((const char *)(&sorted_array[bm.first + (bm.m >> 1)]), _MM_HINT_T0);
            }
            else
            {
                _mm_prefetch((const char *)(&sorted_array[bm.first]), _MM_HINT_T0);
                _mm_prefetch((const char *)(&sorted_array[bm.first + 7]), _MM_HINT_T0);
            }
#endif
        }
#if 1
        else
        {
            last_mile_search_vectorized_step(bm.key, bm.first, bm.m);
            pos = bm.first;
            return 0;
        }
#else
        if(bm.m == 1)
        {
            pos = bm.first;
            return 0;
        }
#endif
    }
    return 1;
}
#endif

template<typename rmi_key_t>
void RMI<rmi_key_t>::lookup_batched(rmi_key_t *key_array, int64_t num_queries, int64_t *pos_array, long thread_index, mm_metrics *mm_info)
{
    int64_t next_query;
    BatchMetadata metadata[BATCH_SIZE];

    for(next_query = 0; next_query < BATCH_SIZE && next_query < num_queries; next_query++)
    {
        BatchMetadata bm;
        memset(&bm, 0, 1*sizeof(BatchMetadata));
        bm.qid = next_query;
        bm.state = GUESS_RMI_ROOT;
        bm.key = key_array[next_query];
#ifdef COLLECT_MMINFO
        bm.range_adjustments = 0;
#endif
        metadata[next_query] = bm;
    }

    int64_t num_active = next_query;

    while(num_active > 0)
    {
        for(int64_t i = 0; i < num_active; i++)
        {
            int64_t pos;
            int status = process_query_one_step(metadata[i], pos, thread_index);
            if(!status) // Query is done (status=0)
            {
                BatchMetadata bm = metadata[i];
                pos_array[bm.qid] = pos;

#ifdef COLLECT_MMINFO
                if (mm_info != nullptr){
                    mm_info[bm.qid].range_adjustments = bm.range_adjustments;
                    mm_info[bm.qid].lookup_accesses = bm.range_adjustments + 1;
                }
#endif
                if(next_query < num_queries)    // there are more queries to process
                {
                    BatchMetadata bm;
                    memset(&bm, 0, 1*sizeof(BatchMetadata));
                    bm.qid = next_query;
                    bm.state = GUESS_RMI_ROOT;
                    bm.key = key_array[next_query];
                    metadata[i] = bm;
                    next_query++;
                }
                else
                {
                    metadata[i] = metadata[num_active - 1];
                    num_active--;
                }
            }
        }
    }
}

template<typename rmi_key_t>
inline int64_t RMI<rmi_key_t>::FCLAMP(double inp, double bound) {
  if (inp < 0.0) return 0;
  return (inp > bound ? bound : (size_t)inp);
}

template<typename rmi_key_t>
void RMI<rmi_key_t>::print_stats()
{
#ifdef STATS
    printf("avg_log2_err = %lf, avg_log2_gap = %lf\n", total_log_err / nq, total_log_gap / nq);
    printf("max_err = %ld, max_gap = %ld\n", max_err, max_gap);
    for(int i = 0; i < 20; i++)
    {
        printf("%ld] log2_err freq = %ld\n", i, err_hist[i]);
    }
#endif
}