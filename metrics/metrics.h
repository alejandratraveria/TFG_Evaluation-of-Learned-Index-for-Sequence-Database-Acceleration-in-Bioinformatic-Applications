#ifndef METRICS_H
#define METRICS_H

#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>


// minimizers information for metrics
typedef struct {
	uint64_t key;
    uint64_t lookup_accesses;
    uint64_t range_adjustments;
    uint64_t num_collisions;
} mm_metrics;

// Get the number of threads
extern int max_num_threads;

// Global arrays for metrics
extern bool* accessed;   // Array to store if a thread has accessed the index

extern uint64_t* num_index_accesses; // Array to store the number of accesses to the index per thread
extern uint64_t* num_total_index_accesses; // Array to store the total number of accesses to the index (final accesses + range adjustment accesses) per thread

extern uint64_t* accesses_w0_range_adjustments; // Array to store the number of accesses to the index without range adjustments per thread
extern uint64_t* accesses_w1_range_adjustment; // Array to store the number of accesses to the index with 1 range adjustment per thread
extern uint64_t* accesses_w2_range_adjustments; // Array to store the number of accesses to the index with 2 range adjustments per thread
extern uint64_t* accesses_w3_range_adjustments; // Array to store the number of accesses to the index with 3 range adjustments per thread
extern uint64_t* accesses_w4_range_adjustments; // Array to store the number of accesses to the index with 4 range adjustments per thread
extern uint64_t* accesses_w5_range_adjustments; // Array to store the number of accesses to the index with 5 range adjustments per thread
extern uint64_t* accesses_w6_range_adjustments; // Array to store the number of accesses to the index with 6 range adjustments per thread
extern uint64_t* accesses_w7_range_adjustments; // Array to store the number of accesses to the index with 7 range adjustments per thread
extern uint64_t* accesses_w8_range_adjustments; // Array to store the number of accesses to the index with 8 range adjustments per thread
extern uint64_t* accesses_wmore8_range_adjustments; // Array to store the number of accesses to the index with more than 8 range adjustments per thread

extern uint64_t* num_collisions; // Array to store the avg number of collisions per access to the index
extern uint64_t* accesses_wo_collisions; // Array to store the number of accesses to the index without collisions
extern uint64_t* accesses_w1_collision; // Array to store the number of accesses to the index with 1 collision
extern uint64_t* accesses_w2to3_collisions; // Array to store the number of accesses to the index with 2 collisions
extern uint64_t* accesses_w4to7_collisions; // Array to store the number of accesses to the index with 4 collisions
extern uint64_t* accesses_w8to15_collisions; // Array to store the number of accesses to the index with 8 collisions
extern uint64_t* accesses_w16to31_collisions; // Array to store the number of accesses to the index with 16 collisions
extern uint64_t* accesses_w32andmore_collisions; // Array to store the number of accesses to the index with 32 or more collisions

extern uint64_t* num_exact_matches; // Array to store the number of exact matches per thread
extern uint64_t* num_mm_matched; // Array to store the number of queries matched per thread
extern uint64_t* total_num_mm; // Array to store the total number of minimizers per thread
extern uint64_t* exec_real_time_seeding; // Array to store the execution time of the seeding per thread
extern uint64_t* exec_real_time_lookup; // Array to store the execution time of the seeds lookup per thread

// Perf metrics
extern int num_counters;
extern int group_counters;

extern uint64_t* cycles;
extern uint64_t* instructions;
extern uint64_t* all_loads;
extern uint64_t* all_stores;
extern uint64_t* llc_references;
extern uint64_t* llc_misses;
extern uint64_t* cache_references;
extern uint64_t* cache_misses;
extern uint64_t* branch_instructions;
extern uint64_t* branch_misses;
extern uint64_t* retired_branch_instructions;
extern uint64_t* mispredicted_branch_instructions;
extern uint64_t* stlb_miss_loads;
extern uint64_t* stlb_miss_stores;

// Topdown metrics
typedef struct {
    uint64_t cycles;
    uint64_t instructions;
    uint64_t cache_misses;
    uint64_t slots;
    uint64_t uops_retired;
    uint64_t backend_slots;
    uint64_t bad_spec_slots;
    uint64_t memory_slots;
} TopdownCounters;
extern TopdownCounters* td_counters;

typedef struct {
    double backend_bound;
    double bad_spec;
    double retiring;
    double memory_bound;
    double frontend_bound;
    double core_bound;
} TopdownMetrics;

typedef struct {
    uint64_t used_threads;
    uint64_t index_accesses; // Number of index accesses (finals)
    uint64_t total_index_accesses;  // Total number of index accesses (finals + range adjustments)

    // Range adjustments
    uint64_t total_accesses_w0_range_adjustments;  // Total number of index accesses without range adjustments
    uint64_t total_accesses_w1_range_adjustments;  // Total number of index accesses with 1 range adjustment
    uint64_t total_accesses_w2_range_adjustments;  // Total number of index accesses with 2 range adjustments
    uint64_t total_accesses_w3_range_adjustments;  // Total number of index accesses with 3 range adjustments
    uint64_t total_accesses_w4_range_adjustments;  // Total number of index accesses with 4 range adjustments
    uint64_t total_accesses_w5_range_adjustments;  // Total number of index accesses with 5 range adjustments
    uint64_t total_accesses_w6_range_adjustments;  // Total number of index accesses with 6 range adjustments
    uint64_t total_accesses_w7_range_adjustments;  // Total number of index accesses with 7 range adjustments
    uint64_t total_accesses_w8_range_adjustments;  // Total number of index accesses with 8 range adjustments
    uint64_t total_accesses_wmore8_range_adjustments;  // Total number of index accesses with more than 8 range adjustments

    // Collisions
    uint64_t total_collisions; // Array to store the avg number of collisions per access to the index
    uint64_t total_accesses_wo_collisions; // Array to store the number of accesses to the index without collisions
    uint64_t total_accesses_w1_collision; // Array to store the number of accesses to the index with 1 collision
    uint64_t total_accesses_w2to3_collisions; // Array to store the number of accesses to the index with 2 collisions
    uint64_t total_accesses_w4to7_collisions; // Array to store the number of accesses to the index with 4 collisions
    uint64_t total_accesses_w8to15_collisions; // Array to store the number of accesses to the index with 8 collisions
    uint64_t total_accesses_w16to31_collisions; // Array to store the number of accesses to the index with 16 collisions
    uint64_t total_accesses_w32andmore_collisions; // Array to store the number of accesses to the index with 32 or more collisions

    // Matches and minimizers
    uint64_t total_exact_matches;  // Total number of exact matches
    uint64_t total_mm_matched;  // Total number of minimizers matched
    uint64_t total_mm;  // Total number of minimizers
    // Execution time
    uint64_t total_exec_time;
    uint64_t avg_exec_time_seeding;
    uint64_t total_exec_time_lookup;
    uint64_t avg_exec_time_lookup;
    // Profiling
    uint64_t num_cycles;
    uint64_t num_instr;
    uint64_t num_all_loads;
    uint64_t num_all_stores;
    uint64_t num_llc_references;
    uint64_t num_llc_misses;
    uint64_t num_cache_references;
    uint64_t num_cache_misses;
    uint64_t num_branch_instr;
    uint64_t num_branch_misses;
    uint64_t num_retired_branch_instr;
    uint64_t num_mispredicted_branch_instr;
    uint64_t num_stlb_miss_loads;
    uint64_t num_stlb_miss_stores;
    // Topdown
    TopdownCounters td_counters_acc; // All fields initialized to 0
    TopdownMetrics td_metrics;
} MetricsResults;

void set_metrics_values(int num_threads, int n_counters);

void compute_metrics_results(MetricsResults *res);

void print_metrics_results(FILE *metrics_file, const MetricsResults *res, int max_num_threads, uint64_t proc_freq, double total_exec_time_app, uint64_t idx_time, uint64_t extract_time, uint64_t seeding_time);

#endif // METRICS_H