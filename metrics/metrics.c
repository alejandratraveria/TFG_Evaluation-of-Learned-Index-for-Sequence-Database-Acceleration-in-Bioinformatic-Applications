// metrics.c
#include "metrics.h"
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>


// Get the number of threads
int max_num_threads;

bool* accessed;

uint64_t* num_index_accesses; 
uint64_t* num_total_index_accesses; 

uint64_t* accesses_w0_range_adjustments; 
uint64_t* accesses_w1_range_adjustment; 
uint64_t* accesses_w2_range_adjustments; 
uint64_t* accesses_w3_range_adjustments;
uint64_t* accesses_w4_range_adjustments; 
uint64_t* accesses_w5_range_adjustments;
uint64_t* accesses_w6_range_adjustments;
uint64_t* accesses_w7_range_adjustments;
uint64_t* accesses_w8_range_adjustments; 
uint64_t* accesses_wmore8_range_adjustments;

uint64_t* num_collisions;
uint64_t* accesses_wo_collisions; 
uint64_t* accesses_w1_collision; 
uint64_t* accesses_w2to3_collisions; 
uint64_t* accesses_w4to7_collisions;
uint64_t* accesses_w8to15_collisions;
uint64_t* accesses_w16to31_collisions;
uint64_t* accesses_w32andmore_collisions;

uint64_t* num_exact_matches; 
uint64_t* num_mm_matched; 
uint64_t* total_num_mm; 

// Perf metrics
int num_counters; // The number of counters in the group
int group_counters = -1; // The group of counters that we want to analyze
uint64_t* cycles;
uint64_t* instructions;
uint64_t* all_loads;
uint64_t* all_stores;
uint64_t* llc_references;
uint64_t* llc_misses;
uint64_t* cache_references;
uint64_t* cache_misses;
uint64_t* branch_instructions;
uint64_t* branch_misses;
uint64_t* retired_branch_instructions;
uint64_t* mispredicted_branch_instructions;
uint64_t* stlb_miss_loads;
uint64_t* stlb_miss_stores;

// Topdown metrics
TopdownCounters* td_counters;

void set_metrics_values(int num_threads, int n_counters){
    max_num_threads = num_threads;
    num_counters = n_counters;
    accessed = (bool*)malloc(max_num_threads * sizeof(bool));
    num_index_accesses = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    num_total_index_accesses = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));

    accesses_w0_range_adjustments = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    accesses_w1_range_adjustment = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    accesses_w2_range_adjustments = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    accesses_w3_range_adjustments = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    accesses_w4_range_adjustments = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    accesses_w5_range_adjustments = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    accesses_w6_range_adjustments = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    accesses_w7_range_adjustments = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    accesses_w8_range_adjustments = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    accesses_wmore8_range_adjustments = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));

    num_collisions = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    accesses_wo_collisions = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    accesses_w1_collision = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    accesses_w2to3_collisions = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    accesses_w4to7_collisions = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    accesses_w8to15_collisions = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    accesses_w16to31_collisions = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    accesses_w32andmore_collisions = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));

    num_exact_matches = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    num_mm_matched = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    total_num_mm = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    cycles = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    instructions = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    all_loads = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    all_stores = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    llc_references = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    llc_misses = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    cache_references = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    cache_misses = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    branch_instructions = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    branch_misses = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    retired_branch_instructions = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    mispredicted_branch_instructions = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    stlb_miss_loads = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));
    stlb_miss_stores = (uint64_t*)malloc(max_num_threads * sizeof(uint64_t));

    // Initialize the topdown metrics arrays
    td_counters = (TopdownCounters*)malloc(max_num_threads * sizeof(TopdownCounters));

    // Initialize all arrays to zero
    for (int i = 0; i < max_num_threads; i++) {
        num_index_accesses[i] = 0;
        num_total_index_accesses[i] = 0;

        accesses_w0_range_adjustments[i] = 0;
        accesses_w1_range_adjustment[i] = 0;
        accesses_w2_range_adjustments[i] = 0;
        accesses_w3_range_adjustments[i] = 0;
        accesses_w4_range_adjustments[i] = 0;
        accesses_w5_range_adjustments[i] = 0;
        accesses_w6_range_adjustments[i] = 0;
        accesses_w7_range_adjustments[i] = 0;
        accesses_w8_range_adjustments[i] = 0;
        accesses_wmore8_range_adjustments[i] = 0;

        num_collisions[i] = 0;
        accesses_wo_collisions[i] = 0;
        accesses_w1_collision[i] = 0;
        accesses_w2to3_collisions[i] = 0;
        accesses_w4to7_collisions[i] = 0;
        accesses_w8to15_collisions[i] = 0;
        accesses_w16to31_collisions[i] = 0;
        accesses_w32andmore_collisions[i] = 0;

        num_exact_matches[i] = 0;
        num_mm_matched[i] = 0;
        total_num_mm[i] = 0;
        accessed[i] = false;
        cycles[i] = 0;
        instructions[i] = 0;
        all_loads[i] = 0;
        all_stores[i] = 0;
        llc_references[i] = 0;
        llc_misses[i] = 0;
        cache_references[i] = 0;
        cache_misses[i] = 0;
        branch_instructions[i] = 0;
        branch_misses[i] = 0;
        retired_branch_instructions[i] = 0;
        mispredicted_branch_instructions[i] = 0;
        stlb_miss_loads[i] = 0;
        stlb_miss_stores[i] = 0;

        td_counters[i].cycles = 0;
        td_counters[i].instructions = 0;
        td_counters[i].cache_misses = 0;
        td_counters[i].slots = 0;
        td_counters[i].uops_retired = 0;
        td_counters[i].backend_slots = 0;
        td_counters[i].bad_spec_slots = 0;
        td_counters[i].memory_slots = 0;
    }
}

void compute_metrics_results(MetricsResults *res) {
    int i;
    for (i = 0; i < max_num_threads; i++) {
        if (accessed[i]) {
            res->used_threads++;
            res->index_accesses += num_index_accesses[i];
            res->total_index_accesses += num_total_index_accesses[i];

            res->total_accesses_w0_range_adjustments += accesses_w0_range_adjustments[i];
            res->total_accesses_w1_range_adjustments += accesses_w1_range_adjustment[i];
            res->total_accesses_w2_range_adjustments += accesses_w2_range_adjustments[i];
            res->total_accesses_w3_range_adjustments += accesses_w3_range_adjustments[i];
            res->total_accesses_w4_range_adjustments += accesses_w4_range_adjustments[i];
            res->total_accesses_w5_range_adjustments += accesses_w5_range_adjustments[i];
            res->total_accesses_w6_range_adjustments += accesses_w6_range_adjustments[i];
            res->total_accesses_w7_range_adjustments += accesses_w7_range_adjustments[i];
            res->total_accesses_w8_range_adjustments += accesses_w8_range_adjustments[i];
            res->total_accesses_wmore8_range_adjustments += accesses_wmore8_range_adjustments[i];

            res->total_collisions += num_collisions[i];
            res->total_accesses_wo_collisions += accesses_wo_collisions[i];
            res->total_accesses_w1_collision += accesses_w1_collision[i];
            res->total_accesses_w2to3_collisions += accesses_w2to3_collisions[i];
            res->total_accesses_w4to7_collisions += accesses_w4to7_collisions[i];
            res->total_accesses_w8to15_collisions += accesses_w8to15_collisions[i];
            res->total_accesses_w16to31_collisions += accesses_w16to31_collisions[i];
            res->total_accesses_w32andmore_collisions += accesses_w32andmore_collisions[i];

            res->total_exact_matches += num_exact_matches[i];
            res->total_mm_matched += num_mm_matched[i];

            res->num_cycles += cycles[i];
            res->num_instr += instructions[i];
            res->num_all_loads += all_loads[i];
            res->num_all_stores += all_stores[i];
            res->num_llc_references += llc_references[i];
            res->num_llc_misses += llc_misses[i];
            res->num_cache_references += cache_references[i];
            res->num_cache_misses += cache_misses[i];
            res->num_branch_instr += branch_instructions[i];
            res->num_branch_misses += branch_misses[i];
            res->num_retired_branch_instr += retired_branch_instructions[i];
            res->num_mispredicted_branch_instr += mispredicted_branch_instructions[i];
            res->num_stlb_miss_loads += stlb_miss_loads[i];
            res->num_stlb_miss_stores += stlb_miss_stores[i];

            res->td_counters_acc.cycles += td_counters[i].cycles;
            res->td_counters_acc.instructions += td_counters[i].instructions;
            res->td_counters_acc.cache_misses += td_counters[i].cache_misses;
            res->td_counters_acc.slots += td_counters[i].slots;
            res->td_counters_acc.uops_retired += td_counters[i].uops_retired;
            res->td_counters_acc.backend_slots += td_counters[i].backend_slots;
            res->td_counters_acc.bad_spec_slots += td_counters[i].bad_spec_slots;
            res->td_counters_acc.memory_slots += td_counters[i].memory_slots;
        }
        res->total_mm += total_num_mm[i];
    }
#ifndef LISA_HASH
    res->total_index_accesses = res->index_accesses + res->total_collisions; // Total accesses = final accesses + accesses due to collisions
#endif

    // Topdown metrics
    res->td_metrics.backend_bound = ((double)res->td_counters_acc.backend_slots) / ((double)res->td_counters_acc.slots) * 100.0;
	res->td_metrics.bad_spec = ((double)res->td_counters_acc.bad_spec_slots) / ((double)res->td_counters_acc.slots) * 100.0;
	res->td_metrics.retiring = ((double)res->td_counters_acc.uops_retired) / ((double)res->td_counters_acc.slots) * 100.0;
	res->td_metrics.memory_bound = ((double)res->td_counters_acc.memory_slots) / ((double)res->td_counters_acc.slots) * 100.0;
	res->td_metrics.frontend_bound = 100.00 - (res->td_metrics.backend_bound + res->td_metrics.bad_spec + res->td_metrics.retiring);
	res->td_metrics.core_bound = res->td_metrics.backend_bound - res->td_metrics.memory_bound;
}

void print_metrics_results(FILE *metrics_file, const MetricsResults *res, int max_num_threads, uint64_t proc_freq, double total_exec_time_app, uint64_t idx_time, uint64_t extract_time, uint64_t seeding_time){
	// Display the metrics
	fprintf(metrics_file, "[M::%s] Metrics results:\n", __func__);
	fprintf(metrics_file, "Num_threads,%u\n", max_num_threads);
	fprintf(metrics_file, "Used threads,%lu\n\n", res->used_threads);
	fprintf(metrics_file, "Number of final accesses to the index,%lu\n", res->index_accesses);
	fprintf(metrics_file, "Total number of accesses to the index,%lu\n\n", res->total_index_accesses);
#ifdef LISA_HASH
	fprintf(metrics_file, "Total number of accesses without range adjustments,%lu\n", res->total_accesses_w0_range_adjustments);
	fprintf(metrics_file, "Total number of accesses with 1 range adjustment,%lu\n", res->total_accesses_w1_range_adjustments);
	fprintf(metrics_file, "Total number of accesses with 2 range adjustments,%lu\n", res->total_accesses_w2_range_adjustments);
	fprintf(metrics_file, "Total number of accesses with 3 range adjustments,%lu\n", res->total_accesses_w3_range_adjustments);
	fprintf(metrics_file, "Total number of accesses with 4 range adjustments,%lu\n", res->total_accesses_w4_range_adjustments);
	fprintf(metrics_file, "Total number of accesses with 5 range adjustments,%lu\n", res->total_accesses_w5_range_adjustments);
	fprintf(metrics_file, "Total number of accesses with 6 range adjustments,%lu\n", res->total_accesses_w6_range_adjustments);
	fprintf(metrics_file, "Total number of accesses with 7 range adjustments,%lu\n", res->total_accesses_w7_range_adjustments);
	fprintf(metrics_file, "Total number of accesses with 8 range adjustments,%lu\n", res->total_accesses_w8_range_adjustments);
	fprintf(metrics_file, "Total number of accesses with more than 8 range adjustments,%lu\n\n", res->total_accesses_wmore8_range_adjustments);
#else
    fprintf(metrics_file, "Total number of collisions,%llu\n", res->total_collisions);
	fprintf(metrics_file, "Total number of accesses without collisions,%llu\n", res->total_accesses_wo_collisions);
	fprintf(metrics_file, "Total number of accesses with 1 collision,%llu\n", res->total_accesses_w1_collision);
	fprintf(metrics_file, "Total number of accesses with 2 to 3 collisions,%llu\n", res->total_accesses_w2to3_collisions);
	fprintf(metrics_file, "Total number of accesses with 4 to 7 collisions,%llu\n", res->total_accesses_w4to7_collisions);
	fprintf(metrics_file, "Total number of accesses with 8 to 15 collisions,%llu\n", res->total_accesses_w8to15_collisions);
	fprintf(metrics_file, "Total number of accesses with 16 to 31 collisions,%llu\n", res->total_accesses_w16to31_collisions);
	fprintf(metrics_file, "Total number of accesses with 32 or more collisions,%llu\n\n", res->total_accesses_w32andmore_collisions);
#endif
    fprintf(metrics_file, "Total number of exact matches,%lu\n", res->total_exact_matches);
	fprintf(metrics_file, "Total number of minimizers,%lu\n", res->total_mm);
	fprintf(metrics_file, "Total number of minimizers matched,%lu\n\n", res->total_mm_matched);

	fprintf(metrics_file, "[M::%s] Profiling results:\n", __func__);
	fprintf(metrics_file, "Number of cycles seeds lookup,%lu\n", res->num_cycles);
	fprintf(metrics_file, "Number of instructions seeds lookup,%lu\n", res->num_instr);
	fprintf(metrics_file, "Loads in seeds lookup,%lu\n", res->num_all_loads);
	fprintf(metrics_file, "Stores in seeds lookup,%lu\n", res->num_all_stores);
	fprintf(metrics_file, "LLC references in seeds lookup,%lu\n", res->num_llc_references);
	fprintf(metrics_file, "LLC misses in seeds lookup,%lu\n", res->num_llc_misses);
	fprintf(metrics_file, "Cache references in seeds lookup,%lu\n", res->num_cache_references);
	fprintf(metrics_file, "Cache misses in seeds lookup,%lu\n", res->num_cache_misses);
	fprintf(metrics_file, "Branch instructions in seeds lookup,%lu\n", res->num_branch_instr);
	fprintf(metrics_file, "Branch misses in seeds lookup,%lu\n", res->num_branch_misses);
	fprintf(metrics_file, "Retired branch instructions,%lu\n", res->num_retired_branch_instr);
	fprintf(metrics_file, "Mispredicted branch instructions,%lu\n", res->num_mispredicted_branch_instr);
	fprintf(metrics_file, "STLB load misses in seeds lookup,%lu\n", res->num_stlb_miss_loads);
	fprintf(metrics_file, "STLB store misses in seeds lookup,%lu\n\n", res->num_stlb_miss_stores);

	fprintf(metrics_file, "IPC,%0.2f\n", res->num_instr*1.0/res->num_cycles);
	fprintf(metrics_file, "LLC miss rate,%0.2f\n", (res->num_llc_misses)*1.0/res->num_llc_references);
	fprintf(metrics_file, "Cache miss rate,%0.2f\n", (res->num_cache_misses)*1.0/res->num_cache_references);
	fprintf(metrics_file, "Misprediction rate,%0.2f\n", (res->num_branch_misses)*1.0/res->num_branch_instr);
	fprintf(metrics_file, "Retired misprediction rate,%0.2f\n\n", (res->num_mispredicted_branch_instr)*1.0/res->num_retired_branch_instr);

	fprintf(metrics_file, "[M::%s] Topdown counters:\n", __func__);
	fprintf(metrics_file, "Cycles,%lu\n", res->td_counters_acc.cycles);
	fprintf(metrics_file, "Instructions,%lu\n", res->td_counters_acc.instructions);
	fprintf(metrics_file, "Cache misses,%lu\n", res->td_counters_acc.cache_misses);
	fprintf(metrics_file, "Slots,%lu\n", res->td_counters_acc.slots);
	fprintf(metrics_file, "Uops retired,%lu\n", res->td_counters_acc.uops_retired);
	fprintf(metrics_file, "Backend slots,%lu\n", res->td_counters_acc.backend_slots);
	fprintf(metrics_file, "Bad spec slots,%lu\n", res->td_counters_acc.bad_spec_slots);
	fprintf(metrics_file, "Memory slots,%lu\n\n", res->td_counters_acc.memory_slots);
	fprintf(metrics_file, "Topdown metrics:\n");
	fprintf(metrics_file, "Useful work,%0.3f\n", res->td_metrics.retiring);
	fprintf(metrics_file, "Bad speculation,%0.3f\n", res->td_metrics.bad_spec);
	fprintf(metrics_file, "Frontend bound,%0.3f\n", res->td_metrics.frontend_bound);
	fprintf(metrics_file, "Backend bound,%0.3f\n", res->td_metrics.backend_bound);
	fprintf(metrics_file, "Memory bound,%0.3f\n", res->td_metrics.memory_bound);
	fprintf(metrics_file, "Core bound,%0.3f\n", res->td_metrics.core_bound);
	fprintf(metrics_file, "\n\n");

    fprintf(metrics_file, "[M::%s] Execution times:\n", __func__);
    fprintf(metrics_file, "Index loading/building time,%0.6f\n", idx_time*1.0/proc_freq);
    fprintf(metrics_file, "Minimizer extraction time,%0.6f\n", extract_time*1.0/proc_freq);
    fprintf(metrics_file, "Seeding time (lookup minimizers in the index),%0.6f\n", seeding_time*1.0/proc_freq);
    fprintf(metrics_file, "Total execution time application,%0.6f\n", total_exec_time_app*1.0/proc_freq);
}