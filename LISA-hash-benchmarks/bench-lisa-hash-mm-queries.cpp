#include <stdio.h>
#include <errno.h>
#include "utils_minimap2.h"
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <immintrin.h>
#include <filesystem>
namespace fs = std::filesystem;

#ifdef LISA_HASH
lisa_hash<uint64_t, uint64_t> *lh;	// Pointer to a lisa_hash object (hash table with key and value of type uint64_t)
#endif

using namespace std;

uint64_t minimizer_lookup_time, alignment_time, dp_time, rmq_time, rmq_t1, rmq_t2, rmq_t3, rmq_t4;
void *km1;
uint64_t km_size = 500000000; // 500 MB
int km_top;
bool enable_vect_dp_chaining = false;

uint64_t proc_freq; // CPU frequency for profiling data


int main(int argc, char *argv[])
{   
    // Set the CPU frequency
	uint64_t tim = __rdtsc();
	sleep(1);
	proc_freq = __rdtsc() - tim;

    fprintf(stderr, "======================================== STARTING benchmark execution: SEEDING phase of genome analysis =======================================\n");

#ifdef LISA_HASH
    fprintf(stderr, "[CONFIG] Running with LISA_HASH (RMI)\n");
#ifdef VECTORIZE
    fprintf(stderr, "[CONFIG] Vectorization: ENABLED\n");
#if (__AVX512BW__)
	fprintf(stderr, "Available AVX512 mode!!\n");	

#elif (__AVX2__)
	fprintf(stderr, "Available AVX2 mode!!\n");

#elif (__AVX__)
	fprintf(stderr, "Available AVX mode!!\n");

#elif (__SSE4_2__)
	fprintf(stderr, "Available SSE4.2 mode!!\n");

#elif (__SSE4_1__)
	fprintf(stderr, "Available SSE4.1 mode!!\n");
#endif
#else
    fprintf(stderr, "[CONFIG] Vectorization: DISABLED\n");
#endif
#ifdef ENABLE_PREFETCH
    fprintf(stderr, "[CONFIG] Prefetching: ENABLED\n");
#else
    fprintf(stderr, "[CONFIG] Prefetching: DISABLED\n");
#endif
#ifdef ENABLE_BATCHING
    fprintf(stderr, "[CONFIG] Batching: ENABLED\n");
#else
    fprintf(stderr, "[CONFIG] Batching: DISABLED\n");
#endif
#else
    fprintf(stderr, "[CONFIG] Running with classic HASH\n");
#endif

	if (remove("minimizers_output.csv") == 0) {
		fprintf(stderr, "[minimizers_output.csv] File deleted successfully\n");
	}

	const char *opt_str = "t:l:";
	ketopt_t o = KETOPT_INIT;
    int c;
	int n_threads = 1;	// By default, use 1 thread and group of counters -1 (NO counters)

	while ((c = ketopt(&o, argc, argv, 1, opt_str, NULL)) >= 0) {
        if (c == 't') n_threads = atoi(o.arg);
        else if (c == 'l') group_counters = atoi(o.arg);
    }

    // Arguments verification: argv[1] -> index file, argv[2] -> query file
    if(argc - o.ind < 1){
        fprintf(stderr, "Invalid number of arguments: 1. index file [2. query file]\n");
		return 1;
	}

	// Perf profiling
	switch (group_counters)
	{
	case 0:
		set_metrics_values(n_threads, 4);	// Set the number of threads to be used and the number of counters 
		perf_init(num_counters, EV_CYCLES, EV_INSTR, MEM_INST_RETIRED_ALL_LOADS, MEM_INST_RETIRED_ALL_STORES);
		break;
	case 1:
		set_metrics_values(n_threads, 4);	// Set the number of threads to be used and the number of counters 
		perf_init(num_counters, LLC_REFERENCES, LLC_MISSES, CACHE_REFERENCES, CACHE_MISSES);
		break;
	case 2:
		set_metrics_values(n_threads, 4);	// Set the number of threads to be used and the number of counters 
		perf_init(num_counters, BRANCH_INSTRUCTIONS, BRANCH_MISSES, BRANCH_INSTRUCTIONS_RETIRED, MISPREDICTED_BRANCH_RETIRED);
		break;
	case 3:
		set_metrics_values(n_threads, 2);	// Set the number of threads to be used and the number of counters 
		perf_init(num_counters, MEM_INST_RETIRED_STLB_MISS_LOADS, MEM_INST_RETIRED_STLB_MISS_STORES);
		break;
	case -1:	// if group_counters == -1 we don't profile
		set_metrics_values(n_threads, 0);
		break;	
	default:
		fprintf(stderr, "[ERROR] \033[1;31m invalid value for the group_counters variable\033[0m\n");
		return 1;
		break;
	}
    
    uint64_t ini_total_time = __rdtsc();	// Start time measurement for the entire execution of the application

    mm_idx_reader_t *idx_rdr;
    mm_idx_t *mi;		// mi : Pointer to the minimap2 index
    
    // Configure options
    char *fnw = 0;
    mm_mapopt_t opt;
    mm_idxopt_t ipt;
    mm_set_opt(0, &ipt, &opt);  // Default options (0)

    size_t last_slash = ((string)argv[o.ind]).find_last_of("/\\");
    string index_filename = (last_slash == string::npos) ? ((string)argv[o.ind]) : ((string)argv[o.ind]).substr(last_slash + 1);

	// --------------------------------------------------------------------------------------------------------------------------------
    // REFERENCE PROCESSING: index loading or building
    // --------------------------------------------------------------------------------------------------------------------------------

    fprintf(stderr, "[REFERENCE]: %s\n", index_filename.c_str());

    uint64_t ini_time = __rdtsc();	// Start time measurement for the index loading/building process

    string input_file(argv[o.ind]);
    string fnw_str;

    int64_t idx = mm_idx_is_idx(argv[o.ind]);
    if (idx < 0) {
        fprintf(stderr, "[ERROR] failed to open file '%s': %s\n", argv[o.ind], strerror(errno));
        return 1;
    } else if (idx) {
        fprintf(stderr, "[INFO] '%s' is a prebuilt index\n", argv[o.ind]);
        fnw = nullptr; // No need to dump the index
    } else {
        fprintf(stderr, "[INFO] '%s' is not an index, it's a reference sequence\n", argv[o.ind]);

        fs::path input_path(input_file);
        fs::path fnw_path = input_path;
        fnw_path += ".mmi";

        fnw_str = fnw_path.string();
        fnw = (char *)fnw_str.c_str();

        fprintf(stderr, "[INFO] the index will be dumped to '%s'\n", fnw);
    }

    idx_rdr = mm_idx_reader_open(argv[o.ind], &ipt, fnw);
    fprintf(stderr, "[INFO - Index] --------------------------------------------------------------------------------------------------------------------------------\n");
    fprintf(stderr, "Index opened successfully: %s\n", argv[o.ind]);

    if (idx_rdr == 0) {
        fprintf(stderr, "[ERROR] failed to open file '%s': %s\n", argv[o.ind], strerror(errno));
         return 1;
    }

    
    mi = mm_idx_reader_read(idx_rdr, n_threads);
    fprintf(stderr, "Index read successfully: %s\n", argv[o.ind]);
        
    // Provide information about the index loading or building process
    if (mm_verbose >= 3)
        fprintf(stderr, "[M::%s::%.3f*%.2f] loaded/built the index for %d target sequence(s)\n",
                __func__, realtime() - mm_realtime0, cputime() / (realtime() - mm_realtime0), mi->n_seq);
            
    // If the verbose level is set to 3, the index statistics will be displayed
    if (mm_verbose >= 3) mm_idx_stat(mi);

    string preset_arg = "map-ont";  // Default preset argument for minimap2
	if (idx_rdr->is_idx != 0){	// If is an index
		size_t pos = input_file.rfind(".mmi");
		input_file = input_file.substr(0, pos); // Remove the .mmi extension
	}
	preset_arg = input_file + "_" + preset_arg + "_minimizers_key_value_sorted"; // Construct the reference minimizers filename

	// Check if the mapping options have been updated
	if (argc > 2) mm_mapopt_update(&opt, mi);

#ifdef LISA_INDEX	// For the case of the LISA index
    mm_idx_dump_hash(preset_arg.c_str(), mi);	// Processes the index and stores the hash table in a file
#endif

    if (argc > 2){
#ifdef LISA_HASH
        fprintf(stderr, "Using LISA_HASH..\n");
        mm_idx_destroy_mm_hash(mi);		// Destruction of any existing hash
        char* prefix = NULL;
        lh = new lisa_hash<uint64_t, uint64_t>(preset_arg, prefix); // Instantiate the LISA hash by providing the built index file (preset_arg)
        fprintf(stderr, "Loading done.\n");
#endif
        uint64_t idx_time = __rdtsc() - ini_time;	// Calculate the total execution time of the index loading/building process

        // --------------------------------------------------------------------------------------------------------------------------------
        // QUERIES PROCESSING: minimizer extraction from the query sequences
        // --------------------------------------------------------------------------------------------------------------------------------

        ini_time = __rdtsc();	// Start time measurement for the minimizer extraction process

        // Single array for all minimizers
        mm128_v all_minimizers = {0,0,0};

        uint32_t *seq_ids = NULL;  // Sequence IDs (optional)
        int *frag_info = NULL;     // Fragment information (optional)
        
        // Extract minimizers
        fprintf(stderr, "[INFO - Queries] ------------------------------------------------------------------------------------------------------------------------------\n");
        fprintf(stderr, "Extracting minimizers from file: %s\n", argv[o.ind + 1]);

        int n_sequences = mm_extract_minimizers_file_simple(mi, argv[o.ind + 1], &opt, &all_minimizers, &seq_ids, &frag_info);

        uint64_t extract_time = __rdtsc() - ini_time;	// Calculate the total execution time of the minimizer extraction process

        uint64_t seeding_time = 0;
        uint64_t num_mm_matched_result = 0;

        if (n_sequences > 0) {
            fprintf(stderr, "Successfully processed %d sequences\n", n_sequences);

            // --------------------------------------------------------------------------------------------------------------------------------
            // SEEDING: lookup minimizers in the index
            // --------------------------------------------------------------------------------------------------------------------------------

            ini_time = __rdtsc();	// Start time measurement for the seeding process
            fprintf(stderr, "[INFO - Seeding] ------------------------------------------------------------------------------------------------------------------------------\n");
            fprintf(stderr, "Looking up minimizers in the index...\n");
            num_mm_matched_result = mm_collect_matches(mi, &all_minimizers);
            fprintf(stderr, "Found %lu matching minimizers\n", num_mm_matched_result);

            // TODO: parallelize all the process

            seeding_time = __rdtsc() - ini_time;	// Calculate the total execution time of the seeding process

            // Free memory
            free_minimizers_simple(&all_minimizers);
            if (seq_ids) free(seq_ids);      
            if (frag_info) free(frag_info);  
        }else {
            fprintf(stderr, "0 sequences processed or an error occurred.\n");
        }
        
    #ifdef LISA_HASH
        mm_idx_destroy_seq(mi);
    #else
        mm_idx_destroy(mi);
    #endif          
        mm_idx_reader_close(idx_rdr); 

        fprintf(stderr, "Execution completed. Resources released successfully.\n");


        uint64_t total_exec_time_app = __rdtsc() - ini_total_time;	// Calculate the total execution time of the mm2-fast application
        // ---------------------------------------------------------------------------------------------------
        // ---------------------------------------------- METRICS --------------------------------------------
        // ---------------------------------------------------------------------------------------------------
        // Open the file to write the metrics information
        char filename[256]; 
        int suffix = 0;
        int i;

        // Generate the metrics file name
        do {
            if (suffix == 0) {
                sprintf(filename, "mm_metrics_g%d_output.csv", group_counters);
            } else {
                sprintf(filename, "mm_metrics_g%d_output_%d.csv", group_counters, suffix);	// If exists, add an incremental suffix
            }
            suffix++;
        } while (access(filename, F_OK) == 0); // Check if the file exists

        FILE *metrics_file = fopen(filename, "w"); 
        if (metrics_file == NULL) {
            fprintf(stderr, "Error: The metrics file could not be opened for writing.\n");
            return 1;
        }

        // Obtain and display the metrics
        MetricsResults metrics = {};  // Initialize all fields to 0
        compute_metrics_results(&metrics);	// Compute the metrics results
        print_metrics_results(metrics_file, &metrics, max_num_threads, proc_freq, total_exec_time_app, idx_time, extract_time, seeding_time);

	    fprintf(stderr, "[INFO - Execution times] ----------------------------------------------------------------------------------------------------------------------\n");
        fprintf(stderr, "Index loading/building time: %0.6f sec\n", idx_time*1.0/proc_freq);
        fprintf(stderr, "Minimizer extraction time: %0.6f sec\n", extract_time*1.0/proc_freq);
        fprintf(stderr, "Seeding time (lookup minimizers in the index): %0.6f sec\n", seeding_time*1.0/proc_freq);
        fprintf(stderr, "Total execution time application: %0.6f sec\n", total_exec_time_app*1.0/proc_freq);

        fclose(metrics_file);

#ifdef LISA_HASH
	    delete lh;
#endif
    } else {
        fprintf(stderr, "No query file provided. EXITING.\n\n");
    }

    return 0;
}