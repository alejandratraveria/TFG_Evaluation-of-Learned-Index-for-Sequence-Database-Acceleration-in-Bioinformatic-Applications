#ifndef PERF_PROFILER_H
#define PERF_PROFILER_H

#include <linux/perf_event.h>
#include <stdint.h>
#include <sys/ioctl.h>
#include <unistd.h>
#include <stdio.h>
//#include <sys/syscall.h>

// Struct to hold file descriptors for counters
typedef struct PerfCounters {
    int cycles_fd;
    int instr_fd;
    int cache_fd;
    int topdown_slots_fd; 
    int topdown_backend_slots_fd;
    int topdown_bad_spec_slots_fd; 
    int topdown_memory_slots_fd; 
    int uops_retired_fd; 
} PerfCounters;

// Function to start profiling a region
PerfCounters start_profiling();

// Function to stop profiling and print results
void stop_profiling(const char* name, PerfCounters counters, uint64_t results[]);

#endif // PERF_PROFILER_H
