#include "perf_profiler.h"
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <linux/perf_event.h>

// Direct syscall for opening performance counters
static long perf_event_open(struct perf_event_attr *hw_event, pid_t pid,
                            int cpu, int group_fd, unsigned long flags) {
    return syscall(__NR_perf_event_open, hw_event, pid, cpu, group_fd, flags);
}

// Setup a performance counter for a specific metric
static int setup_counter(uint32_t type, uint64_t config) {
    struct perf_event_attr pe;
    memset(&pe, 0, sizeof(struct perf_event_attr));
    pe.size = sizeof(struct perf_event_attr);
    pe.type = type;
    pe.config = config;
    pe.disabled = 1;  // Start disabled
    pe.exclude_kernel = 1;
    pe.exclude_hv = 1;
    
    int fd = perf_event_open(&pe, 0, -1, -1, 0);
    if (fd == -1) {
        perror("perf_event_open");
    }
    return fd;
}

static int setup_raw_counter(uint8_t event, uint8_t umask) {
    struct perf_event_attr pe;
    memset(&pe, 0, sizeof(struct perf_event_attr));
    pe.size = sizeof(struct perf_event_attr);
    pe.type = PERF_TYPE_RAW;
    pe.config = (umask << 8) | event;
    pe.disabled = 1;  // Start disabled
    pe.exclude_kernel = 1;
    pe.exclude_hv = 1;
    
    int fd = perf_event_open(&pe, 0, -1, -1, 0);
    if (fd == -1) {
        perror("perf_event_open");
    }
    return fd;

}

// Start profiling: setup and enable counters
PerfCounters start_profiling() {
    PerfCounters counters;
    
    // Setup counters for CPU cycles, instructions, and cache misses
    counters.cycles_fd = setup_counter(PERF_TYPE_HARDWARE, PERF_COUNT_HW_CPU_CYCLES);
    counters.instr_fd = setup_counter(PERF_TYPE_HARDWARE, PERF_COUNT_HW_INSTRUCTIONS);
    counters.cache_fd = setup_counter(PERF_TYPE_HARDWARE, PERF_COUNT_HW_CACHE_MISSES);

    // Setup RAW counters for microarchitecture analysis
    counters.topdown_slots_fd           = setup_raw_counter(0xA4, 0x01);
    counters.uops_retired_fd            = setup_raw_counter(0xC2, 0x02); 
    counters.topdown_backend_slots_fd   = setup_raw_counter(0xA4, 0x02);
    counters.topdown_bad_spec_slots_fd  = setup_raw_counter(0xA4, 0x04);
    counters.topdown_memory_slots_fd    = setup_raw_counter(0xA4, 0x10);

    // Reset and enable the counters
    ioctl(counters.cycles_fd,                   PERF_EVENT_IOC_RESET, 0);
    ioctl(counters.instr_fd,                    PERF_EVENT_IOC_RESET, 0);
    ioctl(counters.cache_fd,                    PERF_EVENT_IOC_RESET, 0);
    ioctl(counters.topdown_slots_fd,            PERF_EVENT_IOC_RESET, 0);
    ioctl(counters.uops_retired_fd,             PERF_EVENT_IOC_RESET, 0);
    ioctl(counters.topdown_backend_slots_fd,    PERF_EVENT_IOC_RESET, 0);
    ioctl(counters.topdown_bad_spec_slots_fd,   PERF_EVENT_IOC_RESET, 0);
    ioctl(counters.topdown_memory_slots_fd,     PERF_EVENT_IOC_RESET, 0);

    ioctl(counters.cycles_fd,                   PERF_EVENT_IOC_ENABLE, 0);
    ioctl(counters.instr_fd,                    PERF_EVENT_IOC_ENABLE, 0);
    ioctl(counters.cache_fd,                    PERF_EVENT_IOC_ENABLE, 0);
    ioctl(counters.topdown_slots_fd,            PERF_EVENT_IOC_ENABLE, 0);
    ioctl(counters.uops_retired_fd,             PERF_EVENT_IOC_ENABLE, 0);
    ioctl(counters.topdown_backend_slots_fd,    PERF_EVENT_IOC_ENABLE, 0);
    ioctl(counters.topdown_bad_spec_slots_fd,   PERF_EVENT_IOC_ENABLE, 0);
    ioctl(counters.topdown_memory_slots_fd,     PERF_EVENT_IOC_ENABLE, 0);
    return counters;
}

// Stop profiling: disable counters, read values, print results, and close FDs
void stop_profiling(const char* name, PerfCounters counters, uint64_t results[]) {
    // Disable counters
    ioctl(counters.cycles_fd,                   PERF_EVENT_IOC_DISABLE, 0);
    ioctl(counters.instr_fd,                    PERF_EVENT_IOC_DISABLE, 0);
    ioctl(counters.cache_fd,                    PERF_EVENT_IOC_DISABLE, 0);
    ioctl(counters.topdown_slots_fd,            PERF_EVENT_IOC_DISABLE, 0);
    ioctl(counters.uops_retired_fd,             PERF_EVENT_IOC_DISABLE, 0);
    ioctl(counters.topdown_backend_slots_fd,    PERF_EVENT_IOC_DISABLE, 0);
    ioctl(counters.topdown_bad_spec_slots_fd,   PERF_EVENT_IOC_DISABLE, 0);
    ioctl(counters.topdown_memory_slots_fd,     PERF_EVENT_IOC_DISABLE, 0);
    // Read counter values
    long long cycles, instructions, cache_misses, slots, uops_retired, backend_slots, bad_spec_slots, memory_slots;
    read(counters.cycles_fd,                   &cycles,         sizeof(long long));
    read(counters.instr_fd,                    &instructions,   sizeof(long long));
    read(counters.cache_fd,                    &cache_misses,   sizeof(long long));
    read(counters.topdown_slots_fd,            &slots,          sizeof(long long));
    read(counters.uops_retired_fd,             &uops_retired,   sizeof(long long));
    read(counters.topdown_backend_slots_fd,    &backend_slots,  sizeof(long long));
    read(counters.topdown_bad_spec_slots_fd,   &bad_spec_slots, sizeof(long long));
    read(counters.topdown_memory_slots_fd,     &memory_slots,   sizeof(long long));
    
    double backend_bound = ((double)backend_slots)  /   ((double)slots) * 100.0; 
    double bad_spec      = ((double)bad_spec_slots) /   ((double)slots) * 100.0; 
    double retiring      = ((double)uops_retired)   /   ((double)slots) * 100.0; 
    double mem_bound     = ((double)memory_slots)   /   ((double)slots) * 100.0; 
    double frontend_bound= 100.00 - (backend_bound + bad_spec + retiring);
    double core_bound    = backend_bound - mem_bound;
    
    // Store results in the array
    results[0] = cycles;
    results[1] = instructions;
    results[2] = cache_misses;
    results[3] = slots;
    results[4] = uops_retired;
    results[5] = backend_slots;
    results[6] = bad_spec_slots;
    results[7] = memory_slots;
    
    // Print results
    //printf("Region %s:\n",name);
    //printf("  Cycles: %lld\n", cycles);
    //printf("  Instructions: %lld\n", instructions);
    //printf("  Cache Misses: %lld\n", cache_misses);
    //printf("  IPC: %.2f\n", (double)instructions / cycles);
    //printf("  Retiring: %.2f\n", retiring);
    //printf("  Bad spec:  %.2f\n", bad_spec); 
    //printf("  Frontend bound: %.2f\n", frontend_bound);
    //printf("  Backend bound: %.2f\n", backend_bound);
    //printf("     Memory bound: %.2f\n", mem_bound);
    //printf("     Core bound: %.2f\n", core_bound);
    
    // Close FDs
    close(counters.cycles_fd                 );
    close(counters.instr_fd                  );
    close(counters.cache_fd                  );
    close(counters.topdown_slots_fd          );
    close(counters.uops_retired_fd           );
    close(counters.topdown_backend_slots_fd  );
    close(counters.topdown_bad_spec_slots_fd );
    close(counters.topdown_memory_slots_fd   );
}
