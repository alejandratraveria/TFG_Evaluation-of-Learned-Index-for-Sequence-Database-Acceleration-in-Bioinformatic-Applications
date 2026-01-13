# Genome Analysis Benchmark: Seeding Phase

## Overview

This benchmark is designed to evaluate the seeding phase's performance of the genome analysis, focusing on minimizer lookup and index access. It supports two indexing strategies:
- **Classic Hash Table:** Traditional hash-based lookup.
- **RMI (Recursive Model Index):** Learned index for fast minimizer search.

The benchmark allows you to test different optimizations, including vectorization, batching, and prefetching, and collects detailed metrics about index access and microarchitectural behaviour.

### Main Steps

1. **Index and RMI Generation:**  
	For a new dataset, you must first build the minimizer index and the RMI structure.
2. **Benchmark Execution:**  
	You can run the benchmark using either the raw sequence (which will rebuild the index each time) or a precomputed index file (`.mmi`).
3. **Metrics Collection:**  
	The benchmark collects various performance metrics and outputs detailed CSV files for further analysis.

## Usage

### 1. Build Index and RMI (First Time for a Dataset)

```bash
source scripts/build_rmi.sh
./scripts/create_index_rmi.sh <reference_file> <preset>
```
- `<reference_file>`: Path to your genome FASTA file.
- `<preset>`: Preset for minimap2 (e.g., `map-ont`).

This step generates the minimizer index and the RMI files.

### 2. Compile the Benchmark

```bash
make -f LISAHash_Makefile clean
make -f LISAHash_Makefile bench_mm_queries lhash=1/0 vect=1/0 pref=1/0 batch=1/0 key_type=uint64 metrics=1/0 mm_info=1/0
```

**Options:**
- `lhash=0`: Use classic hash table.
- `lhash=1`: Use RMI for minimizer lookup.
RMI case optimizations:
- `vect=1`: Enable vectorization for binary search in RMI.
- `batch=1`: Enable batching for search queries.
- `pref=1`: Enable prefetching for search queries.
- `key_type=uint64`: Key type for minimizers.
Collection of metrics:
- `metrics=1`: Enable the collection of metrics and perf counters (if enabled).
- `mm_info=1`: Collect index access metrics and generate output files.
- `topdown=1`: Collect the topdown metrics to obtain a performance stack of the execution.

### 3. Run the Benchmark
```bash
./bench-mm-queries.exe -t <num_threads> -l <counter_group> <index_file> <queries_file>
```
- `<index_file>`: Either the reference FASTA or the `.mmi` index file.
- `<queries_file>`: FASTA containing query sequences.
- `-t <num_threads>`: Number of threads to use (default: 1).
- `-l <counter_group>`: Microarchitectural counter group to measure. Options:
	- `-1`: Disable measurement
	- `0`: Group 0
	- `1`: Group 1
	- `2`: Group 2
	- `3`: Group 3

You can run with either:
- The raw sequence file (will rebuild the index each time), or
- The precomputed index file (`<reference_file>.mmi`).

```bash
./bench-mm-queries.exe <index_file> <queries_file>
```
- `<index_file>`: Either the reference FASTA or the `.mmi` index file.
- `<queries_file>`: FASTA containing query sequences.

## Output Files

- `minimizers_output.csv`: Metrics for each minimizer searched.
- `L1_PARAMETERS_errors.csv`: Error for each RMI leaf node.
- `mm_metrics_g<group>_output.csv`: Perf counters of the selected group and metrics.

## Notes

- For best performance, use the `.mmi` index file for repeated runs.
- The benchmark prints the configuration used at startup, including which index type and optimizations are enabled.

## Getting started and execution example
```bash
git clone <REPO URL>
git submodule update --init --recursive
sudo apt install curl
source scripts/build_rmi.sh
./scripts/create_index_rmi.sh <absolute path reference_file> map-ont
make -f LISAHash_Makefile clean
make -f LISAHash_Makefile bench_mm_queries lhash=1 vect=1 pref=1 batch=1 key_type=uint64 metrics=1 mm_info=1
./bench-mm-queries.exe <index_file> <queries_file>
```
## Switching model
In order to change the used model you have to change the file named `two_layer.rs`, you have different `two_layer_name_of_model.rs` files to choose from.
Also you need to select the `modify_generated_code.sh` that matches with the model you want to use.

Appart from the rust files you also need to change the code for the `RMI.h` file to match the number of parameters the selected model has.