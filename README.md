PBBCache is a cache-partitioning simulator that relies on offline-collected application performance data (e.g., instructions per 
cycle, memory bandwidth consumption, etc.) to determine the degree of throughput, fairness or other relevant metrics for a workload
under a particular partitioning algorithm. A key feature of the simulator is its ability to efficiently determine the optimal
cache-partitioning solution for different optimizations objectives (e.g. throughput or fairness) using a parallel branch-and-bound algorithm.

The main purpose of the simulator is twofold. Firstly, it has been designed as a tool to efficiently assess the
real potential of existing cache-partitioning approaches and to identify their limitations, 
by comparing their effectiveness with that of the optimal solution provided by the simulator. 
Secondly, the ability to study the optimal solution for various optimization metrics makes the simulator
a powerful tool to aid in the design of new cache-partitioning techniques.

## Project dependencies

The project relies on Python 2.7.9 or greater and uses the following Python modules:

* matplotlib
* pandas
* numpy
* sympy
* ipyparallel
* sklearn
* datetime
* pytz

Note: After installing ipyparallel, the binaries that start the cluster are not included in the path sometimes. To do so you must find where are they installed and add them to the PATH:

	find ~ -name 'ipcluster'
	export PATH=$PATH:$HOME/.local/bin

For more information on this library, recur to the [official documentation](https://ipyparallel.readthedocs.io/en/latest/)

## Getting started

A regular use of the simulator would include the following steps:

Include the project common functions module:

```
$ cd simulator; . shrc
```

Start 4 slave (engine) processes:

```
$ ipcluster start -n 4 --daemonize
```

Launch the simulator:

```
$ ./test/sim.py -s input/metrics.csv input/workloads.csv -a opt-stp-bf,ucp -f table -b 35000 -p -m stalls -r 1
W#   Algorithm           BenchID Name                 Cluster       Mask/NR_WAYS    SlowdownNB/ANTT   STP               Slowdown         
W1   opt-stp-bf          1       sphinx306            1             0x700(3)        1.07895           0.91197           1.09653          
W1   opt-stp-bf          2       lbm06                2             0x80(1)         1.0097            0.957             1.04493          
W1   opt-stp-bf          3       libquantum06         3             0x40(1)         1.16813           0.807             1.23916          
W1   opt-stp-bf          4       applu00              4             0x20(1)         1.04627           0.89354           1.11914          
W1   opt-stp-bf          5       soplex06             5             0x1e(4)         1.18749           0.81842           1.22186          
W1   opt-stp-bf          6       milc06               6             0x1(1)          1.02608           0.88443           1.13067          
W1   opt-stp-bf          7       OVERALL              #88           3.087099s       1.14205           5.27237           1.18588          
W1   ucp                 1       sphinx306            1             0x600(2)        1.11045           0.8617            1.1605           
W1   ucp                 2       lbm06                2             0x180(2)        1.00697           0.96071           1.0409           
W1   ucp                 3       libquantum06         3             0x40(1)         1.16813           0.80426           1.24338          
W1   ucp                 4       applu00              4             0x20(1)         1.04627           0.89029           1.12323          
W1   ucp                 5       soplex06             5             0x1e(4)         1.18749           0.81715           1.22376          
W1   ucp                 6       milc06               6             0x1(1)          1.02608           0.87991           1.13648          
W1   ucp                 7       OVERALL              #1            0.001563s       1.15471           5.21402           1.19452 
```          


There are many available command line options, in this example we have used the most basic ones:

* -s specifies a metrics file, that contains a table where each row store the values of various runtime metrics
gathered with performance monitoring counters when a specific application runs alone with a fixed number of cache ways.
* The following parameter (required always) is the workloads file, which specifies the different workloads (one per line) 
that will be considered for the simulation; each workload is encoded as a comma-separated list of program names.
* -a defines the set of partitioning schemes to simulate.
* -f is the output mode
* -p activates the parallel execution of the branch and bound algorithm to calculate the optimal partitioning.
* -m specifies the bandwidth model to be applied, in this case the one based on the stall metric.
* -r allows to only execute a defined range of workloads. It is defined as a comma-separated list and allows to define ranges with - (e.g., if we have 10 workloads and define the flag -r 1,3,5-7,9- we will end up running the workloads 1,3,5,6,7,9,10).

To find all of them you can run the following command:

```
$ ./test/sim.py -h
usage: sim.py [-h] [-s SUMMARY_FILE] [-b MAX_BANDWIDTH] [-p] [-P] [-H] [-C]
              [-W] [-L] [-m BW_MODEL] [-a ALGORITHMS] [-f FORMAT] [-d]
              [-r USE_RANGE] [-w FORCE_WAYS] [-O key=val]
              WORKLOAD_FILE

Test main file for simulator (UCP, Whirlpool, etc.)

positional arguments:
  WORKLOAD_FILE         a workload file

optional arguments:
  -h, --help            show this help message and exit
  -s SUMMARY_FILE, --summary-file SUMMARY_FILE
                        File with the data collected offline for every
                        benchmark
  -b MAX_BANDWIDTH, --max-bandwidth MAX_BANDWIDTH
                        Max bandwidth observed for the target machine
  -p, --parallel        Enable parallel search for the optimal
  -P, --parsim          Launch all simulations in parallel
  -H, --harness         Use harness workload file as input
  -C, --generate-chart  Enable the STP-vs-UNF chart generator
  -W, --show-window     Show matplotlib window with chart (must be used along
                        with -C)
  -L, --list-algorithms
                        List algorithms implemented in the simulator.
  -m BW_MODEL, --bw-model BW_MODEL
                        Select BW model to be applied (['simple', 'stalls'])
  -a ALGORITHMS, --algorithms ALGORITHMS
                        Comma-separated algorithms (['opt-stp', 'yu-petrov',
                        'opt-unf', 'ucp', 'ucp-slowdown', 'kpart',
                        'whirlpool', 'whirlpool-c', 'equal-part', 'optc-stp',
                        'optc-unf', 'opt-unf-np', 'opt-stp-np', 'yu-petrov',
                        'user', 'kpart-opt', 'opt-stp-bf', 'opt-unf-bf',
                        'optc-stp-bf', 'optc-unf-bf'])
  -f FORMAT, --format FORMAT
                        Format style for output
  -d, --debugging       Enable debugging mode (assume normal input files) and
                        keep workload 0
  -r USE_RANGE, --use-range USE_RANGE
                        Pick selected workloads only by specifying a range
  -w FORCE_WAYS, --force-ways FORCE_WAYS
                        Force number of ways of the cache (it must be smaller
                        or equal than the maximum number of ways inferred from
                        input file)
  -O key=val, --option key=val
                        Use generic algorithm-specific options
```


## Project contributors

* Jorge Casas Hernan (<jorcasas@ucm.es>)
* Adrian Garcia Garcia (<adriagar@ucm.es>)
* Juan Carlos Saez Alcaide (<jcsaezal@ucm.es>)

