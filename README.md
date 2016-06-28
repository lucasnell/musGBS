`musGBS`
=====

Pipelines for genotyping by sequencing (GBS) data on mice
-----


These scripts were designed solely for my drive on the HPC cluster at UGA.


### Shell scripts

The order of scripts to be run is as follows:

```
initial_processing.sh –> align_sort_split.sh ─┬─> merge_runs.sh ─> make_stacks.sh ─┐
                                         ... ─┤     [catalog] ┄┄> call_stacks.sh <─┘
        [same as above, for additional runs] ─┘               { ?? } <─┘
```

### Python script

The file `parallel_bash.py` uses Python's `multiprocessing` package to run bash commands
in parallel. I created this file because Stacks does not appear to perform very well
using multiple cores. The Python file allows me to run multiple files at a time, rather
than using multiple cores per file.
