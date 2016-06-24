# Parallel bash


## Python script to run bash script in parallel.


> **Note**:
> It is assumed that the bash script uses named options (e.g., `<script> -c 4`), but
> not positional arguments. If you do not know how to pass named options to bash
> scripts, google "bash getopts".



### Arguments for this script

#### Required

`-s` or `--shellScript` *`FILE`*
- *`FILE`*: Shell script to run in parallel.

`-t` or `--totalCores` *`INT`*
- *`INT`*: How many total cores are available?


#### Optional

`-p` or `--perSampCores` *`INT`*
- *`INT`*: How many cores should the script use per sample? Must be a factor of
  `totalCores`.
- Defaults to 1.


`shellOpts` *`STRING`* `[`*`STRING`* `...]`
- *`STRING`*: Options to input to shell script. The letter associated with the
  option should be in the beginning and separated from the arguments by a ':'; if no
  arguments are required for that option, still include the ':'. Arguments must meet
  one of the following criteria:

    1. Length = 1, in which case they will be the same for all iterations of script.
    2. Length = the number of iterations. If this is the case, items must be
       separated by ';'.

    For multiple options, provide multiple separate strings.

    For example, if you wanted to run a parallel version of the following:

    ```bash
    ./my_script.sh -x A -y 1 -z && ./my_script.sh -x A -y 2 -z
    ```
    ... you would run the following command for this script (assuming each iteration
    of the `my_script.sh` uses 4 cores):
    ```bash
    parallel_bash.py -s my_script.sh -t 8 -p 4 'x:A' 'y:1;2' 'z:'
    ```

- Defaults to "" (no options).
