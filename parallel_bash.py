#!/usr/local/apps/anaconda/3-2.2.0/bin/python

"""Python script to run bash script in parallel."""


'''Shabang for non-cluster:
#!/usr/bin/env python3
'''

import subprocess as sp
import os
import stat
import argparse as ap
from multiprocessing import Pool


# =================================
# Function to run one iteration
# =================================

def runOne(script):
    """Run one iteration.
    Args:
        script (str): String of bash script to run.
    """
    sp.call(script, shell = True)
    return


if __name__ == '__main__':

    # =================================
    # Setting up parser
    # =================================

    Parser_help = \
        """Run shell script in parallel.
        It is assumed that the bash script uses named options (e.g., `<script> -c 4`),
        but not positional arguments.
        If you do not know how to pass named options to bash scripts, search for
        'getopts'."""

    Parser = ap.ArgumentParser(description = ' '.join(Parser_help.split()))

    Parser.add_argument('-s', '--shellScript', required = True,
                        help = 'Shell script to run in parallel.')

    perSamp_help = \
        """How many cores should the script use per sample? Must be a factor of
        `totalCores`. Will be provided to shell script as the input for `-c` (e.g.,
        `<your_script> -c 4`). Defaults to 1."""
    Parser.add_argument('-p', '--perSampCores', type = int,
                        required = False, default = 1,
                        help = ' '.join(perSamp_help.split()))

    Parser.add_argument('-t', '--totalCores', type = int,
                        required = False, default = 1,
                        help = "How many total cores are available? Defaults to 1.")

    shellOpts_help = \
        """Options to input to shell script.
        The letter associated with the option should be in the beginning and separated
        from the arguments by a ':'; if no arguments are required for
        that option, still include the ':'.
        To separate multiple options, use ';'.
        Arguments must be (1) of length 1, in which case they will be the same for all
        iterations of script, or (2) of the same length as the number of iterations
        and separated by ','. For example, if you wanted to run a parallel version of
        `./my_script.sh -x A -y 1 -z && ./my_script.sh -x A -y 2 -z`, you would provide
        the following string here: 'x:A;y:1,2;z:'. Defaults to '' (no options)."""
    Parser.add_argument('-o', '--shellOpts', default = '',
                        help = ' '.join(shellOpts_help.split()))

    # =================================
    # Parsing command line inputs
    # =================================

    args = vars(Parser.parse_args())
    cores = args['totalCores']
    cores_i = args['perSampCores']
    shellOpts = args['shellOpts']
    if args['shellScript'].__class__ == list:
        shellScript = os.path.expanduser(args['shellScript'][0])
    else:
        shellScript = os.path.expanduser(args['shellScript'])

    assert cores % cores_i == 0, '`perSampCores` be a factor of `totalCores`.'

    # Separate `shellOpts`
    if shellOpts != '':
        assert not (shellOpts.count(':') - shellOpts.count(';')) < 1, \
            'Not all options include a ":".'
        assert not (shellOpts.count(':') - shellOpts.count(';')) > 1, \
            'Not all options are separated by a ";".'

    try:
        optDict = {'-' + x.split(':')[0]: x.split(':')[1].split(',')
                   for x in shellOpts.split(';')}
    except IndexError:
        bashScripts = ['%s -c %i' % (shellScript, cores_i)]
    else:
        maxLen = max([len(x) for x in optDict.values()])
        # Expand all option arguments to be the same length
        for i in optDict.keys():
            if len(optDict[i]) == 1:
                optDict[i] = optDict[i] * maxLen
            elif len(optDict[i]) == maxLen:
                continue
            else:
                raise AttributeError(
                    'Argument lengths must be 1 or the same as all others > 1.')
        # Iterate from 1 to `maxLen`, creating strings of script with options
        bashScripts = []
        for j in range(maxLen):
            args_i = [' '.join([' '.join([i, optDict[i][j]])
                                for i in optDict.keys()])][0]
            bashScripts += ['%s -c %i %s' % (shellScript, cores_i, args_i)]

    # Make sure shell script is executable
    perms = os.stat(shellScript)
    os.chmod(shellScript, perms.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    # =================================
    # Run all iterations in parallel
    # =================================

    numThreads = int(cores / cores_i)
    with Pool(processes = numThreads) as pool:
        pool.map(runOne, bashScripts)
