#!/usr/bin/env python

"""Usage:
chaining.py -i INPUT -f FASTA [-ho FILE] [--scaffold=<sn>] [--stranded] [--min_scaffold_size=<sn>] [--min_align_length=<an>] [--format=<format>] [--log=<logfile>] [--debug] [--quiet | --verbose] [--head=<hn>]

-h --help    show this
-i --input input tab alignment file (from lastal)
-f --fasta scaffolds fastafile
--min_scaffold_size=<sn>  min scaffold size [default: 20000]
--min_align_length=<an>  min align size size [default: 200]
--format=<format>   lasttab, bedpe [default: bedpe]
-o FILE  specify output file [default: ./testchains.out]
--log FILE    log file [default: ./testchains.log]
--quiet      print less text
--verbose    print more text
--debug      debug mode
--head=<hn>  print only first hn scaffolds
--stranded  requires chains to respect the orientations
--scaffold=<sca>

"""
# -o FILE  specify output file [default: ./testchains.txt]

import sys
from docopt import docopt
from Alignments import Alignment, AlignSet
from Alignments import read_alignfile, get_scaffold_N_cumsum

verbose = False


# A small wrapper to print to stderr
def eprint(*args, **kwargs):
    if verbose:
        print(*args, file=sys.stderr, **kwargs)


if __name__ == '__main__':
    # __doc__ contient automatiquement la docstring du module
    # en cours
    args = docopt(__doc__, version='0.1')

    verbose = args['--verbose']
    debug = args['--debug']
    head = args['--head']

    output_file = args['-o']
    log_file = args['--log']
    format = args['--format']

    stranded = args['--stranded']
    sel_scaffold = args['--scaffold']

    # if verbose:
    #     eprint(args)

    mapping, scaffolds = \
        read_alignfile(args['--input'],
                       min_scaffold_size=int(args['--min_scaffold_size']),
                       min_align_length=int(args['--min_align_length']))
    fastafile = args['--fasta']
    eprint("#Now computing chains")
    eprint("%d scaffolds" % (len(mapping)))

    # The output file and log file
    with open(output_file, "w") as out, open(log_file, "w") as log:
        i = 1
        for scaffold in sorted(mapping, key=lambda x: scaffolds[x],
                               reverse=True):
            if sel_scaffold is not None and scaffold != sel_scaffold:
                continue
            hsps = AlignSet(mapping[scaffold], debug=debug)
            cumsum = get_scaffold_N_cumsum(scaffold, fastafile)

            sequences = AlignSet(hsps.strict_chaining(qgap_cutoff=5000,
                                                      sgap_cutoff=5000,
                                                      cum_N_sum=cumsum,
                                                      stranded=True
                                                      ),
                                 score_cutoff=100,
                                 debug=debug)
            lischains = AlignSet(sequences.chaining(qgap_cutoff=50000,
                                                    sgap_cutoff=50000,
                                                    cum_N_sum=cumsum,
                                                    stranded=True,
                                                    dropoff=-30000,
                                                    score_cutoff=10000
                                                    ),
                                 score_cutoff=15000,
                                 debug=debug)

            if debug:
                lischains.qdisplay()

            if stranded:
                topchains = lischains
            else:
                chains = AlignSet(lischains.strict_chaining(qgap_cutoff=350000,
                                                            sgap_cutoff=350000,
                                                            cum_N_sum=cumsum,
                                                            stranded=False),
                                  score_cutoff=20000,
                                  size_cutoff=20000)

                topchains = AlignSet(chains.strict_chaining(qgap_cutoff=350000,
                                                            sgap_cutoff=350000,
                                                            cum_N_sum=cumsum,
                                                            stranded=False),
                                     score_cutoff=20000)

            num_chain, coverage = topchains.filter(scaffolds[scaffold],
                                                   coverage=0.99,
                                                   score=10000,
                                                   length_cutoff=20000)
            print("%s\t%d\t%3.3f" % (scaffold, num_chain, coverage), file=log)
            eprint("#%s\t%d\t%d\t%3.3f" % (scaffold, scaffolds[scaffold],
                                           num_chain, coverage))
            topchains.qdisplay(format=format, filter=True, file=out)
            if head and i >= int(head):
                break
            i += 1
