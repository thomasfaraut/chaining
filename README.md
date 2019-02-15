
# chaining
A lis-like chaining algorithm.

```
Usage:
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
```


