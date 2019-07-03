#!/usr/bin/env python
"""
  Usage:
    backtranslate.py <SEQUENCE> <GC_LOWER_BOUND> <GC_UPPER_BOUND>
"""
from docopt import docopt
from backtranslate import backTranslate

options = docopt(__doc__)
SEQUENCE = options["<SEQUENCE>"]
GC_LOWER_BOUND = float(options["<GC_LOWER_BOUND>"])
GC_UPPER_BOUND = float(options["<GC_UPPER_BOUND>"])

codonSeq = backTranslate(SEQUENCE,GC_LOWER_BOUND,GC_UPPER_BOUND)

print(' '.join(codonSeq))
