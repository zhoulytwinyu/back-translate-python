# back-translate python packakge

Back translate peptide sequences into nucleotide sequences within a given GC content range.

## Install
```
pip install git+https://github.com/zhoulytwinyu/back-translate-python.git
```

## Uninstall
```
pip unstall back-translate 
```

## Usage
This package includes a script which can be invoked in bash.
```
#!/usr/bin/bash
back-translate.py <SEQUENCE> <GC_LOWER_BOUND> <GC_UPPER_BOUND>
```
* <SEQUENCE> must be a 1-letter coded, non-degererative amino acid sequence. "*" may be used as the stop.
* <GC_LOWER_BOUND> <GC_UPPER_BOUND> are floating point number within 0.0 - 1.0.

## Algorithm highlight (Greedy algorithm)

### Pick codons using a out-of-order strategy
* Collect and count every amino acid in given peptide sequence. 
* Split the all the amino acids into two groups (being strategic here).
  * Unflexible group: Those whose codon GC contents are all too high or too low for the target GC range.
  * Flexible group: Those whose codong GC contents spans the target GC range and may help tune GC content to our target GC.
* We pick codon for those unflexible amino acids first and then pick codon for flexible amoni acids. We note down the picked codons for each amino acid.
  * We shuffle the amino acids within the two groups first to give them some randomness. It help avoid some corner cases where the algorithm does not work well.
  * When picking codon for amino acids, we note down the GC content of all the picked amino acids and try to choose a codon that brings GC content closer to our target GC range (being greedy here).
### Put the codons back to align with the peptide sequence
* Iterate through the given peptide sequence, choose at random (smoother local GC content) without replacement from our picked codons for the amino acids.

## Caution
This algorithm is essentially probabilistic and may under certain circumstances not able to generate nucleotide sequences within given GC range. In fact, the algorithm only try to optimize sequence toward target GC range; it does not double check whether the generated sequence is indeed within that range. In most cases, it works fine.
