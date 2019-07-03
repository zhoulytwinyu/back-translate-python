#!/usr/bin/env python
from collections import Counter
from bisect import bisect_left,bisect_right
from random import sample, shuffle
import pkg_resources
CODON_FILE = pkg_resources.resource_filename(__name__, 'data/standard_codon.tsv')

def memoize_one(func):
  memory = {"args":None,
            "result":None}
  def new_func(*args):
    if memory["args"] != args:
      memory["args"] = args
      memory["result"] = func(*args)
    return memory["result"]
  return new_func

@memoize_one
def getAminoAcidToCodonMap():
  ret = {}
  fh = open(CODON_FILE)
  for l in fh:
    aa,cds = l.strip().split('\t')
    ret[aa]=cds.split(',')
  return ret

@memoize_one
def getAllAminoAcidSet():
  codonMap = getAminoAcidToCodonMap()
  ret = set(codonMap.keys())
  return ret

@memoize_one
def getCodonGCMap():
  ret = {}
  codonMap = getAminoAcidToCodonMap()
  for cds in codonMap.values():
    for cd in cds:
      ret[cd] = getGCContent(cd)
  return ret

def getGCContent(NTSequence):
  gcCount = 0
  GC = {'G','C'}
  for i in NTSequence:
    if i in GC:
      gcCount+=1
  return gcCount/len(NTSequence)

@memoize_one
def getAminoAcidInfoMap():
  ret = {}
  codonMap = getAminoAcidToCodonMap()
  for aa, cds in codonMap.items():
    codonGCPair = [(cd,getGCContent(cd)) for cd in cds]
    possibleGCs = list(set([p[1] for p in codonGCPair]))
    highGC = max(possibleGCs)
    lowGC = min(possibleGCs)
    GCToCodonMap = {}
    for cd,gc in codonGCPair:
      tmpList = GCToCodonMap.get(cd,[])
      tmpList.append(cd)
      GCToCodonMap[gc] = tmpList
    ret[aa] = { "GCToCodonMap": GCToCodonMap,
                "possibleGCs": possibleGCs,
                "highGC": highGC,
                "lowGC": lowGC
                }
  return ret

def validateAminoAcidSequence(AASeq):
  allAA = getAllAminoAcidSet()
  if set(AASeq).issubset(allAA):
    return true
  else:
    return false

def pickCodon(AA,targetGC):
  infoMap = getAminoAcidInfoMap()[AA]
  possibleGCs = infoMap["possibleGCs"]
  GCToCodonMap = infoMap["GCToCodonMap"]
  diff = [abs(gc-targetGC) for gc in infoMap["possibleGCs"]]
  minIdx = diff.index(min(diff))
  chosenGC = infoMap["possibleGCs"][minIdx]
  chosenCodon = sample(GCToCodonMap[chosenGC],1)[0]
  return chosenCodon

def backTranslate_singleTarget(AASeq,targetGC):
  condomGCMap = getCodonGCMap()
  infoMap = getAminoAcidInfoMap()
  seqLength = len(AASeq)
  AACount = Counter(AASeq)
  # Split amino acid into flexible and unflexible group
  # Unflexible group contains those AA that cannot by adjust GC around targetGC
  AACount_unflexible = {}
  AACount_flexible = {}
  for aa,count in AACount.items():
    if  infoMap[aa]["highGC"]<targetGC or \
        infoMap[aa]["lowGC"]>targetGC:
      AACount_unflexible[aa] = count
    else:
      AACount_flexible[aa] = count
  # Set up
  curTargetGC = targetGC
  sumGC=0
  backTranslationMap = {}
  i=0
  # First, greedily choose codon for unflexible amino acid at random
  pickSequence = []
  for aa,count in AACount_unflexible.items():
    pickSequence.extend([aa]*count)
  shuffle(pickSequence)
  for aa in pickSequence:
    curTargetGC = (targetGC*seqLength-sumGC)/(seqLength-i)
    # Pick codon for AA with the given target GC
    codon = pickCodon(aa,curTargetGC)
    tmpList = backTranslationMap.get(aa,[])
    tmpList.append(codon)
    backTranslationMap[aa] = tmpList
    # Update target GC
    gc = condomGCMap[codon]
    sumGC = gc+sumGC
    i+=1
  # Second, greedily choose codon for flexible amino acid at random
  pickSequence = []
  for aa,count in AACount_flexible.items():
    pickSequence.extend([aa]*count)
  shuffle(pickSequence)
  for aa in pickSequence:
    curTargetGC = (targetGC*seqLength-sumGC)/(seqLength-i)
    # Pick codon for AA with the given target GC
    codon = pickCodon(aa,curTargetGC)
    tmpList = backTranslationMap.get(aa,[])
    tmpList.append(codon)
    backTranslationMap[aa] = tmpList
    # Update target GC
    gc = condomGCMap[codon]
    sumGC = gc+sumGC
    i+=1
  # Shuffle codon to provide some randomness
  for codons in backTranslationMap.values():
    shuffle(codons)
  # Go through AASeq and backtranslate using backTranslationMap
  ret = []
  for aa in AASeq:
    ret.append(backTranslationMap[aa].pop())
  return ret

def backTranslate(AASeq,minGC,maxGC):
  targetGC = (minGC+maxGC)/2
  return backTranslate_singleTarget(AASeq,targetGC)

if __name__=="__main__":
  codonSeq = backTranslate("MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*",0.5,0.55)
  print(' '.join(codonSeq))
