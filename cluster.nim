import algorithm
import math
import itertools
import hts/bam

# Data structure storing information about each read that looks like an STR
type tread* = object
  tid*: int32
  position*: uint32
  repeat*: array[6, char]
  flag*: Flag
  split*: int8
  mapping_quality*: uint8
  repeat_count*: uint8
  read_length*: uint8

type Cluster* = object
  reads*: seq[tread]

proc posmed(cl:Cluster, n:int=5): uint32 =
  ## posmed is the median of the first n positions in the cluster
  let mid = int(min(5, cl.reads.len) / 2 + 0.5)
  return cl.reads[mid].position

proc bytidrep(t:tread): tuple[repeat:array[6, char], tid:int32] =
  return (t.repeat, t.tid)

# Sorts the reads by chromosome (tid) then repeat unit, then by position
proc tread_cmp(a: tread, b:tread): int =
  if a.tid != b.tid: return cmp(a.tid, b.tid)
  for i in 0..<6:
    if a.repeat[i] != b.repeat[i]:
      return cmp(a.repeat[i], b.repeat[i])
  return cmp(a.position, b.position)

proc trim(cl:var Cluster, max_dist:uint32) =
  # drop stuff from start of cluster that is now outside the expected distance
  var lo = cl.posmed(5) - max_dist
  while len(cl.reads) > 0 and cl.reads[0].position < lo:
    cl.reads = cl.reads[1..cl.reads.high]

iterator cluster*(tandems: var seq[tread], max_dist:uint32, min_supporting_reads:int=5): Cluster =
  tandems.sort(tread_cmp)

  for group in groupby(tandems, bytidrep):
    # reps are on same chromosome and have same repeat unit
    var reps: seq[tread] = group.v
    var i = 0
    var c:Cluster
    while i < reps.len:
      # start a new cluster
      var it = reps[i]
      c = Cluster(reads: @[it])
      for j in i+1..reps.high:
        # add any tread that's close enough
        if reps[j].position <= c.posmed(5) + max_dist:
          c.reads.add(reps[j])
          continue

        # remove stuff (at start of cluster) that's now too far away.
        c.trim(max_dist)
        if c.reads.len >= min_supporting_reads:
          yield c
        # increment i to past last j and break out of this cluster
        i = j + 1
        break

    c.trim(max_dist)
    if c.reads.len >= min_supporting_reads:
      yield c
