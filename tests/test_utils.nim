import unittest
include strpkg/utils

suite "utils suite":

  test "test get_tid":
    var targets = @[Target(name: "chr1", tid: 0, length: 10000)]
    check get_tid("chr1", targets) == 0

  test "test median":
    check median_depth(@[1,2,2]) == 2
    check median_depth(@[2000]) == 1047 #Values greater than 1047 will be set to 1047 to avoid overflow
    #echo median_depth(@[-11,2,2]) # This will throw an inded out of range error

  test "test get_chrom":
    var targets = @[Target(name: "chr1", tid: 0, length: 10000)]
    check get_chrom(0, targets) == "chr1"

  test "test mode":
    check mode([1,2,3,3,3,2]) == 3
    check mode([1.1,2.2,3.1,3.1,3.4,2.5]) == 3.1

  test "test get most frequent in CountTable":
    var vals = [1,2,3,3,3,2]
    var t: CountTable[int]
    for x in vals:
      t.inc(x)
    check most_frequent(t, 2) == [3, 2]

  test "test isNaN":
    var x: float = NaN
    var y = 1.5
    check x.isNaN
    check not y.isNaN
