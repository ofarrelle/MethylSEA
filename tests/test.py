from subprocess import check_call, CalledProcessError
import os
import os.path as op
wd = op.dirname(op.realpath(__file__))
os.chdir(wd)

def rm(f):
    try:
        os.unlink(f)
    except OSError:
        pass

assert op.exists('../MethylSEA')
MPath = os.path.abspath('../MethylSEA')
check_call([MPath, '--version'])

rm('test1_CpG.bedGraph')
check_call([MPath, 'extract', 'ct100.fa', 'ct_aln.bam', '-q', '2', '-o', 'test1'])
assert op.exists('test1_CpG.bedGraph')
lines = sum(1 for _ in open('test1_CpG.bedGraph'))
assert lines == 1
rm('test1_CpG.bedGraph')

rm('test2_CpG.bedGraph')
check_call([MPath, 'extract', 'cg100.fa', 'cg_aln.bam', '-q', '2', '-o', 'test2'])
assert op.exists('test2_CpG.bedGraph')
for line in open('test2_CpG.bedGraph'):
    print(line)
lines = sum(1 for _ in open('test2_CpG.bedGraph'))
assert lines > 1
rm('test2_CpG.bedGraph')

# should be none with q > 10
rm('test3_CpG.bedGraph')
check_call([MPath, 'extract', 'cg100.fa', 'cg_aln.bam', '-q', '10', '-o', 'test3'])
assert op.exists('test3_CpG.bedGraph')
lines = sum(1 for _ in open('test3_CpG.bedGraph'))
assert lines == 1
rm('test3_CpG.bedGraph')

# Test the new methylKit option
rm('test4_CpG.methylKit')
rm('test4_CHG.methylKit')
rm('test4_CHH.methylKit')
check_call([MPath, 'extract', '--methylKit', '--CHH', '--CHG', 'cg100.fa', 'cg_aln.bam', '-q', '2', '-o', 'test4'])
assert op.exists('test4_CpG.methylKit')
lines = sum(1 for _ in open('test4_CpG.methylKit'))
assert lines > 1
rm('test4_CpG.methylKit')
lines = sum(1 for _ in open('test4_CHG.methylKit'))
assert lines == 1
rm('test4_CHG.methylKit')
assert op.exists('test4_CHH.methylKit')
lines = sum(1 for _ in open('test4_CHH.methylKit'))
assert lines == 2
rm('test4_CHH.methylKit')

# Check that --minDepth is working, which means there should be no called sites
rm('test5_CpG.bedGraph')
check_call([MPath, 'extract', '--minDepth', '2', 'cg100.fa', 'cg_aln.bam', '-q', '2', '-o', 'test5'])
assert op.exists('test5_CpG.bedGraph')
lines = sum(1 for _ in open('test5_CpG.bedGraph'))
assert lines == 1
rm('test5_CpG.bedGraph')

# Check that --ignoreFlags is working, which means that there are now called sites
rm('test6_CpG.bedGraph')
check_call([MPath, 'extract', '--ignoreFlags', '0xD00', 'cg100.fa', 'cg_aln.bam', '-q', '2', '-o', 'test6'])
assert op.exists('test6_CpG.bedGraph')
lines = sum(1 for _ in open('test6_CpG.bedGraph'))
assert lines == 49
rm('test6_CpG.bedGraph')

# Check that --requireFlags is working
rm('test7_CpG.bedGraph')
check_call([MPath, 'extract', '--requireFlags', '0xD00', 'cg100.fa', 'cg_aln.bam', '-q', '2', '-o', 'test7'])
assert op.exists('test7_CpG.bedGraph')
lines = sum(1 for _ in open('test7_CpG.bedGraph'))
assert lines == 49
rm('test7_CpG.bedGraph')

# Check absolute trimming bounds.
# --nOT 50,50,40,40 on the single kept OT pair (both 100bp, pos 0): read1 is
# trimmed away entirely (50+50), read2 keeps indices 40-59, giving CpG calls at
# even positions 40,42,...,58 = 10 data lines + 1 header = 11. The prior value of
# 12 reflected the pre-f3f16b5 off-by-one that failed to trim the last base.
rm('test8_CpG.bedGraph')
check_call([MPath, 'extract', '--nOT', '50,50,40,40', 'cg100.fa', 'cg_aln.bam', '-q', '2', '-o', 'test8'])
assert op.exists('test8_CpG.bedGraph')
lines = sum(1 for _ in open('test8_CpG.bedGraph'))
assert lines == 11
rm('test8_CpG.bedGraph')

# Check variant filtering (there are 49 lines otherwise)
rm('test9_CpG.bedGraph')
check_call([MPath, 'extract', '-p', '1', '-q', '0', '-o', 'test9', '--minOppositeDepth', '3', '--maxVariantFrac', '0.25', 'cg100.fa', 'cg_with_variants.bam'])
assert op.exists('test9_CpG.bedGraph')
lines = sum(1 for _ in open('test9_CpG.bedGraph'))
assert lines == 48
rm('test9_CpG.bedGraph')

# Check conversion efficiency. 2 read pairs, one mostly converted
# By default, 1 read is MAPQ filtered, another is kept
rm('test10_CpG.bedGraph')
check_call([MPath, 'extract', '-o', 'test10', 'chgchh.fa', 'chgchh_aln.bam'])
assert op.exists('test10_CpG.bedGraph')
lines = sum(1 for _ in open('test10_CpG.bedGraph'))
assert lines == 2
rm('test10_CpG.bedGraph')

# Ensure 2 reads/positions are covered by changing MAPQ
rm('test11_CpG.bedGraph')
check_call([MPath, 'extract', '-o', 'test11', '-q', '5', 'chgchh.fa', 'chgchh_aln.bam'])
assert op.exists('test11_CpG.bedGraph')
lines = sum(1 for _ in open('test11_CpG.bedGraph'))
assert lines == 3
rm('test11_CpG.bedGraph')

# Only 1 read has a conversion efficiency >=0.9
rm('test12_CpG.bedGraph')
check_call([MPath, 'extract', '-o', 'test12', '-q', '5', '--minConversionEfficiency', '0.9', 'chgchh.fa', 'chgchh_aln.bam'])
assert op.exists('test12_CpG.bedGraph')
lines = sum(1 for _ in open('test12_CpG.bedGraph'))
assert lines == 2
rm('test12_CpG.bedGraph')

# No perfectly converted reads
rm('test13_CpG.bedGraph')
check_call([MPath, 'extract', '-o', 'test13', '-q', '5', '--minConversionEfficiency', '1.0', 'chgchh.fa', 'chgchh_aln.bam'])
assert op.exists('test13_CpG.bedGraph')
lines = sum(1 for _ in open('test13_CpG.bedGraph'))
assert lines == 1
rm('test13_CpG.bedGraph')

# Test ignoreNH
rm('test14_CpG.bedGraph')
check_call([MPath, 'extract', '-o', 'test14', '-q', '1', 'cg100.fa', 'NH.bam'])
assert op.exists('test14_CpG.bedGraph')
lines = sum(1 for _ in open('test14_CpG.bedGraph'))
assert lines == 1
rm('test14_CpG.bedGraph')

# Test ignoreNH
rm('test15_CpG.bedGraph')
check_call([MPath, 'extract', '-o', 'test15', '--ignoreNH', '-q', '1', 'cg100.fa', 'NH.bam'])
assert op.exists('test15_CpG.bedGraph')
lines = sum(1 for _ in open('test15_CpG.bedGraph'))
assert lines == 49
rm('test15_CpG.bedGraph')

def bedGraphPositions(fname):
    # Returns the set of 0-based start positions from a bedGraph, skipping the
    # "track ..." header line.
    positions = set()
    for line in open(fname):
        if line.startswith('track'):
            continue
        positions.add(int(line.split()[1]))
    return positions

# biopos_aln.bam holds two synthetic single-end reads mapped to cg100.fa,
# both matching the reference exactly (fully methylated) over their full
# 100bp: "readOT" (flag 0, OT/forward) calls the C side of each CpG at even
# positions 0,2,...,96; "readOB" (flag 16, OB/reverse-complemented in the
# BAM) calls the G side at odd positions 1,3,...,97. Because they're
# single-end, unpaired, and non-overlapping calls, every surviving call is
# reported on its own line, making --five-prime-trim/--three-prime-trim/
# --max-length masking directly verifiable by position.

# Sanity check the fixture itself: with no masking, every position from both
# reads should be called.
rm('test16_CpG.bedGraph')
check_call([MPath, 'extract', 'cg100.fa', 'biopos_aln.bam', '-o', 'test16'])
assert op.exists('test16_CpG.bedGraph')
expected = set(range(0, 98, 2)) | set(range(1, 98, 2))
assert bedGraphPositions('test16_CpG.bedGraph') == expected
rm('test16_CpG.bedGraph')

# --five-prime-trim 10 --three-prime-trim 5 masks 10 bases from the
# biological 5' end and 5 from the biological 3' end of every read. For the
# OT read (not reverse-complemented in the BAM) that's a literal [0,10) /
# [95,100) mask, retaining indices 10-94: even positions 10,12,...,94. For
# the OB read (reverse-complemented in the BAM, so biological 5'/3' are
# swapped in BAM-index order) the mask must flip sides: 5 bases masked at
# BAM-index 0 (the biological 3' end) and 10 bases masked at BAM-index 99
# (the biological 5' end), retaining indices 5-89: odd positions 5,7,...,89.
# If the direction failed to flip for OB (the exact bug --nOB couldn't
# avoid), this would instead retain indices 10-94 and produce odd positions
# 11,13,...,93 -- a different set, so this test directly catches that
# failure mode.
rm('test17_CpG.bedGraph')
check_call([MPath, 'extract', '--five-prime-trim', '10', '--three-prime-trim', '5', 'cg100.fa', 'biopos_aln.bam', '-o', 'test17'])
assert op.exists('test17_CpG.bedGraph')
expected = set(range(10, 95, 2)) | set(range(5, 90, 2))
assert bedGraphPositions('test17_CpG.bedGraph') == expected
rm('test17_CpG.bedGraph')

# --max-length 50 (with trims at 0) isolates the length-aware absolute-position
# cutoff: only the first 50 biological bases of each read are retained. For OT
# that's BAM indices [0,50); for OB (flipped) it's BAM indices [50,100).
rm('test18_CpG.bedGraph')
check_call([MPath, 'extract', '--max-length', '50', 'cg100.fa', 'biopos_aln.bam', '-o', 'test18'])
assert op.exists('test18_CpG.bedGraph')
expected = set(range(0, 50, 2)) | set(range(51, 98, 2))
assert bedGraphPositions('test18_CpG.bedGraph') == expected
rm('test18_CpG.bedGraph')

# --five-prime-trim/--three-prime-trim/--max-length trim relative to a read's
# biological 5'/3' orientation, while --OT/--OB/--CTOT/--CTOB/--nOT/--nOB/
# --nCTOT/--nCTOB trim relative to BAM storage order/strand. These are two
# incompatible schools of thought and must be rejected as mutually exclusive.
rm('test19_CpG.bedGraph')
try:
    check_call([MPath, 'extract', '--five-prime-trim', '10', '--nOT', '5,0,0,0', 'cg100.fa', 'biopos_aln.bam', '-o', 'test19'])
    raise AssertionError("extract should have failed: --five-prime-trim and --nOT are mutually exclusive")
except CalledProcessError:
    pass
assert not op.exists('test19_CpG.bedGraph')

try:
    check_call([MPath, 'extract', '--max-length', '50', '--OT', '5,0,0,0', 'cg100.fa', 'biopos_aln.bam', '-o', 'test19'])
    raise AssertionError("extract should have failed: --max-length and --OT are mutually exclusive")
except CalledProcessError:
    pass
assert not op.exists('test19_CpG.bedGraph')

try:
    check_call([MPath, 'mbias', '--max-length', '50', '--nOB', '5,0,0,0', 'cg100.fa', 'biopos_aln.bam', 'test19_mbias', '--noSVG'])
    raise AssertionError("mbias should have failed: --max-length and --nOB are mutually exclusive")
except CalledProcessError:
    pass

print("Finished correctly")

