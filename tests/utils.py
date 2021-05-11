import pysam
import filecmp


def assert_bam_files_equal(file1, file2, order_matters=False, compare_header=False):
    """Assert that the contents of the two given bam files are equivalent."""

    if order_matters and compare_header:
        assert filecmp.cmp(file1, file2)
        return

    with pysam.AlignmentFile(file1, "rb", check_sq=False, require_index=False) as bam1, \
        pysam.AlignmentFile(file2, "rb", check_sq=False, require_index=False) as bam2:

        if compare_header:
            assert bam1.header == bam2.header

        if order_matters:
            for read1, read2 in zip(bam1, bam2):
                assert read1 == read2

        else:
            f1_reads = set()

            for r in bam1:
                f1_reads.add(r)

            nreads_2 = 0
            for r in bam2:
                assert r in f1_reads
                nreads_2 += 1

            assert nreads_2 == len(f1_reads)


