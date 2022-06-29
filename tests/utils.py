import sys
import pysam
import filecmp

from pathlib import Path


def convert_sam_to_bam(sam_path, out_bam_path):
    with pysam.AlignmentFile(sam_path, "r", check_sq=False, require_index=False) as input_file:
        with pysam.AlignmentFile(out_bam_path, "wb", header=input_file.header) as out_bam_file:
            for r in input_file:
                out_bam_file.write(r)


def assert_read_tags_are_equal(actual_read, expected_read):

    actual_tag_val_dict = dict()
    actual_tag_type_dict = dict()
    for tag, val, tp in actual_read.get_tags(with_value_type=True):
        actual_tag_val_dict[tag] = val
        actual_tag_type_dict[tag] = tp

    expected_tag_val_dict = dict()
    expected_tag_type_dict = dict()
    for tag, val, tp in expected_read.get_tags(with_value_type=True):
        expected_tag_val_dict[tag] = val
        expected_tag_type_dict[tag] = tp

    assert len(actual_tag_val_dict) == len(expected_tag_val_dict), f"Read {actual_read.query_name}: Number of tags not equal: {len(actual_tag_val_dict)} != {len(expected_tag_val_dict)} || {actual_tag_val_dict} | {expected_tag_val_dict}"

    for tag in actual_tag_val_dict.keys():

        assert tag in expected_tag_val_dict, f"Read {actual_read.query_name}: Tag {tag} not in expected tags!"
        assert tag in actual_tag_type_dict, f"Read {actual_read.query_name}: Tag {tag} not in expected tags!"

        actual_val = actual_tag_val_dict[tag]
        actual_tp = actual_tag_type_dict[tag]

        expected_val = expected_tag_val_dict[tag]
        expected_tp = expected_tag_type_dict[tag]

        assert actual_val == expected_val, f"Read {actual_read.query_name}: Actual and expected tag values are not equal for tag {tag}: {actual_val} != {expected_val}"
        assert actual_tp == expected_tp, f"Read {actual_read.query_name}: Actual and expected tag types are not equal for tag {tag}: {actual_tp} != {expected_tp}"


def assert_reads_are_equal(actual_read, expected_read):

    # Go through most fields until tags:
    assert actual_read.query_name == expected_read.query_name, f"Read {actual_read.query_name}: Read names not equal: {actual_read.query_name} != {expected_read.query_name}"
    assert actual_read.flag == expected_read.flag, f"Read {actual_read.query_name}: Read flags not equal"
    assert actual_read.reference_name == expected_read.reference_name, f"Read {actual_read.query_name}: Contig names not equal:  {actual_read.reference_name} != {expected_read.reference_name}"
    assert actual_read.reference_start == expected_read.reference_start, f"Read {actual_read.query_name}: Start position not equal:  {actual_read.reference_start} != {expected_read.reference_start}"
    assert actual_read.mapping_quality == expected_read.mapping_quality, f"Read {actual_read.query_name}: Mapping qualities not equal:  {actual_read.mapping_quality} != {expected_read.mapping_quality}"
    assert actual_read.cigarstring == expected_read.cigarstring, f"Read {actual_read.query_name}: CIGAR strings not equal:  {actual_read.cigarstring} != {expected_read.cigarstring}"
    assert actual_read.next_reference_id == expected_read.next_reference_id, f"Read {actual_read.query_name}: RNEXT fields not equal:  {actual_read.next_reference_id} != {expected_read.next_reference_id}"
    assert actual_read.next_reference_start == expected_read.next_reference_start, f"Read {actual_read.query_name}: PNEXT fields not equal:  {actual_read.next_reference_start} != {expected_read.next_reference_start}"
    assert actual_read.template_length == expected_read.template_length, f"Read {actual_read.query_name}: TLEN fields not equal:  {actual_read.template_length} != {expected_read.template_length}"
    assert actual_read.query_sequence == expected_read.query_sequence, f"Read {actual_read.query_name}: Base sequences not equal"
    assert actual_read.query_qualities == expected_read.query_qualities, f"Read {actual_read.query_name}: Base qualities not equal"

    # Tags should be viewed as a collection, rather than ordered list.
    # As long as the values between the two reads are all present, we should call the tags equal:
    assert_read_tags_are_equal(actual_read, expected_read)


def assert_reads_files_equal(actual_file, expected_file, order_matters=False, compare_header=False):
    """Assert that the contents of the two given bam files are equivalent."""

    actual_file = Path(actual_file)
    expected_file = Path(expected_file)

    actual_file_flags = "rb" if actual_file.resolve().name.endswith(".bam") else "r"
    expected_file_flags = "rb" if expected_file.resolve().name.endswith(".bam") else "r"

    if order_matters and compare_header:
        assert filecmp.cmp(actual_file, expected_file)
        return

    with pysam.AlignmentFile(actual_file, actual_file_flags, check_sq=False, require_index=False) as actual_bam, \
        pysam.AlignmentFile(expected_file, expected_file_flags, check_sq=False, require_index=False) as expected_bam:

        if compare_header:
            assert actual_bam.header == expected_bam.header, "Headers are not equal."

        if order_matters:
            for read1, read2 in zip(actual_bam, expected_bam):
                assert_reads_are_equal(read1, read2)
        else:
            expected_reads = dict()

            for r in expected_bam:
                expected_reads[r.query_name] = r

            nreads_2 = 0
            for r in actual_bam:
                assert r.query_name in expected_reads, f"Read not in expected file: {r.query_name}"
                assert_reads_are_equal(r, expected_reads[r.query_name])
                nreads_2 += 1

            assert nreads_2 == len(expected_reads)


def assert_text_files_equal(actual_file, expected_file, header_prefix="#", order_matters=False, compare_header=False):
    """Assert that the contents of the two given bam files are equivalent."""

    line_no = 1

    if order_matters:
        if compare_header:
            assert filecmp.cmp(actual_file, expected_file)
            return
        else:
            with open(actual_file, 'r') as f_actual:
                with open(expected_file, 'r') as f_expected:


                    # Skip headers:
                    actual_line = f_actual.readline()
                    line_no += 1
                    expected_line = f_expected.readline()
                    while actual_line.startswith(header_prefix):
                        actual_line = f_actual.readline()
                        line_no += 1
                    while expected_line.startswith(header_prefix):
                        expected_line = f_expected.readline()

                    # Get the rest of the lines:

                    while actual_line and expected_line:
                        line_is_ok = actual_line == expected_line
                        if not line_is_ok:
                            print(f"ERROR: Line {line_no} not equal: Actual: {actual_line} || "
                                  f"Expected: {expected_line}", file=sys.stderr)

                        assert line_is_ok

                        actual_line = f_actual.readline()
                        expected_line = f_expected.readline()

                        line_no += 1

                    assert not actual_line and not expected_line
    else:
        if compare_header:
            print("ERROR: This type of comparison is invalid.")
            assert False
        else:
            with open(actual_file, 'r') as f_actual:
                with open(expected_file, 'r') as f_expected:

                    # Skip headers:
                    actual_line = f_actual.readline()
                    line_no += 1
                    expected_line = f_expected.readline()
                    while actual_line.startswith(header_prefix):
                        actual_line = f_actual.readline()
                        line_no += 1
                    while expected_line.startswith(header_prefix):
                        expected_line = f_expected.readline()

                    expected_lines = f_expected.readlines()

                    for i, line in enumerate(f_actual):
                        line_is_ok = line in expected_lines
                        if not line_is_ok:
                            print(f"ERROR: Line {line_no + i} found in expected file: Actual: {actual_line}",
                                  file=sys.stderr)
                        assert line_is_ok

                    has_lines_remaining = len(expected_lines) != i
                    if has_lines_remaining:
                        print(f"ERROR: number of lines is different in actual and expected files.",
                              file=sys.stderr)
                    assert not has_lines_remaining


def cat_file_to_pipe(filename, proc):
    with open(filename, "rb") as input_file:
        input_bytes = input_file.read()

        proc.stdin.write(input_bytes)
        proc.stdin.flush()

    proc.stdin.close()
    proc.wait()
