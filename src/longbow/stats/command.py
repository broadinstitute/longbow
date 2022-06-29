import logging
import re
import sys
import os

import time
import datetime

import click
import click_log
import tqdm

import pysam

import numpy as np
import matplotlib.pyplot as plt

import longbow.utils.constants
from ..utils import bam_utils
from ..utils import plot_utils
from ..utils import model as LongbowModel
from ..utils.model import LibraryModel

from ..segment import command as segment

from ..annotate.command import get_segments


plot_title_path_regex = re.compile(r".*/([^/].*?)/*$")


logging.basicConfig(stream=sys.stderr)
logger = logging.getLogger("stats")
click_log.basic_config(logger)


@click.command(name=logger.name)
@click_log.simple_verbosity_option(logger)
@click.option(
    "-o",
    "--output-prefix",
    default="longbow_stats",
    type=str,
    help="prefix to give to output files",
)
@click.option(
    "-p",
    "--pbi",
    required=False,
    type=click.Path(),
    help="BAM .pbi index file",
)
@click.option(
    "-m",
    "--model",
    help="The model to use for annotation.  If the given value is a pre-configured model name, then that "
         "model will be used.  Otherwise, the given value will be treated as a file name and Longbow will attempt to "
         "read in the file and create a LibraryModel from it.  Longbow will assume the contents are the configuration "
         "of a LibraryModel as per LibraryModel.to_json()."
)
@click.option(
    "-s",
    "--do-simple-splitting",
    required=False,
    is_flag=True,
    default=False,
    help="DEPRECATED.  Do splitting of reads based on splitter delimiters, rather than whole array structure. "
    "This splitting will cause delimiter sequences to be repeated in each read they bound.  "
    "This is now the default setting, and this flag has been DEPRECATED.",
)
@click.argument("input-bam", default="-" if not sys.stdin.isatty() else None, type=click.File("rb"))
def main(pbi, output_prefix, model, do_simple_splitting, input_bam):
    """Calculate and produce stats on the given input bam file."""

    t_start = time.time()

    logger.info("Invoked via: longbow %s", " ".join(sys.argv[1:]))

    # Try to read in a pac bio index file so we can have a good progress bar:
    pbi = f"{input_bam.name}.pbi" if pbi is None else pbi
    read_count = None
    if os.path.exists(pbi):
        read_count = bam_utils.load_read_count(pbi)
        logger.info("Annotating %d reads", read_count)

    if do_simple_splitting:
        logger.warning("Simple splitting is now the default.  \"-s\" / \"--do-simple-splitting\" is now DEPRECATED.")
    do_simple_splitting = True
    logger.info("Using simple splitting mode.")

    # Flush the loggers before the pbar is created.
    # NOTE: this is considered bad form (to access the logger handlers directly)
    for h in logger.handlers:
        h.flush()

    # Track the lengths of our array elements:
    array_element_lengths = []

    pysam.set_verbosity(0)  # silence message about the .bai file not being found
    with pysam.AlignmentFile(
        input_bam, "rb", check_sq=False, require_index=False
    ) as bam_file, tqdm.tqdm(
        desc="Progress",
        unit=" read",
        colour="green",
        total=read_count,
        file=sys.stderr,
        leave=False,
        disable=not sys.stdin.isatty(),
    ) as pbar:

        # Get our model:
        if model is None:
            lb_model = LibraryModel.from_json_obj(bam_utils.get_model_from_bam_header(bam_file.header))
        elif model is not None and LibraryModel.has_prebuilt_model(model):
            lb_model = LibraryModel.build_pre_configured_model(model)
        else:
            lb_model = LibraryModel.from_json_file(model)

        logger.info(f"Using %s: %s", lb_model.name, lb_model.description)

        # Prepare our delimiters for segmentation below:
        delimiters = segment.create_simple_delimiters(lb_model)
        logger.debug(f"Splitting delimiters: %s", str(delimiters))

        # Get our adapter names:
        mas_adapter_names = [array_element_adapters[0] for array_element_adapters in lb_model.array_element_structure]
        # Manually add in the last adapter which is the array end marker:
        mas_adapter_names.append(lb_model.array_element_structure[-1][-1])
        mas_adapter_name_set = set(mas_adapter_names)

        # Create storage point for heatmap data:
        ligation_heat_matrix = np.zeros((len(mas_adapter_names) * 2, len(mas_adapter_names) * 2), dtype=int)
        index_map = dict()
        for i, name in enumerate(mas_adapter_names):
            index_map[name] = i
        for i, name in enumerate(mas_adapter_names):
            index_map[name + "'"] = i + len(mas_adapter_names)

        # Keep track of how many arrays start and end with each element:
        array_start_adapter_counts = {a: 0 for a in mas_adapter_names}
        array_end_adapter_counts = {a: 0 for a in mas_adapter_names}
        array_lengths = []
        ligation_profile_count_dict = dict()

        rc_decorator = "'"

        for read in bam_file:
            # Get our read segments:
            try:
                _, segments = get_segments(read)
            except KeyError:
                logger.error(f"Input bam file does not contain longbow segmented reads!  "
                             f"No {longbow.utils.constants.SEGMENTS_TAG} tag detected on read {read.query_name} !")
                sys.exit(1)

            # Get the list of MAS-seq adapters in the segments and adjust for missing first adapters:
            read_mas_seq_adapters = [s.name for s in segments if s.name in mas_adapter_name_set]
            if segments[0].name not in mas_adapter_name_set and \
                    segments[0].name == lb_model.array_element_structure[0][1] and \
                    len(read_mas_seq_adapters) > 1 and \
                    read_mas_seq_adapters[0] == lb_model.array_element_structure[1][0]:
                read_mas_seq_adapters.insert(0, lb_model.array_element_structure[0][0])
            # NOTE: here model.array_element_structure[0][1] corresponds to "5p_Adapter"
            #       that we use as our second segment in each array element.

            # Segment the array into segments using our actual segmentation algorithm so we have accurate counts:
            if do_simple_splitting:
                segment_tuples = segment.segment_read_with_simple_splitting(read, delimiters, segments)
                array_len = len(segment_tuples)
                logger.debug(f"Split for {read.query_name}: {array_len} segments:")
                # Track segment lengths:
                for seg_tup in segment_tuples:
                    array_element_lengths.append(seg_tup[3] - seg_tup[2] + 1)
            else:
                delimiter_found, delimiter_segments = segment.segment_read_with_bounded_region_algorithm(read, lb_model, segments)
                array_len = sum(delimiter_found)
                # Track segment lengths:
                segment_tuples = []
                for i, seg_list in enumerate(delimiter_segments):
                    if delimiter_found[i]:
                        array_element_lengths.append(seg_list[-1].end - seg_list[0].start + 1)
                        segment_tuples.append(
                            tuple([seg_list[0].name, seg_list[-1].name, seg_list[0].start, seg_list[-1].end])
                        )

            # Here we need to adjust for our segmentation array count.
            # without this step, we will allow array elements consisting of only `random` sections
            # to count as length 1 arrays:
            if len(read_mas_seq_adapters) == 0:
                array_len = 0

            # Increment our start and end adapter counts:
            if len(read_mas_seq_adapters) > 0:
                array_start_adapter_counts[read_mas_seq_adapters[0]] += 1
                array_end_adapter_counts[read_mas_seq_adapters[-1]] += 1

            # Add our length to our list of lengths:
            array_lengths.append(array_len)

            # Adjust names for segmented direction:
            if read.get_tag(longbow.utils.constants.SEGMENTS_RC_TAG) == 1:
                read_mas_seq_adapters = [a + rc_decorator for a in read_mas_seq_adapters]

            # Add our tally for our heatmap:
            if len(read_mas_seq_adapters) > 1:
                cur_adapter = read_mas_seq_adapters[0]
                for next_adapter in read_mas_seq_adapters[1:]:
                    ligation_heat_matrix[index_map[cur_adapter]][index_map[next_adapter]] += 1
                    cur_adapter = next_adapter

            # Track our ligation profile:
            if array_len == 0:
                # Create a string that is descriptive for the non-marker bases:
                ligation_profile_string = "EMPTY (" + " ".join([s.name + rc_decorator for s in segments]) + ")"
            else:
                ligation_profile_string = " ".join(read_mas_seq_adapters)
                
            try:
                ligation_profile_count_dict[ligation_profile_string] += 1
            except KeyError:
                ligation_profile_count_dict[ligation_profile_string] = 1

            logger.debug(f"\tProfile: {ligation_profile_string}")
            logger.debug(f"\tSegments: {segment_tuples}")

            # Update progress:
            pbar.update(1)

    logger.info("Processing statistics...")
    array_lengths = np.array(array_lengths)
    array_element_lengths = np.array(array_element_lengths)

    # Write our stats out to the appropriate files:
    _write_stats(input_bam.name, array_lengths, array_element_lengths, ligation_profile_count_dict, ligation_heat_matrix, lb_model, output_prefix, do_simple_splitting)

    logger.info("Writing complete ligation matrix...")
    _create_ligation_heatmap(output_prefix, ligation_heat_matrix, index_map, f"MAS-seq Ligations\n({lb_model.name})")

    logger.info("Writing reduced ligation matrix...")
    _create_ligation_heatmap_reduced(output_prefix, ligation_heat_matrix, index_map, f"MAS-seq Ligations\n({lb_model.name})")

    logger.info(f"Done. Elapsed time: %2.2fs.", time.time() - t_start)


def _write_stats(input_bam, array_lengths, array_element_lengths, ligation_profile_count_dict, ligation_heat_matrix,
                 lb_model, output_prefix, do_simple_splitting):
    """Write out all basic statistics for the data in the input file."""

    # Calculate histogram of array lengths.  Bins are created around integer values.
    array_length_hist_bins = np.arange(int(np.max(array_lengths))+2)
    array_length_count_hist, _ = np.histogram(array_lengths, bins=array_length_hist_bins-.5)
    array_length_hist_bins = array_length_hist_bins[:-1]

    # Calculate some necessary summary stats on array lengths:
    array_length_mean = np.mean(array_lengths)
    array_length_median = np.median(array_lengths)

    # Calculate some necessary summary stats on array element lengths:
    array_element_length_mean = np.mean(array_element_lengths)
    array_element_length_median = np.median(array_element_lengths)

    num_array_element_length_bins = 150
    try:
        array_element_length_bin_max = int(np.max(array_element_lengths))

        # Create a histogram centered around the mass of the array lengths
        # Track outliers as the last bin of the histogram:
        array_element_length_std = np.std(array_element_lengths)
        sfactor = 4
        range_scaler = np.abs(array_element_length_mean - array_element_length_median)
        if range_scaler < array_element_length_mean * 0.05:
            range_scaler = array_element_length_std
        min_bin = array_element_length_mean - (range_scaler * sfactor)
        if min_bin < 0:
            min_bin = 0
        max_bin = array_element_length_mean + (range_scaler * sfactor)
        bin_width = (max_bin - min_bin) / num_array_element_length_bins

        # Make sure we don't lose the left-hand tail:
        min_bin = 0

        array_element_length_hist_bins = np.arange(min_bin, max_bin, bin_width, dtype=int)
        array_element_length_hist_bins[-1] = array_element_length_bin_max + 1
        array_element_length_count_hist, _ = np.histogram(array_element_lengths,
                                                          bins=array_element_length_hist_bins-(bin_width * .5))

        array_element_length_hist_bins = array_element_length_hist_bins[:-1]
        array_element_length_hist_bins = array_element_length_hist_bins.astype(int)
    except ValueError:
        # In this case, we had no array elements!
        # This is probably because the wrong model was used.
        array_element_length_hist_bins = []
        array_element_length_count_hist = []

    logger.info("Writing summary stats file...")
    _write_summary_stats_file(input_bam,
                              lb_model,
                              output_prefix,
                              do_simple_splitting,
                              array_lengths,
                              array_length_count_hist,
                              array_length_hist_bins,
                              array_element_lengths,
                              array_element_length_count_hist,
                              array_element_length_hist_bins,
                              ligation_profile_count_dict,
                              ligation_heat_matrix)

    logger.info("Writing read stat histograms...")
    _create_array_length_histogram(output_prefix,
                                   "Array Length",
                                   len(array_lengths),
                                   array_length_count_hist,
                                   array_length_hist_bins,
                                   array_length_mean,
                                   array_length_median,
                                   lb_model.name)

    logger.info("Writing array element stat histograms...")
    _create_array_element_length_histogram(output_prefix,
                                           "Array Element Length",
                                           array_element_length_count_hist,
                                           array_element_length_hist_bins,
                                           array_element_length_mean,
                                           array_element_length_median,
                                           lb_model.name)


def _write_summary_stats_file(input_bam,
                              lb_model,
                              output_prefix,
                              do_simple_splitting,
                              array_lengths,
                              array_length_count_hist,
                              array_length_hist_bins,
                              array_element_lengths,
                              array_element_length_count_hist,
                              array_element_length_hist_bins,
                              ligation_profile_count_dict,
                              ligation_heat_matrix,
                              num_ligation_profiles_to_show=40):
    """Write summary statistics for the given input bam to a file."""

    # Calculate summary stats on array lengths:
    array_length_min = np.min(array_lengths)
    array_length_max = np.max(array_lengths)
    array_length_mean = np.mean(array_lengths)
    array_length_median = np.median(array_lengths)
    array_length_std = np.std(array_lengths)

    # Calculate summary stats on array element lengths:
    num_array_elements = len(array_element_lengths)
    if num_array_elements == 0:
        array_element_length_min = 0
        array_element_length_max = 0
        array_element_length_mean = 0
        array_element_length_median = 0
        array_element_length_std = 0
    else:
        array_element_length_min = np.min(array_element_lengths)
        array_element_length_max = np.max(array_element_lengths)
        array_element_length_mean = np.mean(array_element_lengths)
        array_element_length_median = np.median(array_element_lengths)
        array_element_length_std = np.std(array_element_lengths)

    logger.debug("Calculating ligation profiles...")
    ligation_profile_data = _calculate_top_ligation_profiles(len(array_lengths),
                                                             ligation_profile_count_dict,
                                                             num_ligation_profiles_to_show)

    # We need to rotate and flip our count matrix:
    ligation_heat_matrix = np.rot90(np.fliplr(ligation_heat_matrix))

    current_timestamp = time.time()
    timezone = datetime.datetime.now(datetime.timezone.utc).astimezone().tzinfo
    with open(output_prefix + "_summary_stats.txt", 'w') as f:
        f.write("#")
        f.write("=" * 80)
        f.write("\n")
        f.write(f"#Time: {datetime.datetime.fromtimestamp(current_timestamp)} {timezone} ")
        f.write(f"({current_timestamp})\n")
        f.write(f"#Input file: {input_bam}\n")
        f.write(f"#Splitting algorithm: ")
        if do_simple_splitting:
            f.write(f"Simple Splitting")
        else:
            f.write(f"Bounded Region")
        f.write("\n")
        f.write("#")
        f.write("=" * 80)
        f.write("\n")
        f.write("\n")

        logger.debug("Number of reads processed: %d", len(array_lengths))
        logger.debug("Number of array elements found: %d", len(array_element_lengths))

        if len(array_element_lengths) == 0:
            logger.warning("No array elements were found.  Either something went horribly wrong in library prep or "
                           "you used the wrong model when running `longbow stats` (and the latter is more likely - "
                           "CHECK YOUR INPUTS!).")
            f.write("\n")
            f.write("#" + ("=" * 80) + "\n")
            f.write("""#                 __        ___    ____  _   _ ___ _   _  ____\n""")
            f.write("""#                 \ \      / / \  |  _ \| \ | |_ _| \ | |/ ___|\n""")
            f.write("""#                  \ \ /\ / / _ \ | |_) |  \| || ||  \| | |  _\n""")
            f.write("""#                   \ V  V / ___ \|  _ <| |\  || || |\  | |_| |\n""")
            f.write("""#                    \_/\_/_/   \_\_| \_\_| \_|___|_| \_|\____|\n""")
            f.write("#\n")
            f.write("WARNING:                  No array elements were found.\n")
            f.write("WARNING:        Either something went horribly wrong in library prep\n")
            f.write("WARNING:                               or\n")
            f.write("WARNING:        You used the wrong model when running `longbow stats`.\n")
            f.write("WARNING:        The latter is more likely. CHECK YOUR INPUTS!!!!!!!!!!\n")
            f.write("#" + ("=" * 80) + "\n")
            f.write("\n")

        if len(array_lengths) == len(array_element_lengths):
            logger.warning("The number of array elements found was the same as the number of array reads processed.  "
                           "Either something went horribly wrong in library prep or you used the wrong model when "
                           "running `longbow stats` (and the latter is more likely - CHECK YOUR INPUTS!).")
            f.write("\n")
            f.write("#" + ("=" * 80) + "\n")
            f.write("""#                 __        ___    ____  _   _ ___ _   _  ____\n""")
            f.write("""#                 \ \      / / \  |  _ \| \ | |_ _| \ | |/ ___|\n""")
            f.write("""#                  \ \ /\ / / _ \ | |_) |  \| || ||  \| | |  _\n""")
            f.write("""#                   \ V  V / ___ \|  _ <| |\  || || |\  | |_| |\n""")
            f.write("""#                    \_/\_/_/   \_\_| \_\_| \_|___|_| \_|\____|\n""")
            f.write("#\n")
            f.write("WARNING:   The number of array elements found was the same as the number of\n")
            f.write("WARNING:                        array reads processed.\n")
            f.write("WARNING:        Either something went horribly wrong in library prep\n")
            f.write("WARNING:                               or\n")
            f.write("WARNING:        You used the wrong model when running `longbow stats`.\n")
            f.write("WARNING:        The latter is more likely. CHECK YOUR INPUTS!!!!!!!!!!\n")
            f.write("#" + ("=" * 80) + "\n")
            f.write("\n")

        f.write("#" + ("-" * 80) + "\n")
        f.write(f"MAS-seq / Longbow Model:\t{lb_model.name}\n")
        f.write(f"Total Num Reads (Arrays):\t{len(array_lengths)}\n")
        f.write(f"Total Num Array Elements (Segmented Arrays):\t{num_array_elements}\n")
        f.write(f"Output yield gain:\t{num_array_elements/len(array_lengths):.2f}x\n")
        f.write(f"Num unique ligation profiles: {len(ligation_profile_count_dict)}\n")
        f.write("\n")
        f.write("#" + ("-" * 80) + "\n")
        f.write(f"Model structure ({lb_model.name}):\n")
        for array_element_structure in lb_model.array_element_structure:
            s = " ".join(array_element_structure)
            f.write(f"\t{s}\n")
        f.write("\n")

        _write_length_stats_to_file(f, "Array Length", array_length_count_hist,
                                    array_length_hist_bins, array_length_min, array_length_max,
                                    array_length_mean, array_length_median, array_length_std)

        f.write("#" + ("-" * 80) + "\n")
        if do_simple_splitting:
            f.write("#REMINDER: Array splitting performed with Simple Splitting algorithm.\n")
        else:
            f.write("#REMINDER: Array splitting performed with Bounded Region algorithm.\n")
        _write_length_stats_to_file(f, "Array Element Length", array_element_length_count_hist,
                                    array_element_length_hist_bins, array_element_length_min, array_element_length_max,
                                    array_element_length_mean, array_element_length_median, array_element_length_std,
                                    do_outliers=True)

        f.write("#" + ("-" * 80) + "\n")
        f.write("Ligation Matrix Statistics:\n")

        total_count = np.sum(ligation_heat_matrix)

        sub_diagonal_count = 0
        for i in range(1, ligation_heat_matrix.shape[0]):
            sub_diagonal_count += ligation_heat_matrix[i, i-1]
        off_sub_diagonal_count = total_count - sub_diagonal_count

        sub_sub_diagonal_count = 0
        for i in range(2, ligation_heat_matrix.shape[0]):
            sub_sub_diagonal_count += ligation_heat_matrix[i, i-2]

        logger.debug("Heat Matrix Counts:")
        logger.debug("Total Count: %d", total_count)
        logger.debug("Sub Diagonal Count: %d", sub_diagonal_count)
        logger.debug("Off Diagonal Count: %d", off_sub_diagonal_count)
        logger.debug("Sub Sub Diagonal Count: %d", sub_sub_diagonal_count)

        f.write(_get_stat_and_percent_string_if_not_zero("Subdiagonal Count Total (correct segments)", sub_diagonal_count, total_count))
        f.write(_get_stat_and_percent_string_if_not_zero("Off-Subdiagonal Count Total (segmentation / ligation errors)", off_sub_diagonal_count, total_count))
        f.write(_get_stat_and_percent_string_if_not_zero("Sub-Subdiagonal Count Total (missed MAS-seq adapters)", sub_sub_diagonal_count, total_count))
        f.write("\n")

        f.write("#" + ("-" * 80) + "\n")
        f.write("Raw Ligation Matrix:\n")
        _write_heat_matrix(f, ligation_heat_matrix)
        f.write("\n")

        f.write("#" + ("-" * 80) + "\n")
        f.write("Reduced Ligation Matrix:\n")

        reduced_heat_matrix, _ = reduce_heatmap(
            ligation_heat_matrix,
            {i: i for i in range(int(ligation_heat_matrix.shape[0]/2))}
        )
        _write_heat_matrix(f, reduced_heat_matrix)
        f.write("\n")

        f.write("#" + ("-" * 80) + "\n")
        f.write(f"Top {len(ligation_profile_data)} Ligation Profiles:\n")

        field_widths = []
        columns = ["Profile", "Count", "Percent of All Ligations"]
        for i in range(len(ligation_profile_data[0])):
            field_widths.append(max([len(str(p[i])) for p in ligation_profile_data]))
            if field_widths[-1] < len(columns[i]):
                field_widths[-1] = len(columns[i])

        for i, c in enumerate(columns):
            f.write(f"{c:{field_widths[i]}s}\t")
        f.write("\n")

        for profile in ligation_profile_data:
            for i, c in enumerate(profile):
                f.write(f"{c:{field_widths[i]}}\t")
            f.write(f"\n")


def _write_length_stats_to_file(f, name, count_hist, hist_bins, stat_min, stat_max, stat_mean, stat_median, stat_std, do_outliers=False):
    f.write("#" + ("-" * 80) + "\n")
    f.write(f"{name} Stats:\n")
    f.write(f"min:\t{stat_min}\n")
    f.write(f"max:\t{stat_max}\n")
    f.write(f"mean:\t{stat_mean}\n")
    f.write(f"median:\t{stat_median}\n")
    f.write(f"std:\t{stat_std}\n")
    f.write("\n")
    if len(count_hist) > 0:
        f.write("#" + ("-" * 80) + "\n")
        f.write(f"{name} Hist:\n")
        f.write(f"Length   Count\n")
        fstring = 'd' if (isinstance(hist_bins[0], np.signedinteger) or isinstance(hist_bins[0], np.unsignedinteger)) else 'f'
        for i, h in enumerate(count_hist[:-1]):
            f.write(f"{hist_bins[i]:2{fstring}}:\t{h}\n")
        if do_outliers:
            f.write(f"Outliers in higher bins:\t{count_hist[-1]}\n")
        else:
            f.write(f"{hist_bins[-1]:2{fstring}}:\t{count_hist[-1]}\n")
        f.write("\n")
    f.write("\n")


def _get_stat_and_percent_string_if_not_zero(description, stat, total):
    if total == 0:
        logger.warning(f"WARNING: total is zero for {description}: stat: {stat}, total: {total}")
        return f"{description}: {stat} (0.0%)\n"
    elif stat == 0:
        logger.warning(f"WARNING: stat is zero for {description}: stat: {stat}, total: {total}")
        return f"{description}: {stat} (0.0%)\n"
    else:
        return f"{description}: {stat} ({100 * stat / total:.2f}%)\n"


def _write_heat_matrix(f, heat_matrix):
    """Write the given heat matrix to the given file"""

    max_length = _get_num_digits_for_heat_matrix(heat_matrix)

    for i in range(0, heat_matrix.shape[0]):
        for j in range(0, heat_matrix.shape[1]):
            f.write(f"{heat_matrix[i, j]:{max_length}d} ")
        f.write("\n")
    f.write("\n")


def _calculate_top_ligation_profiles(num_reads, ligation_profile_count_dict, num_ligation_profiles_to_show):

    ligation_profiles, profile_counts = zip(*ligation_profile_count_dict.items())
    ligation_profiles = np.array(ligation_profiles)
    profile_counts = np.array(profile_counts)

    # Show the top molecular arrangements:
    if num_ligation_profiles_to_show > len(profile_counts):
        num_ligation_profiles_to_show = len(profile_counts)
    indices = np.argpartition(profile_counts, -num_ligation_profiles_to_show)[-num_ligation_profiles_to_show:]

    # Sort the indices:
    top_counts = profile_counts[indices]
    top_ligations = ligation_profiles[indices]
    idx = np.argsort(top_counts)
    top_counts = top_counts[idx]
    top_ligations = top_ligations[idx]

    # Package up the data:
    data = [
        [top_ligations[i], top_counts[i], f"{(top_counts[i] / num_reads) * 100:2.02f}%"]
        for i in range(len(top_counts) - 1, -1, -1)
    ]
    return data


def _create_array_length_histogram(output_prefix,
                                   stat_name,
                                   num_data_points,
                                   count_hist,
                                   hist_bins,
                                   stat_mean,
                                   stat_median,
                                   model_name):
    """Create a histogram displaying the length distribution of all arrays."""

    handles = []
    fig = plt.figure(figsize=plot_utils.gFIG_SIZE_in)
    ax = fig.add_subplot()

    # Plot the bar graph:
    if len(hist_bins) > 1:
        tic_width = hist_bins[-1] - hist_bins[-2]
    else:
        tic_width = 1

    bar_width = tic_width*0.8
    h = ax.bar(hist_bins, count_hist, width=bar_width, label=f"{stat_name}s")
    handles.append(h)

    # Plot the mean line:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    yrange = ymax - ymin

    h, = ax.plot(
        [stat_mean, stat_mean],
        [ymin - (yrange * .1), ymax + (yrange * .2)],
        "--",
        color=[.4] * 3,
        label=f"Mean {stat_name} = {stat_mean:.02f}"
    )
    handles.append(h)

    h, = ax.plot(
        [stat_median, stat_median],
        [ymin - (yrange * .1), ymax + (yrange * .2)],
        "--",
        color=[0] * 3,
        label=f"Median {stat_name} = {stat_median}"
    )
    handles.append(h)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax + (yrange * .1))

    # Add the labels for the bars:
    yoff = (ymax - ymin) * .025
    for c, l in zip(count_hist, hist_bins):
        h = ax.text(l, c + yoff, f"{c}\n{c / num_data_points * 100:.02f}%", horizontalalignment='center')

    prefix_title = output_prefix
    if "/" in output_prefix:
        while prefix_title.endswith("/"):
            prefix_title = prefix_title[:-1]
        prefix_title = plot_title_path_regex.match(prefix_title).group(1)

    t = f"MAS-seq {stat_name} Counts\n({model_name})"
    ax.set_title(f"{prefix_title}\n{t}")
    ax.set_xlabel(stat_name)
    ax.set_ylabel("Number of Reads")

    ax.set_xticks([x for x in range(0, 1+int(max(hist_bins)), tic_width)])

    fig.legend(handles=handles, loc="upper right")
    plot_utils.fix_plot_visuals(fig)

    plot_utils.save_figure(fig, name=t, prefix=output_prefix)


def _create_array_element_length_histogram(output_prefix,
                                           stat_name,
                                           count_hist,
                                           hist_bins,
                                           stat_mean,
                                           stat_median,
                                           model_name):
    """Create a histogram displaying the length distribution of all arrays."""

    handles = []
    fig = plt.figure(figsize=plot_utils.gFIG_SIZE_in)
    ax = fig.add_subplot()

    # Plot the bar graph:
    ax.plot(hist_bins, count_hist, "-b", label=f"{stat_name}s")
    h, = ax.plot(hist_bins, count_hist, "ob", label=f"{stat_name}s")
    handles.append(h)

    # Plot the mean line:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    yrange = ymax - ymin

    h, = ax.plot(
        [stat_mean, stat_mean],
        [ymin - (yrange * .1), ymax + (yrange * .2)],
        "--",
        color=[.4] * 3,
        label=f"Mean {stat_name} = {stat_mean:.02f}"
    )
    handles.append(h)

    h, = ax.plot(
        [stat_median, stat_median],
        [ymin - (yrange * .1), ymax + (yrange * .2)],
        "--",
        color=[0] * 3,
        label=f"Median {stat_name} = {stat_median}"
    )
    handles.append(h)

    if len(hist_bins) > 1:
        h, = ax.plot(
            [hist_bins[-2] + (hist_bins[-1] - hist_bins[-2])]*2,
            [ymin - (yrange * .1), ymax + (yrange * .2)],
            "-",
            color="red",
            label="Outlier Break"
        )
        handles.append(h)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax + (yrange * .1))

    prefix_title = output_prefix
    if "/" in output_prefix:
        while prefix_title.endswith("/"):
            prefix_title = prefix_title[:-1]
        prefix_title = plot_title_path_regex.match(prefix_title).group(1)

    t = f"MAS-seq {stat_name} Counts\n({model_name})"
    ax.set_title(f"{prefix_title}\n{t}")
    ax.set_xlabel(stat_name)
    ax.set_ylabel("Number of Reads")

    fig.legend(handles=handles, loc="upper right")
    plot_utils.fix_plot_visuals(fig)

    plot_utils.save_figure(fig, name=t, prefix=output_prefix)


def _create_ligation_heatmap(output_prefix, heat_matrix, index_map, title):
    """Plot the given heatmap which represents the ligations between different MAS-seq adapters.
    The resulting plot will represent the forward and reverse complemented ligations separately."""

    # Get our colormap:
    heat_cmap = plot_utils.get_zero_white_cmap(base_cmap=plot_utils.get_heat_cmap("jet", False))

    # Make the figure huge so we can have our numbers fit:
    num_digits = _get_num_digits_for_heat_matrix(heat_matrix, 0)

    # Heuristic here.  We know that 5 digits fit well in a 3x scaled matrix.
    fig_scale = num_digits - 3

    if fig_scale < 3:
        fig_scale = 3

    logger.debug(f"Heatmap Scale: Digits: {num_digits} | Scale: {fig_scale}")
    fig = plt.figure(figsize=[x * fig_scale for x in plot_utils.gFIG_SIZE_in])
    ax = fig.add_subplot()

    # Make data equivalent to bokeh data with some rotations:
    heat_matrix = np.rot90(np.fliplr(heat_matrix))

    heatmap = ax.imshow(heat_matrix, cmap=heat_cmap)
    plt.colorbar(heatmap)

    ligation_overhang_list = list(index_map.keys())
    ax.set_xticks(np.arange(len(ligation_overhang_list)))
    ax.set_yticks(np.arange(len(ligation_overhang_list)))
    ax.set_xticklabels(ligation_overhang_list)
    ax.set_yticklabels(ligation_overhang_list)
    ax.xaxis.tick_top()
    plt.setp(ax.get_xticklabels(), rotation=45)

    ax.set_title(f"{output_prefix} {title}", pad=20)
    ax.set_xlabel("Adapter 1")
    ax.set_ylabel("Adapter 2")

    plot_utils.fix_plot_visuals(fig)

    # Save the figure without numbers first.
    plot_utils.save_figure(fig, name=title, prefix=output_prefix, suffix="no_numbers")

    # Get the color matrix so we can use it to display the counts
    # in an appropriately readable color:
    color_matrix = heatmap.cmap(heatmap.norm(heatmap.get_array()))

    # Add counts in the heatmap:
    for i in range(len(heat_matrix)):
        for j in range(len(heat_matrix[0])):

            if np.mean(color_matrix[i, j][:3]) > 0.5:
                text_color = [0] * 3
            else:
                text_color = [1] * 3

            ax.text(j, i, heat_matrix[i, j],
                    ha="center", va="center", color=text_color, size="xx-small")

    plot_utils.fix_plot_visuals(fig)

    # Save the figure with numbers as well:
    plot_utils.save_figure(fig, name=title, prefix=output_prefix)


# Plot the heatmap we created above:
def _create_ligation_heatmap_reduced(output_prefix, heat_matrix, index_map, title, count_divisor=None, significant_digits=3):
    """Plot the given heatmap which represents the ligations between different MAS-seq adapters.
    The resulting plot will represent both the forward and reverse complemented ligations together
    in the same size-reduced heatmap and does not distinguish between read directions."""

    # Get our colormap:
    heat_cmap = plot_utils.get_zero_white_cmap(base_cmap=plot_utils.get_heat_cmap("jet", False))

    # This heatmap should be the upper left and lower right quadrants of the normal heatmap.
    # We should produce a warning if there are any non-zero counts in the forward->RC indicating
    # squares.

    num_digits = _get_num_digits_for_heat_matrix(heat_matrix, 0)

    # Heuristic here.  We know that 8 digits fit well in a 2x scaled matrix.
    fig_scale = num_digits - 6

    if fig_scale < 2:
        fig_scale = 2

    logger.debug(f"Reduced Heatmap Scale: Digits: {num_digits} | Scale: {fig_scale}")

    reduced_heat_mat, reduced_index_map = reduce_heatmap(heat_matrix, index_map)

    if count_divisor:
        reduced_heat_mat = plot_utils.signif(reduced_heat_mat / count_divisor, significant_digits)

    fig = plt.figure(figsize=[x * fig_scale for x in plot_utils.gFIG_SIZE_in])
    ax = fig.add_subplot()

    # Make data equivalent to bokeh data with some rotations:
    reduced_heat_mat = np.rot90(np.fliplr(reduced_heat_mat))

    heatmap = ax.imshow(reduced_heat_mat, cmap=heat_cmap)
    plt.colorbar(heatmap)

    ligation_overhang_list = list(reduced_index_map.keys())
    ax.set_xticks(np.arange(len(ligation_overhang_list)))
    ax.set_yticks(np.arange(len(ligation_overhang_list)))
    ax.set_xticklabels(ligation_overhang_list)
    ax.set_yticklabels(ligation_overhang_list)
    ax.xaxis.tick_top()
    plt.setp(ax.get_xticklabels(), rotation=45)

    if count_divisor:
        ax.set_title(f"{output_prefix} {title} (count x {count_divisor})", pad=20)
    else:
        ax.set_title(f"{output_prefix} {title}", pad=20)

    ax.set_xlabel("Adapter 1")
    ax.set_ylabel("Adapter 2")

    plot_utils.fix_plot_visuals(fig)

    # Save the figure without numbers first.
    plot_utils.save_figure(fig, name=title, prefix=output_prefix, suffix="reduced_no_numbers")

    # Get the color matrix so we can use it to display the counts
    # in an appropriately readable color:
    color_matrix = heatmap.cmap(heatmap.norm(heatmap.get_array()))

    # Add counts in the heatmap:
    for i in range(len(reduced_heat_mat)):
        for j in range(len(reduced_heat_mat[0])):
            if np.mean(color_matrix[i, j][:3]) > 0.5:
                text_color = [0] * 3
            else:
                text_color = [1] * 3

            # Old color:
            ax.text(j, i, reduced_heat_mat[i, j],
                    ha="center", va="center", color=text_color, size="medium", weight=1000)

    plot_utils.fix_plot_visuals(fig)

    # Save the figure with numbers as well:
    plot_utils.save_figure(fig, name=title, prefix=output_prefix, suffix="reduced")


def _get_num_digits_for_heat_matrix(heat_matrix, offset=1):
    # Make the figure huge so we can have our numbers fit:
    max_heat_val = np.max(heat_matrix)
    if max_heat_val == 0:
        logger.error("ERROR: Heat matrix is basically all zeros.  This is a problem!")
        num_digits = 2
    else:
        num_digits = int(np.ceil(np.log10(max_heat_val))) + offset
    return num_digits


def reduce_heatmap(heat_matrix, index_map):
    """Reduce the given count heat matrix to be 1/4 the size by summing the upper left and lower right quadrants."""

    # Check upper right and lower left quadrants to see if there are any non-zero
    # entries and warn the user if there are:
    half_point = int(heat_matrix.shape[0] / 2)
    ur_sum = heat_matrix[0:half_point, half_point:].sum()
    ll_sum = heat_matrix[half_point:, 0:half_point].sum()
    if ur_sum != 0:
        logger.warning("WARNING: "
                       "Upper right quadrant of heat matrix has nonzero values that are ignored by this method!")
    if ll_sum != 0:
        logger.warning("WARNING: "
                       "Lower left quadrant of heat matrix has nonzero values that are ignored by this method!")
    # Convert given heat matrix into the reduced form for the forward direction only:
    reduced_heat_mat = np.copy(heat_matrix[0:half_point, 0:half_point])
    reduced_heat_mat += np.copy(heat_matrix[half_point:, half_point:])
    reduced_index_map = {k: v for k, v in list(index_map.items())[0:half_point]}

    return reduced_heat_mat, reduced_index_map
