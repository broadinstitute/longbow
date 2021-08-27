import logging
import sys
import os

import time

import click
import click_log
import tqdm

import pysam

import numpy as np
import matplotlib.pyplot as plt

from ..utils import bam_utils
from ..utils import plot_utils
from ..utils.model import LibraryModel

from ..segment import command as segment

from ..annotate.command import get_segments


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
    default="mas15",
    show_default=True,
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
    help="Do splitting of reads based on splitter delimiters, rather than whole array structure. "
    "This splitting will cause delimiter sequences to be repeated in each read they bound.",
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

    # Get our model:
    if LibraryModel.has_prebuilt_model(model):
        logger.info(f"Using %s", LibraryModel.pre_configured_models[model]["description"])
        model = LibraryModel.build_pre_configured_model(model)
    else:
        logger.info(f"Loading model from json file: %s", model)
        model = LibraryModel.from_json_file(model)

    if do_simple_splitting:
        logger.info("Splitting algorithm: Simple Splitting")
    else:
        logger.info("Splitting algorithm: Bounded Region")

    # Prepare our delimiters for segmentation below:
    delimiters = segment.create_simple_delimiters(model)

    # Flush the loggers before the pbar is created.
    # NOTE: this is considered bad form (to access the logger handlers directly)
    for h in logger.handlers:
        h.flush()

    print(delimiters)

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

        # Get our adapter names:
        mas_adapter_names = [array_element_adapters[0] for array_element_adapters in model.array_element_structure]
        # Manually add in the last adapter which is the array end marker:
        mas_adapter_names.append(model.array_element_structure[-1][-1])
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
                             f"No {bam_utils.SEGMENTS_TAG} tag detected on read {read.query_name} !")
                sys.exit(1)

            # Get the list of MAS-seq adapters in the segments and adjust for missing first adapters:
            read_mas_seq_adapters = [s.name for s in segments if s.name in mas_adapter_name_set]
            if segments[0].name not in mas_adapter_name_set and \
                    segments[0].name == "10x_Adapter" and \
                    len(read_mas_seq_adapters) > 1 and \
                    read_mas_seq_adapters[0] == model.array_element_structure[1][0]:
                read_mas_seq_adapters.insert(0, model.array_element_structure[0][0])
            # TODO: Change the "10x_Adapter" above to be dependent on the model itself.

            # Segment the array into segments using our actual segmentation algorithm so we have accurate counts:
            if do_simple_splitting:
                segment_tuples = segment.segment_read_with_simple_splitting(read, delimiters, segments)
                array_len = len(segment_tuples)
                print(f"{read.query_name}\t{array_len}", end="")
            else:
                found_tuple, _ = segment.segment_read_with_bounded_region_algorithm(read, model, segments)
                array_len = sum(found_tuple)

            # Increment our start and end adapter counts:
            if array_len > 0:
                array_start_adapter_counts[read_mas_seq_adapters[0]] += 1
                array_end_adapter_counts[read_mas_seq_adapters[-1]] += 1

            # Add our length to our list of lengths:
            array_lengths.append(array_len)

            # Adjust names for segmented direction:
            if read.get_tag(bam_utils.SEGMENTS_RC_TAG) == 1:
                read_mas_seq_adapters = [a + rc_decorator for a in read_mas_seq_adapters]

            # Add our tally for our heatmap:
            if array_len > 1:
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

            print(f"\t{ligation_profile_string}")
            print(f"\t{segment_tuples}")

            # Update progress:
            pbar.update(1)

    logger.info("Processing statistics...")
    array_lengths = np.array(array_lengths)

    # Write our stats out to the appropriate files:
    _write_stats(array_lengths, ligation_profile_count_dict, model.name, output_prefix)

    logger.info("Writing complete ligation matrix...")
    _create_ligation_heatmap(output_prefix, ligation_heat_matrix, index_map, f"MAS-seq Ligations\n({model.name})")

    logger.info("Writing reduced ligation matrix...")
    _create_ligation_heatmap_reduced(output_prefix, ligation_heat_matrix, index_map, f"MAS-seq Ligations\n({model.name})")

    logger.info(f"Done. Elapsed time: %2.2fs.", time.time() - t_start)


def _write_stats(array_lengths, ligation_profile_count_dict, model_name, output_prefix):
    """Write out all basic statistics for the data in the input file."""

    # Calculate histogram of array lengths.  Bins are created around integer values.
    hist_bins = np.arange(int(np.max(array_lengths))+2)
    count_hist, _ = np.histogram(array_lengths, bins=hist_bins-.5)
    hist_bins = hist_bins[:-1]

    # Calculate summary stats:
    array_length_min = np.min(array_lengths)
    array_length_max = np.max(array_lengths)
    array_length_mean = np.mean(array_lengths)
    array_length_median = np.median(array_lengths)
    array_length_std = np.std(array_lengths)

    logger.info("Writing summary stats file...")
    _write_summary_stats_file(output_prefix,
                              array_lengths,
                              array_length_min,
                              array_length_max,
                              array_length_mean,
                              array_length_median,
                              array_length_std,
                              count_hist,
                              hist_bins,
                              ligation_profile_count_dict)

    logger.info("Writing read stat histograms...")
    _create_array_length_histogram(output_prefix,
                                   len(array_lengths),
                                   count_hist,
                                   hist_bins,
                                   array_length_mean,
                                   array_length_median,
                                   model_name)


def _write_summary_stats_file(output_prefix,
                              array_lengths,
                              array_length_min,
                              array_length_max,
                              array_length_mean,
                              array_length_median,
                              array_length_std,
                              count_hist,
                              hist_bins,
                              ligation_profile_count_dict,
                              num_ligation_profiles_to_show=20):

    """Write summary statistics for the given input bam to a file."""

    logger.debug("Calculating ligation profiles...")
    ligation_profile_data = _calculate_top_ligation_profiles(len(array_lengths),
                                                             ligation_profile_count_dict,
                                                             num_ligation_profiles_to_show)

    with open(output_prefix + "_summary_stats.txt", 'w') as f:
        f.write(f"Total Num Reads (Arrays):\t{len(array_lengths)}\n")
        f.write(f"Total Num Array Elements (Segmented Arrays):\tTBD\n")
        f.write(f"Num unique ligation profiles: {len(ligation_profile_count_dict)}\n")
        f.write("\n")

        f.write(f"Array Length Stats:\n")
        f.write(f"min:\t{array_length_min}\n")
        f.write(f"max:\t{array_length_max}\n")
        f.write(f"mean:\t{array_length_mean}\n")
        f.write(f"median:\t{array_length_median}\n")
        f.write(f"std:\t{array_length_std}\n")

        f.write("\n")
        f.write(f"Array Length Hist:\n")
        for i, h in enumerate(count_hist):
            f.write(f"{hist_bins[i]:2d}:\t{h}\n")

        f.write("\n")
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
                                   num_reads,
                                   count_hist,
                                   hist_bins,
                                   array_length_mean,
                                   array_length_median,
                                   model_name):
    """Create a histogram displaying the length distribution of all arrays."""

    handles = []
    fig = plt.figure(figsize=plot_utils.gFIG_SIZE_in)
    ax = fig.add_subplot()

    # Plot the bar graph:
    h = ax.bar(hist_bins, count_hist, label="Array Lengths")
    handles.append(h)

    # Plot the mean line:
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    yrange = ymax - ymin

    h, = ax.plot(
        [array_length_mean, array_length_mean],
        [ymin - (yrange * .1), ymax + (yrange * .2)],
        "--",
        color=[.4] * 3,
        label=f"Mean Array Length = {array_length_mean:.02f}"
    )
    handles.append(h)

    h, = ax.plot(
        [array_length_median, array_length_median],
        [ymin - (yrange * .1), ymax + (yrange * .2)],
        "--",
        color=[0] * 3,
        label=f"Median Array Length = {array_length_median}"
    )
    handles.append(h)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax + (yrange * .1))

    # Add the labels for the bars:
    yoff = (ymax - ymin) * .025
    for c, l in zip(count_hist, hist_bins):
        h = ax.text(l, c + yoff, f"{c}\n{c / num_reads * 100:.02f}%", horizontalalignment='center')

    t = f"MAS-seq Array Length Counts\n({model_name})"
    ax.set_title(f"{output_prefix} {t}")
    ax.set_xlabel("Array Length")
    ax.set_ylabel("Number of Reads")

    ax.set_xticks([x for x in range(int(max(hist_bins)))])

    fig.legend(handles=handles, loc="upper right")
    plot_utils.fix_plot_visuals(fig)

    plot_utils.save_figure(fig, name=t, prefix=output_prefix)


def _create_ligation_heatmap(output_prefix, heat_matrix, index_map, title):
    """Plot the given heatmap which represents the ligations between different MAS-seq adapters.
    The resulting plot will represent the forward and reverse complemented ligations separately."""

    # Get our colormap:
    heat_cmap = plot_utils.get_zero_white_cmap(base_cmap=plot_utils.get_heat_cmap("jet", False))

    # Make the figure huge so we can have our numbers fit:
    num_digits = int(np.ceil(np.log10(np.max(heat_matrix))))

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

    # Make the figure huge so we can have our numbers fit:
    num_digits = int(np.ceil(np.log10(np.max(heat_matrix))))

    # Heuristic here.  We know that 8 digits fit well in a 2x scaled matrix.
    fig_scale = num_digits - 6

    if fig_scale < 2:
        fig_scale = 2

    logger.debug(f"Reduced Heatmap Scale: Digits: {num_digits} | Scale: {fig_scale}")

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
