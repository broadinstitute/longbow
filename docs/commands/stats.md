---
layout: default
title: stats 
description: "Calculate stats on input file."
nav_order: 9
parent: Commands
---

# Stats

## Description

Stats produces several output files that are useful for quality control.  Examples of files produced are:

- summary statistics
- ligation heatmap plots

All figures are created as both `.png` and `.svg` files.

## Command help

```shell
$ longbow stats --help
Usage: longbow stats [OPTIONS] INPUT_BAM

  Calculate and produce stats on the given input bam file.

Options:
  -v, --verbosity LVL        Either CRITICAL, ERROR, WARNING, INFO or DEBUG
  -o, --output-prefix TEXT   prefix to give to output files
  -p, --pbi PATH             BAM .pbi index file
  -m, --model TEXT           The model to use for annotation.  If the given
                             value is a pre-configured model name, then that
                             model will be used.  Otherwise, the given value
                             will be treated as a file name and Longbow will
                             attempt to read in the file and create a
                             LibraryModel from it.  Longbow will assume the
                             contents are the configuration of a LibraryModel
                             as per LibraryModel.to_json().  [default: mas15]
  -s, --do-simple-splitting  Do splitting of reads based on splitter
                             delimiters, rather than whole array structure.
                             This splitting will cause delimiter sequences to
                             be repeated in each read they bound.
  --help                     Show this message and exit.
```

## Example

```shell
$ longbow stats input.bam
[INFO 2021-08-25 15:01:37    stats] Invoked via: longbow stats input.bam
[INFO 2021-08-25 15:01:37    stats] Using The standard MAS-seq 15 array element model.
Progress: 0 read [00:00, ? read/s]
[INFO 2021-08-25 15:01:40    stats] Processing statistics...
[INFO 2021-08-25 15:01:40    stats] Writing summary stats file...
[INFO 2021-08-25 15:01:40    stats] Writing complete ligation matrix...
[INFO 2021-08-25 15:01:50    stats] Writing reduced ligation matrix...
[INFO 2021-08-25 15:01:53    stats] Done. Elapsed time: 16.36s.
```

```shell
$ ls longbow_stats*
longbow_stats_00_MAS-seq_Array_Length_Counts_mas15.png
longbow_stats_00_MAS-seq_Array_Length_Counts_mas15.svg
longbow_stats_01_MAS-seq_Ligations_mas15_no_numbers.png
longbow_stats_01_MAS-seq_Ligations_mas15_no_numbers.svg
longbow_stats_02_MAS-seq_Ligations_mas15.png
longbow_stats_02_MAS-seq_Ligations_mas15.svg
longbow_stats_03_MAS-seq_Ligations_mas15_reduced_no_numbers.png
longbow_stats_03_MAS-seq_Ligations_mas15_reduced_no_numbers.svg
longbow_stats_04_MAS-seq_Ligations_mas15_reduced.png
longbow_stats_04_MAS-seq_Ligations_mas15_reduced.svg
longbow_stats_summary_stats.txt
```

```shell
$ cat longbow_stats_summary_stats.txt
#================================================================================
#Time: 2021-09-10 15:59:37.763512 PDT (1631314777.7635121)
#Input file: input.bam
#================================================================================

#--------------------------------------------------------------------------------
Total Num Reads (Arrays):	109
Total Num Array Elements (Segmented Arrays):	1635
Output yield gain:	15.00x
Num unique ligation profiles: 2

#--------------------------------------------------------------------------------
Array Length Stats:
min:	15
max:	15
mean:	15.0
median:	15.0
std:	0.0

#--------------------------------------------------------------------------------
Array Length Hist:
Length   Count
 0:	0
 1:	0
 2:	0
 3:	0
 4:	0
 5:	0
 6:	0
 7:	0
 8:	0
 9:	0
10:	0
11:	0
12:	0
13:	0
14:	0
15:	109

#--------------------------------------------------------------------------------
Ligation Matrix Statistics:
Subdiagonal Count Total (correct segments): 1635
Off-Subdiagonal Count Total (segmentation / ligation errors): 0
Sub-Subdiagonal Count Total (missed MAS-seq adapters): 0
Correctly placed segments percentage: 100.00%

#--------------------------------------------------------------------------------
Raw Ligation Matrix:
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
62  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
 0 62  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
 0  0 62  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
 0  0  0 62  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
 0  0  0  0 62  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
 0  0  0  0  0 62  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
 0  0  0  0  0  0 62  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
 0  0  0  0  0  0  0 62  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
 0  0  0  0  0  0  0  0 62  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
 0  0  0  0  0  0  0  0  0 62  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
 0  0  0  0  0  0  0  0  0  0 62  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
 0  0  0  0  0  0  0  0  0  0  0 62  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
 0  0  0  0  0  0  0  0  0  0  0  0 62  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
 0  0  0  0  0  0  0  0  0  0  0  0  0 62  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
 0  0  0  0  0  0  0  0  0  0  0  0  0  0 62  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 47  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 47  0  0  0  0  0  0  0  0  0  0  0  0  0  0 
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 47  0  0  0  0  0  0  0  0  0  0  0  0  0 
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 47  0  0  0  0  0  0  0  0  0  0  0  0 
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 47  0  0  0  0  0  0  0  0  0  0  0 
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 47  0  0  0  0  0  0  0  0  0  0 
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 47  0  0  0  0  0  0  0  0  0 
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 47  0  0  0  0  0  0  0  0 
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 47  0  0  0  0  0  0  0 
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 47  0  0  0  0  0  0 
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 47  0  0  0  0  0 
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 47  0  0  0  0 
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 47  0  0  0 
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 47  0  0 
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 47  0 


#--------------------------------------------------------------------------------
Reduced Ligation Matrix:
  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0 
109   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0 
  0 109   0   0   0   0   0   0   0   0   0   0   0   0   0   0 
  0   0 109   0   0   0   0   0   0   0   0   0   0   0   0   0 
  0   0   0 109   0   0   0   0   0   0   0   0   0   0   0   0 
  0   0   0   0 109   0   0   0   0   0   0   0   0   0   0   0 
  0   0   0   0   0 109   0   0   0   0   0   0   0   0   0   0 
  0   0   0   0   0   0 109   0   0   0   0   0   0   0   0   0 
  0   0   0   0   0   0   0 109   0   0   0   0   0   0   0   0 
  0   0   0   0   0   0   0   0 109   0   0   0   0   0   0   0 
  0   0   0   0   0   0   0   0   0 109   0   0   0   0   0   0 
  0   0   0   0   0   0   0   0   0   0 109   0   0   0   0   0 
  0   0   0   0   0   0   0   0   0   0   0 109   0   0   0   0 
  0   0   0   0   0   0   0   0   0   0   0   0 109   0   0   0 
  0   0   0   0   0   0   0   0   0   0   0   0   0 109   0   0 
  0   0   0   0   0   0   0   0   0   0   0   0   0   0 109   0 


#--------------------------------------------------------------------------------
Top 2 Ligation Profiles:
Profile                                        	Count	Percent of All Ligations	
A B C D E F G H I J K L M N O P                	   62	56.88%                  	
A' B' C' D' E' F' G' H' I' J' K' L' M' N' O' P'	   47	43.12%                  	

```
