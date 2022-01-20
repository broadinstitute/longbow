---
layout: default
title: FAQ
nav_order: 5
description: "Longbow frequently asked questions."
---

# FAQ

## Read Tags

Longbow adds several tags to output reads to indicate various metadata features.  The following table describes 
the read tags added by each tool and what data they contain:

| Read Tag | Type | Tool(s) | Description |
|:---|:---|:---|:---|
| YS | f | Annotate, Demultiplex | Longbow HMM model score (log probability) |
| SG | Z | Annotate, Demultiplex | Read segment information indicating the boundaries and labels of each segment in the read.   For example:   ```  SG:Z:random:0-53,Poly_T:54-85,cDNA:86-823,MARS:824-853,N:854-869,VENUS:870-894,CBC:895-910,UMI:911-920,Poly_T:921-951,cDNA:952-2030,MARS:2031-2060,O:2061-2076,VENUS:2077-2100,CBC:2101-2116,UMI:2117-2126,Poly_T:2127-2157,cDNA:2158-3747,MARS:3748-3777,P:3778-3794 ``` |
| YN | Z | Annotate, Demultiplex | Name of the model used to annotate the reads |
| RC | i | Annotate, Demultiplex | 1 if the read has been reverse-complemented; 0 otherwise |
| XQ | Z | Annotate, Demultiplex | Comma-delimited quality scores ratios for each segment in a given read.  Scores for random segments are set to 0/0.   Scores for other segments are the optimal alignment score (Smith-Waterman) for the expected sequence over 2x the sequence length.   For example:   ```XQ:Z:0/0,12/60,0/0,58/60,30/32,46/46,0/0,0/0,60/60,0/0,58/60,30/32,46/46,0/0,0/0,60/60,0/0,58/60,32/32``` |
| YQ | f | Annotate, Demultiplex | Approximate read quality.  Read quality is approximated by summing the individual segment scores and dividing that sum by the sum of the total possible (best) scores (as defined by the `XQ` tag). |
| ZS | i | Segment, Extract | 1 if the given read is segmented; 0 otherwise |
| YV | i | Filter | 1 if the given read is valid according to the expected order of the model segments; 0 otherwise |
| YK | Z | Filter | Name of the first valid adapter / named non-random region in the given read |
| YG | i | Filter | Number of valid adapters / named non-random regions in the given read |
| ZU | Z | Segment | Unique Molecular Identifier (UMI) sequence for the given segmented read |
| XU | i | Segment | Position (0-based) in the read of the first base in the annotated Unique Molecular Identifier (UMI) sequence for the given segmented read |
| CR | Z | Segment | Cell Barcode (CBC) sequence for the given segmented read |
| XB | i | Segment | Position (0-based) in the read of the first base in the Cell Barcode (CBC) sequence in the given segmented read |
| XF | f | Segment | Confidence factor for the Cell Barcode (CBC) in the given segmented read.  The confidence factor is based on the quality of each base in the CBC and is given by the following:  ``` scale_factor = 100 scale_factor * reduce(operator.mul, map(lambda q: 1. - 10 ** (-(ord(q) - 33.) / 10), qual_string))  ``` |
| XA | Z | Segment | String containing the name of the tag containing the Cell Barcode (CBC) tag followed by the Unique Molecular Identifier (UMI) tag:  `XC-XM` |
| X1 | Z | Segment | Sequence for Spatial Barcode 1 for the given segmented read |
| XP | i | Segment | Position (0-based) in the read of the first base in the Spatial Barcode 1 sequence in the given segmented read |
| X2 | Z | Segment | Sequence for Spatial Barcode 2 for the given segmented read |
| XQ | i | Segment | Position (0-based) in the read of the first base in the Spatial Barcode 2 sequence in the given segmented read |
| XM | Z | Segment | Raw Unique Molecular Identifier (UMI) sequence for the given segmented read (for IsoSeq3 compatibility) |
| XC | Z | Segment | Raw Cell Barcode (CBC) sequence for the given segmented read (for IsoSeq3 compatibility) |
| ic | i | Segment | Sum of number of passes from all ZMWs used to create consensus.  Always set to 1.  (for IsoSeq3 compatibility) |
| im | Z | Segment | ZMW names associated with a given segmented read.  Set to the name of the parent read for a segmented read. (e.g. `m64013e_211031_055434/1/ccs`).  (for IsoSeq3 compatibility) |
| is | i | Segment | Number of ZMWs associated with a given segmented read.  Always set to 1.  (for IsoSeq3 compatibility) |
| it | Z | Segment | List of barcodes / UMIs tagged/clipped during segmentation (e.g. `it:Z:CATTAGGTCATCCCTA,AAATTTTGGA`) (for IsoSeq3 compatibility) |
| zm | i | Segment | ZMW number from which the given read originates.  (for IsoSeq3 compatibility) |
| XN | Z | Segment, Extract | Altered read name given by Longbow to a segmented read (used for debugging).  This name consists of the original read name followed by the start and end positions of the segment on the original read, then the names of the bounding adapters / known regions.  For example:  `XN:Z:m64013e_211029_235558/24/ccs/0_1960/START-MARS` |
### Read Tag Types
The SAM spec defines the following as types for read tags:

| Type | Description |
|:------|:-------------|
| A | printable char |
| i | signed int |
| f | float |
| Z | printable string |
| H | Byte array in hex format |
| B | Integer or numeric array |