OBJECTS=annotate/mas_10+sc_10x5p.expected.bam annotate/mas_15+sc_10x5p.expected.bam filter/mas_10+sc_10x5p.expected.bam filter/mas_15+sc_10x5p.expected.bam segment/mas_10+sc_10x5p.expected.bam segment/mas_15+sc_10x5p.expected.bam sift/mas_10+sc_10x5p.expected.bam sift/mas_15+sc_10x5p.expected.bam pad/mas_10+sc_10x5p.expected.bam pad/mas_15+sc_10x5p.expected.bam extract/mas_10+sc_10x5p.expected.bam extract/mas_15+sc_10x5p.expected.bam

all: $(OBJECTS)

annotate/mas_%+sc_10x5p.expected.bam: mas%_test_input.bam
	mkdir -p annotate
	longbow annotate -m mas_$*+sc_10x5p -t 1 -f -o $@ $<

filter/mas_%+sc_10x5p.expected.bam: annotate/mas_%+sc_10x5p.expected.bam
	mkdir -p filter
	longbow filter -f -o $@ $<

segment/mas_%+sc_10x5p.expected.bam: filter/mas_%+sc_10x5p.expected.bam
	mkdir -p segment
	longbow segment -t 1 -f -o $@ $<

sift/mas_%+sc_10x5p.expected.bam: segment/mas_%+sc_10x5p.expected.bam
	mkdir -p sift
	longbow sift -f -o $@ $<

pad/mas_%+sc_10x5p.expected.bam: sift/mas_%+sc_10x5p.expected.bam
	mkdir -p pad
	longbow pad -t 1 --barcode-tag ZU -n XM -e 2 -f -o $@ $<

extract/mas_%+sc_10x5p.expected.bam: pad/mas_%+sc_10x5p.expected.bam
	mkdir -p extract
	longbow extract -f -o $@ $<
