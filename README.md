# gi2015
Scripts and processed data for reproducing Genome Informatics 2015 talk on junctions found across ~21,500 SRA samples

The presentation itself is in `gi2015.key` and `gi2015.pdf`. The Python script `gi2015.py` generates all the data used in the talk (and much more), but it depends on a list of junctions that's currently unreleased. Its results are contained in the following files whose formats are described below. See `gi2015.py`'s docstring for still more information.

### gi2015.venn.txt
Intersections between junctions obtained from SEQC protocol and junctions
obtained from Rail for the 1720 samples studied by both SEQC and Rail. See
file for details.

### gi2015.sim.txt
Results of running HISAT2 2.0.0-beta and STAR 2.4.2a on the two simulated
human samples from the RGASP spliced alignment paper
"Systematic evaluation of spliced alignment programs for RNA-seq data"
(40 million read pairs each). Junction and overlap accuracy as defined
in the Rail-RNA preprint are given for two protocols: in one, the junctions
that are actually in the sample are provided as annotation, and in the other,
the union of the three gene annotations considered here
(Gencode v19, Ensembl v75, and refGene) are provided.

### gi2015.sex.txt
#### Tab-separated fields, in order of descending field 9

1. project accession number
2. sample accession number
3. experiment accession number
4. run accession number
2. total Y chromosome junction overlaps
3. total junction overlaps in XIST (chrX:73040486-73072588)
4. total overlaps
5. total junctions on chrY
6. total junctions in XIST (chrX:73040486-73072588)
7. total junctions in sample
8. 1 if sample is annotated as male on SRA; 0 if sample is annotated as female
9. field 2 / field 4

### gi2015.[type].tsv, where [type] is in {bottom_10_pct, top_10_pct,
                                    all_words}
Gives common words among samples in bottom and top 10 percent in terms of
proportion of junctions that are annotated -- as well as for all samples.
Only samples with >= 10k junctions are considered.
### Tab-separated fields, in order of descending field 2

1. word
2. number of samples (in which >=10k junctions were found) in which word
    appears

### gi2015.sample_count_submission_date_overlap_geq_20.tsv
#### Tab-separated fields

1. count of samples in which a given junction was found
2. count of projects in which a given junction was found
3. earliest known discovery date (in units of days after February 27, 2009)
    -- this is the earliest known submission date of a sample associated with a
    junction

Above, each junction is covered by at least 20 reads per sample.

### gi2015.[type].stats.tsv, where [type] is in {project, sample}
#### Tab-separated fields

1. [type] count
2. Number of junctions found in >= field 1 [type]s
3. Number of annotated junctions found in >= field 1 [type]s
4. Number of exonskips found in >= field 1 [type]s (exon skip: both 5' and 3'
    splice sites are annotated, but not in the same exon-exon junction)
5. Number of altstartends found in >= field 1 [type]s (altstartend: either 5'
    or 3' splice site is annotated, but not both)
6. Number of novel junctions found in >= field 1 [type]s (novel: both 5' and 
    3' splice sites are unannotated)
7. Number of GT-AG junctions found in >= field 1 [type]s
8. Number of annotated GT-AG junctions found in >= field 1 [type]s
9. Number of GC-AG junctions found in >= field 1 [type]s
10. Number of annotated GC-AG junctions found in >= field 1 [type]s
11. Number of AT-AC junctions found in >= field 1 [type]s
12. Number of annotated AT-AC junctions found in >= field 1 [type]s

### gi2015.stats_by_sample.tsv
#### Tab-separated fields

1. sample index
2. project accession number
3. sample accession number
4. experiment accession number
5. run accession number
6. junction count
7. annotated junction count
8. count of junctions overlapped by at least 5 reads
9. count of annotated junctions overlapped by at least 5 reads
10. total overlap instances
11. total annotated overlap instances

Mathematica 10 was used to make all plots. See the notebook `gi2015.nb`.
