#!/usr/bin/env python
"""
gi2015.py
Abhi Nellore / October 23, 2015

Reproduces data used in Mathematica 10 notebook gi2015.nb and Keynote
presentation gi2015.key for Abhi Nellore's Genome Informatics talk on
October 29, 2015.

PREREQUISITES:
1. Create the hg19 FASTA. Download
    http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit, which
    should be converted to hg19.fa with
    http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
    (or another compile)
    commands:
    chmod +x /path/to/twoBitToFa
    /path/to/twoBitToFa hg19.2bit hg19.fa
2. Get HISAT 2.0.0-beta from https://ccb.jhu.edu/software/hisat2/index.shtml .
3. Get STAR 2.4.2a from https://github.com/alexdobin/STAR .
4. Use create_indexes.sh to generate HISAT and STAR indexes for hg19. The
    script must be edited to suit your configuration; look at the comments
    in it for more details.
5. Get SAMTools at http://samtools.sourceforge.net/ . We used SAMTools 1.2 .

File requirements:
1. all_SRA_introns.tsv.gz: database of introns found across ~21,500 SRA samples
    NOT PROVIDED IN THIS REPO BUT WILL BE RELEASED WITH A PREPRINT SOON.
2. index_to_SRA_accession.tsv: maps sample indexes from all_SRA_introns.tsv.gz
    to SRA run accession numbers (regex: [SED]RR\d+) . (In this repo.)
3. gencode.v19.annotation.gtf.gz
    (from ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/)
4. Homo_sapiens.GRCh37.75.gtf.gz
    (from ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/)
5. refGene.gtf
    (downloaded according to instructions at 
     https://groups.google.com/a/soe.ucsc.edu/forum/#!msg/genome/
     bYEoa_hrSiI/cJ8WjnqXhlIJ ; uses command
     mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N \ 
      -e "select * from refGene;" hg19 | cut -f2- | genePredToGtf
      file stdin refGene.gtf
    Since this isn't necessarily frozen, here's a link:
    http://verve.webfactional.com/misc/refGene.gtf
6. The files http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-1728/
                Human_simulated_reads_1_PE1.fq.gz,
             http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-1728/
                Human_simulated_reads_1_PE2.fq.gz,
             http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-1728/
                Human_simulated_reads_2_PE1.fq.gz,
             http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-1728/
                Human_simulated_reads_2_PE2.fq.gz,
            http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-1728/Sim1_truth.bam,
            http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-1728/Sim2_truth.bam
    which are two hg19 simulations from the paper entitled
    "Systematic evaluation of spliced alignment programs for RNA-seq data"
    by Engstrom et al. (RGASP consortium) in Nature Methods. These files are
    from ArrayExpress and have accession code E-MTAB-1728.
7. http://www.nature.com/nbt/journal/v32/n9/extref/nbt.2957-S4.zip, which
    is Supplementary Data 3 from the paper "A comprehensive assessment of
    RNA-seq accuracy, reproducibility and information content by the
    Sequencing Quality Control Consortium" by SEQC/MAQC-III Consortium
    in Nature Biotech. The junctions on this list are used to make a
    Venn Diagram.
8. The file
    http://verve.webfactional.com/misc/all_illumina_sra_for_human.tsv.gz,
    which has metadata grabbed from the SRA.
9. biosample_tags.tsv, which is in this repo and was generated using
    get_biosample_data.sh . It contains metadata from the NCBI Biosample
    database, including sample submission dates.

all_SRA_introns.tsv.gz is specified as argument of --junctions. Annotations
are read from arguments of command-line parameter --annotations that specify
paths to the GTFs above.

Each line of all_SRA_introns.tsv.gz specifies a different junction and has the
following tab-separated fields.
1. chromosome
2. start position (1-based inclusive)
3. end position (1-based inclusive)
4. strand (+ or -)
5. 5' motif (GT, GC, or AT)
6. 3' motif (AG or AC)
7. comma-separated list of indexes of samples in which junction was found
8. comma-separated list of counts of reads overlapping junctions in
    corresponding sample from field 7. So if field 7 is 4,5,6 and field 8 is
    9,10,11 there are 9 reads overlapping the junction in the sample with
    index 4, 10 reads overlapping the junction in the sample with index 5, and
    11 reads overlapping the junction in the sample with index 6.

Each line of index_to_SRA_accession.tsv specifies a different sample
(specifically, run) on SRA and has the following tab-separated fields.
1. sample index
2. project accession number (regex: [SED]RP\d+)
3. sample accession number (regex: [SED]RS\d+)
4. experiment accession number (regex: [SED]RX\d+)
5. run accession number (regex: [SED]RR\d+)

We used PyPy 2.5.0 with GCC 4.9.2 for our Python implementation and ran:
pypy gi2015.py
    --hisat2-dir /path/to/hisat2-2.0.0-beta
    --star /path/to/STAR
    --annotations /path/to/Homo_sapiens.GRCh37.75.gtf
                  /path/to/gencode.v19.annotation.gtf
                  /path/to/refGene.gtf
    --junctions /path/to/all_SRA_introns.tsv.gz
    --index-to-sra /path/to/index_to_SRA_accession.tsv
    --tmp /path/to/temp_dir_with_200_GB_free_space
    --sort /path/to/coreutils_v8.23/sort
    --rgasp-dir /path/to/dir_with_4_RGASP_FASTQ.GZs_and_2_BAMs
    --cores 8
    --hisat2-idx /path/to/hisat2_idx
    --star-idx /path/to/star_idx
    --seqc /path/to/nbt.2957-S4.zip
    --sra-metadata /path/to/all_illumina_sra_for_human.tsv.gz
    --biosample-metadata /path/to/biosample_tags.tsv

Note that the argument of --hisat2-dir is the directory containing the HISAT 2
binary, but the argument of --STAR is precisely the STAR binary.

The following output was obtained. It is included in this repo because this 
script cannot currently be rerun to obtain results; the input file
all_SRA_introns.tsv.gz is not yet public. Note that an "overlap" below is
an instance where a junction is overlapped by a read. A read that overlaps
two exon-exon junctions contributes two overlaps (or overlap instances).

[basename].venn.txt
Intersections between junctions obtained from SEQC protocol and junctions
obtained from Rail for the 1720 samples studied by both SEQC and Rail. See
file for details.

[basename].sim.txt
Results of running HISAT2 2.0.0-beta and STAR 2.4.2a on the two simulated
human samples from the RGASP spliced alignment paper
"Systematic evaluation of spliced alignment programs for RNA-seq data"
(40 million read pairs each). Junction and overlap accuracy as defined
in the Rail-RNA preprint are given for two protocols: in one, the junctions
that are actually in the sample are provided as annotation, and in the other,
the union of the three gene annotations considered here
(Gencode v19, Ensembl v75, and refGene) are provided.

[basename].sex.txt
Tab-separated fields, in order of descending field 9
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

[basename].[type].tsv, where [type] is in [bottom_10_pct, top_10_pct,
                                    all_words]
Gives common words among samples in bottom and top 10 percent in terms of
proportion of junctions that are annotated -- as well as for all samples.
Only samples with >= 10k junctions are considered.
Tab-separated fields, in order of descending field 2
1. word
2. number of samples (in which >=10k junctions were found) in which word
    appears

[basename].sample_count_submission_date_overlap_geq_20.tsv
Tab-separated fields
1. count of samples in which a given junction was found
2. count of projects in which a given junction was found
3. earliest known discovery date (in units of days after February 27, 2009)
    -- this is the earliest known submission date of a sample associated with a
    junction

Above, each junction is covered by at least 20 reads per sample.

[basename].[type].stats.tsv, where [type] is in [project, sample]
Tab-separated fields
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

[basename].stats_by_sample.tsv
Tab-separated fields
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
"""
import sys
import gzip
import zipfile
import re
import os
import subprocess
from contextlib import contextmanager

def parsed_md(md):
    """ Divides an MD string up by boundaries between ^, letters, and numbers

        md: an MD string (example: 33A^CC).

        Return value: MD string split by boundaries described above.
    """
    md_to_parse = []
    md_group = [md[0]]
    for i, char in enumerate(md):
        if i == 0: continue
        if (re.match('[A-Za-z]', char) is not None) \
            != (re.match('[A-Za-z]', md[i-1]) is not None) or \
            (re.match('[0-9]', char) is not None) \
            != (re.match('[0-9]', md[i-1]) is not None):
            if md_group:
                md_to_parse.append(''.join(md_group))
            md_group = [char]
        else:
            md_group.append(char)
    if md_group:
        md_to_parse.append(''.join(md_group))
    return [char for char in md_to_parse if char != '0']

def indels_introns_and_exons(cigar, md, pos, seq):
    """ Computes indels, introns, and exons from CIGAR, MD string,
        and POS of a given alignment.

        cigar: CIGAR string
        md: MD:Z string
        pos: position of first aligned base
        seq: read sequence

        Return value: tuple (insertions, deletions, introns, exons). Insertions
            is a list of tuples (last genomic position before insertion, 
                                 string of inserted bases). Deletions
            is a list of tuples (first genomic position of deletion,
                                 string of deleted bases). Introns is a list
            of tuples (intron start position (inclusive),
                       intron end position (exclusive),
                       left_diplacement, right_displacement). Exons is a list
            of tuples (exon start position (inclusive),
                       exon end position (exclusive)).
    """
    if cigar == '*':
        return ([], [], [], [])
    insertions, deletions, introns, exons = [], [], [], []
    cigar = re.split(r'([MINDS])', cigar)[:-1]
    md = parsed_md(md)
    seq_size = len(seq)
    cigar_chars, cigar_sizes = [], []
    cigar_index, md_index, seq_index = 0, 0, 0
    max_cigar_index = len(cigar)
    while cigar_index != max_cigar_index:
        if cigar[cigar_index] == 0:
            cigar_index += 2
            continue
        if cigar[cigar_index+1] == 'M':
            aligned_base_cap = int(cigar[cigar_index])
            aligned_bases = 0
            while True:
                try:
                    aligned_bases += int(md[md_index])
                    if aligned_bases <= aligned_base_cap:
                        md_index += 1
                except ValueError:
                    # Not an int, but should not have reached a deletion
                    assert md[md_index] != '^', '\n'.join(
                                                ['cigar and md:',
                                                 ''.join(cigar), ''.join(md)]
                                            )
                    if aligned_bases + len(md[md_index]) > aligned_base_cap:
                        md[md_index] = md[md_index][
                                            :aligned_base_cap-aligned_bases
                                        ]
                        aligned_bases = aligned_base_cap
                    else:
                        aligned_bases += len(md[md_index])
                        md_index += 1
                if aligned_bases > aligned_base_cap:
                    md[md_index] = aligned_bases - aligned_base_cap
                    break
                elif aligned_bases == aligned_base_cap:
                    break
            # Add exon
            exons.append((pos, pos + aligned_base_cap))
            pos += aligned_base_cap
            seq_index += aligned_base_cap
        elif cigar[cigar_index+1] == 'N':
            skip_increment = int(cigar[cigar_index])
            # Add intron
            introns.append((pos, pos + skip_increment,
                            seq_index, seq_size - seq_index))
            # Skip region of reference
            pos += skip_increment
        elif cigar[cigar_index+1] == 'I':
            # Insertion
            insert_size = int(cigar[cigar_index])
            insertions.append(
                    (pos - 1, seq[seq_index:seq_index+insert_size])
                )
            seq_index += insert_size
        elif cigar[cigar_index+1] == 'D':
            assert md[md_index] == '^', '\n'.join(
                                                ['cigar and md:',
                                                 ''.join(cigar), ''.join(md)]
                                            )
            # Deletion
            delete_size = int(cigar[cigar_index])
            md_delete_size = len(md[md_index+1])
            assert md_delete_size >= delete_size
            deletions.append((pos, md[md_index+1][:delete_size]))
            if md_delete_size > delete_size:
                # Deletion contains an intron
                md[md_index+1] = md[md_index+1][delete_size:]
            else:
                md_index += 2
            # Skip deleted part of reference
            pos += delete_size
        else:
            # Soft clip
            assert cigar[cigar_index+1] == 'S'
            # Advance seq_index
            seq_index += int(cigar[cigar_index])
        cigar_index += 2
    '''Merge exonic chunks/deletions; insertions/introns could have chopped
    them up.'''
    new_exons = []
    last_exon = exons[0]
    for exon in exons[1:]:
        if exon[0] == last_exon[1]:
            # Merge ECs
            last_exon = (last_exon[0], exon[1])
        else:
            # Push last exon to new exon list
            new_exons.append(last_exon)
            last_exon = exon
    new_exons.append(last_exon)
    return insertions, deletions, introns, new_exons

def dummy_md_index(cigar):
    """ Creates dummy MD string from CIGAR in case of missing MD.

        cigar: cigar string

        Return value: dummy MD string
    """
    cigar = re.split(r'([MINDS])', cigar)[:-1]
    cigar_index = 0
    max_cigar_index = len(cigar)
    md = []
    while cigar_index != max_cigar_index:
        if cigar[cigar_index] == 0:
            cigar_index += 2
            continue
        if cigar[cigar_index+1] == 'M':
            try:
                if type(md[-1]) is int:
                    md[-1] += int(cigar[cigar_index])
                else:
                    md.append(int(cigar[cigar_index]))
            except IndexError:
                md.append(int(cigar[cigar_index]))
            cigar_index += 2
        elif cigar[cigar_index+1] in 'SIN':
            cigar_index += 2
        elif cigar[cigar_index+1] == 'D':
            md.extend(['^', 'A'*int(cigar[cigar_index])])
            cigar_index += 2
        else:
            raise RuntimeError(
                        'Accepted CIGAR characters are only in [MINDS].'
                    )
    return ''.join(str(el) for el in md)

def junctions_from_sam_stream(sam_stream, include_strand=False):
    """ Obtains set of junctions from SAM file

        ONLY PRIMARY ALIGNMENTS ARE CONSIDERED.

        sam_stream: where to find retrieved alignments in SAM form
        include_strand: True iff strand should be included in returned set

        Return value: set of tuples (chromosome, 1-based junction start
            position (inclusive), 1-based junction end position (exclusive),
            strand (+ or -) iff strand is True)
    """
    junctions_to_return = set()
    for line in sam_stream:
        if line[0] == '@': continue
        try:
            tokens = line.strip().split('\t')
            flag = int(tokens[1])
            if flag & 4:
                continue
            cigar = tokens[5]
            if 'N' not in cigar or flag & 256:
                continue
            rname = tokens[2]
            pos = int(tokens[3])
            seq = tokens[9]
            _, _, junctions, _ = indels_introns_and_exons(
                                        cigar,
                                        dummy_md_index(cigar), pos, seq
                                    )
            if include_strand:
                strand = [
                        field for field in tokens if field[:5] == 'XS:A:'
                    ][0][-1]
                junctions_to_return.update([(rname,) + junction[:2] + (strand,)
                                            for junction in junctions])
            else:
                junctions_to_return.update([(rname,) + junction[:2]
                                            for junction in junctions])
        except IndexError:
            print >>sys.stderr, ('Error found on line: ' + line)
            raise
    return junctions_to_return

def is_gzipped(filename):
    """ Uses gzip magic number to determine whether a file is compressed.

        filename: path to file

        Return value: True iff file filename is gzipped.
    """
    with open(filename, 'rb') as binary_input_stream:
        # Check for magic number
        if binary_input_stream.read(2) == '\x1f\x8b':
            return True
        else:
            return False

@contextmanager
def xopen(filename):
    """ Opens both gzipped and uncompressed files for contextual reading.

        filename: path to file to open

        Yield value: a gzip.open or open object
    """
    if is_gzipped(filename):
        f = gzip.open(filename)
    else:
        f = open(filename)
    try:
        yield f
    finally:
        f.close()

if __name__ == '__main__':
    import argparse
    # Print file's docstring if -h is invoked
    parser = argparse.ArgumentParser(description=__doc__, 
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # Add command-line arguments
    parser.add_argument('--hisat2-dir', type=str, required=True,
            help=('path to directory containing contents of HISAT2; we '
                  'unpacked ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/'
                  'downloads/hisat2-2.0.0-beta-Linux_x86_64.zip to get this')
        )
    parser.add_argument('--star', type=str, required=True,
            help='path to STAR binary'
        )
    parser.add_argument('--star-idx', type=str, required=True,
            help='path to directory containing STAR index for hg19; see '
                 'docstring for more information'
        )
    parser.add_argument('--hisat2-idx', type=str, required=True,
            help='path to basename of HISAT2 index for hg19; see '
                 'docstring for more information'
        )
    parser.add_argument('--samtools', type=str, required=False,
            default='samtools',
            help='path to SAMTools binary'
        )
    parser.add_argument('--annotations', type=str, required=True, nargs='+',
            help='space-separated paths to GTF files encoding known junctions'
        )
    parser.add_argument('--junctions', type=str, required=True,
            help='junctions file; this should be all_SRA_introns.tsv.gz'
        )
    parser.add_argument('--index-to-sra', type=str, required=True,
            help='index to SRA accession numbers file; this should be '
                 'index_to_SRA_accession.tsv'
        )
    parser.add_argument('--sra-metadata', type=str, required=True,
            help='path to SRA metadata file; this should be '
                 'all_illumina_sra_for_human.tsv.gz'
        )
    parser.add_argument('--biosample-metadata', type=str, required=True,
            help='path to Biosample metadata file; this should be '
                 'biosample_tags.tsv'
        )
    parser.add_argument('--seqc', type=str, required=True,
            help='path to SEQC junctions; this should be nbt.2957-S4.zip')
    parser.add_argument('--rgasp-dir', type=str, required=True,
            help=('path to directory containing 6 files from RGASP paper; see '
                  'docstring for which files this should be')
        )
    parser.add_argument('--tmp', type=str, required=False,
            default='/var/tmp',
            help='where to write temporary files; make sure to have ~50GB free'
        )
    parser.add_argument('--cores', type=int, required=False,
            default=8,
            help='number of HISAT2/STAR threads to use'
        )
    parser.add_argument('--basename', type=str, required=False,
            default='gi2015',
            help='basename for output files'
        )
    parser.add_argument('--sort', type=str, required=False,
            default='sort',
            help='path to GNU coreutils (standard UNIX) sort'
        )
    parser.add_argument('--no-sims', action='store_const', const=True,
            default=False,
            help='skip the HISAT2/STAR simulations'
        )
    args = parser.parse_args()

    from collections import defaultdict

    # Map sample indexes to accession number lines
    index_to_sra, index_to_srp, srr_to_index = {}, {}, {}
    srs_to_srr = defaultdict(list)
    # Get sample indexes for all Illumina RNA-seq from SEQC for Venn diagram
    seqc_indexes = set()
    with xopen(args.index_to_sra) as index_stream:
        for line in index_stream:
            partitioned = line.partition('\t')
            sample_index = int(partitioned[0])
            index_to_sra[sample_index] = partitioned[2].strip()
            srp, srs, srx, srr = partitioned[2].strip().split('\t')
            srs_to_srr[srs].append(srr)
            srr_to_index[srr] = sample_index
            index_to_srp[sample_index] = srp
            if srp == 'SRP025982':
                # SEQC hit!
                seqc_indexes.add(sample_index)
    print >>sys.stderr, 'Done mapping sample indexes to samples.'

    from datetime import date
    '''For getting junctions by "earliest detection date"; use units of number
    of days after earliest date. Map sample indexes to submission dates.'''
    all_dates = {}
    with xopen(args.biosample_metadata) as biosample_stream:
        biosample_stream.readline() # header
        for line in biosample_stream:
            tokens = line.strip().split('\t')
            current_date = date(
                    *tuple(
                        [int(el.strip())
                            for el in tokens[10].split('T')[0].split('-')]
                    )
                )
            for srr in srs_to_srr[tokens[9].upper()]:
                all_dates[srr_to_index[srr]] = current_date
    earliest_date = min(all_dates.values())
    for sample_index in all_dates:
        all_dates[sample_index] = (
                all_dates[sample_index] - earliest_date
            ).days
    date_indexes = set(all_dates.keys())
    print >>sys.stderr, 'Done grabbing submission dates from Biosample DB.'

    # Grab all annotated junctions
    annotated_junctions = set()
    annotated_junctions_for_star = set()
    annotated_junctions_for_hisat = set()
    annotated_5p = set()
    annotated_3p = set()
    refs = set(
            ['chr' + str(i) for i in xrange(1, 23)] + ['chrM', 'chrX', 'chrY']
        )
    extended_refs = refs | set([
            'chr6_ssto_hap7', 'chr6_mcf_hap5', 'chr6_cox_hap2',
            'chr6_mann_hap4', 'chr6_apd_hap1', 'chr6_qbl_hap6'
            'chr6_dbb_hap3', 'chr17_ctg5_hap1', 'chr4_ctg9_hap1',
            'chr1_gl000192_random', 'chrUn_gl000225', 'chr4_gl000194_random',
            'chr4_gl000193_random', 'chr9_gl000200_random', 'chrUn_gl000222',
            'chrUn_gl000212', 'chr7_gl000195_random', 'chrUn_gl000223',
            'chrUn_gl000224', 'chrUn_gl000219', 'chr17_gl000205_random',
            'chrUn_gl000215', 'chrUn_gl000216', 'chrUn_gl000217',
            'chr9_gl000199_random', 'chrUn_gl000211', 'chrUn_gl000213',
            'chrUn_gl000220', 'chrUn_gl000218', 'chr19_gl000209_random',
            'chrUn_gl000221', 'chrUn_gl000214', 'chrUn_gl000228',
            'chrUn_gl000227', 'chr1_gl000191_random', 'chr19_gl000208_random',
            'chr9_gl000198_random', 'chr17_gl000204_random', 'chrUn_gl000233',
            'chrUn_gl000237', 'chrUn_gl000230', 'chrUn_gl000242',
            'chrUn_gl000243', 'chrUn_gl000241', 'chrUn_gl000236',
            'chrUn_gl000240', 'chr17_gl000206_random', 'chrUn_gl000232',
            'chrUn_gl000234', 'chr11_gl000202_random', 'chrUn_gl000238',
            'chrUn_gl000244', 'chrUn_gl000248', 'chr8_gl000196_random',
            'chrUn_gl000249', 'chrUn_gl000246', 'chr17_gl000203_random',
            'chr8_gl000197_random', 'chrUn_gl000245', 'chrUn_gl000247',
            'chr9_gl000201_random', 'chrUn_gl000235', 'chrUn_gl000239',
            'chr21_gl000210_random', 'chrUn_gl000231', 'chrUn_gl000229',
            'chrUn_gl000226', 'chr18_gl000207_random'
        ])
    extract_splice_sites_path = os.path.join(args.hisat2_dir,
                                                'extract_splice_sites.py')
    for annotation in args.annotations:
        extract_process = subprocess.Popen(' '.join([
                                            sys.executable,
                                            extract_splice_sites_path,
                                            annotation
                                            if not is_gzipped(
                                               annotation
                                            ) else ('<(gzip -cd %s)'
                                               % annotation
                                            )
                                        ]),
                                        shell=True,
                                        executable='/bin/bash',
                                        stdout=subprocess.PIPE
                                    )
        for line in extract_process.stdout:
            tokens = line.strip().split('\t')
            tokens[1] = int(tokens[1]) + 2
            tokens[2] = int(tokens[2])
            if not tokens[0].startswith('chr'):
                tokens[0] = 'chr' + tokens[0]
            if tokens[0] in refs:
                annotated_junctions.add(tuple(tokens[:-1]))
                if tokens[3] == '+':
                    annotated_5p.add((tokens[0], tokens[1]))
                    annotated_3p.add((tokens[0], tokens[2]))
                else:
                    assert tokens[3] == '-'
                    annotated_3p.add((tokens[0], tokens[1]))
                    annotated_5p.add((tokens[0], tokens[2]))
            if tokens[0] in extended_refs:
                annotated_junctions_for_hisat.add(
                        (tokens[0], tokens[1] - 2, tokens[2], tokens[3])
                    )
                annotated_junctions_for_star.add(tuple(tokens))
        extract_process.stdout.close()
        exit_code = extract_process.wait()
        if exit_code != 0:
            raise RuntimeError(
                'extract_splice_sites.py had nonzero exit code {}.'.format(
                                                                    exit_code
                                                                )
            )
    print >>sys.stderr, 'Found %d annotated junctions.' % (
            len(annotated_junctions)
        )

    if not args.no_sims:
        # Schedule a temporary directory for deletion on script termination
        import tempfile
        import atexit
        import shutil
        try:
            temp_dir = tempfile.mkdtemp(dir=args.tmp)
        except:
            temp_dir = tempfile.mkdtemp()
        atexit.register(shutil.rmtree, temp_dir)

        # Get the truth for both simulations
        true_junctions = []
        for i in range(2):
            samtools_process = subprocess.Popen(
                                [args.samtools, 'view',
                                    os.path.join(args.rgasp_dir,
                                                'Sim%d_truth.bam' % (i+1))],
                                        stdout=subprocess.PIPE
                                    )
            true_junctions.append(
                    junctions_from_sam_stream(samtools_process.stdout,
                                                include_strand=True)
                )
            samtools_process.terminate()
            samtools_process.wait()
            # Sort by QNAME
            subprocess.check_call(' '.join(
                                    [args.samtools, 'view',
                                        os.path.join(args.rgasp_dir,
                                            'Sim%d_truth.bam' % (i+1)), '|',
                                        args.sort, '-T %s -k1,1' % temp_dir,
                                        '>%s' % os.path.join(temp_dir,
                                                    'Sim%d_truth.sorted.sam'
                                                     % (i+1))]),
                                shell=True,
                                executable='/bin/bash')

        print >>sys.stderr, 'Done reading true junctions from simulations.'
        '''Create 2 junction databases for each aligner: the first has the 
        union of all annotated junctions provided as the argument of
        --annotations, the second has exactly the true junctions.'''
        hisat2_annotated_junctions_file = os.path.join(temp_dir,
                                                        'hisat2_junctions.tsv')
        hisat2_true_junctions_files = [os.path.join(temp_dir,
                                                'hisat2_true_junctions%d.tsv'
                                                    % i) for i in range(2)]
        star_annotated_junctions_file = os.path.join(temp_dir,
                                                        'star_junctions.tsv')
        star_true_junctions_files = [os.path.join(temp_dir,
                                                    'star_true_junctions%d.tsv'
                                                        % i) for i in range(2)]
        with open(hisat2_annotated_junctions_file, 'w') as junction_stream:
            print >>junction_stream, '\n'.join(
                    ['\t'.join([str(el) for el in junction])
                        for junction in annotated_junctions_for_hisat]
                )
        with open(star_annotated_junctions_file, 'w') as junction_stream:
            print >>junction_stream, '\n'.join(
                    ['\t'.join([str(el) for el in junction])
                        for junction in annotated_junctions_for_star]
                )
        for i in range(2):
            with open(hisat2_true_junctions_files[i], 'w') as junction_stream:
                for junction in true_junctions[i]:
                    print >>junction_stream, (
                            '%s\t%d\t%d\t%s' % (junction[0],
                                                    int(junction[1]) - 2,
                                                    int(junction[2]) - 1,
                                                    junction[3])
                        )
            with open(star_true_junctions_files[i], 'w') as junction_stream:
                for junction in true_junctions[i]:
                    print >>junction_stream, (
                            '%s\t%s\t%d\t%s' % (junction[0],
                                                    junction[1],
                                                    int(junction[2]) - 1,
                                                    junction[3])
                        )

        '''Run two kinds of simulations: one provides annotated junctions to 
        the aligner; the other provides the true junctions to the aligner.'''
        for i in range(2):
            sim_dir = os.path.join(temp_dir, str(i))
            star_true = os.path.join(sim_dir, 'star_true')
            hisat2_true = os.path.join(sim_dir, 'hisat2_true')
            star_annotated = os.path.join(sim_dir, 'star_annotated')
            hisat2_annotated = os.path.join(sim_dir, 'hisat2_annotated')
            os.makedirs(star_true)
            os.makedirs(hisat2_true)
            os.makedirs(star_annotated)
            os.makedirs(hisat2_annotated)
            star_command = ('set -exo pipefail; cd {working_dir}; '
                            '{star_exe} --genomeDir {star_idx} --readFilesIn '
                            '{left_reads} {right_reads} --runThreadN {cores} '
                            '--sjdbFileChrStartEnd {junction_db} '
                            '--readFilesCommand zcat --outSAMunmapped Within; '
                            'cat Aligned.out.sam | awk '
                            '\'substr($0,1,1)!="@"\' | '
                            '{sort_exe} -T {temp_dir} '
                            '-k1,1 >Aligned.out.sorted.sam')
            hisat2_command = ('set -exo pipefail; cd {working_dir}; '
                              '{hisat2_exe} -x {hisat2_idx} '
                              '-1 {left_reads} -2 {right_reads} '
                              '-S Aligned.out.sam -p {cores} --dta-cufflinks '
                              '--phred64 '
                              '--novel-splicesite-infile {junction_db}; '
                              'cat Aligned.out.sam | awk '
                              '\'substr($0,1,1)!="@"\' '
                              '| {sort_exe} -T {temp_dir} -k1,1 '
                              '>Aligned.out.sorted.sam')
            star_true_command = star_command.format(
                                    working_dir=star_true,
                                    star_exe=args.star,
                                    star_idx=args.star_idx,
                                    left_reads=os.path.join(
                                        args.rgasp_dir,
                                        'Human_simulated_reads_%d_PE1.fq.gz'
                                        % (i+1)
                                    ),
                                    right_reads=os.path.join(
                                        args.rgasp_dir,
                                        'Human_simulated_reads_%d_PE2.fq.gz'
                                        % (i+1)
                                    ),
                                    cores=args.cores,
                                    junction_db=star_true_junctions_files[i],
                                    sort_exe=args.sort,
                                    temp_dir=temp_dir
                                )
            hisat2_true_command = hisat2_command.format(
                                    working_dir=hisat2_true,
                                    hisat2_exe=os.path.join(args.hisat2_dir,
                                                            'hisat2'),
                                    hisat2_idx=args.hisat2_idx,
                                    left_reads=os.path.join(
                                        args.rgasp_dir,
                                        'Human_simulated_reads_%d_PE1.fq.gz'
                                        % (i+1)
                                    ),
                                    right_reads=os.path.join(
                                        args.rgasp_dir,
                                        'Human_simulated_reads_%d_PE2.fq.gz'
                                        % (i+1)
                                    ),
                                    cores=args.cores,
                                    junction_db=hisat2_true_junctions_files[i],
                                    sort_exe=args.sort,
                                    temp_dir=temp_dir
                                )
            subprocess.check_call(star_true_command,
                                        shell=True,
                                        executable='/bin/bash')
            subprocess.check_call(hisat2_true_command,
                                        shell=True,
                                        executable='/bin/bash')
            star_annotated_command = star_command.format(
                                    working_dir=star_annotated,
                                    star_exe=args.star,
                                    star_idx=args.star_idx,
                                    left_reads=os.path.join(
                                        args.rgasp_dir,
                                        'Human_simulated_reads_%d_PE1.fq.gz'
                                        % (i+1)
                                    ),
                                    right_reads=os.path.join(
                                        args.rgasp_dir,
                                        'Human_simulated_reads_%d_PE2.fq.gz'
                                        % (i+1)
                                    ),
                                    cores=args.cores,
                                    junction_db=star_annotated_junctions_file,
                                    sort_exe=args.sort,
                                    temp_dir=temp_dir
                                )
            hisat2_annotated_command = hisat2_command.format(
                                working_dir=hisat2_annotated,
                                hisat2_exe=os.path.join(args.hisat2_dir,
                                                        'hisat2'),
                                hisat2_idx=args.hisat2_idx,
                                left_reads=os.path.join(
                                        args.rgasp_dir,
                                        'Human_simulated_reads_%d_PE1.fq.gz'
                                        % (i+1)
                                    ),
                                right_reads=os.path.join(
                                        args.rgasp_dir,
                                        'Human_simulated_reads_%d_PE2.fq.gz'
                                        % (i+1)
                                    ),
                                cores=args.cores,
                                junction_db=hisat2_annotated_junctions_file,
                                sort_exe=args.sort,
                                temp_dir=temp_dir
                            )
            subprocess.check_call(star_annotated_command,
                                    shell=True,
                                    executable='/bin/bash')
            subprocess.check_call(hisat2_annotated_command,
                                        shell=True,
                                        executable='/bin/bash')

        '''Compute accuracy measures: junction precision and overlap precision.
        See Rail-RNA preprint at
        http://biorxiv.org/content/early/2015/08/11/019067 for appropriate
        definitions.'''
        with open(args.basename + '.sim.txt', 'w') as output_stream:
            for i in range(2):
                for aligner in ['hisat2', 'star']:
                    for sim_type in ['true', 'annotated']:
                        print >>output_stream, ('sim %d, %s junctions, %s'
                                                    % (i+1, sim_type, aligner))
                        print >>output_stream, '\njunction accuracy'
                        with open(os.path.join(temp_dir, str(i),
                                                    aligner + '_' + sim_type,
                                                    'Aligned.out.sorted.sam'
                                                )) as retrieved_stream:
                            retrieved = junctions_from_sam_stream(
                                                        retrieved_stream,
                                                        include_strand=False
                                                    )
                        relevant_instances = len(true_junctions[i])
                        retrieved_instances = len(retrieved)
                        relevant_and_retrieved_instances = len(
                                set(
                                        [el[:3] for el in true_junctions[i]]
                                ).intersection(retrieved)
                            )
                        precision = (float(relevant_and_retrieved_instances)
                                        / retrieved_instances)
                        recall = (float(relevant_and_retrieved_instances)
                                        / relevant_instances)
                        fscore = 2 * precision * recall / (precision + recall)
                        print >>output_stream, 'precision: %.08f' % precision
                        print >>output_stream, 'recall: %.08f' % recall
                        print >>output_stream, 'fscore: %.08f\n' % fscore
                        with open(
                                os.path.join(temp_dir,
                                             'Sim%d_truth.sorted.sam' % (i+1))
                            ) as true_stream, open(
                                os.path.join(temp_dir, str(i),
                                             aligner + '_' + sim_type,
                                             'Aligned.out.sorted.sam')
                            ) as retrieved_stream:
                            (retrieved_instances, relevant_instances,
                                relevant_and_retrieved_instances) = 0, 0, 0
                            while True:
                                t1 = true_stream.readline().strip().split('\t')
                                if not t1[0]: break
                                t2 = true_stream.readline().strip().split('\t')
                                r1 = retrieved_stream.readline(
                                            ).strip().split('\t')
                                while int(r1[1]) & 256:
                                    r1 = retrieved_stream.readline(
                                            ).strip().split('\t')
                                r2 = retrieved_stream.readline(
                                            ).strip().split('\t')
                                while int(r2[1]) & 256:
                                    r2 = retrieved_stream.readline(
                                            ).strip().split('\t')
                                assert t1[0] == t2[0] == r1[0] == r2[0], (
                                        '%s\n%s\n%s\n%s' % (t1, t2, r1, r2)
                                    )
                                _, _, t1_introns, _ = indels_introns_and_exons(
                                                        t1[5],
                                                        dummy_md_index(t1[5]),
                                                        int(t1[3]), t1[9]
                                                    )
                                _, _, t2_introns, _ = indels_introns_and_exons(
                                                        t2[5],
                                                        dummy_md_index(t2[5]),
                                                        int(t2[3]), t2[9]
                                                    )
                                _, _, r1_introns, _ = indels_introns_and_exons(
                                                        r1[5],
                                                        dummy_md_index(r1[5]),
                                                        int(r1[3]), r1[9]
                                                    )
                                _, _, r2_introns, _ = indels_introns_and_exons(
                                                        r2[5],
                                                        dummy_md_index(r2[5]),
                                                        int(r2[3]), r2[9]
                                                    )
                                t1_introns = [
                                        (t1[2],) + el[:2] for el in t1_introns
                                    ]
                                t2_introns = [
                                        (t2[2],) + el[:2] for el in t2_introns
                                    ]
                                r1_introns = [
                                        (r1[2],) + el[:2] for el in r1_introns
                                    ]
                                r2_introns = [
                                        (r2[2],) + el[:2] for el in r2_introns
                                    ]
                                relevant_instances += (
                                        len(t1_introns) + len(t2_introns)
                                    )
                                retrieved_instances += (
                                        len(r1_introns) + len(r2_introns)
                                    )
                                relevant_and_retrieved_instances += max(
                                        len([intron for intron in t1_introns
                                             if intron in r1_introns]) + 
                                        len([intron for intron in t2_introns
                                             if intron in r2_introns]),
                                        len([intron for intron in t1_introns
                                             if intron in r2_introns]) +
                                        len([intron for intron in t2_introns
                                             if intron in r1_introns])
                                    )
                        precision = (float(relevant_and_retrieved_instances)
                                        / retrieved_instances)
                        recall = (float(relevant_and_retrieved_instances)
                                        / relevant_instances)
                        fscore = 2 * precision * recall / (precision + recall)
                        print >>output_stream, '\noverlap accuracy'
                        print >>output_stream, 'precision: %.08f' % precision
                        print >>output_stream, 'recall: %.08f' % recall
                        print >>output_stream, 'fscore: %.08f\n' % fscore
        print >>sys.stderr, 'Done all HISAT2 and STAR simulations.'

    '''Grab SEQC junctions. Three protocols were used: Subread, r-make, and
    NCBI Magic.''' 
    magic_junctions, rmake_junctions, subread_junctions = set(), set(), set()
    seqc_junctions = set()
    with zipfile.ZipFile(args.seqc).open('SupplementaryData3.tab') \
        as seqc_stream:
        seqc_stream.readline() # header
        for line in seqc_stream:
            tokens = line.strip().split('\t')
            tokens[0] = tokens[0].split('.')
            junction = (tokens[0][0], int(tokens[0][1]), int(tokens[0][2]))
            add_junc = False
            if tokens[1] == '1':
                subread_junctions.add(junction)
                add_junc = True
            if tokens[2] == '1':
                rmake_junctions.add(junction)
                add_junc = True
            if tokens[3] == '1':
                magic_junctions.add(junction)
                add_junc = True
            if add_junc:
                seqc_junctions.add(junction)
    print >>sys.stderr, 'Done reading SEQC junctions.'

    # Key: sample index; value: number of junctions found in sample
    junction_counts = defaultdict(int)
    # For junctions in union of annotations specified at command line
    annotated_junction_counts = defaultdict(int)
    # Count total overlap instances and annotated overlap instances
    '''Same as above, but including only junctions covered by at least 5 reads
    in the sample.'''
    junction_counts_geq_5 = defaultdict(int)
    annotated_junction_counts_geq_5 = defaultdict(int)

    overlap_counts = defaultdict(int)
    annotated_overlap_counts = defaultdict(int)
    # Want Y chromosome overlaps to do sex classification
    Y_overlap_counts = defaultdict(int)
    Y_junction_counts = defaultdict(int)
    Xist_junction_counts = defaultdict(int)
    Xist_overlap_counts = defaultdict(int)
    # Mapping counts of samples to junction counts
    sample_count_to_junction_count = defaultdict(int)
    project_count_to_junction_count = defaultdict(int)
    sample_count_to_GTAG_junction_count = defaultdict(int)
    project_count_to_GTAG_junction_count = defaultdict(int)
    sample_count_to_GCAG_junction_count = defaultdict(int)
    project_count_to_GCAG_junction_count = defaultdict(int)
    sample_count_to_ATAC_junction_count = defaultdict(int)
    project_count_to_ATAC_junction_count = defaultdict(int)
    sample_count_to_GTAG_ann_count = defaultdict(int)
    project_count_to_GTAG_ann_count = defaultdict(int)
    sample_count_to_GCAG_ann_count = defaultdict(int)
    project_count_to_GCAG_ann_count = defaultdict(int)
    sample_count_to_ATAC_ann_count = defaultdict(int)
    project_count_to_ATAC_ann_count = defaultdict(int)
    # One of 5' or 3' splice site is in annotation, one isn't
    sample_count_to_altstartend_junction_count = defaultdict(int)
    project_count_to_altstartend_junction_count = defaultdict(int)
    # Both 5' and 3' splice sites are in annotation, but junction is not
    sample_count_to_exonskip_junction_count = defaultdict(int)
    project_count_to_exonskip_junction_count = defaultdict(int)
    # Full junction is in annotation
    sample_count_to_annotated_junction_count = defaultdict(int)
    project_count_to_annotated_junction_count = defaultdict(int)
    # Neither 5' nor 3' is in annotation
    sample_count_to_novel_junction_count = defaultdict(int)
    project_count_to_novel_junction_count = defaultdict(int)
    rail_seqc_junctions = set()
    rail_seqc_junctions_geq_5_samples = set()
    rail_seqc_junctions_geq_10_samples = set()
    rail_seqc_junctions_geq_15_samples = set()
    rail_seqc_junctions_geq_20_samples = set()
    rail_seqc_junctions_geq_25_samples = set()
    rail_seqc_junctions_geq_30_samples = set()
    # For junction-date analyses
    date_to_junction_count = defaultdict(int)
    date_to_junction_count_overlap_geq_20 = defaultdict(int)

    with xopen(args.junctions) as junction_stream, open(
            args.basename
            + '.sample_count_submission_date_overlap_geq_20.tsv', 'w'
        ) as junction_date_stream:
        print >>junction_date_stream, (
                               '# samples in which junction was found'
                               '\t# projects in which junction was found'
                               '\tearliest known discovery date in '
                               'days after %s; format Y-M-D') % (
                                        earliest_date.strftime(
                                            '%Y-%m-%d'
                                        )
                                    )
        for line in junction_stream:
            tokens = line.strip().split('\t')
            junction = (tokens[0], int(tokens[1]), int(tokens[2]))
            if tokens[3] == '+':
                fivep = junction[:2]
                threep = (junction[0], junction[2])
            elif tokens[3] == '-':
                threep = junction[:2]
                fivep = (junction[0], junction[2])
            else:
                raise RuntimeError('Bad strand in line "%s"' % line)
            samples = [int(el) for el in tokens[-2].split(',')]
            coverages = [int(el) for el in tokens[-1].split(',')]
            sample_count = len(samples)
            project_count = len(set([index_to_srp[sample]
                                        for sample in samples]))
            try:
                discovery_date = min(
                        [all_dates[sample] for sample in samples
                            if sample in date_indexes]
                    )
            except ValueError:
                # No discovery date available!
                pass
            else:
                date_to_junction_count[discovery_date] += 1
                if sum(coverages) >= 20:
                    date_to_junction_count_overlap_geq_20[discovery_date] += 1
                    print >>junction_date_stream, '%d\t%d\t%d' % (
                                                            sample_count,
                                                            project_count,
                                                            discovery_date
                                                        )
            samples_and_coverages = zip(samples, coverages)
            sample_count_to_junction_count[sample_count] += 1
            project_count_to_junction_count[project_count] += 1
            if tokens[5] == 'AG':
                if tokens[4] == 'GT':
                    sample_count_to_GTAG_junction_count[sample_count] += 1
                    project_count_to_GTAG_junction_count[project_count] += 1
                elif tokens[4] == 'GC':
                    sample_count_to_GCAG_junction_count[sample_count] += 1
                    project_count_to_GCAG_junction_count[project_count] += 1
                else:
                    raise RuntimeError('Bad motif in line "%s"' % line)
            elif tokens[5] == 'AC':
                if tokens[4] == 'AT':
                    sample_count_to_ATAC_junction_count[sample_count] += 1
                    project_count_to_ATAC_junction_count[project_count] += 1
                else:
                    raise RuntimeError('Bad motif in line "%s"' % line)
            if junction in annotated_junctions:
                sample_count_to_annotated_junction_count[sample_count] += 1
                project_count_to_annotated_junction_count[project_count] += 1
                if tokens[5] == 'AG':
                    if tokens[4] == 'GT':
                        sample_count_to_GTAG_ann_count[sample_count] += 1
                        project_count_to_GTAG_ann_count[project_count] += 1
                    elif tokens[4] == 'GC':
                        sample_count_to_GCAG_ann_count[sample_count] += 1
                        project_count_to_GCAG_ann_count[project_count] += 1
                elif tokens[5] == 'AC':
                    sample_count_to_ATAC_ann_count[sample_count] += 1
                    project_count_to_ATAC_ann_count[project_count] += 1
                for sample, coverage in samples_and_coverages:
                    annotated_junction_counts[sample] += 1
                    annotated_overlap_counts[sample] += coverage
                    if coverage >= 5:
                        annotated_junction_counts_geq_5[sample] += 1
            elif threep in annotated_3p:
                if fivep in annotated_5p:
                    sample_count_to_exonskip_junction_count[sample_count] += 1
                    project_count_to_exonskip_junction_count[
                            project_count
                        ] += 1
                else:
                    sample_count_to_altstartend_junction_count[sample_count] \
                        += 1
                    project_count_to_altstartend_junction_count[
                            project_count
                        ] += 1
            elif fivep in annotated_5p:
                sample_count_to_altstartend_junction_count[sample_count] += 1
                project_count_to_altstartend_junction_count[project_count] += 1
            else:
                sample_count_to_novel_junction_count[sample_count] += 1
                project_count_to_novel_junction_count[project_count] += 1
            seqc_intersect = set(samples).intersection(seqc_indexes)
            if seqc_intersect:
                in_seqc_sample_count = len(seqc_intersect)
                rail_seqc_junctions.add(junction)
                if in_seqc_sample_count >= 5:
                    rail_seqc_junctions_geq_5_samples.add(junction)
                    if in_seqc_sample_count >= 10:
                        rail_seqc_junctions_geq_10_samples.add(junction)
                        if in_seqc_sample_count >= 15:
                            rail_seqc_junctions_geq_15_samples.add(junction)
                            if in_seqc_sample_count >= 20:
                                rail_seqc_junctions_geq_20_samples.add(
                                        junction
                                    )
                                if in_seqc_sample_count >= 25:
                                    rail_seqc_junctions_geq_25_samples.add(
                                        junction
                                    )
                                    if in_seqc_sample_count >= 30:
                                        rail_seqc_junctions_geq_30_samples.add(
                                            junction
                                        )
            for sample, coverage in samples_and_coverages:
                junction_counts[sample] += 1
                overlap_counts[sample] += coverage
                if coverage >= 5:
                    junction_counts_geq_5[sample] += 1
            if junction[0] == 'chrY':
                for sample, coverage in samples_and_coverages:
                    Y_junction_counts[sample] += 1
                    Y_overlap_counts[sample] += coverage
            elif junction[0] == 'chrX':
                if junction[1] >= 73040486 and junction[2] <= 73072588:
                    # In XIST zone
                    for sample, coverage in samples_and_coverages:
                        Xist_junction_counts[sample] += 1
                        Xist_overlap_counts[sample] += coverage
    print >>sys.stderr, 'Done reading junction file.'

    import re
    # For sex analysis
    '''Grab male/female SRRs from SRA metadata; study only samples w/ >= 10k
    junctions'''
    to_print = []
    with xopen(args.sra_metadata) as sra_stream:
        sra_stream.readline() # header
        for line in sra_stream:
            try:
                sample_index = srr_to_index[line.partition('\t')[0]]
            except KeyError:
                # Not in there
                continue
            if junction_counts[sample_index] < 10000: continue
            male_present = (re.search(r"\bmale\b", line) is not None)
            female_present = (re.search(r"\bfemale\b", line) is not None)
            if male_present and not female_present:
                to_print.append((index_to_sra[sample_index],
                              Y_overlap_counts[sample_index],
                              Xist_overlap_counts[sample_index],
                              overlap_counts[sample_index],
                              Y_junction_counts[sample_index],
                              Xist_junction_counts[sample_index],
                              junction_counts[sample_index], 1,
                              float(Y_overlap_counts[sample_index])
                                / overlap_counts[sample_index]))
            elif female_present and not male_present:
                to_print.append((index_to_sra[sample_index],
                              Y_overlap_counts[sample_index],
                              Xist_overlap_counts[sample_index],
                              overlap_counts[sample_index],
                              Y_junction_counts[sample_index],
                              Xist_junction_counts[sample_index],
                              junction_counts[sample_index], 0,
                              float(Y_overlap_counts[sample_index])
                                / overlap_counts[sample_index]))
    to_print.sort(key=lambda x: x[-1], reverse=True)
    with open(args.basename + '.sex.tsv', 'w') as sex_stream:
        print >>sex_stream, ('project\t'
                             'sample\t'
                             'experiment\t'
                             'run\t'
                             'Y chr overlaps\t'
                             'XIST overlaps\t'
                             'total overlaps\t'
                             'Y chr junctions\t'
                             'XIST junctions\t'
                             'total junctions\t'
                             'is male?\t'
                             '(Y chr overlaps / total overlaps)')
        print >>sex_stream, '\n'.join(['%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.08f'
                                        % line for line in to_print])

    '''What kinds of samples have unannotated junctions? Let's study samples
    in the bottom ~10 percent in terms of proportion of junctions covered by 
    >= 5 reads that are annotated. Consider only samples with >= 10k 
    junctions.'''
    proportions = [(sample, float(annotated_junction_counts_geq_5[sample])
                                / junction_counts_geq_5[sample])
                        for sample in junction_counts_geq_5
                        if junction_counts[sample] >= 10000]
    proportions.sort(key=lambda x: x[1])
    proportion_count = len(proportions)
    sample_limit = int(round(proportion_count * .1))
    # Resolve ties, if present
    while (proportions[sample_limit - 1][1] == proportions[sample_limit][1]
                and sample_limit < proportion_count):
        sample_limit += 1
    bottom_10 = set([el[0] for el in proportions[:sample_limit]])
    sample_limit = int(round(proportion_count * .9))
    # Resolve ties, if present
    while (proportions[sample_limit - 1][1] == proportions[sample_limit][1]
                and sample_limit > 1):
        sample_limit -= 1
    top_10 = set([el[0] for el in proportions[sample_limit:]])

    # Get 500 commonest words
    word_incidences = defaultdict(int)
    bottom_10_word_incidences = defaultdict(int)
    top_10_word_incidences = defaultdict(int)
    with xopen(args.sra_metadata) as sra_stream:
        sra_stream.readline() # header
        total_sample_count, bottom_10_sample_count, top_10_sample_count = (
                0, 0, 0
            )
        for line in sra_stream:
            words = set([term.strip().strip(',.@#$%^&*!;:').lower()
                            for term in
                            re.findall(r"[\w']+", line)])
            total_sample_count += 1
            for word in words:
                word_incidences[word] += 1
            try:
                if srr_to_index[line.partition('\t')[0]] in bottom_10:
                    bottom_10_sample_count += 1
                    for word in words:
                        bottom_10_word_incidences[word] += 1
                elif srr_to_index[line.partition('\t')[0]] in top_10:
                    top_10_sample_count += 1
                    for word in words:
                        top_10_word_incidences[word] += 1
            except KeyError:
                continue

    for word_dict, descriptor, sample_count in [(word_incidences, 'all_words',
                                                    total_sample_count),
                                                    (bottom_10_word_incidences,
                                                        'bottom_10_pct',
                                                        bottom_10_sample_count),
                                                    (top_10_word_incidences,
                                                        'top_10_pct',
                                                        top_10_sample_count)]:
        word_dict = sorted(
                word_dict.items(), key=lambda x: x[1], reverse=True
            )[:500]
        with open(args.basename
                    + '.' + descriptor + '.tsv', 'w') as word_stream:
            print >>word_stream, 'word\tsample count (total samples: %d)' % (
                    sample_count
                )
            print >>word_stream, \
                '\n'.join([('%s\t%d' % el) for el in word_dict])

    print >>sys.stderr, 'Done counting word incidences in SRA metadata.'

    '''Aggregate junction stats: how many junctions/overlaps of given type
    are found in >= K samples/projects?'''
    sample_stats_to_aggregate = [sample_count_to_junction_count,
                                 sample_count_to_annotated_junction_count,
                                 sample_count_to_exonskip_junction_count,
                                 sample_count_to_altstartend_junction_count,
                                 sample_count_to_novel_junction_count,
                                 sample_count_to_GTAG_junction_count,
                                 sample_count_to_GTAG_ann_count,
                                 sample_count_to_GCAG_junction_count,
                                 sample_count_to_GCAG_ann_count,
                                 sample_count_to_ATAC_junction_count,
                                 sample_count_to_ATAC_ann_count]
    project_stats_to_aggregate = [project_count_to_junction_count,
                                  project_count_to_annotated_junction_count,
                                  project_count_to_exonskip_junction_count,
                                  project_count_to_altstartend_junction_count,
                                  project_count_to_novel_junction_count,
                                  project_count_to_GTAG_junction_count,
                                  project_count_to_GTAG_ann_count,
                                  project_count_to_GCAG_junction_count,
                                  project_count_to_GCAG_ann_count,
                                  project_count_to_ATAC_junction_count,
                                  project_count_to_ATAC_ann_count]
    for stats, descriptor in [(sample_stats_to_aggregate, 'sample'),
                                (project_stats_to_aggregate, 'project')]:
        max_count, min_count = 0, 1000000000 # way larger than max # samples
        for stat in stats:
            max_count = max(stat.keys() + [max_count])
            min_count = min(stat.keys() + [min_count])
        stat_count = len(stats)
        stat_aggregators = [0 for _ in xrange(stat_count)]
        with open(args.basename + '.' + descriptor + '.stats.tsv', 'w') \
            as stat_stream:
            print >>stat_stream, ('min {descriptor}s\t'
                                  'junctions\t'
                                  'annotated\t'
                                  'exonskips\t'
                                  'altstartend\t'
                                  'novel\t'
                                  'GTAG\t'
                                  'annotated GTAG\t'
                                  'GCAG\t'
                                  'annotated GCAG\t'
                                  'ATAC\t'
                                  'annotated ATAC').format(
                                            descriptor=descriptor
                                        )
            for descriptor_count in xrange(max_count, min_count - 1, -1):
                for i in xrange(stat_count):
                    stat_aggregators[i] += stats[i][descriptor_count]
                print >>stat_stream, '\t'.join(
                                    [str(descriptor_count)]
                                    + [str(el) for el in stat_aggregators]
                                )

    print >>sys.stderr, 'Dumped sample/project-level aggregate junction stats.'

    # Dump junction information by sample
    with open(args.basename + '.stats_by_sample.tsv', 'w') as stat_stream:
        print >>stat_stream, ('sample index\tproject\tsample\texperiment\trun'
                              '\tjunctions\tannotated_junctions'
                              '\tjunctions_geq_5\tannotated_junctions_geq_5'
                              '\toverlaps\tannotated_overlaps')
        for sample_index in sorted(index_to_sra.keys()):
            print >>stat_stream, '\t'.join(
                        [str(el) for el in 
                            [sample_index, index_to_sra[sample_index],
                                junction_counts[sample_index],
                                annotated_junction_counts[sample_index],
                                junction_counts_geq_5[sample_index],
                                annotated_junction_counts_geq_5[sample_index],
                                overlap_counts[sample_index],
                                annotated_overlap_counts[sample_index]]]
                    )
    print >>sys.stderr, 'Dumped junction info by sample.'

    # Dump overlap info for Venn diagram
    with open(args.basename + '.venn.txt', 'w') as venn_stream:
        in_all = set.intersection(
                        magic_junctions, rmake_junctions, subread_junctions
                    )
        in_one = set.union(
                        magic_junctions, rmake_junctions, subread_junctions
                    )
        in_two = set.union(
                    set.intersection(magic_junctions, rmake_junctions),
                    set.intersection(magic_junctions, subread_junctions),
                    set.intersection(rmake_junctions, subread_junctions)
                )
        print >>venn_stream, (
                'total samples studied by SEQC consortium and Rail: %d'
                    % len(seqc_indexes)
            )
        print >>venn_stream, (
                'junctions found by magic, rmake, and subread: %d'
                    % len(in_all)
            )
        print >>venn_stream, (
                'junctions found by magic, rmake, or subread: %d'
                    % len(in_one)
            )
        print >>venn_stream, (
                'junctions found by at least two of '
                '[magic, rmake, subread]: %d'
            ) % len(in_two)
        print >>venn_stream, (
                'junctions found by Rail: %d' % len(rail_seqc_junctions)
            )
        print >>venn_stream, (
                'junctions found by Rail in at least 5 SEQC samples: %d'
                    % len(rail_seqc_junctions_geq_5_samples)
            )
        print >>venn_stream, (
                'junctions found by Rail in at least 10 SEQC samples: %d'
                    % len(rail_seqc_junctions_geq_10_samples)
            )
        print >>venn_stream, (
                'junctions found by Rail in at least 15 SEQC samples: %d'
                    % len(rail_seqc_junctions_geq_15_samples)
            )
        print >>venn_stream, (
                'junctions found by Rail in at least 20 SEQC samples: %d'
                    % len(rail_seqc_junctions_geq_20_samples)
            )
        print >>venn_stream, (
                'junctions found by Rail in at least 25 SEQC samples: %d'
                    % len(rail_seqc_junctions_geq_25_samples)
            )
        print >>venn_stream, (
                'junctions found by Rail in at least 30 SEQC samples: %d'
                    % len(rail_seqc_junctions_geq_30_samples)
            )
        print >>venn_stream, (
                'junctions found by Rail and all of '
                '[magic, rmake, subread]: %d'
            ) % len(in_all.intersection(rail_seqc_junctions))
        print >>venn_stream, (
                'junctions found by Rail and none of '
                '[magic, rmake, subread]: %d'
            ) % len(rail_seqc_junctions - seqc_junctions)
        print >>venn_stream, (
                'junctions found by Rail in at least 5 SEQC samples '
                'and none of [magic, rmake, subread]: %d'
            ) % len(rail_seqc_junctions_geq_5_samples - seqc_junctions)
        print >>venn_stream, (
                'junctions found by Rail in at least 10 SEQC samples '
                'and none of [magic, rmake, subread]: %d'
            ) % len(rail_seqc_junctions_geq_10_samples - seqc_junctions)
        print >>venn_stream, (
                'junctions found by Rail in at least 15 SEQC samples '
                'and none of [magic, rmake, subread]: %d'
            ) % len(rail_seqc_junctions_geq_15_samples - seqc_junctions)
        print >>venn_stream, (
                'junctions found by Rail in at least 20 SEQC samples '
                'and none of [magic, rmake, subread]: %d'
            ) % len(rail_seqc_junctions_geq_20_samples - seqc_junctions)
        print >>venn_stream, (
                'junctions found by Rail in at least 25 SEQC samples '
                'and none of [magic, rmake, subread]: %d'
            ) % len(rail_seqc_junctions_geq_25_samples - seqc_junctions)
        print >>venn_stream, (
                'junctions found by Rail in at least 30 SEQC samples '
                'and none of [magic, rmake, subread]: %d'
            ) % len(rail_seqc_junctions_geq_30_samples - seqc_junctions)
        print >>venn_stream, (
                'junctions found by Rail and one of '
                '[magic, rmake, subread]: %d'
            ) % len(in_one.intersection(rail_seqc_junctions))
        print >>venn_stream, (
                'junctions found by Rail and at least two of '
                '[magic, rmake, subread]: %d'
            ) % len(in_two.intersection(rail_seqc_junctions))

    print >>sys.stderr, 'Dumped Venn info.'