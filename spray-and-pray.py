#!/usr/bin/env python3
from collections import defaultdict
import re
import os
import textwrap
import argparse
import sys
import statistics
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn import metrics
from sklearn.datasets import make_blobs
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram, linkage
from Bio.SeqUtils import gc_fraction
from Bio.Seq import Seq


def remove(stringOrlist, list):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            pass
    outString = "".join(emptyList)
    return outString


def reject_outliers(data, m):
    u = np.mean(data)
    s = np.std(data)
    filtered = [e for e in data if (u - m * s < e < u + m * s)]
    return filtered


def filter(list, items):
    outLS = []
    for i in list:
        if i not in items:
            outLS.append(i)
    return outLS


def delim(line):
    ls = []
    string = ''
    for i in line:
        if i != " ":
            string += i
        else:
            ls.append(string)
            string = ''
    ls.append(string)
    ls = filter(ls, [""])
    return ls


def allButTheLast(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(0, length-1):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x)-1]


def ribosome(seq):
    Dict = defaultdict(lambda: defaultdict(list))
    NTs = ['T', 'C', 'A', 'G']
    stopCodons = ['TAA', 'TAG', 'TGA']
    Codons = []
    for i in range(4):
        for j in range(4):
            for k in range(4):
                codon = NTs[i] + NTs[j] + NTs[k]
                # if not codon in stopCodons:
                Codons.append(codon)

    CodonTable = {}
    AAz = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    AAs = list(AAz)
    k = 0
    for base1 in NTs:
        for base2 in NTs:
            for base3 in NTs:
                codon = base1 + base2 + base3
                CodonTable[codon] = AAs[k]
                k += 1

    prot = []
    for j in range(0, len(seq), 3):
        codon = seq[j:j + 3]
        try:
            prot.append(CodonTable[codon])
        except KeyError:
            prot.append("X")
    protein = ("".join(prot))
    return protein


def codonTable(seq):
    Dict = defaultdict(lambda: defaultdict(list))
    NTs = ['T', 'C', 'A', 'G']
    stopCodons = ['TAA', 'TAG', 'TGA']
    Codons = []
    for i in range(4):
        for j in range(4):
            for k in range(4):
                codon = NTs[i] + NTs[j] + NTs[k]
                # if not codon in stopCodons:
                Codons.append(codon)

    CodonTable = {}
    AAz = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    AAs = list(AAz)
    k = 0
    for base1 in NTs:
        for base2 in NTs:
            for base3 in NTs:
                codon = base1 + base2 + base3
                CodonTable[codon] = AAs[k]
                k += 1

    prot = []
    for j in range(0, len(seq), 3):
        codon = seq[j:j + 3]
        try:
            Dict[CodonTable[codon]][codon].append(codon)
            prot.append(CodonTable[codon])
        except KeyError:
            prot.append("X")
    protein = ("".join(prot))
    return Dict


def tet(seq):
    Dict = defaultdict(list)
    NTs = ['T', 'C', 'A', 'G']
    for i in range(4):
        for j in range(4):
            for k in range(4):
                for l in range(4):
                    if NTs[i] in ["A", "G", "C", "T"] and NTs[j] in ["A", "G", "C", "T"] and NTs[k] in ["A", "G", "C", "T"] and NTs[l] in ["A", "G", "C", "T"]:
                        tet = NTs[i] + NTs[j] + NTs[k] + NTs[l]
                        Dict[tet] = []

    totalKmers = 0
    for m in range(len(seq)):
        TET = (seq[m:m+4])
        if len(TET) == 4:
            if TET[0] in ["A", "G", "C", "T"] and TET[1] in ["A", "G", "C", "T"] and TET[2] in ["A", "G", "C", "T"] and TET[3] in ["A", "G", "C", "T"]:
                Dict[TET].append(TET)
                totalKmers += 1

    return Dict, totalKmers


def Dictparser(Dictionary):
    lowest = float(1000)
    for i in Dictionary:
        if float(Dictionary[i]) < float(lowest):
            lowest = Dictionary[i]
            key = i
    return [i, lowest]


def SeqCoord(seq, start, end):
    return seq[start:end]


def howMany(ls, exclude):
    counter = 0
    for i in ls:
        if i != exclude:
            counter += 1
    return counter


def sum(ls):
    count = 0
    for i in ls:
        count += float(i)
    return count


def cluster(data, maxgap):
    '''Arrange data into groups where successive elements
       differ by no more than *maxgap*

        #->>> cluster([1, 6, 9, 100, 102, 105, 109, 134, 139], maxgap=10)
        [[1, 6, 9], [100, 102, 105, 109], [134, 139]]

        #->>> cluster([1, 6, 9, 99, 100, 102, 105, 134, 139, 141], maxgap=10)
        [[1, 6, 9], [99, 100, 102, 105], [134, 139, 141]]

    '''
    # data = sorted(data)
    data.sort(key=int)
    groups = [[data[0]]]
    for x in data[1:]:
        if abs(x - groups[-1][-1]) <= maxgap:
            groups[-1].append(x)
        else:
            groups.append([x])
    return groups


def lastItem(ls):
    x = ''
    for i in ls:
        if i != "":
            x = i
    return x


def replace(stringOrlist, list, item):
    emptyList = []
    for i in stringOrlist:
        if i not in list:
            emptyList.append(i)
        else:
            emptyList.append(item)
    outString = "".join(emptyList)
    return outString


def fasta(fasta_file):
    count = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            count += 1
            if count % 1000000 == 0:
                print(count)

            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                # header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                # header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    # print(count)
    return Dict


def fasta2(fasta_file):
    count = 0
    seq = ''
    header = ''
    Dict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    for i in fasta_file:
        i = i.rstrip()
        if re.match(r'^>', i):
            count += 1
            if count % 1000000 == 0:
                print(count)

            if len(seq) > 0:
                Dict[header] = seq
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
            else:
                header = i[1:]
                header = header.split(" ")[0]
                seq = ''
        else:
            seq += i
    Dict[header] = seq
    # print(count)
    return Dict


def allButTheFirst(iterable, delim):
    x = ''
    length = len(iterable.split(delim))
    for i in range(1, length):
        x += iterable.split(delim)[i]
        x += delim
    return x[0:len(x)]


parser = argparse.ArgumentParser(
    prog="spray-and-pray.py",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''
    ************************************************************************
    

    
    Developed by Arkadiy Garber; University of Montana, Biological Sciences
    Please send comments and inquiries to arkadiy.garber@mso.umt.edu
    ************************************************************************
    '''))

parser.add_argument('-g', type=str, help="Input bin/assembly in FASTA format", default="NA")

# parser.add_argument('-o', type=str, help="Input ORFs in FASTA amino acid format", default="NA")

parser.add_argument('-ref', type=str, help="Input reference protein database (recommended: nr). Could be FASTA file or "
                                           "DIAMOND database file (with extension .dmnd)", default="NA")

parser.add_argument('-bam', type=str, help="Input sorted BAM file with coverage info (optional)", default="NA")

parser.add_argument('-depth', type=str, help="Input depth/coverage information (optional). "
                                             "Must in the format that is produced from the "
                                             "jgi_summarize_bam_contig_depths script from the Metabat package", default="NA")

parser.add_argument('-out', type=str, help="Basename for output files", default="NA")

parser.add_argument('-dir', type=str, help="output directory", default="NA")

parser.add_argument('-lvl', type=str, help="Level of the taxonomic hierarchy to include in the summary file (Domain, Phylum, Class, Genus, species)", default="NA")

parser.add_argument('-t', type=int, help="number of threads to use for DIAMOND BLAST", default=1)

parser.add_argument('--makedb', type=str, help="if the DIAMOND database does not already exist "
                                                    "(i.e. file with extension .dmnd), and you would like the program t"
                                               "o run  diamond makedb, provide this flag", const=True, nargs="?")

parser.add_argument('--include_id', type=str, help="", const=True, nargs="?")

parser.add_argument('--bin', type=str, help="Including this flag will direct SprayNPray to perform hierarchical "
                                            "clustering based on 1) tetranucleotide frequency, 2) GC-content, 3) codon usage bias, "
                                            "and 4) read coverage (if BAM file is provided). SprayNPray will then split the input contigs into multiple FASTA files, "
                                            "each with its own summary file. This file is incompatible with the --fa flag.", const=True, nargs="?")

parser.add_argument('--spades', type=str, help="is this a SPAdes assembly, with the original SPAdes headers? If so, "
                                               "then you can provide this flag for coverage info.", const=True, nargs="?")

parser.add_argument('--megahit', type=str, help="is this a Megahit assembly, with the original Megahit headers? If so, "
                                               "then you can provide this flag for coverage info.", const=True, nargs="?")

parser.add_argument('--custom_ref', type=str, help="the reference database is not nr", const=True, nargs="?")

parser.add_argument('--meta', type=str, help="contigs are from a mixed community of organisms", const=True, nargs="?")

parser.add_argument('--c', type=str, help="Use complete genes only (no genes that run off edges of contigs)", const=True, nargs="?")

parser.add_argument('--hgt', type=str, help="provide this flag if you'd like the program to output potential HGTs into a separate file. "
                                            "This feature is designed for eukaryotic contigs expected to have HGTs of bacterial origin.", const=True, nargs="?")

parser.add_argument('--euk', type=str, help="provide this flag if you'd like the program to split away eukaryotic contigs from prokaryotic ones.", const=True, nargs="?")

parser.add_argument('--fa', type=str, help="write subset of contigs that match user-specified parameters to a separate FASTA file", const=True, nargs="?")

# parser.add_argument('--include_zero_hits', type=str, help="write subset of contigs that match user-specified parameters to a separate FASTA file", const=True, nargs="?")

parser.add_argument('-blast', type=str, help="DIAMOND BLAST output file from previous run", default="NA")

parser.add_argument('-hits', type=str, help="total number of DIAMOND hits to report in DIAMOND output file (default=100)", default="100")

parser.add_argument('-domain', type=str, help="domain expected among hits to provided contigs, to be written to FASTA file (e.g. Bacteria, Archaea, Eukaryota)", default="NA")

parser.add_argument('-phylum', type=str, help="phylum expected among hits to provided contigs, to be written to FASTA file (e.g. Proteobacteria). "
                                              "If you provide this name, please be sure to also provide the domain name via -domain", default="NA")

parser.add_argument('-Class', type=str, help="class name expected among hits to provided contigs, to be written to FASTA file (e.g. Gammaproteobacteria). "
                                             "If you provide this name, please be sure to also provide the domain and phylum names", default="NA")

parser.add_argument('-genus', type=str, help="genus name expected among hits to provided contigs, to be written to FASTA file (e.g. Shewanella). "
                                             "If you provide this name, please be sure to also provide the domain, phylum, and class names", default="NA")

parser.add_argument('-species', type=str, help="species name expected among hits to provided contigs, to be written to FASTA file (e.g. oneidensis, coli, etc.). "
                                               "If you provide this name, please be sure to also provide the domain, phylum, class, and genus names", default="NA")

parser.add_argument('--phage', type=str, help="add this flag if what you are interested in is phage contigs", const=True, nargs="?")

parser.add_argument('-perc', type=float, help="percentage of total hits to the contig that must be to the specified genus/species for writing to FASTA", default=0)

parser.add_argument('-gc', type=float, help="minimum GC-content of contigs to write to FASTA (default = 0)", default=0)

parser.add_argument('-GC', type=float, help="maximum GC-content of contigs to write to FASTA (default = 100)", default=100)

parser.add_argument('-cov', type=float, help="minimum coverage of contigs to write to FASTA (default = 0)", default=0)

parser.add_argument('-COV', type=float, help="maximum coverage of contigs to write to FASTA (default = 100000000)", default=100000000)

parser.add_argument('-cd', type=float, help="minimum coding density (in hits/kb) to write to FASTA (default = 0)", default=0)

parser.add_argument('-CD', type=float, help="maximum coding density (in hits/kb) to write to FASTA (default = 5)", default=5)

parser.add_argument('-l', type=float, help="minimum length of contig to write to FASTA (default = 300)", default=300)

parser.add_argument('-L', type=float, help="maximum length of contig to write to FASTA (default = 100000000)", default=100000000)

parser.add_argument('-aai', type=float, help="minimum average amino acid identity (percent) to reference proteins (default = 30)", default=30)

parser.add_argument('-minGenes', type=float, help="minimum number of genes that must be on a contig to write to FASTA (default = 1)", default=1)

parser.add_argument('-minLength', type=float, help="minimum length of gene to include in the BLAST analysis (default = 90)", default=90)

parser.add_argument('--test', type=str, help="add this flag during testing of this program's dependencies", const=True, nargs="?")

parser.add_argument('--debug', type=str, help="add this flag to debug an issue caused by a dependency", const=True, nargs="?")

# parser.add_argument('-key', type=str, help="Path to the taxmap_slv_ssu_ref_nr_138.1.txt file, which should be in the repository containing this program", default="NA")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_known_args()[0]


print(".")
# checking paramters:
if args.ref == "NA":
    print("Please provide a reference file via -ref")
    raise SystemExit
else:
    print("Reference database: " + args.ref)

if args.g == "NA":
    print("Please provide an input genome via -g")
    raise SystemExit
else:
    print("Input genome: " + args.g)
    genome = args.g
    if lastItem(genome.split(".")) in ["faa", "ffn"]:
        print("From the file extension of input file, it looks like you provided proteins or gene sequences. "
              "Currently, SprayNPray only takes contigs.")
        answer = input("Did you provide contigs and would like to proceed with analysis? (y/n): ")
        if answer == "y":
            pass
        else:
            raise SystemExit

if args.lvl != "NA":
    if args.lvl in ["Domain", "Phylum", "Class", "Genus", "species"]:
        print("Taxonomic level included in the summary file: " + args.lvl)
    else:
        print("Invalud taxonomic level provided via -lvl")
        raise SystemExit

if args.blast != "NA":
    print("Provided BLAST output file: " + args.blast)

os.system("which spray-and-pray.py > mainDir.txt")

file = open("mainDir.txt")
location = os.getcwd()
for i in file:
    location = i.rstrip()
location = allButTheLast(location, "/")

silvaFile = "/home/ark/MAB/bin/SprayNPray/taxmap_slv_ssu_ref_nr_138.1.txt"
os.system("rm mainDir.txt")

if args.out == "NA" or args.dir == "NA":
    print("Provide output name!")
    raise SystemExit
else:
    outfilename = args.out
    outdir = args.dir

if args.fa:
    print("SprayNPray will write a FASTA file with contigs matching user-specified metrics: " + outfilename + "-contigs.fa")
    print("SprayNPray will write a FASTA file with contigs not matching user-specified metrics: " + outfilename + "-unmatched.contigs.fa\n")
    if args.species != "NA":
        if "NA" in [args.genus, args.Class, args.phylum, args.domain]:
            print("If species name is provided, please provide also the Genus, Class, Phylum, and Domain names")
            raise SystemExit
        else:
            print("species restriction: " + args.species)

    if args.genus != "NA":
        if "NA" in [args.Class, args.phylum, args.domain]:
            print("If Genus name is provided, please provide also the Class, phylum, and domain names")
            raise SystemExit
        else:
            print("Genus restriction: " + args.genus)

    if args.Class != "NA":
        if "NA" in [args.phylum, args.domain]:
            print("If Class name is provided, please provide also the phylum and domain names")
            raise SystemExit
        else:
            print("Class restriction: " + args.Class)

    if args.phylum != "NA":
        if "NA" in [args.domain]:
            print("If phylum name is provided, please provide also the domain name")
            raise SystemExit
        else:
            print("Phylum restriction: " + args.phylum)

    if args.domain:
        print("Domain restriction: " + args.domain)

file = open(args.g)
file = fasta(file)
total = 0
for i in file.keys():
    total += len(file[i])

if total < 20000:
    if args.meta:
        pass
    else:
        print("Looks like there are less than 20000 characters in your provided sequences file. Please re-run the script with the --meta flag")
        raise SystemExit

if len(file.keys()) == 0:
    print("SprayNPray did not detect any sequences in your provided contigs file. Please check your input FASTA file.")
    raise SystemExit

if args.blast == "NA":

    print("Running Prodigal: calling ORFs from provided contigs")

    if args.debug:
        if args.c:

            if args.meta:
                os.system("prodigal -i %s -a %s-proteins.faa -d %s-cds.ffn -p meta -c" % (
                args.g, args.g, args.g))
            else:
                os.system(
                    "prodigal -i %s -a %s-proteins.faa -d %s-cds.ffn -c" % (args.g, args.g, args.g))
        else:
            if args.meta:
                os.system("prodigal -i %s -a %s-proteins.faa -d %s-cds.ffn -p meta" % (
                args.g, args.g, args.g))
            else:
                os.system("prodigal -i %s -a %s-proteins.faa -d %s-cds.ffn" % (args.g, args.g, args.g))
    else:
        if args.c:

            if args.meta:
                print("prodigal -i %s -a %s-proteins.faa -d %s-cds.ffn -p meta -c" % (args.g, args.g, args.g))
                os.system("prodigal -i %s -a %s-proteins.faa -d %s-cds.ffn -p meta -c -q > /dev/null 2>&1" % (args.g, args.g, args.g))
            else:
                os.system("prodigal -i %s -a %s-proteins.faa -d %s-cds.ffn -c -q > /dev/null 2>&1" % (args.g, args.g, args.g))
        else:
            if args.meta:
                os.system("prodigal -i %s -a %s-proteins.faa -d %s-cds.ffn -p meta -q > /dev/null 2>&1" % (args.g, args.g, args.g))
            else:
                os.system("prodigal -i %s -a %s-proteins.faa -d %s-cds.ffn -q > /dev/null 2>&1" % (args.g, args.g, args.g))

    # checking the prodigal-produced .faa file that will be used for downstream analysis
    faa = open("%s-proteins.faa" % args.g)
    faa = fasta(faa)
    out = open("%s-proteins-new.faa" % args.g, "w")
    count = 0
    residues = 0
    Xs = 0
    for i in faa.keys():
        print(i)
        print(faa[i])
        if len(faa[i]) >= args.minLength:
            out.write(">" + i + "\n")
            out.write(faa[i] + "\n")
            count += 1
            residues += len(faa[i])
            # residues = "1000"
            Xs += list(faa[i]).count("X")

    if Xs == residues:
        print(Xs)
        print(residues)
        print("")
        print("Looks like you have provided protein sequences instead of contigs to SprayNPray. "
              "Currently, this program only accepts as input a FASTA-formatted file of contigs, "
              "containing nucleotide sequences (e.g. files with extensions fna, fa, or fasta)."
              "If you believe this message is in error, please let the developers know by starting a GitHub Issue.")
        raise SystemExit

    if count == 0:
        print("SprayNPray did not detect any protein sequences in the Prodigal output. Looks like Prodigal did not "
              "perform as expected. Please check your input contigs, as well the Prodigal installation within the SprayNPray "
              "conda environment. You can type \'prodigal -h\' to check if prodigal is installed properly.")
        raise SystemExit

    out.close()
    os.system("mv %s-proteins-new.faa %s-proteins.faa" % (args.g, args.g))

    db = args.ref

    print("Running Diamond BLAST")
    print(
        "diamond blastp --db %s --query %s-proteins.faa "
        "--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle "
        "--out %s.blast --max-target-seqs %s --evalue 1E-15 --threads %s --query-cover 50 --subject-cover 50"
        % (db, args.g, args.g, args.hits, args.t))
    os.system(
        "diamond blastp --db %s --query %s-proteins.faa "
        "--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle "
        "--out %s.blast --max-target-seqs %s --evalue 1E-15 --threads %s --query-cover 50 --subject-cover 50"
        % (db, args.g, args.g, args.hits, args.t))

    blastFile = "%s.blast" % (args.g)
else:
    blastFile = args.blast

blast = open(blastFile)
count = 0
for j in blast:
    count += 1

if count == 0:
    print("It looks like there was an issue running DIAMOND. The file that was created did not have any hits. "
          "Please check your input file, as well as the DIAMOND installation")
    raise SystemExit

try:

    blast = open(blastFile)
except FileNotFoundError:
    print("no blast file")
    raise SystemExit

out = open("%s/%s-top%s.csv" % (outdir, outfilename, args.hits), "w")
out.write("orf,taxa,top_hit\n")
for i in blast:
    ls = i.rstrip().split("\t")
    try:
        name = (ls[12].split("[")[1])
        name = name[0:len(name) - 1]
        out.write(ls[0] + "," + replace(name, [","], ";") + "," + replace(ls[12], [","], ";") + "\n")
    except IndexError:
        pass
out.close()
os.system("mv %s-top%s.csv %s/" % (outfilename, args.hits, outdir))

if args.bam != "NA":
    print("Extracting coverage information from the provided BAM files")
    if args.debug:
        os.system("jgi_summarize_bam_contig_depths --outputDepth %s.depth %s" % (args.g, args.bam))
    else:
        os.system("jgi_summarize_bam_contig_depths --outputDepth %s.depth %s > /dev/null 2>&1" % (args.g, args.bam))

print("Calculating GC-content")
gcDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
total = 0
for i in file.keys():
    seq = file[i]
    gcDict[i] = str(gc_fraction(seq))

print("Preparing summary: %s.csv" % outfilename)

aaiDict = defaultdict(list)
blastDict = defaultdict(list)
redunDict = defaultdict(list)
blast = open(blastFile)
for i in blast:
    ls = i.rstrip().split("\t")
    contig = allButTheLast(ls[0], "_")
    name = ls[12]
    try:
        name = name.split("]")[0]
        name = name.split("[")[1]
    except IndexError:
        if args.custom_ref:
            name = ls[1]
        else:
            name = "NA"
    aai = ls[2]
    if ls[0] not in redunDict.keys():
        redunDict[ls[0]].append(name)

        if args.include_id:
            blastDict[contig].append(name + " " + ls[2])
        else:
            blastDict[contig].append(name)
        aaiDict[contig].append(float(aai))
blast.close()

if args.bam != "NA":
    depthDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    depth = open("%s.depth" % args.g)
    for i in depth:
        ls = i.rstrip().split("\t")
        if ls[1] != "contigLen":
            depthDict[ls[0]]["length"] = float(ls[1])
            depthDict[ls[0]]["depth"] = ls[2]
elif args.depth != "NA":
    depthDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    depth = open("%s" % args.depth)
    for i in depth:
        ls = i.rstrip().split("\t")
        if ls[1] != "contigLen":
            depthDict[ls[0]]["length"] = float(ls[1])
            depthDict[ls[0]]["depth"] = ls[2]

# reading silva headers
silvaDict = defaultdict(lambda: defaultdict(lambda: 'unclassified'))
try:
    silva = open(silvaFile)
except FileNotFoundError:
    print(silvaFile)
    os.system("gunzip %s.gz" % silvaFile)
    silva = open(silvaFile)
    # print("SprayNPray cannot find the following file: taxmap_slv_ssu_ref_nr_138.1.txt. \n"
    #       "There is a good chance that it is present in its gzipped form in the SprayNPray \n"
    #       "directory/folder on your system. Please unzip this file and try running the program \n"
    #       "again. If you just waited for a length DIAMOND run to finish, you can provide the \n"
    #       "DIAMOND BLAST output (%s.blast) to the command when you re-run using the -blast argument" % args.g)
    # raise SystemExit

for i in silva:
    ls = i.rstrip().split("\t")
    if ls[0] != "primaryAccession":

        Domain = (ls[3].split(";")[0])
        Phylum = (ls[3].split(";")[1])
        if Phylum == "":
            Phylum = Domain

        try:
            Class = (ls[3].split(";")[2])
        except IndexError:
            Class = (ls[3].split(";")[0])

        if Domain in ["Bacteria", "Archaea"]:
            Genus = lastItem(ls[3].split(";"))
            if Genus.split(" ") == "Candidatus":
                Genus = Genus.split(" ")[1]

            if Genus == "Oikopleura":
                Domain = "Eukaryota"
            if Genus == "Diplosphaera":
                Domain = "Eukaryota"
            if Genus == "Planococcus":
                Domain = "Eukaryota"

        elif Domain in ["Eukaryota"]:
            Genus = ls[4].split(" ")[0]
            if Genus == "Candidatus":
                Genus = ls[4].split(" ")[1]

            if lastItem(ls[3].split(";")) == "Chloroplast":
                Genus = "Chloroplast" + "_" + ls[4].split(" ")[0]

            if lastItem(ls[3].split(";")) == "Mitochondria":
                Genus = "Mitochondria" + "_" + ls[4].split(" ")[0]

            if Genus == "uncultured":
                Genus = "uncultured_" + lastItem(ls[3].split(";"))

            if Genus == "Labrys":
                Domain = "Bacteria"
            if Genus == "Halofilum":
                Domain = "Bacteria"
            if Genus == "Bacillus":
                Domain = "Bacteria"
            if Genus == "Lactobacillus":
                Domain = "Bacteria"
            if Genus == "Pseudomonas":
                Domain = "Bacteria"
            if Genus == "Arthrobacter":
                Domain = "Bacteria"
            if Genus == "Paracoccus":
                Domain = "Bacteria"
            if Genus == "Ensifer":
                Domain = "Bacteria"
            if Genus == "Arthrobacter":
                Domain = "Bacteria"
            if Genus == "Aeromonas":
                Domain = "Bacteria"
            if Genus == "Acinetobacter":
                Domain = "Bacteria"
            if Genus == "Edwardsiella":
                Domain = "Bacteria"
            if Genus == "Mesorhizobium":
                Domain = "Bacteria"
            if Genus == "Kitasatospora":
                Domain = "Bacteria"
            if Genus == "Clostridium":
                Domain = "Bacteria"
            if Genus == "Rhodocista":
                Domain = "Bacteria"
            if Genus == "Actinomyces":
                Domain = "Bacteria"

        silvaDict[Genus]["Domain"] = Domain
        silvaDict[Genus]["Phylum"] = Phylum
        silvaDict[Genus]["Class"] = Class
        silvaDict[Genus]["Genus"] = Genus

        silvaDict[Class]["Phylum"] = Phylum
        silvaDict[Class]["Domain"] = Domain
        silvaDict[Class]["Class"] = Class

        silvaDict[Phylum]["Phylum"] = Phylum
        silvaDict[Phylum]["Domain"] = Domain

        silvaDict[Domain]["Domain"] = Domain

out = open(outdir + "/" + outfilename + ".csv", "w")
out.write("contig" + "," + "contig_length" + "," + "hits_per_kb" + "," + "cov" + "," + "GC-content" + "," + "Average_AAI" + "," + "closest_blast_hits" + "\n")
for i in file.keys():
    if args.bam != "NA":
        depth = depthDict[i]["depth"]
        length = depthDict[i]["length"]
    elif args.depth != "NA":
        depth = depthDict[i]["depth"]
        length = depthDict[i]["length"]
    elif args.spades:
        depth = lastItem(i.split("_"))
        length = len(file[i])
    elif args.megahit:
        depth = i.split(" ")[2].split("=")[1]
        length = len(file[i])
    else:
        depth = "Unknown"
        length = len(file[i])
    gc = gcDict[i]
    hitsList = blastDict[i.split(" ")[0]]
    if len(hitsList) > 0:
        print(i.split(" ")[0])
        print(hitsList)
        print("")
    try:
        AAI = statistics.mean(aaiDict[i.split(" ")[0]])
    except statistics.StatisticsError:
        AAI = "NA"

    if len(hitsList) > 0:
        out.write(i + "," + str(length) + "," + str(len(hitsList) / (length / 1000)) + "," + str(depth) + "," + str(gc) + "," + str(AAI) + ",")

        if args.lvl == "NA":
            for j in hitsList:
                try:
                    out.write(j + "; ")
                except TypeError:
                    pass
        else:
            for j in hitsList:
                try:
                    Genus = j.split(" ")[0]
                    if Genus in ["Candidatus", "uncultured"]:
                        Genus = j.split(" ")[1]

                    try:
                        species = j.split(" ")[1]
                        if species == "sp.":
                            species = j.split(" ")[2]
                    except IndexError:
                        species = "unclassified"

                    Domain = silvaDict[Genus]["Domain"]
                    Phylum = silvaDict[Genus]["Phylum"]
                    Class = silvaDict[Genus]["Class"]

                    if len(silvaDict[Genus]) == 0:
                        Domain = "unclassifed"
                        Phylum = "unclassifed"
                        Class = "unclassifed"

                    if re.findall(r'symbiont', j):
                        Domain = j
                        Phylum = j
                        Class = j

                    if args.lvl == "Domain":
                        out.write(str(Domain) + "; ")

                    elif args.lvl == "Phylum":
                        out.write(str(Phylum) + "; ")

                    elif args.lvl == "Class":
                        out.write(str(Class) + "; ")

                    elif args.lvl == "Genus":
                        out.write(str(Genus) + "; ")

                    elif args.lvl == "species":
                        out.write(str(species) + "; ")

                    else:
                        break

                except TypeError:
                    pass

        out.write("\n")
out.close()
os.system("mv %s.csv %s" % (outfilename, outdir))
# os.system("mv %s-proteins.faa %s/" % (args.g, outdir))
# os.system("mv %s-cds.ffn %s/" % (args.g, outdir))
# if args.blast == "NA":
    # os.system("mv %s.blast %s/" % (args.g, outdir))
    # blastFile = "%s/%s.blast" % (outdir, args.g)


############## WORDCLOUD LOOP ####################
out = open(outdir + "/" + outfilename + ".words.csv", "w")
summary = open("%s/%s.csv" % (outdir, outfilename))
for i in summary:
    ls = i.rstrip().split(",")
    if ls[6] != "closest_blast_hits":
        ls2 = (ls[6].split("; "))
        contig = ls[0]
        for j in ls2:
            if (j.split(" ")[0]) == "Candidatus":
                word = (j.split(" ")[1])
            else:
                word = (j.split(" ")[0])

            if args.bam != "NA":
                cov = int(round(float(ls[3])))

                for k in range(0, cov):
                    if not re.findall(r'unclass', word):
                        out.write(word + "\n")

            else:
                if not re.findall(r'unclass', word):
                    out.write(word + "\n")
out.close()
os.system("mv %s.words.csv %s/" % (outfilename, outdir))

os.system("echo ${rscripts} > r.txt")
Rfile = open("r.txt")
for i in Rfile:
    Rdir = (i.rstrip())
os.system("rm r.txt")

try:
    test = open(Rdir + "/wordcloud.R")
except:
    os.system("which spray-and-pray.py > r.txt")
    Rfile = open("r.txt")
    for i in Rfile:
        Rdir = (i.rstrip())
    Rdir = allButTheLast(Rdir, "/")
    os.system("rm r.txt")

if args.debug:
    os.system("Rscript --vanilla %s/wordcloud.R %s/%s.words.csv %s/%s.words.tiff" % (Rdir, outdir, outfilename, outdir, outfilename))
else:
    os.system("Rscript --vanilla %s/wordcloud.R %s/%s.words.csv %s/%s.words.tiff > /dev/null 2>&1" % (Rdir, outdir, outfilename, outdir, outfilename))
#################################

if args.fa:
    print("Writing contigs based on user-specified metrics")
    summary = open("%s/%s.csv" % (outdir, outfilename))
    out = open('%s/%s-contigs.fa' % (outdir, outfilename), "w")
    out2 = open('%s/%s-unmatched.contigs.fa' % (outdir, outfilename), "w")

    for i in summary:
        ls = i.rstrip().split(",")

        if ls[1] != "contig_length":
            length = float(ls[1])
            hitsperkb = float(ls[2])
            gc = float(ls[4])
            hits = ls[6].split("; ")
            totalHits = len(hits)

            try:
                aai = float(ls[5])
            except ValueError:
                aai = 100

            if ls[3] != "Unknown":
                cov = float(ls[3])
            else:
                cov = 0

            if args.domain != "NA":
                # doing the math
                matches = 0

                for j in hits:
                    Genus = j.split(" ")[0]

                    try:
                        species = j.split(" ")[1]
                        if species == "sp.":
                            species = j.split(" ")[2]
                    except IndexError:
                        species = "unclassified"

                    Domain = silvaDict[Genus]["Domain"]
                    Phylum = silvaDict[Genus]["Phylum"]
                    Class = silvaDict[Genus]["Class"]
                    if len(silvaDict[Genus]) > 0:
                        domain = "unclassifed"
                        phylum = "unclassifed"
                        Class = "unclassifed"

                    if args.domain != "NA":

                        if args.phylum != "NA":

                            if args.Class != "NA":

                                if args.genus != "NA":

                                    if args.species != "NA":

                                        if species == args.species:
                                            genusChoices = args.genus

                                            if Genus in genusChoices.split(","):

                                                if args.phage:

                                                    if j.split(" ")[1] == "phage":
                                                        matches += 1
                                                else:
                                                    if j.split(" ")[1] == "phage":
                                                        pass
                                                    else:
                                                        matches += 1
                                    else:
                                        genusChoices = args.genus
                                        if Genus in genusChoices.split(","):
                                            matches += 1
                                else:
                                    if Class == args.Class:
                                        matches += 1
                            else:
                                if Phylum == args.phylum:
                                    matches += 1
                        else:
                            if Domain == args.domain:
                                matches += 1
                    else:
                        if args.phage:
                            if j.split(" ")[1] == "phage":
                                matches += 1
                            else:
                                matches += 1

                perc = (matches / totalHits) * 100

                if hits[0] == '' and totalHits == 1:
                    if args.minGenes == 0:
                        perc = 100

            else:
                perc = 100

            # if args.euk:
            #     if gc >= 0 and gc <= 40 and hitsperkb >= 0 and hitsperkb <= 0.5 and totalHits >= args.minGenes:
            #         out.write(">" + ls[0] + "\n")
            #         out.write(file[ls[0]] + "\n")
            #     else:
            #         out2.write(">" + ls[0] + "\n")
            #         out2.write(file[ls[0]] + "\n")
            #
            # else:

            if perc >= args.perc and gc >= args.gc and gc <= args.GC and length >= args.l and length <= args.L and \
                            cov >= args.cov and cov <= args.COV and aai >= args.aai and hitsperkb >= args.cd and \
                            hitsperkb <= args.CD and totalHits >= args.minGenes:
                out.write(">" + ls[0] + "\n")
                out.write(file[ls[0]] + "\n")
            else:
                out2.write(">" + ls[0] + "\n")
                out2.write(file[ls[0]] + "\n")

    out.close()


if args.euk:
    print("Writing contigs based on user-specified metrics")
    summary = open("%s/%s.csv" % (outdir, outfilename))
    out = open('%s/%s-eukaryotic.contigs.fa' % (outdir, outfilename), "w")
    out2 = open('%s/%s-prokaryotic.contigs.fa' % (outdir, outfilename), "w")

    for i in summary:
        ls = i.rstrip().split(",")

        if ls[1] != "contig_length":
            length = float(ls[1])
            hitsperkb = float(ls[2])
            gc = float(ls[4])
            hits = ls[6].split("; ")
            totalHits = len(hits)

            try:
                aai = float(ls[5])
            except ValueError:
                aai = 100

            if ls[3] != "Unknown":
                cov = float(ls[3])
            else:
                cov = 0

            if args.domain != "NA":
                # doing the math
                matches = 0

                for j in hits:
                    Genus = j.split(" ")[0]

                    try:
                        species = j.split(" ")[1]
                        if species == "sp.":
                            species = j.split(" ")[2]
                    except IndexError:
                        species = "unclassified"

                    Domain = silvaDict[Genus]["Domain"]
                    Phylum = silvaDict[Genus]["Phylum"]
                    Class = silvaDict[Genus]["Class"]
                    if len(silvaDict[Genus]) > 0:
                        domain = "unclassifed"
                        phylum = "unclassifed"
                        Class = "unclassifed"

                    if args.domain != "NA":

                        if args.phylum != "NA":

                            if args.Class != "NA":

                                if args.genus != "NA":

                                    if args.species != "NA":

                                        if species == args.species:
                                            genusChoices = args.genus

                                            if Genus in genusChoices.split(","):

                                                if args.phage:

                                                    if j.split(" ")[1] == "phage":
                                                        matches += 1
                                                else:
                                                    if j.split(" ")[1] == "phage":
                                                        pass
                                                    else:
                                                        matches += 1
                                    else:
                                        genusChoices = args.genus
                                        if Genus in genusChoices.split(","):
                                            matches += 1
                                else:
                                    if Class == args.Class:
                                        matches += 1
                            else:
                                if Phylum == args.phylum:
                                    matches += 1
                        else:
                            if Domain == args.domain:
                                matches += 1
                    else:
                        if args.phage:
                            if j.split(" ")[1] == "phage":
                                matches += 1
                            else:
                                matches += 1

                perc = (matches / totalHits) * 100

                if hits[0] == '' and totalHits == 1:
                    if args.minGenes == 0:
                        perc = 100

            else:
                perc = 100

            if args.euk:
                if gc >= 0 and gc <= 40 and hitsperkb >= 0 and hitsperkb <= 0.5 and totalHits >= args.minGenes:
                    out.write(">" + ls[0] + "\n")
                    out.write(file[ls[0]] + "\n")
                else:
                    out2.write(">" + ls[0] + "\n")
                    out2.write(file[ls[0]] + "\n")

            else:

                if perc >= args.perc and gc >= args.gc and gc <= args.GC and length >= args.l and length <= args.L and \
                                cov >= args.cov and cov <= args.COV and aai >= args.aai and hitsperkb >= args.cd and \
                                hitsperkb <= args.CD and totalHits >= args.minGenes:
                    out.write(">" + ls[0] + "\n")
                    out.write(file[ls[0]] + "\n")
                else:
                    out2.write(">" + ls[0] + "\n")
                    out2.write(file[ls[0]] + "\n")

    out.close()

if args.hgt:
    print("Looking for HGT candidates")
    summaryDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    summary = open("%s/%s.csv" % (outdir, outfilename))
    for i in summary:
        ls = i.rstrip().split(",")
        if ls[2] != "hits_per_kb":
            counter = 0
            for j in ls[6].split("; "):
                if re.findall(r'Bacteria', j):
                    counter += 1

            if counter > -1:
                if float(ls[2]) < 0.25:
                    summaryDict[ls[0]] = ls

    prots = open("%s-proteins.faa" % args.g)
    prots = fasta2(prots)

    outSeq = open("%s/%s.hgt.fasta" % (outdir, outfilename), "w")
    out = open("%s/%s.hgt.csv" % (outdir, outfilename), "w")
    out.write("contig,orf,blastHit,seq\n")
    blastDict = defaultdict(lambda: defaultdict(list))
    blast = open(blastFile)
    for i in blast:
        ls = i.rstrip().split("\t")
        contig = allButTheLast(ls[0], "_")
        if contig in summaryDict.keys():
            try:
                name = (ls[12].split("[")[1])
                name = name[0:len(name) - 1]
                Genus = (name.split(" ")[0])
                Domain = (silvaDict[Genus]["Domain"])
                blastDict[ls[0]][Domain].append(replace(ls[12], [","], ";"))

                # if Domain == "Bacteria":
                #     outSeq.write(">" + ls[0] + "\n")
                #     outSeq.write(prots[ls[0]] + "\n")
                #     out.write(
                #         contig + "," + ls[0] + "," + str(replace(ls[12], [","], ";")) + "," + prots[ls[0]] + "\n")

            except IndexError:
                pass

    for i in blastDict.keys():
        Bacteria = blastDict[i]["Bacteria"]
        Eukaryota = blastDict[i]["Eukaryota"]
        if len(Bacteria) > len(Eukaryota):
            contig = allButTheLast(i, "_")
            outSeq.write(">" + i + "\n")
            outSeq.write(prots[i] + "\n")
            out.write(contig + "," + i + "," + str(Bacteria[0]) + "," + prots[i] + "\n")

    outSeq.close()
    out.close()

if args.bin:
    print("Starting the binning algorithm")

    # os.system("echo ${hmms} > wd.txt")
    # WD = open("wd.txt")
    # wd = ''
    # for i in WD:
    #     wd = (i.rstrip())
    # os.system("rm wd.txt")
    #
    # try:
    #     test = open(wd + "/Universal_Hug_et_al.hmm")
    #     hmms = wd + "/Universal_Hug_et_al.hmm"
    # except FileNotFoundError:

    os.system("which spray-and-pray.py > wd.txt")

    WD = open("wd.txt")
    wd = ''

    for i in WD:
        wd = (i.rstrip())
    wd = allButTheLast(wd, "/")

    hmms = wd + "/Universal_Hug_et_al.hmm"
    os.system("rm wd.txt")

    cds = open("%s-cds.ffn" % (args.g))
    cds = fasta2(cds)

    print("Estimating the number of genomes")
    if args.debug:
        os.system(
            "hmmsearch --tblout universal.tblout --cpu %s %s %s-proteins.faa" % (args.t, hmms, args.g))
    else:
        os.system("hmmsearch --tblout universal.tblout --cpu %s %s %s-proteins.faa > /dev/null 2>&1" % (args.t, hmms, args.g))
    tblout = open("universal.tblout")
    tbloutDict = defaultdict(list)
    for i in tblout:
        if not re.match(r'#', i):
            ls = delim(i)
            tbloutDict[ls[2]].append(ls[0])

    hits = []
    for i in tbloutDict.keys():
        hits.append(len(tbloutDict[i]))

    filteredHits = reject_outliers(hits, 2)
    numClusters = round(statistics.mean(filteredHits))
    print("Predicting " + str(numClusters) + " clusters")

    # ************************************************
    print("Calculating tetranucleotide frequency")
    gcDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    tetDict2 = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    GC = 0
    total = 0
    count = 0
    for i in file.keys():
        count += 1
        perc = count/len(file.keys())
        perc = perc*100
        sys.stdout.write("calculating: %d%%   \r" % (perc))
        sys.stdout.flush()
        seq = file[i]
        TET = tet(seq)
        tetDict = TET[0]
        totalKmers = TET[1]
        for j in tetDict.keys():
            tetDict2[i][j] = len(tetDict[j]) / totalKmers

    # ************************************************
    print("Extracting GC and gene density from summary file")
    summaryDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    summary = open("%s/%s.csv" % (outdir, outfilename))
    for i in summary:
        ls = i.rstrip().split(",")
        if ls[0] != "contig":
            summaryDict[ls[0]]["den"] = float(ls[2])
            summaryDict[ls[0]]["gc"] = float(ls[4])

    # ************************************************
    print("Calculating codon usage bias")
    codonDict2 = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for i in cds.keys():
        contig = allButTheLast(i, "_")
        seq = (cds[i])
        Dict = (codonTable(seq))
        for j in Dict.keys():
            for k in Dict[j]:
                codonDict2[contig][j][k].append(len(Dict[j][k]))

    CodonList = []
    codonDict3 = defaultdict(lambda: defaultdict(lambda: 0))
    for i in sorted(codonDict2.keys()):
        contig = i
        for j in sorted(codonDict2[i]):
            localtotal = 0
            for k in sorted(codonDict2[i][j]):
                Codons = sum(codonDict2[i][j][k])
                localtotal += Codons

            for k in sorted(codonDict2[i][j]):
                Codons = sum(codonDict2[i][j][k])
                codonDict3[i][k] = Codons / localtotal
                if k not in CodonList:
                    CodonList.append(k)

    # ************************************************
    print("Consolidating")
    for i in tetDict2.keys():
        tetDict2[i]["den"] = summaryDict[i]["den"]
        tetDict2[i]["gc"] = summaryDict[i]["gc"]
        for j in CodonList:
            tetDict2[i][j] = codonDict3[i][j]
        if args.bam != "NA":
            tetDict2[i]["depth1"] = depthDict[i]["depth"]
            tetDict2[i]["depth2"] = depthDict[i]["depth"]
            tetDict2[i]["depth3"] = depthDict[i]["depth"]
            tetDict2[i]["depth4"] = depthDict[i]["depth"]
            tetDict2[i]["depth5"] = depthDict[i]["depth"]
            tetDict2[i]["depth6"] = depthDict[i]["depth"]
            tetDict2[i]["depth7"] = depthDict[i]["depth"]
            tetDict2[i]["depth8"] = depthDict[i]["depth"]
            tetDict2[i]["depth9"] = depthDict[i]["depth"]
            tetDict2[i]["depth10"] = depthDict[i]["depth"]

    print("Transforming data")
    tetDict3 = defaultdict(list)
    for i in tetDict2.keys():
        for j in tetDict2[i]:
            tetDict3[i].append(tetDict2[i][j])

    listOfTet = []
    listOfContigs = []
    for i in sorted(tetDict3.keys()):
        listOfContigs.append(i)
        listOfTet.append(tetDict3[i])

    X = np.array(listOfTet, dtype=float)

    linked = linkage(X, 'single')

    labelList = range(1, len(file.keys()) + 1)

    cluster = AgglomerativeClustering(n_clusters=numClusters, affinity='euclidean', linkage='ward')
    cluster.fit_predict(X)

    contigClusterDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    clusterDict = defaultdict(list)
    for i in range(0, len(listOfContigs)):
        clusterDict[cluster.labels_[i]].append(listOfContigs[i])
        contigClusterDict[listOfContigs[i]] = cluster.labels_[i]

    header = ''
    summaryDict = defaultdict(lambda: defaultdict(lambda: 'EMPTY'))
    summary = open("%s/%s.csv" % (outdir, outfilename))
    out = open("%s/%s-wClusters.csv" % (outdir, outfilename), "w")
    for i in summary:
        ls = i.rstrip().split(",")
        if ls[1] != "contig_length":
            summaryDict[ls[0]] = i.rstrip()
            out.write(str(contigClusterDict[ls[0]]) + "," + i.rstrip() + '\n')

        else:
            header = i.rstrip()
            out.write("cluster," + i.rstrip() + "\n")
    out.close()

    os.system("mkdir -p %s/bins" % (outdir))
    print("Writing predicted clusters to bins:")
    for i in clusterDict.keys():
        size = 0
        out = open("%s/bins/cluster_%s.csv" % (outdir, i), "w")
        out.write(header + "\n")
        outFASTA = open("%s/bins/cluster_%s.fa" % (outdir, i), "w")
        for j in clusterDict[i]:
            out.write(summaryDict[j] + "\n")
            outFASTA.write(">" + j + "\n")
            outFASTA.write(file[j] + '\n')
            size += int(summaryDict[j].split(",")[1])
        print("--cluster_%s.fa: %s Mb" % (i, str(size/1000000)))
        outFASTA.close()
        out.close()

    os.system("mv %s-proteins.faa %s/" % (args.g, outdir))
    os.system("mv %s-cds.ffn %s/" % (args.g, outdir))

    if args.test:
        print("\nSprayNPray finished without errors. Checking the output files to make sure everything is in check...")
    else:
        print("\nSprayNPray finished successfully. Thank you for using.")
    os.system("sleep 1")
    # inputGenome = args.g
    # os.system("mv %s.* %s/" % (args.g, outdir))
    # os.system("mv %s.spraynpray* %s/" % (allButTheLast(inputGenome, "."), outdir))
    # os.system("mv %s-* %s/" % (args.g, outdir))
    # if args.makedb:
    #     os.system("mv %s.dmnd %s/" % (args.ref, outdir))
    # os.system("mv universal.tblout %s/" % (outdir))











