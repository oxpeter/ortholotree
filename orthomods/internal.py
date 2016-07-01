"""
Functions to help identify, extract, align and construct a phylogeny for a target gene.
The primary purpose of this module is for a user to confirm orthology of their gene, not to understand in detail the phylogenetic relationships between all closely related proteins.
"""

from math import floor
import os
import re

from Bio import AlignIO
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import numpy as np


####### Classes ################################

class Consensus():
    def __init__(self, fastafile):
        self.fastafile = fastafile

        # collect all sequences and find longest:
        self.all_seqs = {}
        maxlen = 0
        for defline, seq, species in get_gene_fastas(fastafile=self.fastafile):
            self.all_seqs[defline] = seq
            if len(seq) > maxlen:
                maxlen = len(seq)

        self.maxlen = maxlen

        # initiate creation of consensus sequences
        self.consensus_pc()

    def __rep__(self):
        return '%r' % self.fastafile

    def __str__(self):
        return  "Number of sequences: %d\nMax alignment sequence length = %d" % (len(self.all_seqs), self.maxlen)

    def consensus_pc(self, keep_gaps=False):
        # iterate through each position and get percentage of highest represented marker
        self.consensus = {}
        self.consensus_pc = {}

        for i in range(self.maxlen):
            counts = {'A':0, 'B':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0,
                        'H':0, 'I':0, 'J':0, 'K':0, 'L':0, 'M':0, 'N':0,
                        'O':0, 'P':0, 'Q':0, 'R':0, 'S':0, 'T':0, 'U':0,
                        'V':0, 'W':0, 'X':0, 'Y':0, 'Z':0, '-':0, 'null':0,
                        }
            for d in self.all_seqs:
                try:
                    counts[self.all_seqs[d][i].upper()] += 1
                except KeyError:
                    counts['null'] += 1
                except IndexError:
                    counts['null'] += 1

            if not keep_gaps:
                gaps = counts['-']
                del counts['-']
            nulls = counts['null']
            del counts['null']
            max_aa_count = max(counts.values())
            self.consensus[i]    = [ aa for aa in counts if counts[aa] == max_aa_count ]
            self.consensus_pc[i] = (1.0 * max_aa_count / (sum(counts.values())))

    def make_sliding_consensus(self, window=20, window_pc=False):
        """
        Creates the global average consensus across a specified sliding window distance
        """
        if window_pc:
            window = len(self.consensus_pc.values()) * window / 100.0
        self.sliding_cons = {}
        for i in range(len(self.consensus_pc.values())):
            start = int(i - window / 2.0)
            if start < 0 :
                start = 0 # negative values will mess up the slicing
            end = int(i + window / 2.0)
            self.sliding_cons[i] = np.mean( self.consensus_pc.values()[start:end])

    def make_local_sliders(self, window=20, window_pc=False):
        """
        Creates the average level of consensus for each sequence across the specified
        window size.

        NB: if keep_gaps is false for the consensus, then a sequence with a perfect
        consensus match, but a few gaps, will have lower scores in those positions
        flanking the gaps, as a gap will not be considered a consensus sequence, and
        therefore score 0 when calculating the percentage.
        """
        if window_pc:
            window = len(self.consensus_pc.values()) * window / 100.0

        self.sliding_local = { seq:[] for seq in self.all_seqs }
        for seq in self.all_seqs:
            for i,bp in enumerate(self.all_seqs[seq]):
                start = int(i - window / 2.0)
                if start < 0 :
                    start = 0  # negative values will mess up the slicing
                end = int(i + window / 2.0)

                # calculate what percentage of local sites match the consensus
                # (note the consideration for multiple sequences --> "l in c"
                idx_pc = sum( 1 for l,c in zip(
                                        self.all_seqs[seq][start:end],
                                        self.consensus.values()[start:end],
                                            ) if l in c
                            ) / float(end - start)

                self.sliding_local[seq].append(idx_pc)

class HMMer():
    """
    A class for parsing HMMer results files.

    this class efficiently parses the hmmer result file, making all
    the different elements available for use.

    INPUT:
    handle to hmmer result file.

    ATTRIBUTES:
    query: query name
    target: list of significant target matches
    domain_seq: dictionary of domain sequences for each domain for each target
    domain_prb: domain probability scores for each alignment
    stats: indexed by target -->
           {'eval'   : Evalue for target
            'score'  : HMM score for target
            'bias'   : target bias
            'dom_e'  : e-value for best domain
            'dom_s'  : score for best domain
            'dom_b'  : bias for best domain
            'dom_exp': number of expected domains
            'dom_no' : actual number of domains found
            'desc'   : description of target (from original defline)
            }
    domain_stats: indexed by target then by domain (starting at 1) -->
           {'score'
            'bias'
            'c-eval'
            'i-eval'
            'hmmfrom'
            'hmmto'
            'alifrom'
            'alito'
            'envfrom'
            'envto'
            'acc'
            }


    """
    def __init__(self, hmmer_handle):
        self.query = None
        self.stats = {}
        self.targets = []

        self.domain_seq = {}   # the sequence of each domain
        self.domain_aln = {}   # the alignment summary for each domain
        self.domain_prb = {}   # the probability score for each domain
        self.domain_stats = {} # the result scores for each domain (indexed by target then dom)

        # initialise variables:
        current_target = None

        for line in hmmer_handle:
            if len(line) == 0 or line[0] == '#':
                continue
            else:
                cols = line.split()
                if len(cols) == 0:
                    continue


            if cols[0:5] == ['Domain', 'annotation', 'for', 'each', 'sequence']:
                in_complete = False
                in_annotation = True
            elif cols[1:3] == ['inclusion', 'threshold']:
                in_complete = False
                in_annotation = False
            elif not self.query and cols[0] == 'Query:':
                self.query = cols[1]
                in_complete = True
                in_annotation = False

            if in_complete:
                if is_number(cols[0]):
                    # load the search scores into the stats dic
                    self.stats[cols[8]] = { 'eval'   :float(cols[0]),
                                            'score'  :float(cols[1]),
                                            'bias'   :float(cols[2]),
                                            'dom_e'  :float(cols[3]),
                                            'dom_s'  :float(cols[4]),
                                            'dom_b'  :float(cols[5]),
                                            'dom_exp':float(cols[6]),
                                            'dom_no' :int(cols[7]),
                                            'desc'   :" ".join(cols[9:])}
                    self.targets.append(cols[8])
                else:
                    continue

            elif in_annotation:
                if cols[0] == '>>': # ie, new target reached
                    # update stats for last domain
                    if current_target and domain_counter > 0:
                        self.domain_seq[current_target][domain_counter] = tseq
                        self.domain_prb[current_target][domain_counter] = prob

                    # initialise variables and libraries
                    current_target = cols[1]
                    domain_counter = 0
                    self.domain_stats[current_target] = {}
                    self.domain_seq[current_target] = {}
                    self.domain_prb[current_target] = {}

                elif cols[0] == '==':
                    # update stats for last domain
                    if current_target and domain_counter > 0:
                        self.domain_seq[current_target][domain_counter] = tseq
                        self.domain_prb[current_target][domain_counter] = prob

                    # update counter
                    domain_counter += 1

                    # reset variables
                    tseq = ""
                    algn = ""
                    prob = ""

                elif cols[:4] == ['Internal', 'pipeline', 'statistics', 'summary:']:
                    # end of file wrap up #
                    # update stats for last domain
                    if current_target and domain_counter > 0:
                        self.domain_seq[current_target][domain_counter] = tseq
                        self.domain_prb[current_target][domain_counter] = prob
                    break

                elif is_number(cols[0]) and domain_counter == 0: # parse domain result table
                    self.domain_stats[current_target][int(cols[0])] = {'score':float(cols[2]),
                                                        'bias'      :float(cols[3]),
                                                        'c-eval'    :float(cols[4]),
                                                        'i-eval'    :float(cols[5]),
                                                        'hmmfrom'   :int(cols[6]),
                                                        'hmmto'     :int(cols[7]),
                                                        'alifrom'   :int(cols[9]),
                                                        'alito'     :int(cols[10]),
                                                        'envfrom'   :int(cols[12]),
                                                        'envto'     :int(cols[13]),
                                                        'acc'       :float(cols[15])}
                else:
                    if cols[0] == self.query: # query sequence
                        continue
                    elif cols[0] == current_target:
                        tseq += cols[2]
                    elif cols[-1] == 'PP': # probability scores
                        prob += cols[0]
                    ## TODO: Figure out how to efficiently parse the alignment string
                    else:
                        continue

        hmmer_handle.close()

    def __rep__(self):
        return "%r" % self.query

    def __str__(self):
        return "%s:\n%s targets" % (self.query, len(self.targets))






####### File conversion ########################
def phylipise(species, number, size=8):
    padding = size - len(species) - len(str(number))
    if padding < 1:
        unpaddedname = "%s%s" % (species, number)
        shortname = unpaddedname[:10]
    else:
        shortname = "%s%s%s" % (species, "0" * padding, number)
    return shortname

def make_phylip(fastaalignment, logfile):
    "Convert a fasta file alignment to phylip format"
    phylip_alignment = logfile[:-3] + 'phylip'

    input_handle = open(fastaalignment, 'rb')
    output_handle = open(phylip_alignment, 'w')

    alignment = AlignIO.read( input_handle, "fasta")
    AlignIO.write(alignment, output_handle, "phylip")

    input_handle.close()
    output_handle.close()

    return phylip_alignment

def fix_leaky_pipes(genename):
    return genename.replace("|", "\|")

def remove_illegal_characters(defline):
    """RaXML does not allow certain characters in the taxon name (as determined by the
    defline in the fasta file). This function takes the first whitespace separated
    'word', and converts all remaining illegal characters into underscores."""
    illegal_chars = [":", ",", ")", "(", ";", "]", "[", "'" ]
    if re.match(">", defline):
        newname = defline[1:].split()[0]
    else:
        newname = defline.split()[0]

    for char in illegal_chars:
        newname = newname.replace(char, "_")
    return newname

####### fasta file operations ##################
def count_genes(genes=[], fastafile=None):
    "evaluates the number of genes provided between a gene list and a fasta file"
    genenum = 0
    if not isinstance(genes,list):
        genes = [genes]
    genes = [ g for g in genes if g != '' ]

    if fastafile:
        # count number of genes provided:
        handle = os.popen("grep -c '^>' " + fastafile)
        result = re.search("(\d*)", handle.readline())
        handle.close()
        if result:
            try:
                genenum = int(result.group(1))
            except ValueError:
                genenum = 2
                print "ValueError calculating genenum"
                """
                putting one will ensure hmmer model is built if there is an error
                counting the number of genes in the fasta file
                """
        else:
            genenum = 2
            print "No result found for genenum"
    return len(genes), genenum

def parsefasta(fastafile, verbalise=lambda *a: None):
    """
    creates a generator that yields each successive sequence and defline.
    """
    handle = open(fastafile, 'rb')
    seq = ""
    for line in handle:
        if line[0] == '>':
            seq = seq.replace(" ","")
            if is_validfasta(seq, verbalise=verbalise):
                yield defline, seq
            defline = line.strip()
            seq = ""
        else:
            seq += line.strip()
    else:
        seq = seq.replace(" ","")
        if is_validfasta(seq, verbalise=verbalise):
            yield defline, seq

def is_validfasta(seq, verbalise=lambda *a: None):
    if len(seq) == 0:
        return False
    badcharacters = re.findall('[\(\)\!\@\#\$\%\^\&\>\<\\\|\/\:\;]',seq)
    if badcharacters:
        verbalise("R", "Invalid characters found in fastafile: %s" % " ".join(set(badcharacters)))
        return False
    else:
        return True

def find_gene(fastafile, gene, verbalise=lambda *a: None):
    genename = fix_leaky_pipes(gene)
    for defline, seq in parsefasta(fastafile):
        if re.search( '(\W)?' +  genename + '([\W\s].*)?$', defline):  #  '[\s\|\$]'
            return defline, seq
    else:
        return None, None

def find_genes(fastafiles, genes, verbalise=lambda *a: None):
    """
    This function allows extraction of multiple genes, hopefully to speed up large
    iterative searches that relied on the older find_gene function.
    NB: this function must return a different format than the old one, as it is dealing
    with multiple genes. It therefore returns a dictionary keyed by defline.

    The iterations will stop as soon as the number of entries saved in the dictionary
    is the same as the number of entries provided. If it reaches the end of the file, then
    it will return the dictionary as it stands (which may be empty, if no matches were
    found).

    the duplicates feature has only limited insurance - it will only flag duplicates found
    while still searching through the fasta files, but because the search ends once the
    size of the dictionary matches the number of genes requested, any duplicates that
    exist AFTER the last sequence parsed will not be identified as duplicates.
    """



    if isinstance(genes, str):
        genes = [genes]
    genenames = [ fix_leaky_pipes(gene) for gene in genes ]

    if isinstance(fastafiles, str):
        fastafiles = [fastafiles]

    dup_idx = {}
    duplicates = []  # to store genenames that are found more than once
    genedic = {}
    for fastafile in fastafiles:
        for defline, seq in parsefasta(fastafile):
            for g in genenames:
                if re.search( '(\W)?' +  g + '([\W\s].*)?$', defline):  #  '[\s\|\$]'

                    # duplicate insurance check:
                    if g in dup_idx:
                        duplicates.append(g)
                        dup_idx[g].append(defline)
                    else:
                        dup_idx[g] = [defline]

                    genedic[defline] = seq
                    if len(genedic) == len(genenames):
                        return genedic
    else:
        # remove duplicates that were found:
        for g in duplicates:
            for defline in dup_idx[g]:
                del genedic[defline]
        return genedic

def get_gene_fastas(genes=None, fastafile=None,
                    startpos=0, endpos=None,
                    specieslist = [], species=None,
                    comment=None, short=False,
                    dbpaths={}, verbalise=(lambda *a: None)):
    """
    Can either be given as a transcript name to be searched within the peptide databases,
    or can be a fasta file.

    Function used to return the species as the third element in the tuple, but no longer
    does so. The third None still exists so the dependent functions don't break. This will
    hopefully be removed in a future upgrade.
    """

    if genes:
        if dbpaths=={} or specieslist == []:
            verbalise("R", "No database or specieslist provided :(")
            yield None, None, None
            raise StopIteration

        # extract all sequences:
        if species in specieslist:
            reportedspecies = species
            try:
                seqdic = find_genes(dbpaths[species + '_lpep'], genes, verbalise=verbalise)
            except KeyError:
                seqdic = find_genes(species, genes, verbalise=verbalise)

        else:
            seqdic = find_genes([dbpaths[sp + '_lpep'] for sp in specieslist],
                                genes,
                                verbalise=verbalise)

        for defline, seq in seqdic.items():
            # create fasta file from extracted sequence:
            if seq:
                if short:
                    name = phylipise(defline, short)
                    defline = ">%s %s" % (name, comment)
                elif comment:
                    defline = "%s %s" % (defline, comment)
                else:
                    defline = "%s" % (defline)

                yield defline, seq[startpos:endpos], None

    if fastafile:
        for defline, seq in parsefasta(fastafile, verbalise=verbalise):
            yield defline, seq[startpos:endpos], None

    if not genes and not fastafile:
        yield None, None, None

def rank_scores(homologlist, thresh1=0, thresh2=None, genename=None, outfile=None, showplot=False):
    yvalues = sorted([val[1] for val in homologlist.values()], reverse=True)
    plt.plot(yvalues)
    score_cutoff = thresh1 * max(yvalues)
    sample_cutoff = sum(1 for s in yvalues if s >= thresh1 * max(yvalues))
    plt.axhline( score_cutoff , color='r' )
    if thresh2:
        plt.axhline( thresh2 * max(yvalues) , color='r' )
    plt.axvline( sample_cutoff -1 , color='g' )
    plt.text(sample_cutoff + 1,score_cutoff + 10 , "(%d,%d)" % (sample_cutoff,score_cutoff) )
    plt.xlabel("Gene rank")
    plt.ylabel("Phmmer score")
    plt.title("Ranking of phmmer scores for alignment with %s" % genename)
    if outfile:
        plt.savefig(outfile, format='png')
    if showplot:
        plt.show()
    else:
        plt.close()

def find_holes(seq):
    allholes = re.findall('([A-Za-z]+)(-{50,}[A-Za-z])', seq)
    dists = []
    for hole in allholes:
        dists.append(len(hole[0]))
        dists.append(len(hole[1]))
    return dists

def find_biggest_hole(seq):
    allholes = re.findall('-+', seq)
    if len(allholes) > 0:
        biggesthole = len(sorted(allholes)[-1])
        pattern = '(.+)-{' + str(biggesthole) + '}(.+)'
        bigsearch = re.search(pattern, seq)
        return len(bigsearch.group(1)), biggesthole, len(bigsearch.group(2))
    else:
        return 0,0,len(seq)

def get_pcmatch(seq):
    if len(seq) == 0:
        return 0, 0
    minigaps = len(re.findall('-', seq))
    width = len(seq)
    matches = width - minigaps
    assert matches >= 0
    pcmatch = 1.0 * matches / width
    return width, pcmatch

def display_alignment(fastafile, conversiondic={}, outfile=None, showplot=True,
                        gapthresh=0.05):
    fig = build_alignment(fastafile, conversiondic, gapthresh=gapthresh)
    if outfile:
        fig.savefig(outfile, format='png')
    if showplot:
        fig.show()
    else:
        plt.close()

def get_graphing_name(defline, conversiondic={}, truncate_name=False):
    namesearch = re.search('^>*(\S+)', defline)
    if namesearch:
        genename = namesearch.group(1)
    else:
        genename = defline.strip()
    if genename in conversiondic:
        fullname = conversiondic[genename][0]
    else:
        fullname = genename

    if truncate_name and len(fullname) > 11:
        graphingname = "...".join([genename[:6],genename[-5:]])
    else:
        graphingname = fullname

    return graphingname

def build_alignment(fastafile, conversiondic={}, img_width=10, gapthresh=0.05,
                    truncate_name=False, graph_style='consensus'):
    """
    Draw an alignment graph in the vein of BLAST alignment results on NCBI.
    colour scale represents the percentage of alignment positions filled with actual
    sequence, but does not represent the fit of that alignment. This is indicated by
    adding a consensus bar at the bottom - high consensus meaning most amino acids/base
    pairs are identical in a given sliding window.

    graph_style can be 'consensus', 'amino', or 'block'
    """
    # set similarity-based color scheme:
    nrml = mpl.colors.Normalize(vmin=0, vmax=1)
    sm = plt.cm.ScalarMappable(cmap=cm.jet, norm=nrml)
    sm._A = []

    #get consensus percentages for each position:
    consensus = Consensus(fastafile)

    # calculate sliding average:
    consensus.make_sliding_consensus(20)
    sliding_colors = sm.to_rgba(consensus.sliding_cons.values())
    consensus.make_local_sliders(20)

    # setting color maps:
    if graph_style == 'amino':
        # color based on the peptide sequence
        acma = {'A':[200,200,200,256], 'B':[0,0,0,256],      'C':[230,230,0,256],
                'D':[230,10,10,256],
                'E':[230,10,10,256],  'F':[50,50,170,256],  'G':[235,235,235,256],
                'H':[130,130,210,256],'I':[15,130,15,256],  'J':[0,0,0,256],
                'K':[20,90,255,256],  'L':[15,130,15,256],  'M':[230,230,0,256],
                'N':[0,220,220,256],  'O':[0,0,0,256],      'P':[220,150,130,256],
                'Q':[0,220,220,256],  'R':[20,90,255,256],  'S':[250,150,0,256],
                'T':[250,150,0,256],  'U':[0,0,0,256],      'V':[15,130,15,256],
                'W':[180,90,180,256], 'X':[0,0,0,256],      'Y':[50,50,170,256],
                'Z':[0,0,0,256],      '-':[256,256,256,0],  'null':[256,256,256,0],
                }
        acm = {}
        for aa in acma:
            acm[aa] = [ n/256.0 for n in acma[aa] ]
        # assign colors to each sequence based on percentage consensus:
        colorme = { k:[] for k in consensus.all_seqs }
        for defline, seq in consensus.all_seqs.items():
            for aa in seq:
                colorme[defline].append(acm[aa])

    elif graph_style == 'consensus':
        # assign colors to each sequence based on percentage consensus:
        colorme = { k:[] for k in consensus.all_seqs }

        for defline, seq in consensus.all_seqs.items():
            for i,pc in enumerate(consensus.sliding_local[defline]):
                if consensus.all_seqs[defline][i] == '-':
                    colorme[defline].append((1.0,1.0,1.0,0.0))
                else:
                    colorme[defline].append(sm.to_rgba(pc))

    if graph_style in ['consensus', 'amino']:
        # get coords for alignment (also sort sequences alphabetically):
        graphingnames = {
                defline:get_graphing_name(
                                    defline,
                                    conversiondic,
                                    True
                                        ) for defline in consensus.all_seqs }
        keynames = sorted([ (graphingnames[d],d) for d in colorme ],
                            reverse=True,
                            key=lambda x: x[0])
        name_pos = np.arange(len(colorme)) + 0.5
        y_frame = { k:y for k,y in zip(keynames, name_pos)}

        # set all plotting values into lists for loading into plt.barh:
        y_pos, lefts, widths, colors, bh_lefts, bh_widths = [], [], [], [], [], []
        for ((gname, defline), namepos) in y_frame.items():
            y_pos += [namepos] * len(consensus.all_seqs[defline])
            lefts += range(len(consensus.all_seqs[defline]))
            widths += [1] * len(consensus.all_seqs[defline])
            colors += colorme[defline]


    elif graph_style == 'block':
        # the original gap-based color scheme:
        # find gaps:
        graph_points = {}
        for defline, seq in consensus.all_seqs.items():
            # determine the smallest reportable gap size is:
            repgap = int(gapthresh * len(seq))

            # get distances and coverage percentages
            points = re.search('^(-*)(\S+[A-Za-z])(-*)$', seq)
            if points:
                pattern = '-{' + str(repgap) + ',}'
                fragments = re.findall(pattern, points.group(2))
                # set starting pos to beginning of matching sequence:
                spos = len(points.group(1))
                if len(fragments) > 0:
                    """
                    dists is a list of tuples, each tuple containing the start position of
                    a large gap,  the length of the gap, the start of the preceding non-gap
                    fragment, its width and the % match.
                    """
                    dists = []
                    for frag in fragments:
                        nextgap = seq.find(frag, spos)
                        width, pcmatch = get_pcmatch(seq[spos:nextgap])
                        dists.append((nextgap,len(frag), spos, width, pcmatch))
                        spos = nextgap + len(frag)

                    else:
                        lastfrag = points.group(3)
                        nextgap = len(seq) - len(lastfrag)
                        width, pcmatch = get_pcmatch(seq[spos:nextgap])
                        dists.append((0,0, spos, width, pcmatch))

                else:
                    width, pcmatch = get_pcmatch(points.group(2))
                    dists = [(0,0,spos, width, pcmatch)]

            else:
                dists = [(0,0,0,1,0)]


            # get name (convert if possible):
            graphingname = get_graphing_name(defline, conversiondic, True)
            graph_points[graphingname] = dists

        # get coords for alignment:
        keynames = sorted(graph_points.keys(), reverse=True)
        name_pos = np.arange(len(graph_points)) + 0.5
        y_frame = { k:y for k,y in zip(keynames, name_pos)}

        y_pos, lefts, widths, colors, bh_lefts, bh_widths = [], [], [], [], [], []
        for k in keynames:
            for dists in graph_points[k]:
                y_pos.append(y_frame[k])
                lefts.append(dists[2])
                widths.append(dists[3])
                colors.append(sm.to_rgba(dists[4]))
                #bh_lefts.append(dists[0])
                #bh_widths.append(dists[1])
        keynames = [ (n,n) for n in keynames ] # to make compatible with amino and consensus

    # plot graph:
    """
    if 30 > len(keynames) :
        fig = plt.figure(figsize=(img_width,img_width*1))
    elif 60 > len(keynames) >= 30:
        fig = plt.figure(figsize=(img_width,img_width*2))
    elif 90 > len(keynames) >= 60:
        fig = plt.figure(figsize=(img_width,img_width*3))
    else:
    """
    fig = plt.figure(figsize=(img_width,int(len(keynames)/3) + 2))

    # plot alignments:
    ax1 = plt.subplot2grid((12,10),(0,0), colspan=9, rowspan=9)
    plt.barh(left=lefts,
                width=widths,
                bottom=y_pos,
                height=0.8,
                color=colors,
                edgecolor=colors)

    #if graph_style == 'block':
    #    plt.barh(left=bh_lefts, width=bh_widths, bottom=y_pos, height=0.8, color='white',
    #            alpha=0.5)
    plt.yticks(name_pos + 0.4, [k[0] for k in keynames])
    plt.xlabel("position (aa)")
    plt.title("Peptide alignment and amino acid sequence")
    plt.tight_layout()

    # plot legend:
    ax2 = plt.subplot2grid((12,10),(0,9), colspan=1,rowspan=5)
    cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cm.jet, norm=nrml, orientation='vertical')
    plt.tight_layout()

    # plot consensus colors:
    ax3 = plt.subplot2grid((12,10),(9,0), colspan=9,rowspan=2)
    size = consensus.maxlen
    plt.barh(left=range(size), bottom=[1]*size,
                height=[0.8]*size, width=[1]*size,
                color=sliding_colors,
                edgecolor=sliding_colors,)
    plt.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
    plt.xlabel("20 aa sliding global average consensus (%)")
    plt.tight_layout()
    return fig


####### hmmer functions ########################

def parse_the_hmmer(handle):
    """
    ####### DEPRECATED ######
    parses the protein matches from a hmmer search and returns a dictionary of peptides
    and their associated score and p-value.
    """
    parse_dic = {}
    lcount = 0
    collected = 0
    for line in handle:
        lcount += 1
        if line[0] in ['#', 'Q', 'D', 'S']:
            continue
        elif len(line) < 2:
            continue
        elif line[0] == '>':
            break
        elif line.split()[1] in ['hits', 'inclusion', 'annotation']:
            break
        else:
            try:
                score = float(line.split()[1])
                pvalue = eval(line.split()[0])
            except ValueError:
                continue
            else:
                parse_dic[line.split()[8]] = (pvalue, score)

    handle.close()
    return parse_dic

####### Phylogeny creation/manipulation ########


def rename_newick(raxml_final, conversiondic={}):
    #replace short names in newick file with full names
    if os.path.exists(raxml_final):
        handle = open(raxml_final, 'rb')
        newfile = raxml_final[:-3] + "scored.nwk"
        newhandle = open(newfile, 'w')
        for line in handle: # should only be one line in file
            for shortname in conversiondic:
                line = re.sub(shortname,
                            conversiondic[shortname][0] + "_" + str(conversiondic[shortname][1]),
                            line )
            newhandle.write(line)
        handle.close()
        newhandle.close()
    else:
        newfile = None
    return newfile

################################################
def is_number(s):
    try:
        float(s)
    except ValueError:
        return False
    else:
        return True

def main():
    pass

if __name__ == '__main__':
    main()





