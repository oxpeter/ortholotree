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

def trim_name_dross(genename):
    if genename.find('|') >= 0:
        return genename[genename.find('|')+1:]
    else:
        return genename

####### fasta file operations ##################
def count_genes(genes=[], fastafile=None):
    "evaluates the number of genes provided between a gene list and a fasta file"

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
    badcharacters = re.findall('[\(\)\!\@\#\$\%\^\&\*\>\<\\\|\/\:\;]',seq)
    if badcharacters:
        verbalise("R", "Invalid characters found in fastafile: %s" % " ".join(set(badcharacters)))
        return False
    else:
        return True

def find_gene(fastafile, gene, verbalise=lambda *a: None):
    genename = trim_name_dross(gene)
    for defline, seq in parsefasta(fastafile):
        if re.search( '[^\w]' + genename + '(\s.*)?$', defline):  #  '[\s\|\$]'
            return defline, seq
    else:
        return None, None

def get_gene_fastas(genes=None, fastafile=None,
                    startpos=0, endpos=None,
                    specieslist = [], species=None,
                    comment=None, short=False,
                    dbpaths={}, verbalise=(lambda *a: None)):
    """
    Can either be given as a transcript name to be searched within the peptide databases,
    or can be a fasta file.
    """

    if genes:
        if dbpaths=={} or specieslist == []:
            yield None, None, None
            raise StopIteration
        if isinstance(genes, str):
            genes = [genes]
        for gene in genes:
            if species in specieslist:
                reportedspecies = species
                try:
                    defline, seq = find_gene(dbpaths[species + '_lpep'], gene, verbalise=verbalise)
                except KeyError:
                    defline, seq = find_gene(species, gene, verbalise=verbalise)
                if seq:
                    seq = seq[startpos:endpos]

                    if len(seq) == 0:
                        verbalise("R",
                            "Transcript %s (%s) could not be found in the database." % (gene, species))
                        defline, seq, reportedspecies = None, None, None

            else:   # if no species is given, check all peptide files in dbpaths
                for sp in specieslist:
                    seq = ""
                    defline, seq = find_gene(dbpaths[sp + '_lpep'], gene, verbalise=verbalise)
                    if seq:
                        seq = seq[startpos:endpos]
                        if len(seq) > 0:   # found a match!
                            reportedspecies = sp
                            break
                else:
                    verbalise("R",
                        "Transcript %s (%s) could not be found in the database." % (gene, species))
                    defline, seq, reportedspecies = None, None, None

            # create fasta file from extracted sequence:
            if seq:
                if short:
                    name = phylipise(reportedspecies, short)
                    defline = ">%s (%s) %s" % (name, reportedspecies, comment)
                else:
                    defline = ">%s (%s) %s" % (gene, reportedspecies, comment)

            yield defline, seq, reportedspecies

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
    pcmatch = 10.0 * matches / width
    return width, pcmatch

def consensus_pc(fastafile):
    all_seqs = {}
    maxlen = 0
    for defline, seq, species in get_gene_fastas(fastafile=fastafile):
        all_seqs[defline] = seq
        if len(seq) > maxlen:
            maxlen = len(seq)

    consensus = {}
    print "Max alignment sequence length = %d" % (maxlen)
    for i in range(maxlen):
        counts = {'A':0, 'B':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0,
                    'H':0, 'I':0, 'J':0, 'K':0, 'L':0, 'M':0, 'N':0,
                    'O':0, 'P':0, 'Q':0, 'R':0, 'S':0, 'T':0, 'U':0,
                    'V':0, 'W':0, 'X':0, 'Y':0, 'Z':0, '-':0, 'null':0,
                    }
        for d in all_seqs:
            try:
                counts[all_seqs[d][i].upper()] += 1
            except KeyError:
                counts['null'] += 1
            except IndexError:
                counts['null'] += 1

        consensus[i] = 1.0 * max(counts.values()) / (sum(counts.values()) - counts['null'])
    return consensus

def sliding_average(float_list, window=20, window_pc=False):
    if window_pc:
        window = len(float_list) * window / 100
    sliding_averages = []
    for i in range(len(float_list)):
        sliding_averages.append(np.mean(float_list[i:i+window]))
    return sliding_averages

def display_alignment(fastafile, conversiondic={}, outfile=None, showplot=True,
                        gapthresh=0.05):
    fig = build_alignment(fastafile, conversiondic, gapthresh=gapthresh)
    if outfile:
        fig.savefig(outfile, format='png')
    if showplot:
        fig.show()
    else:
        fig.close()

def build_alignment(fastafile, conversiondic={}, img_width=10, gapthresh=0.05):
    """
    Draw an alignment graph in the vein of BLAST alignment results on NCBI.
    colour scale represents the % match as base 10, to allow flooring of actual percentages
    to find appropriate colour group. Key number represents the minimum value the % match
    must be to join that colour group.
    """
    #cmap = {0:'white', 1:'silver', 2:'tan' ,3:'cornflowerblue' ,4:'blue' ,5:'darkcyan' ,6:'green', 7:'gold' ,8:'orangered' ,9:'red' ,10:'maroon'}
    #cmap = cm.cool
    print "setting color maps"
    nrml = mpl.colors.Normalize(vmin=0, vmax=1)
    sm = plt.cm.ScalarMappable(cmap=cm.jet, norm=nrml)
    sm._A = []

    graph_points = {}
    hole_points = {}

    #get consensus percentages for each position:
    cons = consensus_pc(fastafile)
    # calculate sliding average:
    slave = sliding_average(cons.values())
    sliding_colors = sm.to_rgba(slave)

    # find gaps:
    for defline, seq, species in get_gene_fastas(fastafile=fastafile):
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
        namesearch = re.search('>(\S+)', defline)
        if namesearch:
            genename = namesearch.group(1)
        else:
            genename = defline.strip()
        if genename in conversiondic:
            fullname = conversiondic[genename][0]
        else:
            fullname = genename

        graph_points[fullname] = dists

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
            colors.append(sm.to_rgba(dists[4]/10.0))
            bh_lefts.append(dists[0])
            bh_widths.append(dists[1])

    # plot graph:
    if 30 > len(keynames) :
        fig = plt.figure(figsize=(img_width,img_width*1))
    elif 60 > len(keynames) >= 30:
        fig = plt.figure(figsize=(img_width,img_width*2))
    elif 90 > len(keynames) >= 60:
        fig = plt.figure(figsize=(img_width,img_width*3))
    else:
        fig = plt.figure(figsize=(img_width,int(len(keynames)/2.8)))

    # plot alignments:
    ax1 = plt.subplot2grid((12,10),(0,0), colspan=9, rowspan=9)
    plt.barh(left=lefts,    width=widths,    bottom=y_pos, height=0.8, color=colors)
    plt.barh(left=bh_lefts, width=bh_widths, bottom=y_pos, height=0.8, color='white',
            alpha=0.5)
    plt.yticks(name_pos + 0.4, keynames)
    plt.xlabel("position (aa)")
    plt.title("Gaps in alignment (%)")
    plt.tight_layout()

    # plot legend:
    ax2 = plt.subplot2grid((12,10),(0,9), colspan=1,rowspan=5)
    cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cm.jet, norm=nrml, orientation='vertical')
    plt.tight_layout()

    # plot consensus colors:
    ax3 = plt.subplot2grid((12,10),(9,0), colspan=9,rowspan=2)
    size = len(slave)
    plt.barh(left=range(size), bottom=[1]*size,
                height=[0.8]*size, width=[1]*size,
                color=sliding_colors,
                edgecolor=sliding_colors,)
    plt.tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
    plt.xlabel("25 aa sliding average consensus (%)")
    plt.tight_layout()
    return fig


####### hmmer functions ########################

def parse_the_hmmer(handle):
    """
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

def main():
    pass

if __name__ == '__main__':
    main()





