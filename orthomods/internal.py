"""
Functions to help identify, extract, align and construct a phylogeny for a target gene.
The primary purpose of this module is for a user to confirm orthology of their gene, not to understand in detail the phylogenetic relationships between all closely related proteins.
"""

from math import floor
import os
import re

from Bio import AlignIO
import matplotlib.pyplot as plt
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
def extractseq(geneID, db="", startpos=0, endpos=-1):
    """ extracts sequence of geneID from the current annotations. type is cds,  pep or fasta.
    """
    geneseq = ""
    fobj = open(db, 'rb')
    for line in fobj:
        if line[0] == '>':
            query = re.search( geneID + '[\s]', line)
        if query:
            thisline = fobj.next()

            while thisline[0] != '>':
                geneseq += thisline.strip()
                try:
                    thisline = fobj.next()
                except StopIteration:
                    break
            else:
                break
    fobj.close()
    return geneseq[startpos:endpos]

def get_gene_fastas(genes=None, species=None, fastafile=None,
                    specieslist = [], comment=None, short=False, dbpaths={}):
    """
    Can either be given as a transcript name to be searched within the peptide databases,
    or can be a fasta file.
    """
    if genes:
        for gene in genes:
            if species in specieslist:
                reportedspecies = species
                seq = extractseq(gene, db=dbpaths[species + '_lpep'])

                if len(seq) == 0:
                    print "Transcript %s could not be extracted from the LNRP database for species %s" % (gene, species)
                    exit()
            else:   # if no species is given, check all LNRP files
                for sp in specieslist:
                    seq = extractseq(gene, db=dbpaths[sp + '_lpep'])
                    if len(seq) > 0:   # found a match!
                        reportedspecies = sp
                        break
                else:
                    print "Transcript %s could not be extracted from the LNRP database for species %s" % (gene, species)
                    defline, seq, reportedspecies = None, None, None

            # create fasta file from extracted sequence:
            if short:
                name = phylipise(reportedspecies, short)
                defline = ">%s (%s) %s" % (name, reportedspecies, comment)
            else:
                defline = ">%s (%s) %s" % (gene, reportedspecies, comment)

            yield defline, seq, reportedspecies

    elif fastafile:
        handle = open(fastafile, 'rb')
        seq = ""
        for line in handle:
            if line[0] == '>':
                if seq != "":
                    yield defline, seq, None
                defline = line.strip()
                seq = ""
            else:
                seq += line.strip()
        else:
            yield defline, seq, None

    else:
        yield None, None, None

def rank_scores(homologlist, thresh1, thresh2=None, genename=None, outfile=None, showplot=False):
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

def display_alignment(fastafile, conversiondic={}, outfile=None, showplot=True):
    """
    Draw an alignment graph in the vein of BLAST alignment results on NCBI.
    colour scale represents the % match as base 10, to allow flooring of actual percentages
    to find appropriate colour group. Key number represents the minimum value the % match
    must be to join that colour group.
    """
    cmap = {0:'white', 1:'silver', 2:'tan' ,3:'cornflowerblue' ,4:'blue' ,5:'darkcyan' ,6:'green', 7:'gold' ,8:'orangered' ,9:'red' ,10:'maroon'}


    graph_points = {}
    hole_points = {}
    for defline, seq, species in get_gene_fastas(fastafile=fastafile):
        # determine the smallest reportable gap size is:
        repgap = int(len(seq)/10)

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
            colors.append(cmap[floor(dists[4])])
            bh_lefts.append(dists[0])
            bh_widths.append(dists[1])

    # plot graph:
    if 75 > len(keynames) > 25:
        plt.figure(figsize=(10,10))
    elif len(keynames) >= 75:
        plt.figure(figsize=(10,20))
    plt.barh(left=lefts,    width=widths,    bottom=y_pos, height=0.8, color=colors)
    plt.barh(left=bh_lefts, width=bh_widths, bottom=y_pos, height=0.8, color='white',
            alpha=0.5)
    plt.yticks(name_pos + 0.4, keynames)
    plt.xlabel("position (aa)")
    plt.title("Alignment of genes")
    plt.tight_layout()

    if outfile:
        plt.savefig(outfile, format='png')
    if showplot:
        plt.show()
    else:
        plt.close()

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
        if len(line) < 2:
            continue
        if line.split()[1] in ['hits', 'inclusion', 'annotation']:
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





