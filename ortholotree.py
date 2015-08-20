#! /usr/bin/env python
"""
A wrapper to identify, extract, align and phylogenise a target gene. The primary purpose
of this module is for a user to confirm orthology of their gene, not to understand in
detail the phylogenetic relationships between all closely related proteins.
"""

import argparse
import os
import tempfile

import config
from orthomods import internal, external

############################################################################

dbpaths = config.import_paths()

############################################################################

def define_arguments():
    parser = argparse.ArgumentParser(description=
            "A module to perform a variety of gene term related analyses")
    ### input options ###
    # logging options:
    parser.add_argument("-q", "--quiet", action='store_true',default=False,
                        help="print fewer messages and output details")
    parser.add_argument("-o", "--output", type=str, default='genematch.out',
                        help="specify the filename to save results to")
    parser.add_argument("-d", "--directory", type=str,
                        help="specify the directory to save results to")
    parser.add_argument("-D", "--display_on", action='store_true',default=False,
                        help="display graph results (eg for p value calculation)")


    # data file options:
    parser.add_argument("-f", "--fasta", type=str,
                        help="Fasta file of gene to analyse (must be amino acids)")
    parser.add_argument("-g", "--gene", type=str,
                        help="""Gene or transcript name. This name will be checked
                        against the peptide database to see if there is a match.""")
    parser.add_argument("-s", "--species", type=str, default="",
                        help="""Four letter abbreviation of species, to speed up
                        searching for genes. E.g Cerapachys biroi = Cbir""")
    parser.add_argument("-n", "--name_conversion", type=str,
                        help="""Name conversion file for RAxML only analysis. Format is
                        the same as that outputted by a standard run. If rerunning an
                        analysis, it is therefore possible to supply the mafft.out
                        file and the name_conversion.txt file to create a phylogeny with
                        the original peptide names.""")

    # analysis options:
    parser.add_argument("-B", "--buildhmmer", action='store_true',
                        help="""Use sequences supplied in fasta file or gene to build a
                        hmmer model and use hmmer to extract sequence matches (more
                        sensitive than using the default phmmer setting, but requires
                        knowledge of the genes to input). This may be useful after a
                        preliminary run, using a resulting orthologous cluster as the
                        new input.""")
    parser.add_argument("-c", "--mincollect", type=int, default=2,
                        help="""minimum number of homologs to collect from each species,
                        regardless of how poor the similarity""")
    parser.add_argument("-t", "--scorethresh", type=float, default=0.5,
                        help="""minimum fraction of the score of the best matching
                        gene for a species for including additional matches""")
    parser.add_argument("-p", "--threads", type=int, default=2,
                        help="number of threads to use for RAxML calculations")
    parser.add_argument("-b", "--bootstrap", type=int,
                        help="""Perform rapid bootstrap analysis. Requires a random
                         integer for the bootstrap seed.""")
    parser.add_argument("-a", "--globalthresh", type=float, default=0.2,
                        help="""Identify the top proportion of all results globally.""")
    parser.add_argument("-R", "--raxml_only", action='store_true',
                        help="""Construct phylogeny using RAxML using fasta alignment
                        provided in fasta file """)
    parser.add_argument("-x", "--exclude_genes", type=str,
                        help="""A comma-separated list of genes to exclude from the
                        alignment and phylogeny""")
    parser.add_argument("-e", "--exclude_species", type=str,
                        help="""A comma-separated list of species to exclude from the
                        alignment and phylogeny. Use the four-letter abbreviation.""")
    parser.add_argument("-l", "--maxlength", type=int,
                        help="""If provided, will remove all genes longer than this
                        size. Useful for removing concatenated genes that otherwise
                        are orthologous. It is not recommended that you use this
                        flag until you have looked at your results without it, and
                        preferably tried to eliminate long genes by using better search
                        models.""")
    parser.add_argument("-m", "--minlength", type=int,
                        help="""If provided, will remove all genes shorter than this
                        size. It is not recommended that you use this
                        flag until you have looked at your results without it, and
                        preferably tried to eliminate short genes by using better search
                        models.""")

    return parser



################################################

if __name__ == '__main__':
    parser = define_arguments()
    args = parser.parse_args()

    verbalise = config.check_verbose(not(args.quiet))
    logfile = config.create_log(args, outdir=args.directory, outname=args.output)

    temp_dir = tempfile.mkdtemp()

    ############ if generating phylogeny only ###########
    if args.raxml_only:
        if not args.fasta:
            verbalise("R", "No fasta alignment was supplied!")
            exit()
        verbalise("B", "Running RAxML analysis to construct phylogeny using supplied alignment")
        phylip_alignment = internal.make_phylip(args.fasta, logfile)
        raxml_final = external.raxml_phylogeny(phylip_alignment, logfile, bootstrap=args.bootstrap)
        if args.name_conversion:
            handle = open(args.name_conversion, 'rb')
            conv_dic = { line.split()[0]:(line.split()[2], line.split()[1]) for line in handle }
            handle.close()
            internal.rename_newick(raxml_final, conversiondic=conv_dic)
        exit()

    # initialise dictionary of all accessible longest non-redundant peptide fasta files
    specieslist = [ "Ador", "Aech", "Aflo", "Amel", "Apis", "Aros", "Bimp", "Bmor",
                "Bter", "Cele", "Cflo", "Csol", "Dcit", "Fari", "Hsal", "Lhum",
                "Mdem", "Mpha", "Mrot", "Nvit", "Oabi", "Pbar", "Pcan", "Sinv",
                "Tcas", "Waur", "Cbir", "Ebur", "Dmel" ]
    lpep_paths = { s:dbpaths[s+'_lpep'] for s in specieslist }

    ######### Get protein sequences #########
    genes = config.make_a_list(args.gene)
    homologlist = external.get_similar_sequences(temp_dir,
                                        buildhmmer=args.buildhmmer,
                                        fastafile=args.fasta,
                                        specieslist=specieslist,
                                        species=args.species,
                                        mincollect=args.mincollect,
                                        globalthresh=args.globalthresh,
                                        localthresh=args.scorethresh)

    ######### Extract identified sequences from LNRP fasta files #########
    conv_handle = open(logfile[:-3] + 'name_conversion.txt', 'w')
    conv_dic = {}
    itercount = 0
    previousseq = ""
    seqdic = {}         # loaded up to remove duplicate sequences
    excluded_genes = config.make_a_list(args.exclude_genes)
    excluded_species = config.make_a_list(args.exclude_species)

    for homolog in sorted(homologlist):
        # remove excluded genes before bothering to look up their sequence:
        searchname = internal.trim_name_dross(homolog)
        if searchname in excluded_genes:
            continue
        if homologlist[homolog][0] in excluded_species:
            continue

        # extract sequences of remaining genes and add to conversion dictionary
        itercount += 1

        for defline, seq, spec in internal.get_gene_fastas(genes=[searchname],
                                    species=homologlist[homolog][0],
                                    fastafile=None,
                                    specieslist = specieslist,
                                    comment=homologlist[homolog][1],
                                    short=itercount):

            # filter based on size:
            if args.maxlength and len(seq) > args.maxlength:
                continue
            if args.minlength and len(seq) < args.minlength:
                continue

            seqdic[seq] = defline

        shortname = internal.phylipise(homologlist[homolog][0], itercount)
        conv_handle.write("%s %-5d %s\n" % (shortname,
                                          homologlist[homolog][1],
                                          homolog))
        conv_dic[shortname] = (homolog, homologlist[homolog][1])
    conv_handle.close()

    verbalise("G", "%d non-excluded homologous sequences found" % len(seqdic))

    # write multifasta file:
    homolog_fasta = os.path.join(temp_dir,"homolog.fasta")
    handle = open(homolog_fasta, 'w')
    for seq in seqdic:
        fastaseq = "%s\n%s\n" % (seqdic[seq], seq)
        handle.write(fastaseq)
    handle.close()

    # show score curve for identified genes:
    internal.rank_scores(homologlist, args.scorethresh, args.globalthresh, args.gene,
                logfile[:-3] + "ranking.png", showplot=args.display_on)
    # show raw protein sizes:
    internal.display_alignment(  homolog_fasta,
                        conversiondic=conv_dic,
                        outfile=logfile[:-3] + 'homologs.png',
                        showplot=args.display_on)

    ######### MAFFT alignment of extracted sequences #########
    mafft_alignment = logfile[:-3] + 'mafft.fa'
    external.mafft_align(homolog_fasta, mafft_alignment)
    internal.display_alignment(mafft_alignment,
                        conversiondic=conv_dic,
                        outfile=logfile[:-3] + 'mafft.png',
                        showplot=args.display_on)

    ######### RaXML phylogenetic analysis of alignment #########
    """
    Using PROT-LG-GAMMA model, a single tree and no bootstrapping (though all of these
    could be setup to allow overriding for fringe case analyses).
    """
    verbalise("B", "Running RAxML analysis to construct phylogeny")
    phylip_alignment = internal.make_phylip(mafft_alignment, logfile)
    if args.bootstrap:
                raxml_best = external.raxml_phylogeny(phylip_alignment,
                                        logfile,
                                        bootstrap=False)
                best_renamed = internal.rename_newick(raxml_best, conversiondic=conv_dic)
                raxml_bstrap = external.raxml_phylogeny(phylip_alignment,
                                        logfile,
                                        bootstrap=args.bootstrap)
                bstrap_renamed = internal.rename_newick(raxml_bstrap, conversiondic=conv_dic)
                final_tree = external.apply_boostrap(best_renamed, bstrap_renamed, logfile)
                verbalise("Y",
                    "Best tree with bootstrap support can be found at %s" % final_tree)
    else:
        raxml_final = external.raxml_phylogeny(phylip_alignment, logfile, bootstrap=args.bootstrap)
        raxml_renamed = internal.rename_newick(raxml_final, conversiondic=conv_dic)

    # clean up temp files and directory
    for file in [ homolog_fasta, ]:
        if os.path.exists(file):
            os.remove(file)

    os.rmdir(temp_dir)  # dir must be empty!