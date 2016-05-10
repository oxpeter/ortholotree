"""
Wrappers for the external programs. MAFFT is a sequence aligner, RAxML constructs
phylogenies, and HMMER performs protein searches.
"""
# TODO: remove all dependencies on the internal module
# TODO: push as much of this code into the internal module, making this only
#       a set of bare-bones wrappers.

import os

from orthomods import internal


####### MAFFT functions ########################
def mafft_align(inputfasta, outputfasta):
    mafft = os.popen( 'mafft --quiet ' + inputfasta )
    handle = open( outputfasta, 'w')
    for line in mafft:
        handle.write(line)
    handle.close()
    mafft.close()

####### HMMer functions ########################
def get_similar_sequences(temp_dir, buildhmmer=False, fastafile=None,
                        specieslist={}, species=None, genes=[], dbpaths={},
                        mincollect=2, globalthresh=0.2, localthresh=0.8,
                        verbalise=lambda *a: None):

    # clean gene list type and content:
    if not isinstance(genes, list):
        genes = [genes]
    genes = [ g for g in genes if g != '' ]

    # count genes provided:
    genelist_num, fasta_num = internal.count_genes(genes, fastafile)
    verbalise("Y", "Genelist size:%d\nFasta size:%d" % (genelist_num, fasta_num ))

    # if fasta files are provided, create a temp fastafile to search against with hmmer:
    if fastafile:
        extra_file = os.path.join(temp_dir, "query_fasta")
        handle = open(extra_file, 'w')
        for defline, seq in internal.parsefasta(fastafile):
            handle.write(">%s\n%s\n" % (defline, seq))
        handle.close()

        extra_file_search = extra_file
    else:
        extra_file_search = None

    if genelist_num + fasta_num > 1:
        buildhmmer = True

    if buildhmmer:
        hmminput = os.path.join(temp_dir, "hmminput.fa")
        handle = open(hmminput, 'w')
        seqcount = 0
        verbalise("B", "Extracting sequence data from %d peptides" % len(genes))

        for defline, seq, species in internal.get_gene_fastas(genes=genes,
                                                    species=None,
                                                    fastafile=fastafile,
                                                    specieslist=specieslist,
                                                    dbpaths=dbpaths):

            if seq:
                seqcount += 1
                fasta_seq = "%s\n%s\n" % (defline, seq)
                handle.write(fasta_seq)
        handle.close()
        if seqcount == 0:
            verbalise("R", "No gene sequences were found.")
            return {}
        # create alignment of input sequences:
        mafft_align1 = os.path.join(temp_dir, "mafft_align_input.fa")
        mafft_align(hmminput, mafft_align1)

        verbalise("B", "Creating hidden markov model from %d sequences" % seqcount)
        # create hmmbuild model of alignment:
        hmmmodel = os.path.join(temp_dir, "hmmmodel.fa")
        open(hmmmodel, 'a').close()
        handle = os.popen(" ".join(['hmmbuild --informat afa', hmmmodel, mafft_align1]))
        handle.close()

        homologlist = hmmer_search(None,
                                    specieslist,
                                    query_species=species,
                                    minthresh=localthresh,
                                    temp_dir=temp_dir,
                                    dbpaths=dbpaths,
                                    mincollect=mincollect,
                                    globalthresh=globalthresh,
                                    hmmfile=hmmmodel,
                                    verbalise=verbalise,
                                    extra_file_search=extra_file_search)

        os.remove(mafft_align1)
        os.remove(hmminput)

    else:
        verbalise("B", "Extracting sequence from %s" % genes)
        if not isinstance(genes, list):
            genes = [genes]
        # run phmmer on a single input gene/sequence:
        seq = ""
        for defline, seq, species in internal.get_gene_fastas(genes=genes,
                                                    species=species,
                                                    fastafile=fastafile,
                                                    specieslist=specieslist,
                                                    dbpaths=dbpaths):

            fasta_seq = "%s\n%s\n" % (defline, seq)
            verbalise("C", fasta_seq)

        if not seq:
            verbalise("R", "No genes sequences were found.")
            return {}
        ## phmmer all lpep files
        homologlist = hmmer_search(fasta_seq,
                                    specieslist,
                                    query_species=species,
                                    minthresh=localthresh,
                                    dbpaths=dbpaths,
                                    temp_dir=temp_dir,
                                    mincollect=mincollect,
                                    globalthresh=globalthresh,
                                    hmmfile=None,
                                    verbalise=verbalise,
                                    extra_file_search=extra_file_search)

    return homologlist

def hmmer_search(fasta_seq, specieslist, query_species,  temp_dir, dbpaths={},
                    minthresh=0.8, mincollect=2, globalthresh=0.01, hmmfile=None,
                    verbalise=lambda *a: None, extra_file_search=None):
    """
    HMMER search of longest non-redundant peptide fasta files using either phmmer (with
    query protein as input) or hmmersearch (using constructed hmmfile as input).
    Finds best score, and collects all proteins with score > PC% (default 80%) of the
    best score, but limited to no less than X (default 2) proteins per species.
    """

    # set search type and input files:
    if hmmfile:
        searchcmd = 'hmmsearch'
        hmminput = hmmfile
    else:
        searchcmd = 'phmmer'
        hmminput = os.path.join(temp_dir,"seq.fasta")
        handle = open(hmminput, 'w')
        handle.write(fasta_seq)
        handle.close()

    verbalise("B", "Finding similar sequences using %s" % searchcmd)
    homologlist = {}
    all_results = {}
    filtered_results = {}
    has_bestscore = False
    if query_species:
        search_space = [query_species] + [ s for s in specieslist if s != query_species ]
    else:
        search_space = specieslist

    if extra_file_search:
        search_space += [ extra_file_search ]

    for sp in search_space:
        try:
            phandle = os.popen( " ".join([searchcmd, hmminput, dbpaths[sp + '_lpep']]) )
        except KeyError:
            try:
                phandle = os.popen( " ".join([ searchcmd, hmminput, sp ]) )
            except:
                print "error:"
                print sp + '_lpep', "not found"
                continue

        # parse phmmer results:
        all_results[sp] = internal.parse_the_hmmer(phandle)

        # determine cutoff threshold for filtering results:
        if sp == query_species:
            bestscore = max( v[1] for v in all_results[sp].values())
            has_bestscore = True
        elif has_bestscore:
            pass
        else:
            score_list = [ v[1] for v in all_results[sp].values() ]
            if len(score_list) > 0:
                bestscore = max( score_list )
            else:
                bestscore = 0
        cutoff_thresh = minthresh * bestscore

        # filter local file based on parameters given:
        ars = all_results[sp]
        filtered_results[sp] = { gene:(sp, ars[gene][1], ars[gene][0]) for i,gene in enumerate(ars) if ars[gene][1] >= cutoff_thresh or i < mincollect}

    # filter for global threshold (ie, based on % of best match). Most useful if no
    # species has been specified for fasta file.
    if filtered_results == {}:
        return {}
    try:
        bestscore = max( v[1] for s in filtered_results for v in filtered_results[s].values())
    except ValueError:
        bestscore = 0
    global_thresh = globalthresh * bestscore
    homologlist = { k:v for nd in filtered_results.values() for (k,v) in nd.items() if v[1] >= global_thresh }

    # clean up temporary files
    os.remove(hmminput)

    return homologlist

####### RAxML functions ########
def raxml(logfile, fastafile, bootstrap=False, threads=2,
            name_conversion=None, phylipize=False,
            verbalise=lambda *a: None, ):

    if phylipize:
        alignment = internal.make_phylip(fastafile, logfile)
    else:
        alignment = fastafile

    print "phylip done"
    if bootstrap:
        raxml_best = raxml_phylogeny(alignment,
                                                logfile,
                                                bootstrap=False,
                                                threads=threads)
        raxml_bstrap = raxml_phylogeny(alignment,
                                                logfile,
                                                bootstrap=bootstrap,
                                                threads=threads)
        raxml_final = apply_boostrap(raxml_best, raxml_bstrap, logfile)

        verbalise("Y",
                "Best tree with bootstrap support can be found at %s" % raxml_final)

    else:
        print "Starting raxml"
        raxml_final = raxml_phylogeny(alignment,
                                        logfile,
                                        bootstrap=bootstrap,
                                        threads=threads)
        verbalise("Y", "Best tree can be found at %s" % raxml_final)
    if name_conversion:
        handle = open(name_conversion, 'rb')
        conv_dic = { line.split()[0]:(line.split()[2], line.split()[1]) for line in handle }
        handle.close()
        internal.rename_newick(raxml_final, conversiondic=conv_dic)
    return raxml_final


def raxml_phylogeny(phylip_alignment, logfile, bootstrap=False, threads=2):
    print "raxml_phylogeny"
    if threads < 2:
        verbalise("R", "Minimum number of threads allowed is 2")
        exit()
    if bootstrap:
        bootstrapopt = '-N 100 -x ' + str(bootstrap)
        bs_str = 'bs.'
    else:
        bootstrapopt = '-N 1'
        bs_str = ""
    print "Bootstrap: ", bootstrap
    # determine path of final RAxML output file
    raxml_outfile = os.path.basename(logfile[:-3] + bs_str + 'raxml.out')
    if bootstrap:
        prefix = "RAxML_bootstrap."
    else:
        prefix = "RAxML_bestTree."
    raxml_final = os.path.join(os.path.dirname(logfile), prefix + raxml_outfile)

    cmd = " ".join(["RAxML", "-s", phylip_alignment,
                                    '-T', str(threads),
                                    '-w', os.path.dirname(logfile),
                                    "-n", raxml_outfile,
                                    '-p', "12345",
                                    "-m", 'PROTGAMMALG',
                                    bootstrapopt,
                                    ])
    os.system( cmd )
    return raxml_final

def apply_boostrap(besttree, bootstraps, logfile):
    """
    Takes the bootstrap trees and the best tree and puts the bootstrap support onto the
    best tree.
    """
    outfile = os.path.basename(logfile[:-3] + 'raxml_merged.nwk')
    cmd = " ".join(["raxmlHPC",     '-z', bootstraps,
                                    '-w', os.path.dirname(logfile),
                                    '-n', outfile,
                                    '-m', 'PROTGAMMALG',
                                    '-t', besttree,
                                    '-f', 'b',
                                    ])
    os.system(cmd)
    return os.path.join(os.path.dirname(logfile), "RAxML_bipartitions." + outfile)

def main():
    pass

if __name__ == '__main__':
    main()

