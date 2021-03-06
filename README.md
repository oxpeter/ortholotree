# ortholotree
search and evaluate gene orthologs using hmm and phylogeny

ortholotree was written to first and foremost look for gene orthologs. Its functionality
would suggest that it can do more than this, but this is not its intent by design. It
can be run directly through CLI, and the options are hopefully explained sufficiently
in the --help option to allow relatively straightforward usage. More complete documents
may follow in time, if there is enough demand...

To install, change into the downloaded directory and type:

python setup.py install

You may need sudo privileges depending on your system.

This will install a config file that can be automatically generated by running config.py,
or by manually populating it with file paths. The format is simply two columns: KEY  PATH
The necessary keys are determined by the modules, and therefore cannot be arbitrarily set.
However, if you initialise the file with config.py, it will put all the necessary KEYs in
the file, even if it cannot find the appropriate file path.
# prerequisites
orthotree relies on the
*hmmer (for protein searching - http://hmmer.janelia.org)
*mafft (for sequence alignment - http://mafft.cbrc.jp/alignment/software/)
*RAxML (for phylogeny construction - http://sco.h-its.org/exelixis/web/software/raxml)

you will also need to install the following python libraries:
*Biopython
*Matplotlib
*Numpy

# usage
ortholotree will work best with a couple of iterations of usage, depending on your
starting point.

__1st pass__

If you wish to find the ortholog of a particular gene, then you will
start by providing the peptide ID to the -g flag, or the link to a fasta file
containing a single sequence using the -f flag. This will run the phmmer search on all
protein databases using the query peptide. It will then filter the results according to
three parameters:

*Local threshold (-t) For each species searched, collect all matches that are above this
value. If the species of the query sequence is provided, and a peptide database exists
for this species, then the best match is the score of the query matched to itself. All
other scores are therefore relative to the best score. If the species is either unknown
or not in a db, then the best score is the highest score returned for each peptide db
searched.

*Minimum collection (-c) For each species, the best n results will be collected, regardless
of whether they pass the local threshold. If you only want genes that pass the local
threshold, then pass a value of 0 to this flag.

*Global threshold (-a) Once all genes passing the previous two filters are collated, only
genes with a score that is above the specified percentage of the best score will be kept.

__2nd pass__

Using the results of the first pass, select likely ortholog candidates (this may simply
be the gene in your target species, as well as your initial query gene). Pass all genes
as a comma-separated list to the -g flag, and specify -B to tell the program to build a
hidden markov model using an alignment of the specified proteins.

The various peptide databases will then be searched using the hmm, providing a more
sensitive search. You can filter the hits using the same parameter flags as in the first
pass, but typically you can use higher percentage values (as the orthologs will generally
be a better match to the model than to a single peptide, particularly from highly
divergent species.

BE WARNED! If you provide non-orthologous genes to the model, you will probably pull out
multiple gene families, giving too many genes, or not enough.


__3rd pass__

This final pass will often not be necessary, but if you want to tidy up your results,
then use the tree built by the 2nd pass to identify the orthologs from all species
available, and use all identified orthologs to build the hmm model as specified above.
You can also use some additional flags to help remove unwanted hits:

*Exclude gene (-x) a comma-separated list of genes that you do not want in your alignment
or phylogeny.

*Exclude species (-e) like the -x flag, only all genes from a species will be excluded.
Species names are specified using the four-letter key used to construct the peptide
database paths (E.g. Cerapachys biroi is Cbir).

*Maximum length (-l) to prevent gene fusions from being picked up and causing havoc with
your alignment and phylogeny, you can specify a length (in amino acids), for which all
proteins must be smaller than to remain.

If you would like more information on how robust your phylogeny is, use the -b flag to
run a bootstrap analysis. This will produce 100 bootstrapped trees, then show the
bootstrap support on the (independently calculated) most likely tree.





