#!/usr/bin/env python
"""
This is the script containing all the test code for the ortholotree package.
"""

import unittest
from pkg_resources import resource_filename, resource_exists, resource_stream

from orthomods import internal, external

class TestUnitTest(unittest.TestCase):
    def test_truepasses(self):
        self.assertTrue(True)

    def test_falsepasses(self):
        self.assertFalse(False)

    def test_equalpasses(self):
        self.assertEqual(1,1)
        self.assertEqual('a', 'a')
        self.assertEqual([1,'a'], [1,'a'])

    def test_typeerror(self):
        with self.assertRaises(TypeError):
            " ".split(2)

class FastaTestCase(unittest.TestCase):
    def setUp(self):
        self.fastafile = resource_filename('testcode', 'fasta_eg1.txt')
        self.fastafileillegal = resource_filename('testcode', 'fasta_eg3.txt')
        self.fastafiletwo = resource_filename('testcode', 'fasta_eg2.txt')
        self.fastaalignment = resource_filename('testcode', 'fasta_alignment_eg5.txt')
        self.dbpaths = {'Amel_lpep':resource_filename('testcode', 'Amel_longest_pep.txt'),
                        'Cbir_lpep':resource_filename('testcode', 'Cbir_longest_pep.txt')}
        self.specieslist = ['Cbir', 'Amel']

class TestGetGeneFastas(FastaTestCase):
    def test_no_argsgiven(self):
        for defline, seq, species in internal.get_gene_fastas():
            self.assertEqual(defline, None)
            self.assertEqual(seq, None)
            self.assertEqual(species, None)

    def test_onlyfastafilegiven(self):
        for defline, seq, species in internal.get_gene_fastas(fastafile=self.fastafile):
            self.assertEqual(defline, '>Cbir|LOC12345 testgene')
            self.assertEqual(seq,
             'ABCDEFGABCDEFGABCDEFGABCDEFGABCDEFGhijklmnopHIJKLMhijklmnopHIJKLMNOPQRSTUV')
            self.assertEqual(species, None)

    def test_onlygenegiven(self):
        for defline, seq, species in internal.get_gene_fastas(genes=['NP_001035293.1']):
            self.assertEqual(defline, None)
            self.assertEqual(seq, None)
            self.assertEqual(species, None)

    def test_dbpath_no_specieslist(self):
        for defline, seq, species in internal.get_gene_fastas(genes=['NP_001035293.1'],
                                                        dbpaths=self.dbpaths):
            self.assertEqual(defline, None)
            self.assertEqual(seq, None)
            self.assertEqual(species, None)

    def test_geneand_dbpaths(self):
        for defline, seq, species in internal.get_gene_fastas(genes=['NP_001035293.1'],
                                                        dbpaths=self.dbpaths,
                                                        specieslist=self.specieslist):
            self.assertEqual(defline, '>Amel|NP_001035293.1')
            self.assertEqual(seq,
'MPILIPHRNPASANYYENKDGARIVKASHFELDYMLGRKITFFCMATGFPRPEITWLKDGIELYHHKFFQVHEWPVGNDTLKSKMEIDPATQKDAGYYECQADNQYAVDRRGFRTDYVMISY')
            self.assertEqual(species, None)

    def test_nogenelist(self):
        for defline, seq, species in internal.get_gene_fastas(genes='NP_001035293.1',
                                                        dbpaths=self.dbpaths,
                                                        specieslist=self.specieslist):
            self.assertEqual(defline, '>Amel|NP_001035293.1')
            self.assertEqual(seq,
'MPILIPHRNPASANYYENKDGARIVKASHFELDYMLGRKITFFCMATGFPRPEITWLKDGIELYHHKFFQVHEWPVGNDTLKSKMEIDPATQKDAGYYECQADNQYAVDRRGFRTDYVMISY')
            self.assertEqual(species, None)

    def test_duplicate_matches(self):
        for defline, seq, species in internal.get_gene_fastas(genes=['XP','XP_006570708.1'],
                                                        dbpaths=self.dbpaths,
                                                        specieslist=self.specieslist):
            self.assertEqual(defline, '>Amel|XP_006570708.1')
            self.assertEqual(seq,
            'MEIAASAMLDGLKNNRISKLALSRFLSQSLVSCILLGLLLEFRAQLETTGSPANKPASSASGSGTGGTSGTSVGANNLTTSGIATGSSGSGSGATVSGIGTVNAGSSGINIGANVGTGNVASGTVESRTSTIGVQNKQLQNVKGEHPSKAFLQNRSMSLVDMYIDNSEPSENVGQIHFSLEYDFQNTTLILRIIQGKDLPAKDLSGTSDPYVRVTLLPDKKHRLETKIKRRTLNPRWNETFYFEGFPIQKLQSRVLHLHVFDYDRFSRDDSIGEMFLPLCQVDFSDKPSFWKALKPPAKDKCGELLCSLCYHPSNSVLTLTLLKARNLKAKDINGKSDPYVKVWLQFGDKRIEKRKTPIFKCTLNPVFNEAFSFNVPWEKIRECSLDVMVMDFDNIGRNELIGRIQLAGKNGSGASETKHWQDMITKPRQTIVQWHRLKPE' )
            self.assertEqual(species, None)

    def test_incompletegenenameend(self):
        for defline, seq, species in internal.get_gene_fastas(genes=['NP_0010352'],
                                                        dbpaths=self.dbpaths,
                                                        specieslist=self.specieslist):
            self.assertEqual(defline, None)
            self.assertEqual(seq, None)
            self.assertEqual(species, None)

    def test_badfastaseq(self):
        for defline, seq, species in internal.get_gene_fastas(fastafile=self.fastafileillegal):
            self.assertEqual(defline, None)
            self.assertEqual(seq, None)
            self.assertEqual(species, None)

    def test_genewithspecies(self):
        for defline, seq, species in internal.get_gene_fastas(genes='NP_001035293.1',
                                                        dbpaths=self.dbpaths,
                                                        specieslist=self.specieslist,
                                                        species='Amel'):
            self.assertEqual(defline, '>Amel|NP_001035293.1')
            self.assertEqual(seq,
'MPILIPHRNPASANYYENKDGARIVKASHFELDYMLGRKITFFCMATGFPRPEITWLKDGIELYHHKFFQVHEWPVGNDTLKSKMEIDPATQKDAGYYECQADNQYAVDRRGFRTDYVMISY')
            self.assertEqual(species, None)

class TestAlignments(FastaTestCase):
    def test_consensus_creation(self):
        consensus = internal.Consensus(self.fastaalignment)
        self.assertEqual(consensus.maxlen, 7)
        self.assertEqual(consensus.consensus_pc[0], 0.5)
        self.assertEqual(consensus.consensus_pc[1], 1.0)
        self.assertEqual(consensus.consensus_pc[2], 0.5)
        self.assertEqual(consensus.consensus_pc[4], 0.1)
        self.assertTrue('A' in consensus.consensus[0])
        self.assertFalse('H' in consensus.consensus[0])
        self.assertFalse('I' in consensus.consensus[0])
        self.assertTrue('C' in consensus.consensus[2])
        self.assertTrue('D' in consensus.consensus[2])
        self.assertEqual(len(consensus.consensus[2]),2)
        self.assertEqual(len(consensus.consensus[0]),1)

class TestFixLeakyPipes(unittest.TestCase):
    def test_namewithpipe_fixed(self):
        self.assertEqual(internal.fix_leaky_pipes('Cbir|LOC12345'),
                                'Cbir\|LOC12345')

    def test_namewithoutpipe_maintained(self):
        self.assertEqual(internal.fix_leaky_pipes('CbirLOC12345'),
                                'CbirLOC12345')

class TestPCmatch(unittest.TestCase):
    def test_pcmatch100(self):
        self.assertEqual(internal.get_pcmatch('ABCDEFGHIJ'),
                            (10,1))
    def test_pcmatch0(self):
        self.assertEqual(internal.get_pcmatch('----------'),
                            (10,0))
    def test_pcmatch50(self):
        self.assertEqual(internal.get_pcmatch('ABCDE-----'),
                            (10,0.5))
        self.assertEqual(internal.get_pcmatch('A-B-C-D-E-'),
                            (10,0.5))

class HMMerTestCase(unittest.TestCase):
    def setUp(self):
        self.hmmerfile = resource_filename('testcode', 'phmmer_out.txt')

class TestHmmerParsing(HMMerTestCase):
    def test_make_instance(self):
        with open(self.hmmerfile) as handle:
            hmmer = internal.HMMer(handle)
        self.assertTrue(hmmer)

    def test_hmmer_stats(self):
        with open(self.hmmerfile) as handle:
            hmmer = internal.HMMer(handle)

        self.assertEqual(hmmer.query, 'Amel|AST')
        self.assertEqual(len(hmmer.targets), 6)
        self.assertEqual(hmmer.targets[0], 'sp|P85797|ALLS_APIME')
        self.assertEqual(len(hmmer.domain_seq), 10)
        self.assertEqual(len(hmmer.domain_prb), 10)
        self.assertEqual(len(hmmer.domain_stats), 10)
        self.assertEqual(len(hmmer.stats), 6)

    def test_hmmer_sequences(self):
        with open(self.hmmerfile) as handle:
            hmmer = internal.HMMer(handle)

        dmel_dom1 = "GGDNIDK-R-VERYAFGLGRRAYMYTNGGP-GMK"
        dmel_dom2 = "KRVERYAFGLGRRaYMYTNGGpgmKRLPVYNFGLGKRsRPYSFGLGKRSD-YDYDQDNEIDYRvPPAN-YLAAERAVRPGRQNKRTTrpQPFNFGLGRR"
        self.assertEqual(hmmer.domain_seq['sp|Q9VC44|ALLS_DROME'][1],dmel_dom1)
        self.assertEqual(hmmer.domain_seq['sp|Q9VC44|ALLS_DROME'][2],dmel_dom2)

        dmel_prob1 = "4454443.3.3569***********99998.444"
        dmel_prob2 = "78999******985678886433478889******99678********98.59999999**9624444.3445555566667776532389******97"
        self.assertEqual(hmmer.domain_prb['sp|Q9VC44|ALLS_DROME'][1],dmel_prob1)
        self.assertEqual(hmmer.domain_prb['sp|Q9VC44|ALLS_DROME'][2],dmel_prob2)

        dmel_desc = "Allatostatin-A OS=Drosophila melanogas"
        self.assertEqual(hmmer.stats['sp|Q9VC44|ALLS_DROME']['eval'], 1.5e-07)
        self.assertEqual(hmmer.stats['sp|Q9VC44|ALLS_DROME']['score'], 36.7)
        self.assertEqual(hmmer.stats['sp|Q9VC44|ALLS_DROME']['bias'], 7.3)
        self.assertEqual(hmmer.stats['sp|Q9VC44|ALLS_DROME']['dom_e'], 1.5e-06)
        self.assertEqual(hmmer.stats['sp|Q9VC44|ALLS_DROME']['dom_s'], 33.4)
        self.assertEqual(hmmer.stats['sp|Q9VC44|ALLS_DROME']['dom_b'], 6.6)
        self.assertEqual(hmmer.stats['sp|Q9VC44|ALLS_DROME']['dom_exp'], 2.2)
        self.assertEqual(hmmer.stats['sp|Q9VC44|ALLS_DROME']['dom_no'], 2)
        self.assertEqual(hmmer.stats['sp|Q9VC44|ALLS_DROME']['desc'], dmel_desc)

if __name__ == '__main__':
    unittest.main()