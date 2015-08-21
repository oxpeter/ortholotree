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
            self.assertEqual(defline, '>NP_001035293.1 (Amel) None')
            self.assertEqual(seq,
'MPILIPHRNPASANYYENKDGARIVKASHFELDYMLGRKITFFCMATGFPRPEITWLKDGIELYHHKFFQVHEWPVGNDTLKSKMEIDPATQKDAGYYECQADNQYAVDRRGFRTDYVMISY')
            self.assertEqual(species, 'Amel')

    def test_nogenelist(self):
        for defline, seq, species in internal.get_gene_fastas(genes='NP_001035293.1',
                                                        dbpaths=self.dbpaths,
                                                        specieslist=self.specieslist):
            self.assertEqual(defline, '>NP_001035293.1 (Amel) None')
            self.assertEqual(seq,
'MPILIPHRNPASANYYENKDGARIVKASHFELDYMLGRKITFFCMATGFPRPEITWLKDGIELYHHKFFQVHEWPVGNDTLKSKMEIDPATQKDAGYYECQADNQYAVDRRGFRTDYVMISY')
            self.assertEqual(species, 'Amel')

    def test_incompletegenenamebeginning(self):
        for defline, seq, species in internal.get_gene_fastas(genes=['001035293.1'],
                                                        dbpaths=self.dbpaths,
                                                        specieslist=self.specieslist):
            self.assertEqual(defline, None)
            self.assertEqual(seq, None)
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
            self.assertEqual(defline, '>NP_001035293.1 (Amel) None')
            self.assertEqual(seq,
'MPILIPHRNPASANYYENKDGARIVKASHFELDYMLGRKITFFCMATGFPRPEITWLKDGIELYHHKFFQVHEWPVGNDTLKSKMEIDPATQKDAGYYECQADNQYAVDRRGFRTDYVMISY')
            self.assertEqual(species, 'Amel')


class TestTrimNameDross(unittest.TestCase):
    def test_namewithpipe_trimmed(self):
        self.assertEqual(internal.trim_name_dross('Cbir|LOC12345'),
                                'LOC12345')

    def test_namewithoutpipe_maintained(self):
        self.assertEqual(internal.trim_name_dross('CbirLOC12345'),
                                'CbirLOC12345')


class TestPCmatch(unittest.TestCase):
    def test_pcmatch100(self):
        self.assertEqual(internal.get_pcmatch('ABCDEFGHIJ'),
                            (10,10))
    def test_pcmatch0(self):
        self.assertEqual(internal.get_pcmatch('----------'),
                            (10,0))
    def test_pcmatch50(self):
        self.assertEqual(internal.get_pcmatch('ABCDE-----'),
                            (10,5))
        self.assertEqual(internal.get_pcmatch('A-B-C-D-E-'),
                            (10,5))



if __name__ == '__main__':
    unittest.main()