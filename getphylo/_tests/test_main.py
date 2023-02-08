import os
import unittest
from unittest.mock import patch
from io import StringIO
#from tempfile import TemporaryDirectory

from getphylo.utils import io as gp_io
from getphylo import screen
from getphylo.align import (
    get_locus_length,
    make_fasta_for_alignments,
    get_locus_from_tsv
)

class TestAlignmentLength(unittest.TestCase):
    def test_get_locus_from_tsv(self):
        mock_tsv = [
            ['locus1','target_locus2', 100],
            ['locus2','off_target_locus2', 100],
            ['locus3','off_target_locus2', 100]
        ]
        mock_tsv2 = [
            ['locus1','target_locus8', 100],
            ['locus2','off_target_locus7', 100],
            ['locus3','off_target_locus9', 100]
        ]
        with patch.object(screen, 'get_locus', return_value='AAAAAA') as patched_read:
            with patch.object(gp_io, 'read_tsv', return_value= mock_tsv):
                query = 'locus1'
                assert get_locus_from_tsv(query,'dummy.fasta') == ('>dummy_target_locus2', 'AAAAAA')
                patched_read.assert_called_once_with('dummy.fasta', 'target_locus2')

            with patch.object(gp_io, 'read_tsv', return_value=mock_tsv2):
                assert get_locus_from_tsv('locus2','dummy.fasta') == ('>dummy_off_target_locus7', 'AAAAAA')
                with self.assertRaisesRegex(ValueError, 'not found'):
                    assert get_locus_from_tsv('bad','dummy.fasta')

    def test_get_locus_length(self):
        assert get_locus_length(['xxxxx','yyyy']) == 4
        assert get_locus_length(['>my_sequence','MYSEQENCE']) == 9
        assert get_locus_length(['>my_sequence','MYSEQENCE\n']) == 9
        with self.assertRaisesRegex(ValueError, 'empty'):
                get_locus_length([]) 
        #assert get_locus_length(['MYSEQENCE', '>my_sequence']) == 9
        #assert get_locus_length('>my_sequence') == 9
        #with patch_open(return_value=StringIO(">bob\nstu\nff\n")):
            #read_fasta("thing")

        ###THINK OF TESTS!

#def test_do_alignments(self):
    
#def test_get_locus_alignment(self):

#def test_make_combined_alignment(self):
    