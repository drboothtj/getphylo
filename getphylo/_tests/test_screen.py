import glob
import os
import unittest
from unittest.mock import patch
from io import StringIO

from getphylo.screen import write_pa_table
from getphylo.utils import io as gp_io


class TestWritePATable(unittest.TestCase):
    def test_write_get_pa_table(self):
        pa_data = [['strain1', 'strain2','strain3','strain4'],[0,1,1,0],[2,1,1,1],[0,1,0,1]]
        loci = ['locus1','locus2','locus3']
        output = 'output'
        expected_filename = 'output/presence_absence_table.csv'
        expected_result = (
            ['strain;locus1;locus2;locus3','strain1;0;2;0','strain2;1;1;1','strain3;1;1;0','strain4;0;1;1']
        )

        with patch.object(gp_io, 'write_to_file') as patched_read:
            write_pa_table(pa_data, loci, output)
            patched_read.assert_called_once_with(expected_filename, expected_result)

#test (main) checkpoint is correct -> add to checkpoint check to io?

# assert the presence of required files and suggest a different checkpoin







