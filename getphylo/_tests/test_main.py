import glob
import os
import unittest
from unittest.mock import patch
from io import StringIO

from getphylo.main import check_seed
from getphylo.utils.checkpoint import Checkpoint
from getphylo.utils.errors import BadSeedError

class TestCheckSeed(unittest.TestCase):
    def test_check_seed(self):
        '''Basic test for check_seed()
            Arguments: Self
            Returns: None
        '''
        checkpoint1 = 0
        search_string1 = 'dummy'
        glob1 = ['one.gbk','two.gbk','three.gbk']
        with patch.object(glob, 'glob', return_value=glob1):
            assert check_seed(checkpoint1, search_string1) == 'one.gbk'
        
    def test_check_seed_checkpoints(self):
        '''
        Test that checkpoint inputs are triggering correct errors or results
            Arguments: Self
            Returns: None
        '''
        checkpoint1 = 0
        checkpoint2 = Checkpoint.START
        checkpoint3 = Checkpoint.SINGLETONS_THRESHOLDED
        checkpoint4 = 500
        checkpoint5 = 'start'
        checkpoint6 = [Checkpoint.START, Checkpoint.SINGLETONS_THRESHOLDED]
        search_string1 = 'dummy'
        glob1 = ['one.gbk','two.gbk','three.gbk']
        with patch.object(glob, 'glob', return_value=glob1):
            assert check_seed(checkpoint1, search_string1) == 'one.gbk'
        with patch.object(glob, 'glob', return_value=glob1):
            assert check_seed(checkpoint2, search_string1) == 'one.gbk'
        with patch.object(glob, 'glob', return_value=glob1):
            with self.assertRaisesRegex(BadSeedError, 'checkpoint has been set'):
                check_seed(checkpoint3, search_string1) == 'one.gbk'
        with patch.object(glob, 'glob', return_value=glob1):
            with self.assertRaisesRegex(BadSeedError, 'checkpoint has been set'):
                check_seed(checkpoint4, search_string1) == 'one.gbk'
        with patch.object(glob, 'glob', return_value=glob1):
            with self.assertRaises(TypeError):
                check_seed(checkpoint5, search_string1) == 'one.gbk'
        with patch.object(glob, 'glob', return_value=glob1):
            with self.assertRaises(TypeError):
                check_seed(checkpoint6, search_string1) == 'one.gbk'

    def test_check_seed_glob(self):
        '''
        Test that an empty glob will trigger a BadSeedError
            Arguments: Self
            Returns: None
        '''
        checkpoint1 = 0
        search_string1 = 'dummy'
        glob2 = []
        with patch.object(glob, 'glob', return_value=glob2):
            with self.assertRaisesRegex(BadSeedError, 'No files found'):
                check_seed(checkpoint1, search_string1) 
