##TODO##
#add a test for extension2 error in test_change_extension!
#combine into a single utils folder and script
#add test for get_records_from_genbank
#does read_file need a test?

import os
import unittest
import glob
from unittest.mock import patch

from getphylo.utils.io import(
    count_files,
    change_extension,
    get_records_from_genbank,
    make_folder,
    )

class Test_io(unittest.TestCase):
    def test_count_files(self):
        list1 = ['a.gb', 'y.gb', 'z.gb']
        no_list = []
        with patch.object(glob, 'glob', return_value = list1):
            assert count_files('dummy') == 3        
        with patch.object(glob, 'glob', return_value = []):
            assert count_files('dummy') == 0
    
    def test_change_extension(self):
        filename1 = 'test.gbk'
        filename2 = 'test.test.gbk'
        filename3 = 123
        extension1 = 'tsv'
        extension2 = '.tsv'
        #add a test for extension2 error!
        extension3 = 7
        assert change_extension(filename1, extension1) == 'test.tsv'
        assert change_extension(filename2, extension1) == 'test.test.tsv'
        with self.assertRaises(TypeError):
            assert change_extension(filename3, extension1)
        with self.assertRaises(TypeError):
            assert change_extension(filename1, extension3)

    def test_get_records_from_genbank(self):
        #add tests
        print('add me later')

    def test_make_folder(self):
        with patch.object(os.path, 'exists', return_value = True):
            with self.assertRaisesRegex(Exception, 'already exists'):
                make_folder('dummy')
        with patch.object(os.path, 'exists', return_value = False):
            with patch.object(os, 'mkdir') as mkdir:
                make_folder('dummy')
                mkdir.assert_called_once_with('dummy')

    def test_read_file(self):
        #check how to test
        print('do I need a test?')