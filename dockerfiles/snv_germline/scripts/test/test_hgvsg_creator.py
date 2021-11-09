#################################################################
#   Libraries
#################################################################
import sys, os
import pytest
from subprocess import Popen

from hgvsg_creator import (
                            main as main_hgvsg_creator
                           )


#################################################################
#   Tests
#################################################################

def test_full_process():
    # Variables and Run
    args = {'inputfile': 'test/files/test_reduced_sorted','outputfile':'output.vcf','chainfile':'test/files/test.chain'}
    # Test
    main_hgvsg_creator(args)
    a = os.popen('bgzip -c -d output.vcf.gz')
    b = os.popen('bgzip -c -d test/files/test_reduced_sorted_hgvsg.vcf.gz')

    assert [row for row in a.read()] == [row for row in b.read()]
    # Clean
    os.remove('output.vcf.gz')
    os.remove('output.vcf.gz.tbi')


#################################################################
#   Errors
#################################################################

def test_nonstandard_variant():
    # Variables
    args = {'inputfile': 'test/files/test_reduced_sorted_fail','outputfile':'output.vcf'}
    # Run and Tests
    with pytest.raises(Exception, match="Unexpected variant format found. Quitting and deleting output. Variant: chr21\t1000\tTTT\tTT"):
        main_hgvsg_creator(args)
