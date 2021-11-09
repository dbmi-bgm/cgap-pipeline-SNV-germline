#################################################################
#   Libraries
#################################################################
import sys, os
import pytest
from subprocess import Popen

from liftover_hg19 import (
                            main as main_liftover_hg19
                           )


#################################################################
#   Tests
#################################################################

def test_full_process():
    # Variables and Run
    args = {'inputfile': 'test/files/test_reduced_sorted','outputfile':'output.vcf','chainfile':'test/files/test.chain'}
    # Test
    main_liftover_hg19(args)
    a = os.popen('bgzip -c -d output.vcf.gz')
    b = os.popen('bgzip -c -d test/files/test_reduced_lo_sorted.vcf.gz')

    assert [row for row in a.read()] == [row for row in b.read()]
    # Clean
    os.remove('output.vcf.gz')
    os.remove('output.vcf.gz.tbi')
