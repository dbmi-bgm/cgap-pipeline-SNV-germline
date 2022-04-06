#################################################################
#   Libraries
#################################################################
import sys, os
import pytest

from dbSNP_ID_fixer import (
                            main as main_dbSNP_ID_fixer
                           )

#################################################################
#   Tests
#################################################################

def test_dbSNP_ID_fixer_region1():
    # Variables
    args = {'dbSNPvcf': 'test/files/dbSNPfix_dbSNP151_data_test.vcf.gz', 'inputvcf': 'test/files/dbSNPfix_sample_data_test.vcf.gz', 'regionfile': 'chr19:1-350000'}
    # Run
    main_dbSNP_ID_fixer(args)
    # Tests
    assert [row for row in open('chr19:1-350000')] == [row for row in open('test/files/dbSNPfix_region1')]
    # Clean
    os.remove('chr19:1-350000')

def test_dbSNP_ID_fixer_region2():
    # Variables
    args = {'dbSNPvcf': 'test/files/dbSNPfix_dbSNP151_data_test.vcf.gz', 'inputvcf': 'test/files/dbSNPfix_sample_data_test.vcf.gz', 'regionfile': 'chr19:350001-20000000'}
    # Run
    main_dbSNP_ID_fixer(args)
    # Tests
    assert [row for row in open('chr19:350001-20000000')] == [row for row in open('test/files/dbSNPfix_region2')]
    # Clean
    os.remove('chr19:350001-20000000')
