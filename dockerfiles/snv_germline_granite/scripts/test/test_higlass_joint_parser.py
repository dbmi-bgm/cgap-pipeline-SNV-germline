
#################################################################
#   Libraries
#################################################################
import sys, os
import pytest

from higlass_joint_parser import (
                            main as main_higlass_joint_parser
                           )


#################################################################
#   Tests
#################################################################
probandlist = ['test_1', 'test_2', 'test_3', 'test_4', 'test_5', 'test_6', 'test_7', 'test_8', 
'test_9', 'test_10', 'test_11', 'test_12', 'test_13', 'test_14', 'test_15', 'test_16', 
'test_17', 'test_18', 'test_19', 'test_20', 'test_21', 'test_22', 'test_23', 'test_24', 
'test_25', 'test_26', 'test_27', 'test_28', 'test_32', 'test_36', 'test_40', 'test_51', 
'test_52', 'test_58', 'test_60', 'test_64', 'test_68', 'test_72']

def test_joint_parser_v3_v2_together():
    # Variables and Run
    args = {'inputfile': 'test/files/test_in_joint.vcf.gz', 'probandlist': probandlist, 'outputfile':'output.vcf', 'gnomAD':['v3', 'v2']}
    # Test
    main_higlass_joint_parser(args)
    a = os.popen('bgzip -c -d output.vcf.gz')
    b = os.popen('bgzip -c -d test/files/test_out_joint.vcf.gz')

    assert [row for row in a.read()] == [row for row in b.read()]

    # Clean
    os.remove('output.vcf.gz')
    os.remove('output.vcf.gz.tbi')

def test_joint_parser_v3():
    # Variables and Run
    args = {'inputfile': 'test/files/test_in_joint.vcf.gz', 'probandlist': probandlist, 'outputfile':'output.vcf', 'gnomAD':['v3']}
    # Test
    main_higlass_joint_parser(args)
    a = os.popen('bgzip -c -d output.vcf.gz')
    b = os.popen('bgzip -c -d test/files/test_out_joint_v3.vcf.gz')

    assert [row for row in a.read()] == [row for row in b.read()]

    # Clean
    os.remove('output.vcf.gz')
    os.remove('output.vcf.gz.tbi')

def test_joint_parser_v2():
    # Variables and Run
    args = {'inputfile': 'test/files/test_in_joint.vcf.gz', 'probandlist': probandlist, 'outputfile':'output.vcf', 'gnomAD':['v2']}
    # Test
    main_higlass_joint_parser(args)
    a = os.popen('bgzip -c -d output.vcf.gz')
    b = os.popen('bgzip -c -d test/files/test_out_joint_v2.vcf.gz')

    assert [row for row in a.read()] == [row for row in b.read()]

    # Clean
    os.remove('output.vcf.gz')
    os.remove('output.vcf.gz.tbi')
