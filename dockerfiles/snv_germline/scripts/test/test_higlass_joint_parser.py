
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

def test_joint_parser_v3_v2_together():
    # Variables and Run
    args = {'inputfile': 'test/files/test_in_joint.vcf.gz', 'probandlist': 'test/files/test_id_list.txt', 'outputfile':'output.vcf', 'gnomAD':['v3', 'v2']}
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
    args = {'inputfile': 'test/files/test_in_joint.vcf.gz', 'probandlist': 'test/files/test_id_list.txt', 'outputfile':'output.vcf', 'gnomAD':['v3']}
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
    args = {'inputfile': 'test/files/test_in_joint.vcf.gz', 'probandlist': 'test/files/test_id_list.txt', 'outputfile':'output.vcf', 'gnomAD':['v2']}
    # Test
    main_higlass_joint_parser(args)
    a = os.popen('bgzip -c -d output.vcf.gz')
    b = os.popen('bgzip -c -d test/files/test_out_joint_v2.vcf.gz')

    assert [row for row in a.read()] == [row for row in b.read()]

    # Clean
    os.remove('output.vcf.gz')
    os.remove('output.vcf.gz.tbi')
