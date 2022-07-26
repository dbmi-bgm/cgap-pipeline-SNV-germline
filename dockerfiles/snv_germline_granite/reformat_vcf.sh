#!/bin/bash

# variables from command line
input_vcf=$1
output_vcf=$2

# run
portal_reformat_vcf.py -i $input_vcf -o $output_vcf || exit 1

# compress and index output vcf
bgzip $output_vcf || exit 1
tabix -p vcf ${output_vcf}.gz || exit 1

