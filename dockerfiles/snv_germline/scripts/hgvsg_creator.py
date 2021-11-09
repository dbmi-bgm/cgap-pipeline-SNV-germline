#!/usr/bin/env python3

################################################
#
#  Script to add HGVSG field for all variants
#   into the INFO column of sample vcf files
#  Note: Should be run on filtered variant list
#
################################################

################################################
#   Libraries
################################################
from granite.lib import vcf_parser
from os import remove
import sys, argparse, subprocess

################################################
#Dictionary for converting chr to NCBI accession
################################################

# accession numbers for the 23 chromosomes haven't changed between GRCh38 (https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/) and the GRCh38.p13 (https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39/)
NCBI_CONVERT = {}
NCBI_CONVERT["chr1"] = "NC_000001.11"
NCBI_CONVERT["chr2"] = "NC_000002.12"
NCBI_CONVERT["chr3"] = "NC_000003.12"
NCBI_CONVERT["chr4"] = "NC_000004.12"
NCBI_CONVERT["chr5"] = "NC_000005.10"
NCBI_CONVERT["chr6"] = "NC_000006.12"
NCBI_CONVERT["chr7"] = "NC_000007.14"
NCBI_CONVERT["chr8"] = "NC_000008.11"
NCBI_CONVERT["chr9"] = "NC_000009.12"
NCBI_CONVERT["chr10"] = "NC_000010.11"
NCBI_CONVERT["chr11"] = "NC_000011.10"
NCBI_CONVERT["chr12"] = "NC_000012.12"
NCBI_CONVERT["chr13"] = "NC_000013.11"
NCBI_CONVERT["chr14"] = "NC_000014.9"
NCBI_CONVERT["chr15"] = "NC_000015.10"
NCBI_CONVERT["chr16"] = "NC_000016.10"
NCBI_CONVERT["chr17"] = "NC_000017.11"
NCBI_CONVERT["chr18"] = "NC_000018.10"
NCBI_CONVERT["chr19"] = "NC_000019.10"
NCBI_CONVERT["chr20"] = "NC_000020.11"
NCBI_CONVERT["chr21"] = "NC_000021.9"
NCBI_CONVERT["chr22"] = "NC_000022.11"
NCBI_CONVERT["chrX"] = "NC_000023.11"
NCBI_CONVERT["chrY"] = "NC_000024.10"
NCBI_CONVERT["chrM"] = "NC_012920.1"

################################################
#   Functions
################################################

# designed given best practices at http://varnomen.hgvs.org/recommendations/DNA/
def create_hgvsg(CHROM, POS, REF, ALT):
    try:
        # first check that there are no non-standard variants:
        if len(REF) > 1 and len(ALT) > 1 and REF[0] == ALT[0] and REF[1] == ALT[1]:
            HGVSG = 'FAIL'
        else:
            # start by identifying genomic versus mitochondrial
            # note that mito requires some code below for variants that span the end to beginning of the circular DNA (e.g., base 16569 - base 2)
            HGVSG = NCBI_CONVERT[CHROM]
            if HGVSG != "NC_012920.1":
                HGVSG += ':g.'
            else:
                HGVSG += ':m.'
            # deal with dashes as a first case - they only appear in ALT and are all deletions
            # deletions - http://varnomen.hgvs.org/recommendations/DNA/variant/deletion/
            if ALT == "-":
                if len(REF) == 1:
                    HGVSG += str(POS) + 'del'
                else:
                    end_position = POS + len(REF) - 1
                    if ':m.' in HGVSG and end_position > 16569:
                        circular_position = end_position - 16569
                        HGVSG += str(POS) + '_' + str(circular_position) + 'del'
                    else:
                        HGVSG += str(POS) + '_' + str(end_position) + 'del'
            # for all other cases
            else:
                # delins first because they are the weirdest - http://varnomen.hgvs.org/recommendations/DNA/variant/delins/
                if (REF[0] != ALT[0] and len(REF) > 1) or (REF[0] != ALT[0] and len(ALT) > 1):
                    # single-base deletion with multi-base insertion (single to single is just a substitution)
                    if len(REF) == 1:
                        HGVSG += str(POS)+'delins'+ALT
                    # multi-base deletion and any insertion
                    else:
                        end_position = POS + len(REF) - 1
                        if ':m.' in HGVSG and end_position > 16569:
                            circular_position = end_position - 16569
                            HGVSG += str(POS) + '_' + str(circular_position) + 'delins' + ALT
                        else:
                            HGVSG += str(POS) + '_' + str(end_position) + 'delins' + ALT
                else:
                    # substitution - http://varnomen.hgvs.org/recommendations/DNA/variant/substitution/
                    if len(REF) == len(ALT):
                        HGVSG += str(POS) + REF + '>' + ALT
                    # normal deletions have a REF e.g. ATTT and ALT e.g., A so position needs +1
                    elif len(REF) > len(ALT):
                        start_position = POS + 1
                        if len(REF) == 2:
                            if ':m.' in HGVSG and start_position > 16569:
                                circular_position = start_position - 16569
                                HGVSG += str(circular_position) + 'del'
                            else:
                                HGVSG += str(start_position) + 'del'
                        else:
                            end_position = POS + len(REF) - 1
                            if ':m.' in HGVSG and end_position > 16569:
                                if start_position > 16569:
                                    start_position = start_position - 16569
                                    circular_position = start_position + len(REF) - 2
                                    HGVSG += str(start_position) + '_' + str(circular_position) + 'del'
                                else:
                                    circular_position = end_position - 16569
                                    HGVSG += str(start_position) + '_' + str(circular_position) + 'del'
                            else:
                                HGVSG += str(start_position) + '_' + str(end_position) + 'del'
                    # insertions have REF e.g., A and ALT e.g., ATTTTTTT and they want the inserted bases and the positions that they are between
                    # http://varnomen.hgvs.org/recommendations/DNA/variant/insertion/
                    elif len(REF) < len(ALT):
                        if ':m.' in HGVSG and POS == 16569:
                            HGVSG += str(POS) + '_' + str(1) + 'ins' + ALT[1:]
                        else:
                            HGVSG += str(POS) + '_' + str(POS + 1) + 'ins' + ALT[1:]
    except:
        HGVSG = ''
    return HGVSG

def main(args):
    vcf = vcf_parser.Vcf(args['inputfile'])
    hgvsg_definition = '##INFO=<ID=hgvsg,Number=1,Type=String,Description="hgvsg created from variant following best practices - http://varnomen.hgvs.org/recommendations/DNA/">'
    vcf.header.add_tag_definition(hgvsg_definition)

    with open(args['outputfile'], 'w') as fo:
        vcf.write_header(fo)
        for vnt_obj in vcf.parse_variants():
            entry = create_hgvsg(vnt_obj.CHROM,vnt_obj.POS,vnt_obj.REF,vnt_obj.ALT)
            if entry:
                if entry != 'FAIL':
                    hgvsg_entry = 'hgvsg='+entry
                    vnt_obj.add_tag_info(hgvsg_entry)
                else:
                    remove(args['outputfile'])
                    raise Exception('Unexpected variant format found. Quitting and deleting output. Variant: '+vnt_obj.CHROM+"\t"+str(vnt_obj.POS)+"\t"+vnt_obj.REF+"\t"+vnt_obj.ALT)
                    #print("ERROR: Non-standard / Unexpected Variant Found")
                    #print("ERROR: Quitting and Deleting Output")
                    #print("ERROR: Variant:",vnt_obj.CHROM,str(vnt_obj.POS),vnt_obj.REF,vnt_obj.ALT)
                    #remove(args['outputfile'])
                    #sys.exit()
            vcf.write_variant(fo, vnt_obj)

    subprocess.run(["bgzip", args['outputfile']])
    subprocess.run(["tabix",args['outputfile']+".gz"])

################################################
#   Main
################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Add hgvsg INFO field for each qualified variant')

    parser.add_argument('-i','--inputfile', help='input VCF file', required=True)
    parser.add_argument('-o','--outputfile', help='output VCF file', required=True)

    args = vars(parser.parse_args())

    main(args)
