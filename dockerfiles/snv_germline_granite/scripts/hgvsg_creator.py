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

##################################################
#Dictionaries for converting chr to NCBI accession
##################################################

# GRCh38 - accession numbers for the 23 chromosomes haven't changed between GRCh38 (https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/) and the GRCh38.p13 (https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.39/)
NCBI_CONVERT =  {
                "chr1" : "NC_000001.11",
                "chr2" : "NC_000002.12",
                "chr3" : "NC_000003.12",
                "chr4" : "NC_000004.12",
                "chr5" : "NC_000005.10",
                "chr6" : "NC_000006.12",
                "chr7" : "NC_000007.14",
                "chr8" : "NC_000008.11",
                "chr9" : "NC_000009.12",
                "chr10" : "NC_000010.11",
                "chr11" : "NC_000011.10",
                "chr12" : "NC_000012.12",
                "chr13" : "NC_000013.11",
                "chr14" : "NC_000014.9",
                "chr15" : "NC_000015.10",
                "chr16" : "NC_000016.10",
                "chr17" : "NC_000017.11",
                "chr18" : "NC_000018.10",
                "chr19" : "NC_000019.10",
                "chr20" : "NC_000020.11",
                "chr21" : "NC_000021.9",
                "chr22" : "NC_000022.11",
                "chrX" : "NC_000023.11",
                "chrY" : "NC_000024.10",
                "chrM" : "NC_012920.1"
                }
#GRCh37 (hg19) - https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/
NCBI_CONVERT_hg19 = {
                    "chr1" : "NC_000001.10",
                    "chr2" : "NC_000002.11",
                    "chr3" : "NC_000003.11",
                    "chr4" : "NC_000004.11",
                    "chr5" : "NC_000005.9",
                    "chr6" : "NC_000006.11",
                    "chr7" : "NC_000007.13",
                    "chr8" : "NC_000008.10",
                    "chr9" : "NC_000009.11",
                    "chr10" : "NC_000010.10",
                    "chr11" : "NC_000011.9",
                    "chr12" : "NC_000012.11",
                    "chr13" : "NC_000013.10",
                    "chr14" : "NC_000014.8",
                    "chr15" : "NC_000015.9",
                    "chr16" : "NC_000016.9",
                    "chr17" : "NC_000017.10",
                    "chr18" : "NC_000018.9",
                    "chr19" : "NC_000019.9",
                    "chr20" : "NC_000020.10",
                    "chr21" : "NC_000021.8",
                    "chr22" : "NC_000022.10",
                    "chrX" : "NC_000023.10",
                    "chrY" : "NC_000024.9",
                    "chrM" : "NC_001807.4" # https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/ hg19 has this mito in main assembly (used for LO)
                    }

################################################
#   Functions
################################################

# designed given best practices at http://varnomen.hgvs.org/recommendations/DNA/
def create_hgvsg(CHROM, POS, REF, ALT, DICT):
    try:
        # first check that there are no non-standard variants:
        if len(REF) > 1 and len(ALT) > 1 and REF[0] == ALT[0] and REF[1] == ALT[1]:
            HGVSG = 'FAIL'
        else:
            # start by identifying genomic versus mitochondrial
            # note that mito requires some code below for variants that span the end to beginning of the circular DNA (e.g., a deletion of 16568-3)
            HGVSG = DICT[CHROM]
            if HGVSG != "NC_012920.1" and HGVSG != "NC_001807.4":
                HGVSG += ':g.'
            else:
                HGVSG += ':m.'
                # hg38 and hg19 mito have slightly different lengths
                if HGVSG == "NC_012920.1:m.":
                    mito_length = 16569
                if HGVSG == "NC_001807.4:m.":
                    mito_length = 16571

            # now on to the different types of variants
            # deal with dashes as a first case - they only appear in ALT and are all deletions
            # deletions - http://varnomen.hgvs.org/recommendations/DNA/variant/deletion/
            if ALT == "-":
                if len(REF) == 1:
                    HGVSG += str(POS) + 'del'
                else:
                    end_position = POS + len(REF) - 1
                    if ':m.' in HGVSG and end_position > mito_length:
                        circular_position = end_position - mito_length
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
                        if ':m.' in HGVSG and end_position > mito_length:
                            circular_position = end_position - mito_length
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
                            if ':m.' in HGVSG and start_position > mito_length:
                                circular_position = start_position - mito_length
                                HGVSG += str(circular_position) + 'del'
                            else:
                                HGVSG += str(start_position) + 'del'
                        else:
                            end_position = POS + len(REF) - 1
                            if ':m.' in HGVSG and end_position > mito_length:
                                if start_position > mito_length:
                                    start_position = start_position - mito_length
                                    circular_position = start_position + len(REF) - 2
                                    HGVSG += str(start_position) + '_' + str(circular_position) + 'del'
                                else:
                                    circular_position = end_position - mito_length
                                    HGVSG += str(start_position) + '_' + str(circular_position) + 'del'
                            else:
                                HGVSG += str(start_position) + '_' + str(end_position) + 'del'
                    # insertions have REF e.g., A and ALT e.g., ATTTTTTT and they want the inserted bases and the positions that they are between
                    # http://varnomen.hgvs.org/recommendations/DNA/variant/insertion/
                    elif len(REF) < len(ALT):
                        if ':m.' in HGVSG and POS == mito_length:
                            HGVSG += str(POS) + '_' + str(1) + 'ins' + ALT[1:]
                        else:
                            HGVSG += str(POS) + '_' + str(POS + 1) + 'ins' + ALT[1:]
    except:
        HGVSG = ''
    return HGVSG

def main(args):
    vcf = vcf_parser.Vcf(args['inputfile'])
    hgvsg_definition = '##INFO=<ID=hgvsg,Number=1,Type=String,Description="hgvsg created from variant following best practices - http://varnomen.hgvs.org/recommendations/DNA/">'
    hgvsg_19_definition = '##INFO=<ID=hgvsg_hg19,Number=1,Type=String,Description="hgvsg for liftover coordinates in hg19 created from variant following best practices - http://varnomen.hgvs.org/recommendations/DNA/">'
    vcf.header.add_tag_definition(hgvsg_19_definition)
    vcf.header.add_tag_definition(hgvsg_definition)

    with open(args['outputfile'], 'w') as fo:
        vcf.write_header(fo)
        for vnt_obj in vcf.parse_variants():
            entry = create_hgvsg(vnt_obj.CHROM,vnt_obj.POS,vnt_obj.REF,vnt_obj.ALT,NCBI_CONVERT)
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
            try:
                entry_hg19 = create_hgvsg(vnt_obj.CHROM,int(vnt_obj.INFO.split("hg19_pos=")[1].split(";")[0]),vnt_obj.REF,vnt_obj.ALT,NCBI_CONVERT_hg19)
            except:
                entry_hg19 = ''
            if entry_hg19:
                if entry_hg19 != 'FAIL':
                    hgvsg_hg19_entry = 'hgvsg_hg19='+entry_hg19
                    vnt_obj.add_tag_info(hgvsg_hg19_entry)
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
