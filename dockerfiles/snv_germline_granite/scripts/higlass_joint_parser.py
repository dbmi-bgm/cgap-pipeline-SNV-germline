#!/usr/bin/env python3

################################################
#
#   Script to parse and clean joint-called VCF
#        and carry out Fisher exact tests
#         for cohort viewing in HiGlass
#
#   Following Sentieon Joint Calling, the VCF
#   must be run through the following scripts:
#   1. bcftools-norm-multiallelics.sh
#   2. vep-annot.sh
#   3. portal_reformat_vcf.py
#
#   Step 3's output is input for this parser
#
#   A txt list of case IDs is also required
#
#   This parser calculates AN, AC, and AF for
#   probands (cases) and reports that in the
#   INFO field, alongside gnomAD v2 and v3
#   AN, AC, and AF and VEP annotations for the
#   most severe consequence in any transcript
#   for each variant
#
#   The most severe conseqence level (HIGH,
#   MODERATE, LOW, and MODIFIER) is also
#   calculated, given the VEP_order dictionary,
#   and added to the new INFO field
#
#   Fisher exact can be carried out against
#   gnomAD v2 and/or gnomAD v3 to create
#   odds ratios and p-vales for case vs controls
#   which are added to the new INFO field
#
#   All other entries are removed from the
#   INFO field. All individual IDs are
#   removed from the header columns and
#   genotype information is dropped

################################################

################################################
#   Libraries
################################################

from granite.lib import vcf_parser
import math
import sys, argparse, subprocess
from scipy.stats import fisher_exact

################################################
#   Top level variables
################################################

# AF is reported as 3 siginificant digits
# from Sentieon/GATK so we will use that for
# AF and p-values

# significant digits when calculated below 1
# e.g., AF and p
significant_digits = 3

#significant digits when calculated above 1
# e.g., OR and log10
round_digits = 4

################################################
#   Functions
################################################

def fisher_calculation(dict, gnomAD_version, proband_alt, proband_ref, gnomAD_alt, gnomAD_ref):
    '''
    This function is called within fisher_exact_gnomAD
    and runs the actual Fisher exact test once the values
    have been generated by the main function. It returns an updated
    info_dict.
    '''
    oddsratio, pvalue = fisher_exact([[proband_alt, proband_ref], [gnomAD_alt, gnomAD_ref]], alternative='greater')
    dict["fisher_gnomAD"+gnomAD_version+"_p"] = round(pvalue, significant_digits - int(math.floor(math.log10(abs(pvalue)))) - 1)
    dict["fisher_gnomAD"+gnomAD_version+"_OR"] = round(oddsratio, round_digits)
    if pvalue == 1:
        dict["fisher_gnomAD"+gnomAD_version+"_minuslog10p"] = 0
    else:
        dict["fisher_gnomAD"+gnomAD_version+"_minuslog10p"] = round(-math.log10(pvalue), round_digits)
    return dict

def fisher_exact_gnomAD(dict, gnomAD_version):
    '''
    This generated the necessary inputs for a  Fisher
    exact test (2 by 2) on the probands (cases)
    versus gnomAD version(s) supplied by the user
    and calls the fisher_calculation function to
    carry out the test. It returns an updated
    info_dict.
    '''
    gnomAD_dict = {"v2": "gnomADe2_", "v3": "gnomADg_"}
    gnomAD_base = gnomAD_dict[gnomAD_version]

    #can only run the test if there are values for gnomAD
    if dict[gnomAD_base+"AC"] != '' and dict[gnomAD_base+"AN"] != '':
        #store the proband info
        proband_alt = int(dict["AC_proband"])
        proband_ref = int(dict["AN_proband"]) - proband_alt

        #then access the gnomAD info
        #gnomAD v2 has some entries with multiple values so we need to condition on their absence first
        if "&" not in dict[gnomAD_base+"AC"] and "&" not in dict[gnomAD_base+"AN"]:
            gnomAD_alt = int(dict[gnomAD_base+"AC"])
            gnomAD_ref = int(dict[gnomAD_base+"AN"]) - gnomAD_alt

            #carry out the test, store the p value and OR value, then generate and store the -log10(p) values as well
            dict = fisher_calculation(dict, gnomAD_version, proband_alt, proband_ref, gnomAD_alt, gnomAD_ref)

        # if the entry does have an "&" we want to select the most rare to compare to
        else:
            # found that some variants with "&" in AC and AN do not have "&" in AF, which means we can't use AF index to pull AC and AN.
            # need to carry out the calculation
            gnomAD_alt_list = [int(x) for x in dict[gnomAD_base+"AC"].split("&")]
            gnomAD_AN_list = [int(x) for x in dict[gnomAD_base+"AN"].split("&")]

            # if the lists are of different lengths, we can't use a shared index, so we have to ignore variant
            if len(gnomAD_alt_list) != len(gnomAD_AN_list):
                gnomAD_alt = 0
                gnomAD_ref = 0
            else:
                # want to select the most rare from the list of possibilites
                temp_idx = 0
                #temp_AF = gnomAD_alt_list[temp_idx]/gnomAD_AN_list[temp_idx]
                for idx, value in enumerate(gnomAD_alt_list):
                    try:
                        if idx == 0:
                            temp_AF = gnomAD_alt_list[temp_idx]/gnomAD_AN_list[temp_idx]
                        else:
                            if gnomAD_alt_list[idx]/gnomAD_AN_list[idx] < temp_AF:
                                temp_AF = gnomAD_alt_list[idx]/gnomAD_AN_list[idx]
                                temp_idx = idx
                    except Exception:
                        #occurs with 0 in denominator
                        gnomAD_alt = 0
                        gnomAD_ref = 0
                        temp_idx = "FAIL"
                if temp_idx != "FAIL":
                    gnomAD_alt = gnomAD_alt_list[temp_idx]
                    gnomAD_ref = gnomAD_AN_list[temp_idx] - gnomAD_alt_list[temp_idx]

            # finally, we can calculate the same as above
            dict = fisher_calculation(dict, gnomAD_version, proband_alt, proband_ref, gnomAD_alt, gnomAD_ref)

    #if there aren't values for gnomAD, store empty strings, which will be replaced by NA downstream
    else:
        dict["fisher_gnomAD"+gnomAD_version+"_p"] = ''
        dict["fisher_gnomAD"+gnomAD_version+"_OR"] = ''
        dict["fisher_gnomAD"+gnomAD_version+"_minuslog10p"] = ''
    return dict

def main(args):
    vcf = vcf_parser.Vcf(args['inputfile'])
    VEPtag = 'CSQ'
    GENEtag = 'GENES'

    #we want gnomAD v2 and v3 allele frequencies etc.
    gnomADg_AC_idx = vcf.header.get_tag_field_idx(VEPtag, 'gnomADg_AC')
    gnomADg_AN_idx = vcf.header.get_tag_field_idx(VEPtag, 'gnomADg_AN')
    gnomADg_AF_idx = vcf.header.get_tag_field_idx(VEPtag, 'gnomADg_AF')
    gnomADe2_AC_idx = vcf.header.get_tag_field_idx(VEPtag, 'gnomADe2_AC')
    gnomADe2_AN_idx = vcf.header.get_tag_field_idx(VEPtag, 'gnomADe2_AN')
    gnomADe2_AF_idx = vcf.header.get_tag_field_idx(VEPtag, 'gnomADe2_AF')
    worst_trscrpt_idx = vcf.header.get_tag_field_idx(VEPtag, 'most_severe')

    #we want worst consequence to generate the level of the worst consequence (below)
    most_severe_consequence_idx = vcf.header.get_tag_field_idx(GENEtag, 'most_severe_consequence')

    #VEP_order taken from portal_reformat_vcf.py. If this is updated, you will also need to update the if/elif/else for level calculation
    VEP_order = {
                    # HIGH
                    'transcript_ablation': 1,
                    'splice_acceptor_variant': 2,
                    'splice_donor_variant': 3,
                    'stop_gained': 4,
                    'frameshift_variant': 5,
                    'stop_lost': 6,
                    'start_lost': 7,
                    'transcript_amplification': 8,
                    # MODERATE
                    'inframe_insertion': 9,
                    'inframe_deletion': 10,
                    'missense_variant': 11,
                    'protein_altering_variant': 12,
                    # LOW
                    'splice_region_variant': 13,
                    'incomplete_terminal_codon_variant': 14,
                    'start_retained_variant': 15,
                    'stop_retained_variant': 16,
                    'synonymous_variant': 17,
                    # MODIFIER
                    'coding_sequence_variant': 18,
                    'mature_miRNA_variant': 19,
                    '5_prime_UTR_variant': 20,
                    '3_prime_UTR_variant': 21,
                    'intron_variant': 22,
                    'MODIFIER': 23
                }

    #the proband list contains the cases (or probands) that we want to use for case AN, AC, and AF calculations. currently, gnomAD is going to be used for controls
    #note that these values will differ from the AN, AC, and AF calculated by Sentieon which represent both the probands and family members (all input files)
    proband_list = []
    with open(args['probandlist'], 'r') as pl:
        for p in pl:
            proband_list.append(p.strip("\n"))

    #below, we do the calculations and parsing
    with open(args['outputfile'], 'w') as fo:
        #remove all sample IDs from the header columns - not needed at this time and also makes the data more anonymous
        vcf.header.columns='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        vcf.write_header(fo)
        for vnt_obj in vcf.parse_variants():
            #create an empty dictionary to store the new INFO field keys and values
            info_dict = {}
            #also create an ordered list of the keys we will be including so that we can loop this to write them in the same order for each variant
            info_list = ["AC_proband", "AN_proband", "AF_proband", "gnomADg_AC", "gnomADg_AN", "gnomADg_AF", "gnomADe2_AC", "gnomADe2_AN", "gnomADe2_AF",  "most_severe_consequence", "level_most_severe_consequence"]
            for version in args['gnomAD']:
                info_list.append("fisher_gnomAD"+version+"_OR")
                #info_list.append("fisher_gnomAD"+version+"_p")
                info_list.append("fisher_gnomAD"+version+"_minuslog10p")

            #get the index for genotype (GT) and pull genotypes for all probands
            GT_idx = vnt_obj.FORMAT.split(":").index("GT")
            proband_genotypes = {}
            for proband in proband_list:
                proband_genotypes[proband] = vnt_obj.GENOTYPES[proband].split(":")[GT_idx]

            # if we want the raw (cohort-wide) AC, AN, and AF values, bring this back and add them to info_list
            # info_dict["AC"] = vnt_obj.get_tag_value("AC")
            # info_dict["AN"] = vnt_obj.get_tag_value("AN")
            # info_dict["AF"] = vnt_obj.get_tag_value("AF")

            #for the probands, calculate AC, AN and AF
            AC_proband = 0
            AN_proband = 0
            for proband in proband_genotypes:
                # add 2 to total count if there is a genotype called
                if proband_genotypes[proband] != "./.":
                    AN_proband += 2
                    # then count the alternative alleles
                    if proband_genotypes[proband] == "0/0":
                        pass
                    elif proband_genotypes[proband] == "1/0" or proband_genotypes[proband] == "0/1":
                        AC_proband += 1
                    elif proband_genotypes[proband] == "1/1":
                        AC_proband += 2
                    else:
                        raise Exception("Unexpected genotype found for variant "+vnt_obj.CHROM+"\t"+str(vnt_obj.POS)+"\t"+proband+"\t"+proband_genotypes[proband]+"\n"+"Did you run bcftools norm multiallelics?")

            #with counting done, we need to calculate AF and store all 3 values in the info_dict
            if AC_proband == 0:
                info_dict["AF_proband"] = 0
            else:
                AF_proband = AC_proband / AN_proband
                #AF is reported as 3 siginificant digits from Sentieon/GATK, so this bit of math replicates that
                info_dict["AF_proband"] = round(AF_proband, significant_digits - int(math.floor(math.log10(abs(AF_proband)))) - 1)


            info_dict["AC_proband"] = AC_proband
            info_dict["AN_proband"] = AN_proband

            #now we need to get the worst consequence and generate the level
            try:
                val_get = vnt_obj.get_tag_value(VEPtag)
            except Exception:
                print("VEP not run?")
            trscrpt_list = val_get.split(',')
            for t in trscrpt_list:
                if t.split('|')[worst_trscrpt_idx] == "1":
                    # if it is the worst transcript (1 per variant)
                    # we want to get all relevant gnomAD annotations
                    info_dict["gnomADg_AC"] = t.split('|')[gnomADg_AC_idx]
                    info_dict["gnomADg_AN"] = t.split('|')[gnomADg_AN_idx]
                    info_dict["gnomADg_AF"] = t.split('|')[gnomADg_AF_idx]
                    info_dict["gnomADe2_AC"] = t.split('|')[gnomADe2_AC_idx]
                    info_dict["gnomADe2_AF"] = t.split('|')[gnomADe2_AF_idx]
                    info_dict["gnomADe2_AN"] = t.split('|')[gnomADe2_AN_idx]
                    break
            info_dict["most_severe_consequence"] = vnt_obj.get_tag_value(GENEtag).split('|')[most_severe_consequence_idx]
            try: i = VEP_order[info_dict["most_severe_consequence"]]
            except Exception:
                i = 23
            if i >= 1 and i <= 8:
                info_dict["level_most_severe_consequence"] = "HIGH"
            elif i > 8 and i <= 12:
                info_dict["level_most_severe_consequence"] = "MODERATE"
            elif i > 12 and i <= 17:
                info_dict["level_most_severe_consequence"] = "LOW"
            elif i > 17 and i <= 23:
                info_dict["level_most_severe_consequence"] = "MODIFIER"

            #run the Fisher exact test
            for version in args['gnomAD']:
                info_dict = fisher_exact_gnomAD(info_dict, version)

            #create new INFO field and replace existing one
            INFO = ""
            for field in info_list:
                INFO+=field+"="
                if info_dict[field] == '':
                    INFO+="NA"
                else:
                    INFO+=str(info_dict[field])
                INFO+=";"
            vnt_obj.INFO = INFO.strip(";")
            vnt_obj.IDs_genotypes = None
            vcf.write_variant(fo, vnt_obj)


    subprocess.run(["bgzip", args['outputfile']])
    subprocess.run(["tabix",args['outputfile']+".gz"])

################################################
#   Main
################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Add hgvsg INFO field for each qualified variant')

    parser.add_argument('-i','--inputfile', help='input VCF file', required=True)
    parser.add_argument('-p','--probandlist', help='txt file with all cases (proband) IDs on individual lines', required=True)
    parser.add_argument('-o','--outputfile', help='output VCF file', required=True)
    parser.add_argument('-g','--gnomAD', nargs='*', help='gnomAD version(s) to use as control (v2 and/or v3)', choices=['v2', 'v3'], required=True)

    args = vars(parser.parse_args())

    main(args)
