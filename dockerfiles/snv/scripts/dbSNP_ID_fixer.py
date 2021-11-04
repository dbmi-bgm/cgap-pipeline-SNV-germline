#!/usr/bin/env python3

################################################
#
#  Script to add dbSNP rsIDs from dbSNP vcf
#   into the ID column of sample vcf files
#  Note: Must be run after multiallelic split
#
################################################

################################################
#   Libraries
################################################
import sys,argparse
import tabix

################################################
#   Classes
################################################

# for the dbSNP class, a list is created for ALT to allow comma-separated alleles (e.g., A -> AT,ATT,ATTT,etc.)
class dbSNP_entry:
    def __init__(self, CHROM, POS, ID, REF, ALT):
        self.CHROM = CHROM
        self.POS = POS
        self.ID = ID
        self.REF = REF
        self.ALT = ALT.split(",") #comment here and in documentation as well

# for the sample_entry class, we have already carried out a step to place multiallelic variants onto separate lines
# note that ALT here does not create a list
# this code should not be run prior to running bcftools norm
class sample_entry:
    def __init__(self, CHROM, POS, ID, REF, ALT):
        self.CHROM = CHROM
        self.POS = POS
        self.ID = ID
        self.REF = REF
        self.ALT = ALT

################################################
#   Functions
################################################

def main(args):
    # open dbSNP_vcf and sample_vcf with tabix
    dbSNP_vcf = args['dbSNPvcf']
    dbSNP_tabix = tabix.open(dbSNP_vcf)

    sample_vcf = args['inputvcf']
    sample_tabix = tabix.open(sample_vcf)

    # regionfile is used for both the region query from the sample_vcf and also the output file name
    in_region = args['regionfile']
    out_file_name = args['regionfile']

    with open(out_file_name, 'w') as fo:
        # create tabix iterator for full region in sample_vcf
        try:
            sample_region = sample_tabix.querys(in_region)
        except tabix.TabixError:
            return
        for sample_record in sample_region:
            # create variant object and corresponding region to query from dbSNP_vcf
            vnt_obj = sample_entry(sample_record[0],sample_record[1],sample_record[2],sample_record[3],sample_record[4])
            dbSNP_region = vnt_obj.CHROM+":"+str(vnt_obj.POS)+"-"+str(vnt_obj.POS)
            overlapping_rsIDs = dbSNP_tabix.querys(dbSNP_region)

            # create empty list for replacement rsIDs
            rs_IDs = []

            # iterate through overlapping dbSNP entries to idenitfy matches through position, and reference / alternative alleles
            for db_record in overlapping_rsIDs:
                db_obj = dbSNP_entry(db_record[0],db_record[1],db_record[2],db_record[3],db_record[4])
                if vnt_obj.REF == db_obj.REF and str(vnt_obj.POS) == str(db_obj.POS):
                    for alt_allele in db_obj.ALT:
                        if vnt_obj.ALT == alt_allele:
                            rs_IDs.append(db_obj.ID) # save those that correspond to the variant

            # change rs_IDs list into a string depending on results of search
            if len(rs_IDs) > 0:
                replacement_ID = ';'.join([str(elem) for elem in rs_IDs])
            else:
                replacement_ID = "." #important to replace incorrect rsIDs with "." due to carryover of incorrect rsIDs in multiallelic splitting

            # replace the existing ID column with the replacement_ID column, depending on results of our search and existing IDs
            if vnt_obj.ID == ".":
                sample_record[2] = replacement_ID
            else:
                current_IDs = vnt_obj.ID.split(";")
                if len(current_IDs) == 1:
                    if 'rs' in current_IDs[0]:
                        sample_record[2] = replacement_ID
                    else:
                        if replacement_ID != ".":
                            sample_record[2] = sample_record[2]+";"+replacement_ID
                        else:
                            sample_record[2] = sample_record[2]
                if len(current_IDs) > 1:
                    keep_ID_list = []
                    for current_ID in current_IDs:
                        if 'rs' not in current_ID:
                            keep_ID_list.append(current_ID)
                    if len(keep_ID_list) > 0:
                        keep_IDs = ';'.join([str(elem) for elem in keep_ID_list])
                        if replacement_ID != ".":
                            sample_record[2] = keep_IDs+";"+replacement_ID
                        else:
                            sample_record[2] = keep_IDs
                    else:
                        sample_record[2] = replacement_ID
            fo.write('\t'.join(sample_record)+'\n')

################################################
#   MAIN
################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-db', '--dbSNPvcf', help='input dbSNP vcf', required=True)
    parser.add_argument('-i', '--inputvcf',  help='input sample vcf', required=True)
    parser.add_argument('-r', '--regionfile', help='regionfile', required=True)

    args = vars(parser.parse_args())

    main(args)
