'''
vcf2tab
A tool to convert vcf data to tabular format.
Author: Khalid Mahmood
Contact: khalid.mahmood@unimelb.edu.au
Copyright: 2016
'''

#!/usr/bin/python

# from utils import findlist
# from annotations import getTabixVal,getTabixValCondel,getTabixBool
# from annotations import getfathmm,adjust_scores

import sys
import os
import argparse
import getopt
import vcf
import re
import array

# MAIN

def main(argv):

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--vcf", type=str, dest="vcf", help="Input variant file (vcf)", required=True)
    parser.add_argument("-o", "--output", type=str, dest="out", help="Output file (tabular)", required=True)
    parser.add_argument("-v", "--verbosity", action="count", default=0)

    args = parser.parse_args()
    outputfile = open(args.out, "w")

    if args.verbosity >= 2:
        print "{} to the power {} equals {}".format(args.v, args.o, answer)
    elif args.verbosity >= 1:
        print "{}^{} == {}".format(args.x, args.y, answer)


    # Output header
    outputfile.write("chr\tpos\tid\tref\talt\tgene\tfeature\tfeature_type\tconsequence\n")

    vcf_row = {}

    vcf_reader = vcf.Reader(open(args.vcf, 'r'))
    for record in vcf_reader:
        current_chr = record.CHROM
        current_id = record.ID
        current_pos = record.POS
        current_ref = record.REF
        current_alt = ','.join(str(v) for v in record.ALT)

        # CHECK INDEL AND MNP
        indel = True if ((len(current_ref) > 1 or len(current_alt) > 1) and \
                ("," not in current_ref and "," not in current_alt)) else False
        mnp = True if len(record.ALT) > 1 else False
        mnpflag = "%s" % mnp

        # VEP fields
        current_gene, current_feature = '',''
        current_feature_type, current_consequence = '',''
        if "CSQ" in record.INFO:
            csq = record.INFO['CSQ']

            # BELOW: THERE ARE A COUPLE OF OPTIONS TO PROCEED
            # For going through annotations for all transcript
            for current_csq_element in csq:
                current_csq = current_csq_element.split('|')
                current_gene = current_csq[0]
                current_feature = current_csq[1]
                current_feature_type = current_csq[2]
                current_consequence = current_csq[3]
                out_str = [ current_chr, str(current_pos), str(current_id), current_ref, current_alt,
                            current_gene, current_feature, current_feature_type, current_consequence]
                out_str = [x or 'None' for x in out_str]
                outputfile.write("\t".join(out_str))
                outputfile.write("\n")

        else:
            current_gene, current_feature = '',''
            current_feature_type, current_consequence = '',''
            out_str = [ current_chr, str(current_pos), str(current_id), current_ref, current_alt,
                        current_gene, current_feature, current_feature_type, current_consequence]
            out_str = [x or 'None' for x in out_str]
            outputfile.write("\t".join(out_str))
            outputfile.write("\n")

    outputfile.close()

if __name__ == "__main__":
    main(sys.argv)
