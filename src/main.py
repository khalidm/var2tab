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
# import pysam

# def getcadd(cadd_tbx, current_chr, current_pos, current_ref, current_alt):
#     current_chr = current_chr.translate(None, 'chr')
#     data = cadd_tbx.fetch(current_chr, current_pos-1, current_pos)
#     cadd_phred, cadd_priPhCons, cadd_GerpRS = '','',''
#     cadd_polysift, cadd_test1, cadd_test2 = '','',''
#
#     if data is not None:
#         for row in data:
#             row_info = row.split("\t")
#             cadd_ref = row_info[2]
#             cadd_alt = row_info[4]
#             if(cadd_ref == current_ref and cadd_alt == current_alt):
#                 cadd_phred = row_info[115]
#                 cadd_priPhCons = row_info[18]
#                 cadd_GerpRS = row_info[26]
#                 if "damaging" in row_info[110] or "deleterious" in row_info[112]:
#                     cadd_polysift = "del"
#                 break
#     else:
#         cadd_phred = 'NA'
#
#     return cadd_phred, cadd_priPhCons, cadd_GerpRS, \
#             cadd_polysift
#
# # return allele frequency given the allele count and assuming allele number = (total allele number/2)
# def getAF(ac, an):
#     if(float(an)>0.0):
#         af_temp = float(ac) / an
#         #newlist = round(af_temp, 8)
#         newlist = af_temp
#     else:
#         newlist = 0.0
#     return str(newlist)
#
# # return index of the current alt allele from exac multiallelic data
# def getexacallele(exac_tbx, current_chr, current_pos, current_ref, current_alt):
#     current_chr = current_chr.translate(None, 'chr')
#     data = exac_tbx.fetch(current_chr, current_pos-1, current_pos)
#     index = -2
#     row = 0
#     found = False
#     exac_filter_return = ""
#     if data:
#         for exac_row in data:
#             exac_pos = int(exac_row.split("\t")[1])
#
#             exac_ref_temp = exac_row.split("\t")[3]
#             exac_alt_temp = exac_row.split("\t")[4]
#             exac_filter = exac_row.split("\t")[6]
#
#             exac_alt_row = exac_alt_temp.split(",")
#             exac_ref_row = exac_ref_temp.split(",")
#
#             #if(current_pos == exac_pos and current_ref in exac_ref_row and \
#             if(current_pos == exac_pos and current_alt in exac_alt_row ):
#                 index = exac_alt_row.index(current_alt)
#                 exac_filter_return = exac_filter
#                 row += 1
#                 break
#     else:
#         index = -2
#
#     #print "Row = " + str(row) + " " + str(found)
#     return index, exac_filter_return

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
    outputfile.write("chr\tpos\tid\tref\talt\tannotation\tgene_name\tlof" \
            "\texon\taa_pos\tSIFT\tPOLYPHEN\tCONDEL\tAF\tGMAF\t1kgEMAF\tESPEMAF\t" \
            "ExAC_AF\tExAC_EAS\tExAC_NFE\tExAC_FIN\tExAC_SAS\tExAC_AFR\tExAC_AMR\tExAC_OTH\t" \
            "CADD\tmaxCADD\tpriPhCons\tGerpRS\tFATHMM\t" \
            "Mapability\tPromoter\tEnhancer\tRepeat\tPfam\t" \
            "CPG\tClinVar\tGWAS\tMNP_FLAG\tExAC_FLAG\n")

    vcf_row = {}

    vcf_reader = vcf.Reader(open(args.vcf, 'r'))
    for record in vcf_reader:
        current_chr = record.CHROM
        current_id = record.ID
        current_pos = record.POS
        current_ref = record.REF
        current_alt = ','.join(str(v) for v in record.ALT)
        current_AF = record.INFO['AF']
        current_het_nfe = ''
        current_hom_nfe = ''

        # CHECK INDEL AND MNP
        indel = True if ((len(current_ref) > 1 or len(current_alt) > 1) and \
                ("," not in current_ref and "," not in current_alt)) else False
        mnp = True if len(record.ALT) > 1 else False
        mnpflag = "%s" % mnp

        # VEP
        current_sift, current_polyphen, current_consequence, current_LOF = '','','',''
        current_sift_score, current_polyphen_score = 0.9999, 0.0001
        current_gmaf, current_eur_maf, current_ea_maf = '','',''
        current_feature, current_feature_type = '',''
        if "CSQ" in record.INFO:
            csq = record.INFO['CSQ'][0].split('|')
            current_feature, current_feature_type = csq[2], csq[3]
            current_consequence = csq[4]

            #print csq[24] + "-" + csq[25]
            current_sift = csq[23]
            current_polyphen = csq[24]
            if ( len(current_sift) > 0):
                current_sift_score = re.findall(r'[0-9.]+', current_sift)[0]
            if ( len(current_polyphen) > 0):
                current_polyphen_score = re.findall(r'[0-9.]+', current_polyphen)[0]


            current_gmaf, current_eur_maf, current_ea_maf =  csq[31], csq[35], csq[37]
            #current_LOF = csq[48]
        else:
            current_feature, current_feature_type, current_consequence = '','',''
            current_sift, current_polyphen, current_eur_maf = '','',''
            current_ea_maf, current_LOF, current_gmaf = '','',''

        # # SnpEff
        # ann = record.INFO['ANN'][0].split('|')
        # annotation = ann[1]
        # #   GENE INFORMATION
        # current_gene, current_exon, current_aa_pos = ann[3], ann[8], ann[10]

        #CADD SNP
        # cadd_phred_temp = ''
        # cadd_phred = ''
        # indel_str= ''
        # mnp_cadds = []
        # cadd_scores = []
        # fathmm_score = 0.0
        # for alt in record.ALT:
        #     if(len(current_ref) == 1 and len(alt) == 1):
        #         (cadd_phred_temp, cadd_priPhCons, cadd_GerpRS, cadd_polysift) = \
        #                 getcadd(cadd_tbx, current_chr, current_pos, current_ref, alt)
        #         mnp_cadds.append(str(alt) + ":" + cadd_phred_temp)
        #         cadd_scores.append(cadd_phred_temp)
        #         # GET FATHMM SCORE
        #         fathmm_score = getfathmm(fathmm_tbx, current_chr, current_pos, current_ref, alt)
        #     else: # IF VAR IS AN INDEL
        #         (cadd_phred_temp, cadd_priPhCons, cadd_GerpRS, cadd_polysift) = \
        #                 getcadd(cadd_indel_tbx, current_chr, current_pos, current_ref, alt)
        #         mnp_cadds.append(str(alt) + ":" + cadd_phred_temp)
        #         cadd_scores.append(cadd_phred_temp)
        # cadd_phred = ",".join(mnp_cadds)
        # indel_str = "."

        # INSERT OTHER TABIX BASED ANNOTATORS BELOW
        current_AF = record.INFO['AF']

        # RESCORE SCORES FOR PROTEIN TRUNCATING MUTATIONS
        # (current_condel, current_sift, current_polyphen, fathmm_score) = adjust_scores(current_condel, current_sift, \
        #         current_polyphen, fathmm_score, annotation)

        out_str = [ current_chr, str(current_pos), str(current_id), current_ref, current_alt,
                annotation, current_gene, current_LOF, current_exon,
                current_aa_pos, str(current_sift_score), str(current_polyphen_score), str(current_condel), str(current_AF),
                current_gmaf, current_eur_maf, current_ea_maf,
                str(current_exac_af), str(current_exac_eas), str(current_exac_nfe), str(current_exac_fin),
                str(current_exac_sas), str(current_exac_afr), str(current_exac_amr), str(current_exac_oth),
                cadd_phred, str(max(cadd_scores)), cadd_priPhCons, cadd_GerpRS,
                str(fathmm_score), str(current_mapability), current_promoter, current_enhancer,
                current_rmsk, current_pfam, current_cpg, current_clinvar, current_gwas,
                mnpflag, exac_flag]
        out_str = [x or '.' for x in out_str]

        # filter away!


        outputfile.write("\t".join(out_str))
        outputfile.write("\n")
        outputfile.close()

if __name__ == "__main__":
    main(sys.argv)
