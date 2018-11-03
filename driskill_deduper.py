#!/usr/bin/env python
import argparse

def get_arguments():
    parser = argparse.ArgumentParser(description='k-mer size program')
    parser.add_argument("-f", "--file", help="input string", required=True, type=str)
    parser.add_argument("-o", "--out_file", help="input string", required=True, type=str)
    parser.add_argument("-u", "--UMI_file", help="input string", required=True, type=str)



    return parser.parse_args()
    
args = get_arguments()
file = args.file
out_file = args.out_file
UMI_file = args.UMI_file

###################################################################################################################
help="\n -f is in file \n -o is out file name \n -u is UMI file for dictionary \n -p is for paired end file"
paired="This program does not take paired end reads"

###################################################################################################################


import re
import sys
import os

UMI_dict={}
for line in open(UMI_file).readlines():
    umi=line.strip()
    UMI_dict[umi]=None


ctr=1

unique_dict={}
current_rname=""
previous_rname=""

ctr_umi_notin_dict=0
forward_reads=0
reverse_reads=0
unique_reads=0
###################################################################################################################

def pair_end_reads(flag):
    if((int(flag) & 1) == 1):
        sys.exit("Error, file contains pair-end reads")
    else:
        return True 

def forward_adjust(pos):
    softclip=re.search(r"^(\d+)S", cigar)
    if softclip is None:
        softclip=softclip
    else:
        s1=int(softclip.group(1))
        pos=pos-s1

    string=current_rname+UMI+str(pos)+"Plus"
    return string

def reverse_adjust(pos):
    pos=pos+total_S+total_M+total_D+total_N-1

    string=current_rname+UMI+str(pos)+"Minus"
    return string

###################################################################################################################

with open (file, "r") as fh, \
open(out_file + '_deduped.sam', 'w') as PCR_removed, \
open(out_file + '_undetermine.sam', 'w') as undetermined, \
open(out_file + '_count_data', 'w') as count_data:
    
###################################################################################################################
    
    for line in fh:
        if line.startswith("@"):
            PCR_removed.write(line)
            continue
        else:
            
            parts = line.split("\t")#this splits all the colums in the line into lists(parts)
            qname=parts[0]
            flag=parts[1]
            pair_end_reads(flag)
            
            current_rname=parts[2]
            
            if ctr==1:
                previous_rname=current_rname
            ctr+=1
            
            if current_rname != previous_rname:
                previous_rname=current_rname
                del unique_dict
                unique_dict={}
                
            pos=int(parts[3])
            cigar=parts[5]
            rnext=parts[6]            
            
            
            S = re.findall(r'(\d+)S$', cigar)#finds all end S's
            S=list(map(int,S))
            total_S = sum([int(i) for i in S])
            
            M = re.findall(r"(\d+)M", cigar)
            M = list(map(int,M))
            total_M = sum([int(i) for i in M])
            
            D = re.findall(r"(\d+)D", cigar)
            D = list(map(int,D))
            total_D = sum([int(i) for i in D])
            
            N = re.findall(r"(\d+)N", cigar)
            N = list(map(int,N))
            total_N = sum([int(i) for i in N])
            
            qname_parts=qname.split(":")
            UMI=qname_parts[7]
            
###################################################################################################################
                           

            if UMI not in UMI_dict.keys():
                ctr_umi_notin_dict+=1
                undetermined.write(line)
            else:
                if UMI in UMI_dict.keys():

###################################################################################################################

                    if((int(flag) & 16) != 16):
                        forward_reads+=1
                        string=forward_adjust(pos)
                        #print(string)

###################################################################################################################

                    
                    if((int(flag) & 16) == 16):
                        reverse_reads+=1
                        string=reverse_adjust(pos)
                        #print(string)
                        
###################################################################################################################

                    
                    if string not in unique_dict:
                        unique_reads+=1
                        unique_dict[string]=None
                        #print(unique_dict)
                        PCR_removed.write(line)
                            

###################################################################################################################

    count_data.write("{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n".format("UMI Not Found in Dictionary", ctr_umi_notin_dict, \
                                                               "Number of Forward Reads", forward_reads, \
                                                               "Number of Reverse Reads",  reverse_reads, \
                                                               "Number of Unique Reads Found", unique_reads)) 




