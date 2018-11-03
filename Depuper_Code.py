
# coding: utf-8

# # Depuper Code
# 

# example of see respository for real test_umi_dict.txt 
# CTGTTCAC
# CTGTTCAG
# CTGTTCAT
# CTGTTCAA
# CTGTTCTA
# CTGTTCTC
# CTGTTCTG
# CTGTTCTT

# example of see respository for real test_sam.txt
# changed second entry back and forth from flag 1 to 0 to test for paired end on real file
# 
# @HD	VN:1.0	SO:unsorted
# @PG	ID:GSNAP	PN:gsnap	VN:2017-10-12	CL:gsnap.avx2 --gunzip -t 26 -A sam -m 5 -d mm10_chroms -D /projects/bgmp/coonrod/mmu/INTEL -s /projects/bgmp/coonrod/mmu/INTEL/mm10_chroms/mm10_chroms.maps/Mus_musculus.GRCm38.89.splicesites.iit --split-output=/projects/bgmp/coonrod/deduper/gsnap//Datset1 /projects/bgmp/coonrod/deduper//Dataset1.fastq_dups.gz
# @SQ	SN:1	LN:195471971
# @SQ	SN:2	LN:182113224
# @SQ	SN:3	LN:160039680
# @SQ	SN:4	LN:156508116
# @SQ	SN:5	LN:151834684
# @SQ	SN:6	LN:149736546
# @SQ	SN:7	LN:145441459
# @SQ	SN:8	LN:129401213
# @SQ	SN:9	LN:124595110
# @SQ	SN:10	LN:130694993
# @SQ	SN:11	LN:122082543
# @SQ	SN:12	LN:120129022
# @SQ	SN:13	LN:120421639
# @SQ	SN:14	LN:124902244
# @SQ	SN:15	LN:104043685
# @SQ	SN:16	LN:98207768
# @SQ	SN:17	LN:94987271
# @SQ	SN:18	LN:90702639
# @SQ	SN:19	LN:61431566
# @SQ	SN:X	LN:171031299
# @SQ	SN:Y	LN:91744698
# @SQ	SN:MT	LN:16299
# 1)NS500451:single-end:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC	0	1	100	36	5S10M5I5M10N5D71M40S	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
# 2)NS500451:paired-end:HWKTMBGXX:1:11101:24260:1121:CTGTTCAG	0	1	100	36	71M	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
# 3)NS500451:umi-in-dictionary:HWKTMBGXX:1:11101:24260:1121:CTGTTCTT	0	3	100	36	71M	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
# 4)NS500451:umi-not-in-dictionary:HWKTMBGXX:1:11101:24260:1121:TTTTTTTT	16	4	100	36	71M	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
# 5)NS500451:flag-is-forward:HWKTMBGXX:1:11101:24260:1121:CTGTTCAT	16	5	100	36	71M	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
# 6)NS500451:flag-is-reverse:HWKTMBGXX:1:11101:24260:1121:CTGTTCAA	0	6	100	36	71M	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
# 7)NS500451:forward&has-s:HWKTMBGXX:1:11101:24260:1121:CTGTTCTA	16	6	100	36	40S71M	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
# 8)NS500451:forward&noaction:HWKTMBGXX:1:11101:24260:1121:CTGTTCTC	0	8	100	36	71M40S	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
# 9)NS500451:forward&has-s:HWKTMBGXX:1:11101:24260:1121:CTGTTCTC	16	9	100	36	5S71M40S	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
# 10)NS500451:forward&has-s:HWKTMBGXX:1:11101:24260:1121:CTGTTCTC	0	10	100	36	5S71M40S	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
# 11)NS500451:forward&has-s,dupe:HWKTMBGXX:1:11101:24260:1121:CTGTTCTC	0	10	100	36	5S71M40S	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU

# In[ ]:


#test code
import re
import sys
import os

#################################################################################################################
#creating a dictionary with umis

#this is the path to the test dictionary with umis
test_umi_dict="/Users/mandiedriskill/GitHub/2018_Bi621/deduper-mandiedriskill/test_umi_dict.txt"
#creating an empty dictionary for umis
UMI_dict={}
#read the lines in the file
for line in open(test_umi_dict).readlines():
    #stripping the lines in the file
    umi=line.strip()
    #setting the dictionary. Keys are umis and values are none
    UMI_dict[umi]=None
    #check dictionary to make sure it looks right
    #print(UMI_dict)

#################################################################################################################
#creating counters, empty variables and dictionaries

#creating a counter to keep track of lines in the test file
ctr=1

#creating a dictionary to keep track of unique reads
unique_dict={}
#setting the variable to empty
current_rname=""
#setting the variable to empty
previous_rname=""

#creating counters for unmis not in dictionary, forward reads, reverse reads, and unique/deduped reads
ctr_umi_notin_dict=0
forward_reads=0
reverse_reads=0
unique_reads=0
###################################################################################################################
#creating functions

#creating a function for determining if the read is paired-end
def pair_end_reads(flag):
    #if the reads are paired-ended then do stuff
    if((int(flag) & 1) == 1):
        #exiting the code and outputting an error code
        sys.exit("Error, file contains pair-end reads")
    else:
        #if the reads is not paired-end then keep doing stuff
        return True 
    
#creating a function for adjusting the forward postition and making a unique string
def forward_adjust(pos):
    #setting variable to the left S
    softclip=re.search(r"^(\d+)S", cigar)
    # if not S is found then do stuff 
    if softclip is None:
        #softclip does nothing if it has no value
        softclip=softclip
    else:
        #if softclip has value then it is set to group 1 integer
        s1=int(softclip.group(1))
        #the position is adjusted if there is a integer
        pos=pos-s1
    
    #spring is set to a combination of rname,UMI,POS,and Plus
    string=current_rname+UMI+str(pos)+"Plus"
    #returning the funciton
    return string

#creating a function for adjusting the reverse postition and making a unique string
def reverse_adjust(pos):
    #the position is adjusted with the following values
    pos=pos+total_S+total_M+total_D+total_N-1

    #string is set to a combinaiton of rname,UMI,POS, and Minus
    string=current_rname+UMI+str(pos)+"Minus"
    #returning the function
    return string

###################################################################################################################
#opening the test files and writing to output files.

file = "/Users/mandiedriskill/GitHub/2018_Bi621/deduper-mandiedriskill/test_sam.txt"
o = "test_file"

with open (file, "r") as fh, open(o + '_deduped.sam', 'w') as PCR_removed, open(o + '_undetermine.sam', 'w') as undetermined, open(o + '_count_data', 'w') as count_data:
    
###################################################################################################################
#splitting the line into parts and keeping track of line to erase dictionary when a new chromosome is encountered

    
    for line in fh:
        #if line starts with @ do stuff
        if line.startswith("@"):
            #writing the @ lines to the pcr_removed file
            PCR_removed.write(line)
            continue
        else:
            #splitting the line by tabs so parts can be grabbed
            parts = line.split("\t")#this splits all the colums in the line into lists(parts)
            #making sure parts work
            #print(parts)
            #qname is parts 0
            qname=parts[0]
            #making sure qname works
            print(qname)
            #flag is parts 1
            flag=parts[1]
            #making sure flag works
            #print(flag)
            #using the paired end function
            pair_end_reads(flag)
            #current rname is parts 2
            current_rname=parts[2]
            #making sure current_rname works
            #print(current_rname)
            
            #if the counter is on line 1 then do stuff
            if ctr==1:
                #setting previous name to current name
                previous_rname=current_rname
                #making sure previous_rname is working
                #print(previous_rname)
            #keeping track of ctr lines    
            ctr+=1
            
            #if the current_rname doesn't equal previous_rname then do stuff
            if current_rname != previous_rname:
                #reinitializing previous_rname to current_rname
                previous_rname=current_rname
                #deleting the unique_dict when a new chromosome is encountered
                del unique_dict
                #reinitializing the dictionary
                unique_dict={}
            
            
            pos=int(parts[3])
            #making sure pos is correct
            #print(pos)
            cigar=parts[5]
            #making sure cigar is correct
            #print(cigar)
            rnext=parts[6] 
            #making sure rnext is correct
            #print(rnext)
            
###################################################################################################################
#isolating and suming Ds, Ms, Is, and end Ss.
            
            #finding all end Ss 
            S = re.findall(r'(\d+)S$', cigar)#finds all end S's
            #making sure S is correct
            #print(S)
            #turning Ss into integers
            S=list(map(int,S))
            #making sure S is correct
            #print(S)
            #summing all the Ss
            total_S = sum([int(i) for i in S])
            #making sure total_S is correct
            #print(total_S)
            
            #finding all Ms
            M = re.findall(r"(\d+)M", cigar)
            #making sure M is correct
            #print(M)
            #turning Ms into integers
            M = list(map(int,M))
            #making sure M is correct
            #print(M)
            #summing all the Ms
            total_M = sum([int(i) for i in M])
            #making sure total_M is correct
            #print(total_M)
            
            
            #finding all Ds
            D = re.findall(r"(\d+)D", cigar)
            #making sure D is correct
            #print(D)
            #turning Ds into integers
            D = list(map(int,D))
            #making sure D is correct
            #print(D)
            #summing all the Ds
            total_D = sum([int(i) for i in D])
            #making sure total_D is correct
            #print(total_D)
            
            #finding all Ns
            N = re.findall(r"(\d+)N", cigar)
            #making sure N is correct
            #print(N)
            #turning Ns into integers
            N = list(map(int,N))
            #making sure N is correct
            #print(N)
            #summing all the Ns
            total_N = sum([int(i) for i in N])
            #making sure total_N is correct
            #print(total_N)
            
            #splitting qname into parts to isolate UMI
            qname_parts=qname.split(":")
            #making sure qname is correct
            #print(qname)
            #Setting the umi to a variable
            UMI=qname_parts[7]
            #making sure UMI is correct
            #print(UMI)
            
###################################################################################################################
                           
            #If the UMI is not in the dictionary then do stuff
            if UMI not in UMI_dict.keys():
                #keeping track of how many UMIs are not in the dictionary
                ctr_umi_notin_dict+=1
                #making sure correct lines are going into the dictionary
                #print(line)
                #writing the reads with bad UMIs to a undetermined file
                undetermined.write(line)
            else:
                #If the UMI is in the dictionary then do stuff
                if UMI in UMI_dict.keys():
                    #making sure correct lines are still there
                    #print(line)

###################################################################################################################
# Forward strand     
                    #if forward stand then do stuff
                    if((int(flag) & 16) != 16):
                        #counting the number for forward reads
                        forward_reads+=1
                        #setting string to the adjusted forward position
                        string=forward_adjust(pos)
                        #making sure pos was adjusted and string is correct
                        #print(string)
                        #making sure corerct line are being used for forward function
                        #print(line)

###################################################################################################################
# reverse strand
                    #if reverse stand then do stuff
                    if((int(flag) & 16) == 16):
                        #counting the number of reverse reads
                        reverse_reads+=1
                        #setting string to the adjusted reverse potition
                        string=reverse_adjust(pos)
                        #making sure pos was adjusted and string is correct 
                        #print(string)
                        #making sure corerct line are being used for reverse function
                        #print(line)
                        
###################################################################################################################
#Filling the dictionary with unique reads and writing them to a deduped file
                    
                    #if unique string not in the dictionary then do stuff
                    if string not in unique_dict:
                        #keeping track of the number of unique/dedpued reads
                        unique_reads+=1
                        #filling the dictionary with unique reads (keys) and values as none
                        unique_dict[string]=None
                        #making sure dicitonary is putting in unique reads and dumping dictionary for new chrom.
                        print(unique_dict)
                        #writing lines that are unique
                        PCR_removed.write(line)
                            

###################################################################################################################
    #writing count information to a count_data file.
    count_data.write("{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n".format("UMI Not Found in Dictionary", ctr_umi_notin_dict,                                                                "Number of Forward Reads", forward_reads,                                                                "Number of Reverse Reads",  reverse_reads,                                                                "Number of Unique Reads Found", unique_reads)) 
#making sure count data is correct.
print("UMI Not Found in Dictionary", ctr_umi_notin_dict)
print("Number of Forward Reads", forward_reads)
print("Number of Reverse Reads",  reverse_reads)
print("Number of Unique Reads Found", unique_reads)


# #test output
# cat test_file_count_data 
# UMI Not Found in Dictionary
# 1
# Number of Forward Reads
# 7
# Number of Reverse Reads
# 3
# Number of Unique Reads Found
# 9
# 

# #test output
# cat test_file_deduped.sam 
# @HD	VN:1.0	SO:unsorted
# @PG	ID:GSNAP	PN:gsnap	VN:2017-10-12	CL:gsnap.avx2 --gunzip -t 26 -A sam -m 5 -d mm10_chroms -D /projects/bgmp/coonrod/mmu/INTEL -s /projects/bgmp/coonrod/mmu/INTEL/mm10_chroms/mm10_chroms.maps/Mus_musculus.GRCm38.89.splicesites.iit --split-output=/projects/bgmp/coonrod/deduper/gsnap//Datset1 /projects/bgmp/coonrod/deduper//Dataset1.fastq_dups.gz
# @SQ	SN:1	LN:195471971
# @SQ	SN:2	LN:182113224
# @SQ	SN:3	LN:160039680
# @SQ	SN:4	LN:156508116
# @SQ	SN:5	LN:151834684
# @SQ	SN:6	LN:149736546
# @SQ	SN:7	LN:145441459
# @SQ	SN:8	LN:129401213
# @SQ	SN:9	LN:124595110
# @SQ	SN:10	LN:130694993
# @SQ	SN:11	LN:122082543
# @SQ	SN:12	LN:120129022
# @SQ	SN:13	LN:120421639
# @SQ	SN:14	LN:124902244
# @SQ	SN:15	LN:104043685
# @SQ	SN:16	LN:98207768
# @SQ	SN:17	LN:94987271
# @SQ	SN:18	LN:90702639
# @SQ	SN:19	LN:61431566
# @SQ	SN:X	LN:171031299
# @SQ	SN:Y	LN:91744698
# @SQ	SN:MT	LN:16299
# 1)NS500451:single-end:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC	0	1	100	36	5S10M5I5M10N5D71M40S	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
# 2)NS500451:paired-end:HWKTMBGXX:1:11101:24260:1121:CTGTTCAG	0	1	100	36	71M	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
# 3)NS500451:umi-in-dictionary:HWKTMBGXX:1:11101:24260:1121:CTGTTCTT	0	3	100	36	71M	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
# 5)NS500451:flag-is-forward:HWKTMBGXX:1:11101:24260:1121:CTGTTCAT	16	5	100	36	71M	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
# 6)NS500451:flag-is-reverse:HWKTMBGXX:1:11101:24260:1121:CTGTTCAA	0	6	100	36	71M	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
# 7)NS500451:forward&has-s:HWKTMBGXX:1:11101:24260:1121:CTGTTCTA	16	6	100	36	40S71M	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
# 8)NS500451:forward&noaction:HWKTMBGXX:1:11101:24260:1121:CTGTTCTC	0	8	100	36	71M40S	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
# 9)NS500451:forward&has-s:HWKTMBGXX:1:11101:24260:1121:CTGTTCTC	16	9	100	36	5S71M40S	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
# 10)NS500451:forward&has-s:HWKTMBGXX:1:11101:24260:1121:CTGTTCTC	0	10	100	36	5S71M40S	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU

# #test output
# cat test_file_undetermine.sam 
# 4)NS500451:umi-not-in-dictionary:HWKTMBGXX:1:11101:24260:1121:TTTTTTTT	16	4	100	36	71M	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
# 

# Directory with files to run real code on
# 
# /projects/bgmp/shared/deduper
# Dataset1.sam  
# Dataset2.sam  
# Dataset3.sam

# Copied files to own directory
# 
# /projects/bgmp/mdriskil
# cp Data*.sam /projects/bgmp/mdriskil/deduper

# First sorting the following sam files with SamTools: Dataset1.sam  Dataset2.sam  Dataset3.sam
# 
# nano samtools_sort.slurm
# 
# ##################################################################################################################
# #!/usr/bin/env bash
# #SBATCH --partition=short    ### Partition
# #SBATCH --job-name=sam_view.job  ### Job Name
# #SBATCH --output=sam_view.job.out    ### File in which to store job output
# #SBATCH --time=0-23:00:00   ### Wall clock time limit in Days-HH:MM:SS
# #SBATCH --nodes=1           ### Number of nodes needed for the job
# #SBATCH --ntasks-per-node=28 ### Number of tasks to be launcged per Node
# 
# module purge
# module load samtools/1.5
# 
# cd /projects/bgmp/mdriskil/deduper
# 
# files=`ls -1 ./*.sam`
# for file in ${files}; do
#          trimmedfilename=`echo ${file} | cut -d "/" -f2 | cut -d "." -f1`;
# samtools sort -o sorted_${trimmedfilename}.sam -O sam ${file};
# 
# done
# exit
# 
# ##################################################################################################################
# 
# chmod 755 samtools_sort.slurm
# sbatch samtools_sort.slurm
# 
# ##################################################################################################################
# 
# output 
# sorted_Dataset1.sam
# sorted_Dataset2.sam
# sorted_Dataset3.sam
# 

# Making sure files sorted correctly
# 
# cat sorted_Dataset1.sam | grep -v "^@" | cut -f3 | uniq -c
# 1013180 2
# 
# cat sorted_Dataset2.sam | grep -v "^@" | cut -f3 | uniq -c
# 1382109 2
# 
# cat sorted_Dataset3.sam | grep -v "^@" | cut -f3 | uniq -c
#  239922 1
# 1278754 2
#  184971 3
#  208866 4
#  183985 5
#  169727 6
#  577855 7
#  266089 8
#  262789 9
#  191431 10
#  604610 11
#  122285 12
#  160159 13
#  131529 14
#  167429 15
#  129414 16
#  194872 17
#   96470 18
#  229876 19
#  104108 X
#     736 Y

# In[ ]:


#Actual code (driskill_deduper.py) going to run on Dataset1.sam  Dataset2.sam  Dataset3.sam

nano driskill_deduper.py

##################################################################################################################


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

with open (file, "r") as fh, open(out_file + '_deduped.sam', 'w') as PCR_removed, open(out_file + '_undetermine.sam', 'w') as undetermined, open(out_file + '_count_data', 'w') as count_data:
    
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
                        print(unique_dict)
                        PCR_removed.write(line)
                            

###################################################################################################################

    count_data.write("{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n".format("UMI Not Found in Dictionary", ctr_umi_notin_dict,                                                                "Number of Forward Reads", forward_reads,                                                                "Number of Reverse Reads",  reverse_reads,                                                                "Number of Unique Reads Found", unique_reads)) 

##################################################################################################################

chmod 755 driskill_deduper.py


# cat STL96.txt
# AACGCCAT
# AAGGTACG
# AATTCCGG
# ACACAGAG
# ACACTCAG
# ACACTGTG
# ACAGGACA
# ACCTGTAG
# ACGAAGGT
# ACGACTTG
# ACGTCAAC
# ACGTCATG
# ACTGTCAG
# ACTGTGAC
# AGACACTC
# AGAGGAGA
# AGCATCGT
# AGCATGGA
# AGCTACCA
# AGCTCTAG
# AGGACAAC
# AGGACATG
# AGGTTGCT
# AGTCGAGA
# AGTGCTGT
# ATAAGCGG
# ATCCATGG
# ATCGAACC
# ATCGCGTA
# ATCGTTGG
# CAACGATC
# CAACGTTG
# CAACTGGT
# CAAGTCGT
# CACACACA
# CAGTACTG
# CATCAGCA
# CATCGTTC
# CCAAGGTT
# CCTAGCTT
# CGATTACG
# CGCCTATT
# CGTTCCAT
# CGTTGGAT
# CTACGTTC
# CTACTCGT
# CTAGAGGA
# CTAGGAAG
# CTAGGTAC
# CTCAGTCT
# CTGACTGA
# CTGAGTGT
# CTGATGTG
# CTGTTCAC
# CTTCGTTG
# GAACAGGT
# GAAGACCA
# GAAGTGCA
# GACATGAG
# GAGAAGAG
# GAGAAGTC
# GATCCTAG
# GATGTCGT
# GCCGATAT
# GCCGATTA
# GCGGTATT
# GGAATTGG
# GGATAACG
# GGCCTAAT
# GGCGTATT
# GTCTTGTC
# GTGATGAG
# GTGATGTC
# GTGTACTG
# GTGTAGTC
# GTTCACCT
# GTTCTGCT
# GTTGTCGA
# TACGAACC
# TAGCAAGG
# TAGCTAGC
# TAGGTTCG
# TATAGCGC
# TCAGGACT
# TCCACATC
# TCGACTTC
# TCGTAGGT
# TCGTCATC
# TGAGACTC
# TGAGAGTG
# TGAGTGAG
# TGCTTGGA
# TGGAGTAG
# TGTGTGTG
# TTCGCCTA
# TTCGTTCG

# created slurm to run on talapas to dedupe sorted_Dataset1.sam sorted_Dataset2.sam sorted_Dataset3.sam
# 
# nano dedupe.slurm
# 
# ##################################################################################################################
# 
# #!/usr/bin/env bash
# #SBATCH --partition=short    ### Partition
# #SBATCH --job-name=dedupejob  ### Job Name
# #SBATCH --output=dedupejob.out    ### File in which to store job output
# #SBATCH --time=0-00:30:00   ### Wall clock time limit in Days-HH:MM:SS
# #SBATCH --nodes=1           ### Number of nodes needed for the job
# #SBATCH --ntasks-per-node=14 ### Number of tasks to be launcged per Node
# 
# cd /projects/bgmp/mdriskil/deduper
# 
# files=`ls -1 ./sorted_*`
# for file in ${files}; do
#          trimmedfilename=`echo ${file} | cut -d "/" -f2 | cut -d "_" -f2 | cut -d "." -f1`;
# ./driskill_deduper.py -f ${file} -o ${trimmedfilename} -u STL96.txt;
# done
# exit
# 
# chmod 755 dedupe.slurm
# sbatch dedupe.slurm
# 
# ###################################################################################################################
# output files
# 
# Dataset1_count_data
# Dataset1_deduped.sam
# Dataset1_undetermine.sam
# 
# Dataset2_count_data
# Dataset2_deduped.sam
# Dataset2_undetermine.sam
# 
# Dataset3_count_data
# Dataset3_deduped.sam
# Dataset3_undetermine.sam

# Number of lines that were unique/deduplicated
# cat Dataset1_deduped.sam | grep -v "^@" | wc -l
# 635985
# 
# cat Dataset2_deduped.sam | grep -v "^@" | wc -l
# 744188
# 
# cat Dataset3_deduped.sam | grep -v "^@" | wc -l
# 4053644
# ##################################################################################################################
# 
# 
# Number of lines that that had umis that were not in the dictionary
# cat Dataset1_undetermine.sam | grep -v "^@" | wc -l
# 8190
# 
# cat Dataset2_undetermine.sam | grep -v "^@" | wc -l
# 9319
# 
# cat Dataset3_undetermine.sam | grep -v "^@" | wc -l
# 48939
# ##################################################################################################################
# 
# 
# Count data for all the files
# cat Dataset1_count_data 
# UMI Not Found in Dictionary
# 8190
# Number of Forward Reads
# 1004945
# Number of Reverse Reads
# 45
# Number of Unique Reads Found
# 635985
# 
# cat Dataset2_count_data 
# UMI Not Found in Dictionary
# 9319
# Number of Forward Reads
# 1372737
# Number of Reverse Reads
# 53
# Number of Unique Reads Found
# 744188
# 
# cat Dataset3_count_data 
# UMI Not Found in Dictionary
# 48939
# Number of Forward Reads
# 5658335
# Number of Reverse Reads
# 13876
# Number of Unique Reads Found
# 4053644
# 
# 
# 
