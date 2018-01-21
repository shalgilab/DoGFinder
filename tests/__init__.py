#!/usr/bin/python
import commands
import subprocess
import os

print "Running tests:\n"


path_to_dir = os.getcwd()
bam_file_HS = path_to_dir+'/tests/data/HS.bam'
bam_file_sorted = path_to_dir+'/tests/data/HS.sorted_PE.sorted_DS.bam'
bam_file_sorted_un = path_to_dir+'/tests/data/UN.sorted_PE.sorted_DS.bam'

bam_file_UN = path_to_dir+'/tests/data/UN.bam'
gtf_file=path_to_dir+'/tests/data/test_gtf.gtf'
luci=path_to_dir+'/tests/data/loci_annotation.bed'



path_dog_ann=path_to_dir+'/tests/data/union_dog_annotation.bed'


Get_luci_annotation_cmd="Get_loci_annotation -out "+path_to_dir+"/tests/data -gtf "+gtf_file
print Get_luci_annotation_cmd
os.system(Get_luci_annotation_cmd) 


Paired_End_Mod_cmd="Pre_Process -bam "+bam_file_HS+","+bam_file_UN+" -ref "+luci
print Paired_End_Mod_cmd
os.system(Paired_End_Mod_cmd) 
 

Get_DoGs_cmd="Get_DoGs -out "+path_to_dir+"/tests/data -bam "+bam_file_sorted+" -a "+luci+" -suff HS "
print Get_DoGs_cmd
os.system(Get_DoGs_cmd) 

Get_DoGs_cmd="Get_DoGs -out "+path_to_dir+"/tests/data -bam "+bam_file_sorted_un+" -a "+luci+" -suff UN "
print Get_DoGs_cmd
os.system(Get_DoGs_cmd)


Get_DoGs_cmd="Union_DoGs_annotation -dog "+path_to_dir+"/tests/data/Final_Dog_annotation_HS.bed,"+path_to_dir+"/tests/data/Final_Dog_annotation_UN.bed -out "+path_to_dir+"/tests/data"
print Get_DoGs_cmd
os.system(Get_DoGs_cmd)

Get_DoGs_cmd="Common_DoGs_annotation -comm "+path_to_dir+"/tests/data/Final_Dog_annotation_HS.bed,"+path_to_dir+"/tests/data/Final_Dog_annotation_UN.bed -out "+path_to_dir+"/tests/data"
print Get_DoGs_cmd
os.system(Get_DoGs_cmd)


Get_DoGs_rpkm_cmd="Get_DoGs_rpkm -out "+path_to_dir+"/tests/data -bam "+bam_file_sorted+" -dog "+path_dog_ann+" -s "
print Get_DoGs_rpkm_cmd
os.system(Get_DoGs_rpkm_cmd) 


