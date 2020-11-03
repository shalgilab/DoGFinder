#!/usr/bin/env python2
from __future__ import division
import pybedtools
import numpy as np
import numpy 
import sys
import pysam
import pandas as pd
import itertools
import os.path
import argparse
import subprocess
import shutil
import commands
from multiprocessing import Lock, Process, Queue, current_process, cpu_count
# Function to create annotation file later to be used for bedtools coverage
# create annotation for coverage: for window size coverage:
def Create_cov_annotation(window, new_ann,zero_ref,max_dog_length ):

	#Downstream of gene annotation:
	df_DoG_annotation_window=pd.DataFrame(data=new_ann);
	df_DoG_annotation_window=df_DoG_annotation_window.set_index([new_ann[3]]);
	df_DoG_annotation_combine_windows=pd.DataFrame();

	#+strand
	df_DoG_annotation_window.loc[df_DoG_annotation_window[5]=='+',1]=df_DoG_annotation_window.loc[df_DoG_annotation_window[5]=='+',2]+zero_ref
	df_DoG_annotation_window.loc[df_DoG_annotation_window[5]=='+',2]=df_DoG_annotation_window.loc[df_DoG_annotation_window[5]=='+',1]+window

	#-strand
	df_DoG_annotation_window.loc[df_DoG_annotation_window[5]=='-',2]=df_DoG_annotation_window.loc[df_DoG_annotation_window[5]=='-',1]-zero_ref
	df_DoG_annotation_window.loc[df_DoG_annotation_window[5]=='-',1]=df_DoG_annotation_window.loc[df_DoG_annotation_window[5]=='-',2]-window
	df_DoG_annotation_window[3] =df_DoG_annotation_window[3]+'@'+str(1)
	df_DoG_annotation_combine_windows=[df_DoG_annotation_combine_windows,df_DoG_annotation_window];
	df_DoG_annotation_combine_windows= pd.concat(df_DoG_annotation_combine_windows)
	#loop to create all annotations: with overlap:
	move_end=window+int(window/2);
	k=2;

	while move_end<=max_dog_length:
		df_DoG_annotation_window.loc[df_DoG_annotation_window[5]=='+',1]=df_DoG_annotation_window.loc[df_DoG_annotation_window[5]=='+',1]+int(window/2)
		df_DoG_annotation_window.loc[df_DoG_annotation_window[5]=='+',2]=df_DoG_annotation_window.loc[df_DoG_annotation_window[5]=='+',2]+int(window/2)
		df_DoG_annotation_window.loc[df_DoG_annotation_window[5]=='-',1]=df_DoG_annotation_window.loc[df_DoG_annotation_window[5]=='-',1]-int(window/2)
		df_DoG_annotation_window.loc[df_DoG_annotation_window[5]=='-',2]=df_DoG_annotation_window.loc[df_DoG_annotation_window[5]=='-',2]-int(window/2)
		df_DoG_annotation_window[3] = df_DoG_annotation_window[3].replace(to_replace='@'+str(k-1), value='@'+str(k), regex=True)
		df_DoG_annotation_combine_windows=[df_DoG_annotation_combine_windows,df_DoG_annotation_window];
		df_DoG_annotation_combine_windows = pd.concat(df_DoG_annotation_combine_windows)
		k=k+1;
		move_end=move_end+int(window/2);
	df_DoG_annotation_combine_windows.loc[df_DoG_annotation_combine_windows[1]<0,1]=0
	df_DoG_annotation_combine_windows.loc[df_DoG_annotation_combine_windows[2]<0,2]=0
	t1=df_DoG_annotation_combine_windows[2]==0
	t2=df_DoG_annotation_combine_windows[1]==0
	df_DoG_annotation_combine_windows.loc[t1 &t2 ,2]=1

	return df_DoG_annotation_combine_windows;

# create annotation for coverage: for tot coverage
def Create_cov_annotation_tot(window, new_ann,zero_ref,max_dog_length ):

	#Downstream of gene annotation:
	df_DoG_annotation_window=pd.DataFrame(data=new_ann);
	df_DoG_annotation_window=df_DoG_annotation_window.set_index([new_ann[3]]);
	df_DoG_annotation_combine_windows=pd.DataFrame();

	#+strand
	df_DoG_annotation_window.loc[df_DoG_annotation_window[5]=='+',2]=df_DoG_annotation_window.loc[df_DoG_annotation_window[5]=='+',2]+window+zero_ref

	#-strand
	df_DoG_annotation_window.loc[df_DoG_annotation_window[5]=='-',1]=df_DoG_annotation_window.loc[df_DoG_annotation_window[5]=='-',1]-window-zero_ref

	df_DoG_annotation_window[3] =df_DoG_annotation_window[3]+'@'+str(1)
	df_DoG_annotation_combine_windows=[df_DoG_annotation_combine_windows,df_DoG_annotation_window];
	df_DoG_annotation_combine_windows= pd.concat(df_DoG_annotation_combine_windows)
	#loop to create all annotations: with overlap:
	move_end=window+window;
	k=2;

	while move_end<=max_dog_length:
		df_DoG_annotation_window.loc[df_DoG_annotation_window[5]=='+',2]=df_DoG_annotation_window.loc[df_DoG_annotation_window[5]=='+',2]+int(window)
		df_DoG_annotation_window.loc[df_DoG_annotation_window[5]=='-',1]=df_DoG_annotation_window.loc[df_DoG_annotation_window[5]=='-',1]-int(window)
		df_DoG_annotation_window[3] = df_DoG_annotation_window[3].replace(to_replace='@'+str(k-1), value='@'+str(k), regex=True)
		df_DoG_annotation_combine_windows=[df_DoG_annotation_combine_windows,df_DoG_annotation_window];
		df_DoG_annotation_combine_windows = pd.concat(df_DoG_annotation_combine_windows)
		k=k+1;
		move_end=move_end+int(window);
	df_DoG_annotation_combine_windows.loc[df_DoG_annotation_combine_windows[1]<0,1]=0
	df_DoG_annotation_combine_windows.loc[df_DoG_annotation_combine_windows[2]<0,2]=0
	t1=df_DoG_annotation_combine_windows[2]==0
	t2=df_DoG_annotation_combine_windows[1]==0
	df_DoG_annotation_combine_windows.loc[t1 &t2 ,2]=1

	return df_DoG_annotation_combine_windows;


#This function return the DoGs that have ended and the ones for next loop:
def DoGs_annotation( df_annotation_combine_windows_coverage,window,wind_cut,df_annotation,moving_window_size_ds,new_ann,zero_ref,df_max_dog_length,new_old,max_id_k,max_id):

	df_annotation_combine_windows_coverage['names'], df_annotation_combine_windows_coverage['id'] = zip(*df_annotation_combine_windows_coverage[3].map(lambda x: x.split('@')))
	df_annotation_combine_windows_coverage.id = df_annotation_combine_windows_coverage.id.astype(int)
	df_annotation_cov = pd.DataFrame({'names': df_annotation_combine_windows_coverage['names'], 'id' : df_annotation_combine_windows_coverage['id'],'coverage' : 	df_annotation_combine_windows_coverage[9]})
	df_DoG_annotation_trans=df_annotation_cov.pivot(index='names', columns='id',values='coverage');

	numpyMatrix = df_DoG_annotation_trans.as_matrix()
	location_DoG_stop = np.empty((len(numpyMatrix)))
	location_DoG_stop[:] = np.NAN

	for k in range(0, len(numpyMatrix)):
		if sum(numpyMatrix[k,:]<wind_cut)>0:
			location_DoG_stop[k]=np.argmax(numpyMatrix[k,:]<wind_cut)
		if max_id_k>=0:
			location_DoG_stop[k]=max_id_k
	#DoGs end coordination:
	
	DoGs_df=pd.DataFrame(location_DoG_stop,index=df_DoG_annotation_trans.index)

	idx=DoGs_df.isnull().values
	DoGs_df['size_f']=moving_window_size_ds.loc[DoGs_df.index,12]+max_id_k*int(window)+zero_ref
	DoGs_df.loc[DoGs_df['size_f']>df_max_dog_length.loc[DoGs_df.index.values,'size'],0]=max_id
	DoGs_df=pd.DataFrame(DoGs_df[0],index=DoGs_df.index)
	DoGs_end_df=pd.DataFrame(DoGs_df[DoGs_df[0].notnull()])
	DoGs_next_loop_df=DoGs_df[DoGs_df[0].isnull()].index.values

	if len(DoGs_end_df)>0:
		DoGs_end_df=pd.DataFrame(DoGs_end_df,index=DoGs_end_df.index)
		DoGs_end_df.columns=['size'];
		if not new_old:
			DoGs_end_df['size']=moving_window_size_ds.loc[DoGs_end_df.index,12]+DoGs_end_df['size']*int(window)+zero_ref
			DoGs_end_df.loc[DoGs_end_df['size']>df_max_dog_length.loc[DoGs_end_df.index.values,'size'],'size']=df_max_dog_length.loc[DoGs_end_df.index.values,'size']
			
		else :
			DoGs_end_df['size']=moving_window_size_ds.loc[DoGs_end_df.index,12]+DoGs_end_df['size']*int(window/2)+zero_ref

		df_DoGs_end_annotations=pd.concat([df_annotation.loc[DoGs_end_df.index],DoGs_end_df],axis=1);
		df_DoGs_end_annotations.loc[df_DoGs_end_annotations['strand']=='+','start']=df_DoGs_end_annotations.loc[df_DoGs_end_annotations['strand']=='+','end']
		df_DoGs_end_annotations.loc[df_DoGs_end_annotations['strand']=='+','end']=df_DoGs_end_annotations.loc[df_DoGs_end_annotations['strand']=='+','start']+df_DoGs_end_annotations.loc	[df_DoGs_end_annotations['strand']=='+','size']

		df_DoGs_end_annotations.loc[df_DoGs_end_annotations['strand']=='-','end']=df_DoGs_end_annotations.loc[df_DoGs_end_annotations['strand']=='-','start']
		df_DoGs_end_annotations.loc[df_DoGs_end_annotations['strand']=='-','start']=df_DoGs_end_annotations.loc[df_DoGs_end_annotations['strand']=='-','end']-df_DoGs_end_annotations.loc[df_DoGs_end_annotations['strand']=='-','size']

		df_DoGs_end_annotations=pd.DataFrame(df_DoGs_end_annotations[df_DoGs_end_annotations.columns[:-1]])
	else:
		df_DoGs_end_annotations=pd.DataFrame()
	if len(DoGs_next_loop_df)>0:
		df_annot_next=pd.DataFrame(new_ann.loc[DoGs_next_loop_df])
	else:
		df_annot_next=pd.DataFrame()

	return df_annot_next,df_DoGs_end_annotations;



# Examine bam file
def Examine_bam(bam_files_names,ref_ann):
	#print "Examine bam file: %s "%(bam_files_names)
	cmd = 'infer_experiment.py -r '+ref_ann+' -i '+bam_files_names
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
	out, err = p.communicate()
	result = out.split('\n')
	count=0;
	S= numpy.zeros(shape=(3,1))
	i=1;
	PE='False'
	for lin in result:
		if i==3 and 'PairEnd' in lin:
			PE='True'
		if not lin.startswith('#'):
			if len(lin.split(":",1))>1:
				S[count]= float(lin.split(":",1)[1])
				count=count+1
		i=i+1
	if S[1]> 0.6:
		strand='first'
	#	print "Leading strand is: %s strand"%(strand)
	elif S[2]> 0.6:
		strand='second'
	#	print "Leading strand is: %s strand"%(strand)
	else:
		strand='non'
	#	print 'Not strand specific'

	if PE=='False':
		PE=='False'
	#	print 'Not Paired End file'
    

	return PE,strand;


# PE modification
def PE_Mod(PE,strand,path_bam,output_path):
	if PE=='True':
		#print "\nProcessing PE bam file...may take some time" 
		base_name=os.path.basename(path_bam)    
		out_name=output_path+"/"+base_name.split(".bam")[0]+".sorted_PE.bam";
		infile = pysam.AlignmentFile(path_bam, "rb")
		sam_strng=commands.getstatusoutput('samtools')
		temp=sam_strng[1].split("Version:")[1][:7]
		sam_version=int(filter(str.isdigit, temp))
		if sam_version==119:
			outfile = pysam.AlignmentFile(output_path+"/"+base_name.split(".bam")[0]+"_temp.bam", "wb", template=infile)
		else:
			outfile = pysam.AlignmentFile(output_path+"/"+base_name.split(".bam")[0]+"_temp.bam", "wb",add_sam_header=True, template=infile)



		for read in infile:
			Flag=read.flag;
			if strand=='second':
				if (Flag==145 or Flag==147 or Flag==153 or Flag==97 or Flag==99  or Flag==73 ):
					if Flag<100 :
						read.flag=145
					read.set_tag("XS",'-', replace=True)
					outfile.write(read)
				if (Flag==161 or Flag==163 or Flag==81 or Flag==83 or Flag==89  or Flag==137 ):
					if Flag<100 :
						read.flag=161
					read.set_tag("XS",'+', replace=True)
					outfile.write(read)
			if strand=='first' or 'non':      
				if (Flag==145 or Flag==147 or Flag==153 or Flag==97 or Flag==99  or Flag==73 ):
					if Flag>100 :
						read.flag=97
					read.set_tag("XS",'+', replace=True)
					outfile.write(read)
				if (Flag==161 or Flag==163 or Flag==81 or Flag==83 or Flag==89  or Flag==137 ):
					if Flag>100 :
						read.flag=81
					read.set_tag("XS",'-', replace=True)       
					outfile.write(read)

		infile.close()
		outfile.close()		

		#print "Sorting PE bam file" 
		if sam_version==119:
			pysam.sort("-o", out_name, output_path+"/"+base_name.split(".bam")[0]+"_temp.bam")
		else:
			pysam.sort("-O", "BAM","-T",output_path+"/"+base_name.split(".bam")[0]+"_temp.sorted","-o", out_name, output_path+"/"+base_name.split(".bam")[0]+"_temp.bam")

		pysam.index(out_name)
		os.remove(output_path+"/"+base_name.split(".bam")[0]+"_temp.bam")

	if PE=='False':
		out_name=path_bam

	return out_name

def count_reads(path_bam,output_path):
	base_name=os.path.basename(path_bam)    
	path_temp=output_path+"/"+base_name.split(".bam")[0]+".temp.sorted";
	path_new_sort=output_path+"/"+base_name.split(".bam")[0]+".sorted.bam";

	if  not(os.path.isfile(path_bam+".bai")):
	        print "No index file at : ",path_bam+".bai... Sorting and indexing..."
	        pysam.sort("-O", "BAM","-T",path_temp,"-o",path_new_sort,path_bam )
	        pysam.index(path_new_sort)
	        path_bam=path_new_sort
	idx_stat = pysam.idxstats(path_bam)
	stat_string=idx_stat.split("\n")
	count=0;
	for k in range(0,len(stat_string)):
		temp=stat_string[k].split("\t")
		if len(stat_string[k].split("\t")) > 2:
			count=count+int(temp[2])	
	number_lines = count;
	tot_reads_m=number_lines/1000000;
#	print "Bam file: %s has %s mapped Mreads" %(path_bam,str(tot_reads_m))
	return tot_reads_m


def Downsample(k,path_bam,output_path,reads_min_count):

	base_name=os.path.basename(path_bam)  
	path_new=output_path+"/"+base_name.split(".bam")[0]+".DS.bam";
	path_temp=output_path+"/"+base_name.split(".bam")[0]+".temp.sorted";
	path_new_sort=output_path+"/"+base_name.split(".bam")[0]+".sorted_DS.bam";
	reads_count=count_reads(path_bam,output_path)
	s=reads_min_count/reads_count
	string=str(s)
	s=string.split('.')[1]

	if reads_count==reads_min_count:
		shutil.copy2(path_bam, path_new_sort)
		#os.rename(path_bam, path_new_sort)
		pysam.index(path_new_sort)

	else :
#		pysam.view("-b","-s",str(k)+'.'+str(s),"-O","BAM","-o",path_new,path_bam,catch_stdout=False)
		pysam.view("-b","-s",str(k)+'.'+str(s),"-O","BAM","-o",path_new,path_bam,catch_stdout=False)
		pysam.sort("-O", "BAM","-T",path_temp,"-o",path_new_sort,path_new )
		pysam.index(path_new_sort)
		os.remove(path_new)
	

	#print "Downsampled Bam file at: %s " %(path_new_sort)



