# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 21:28:10 2016

downsample pacbio data by assign each read a probablity based on the read length

@author: Nan
"""
from Bio import SeqIO
import numpy as np
from scipy.stats import lognorm
import matplotlib.pyplot as plt
import random

random.seed(0)
target_coverage = 15
genome_length = 4641652
target_read_length = genome_length*target_coverage


records=SeqIO.parse("D:/Data/20161125/filtered_subreads_first1k.fastq", "fastq")

read_length = []
read_dict = {}

for record in records:
    read_length.append(len(record.seq))
    read_dict[record.id] = False
    
data = np.array(read_length)
read_length_sum = np.sum(data)
print read_length_sum

ratio = float(target_read_length)/read_length_sum
    
sigma, loc, scale = lognorm.fit(data, floc=0)

# print sigma, loc, scale
# print lognorm.mean(sigma, loc=loc, scale=scale)


read_length_count = 0
"""
y_value = lognorm.pdf(data, sigma, loc, scale)
background = np.median(y_value)
"""



end_point = lognorm.interval(0.5, sigma, loc, scale)
print end_point

# calculate the homogenesous distribution as a comparable reference
background = 0.5/(end_point[1] - end_point[0])
print background


record_dict = SeqIO.index("D:/Data/20161125/filtered_subreads_first1k.fastq", "fastq")
target_seq = []
i= 0
id_list = list(record_dict.keys())
seq_num = len(id_list)
while read_length_count <= target_read_length:
    print read_length_count
    if i == seq_num:
        i = 0
    record_id = id_list[i]
    record = record_dict[record_id]
    rand = random.random()
    
    if (not(read_dict[record.id])):
        #print "haha"
        print record.id
        dist_value = lognorm.pdf(len(record.seq), sigma, loc, scale)
        if 0<=rand<= ratio*dist_value/background:
            # print ratio*dist_value/background
            target_seq.append(record)
            read_dict[record.id] = True
            read_length_count += len(record.seq)
        
    i += 1

        
SeqIO.write(target_seq, "D:/Data/20161229/target.fastq", "fastq")
            

"""
x_fit = np.linspace(data.min(),data.max(),100)
pdf_fitted = lognorm.pdf(x_fit, sigma, loc, scale)

print lognorm.pdf(10000, sigma, loc, scale)

plt.plot(x_fit, pdf_fitted)

plt.show()
"""