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

initial_coverage = 100
traget_coverage = 15

records=SeqIO.parse("D:/Data/20161125/filtered_subreads_first1k.fastq", "fastq")

read_length = []

for record in records:
    read_length.append(len(record.seq))
    
data = np.array(read_length)
    
sigma, loc, scale = lognorm.fit(data, floc=0)

print sigma, loc, scale
print lognorm.mean(sigma, loc=loc, scale=scale)
"""
x_fit = np.linspace(data.min(),data.max(),100)
pdf_fitted = lognorm.pdf(x_fit, sigma, loc, scale)

print lognorm.pdf(10000, sigma, loc, scale)

plt.plot(x_fit, pdf_fitted)

plt.show()
"""