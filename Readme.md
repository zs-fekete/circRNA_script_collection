#Scripts used for detecting circRNAs

This will be replacing my "circRNA-tools". That one stays on for a while but won't be updated


##CircRNA detection

Added scripts:

###circSTAR_v3

circRNA detector script adapted for rabbit, with added function of returning circRNA supporting reads. I recommend reading Smid et al. (2019) to understand this one.

+ Original script: https://bitbucket.org/snippets/MSmid/Le949d/identify-circularrna-reads
+ Original publication: Smid M, Wilting SM, Uhr K, et al. The circular RNome of primary breast cancer. Genome Res. 2019;29(3):356-366. doi:10.1101/gr.238121.118


###Quantification

####CircSTAR_readcounts

Collects and counts circRNA supporting reads across samples as well as counting how many samples in total contained at least one supporting read trio.
Takes an arbitrary amount of circSTAR output files as input
