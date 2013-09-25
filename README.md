# BamSlider

An **experimental** sliding window for BAM/SAM files using Python's
[deques](http://docs.python.org/2/library/collections.html#deque-objects). This
is all in alpha, but I hope to develop it as I get more time.

`BamSlider` is a iterator class. You can use it like so:


    slider = BamSlider("test.sorted.sam", size=10000, step=5000)
    for window, reads in slider:
        print window, [r.mapq for r in reads]

Each `window` here is a named tuple of seqname (i.e. chromosome or
contig name), start position, and end position. `reads` is a list of
all reads with **start positions** in this window. Note that pairing
is not explicitly handled.

## Todo

 - Unit testing -- **high priorty**
 - Maybe add option to handle pairing of properly-paired reads. This
   would require using a dictionary to group reads, for efficiency.
   
