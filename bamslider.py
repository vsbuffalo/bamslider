import pysam
from collections import namedtuple, deque, counter



class BamSlider(object):
    """BamSlider is a simple class that uses PySam to create a windowed
    iterator (slider) across a sorted SAM/BAM file. This iterator
    yields a tuples (window_info, reads) where window_info is a named
    tuple indicating the start and end positions, and chromosome name.

    Groups of reads in
    a window will be yielded; note that this is not a pileup, which
    are all reads covering a particular base. Also, this is unlike a
    moving average in that reads will not be processed twice, even if
    they overlap the window.

    Note that there are some differing options in how we process
    windows.

    # Read Pairs and Improperly Paired Reads
    Read pairs are processed by the left-most mate, as long as they
    are proper pairs. The motivation here is that for most statistics,
    we care about processing reads in discrete fragments. If a
    left-most read's pair falls outside the window, it is still
    treated as belonging to this window (note that improper pairs are
    the exception here; improper pairs are processed as two single
    reads).

    # Overlapping Alignments    
    In some cases, our shift length could be less than our window
    size. In these cases, we could end up with cases where our window
    leads to redundancy, as a read overlaps both windows. To prevent
    this, we simply require that the left-most start position of the
    read by within the shift length away from the window start only.

    """
    
    def __init__(self, filename, size, mode=None):
        """
        Initialize a BamSlider by opening a SAM/BAM file with PySam.
        """            
        if size < 1:
            raise ValueError, "shift size must be >= 1"
        self.filename = filename
        self.size = size
        self.shift = shift
        self._mode = mode
        self._samfile = pysam.Samfile(filename, mode)
        self._deque = deque()
        self._type_counter = counter()
    
    def __iter__(self):
        """
        
        """
        return self

    def next(self):
        """
        
        """
        last_pos = -1
        seen_qnames = set()
        for read in self._samfile:
            if read.is_unmapped:
                continue
            if read.is_proper_pair:
                
        
        

if __name__ == "__main__":
    pass
