import pysam
from collections import namedtuple, deque, Counter
import pdb

Window = namedtuple('Window', ("seqname", "start", "end"))

class BamSlider(object):
    """BamSlider is a simple class that uses PySam to create a windowed
    iterator (slider) across a sorted SAM/BAM file. This iterator
    yields a tuples (window_info, reads) where window_info is a named
    tuple indicating the start and end positions, and chromosome name,
    and reads is a list of all reads with a start position in this
    range.

    It makes sense to have BamSlider's window size longer than the
    longest expected paired-end fragment size so that both pairs are
    included in window. PySam's Samfile.mate() is not particulary
    fast, so we hash all sequences in window by readname in an
    OrderedDict, to group pairs together.
    """
    
    def __init__(self, filename, size, step, mode=None):
        """Initialize a BamSlider with window length `size` and step size
        `step`.
        """            
        if size < 1:
            raise ValueError, "shift size must be >= 1"
        self.filename = filename
        self.size = size
        self.step = step
        self._mode = mode
        self._samfile = pysam.Samfile(filename, mode)
        self._deque = deque()
        self._type_counts = Counter()
        
        # for storing state in iterator
        self._last_pos = None
        self._last_tid = -1
        self._window_starts = None
        self._wstart = None
        self._wend = None

    
    def __iter__(self):
        """
        Iterator method. 
        """
        return self._windows()

    def _windows(self):
        """This generator method generates aligned reads in window sizes
        specified by `size` and `step`.

        Caution: beware of edge effects with large windows. If you
        have a 10kb window that jumps in 10kb chunks, the last chunk
        may less than 10kb wide, leadering to noise.

        """
        for read in self._samfile:
            new_chromosome = (self._last_pos < 0 and self._last_tid is None) or self._last_tid != read.tid
            if new_chromosome:
                # if the deque isn't empty, dump it via yield
                if len(self._deque):
                    this_window = list(self._deque)
                    self._deque = deque()
                    yield (Window(self._samfile.references[read.tid], self._wstart, self._wend), this_window)

                # initialize sequence window for the length of this
                # chromosome
                self._window_starts = iter(xrange(0, self._samfile.lengths[read.tid]-self.size, self.step))
                self._wstart = next(self._window_starts)
                self._wend = self._wstart + self.size
                self._last_pos = read.pos
                self._last_tid = read.tid
            
            if read.is_unmapped:
                self._type_counts["unmapped"] += 1
                continue

            # always check things are sorted
            if self._last_pos > read.pos:
                raise ValueError, "BAM/SAM file not sorted"
            self._last_pos = read.pos

            # look at where we are, and if we see a read that has a
            # start position that's past this window, yield the
            # current results.
            #
            # Note that we can only left pop off the deque up to the
            # start of the next window.
            while read.pos > self._wend:
                # keep yielding windows until we hit a read where this
                # isn't the case.

                # build list of reads for this window
                this_window = list()
                # we can't modify deque with popleft while iterating it 
                num_remove = 0
                for i in xrange(len(self._deque)):
                    dequed_read = self._deque[i]
                    if dequed_read.pos < self._wstart + self.step:
                        # remove reads with start positions that are before next start
                        num_remove += 1
                    this_window.append(dequed_read)
                for i in range(num_remove):
                    self._deque.popleft()
                    
                last_wstart, last_wend = self._wstart, self._wend
                self._wstart = next(self._window_starts)
                self._wend = self._wstart + self.size
                yield (Window(self._samfile.references[read.tid], last_wstart, last_wend), this_window)

            # append this read to the end of the deque
            self._deque.append(read)

if __name__ == "__main__":
    slider = BamSlider("test.sorted.sam", 10000, 5000)
    for window, reads in slider:
        print window, [r.mapq for r in reads]
    
