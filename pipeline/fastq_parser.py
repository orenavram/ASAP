import gzip

class FastqParser(object):
    """Returns a read-by-read fastQ parser analogous to file.readline()"""

    def __init__(self, filePath, headerSymbols=['@', '+']):
        """Returns a read-by-read fastQ parser analogous to file.readline().
        Exmpl: parser.next()
        -OR-
        Its an iterator so you can do:
        for rec in parser:
            ... do something with rec ...

        rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
        """
        # if filePath.endswith('.gz'):
        #     self._file = gzip.open(filePath)
        # else:
        self._file = open(filePath)
        self._currentLineNumber = 0
        self._hdSyms = headerSymbols

    def __iter__(self):
        return self

    def __next__(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: tuple: (seqHeader,seqStr,qualHeader,qualStr)"""
        # ++++ Get Next Four Lines ++++
        elemList = []
        for i in range(4):
            line = self._file.readline()
            self._currentLineNumber += 1  ## increment file position
            if line:
                elemList.append(line.strip())
            else:
                elemList.append(None)

        # -- Check for acceptable end of file --
        if elemList.count(None) == 4:
            raise StopIteration
        # ++++ Check Lines For Expected Form ++++
        # -- Make sure we got 4 full lines of data --
        assert all(elemList), \
            f'It looks like I encountered a premature EOF or empty line.\nPlease check FastQ file near line number {self._currentLineNumber} (plus or minus ~4 lines) and try again**'
        # -- Make sure we are in the correct "register" --
        assert elemList[0].startswith(self._hdSyms[0]), \
            f"The 1st line in fastq element does not start with '{self._hdSyms[0]}'.\nPlease check FastQ file near line number {self._currentLineNumber} (plus or minus ~4 lines) and try again**"
        assert elemList[2].startswith(self._hdSyms[1]), \
            f"The 3rd line in fastq element does not start with '{self._hdSyms[1]}'.\nPlease check FastQ file near line number {self._currentLineNumber} (plus or minus ~4 lines) and try again**"
        # -- Make sure the seq line and qual line have equal lengths --
        assert len(elemList[1]) == len(elemList[3]), \
            f"The length of Sequence data and Quality data of the last record aren't equal.\nPlease check FastQ file near line number {self._currentLineNumber} (plus or minus ~4 lines) and try again**"

        # ++++ Return fatsQ data as tuple ++++
        return tuple(elemList)