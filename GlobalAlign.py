SD = {"A": 0, "C": 1, "G": 2, "T": 3, "GAP": 4}  # scores dictionary


class GlobalAlign:
    """
    Implements global alignment. For my considerations, picked to be the
    father class of alignments (mostly in order to prevent code duplication)
    """

    def __init__(self, s, t, align_type, scores):
        # sequence data
        self._s = s   # rows
        self._t = t   # cols
        self._type = align_type
        # table data
        self._rowSize = len(s) + 1  # +1 for the starting gap
        self._colSize = len(t) + 1
        self._scores = scores     # pandas' data frame
        self._V = []  # will be a 2D table
        # alignment printing data
        self._bestScore = 0
        self._Ptr = []

    def doAlign(self):
        self._align()
        self._buildSol()

    def _initialVandPtr(self):
        """
        Initialize an zero 2D-table (the rows arrays) and than initialize
        the first row and column of the matrix with gap penalties.
        """
        self._V = [[0 for x in range(self._colSize)] for y in range(self._rowSize)]
        self._Ptr = [["" for x in range(self._colSize)] for y in range(self._rowSize)]
        for j in range(1, self._colSize):
            self._V[0][j] = int(self._scores.at[SD["GAP"], self._t[
                j - 1]]) + self._V[0][j - 1]
            self._Ptr[0][j] = "L"
        for i in range(1, self._rowSize):
            self._V[i][0] = int(self._scores.at[SD["GAP"], self._s[
                i - 1]]) + self._V[i - 1][0]
            self._Ptr[i][0] = "U"

    def _align(self):
        """
        Implements dynamic algorithm and fills matrix V in order to find the
         best score and trace at position i,j of sequences in length n, m
        """
        self._initialVandPtr()
        for i in range(1, self._rowSize):
            sChar = self._s[i - 1]
            for j in range(1, self._colSize):
                tChar = self._t[j - 1]
                # Each variable is a tuple of (highest value, The direction
                # which it came from). D for diagonal, U for up, L for left.
                match = (self._V[i - 1][j - 1] + self._scores.at[SD[
                                                                   sChar], tChar], "D")
                sAndGap = (self._V[i - 1][j] + self._scores.at[4, sChar], "U")
                tAndGap = (self._V[i][j - 1] + self._scores.at[4, tChar], "L")
                maxVal = max(match, sAndGap, tAndGap, key=lambda tup:tup[0])
                self._V[i][j] = maxVal[0]
                self._Ptr[i][j] = maxVal[1]


    def _buildSol(self):
        """
        Extracts the solution from the score matrix V with the help of
        the trace matrix Ptr. Creates two strings which represents the two
        original sequences as they put in the alignment, including gaps.
        """
        resS, resT = "", ""
        i, j = self._rowSize - 1, self._colSize - 1
        self._bestScore = self._V[i][j] # Bottom  right corner of the matrix.
        # Trace the path back and creates the solutions from end do start.
        while i > 0 and j > 0:
            if self._Ptr[i][j] == "D":
                resS += self._s[i - 1]
                resT += self._t[j - 1]
                i, j = i-1, j-1
            elif self._Ptr[i][j] == "U":
                resS += self._s[i - 1]
                resT += "-"
                i = i - 1
            else:
                resT += self._t[j - 1]
                resS += "-"
                j = j - 1
        # complete the other string
        while i > 0:
            resS += self._s[i - 1]
            resT += '-'
            i -= 1
        while j > 0:
            resS += '-'
            resT += self._t[j - 1]
            j -= 1

        resS = resS[::-1]
        resT = resT[::-1]
        self.printAlign(resS, resT, self._bestScore, self._type)

    def printAlign(self, s, t, score, atype):
        for i in range(0, len(s), 50):
            print(''.join(s[i:i + 50]))
            print(''.join(t[i:i + 50]))
            print()
        print(atype + ": ", score)
