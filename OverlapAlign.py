from GlobalAlign import GlobalAlign
SD = {"A": 0, "C": 1, "G": 2, "T": 3, "GAP": 4}  # scores dictionary


class OverlapAlign(GlobalAlign):
    """
    Successor of GlobalAlign. Differ from global alignment by the
    initialization of first column, and the fill of the last row in while
    calculate the matrix V. In both of them, penalties for gaps are ignores,
     because we do allow gaps before and after the alignment.
    """
    def __init__(self, s, t, align_type, scores):
        super().__init__(s, t, align_type, scores)

    def _initialVandPtr(self):
        """
        The initialization at overlap requires non-penalty score if we align
        the prefix of s with gaps. So at the first column of the matrix we
        initialize without penalty.
        """
        self._V = [[0 for x in range(self._colSize)] for y in range(self._rowSize)]
        self._Ptr = [["" for x in range(self._colSize)] for y in range(self._rowSize)]
        for j in range(1, self._colSize):
            self._V[0][j] = int(self._scores.at[SD["GAP"], self._t[
                j - 1]]) + self._V[0][j - 1]
            self._Ptr[0][j] = "L"

    def _align(self):
        """
        Here the difference from global alignment is at the last row. we do let
        the t suffix be aligned with gaps, so while we fill the last row,
        we don't add gap penalty to the value.
        :return:
        """
        self._initialVandPtr()
        for i in range(1, self._rowSize):
            sChar = self._s[i - 1]
            for j in range(1, self._colSize):
                tChar = self._t[j - 1]
                match = (self._V[i - 1][j - 1] + self._scores.at[SD[
                                                                   sChar], tChar], "D")
                sAndGap = (self._V[i - 1][j] + self._scores.at[4, sChar], "U")
                tAndGap = (self._V[i][j - 1] + self._scores.at[4, tChar], "L")
                if i == self._rowSize-1:
                    tAndGap = (self._V[i][j-1], "L")
                maxVal = max(match, sAndGap, tAndGap, key=lambda tup:tup[0])
                self._V[i][j] = maxVal[0]
                self._Ptr[i][j] = maxVal[1]
