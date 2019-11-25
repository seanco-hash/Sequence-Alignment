from GlobalAlign import GlobalAlign
SD = {"A": 0, "C": 1, "G": 2, "T": 3, "GAP": 4}  # scores dictionary


class LocalAlign(GlobalAlign):
    """
    Successor of GlobalAlign.. In this alignment type we use additional
    possible values to the matrix cells while we fill it. It is possible
    that the best subsequent will start at each single cell.
    """
    def __init__(self, s, t, align_type, scores):
        super().__init__(s, t, align_type, scores)
        self._maxi, self._maxj = 0, 0

    def _align(self):
        """
        Note that while filling the matrix, we actually find the best
        alignment in each state - cell i,j. So we add an option of calculate
        the cell score alone, and we sign it with the letter F for First -
        if it is the max score at the point, it mean that this cell is the
        first in the sequence.
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
                # New state, single cell may be scored higher alone than as
                # part of other subsequent
                alone = (self._scores.at[SD[sChar], tChar], "F")
                maxVal = max(match, sAndGap, tAndGap, alone, key=lambda \
                        tup:tup[0])
                self._V[i][j] = maxVal[0]
                self._Ptr[i][j] = maxVal[1]
                # Saves the best score and its indexes at the matrix.
                if maxVal[0] > self._bestScore:
                    self._bestScore = maxVal[0]
                    self._maxi = i
                    self._maxj = j

    def _buildSol(self):
        """
        # The difference in local alignment is that we find the max value at
        #  the matrix and trace our way back from it. It will give us the
        # best subsequent alignment
        :return:
        """
        i, j = self._maxi, self._maxj  # the max value's coordinated in the matrix.
        resS, resT = "", ""
        while i > 0 and j > 0:
            if self._Ptr[i][j] == "D":
                resS += self._s[i - 1]
                resT += self._t[j - 1]
                i, j = i-1, j-1
            elif self._Ptr[i][j] == "U":
                resS += self._s[i - 1]
                resT += "-"
                i = i - 1
            elif self._Ptr[i][j] == "F":
                resS += self._s[i - 1]
                resT += self._t[j - 1]
                i, j = 0, 0
            else:
                resT += self._t[j - 1]
                resS += "-"
                j = j - 1

        resS = resS[::-1]
        resT = resT[::-1]
        self.printAlign(resS, resT, self._bestScore, self._type)
