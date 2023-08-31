from source.Matrix import *
from abc import ABC, abstractmethod


class Alignment(ABC):
    def __init__(self, seq1, seq2, I, E, submat, k):
        self.__maxSols = k
        self.__subMat = submat
        self.__solutions = []
        self.__I = I
        self.__E = E
        self.__seq1 = seq1
        self.__seq2 = seq2

    def __getitem__(self, item):
        """
        Return score at specified position
        """
        return self.__S[item]

    @abstractmethod
    def calculate_score(self):
        """
        Calculate the scores of Score Matrix S, according to the values in S, V and W
        """
        pass

    def backtrack(self, i, j, h_seq, v_seq):
        """
        Use the Score Matrices to find the path taken and form the two strings of the alignment
        :param i: int, row index of current position
        :param j: int, col index of current position
        :return: list of tuples of size 5, (str, str, float, (int,int), (int,int)) for seq1 with gaps, seq2 with the
         inserted gaps and score of the alignment and the positions of beginning and end of the alignements sequences
        """
        pass

    def get_solution(self):
        """
        :return: list of tuples of size 5, (str, str, float, (int,int), (int,int)) for seq1 with gaps, seq2 with the
         inserted gaps and score of the alignment and the positions of beginning and end of the alignements sequences
        """
        if self.__solutions is not None:
            return self.__solutions
        self.run()
        return self.__solutions

    def compute_scores(self, alignment):
        """
        Compute the identity, similarity and number of gaps of an alignment
        :param alignment: list of tuples of size 5, (str, str, float, (int,int), (int,int)) for seq1 with gaps, seq2
        with the inserted gaps and score of the alignment and the positions of beginning and end of the alignments sequences
        :return: list of tuples of 3 floats (rounded two decimal places) respectively identity, similarity and gaps rates (in %)
        """
        results = list()
        for Tuple in alignment:
            seq1, seq2 = Tuple[0], Tuple[1]
            length = len(seq1)

            gaps_seq1 = 0
            gaps_seq2 = 0
            matches = 0
            similarities = 0
            for i in range(length):
                # A bit ugly but whatever
                if seq1[i] == seq2[i]:
                    matches += 1
                    similarities += 1
                elif seq1[i] == "-":
                    gaps_seq1 += 1
                elif seq2[i] == "-":
                    gaps_seq2 += 1
                else:
                    # Get value in S. Similar only if subMat(i, j) > 0
                    similar = self.__subMat[seq1[i], seq2[i]]
                    similarities += 1 if similar > 0 else 0

            results.append((round(matches / length * 100, 2), round(
                similarities / length * 100, 2), round((gaps_seq1 + gaps_seq2) / (2 * length) * 100, 2)))

        return results

    @abstractmethod
    def run(self):
        """
        Run the alignment algorithm according to the parameters
        :return:
        """
        pass


class NeedlemanWunsch(Alignment):
    def __init__(self, seq1, seq2, I, E, submat, k):
        """
        Global alignment algorithm
        :param seq1: str, first sequence of amino acids to align
        :param seq2: str, second sequence of amino acids to align
        :param I: float, initial gap value
        :param E: float, extension gap value
        :param submat: SubstitutionMatrix object
        :param k: int, maximum number of solutions to find
        """
        super().__init__(seq1, seq2, I, E, submat, k)
        self.__S = ScoreMatrix(seq1, seq2, I, E, "S")
        self.__V = ScoreMatrix(seq1, seq2, I, E, "V")
        self.__W = ScoreMatrix(seq1, seq2, I, E, "W")

    def calculate_score(self):
        for i in range(1, self.__S.get_num_rows()):
            for j in range(1, self.__S.get_num_cols()):
                self.__V[i, j] = max(
                    [self.__S[i - 1, j] - self._Alignment__I,
                     self.__V[i - 1, j] - self._Alignment__E])

                self.__W[i, j] = max(
                    [self.__S[i, j - 1] - self._Alignment__I,
                     self.__W[i, j - 1] - self._Alignment__E])

                self.__S[i, j] = max([self.__S[i - 1, j - 1] + self._Alignment__subMat[self._Alignment__seq1[i - 1], self._Alignment__seq2[j - 1]],
                                                self.__W[i, j], self.__V[i, j]])

    def run(self):
        self.calculate_score()
        self.backtrack(self.__S.get_num_rows() - 1,
                       self.__S.get_num_cols() - 1, "", "")

    def __deconstruct(self, h_seq, v_seq):
        """
        'pop' character in first place from string
        """
        return (h_seq[1:], v_seq[1:])

    def backtrack(self, i, j, h_seq, v_seq):
        if self._Alignment__maxSols == 0:
            return
        # End of alignment
        if i == 0 and j == 0:
            self._Alignment__maxSols -= 1
            self._Alignment__solutions.append(
                (v_seq, h_seq, float(self.__S[len(self._Alignment__seq1), len(self._Alignment__seq2)]), (0, 0),
                 (len(self._Alignment__seq1), len(self._Alignment__seq2))))
        else:
            # Fill the remaining alignment with gaps
            if i == 0:
                v_seq = "-" + v_seq
                h_seq = self._Alignment__seq2[j - 1] + h_seq
                self.backtrack(i, j - 1, h_seq, v_seq)
            elif j == 0:
                h_seq = "-" + h_seq
                v_seq = self._Alignment__seq1[i - 1] + v_seq
                self.backtrack(i - 1, j, h_seq, v_seq)
            else:
                # Match
                if self.__S[i, j] == self._Alignment__subMat[self._Alignment__seq1[i - 1], self._Alignment__seq2[j - 1]] \
                        + self.__S[i - 1][j - 1]:
                    v_seq = self._Alignment__seq1[i - 1] + v_seq
                    h_seq = self._Alignment__seq2[j - 1] + h_seq
                    self.backtrack(i - 1, j - 1, h_seq, v_seq)
                    h_seq, v_seq = self.__deconstruct(h_seq, v_seq)
                # Horizontal gap
                if self.__S[i, j] == self.__V[i, j]:
                    h_seq = "-" + h_seq
                    v_seq = self._Alignment__seq1[i - 1] + v_seq
                    self.backtrack(i - 1, j, h_seq, v_seq)
                    h_seq, v_seq = self.__deconstruct(h_seq, v_seq)
                # Vertical gap
                if self.__S[i, j] == self.__W[i, j]:
                    v_seq = "-" + v_seq
                    h_seq = self._Alignment__seq2[j - 1] + h_seq
                    self.backtrack(i, j - 1, h_seq, v_seq)
                    h_seq, v_seq = self.__deconstruct(h_seq, v_seq)


class SmithWaterman(Alignment):
    def __init__(self, seq1, seq2, I, E, submat, k):
        """
        Local alignment algorithm
        :param seq1: str, first sequence of amino acids to align
        :param seq2: str, second sequence of amino acids to align
        :param I: float, initial gap value
        :param E: float, extension gap value
        :param submat: SubstitutionMatrix object
        :param k: int, maximum number of solutions to find
        """
        super().__init__(seq1, seq2, I, E, submat, k)
        self.__S = ScoreMatrix(seq1, seq2, I, E, "S", False)
        self.__V = ScoreMatrix(seq1, seq2, I, E, "V", False)
        self.__W = ScoreMatrix(seq1, seq2, I, E, "W", False)

    def __erase_path(self, i, j):
        """
        Called by backtrack for better legibility
        """
        self.__S[i, j] = 0
        self.__V[i, j] = 0
        self.__W[i, j] = 0

    def backtrack(self, i, j, h_seq, v_seq):
        """
        i: Vertical sequence index
        j: Horizontal sequence index
        h_seq, v_seq: Horizontal, Vertical sequences respectively being constructed

        Return: [seq1 aligned, seq2 aligned, Score placeholder, (start alignment seq1, start alignment seq2)]
        """

        # End of alignment. Return
        if self.__S[i, j] == 0:
            self._Alignment__maxSols -= 1
            return [v_seq, None, i+1]
        else:
            # Match
            if self.__S[i, j] == self._Alignment__subMat[self._Alignment__seq1[i - 1], j - 1] \
                    + self.__S[i - 1][j - 1]:
                self.__erase_path(i, j)
                i -= 1
                j -= 1
                v_seq = self._Alignment__seq1[i] + v_seq
            # Horizontal gap
            elif self.__S[i, j] == self.__V[i, j]:
                self.__erase_path(i, j)
                i -= 1
                v_seq = self._Alignment__seq1[i] + v_seq
            # Vertical gap
            elif self.__S[i, j] == self.__W[i, j]:
                self.__erase_path(i, j)
                j -= 1
                v_seq = "-" + v_seq
        return self.backtrack(i, j, h_seq, v_seq)

    def recalculate(self):
        """
        The path taken must be erased (values put to 0 and all the values in the matrix below the last cell of the path
        must be recomputed. The values at 0 must stay at 0.
        :return: None, but the ScoreMatrix S has been modified accordingly
        """
        # Felt lazy
        # Recalculates the entire matrix
        # Might change later, I don't know
        # Also I have a fast enough computer so... ¯\_(ツ)_/¯

        for i in range(1, self.__S.get_num_rows()):
            for j in range(1, self.__S.get_num_cols()):
                if self.__S[i, j] == 0:
                    continue
                self.__V[i, j] = max(
                    [self.__S[i - 1, j] - self._Alignment__I, self.__V[i - 1, j] - self._Alignment__E, 0])

                self.__W[i, j] = max(
                    [self.__S[i, j - 1] - self._Alignment__I, self.__W[i, j - 1] - self._Alignment__E, 0])

                self.__S[i, j] = max([self.__S[i - 1, j - 1] + self._Alignment__subMat[self._Alignment__seq1[i - 1], j-1],
                                      self.__W[i, j], self.__V[i, j], 0])

    def calculate_score(self):
        """
        Negative values are replaced by 0
        """
        for i in range(1, self.__S.get_num_rows()):
            for j in range(1, self.__S.get_num_cols()):
                self.__V[i, j] = max(
                    [self.__S[i - 1, j] - self._Alignment__I, self.__V[i - 1, j] - self._Alignment__E, 0])

                self.__W[i, j] = max(
                    [self.__S[i, j - 1] - self._Alignment__I, self.__W[i, j - 1] - self._Alignment__E, 0])

                self.__S[i, j] = max([self.__S[i - 1, j - 1] + self._Alignment__subMat[self._Alignment__seq1[i - 1], j-1],
                                                self.__W[i, j], self.__V[i, j], 0])

    def run(self):
        self.calculate_score()
        for num in range(self._Alignment__maxSols):
            i, j, m = self.__S.get_max()
            sol = self.backtrack(i, j, "", "")
            sol[1] = round(float(m), 2)  # Place max
            sol[2] = (sol[2], i)
            self._Alignment__solutions.append(tuple(sol))
            self.recalculate()
