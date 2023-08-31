from math import log2, sqrt
from source.Matrix_P1 import *


class PSSM(Matrix):
    global PROB_AA
    PROB_AA = {'A': 0.0828, 'Q': 0.0394, 'L': 0.0967, 'S': 0.065,
               'R': 0.0553, 'E': 0.0676, 'K': 0.0585, 'T': 0.0532, 'N': 0.0405,
               'G': 0.0709, 'M': 0.0243, 'W': 0.0107, 'D': 0.0545, 'H': 0.0227,
               'F': 0.0386, 'Y': 0.0291, 'C': 0.0136, 'I': 0.0599, 'P': 0.0468,
               'V': 0.0687}

    def __init__(self, multiple_alignment):
        self.alignments = []
        self.row_labels = list(PROB_AA.keys())

        # Extract ailgnments from file
        with open(multiple_alignment, "r") as alignments:
            line = alignments.readline()
            count = 0
            while line:
                if count % 2 != 0:
                    self.alignments.append(line.strip("\n"))

                count += 1
                line = alignments.readline()

        # self.occurrences = list containing a dict of occurrences where self.occurrences[i]
        # is the number of occurrences of each amino acid at position i of the sequences
        self.occurences = [dict.fromkeys(PROB_AA, 0)
                           for i in range(len(self.alignments[0]))]

        super().__init__(len(PROB_AA), len(self.alignments[0]), 0)
        self.count_occurences()
        self.PWM()
        self.calculate_probabilities()

    def count_occurences(self):
        """
        Count the occurences of amino acids at each position
        """
        for i in range(len(self.alignments[0])):
            for alignment in self.alignments:
                key = alignment[i]
                if key != '-':
                    self.occurences[i][key] += 1

    def PWM(self):
        """
        Position-weight matrix is a matrix using the occurences of amino acids at specific position and
        normalize the counts according to the number of sequences to obtain the frequency of the amino acid at
        a position. We use here also pseudocount (beta) to avoid a null frequency at a position. The gaps are ignored.
        Corresponds to the q(u,a) formula
        """

        for i in range(self.get_num_rows()):
            label = self.row_labels[i]
            for j in range(self.get_num_cols()):
                self[label, j] = self.q(j, label)

    def f(self, u, b):
        return self.occurences[u][b] / len(self.alignments)

    def q(self, u, a):
        ALPHA = len(self.alignments) - 1
        BETA = sqrt(len(self.alignments))
        return ((ALPHA * self.f(u, a)) + (BETA * PROB_AA[a])) / (ALPHA + BETA)

    def calculate_probabilities(self):
        """
        To do the transition from PWM to PSSM, you have to transform the normalized frequency to a probability, with
        the help of the log-odds. This is the logarithm between the probability of an event and
        of this event occurring randomly. The random event here is the probability of finding the amino acid anywhere,
        the probability of the event is our observation of having this amino acid at this specific position.
        Corresponds to m(u,a) formula
        """
        for i in range(self.get_num_rows()):
            label = self.row_labels[i]
            for j in range(self.get_num_cols()):

                self[label, j] = round(
                    log2(self[label, j] / PROB_AA[label]), 2)

    def __setitem__(self, index, value):
        """
        index = (letter, column)
        """

        i = self.row_labels.index(index[0])
        super().__setitem__((i, index[1]), value)

    def __getitem__(self, item):
        i = self.row_labels.index(item[0])
        return super().__getitem__((i, item[1]))
