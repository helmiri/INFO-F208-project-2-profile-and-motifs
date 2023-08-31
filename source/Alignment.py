from source.Matrix import *
from source.Alignment_P1 import *
from abc import ABC, abstractmethod

# Just modified SmithWaterman a little bit to exclude the second sequence to compare
# I know, very bad OOP but I'm lazy and it works so yeah...


class SmithWatermanProfile(SmithWaterman):
    def __init__(self, seq, I, E, PSSM, k):
        """
        Local alignment algorithm
        :param seq: str, sequence of amino acids to find matches for the profile matrix
        :param I: float, initial gap value
        :param E: float, extension gap value
        :param PSSM: PSSM object
        :param k: int, maximum number of solutions to find
        """
        super().__init__(seq, PSSM.get_num_cols(), I, E,
                         PSSM, k)
