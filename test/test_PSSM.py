import unittest

from source.Matrix import *

DIR = "ressources" # must be changed when testing locally (annoying relative path in Github workflow)


class TestPSSM(unittest.TestCase):

    def test_score(self):
        pssm = PSSM(f"{DIR}/msaresults-MUSCLE.fasta")
        self.assertEqual(-0.81, pssm['L', 0])
        self.assertEqual(0.78, pssm['I', 17])
        self.assertEqual(-0.51, pssm['S', 56])
        self.assertEqual(6.15, pssm['W', 54])
        self.assertEqual(6.43, pssm['W', 5])


if __name__ == '__main__':
    unittest.main()
