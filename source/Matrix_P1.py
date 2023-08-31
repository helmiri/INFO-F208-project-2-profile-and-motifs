import numpy as np


class Matrix:
    def __init__(self, nrows=0, ncols=0, value=None):
        self.__matrix = np.full([nrows, ncols], value, float)

    def __getitem__(self, item):
        return self.__matrix[item]

    def __setitem__(self, index, value):
        """
        The fact that we can use matrix[tuple] or matrix[int] interchangably
        is so convenient
        """
        self.__matrix[index] = value

    def get_num_cols(self):
        return self.__matrix.shape[1]

    def get_num_rows(self):
        return self.__matrix.shape[0]

    def get_max(self):
        """
        :return: tuple of size 3 of the row index, col index and value of the cell with the maximum value in the matrix
        If several occurences of the maximum value are found, the lowest indices are picked.
        """
        result = [0, 0, self.__matrix[0][0]]

        for i in range(self.get_num_rows()):
            for j in range(self.get_num_cols()):
                if self.__matrix[i, j] > result[2]:
                    result[0], result[1], result[2] = i, j, self.__matrix[i][j]

        return tuple(result)

    def set_value(self, i, j, value):
        self.__matrix[i, j] = value


class SubstitutionMatrix(Matrix):
    """
    Format : free
    """

    def __init__(self, file):
        self.parse_file(file)

    def parse_file(self, file):
        bufferMat = list()
        i = 0

        # Save labels for index search
        self.__labels_order = np.empty([0, 0])
        with open(file) as source:
            for line in source:
                buffer = line.split()
                if buffer == []:
                    break
                # Skip commented lines
                if buffer[0] == '#':
                    continue

                # Save label order for index search
                # Assumed labels_cols = labels_rows in a NxN matrix
                if self.__labels_order.size == 0:
                    self.__labels_order = np.array(buffer[:-4])
                    continue

                # Strip col labels + [B, Z, X, *] and construct temp matrix
                bufferMat.append(buffer[1:-4])
                i += 1

        # Strip row labels [B, Z, X, *]
        bufferMat = bufferMat[:-4]

        # Initialize target matrix
        super().__init__(len(bufferMat), len(bufferMat[0]))

        # Convert to int and assign
        for i in range(len(bufferMat)):
            super().__setitem__(i, np.array(bufferMat[i], int))

    def __getitem__(self, item):
        """
        Search by label
        """
        i = np.where(self.__labels_order == item[0])[0][0]
        j = np.where(self.__labels_order == item[1])[0][0]
        return super().__getitem__((i, j))


# Proposition of implementation, but free for you to use it or not, not included in tests
class ScoreMatrix(Matrix):
    """
    Format : list of list of float values
    """

    def __init__(self, seq1, seq2, I, E, matType, is_global=True):
        self.__seq1 = seq1
        self.__seq2 = seq2
        self.__type = matType

        super().__init__(len(seq1) + 1, seq2 + 1, 0)
        if matType == "S":

            if is_global:
                self.set_value(1, 0, -I)
                self.set_value(0, 1, -I)

                for j in range(2, self.get_num_cols()):
                    self.set_value(0, j, super().__getitem__((0, j - 1)) - E)

                for i in range(2, self.get_num_rows()):
                    self.set_value(i, 0, super().__getitem__((i - 1, 0)) - E)

        elif matType == "V":
            for j in range(self.get_num_cols()):
                self.set_value(0, j,  -np.inf)

        elif matType == "W":
            for i in range(self.get_num_rows()):
                self.set_value(i, 0, -np.inf)
