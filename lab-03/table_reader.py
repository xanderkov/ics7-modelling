from math import log

class TableReader:
    def __init__(self, filename=None, arr=None):
        self.filename = filename
        self.dot_arr = arr

    def print_newt_table(self, table):
        n = len(table)
        for i in range(n):
            for j in range(n - i + 1):
                print("{:<15.8f}".format(table[i][j]), end=' ')
            print()

    def load_dots_from_file(self, ind=1):
        try:
            f = open(self.filename)
        except:
            print("File doesn't exist")
            return []
        dots = []
        line = f.readline()
        if ind == 1:
            while line:
                try:
                    x, y = map(float, line.split())
                    dots.append([log(x), log(y), x, y])
                    # dots.append([x, y, x, y])
                except:
                    print("File wasn't properly read")
                    break
                line = f.readline()
        f.close()
        self.dot_arr = dots

    def make_newton_table(self):
        n = len(self.dot_arr)
        m = len(self.dot_arr) + 1
        table = [0] * n
        for i in range(n):
            table[i] = [0] * m
        for i in range(n):
            table[i][0] = self.dot_arr[i][0]
            table[i][1] = self.dot_arr[i][1]
        for j in range(2, n + 1):
            for i in range(n - j + 1):
                table[i][j] = (table[i][j - 1] - table[i + 1][j - 1]) / (table[i][0] - table[i + (j - 1)][0])
        return table

    def make_table_from_file(self):
        self.load_dots_from_file()
        table = self.make_newton_table()
        return table

