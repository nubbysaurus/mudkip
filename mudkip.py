"""
@file   mudkip.py

Based on input a given field configuration, optimize the input of water by an
irrigation system such that the total available water volume remains greater
than a given level.

@author Jake Lee/nubby  (jlee211@ucsc.edu)
@date   3 Dec 2024
"""
import sys
#import cvx


# Module-level variables.
_DEFAULT_FARM_DATA_PATH = "./data/sample.csv"

# Datatypes.
class Farm(object):
    """
    Each Farm can be represented by an n-sized array of constituent Field nodes,
    F, and a corresponding node-adjacency matrix of size n x n, N. Each element
    of N at row i, column j (== (i,j)) is equal to the following:

    N(i,j) = 1 if F(i) is adjacent to F(j), and
    N(i,j) = 0 otherwise, which includes i = j, and diagonals!

    [For now,] Input data must be structured so that it is fit to a square N
    matrix, i.e. with input Field structure of [1, 2, 3, 4, 5, 6, 7, 8, 9], the
    Farm can be represented by the matrix:
        
        [
            [1, 2, 3],
            [4, 5, 6],
            [7, 8, 9]
        ], and a corresponding N matrix of:

        #    1  2  3  4  5  6  7  8  9
        [
            [0, 1, 0, 1, 0, 0, 0, 0, 0],    # 1
            [1, 0, 1, 0, 1, 0, 0, 0, 0],    # 2
            [0, 1, 0, 0, 0, 1, 0, 0, 0],    # 3
            [1, 0, 0, 0, 1, 0, 1, 0, 0],    # 4
            [0, 1, 0, 0, 0, 1, 0, 1, 0],    # 5
            [0, 0, 1, 0, 1, 0, 0, 0, 1],    # 6
            [0, 0, 0, 1, 0, 0, 0, 1, 0],    # 7
            [0, 0, 0, 0, 1, 0, 1, 0, 1],    # 8
            [0, 0, 0, 0, 0, 1, 0, 1, 0],    # 9
        ]

    @param  data    (list[dict])  Farm data [that we care about for now] in dict
                                  format:

        [
            {
                "Elevation": float,
                "Name": str,
                "WaterVolume": float
            }, ...
        ]

    TODO(nubby): Add initial water content.
    TODO(nubby): Add Field coordinates and size, then describe adjacency with
                 more nuance.
    """
    def __init__(self, data: list[dict]):
        self.F = [] 
        self.N = []
        self.size = len(size.keys())
        self._build_farm(data)

    def _build_farm(self, data: list[dict]):
        """_build_farm(size)
        Initialize a new array of Fields of a desired size and build their
        adjacency matrix.
        """
        if len(data) == 9:
            #    1  2  3  4  5  6  7  8  9
            self.N = [
                [0, 1, 0, 1, 0, 0, 0, 0, 0],    # 1
                [1, 0, 1, 0, 1, 0, 0, 0, 0],    # 2
                [0, 1, 0, 0, 0, 1, 0, 0, 0],    # 3
                [1, 0, 0, 0, 1, 0, 1, 0, 0],    # 4
                [0, 1, 0, 0, 0, 1, 0, 1, 0],    # 5
                [0, 0, 1, 0, 1, 0, 0, 0, 1],    # 6
                [0, 0, 0, 1, 0, 0, 0, 1, 0],    # 7
                [0, 0, 0, 0, 1, 0, 1, 0, 1],    # 8
                [0, 0, 0, 0, 0, 1, 0, 1, 0],    # 9
            ]
        else:
            raise ValueError(
                "WIP: Currently only supports Farms of exactly 9 Field nodes!"
            )
        [self.F.append(Field(field_data)) for field_data in data]


class Field(object):
    """Field()
    Container for information and dynamics relating to each Field node.
    NOTE: These aim to include, if not add to, all of the properties of Field
          nodes used in APSIM.

    A Field imports the :

        {
            "Elevation": float,
            "Name": str,
            "WaterVolume": float
        }
    """
    def __init__(self, data: dict):
        pass


# Helpers.
def _load_csv(fpath: str) -> list[dict]:
    data = []
    print(f"Reading from {fpath}...")
    with open(fpath, "r+") as csvs:
        reader = csv.DictReader(csvs)
        for row in reader:
            data.append(row)
    print(f"    DONE")
    if (not data):
        raise OSError(
            f"WARNING!! {fpath} is an empty file!"
        )
    return data


# Main functionality.
def run_mudkip(farm: Farm):
    pass

def clean_up(e: Exception):
    print(f"{e}; good bye.")
    sys.exit(0)

if __name__ == "__main__":
    try:
        farm_data = _load_csv(_DEFAULT_FARM_DATA_PATH)
    except OSError as e:
        clean_up(e)

    try:
        farm = Farm(farm_data)
    except ValueError as e:
        clean_up(e)
    watering_schedule = run_mudkip(farm)
