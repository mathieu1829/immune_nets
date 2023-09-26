# generate all types networks and save it in specified directory
import argparse
from src.creation.algorithms.common_methods import *
from src.creation.algorithms.simple_distance import *
from src.creation.algorithms.simple_vector_distance_v2 import *
from src.creation.algorithms.simple_vector_distance import *

from src.creation.inputStrategies.db_strategy import *
from src.creation.outputStrategies.db_strategy import db_strategy as db_output_strategy
from src.creation.inputStrategies.csv_strategy import *
from src.creation.enums.matrices import *


parser = argparse.ArgumentParser()
parser.add_argument('-i','--input', help='Provide path to file with clonotypes')
parser.add_argument('-s','--input_strategy', help='Provide strategy to process input data (csv, db)')
parser.add_argument('-g','--output_strategy', help='Provide strategy to process output data ')
parser.add_argument('-o','--output', help='Provide path to location where the results shall be saved.')

def main():
    args = parser.parse_args()
    path = args.input #"..\\..\\tests\\test_data\\test_clonotypes.csv"
    input_strategy = args.input_strategy
    match input_strategy:
        case "db":
            df = db_strategy("data")
        case "csv":
            df = csv_strategy(path)
    if df is None:
        print("ERROR: invalid strategy")

    SimpleDistance(db_output_strategy).createGraph(clonotypes=df,matrix=Matrices.BLOSUM62)
    SimpleDistance(db_output_strategy).createGraph(clonotypes=df,matrix=Matrices.PAM250)
    simple_vector_distance(db_output_strategy).createGraph(clonotypes=df,matrix=Matrices.BLOSUM62)
    simple_vector_distance(db_output_strategy).createGraph(clonotypes=df,matrix=Matrices.PAM250)
    simple_vector_distance_v2(db_output_strategy).createGraph(clonotypes=df,matrix=Matrices.BLOSUM62)
    simple_vector_distance_v2(db_output_strategy).createGraph(clonotypes=df,matrix=Matrices.PAM250)



if __name__ == "__main__":

    main()
