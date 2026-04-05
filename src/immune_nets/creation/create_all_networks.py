# generate all types networks and save it in specified directory
import argparse
from immune_nets.creation.algorithms.common_methods import *
from immune_nets.creation.algorithms.simpleDistance import *
from immune_nets.creation.distance.alignment import sequenceAligner
from immune_nets.creation.algorithms.simpleVectorDistanceV2 import *
from immune_nets.creation.algorithms.simpleVectorDistance import *
from pathlib import Path
import importlib
import os
import inspect 
from immune_nets.creation.distance.alignment import sequenceAligner
from immune_nets.creation.distance.hamming import hammingDistance
from immune_nets.creation.distance.negativeHamming import negativeHammingDistance
from immune_nets.entities import ImmuneRepertoire

from immune_nets.creation.io_strategies.db_strategy import db_strategy 
from immune_nets.creation.io_strategies.test_csv_strategy import test_csv_strategy
from immune_nets.creation.enums.matrices import *
import uuid

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input', help='Provide path to file with clonotypes')
parser.add_argument('-s','--input_strategy', help='Provide strategy to process input data (csv, db)')
parser.add_argument('-g','--output_strategy', help='Provide strategy to process output data ')
parser.add_argument('-o','--output', help='Provide path to location where the results shall be saved.')

def getAllAlgorithms():
    path = Path(__file__).parent / "algorithms" 
    argList = [ algo.rstrip(".py") for algo in os.listdir(path) if algo != "algorithm.py" and algo != "__pycache__" and algo != "common_methods.py" ]
    algoList = []
    for algo in argList:
        some_algorithm = importlib.import_module(f'immune_nets.creation.algorithms.{algo}', package=None)
        algoList.append(dict(inspect.getmembers(some_algorithm,predicate=inspect.isfunction))[algo])
        # print(dict(inspect.getmembers(some_algorithm,predicate=inspect.isfunction)))
    return algoList
    

       
    


def main():
    args = parser.parse_args()
    # path = args.input #"..\\..\\tests\\test_data\\test_clonotypes.csv"
    path = "tests/test_data/test_clonotypes.csv"
    # input_strategy = args.input_strategy
    input_strategy = "csv"
    match input_strategy:
        case "db":
            df = db_strategy().input()
        case "csv":
            df = test_csv_strategy().input(path)
            df.clones.name = "testData"
    if df is None:
        print("ERROR: invalid strategy")

    for algo in getAllAlgorithms():
        for dist in [hammingDistance(), negativeHammingDistance(), sequenceAligner('BLOSUM62'), sequenceAligner('PAM250')]:
            algo(repertoire=df, distance=dist, strategy = db_strategy().output)

    # SimpleDistance(db_strategy().output).createGraph(clonotypes=df,matrix=Matrices.BLOSUM62)
    # SimpleDistance(db_strategy().output).createGraph(clonotypes=df,matrix=Matrices.PAM250)
    # simpleVectorDistance(db_strategy().output).createGraph(clonotypes=df,matrix=Matrices.BLOSUM62)
    # simpleVectorDistance(db_strategy().output).createGraph(clonotypes=df,matrix=Matrices.PAM250)
    # simpleVectorDistanceV2(db_strategy().output).createGraph(clonotypes=df,matrix=Matrices.BLOSUM62)
    # simpleVectorDistanceV2(db_strategy().output).createGraph(clonotypes=df,matrix=Matrices.PAM250)

if __name__ == "__main__":
    main()
    # getAllAlgorithms()
