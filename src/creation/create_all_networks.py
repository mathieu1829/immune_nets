# generate all types networks and save it in specified directory
import argparse
from algoritms.common_methods import *
from algoritms.simple_distance import *

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input', help='Provide path to file with clonotypes')
parser.add_argument('-o','--output', help='Provide path to location where the results shall be saved.')

def main():
    args = parser.parse_args()
    path = args.input #"..\\..\\tests\\test_data\\test_clonotypes.csv"
    df = pd.read_csv(path)
    df['tcra_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRA"))
    df['tcrb_aa'] = df['cdr3s_aa'].apply(lambda x: split_tcr_column(x, subunit="TRB"))

    df_net = create_BLOSUM62_network(df)
    df_net.to_csv(args.output + "results.csv")



if __name__ == "__main__":

    main()
