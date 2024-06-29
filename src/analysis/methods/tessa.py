from pathlib import Path
import os

def Tessa (tessaPath, inputTCRPath):
    parent = Path(__file__).parent 
    if tessaPath[0] != "/":
        path = Path(__file__).parent / tessaPath
    else:
        path = tessaPath
    # cmd = f"python3 {path} -tcr {inputTCRPath} -model {parent}/BriseisEncoder/TrainedEncoder.h5 -embedding_vectors {parent}/BriseisEncoder/Atchley_factors.csv -output_TCR test.csv -output_VJ testVJ.csv -output_log test.log -exp {inputGEXPath} -output_tessa {outputPath} -within_sample_networks FALSE"
    cmd = f"python {path}/BriseisEncoder/BriseisEncoder.py -tcr {inputTCRPath} -model {path}/BriseisEncoder/TrainedEncoder.h5 -embeding_vectors {path}/BriseisEncoder/Atchley_factors.csv -output_TCR test.csv -output_VJ testVJ.csv -output_log test.log"
    print(cmd)
    os.system(cmd)

# Tessa("/home/myc0plasmus/Documents/python/TESSA","/home/myc0plasmus/Documents/python/TESSA/example_data/example_TCRmeta.csv")
