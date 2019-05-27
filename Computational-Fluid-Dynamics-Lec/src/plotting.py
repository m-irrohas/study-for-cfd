import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
import pathlib

def visualize(filename, all_plot=True, save_fig=True):
    """結果をプロットするだけ
    Arg:
        filename(str) <= .csv
    Return:
        filename(のstem).png <= (save_fig==True)
        None <= (save_fig==False)
    """
    rows = []
    with open(filename, mode="r") as f:
        reader = csv.reader(f, delimiter=' ', skipinitialspace=True)
        for row in reader:
            row = row[:-1]
            row = list(map(float, row))
            rows.append(row)
    
    plt.clf()
    X = rows[0]
    U = rows[1:]
    if all_plot:
        for u in U:
            plt.plot(X,u)
    else:
        pass #とりあえず
        
    file_path = pathlib.Path(filename)
    file_stem = file_path.stem
    if save_fig:
        plt.savefig(file_stem+".png")
    plt.show()
    return None

ftcs1 = "../output/ftcs1.csv"
ftcs2 = "../output/ftcs2.csv"
visualize(ftcs1)
visualize(ftcs2)