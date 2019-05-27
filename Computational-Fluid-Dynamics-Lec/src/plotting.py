import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
import pathlib
import glob

def visualize(filename, all_plot=False, save_fig=True):
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
        u_start = U[0]
        u_end = U[-1]
        plt.plot(X,u_start, label="initial")
        plt.plot(X, u_end, label="end")

    file_path = pathlib.Path(filename)
    file_stem = file_path.stem
    plt.legend()
    if save_fig:
        plt.savefig(file_stem+".png")
    plt.show()
    return None

now_path = pathlib.Path.cwd()
ftcs1_path = now_path/"output"/"ftcs1.csv"
visualize(ftcs1_path)
ftcs2_path = now_path/"output"/"ftcs2.csv"
visualize(ftcs2_path)