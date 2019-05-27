import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
import pathlib
import glob

def visualize(filename, correct_filename, all_plot=False, save_fig=True):
    """結果をプロットするだけ
    Arg:
        filename(str) <= .csv
    Return:
        filename(のstem).png <= (save_fig==True)
        None <= (save_fig==False)
    """
    rows = []
    corrects = []
    with open(filename, mode="r") as f:
        reader = csv.reader(f, delimiter=' ', skipinitialspace=True)
        for row in reader:
            row = row[:-1]
            row = list(map(float, row))
            rows.append(row)

    with open(correct_filename, mode="r") as f:
        reader = csv.reader(f, delimiter=' ', skipinitialspace=True)
        for row in reader:
            row = row[:-1]
            row = list(map(float, row))
            corrects.append(row)
    plt.clf()
    X = rows[0]
    U = rows[1:]
    u_right = corrects[1]
    if all_plot:
        for u in U:
            plt.plot(X,u)
    else:
        u_start = U[0]
        u_end = U[-1]
        plt.plot(X,u_start, label="initial")
        plt.plot(X, u_end, label="end")
        plt.plot(X, u_right, label="correct")

    file_path = pathlib.Path(filename)
    file_stem = file_path.stem
    plt.legend()
    if save_fig:
        plt.savefig("output/"+file_stem+".png")
    plt.show()
    return None

now_path = pathlib.Path.cwd()
correct1_path = now_path/"output"/"correct1.csv"
correct2_path = now_path/"output"/"correct2.csv"
ftcs1_path = now_path/"output"/"ftcs1.csv"
ftcs2_path = now_path/"output"/"ftcs2.csv"
LW1_path = now_path/"output"/"LW1.csv"
LW2_path = now_path/"output"/"LW2.csv"
visualize(ftcs1_path, correct1_path)
visualize(ftcs2_path, correct2_path)
visualize(LW1_path, correct1_path)
visualize(LW2_path, correct2_path)