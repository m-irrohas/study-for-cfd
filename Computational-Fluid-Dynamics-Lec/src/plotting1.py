import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd

def visualize(filename, all_plot=True, save_fig=True):
    """結果をプロットするだけ
    Arg:
        filename(str)
    Return:
        None
    """
    rows = []
    with open(filename, mode="r") as f:
        reader = csv.reader(f, delimiter=' ', skipinitialspace=True)
        for row in reader:
            row = row[:-1]
            row = list(map(float, row))
            rows.append(row)
    
    print(rows)
    X = rows[0]
    U = rows[1:]
    if all_plot:
        for u in U:
            plt.plot(X,u)
    else:
        pass #とりあえず
    if save_fig:
        plt.savefig("./output/ftcs_v1.png")
    plt.show()
    
    return None

visualize("test.csv")