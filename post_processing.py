import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import os
from matplotlib import rcParams
import numpy as np

# image_dir = "C:\\Users\\Nikhil Melgiri\\ME559\\Abaqus_CAE\\images"
# data_dir = "C:\\Users\\Nikhil Melgiri\\ME559\\Abaqus_CAE\\data"


def segment_dataframe(df: pd.DataFrame):
    dfs = dict(tuple(df.groupby('category')))
    for key, large_df in dfs.items():
        fraction_dfs = dict(tuple(large_df.groupby('fraction')))
        # order each dataframe in the dictionary by the "number of elements" column
        for fraction, small_df in fraction_dfs.items():
            fraction_dfs[fraction] = small_df.sort_values(by="number of elements")
        dfs[key] = fraction_dfs
    return dfs


def mesh_analysis(df: pd.DataFrame, tolerance: float = 0.01)->dict:
    dfs = segment_dataframe(df)
    
    convergence_dict = {}
    for key, df in dfs.items():
        convergence_dict[key] = {}
        for fraction, small_df in df.items():
            dstress_delements = []
            for i in range(1, len(small_df)):
                dstress = small_df["max stress"].iloc[i] - small_df["max stress"].iloc[i-1]
                average_stress = (small_df["max stress"].iloc[i] + small_df["max stress"].iloc[i-1])/2
                true_average_stress = dstress / average_stress
                delements = small_df["number of elements"].iloc[i] - small_df["number of elements"].iloc[i-1]
                dstress_delements.append(abs(true_average_stress / delements))
            convergence_dict[key][fraction] = all([abs(x) < tolerance for x in dstress_delements])
    return convergence_dict


def plot_convergence(df: pd.DataFrame):
    segment_df = segment_dataframe(df)
    for key, df in segment_df.items():
        for fraction, small_df in df.items():
            x_values = []
            y_values = []
            for i in range(1, len(small_df)):
                dstress = small_df["max stress"].iloc[i] - small_df["max stress"].iloc[i-1]
                average_stress = (small_df["max stress"].iloc[i] + small_df["max stress"].iloc[i-1])/2
                true_average_stress = dstress / average_stress
                delements = small_df["number of elements"].iloc[i] - small_df["number of elements"].iloc[i-1]
                x_values.append(small_df["number of elements"].iloc[i])
                y_values.append(abs(true_average_stress / delements))
            plt.plot(x_values, y_values, label=fraction)
        plt.legend()
        plt.title(key)
        plt.xlabel("Number of Elements")
        plt.ylabel("Stress Change per Element")
        # plt.show()
        key = key.replace(" ", "") + "Graph.png"
        plt.savefig(key)
        plt.clf()


def get_x_y_for_correlations(df: pd.DataFrame, applied_stress: float=100000):
    
    segment_df = segment_dataframe(df)
    x_values_pin = []
    y_values_pin = []
    x_values_no_pin = []
    y_values_no_pin = []
    pin_df = segment_df["Pin"]
    no_pin_df = segment_df["No Pin"]
    for fraction, small_df in pin_df.items():
        x_values_pin.append(fraction)
        y_values_pin.append(small_df["max stress"].sum() / len(small_df) / applied_stress)
    for fraction, small_df in no_pin_df.items():
        x_values_no_pin.append(fraction)
        y_values_no_pin.append(small_df["max stress"].sum() / len(small_df) / applied_stress)
    return x_values_pin, y_values_pin, x_values_no_pin, y_values_no_pin

def create_dataframe(data_dir: str):
    df = pd.DataFrame()
    for file in os.listdir(data_dir):
        data = pd.read_csv(os.path.join(data_dir, file))
        data["image"] = file.replace(".csv", ".png")
        df = pd.concat([df, data], ignore_index=True)
    df.columns = df.columns.str.strip()
    return df

if __name__ == "__main__":
    data_dir = "C:\\Users\\Nikhil Melgiri\\ME559\\Abaqus_CAE\\data"
    df = create_dataframe(data_dir)
    dataframes = mesh_analysis(df)
    # plot_convergence(df)
    x_values_pin, y_values_pin, x_values_no_pin, y_values_no_pin = get_x_y_for_correlations(df)
    pin_polyfit = np.polyfit(x_values_pin, y_values_pin, 2)
    no_pin_polyfit = np.polyfit(x_values_no_pin, y_values_no_pin, 2)
    pin_poly = np.poly1d(pin_polyfit)
    no_pin_poly = np.poly1d(no_pin_polyfit)
    residuals_pin = y_values_pin - pin_poly(x_values_pin)
    residuals_no_pin = y_values_no_pin - no_pin_poly(x_values_no_pin)
    ss_res_pin = np.sum(residuals_pin**2)
    ss_res_no_pin = np.sum(residuals_no_pin**2)
    # plt.plot(x_values_pin, y_values_pin, label="Pin")
    # plt.plot(x_values_no_pin, y_values_no_pin, label="No Pin")
    # plt.plot(x_values_pin, pin_poly(x_values_pin), label="Pin Fit", linestyle="--")
    # plt.plot(x_values_no_pin, no_pin_poly(x_values_no_pin), label="No Pin Fit", linestyle="--")
    # plt.legend()
    # plt.title("K_f vs. Fraction")
    # plt.xlabel("Fraction (Diameter / Width)")
    # plt.ylabel("K_f")
    # plt.savefig("K_f vs. Fraction Fit.png")
    print(pin_poly)
    print(no_pin_poly)
    
    
