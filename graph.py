import numpy
import pandas as pd
import matplotlib.pyplot as plt

def build_graph(file, rows_x, rows_y, titol, x_label, y_label, x_selection_rows, x_selection_value, y_selection_rows, y_selection_value):
    df = pd.read_csv(file, header = None, names = ["Points", "Num Clusters", "Num Thread", "Total RunTime (36 execution)", "Mean runtime", "Standard Deviation", "Best case", "Worst case"])
    i = 0
    df_speedup = []
    while i < len(x_selection_value):
        df_speedup.append(df.loc[(df[x_selection_rows] == x_selection_value[i]) & (df[y_selection_rows] == y_selection_value[i])])
        i += 1
    
    for dfi in df_speedup:
        dfi["Speedup"] = dfi['Total RunTime (36 execution)'].iloc[0] / dfi['Total RunTime (36 execution)']
        #print(dfi)
    i = 0
    for i, dfi in enumerate(df_speedup):
        plt.plot(dfi[rows_x], dfi[rows_y], label="Num. Clusters: " + str(y_selection_value[i]))
        i += 1
    plt.xlabel(rows_x)
    plt.ylabel(rows_y)
    plt.title(titol)
    plt.legend()
    plt.grid(True)
    #plt.show()
    plt.savefig("images/SpeedUpGraph_Points_" + str(x_selection_value[0]) + "_V2.png")
    
def build_graph_comparison(file1, file2, rows_x, rows_y, titol, x_selection_rows, x_selection_value, y_selection_rows, y_selection_value):
    df1 = pd.read_csv(file1, header = None, names = ["Points", "Num Clusters", "Num Thread", "Total RunTime (36 execution)", "Mean runtime", "Standard Deviation", "Best case", "Worst case"])
    df2 = pd.read_csv(file2, header = None, names = ["Points", "Num Clusters", "Num Thread", "Total RunTime (36 execution)", "Mean runtime", "Standard Deviation", "Best case", "Worst case"])

    df_selection1 = df1.loc[(df1[x_selection_rows] == x_selection_value) & (df1[y_selection_rows] == y_selection_value)]
    df_selection2 = df2.loc[(df2[x_selection_rows] == x_selection_value) & (df2[y_selection_rows] == y_selection_value)]

    df_selection1 = df_selection1.add_suffix('_v1')
    df_selection2 = df_selection2.add_suffix('_v2')
    
    df = pd.concat([df_selection1, df_selection2], axis = 1)
    df["Speedup"] = df['Total RunTime (36 execution)_v1'] / df['Total RunTime (36 execution)_v2']

    plt.plot(df[rows_x], df[rows_y], label = "n Points: " + str(x_selection_value) + "n Clusters: " + str(y_selection_value))

    plt.xlabel("Num. Threads")
    plt.ylabel(rows_y)
    plt.title(titol)
    plt.legend()
    plt.grid(True)
    #plt.show()
    plt.savefig("images/SpeedUpGraph_Points_" + str(x_selection_value) + "_Clusters_" + str(y_selection_value)+ "_V2.png")
    
def build_graph_comparison_cumulative(file1, file2, rows_x, rows_y, titol, y_selection_rows, y_selection_value):
    df1 = pd.read_csv(file1, header = None, names = ["Points", "Num Clusters", "Num Thread", "Total RunTime (36 execution)", "Mean runtime", "Standard Deviation", "Best case", "Worst case"])
    df2 = pd.read_csv(file2, header = None, names = ["Points", "Num Clusters", "Num Thread", "Total RunTime (36 execution)", "Mean runtime", "Standard Deviation", "Best case", "Worst case"])

    df_selection1 = df1.loc[(df1[y_selection_rows] == y_selection_value)]
    df_selection2 = df2.loc[(df2[y_selection_rows] == y_selection_value)]
    
    normalized_factor = (df_selection1["Points"] ** 2)
    
    df_selection1['Total RunTime (36 execution)'] /= normalized_factor
    df_selection2['Total RunTime (36 execution)'] /= normalized_factor

    df_selection1 = df_selection1.add_suffix('_v1')
    df_selection2 = df_selection2.add_suffix('_v2')
    
    df = pd.concat([df_selection1, df_selection2], axis = 1)
    df["Speedup"] = df['Total RunTime (36 execution)_v1'] / df['Total RunTime (36 execution)_v2']
    df_speedup_mean = df.groupby('Num Thread_v1')['Speedup'].mean().reset_index()
    print(df)

    plt.plot(df_speedup_mean[rows_x], df_speedup_mean[rows_y], label = "n Clusters: " + str(y_selection_value))
    plt.xlabel("Num. Threads")
    plt.ylabel(rows_y)
    plt.title(titol)
    plt.legend()
    plt.grid(True)
    #plt.show()
    plt.savefig("images/SpeedUpGraph_Points_" + str(y_selection_value) + "_V2.png")
        
#build_graph("results/Results_Medoids_Windows_Version2.csv", "Num Thread", "Speedup", "Speedup versus Number of Threads (Num. Points: 500) Version 2", "num Points", "num clusters", "Points", [500,500,500,500,500], "Num Clusters", [3,5,8,10,16])        
build_graph_comparison("results/Results_Medoids_Windows_Version1.csv","results/Results_Medoids_Windows_Version2.csv", "Num Thread_v1", "Speedup", "Speedup Comparison between Program Version 1 and Version 2", "Points", 5000, "Num Clusters", 16)        
#build_graph_comparison_cumulative("results/Results_Medoids_Windows_Version1.csv","results/Results_Medoids_Windows_Version2.csv", "Num Thread_v1", "Speedup", "Speedup Comparison between Program Version 1 and Version 2","Num Clusters", 3)        

#file_path = "results/Results_Medoids_Windows_Version1.csv"

#df = pd.read_csv(file_path, header = None, names = ["Points", "Num Clusters", "Num Thread", "Total RunTime (36 execution)", "Mean runtime", "Standard Deviation", "Best case", "Worst case"])
#df_speedUp_10000 = df.loc[(df["Points"] == 10000) & (df["Num Clusters"] == 10)]
#df_speedUp_5000 = df.loc[(df["Points"] == 5000) & (df["Num Clusters"] == 5)]
#df_speedUp_2500 = df.loc[(df["Points"] == 2500) & (df["Num Clusters"] == 5)]
#df_speedUp_1000 = df.loc[(df["Points"] == 1000) & (df["Num Clusters"] == 5)]
#df_speedUp_500 = df.loc[(df["Points"] == 500) & (df["Num Clusters"] == 5)]
#df_speedUp_3 = df.loc[(df["Points"] == 5000) & (df["Num Clusters"] == 3)]
#df_speedUp_5 = df.loc[(df["Points"] == 5000) & (df["Num Clusters"] == 5)]
#df_speedUp_10 = df.loc[(df["Points"] == 5000) & (df["Num Clusters"] == 10)]
#df_speedUp_10000['Speedup'] = df_speedUp_10000['Total RunTime (36 execution)'].iloc[0] / df_speedUp_10000['Total RunTime (36 execution)']
#df_speedUp_5000['Speedup'] = df_speedUp_5000['Total RunTime (36 execution)'].iloc[0] / df_speedUp_5000['Total RunTime (36 execution)']
#df_speedUp_2500['Speedup'] = df_speedUp_2500['Total RunTime (36 execution)'].iloc[0] / df_speedUp_2500['Total RunTime (36 execution)']
#df_speedUp_1000['Speedup'] = df_speedUp_1000['Total RunTime (36 execution)'].iloc[0] / df_speedUp_1000['Total RunTime (36 execution)']
#df_speedUp_500['Speedup'] = df_speedUp_500['Total RunTime (36 execution)'].iloc[0] / df_speedUp_500['Total RunTime (36 execution)']
#df_speedUp_3['Speedup'] = df_speedUp_3['Total RunTime (36 execution)'].iloc[0] / df_speedUp_3['Total RunTime (36 execution)']
#df_speedUp_5['Speedup'] = df_speedUp_5['Total RunTime (36 execution)'].iloc[0] / df_speedUp_5['Total RunTime (36 execution)']
#df_speedUp_10['Speedup'] = df_speedUp_10['Total RunTime (36 execution)'].iloc[0] / df_speedUp_10['Total RunTime (36 execution)']
#plt.plot(df_speedUp_10000["Num Thread"], df_speedUp_10000["Speedup"], label="Curva10000")
#plt.plot(df_speedUp_3["Num Thread"], df_speedUp_3["Speedup"], label="Curva3")
#plt.plot(df_speedUp_5["Num Thread"], df_speedUp_5["Speedup"], label="Curva5")
#plt.plot(df_speedUp_10["Num Thread"], df_speedUp_10["Speedup"], label="Curva10")
#plt.plot(df_speedUp_500["Num Thread"], df_speedUp_500["Speedup"], label="Curva500")
#plt.xlabel("Num. Thread")
#plt.ylabel("SpeedUp")
#plt.title("k-medoids_V1 speedup")
#plt.legend()
#plt.grid(True)
#plt.show()

