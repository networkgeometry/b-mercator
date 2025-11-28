import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

MARKERS = ['o']
# COLORS = ['#2D1E2F', '#F7B32B', '#21A179', '#8E7DBE', '#F72C25']
# COLORS = ["#2D1E2F", "#e76254", "#f7aa58", "#ffe6b7", 
#           "#72bcd5", "#376795","#1e466e"]
COLORS = ["#2D1E2F", "#e76254", "#f7aa58", "#72bcd5", "#376795","#1e466e"]



def load_files(path, filename, ext):
    df = pd.read_csv(f"{path}/{filename}.{ext}", sep="\s+")
    return df


def load_vprops(path, filename):
    return load_files(path, filename, 'inf_vprop')


def load_vprops_bipartite_nodes(path, filename):
    return load_files(path, filename, 'inf_vprop_bipartite_nodes')


def load_vprops_bipartite_features(path, filename):
    return load_files(path, filename, 'inf_vprop_bipartite_features')


def load_vstats(path, filename):
    return load_files(path, filename, 'inf_vstat')


def load_vstats_bipartite_nodes(path, filename):
    return load_files(path, filename, 'inf_vstat_bipartite_nodes')


def load_vstats_bipartite_features(path, filename):
    return load_files(path, filename, 'inf_vstat_bipartite_features')


def load_vstat_obs(path, filename):
    return load_files(path, filename, 'obs_vstat')


def load_vstat_obs_bipartite_nodes(path, filename):
    return load_files(path, filename, 'obs_vstat_bipartite_nodes')


def load_vstat_obs_bipartite_features(path, filename):
    return load_files(path, filename, 'obs_vstat_bipartite_features')


def load_pcon(path, filename):
    df = pd.read_csv(f"{path}/{filename}.inf_pconn", sep="\s+", comment="#")
    df.columns = ["RescaledDist", "InfConnProb", "ThConnProb"]
    return df


def load_bipartite_pcon(path, filename):
    df = pd.read_csv(f"{path}/{filename}.inf_pconn_bipartite",
                     sep="\s+", comment="#")
    df.columns = ["RescaledDist", "InfConnProb", "ThConnProb"]
    return df


def plot_vstats(df_obs, df_vstats, x1_col, y1_col, x2_col, y2_col, y2err_col, xlabel, ylabel, labels, colors=COLORS):
    x1, y1 = df_obs[x1_col], df_obs[y1_col]
    plt.plot(x1, y1, color=colors[0], label="Data",
             marker=MARKERS[0], ms=6, linewidth=1, alpha=0.7)

    for i, df_vstat in enumerate(df_vstats):
        x2, y2, y2err = df_vstat[x2_col], df_vstat[y2_col], df_vstat[y2err_col]
        plt.plot(x2, y2, color=colors[i+1], alpha=0.7, label=labels[i])
        plt.fill_between(x2, y2 - y2err, y2 + y2err,
                         color=colors[i+1], alpha=0.2)

    plt.yscale('log')
    plt.xscale('log')
    #plt.legend(fontsize=18, frameon=False, loc='lower left')
    plt.xlabel(xlabel, fontsize=24)
    plt.ylabel(ylabel, fontsize=24)


def plot_degree_distribution(df_obs, df_vstats, labels, xlabel=r'$k$', ylabel=r'$P_c(k)$', colors=COLORS):
    plot_vstats(df_obs, df_vstats, "#", "DegDist", "#", "DegDistEnsStd", "CDegDistEns",
                xlabel, ylabel, labels, colors)


def plot_clustering_coefficient(df_obs, df_vstats, labels, xlabel=r'$k$', ylabel=r'$\bar{c}(k)$', colors=COLORS):
    plot_vstats(df_obs, df_vstats, "#", "NbTriang", "#", "NbTriangEnsStd", "ClustEns",
                xlabel, ylabel, labels, colors)
    plt.xlim(1.9, 200)
    plt.yscale('linear')


def plot_number_of_triangles(df_vprops, labels, log_scale=True, colors=COLORS, additional_label=''):
    for i, df_vprop in enumerate(df_vprops):
        obs = df_vprop.iloc[:, 10]
        inf = df_vprop.iloc[:, 11]
        plt.scatter(obs, inf, color=colors[i+1],
                    label=labels[i], alpha=0.3, marker='.')

    xx = range(1, int(max(obs)))
    plt.plot(xx, xx, linewidth=2, linestyle='--', color='black')
    if log_scale:
        plt.xscale('log')
        plt.yscale('log')

    plt.ylabel(r'inferred $\# \triangle$' + f'{additional_label}', fontsize=20)
    plt.xlabel(r'original $\# \triangle$' + f'{additional_label}', fontsize=20)


def plot_sum_degree_neighbours(df_vprops, labels, log_scale=True, colors=COLORS, additional_label=''):
    for i, df_vprop in enumerate(df_vprops):
        obs = df_vprop.iloc[:, 7]
        inf = df_vprop.iloc[:, 8]
        plt.scatter(obs, inf, color=colors[i+1],
                    label=labels[i], alpha=0.3, marker='.')

    xx = range(int(min(obs)), int(max(obs)))
    plt.plot(xx, xx, linewidth=2, linestyle='--', color='black')
    if log_scale:
        plt.xscale('log')
        plt.yscale('log')

    plt.ylabel(r'inferred $\bar{k}_{nn}$' + f'{additional_label}', fontsize=20)
    plt.xlabel(r'original $\bar{k}_{nn}$' + f'{additional_label}', fontsize=20)


def plot_connection_probabilities(df_pconns, labels, ylabel=r'$p(\chi)$', ylim_min=1e-5):
    for i, df in enumerate(df_pconns):
        non_zeros = np.where(df.InfConnProb.values != 0)
        x = df.RescaledDist.values[non_zeros]
        y = df.InfConnProb.values[non_zeros]
        plt.plot(x, y, color=COLORS[i+1], label=labels[i])
        plt.plot(df.RescaledDist, df.ThConnProb,
                 color=COLORS[i+1], linestyle='--', label='expected')

    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(ylim_min, 3)
    plt.ylabel(ylabel, fontsize=24)
    plt.xlabel(r"$\chi$", fontsize=24)


def plot_average_degree_neighbour(df_obs, df_vstats, labels, xlabel='$k$',
                                  ylabel=r"$\bar{k}_{nn}(k)$", colors=COLORS):
    plot_vstats(df_obs, df_vstats, "#", "SumDegN", "#", "SumDegNEnsStd", "AvgDegNEns",
                xlabel, ylabel, labels, colors=colors)
    plt.yscale('linear')
