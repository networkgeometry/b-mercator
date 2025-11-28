from yaml import load, Loader
import os
import argparse
import textwrap
import itertools
from tqdm import tqdm

def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(
            """    
    Example:

    > python test/generate_bipartiteSD_networks.py \
        -o [output_folder]
        -c [config.yaml]
    """
        ),
    )
    parser.add_argument(
        "-o", "--output_folder", type=str, required=True, help="Path to output folder"
    )
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        required=True,
        help="Path to config file with network parameters",
    )
    args = parser.parse_args()
    return args


def construct_folder(
    dim, Beta_s, gamma_s, Ns_obs, kmean_s, gamma_n, kmean_n, gamma_f, N_f, Beta_bi, correlation, i
):
    return f"dim_{dim}_B_s_{Beta_s}_g_s_{gamma_s}_Ns_obs_{Ns_obs}_k_s_{kmean_s}_g_n_{gamma_n}_k_n_{kmean_n}_g_f_{gamma_f}_N_f_{N_f}_B_bi_{Beta_bi}_c_{correlation}_i_{i}/"


if __name__ == "__main__":
    args = parse_args()
    config = load(open(args.config, "r"), Loader=Loader)

    for parameter in config:
        if "Beta_s" in parameter.keys():
            Beta_s = parameter["Beta_s"]
        elif "gamma_s" in parameter.keys():
            gamma_s = parameter["gamma_s"]
        elif "Ns_obs" in parameter.keys():
            Ns_obs = parameter["Ns_obs"]
        elif "kmean_s" in parameter.keys():
            kmean_s = parameter["kmean_s"]
        elif "gamma_n" in parameter.keys():
            gamma_n = parameter["gamma_n"]
        elif "kmean_n" in parameter.keys():
            kmean_n = parameter["kmean_n"]
        elif "gamma_f" in parameter.keys():
            gamma_f = parameter["gamma_f"]
        elif "N_f" in parameter.keys():
            N_f = parameter["N_f"]
        elif "Beta_b" in parameter.keys():
            Beta_b = parameter["Beta_b"]
        elif "ntimes" in parameter.keys():
            ntimes = parameter["ntimes"]
        elif "dimension" in parameter.keys():
            dimension = parameter["dimension"]
        elif "correlation" in parameter.keys():
            correlation = parameter["correlation"]

    os.system("g++ -O3 --std=c++17 -o gen_SD src/generatingSDBipartiteSD_unix.cpp")

    all_params = [
        Beta_s,
        gamma_s,
        Ns_obs,
        kmean_s,
        gamma_n,
        kmean_n,
        gamma_f,
        N_f,
        Beta_b,
        dimension,
        correlation,
    ]

    for i in tqdm(range(ntimes)):
        for (
            Beta_s_val,
            gamma_s_val,
            Ns_obs_val,
            kmean_s_val,
            gamma_n_val,
            kmean_n_val,
            gamma_f_val,
            N_f_val,
            Beta_b_val,
            dimension_val,
            correlation_val
        ) in itertools.product(*all_params):

            output_folder = construct_folder(
                dimension_val,
                Beta_s_val,
                gamma_s_val,
                Ns_obs_val,
                kmean_s_val,
                gamma_n_val,
                kmean_n_val,
                gamma_f_val,
                N_f_val,
                Beta_b_val,
                correlation_val,
                i,
            )

            folder = f'{args.output_folder}/{output_folder}'
            os.makedirs(folder)
            command = f"./gen_SD -d {dimension_val} -b {Beta_s_val*dimension_val} -g {gamma_s_val} -n {Ns_obs_val} -k {kmean_s_val} -e {gamma_n_val} -q {kmean_n_val} -t {gamma_f_val} -f {N_f_val} -p {Beta_b_val*dimension_val} -c {correlation_val} -r {folder} -v"
            os.system(command)
            