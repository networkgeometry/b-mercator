from numba import jit
import numpy as np
import pandas as pd
from tqdm import tqdm


# @jit(nopython=True)
def euclidean_to_hyperspherical_coordinates(vec):
    # From: https://en.wikipedia.org/wiki/N-sphere
    # vec -- coordinates of node with size D+1
    r = np.linalg.norm(vec)
    angles = [r]
    for i in range(len(vec) - 2):
        bottom = 0
        for j in range(i, len(vec)):
            bottom += vec[j] * vec[j]
        bottom = np.sqrt(bottom)
        angles.append(np.arccos(vec[i] / bottom))

    denominator = np.sqrt(vec[-1] * vec[-1] + vec[-2] * vec[-2])
    if denominator < 1e-15:
        theta = 0
    else:
        theta = np.arccos(vec[-2] / denominator)
    
    if vec[-1] < 0:
        theta = 2 * np.pi - theta

    angles.append(theta)
    return angles


# @jit(nopython=True)
def hyperspherical_to_euclidean_coordinates(v):
    positions = []
    angles = v[1:]
    r = v[0]
    for i in range(len(angles)):
        val = np.cos(angles[i])
        for j in range(i):
            val *= np.sin(angles[j])
        positions.append(r * val)

        if i == (len(angles) - 1):
            val = np.sin(angles[i])
            for j in range(i):
                val *= np.sin(angles[j])
            positions.append(r * val)
    return positions


@jit(nopython=True)
def compute_angular_distances(x, y):
    angular_distances = []
    for v, u in zip(x, y):
        angular_distances.append(
            np.arccos(np.dot(v, u) / (np.linalg.norm(v) * np.linalg.norm(u))))
    return angular_distances


def get_rotation_matrix_SD(u, v=np.matrix([1, 0, 0, 0], dtype=float)):
    # From https://math.stackexchange.com/a/4167838
    u = np.matrix(u)
    v = np.matrix(v)    
    uv = u + v
    
    R = np.eye(4) - uv.T / (1 + u @ v.T) * uv + 2 * np.multiply(u, v.T)
    return R


def rotation_matrix_XY(theta):
    m = np.matrix([[1, 0, 0, 0],
                   [0, 1, 0, 0],
                   [0, 0, np.cos(theta), np.sin(theta)],
                   [0, 0, -np.sin(theta), np.cos(theta)]])
    return m


def rotation_matrix_XW(theta):
    m = np.matrix([[1, 0, 0, 0],
                   [0, np.cos(theta), -np.sin(theta), 0],
                   [0, np.sin(theta), np.cos(theta), 0],
                   [0, 0, 0, 1]])
    return m


def apply_pipeline_matrix_S3(pos1, pos2, theta_num=50):
    mean_distances = []
    pos1_rotation_matrices = []
    pos2_rotation_matrices = []
    thetas1 = []
    thetas2 = []

    for i in (range(len(pos1))):
        node_i = pos1[i]
        node_i_prime = pos2[i]

        m_i = get_rotation_matrix_SD(node_i)
        m_j = get_rotation_matrix_SD(node_i_prime)

        pos1_axis = np.array([(m_i @ v).A1 for v in pos1])
        pos2_axis = np.array([(m_j @ v).A1 for v in pos2])

        for theta1 in np.linspace(0, 2*np.pi, num=theta_num):
            pos2_XY = np.matmul(rotation_matrix_XY(theta1), pos2_axis.transpose()).T

            for theta2 in np.linspace(0, 2*np.pi, num=theta_num):
                pos2_XY_XW = np.matmul(rotation_matrix_XW(theta2), pos2_XY.transpose()).T
                
                mean_distance = compute_angular_distances(
                    pos1_axis, pos2_XY_XW)
                mean_distances.append(mean_distance)

                pos1_rotation_matrices.append(m_i)
                pos2_rotation_matrices.append(m_j)
                thetas1.append(theta1)
                thetas2.append(theta2)

    results = pd.DataFrame(mean_distances).T
    min_distance_idx = np.argmin(results.mean(axis=0).values)
    min_pos1_rotation_matrix = pos1_rotation_matrices[min_distance_idx]
    min_pos2_rotation_matrix = pos2_rotation_matrices[min_distance_idx]
    min_theta1 = thetas1[min_distance_idx]
    min_theta2 = thetas2[min_distance_idx]

    out_pos1 = np.array([(min_pos1_rotation_matrix @ v).A1 for v in pos1])
    out_pos2 = np.array([(min_pos2_rotation_matrix @ v).A1 for v in pos2])
    out_pos2 = np.matmul(rotation_matrix_XY(min_theta1), out_pos2.transpose()).T
    out_pos2 = np.matmul(rotation_matrix_XW(min_theta2), out_pos2.transpose()).T
    out_pos2 = np.array(out_pos2)
    
    out_pos1_spherical = np.array([euclidean_to_hyperspherical_coordinates(v) for v in out_pos1])
    out_pos2_spherical = np.array([euclidean_to_hyperspherical_coordinates(v) for v in out_pos2])
    return out_pos1, out_pos2, out_pos1_spherical, out_pos2_spherical


def apply_pipeline_matrix_with_loading_S3(real_coords_path, inf_coords_path, theta_num=20, new_version=False):
    # Load nodes' positions
    real_coords = pd.read_csv(
        real_coords_path, sep="\s+", comment="#", header=None)
    if new_version:
        real_coords.columns = ['index', 'kappa', 'radius', 'pos0',
                            'pos1', 'pos2', 'pos3', 'realdeg', 'expdegree']
    else:
        real_coords.columns = ['index', 'kappa', 'pos0',
                            'pos1', 'pos2', 'pos3', 'realdeg', 'expdegree']
    inf_coords = pd.read_csv(
        inf_coords_path, comment="#", header=None, sep="\s+")
    if new_version:
        inf_coords.columns = ['index', 'inf_kappa', 'inf_hyp_radius',
                            'inf_pos0', 'inf_pos1', 'inf_pos2', 'inf_pos3']
    else:        
        inf_coords.columns = ['index', 'inf_kappa',
                            'inf_pos0', 'inf_pos1', 'inf_pos2', 'inf_pos3']
    df = inf_coords.merge(real_coords, on="index")
    real_coords_all = df[["pos0", "pos1", "pos2", "pos3"]].values
    inf_coords_all = df[["inf_pos0", "inf_pos1", "inf_pos2", "inf_pos3"]].values

    r1 = np.mean([np.linalg.norm(r) for r in real_coords_all])
    r2 = np.mean([np.linalg.norm(r) for r in inf_coords_all])
    real_coords_all /= r1
    inf_coords_all /= r2
    return apply_pipeline_matrix_S3(real_coords_all, inf_coords_all, theta_num=theta_num)
