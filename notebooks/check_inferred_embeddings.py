from numba import jit
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy import stats


@jit(nopython=True)
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


@jit(nopython=True)
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


def get_rotation_matrix(vec, axis=np.array([1, 0, 0])):
    # From: https://math.stackexchange.com/a/476311
    a = vec / np.linalg.norm(vec)
    b = axis
    v = np.cross(a, b)
    s = np.linalg.norm(v)
    c = np.dot(a, b)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    r = np.eye(3) + vx + np.dot(vx, vx) * (1-c)/(s**2)
    return np.matrix(r)


def rotation_along_Xaxis(vec, theta):
    m = np.matrix([[1, 0, 0],
                   [0, np.cos(theta), np.sin(theta)],
                   [0, -np.sin(theta), np.cos(theta)]])
    return (m @ vec).A1


def rotation_matrix_XY(theta):
    m = np.matrix([[1, 0, 0],
                   [0, np.cos(theta), np.sin(theta)],
                   [0, -np.sin(theta), np.cos(theta)]])
    return m


def rotation_matrix_Z(theta):
    m = np.matrix([[np.cos(theta), -np.sin(theta), 0],
                   [np.sin(theta), np.cos(theta), 0],
                   [0, 0, 1]])
    return m


def rotation_matrix_XYZ(gamma, beta, alpha):
    X11 = np.cos(alpha) * np.cos(beta)
    X12 = np.cos(alpha) * np.sin(beta) * np.sin(gamma) - np.sin(alpha) * np.cos(gamma)
    X13 = np.cos(alpha) * np.sin(beta) * np.cos(gamma) + np.sin(alpha) * np.sin(gamma)
    
    X21 = np.sin(alpha) * np.cos(beta)
    X22 = np.sin(alpha) * np.sin(beta) * np.sin(gamma) + np.cos(alpha) * np.cos(gamma)
    X23 = np.sin(alpha) * np.sin(beta) * np.cos(gamma) - np.cos(alpha) * np.sin(gamma)

    X31 = -np.sin(beta)
    X32 = np.cos(beta) * np.sin(gamma)
    X33 = np.cos(beta) * np.cos(gamma)
    
    m = np.matrix([[X11, X12, X13],
                   [X21, X22, X23],
                   [X31, X32, X33]])
    return m



def apply_pipeline_matrix(pos1, pos2, theta_num=50):
    mean_distances = []
    pos1_rotation_matrices = []
    pos2_rotation_matrices = []
    thetas = []

    for i in (range(len(pos1))):
        node_i = pos1[i]
        node_i_prime = pos2[i]

        m_i = get_rotation_matrix(node_i)
        m_j = get_rotation_matrix(node_i_prime)

        pos1_axis = np.array([(m_i @ v).A1 for v in pos1])
        pos2_axis = np.array([(m_j @ v).A1 for v in pos2])

        for theta in np.linspace(0, 2*np.pi, num=theta_num):
            # rotate only the second coordinates
            #pos2_axis_theta = np.array(
            #    [rotation_along_Xaxis(v, theta) for v in pos2_axis])
            pos2_axis_theta = np.matmul(rotation_matrix_XY(theta), pos2_axis.transpose()).T

            mean_distance = compute_angular_distances(pos1_axis, pos2_axis_theta)
            mean_distances.append(np.mean(mean_distance))

            pos1_rotation_matrices.append(m_i)
            pos2_rotation_matrices.append(m_j)
            thetas.append(theta)


    min_distance_idx = np.argmin(np.array(mean_distances))
    min_pos1_rotation_matrix = pos1_rotation_matrices[min_distance_idx]
    min_pos2_rotation_matrix = pos2_rotation_matrices[min_distance_idx]
    min_theta = thetas[min_distance_idx]

    out_pos1 = np.array([(min_pos1_rotation_matrix @ v).A1 for v in pos1])
    out_pos2 = np.array([(min_pos2_rotation_matrix @ v).A1 for v in pos2])
    out_pos2 = np.matmul(rotation_matrix_XY(min_theta), out_pos2.transpose()).T
    out_pos2 = np.array(out_pos2)

    out_pos1_spherical = np.array([euclidean_to_hyperspherical_coordinates(v) for v in out_pos1])
    out_pos2_spherical = np.array([euclidean_to_hyperspherical_coordinates(v) for v in out_pos2])

    rotation_dict = {'min_pos1_rotation_matrix': min_pos1_rotation_matrix, 
                     'min_pos2_rotation_matrix': min_pos2_rotation_matrix,
                     'min_theta': min_theta}
    return out_pos1, out_pos2, out_pos1_spherical, out_pos2_spherical, rotation_dict


def load_coordinates(real_coords_path, inf_coords_path, cutoff_percent=0, min_degree=None, is_gen_coord=False):
    # INFO: Real coords and inf coords are just two .inf_coord files or .gen_coord files (change optional parameter) 
    real_coords = pd.read_csv(real_coords_path, sep="\s+", comment="#", header=None)
    inf_coords = pd.read_csv(inf_coords_path, sep="\s+", comment="#", header=None)

    if is_gen_coord:
        real_coords.columns = ['index', 'kappa', 'radius', 'pos0', 'pos1', 'pos2', 'real_deg', 'exp_deg']
        inf_coords.columns = ['index', 'inf_kappa', 'inf_hyp_rad', 'inf_pos0', 'inf_pos1', 'inf_pos2', 'real_deg', 'exp_deg']
    else:
        real_coords.columns = ['index', 'kappa', 'hyp_rad', 'pos0', 'pos1', 'pos2']
        inf_coords.columns = ['index', 'inf_kappa', 'inf_hyp_rad', 'inf_pos0', 'inf_pos1', 'inf_pos2']
    
    df = inf_coords.merge(real_coords, on="index")
    
    real_coords_all = df[["pos0", "pos1", "pos2"]].values
    inf_coords_all = df[["inf_pos0", "inf_pos1", "inf_pos2"]].values

    r1 = np.linalg.norm(real_coords_all[0])
    r2 = np.linalg.norm(inf_coords_all[0])
    real_coords_all /= r1
    inf_coords_all /= r2
    return real_coords_all, inf_coords_all


def apply_pipeline_matrix_with_loading(real_coords_path, inf_coords_path, cutoff_percent=0.005, theta_num=50, is_gen_coord=False):
    real_coords_all, inf_coords_all = load_coordinates(real_coords_path, inf_coords_path, cutoff_percent, is_gen_coord=is_gen_coord)
    return apply_pipeline_matrix(real_coords_all, inf_coords_all, theta_num=theta_num)


def apply_pipeline_matrix_with_loading_and_fix_rotation(real_coords_path, inf_coords_path, cutoff_percent=0.005, theta_num=50, final_theta_z=None, is_gen_coord=False):
    real_coords_all, best_inf_coords_euclidean, real_coords_spherical, best_inf_coords_spherical, rotation_dict =  \
        apply_pipeline_matrix_with_loading(real_coords_path, inf_coords_path, cutoff_percent, theta_num, is_gen_coord=is_gen_coord)

    # Rotate on Z-axis and find the angle corresponding to the maximum value of pearson correlation for the second angle (phi2)
    all_thetas = np.linspace(0, 2*np.pi, theta_num)
    all_pearson_phi1 = []
    all_pearson_phi2 = []

    for theta_z in all_thetas:
        rotate_best_inf_coords_euclidean = np.matmul(rotation_matrix_XY(theta_z), 
                                                    best_inf_coords_euclidean.transpose()).T
        rotate_best_inf_coords_euclidean = np.array(rotate_best_inf_coords_euclidean)

        rotate_best_inf_coords_spherical = np.array(
            [euclidean_to_hyperspherical_coordinates(v) for v in rotate_best_inf_coords_euclidean])

        x = real_coords_spherical[:, 1]
        y = rotate_best_inf_coords_spherical[:, 1]
        pearson_phi1 = stats.pearsonr(x, y)[0]
        all_pearson_phi1.append(pearson_phi1)
        
        x = real_coords_spherical[:, 2]
        y = rotate_best_inf_coords_spherical[:, 2]
        pearson_phi2 = stats.pearsonr(x, y)[0]
        all_pearson_phi2.append(pearson_phi2)

    idx = np.argmax(np.abs(np.array(all_pearson_phi2)))
    if final_theta_z is None:
        theta_z = all_thetas[idx]
    else:
        theta_z = final_theta_z
    #print('The best Z-axis rotation angle: ', theta_z)
    rotate_best_inf_coords_euclidean = np.matmul(rotation_matrix_XY(theta_z), 
                                                best_inf_coords_euclidean.transpose()).T
    rotate_best_inf_coords_spherical = np.array(
        [euclidean_to_hyperspherical_coordinates(v) for v in np.array(rotate_best_inf_coords_euclidean)])
    if all_pearson_phi1[idx] < 0:
        rotate_best_inf_coords_spherical[:, 1] = np.pi - rotate_best_inf_coords_spherical[:, 1]
    if all_pearson_phi2[idx] < 0:
        rotate_best_inf_coords_spherical[:, 2] = 2*np.pi - rotate_best_inf_coords_spherical[:, 2]

    rotate_best_inf_coords_euclidean = np.array(
        [hyperspherical_to_euclidean_coordinates(v) for v in rotate_best_inf_coords_spherical])

    rotation_dict['theta_z'] = theta_z
    rotation_dict['all_pearson_phi1[idx] < 0'] = all_pearson_phi1[idx] < 0
    rotation_dict['all_pearson_phi2[idx] < 0'] = all_pearson_phi2[idx] < 0
    return real_coords_all, rotate_best_inf_coords_euclidean, real_coords_spherical, rotate_best_inf_coords_spherical, rotation_dict


def apply_pipeline_matrix_with_loading_and_rotate_all_euclidean_use_all_nodes(real_coords_path, inf_coords_path, cutoff_percent=0.005, theta_num=50, is_gen_coord=False):
    _, _, _, _, rotation_dict =  \
        apply_pipeline_matrix_with_loading_and_fix_rotation(real_coords_path, inf_coords_path, cutoff_percent, theta_num, is_gen_coord=is_gen_coord)

    # Load all nodes
    real_coords_all, inf_coords_all = load_coordinates(real_coords_path, inf_coords_path, min_degree=0, is_gen_coord=is_gen_coord)
    
    # First rotation (aligning the nodes)
    real_coords_all = np.array([(rotation_dict['min_pos1_rotation_matrix'] @ v).A1 for v in real_coords_all])
    inf_coords_all = np.array([(rotation_dict['min_pos2_rotation_matrix'] @ v).A1 for v in inf_coords_all])
    inf_coords_all = np.array(np.matmul(rotation_matrix_XY(rotation_dict['min_theta']), inf_coords_all.transpose()).T)

    real_coords_spherical = np.array([euclidean_to_hyperspherical_coordinates(v) for v in real_coords_all])
    inf_coords_spherical = np.array([euclidean_to_hyperspherical_coordinates(v) for v in inf_coords_all])

    # Second rotation (along X-axis)
    inf_coords_all = np.matmul(rotation_matrix_XY(rotation_dict['theta_z']), inf_coords_all.transpose()).T
    inf_coords_spherical = np.array([euclidean_to_hyperspherical_coordinates(v) for v in np.array(inf_coords_all)])
    if rotation_dict['all_pearson_phi1[idx] < 0']:
        inf_coords_spherical[:, 1] = np.pi - inf_coords_spherical[:, 1]
    if rotation_dict['all_pearson_phi2[idx] < 0']:
        inf_coords_spherical[:, 2] = 2*np.pi - inf_coords_spherical[:, 2]
    inf_coords_all = np.array([hyperspherical_to_euclidean_coordinates(v) for v in inf_coords_spherical])
   
    return real_coords_all, inf_coords_all, real_coords_spherical, inf_coords_spherical


def apply_pipeline_matrix_S2(pos1, pos2, theta_num=20):
    out_pos1, out_pos2, out_pos1_spherical, out_pos2_spherical, _ = apply_pipeline_matrix(pos1, pos2, theta_num=theta_num)

    all_thetas = np.linspace(0, 2*np.pi, theta_num)
    all_pearson_phi1 = []
    all_pearson_phi2 = []

    for theta_z in all_thetas:
        rotate_out_pos2 = np.matmul(rotation_matrix_XY(theta_z), 
                                                    out_pos2.transpose()).T
        rotate_out_pos2 = np.array(rotate_out_pos2)

        rotate_out_pos2_spherical = np.array(
            [euclidean_to_hyperspherical_coordinates(v) for v in rotate_out_pos2])

        x = out_pos1_spherical[:, 1]
        y = rotate_out_pos2_spherical[:, 1]
        pearson_phi1 = stats.pearsonr(x, y)[0]
        all_pearson_phi1.append(pearson_phi1)
        
        x = out_pos1_spherical[:, 2]
        y = rotate_out_pos2_spherical[:, 2]
        pearson_phi2 = stats.pearsonr(x, y)[0]
        all_pearson_phi2.append(pearson_phi2)

    idx = np.argmax(np.abs(np.array(all_pearson_phi2)))
    theta_z = all_thetas[idx]

    #print('The best Z-axis rotation angle: ', theta_z)
    rotate_out_pos2 = np.matmul(rotation_matrix_XY(theta_z), out_pos2.transpose()).T
    rotate_out_pos2_spherical = np.array(
        [euclidean_to_hyperspherical_coordinates(v) for v in np.array(rotate_out_pos2)])
    
    if all_pearson_phi1[idx] < 0:
        rotate_out_pos2_spherical[:, 1] = np.pi - rotate_out_pos2_spherical[:, 1]
    if all_pearson_phi2[idx] < 0:
        rotate_out_pos2_spherical[:, 2] = 2*np.pi - rotate_out_pos2_spherical[:, 2]

    rotate_out_pos2 = np.array(
        [hyperspherical_to_euclidean_coordinates(v) for v in rotate_out_pos2_spherical])

    return out_pos1, rotate_out_pos2, out_pos1_spherical, rotate_out_pos2_spherical


