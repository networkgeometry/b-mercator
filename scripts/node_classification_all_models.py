from karateclub import DeepWalk, Role2Vec, LaplacianEigenmaps, NetMF
from karateclub import FeatherNode, MUSAE
import umap
from tqdm import tqdm
import numpy as np
import scipy as sp
# shim scipy.errstate -> numpy.errstate
if not hasattr(sp, "errstate"):
    sp.errstate = np.errstate

from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score, accuracy_score
from sklearn.neighbors import KNeighborsClassifier


def node_classification(dim, df, df_nodes, test_size=0.8, k_neighbours=10, n_times=5):
    if dim == 1:
        cols = "inf_theta"
    else:
        cols = [f"inf_p{i}" for i in range(dim + 1)]

    f1_unipartite, accuracy_unipartite = [], []
    f1_bipartite, accuracy_bipartite = [], []

    for _ in tqdm(range(n_times)):
        y_pred, y_true = predict_labels_kneighbours(
            dim, df[cols].values, df["class"].values, test_size, k_neighbours
        )
        f1_unipartite.append(f1_score(y_true, y_pred, average="micro"))
        accuracy_unipartite.append(accuracy_score(y_true, y_pred))

        y_pred, y_true = predict_labels_kneighbours(
            dim,
            df_nodes[cols].values,
            df_nodes["class"].values,
            test_size,
            k_neighbours,
        )
        f1_bipartite.append(f1_score(y_true, y_pred, average="micro"))
        accuracy_bipartite.append(accuracy_score(y_true, y_pred))

    return {
        "f1_unipartite": f1_unipartite,
        "f1_bipartite": f1_bipartite,
        "accuracy_unipartite": accuracy_unipartite,
        "accuracy_bipartite": accuracy_bipartite,
    }


def predict_labels_kneighbours(dim, pos, labels, test_size, k_neighbours):
    if dim == 1:
        pos = pos.reshape(-1, 1)
    pos_train, pos_test, labels_train, labels_test = train_test_split(
        pos, labels, test_size=test_size
    )

    metric = compute_angle_S1 if dim == 1 else compute_angle_SD
    neigh = KNeighborsClassifier(
        n_neighbors=k_neighbours, metric=metric, weights="distance"
    )
    neigh.fit(pos_train, labels_train)
    predicted_labels = neigh.predict(pos_test)
    return predicted_labels, labels_test


def compute_angle_SD(p1, p2):
    return np.arccos(
        np.clip(np.dot(p1, p2) / (np.linalg.norm(p1) * np.linalg.norm(p2)), -1, 1)
    )


def compute_angle_S1(t1, t2):
    return np.pi - np.fabs(np.pi - np.fabs(t1 - t2))


def train_embeddings(method, g_relabeled, dimensions, node_feature_matrix=None):
    if method == "DeepWalk":
        model = DeepWalk(dimensions=dimensions)
    elif method == "Role2Vec":
        model = Role2Vec(dimensions=dimensions)
    elif method == "LaplacianEigenmaps":
        model = LaplacianEigenmaps(dimensions=dimensions)
    elif method == "NetMF":
        model = NetMF(dimensions=dimensions)
    elif method == "FeatherNode":
        model = FeatherNode(reduction_dimensions=dimensions)
    elif method == "MUSAE":
        model = MUSAE(dimensions=dimensions)
    elif method == "UMAP":
        model = umap.UMAP(n_components=dimensions)

    if method in ["FeatherNode", "MUSAE"]:
        model.fit(g_relabeled, node_feature_matrix)
    elif method in ["UMAP"]:
        return model.fit_transform(node_feature_matrix)
    else:
        model.fit(g_relabeled)

    return model.get_embedding()


def node_classification_external(
    method,
    dim,
    g_relabeled,
    y_values,
    node_feature_matrix=None,
    test_size=0.8,
    k_neighbours=10,
    n_times=5,
):
    emb = train_embeddings(
        method, g_relabeled, dim, node_feature_matrix=node_feature_matrix
    )
    f1_vals, accuracy_vals = [], []

    for _ in range(n_times):
        y_pred, y_true = predict_labels_kneighbours_external(
            emb, y_values, test_size, k_neighbours
        )
        f1_vals.append(f1_score(y_true, y_pred, average="micro"))
        accuracy_vals.append(accuracy_score(y_true, y_pred))

    return {
        "f1_vals": f1_vals,
        "accuracy_vals": accuracy_vals,
    }


def predict_labels_kneighbours_external(pos, labels, test_size, k_neighbours):
    pos_train, pos_test, labels_train, labels_test = train_test_split(
        pos, labels, test_size=test_size
    )
    neigh = KNeighborsClassifier(
        n_neighbors=k_neighbours, 
        metric="cosine",
        algorithm="brute",  # required for cosine in sklearn
        weights="distance"
    )
    neigh.fit(pos_train, labels_train)
    predicted_labels = neigh.predict(pos_test)
    return predicted_labels, labels_test
