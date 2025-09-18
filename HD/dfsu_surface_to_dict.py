import numpy as np
import pandas as pd
from mikeio import Dfsu
import datetime
import os
import scipy.io as sio

def _to_matlab_datenum(ts):
    if isinstance(ts, (np.datetime64, pd.Timestamp)):
        dt = pd.to_datetime(ts).to_pydatetime()
    elif isinstance(ts, datetime.datetime):
        dt = ts
    else:
        raise TypeError("Unsupported timestamp type for conversion")
    day_frac = (dt.hour * 3600 + dt.minute * 60 + dt.second + dt.microsecond / 1e6) / 86400.0
    return dt.toordinal() + 366 + day_frac

def dfsu_to_surface_dict(fname: str, dt: int = 1, matlab_index: bool = True, overwrite: bool = False):
    """
    Convert a DFSU file to a MATLAB-compatible dict and save as dfsu.mat in the same folder.
    Automatically loads dfsu.mat if it exists (unless overwrite=True).
    """
    out_dir = os.path.dirname(fname)
    mat_path = os.path.join(out_dir, "dfsu.mat")

    # Load existing file if available
    if os.path.exists(mat_path) and not overwrite:
        print(f"Loading existing {mat_path}")
        return sio.loadmat(mat_path)

    # Read DFSU
    dfs = Dfsu(fname)
    geom = dfs.geometry

    # Select time indices
    Nt_total = dfs.n_timesteps
    time_indices = list(range(0, Nt_total, int(dt)))
    if not time_indices:
        raise ValueError("No timesteps selected")
    ds = dfs.read(items=0, time=time_indices, keepdims=True)

    # Surface zn
    top_el = geom.top_elements
    node_ids_surf, _ = geom._get_nodes_and_table_for_elements(top_el, node_layers="top")
    zn_all = ds[0]._zn
    zn_surf = zn_all[:, node_ids_surf].T
    zn_below = zn_all[:, node_ids_surf - 1].T
    dzVertex = zn_surf - zn_below

    # Mesh nodes
    geom2d = geom.geometry2d
    node_coords_2d = np.asarray(geom2d.node_coordinates)
    xMesh = node_coords_2d[:, 0].astype(np.float64).reshape(-1, 1)
    yMesh = node_coords_2d[:, 1].astype(np.float64).reshape(-1, 1)
    zMesh = node_coords_2d[:, 2].astype(np.float64).reshape(-1, 1)

    # Element coordinates (collapsed in vertical)
    ele_coords = np.asarray(geom.element_coordinates)
    Nlayers = geom.n_layers
    Ne_total = ele_coords.shape[0]
    Ne_horiz = Ne_total // Nlayers

    elem3d = ele_coords.reshape(Ne_horiz, Nlayers, 3)
    x_elem = elem3d[:, 0, 0].astype(np.float64).reshape(-1, 1)
    y_elem = elem3d[:, 0, 1].astype(np.float64).reshape(-1, 1)
    z_elem_layers = elem3d[:, :, 2].astype(np.float64)

    # Time
    times = dfs.time[time_indices]
    dateTime = np.array([_to_matlab_datenum(t) for t in times], dtype=np.float64).reshape(-1, 1)

    # Surface connectivity
    elem_table_2d = geom2d.element_table
    max_nodes = max(len(row) for row in elem_table_2d)
    meshIndices = np.full((len(elem_table_2d), max_nodes), -1, dtype=np.int32)
    for i, row in enumerate(elem_table_2d):
        meshIndices[i, :len(row)] = row
    if matlab_index:
        meshIndices += 1

    # Assemble output dict
    out = {
        "name": os.path.basename(fname),
        "dateTime": dateTime,
        "x": x_elem,
        "y": y_elem,
        "z": z_elem_layers,
        "xMesh": xMesh,
        "yMesh": yMesh,
        "zMesh": zMesh,
        "zCoordinate": zn_surf.astype(np.float64),
        "dzVertex": dzVertex.astype(np.float64),
        "meshIndices": meshIndices,
    }

    # Save as MATLAB binary file (v7) — works with MATLAB load
    sio.savemat(mat_path, out)
    print(f"Saved results to {mat_path}")

    print("Finished. Sizes:")
    for k, v in out.items():
        if isinstance(v, np.ndarray):
            print(f"  {k}: {v.shape}")
        else:
            print(f"  {k}: {type(v)}")

    return out
