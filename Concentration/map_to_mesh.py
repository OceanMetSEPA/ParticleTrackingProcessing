import numpy as np
from scipy.spatial import cKDTree
from concurrent.futures import ThreadPoolExecutor
import os
import warnings
import scipy.io as sio

def _load_particle_dict(input_arg):
    """
    Load particle_dict from input:
    - dict: returned as is
    - file: load MAT-file
    - folder: find [siteName]_trackStruct.mat inside folder
    """
    if isinstance(input_arg, dict):
        return input_arg, None  # no siteName
    elif isinstance(input_arg, str):
        if os.path.isfile(input_arg):
            matfile = input_arg
        elif os.path.isdir(input_arg):
            # find [foldername]_trackStruct.mat
            site_name = os.path.basename(os.path.normpath(input_arg))
            matfile = os.path.join(input_arg, f"{site_name}_trackStruct.mat")
            if not os.path.exists(matfile):
                warnings.warn(f"{matfile} not found.")
                return None, site_name
        else:
            raise ValueError("Input must be a dict, a MAT file path, or a folder path.")
        data = sio.loadmat(matfile)
        site_name = os.path.splitext(os.path.basename(matfile))[0].replace("_trackStruct","")
        return data, site_name
    else:
        raise TypeError("Input must be a dict, a MAT file path, or a folder path.")


def map_particle_tracks_to_mesh(
    particle_input,
    dfsu_dict,
    matlab_indexing=False,
    print_every=1000,
    n_workers=4,
    save_result=True
):
    """
    Map particle tracks to mesh. particle_input can be:
        - dict
        - path to MAT-file
        - folder containing [siteName]_trackStruct.mat
    
    The resulting mapped struct is optionally saved to [siteName]_trackStruct_dfsu.mat
    in the same folder as the input MAT-file or folder.
    """
    # Load particle_dict
    particle_dict, site_name = _load_particle_dict(particle_input)
    if particle_dict is None:
        return None  # couldn't find MAT-file

    if site_name is None:
        site_name = "trackStruct"

    out_dir = None
    if isinstance(particle_input, str):
        out_dir = os.path.dirname(particle_input)
        if out_dir == "":
            out_dir = os.getcwd()
    else:
        out_dir = os.getcwd()

    save_path = os.path.join(out_dir, f"{site_name}_trackStruct_dfsu.mat")

    print("Starting parallel mapping of particle tracks to mesh...")

    xMesh = np.ravel(dfsu_dict["xMesh"])
    yMesh = np.ravel(dfsu_dict["yMesh"])
    meshIndices = np.array(dfsu_dict["meshIndices"], dtype=int)  # [N_cells, 3]

    # Build KD-tree
    tree = cKDTree(np.c_[xMesh, yMesh])
    print("KD-tree built.")

    # Precompute vertex -> triangle mapping
    vertex_to_tri = [[] for _ in range(len(xMesh))]
    for tri_idx, verts in enumerate(meshIndices):
        for v in verts:
            vertex_to_tri[v].append(tri_idx)

    # Timesteps
    vol_times = dfsu_dict["dateTime"].ravel()
    part_times = particle_dict["dateTime"].ravel()
    intersect_times = np.intersect1d(vol_times, part_times)

    print(f"Volume timesteps: {len(vol_times)}, Particle timesteps: {len(part_times)}")
    print(f"Intersecting timesteps: {len(intersect_times)}")

    vol_idx_map = {t: i for i, t in enumerate(vol_times)}
    part_idx_map = {t: i for i, t in enumerate(part_times)}

    N_particles = particle_dict["x"].shape[0]
    N_timesteps = len(intersect_times)

    meshIndex = np.full((N_particles, N_timesteps), np.nan)
    waterSurface = np.full((N_particles, N_timesteps), np.nan)
    depthBelowSurface = np.full((N_particles, N_timesteps), np.nan)

    def process_timestep(ti_tup):
        ti, t = ti_tup
        ti_vol = vol_idx_map[t]
        ti_part = part_idx_map[t]

        meshIndex_t = np.full(N_particles, np.nan)
        waterSurface_t = np.full(N_particles, np.nan)
        depthBelowSurface_t = np.full(N_particles, np.nan)

        xp = particle_dict["x"][:, ti_part]
        yp = particle_dict["y"][:, ti_part]
        zp = particle_dict["z"][:, ti_part]

        _, node_ids = tree.query(np.c_[xp, yp])

        for pi, nid in enumerate(node_ids):
            if pi % print_every == 0:
                print(f"Timestep {ti+1}, particle {pi+1}/{N_particles}")

            candidate_tri = vertex_to_tri[nid]
            if not candidate_tri:
                continue

            # Check which candidate triangle contains particle
            for tri_idx in candidate_tri:
                verts = meshIndices[tri_idx]
                xverts = xMesh[verts]
                yverts = yMesh[verts]
                zverts = dfsu_dict["zCoordinate"][verts, ti_vol]

                # Barycentric coords
                x1, x2, x3 = xverts
                y1, y2, y3 = yverts
                detT = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1)
                if detT == 0:
                    inside = False
                else:
                    l1 = ((y2 - y3)*(xp[pi]-x3) + (x3 - x2)*(yp[pi]-y3)) / detT
                    l2 = ((y3 - y1)*(xp[pi]-x3) + (x1 - x3)*(yp[pi]-y3)) / detT
                    l3 = 1 - l1 - l2
                    inside = (l1 >= -1e-12) & (l2 >= -1e-12) & (l3 >= -1e-12)
                if inside:
                    meshIndex_t[pi] = tri_idx + 1 if matlab_indexing else tri_idx
                    top = np.mean(zverts) if detT==0 else l1*zverts[0]+l2*zverts[1]+l3*zverts[2]
                    waterSurface_t[pi] = top
                    depthBelowSurface_t[pi] = top - zp[pi]
                    break  # stop at first triangle
        return ti, meshIndex_t, waterSurface_t, depthBelowSurface_t

    # Parallel execution
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        results = list(executor.map(process_timestep, enumerate(intersect_times)))

    # Collect results
    for ti, mi_t, ws_t, dbs_t in results:
        meshIndex[:, ti] = mi_t
        waterSurface[:, ti] = ws_t
        depthBelowSurface[:, ti] = dbs_t

    # Build mapped_dict
    mapped_dict = particle_dict.copy()
    mapped_dict["dateTime"] = intersect_times.reshape(-1,1)
    mapped_dict["index"] = np.arange(1, N_timesteps+1).reshape(-1,1)
    mapped_dict["meshIndex"] = meshIndex
    mapped_dict["waterSurface"] = waterSurface
    mapped_dict["depthBelowSurface"] = depthBelowSurface

    # Filter particle fields
    for k, v in particle_dict.items():
        if isinstance(v, np.ndarray) and v.ndim==2 and v.shape[1]==len(part_times):
            idxs = [part_idx_map[t] for t in intersect_times]
            mapped_dict[k] = v[:, idxs]

    # Save result if requested
    if save_result:
        sio.savemat(save_path, mapped_dict)
        print(f"Saved mapped struct to {save_path}")

    print("Mapping complete!")
    return mapped_dict
