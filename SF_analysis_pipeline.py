"""
SF_analysis_pipeline.py

Pipeline for single-fiber myonuclei analysis:
- Fiji macro produces STDIP, Skel, and Z-stack images in folders:
  [Subject]/[Timepoint]/[Leg]/STDIP|Skel|TIFs
- This script segments nuclei, computes morphometrics, orientation vs. skeleton,
  assigns Z positions, clusters nuclei in 3D, measures fiber diameter along the skeleton,
  computes nucleus-to-skeleton distances, and produces biopsy-level summaries merged with Imaris.

Outputs:
  * Per-fiber CSVs: included nuclei, excluded nuclei (with reasons), fiber width profile
  * Overlay PNGs
  * Subject workbook with one sheet per timepoint (included nuclei only)
  * Biopsy-level summary CSV per leg (integrated with Imaris master)
  * Excluded fibers CSV per leg (no Imaris match)

Dependencies: numpy, pandas, scikit-image, scikit-learn, scipy, openpyxl, pillow
"""

# Suggested CLI input:
# python SF_analysis_pipeline.py "J:\MacroOutput Test" --imaris_master "C:\Users\Dreyer Lab\Desktop\23_Wu_SF_ImarisResults.xlsx" --min_area_px 200 --max_area_px 0 --pixel_size_xy_um 0.3291532 --z_scale_um_per_index 2.7 --z_std_threshold 2.0 --dbscan_eps_um 20 --skeleton_radius_px 20


import os
import re
import warnings
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from skimage import io, morphology, measure
from skimage.segmentation import find_boundaries
from sklearn.cluster import DBSCAN
from PIL import Image, ImageDraw, ImageFont
from scipy.ndimage import distance_transform_edt
from skimage.graph import route_through_array
from scipy.spatial import distance

# ----------------------------- Filename parsing ------------------------------
FILENAME_REGEX = re.compile(
    r'(?:.*_)?23_Wu_(?P<Subject>\d+)_'           # Subject numeric
    r'(?P<Day>[A-Za-z0-9]+)_'                    # Timepoint token
    r'(?P<Side>[LR])_?(?P<Leg>\d+)?\.lif'        # Side L/R, optional leg code
    r'(?:\s*[-_]\s*)+Fiber(?P<Fiber>\d+-\d+)_.*\.(?:tif|tiff)$',
    re.IGNORECASE
)

# ----------------------------- Utilities ------------------------------------
def _ensure_dir(p: str):
    os.makedirs(p, exist_ok=True)

def _imread(path: str):
    return io.imread(path)

def _to_bool(img: np.ndarray) -> np.ndarray:
    return img if img.dtype == bool else (img > 0)

def _label_mask(bin_mask: np.ndarray, min_area_px: int = 0) -> np.ndarray:
    lab = measure.label(bin_mask, connectivity=2)
    if min_area_px and min_area_px > 0:
        lab = morphology.remove_small_objects(lab, min_size=min_area_px)
    return lab

def _classify_shape(aspect_ratio: float) -> str:
    if aspect_ratio is None or np.isnan(aspect_ratio):
        return "Unknown"
    if aspect_ratio > 0.8:
        return "Spherical"
    elif aspect_ratio < 0.55:
        return "Ellipsoid"
    else:
        return "Intermediate"

def _orientation_deg_from_rad(theta: float) -> float:
    if np.isnan(theta):
        return np.nan
    return float(np.degrees(theta) % 180.0)

def _principal_angle_from_local_points(xy: np.ndarray) -> float:
    if xy.shape[0] < 2:
        return np.nan
    mean = xy.mean(axis=0)
    cov = np.cov((xy - mean).T)
    w, v = np.linalg.eig(cov)
    vec = v[:, np.argmax(w)]
    return float(np.degrees(np.arctan2(vec[1], vec[0])) % 180.0)

def _relative_angle(a: float, b: float) -> float:
    if np.isnan(a) or np.isnan(b):
        return np.nan
    d = abs(a - b)
    return float(180.0 - d if d > 90.0 else d)

def _parse_fiber_num_and_type(fiber_code: str) -> Tuple[Optional[int], Optional[int]]:
    """
    Parse codes like '2-1' into (fiber_num=2, fiber_type=1).
    """
    m = re.match(r'^\s*(\d+)-(\d+)\s*$', str(fiber_code))
    if not m:
        return None, None
    return int(m.group(1)), int(m.group(2))

def _is_endpoint(y: int, x: int, skel: np.ndarray) -> bool:
    """Pixel with exactly one 8-connected skeleton neighbor."""
    y0, y1 = max(0, y - 1), min(skel.shape[0], y + 2)
    x0, x1 = max(0, x - 1), min(skel.shape[1], x + 2)
    window = skel[y0:y1, x0:x1]
    return (window.sum() - 1) == 1

def _find_skeleton_endpoints(skel: np.ndarray) -> List[Tuple[int, int]]:
    """Return list of (row, col) endpoints in a skeleton."""
    ys, xs = np.where(skel > 0)
    coords = np.column_stack((ys, xs))
    endpoints = [(int(y), int(x)) for y, x in coords if _is_endpoint(int(y), int(x), skel)]
    return endpoints

# --------------------------- Core computations -------------------------------
@dataclass
class Params:
    pixel_size_xy_um: float
    z_scale_um_per_index: float
    z_std_threshold: float
    dbscan_eps_um: float 
    dbscan_min_samples: int
    skeleton_radius_px: float 
    min_area_px: int 
    max_area_px: int
    fiber_width_step_um: float
    fiber_width_max_radius_um: float

def measure_nuclei_from_binary(mask_bin: np.ndarray,
                               params: Params) -> Tuple[pd.DataFrame, np.ndarray]:
    lab = _label_mask(mask_bin, min_area_px=0)
    props = measure.regionprops_table(
        lab,
        properties=('label', 'area', 'perimeter', 'centroid',
                    'major_axis_length', 'minor_axis_length', 'orientation')
    )
    df = pd.DataFrame(props).rename(columns={
        'label': 'Label',
        'area': 'Area_px',
        'perimeter': 'Perimeter_px',
        'centroid-0': 'Centroid_Y',
        'centroid-1': 'Centroid_X',
        'major_axis_length': 'MajorAxis_px',
        'minor_axis_length': 'MinorAxis_px',
        'orientation': 'Orientation_rad'
    })

    df['Aspect_Ratio'] = df['MinorAxis_px'] / df['MajorAxis_px']
    df['Orientation_deg'] = df['Orientation_rad'].apply(_orientation_deg_from_rad)
    df['Shape_Class'] = df['Aspect_Ratio'].apply(_classify_shape)

    px_um = float(params.pixel_size_xy_um)
    df['Area_um2'] = df['Area_px'] * (px_um ** 2)
    df['Perimeter_um'] = df['Perimeter_px'] * px_um
    df['MajorAxis_um'] = df['MajorAxis_px'] * px_um
    df['MinorAxis_um'] = df['MinorAxis_px'] * px_um
    df['Centroid_X_um'] = df['Centroid_X'] * px_um
    df['Centroid_Y_um'] = df['Centroid_Y'] * px_um

    # Area filtering flag
    min_px = int(params.min_area_px) if params.min_area_px is not None else 0
    max_px = int(params.max_area_px) if params.max_area_px is not None else 0
    df['IncludedByArea'] = True
    if min_px > 0:
        df.loc[df['Area_px'] < min_px, 'IncludedByArea'] = False
    if max_px > 0:
        df.loc[df['Area_px'] > max_px, 'IncludedByArea'] = False

    return df, lab

def compute_mean_z_for_labels(z_stack: np.ndarray,
                              lab_mask: np.ndarray,
                              params: Params) -> pd.DataFrame:
    """
    Calculates Mean Z and Z-Standard Deviation (index units) for each label.
    High Z-STD indicates a vertical overlap (stacking), Low Z-STD indicates a planar object (e.g., Rouleaux).
    """
    if z_stack.ndim != 3:
        raise ValueError("Z-stack must be 3D array (Z, Y, X).")
    
    # Map max intensity Z-index
    z_idx_map = np.argmax(z_stack, axis=0)

    recs = []
    for region in measure.regionprops(lab_mask):
        coords = region.coords
        if coords.size == 0:
            continue
        
        # Get Z indices for all pixels in this label
        z_vals = z_idx_map[coords[:, 0], coords[:, 1]]
        
        z_mean_idx = np.mean(z_vals)
        z_std_idx = np.std(z_vals)
        
        mean_z_um = float(z_mean_idx) * float(params.z_scale_um_per_index)
        
        recs.append({
            'Label': int(region.label), 
            'Mean_Z_um': mean_z_um,
            'Z_Std_Index': float(z_std_idx)
        })
        
    return pd.DataFrame(recs)

def skeleton_angles_near_nuclei(nuc_df: pd.DataFrame,
                                skel_img: np.ndarray,
                                params: Params) -> pd.Series:
    skel_pts = np.column_stack(np.where(skel_img > 0))
    if skel_pts.shape[0] == 0 or nuc_df.empty:
        return pd.Series(np.nan, index=nuc_df.index, name='Skeleton_Angle_deg')

    Ys = nuc_df['Centroid_Y'].values
    Xs = nuc_df['Centroid_X'].values
    out = np.full((len(nuc_df),), np.nan, dtype=float)
    R = float(params.skeleton_radius_px)

    for i, (cx, cy) in enumerate(zip(Xs, Ys)):
        d = np.sqrt((skel_pts[:, 1] - cx) ** 2 + (skel_pts[:, 0] - cy) ** 2)
        idx = np.where(d < R)[0]
        if idx.size < 5:
            k = min(10, skel_pts.shape[0])
            idx = np.argsort(d)[:k]
        ang = _principal_angle_from_local_points(
            np.column_stack([skel_pts[idx, 1], skel_pts[idx, 0]]).astype(float)
        )
        out[i] = ang
    return pd.Series(out, index=nuc_df.index, name='Skeleton_Angle_deg')

def dbscan_3d_labels(nuc_df: pd.DataFrame, params: Params) -> pd.Series:
    x = nuc_df['Centroid_X_um'].astype(float)
    y = nuc_df['Centroid_Y_um'].astype(float)
    z = nuc_df['Mean_Z_um'].astype(float)
    valid = x.notna() & y.notna() & z.notna()
    labels = pd.Series(np.nan, index=nuc_df.index, name='Cluster3D')
    if valid.sum() == 0:
        return labels

    coords = np.column_stack([x[valid].values, y[valid].values, z[valid].values])
    db = DBSCAN(eps=float(params.dbscan_eps_um), min_samples=int(params.dbscan_min_samples))
    lab = db.fit_predict(coords)
    labels.loc[valid] = lab.astype(float)
    return labels

def compute_neighbor_features(nuc_df: pd.DataFrame,
                              ks: Tuple[int, ...] = (3, 5)) -> pd.DataFrame:
    """
    Compute k-nearest-neighbor features for each nucleus using 3D (X_um, Y_um, Mean_Z_um).
    Returns a dataframe indexed like nuc_df with additional columns:

      - NN{k}_MeanDist_um
      - NN{k}_MeanAspectRatio
      - NN{k}_Frac_Spherical
      - NN{k}_Frac_Ellipsoid
      - NN{k}_Frac_Intermediate
      - NN{k}_MeanRelAngle_deg
      - NN{k}_Frac_Parallel
      - NN{k}_Frac_Perpendicular
      - NN{k}_Frac_IM_Orientation

    where 'Parallel', 'Perpendicular', 'IM_Orientation' are based on Relative_Angle_deg bins.
    """
    # Require 3D coords and core nuclear features
    coord_cols = ['Centroid_X_um', 'Centroid_Y_um', 'Mean_Z_um']
    for c in coord_cols:
        if c not in nuc_df.columns:
            return pd.DataFrame(index=nuc_df.index)

    coords = nuc_df[coord_cols].astype(float)
    valid_mask = coords.notna().all(axis=1)

    if valid_mask.sum() < 2:
        # Not enough nuclei to define neighbors
        return pd.DataFrame(index=nuc_df.index)

    coords_valid = coords[valid_mask].to_numpy()
    indices_valid = coords[valid_mask].index.to_numpy()

    # Pairwise distances (µm)
    dist_mat = distance.cdist(coords_valid, coords_valid)
    np.fill_diagonal(dist_mat, np.inf)

    # Prepare output
    out = pd.DataFrame(index=nuc_df.index, dtype=float)

    # Helper views for neighbor statistics
    shape_series = nuc_df['Shape_Class'] if 'Shape_Class' in nuc_df.columns else pd.Series(index=nuc_df.index, dtype=object)
    ar_series = nuc_df['Aspect_Ratio'] if 'Aspect_Ratio' in nuc_df.columns else pd.Series(index=nuc_df.index, dtype=float)
    rel_angle = nuc_df['Relative_Angle_deg'] if 'Relative_Angle_deg' in nuc_df.columns else pd.Series(index=nuc_df.index, dtype=float)

    # Orientation class for each nucleus (same bins as summarize_fiber)
    parallel_mask = (rel_angle >= 0) & (rel_angle <= 10)
    perp_mask = (rel_angle >= 80) & (rel_angle <= 90)
    im_mask = (rel_angle > 10) & (rel_angle < 80)

    for row_i, idx in enumerate(indices_valid):
        d_row = dist_mat[row_i]

        # Skip pathological rows
        if not np.isfinite(d_row).any():
            continue

        # Sort other nuclei by distance
        order = np.argsort(d_row)

        for k in ks:
            if k > len(order):
                continue  # not enough neighbors for this k

            neigh_pos = order[:k]
            neigh_indices = indices_valid[neigh_pos]

            # Distances
            d_k = d_row[neigh_pos]
            out.loc[idx, f'NN{k}_MeanDist_um'] = float(np.nanmean(d_k)) if d_k.size else np.nan

            # Aspect ratio
            ar_k = ar_series.loc[neigh_indices].astype(float)
            out.loc[idx, f'NN{k}_MeanAspectRatio'] = float(ar_k.mean()) if not ar_k.empty else np.nan

            # Shape class composition
            shapes_k = shape_series.loc[neigh_indices].astype(str)
            n_k = len(shapes_k)
            if n_k > 0:
                for cls in ['Spherical', 'Ellipsoid', 'Intermediate']:
                    frac = (shapes_k == cls).sum() / n_k
                    out.loc[idx, f'NN{k}_Frac_{cls}'] = float(frac)
            else:
                for cls in ['Spherical', 'Ellipsoid', 'Intermediate']:
                    out.loc[idx, f'NN{k}_Frac_{cls}'] = np.nan

            # Orientation-relative-to-skeleton distribution in neighbors
            par_k = parallel_mask.loc[neigh_indices]
            perp_k = perp_mask.loc[neigh_indices]
            im_k = im_mask.loc[neigh_indices]
            n_k_orient = float(n_k)

            if n_k_orient > 0:
                out.loc[idx, f'NN{k}_Frac_Parallel'] = float(par_k.sum() / n_k_orient)
                out.loc[idx, f'NN{k}_Frac_Perpendicular'] = float(perp_k.sum() / n_k_orient)
                out.loc[idx, f'NN{k}_Frac_IM_Orientation'] = float(im_k.sum() / n_k_orient)
            else:
                out.loc[idx, f'NN{k}_Frac_Parallel'] = np.nan
                out.loc[idx, f'NN{k}_Frac_Perpendicular'] = np.nan
                out.loc[idx, f'NN{k}_Frac_IM_Orientation'] = np.nan

            # Mean Relative_Angle_deg of neighbors
            rel_k = rel_angle.loc[neigh_indices].astype(float)
            out.loc[idx, f'NN{k}_MeanRelAngle_deg'] = float(rel_k.mean()) if not rel_k.empty else np.nan

    return out

def compute_fiber_width_profile(nuc_df: pd.DataFrame,
                                skel_img: np.ndarray,
                                params: Params) -> pd.DataFrame:
    """Measure fiber diameter at regular intervals (e.g., every 100 µm) along the skeleton."""
    if nuc_df.empty or 'Centroid_Y' not in nuc_df.columns or 'Centroid_X' not in nuc_df.columns:
        return pd.DataFrame()

    centroids = nuc_df[['Centroid_Y', 'Centroid_X']].to_numpy().astype(float)
    valid = np.isfinite(centroids).all(axis=1)
    centroids = centroids[valid]
    if centroids.size == 0:
        return pd.DataFrame()

    skel = _to_bool(skel_img)
    if skel.sum() == 0:
        return pd.DataFrame()

    endpoints = _find_skeleton_endpoints(skel)
    if len(endpoints) < 2:
        return pd.DataFrame()

    if len(endpoints) > 2:
        endpoints_arr = np.array(endpoints)
        dmat = distance.squareform(distance.pdist(endpoints_arr))
        i0, i1 = np.unravel_index(np.argmax(dmat), dmat.shape)
        start = tuple(endpoints_arr[i0])
        end = tuple(endpoints_arr[i1])
    else:
        start, end = endpoints[0], endpoints[1]

    dist_map = distance_transform_edt(~skel)
    try:
        path, _ = route_through_array(dist_map, start, end, fully_connected=True)
    except Exception:
        return pd.DataFrame()

    path = np.asarray(path, dtype=float)  # (N, 2) rows = (y, x)
    if path.ndim != 2 or path.shape[0] < 2:
        return pd.DataFrame()

    diffs = np.diff(path, axis=0)
    step_lengths_px = np.linalg.norm(diffs, axis=1)
    if step_lengths_px.size == 0:
        return pd.DataFrame()

    cumdist_px = np.concatenate([[0.0], np.cumsum(step_lengths_px)])
    px_to_um = float(params.pixel_size_xy_um)
    cumdist_um = cumdist_px * px_to_um
    total_length_um = cumdist_um[-1]

    step_um = float(getattr(params, "fiber_width_step_um", 100.0))
    if step_um <= 0:
        step_um = 100.0

    target_positions_um = np.arange(0.0, total_length_um + step_um * 0.5, step_um)
    if target_positions_um.size == 0:
        return pd.DataFrame()

    sample_indices = np.searchsorted(cumdist_um, target_positions_um)
    sample_indices = np.clip(sample_indices, 0, len(path) - 1)
    sample_indices = np.unique(sample_indices)

    max_radius_um = float(getattr(params, "fiber_width_max_radius_um", 100.0))
    max_radius_px = max(1, int(round(max_radius_um / px_to_um)))

    results = []

    for pos_idx in sample_indices:
        center = path[pos_idx]  # (y, x)
        y, x = center.astype(float)

        k = 5
        i0 = max(0, pos_idx - k)
        i1 = min(len(path) - 1, pos_idx + k)
        prev_pt = path[i0].astype(float)
        next_pt = path[i1].astype(float)

        dy, dx = next_pt[0] - prev_pt[0], next_pt[1] - prev_pt[1]
        norm = np.hypot(dy, dx)
        if norm == 0:
            continue

        dir_vec = np.array([dy, dx]) / norm
        perp_vec = np.array([-dx, dy]) / norm  # +90° perpendicular

        disp = centroids - np.array([[y, x]])  # (N, 2), (Δy, Δx)
        dists_px = np.linalg.norm(disp, axis=1)
        within_radius = dists_px <= max_radius_px
        if not within_radius.any():
            continue

        disp_local = disp[within_radius]
        along = disp_local @ dir_vec
        perp_proj = disp_local @ perp_vec

        local_mask = np.abs(along) <= max_radius_px
        if not local_mask.any():
            continue

        perp_proj = perp_proj[local_mask]
        disp_local = disp_local[local_mask]

        left_mask = perp_proj < 0
        right_mask = perp_proj > 0
        if not (left_mask.any() and right_mask.any()):
            continue

        left_indices = np.where(left_mask)[0]
        right_indices = np.where(right_mask)[0]

        left_idx_local = left_indices[np.argmin(perp_proj[left_indices])]    # most negative
        right_idx_local = right_indices[np.argmax(perp_proj[right_indices])]  # most positive

        left_centroid = np.array([y, x]) + disp_local[left_idx_local]
        right_centroid = np.array([y, x]) + disp_local[right_idx_local]

        diameter_px = np.linalg.norm(left_centroid - right_centroid)
        diameter_um = diameter_px * px_to_um
        position_um = cumdist_um[pos_idx]

        results.append({
            "Position_um": position_um,
            "Diameter_um": diameter_um,
            "Skel_Y": y, "Skel_X": x,
            "Left_Y": left_centroid[0], "Left_X": left_centroid[1],
            "Right_Y": right_centroid[0], "Right_X": right_centroid[1],
        })

    if not results:
        return pd.DataFrame()

    return pd.DataFrame(results)

def save_labeled_overlay(stdip_mask: np.ndarray,
                         lab_mask: np.ndarray,
                         nuc_df: pd.DataFrame,
                         out_path: str):
    """Create RGB overlay with STDIP background, red outlines, and white label IDs."""
    bg = (stdip_mask.astype(np.uint8) * 255)
    if bg.ndim != 2:
        bg = bg.squeeze()
        if bg.ndim != 2:
            bg = bg[..., 0]
    rgb = np.dstack([bg, bg, bg]).astype(np.uint8)

    bounds = find_boundaries(lab_mask, mode='outer')
    rgb[bounds, 0] = 255; rgb[bounds, 1] = 0; rgb[bounds, 2] = 0

    try:
        img = Image.fromarray(rgb)
        draw = ImageDraw.Draw(img)
        try:
            font = ImageFont.truetype("arial.ttf", size=12)
        except Exception:
            font = ImageFont.load_default()
        for _, row in nuc_df.iterrows():
            x = float(row.get('Centroid_X', np.nan))
            y = float(row.get('Centroid_Y', np.nan))
            lab_id = int(row.get('Label', -1))
            if not np.isnan(x) and not np.isnan(y) and lab_id >= 0:
                draw.text((x + 1, y + 1), str(lab_id), fill=(255, 255, 255),
                          stroke_width=1, stroke_fill=(0, 0, 0), font=font)
        img.save(out_path)
    except Exception:
        io.imsave(out_path, rgb)

# ---------------------------- File matching ----------------------------------
def parse_file_metadata(path: str) -> Optional[Tuple[str, str, str, str]]:
    name = os.path.basename(path)
    m = FILENAME_REGEX.match(name)
    if not m:
        return None
    return (m.group('Subject'), m.group('Day'), m.group('Side'), m.group('Fiber'))

def _gather_files(folder: str, exts=('.tif', '.tiff')) -> List[str]:
    if not os.path.isdir(folder):
        return []
    return sorted([os.path.join(folder, f) for f in os.listdir(folder)
                   if f.lower().endswith(exts)])

def _build_fiber_map(leg_dir: str) -> Dict[Tuple[str, str, str, str], Tuple[str, str, str]]:
    stdip = _gather_files(os.path.join(leg_dir, 'STDIP'))
    skel = _gather_files(os.path.join(leg_dir, 'Skel'))
    zstk = _gather_files(os.path.join(leg_dir, 'TIFs'))

    maps = {}
    for group, files in [('STDIP', stdip), ('Skel', skel), ('TIFs', zstk)]:
        for p in files:
            meta = parse_file_metadata(p)
            if not meta:
                continue
            maps.setdefault((group, meta), []).append(p)

    meta_keys = {meta for (grp, meta) in maps.keys()}
    fiber_map = {}
    for meta in sorted(meta_keys):
        s_list = maps.get(('STDIP', meta), [])
        k_list = maps.get(('Skel', meta), [])
        z_list = maps.get(('TIFs', meta), [])
        s_choice = s_list[0] if s_list else None
        k_choice = k_list[0] if k_list else None
        z_choice = z_list[0] if z_list else None
        if s_choice and k_choice and z_choice:
            fiber_map[meta] = (s_choice, k_choice, z_choice)
    return fiber_map

# ----------------------------- Distance to skeleton --------------------------
def compute_distance_to_skeleton(nuc_df: pd.DataFrame,
                                 skel_img: np.ndarray,
                                 params: Params) -> pd.Series:
    """Compute distance (µm) from each nucleus centroid to nearest skeleton pixel."""
    skel_bool = _to_bool(skel_img)
    dist_map = distance_transform_edt(~skel_bool)
    px_um = float(params.pixel_size_xy_um)

    distances = []
    for _, row in nuc_df.iterrows():
        y = int(round(row['Centroid_Y']))
        x = int(round(row['Centroid_X']))
        if 0 <= y < dist_map.shape[0] and 0 <= x < dist_map.shape[1]:
            d_px = dist_map[y, x]
            distances.append(d_px * px_um)
        else:
            distances.append(np.nan)
    return pd.Series(distances, index=nuc_df.index, name='DistanceToSkel')

# ----------------------------- Biopsy summary row ----------------------------
def summarize_fiber(nuc_df: pd.DataFrame,
                    width_df: pd.DataFrame,
                    imaris_row: Optional[pd.Series],
                    excluded_df: Optional[pd.DataFrame] = None) -> Dict:
    """
    Summarize statistics for one fiber.
    nuc_df: The INCLUDED nuclei.
    excluded_df: The EXCLUDED nuclei (used for rejection counts).
    """
    summary = {}
    summary['Subject'] = nuc_df['Subject'].iloc[0]
    summary['Timepoint'] = nuc_df['Timepoint'].iloc[0]
    summary['Side'] = nuc_df['Side'].iloc[0]
    summary['Fiber'] = nuc_df['Fiber'].iloc[0]
    summary['FiberType'] = nuc_df['FiberType'].iloc[0]

    # Imaris metrics
    if imaris_row is not None:
        summary['ImarisMyonuclei'] = imaris_row.get('Imaris Myonuclei', np.nan)
        summary['ImarisFiberLength'] = imaris_row.get('ImarisFiberLength', np.nan)
        summary['ImarisMyonuclei/mm'] = imaris_row.get('Nuclei/mm', np.nan)
    else:
        summary['ImarisMyonuclei'] = np.nan
        summary['ImarisFiberLength'] = np.nan
        summary['ImarisMyonuclei/mm'] = np.nan

    # In this logic, 'nuc_df' already contains only the included ones
    included = nuc_df
    summary['PythonMyonuclei'] = int(len(included))

    # Excluded Counts
    if excluded_df is not None and not excluded_df.empty:
        summary['Excluded_Total'] = len(excluded_df)
        if 'ExclusionReason' in excluded_df.columns:
            summary['Excluded_By_Area'] = int((excluded_df['ExclusionReason'] == 'Area').sum())
            summary['Excluded_By_Z'] = int((excluded_df['ExclusionReason'] == 'Z-Depth').sum())
        else:
            summary['Excluded_By_Area'] = np.nan
            summary['Excluded_By_Z'] = np.nan
    else:
        summary['Excluded_Total'] = 0
        summary['Excluded_By_Area'] = 0
        summary['Excluded_By_Z'] = 0

    # Diameter stats
    if not width_df.empty:
        summary['MeanDiameter'] = float(width_df['Diameter_um'].mean())
        summary['MinDiameter'] = float(width_df['Diameter_um'].min())
        summary['MaxDiameter'] = float(width_df['Diameter_um'].max())
        summary['DiameterSD'] = float(width_df['Diameter_um'].std())
    else:
        summary['MeanDiameter'] = summary['MinDiameter'] = summary['MaxDiameter'] = summary['DiameterSD'] = np.nan

    # Fiber volume (µm^3) cylinder approximation: pi * r^2 * h
    if (not np.isnan(summary['MeanDiameter'])) and (not np.isnan(summary['ImarisFiberLength'])):
        r = summary['MeanDiameter'] / 2.0
        h = summary['ImarisFiberLength']
        summary['FiberVolume'] = float(np.pi * (r ** 2) * h)
    else:
        summary['FiberVolume'] = np.nan

    # Myonuclear domain volume (µm^3 per nucleus)
    # Using Python myonuclei count (segmented nuclei)
    if summary['PythonMyonuclei'] > 0 and not np.isnan(summary['FiberVolume']):
        summary['DomainVolume_um3_Python'] = float(summary['FiberVolume'] / summary['PythonMyonuclei'])
    else:
        summary['DomainVolume_um3_Python'] = np.nan

    # Optional: same concept using Imaris myonuclei count
    if (not np.isnan(summary['FiberVolume'])) and (not np.isnan(summary['ImarisMyonuclei'])) and summary['ImarisMyonuclei'] > 0:
        summary['DomainVolume_um3_Imaris'] = float(summary['FiberVolume'] / summary['ImarisMyonuclei'])
    else:
        summary['DomainVolume_um3_Imaris'] = np.nan

    # --- Neighbor shape and orientation context (fiber-level means of NN3) ---

    # Mean NN3 distance (how tightly packed neighbors are)
    if 'NN3_MeanDist_um' in included.columns:
        summary['Mean_NN3_MeanDist_um'] = float(included['NN3_MeanDist_um'].mean())
        summary['SD_NN3_MeanDist_um'] = float(included['NN3_MeanDist_um'].std())
    else:
        summary['Mean_NN3_MeanDist_um'] = np.nan
        summary['SD_NN3_MeanDist_um'] = np.nan

    # Shape context of neighbors (fractions among nearest neighbors)
    for col in ['NN3_Frac_Spherical', 'NN3_Frac_Ellipsoid', 'NN3_Frac_Intermediate']:
        if col in included.columns:
            summary[f'Mean_{col}'] = float(included[col].mean())
        else:
            summary[f'Mean_{col}'] = np.nan

    # Orientation context of neighbors (fraction parallel / perpendicular / intermediate)
    for col in ['NN3_Frac_Parallel', 'NN3_Frac_Perpendicular', 'NN3_Frac_IM_Orientation']:
        if col in included.columns:
            summary[f'Mean_{col}'] = float(included[col].mean())
        else:
            summary[f'Mean_{col}'] = np.nan


    # Shape classes
    for cls in ['Ellipsoid', 'Intermediate', 'Spherical']:
        count = int((included['Shape_Class'] == cls).sum())
        summary[f'#_{cls}'] = count
        summary[f'%_{cls}'] = (count / summary['PythonMyonuclei'] * 100) if summary['PythonMyonuclei'] else np.nan

    # Orientation bins
    rel = included['Relative_Angle_deg']
    parallel = int(((rel >= 0) & (rel <= 10)).sum())
    perpendicular = int(((rel >= 80) & (rel <= 90)).sum())
    intermediate = int(((rel > 10) & (rel < 80)).sum())

    summary['#_Parallel'] = parallel
    summary['%_Parallel'] = (parallel / summary['PythonMyonuclei'] * 100) if summary['PythonMyonuclei'] else np.nan
    summary['#_Perpendicular'] = perpendicular
    summary['%_Perpendicular'] = (perpendicular / summary['PythonMyonuclei'] * 100) if summary['PythonMyonuclei'] else np.nan
    summary['#_IM_Orientation'] = intermediate
    summary['%_IM_Orientation'] = (intermediate / summary['PythonMyonuclei'] * 100) if summary['PythonMyonuclei'] else np.nan

    # Clusters (exclude noise label -1)
    if 'Cluster3D' in included.columns:
        valid_clusters = included['Cluster3D'].dropna()
        summary['Cluster#'] = int((valid_clusters[valid_clusters != -1]).nunique())
    else:
        summary['Cluster#'] = np.nan

    # Rouleaux and Central Rouleaux
    rouleaux_mask = included['Aspect_Ratio'] > 4
    summary['Rouleaux_#'] = int(rouleaux_mask.sum())

    central_mask = (included['Aspect_Ratio'] > 4) & (included['DistanceToSkel'] <= 15.0)
    summary['CentralRouleaux_#'] = int(central_mask.sum())

    return summary

# ----------------------------- Process fiber ---------------------------------
def process_fiber(meta_key: Tuple[str, str, str, str],
                  stdip_path: str,
                  skel_path: str,
                  zstack_path: str,
                  out_dir: str,
                  params: Params) -> pd.DataFrame:
    subject, day, side, fiber_code = meta_key
    fiber_num, fiber_type_num = _parse_fiber_num_and_type(fiber_code)
    fiber_col_value = fiber_num if fiber_num is not None else fiber_code
    fiber_type_col_value = fiber_type_num

    base = f"{subject}_{day}_{side}_Fiber{fiber_code}"
    _ensure_dir(out_dir)

    # Read inputs
    mask = _to_bool(_imread(stdip_path))
    skel = _imread(skel_path)
    zstk = _imread(zstack_path)

    # Measurements & Z calculations
    nuc_df, lab = measure_nuclei_from_binary(mask, params)
    z_df = compute_mean_z_for_labels(zstk, lab, params)
    
    # Merge Z info (Mean Z, Z_Std_Index)
    nuc_df = nuc_df.merge(z_df, on='Label', how='left')

    # --- FILTERING LOGIC ---
    # 1. Z-Consistency Flag: Exclude if vertical spread (std dev) is too high
    nuc_df['IncludedByZ'] = nuc_df['Z_Std_Index'] <= params.z_std_threshold
    
    # 2. Combined Flag: Must pass both Area check and Z check
    # Note: 'IncludedByArea' was computed in measure_nuclei_from_binary
    nuc_df['Included'] = nuc_df['IncludedByArea'] & nuc_df['IncludedByZ']
    # -----------------------

    # Angles
    nuc_df['Skeleton_Angle_deg'] = skeleton_angles_near_nuclei(nuc_df, skel, params)
    nuc_df['Relative_Angle_deg'] = [
        _relative_angle(a, b) for a, b in zip(nuc_df['Orientation_deg'], nuc_df['Skeleton_Angle_deg'])
    ]

    # Distances to skeleton (µm)
    nuc_df['DistanceToSkel'] = compute_distance_to_skeleton(nuc_df, skel, params)

    # Clustering (only meaningful for those with valid Z, but run on all for completeness)
    nuc_df['Cluster3D'] = dbscan_3d_labels(nuc_df, params)

    # Nearest-neighbor features (compute on INCLUDED nuclei only; fill NaN for excluded)
    try:
        nn_features = compute_neighbor_features(nuc_df[nuc_df['Included']])
        for col in nn_features.columns:
            nuc_df[col] = np.nan  # initialize
        nuc_df.loc[nn_features.index, nn_features.columns] = nn_features
    except Exception as e:
        warnings.warn(f"    Failed NN computation for {base}: {e}")


    # Metadata
    nuc_df.insert(0, 'Subject', subject)
    nuc_df.insert(1, 'Timepoint', day)
    nuc_df.insert(2, 'Side', side)
    nuc_df.insert(3, 'Fiber', fiber_col_value)
    nuc_df.insert(4, 'FiberType', fiber_type_col_value)

    # Split included/excluded based on the master 'Included' flag
    included_df = nuc_df[nuc_df['Included']].copy()
    excluded_df = nuc_df[~nuc_df['Included']].copy()

    # Determine exclusion reason for the excluded dataframe
    if not excluded_df.empty:
        # Priority: Check Area first, then Z. 
        # If it failed both, it will be labeled as Area (can adjust priority order if needed)
        conditions = [
            (~excluded_df['IncludedByArea']),
            (~excluded_df['IncludedByZ'])
        ]
        choices = ['Area', 'Z-Depth']
        excluded_df['ExclusionReason'] = np.select(conditions, choices, default='Unknown')

    # Save CSVs
    out_csv = os.path.join(out_dir, f"{base}_nuclei_results.csv")
    included_df.to_csv(out_csv, index=False)

    out_csv_ex = os.path.join(out_dir, f"{base}_excluded_nuclei.csv")
    excluded_df.to_csv(out_csv_ex, index=False)

    # Fiber-width profile
    width_df = compute_fiber_width_profile(included_df, skel, params)
    if not width_df.empty:
        width_df.insert(0, 'Subject', subject)
        width_df.insert(1, 'Timepoint', day)
        width_df.insert(2, 'Side', side)
        width_df.insert(3, 'Fiber', fiber_col_value)
        width_df.insert(4, 'FiberType', fiber_type_col_value)
        width_csv = os.path.join(out_dir, f"{base}_fiber_width_profile.csv")
        width_df.to_csv(width_csv, index=False)

    # Overlay (Draw included only, or all? Usually draw all but maybe different colors? 
    # Current function draws all labels in red/white. Keeping as is.)
    overlay_path = os.path.join(out_dir, f"{base}_overlay.png")
    save_labeled_overlay(mask, lab, nuc_df, overlay_path)

    return included_df, excluded_df  # return both for summary

# ----------------------------- Imaris master ---------------------------------
def load_imaris_master(master_path: str) -> pd.DataFrame:
    # Supports .xlsx or .csv; infer by extension
    if master_path.lower().endswith(('.xlsx', '.xls')):
        df = pd.read_excel(master_path)
    else:
        df = pd.read_csv(master_path)
    # Normalize expected column names
    df = df.rename(columns={
        'Subject': 'Subject',
        'Timepoint': 'Timepoint',
        'Leg': 'Side',  # Map Imaris 'Leg' to pipeline 'Side' (R/L)
        'Fiber#': 'Fiber',
        'FiberType': 'FiberType',
        'Imaris Myonuclei': 'Imaris Myonuclei',
        'ImarisFiberLength': 'ImarisFiberLength',
        'Nuclei/mm': 'Nuclei/mm'
    })
    # Ensure numeric Fiber and FiberType if stored as text
    df['Fiber'] = pd.to_numeric(df['Fiber'], errors='coerce')
    df['FiberType'] = pd.to_numeric(df['FiberType'], errors='coerce')
    return df

def find_imaris_row(imaris_df: pd.DataFrame,
                    subject: str, timepoint: str, side: str,
                    fiber: Optional[int], fibertype: Optional[int]) -> Optional[pd.Series]:
    q = (imaris_df['Subject'].astype(str) == str(subject)) & \
        (imaris_df['Timepoint'].astype(str) == str(timepoint)) & \
        (imaris_df['Side'].astype(str) == str(side))
    if fiber is not None:
        q = q & (imaris_df['Fiber'] == fiber)
    if fibertype is not None and 'FiberType' in imaris_df.columns:
        q = q & (imaris_df['FiberType'] == fibertype)
    match = imaris_df[q]
    if match.empty:
        return None
    return match.iloc[0]

# ----------------------------- Batch orchestration ---------------------------
def batch_subject(base_dir: str, master_path: str, params: Params):
    # Process each subject directory
    for subject in sorted(os.listdir(base_dir)):
        subj_dir = os.path.join(base_dir, subject)
        if not os.path.isdir(subj_dir):
            continue

        # Load Imaris master once per subject
        imaris_df = load_imaris_master(master_path)

        timepoint_to_rows: Dict[str, List[pd.DataFrame]] = {}

        for timepoint in sorted(os.listdir(subj_dir)):
            tp_dir = os.path.join(subj_dir, timepoint)
            if not os.path.isdir(tp_dir):
                continue

            for leg in sorted(os.listdir(tp_dir)):
                leg_dir = os.path.join(tp_dir, leg)
                if not os.path.isdir(leg_dir):
                    continue

                print(f"\n== {subject}/{timepoint}/{leg} ==")
                fiber_map = _build_fiber_map(leg_dir)
                if not fiber_map:
                    all_candidates = (
                        _gather_files(os.path.join(leg_dir, 'STDIP')) +
                        _gather_files(os.path.join(leg_dir, 'Skel')) +
                        _gather_files(os.path.join(leg_dir, 'TIFs'))
                    )
                    print("  Debug: parsed metadata from files in this leg:")
                    for p in all_candidates:
                        m = FILENAME_REGEX.match(os.path.basename(p))
                        print("   -", os.path.basename(p), "->", (m.groupdict() if m else "NO MATCH"))

                # Collect biopsy-level summaries and excluded fibers
                biopsy_summaries: List[Dict] = []
                excluded_fibers: List[Dict] = []

                for meta_key, (stdip, skel, zstk) in fiber_map.items():
                    m_subject, m_day, m_side, m_fiber_code = meta_key
                    if m_subject != subject or m_day != timepoint or m_side != leg:
                        warnings.warn(f"Metadata mismatch: parsed {meta_key} but folder is {subject}/{timepoint}/{leg}.")

                    fiber_num, fiber_type_num = _parse_fiber_num_and_type(m_fiber_code)
                    fiber_tag = f"{m_subject}_{m_day}_{m_side}_Fiber{m_fiber_code}"
                    out_dir = os.path.join(leg_dir, f"{fiber_tag}_output")

                    try:
                        included_df, excluded_df = process_fiber(meta_key, stdip, skel, zstk, out_dir, params)
                        timepoint_to_rows.setdefault(m_day, []).append(included_df.copy())

                        # Load width profile if exists
                        width_csv = os.path.join(out_dir, f"{fiber_tag}_fiber_width_profile.csv")
                        width_df = pd.read_csv(width_csv) if os.path.exists(width_csv) else pd.DataFrame()

                        # Attempt to find Imaris row
                        im_row = find_imaris_row(imaris_df, m_subject, m_day, m_side, fiber_num, fiber_type_num)

                        # Summarize fiber
                        summary_row = summarize_fiber(included_df, width_df, im_row, excluded_df)

                        if im_row is None:
                            summary_row['ExcludedReason'] = 'No Imaris match'
                            excluded_fibers.append(summary_row)
                        else:
                            biopsy_summaries.append(summary_row)

                    except Exception as e:
                        warnings.warn(f"    Failed {fiber_tag}: {e}")

                # Write biopsy-level summaries for this leg
                if biopsy_summaries or excluded_fibers:
                    biopsy_df = pd.DataFrame(biopsy_summaries)
                    excl_df = pd.DataFrame(excluded_fibers)

                    # Order columns per your schema
                    desired_cols = [
                        'Subject','Timepoint','Side','Fiber','FiberType',
                        'ImarisMyonuclei','ImarisFiberLength','ImarisMyonuclei/mm',
                        'PythonMyonuclei','Excluded_Total','Excluded_By_Area','Excluded_By_Z',
                        'MeanDiameter','MinDiameter','MaxDiameter','DiameterSD',
                        'FiberVolume',
                        # Volume-based myonuclear domain (footprint)
                        'DomainVolume_um3_Python','DomainVolume_um3_Imaris',
                        # Neighbor distance context
                        'Mean_NN3_MeanDist_um','SD_NN3_MeanDist_um',
                        # Neighbor shape context
                        'Mean_NN3_Frac_Spherical','Mean_NN3_Frac_Ellipsoid','Mean_NN3_Frac_Intermediate',
                        # Neighbor orientation context
                        'Mean_NN3_Frac_Parallel','Mean_NN3_Frac_Perpendicular','Mean_NN3_Frac_IM_Orientation',
                        # Global shape/orientation counts
                        '#_Ellipsoid','%_Ellipsoid','#_Intermediate','%_Intermediate','#_Spherical','%_Spherical',
                        '#_Parallel','%_Parallel','#_IM_Orientation','%_IM_Orientation','#_Perpendicular','%_Perpendicular',
                        'Cluster#','Rouleaux_#','CentralRouleaux_#'
                    ]
                    # Filter for columns that actually exist
                    biopsy_df = biopsy_df[[c for c in desired_cols if c in biopsy_df.columns]]
                    excl_df = excl_df[[c for c in desired_cols if c in excl_df.columns] + ['ExcludedReason']]

                    out_biopsy = os.path.join(leg_dir, f"{subject}_{timepoint}_{leg}_biopsy_summary.csv")
                    biopsy_df.to_csv(out_biopsy, index=False)

                    out_excl = os.path.join(leg_dir, f"{subject}_{timepoint}_{leg}_excluded_fibers.csv")
                    excl_df.to_csv(out_excl, index=False)

                    print(f"  Biopsy summary written: {out_biopsy}")
                    if not excl_df.empty:
                        print(f"  Excluded fibers written: {out_excl}")

        # Subject workbook (included nuclei only)
        if timepoint_to_rows:
            xlsx_path = os.path.join(subj_dir, f"{subject}_results.xlsx")
            with pd.ExcelWriter(xlsx_path, engine='openpyxl') as writer:
                for tp, df_list in sorted(timepoint_to_rows.items()):
                    all_tp = pd.concat(df_list, ignore_index=True) if df_list else pd.DataFrame()

                    base_cols = [
                        'Subject','Timepoint','Side','Fiber','FiberType',
                        'Label','Centroid_X','Centroid_Y','Centroid_X_um','Centroid_Y_um',
                        'Area_px','Area_um2','Perimeter_px','Perimeter_um',
                        'MajorAxis_px','MinorAxis_px','MajorAxis_um','MinorAxis_um',
                        'Aspect_Ratio','Shape_Class','Orientation_deg','Mean_Z_um','Z_Std_Index',
                        'Skeleton_Angle_deg','Relative_Angle_deg','Cluster3D','Included','IncludedByArea','IncludedByZ','DistanceToSkel'
                    ]
                    col_order = [c for c in base_cols if c in all_tp.columns] + \
                                [c for c in all_tp.columns if c not in set(base_cols)]

                    out_df = all_tp[col_order] if col_order else all_tp
                    out_df.to_excel(writer, sheet_name=str(tp), index=False)

            print(f"Subject workbook written: {xlsx_path}")
        else:
            print(f"No results aggregated for subject {subject}.")

# ----------------------------- CLI entrypoint --------------------------------
def main():
    import argparse
    p = argparse.ArgumentParser(description="CellProfiler-free single-fiber nuclei analysis with biopsy-level summaries.")
    p.add_argument("base_dir", help="Base directory containing subject/timepoint/leg folders.")
    p.add_argument("--imaris_master", required=True, help="Path to Imaris master .xlsx or .csv.")
    p.add_argument("--pixel_size_xy_um", type=float, default=1.0, help="Microns per pixel in XY.")
    p.add_argument("--z_scale_um_per_index", type=float, default=2.7, help="Microns per Z index step.")
    p.add_argument("--z_std_threshold", type=float, default=2.0, help="Max std of per-pixel Z indices. > this = vertical overlap (Excluded).")
    p.add_argument("--dbscan_eps_um", type=float, default=20.0, help="DBSCAN eps in microns (after XY conversion).")
    p.add_argument("--dbscan_min_samples", type=int, default=2, help="DBSCAN min_samples.")
    p.add_argument("--skeleton_radius_px", type=float, default=20.0, help="Neighborhood radius (px) for skeleton PCA.")
    p.add_argument("--min_area_px", type=int, default=0, help="Minimum area filter (px^2). 0 disables.")
    p.add_argument("--max_area_px", type=int, default=0, help="Maximum area filter (px^2). 0 disables.")
    p.add_argument("--fiber_width_step_um", type=float, default=100.0, help="Step size along fiber for diameter sampling (µm).")
    p.add_argument("--fiber_width_max_radius_um", type=float, default=100.0, help="Max search radius from skeleton for nuclei (µm).")

    args = p.parse_args()
    params = Params(
        pixel_size_xy_um=args.pixel_size_xy_um,
        z_scale_um_per_index=args.z_scale_um_per_index,
        z_std_threshold=args.z_std_threshold,
        dbscan_eps_um=args.dbscan_eps_um,
        dbscan_min_samples=args.dbscan_min_samples,
        skeleton_radius_px=args.skeleton_radius_px,
        min_area_px=args.min_area_px,
        max_area_px=args.max_area_px,
        fiber_width_step_um=args.fiber_width_step_um,
        fiber_width_max_radius_um=args.fiber_width_max_radius_um,
    )
    batch_subject(args.base_dir, args.imaris_master, params)

if __name__ == "__main__":
    main()