import os
import cv2
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from skimage import io, measure
from scipy.ndimage import label

# === User Inputs ===
nuclei_mask_path = r"SDTIP_path"
skeleton_mask_path = r"skeleton_path"
output_dir = r"output_dir"

os.makedirs(output_dir, exist_ok=True)

# === Load and binarize nuclei mask ===
nuclei_mask = io.imread(nuclei_mask_path)
nuclei_bin = (nuclei_mask > 0).astype(np.uint8)
if nuclei_bin.ndim > 2:
    nuclei_bin = np.any(nuclei_bin, axis=2).astype(np.uint8)

# === Load and binarize skeleton mask ===
skeleton_mask = io.imread(skeleton_mask_path)
skel_bin = (skeleton_mask > 0).astype(np.uint8)
if skel_bin.ndim > 2:
    skel_bin = np.any(skel_bin, axis=2).astype(np.uint8)

# Skeleton coordinates (row, col)
skel_coords = np.column_stack(np.where(skel_bin))

# === Measure orientations and relative angle ===
labeled_nuclei, _ = label(nuclei_bin)

results = []
for region in measure.regionprops(labeled_nuclei):
    if region.area < 5:
        continue

    # Single nucleus mask
    mask_single = (labeled_nuclei == region.label).astype(np.uint8)
    contours, _ = cv2.findContours(mask_single, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if not contours or len(contours[0]) < 5:
        continue

    (cx, cy), (w, h), nuc_ang = cv2.fitEllipse(contours[0])

    # Make sure w is major axis
    if w < h:
        w, h = h, w
        nuc_ang = (nuc_ang + 90) % 180

    # Local skeleton orientation via PCA
    dists = np.sqrt((skel_coords[:, 0] - int(cy))**2 +
                    (skel_coords[:, 1] - int(cx))**2)
    radius = 20
    idx = np.where(dists < radius)[0]
    if len(idx) < 5:
        idx = np.argsort(dists)[:10]

    local_pts = skel_coords[idx]
    pts = np.column_stack((local_pts[:, 1], local_pts[:, 0])).astype(float)
    mean_pts = pts.mean(axis=0)
    cov = np.cov((pts - mean_pts).T)
    eigvals, eigvecs = np.linalg.eig(cov)
    principal = eigvecs[:, np.argmax(eigvals)]
    skel_ang = np.degrees(np.arctan2(principal[1], principal[0])) % 180

    # Relative angle 0–90°
    rel_ang = abs(nuc_ang - skel_ang)
    if rel_ang > 90:
        rel_ang = 180 - rel_ang

    results.append({
        "Centroid_X": cx,
        "Centroid_Y": cy,
        "Nucleus_Angle": nuc_ang,
        "Skeleton_Angle": skel_ang,
        "Relative_Angle": rel_ang
    })

df = pd.DataFrame(results)
if df.empty:
    raise RuntimeError("No valid nuclei found for orientation plotting.")

# === Build overlay image ===
overlay = cv2.cvtColor((nuclei_bin * 255).astype(np.uint8), cv2.COLOR_GRAY2BGR)
for (y, x) in skel_coords:
    overlay[int(y), int(x)] = (255, 0, 0)  # blue skeleton

# === Plot overlay + orientation arrows ===
fig, ax = plt.subplots(figsize=(8, 8))
ax.imshow(overlay)

for _, r in df.iterrows():
    x0, y0 = r['Centroid_X'], r['Centroid_Y']
    ang_rad = np.deg2rad(r['Nucleus_Angle'])
    dx, dy = np.cos(ang_rad), np.sin(ang_rad)
    length = 15

    cmap = plt.get_cmap('cool')  # magenta→cyan
    color = cmap(r['Relative_Angle'] / 90.0)

    ax.arrow(x0, y0, length*dx, length*dy,
             color=color, linewidth=3,
             head_width=2.5, head_length=3)

ax.axis('off')

# Colorbar for relative angle
sm = plt.cm.ScalarMappable(cmap='cool', norm=plt.Normalize(0, 90))
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, orientation='horizontal', pad=0.02)
cbar.set_label('Relative Angle (°)', fontsize=12)
cbar.ax.tick_params(labelsize=12)

# Save + show
output_fig = os.path.join(output_dir, 'nuclei_orientation_overlay.tiff')
plt.savefig(output_fig, dpi=300, bbox_inches='tight')
plt.show()

print(f"Saved overlay visualization to {output_fig}")
