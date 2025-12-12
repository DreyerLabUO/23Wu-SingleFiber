import os
import cv2
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from skimage import io

# === User Inputs ===
nuclei_csv = r"nuclei_measurements_path"        
skeleton_img = r"skeleton_path"     
output_dir = r"desired_output_dir"               # Directory to save outputs


radius = 20  # neighborhood radius for PCA (pixels)
nucleus_index = 0  # index of nucleus to visualize (choose example nucleus)

# --- Load data ---
df = pd.read_csv(nuclei_csv)
nucleus = df.iloc[nucleus_index]
cx, cy = nucleus["Location_Center_X"], nucleus["Location_Center_Y"]
nuc_ang = np.deg2rad(nucleus["AreaShape_Orientation"])

skel = io.imread(skeleton_img)
skel_coords = np.column_stack(np.where(skel > 0))

# --- Find local neighborhood points ---
dists = np.sqrt((skel_coords[:,0] - cy)**2 + (skel_coords[:,1] - cx)**2)
idx = np.where(dists < radius)[0]
if len(idx) < 5:
    idx = np.argsort(dists)[:10]
local_pts = skel_coords[idx]
pts = np.column_stack((local_pts[:,1], local_pts[:,0])).astype(float)

# --- PCA ---
mean_pts = pts.mean(axis=0)
cov = np.cov((pts - mean_pts).T)
eigvals, eigvecs = np.linalg.eig(cov)
principal = eigvecs[:, np.argmax(eigvals)]
secondary = eigvecs[:, np.argmin(eigvals)]

# --- Visualization ---
fig, ax = plt.subplots(figsize=(6, 6), dpi=300)
ax.imshow(skel, cmap="gray")
ax.scatter(pts[:,0], pts[:,1], c="blue", s=20, label="Skeleton Points")
ax.scatter(cx, cy, marker="*", c="yellow", s=120, label="Nucleus Centroid")

# Principal (red) and secondary (green) PCA axes
scale = 30
ax.arrow(mean_pts[0], mean_pts[1], scale*principal[0], scale*principal[1],
         color="red", width=1, label="Fiber Orientation (PCA)")
ax.arrow(mean_pts[0], mean_pts[1], scale*secondary[0], scale*secondary[1],
         color="green", width=1, label="Secondary PCA Axis")

# Nucleus orientation arrow (purple)
ax.arrow(cx, cy, 20*np.cos(nuc_ang), 20*np.sin(nuc_ang), color="purple", width=1, label="Nucleus Axis")

ax.set_title("PCA-Derived Local Fiber Orientation", fontsize=16, weight="bold")
ax.axis("off")
ax.legend(fontsize=10, loc="upper right")

plt.tight_layout()
plt.savefig(os.path.join(output_dir, "pca_orientation_diagram.png"), dpi=300, bbox_inches="tight")
plt.show()