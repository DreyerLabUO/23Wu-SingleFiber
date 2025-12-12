# **Parameter Guide — Single-Fiber Myonuclear Analysis Pipeline**

*A practical reference for selecting and tuning parameters in `SF_analysis_pipeline.py`.*

This guide describes all user-adjustable parameters in the pipeline, explains what they control, how they affect downstream calculations, and provides recommended ranges for most imaging scenarios.

---

# **1. Imaging Geometry Parameters**

These parameters must match your microscope acquisition settings.

---

## **`--pixel_size_xy_um`**

**Default:** `1.0`
**Meaning:** Microns per pixel in the XY plane.
**Used for:**

* Converting nuclei centroids to µm
* Computing distances, DBSCAN radius, NN3/NN5 evaluation
* Diameter estimation
* Volume calculation

### **When to adjust:**

Always. Set this to the correct value for your imaging objective + camera pixel pitch.

---

## **`--z_scale_um_per_index`**

**Default:** `2.7`
**Meaning:** Distance (microns) between Z-slices.
**Used for:**

* Computing 3D nuclear positions
* Z-consistency filtering
* DBSCAN 3D clustering

### **When to adjust:**

Always—must match the Z-step used in acquisition.

---

# **2. Segmentation Filters**

These parameters control which nuclei are included in the analysis.

---

## **`--min_area_px`**

**Default:** `0`
**Meaning:** Minimum nuclear area (in pixels).
**Purpose:** Filter out debris, noise, small fragments.

### **When to adjust:**

* Segmentation produces small off-target objects
* Low background images → smaller threshold
* Noisy images → larger threshold

---

## **`--max_area_px`**

**Default:** `0` (disabled)
**Meaning:** Maximum allowed nuclear area.
**Purpose:** Remove merged nuclei or segmentation artifacts.

### **When to adjust:**

* Multiple nuclei merge into a single large blob
* Hyperintense background causes mask expansion

---

## **`--z_std_threshold`**

**Default:** `2.0`
**Meaning:** Maximum allowed standard deviation of per-pixel Z indices.
**Purpose:** Filter out *vertically overlapping* nuclei (e.g., stacked objects or ambiguous 3D merges).

### **Effects of tuning:**

* **Lower threshold (1.5–1.8):** More strict → fewer nuclei, higher confidence
* **Higher threshold (2.5–3.0):** More tolerant → better for noisy Z-stacks

### **Troubleshooting hint:**

If *all* nuclei get excluded → Z-stack likely mislabeled or saturated.

---

# **3. Orientation & Geometry Parameters**

---

## **`--skeleton_radius_px`**

**Default:** `20 px`
**Meaning:** Pixel radius around each nucleus used to extract local skeleton points for PCA estimation of fiber direction.

### **When to adjust:**

* Skeleton is broken → **increase radius**
* Fiber is very curved → **decrease radius**
* Orientation appears noisy → **increase mildly**

### **Typical:** 10–40 px.

---

## **`--fiber_width_step_um`**

**Default:** `100 µm`
**Meaning:** Sampling interval for estimating fiber diameter along the skeleton.

### **When to adjust:**

* Need more detailed width profile → lower to 50 µm
* Long fibers, large datasets → keep at 100–150 µm

---

## **`--fiber_width_max_radius_um`**

**Default:** `100 µm`
**Meaning:** Maximum distance from the skeleton to search for nuclei representing fiber edges.

### **When to adjust:**

* **Thick fibers** → increase radius
* **Sparse nuclei** → increase radius
* **Crowded regions** or central nuclei → may need lower radius

### **Typical:** 50–150 µm.

---

# **4. Clustering Parameters (3D Nuclear Spatial Organization)**

---

## **`--dbscan_eps_um`**

**Default:** `20 µm`
**Meaning:** Maximum distance between nuclei for them to be considered “neighbors” in DBSCAN.

### **When to adjust:**

* If every nucleus forms its own cluster → increase
* If clusters merge too easily → decrease

### Perform the [**elbow method for DBSCAN**](https://medium.com/@tarammullin/dbscan-parameter-estimation-ff8330e3a3bd) to validate optimal epsilon value for each project.

---

## **`--dbscan_min_samples`**

**Default:** `2`
**Meaning:** Minimum number of nuclei required to form a cluster.

### **Guidelines:**

* Set to **2** for permissive clustering (detect small local groups)
* Increase to **3–4** for stricter, larger clusters

---

# **5. Nearest-Neighbor Parameters (NN3 / NN5 Analysis)**

Nearest-neighbor metrics quantify **local nuclear environments**:

* Mean NN distance
* Proportion of neighboring shapes (spherical / ellipsoid / intermediate)
* Orientation context (parallel / perpendicular / intermediate relative to fiber)
* Mean relative angle of neighbors

### **Important:**

**There are no user-adjustable parameters.**
These computations are entirely data-driven and scale naturally with nuclear density.

If NN metrics output NaN → too few nuclei passed filtering.

---

# **6. Additional Internal Behavior**

### **Filename Parsing**

To be compatible with the published pipeline, all images must follow Dreyer Lab’s naming convention:

```
23_Wu_<Subject>_<Timepoint>_<Side>_-_FiberX-Y_*.tif
```

If mismatched:

* Files will be skipped
* The troubleshooting page provides filename debugging steps

To adapt to differing naming convention:

* Modify Regex patterns accordingly throughout script


---

# **7. How to Diagnose Parameter Issues**

| Symptom                    | Likely Cause                | Parameter to Adjust                    |
| -------------------------- | --------------------------- | -------------------------------------- |
| Too many nuclei excluded   | Z filter too strict         | Increase `--z_std_threshold`           |
| Orientation looks unstable | Not enough skeleton context | Increase `--skeleton_radius_px`        |
| Diameter profile empty     | Too few nuclei near edges   | Increase `--fiber_width_max_radius_um` |
| Clusters too small         | DBSCAN too strict           | Increase `--dbscan_eps_um`             |
| All nuclei in one cluster  | DBSCAN too permissive       | Decrease `--dbscan_eps_um`             |

---

# **8. Parameter Interaction Notes**

* **Z-STD filtering** strongly influences NN metrics—fewer nuclei → weaker microenvironment representation
* Proper **pixel size** is crucial; wrong scaling propagates error into all distance-based analyses. Measure representative objects in FIJI to gauge scaling
* **Skeleton radius** and **step size** interact—smaller radius + sharp curvature may distort PCA-based orientation
* **DBSCAN** is sensitive to both XY and Z scaling; ensure correct calibration

---

# **9. Summary**

This parameter guide assists users in tuning the pipeline for:

* Variable imaging modalities
* Different magnifications
* Complex nuclear morphologies
* Noisy z-stacks
* High-density vs sparse fibers

Proper parameter selection ensures accurate, biologically meaningful quantification of myonuclear organization and microenvironmental structure.
