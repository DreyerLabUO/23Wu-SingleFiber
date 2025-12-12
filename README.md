# **Single-Fiber Myonuclear Analysis Pipeline**

**A reproducible Python framework for quantifying nuclear morphology, spatial organization, microenvironment structure, and fiber geometry in single skeletal muscle fibers.**

This repository contains the full analysis workflow used in the Dreyer Lab’s single-fiber nuclear organization studies.
The pipeline integrates output from a FIJI/ImageJ segmentation macro with automated Python-based quantification, and is designed for high-throughput, auditable analysis across subjects and timepoints.

---

## **Features**

* Nuclear morphometrics (area, axis lengths, shape classification)
* Z-position inference & **Z-consistency filtering**
* Orientation relative to the fiber’s local axis (PCA-based)
* Distance to fiber skeleton
* **3D DBSCAN clustering** of nuclei
* **Nearest-neighbor microenvironment analysis (NN3 & NN5)**
* Fiber diameter estimation along the length of the fiber
* Fiber volume + myonuclear domain estimation
* Automatic per-fiber, per-leg/biopsy, and per-subject summary outputs
* Fully batch-automated processing across folders

---

## **Repository Structure**

```
23Wu-SingleFiber/
│
├── SF_analysis_pipeline.py         # Main analysis script
├── README.md                       # You are here
├── requirements.txt                # Python dependencies
├── docs/
│     ├── troubleshooting.md
│     └── parameter_guide.md
├── examples/
│     ├── example_macro_output/
│     └── example_results/
├── visualizations/                 # Manuscript figures 8-10
└── LICENSE
```

---

# **Quickstart**

### **1. Create and activate a Conda environment**
Ensure you have Python 3.10+ installed. It is recommended to use a virtual environment to manage dependencies.
```bash
conda create --name sfpipeline python=3.10
conda activate sfpipeline
```

### **2. Install dependencies**

**Option A: Core Analysis (Lightweight)** If you only need to run the main script.
```bash
pip install -r requirements.txt
```

Or manually:

```
python -m pip install numpy pandas scikit-image scikit-learn scipy openpyxl Pillow
```

**Option B: Full Development (Analysis + Representative Figures)** If you are interested in also running Jupyter Notebooks or generating figures (as seen in published manuscript), install the full suite.
```bash
pip install -r requirements-vis.txt
```

Or manually:

```
python -m pip install numpy pandas scikit-image scikit-learn scipy openpyxl Pillow matplotlib seaborn plotly jupyterlab notebook
```


### **3. Prepare your data directory**

Your macro output should follow:

```
<base_dir>/
    20/
        BL/
            L/
                STDIP/
                Skel/
                TIFs/
        D14/
            R/
                STDIP/
                Skel/
                TIFs/
```

With appropriate naming convention, each fiber must have:

* 1 STDIP mask (eg. ```STDIP_23_Wu_20_D14_R_02.lif_-_Fiber1-1_Merged.tif```)
* 1 Skel image (eg. ```Skel_23_Wu_20_D14_R_02.lif_-_Fiber1-1_Merged.tif```)
* 1 TIF z-stack (eg. ```23_Wu_20_D14_R_02.lif_-_Fiber1-1_Merged.tif```)

---

### **4. Run the pipeline**

```bash
python SF_analysis_pipeline.py "D:\MacroOutput" --imaris_master "C:\path\to\Imaris.xlsx" \
    --pixel_size_xy_um 0.329 \
    --z_scale_um_per_index 2.7 \
    --min_area_px 200 \
    --max_area_px 0 \
    --z_std_threshold 2.0 \
    --dbscan_eps_um 20 \
    --skeleton_radius_px 20 \
    --fiber_width_step_um 100 \
    --fiber_width_max_radius_um 100
```

---

# **Major Processing Steps**

1. Load masks, skeletons, and z-stacks
2. Segment & filter nuclei (area + Z-consistency)
3. Extract nuclear geometry & orientation
4. Compute relative orientation to local skeleton (PCA)
5. Compute 3D nuclear positions
6. DBSCAN clustering
7. Compute NN3 & NN5 microenvironment metrics
8. Estimate fiber diameter along path
9. Write per-fiber CSV and visualization
10. Aggregate biopsy and subject-level metrics


---

# **Key Parameters**

For typical imaging setups:

| Parameter                     | Controls               | When to adjust                |
| ----------------------------- | ---------------------- | ----------------------------- |
| [`--pixel_size_xy_um`](./docs/parameter_guide.md#--pixel_size_xy_um)          | XY scaling             | Always match your microscope  |
| [`--z_scale_um_per_index`](./docs/parameter_guide.md#--z_scale_um_per_index)      | Z-step size            | Match acquisition settings    |
| [`--min_area_px`](./docs/parameter_guide.md#--min_area_px)               | Small-object filtering | Remove debris                 |
| [`--max_area_px`](./docs/parameter_guide.md#--max_area_px)               | Large-object filtering | Remove merged nuclei          |
| [`--z_std_threshold`](./docs/parameter_guide.md#--z_std_threshold)           | Z-consistency filter   | Adjust based on stack quality |
| [`--dbscan_eps_um`](./docs/parameter_guide.md#--dbscan_eps_um)             | 3D cluster sensitivity | Depends on nuclear packing    |
| [`--skeleton_radius_px`](./docs/parameter_guide.md#--skeleton_radius_px)        | PCA neighborhood       | Curved/broken skeletons       |
| [`--fiber_width_step_um`](./docs/parameter_guide.md#--fiber_width_step_um)       | Diameter sampling      | Resolution of width profile   |
| [`--fiber_width_max_radius_um`](./docs/parameter_guide.md#--fiber_width_max_radius_um) | Search radius          | Fiber thickness               |

Full descriptions appear in [**docs/parameter_guide.md**](./docs/parameter_guide.md).

---

# **Output Files**

## **Per-Fiber (`*_output/`)**

* `nuclei_results.csv` — morphometrics, Z-metrics, NN3/NN5, orientation, DBSCAN labels
* `excluded_nuclei.csv` — nuclei removed by filters
* `fiber_width_profile.csv` — diameter along the fiber
* `overlay.png` — quality control visualization

## **Per-Leg Biopsy Summary**

* `<Subject>_<Timepoint>_<Side>_biopsy_summary.csv`
  Includes fiber volume, myonuclear domain, neighbor metrics, clustering, shape distribution, etc.

## **Per-Subject Workbook**

* `<Subject>_results.xlsx`
  One sheet per timepoint with all nuclei-level metrics.

---

# **Citation**

If this pipeline contributes to your research, please cite the corresponding Dreyer Lab manuscript.

---

# **Troubleshooting**

See [`docs/troubleshooting.md`](./docs/troubleshooting.md) for the included troubleshooting guide.
