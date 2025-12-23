# **Single-Fiber Myonuclear Analysis Pipeline**

**A reproducible framework for quantifying nuclear morphology, spatial organization, microenvironment structure, and fiber geometry in single skeletal muscle fibers.**

This repository contains the full analysis workflow used in the Dreyer Lab’s single-fiber nuclear organization studies.

The workflow consists of **two main parts**:

1. **Part I**: Automated 2D image processing using a custom FIJI/ImageJ Macro.

2. **Part II**: Quantitative analysis using a Python Pipeline.

***Note**: Sample data is provided for testing either stage independently. Users wishing to complete a full analysis should run Part I first to generate the necessary inputs for Part II.*

---

## **Features**

* **Automated Processing**: Batch conversion of raw confocal stacks into analyzable 2D projections and skeletons.
* **Nuclear Morphometrics**: Area, axis lengths, shape classification.
* **Spatial Organization**: Z-position inference, Z-consistency filtering, and 3D DBSCAN clustering.
* **Microenvironment**: Nearest-neighbor analysis (NN3 & NN5) and myonuclear domain estimation.
* **Geometry**: Fiber diameter estimation along the length of the fiber.
* **Reporting**: Automatic per-fiber, per-leg/biopsy, and per-subject summary outputs.

---

## **Repository Structure**

```
23Wu-SingleFiber/
│
├── SF_analysis_pipeline.py         # Python analysis script
├── 23_Wu_SF_Analysis.ijm           # FIJI/ImageJ macro script
├── README.md                       # You are here
├── requirements.txt                # Python dependencies (lightweight)
├── requirements.txt                # Python dependencies (with figures)
├── docs/
│     ├── troubleshooting.md
│     └── parameter_guide.md
├── examples/
│     ├── example_inputs/           # Example inputs, hosted on Dryad
│     └── example_results/
├── visualizations/                 # Manuscript figures 8-10
└── LICENSE
```

---
# **Part I: Automated 2D Image Processing (FIJI Macro)**

This step processes raw ```.lif``` (Leica) project files into the standardized directory structure required for the Python analysis. An example Leica project is available for download on Dryad. 

### A. **Setup and Execution**

1. **Download Fiji/ImageJ** and enable the Bio-Formats update site.

2. **Import the Macro**: Launch Fiji and open the ```23_Wu_SF_Analysis.ijm``` file from this repository.

3. **Prepare Directories**: Create an *Input* folder containing your project (```.lif```) files and an empty *Output* folder for the processed images.

4. **Run**: Click "Run" on the macro editor window. Select the input and output directories when prompted.

5. **Verify**: Once finished, check the output folder. The macro automatically organizes files into ```STDIP```, ```Skel```, and ```TIFs``` directories.

### B. **File Naming & Regex Adaptations**

The macro uses Regular Expressions (Regex) to parse metadata from filenames. Our standardized convention is:
```[StudyCode]_[SubjectID]_[Timepoint]_[Leg]_Merged.lif```

* **Standard Regex Pattern**: ```([0-9]+)_([A-Za-z]+)_([0-9]+)_([A-Za-z0-9]+)_([A-Za-z]+)_Merged```

* **Example**: ```23_Wu_01_D14_L_Merged.lif``` (Study 23 Wu, Subject 01, Day 14, Left Leg).

**Adaptation**: If your naming convention differs, you must adjust the Regex pattern in the ```.ijm``` file (lines responsible for file parsing) to match your structure.

### **Automated Processing Steps**

The macro performs the following operations automatically:

1. **Raw TIF Saving**: Extracts and saves the raw z-stack. Used downstream to identify the 3D Z-position of nuclei.

2. **STD Z-Projection (STDIP)**:
* Applies rolling ball background subtraction (radius = 30 px).
* Collapses z-stack into a 2D Standard Deviation Intensity Projection (STDIP).
* Applies median filter and Otsu thresholding to create a binary mask of nuclei.

3. **Skeleton Generation (Skel)**:
* Creates a fiber mask by blurring and thresholding the projection.
* Isolates large objects (>70,000 px²) to remove artifacts.
* Skeletonizes the mask to represent fiber orientation and length.


# **Part II: Analysis with Python**

Once the images are processed (or if using the provided ```example_macro_output``` data, accessible via Dryad), use the Python pipeline to quantify morphology and organization.


### 1. **Create and activate a Conda environment**
Ensure you have Python 3.10+ installed. It is recommended to use a virtual environment to manage dependencies.
```bash
conda create --name sfpipeline python=3.10
conda activate sfpipeline
```

### 2. **Install dependencies**

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


### 3. **Prepare your data directory**

If you ran **Part I**, your output folder is already structured correctly. If you are testing **Part II only**, use the provided ```example_macro_ouput```.

The required structure is:
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

### 4. **Run the pipeline**

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

### **Automated Processing Steps**

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

* `nuclei_results.csv` — Morphometrics, Z-metrics, NN3/NN5, Orientation, DBSCAN labels
* `excluded_nuclei.csv` — Nuclei removed by filters
* `fiber_width_profile.csv` — Diameter along the fiber
* `overlay.png` — Quality control visualization

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
