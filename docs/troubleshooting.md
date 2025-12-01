
# **Troubleshooting - Single-Fiber Pipeline**

This page lists common issues, causes, and solutions.

---

## **1. No fibers detected / empty fiber map**

**Cause:** Filename pattern not matched.
**Fix:** Ensure filenames follow the required format:

```
..._23_Wu_<Subject>_<Timepoint>_<Side>.lif_-_FiberX-Y_...
```

---

## **2. nuclei_results.csv is empty**

**Possible causes & fixes:**

* **Filters too strict**

  * Lower `--min_area_px`, `--z_std_threshold`
* **Mask file is incorrect (all black)**

  * Check STDIP output
  * Open in FIJI to verify segmentation

---

## **3. Z-STD too high for all nuclei**

**Cause:** Z-stack is incorrectly exported or saturated.
**Fix:**

* Verify Z intensity profile
* Lower `--z_scale_um_per_index` only if slice spacing was incorrect
* Increase `--z_std_threshold` slightly (2.5–3.0)

---

## **4. Diameter profile is empty**

**Causes:**

* No skeleton detected → check Skel/ folder
* Nuclei too sparse within radius → increase `--fiber_width_max_radius_um`

---

## **5. Orientation values look random**

**Fix:**

* Increase `--skeleton_radius_px` to give PCA more context
* Confirm that the skeleton correctly follows the fiber outline

---

## **6. DBSCAN creates only noise clusters**

**Fix:**

* Raise `--dbscan_eps_um` (Perform elbow method to determine optimal epsilon; eg. 20 → 25–30 µm)
* Lower `--dbscan_min_samples` (2 → 1–2)

---

## **7. NN3/NN5 values are NaN**

**Cause:** Too few included nuclei.
**Fix:**

* Inspect excluded_nuclei.csv
* Confirm area/Z filters aren’t removing everything

---

## **8. Excel converts fiber IDs into dates**

**Solution:**
The script already enforces text formatting in `<Subject>_results.xlsx`, but if manually editing:

* Prepend an apostrophe: `'1-1'`
* Or set the column to “Text” in Excel

---

## **9. Out-of-memory errors**

**Fix:**

* Run per-subject instead of full batch
* Use a smaller dataset to validate workflow
* Ensure >16 GB RAM for large cohorts (64+ GB RAM recommended)

---

## **10. Overlay.png appears blank**

**Cause:** Mask is all zeros or misaligned with lab image.
**Fix:**

* Check STDIP mask generation in FIJI
* Open overlay with contrast adjusted

