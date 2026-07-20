def append_sample_metadata(df,
                           sample_id=None,
                           image_id=None,
                           group=None,
                           pixel_size_um=None,
                           magnification=None):
    """
    Add sample/image metadata columns to a dataframe.

    Aim:
    - keep every exported table self-describing across multi-image runs
    - avoid losing sample identity when tables are pooled downstream
    """
    df = df.copy()
    df["sample_id"] = sample_id
    df["image_id"] = image_id
    df["group"] = group
    df["pixel_size_um"] = pixel_size_um
    df["magnification"] = magnification
    return df

def safe_fraction(numerator, denominator, fallback=np.nan):
    """
    Safely divide two numbers and return a fallback if the denominator is zero.

    Aim:
    - avoid division-by-zero failures in image-level QC summaries
    - keep batch outputs numerically stable on empty or sparse images
    """
    try:
        denominator = float(denominator)
        if denominator == 0:
            return fallback
        return float(numerator) / denominator
    except Exception:
        return fallback

def robust_quantile(x, q, fallback=None):
    x = np.asarray(x)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return fallback
    return float(np.quantile(x, q))

def summarize_components(mask, connectivity=2):
    """
    Return connected-component summary statistics.
    """
    lbl = measure.label(mask, connectivity=connectivity)
    props = measure.regionprops_table(
        lbl,
        properties=("label", "area", "bbox", "eccentricity", "solidity", "extent")
    )
    df = pd.DataFrame(props)
    return lbl, df

def estimate_cleanup_parameters(component_df,
                                min_component_quantile=0.10,
                                area_threshold_quantile=0.25,
                                large_object_quantile=0.99,
                                max_reasonable_area_multiplier=4.0,
                                min_plausible_object_area_px=16,
                                max_plausible_object_area_px=50000,
                                min_min_size_px=4,
                                max_min_size_px=5000,
                                min_hole_area_px=4,
                                max_hole_area_px=10000):
    """
    Estimate morphology cleanup parameters from connected-component areas.

    Aim of this update:
    - keep cleanup adaptive to the current image
    - reduce instability from extreme debris or giant artifacts
    - return a few extra diagnostics for downstream QC

    Strategy:
    - use plausible component sizes to estimate cleanup thresholds
    - use the full distribution only for upper-tail large-object flagging
    - cap thresholds to avoid pathological values on unusual images
    """
    # Return safe defaults if no components exist
    if component_df.shape[0] == 0:
        return {
            "min_size": 0,
            "area_threshold": 0,
            "large_object_cutoff": np.inf,
            "median_area": np.nan,
            "n_components": 0,
            "n_plausible_components": 0,
            "fraction_plausible_components": np.nan
        }

    # Extract positive finite component areas
    areas = np.asarray(component_df["area"].values, dtype=float)
    areas = areas[np.isfinite(areas)]
    areas = areas[areas > 0]

    # Return safe defaults if usable areas disappear after filtering
    if areas.size == 0:
        return {
            "min_size": 0,
            "area_threshold": 0,
            "large_object_cutoff": np.inf,
            "median_area": np.nan,
            "n_components": 0,
            "n_plausible_components": 0,
            "fraction_plausible_components": np.nan
        }

    # ------------------------------------------------------------
    # Background:
    # Use a plausibility-filtered subset to estimate cleanup thresholds
    # so tiny debris and huge artifacts do not dominate the lower-tail summaries.
    # ------------------------------------------------------------
    # Keep only biologically plausible objects for threshold estimation
    plausible_areas = areas[
        (areas >= min_plausible_object_area_px) &
        (areas <= max_plausible_object_area_px)
    ]

    # Fallback if the plausible filter is too strict for a given image
    if plausible_areas.size == 0:
        plausible_areas = areas.copy()

    # Estimate lower-tail cleanup thresholds from plausible objects
    q_min = robust_quantile(plausible_areas, min_component_quantile, fallback=1)
    q_hole = robust_quantile(plausible_areas, area_threshold_quantile, fallback=q_min)
    # Estimate upper-tail large-object threshold from all objects
    q_large = robust_quantile(areas, large_object_quantile, fallback=np.max(areas))
    median_area = robust_quantile(plausible_areas, 0.50, fallback=q_hole)

    # Bound minimum object removal threshold to a safe range
    min_size = max(min_min_size_px, int(round(q_min)))
    min_size = min(min_size, max_min_size_px)
    # Bound small-hole filling threshold to a safe range
    area_threshold = max(min_hole_area_px, int(round(q_hole)))
    area_threshold = min(area_threshold, max_hole_area_px)

    # Constrain the extreme upper threshold so one bizarre giant object does not dominate.
    capped_large = q_large
    if median_area is not None and np.isfinite(median_area):
        capped_large = min(q_large, median_area * max_reasonable_area_multiplier)

    # Cap extreme large-object cutoff so one giant artifact does not dominate
    large_object_cutoff = max(area_threshold, int(round(capped_large)))

    # Return cleanup settings plus a few QC diagnostics
    return {
        "min_size": int(min_size),
        "area_threshold": int(area_threshold),
        "large_object_cutoff": int(large_object_cutoff),
        "median_area": float(median_area),
        "n_components": int(areas.size),
        "n_plausible_components": int(plausible_areas.size),
        "fraction_plausible_components": float(plausible_areas.size / areas.size)
    }

def clean_input_mask(mask,
                     min_size,
                     area_threshold,
                     connectivity=2,
                     fill_holes=True,
                     border_clear=False):
    """
    Conservative cleanup for binary nuclear masks.

    Aim:
    - remove tiny disconnected bright debris
    - fill small internal holes inside likely nuclei
    - optionally remove border-touching objects

    Notes:
    - this function is intentionally conservative
    - it assumes the DAPI input is already binarized
    """
    out = np.asarray(mask).astype(bool).copy()

    # Remove very small bright debris
    out = morphology.remove_small_objects(
        out,
        min_size=min_size,
        connectivity=connectivity
    )

    # Fill small internal holes
    if fill_holes:
        out = morphology.remove_small_holes(
            out,
            area_threshold=area_threshold,
            connectivity=connectivity
        )

    # Optionally remove border-touching objects
    if border_clear:
        out = segmentation.clear_border(out)

    return out

def estimate_nuclear_scale_parameters(component_df,
                                      pixel_size_um=None,
                                      sigma_as_fraction_of_radius=0.20,
                                      peak_footprint_as_fraction_of_diameter=0.50,
                                      min_dt_smoothing_sigma=0.8,
                                      max_dt_smoothing_sigma=3.0,
                                      min_peak_footprint=5,
                                      max_peak_footprint=21,
                                      min_plausible_nucleus_area_px=16,
                                      max_plausible_nucleus_area_px=50000):
    """
    Estimate scale-aware watershed parameters for splitting merged nuclei.

    Aim:
    - reduce dependence on demo-specific raw-pixel settings
    - derive watershed parameters from the observed nucleus-size distribution
    - preserve a simple fallback if image scale metadata are unavailable
    """
    # Return conservative fallback values if no components exist
    if component_df.shape[0] == 0:
        return {
            "estimated_typical_nucleus_area_px": np.nan,
            "estimated_nucleus_radius_px": np.nan,
            "estimated_nucleus_diameter_px": np.nan,
            "dt_smoothing_sigma": float(min_dt_smoothing_sigma),
            "peak_footprint": int(min_peak_footprint),
            "pixel_size_um": pixel_size_um
        }

    # Extract positive finite component areas
    areas = np.asarray(component_df["area"].values, dtype=float)
    areas = areas[np.isfinite(areas)]
    areas = areas[areas > 0]

    # Focus on plausible nuclei when estimating typical scale
    plausible_areas = areas[
        (areas >= min_plausible_nucleus_area_px) &
        (areas <= max_plausible_nucleus_area_px)
    ]

    # Fall back to all areas if plausibility bounds are too strict
    if plausible_areas.size == 0:
        plausible_areas = areas.copy()

    # Use median plausible area as a stable estimate of typical nucleus size
    typical_area = robust_quantile(plausible_areas, 0.50, fallback=np.median(plausible_areas))
    # Convert area to approximate radius and diameter
    radius_px = np.sqrt(typical_area / np.pi)
    diameter_px = 2.0 * radius_px

    # Set smoothing relative to estimated nucleus radius
    sigma = radius_px * sigma_as_fraction_of_radius
    sigma = float(np.clip(sigma, min_dt_smoothing_sigma, max_dt_smoothing_sigma))

    # Set local-max footprint relative to estimated nucleus diameter
    peak_footprint = int(round(diameter_px * peak_footprint_as_fraction_of_diameter))
    peak_footprint = int(np.clip(peak_footprint, min_peak_footprint, max_peak_footprint))

    # Use an odd footprint for a symmetric local-max neighborhood
    if peak_footprint % 2 == 0:
        peak_footprint += 1
        if peak_footprint > max_peak_footprint:
            peak_footprint -= 2

    # Return estimated nucleus scale plus resolved watershed settings
    return {
        "estimated_typical_nucleus_area_px": float(typical_area),
        "estimated_nucleus_radius_px": float(radius_px),
        "estimated_nucleus_diameter_px": float(diameter_px),
        "dt_smoothing_sigma": float(sigma),
        "peak_footprint": int(peak_footprint),
        "pixel_size_um": pixel_size_um
    }

def split_merged_nuclei(mask, sigma=1.0, peak_footprint=9, connectivity=2):
    """
    Attempt to split merged nuclei using distance-transform watershed.
    """
    if np.sum(mask) == 0:
        return np.zeros_like(mask, dtype=np.int32)

    distance = ndi.distance_transform_edt(mask)

    if sigma is not None and sigma > 0:
        distance_smooth = ndi.gaussian_filter(distance, sigma=sigma)
    else:
        distance_smooth = distance

    coords = feature.peak_local_max(
        distance_smooth,
        labels=mask,
        footprint=np.ones((peak_footprint, peak_footprint), dtype=bool),
        exclude_border=False
    )

    markers = np.zeros(mask.shape, dtype=np.int32)
    if coords.shape[0] == 0:
        # fallback: plain connected components
        return measure.label(mask, connectivity=connectivity)

    for i, (r, c) in enumerate(coords, start=1):
        markers[r, c] = i

    markers = ndi.label(markers > 0)[0]
    labels = watershed(-distance_smooth, markers, mask=mask)
    return labels.astype(np.int32)


def regionprops_dataframe(label_img):
    if np.max(label_img) == 0:
        return pd.DataFrame(columns=[
            "label", "area", "centroid-0", "centroid-1",
            "bbox-0", "bbox-1", "bbox-2", "bbox-3",
            "eccentricity", "solidity", "extent",
            "major_axis_length", "minor_axis_length",
            "equivalent_diameter_area", "perimeter"
        ])

    props = measure.regionprops_table(
        label_img,
        properties=(
            "label", "area", "centroid", "bbox",
            "eccentricity", "solidity", "extent",
            "major_axis_length", "minor_axis_length",
            "equivalent_diameter_area", "perimeter"
        )
    )
    return pd.DataFrame(props)

def add_border_flag(df, label_img):
    """
    Add a touches_border flag to a regionprops dataframe.

    Aim:
    - keep partial edge nuclei by default
    - allow later QC or filtering without deleting them up front
    """
    df = df.copy()

    # Return an empty boolean column if no rows exist
    if df.shape[0] == 0:
        df["touches_border"] = pd.Series(dtype=bool)
        return df

    # Store image bounds for bbox checks
    nrows, ncols = label_img.shape
    touches = []

    # Mark objects whose bounding box touches any image edge
    for _, row in df.iterrows():
        r0 = int(row["bbox-0"])
        c0 = int(row["bbox-1"])
        r1 = int(row["bbox-2"])
        c1 = int(row["bbox-3"])

        is_border = (r0 <= 0) or (c0 <= 0) or (r1 >= nrows) or (c1 >= ncols)
        touches.append(bool(is_border))

    # Append the border-touching QC flag
    df["touches_border"] = touches
    return df

def build_tissue_mask(
    nuclear_mask,
    cytoplasmic_mask,
    target_mask=None,
    mode="conservative",
    source="nuc_plus_cyt",
    min_size=None,
    hole_area=None,
    closing_radius=None,
    dilation_radius=None,
    component_policy=None,
    keep_n_largest=None,
    min_component_area=None,
    nuclear_rescue_radius=None
):
    """Build a broad tissue-support mask for removing out-of-sample artifacts."""

    # Convert inputs to boolean arrays so logical and morphology operations behave consistently.
    nuclear_mask = np.asarray(nuclear_mask).astype(bool)
    cytoplasmic_mask = np.asarray(cytoplasmic_mask).astype(bool)

    # Convert the optional target mask only if it is provided.
    if target_mask is not None:
        target_mask = np.asarray(target_mask).astype(bool)

    # Set mode-specific defaults, while allowing explicit arguments to override them.
    if mode == "regular":
        if min_size is None:
            min_size = 5000
        if hole_area is None:
            hole_area = 20000
        if closing_radius is None:
            closing_radius = 15
        if dilation_radius is None:
            dilation_radius = 20
        if component_policy is None:
            component_policy = "largest_n"
        if keep_n_largest is None:
            keep_n_largest = 1
        if min_component_area is None:
            min_component_area = 50000
        if nuclear_rescue_radius is None:
            nuclear_rescue_radius = 0

    elif mode == "conservative":
        if min_size is None:
            min_size = 1000
        if hole_area is None:
            hole_area = 20000
        if closing_radius is None:
            closing_radius = 9
        if dilation_radius is None:
            dilation_radius = 12
        if component_policy is None:
            component_policy = "min_area"
        if keep_n_largest is None:
            keep_n_largest = 1
        if min_component_area is None:
            min_component_area = 20000
        if nuclear_rescue_radius is None:
            nuclear_rescue_radius = 5

    else:
        raise ValueError(f"Unsupported mode: {mode}")

    # Build the initial tissue-support seed from the requested channel set.
    # The nucleus + cytoplasm option is usually the safest default.
    if source == "nuc_plus_cyt":
        tissue_mask = np.logical_or(nuclear_mask, cytoplasmic_mask)
    elif source == "all_channels":
        if target_mask is None:
            raise ValueError("target_mask must be provided when source='all_channels'")
        tissue_mask = nuclear_mask | cytoplasmic_mask | target_mask
    else:
        raise ValueError(f"Unsupported tissue source: {source}")

    # Bridge small gaps and fractures in the tissue-support mask.
    if closing_radius > 0:
        tissue_mask = morphology.binary_closing(
            tissue_mask,
            footprint=morphology.disk(closing_radius)
        )

    # Fill small holes inside the tissue footprint.
    if hole_area > 0:
        tissue_mask = morphology.remove_small_holes(
            tissue_mask,
            area_threshold=hole_area
        )

    # Remove very small disconnected objects that are unlikely to be real tissue.
    if min_size > 0:
        tissue_mask = morphology.remove_small_objects(
            tissue_mask,
            min_size=min_size
        )

    # Slightly expand the tissue mask to preserve edge-adjacent valid signal.
    if dilation_radius > 0:
        tissue_mask = morphology.binary_dilation(
            tissue_mask,
            footprint=morphology.disk(dilation_radius)
        )

    # Label connected tissue regions for component-based filtering.
    tissue_labels = measure.label(tissue_mask)
    props = measure.regionprops(tissue_labels)

    # Keep only the largest N components when using the stricter policy.
    if component_policy == "largest_n":
        props = sorted(props, key=lambda x: x.area, reverse=True)
        keep_labels = [
            p.label for p in props[:keep_n_largest]
            if p.area >= min_component_area
        ]

    # Keep all components above the minimum area when using the more permissive policy.
    elif component_policy == "min_area":
        keep_labels = [
            p.label for p in props
            if p.area >= min_component_area
        ]
    else:
        raise ValueError(f"Unsupported component_policy: {component_policy}")

    # Rebuild the mask from retained components, or return an empty mask if none pass.
    if keep_labels:
        tissue_mask = np.isin(tissue_labels, keep_labels)
    else:
        tissue_mask = np.zeros_like(tissue_mask, dtype=bool)

    # Optionally rescue nucleus-supported regions that may have been trimmed away.
    if nuclear_rescue_radius > 0:
        nuclear_rescue = morphology.binary_dilation(
            nuclear_mask,
            footprint=morphology.disk(nuclear_rescue_radius)
        )
        tissue_mask = np.logical_or(tissue_mask, nuclear_rescue)

    # Return the final cleaned tissue-support mask.
    return tissue_mask

# Save a single binary mask as a grayscale image
def save_mask_plot(ar, out_path, cmap="gray", title=None):
    plt.imshow(ar, cmap=cmap)
    plt.xticks([])
    plt.yticks([])
    plt.savefig(out_path, dpi=600, bbox_inches="tight")
    plt.close()

# Save a three-layer overlay image
def overlay_plot(arr_1, arr_2, out_path, arr_3=None):
    plt.imshow(arr_1, cmap=ListedColormap([(1, 1, 1, 0), "dodgerblue"]), alpha=0.5)
    plt.imshow(arr_2, cmap=ListedColormap([(1, 1, 1, 0), "orange"]), alpha=0.5)
    if arr_3 is not None:
        plt.imshow(arr_3, cmap=ListedColormap([(1, 1, 1, 0), "violet"]), alpha=0.5)
    plt.xticks([])
    plt.yticks([])
    plt.savefig(out_path, dpi=600, bbox_inches="tight")
    plt.close()

def save_labeled_overlay(binary_mask, label_img, out_path, max_objects=300):
    """
    Save a QC image showing the cleaned binary mask with a semi-random label overlay.

    Aim:
    - make object separation visually inspectable after cleanup
    - keep dense images readable by optionally truncating the displayed label range
    """
    binary_mask = np.asarray(binary_mask).astype(bool)
    label_img = np.asarray(label_img)

    plt.figure(figsize=(8, 8))
    plt.imshow(binary_mask, cmap="gray", alpha=0.6)

    if np.max(label_img) > 0:
        shown = label_img.copy()
        if max_objects is not None and max_objects > 0:
            shown[shown > max_objects] = 0
        plt.imshow(shown, cmap="nipy_spectral", alpha=0.5)

    plt.xticks([])
    plt.yticks([])
    plt.savefig(out_path, dpi=600, bbox_inches="tight")
    plt.close()

def build_mask_run_summary(input_mask,
                           cleaned_mask,
                           initial_df,
                           final_df,
                           sample_id=None,
                           image_id=None,
                           group=None,
                           pixel_size_um=None,
                           magnification=None):
    """
    Build an image-level summary for a binary-mask cleanup workflow.

    Metrics:
    - image_area_px: total pixel count in the image
    - input_foreground_pixels: number of positive pixels before cleanup
    - cleaned_foreground_pixels: number of positive pixels after cleanup
    - input_foreground_fraction: fraction of image positive before cleanup
    - cleaned_foreground_fraction: fraction of image positive after cleanup
    - n_initial_components: number of connected components before cleanup
    - n_final_components: number of connected components after cleanup
    - fraction_foreground_retained: cleaned/input foreground ratio
    - fraction_components_retained: final/initial component ratio

    Aim:
    - create one compact per-image QC row
    - support later cross-sample batch review and troubleshooting
    """
    input_mask = np.asarray(input_mask).astype(bool)
    cleaned_mask = np.asarray(cleaned_mask).astype(bool)

    image_area_px = int(input_mask.size)
    input_foreground_pixels = int(input_mask.sum())
    cleaned_foreground_pixels = int(cleaned_mask.sum())

    n_initial_components = int(initial_df.shape[0]) if initial_df is not None else 0
    n_final_components = int(final_df.shape[0]) if final_df is not None else 0

    return {
        "sample_id": sample_id,
        "image_id": image_id,
        "group": group,
        "pixel_size_um": pixel_size_um,
        "magnification": magnification,
        "image_height_px": int(input_mask.shape[0]),
        "image_width_px": int(input_mask.shape[1]),
        "image_area_px": image_area_px,
        "input_foreground_pixels": input_foreground_pixels,
        "cleaned_foreground_pixels": cleaned_foreground_pixels,
        "input_foreground_fraction": safe_fraction(input_foreground_pixels, image_area_px),
        "cleaned_foreground_fraction": safe_fraction(cleaned_foreground_pixels, image_area_px),
        "n_initial_components": n_initial_components,
        "n_final_components": n_final_components,
        "fraction_foreground_retained": safe_fraction(cleaned_foreground_pixels, input_foreground_pixels),
        "fraction_components_retained": safe_fraction(n_final_components, n_initial_components)
    }

def flag_sparse_or_artifact_image(input_foreground_fraction,
                                  cleaned_foreground_fraction,
                                  n_initial_components,
                                  n_final_components,
                                  low_signal_threshold=0.001,
                                  collapse_ratio_threshold=0.10):
    """
    Return lightweight QC flags for sparse or unstable images.

    Aim:
    - provide batch-friendly warning columns without hard-failing the run
    - distinguish low-signal images from images where cleanup drastically changed structure
    """
    flags = {}

    flags["flag_low_signal_image"] = bool(
        np.isfinite(input_foreground_fraction) and input_foreground_fraction < low_signal_threshold
    )

    component_ratio = safe_fraction(n_final_components, n_initial_components)
    flags["flag_component_collapse"] = bool(
        np.isfinite(component_ratio) and component_ratio < collapse_ratio_threshold
    )

    flags["flag_cleanup_expanded_foreground"] = bool(
        np.isfinite(input_foreground_fraction) and
        np.isfinite(cleaned_foreground_fraction) and
        cleaned_foreground_fraction > input_foreground_fraction * 1.5
    )

    return flags
