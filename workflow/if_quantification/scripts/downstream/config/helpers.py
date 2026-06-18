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


def robust_quantile(x, q, fallback=None):
    x = np.asarray(x)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return fallback
    return float(np.quantile(x, q))


def estimate_cleanup_parameters(component_df,
                                min_component_quantile=0.10,
                                area_threshold_quantile=0.25,
                                large_object_quantile=0.99,
                                max_reasonable_area_multiplier=4.0):
    """
    Estimate morphology parameters from connected-component areas.

    Strategy:
    - min_size: lower-tail component area estimate, to remove tiny debris
    - area_threshold: lower-mid component area estimate, used for hole removal
    - large_object_cutoff: upper-tail estimate to identify unusually large merged objects
    """
    if component_df.shape[0] == 0:
        return {
            "min_size": 0,
            "area_threshold": 0,
            "large_object_cutoff": np.inf
        }

    areas = np.asarray(component_df["area"].values, dtype=float)
    areas = areas[np.isfinite(areas)]
    areas = areas[areas > 0]

    if areas.size == 0:
        return {
            "min_size": 0,
            "area_threshold": 0,
            "large_object_cutoff": np.inf
        }

    q10 = robust_quantile(areas, min_component_quantile, fallback=1)
    q25 = robust_quantile(areas, area_threshold_quantile, fallback=q10)
    q99 = robust_quantile(areas, large_object_quantile, fallback=np.max(areas))
    median_area = robust_quantile(areas, 0.50, fallback=q25)

    min_size = max(1, int(round(q10)))
    area_threshold = max(1, int(round(q25)))

    # Constrain extreme upper threshold so one bizarre giant object does not dominate.
    capped_large = min(q99, median_area * max_reasonable_area_multiplier) if median_area is not None else q99
    large_object_cutoff = max(area_threshold, int(round(capped_large)))

    return {
        "min_size": min_size,
        "area_threshold": area_threshold,
        "large_object_cutoff": large_object_cutoff,
        "median_area": float(median_area),
        "n_components": int(areas.size)
    }


def clean_nuclear_mask(mask,
                       min_size,
                       area_threshold,
                       connectivity=2,
                       fill_holes=True,
                       border_clear=False):
    """
    Conservative cleanup for binary nuclear masks.
    """
    out = mask.copy()

    # Remove very small bright debris
    out = morphology.remove_small_objects(out, min_size=min_size, connectivity=connectivity)

    # Fill small internal holes
    if fill_holes:
        out = morphology.remove_small_holes(out, area_threshold=area_threshold, connectivity=connectivity)

    # Optionally remove border-touching objects
    if border_clear:
        out = segmentation.clear_border(out)

    return out


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


def save_mask(path, arr):
    # uint8 is convenient for masks
    io.imsave(path, util.img_as_ubyte(arr > 0), check_contrast=False)


def save_label_tiff(path, lbl):
    # save labels as uint16 when possible
    out = np.asarray(lbl)
    if out.max() <= np.iinfo(np.uint16).max:
        out = out.astype(np.uint16)
    else:
        out = out.astype(np.uint32)
    io.imsave(path, out, check_contrast=False)


def plot_component_area_hist(component_df, out_path=None):
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    if component_df.shape[0] > 0:
        ax.hist(component_df["area"], bins=50)
    ax.set_title("Connected-component area distribution", fontsize=14)
    ax.set_xlabel("Component area (pixels)", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    if out_path is not None:
        fig.savefig(out_path, dpi=200)
    plt.close(fig)


def plot_segmentation_qc(binary_mask, cleaned_mask, label_img, out_path=None, max_objects=500):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)

    axes[0].imshow(binary_mask, cmap="gray")
    axes[0].set_title("Input binary DAPI", fontsize=14)
    axes[0].axis("off")

    axes[1].imshow(cleaned_mask, cmap="gray")
    axes[1].set_title("Cleaned nuclear mask", fontsize=14)
    axes[1].axis("off")

    # optionally simplify very dense overlays
    lbl_show = label_img.copy()
    unique_labels = np.unique(lbl_show)
    unique_labels = unique_labels[unique_labels != 0]
    if unique_labels.size > max_objects:
        keep = set(unique_labels[:max_objects])
        lbl_show = np.where(np.isin(lbl_show, list(keep)), lbl_show, 0)

    overlay = color.label2rgb(lbl_show, image=cleaned_mask, bg_label=0, alpha=0.35)
    axes[2].imshow(overlay)
    axes[2].set_title("Segmentation overlay", fontsize=14)
    axes[2].axis("off")

    if out_path is not None:
        fig.savefig(out_path, dpi=200)
    plt.close(fig)
