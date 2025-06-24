from osgeo import gdal
from scipy import ndimage
import numpy as np
from scipy.interpolate import interp1d, griddata
from scipy.ndimage import gaussian_filter, median_filter, grey_erosion, grey_dilation, uniform_filter
import os
from multiprocessing import Pool
import networkx as nx


def filter_and_custom_close_image(arr, median_kernel=5, gaussian_sigma=3, erosion_size=3):
    # Apply median filter, Gaussian filter, erosion and dilation to smooth the input array
    nan_mask = np.isnan(arr)
    arr_filled = np.nan_to_num(arr, nan=0)
    median_filtered = median_filter(arr_filled, size=median_kernel)
    gaussian_filtered = gaussian_filter(median_filtered, sigma=gaussian_sigma)
    eroded = grey_erosion(gaussian_filtered, size=(erosion_size, erosion_size))
    dilated = grey_dilation(eroded, size=(erosion_size, erosion_size))
    dilated[nan_mask] = np.nan
    return dilated


def apply_grey_erosion(arr, erosion_size=7):
    # Apply grey erosion while preserving NaN values
    nan_mask = np.isnan(arr)
    arr_filled = np.nan_to_num(arr, nan=0)
    eroded = grey_erosion(arr_filled, size=(erosion_size, erosion_size))
    eroded[nan_mask] = np.nan
    return eroded


def masked_mean_filter_fast(input_array, mask, filter_size):
    # Apply mean filtering only within the valid mask
    input_array = np.where(mask, input_array, 0)
    count_mask = uniform_filter(mask.astype(float), size=filter_size)
    sum_array = uniform_filter(input_array, size=filter_size)
    with np.errstate(invalid='ignore', divide='ignore'):
        result = sum_array / count_mask
    result[~mask] = input_array[~mask]
    return result


def interpolate_2d_array(array):
    # Interpolate missing values in 2D array along rows and columns
    rows, cols = array.shape
    x = np.arange(cols)
    interpolated_rows = []
    for row in array:
        if np.isnan(row).all():
            interpolated_rows.append(np.full(cols, np.nan))
        else:
            interpolated_rows.append(interp1d(x[~np.isnan(row)], row[~np.isnan(row)], kind='linear', bounds_error=False)(x))
    interpolated_rows = np.array(interpolated_rows)

    interpolated_cols = []
    for col in array.T:
        if np.isnan(col).all():
            interpolated_cols.append(np.full(rows, np.nan))
        else:
            interpolated_cols.append(interp1d(np.arange(rows)[~np.isnan(col)], col[~np.isnan(col)], kind='linear', bounds_error=False)(np.arange(rows)))
    interpolated_cols = np.array(interpolated_cols).T

    interpolated = np.where(np.isnan(interpolated_rows) & np.isnan(interpolated_cols),
                            np.nan,
                            np.where(np.isnan(interpolated_rows), interpolated_cols,
                                     np.where(np.isnan(interpolated_cols), interpolated_rows,
                                              (interpolated_rows + interpolated_cols) / 2)))
    if np.isnan(interpolated).any():
        # Use nearest-neighbor interpolation for remaining NaNs
        interpolated = np.ma.masked_invalid(interpolated)
        x2 = np.arange(0, interpolated.shape[1])
        y2 = np.arange(0, interpolated.shape[0])
        xx2, yy2 = np.meshgrid(x2, y2)
        x2 = xx2[~interpolated.mask]
        y2 = yy2[~interpolated.mask]
        newarr2 = interpolated[~interpolated.mask].data
        final_interpolated = griddata((x2, y2), newarr2.ravel(), (xx2, yy2), method='nearest')
        return final_interpolated
    else:
        return interpolated


def filter_block(dsm_data, filterd_data, compare_arr, interpolate_func,
                 pixel_threshold=250000, min_value=1, iteration_value=0.25):
    # Perform filtering and interpolation on local blocks to reduce artifacts
    length, width = dsm_data.shape
    if np.all(np.isnan(dsm_data)) or length * width < pixel_threshold:
        return dsm_data
    DSM = dsm_data.copy()
    max_value = np.nanmax(DSM)
    while np.nanmax(filterd_data) <= max_value:
        img_binary = np.where(DSM > filterd_data, 1, 0)
        zeros_arr = np.zeros_like(DSM)
        labeled_array, num_features = ndimage.label(img_binary)
        masks = [(labeled_array == label_value) for label_value in range(1, num_features + 1)]
        for mask in masks:
            int_mask = mask.astype(int)
            count_ones = np.count_nonzero(int_mask == 1)
            if count_ones < pixel_threshold:
                ones_indices = np.argwhere(int_mask == 1)
                values, compare_values = [], []
                for i, j in ones_indices:
                    # Only consider edge pixels
                    if (i > 0 and int_mask[i - 1, j] == 0) or (i < int_mask.shape[0] - 1 and int_mask[i + 1, j] == 0) or \
                       (j > 0 and int_mask[i, j - 1] == 0) or (j < int_mask.shape[1] - 1 and int_mask[i, j + 1] == 0):
                        values.append(DSM[i, j])
                        compare_values.append(compare_arr[i, j])
                if values:
                    mean_value = np.mean(values)
                    mean_compare_value = np.mean(compare_values)
                    diff = mean_value - mean_compare_value
                    if diff > min_value:
                        zeros_arr += int_mask
        interp_mask = np.where(zeros_arr == 1)
        zz_interp = DSM.astype(float)
        zz_interp[interp_mask] = np.nan
        zz_interp_lin = interpolate_func(zz_interp)
        zz_interp_lin_filtered = masked_mean_filter_fast(zz_interp_lin, zeros_arr.astype(bool), 10)
        DSM[interp_mask] = zz_interp_lin_filtered[interp_mask]
        filterd_data += iteration_value
        max_value = np.nanmax(DSM)
    return DSM


def grey_opening(arr, structure_size=5):
    # Apply grey opening (erosion followed by dilation)
    eroded = grey_erosion(arr, size=(structure_size, structure_size))
    opened = grey_dilation(eroded, size=(structure_size, structure_size))
    return opened


def read_data(file_path, nodata_value):
    # Load raster and convert nodata values to NaN
    dsm_dataset = gdal.Open(file_path)
    dsm_data = dsm_dataset.ReadAsArray().astype(float)
    dsm_data[dsm_data == nodata_value] = np.nan
    mask = np.isnan(dsm_data)
    return dsm_dataset, dsm_data, mask


def save_raster(dem_data, mask, dataset, output_path, nodata_value):
    # Save processed array to a new raster file
    dem_data[mask] = np.nan
    driver = gdal.GetDriverByName('GTiff')
    dem_dataset = driver.Create(output_path, dataset.RasterXSize, dataset.RasterYSize, 1, gdal.GDT_Float32)
    dem_dataset.SetGeoTransform(dataset.GetGeoTransform())
    dem_dataset.SetProjection(dataset.GetProjection())
    band = dem_dataset.GetRasterBand(1)
    band.WriteArray(dem_data)
    band.SetNoDataValue(nodata_value)
    band.FlushCache()
    dem_dataset = None


def min_cut_synthesis(array1, array2):
    # Perform graph cut fusion of two overlapping image blocks
    rows, cols = array1.shape
    G = nx.DiGraph()
    source = "source"
    sink = "sink"
    G.add_node(source)
    G.add_node(sink)

    for i in range(rows):
        for j in range(cols):
            node = (i, j)
            G.add_node(node)
            G.add_edge(source, node, capacity=array1[i, j])
            G.add_edge(node, sink, capacity=array2[i, j])

    cut_value, partition = nx.minimum_cut(G, source, sink)
    reachable, non_reachable = partition

    result = np.zeros_like(array1)
    for node in reachable:
        if node != source:
            i, j = node
            result[i, j] = array1[i, j]
    for node in non_reachable:
        if node != sink:
            i, j = node
            result[i, j] = array2[i, j]
    return result


def block_operation_with_overlap(array, block_size, overlap_size, operation_func, interpolate_func,
                                 median_kernel, gaussian_sigma, erosion_size_filter, erosion_size_compare,
                                 pixel_threshold, min_value, iteration_value, use_min_cut_for_overlap=False):
    # Divide image into blocks with overlap, process, and fuse results
    array_shape = array.shape
    result_array = np.full_like(array, np.nan)

    for i in range(0, array_shape[0], block_size - overlap_size):
        for j in range(0, array_shape[1], block_size - overlap_size):
            block_start_row = i
            block_end_row = min(i + block_size, array_shape[0])
            block_start_col = j
            block_end_col = min(j + block_size, array_shape[1])

            current_block = array[block_start_row:block_end_row, block_start_col:block_end_col]
            filterd_data = filter_and_custom_close_image(current_block, median_kernel, gaussian_sigma, erosion_size_filter)
            compare_arr = apply_grey_erosion(current_block, erosion_size_compare)
            result_block = operation_func(current_block, filterd_data, compare_arr, interpolate_func,
                                          pixel_threshold, min_value, iteration_value)

            existing_block = result_array[block_start_row:block_end_row, block_start_col:block_end_col]
            blended_block = np.where(np.isnan(existing_block), result_block, (existing_block + result_block) / 2)
            result_array[block_start_row:block_end_row, block_start_col:block_end_col] = blended_block

            # Fuse overlapping left area using either minimum or graph cut
            if block_start_col != 0:
                overlap_start_col = block_start_col
                overlap_end_col = overlap_start_col + overlap_size

                left_existing = result_array[block_start_row:block_end_row, overlap_start_col:overlap_end_col]
                left_result = result_block[:, 0:overlap_size]

                left_existing_filled = np.where(np.isnan(left_existing), 0, left_existing)
                left_result_filled = np.where(np.isnan(left_result), 0, left_result)

                if use_min_cut_for_overlap:
                    left_blend = min_cut_synthesis(left_existing_filled, left_result_filled)
                else:
                    left_blend = np.minimum(left_existing_filled, left_result_filled)

                result_array[block_start_row:block_end_row, overlap_start_col:overlap_end_col] = left_blend

    return result_array


def process_tif_file(input_path, output_path, nodata_value,
                     median_kernel, gaussian_sigma, erosion_size_filter, erosion_size_compare,
                     pixel_threshold, min_value, iteration_value, grey_open_size,
                     block_size, overlap_size, use_min_cut_for_overlap=False):
    # Load, process, and save a single raster
    print(f"Processing: {input_path}")
    dataset, dsm_data, nan_mask = read_data(input_path, nodata_value)
    dem_data = block_operation_with_overlap(dsm_data, block_size, overlap_size, filter_block, interpolate_2d_array,
                                            median_kernel, gaussian_sigma, erosion_size_filter, erosion_size_compare,
                                            pixel_threshold, min_value, iteration_value, use_min_cut_for_overlap)
    dem_data = grey_opening(dem_data, grey_open_size)
    save_raster(dem_data, nan_mask, dataset, output_path, nodata_value)
    print(f"Saved: {output_path}")


def batch_process_tif_files(input_folder, output_folder, nodata_value,
                            median_kernel, gaussian_sigma, erosion_size_filter, erosion_size_compare,
                            pixel_threshold, min_value, iteration_value, grey_open_size,
                            block_size, overlap_size, use_min_cut_for_overlap=False):
    # Batch process all .tif files in a folder
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    tasks = [(os.path.join(input_folder, f),
              os.path.join(output_folder, f),
              nodata_value,
              median_kernel, gaussian_sigma, erosion_size_filter, erosion_size_compare,
              pixel_threshold, min_value, iteration_value, grey_open_size,
              block_size, overlap_size, use_min_cut_for_overlap)
             for f in os.listdir(input_folder) if f.lower().endswith('.tif')]
    with Pool(processes=4) as pool:
        pool.starmap(process_tif_file, tasks)


if __name__ == "__main__":
    # Set input and output folders
    input_folder = r"D:\fushi\Vaihingen"
    output_folder = r"D:\fushi\Vaihingen\SAVE\test"

    # Parameters
    nodata_value = -99999
    median_kernel = 5
    gaussian_sigma = 2
    erosion_size_filter = 120
    erosion_size_compare = 7
    pixel_threshold = 250000
    min_value = 1
    iteration_value = 0.25
    grey_open_size = 5
    block_size = 1024
    overlap_size = 64

    # Toggle fusion method for overlapping blocks
    use_min_cut_for_overlap = False  # True = graph cut; False = np.minimum

    batch_process_tif_files(input_folder, output_folder, nodata_value,
                            median_kernel, gaussian_sigma, erosion_size_filter, erosion_size_compare,
                            pixel_threshold, min_value, iteration_value, grey_open_size,
                            block_size, overlap_size, use_min_cut_for_overlap)
