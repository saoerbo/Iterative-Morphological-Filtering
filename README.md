A Iterative Morphological Filtering Algorithm for Producing DEM from Gridded DSM
Shaobo Linghu, Wenlong Song, Yizhu Lu, Kaizheng Xiang, Hongjie Liu, Long Chen, Tianshi Feng, Rongjie Gui, Yao Zhao and Haider Abbas

This is an algorithm for producing DEM from gridded DSM

Parameter explanation：
*input_folder is the folder where the image is read
*output_folder is the folder where the processed image is output
*nodata_value is the nodata value in the raster imagery and is used for masking
*median_kernel and gaussian_sigma are filtering parameters used to preprocess the image, remove the small noise in the image, and adjust the size according to the actual situation
*erosion_size_filter is the filter parameter used to generate the threshold basemap for morphological operation, and the value is not less than the number of pixels with the maximum radius of the object in the region
*erosion_size_compare is the filtering parameter used to generate the height difference calculation basemap, which can be appropriately increased for the image with blurred edges, generally 3~7 size
*pixel_threshold is the maximum number of feature pixels in the image area, which should be set according to the actual situation, but it should be noted that it is not possible to set the number of block cells to be processed in chunks, otherwise the processing of the image block will be skipped
*min_value is the segmentation threshold for judging the features, that is, the area where the average edge is how many meters above the surface is divided into feature areas
*iteration_value is an iterative binary iteration value, which can be adjusted according to the requirements, and should not be too high, otherwise the accuracy will be affected
*grey_open_size is a parameter for morphological operation of the interpolated image, which is used to smooth the interpolation area and adjust it according to the actual situation
*block_size is the size of the block of the image block operation, which can be adjusted according to the actual situation, for smaller images, it can be directly set to the block size larger than the image size, which can avoid blocking, if you need to divide the block, please choose a reasonable size to avoid the appearance of narrow strip blocks when blocking
*overlap_size is the size of the overlapping area between blocks, which can be set to 0 to cancel the overlapping area
*use_min_cut_for_overlap is used to set whether to enable the graph cutting algorithm, the overlapping areas without the graph cutting algorithm will be filled with the minimum value, set True on, False off
Original DSM
![image](https://github.com/user-attachments/assets/bb5c685a-4f58-4754-94bb-702b0f44a1e4)
DEM from IMF
![image](https://github.com/user-attachments/assets/a3b474bf-04f0-4327-8a3f-a8d03d159b75)
