#include "enhancement.h"
#include "thresholding.h"

extremeValues findExtremeValues(dataType* imageDataPtr, const size_t height, const size_t width, const BoundingBox2D box) {
	
	dataType minData = INFINITY, maxData = -1 * INFINITY;
	for (size_t i = box.i_min; i <= box.i_max; i++) {
		for (size_t j = box.j_min; j <= box.j_max; j++) {
			size_t xd = x_new(i, j, height);
			if (imageDataPtr[xd] > maxData) {
				maxData = imageDataPtr[xd];
			}
			if (imageDataPtr[xd] < minData) {
				minData = imageDataPtr[xd];
			}
		}
	}
	return (extremeValues) { minData, maxData };
}

dataType getImageMeanValue(dataType* imageDataPtr, const size_t length, const size_t width) {
	if (imageDataPtr == NULL)
		return;

	size_t dim2D = length * width;
	dataType mean_value = 0.0;
	for (size_t i = 0; i < dim2D; i++) {
		mean_value += imageDataPtr[i];
	}
	return (mean_value / (dataType)dim2D);
}

dataType getImageMeanValue3D(dataType** imageDataPtr, const size_t length, const size_t width, const size_t height) {
	if (imageDataPtr == NULL)
		return;

	size_t dim2D = length * width;
	dataType mean_value = 0.0;
	for (size_t k = 0; k < height; k++) {
		for (size_t i = 0; i < dim2D; i++) {
			mean_value += imageDataPtr[k][i];
		}
	}
	return (mean_value / (dataType)(dim2D * height));
}

bool computeImageHistogram(dataType* imageDataPtr, const size_t height, const size_t width, const size_t bins, size_t* histogram) {

	if (imageDataPtr == NULL || histogram == NULL)
		return false;

	BoundingBox2D box = { 0, height - 1, 0, width - 1 };
	extremeValues MinMax = findExtremeValues(imageDataPtr, height, width, box);
	dataType step = (MinMax.max - MinMax.min) / (dataType)bins;

	size_t count_elements;
	for (size_t n = 0; n < bins; n++) {
		dataType min_bin = n * step;
		dataType max_bin = (n + 1) * step;
		count_elements = 0;
		for (size_t i = 0; i < height; i++) {
			for (size_t j = 0; j < width; j++) {
				size_t xd = x_new(i, j, height);
				if (imageDataPtr[xd] >= min_bin && imageDataPtr[xd] <= max_bin) {
					count_elements++;
				}
			}
		}
		histogram[n] = count_elements;
	}
	return true;
}

size_t findHistogramPeaks(dataType* imageDataPtr, const size_t height, const size_t width, const size_t bins, const char* store_file) {
	if (imageDataPtr == NULL || bins == 0 || bins > 256)
		return false;

	size_t* histogram = (size_t*)malloc(bins * sizeof(size_t));
	for (size_t i = 0; i < bins; i++) {
		histogram[i] = 0;
	}

	BoundingBox2D box = { 0, height - 1, 0, width - 1 };
	extremeValues MinMax = findExtremeValues(imageDataPtr, height, width, box);
	dataType step = (MinMax.max - MinMax.min) / (dataType)bins;

	computeImageHistogram(imageDataPtr, height, width, bins, histogram);
	//FILE* file_histogram;
	//if (fopen_s(&file_histogram, store_file, "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(file_histogram, "x,y\n");
	//for (size_t n = 0; n < bins; n++) {
	//	dataType inf_bin = n * step;
	//	fprintf(file_histogram, "%f,%d\n", inf_bin, histogram[n]);
	//}
	//fclose(file_histogram);

	//const char storing_path[] = "C:/Users/Konan Allaly/Documents/Tests/output/histogram/peaks_slice_240.csv";
	//FILE* file_peaks;
	//if (fopen_s(&file_peaks, storing_path, "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(file_peaks, "x,y\n");

	size_t count_peaks = 0;
	for (size_t i = 1; i < bins - 1; i++) {
		dataType frequence = histogram[i] / ((dataType)(height * width));
		if (frequence > 0.005 && (histogram[i] > histogram[i - 1]) && (histogram[i] > histogram[i + 1])) {
			//dataType inf_range = i * step;
			//dataType sup_range = (i + 1) * step;
			//printf("Peak %d: %f to %f\n", index_peak, inf_range, sup_range);

			//dataType inf_bin = i * step;
			//fprintf(file_peaks, "%f,%d\n", inf_bin, histogram[i]);
			count_peaks++;
		}
	}
	//printf("we found : %d peaks\n", count_peaks);
	//fclose(file_peaks);

	free(histogram);
	return count_peaks;
}

size_t histogramPeaks3D(dataType** imageDataPtr, const size_t length, const size_t width, const size_t height, const size_t bins, const char* store_file) {
	if (imageDataPtr == NULL || bins == 0 || bins > 256)
		return false;

	size_t* histogram = (size_t*)malloc(bins * sizeof(size_t));
	for (size_t i = 0; i < bins; i++) {
		histogram[i] = 0;
	}

	//Find minimum and maximum values
	dataType minData = imageDataPtr[0][0];
	dataType maxData = imageDataPtr[0][0];
	for (size_t i = 0; i < height; i++) {
		for (size_t j = 0; j < length * width; j++) {
			if (imageDataPtr[i][j] > maxData) {
				maxData = imageDataPtr[i][j];
			}
			if (imageDataPtr[i][j] < minData) {
				minData = imageDataPtr[i][j];
			}
		}
	}
	dataType step = (maxData - minData) / (dataType)bins;

	computeHistogram(imageDataPtr, histogram, length, width, height, bins);
	//FILE* file_histogram;
	//if (fopen_s(&file_histogram, store_file, "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(file_histogram, "x,y\n");
	//for (size_t n = 0; n < bins; n++) {
	//	dataType inf_bin = n * step;
	//	fprintf(file_histogram, "%f,%d\n", inf_bin, histogram[n]);
	//}
	//fclose(file_histogram);

	//const char storing_path[] = "C:/Users/Konan Allaly/Documents/Tests/output/removed/peaks_modified.csv";
	//FILE* file_peaks;
	//if (fopen_s(&file_peaks, storing_path, "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(file_peaks, "x,y\n");

	size_t count_peaks = 0;
	size_t index_peak = 1;
	for (size_t i = 1; i < bins - 1; i++) {
		dataType frequence = histogram[i] / ((dataType)(height * width * length));
		if (frequence > 0.005 && (histogram[i] > histogram[i - 1]) && (histogram[i] > histogram[i + 1])) {
			dataType inf_bin = i * step;
			dataType sup_bin = (i + 1) * step;
			printf("Peak %d: %f to %f\n", index_peak, inf_bin, sup_bin);

			//dataType inf_bin = i * step;
			//fprintf(file_peaks, "%f,%d\n", inf_bin, histogram[i]);
			count_peaks++;
			index_peak++;
		}
	}
	//printf("we found : %d peaks\n", count_peaks);
	//fclose(file_peaks);

	free(histogram);
	return count_peaks;
}

bool histogramEqualization(dataType* imageDataPtr, const size_t height, const size_t width) {
	
	if (imageDataPtr == NULL)
		return false;

	return true;
}

bool CLACHE(dataType* imageDataPtr, const size_t height, const size_t width, const size_t nbDivision) {
	
	if (imageDataPtr == NULL)
		return false;

	//Divide the image into 64 equal regions
	size_t h = height / nbDivision;
	size_t w = width / nbDivision;

	size_t idRegion = 0;
	for (size_t n = 0; n < nbDivision; n++) {
		for (size_t m = 0; m < nbDivision; m++) {
			size_t i_min = n * h;
			size_t i_max = (n + 1) * h - 1;
			size_t j_min = m * w;
			size_t j_max = (m + 1) * w - 1;
			BoundingBox2D box = {i_min, i_max, j_min, j_max};
			idRegion++;

			//bool computeImageHistogram(dataType * imageDataPtr, const size_t height, const size_t width, const size_t bins, size_t * histogram)
			
			//if (idRegion != 36 && idRegion != 37 && idRegion != 28 && idRegion != 29) {
			//	for (size_t i = i_min; i <= i_max; i++) {
			//		for (size_t j = j_min; j <= j_max; j++) {
			//			imageDataPtr[x_new(i, j, height)] = 0;
			//		}
			//	}
			//}

			extremeValues MinMax = findExtremeValues(imageDataPtr, height, width, box);
			printf("Pixel range region %d : %f\n", idRegion, MinMax.max - MinMax.min);

		}
	}

	return true;
}