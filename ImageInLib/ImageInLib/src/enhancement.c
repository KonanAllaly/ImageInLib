#include "enhancement.h"

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