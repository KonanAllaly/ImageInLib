#pragma warning(disable : 6011)
#pragma warning(disable : 4244)
#pragma warning(disable : 4477)
#pragma warning(disable : 4267)
#pragma warning(disable : 4700)

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string.h> 
#include <math.h>
#include <time.h>

#include "file.h"
#include "common_functions.h"
#include "Labelling.h"
#include "template_functions.h"
#include "../src/vtk_params.h"
#include "common_vtk.h"
#include "src/imageInterpolation.h"
#include "src/thresholding.h"
#include "morphological_change.h"
#include "distanceMaps.h"
#include "filtering.h"

#define _sx 1.171875
#define _sy 1.171875
#define _sz 2.5
#define _s_pet 4

#define thres_min 300
#define thres_max 700


int main() {

	size_t i, j, k;

	//image Dimensions
	size_t Width = 144;
	size_t Length = 144;
	const size_t dim2D = Width * Length;
	const size_t Height = 255;

	std::string inputPath = "C:/Users/Konan Allaly/Documents/Tests/input/";
	std::string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";

	dataType** imageData = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D];
	}
	if (imageData == NULL)
		return false;

	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			imageData[k][i] = 0.0;
		}
	}

	std::string inputImagePath = inputPath + "patient2_pet.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, inputImagePath.c_str(), LOAD_DATA, false);

	dataType minData = 100000.0, maxData = 0.0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (imageData[k][i] > maxData)
				maxData = imageData[k][i];

			if (imageData[k][i] < minData)
				minData = imageData[k][i];
		}
	}

	//std::cout << "max data = " << maxData << ", min data = " << minData << std::endl;

	rescaleNewRange(imageData, Length, Width, Height, 0, 1.0, maxData, minData);

	Image_Data F_Image;
	F_Image.imageDataPtr = imageData;
	F_Image.height = Height; F_Image.length = Length; F_Image.width = Width; 

	Filter_Parameters parameters = {1.0, 1.0, 1e-3, 1, 1.4, 1e-3, 1e-6, 1e-6, 1, 1, 100};

	//filterImage(F_Image, parameters, MEAN_CURVATURE_FILTER);
	filterImage(F_Image, parameters, GEODESIC_MEAN_CURVATURE_FILTER);

	std::string  outputImagePath = outputPath + "filtered_GMC_PET_p2_v2.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	//Image_Data image1;
	//image1.imageDataPtr = imageData;
	//image1.origin = { 286.585938, 216.585938, -1021.400024};
	//image1.height = Height; image1.length = Length; image1.width = Width;
	//image1.spacing = {_s_pet, _s_pet, _s_pet};
	//image1.orientation.v1 = { 1, 0, 0 }; image1.orientation.v2 = { 0, 1, 0 }; image1.orientation.v3 = { 0, 0, -1 };

	//const size_t Height_new = 406, Width_new = 512, Length_new = 512;
	//dataType** image = new dataType * [Height_new];
	//dataType** distanceMap = new dataType * [Height_new];
	//int** ball = new int* [Height_new];
	//int** segmentedImage = new int * [Height_new];
	//bool** statusArray = new bool* [Height_new];
	//for (k = 0; k < Height_new; k++) {
	//	image[k] = new dataType [Length_new * Width_new];
	//	distanceMap[k] = new dataType[Length_new * Width_new];
	//	ball[k] = new int[Length_new * Width_new];
	//	segmentedImage[k] = new int[Length_new * Width_new];
	//	statusArray[k] = new bool[Length_new * Width_new];
	//}
	//if (image == NULL || distanceMap == NULL || ball == NULL || segmentedImage == NULL || statusArray == NULL)
	//	return false;

	//for (k = 0; k < Height_new; k++) {
	//	for (i = 0; i < Length_new * Width_new; i++) {
	//		image[k][i] = 0.0;
	//		distanceMap[k][i] = 0.0;
	//		ball[k][i] = 0;
	//		segmentedImage[k][i] = 0;
	//		statusArray[k][i] = false;
	//	}
	//}

	//Image_Data image2;
	//image2.imageDataPtr = image;
	//image2.origin = {300.0, 230.0, -1022.50 };
	//image2.height = Height_new; image2.width = Width_new; image2.length = Length_new;
	//image2.spacing = { _sx, _sy, _sz };
	//image2.orientation.v1 = { 1, 0, 0 }; image2.orientation.v2 = {0, 1, 0}; image2.orientation.v3 = { 0, 0, -1 };

	//interpolationMethod method = TRILINEAR;
	//imageInterpolation3D(image1, image2, method);

	//rescaleNewRange(image, Length_new, Width_new, Height_new, 0, 4000, 40335, -0.3997);

	//thresholding3dFunctionN(image, Length_new, Width_new, Height_new, thres_min, thres_max, 0.0, 1.0);

	//Distance_Map_Params params_dist = {0.5, 1.0, 0.0, 1000000, 1e-3};
	//computeDistanceMap(distanceMap, image, Length_new, Width_new, Height_new, params_dist, FAST_SWEEP);

	////find the point with the max distance 
	//dataType dist_max = 0.0;
	//int imax = 0, jmax = 0, kmax = 0;
	//for (k = 0; k < Height_new; k++) {
	//	for (i = 0; i < Length_new; i++) {
	//		for (j = 0; j < Width_new; j++) {
	//			if (distanceMap[k][x_new(i, j, Length_new)] > dist_max) {
	//				dist_max = distanceMap[k][x_new(i, j, Length_new)];
	//				imax = i; jmax = j; kmax = k;
	//			}
	//		}
	//	}
	//}

	////Fill the ball around point with the highest distance
	//for (k = 0; k < Height_new; k++) {
	//	for (i = 0; i < Length_new; i++) {
	//		for (j = 0; j < Width_new; j++) {
	//			if (sqrt((i - imax) * (i - imax) + (j - jmax) * (j - jmax) + (k - kmax) * (k - kmax)) < 15) {
	//				ball[k][x_new(i, j, Length_new)] = 1;
	//			}
	//		}
	//	}
	//}

	//dilatation3dHeighteenNeigbours(image, Length_new, Width_new, Height_new, 1.0, 0.0);
	//erosion3dHeighteenNeigbours(image, Length_new, Width_new, Height_new, 1.0, 0.0);
	//dilatation3dHeighteenNeigbours(image, Length_new, Width_new, Height_new, 1.0, 0.0);
	//erosion3dHeighteenNeigbours(image, Length_new, Width_new, Height_new, 1.0, 0.0);
	//labelling3D(image, segmentedImage, statusArray, Length_new, Width_new, Height_new, 1.0);

	/*
	//number of region voxels
	int numberOfRegionsCells = 0;
	for (k = 0; k < Height_new; k++) {
		for (i = 0; i < Length_new * Width_new; i++) {
			if (image[k][i] == 1.0) {
				numberOfRegionsCells++;
			}
		}
	}
	//Counting
	int* countingArray = (int*)malloc(numberOfRegionsCells * sizeof(int));
	if (countingArray == NULL) return false;
	for (i = 0; i < numberOfRegionsCells; i++) 
		countingArray[i] = 0;

	for (k = 0; k < Height_new; k++) {
		for (i = 0; i < Length_new * Width_new; i++) {
			if (segmentedImage[k][i] > 0) {
				countingArray[segmentedImage[k][i]]++;
			}
		}
	}
	//Find the biggest region
	int maxElement = 0;
	for (k = 0; k < numberOfRegionsCells; k++) {
		if (countingArray[k] > maxElement) {
			maxElement = countingArray[k];
		}
	}
	//Keep just the biggest region 
	for (k = 0; k < Height_new; k++) {
		for (i = 0; i < Length_new * Width_new; i++) {
			if (countingArray[segmentedImage[k][i]] != maxElement) {
				segmentedImage[k][i] = 0;
			}
		}
	}
	*/

	//std::string  outputImagePath = outputPath + "distance_map.raw";
	//manageRAWFile3D<dataType>(distanceMap, Length_new, Width_new, Height_new, outputImagePath.c_str(), STORE_DATA, false);

	//outputImagePath = outputPath + "ball.raw";
	//manageRAWFile3D<int>(ball, Length_new, Width_new, Height_new, outputImagePath.c_str(), STORE_DATA, false);
	
	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
	}
	delete[] imageData;

	//for (k = 0; k < Height_new; k++) {
	//	delete[] image[k];
	//	delete[] segmentedImage[k];
	//	delete[] statusArray[k];
	//	delete[] distanceMap[k];
	//}
	//delete[] image;
	//delete[] segmentedImage;
	//delete[] statusArray;
	//delete[] distanceMap;

	//free(countingArray);

	return EXIT_SUCCESS;
}
