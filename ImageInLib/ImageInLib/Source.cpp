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
#include "src/imageInterpolation.h"
#include "segmentation.h"
#include "src/segmentation3D_subsurf.h"
#include "src/segmentation3d_gsubsurf.h"

#define _sx 1.171875
#define _sy 1.171875
#define _sz 2.5
#define _s_pet 4

#define thres_min 900
#define thres_max 1300


int main() {

	size_t i, j, k;

	std::string inputPath = "C:/Users/Konan Allaly/Documents/Tests/input/";
	std::string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";

	/*======== Data container for original ct data ============*/
	//const size_t Height = 406, Width = 512, Length = 512;
	//const size_t dim2D = Width * Length;
	//dataType** imageData = new dataType * [Height];
	//short** image = new short * [Height];
	//for (k = 0; k < Height; k++) {
	//	imageData[k] = new dataType[dim2D];
	//	image[k] = new short[dim2D];
	//}
	//if (imageData == NULL || image == NULL)
	//	return false;

	//Image_Data original_CT;
	//original_CT.imageDataPtr = imageData;
	//original_CT.height = Height; original_CT.length = Length; original_CT.width = Width;
	//original_CT.origin = { -300.0, -230.0, -1339.29999 };
	//original_CT.spacing = { _sx, _sy, _sz };
	//original_CT.orientation.v1 = { 1, 0, 0 }; original_CT.orientation.v2 = { 0, 1, 0 }; original_CT.orientation.v3 = { 0, 0, 1 };

	////Initialization
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		imageData[k][i] = 0.0;
	//		image[k][i] = 0;
	//	}
	//}

	//std::string inputImagePath = inputPath + "patient2_ct.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, inputImagePath.c_str(), LOAD_DATA, false);

	//std::string outputImagePath = outputPath + "orginal_ct.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	/*===================== Data container for orginal PET data =============*/

	//const size_t height_pet = 255, length_pet = 144, width_pet = 144;
	//dataType** image_pet = new dataType * [height_pet];
	//for (k = 0; k < height_pet; k++) {
	//	image_pet[k] = new dataType[length_pet * width_pet];
	//}
	//if (image_pet == NULL)
	//	return false;

	//Image_Data original_PET;
	//original_PET.imageDataPtr = image_pet;
	//original_PET.height = height_pet; original_PET.length = length_pet; original_PET.width = width_pet;
	//original_PET.origin = { -286.586, -216.586, -1338.80 };
	//original_PET.spacing = { _s_pet, _s_pet, _s_pet };
	//original_PET.orientation.v1 = { 1, 0, 0 }; original_PET.orientation.v2 = { 0, 1, 0 }; original_PET.orientation.v3 = { 0, 0, 1 };

	////Initialization
	//for (k = 0; k < height_pet; k++) {
	//	for (i = 0; i < length_pet * width_pet; i++) {
	//		image_pet[k][i] = 0.0;
	//	}
	//}
	//imageInterpolation3D(original_CT, original_PET, TRILINEAR);

	//outputImagePath = outputPath + "ct_interpolated_to_pet_dim.raw";
	//manageRAWFile3D<dataType>(image_pet, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	/*==================== Downsampling ================*/

	////Down sampled PET data structure
	//const size_t height = 127, length = 72, width = 72;
	//dataType** image_down = new dataType * [height];
	//for (k = 0; k < height; k++) {
	//	image_down[k] = new dataType[length * width];
	//}

	//Image_Data ImageD;
	//ImageD.imageDataPtr = image_down;
	//ImageD.height = height; ImageD.length = length; ImageD.width = width;
	//ImageD.origin = { -286.586, -216.586, -1338.80 };
	//ImageD.spacing = { 2.0, 2.0, 2.0 };
	//ImageD.orientation.v1 = { 1, 0, 0 }; ImageD.orientation.v2 = { 0, 1, 0 }; ImageD.orientation.v3 = { 0, 0, 1 };

	//original_PET.spacing = { 1.0, 1.0, 1.0 };
	//imageInterpolation3D(original_PET, ImageD, TRILINEAR);

	//outputImagePath = outputPath + "down_sampled.raw";
	//manageRAWFile3D<dataType>(image_down, length, width, height, outputImagePath.c_str(), STORE_DATA, false);

	//for (k = 0; k < Height; k++) {
	//	delete[] imageData[k];
	//	delete[] image[k];
	//}
	//delete[] imageData;
	//delete[] image;

	//for (k = 0; k < height_pet; k++) {
	//	delete[] image_pet[k];
	//}
	//delete[] image_pet;

	//for (k = 0; k < height; k++) {
	//	delete[] image_down[k];
	//}
	//delete[] image_down;

	/*=============== Work with down sampled ct image ================*/

	const size_t height = 127, length = 72, width = 72, dim2D = length * width;

	dataType** imageData = new dataType * [height];
	dataType** maskThreshold = new dataType * [height];
	dataType** distanceMap = new dataType * [height];
	for (k = 0; k < height; k++) {
		imageData[k] = new dataType[dim2D];
		maskThreshold[k] = new dataType[dim2D];
		distanceMap[k] = new dataType[dim2D];
	}
	if (imageData == NULL || maskThreshold == NULL || distanceMap == NULL)
		return false;

	//Initialization
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			imageData[k][i] = 0.0;
			maskThreshold[k][i] = 0.0;
			distanceMap[k][i] = 0.0;
		}
	}

	std::string inputImagePath = outputPath + "down_sampled.raw";
	manageRAWFile3D<dataType>(imageData, length, width, height, inputImagePath.c_str(), LOAD_DATA, false);

	//Copy and threshold
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (imageData[k][i] >= thres_min && imageData[k][i] <= thres_max) {
				maskThreshold[k][i] = 1.0;
			}
			else {
				maskThreshold[k][i] = 0.0;
			}
		}
	}

	//Distance Map
	Distance_Map_Params params_dist = {0.5, 1.0, 0.0, 1000000, 1e-3};
	computeDistanceMap(distanceMap, maskThreshold, length, width, height, params_dist, FAST_SWEEP);

	//find the point with the max distance 
	dataType dist_max = 0.0;
	int imax = 0, jmax = 0, kmax = 0;
	for (k = 0; k < height; k++) {
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {
				if (distanceMap[k][x_new(i, j, length)] > dist_max) {
					dist_max = distanceMap[k][x_new(i, j, length)];
					imax = i; jmax = j; kmax = k;
				}
			}
		}
	}
	std::cout << "Center ct (" << imax << ", " << jmax << ", " << kmax << ")" << std::endl;

	dataType minData = 100000.0, maxData = 0.0;
	for (k = 0; k < height; k++) {
		for (i = 0; i < length * width; i++) {
			if (imageData[k][i] > maxData)
				maxData = imageData[k][i];
			if (imageData[k][i] < minData)
				minData = imageData[k][i];
		}
	}
	rescaleNewRange(imageData, length, width, height, 0, 1.0, maxData, minData);

	/*===============================================================*/
	//Image_Data Image_CT;
	//Image_CT.imageDataPtr = imageData;
	//Image_CT.height = Height; Image_CT.length = Length; Image_CT.width = Width;
	//Image_CT.origin = { -300.0, -230.0, -1339.29999 };
	//Image_CT.spacing = { 1.171875, 1.171875, 2.5 };
	//Image_CT.orientation.v1 = { 1, 0, 0 }; Image_CT.orientation.v2 = { 0, 1, 0 }; Image_CT.orientation.v3 = { 0, 0, 1 };

	//dataType minData = 100000.0, maxData = 0.0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		if (imageData[k][i] > maxData)
	//			maxData = imageData[k][i];
	//		if (imageData[k][i] < minData)
	//			minData = imageData[k][i];
	//	}
	//}
	//std::cout << "max data = " << maxData << ", min data = " << minData << std::endl;
	//rescaleNewRange(imageData, Length, Width, Height, 0, 1.0, maxData, minData);

	//Image_Data F_Image;
	//F_Image.imageDataPtr = imageData;
	//F_Image.height = Height; F_Image.length = Length; F_Image.width = Width; 
	//Filter_Parameters parameters = {1.2, 1.0, 1e-3, 1, 1.4, 1e-3, 1e-6, 1e-6, 1, 1, 100};
	//filterImage(F_Image, parameters, MEAN_CURVATURE_FILTER);
	//filterImage(F_Image, parameters, GEODESIC_MEAN_CURVATURE_FILTER);
	//std::string  outputImagePath = outputPath + "filteredMC_CT_p2.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

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

	//thresholding3dFunctionN(maskThreshold, Length, Width, Height, thres_min, thres_max, 0.0, 1.0);
	//Distance_Map_Params params_dist = {0.5, 1.0, 0.0, 1000000, 1e-3};
	//computeDistanceMap(distanceMap, maskThreshold, Length, Width, Height, params_dist, FAST_SWEEP);

	////find the point with the max distance 
	//dataType dist_max = 0.0;
	//int imax = 0, jmax = 0, kmax = 0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (distanceMap[k][x_new(i, j, Length)] > dist_max) {
	//				dist_max = distanceMap[k][x_new(i, j, Length)];
	//				imax = i; jmax = j; kmax = k;
	//			}
	//		}
	//	}
	//}
	//std::cout << "Center ct (" << imax << ", " << jmax << ", " << kmax << ")" << std::endl;

	////Fill the ball around point with the highest distance
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (sqrt((i - imax) * (i - imax) + (j - jmax) * (j - jmax) + (k - kmax) * (k - kmax)) < 15) {
	//				ball[k][x_new(i, j, Length)] = 1;
	//			}
	//		}
	//	}
	//}

	////PET data structure
	//const size_t height_pet = 255, length_pet = 144, width_pet = 144;
	//dataType** image_pet = new dataType * [height_pet];
	//for (k = 0; k < height_pet; k++) {
	//	image_pet[k] = new dataType[length_pet * width_pet];
	//}
	//inputImagePath = outputPath + "filtered_MC_PET_p2_v2.raw";
	//manageRAWFile3D<dataType>(image_pet, length_pet, width_pet, height_pet, inputImagePath.c_str(), LOAD_DATA, false);

	//dataType minData = 100000.0, maxData = 0.0;
	//for (k = 0; k < height_pet; k++) {
	//	for (i = 0; i < length_pet * width_pet; i++) {
	//		if (image_pet[k][i] > maxData)
	//			maxData = image_pet[k][i];
	//		if (image_pet[k][i] < minData)
	//			minData = image_pet[k][i];
	//	}
	//}
	//rescaleNewRange(image_pet, length_pet, width_pet, height_pet, 0, 1.0, maxData, minData);

	//dataType minData = 100000.0, maxData = 0.0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length * Width; i++) {
	//		if (imageData[k][i] > maxData)
	//			maxData = imageData[k][i];
	//		if (imageData[k][i] < minData)
	//			minData = imageData[k][i];
	//	}
	//}
	//rescaleNewRange(imageData, Length, Width, Height, 0, 1.0, maxData, minData);

	//Image_Data Image_CT;
	//Image_CT.imageDataPtr = imageData;
	//Image_CT.height = Height; Image_CT.length = Length; Image_CT.width = Width;
	//Filter_Parameters filterParameters; 
	//filterParameters.coef = 1e-4; filterParameters.edge_detector_coefficient = 1000; filterParameters.eps2 = 1e-4;
	//filterParameters.h = 1.0; filterParameters.maxNumberOfSolverIteration = 100; filterParameters.omega_c = 1.5; filterParameters.p = 1;
	//filterParameters.sigma = 1e-3; filterParameters.timeStepSize = 0.25; filterParameters.timeStepsNum = 1; filterParameters.tolerance = 1e-3;
	////Filter_Parameters parameters = { 1.2, 1.0, 1e-3, 1, 1.4, 1e-3, 1e-6, 1e-6, 1, 1, 100 };
	////filterImage(Image_CT, parameters, MEAN_CURVATURE_FILTER);

	////filterImage(Image_CT, filterParameters, GEODESIC_MEAN_CURVATURE_FILTER);
	//std::string  outputImagePath = outputPath + "filteredMC_CT_p2.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, outputImagePath.c_str(), LOAD_DATA, false);

	//Image_Data F_Image;
	//F_Image.imageDataPtr = image_pet;
	//F_Image.height = height_pet; F_Image.length = length_pet; F_Image.width = width_pet;
	//Filter_Parameters parameters = { 1.2, 1.0, 1e-3, 1, 1.4, 1e-3, 1e-6, 1e-6, 1, 1, 100 };
	//filterImage(F_Image, parameters, MEAN_CURVATURE_FILTER);

	/*Point3D centerCT = { imax, jmax, kmax };
	Point3D  centerPET = { 0.0, 0.0, 0.0 };
	VoxelSpacing spacingCT = { _sx, _sy, _sz }, spacingPET = { _s_pet, _s_pet, _s_pet };
	Point3D originCT = { -300.0, -230.0, -1339.29999 }, originPET = { -286.585938, -216.585938, -1338.80 };
	OrientationMatrix orientation; orientation.v1 = { 1, 0, 0 }; orientation.v2 = { 0, 1, 0 }; orientation.v3 = { 0, 0, 1 };*/

	//centerPET = getRealCoordFromImageCoord3D(centerCT, originCT, spacingCT, orientation);
	//centerPET = getImageCoordFromRealCoord3D(centerPET, originPET, spacingPET, orientation);

	//int ipet = (int)centerPET.x, jpet = (int)centerPET.y, kpet = (int)centerPET.z;
	////dataType** ballPET = new dataType * [255];
	////for (k = 0; k < 255; k++) {
	////	ballPET[k] = new dataType[144 * 144];
	////}
	//////Fill the ball around point with the highest distance
	////for (k = 0; k < 255; k++) {
	////	for (i = 0; i < 144; i++) {
	////		for (j = 0; j < 144; j++) {
	////			if (sqrt((i - ipet) * (i - ipet) + (j - jpet) * (j - jpet) + (k - kpet) * (k - kpet)) < 5) {
	////				ballPET[k][x_new(i, j, 144)] = 1;
	////			}
	////			else {
	////				ballPET[k][x_new(i, j, 144)] = 0;
	////			}
	////		}
	////	}
	////}
	//std::cout << "Center pet : (" << ipet << ", " << jpet << ", " << kpet << ")" << std::endl;

	/*================= Interpolation of CT data Edge detector ====================*/

	////CT edge detector
	//dataType** edge_ct = new dataType * [Height];
	//for (k = 0; k < Height; k++) {
	//	edge_ct[k] = new dataType[dim2D];
	//}
	//Image_Data Image_CT;
	//Image_CT.imageDataPtr = edge_ct;
	//Image_CT.height = Height; Image_CT.length = Length; Image_CT.width = Width;
	//Image_CT.origin = originCT;
	//Image_CT.spacing = spacingCT;
	//Image_CT.orientation = orientation;

	//Load edge detector ct
	//std::string  outputImagePath = outputPath + "_edgeDetector_K=100000.raw";
	//manageRAWFile3D<dataType>(edge_ct, Length, Width, Height, outputImagePath.c_str(), LOAD_DATA, false);

	////PET edge detector coming from CT data
	//dataType** edge_pet = new dataType * [height_pet];
	//dataType** edge_ct = new dataType * [height_pet];
	//dataType** ct_50_pet_50 = new dataType * [height_pet];
	//dataType** ct_25_pet_75 = new dataType * [height_pet];
	//dataType** ct_75_pet_25 = new dataType * [height_pet];
	//for (k = 0; k < height_pet; k++) {
	//	edge_pet[k] = new dataType[width_pet * length_pet];
	//	edge_ct[k] = new dataType[width_pet * length_pet];
	//	ct_50_pet_50[k] = new dataType[width_pet * length_pet];
	//	ct_25_pet_75[k] = new dataType[width_pet * length_pet];
	//	ct_75_pet_25[k] = new dataType[width_pet * length_pet];
	//}
	//Image_Data Image_PET;
	//Image_PET.imageDataPtr = edge_pet;
	//Image_PET.height = height_pet; Image_PET.length = length_pet; Image_PET.width = width_pet;
	//Image_PET.origin = originPET;
	//Image_PET.spacing = spacingPET;
	//Image_PET.orientation = orientation;

	//std::string  outputImagePath = outputPath + "_edgeDetector_pet.raw";
	//manageRAWFile3D<dataType>(edge_pet, length_pet, width_pet, height_pet, outputImagePath.c_str(), LOAD_DATA, false);

	//outputImagePath = outputPath + "_edgeDetector_ct.raw";
	//manageRAWFile3D<dataType>(edge_ct, length_pet, width_pet, height_pet, outputImagePath.c_str(), LOAD_DATA, false);

	//for (k = 0; k < height_pet; k++) {
	//	for (i = 0; i < length_pet * width_pet; i++) {
	//		ct_50_pet_50[k][i] = 0.5 * edge_ct[k][i] + 0.5 * edge_pet[k][i];
	//		ct_25_pet_75[k][i] = 0.25 * edge_ct[k][i] + 0.75 * edge_pet[k][i];
	//		ct_75_pet_25[k][i] = 0.75 * edge_ct[k][i] + 0.25 * edge_pet[k][i];
	//	}
	//}

	//outputImagePath = outputPath + "_edge_ct_50_pet_50.raw";
	//manageRAWFile3D<dataType>(ct_50_pet_50, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);
	//outputImagePath = outputPath + "_edge_ct_75_pet_25.raw";
	//manageRAWFile3D<dataType>(ct_75_pet_25, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);
	//outputImagePath = outputPath + "_edge_ct_25_pet_75.raw";
	//manageRAWFile3D<dataType>(ct_50_pet_50, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);
	//imageInterpolation3D(Image_CT, Image_PET, TRILINEAR);

	////imageInterpolation3D(Image_CT, Image_PET, NEAREST_NEIGHBOR);
	//std::string outputImagePath = outputPath + "interpolated_ct_trilinear.raw";
	//manageRAWFile3D<dataType>(edge_pet, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	//for (k = 0; k < height_ct; k++) {
	//	delete[] edge_ct[k];
	//}
	//delete[] edge_ct;

	//for (k = 0; k < height_pet; k++) {
	//	delete[] edge_pet[k];
	//}
	//delete[] edge_pet;

	//======================  Labelling ========================

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
	//manageRAWFile3D<dataType>(distanceMap, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);
	//outputImagePath = outputPath + "ballPET.raw";
	//manageRAWFile3D<dataType>(ballPET, 144, 144, 255, outputImagePath.c_str(), STORE_DATA, false);

	//FilterParameters filtering_parameters;
	//filtering_parameters.tau = 0.25; filtering_parameters.h = 1.0;
	//filtering_parameters.omega_c = 1.4; filtering_parameters.tolerance = 1e-3;
	//filtering_parameters.maxNumberOfSolverIteration = 1000;
	//filtering_parameters.eps2 = 0.000001; filtering_parameters.edge_detector_coefficient = 1000;

	//SegmentationParameters segmentation_parameters;
	//segmentation_parameters.h = 1.0; segmentation_parameters.edge_detector_coefficient = 100; segmentation_parameters.eps2 = 0.0001;
	//segmentation_parameters.maxNumberOfSolverIteration = 100; segmentation_parameters.omega_c = 1.5; segmentation_parameters.tolerance = 1e-3;
	//segmentation_parameters.tau = 10.0; segmentation_parameters.timeStepsNum = 1000;

	//free(countingArray);

	//============= Segmentation ================

	Point3D* center_seg = new Point3D[1];
	//center_seg->x = ipet; center_seg->y = jpet; center_seg->z = kpet;
	center_seg->x = imax; center_seg->y = jmax; center_seg->z = kmax;
	dataType** initial_seg = new dataType * [height];
	for (k = 0; k < height; k++) {
		initial_seg[k] = new dataType[dim2D];
	}
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			initial_seg[k][i] = 0.0;
		}
	}
	generateInitialSegmentationFunctionForMultipleCentres(initial_seg, length, width, height, center_seg, 1.0, 10, 1);

	//Parameters for CT
	//Segmentation_Parameters seg_params; seg_params.coef = 10000; seg_params.eps2 = 1e-6; seg_params.gauss_seidelTolerance = 1e-4;
	//seg_params.h = 1.0; seg_params.maxNoGSIteration = 1000; seg_params.maxNoOfTimeSteps = 500; seg_params.mod = 1;
	//seg_params.numberOfTimeStep = 500; seg_params.omega_c = 1.5; seg_params.segTolerance = 1e-4;
	//dataType Tmax = seg_params.maxNoOfTimeSteps; seg_params.tau = 1.0;
	//Filter_Parameters parameters; parameters.coef = 1e-4; parameters.edge_detector_coefficient = 1000; parameters.eps2 = 1e-4;
	//parameters.h = 1.0; parameters.maxNumberOfSolverIteration = 100; parameters.omega_c = 1.5; parameters.p = 1;
	//parameters.sigma = 1e-4; parameters.timeStepSize = 1.2; parameters.timeStepsNum = 1; parameters.tolerance = 1e-7;

	//Segmentation_Parameters seg_params; 
	//seg_params.coef = 100000; seg_params.eps2 = 1e-6; seg_params.gauss_seidelTolerance = 1e-6;
	//seg_params.h = 1.0; seg_params.maxNoGSIteration = 100; seg_params.maxNoOfTimeSteps = 1; seg_params.mod = 10;
	//seg_params.numberOfTimeStep = 1; seg_params.omega_c = 1.5; seg_params.segTolerance = 1e-10; seg_params.tau = 8.0;
	//seg_params.coef_conv = 1.0; seg_params.coef_dif = 1.0;

	//Filter_Parameters parameters;
	//parameters.coef = 1e-2; parameters.edge_detector_coefficient = 1; parameters.eps2 = 1e-4;
	//parameters.h = 1.0; parameters.maxNumberOfSolverIteration = 100; parameters.omega_c = 1.5; parameters.p = 1;
	//parameters.sigma = 1e-3; parameters.timeStepSize = 0.25; parameters.timeStepsNum = 1; parameters.tolerance = 1e-3;

	//Filter_Parameters parameters;
	//parameters.timeStepSize = 0.25; parameters.h = 1.0;
	//parameters.omega_c = 1.4; parameters.tolerance = 1e-3;
	//parameters.maxNumberOfSolverIteration = 1000;
	//parameters.eps2 = 1e-6; parameters.edge_detector_coefficient = 1000;

	//Segmentation_Parameters seg_params;
	//seg_params.h = 1.0; seg_params.coef = 100; seg_params.eps2 = 0.0001;
	//seg_params.maxNoGSIteration = 10000; seg_params.omega_c = 1.5; seg_params.gauss_seidelTolerance = 1e-6;
	//seg_params.tau = 10.0; seg_params.maxNoOfTimeSteps = 1000; seg_params.mod = 1;
	//seg_params.coef_conv = 1.0; seg_params.coef_dif = 1.0;

	Segmentation_Parameters seg_params;
	seg_params.tau = 1.0; seg_params.h = 1.0; seg_params.omega_c = 1.4; seg_params.mod = 1;
	seg_params.maxNoGSIteration = 100; seg_params.coef = 10000; seg_params.gauss_seidelTolerance = 1e-6;
	seg_params.maxNoOfTimeSteps = 1000; seg_params.segTolerance = 1e-6; seg_params.eps2 = 1e-6;
	seg_params.coef_conv = 1.0; seg_params.coef_dif = 1.0;

	Filter_Parameters parameters; parameters.coef = 1e-4; parameters.edge_detector_coefficient = 1000; parameters.eps2 = 1e-4;
	parameters.h = 1.0; parameters.maxNumberOfSolverIteration = 100; parameters.omega_c = 1.5; parameters.p = 1;
	parameters.sigma = 1e-4; parameters.timeStepSize = 1.2; parameters.timeStepsNum = 1; parameters.tolerance = 1e-7;

	unsigned char segmentPath[] = "C:/Users/Konan Allaly/Documents/Tests/output/segmentation/";

	Image_Data Image;
	Image.imageDataPtr = imageData;
	Image.length = length; Image.width = width; Image.height = height;

	////subsurfSegmentation(F_Image, initial_seg, seg_params, parameters, center_seg, 1, segmentPath);
	generalizedSubsurfSegmentation(Image, initial_seg, seg_params, parameters, center_seg, 1, segmentPath);
	
	//for (k = 0; k < Height; k++) {
	//	delete[] imageData[k];
	//	delete[] image[k];
	//	delete[] distanceMap[k];
	//	delete[] ball[k];
	//	delete[] maskThreshold[k];
	//	delete[] initial[k];
	//	//delete[] edge_ct[k];
	//}
	//delete[] imageData;
	//delete[] image;
	//delete[] distanceMap;
	//delete[] ball;
	//delete[] maskThreshold;
	//delete[] initial;
	////delete[] edge_ct;

	////for (k = 0; k < height_pet; k++) {
	////	delete[] image_pet[k];
	////	delete[] edge_pet[k];
	////	delete[] edge_ct[k];
	////	delete[] ct_50_pet_50[k];
	////	delete[] ct_25_pet_75[k];
	////	delete[] ct_75_pet_25[k];
	////}
	////delete[] image_pet;
	////delete[] edge_pet;
	////delete[] edge_ct;
	////delete[] ct_50_pet_50;
	////delete[] ct_25_pet_75;
	////delete[] ct_75_pet_25;

	//for (k = 0; k < height_pet; k++) {
	//	delete[] initial_seg[k];
	//}
	//delete[] initial_seg;

	for (k = 0; k < height; k++) {
		delete[] imageData[k];
		delete[] maskThreshold[k];
		delete[] distanceMap[k];
		delete[] initial_seg[k];
	}
	delete[] imageData;
	delete[] maskThreshold;
	delete[] distanceMap;
	delete[] initial_seg;
	delete[] center_seg;

	//for (k = 0; k < 255; k++) {
	//	delete[] ballPET[k];
	//}
	//delete[] ballPET;

	return EXIT_SUCCESS;
}
