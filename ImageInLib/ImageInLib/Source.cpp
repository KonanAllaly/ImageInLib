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

	Point3D originCT = { -300.0, -230.0, -1339.29999 };
	Point3D originPET = { -286.586, -216.586, -1338.80 };
	VoxelSpacing spacingCT = { _sx, _sy, _sz };
	VoxelSpacing spacingPET = { _s_pet, _s_pet, _s_pet };
	OrientationMatrix orientation;
	orientation.v1 = { 1.0, 0.0, 0.0 }; orientation.v2 = { 0.0, 1.0, 0.0 }; orientation.v3 = { 0.0, 0.0, 1.0 };

	/*======== Data container for original ct data ============*/
	const size_t Height = 406, Width = 512, Length = 512;
	const size_t dim2D = Width * Length;
	dataType** image_ct = new dataType * [Height];
	dataType** edge_ct = new dataType * [Height];
	dataType** image_ct_filtered = new dataType * [Height];
	dataType** maskThreshold = new dataType * [Height];
	dataType** distanceMap = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		image_ct[k] = new dataType[dim2D];
		edge_ct[k] = new dataType[dim2D];
		image_ct_filtered[k] = new dataType[dim2D];
		maskThreshold[k] = new dataType[dim2D];
		distanceMap[k] = new dataType[dim2D];
	}
	if (image_ct == NULL || edge_ct == NULL || image_ct_filtered == NULL || maskThreshold == NULL || distanceMap == NULL)
		return false;

	//Initialization
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			image_ct[k][i] = 0.0;
			edge_ct[k][i] = 0.0;
			image_ct_filtered[k][i] = 0.0;
			maskThreshold[k][i] = 0.0;
			distanceMap[k][i] = 0.0;
		}
	}

	Image_Data IMAGE_CT;
	IMAGE_CT.imageDataPtr = image_ct_filtered;
	IMAGE_CT.height = Height; IMAGE_CT.length = Length; IMAGE_CT.width = Width;
	IMAGE_CT.origin = originCT;
	IMAGE_CT.spacing = spacingCT;
	IMAGE_CT.orientation = orientation;

	std::string inputImagePath = inputPath + "patient2_ct.raw";
	manageRAWFile3D<dataType>(image_ct, Length, Width, Height, inputImagePath.c_str(), LOAD_DATA, false);

	inputImagePath = inputPath + "filtered_patient2_ct.raw";
	manageRAWFile3D<dataType>(image_ct_filtered, Length, Width, Height, inputImagePath.c_str(), LOAD_DATA, false);

	//std::string outputImagePath = outputPath + "orginal_ct.raw";
	//manageRAWFile3D<dataType>(image_ct, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	/*===================   Find the point with the higest distance ==========*/

	//fill the threshold container
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			maskThreshold[k][i] = image_ct[k][i];
		}
	}

	thresholding3dFunctionN(maskThreshold, Length, Width, Height, thres_min, thres_max, 0.0, 1.0);

	//std::string outputImagePath = outputPath + "threshold.raw";
	//manageRAWFile3D<dataType>(maskThreshold, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	//Distance Map
	Distance_Map_Params params_dist = { 0.5, 1.0, 0.0, 1000000, 1e-3 };
	computeDistanceMap(distanceMap, maskThreshold, Length, Width, Height, params_dist, FAST_SWEEP);

	//outputImagePath = outputPath + "distance_map.raw";
	//manageRAWFile3D<dataType>(distanceMap, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	//find the point with the max distance 
	dataType dist_max = 0.0;
	int imax = 0, jmax = 0, kmax = 0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				if (distanceMap[k][x_new(i, j, Length)] > dist_max) {
					dist_max = distanceMap[k][x_new(i, j, Length)];
					imax = i; jmax = j; kmax = k;
				}
			}
		}
	}
	std::cout << "Center ct (" << imax << ", " << jmax << ", " << kmax << ")" << std::endl;

	/*===================== Data container for orginal PET data =============*/

	const size_t height_pet = 255, length_pet = 144, width_pet = 144;
	dataType** image_pet = new dataType * [height_pet];
	dataType** edge_pet = new dataType * [height_pet];
	for (k = 0; k < height_pet; k++) {
		image_pet[k] = new dataType[length_pet * width_pet];
		edge_pet[k] = new dataType[length_pet * width_pet];
	}
	if (image_pet == NULL || edge_pet == NULL)
		return false;

	for (k = 0; k < height_pet; k++) {
		for (i = 0; i < length_pet * width_pet; i++) {
			image_pet[k][i] = 0.0;
			edge_pet[k][i] = 0.0;
		}
	}

	Image_Data IMAGE_PET;
	IMAGE_PET.imageDataPtr = image_pet;
	IMAGE_PET.height = height_pet; IMAGE_PET.length = length_pet; IMAGE_PET.width = width_pet;
	IMAGE_PET.origin = originPET;
	IMAGE_PET.spacing = spacingPET;
	IMAGE_PET.orientation = orientation;

	inputImagePath = inputPath + "filtered_patient2_pet.raw";
	manageRAWFile3D<dataType>(image_pet, length_pet, width_pet, height_pet, inputImagePath.c_str(), LOAD_DATA, false);

	/*==================  Interpolated Coordinates ========*/

	Point3D centerCT = { imax, jmax, kmax }, centerPET = { 0.0, 0.0, 0.0 };

	centerPET = getRealCoordFromImageCoord3D(centerCT, originCT, spacingCT, orientation);
	centerPET = getImageCoordFromRealCoord3D(centerPET, originPET, spacingPET, orientation);

	int ipet = (int)centerPET.x, jpet = (int)centerPET.y, kpet = (int)centerPET.z;
	std::cout << "Center pet (" << ipet << ", " << jpet << ", " << kpet << ")" << std::endl;

	/*================= Interpolation input Image and Edge detector ====================*/

	//Load edge detector ct
	inputImagePath = inputPath + "edge_detector_ct.raw";
	manageRAWFile3D<dataType>(edge_ct, Length, Width, Height, inputImagePath.c_str(), LOAD_DATA, false);

	//Load edge detector pet
	inputImagePath = inputPath + "edge_detector_pet.raw";
	manageRAWFile3D<dataType>(edge_pet, length_pet, width_pet, height_pet, inputImagePath.c_str(), LOAD_DATA, false);

	Image_Data EDGE_CT;
	EDGE_CT.imageDataPtr = edge_ct;
	EDGE_CT.height = Height; EDGE_CT.length = Length; EDGE_CT.width = Width;
	EDGE_CT.origin = originCT;
	EDGE_CT.spacing = spacingCT;
	EDGE_CT.orientation = orientation;

	Image_Data EDGE_PET;
	EDGE_PET.imageDataPtr = edge_pet;
	EDGE_PET.height = height_pet; EDGE_PET.length = length_pet; EDGE_PET.width = width_pet;
	EDGE_PET.origin = originPET;
	EDGE_PET.spacing = spacingPET;
	EDGE_PET.orientation = orientation;

	//Interpolate From CT dim to PET dim
	dataType** ct_interpolated = new dataType * [height_pet];
	dataType** edge_ct_interpolated = new dataType * [height_pet];
	for (k = 0; k < height_pet; k++) {
		ct_interpolated[k] = new dataType[length_pet * width_pet];
		edge_ct_interpolated[k] = new dataType[length_pet * width_pet];
	}
	if (ct_interpolated == NULL || edge_ct_interpolated == NULL)
		return false;

	//Initialization
	for (k = 0; k < height_pet; k++) {
		for (i = 0; i < length_pet * width_pet; i++) {
			ct_interpolated[k][i] = 0.0;
			edge_ct_interpolated[k][i] = 0.0;
		}
	}

	Image_Data CT_interpolate;
	CT_interpolate.imageDataPtr = ct_interpolated;
	CT_interpolate.height = height_pet; CT_interpolate.length = length_pet; CT_interpolate.width = width_pet;
	CT_interpolate.origin = originPET; CT_interpolate.spacing = spacingPET; CT_interpolate.orientation = orientation;

	Image_Data edge_CT_interpolate;
	edge_CT_interpolate.imageDataPtr = edge_ct_interpolated;
	edge_CT_interpolate.height = height_pet; edge_CT_interpolate.length = length_pet; edge_CT_interpolate.width = width_pet;
	edge_CT_interpolate.origin = originPET; edge_CT_interpolate.spacing = spacingPET; edge_CT_interpolate.orientation = orientation;

	imageInterpolation3D(IMAGE_CT, CT_interpolate, TRILINEAR);
	imageInterpolation3D(EDGE_CT, edge_CT_interpolate, TRILINEAR);

	//outputImagePath = outputPath + "edge_ct_in_pet_dim.raw";
	//manageRAWFile3D<dataType>(edge_ct_interpolated, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	//std::string  outputImagePath = outputPath + "interpolated_ct.raw";
	//manageRAWFile3D<dataType>(ct_interpolated, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	//outputImagePath = outputPath + "_edge_ct_interpolated.raw";
	//manageRAWFile3D<dataType>(edge_ct_interpolated, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	/*=============== Compute new edge detector ===========*/

	dataType** edge_ct_50_pet_50 = new dataType * [height_pet];
	dataType** edge_ct_75_pet_25 = new dataType * [height_pet];
	dataType** edge_ct_25_pet_75 = new dataType * [height_pet];
	for (k = 0; k < height_pet; k++) {
		edge_ct_50_pet_50[k] = new dataType[length_pet * width_pet];
		edge_ct_75_pet_25[k] = new dataType[length_pet * width_pet];
		edge_ct_25_pet_75[k] = new dataType[length_pet * width_pet];
	}
	if (edge_ct_50_pet_50 == NULL || edge_ct_75_pet_25 == NULL || edge_ct_25_pet_75 == NULL)
		return false;

	//Initialization
	for (k = 0; k < height_pet; k++) {
		for (i = 0; i < length_pet * width_pet; i++) {
			edge_ct_50_pet_50[k][i] = 0.0;
			edge_ct_75_pet_25[k][i] = 0.0;
			edge_ct_25_pet_75[k][i] = 0.0;
		}
	}

	for (k = 0; k < height_pet; k++) {
		for (i = 0; i < length_pet * width_pet; i++) {
			edge_ct_50_pet_50[k][i] = 0.5 * edge_ct_interpolated[k][i] + 0.5 * edge_pet[k][i];
			edge_ct_75_pet_25[k][i] = 0.75 * edge_ct_interpolated[k][i] + 0.25 * edge_pet[k][i];
			edge_ct_25_pet_75[k][i] = 0.25 * edge_ct_interpolated[k][i] + 0.75 * edge_pet[k][i];
		}
	}

	//outputImagePath = outputPath + "edge_ct_50_pet_50.raw";
	//manageRAWFile3D<dataType>(edge_ct_50_pet_50, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	//outputImagePath = outputPath + "edge_ct_75_pet_25.raw";
	//manageRAWFile3D<dataType>(edge_ct_75_pet_25, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	//outputImagePath = outputPath + "edge_ct_25_pet_75.raw";
	//manageRAWFile3D<dataType>(edge_ct_25_pet_75, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	/*==================== Downsampling ================*/

	//const size_t height = 127, length = 72, width = 72;
	const size_t height = 170, length = 96, width = 96;
	dataType** image_down = new dataType * [height];
	for (k = 0; k < height; k++) {
		image_down[k] = new dataType[length * width];
	}

	//down sample edge ct
	edge_CT_interpolate.origin = originPET;
	edge_CT_interpolate.spacing = { 1.0, 1.0, 1.0 };

	Image_Data ImageD;
	ImageD.imageDataPtr = image_down;
	ImageD.height = height; ImageD.length = length; ImageD.width = width;
	ImageD.origin = originPET;
	//ImageD.spacing = { 2.0, 2.0, 2.0 };
	ImageD.spacing = { 1.5, 1.5, 1.5 };
	ImageD.orientation = orientation;

	for (k = 0; k < height; k++) {
		for (i = 0; i < length * width; i++) {
			image_down[k][i] = 0.0;
		}
	}

	imageInterpolation3D(edge_CT_interpolate, ImageD, TRILINEAR);

	std::string outputImagePath = outputPath + "edge_ct_down_sampled_72_127.raw";
	manageRAWFile3D<dataType>(image_down, length, width, height, outputImagePath.c_str(), STORE_DATA, false);

	Point3D centerDown = { 0.0, 0.0, 0.0 };

	centerDown = getRealCoordFromImageCoord3D(centerPET, originPET, edge_CT_interpolate.spacing, orientation);
	centerDown = getImageCoordFromRealCoord3D(centerDown, ImageD.origin, ImageD.spacing, orientation);

	int idown = (int)centerDown.x, jdown = (int)centerDown.y, kdown = (int)centerDown.z;

	std::cout << "Center down sampled data (" << idown << ", " << jdown << ", " << kdown << ")" << std::endl;

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
	center_seg->x = idown; center_seg->y = jdown; center_seg->z = kdown;
	//center_seg->x = ipet; center_seg->y = jpet; center_seg->z = kpet;
	//center_seg->x = imax; center_seg->y = jmax; center_seg->z = kmax;
	 
	/*
	dataType** initial_seg = new dataType * [height_pet];
	for (k = 0; k < height_pet; k++) {
		initial_seg[k] = new dataType[length_pet * width_pet];
	}
	for (k = 0; k < height_pet; k++) {
		for (i = 0; i < width_pet * length_pet; i++) {
			initial_seg[k][i] = 0.0;
		}
	}
	*/

	dataType** initial_seg = new dataType * [height];
	for (k = 0; k < height; k++) {
		initial_seg[k] = new dataType[length * width];
	}
	for (k = 0; k < height; k++) {
		for (i = 0; i < width * length; i++) {
			initial_seg[k][i] = 0.0;
		}
	}

	generateInitialSegmentationFunctionForMultipleCentres(initial_seg, length, width, height, center_seg, 1.0, 10.0, 1.0);

	Image_Data ImageToBeSegmented;
	ImageToBeSegmented.imageDataPtr = image_down;
	ImageToBeSegmented.height = height; ImageToBeSegmented.length = length; ImageToBeSegmented.width = width;

	Segmentation_Parameters seg_params;
	seg_params.tau = 1.0; seg_params.h = 6.0; seg_params.omega_c = 1.4; seg_params.mod = 10;
	seg_params.maxNoGSIteration = 100; seg_params.coef = 10000; seg_params.gauss_seidelTolerance = 1e-6;
	seg_params.maxNoOfTimeSteps = 12000; seg_params.segTolerance = 1e-6; seg_params.eps2 = 1e-6;
	seg_params.coef_conv = 1.0; seg_params.coef_dif = 1.0;

	Filter_Parameters parameters; parameters.coef = 1e-4; parameters.edge_detector_coefficient = 1000; parameters.eps2 = 1e-4;
	parameters.h = 1.0; parameters.maxNumberOfSolverIteration = 100; parameters.omega_c = 1.4; parameters.p = 1;
	parameters.sigma = 0.1; parameters.timeStepSize = 1.0; parameters.timeStepsNum = 1; parameters.tolerance = 1e-5;

	unsigned char segmentPath[] = "C:/Users/Konan Allaly/Documents/Tests/output/segmentation/";

	//Image_Data Image;
	//Image.imageDataPtr = image_ct_filtered;
	//Image.length = Length; Image.width = Width; Image.height = Height;

	////subsurfSegmentation(F_Image, initial_seg, seg_params, parameters, center_seg, 1, segmentPath);
	generalizedSubsurfSegmentation(ImageToBeSegmented, initial_seg, seg_params, parameters, center_seg, 1, segmentPath);

	for (k = 0; k < height; k++) {
		delete[] image_down[k];
		delete[] initial_seg[k];
	}
	delete[] image_down;
	delete[] initial_seg;

	for (k = 0; k < Height; k++) {
		delete[] image_ct[k];
		delete[] edge_ct[k];
		delete[] image_ct_filtered[k];
		delete[] maskThreshold[k];
		delete[] distanceMap[k];
	}
	delete[] image_ct;
	delete[] edge_ct;
	delete[] image_ct_filtered;
	delete[] maskThreshold;
	delete[] distanceMap;

	for (k = 0; k < height_pet; k++) {
		delete[] image_pet[k];
		delete[] edge_pet[k];
		delete[] edge_ct_interpolated[k];
		delete[] ct_interpolated[k];
		delete[] edge_ct_50_pet_50[k];
		delete[] edge_ct_75_pet_25[k];
		delete[] edge_ct_25_pet_75[k];
		//delete[] initial_seg[k];
	}
	delete[] image_pet;
	delete[] edge_pet;
	delete[] edge_ct_interpolated;
	delete[] ct_interpolated;
	delete[] edge_ct_50_pet_50;
	delete[] edge_ct_75_pet_25;
	delete[] edge_ct_25_pet_75;
	//delete[] initial_seg;

	delete[] center_seg;

	//for (k = 0; k < 255; k++) {
	//	delete[] ballPET[k];
	//}
	//delete[] ballPET;

	return EXIT_SUCCESS;
}
