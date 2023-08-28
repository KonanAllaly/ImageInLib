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

#define thres_min 995
#define thres_max 1213

int main() {

	int i, j, k;

	std::string inputPath = "C:/Users/Konan Allaly/Documents/Tests/input/";
	std::string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";

	//Point3D originCT = { -300.0, -230.0, -1022.5 }; // patient 2
	//Point3D originPET = { -286.586, -216.586, -1021.40 }; // patient 2
	//VoxelSpacing spacingCT = { _sx, _sy, _sz };
	//VoxelSpacing spacingPET = { _s_pet, _s_pet, _s_pet };

	//Point3D originCT = { -300.0, -230.0, -1339.2999 }; // patient 3
	//Point3D originPET = { -286.586, -216.586, -1338.800049 }; // patient 3
	//VoxelSpacing spacingCT = { _sx, _sy, _sz };
	//VoxelSpacing spacingPET = { _s_pet, _s_pet, _s_pet };

	Point3D originCT = { -249.5117, -403.0117, 675.0 }; // patient 6
	Point3D originPET = { -360.1, -515.421, 677.0 }; // patient 6
	VoxelSpacing spacingCT = { 0.97656, 0.97656, 1.5 };
	VoxelSpacing spacingPET = { 3.3, 3.3, 3.0 };
	
	OrientationMatrix orientation;
	orientation.v1 = { 1.0, 0.0, 0.0 }; orientation.v2 = { 0.0, 1.0, 0.0 }; orientation.v3 = { 0.0, 0.0, 1.0 };

	/*======== Data container for original ct data ============*/
	
	const size_t Height = 330, Width = 512, Length = 512;
	const size_t dim2D = Width * Length;
	
	short** image = new short * [Height];
	dataType** image_ct = new dataType * [Height];
	dataType** maskThreshold = new dataType * [Height];
	dataType** distanceMap = new dataType * [Height];
	dataType** ball = new dataType * [Height];
	int** segment = new int* [Height];
	bool** status = new bool* [Height];
	for (k = 0; k < Height; k++) {
		image[k] = new short[dim2D];
		image_ct[k] = new dataType[dim2D];
		maskThreshold[k] = new dataType[dim2D];
		distanceMap[k] = new dataType[dim2D];
		ball[k] = new dataType[dim2D];
		segment[k] = new int[dim2D];
		status[k] = new bool[dim2D];
	}
	if (image == NULL || image_ct == NULL || maskThreshold == NULL || distanceMap == NULL)
		return false;

	//Initialization
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			image[k][i] = 0.0;
			image_ct[k][i] = 0.0;
			maskThreshold[k][i] = 0.0;
			distanceMap[k][i] = 0.0;
			ball[k][i] = 0.0;
			segment[k][i] = 0;
			status[k][i] = false;
		}
	}

	std::string inputImagePath = inputPath + "patient4.raw";
	manageRAWFile3D<short>(image, Length, Width, Height, inputImagePath.c_str(), LOAD_DATA, false);
	//manageRAWFile3D<dataType>(image_ct, Length, Width, Height, inputImagePath.c_str(), LOAD_DATA, false);

	Image_Data IMAGE_CT;
	IMAGE_CT.imageDataPtr = image_ct;
	IMAGE_CT.height = Height; IMAGE_CT.length = Length; IMAGE_CT.width = Width;
	IMAGE_CT.origin = originCT;
	IMAGE_CT.spacing = spacingCT;
	IMAGE_CT.orientation = orientation;

	/*===================   Find the point with the highest distance ==========*/

	//fill the threshold container
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			image_ct[k][i] = (dataType)image[k][i];
			maskThreshold[k][i] = image_ct[k][i];
		}
	}

	std::string outputImagePath = outputPath + "LoadedP4.raw";
	manageRAWFile3D<dataType>(image_ct, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	//thresholding3dFunctionN(maskThreshold, Length, Width, Height, thres_min, thres_max, 0.0, 1.0);

	//Distance Map
	Distance_Map_Params params_dist = { 0.5, 1.0, 0.0, 1000000, 1e-3 };
	//computeDistanceMap(distanceMap, maskThreshold, Length, Width, Height, params_dist, FAST_SWEEP);

	outputImagePath = outputPath + "distanceP4.raw";
	//manageRAWFile3D<dataType>(distanceMap, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	//Find the point with the highest distance
	dataType dist_max = 0.0;
	int imax = 0, jmax = 0, kmax = 0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				if (distanceMap[k][x_new(i, j, Length)] > dist_max) {
					dist_max = distanceMap[k][x_new(i, j, Length)];
					kmax = k; imax = i; jmax = j;
				}
			}
		}
	}
	//std::cout << "Distance max = " << dist_max << std::endl;
	//std::cout << "Center CT (" << imax << ", " << jmax << ", " << kmax << ")" << std::endl;

	//ball in the Liver
	dataType radius = 5.0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				if (sqrt((imax - i) * (imax - i) + (jmax - j) * (jmax - j) + (kmax - k) * (kmax - k)) < radius) {
					ball[k][x_new(i, j, Length)] = 1;
				}
				else {
					ball[k][x_new(i, j, Length)] = 0;
				}
			}
		}
	}
	outputImagePath = outputPath + "ballP4.raw";
	//manageRAWFile3D<dataType>(ball, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	//Find threshold values around the point with the highest distance
	dataType local_max = 0.0, local_min = 1000000;
	size_t x = 0, voxel_used = 0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				x = x_new(i, j, Length);
				if (sqrt((imax - i) * (imax - i) + (jmax - j) * (jmax - j) + (kmax - k) * (kmax - k)) < radius) {
					voxel_used++;
					if (image_ct[k][x] > local_max) {
						local_max = image_ct[k][x];
					}
					if (image_ct[k][x] < local_min) {
						local_min = image_ct[k][x];
					}
				}
			}
		}
	}
	//std::cout << "Local max = " << local_max << ", " << " Local min = " << local_min << std::endl;
	//std::cout << voxel_used << " voxels were used to find local threshold" << std::endl;

	//fill the threshold container
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			maskThreshold[k][i] = image_ct[k][i];
		}
	}
	//thresholding3dFunctionN(maskThreshold, Length, Width, Height, local_min, local_max, 0.0, 1.0);

	/*===================== Labelling CT ========================================*/

	//erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	//erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	//erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	//erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	//erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	//erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	//labelling3D(maskThreshold, segment, status, Length, Width, Height, 1.0);

	//number of region voxels
	int numberOfRegionsCells = 0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (maskThreshold[k][i] == 1.0) {
				numberOfRegionsCells++;
			}
		}
	}
	//std::cout << "number of cell : " << numberOfRegionsCells << std::endl;

	//Counting
	int* countingArray = (int*)malloc(numberOfRegionsCells * sizeof(int));
	if (countingArray == NULL) return false;
	for (i = 0; i < numberOfRegionsCells; i++)
		countingArray[i] = 0;

	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (segment[k][i] > 0) {
				countingArray[segment[k][i]]++;
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
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (countingArray[segment[k][i]] != maxElement) {
				segment[k][i] = 0;
			}
		}
	}

	int liver_number_of_voxel = 0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (segment[k][i] != 0) {
				maskThreshold[k][i] = 1.0;
				liver_number_of_voxel++;
			}
			else {
				maskThreshold[k][i] = 0.0;
			}
		}
	}
	//std::cout << "The segmented Liver contains : " << liver_number_of_voxel << " voxels" << std::endl;
	//outputImagePath = outputPath + "biggest_regionP7.raw";
	//manageRAWFile3D<dataType>(maskThreshold, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	free(countingArray);

	/*===================== Data container for orginal PET data =============*/

	//const size_t height_pet = 255, length_pet = 144, width_pet = 144; // patient2

	const size_t height_pet = 318, length_pet = 144, width_pet = 144; // patient3
	 
	//const size_t height_pet = 192, length_pet = 81, width_pet = 120; // patient1
	//const size_t height_pet = 195, length_pet = 128, width_pet = 128; // patient4, patient5
	//const size_t height_pet = 440, length_pet = 220, width_pet = 220; // patient6
	//const size_t height_pet = 440, length_pet = 220, width_pet = 220; // patient7
	
	dataType** image_pet = new dataType * [height_pet];
	short** pet = new short * [height_pet];
	dataType** mask_pet = new dataType * [height_pet];
	dataType** ball_pet = new dataType * [height_pet];
	//int** segment = new int * [height_pet];
	//bool** status = new bool* [height_pet];
	for (k = 0; k < height_pet; k++) {
		image_pet[k] = new dataType[length_pet * width_pet];
		pet[k] = new short[length_pet * width_pet];
		mask_pet[k] = new dataType[length_pet * width_pet];
		ball_pet[k] = new dataType[length_pet * width_pet];
		//segment[k] = new int[length_pet * width_pet];
		//status[k] = new bool[length_pet * width_pet];
	}
	if (image_pet == NULL || pet == NULL || mask_pet == NULL)
		return false;

	for (k = 0; k < height_pet; k++) {
		for (i = 0; i < length_pet * width_pet; i++) {
			image_pet[k][i] = 0.0;
			pet[k][i] = 0;
			mask_pet[k][i] = 0.0;
			ball_pet[k][i] = 0.0;
			//segment[k][i] = 0;
			//status[k][i] = false;
		}
	}

	Image_Data IMAGE_PET;
	IMAGE_PET.imageDataPtr = image_pet;
	IMAGE_PET.height = height_pet; IMAGE_PET.length = length_pet; IMAGE_PET.width = width_pet;
	IMAGE_PET.origin = originPET;
	IMAGE_PET.spacing = spacingPET;
	IMAGE_PET.orientation = orientation;

	//inputImagePath = inputPath + "patient3_pet.raw";
	//manageRAWFile3D<dataType>(image_pet, length_pet, width_pet, height_pet, inputImagePath.c_str(), LOAD_DATA, true);

	//outputImagePath = outputPath + "loadedPET.raw";
	//manageRAWFile3D<dataType>(image_pet, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	////copy to mask
	//for (k = 0; k < height_pet; k++) {
	//	for (i = 0; i < length_pet * width_pet; i++) {
	//		mask_pet[k][i] = image_pet[k][i];
	//	}
	//}

	//thresholding3dFunctionN(mask_pet, length_pet, width_pet, height_pet, 5700, 9798, 0.0, 1.0); // patient 2
	//thresholding3dFunctionN(mask_pet, length_pet, width_pet, height_pet, 9000, 10000, 0.0, 1.0); // patient 3
	//thresholding3dFunctionN(mask_pet, length_pet, width_pet, height_pet, 8000, 15000, 0.0, 1.0); //patient 4
	//thresholding3dFunctionN(mask_pet, length_pet, width_pet, height_pet, 4250, 8000, 0.0, 1.0); //patient 5,6
	//thresholding3dFunctionN(mask_pet, length_pet, width_pet, height_pet, 4250, 8000, 0.0, 1.0); //patient 7

	/*==================  Interpolated Coordinates ========*/

	//Point3D centerCT = { imax, jmax, kmax }, centerPET = { 0.0, 0.0, 0.0 };
	//centerPET = getRealCoordFromImageCoord3D(centerCT, originCT, spacingCT, orientation);
	//centerPET = getImageCoordFromRealCoord3D(centerPET, originPET, spacingPET, orientation);

	//int ipet = (int)centerPET.x, jpet = (int)centerPET.y, kpet = (int)centerPET.z;
	////std::cout << "Center PET (" << ipet << ", " << jpet << ", " << kpet << ")" << std::endl;

	/*================= Local threshold for PET ==============*/

	////Find threshold values around the point with the highest distance
	//dataType local_max = 0.0, local_min = 100000000, radius = 7.0;
	//size_t x = 0, nb_local_voxels = 0;
	//for (k = 0; k < height_pet; k++) {
	//	for (i = 0; i < length_pet; i++) {
	//		for (j = 0; j < width_pet; j++) {
	//			x = x_new(i, j, length_pet);
	//			if (sqrt((ipet - i) * (ipet - i) + (jpet - j) * (jpet - j) + (kpet - k) * (kpet - k)) < radius) {
	//				ball_pet[k][x] = 1.0;
	//				nb_local_voxels++;
	//				if (image_pet[k][x] > local_max) {
	//					local_max = image_pet[k][x];
	//				}
	//				if (image_pet[k][x] < local_min) {
	//					local_min = image_pet[k][x];
	//				}
	//			}
	//		}
	//	}
	//}
	//std::cout << "Local max = " << local_max << ", " << " Local min = " << local_min << std::endl;
	//std::cout << nb_local_voxels << " voxels were used to find those local threshold values" << std::endl;

	//outputImagePath = outputPath + "ball_pet.raw";
	//manageRAWFile3D<dataType>(ball_pet, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	//thresholding3dFunctionN(mask_pet, length_pet, width_pet, height_pet, local_min, local_max, 0.0, 1.0);

	//======================  Labelling PET ========================

	//erosion3dHeighteenNeigbours(mask_pet, length_pet, width_pet, height_pet, 1.0, 0.0);
	////erosion3dHeighteenNeigbours(mask_pet, length_pet, width_pet, height_pet, 1.0, 0.0);
	////erosion3dHeighteenNeigbours(mask_pet, length_pet, width_pet, height_pet, 1.0, 0.0);
	////erosion3dHeighteenNeigbours(mask_pet, length_pet, width_pet, height_pet, 1.0, 0.0);
	////erosion3dHeighteenNeigbours(mask_pet, length_pet, width_pet, height_pet, 1.0, 0.0);
	// 
	//labelling3D(mask_pet, segment, status, length_pet, width_pet, height_pet, 1.0);

	////number of region voxels
	//int numberOfRegionsCells = 0;
	//for (k = 0; k < height_pet; k++) {
	//	for (i = 0; i < length_pet * width_pet; i++) {
	//		if (mask_pet[k][i] == 1.0) {
	//			numberOfRegionsCells++;
	//		}
	//	}
	//}
	//std::cout << "number of cell : " << numberOfRegionsCells << std::endl;
	// 
	////Counting
	//int* countingArray = (int*)malloc(numberOfRegionsCells * sizeof(int));
	//if (countingArray == NULL) return false;
	//for (i = 0; i < numberOfRegionsCells; i++) 
	//	countingArray[i] = 0;

	//for (k = 0; k < height_pet; k++) {
	//	for (i = 0; i < length_pet * width_pet; i++) {
	//		if (segment[k][i] > 0) {
	//			countingArray[segment[k][i]]++;
	//		}
	//	}
	//}
	////Find the biggest region
	//int maxElement = 0;
	//for (k = 0; k < numberOfRegionsCells; k++) {
	//	if (countingArray[k] > maxElement) {
	//		maxElement = countingArray[k];
	//	}
	//}
	////Keep just the biggest region 
	//for (k = 0; k < height_pet; k++) {
	//	for (i = 0; i < length_pet * width_pet; i++) {
	//		if (countingArray[segment[k][i]] != maxElement) {
	//			segment[k][i] = 0;
	//		}
	//	}
	//}

	//int liver_volume = 0;
	//for (k = 0; k < height_pet; k++) {
	//	for (i = 0; i < length_pet * width_pet; i++) {
	//		if (segment[k][i] != 0) {
	//			mask_pet[k][i] = 1.0;
	//			liver_volume++;
	//		}
	//		else {
	//			mask_pet[k][i] = 0.0;
	//		}
	//	}
	//}
	//std::cout << "The segmented Liver volume is : " << liver_volume << " voxels" << std::endl;

	//outputImagePath = outputPath + "biggest_region_pet.raw";
	//manageRAWFile3D<dataType>(mask_pet, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	//free(countingArray);

	/*============= Segmentation ================*/

	//Point3D* center_seg = new Point3D[1];
	//center_seg->x = idown; center_seg->y = jdown; center_seg->z = kdown;
	//center_seg->x = ipet; center_seg->y = jpet; center_seg->z = kpet;
	//center_seg->x = imax; center_seg->y = jmax; center_seg->z = kmax;

	//dataType** initial_seg = new dataType * [height_pet];
	//for (k = 0; k < height_pet; k++) {
	//	initial_seg[k] = new dataType[length_pet * width_pet];
	//}
	//for (k = 0; k < height_pet; k++) {
	//	for (i = 0; i < width_pet * length_pet; i++) {
	//		initial_seg[k][i] = 0.0;
	//	}
	//}

	//generateInitialSegmentationFunctionForMultipleCentres(initial_seg, length_pet, width_pet, height_pet, center_seg, 1.0, 10.0, 1.0);

	//Image_Data ImageToBeSegmented;
	//ImageToBeSegmented.imageDataPtr = image_pet;
	//ImageToBeSegmented.height = height_pet; ImageToBeSegmented.length = length_pet; ImageToBeSegmented.width = width_pet;

	//Segmentation_Parameters seg_params;
	//seg_params.tau = 1.0; seg_params.h = 0.1; seg_params.omega_c = 1.4; seg_params.mod = 10;
	//seg_params.maxNoGSIteration = 100; seg_params.coef = 100; seg_params.gauss_seidelTolerance = 1e-6;
	//seg_params.maxNoOfTimeSteps = 4000; seg_params.segTolerance = 1e-6; seg_params.eps2 = 1e-6;
	//seg_params.coef_conv = 10.0; seg_params.coef_dif = 1.0;

	//Smoothing by heat explicit so we don't need all the parameters
	//Filter_Parameters parameters; parameters.coef = 1e-4;
	//parameters.h = 6.0; parameters.timeStepSize = 1.0; parameters.timeStepsNum = 10;
	//Filter_Parameters parameters; parameters.coef = 1e-4; parameters.edge_detector_coefficient = 1000; parameters.eps2 = 1e-4;
	//parameters.h = 6.0; parameters.maxNumberOfSolverIteration = 100; parameters.omega_c = 1.4; parameters.p = 1;
	//parameters.sigma = 0.1; parameters.timeStepSize = 1.0; parameters.timeStepsNum = 1; parameters.tolerance = 1e-6;

	//unsigned char segmentPath[] = "C:/Users/Konan Allaly/Documents/Tests/output/segmentation/";

	////subsurfSegmentation(F_Image, initial_seg, seg_params, parameters, center_seg, 1, segmentPath);
	//generalizedSubsurfSegmentation(ImageToBeSegmented, initial_seg, seg_params, parameters, center_seg, 1, segmentPath);

	for (k = 0; k < Height; k++) {
		delete[] image_ct[k];
		delete[] ball[k];
		delete[] maskThreshold[k];
		delete[] distanceMap[k];
		delete[] segment[k];
		delete[] status[k];
	}
	delete[] image_ct;
	delete[] ball;
	delete[] maskThreshold;
	delete[] distanceMap;
	delete[] segment;
	delete[] status;

	for (k = 0; k < height_pet; k++) {
		delete[] image_pet[k];
		delete[] pet[k];
		delete[] mask_pet[k];
		delete[] ball_pet[k];
		//delete[] segment[k];
		//delete[] status[k];
	}
	delete[] image_pet;
	delete[] pet;
	delete[] mask_pet;
	delete[] ball_pet;
	//delete[] segment;
	//delete[] status;

	return EXIT_SUCCESS;
}
