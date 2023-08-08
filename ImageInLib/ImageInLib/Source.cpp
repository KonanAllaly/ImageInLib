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

#define thres_min 0.225
#define thres_max 0.325

//#define thres_min 900 //984
//#define thres_max 1300 //1274


int main() {

	int i, j, k;

	std::string inputPath = "C:/Users/Konan Allaly/Documents/Tests/input/";
	std::string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";

	Point3D originCT = { -300.0, -230.0, -1339.29999 };
	Point3D originPET = { -286.586, -216.586, -1338.80 };
	VoxelSpacing spacingCT = { _sx, _sy, _sz };
	VoxelSpacing spacingPET = { _s_pet, _s_pet, _s_pet };
	OrientationMatrix orientation;
	orientation.v1 = { 1.0, 0.0, 0.0 }; orientation.v2 = { 0.0, 1.0, 0.0 }; orientation.v3 = { 0.0, 0.0, 1.0 };

	/*======== Data container for original ct data ============*/
	
	//const size_t Height = 406, Width = 512, Length = 512;
	//const size_t dim2D = Width * Length;
	//dataType** image_ct = new dataType * [Height];
	//dataType** edge_ct = new dataType * [Height];
	//dataType** image_ct_filtered = new dataType * [Height];
	//dataType** maskThreshold = new dataType * [Height];
	//dataType** distanceMap = new dataType * [Height];
	//for (k = 0; k < Height; k++) {
	//	image_ct[k] = new dataType[dim2D];
	//	edge_ct[k] = new dataType[dim2D];
	//	image_ct_filtered[k] = new dataType[dim2D];
	//	maskThreshold[k] = new dataType[dim2D];
	//	distanceMap[k] = new dataType[dim2D];
	//}
	//if (image_ct == NULL || edge_ct == NULL || image_ct_filtered == NULL || maskThreshold == NULL || distanceMap == NULL)
	//	return false;

	////Initialization
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		image_ct[k][i] = 0.0;
	//		edge_ct[k][i] = 0.0;
	//		image_ct_filtered[k][i] = 0.0;
	//		maskThreshold[k][i] = 0.0;
	//		distanceMap[k][i] = 0.0;
	//	}
	//}

	//std::string inputImagePath = inputPath + "patient2_ct_rescaled.raw";
	//manageRAWFile3D<dataType>(image_ct, Length, Width, Height, inputImagePath.c_str(), LOAD_DATA, false);

	//Image_Data IMAGE_CT;
	//IMAGE_CT.imageDataPtr = image_ct;
	//IMAGE_CT.height = Height; IMAGE_CT.length = Length; IMAGE_CT.width = Width;
	//IMAGE_CT.origin = originCT;
	//IMAGE_CT.spacing = spacingCT;
	//IMAGE_CT.orientation = orientation;

	//inputImagePath = inputPath + "filtered_patient2_ct.raw";
	//manageRAWFile3D<dataType>(image_ct_filtered, Length, Width, Height, inputImagePath.c_str(), LOAD_DATA, false);

	/*===================   Find the point with the higest distance ==========*/

	////fill the threshold container
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		maskThreshold[k][i] = image_ct[k][i];
	//	}
	//}

	//thresholding3dFunctionN(maskThreshold, Length, Width, Height, thres_min, thres_max, 0.0, 1.0);

	////std::string outputImagePath = outputPath + "threshold.raw";
	////manageRAWFile3D<dataType>(maskThreshold, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	//Distance Map
	//Distance_Map_Params params_dist = { 0.5, 1.0, 0.0, 1000000, 1e-3 };
	//computeDistanceMap(distanceMap, maskThreshold, Length, Width, Height, params_dist, FAST_SWEEP);

	//std::string outputImagePath = outputPath + "distance.raw";
	//manageRAWFile3D<dataType>(distanceMap, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

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

	//dataType ref_value = image_ct_filtered[kmax][x_new(imax, jmax, Length)];

	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			maskThreshold[k][x_new(i, j, Length)] = exp(pow(image_ct_filtered[k][x_new(i, j, Length)] - ref_value,2) / (2 * pow(20,2)));
	//		}
	//	}
	//}

	//std::string outputImagePath = outputPath + "enhanced.raw";
	//manageRAWFile3D<dataType>(maskThreshold, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	/*===================== Data container for orginal PET data =============*/

	const size_t height_pet = 255, length_pet = 144, width_pet = 144; // patient2
	//const size_t height_pet = 318, length_pet = 144, width_pet = 144; // patient3
	//const size_t height_pet = 192, length_pet = 81, width_pet = 120; // patient1
	//const size_t height_pet = 195, length_pet = 128, width_pet = 128; // patient4, patient5
	//const size_t height_pet = 440, length_pet = 220, width_pet = 220; // patient6
	//const size_t height_pet = 440, length_pet = 220, width_pet = 220; // patient7
	dataType** image_pet = new dataType * [height_pet];
	short** pet = new short * [height_pet];
	dataType** mask_pet = new dataType * [height_pet];
	//dataType** distance_pet = new dataType * [height_pet];
	int** segment = new int * [height_pet];
	bool** status = new bool* [height_pet];
	for (k = 0; k < height_pet; k++) {
		image_pet[k] = new dataType[length_pet * width_pet];
		pet[k] = new short[length_pet * width_pet];
		mask_pet[k] = new dataType[length_pet * width_pet];
		//distance_pet[k] = new dataType[length_pet * width_pet];
		segment[k] = new int[length_pet * width_pet];
		status[k] = new bool[length_pet * width_pet];
	}
	if (image_pet == NULL || pet == NULL || mask_pet == NULL)
		return false;

	for (k = 0; k < height_pet; k++) {
		for (i = 0; i < length_pet * width_pet; i++) {
			image_pet[k][i] = 0.0;
			pet[k][i] = 0.0;
			mask_pet[k][i] = 0.0;
			//distance_pet[k][i] = 0.0;
			segment[k][i] = 0;
			status[k][i] = false;
		}
	}

	Image_Data IMAGE_PET;
	IMAGE_PET.imageDataPtr = image_pet;
	IMAGE_PET.height = height_pet; IMAGE_PET.length = length_pet; IMAGE_PET.width = width_pet;
	IMAGE_PET.origin = originPET;
	IMAGE_PET.spacing = spacingPET;
	IMAGE_PET.orientation = orientation;

	//std::string inputImagePath = inputPath + "patient3_pet.raw";
	std::string inputImagePath = inputPath + "patient2_pet.raw";
	//manageRAWFile3D<short>(pet, length_pet, width_pet, height_pet, inputImagePath.c_str(), LOAD_DATA, true);
	manageRAWFile3D<dataType>(image_pet, length_pet, width_pet, height_pet, inputImagePath.c_str(), LOAD_DATA, false);

	//find min, max data
	dataType minData = 10000000, maxData = 0;
	for (k = 0; k < height_pet; k++) {
		for (i = 0; i < length_pet * width_pet; i++) {
			
			if (image_pet[k][i] > maxData) {
				maxData = image_pet[k][i];
			}

			if (image_pet[k][i] < minData) {
				minData = image_pet[k][i];
			}
		}
	}
	std::cout << "min = " << minData << ", max = " << maxData << std::endl;

	//rescaleNewRange(image_pet, length_pet, width_pet, height_pet, 0, 255, maxData, minData);

	std::string outputImagePath = outputPath + "rescale.raw";
	//manageRAWFile3D<dataType>(image_pet, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);
	
	//int length_histo = (int)(255 - 0 + 1);
	int length_histo = (int)(maxData - minData + 1);

	int* histogram = new int[length_histo];
	for (j = 0; j < length_histo; j++) histogram[j] = 0;
	int value = 0, total = 0;
	for (k = 0; k < height_pet; k++) {
		for (i = 0; i < length_pet * width_pet; i++) {
			value = (int)image_pet[k][i];
			histogram[value]++;
		}
	}

	outputImagePath = outputPath + "intensity.csv";
	FILE* intensity;
	if (fopen_s(&intensity, outputImagePath.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}

	outputImagePath = outputPath + "count.csv";
	FILE* count;
	if (fopen_s(&count, outputImagePath.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}

	outputImagePath = outputPath + "histogram.csv";
	FILE* histo;
	if (fopen_s(&histo, outputImagePath.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	
	int cpt;
	//fprintf(intensity, "intensity \n");
	//fprintf(count, "count \n");

	for (k = 0; k < length_histo; k++) {
		//if (histogram[k] > 100) {
		//	fprintf(intensity, "%d, \n", k);
		//	fprintf(count, "%d, \n", histogram[k]);
		//}
		//fprintf(intensity, "%d, \n", k);
		//fprintf(count, "%d, \n", histogram[k]);
		//total = total + histogram[k];
		fprintf(intensity, "%d, \n", k);
		fprintf(count, "%d, \n", histogram[k]);
		fprintf(histo, "%d, %d \n", k, histogram[k]);
	}
	fclose(intensity);
	fclose(count);
	fclose(histo);

	//std::cout << "number of elements : " << cpt << std::endl;
	//std::cout << "voxel number : " << total << std::endl;
	//std::cout << "intensity value 0 has : " << histogram[0] << " elements" << std::endl;

	delete[] histogram;

	//std::string outputImagePath = outputPath + "ct_down_sampled_144_255.raw";
	//manageRAWFile3D<dataType>(image_pet, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	//Copy
	for (k = 0; k < height_pet; k++) {
		for (i = 0; i < length_pet * width_pet; i++) {
			//image_pet[k][i] = (dataType)pet[k][i];
			mask_pet[k][i] = image_pet[k][i];
		}
	}

	//std::string outputImagePath = outputPath + "loaded_pet_p1.raw";
	//manageRAWFile3D<dataType>(mask_pet, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	////Rescalling
	//dataType maxData = 0.0, minData = 10000;
	//for (k = 0; k < height_pet; k++) {
	//	for (i = 0; i < length_pet; i++) {
	//		for (j = 0; j < width_pet; j++) {
	//			if (image_pet[k][x_new(i, j, length_pet)] > maxData) {
	//				maxData = image_pet[k][x_new(i, j, length_pet)];
	//			}
	//			if (image_pet[k][x_new(i, j, length_pet)] < minData) {
	//				minData = image_pet[k][x_new(i, j, length_pet)];
	//			}
	//		}
	//	}
	//}
	//std::cout << "max = " << maxData << " , min = " << minData << std::endl;
	//rescaleNewRange(image_pet, length_pet, width_pet, height_pet, 0.0, 1.0, maxData, minData);

	//thresholding3dFunctionN(mask_pet, length_pet, width_pet, height_pet, 5700, 9798, 0.0, 1.0); // patient 2
	//thresholding3dFunctionN(mask_pet, length_pet, width_pet, height_pet, 9000, 10000, 0.0, 1.0); // patient 3
	//thresholding3dFunctionN(mask_pet, length_pet, width_pet, height_pet, 8000, 15000, 0.0, 1.0); //patient 4
	//thresholding3dFunctionN(mask_pet, length_pet, width_pet, height_pet, 4250, 8000, 0.0, 1.0); //patient 5,6
	//thresholding3dFunctionN(mask_pet, length_pet, width_pet, height_pet, 4250, 8000, 0.0, 1.0); //patient 7

	//outputImagePath = outputPath + "threshold.raw";
	//manageRAWFile3D<dataType>(mask_pet, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	//Distance Map
	//Distance_Map_Params params_dist = { 0.5, 1.0, 0.0, 1000000, 1e-3 };
	//computeDistanceMap(distance_pet, mask_pet, length_pet, width_pet, height_pet, params_dist, FAST_SWEEP);
	//std::string outputImagePath = outputPath + "distance_pet.raw";
	//manageRAWFile3D<dataType>(distance_pet, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	////find the point with the max distance 
	//dataType dist_max = 0.0;
	//int ipet = 0, jpet = 0, kpet = 0;
	//for (k = 0; k < height_pet; k++) {
	//	for (i = 0; i < length_pet; i++) {
	//		for (j = 0; j < width_pet; j++) {
	//			if (distance_pet[k][x_new(i, j, length_pet)] > dist_max) {
	//				dist_max = distance_pet[k][x_new(i, j, length_pet)];
	//				ipet = i; jpet = j; kpet = k;
	//			}
	//		}
	//	}
	//}
	//std::cout << "Coordinates of the starting point : (" << ipet << " ," << jpet << " ," << kpet << ")" << std::endl;

	/*==================  Interpolated Coordinates ========*/

	//Point3D centerCT = { imax, jmax, kmax }, centerPET = { 0.0, 0.0, 0.0 };

	//centerPET = getRealCoordFromImageCoord3D(centerCT, originCT, spacingCT, orientation);
	//centerPET = getImageCoordFromRealCoord3D(centerPET, originPET, spacingPET, orientation);

	//int ipet = (int)centerPET.x, jpet = (int)centerPET.y, kpet = (int)centerPET.z;
	//std::cout << "Center pet (" << ipet << ", " << jpet << ", " << kpet << ")" << std::endl;

	/*================= Interpolation input Image and Edge detector ====================*/

	//Load edge detector ct
	//inputImagePath = inputPath + "edge_detector_ct.raw";
	//manageRAWFile3D<dataType>(edge_ct, Length, Width, Height, inputImagePath.c_str(), LOAD_DATA, false);

	//Load edge detector pet
	//inputImagePath = inputPath + "edge_detector_pet.raw";
	//manageRAWFile3D<dataType>(edge_pet, length_pet, width_pet, height_pet, inputImagePath.c_str(), LOAD_DATA, false);

	//Image_Data EDGE_CT;
	//EDGE_CT.imageDataPtr = edge_ct;
	//EDGE_CT.height = Height; EDGE_CT.length = Length; EDGE_CT.width = Width;
	//EDGE_CT.origin = originCT;
	//EDGE_CT.spacing = spacingCT;
	//EDGE_CT.orientation = orientation;

	//Image_Data EDGE_PET;
	//EDGE_PET.imageDataPtr = edge_pet;
	//EDGE_PET.height = height_pet; EDGE_PET.length = length_pet; EDGE_PET.width = width_pet;
	//EDGE_PET.origin = originPET;
	//EDGE_PET.spacing = spacingPET;
	//EDGE_PET.orientation = orientation;

	//Interpolate From CT dim to PET dim
	//dataType** ct_interpolated = new dataType * [height_pet];
	//dataType** edge_ct_interpolated = new dataType * [height_pet];
	//for (k = 0; k < height_pet; k++) {
	//	ct_interpolated[k] = new dataType[length_pet * width_pet];
	//	edge_ct_interpolated[k] = new dataType[length_pet * width_pet];
	//}
	//if (ct_interpolated == NULL || edge_ct_interpolated == NULL)
	//	return false;

	////Initialization
	//for (k = 0; k < height_pet; k++) {
	//	for (i = 0; i < length_pet * width_pet; i++) {
	//		ct_interpolated[k][i] = 0.0;
	//		edge_ct_interpolated[k][i] = 0.0;
	//	}
	//}

	//Image_Data CT_interpolate;
	//CT_interpolate.imageDataPtr = ct_interpolated;
	//CT_interpolate.height = height_pet; CT_interpolate.length = length_pet; CT_interpolate.width = width_pet;
	//CT_interpolate.origin = originPET; CT_interpolate.spacing = spacingPET; CT_interpolate.orientation = orientation;

	//Image_Data edge_CT_interpolate;
	//edge_CT_interpolate.imageDataPtr = edge_ct_interpolated;
	//edge_CT_interpolate.height = height_pet; edge_CT_interpolate.length = length_pet; edge_CT_interpolate.width = width_pet;
	//edge_CT_interpolate.origin = originPET; edge_CT_interpolate.spacing = spacingPET; edge_CT_interpolate.orientation = orientation;

	//imageInterpolation3D(IMAGE_CT, CT_interpolate, TRILINEAR);
	//imageInterpolation3D(EDGE_CT, edge_CT_interpolate, TRILINEAR);

	//outputImagePath = outputPath + "edge_ct_in_pet_dim.raw";
	//manageRAWFile3D<dataType>(edge_ct_interpolated, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	//std::string  outputImagePath = outputPath + "interpolated_ct.raw";
	//manageRAWFile3D<dataType>(ct_interpolated, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	//outputImagePath = outputPath + "_edge_ct_interpolated.raw";
	//manageRAWFile3D<dataType>(edge_ct_interpolated, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	/*=============== Compute new edge detector ===========*/

	//dataType** edge_ct_50_pet_50 = new dataType * [height_pet];
	//dataType** edge_ct_75_pet_25 = new dataType * [height_pet];
	//dataType** edge_ct_25_pet_75 = new dataType * [height_pet];
	//for (k = 0; k < height_pet; k++) {
	//	edge_ct_50_pet_50[k] = new dataType[length_pet * width_pet];
	//	edge_ct_75_pet_25[k] = new dataType[length_pet * width_pet];
	//	edge_ct_25_pet_75[k] = new dataType[length_pet * width_pet];
	//}
	//if (edge_ct_50_pet_50 == NULL || edge_ct_75_pet_25 == NULL || edge_ct_25_pet_75 == NULL)
	//	return false;

	////Initialization
	//for (k = 0; k < height_pet; k++) {
	//	for (i = 0; i < length_pet * width_pet; i++) {
	//		edge_ct_50_pet_50[k][i] = 0.0;
	//		edge_ct_75_pet_25[k][i] = 0.0;
	//		edge_ct_25_pet_75[k][i] = 0.0;
	//	}
	//}

	//for (k = 0; k < height_pet; k++) {
	//	for (i = 0; i < length_pet * width_pet; i++) {
	//		edge_ct_50_pet_50[k][i] = 0.5 * edge_ct_interpolated[k][i] + 0.5 * edge_pet[k][i];
	//		edge_ct_75_pet_25[k][i] = 0.75 * edge_ct_interpolated[k][i] + 0.25 * edge_pet[k][i];
	//		edge_ct_25_pet_75[k][i] = 0.25 * edge_ct_interpolated[k][i] + 0.75 * edge_pet[k][i];
	//	}
	//}

	//outputImagePath = outputPath + "edge_ct_50_pet_50.raw";
	//manageRAWFile3D<dataType>(edge_ct_50_pet_50, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	//outputImagePath = outputPath + "edge_ct_75_pet_25.raw";
	//manageRAWFile3D<dataType>(edge_ct_75_pet_25, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	//outputImagePath = outputPath + "edge_ct_25_pet_75.raw";
	//manageRAWFile3D<dataType>(edge_ct_25_pet_75, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	/*==================== Downsampling ================*/

	//const size_t height = 127, length = 72, width = 72;
	////const size_t height = 170, length = 96, width = 96;
	//dataType** image_down = new dataType * [height];
	//dataType** mask_down = new dataType * [height];
	//dataType** distance_down = new dataType * [height];
	//for (k = 0; k < height; k++) {
	//	image_down[k] = new dataType[length * width];
	//	mask_down[k] = new dataType[length * width];
	//	distance_down[k] = new dataType[length * width];
	//}
	//if (image_down == NULL || maskThreshold == NULL || distance_down == NULL)
	//	return false;

	////Initialization
	//for (k = 0; k < height; k++) {
	//	for (i = 0; i < length * width; i++) {
	//		image_down[k][i] = 0.0;
	//		mask_down[k][i] = 0.0;
	//		distance_down[k][i] = 0.0;
	//	}
	//}

	//down sample edge ct
	//edge_CT_interpolate.origin = originPET;
	//edge_CT_interpolate.spacing = { 1.0, 1.0, 1.0 };

	//Image_Data ImageD;
	//ImageD.imageDataPtr = image_down;
	//ImageD.height = height; ImageD.length = length; ImageD.width = width;
	//ImageD.origin = originPET;
	//ImageD.spacing = { 8.0, 8.0, 8.0 };
	//ImageD.orientation = orientation;

	//imageInterpolation3D(IMAGE_PET, ImageD, TRILINEAR);
	//std::string outputImagePath = outputPath + "ct_down_sampled_72_127.raw";
	//manageRAWFile3D<dataType>(image_down, length, width, height, outputImagePath.c_str(), STORE_DATA, false);

	////Copy
	//for (k = 0; k < height; k++) {
	//	for (i = 0; i < length * width; i++) {
	//		mask_down[k][i] = image_down[k][i];
	//	}
	//}

	//thresholding3dFunctionN(mask_down, length, width, height, thres_min, thres_max, 0.0, 1.0);

	//outputImagePath = outputPath + "threshold_down.raw";
	//manageRAWFile3D<dataType>(image_down, length, width, height, outputImagePath.c_str(), STORE_DATA, false);

	//Distance Map
	//Distance_Map_Params params_dist = { 0.5, 1.0, 0.0, 1000000, 1e-3 };
	//computeDistanceMap(distance_down, mask_down, length, width, height, params_dist, FAST_SWEEP);
	//std::string outputImagePath = outputPath + "distance_down2.raw";
	//manageRAWFile3D<dataType>(distance_down, length, width, height, outputImagePath.c_str(), STORE_DATA, false);

	////find the point with the max distance 
	//dataType dist_max = 0.0;
	//int idown = 0, jdown = 0, kdown = 0;
	//for (k = 0; k < height; k++) {
	//	for (i = 0; i < length; i++) {
	//		for (j = 0; j < width; j++) {
	//			if (distance_down[k][x_new(i, j, length)] > dist_max) {
	//				dist_max = distance_down[k][x_new(i, j, length)];
	//				idown = i; jdown = j; kdown = k;
	//			}
	//		}
	//	}
	//}
	//std::cout << "Coordinates of the starting point : (" << idown << " ," << jdown << " ," << kdown << ")" << std::endl;

	////use the filtered image
	//IMAGE_CT.imageDataPtr = image_ct_filtered;
	//imageInterpolation3D(IMAGE_CT, IMAGE_PET, TRILINEAR);
	//imageInterpolation3D(IMAGE_PET, ImageD, TRILINEAR);

	////Rescalling
	//dataType maxData = 0.0, minData = 10000;
	//for (k = 0; k < height; k++) {
	//	for (i = 0; i < length; i++) {
	//		for (j = 0; j < width; j++) {
	//			if (image_down[k][x_new(i, j, length)] > maxData) {
	//				maxData = image_down[k][x_new(i, j, length)];
	//			}
	//			if (image_down[k][x_new(i, j, length)] < minData) {
	//				minData = image_down[k][x_new(i, j, length)];
	//			}
	//		}
	//	}
	//}
	//std::cout << "max = " << maxData << " , min = " << minData << std::endl;
	//rescaleNewRange(image_down, length, width, height, 0.0, 1.0, maxData, minData);

	//Point3D centerDown = { 0.0, 0.0, 0.0 };

	//centerDown = getRealCoordFromImageCoord3D(centerPET, originPET, edge_CT_interpolate.spacing, orientation);
	//centerDown = getImageCoordFromRealCoord3D(centerDown, ImageD.origin, ImageD.spacing, orientation);

	//int idown = (int)centerDown.x, jdown = (int)centerDown.y, kdown = (int)centerDown.z;
	//std::cout << "Center down sampled data (" << idown << ", " << jdown << ", " << kdown << ")" << std::endl;

	//======================  Labelling ========================

	//dilatation3dHeighteenNeigbours(mask_pet, length_pet, width_pet, height_pet, 1.0, 0.0);
	
	//thresholding3dFunctionN(mask_down, length, width, height, thres_min, thres_max, 0.0, 1.0);
	//dilatation3dHeighteenNeigbours(image, Length_new, Width_new, Height_new, 1.0, 0.0);

	//erosion3dHeighteenNeigbours(mask_pet, length_pet, width_pet, height_pet, 1.0, 0.0);

	//dilatation3dHeighteenNeigbours(mask_pet, length_pet, width_pet, height_pet, 1.0, 0.0);
	//erosion3dHeighteenNeigbours(mask_pet, length_pet, width_pet, height_pet, 1.0, 0.0);

	//dilatation3dHeighteenNeigbours(image, Length_new, Width_new, Height_new, 1.0, 0.0);
	//erosion3dHeighteenNeigbours(image, Length_new, Width_new, Height_new, 1.0, 0.0);
	//labelling3D(image, segmentedImage, statusArray, Length_new, Width_new, Height_new, 1.0);

	//outputImagePath = outputPath + "eroded.raw";
	//manageRAWFile3D<dataType>(mask_pet, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	//labelling3D(mask_pet, segment, status, length_pet, width_pet, height_pet, 1.0);

	//number of region voxels
	int numberOfRegionsCells = 0;
	for (k = 0; k < height_pet; k++) {
		for (i = 0; i < length_pet * width_pet; i++) {
			if (mask_pet[k][i] == 1.0) {
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

	for (k = 0; k < height_pet; k++) {
		for (i = 0; i < length_pet * width_pet; i++) {
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
	for (k = 0; k < height_pet; k++) {
		for (i = 0; i < length_pet * width_pet; i++) {
			if (countingArray[segment[k][i]] != maxElement) {
				segment[k][i] = 0;
			}
		}
	}

	for (k = 0; k < height_pet; k++) {
		for (i = 0; i < length_pet * width_pet; i++) {
			if (segment[k][i] != 0) {
				mask_pet[k][i] = 1.0;
			}
			else {
				mask_pet[k][i] = 0.0;
			}
		}
	}

	//dilatation3dHeighteenNeigbours(mask_pet, length_pet, width_pet, height_pet, 1.0, 0.0);

	//outputImagePath = outputPath + "biggest_region_p3.raw";
	//manageRAWFile3D<dataType>(mask_pet, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	//free(countingArray);

	/*==================== Filtering ==========================*/
	
	//Heat explicit so we don't need all the parameters
	//Filter_Parameters parameters; parameters.coef = 1e-4;
	//parameters.h = 6.0; parameters.timeStepSize = 1.0; parameters.timeStepsNum = 50;

	//filterImage(IMAGE_CT, parameters, LINEAR_HEATEQUATION_EXPLICIT);
	//std::string outputImagePath = outputPath + "smoothed.raw";
	//manageRAWFile3D<dataType>(image_ct, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

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

	//for (k = 0; k < height; k++) {
	//	delete[] image_down[k];
	//	//delete[] initial_seg[k];
	//	delete[] mask_down[k];
	//	delete[] distance_down[k];
	//}
	//delete[] image_down;
	////delete[] initial_seg;
	//delete[] mask_down;
	//delete[] distance_down;

	//for (k = 0; k < Height; k++) {
	//	delete[] image_ct[k];
	//	delete[] edge_ct[k];
	//	delete[] image_ct_filtered[k];
	//	delete[] maskThreshold[k];
	//	delete[] distanceMap[k];
	//	//delete[] initial_seg[k];
	//}
	//delete[] image_ct;
	//delete[] edge_ct;
	//delete[] image_ct_filtered;
	//delete[] maskThreshold;
	//delete[] distanceMap;
	////delete[] initial_seg;

	for (k = 0; k < height_pet; k++) {
		delete[] image_pet[k];
		delete[] pet[k];
		delete[] mask_pet[k];
		//delete[] distance_pet[k];
		delete[] segment[k];
		delete[] status[k];
		//delete[] initial_seg[k];
		//delete[] edge_ct_interpolated[k];
		//delete[] ct_interpolated[k];
		//delete[] edge_ct_50_pet_50[k];
		//delete[] edge_ct_75_pet_25[k];
		//delete[] edge_ct_25_pet_75[k];
	}
	delete[] image_pet;
	delete[] pet;
	delete[] mask_pet;
	//delete[] distance_pet;
	delete[] segment;
	delete[] status;
	//delete[] initial_seg;
	//delete[] edge_ct_interpolated;
	//delete[] ct_interpolated;
	//delete[] edge_ct_50_pet_50;
	//delete[] edge_ct_75_pet_25;
	//delete[] edge_ct_25_pet_75;

	//delete[] center_seg;

	free(countingArray);

	return EXIT_SUCCESS;
}
