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
#include "src/shape_registration.h"

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

	Point3D originCT = { -300.0, -230.0, -1022.5 }; // patient 2
	Point3D originPET = { -286.586, -216.586, -1021.40 }; // patient 2
	VoxelSpacing spacingCT = { _sx, _sy, _sz };
	VoxelSpacing spacingPET = { _s_pet, _s_pet, _s_pet };

	//Point3D originCT = { -300.0, -230.0, -1339.2999 };
	//Point3D originPET = { -286.586, -216.586, -1338.800049 };
	//VoxelSpacing spacingCT = { _sx, _sy, _sz };
	//VoxelSpacing spacingPET = { _s_pet, _s_pet, _s_pet };

	//Point3D originCT = { -249.5117, -403.0117, 675.0 }; // patient 6
	//Point3D originPET = { -360.1, -515.421, 677.0 }; // patient 6
	//VoxelSpacing spacingCT = { 0.97656, 0.97656, 1.5 };
	//VoxelSpacing spacingPET = { 3.3, 3.3, 3.0 };
	
	OrientationMatrix orientation;
	orientation.v1 = { 1.0, 0.0, 0.0 }; orientation.v2 = { 0.0, 1.0, 0.0 }; orientation.v3 = { 0.0, 0.0, 1.0 };

	/*========== 2D Labeling for the wiki ======================*/

	//std::string Image2DPath = inputPath + "memboa.raw";
	//dataType * image2D = new dataType[260 * 260];
	//dataType* mask2D = new dataType[260 * 260];
	//manageRAWFile2D<dataType>(image2D, 260, 260, Image2DPath.c_str(), LOAD_DATA, false);
	//Image2DPath = outputPath + "memboa_loaded.raw";
	//manageRAWFile2D<dataType>(image2D, 260, 260, Image2DPath.c_str(), STORE_DATA, false);
	//for (i = 0; i < 260 * 260; i++) {
	//	if (image2D[i] <= 125) {
	//		mask2D[i] = 0.0;
	//	}
	//	else {
	//		mask2D[i] = 1.0;
	//	}
	//}
	//Image2DPath = outputPath + "threshold_memboa.raw";
	//manageRAWFile2D<dataType>(mask2D, 260, 260, Image2DPath.c_str(), STORE_DATA, false);
	//delete[] image2D;
	//delete[] mask2D;

	/*======== Data container for original ct data ============*/
	
	const size_t Height = 406, Width = 512, Length = 512;
	const size_t dim2D = Width * Length;
	
	dataType** image_ct = new dataType * [Height];
	dataType** maskThreshold = new dataType * [Height];
	dataType** distanceMap = new dataType * [Height];
	int** segment = new int* [Height];
	bool** status = new bool* [Height];
	for (k = 0; k < Height; k++) {
		image_ct[k] = new dataType[dim2D];
		maskThreshold[k] = new dataType[dim2D];
		distanceMap[k] = new dataType[dim2D];
		segment[k] = new int[dim2D];
		status[k] = new bool[dim2D];
	}
	if (image_ct == NULL || maskThreshold == NULL || distanceMap == NULL)
		return false;

	//Initialization
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			image_ct[k][i] = 0.0;
			maskThreshold[k][i] = 0.0;
			distanceMap[k][i] = 0.0;
			segment[k][i] = 0;
			status[k][i] = false;
		}
	}

	std::string inputImagePath = inputPath + "patient2_ct.raw";
	manageRAWFile3D<dataType>(image_ct, Length, Width, Height, inputImagePath.c_str(), LOAD_DATA, false);

	/*===================   Find the point with the highest distance ==========*/

	//fill the threshold container
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			maskThreshold[k][i] = image_ct[k][i];
		}
	}

	std::string outputImagePath = outputPath + "loaded.raw";
	//manageRAWFile3D<dataType>(maskThreshold, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	thresholding3dFunctionN(maskThreshold, Length, Width, Height, thres_min, thres_max, 0.0, 1.0);

	outputImagePath = outputPath + "threshold.raw";
	//manageRAWFile3D<dataType>(maskThreshold, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	//Distance Map
	Distance_Map_Params params_dist = { 0.5, 1.0, 0.0, 1000000, 1e-3 };
	computeDistanceMap(distanceMap, maskThreshold, Length, Width, Height, params_dist, FAST_SWEEP);

	outputImagePath = outputPath + "distance.raw";
	//manageRAWFile3D<dataType>(distanceMap, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	Point3D point_found = getPointWithTheHighestValue(distanceMap, Length, Width, Height);
	int imax = (int)point_found.x; 
	int jmax = (int)point_found.y;
	int kmax = (int)point_found.z;

	//point_found = getRealCoordFromImageCoord3D(point_found, originCT, spacingCT, orientation);
	//std::cout << "Coordinates of the approximative centroid : (" << point_found.x << "," << point_found.y << "," << point_found.z << ")" << std::endl;

	//std::cout << "Distance max = " << dist_max << std::endl;
	//std::cout << "Center CT (" << imax << ", " << jmax << ", " << kmax << ")" << std::endl;

	//ball in the Liver
	dataType radius = 5.0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (sqrt((imax - i) * (imax - i) + (jmax - j) * (jmax - j) + (kmax - k) * (kmax - k)) < radius) {
	//				//ball[k][x_new(i, j, Length)] = 1;
	//				image_ct[k][x_new(i, j, Length)] = 0;
	//			}
	//			//else {
	//			//	ball[k][x_new(i, j, Length)] = 0;
	//			//}
	//		}
	//	}
	//}
	//outputImagePath = outputPath + "ball_image_ct.raw";
	//manageRAWFile3D<dataType>(image_ct, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	//Find threshold values around the point with the highest distance
	dataType local_max = 0.0, local_min = 1000000;// radius = 10.0;
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
	thresholding3dFunctionN(maskThreshold, Length, Width, Height, local_min, local_max, 0.0, 1.0);

	/*===================== Labelling CT ========================================*/

	erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	//erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);

	//outputImagePath = outputPath + "afterErosion.raw";
	//manageRAWFile3D<dataType>(maskThreshold, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	labelling3D(maskThreshold, segment, status, Length, Width, Height, 1.0);

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
				maskThreshold[k][i] = 0.0;
			}
			else {
				maskThreshold[k][i] = 1.0;
			}
		}
	}

	////Apply dilatation before saving
	//dataType min_liver = 100000, max_liver = 0.0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		if (segment[k][i] != 0) {
	//			maskThreshold[k][i] = 1.0;
	//			if (image_ct[k][i] < min_liver) {
	//				min_liver = image_ct[k][i];
	//			}
	//			if (image_ct[k][i] > max_liver) {
	//				max_liver = image_ct[k][i];
	//			}
	//		}
	//		else {
	//			maskThreshold[k][i] = 0.0;
	//		}
	//	}
	//}
	////dilatation3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	////dilatation3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	////dilatation3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	////dilatation3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);

	//outputImagePath = outputPath + "biggest_region.raw";
	//manageRAWFile3D<dataType>(maskThreshold, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	//dilatation3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);

	//int liver_number_of_voxel = 0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		//if (segment[k][i] != 0) {
	//		//	maskThreshold[k][i] = image_ct[k][i];
	//		//	liver_number_of_voxel++;
	//		//}
	//		//else {
	//		//	maskThreshold[k][i] = 0.0;
	//		//}
	//		if (maskThreshold[k][i] == 1) {
	//			maskThreshold[k][i] = image_ct[k][i];
	//			liver_number_of_voxel++;
	//		}
	//		else {
	//			maskThreshold[k][i] = 0.0;
	//		}
	//	}
	//}
	//std::cout << "The segmented Liver contains : " << liver_number_of_voxel << " voxels" << std::endl;

	//std::cout << "Stat inside the segmented Live : max = " << max_liver << ", min = " << min_liver << std::endl;
	//outputImagePath = outputPath + "biggest_region_with_one_dilatation.raw";
	//manageRAWFile3D<dataType>(maskThreshold, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	//dataType * center = (dataType*)malloc(3 * sizeof(dataType));
	//centroidImage(maskThreshold, center, Height, Length, Width, 0.0);
	//Point3D Liver_centroid = { center[0], center[1], center[2] };
	//Liver_centroid = getRealCoordFromImageCoord3D(Liver_centroid, originCT, spacingCT, orientation);

	//Point3D centroid_slicer = { 76.50, -59.9137, -588.443 };

	//std::cout << "Coordinates of the centroid : (" << center[0] << "," << center[1] << "," << center[2] << ")" << std::endl;
	//std::cout << "Coordinates of the centroid From Slicer : (" << centroid_slicer.x << "," << centroid_slicer.y << "," << centroid_slicer.z << ")" << std::endl;
	//std::cout << "Coordinates of the centroid in real world coordinates : (" << Liver_centroid.x << "," << Liver_centroid.y << "," << Liver_centroid.z << ")" << std::endl;

	//std::cout << "Distance centroid slicer VS centroid code : " << getPoint3DDistance(Liver_centroid, centroid_slicer) << std::endl;
	//std::cout << "Distance approximative centroid VS centroid Slicer : " << getPoint3DDistance(point_found, centroid_slicer) << std::endl;
	//std::cout << "Distance approximative centroid VS centroid code : " << getPoint3DDistance(point_found, Liver_centroid) << std::endl;

	//free(center);
	free(countingArray);

	/*===================== Data container for orginal PET data =============*/

	const size_t height_pet = 255, length_pet = 144, width_pet = 144; // patient2

	//const size_t height_pet = 318, length_pet = 144, width_pet = 144; // patient3
	 
	//const size_t height_pet = 192, length_pet = 81, width_pet = 120; // patient1
	//const size_t height_pet = 195, length_pet = 128, width_pet = 128; // patient4, patient5
	//const size_t height_pet = 440, length_pet = 220, width_pet = 220; // patient6
	//const size_t height_pet = 440, length_pet = 220, width_pet = 220; // patient7
	
	dataType** image_pet = new dataType * [height_pet];
	dataType** mask_pet = new dataType * [height_pet];
	dataType** inproved_liver = new dataType * [height_pet];
	//int** segment = new int * [height_pet];
	//bool** status = new bool* [height_pet];
	for (k = 0; k < height_pet; k++) {
		image_pet[k] = new dataType[length_pet * width_pet];
		mask_pet[k] = new dataType[length_pet * width_pet];
		inproved_liver[k] = new dataType[length_pet * width_pet];
		//segment[k] = new int[length_pet * width_pet];
		//status[k] = new bool[length_pet * width_pet];
	}
	if (image_pet == NULL || mask_pet == NULL)
		return false;

	for (k = 0; k < height_pet; k++) {
		for (i = 0; i < length_pet * width_pet; i++) {
			image_pet[k][i] = 0.0;
			mask_pet[k][i] = 0.0;
			inproved_liver[k][i] = 0;
			//segment[k][i] = 0;
			//status[k][i] = false;
		}
	}

	//inputImagePath = inputPath + "patient3_pet.raw";
	//manageRAWFile3D<dataType>(image_pet, length_pet, width_pet, height_pet, inputImagePath.c_str(), LOAD_DATA, false);

	//std::string outputImagePath = outputPath + "loadedPET.raw";
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

	//inputImagePath = inputPath + "_seg_func_00.raw";
	//manageRAWFile3D<dataType>(mask_pet, length_pet, width_pet, height_pet, inputImagePath.c_str(), LOAD_DATA, false);

	//inputImagePath = inputPath + "_seg_func_300.raw";
	//manageRAWFile3D<dataType>(inproved_liver, length_pet, width_pet, height_pet, inputImagePath.c_str(), LOAD_DATA, false);

	//for (k = 0; k < height_pet; k++) {
	//	for (i = 0; i < length_pet * width_pet; i++) {
	//		if (inproved_liver[k][i] >= 0.75) {
	//			inproved_liver[k][i] = 1.0;
	//		}
	//		else {
	//			inproved_liver[k][i] = 0.0;
	//		}
	//	}
	//}

	//dataType* center1 = (dataType*)malloc(3 * sizeof(dataType));
	//dataType* center2 = (dataType*)malloc(3 * sizeof(dataType));

	//centroidImage(mask_pet, center1, height_pet, length_pet, width_pet, 0.0);
	//centroidImage(inproved_liver, center2, height_pet, length_pet, width_pet, 0.0);

	//Point3D segment_centroid = { center1[0], center1[1], center1[2] };
	//Point3D inproved_centroid = { center2[0], center2[1], center2[2] };
	//Point3D slicer_centroid = { -72.5352, 71.2067, -510.802 };

	//segment_centroid = getRealCoordFromImageCoord3D(segment_centroid, originPET, spacingPET, orientation);
	//inproved_centroid = getRealCoordFromImageCoord3D(inproved_centroid, originPET, spacingPET, orientation);
	//point_found = getRealCoordFromImageCoord3D(point_found, originCT, spacingCT, orientation);

	//double dist1 = getPoint3DDistance(point_found, segment_centroid);
	//std::cout << " The distance segment liver - point found = " << dist1 << std::endl;

	//double dist2 = getPoint3DDistance(point_found, inproved_centroid);
	//std::cout << " The distance inproved segment centroid - point found = " << dist2 << std::endl;

	//double dist3 = getPoint3DDistance(segment_centroid, inproved_centroid);
	//std::cout << " The distance inproved segment centroid - segment centroid = " << dist3 << std::endl;

	//double dist4 = getPoint3DDistance(point_found, slicer_centroid);
	//std::cout << " The distance found point - slicer centroid = " << dist4 << std::endl;

	//dist4 = getPoint3DDistance(inproved_centroid, slicer_centroid);
	//std::cout << " The distance inproved - slicer centroid = " << dist4 << std::endl;

	//dist4 =  getPoint3DDistance(segment_centroid, slicer_centroid);
	//std::cout << " The distance segment - slicer centroid = " << dist4 << std::endl;

	//dataType** ball1 = new dataType * [height_pet];
	//dataType** ball2 = new dataType * [height_pet];
	//dataType** ball3 = new dataType * [height_pet];
	//dataType** ball4 = new dataType * [height_pet];
	//for (k = 0; k < height_pet; k++) {
	//	ball1[k] = new dataType[length_pet * width_pet];
	//	ball2[k] = new dataType[length_pet * width_pet];
	//	ball3[k] = new dataType[length_pet * width_pet];
	//	ball4[k] = new dataType[length_pet * width_pet];
	//}

	//Point3D current_point = {0.0, 0.0, 0.0};
	//for (k = 0; k < height_pet; k++) {
	//	for (i = 0; i < length_pet; i++) {
	//		for (j = 0; j < width_pet; j++) {
	//			current_point.x = (dataType)i; current_point.y = (dataType)j; current_point.z = (dataType)k;
	//			current_point = getRealCoordFromImageCoord3D(current_point, originPET, spacingPET, orientation);
	//			dist1 = getPoint3DDistance(current_point, point_found);
	//			dist2 = getPoint3DDistance(current_point, segment_centroid);
	//			dist3 = getPoint3DDistance(current_point, inproved_centroid);
	//			dist4 = getPoint3DDistance(current_point, slicer_centroid);
	//			x = x_new(i, j, length_pet);
	//			if (dist1 < 10) {
	//				ball1[k][x] = 1.0;
	//			}
	//			else {
	//				ball1[k][x] = 0.0;
	//			}
	//			if (dist2 < 10) {
	//				ball2[k][x] = 1.0;
	//			}
	//			else {
	//				ball2[k][x] = 0.0;
	//			}
	//			if (dist3 < 10) {
	//				ball3[k][x] = 1.0;
	//			}
	//			else {
	//				ball3[k][x] = 0.0;
	//			}
	//			if (dist4 < 10) {
	//				ball4[k][x] = 1.0;
	//			}
	//			else {
	//				ball4[k][x] = 0.0;
	//			}
	//		}
	//	}
	//}

	//outputImagePath = outputPath + "ball_around_found_point.raw";
	//manageRAWFile3D<dataType>(ball1, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	//outputImagePath = outputPath + "ball_around_segment_centroid.raw";
	//manageRAWFile3D<dataType>(ball2, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	//outputImagePath = outputPath + "ball_around_inproved_centroid.raw";
	//manageRAWFile3D<dataType>(ball3, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	//outputImagePath = outputPath + "ball_around_slicer_centroid.raw";
	//manageRAWFile3D<dataType>(ball4, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	//free(center1);
	//free(center2);
	//for (k = 0; k < height_pet; k++) {
	//	delete[] ball1[k];
	//	delete[] ball2[k];
	//	delete[] ball3[k];
	//	delete[] ball4[k];
	//}
	//delete[] ball1;
	//delete[] ball2;
	//delete[] ball3;
	//delete[] ball4;

	/*==================  Interpolation ========*/

	Image_Data IMAGE_CT;
	IMAGE_CT.imageDataPtr = image_ct;
	IMAGE_CT.height = Height; IMAGE_CT.length = Length; IMAGE_CT.width = Width;
	IMAGE_CT.origin = originCT;
	IMAGE_CT.spacing = spacingCT;
	IMAGE_CT.orientation = orientation;

	Image_Data IMAGE_PET;
	IMAGE_PET.imageDataPtr = image_pet;
	IMAGE_PET.height = height_pet; IMAGE_PET.length = length_pet; IMAGE_PET.width = width_pet;
	IMAGE_PET.origin = originPET;
	IMAGE_PET.spacing = spacingPET;
	IMAGE_PET.orientation = orientation;

	Point3D centerCT = { imax, jmax, kmax }, centerPET = { 0.0, 0.0, 0.0 };
	centerPET = getRealCoordFromImageCoord3D(centerCT, originCT, spacingCT, orientation);
	centerPET = getImageCoordFromRealCoord3D(centerPET, originPET, spacingPET, orientation);
	int ipet = (int)centerPET.x, jpet = (int)centerPET.y, kpet = (int)centerPET.z;

	//inputImagePath = inputPath + "patient3_pet.raw";
	//manageRAWFile3D<dataType>(image_pet, length_pet, width_pet, height_pet, inputImagePath.c_str(), LOAD_DATA, false);

	//Statistics stats = getStatistics(IMAGE_PET, centerPET, 5.0);
	////std::cout << " Stats : max = " << stats.max_data << ", min = " << ", mean = " << stats.mean_data << ", sd = " << stats.sd_data << std::endl;

	//centerPET = getImageCoordFromRealCoord3D(centerPET, originPET, spacingPET, orientation);
	//int ipet = (int)centerPET.x, jpet = (int)centerPET.y, kpet = (int)centerPET.z;
	////std::cout << "Center PET (" << ipet << ", " << jpet << ", " << kpet << ")" << std::endl;

	imageInterpolation3D(IMAGE_CT, IMAGE_PET, TRILINEAR);

	//Point3D current_point = { 0.0, 0.0, 0.0 };
	//double distPoint = 0.0;
	//for (k = 0; k < height_pet; k++) {
	//	for (i = 0; i < length_pet; i++) {
	//		for (j = 0; j < width_pet; j++) {
	//			current_point.x = i; current_point.y = j; current_point.z = k;
	//			distPoint = getPoint3DDistance(centerPET, current_point);
	//			if (distPoint < 5) {
	//				mask_pet[k][x_new(i, j, length_pet)] = 1.0;
	//				image_pet[k][x_new(i, j, length_pet)] = 4000;
	//			}
	//			else {
	//				mask_pet[k][x_new(i, j, length_pet)] = 0.0;
	//			}
	//		}
	//	}
	//}
	//outputImagePath = outputPath + "ball_in_image.raw";
	//manageRAWFile3D<dataType>(image_pet, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	//std::string outputImagePath = outputPath + "interpolated_ct.raw";
	//manageRAWFile3D<dataType>(image_pet, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	/*================= Local threshold for PET ==============*/

	////Find threshold values around the point with the highest distance
	//dataType local_max = 0.0, local_min = 100000000, radius = 5.0;
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

	////outputImagePath = outputPath + "ball_pet.raw";
	////manageRAWFile3D<dataType>(ball_pet, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	//thresholding3dFunctionN(mask_pet, length_pet, width_pet, height_pet, local_min, local_max, 0.0, 1.0);

	//======================  Labelling PET ========================

	//erosion3dHeighteenNeigbours(mask_pet, length_pet, width_pet, height_pet, 1.0, 0.0);
	////erosion3dHeighteenNeigbours(mask_pet, length_pet, width_pet, height_pet, 1.0, 0.0);
	////erosion3dHeighteenNeigbours(mask_pet, length_pet, width_pet, height_pet, 1.0, 0.0);
	////erosion3dHeighteenNeigbours(mask_pet, length_pet, width_pet, height_pet, 1.0, 0.0);
	//erosion3dHeighteenNeigbours(mask_pet, length_pet, width_pet, height_pet, 1.0, 0.0); 
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
	//			//mask_pet[k][i] = 1.0;
	//			mask_pet[k][i] = image_pet[k][i];
	//			liver_volume++;
	//		}
	//		else {
	//			mask_pet[k][i] = 0.0;
	//		}
	//	}
	//}
	//std::cout << "The segmented Liver volume is : " << liver_volume << " voxels" << std::endl;

	//std::string outputImagePath = outputPath + "biggest_region_pet_p3.raw";
	//manageRAWFile3D<dataType>(mask_pet, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	//free(countingArray);

	/*============= Find image centroid ============*/

	

	/*============= Segmentation ================*/

	Point3D* center_seg = new Point3D[1];
	//center_seg->x = idown; center_seg->y = jdown; center_seg->z = kdown;
	//center_seg->x = imax; center_seg->y = jmax; center_seg->z = kmax;
	center_seg->x = ipet; center_seg->y = jpet; center_seg->z = kpet;
	
	dataType** initial_seg = new dataType * [height_pet];
	for (k = 0; k < height_pet; k++) {
		initial_seg[k] = new dataType[length_pet * width_pet];
	}
	for (k = 0; k < height_pet; k++) {
		for (i = 0; i < width_pet * length_pet; i++) {
			initial_seg[k][i] = 0.0;
		}
	}

	IMAGE_CT.imageDataPtr = maskThreshold;
	IMAGE_PET.imageDataPtr = initial_seg;
	imageInterpolation3D(IMAGE_CT, IMAGE_PET, TRILINEAR);

	//generateInitialSegmentationFunctionForMultipleCentres(initial_seg, length_pet, width_pet, height_pet, center_seg, 1.0, 10.0, 1.0);

	rescaleNewRange(image_pet, length_pet, width_pet, height_pet, 0.0, 1.0, 4000.0, 0.0);

	outputImagePath = outputPath + "/segmentation/_seg_func_00.raw";
	//manageRAWFile3D<dataType>(initial_seg, length_pet, width_pet, height_pet, outputImagePath.c_str(), STORE_DATA, false);

	Image_Data ImageToBeSegmented;
	ImageToBeSegmented.imageDataPtr = image_pet;
	ImageToBeSegmented.height = height_pet; ImageToBeSegmented.length = length_pet; ImageToBeSegmented.width = width_pet;

	Segmentation_Parameters seg_params;
	seg_params.tau = 0.1; seg_params.h = 1.0; seg_params.omega_c = 1.4; seg_params.mod = 10;
	seg_params.maxNoGSIteration = 100; seg_params.coef = 50000; seg_params.gauss_seidelTolerance = 1e-3;
	seg_params.maxNoOfTimeSteps = 100; seg_params.segTolerance = 1e-6; seg_params.eps2 = 1e-6;
	seg_params.coef_conv = 1.0; seg_params.coef_dif = 0.010;

	//Smoothing by heat explicit so we don't need all the parameters
	Filter_Parameters parameters; parameters.coef = 1e-4;
	parameters.h = seg_params.h; parameters.timeStepSize = parameters.h/6.0; parameters.timeStepsNum = 1; //10
	//Filter_Parameters parameters; parameters.coef = 1e-4; parameters.edge_detector_coefficient = 1000; parameters.eps2 = 1e-4;
	//parameters.h = 6.0; parameters.maxNumberOfSolverIteration = 100; parameters.omega_c = 1.4; parameters.p = 1;
	//parameters.sigma = 0.1; parameters.timeStepSize = 1.0; parameters.timeStepsNum = 1; parameters.tolerance = 1e-6;

	unsigned char segmentPath[] = "C:/Users/Konan Allaly/Documents/Tests/output/segmentation/";

	//subsurfSegmentation(ImageToBeSegmented, initial_seg, seg_params, parameters, center_seg, 1, segmentPath);
	generalizedSubsurfSegmentation(ImageToBeSegmented, initial_seg, seg_params, parameters, center_seg, 1, segmentPath);

	delete[] center_seg;

	for (k = 0; k < Height; k++) {
		delete[] image_ct[k];
		delete[] maskThreshold[k];
		delete[] distanceMap[k];
		delete[] segment[k];
		delete[] status[k];
	}
	delete[] image_ct;
	delete[] maskThreshold;
	delete[] distanceMap;
	delete[] segment;
	delete[] status;

	for (k = 0; k < height_pet; k++) {
		delete[] image_pet[k];
		delete[] mask_pet[k];
		delete[] initial_seg[k];
		delete[] inproved_liver[k];
		//delete[] segment[k];
		//delete[] status[k];
	}
	delete[] image_pet;
	delete[] mask_pet;
	delete[] initial_seg;
	delete[] inproved_liver;
	//delete[] segment;
	//delete[] status;

	return EXIT_SUCCESS;
}
