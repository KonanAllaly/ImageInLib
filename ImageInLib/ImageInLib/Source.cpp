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

//#define thres_min 995
//#define thres_max 1213

#define orig_thres_min -29
#define orig_thres_max 189
#define orig_min -3024

#define thres_min orig_thres_min - orig_min
#define thres_max orig_thres_max - orig_min

int main() {

	int i, j, k;
	
	std::string inputPath = "C:/Users/Konan Allaly/Documents/Tests/input/vtk/";
	std::string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";

	//Point3D originCT = { -300.0, -230.0, -1022.5 }; // patient 2
	//Point3D originPET = { -286.586, -216.586, -1021.40 }; // patient 2
	//VoxelSpacing spacingCT = { _sx, _sy, _sz };
	//VoxelSpacing spacingPET = { _s_pet, _s_pet, _s_pet };

	//Point3D originCT = { -300.0, -230.0, -1339.2999 };
	//Point3D originPET = { -286.586, -216.586, -1338.800049 };
	//VoxelSpacing spacingCT = { _sx, _sy, _sz };
	//VoxelSpacing spacingPET = { _s_pet, _s_pet, _s_pet };

	//Point3D originCT = { -249.5117, -403.0117, 675.0 }; // patient 6
	//Point3D originPET = { -360.1, -515.421, 677.0 }; // patient 6
	//VoxelSpacing spacingCT = { 0.97656, 0.97656, 1.5 };
	//VoxelSpacing spacingPET = { 3.3, 3.3, 3.0 };

	//Point3D originCT = { -249.5117, -446.01, -1123.5 }; // patient 7
	//Point3D originPET = { -362.832, -557.989, -1123.5 }; // patient 7
	//VoxelSpacing spacingCT = { 0.97656, 0.97656, 1.5 };
	//VoxelSpacing spacingPET = { 3.906, 3.906, 5.9895 };
	
	//OrientationMatrix orientation;
	//orientation.v1 = { 1.0, 0.0, 0.0 }; orientation.v2 = { 0.0, 1.0, 0.0 }; orientation.v3 = { 0.0, 0.0, 1.0 };

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
	
	//int Height = 330, Width = 512, Length = 512;
	//dataType** image_ct = new dataType * [Height];
	//for (k = 0; k < Height; k++) {
	//	image_ct[k] = new dataType[Length * Width];
	//}
	//if (image_ct == NULL)
	//	return false;
	////Initialization
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length * Width; i++) {
	//		image_ct[k][i] = 0;
	//	}
	//}
	//dataType** image_ct = new dataType * [Height];
	//dataType** maskThreshold = new dataType * [Height];
	//dataType** distanceMap = new dataType * [Height];
	//int** segment = new int* [Height];
	//bool** status = new bool* [Height];
	//for (k = 0; k < Height; k++) {
	//	image_ct[k] = new dataType[dim2D];
	//	maskThreshold[k] = new dataType[dim2D];
	//	distanceMap[k] = new dataType[dim2D];
	//	segment[k] = new int[dim2D];
	//	status[k] = new bool[dim2D];
	//}
	//if (image_ct == NULL || maskThreshold == NULL || distanceMap == NULL)
	//	return false;
	////Initialization
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		image_ct[k][i] = 0.0;
	//		maskThreshold[k][i] = 0.0;
	//		distanceMap[k][i] = 0.0;
	//		segment[k][i] = 0;
	//		status[k][i] = false;
	//	}
	//
	//typedef struct {
	//	double spacing[3];
	//	double origin[3];
	//	int dimensions[3];
	//	VtkDataType vDataType;
	//	dataType** dataPointer;
	//	vtkOperation operation;
	//} Vtk_File_Info;

	Vtk_File_Info* ctDataCointainer = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info)); 
	ctDataCointainer->operation = copyFrom;
	//infoP4->spacing[0] = 0.976562; infoP4->spacing[1] = 0.976462; infoP4->spacing[2] = 1.5;
	//infoP4->origin[0] = -250; infoP4->origin[1] = -250; infoP4->origin[2] = -876.5;
	//infoP4->dimensions[2] = Height; infoP4->dimensions[1] = Length; infoP4->dimensions[0] = Width;
	//infoP4->dataPointer = image_ct;  
	//infoP4->vDataType = dta_Dbl;

	std::string inputImagePath = inputPath + "Patient4_ct.vtk";
	readVtkFile(inputImagePath.c_str(), ctDataCointainer);

	int Height = ctDataCointainer->dimensions[2];
	int Length = ctDataCointainer->dimensions[1];
	int Width = ctDataCointainer->dimensions[0];
	
	//std::cout << "Height = " << Height << ", Width = " << Width << ", Length = " << Length << std::endl;

	dataType min_data = 100000000.0, max_data = -100000000.0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length * Width; i++) {
			if (ctDataCointainer->dataPointer[k][i] > max_data) {
				max_data = ctDataCointainer->dataPointer[k][i];
			}
			if (ctDataCointainer->dataPointer[k][i] < min_data) {
				min_data = ctDataCointainer->dataPointer[k][i];
			}
		}
	}
	std::cout << "Min = " << min_data << ", Max = " << max_data << std::endl;
	std::cout << "Spacing : (" << ctDataCointainer->spacing[0] << "," << ctDataCointainer->spacing[1] << ","  << ctDataCointainer->spacing[2] << ")" << std::endl;
	std::cout << "Image origin : (" << ctDataCointainer->origin[0] << "," << ctDataCointainer->origin[1] << "," << ctDataCointainer->origin[2] << ")" << std::endl;

	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length * Width; i++) {
			ctDataCointainer->dataPointer[k][i] = ctDataCointainer->dataPointer[k][i] + (-1) * min_data;
		}
	}

	//std::string outputImagePath = outputPath + "patient4_ct_new.raw";
	//manageRAWFile3D<dataType>(ctDataCointainer->dataPointer, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	/*===================   Find the point with the highest distance ==========*/

	int dim2D = Length * Width;
	dataType** maskThreshold = new dataType * [Height];
	dataType** distanceMap = new dataType * [Height];
	dataType** ball = new dataType * [Height];
	int** segment = new int* [Height];
	bool** status = new bool* [Height];
	for (k = 0; k < Height; k++) {
		maskThreshold[k] = new dataType[dim2D];
		distanceMap[k] = new dataType[dim2D];
		ball[k] = new dataType[dim2D];
		segment[k] = new int[dim2D];
		status[k] = new bool[dim2D];

		if (maskThreshold[k] == NULL || distanceMap[k] == NULL || ball[k] == NULL || segment[k] == NULL || status[k] == NULL)
			return false;
	}
	if (maskThreshold == NULL || distanceMap == NULL || ball == NULL || segment == NULL || status == NULL)
		return false;

	//fill the threshold container
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			maskThreshold[k][i] = ctDataCointainer->dataPointer[k][i];
			distanceMap[k][i] = 0.0;
			ball[k][i] = 0.0;
			segment[k][i] = 0;
			status[k][i] = false;
		}
	}

	//std::string outputImagePath = outputPath + "patient5_ct_new.raw";
	//manageRAWFile3D<dataType>(maskThreshold, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	thresholding3dFunctionN(maskThreshold, Length, Width, Height, thres_min, thres_max, 0.0, 1.0);

	//outputImagePath = outputPath + "threshold.raw";
	//manageRAWFile3D<dataType>(maskThreshold, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	//Distance Map
	Distance_Map_Params params_dist = { 0.5, 1.0, 0.0, 1000000, 1e-3 };
	computeDistanceMap(distanceMap, maskThreshold, Length, Width, Height, params_dist, FAST_SWEEP);

	//outputImagePath = outputPath + "distance.raw";
	//manageRAWFile3D<dataType>(distanceMap, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	Point3D point_found = getPointWithTheHighestValue(distanceMap, Length, Width, Height);
	int imax = (int)point_found.x; 
	int jmax = (int)point_found.y;
	int kmax = (int)point_found.z;

	//point_found = getRealCoordFromImageCoord3D(point_found, originCT, spacingCT, orientation);
	//std::cout << "Coordinates of the approximative centroid : (" << point_found.x << "," << point_found.y << "," << point_found.z << ")" << std::endl;

	//ball in the Liver
	int x = 0;
	dataType radius = 10.0, dist = 0.0;
	Point3D current_point;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				x = x_new(i, j, Length);
				current_point.x = (dataType)i; current_point.y = (dataType)j; current_point.z = (dataType)k;
				dist = getPoint3DDistance(current_point, point_found);
				if (dist <= radius) {
					ball[k][x] = 1.0;
				}
				else {
					ball[k][x] = 0.0;
				}
			}
		}
	}

	//outputImagePath = outputPath + "ball_in_liver_ct_p5.raw";
	//manageRAWFile3D<dataType>(ball, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	//Find threshold values around the point with the highest distance
	dataType local_max = 0.0, local_min = 1000000;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				x = x_new(i, j, Length);
				current_point.x = (dataType)i; current_point.y = (dataType)j; current_point.z = (dataType)k;
				dist = getPoint3DDistance(current_point, point_found);
				if ( dist < radius) {
					
					if (ctDataCointainer->dataPointer[k][x] > local_max) {
						local_max = ctDataCointainer->dataPointer[k][x];
					}
					if (ctDataCointainer->dataPointer[k][x] < local_min) {
						local_min = ctDataCointainer->dataPointer[k][x];
					}
				}
			}
		}
	}
	std::cout << "Local max = " << local_max << ", " << " Local min = " << local_min << std::endl;

	//fill the threshold container
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			maskThreshold[k][i] = ctDataCointainer->dataPointer[k][i];
		}
	}
	thresholding3dFunctionN(maskThreshold, Length, Width, Height, local_min, local_max, 0.0, 1.0);

	/*===================== Labelling CT ========================================*/

	erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	//erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	//erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	//erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);

	//outputImagePath = outputPath + "afterErosion_p5.raw";
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

	//outputImagePath = outputPath + "biggest_region_p5.raw";
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
	//free(countingArray);

	/*===================== Data container interpolated data =============*/
	 
	int height = 83, width = 128, length = 128;
	
	dataType** imageDown = new dataType * [height];
	for (k = 0; k < height; k++) {
		imageDown[k] = new dataType[length * width];
	}
	if (imageDown == NULL)
		return false;

	for (k = 0; k < height; k++) {
		for (i = 0; i < length * width; i++) {
			imageDown[k][i] = 0.0;
		}
	}

	/*==================  Interpolation ========*/

	OrientationMatrix orientation;
	orientation.v1 = { 1.0, 0.0, 0.0 }; orientation.v2 = { 0.0, 1.0, 0.0 }; orientation.v3 = { 0.0, 0.0, 1.0 };

	Image_Data IMAGE_CT;
	IMAGE_CT.imageDataPtr = ctDataCointainer->dataPointer;
	IMAGE_CT.height = Height; IMAGE_CT.length = Length; IMAGE_CT.width = Width;
	IMAGE_CT.origin.x = ctDataCointainer->origin[0]; IMAGE_CT.origin.y = ctDataCointainer->origin[1]; IMAGE_CT.origin.z = ctDataCointainer->origin[2];
	IMAGE_CT.spacing.sx = ctDataCointainer->spacing[0]; IMAGE_CT.spacing.sy = ctDataCointainer->spacing[1]; IMAGE_CT.spacing.sz = ctDataCointainer->spacing[2];
	IMAGE_CT.orientation = orientation;

	Image_Data IMAGE_DOWN;
	IMAGE_DOWN.imageDataPtr = imageDown;
	IMAGE_DOWN.height = height; IMAGE_DOWN.length = length; IMAGE_DOWN.width = width;
	IMAGE_DOWN.origin = IMAGE_CT.origin;
	IMAGE_DOWN.spacing = {3.906248, 3.906348, 10.0};
	IMAGE_DOWN.orientation = orientation;

	Point3D centerCT = { imax, jmax, kmax }, centerDown = { 0.0, 0.0, 0.0 };
	 
	centerDown = getRealCoordFromImageCoord3D(centerCT, IMAGE_CT.origin, IMAGE_CT.spacing, orientation);
	centerDown = getImageCoordFromRealCoord3D(centerDown, IMAGE_DOWN.origin, IMAGE_DOWN.spacing, orientation);
	int idown = (int)centerDown.x, jdown = (int)centerDown.y, kdown = (int)centerDown.z;

	imageInterpolation3D(IMAGE_CT, IMAGE_DOWN, TRILINEAR);

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

	//std::string outputImagePath = outputPath + "interpolated_ct_p4.raw";
	//manageRAWFile3D<dataType>(imageDown, length, width, height, outputImagePath.c_str(), STORE_DATA, false);

	/*============= Segmentation ================*/

	Point3D* center_seg = new Point3D[1];
	center_seg->x = idown; center_seg->y = jdown; center_seg->z = kdown;
	
	dataType** initial_seg = new dataType * [height];
	for (k = 0; k < height; k++) {
		initial_seg[k] = new dataType[length * width];
	}
	for (k = 0; k < height; k++) {
		for (i = 0; i < width * length; i++) {
			initial_seg[k][i] = 0.0;
		}
	}

	IMAGE_CT.imageDataPtr = maskThreshold;
	IMAGE_DOWN.imageDataPtr = initial_seg;
	imageInterpolation3D(IMAGE_CT, IMAGE_DOWN, TRILINEAR);

	//generateInitialSegmentationFunctionForMultipleCentres(initial_seg, length_pet, width_pet, height_pet, center_seg, 1.0, 10.0, 1.0);
	rescaleNewRange(imageDown, length, width, height, 0.0, 1.0, 6095.0, 0.0);

	std::string outputImagePath = outputPath + "/segmentation/_seg_func_00.raw";
	manageRAWFile3D<dataType>(initial_seg, length, width, height, outputImagePath.c_str(), STORE_DATA, false);

	Image_Data ImageToBeSegmented;
	ImageToBeSegmented.imageDataPtr = imageDown;
	ImageToBeSegmented.height = height; ImageToBeSegmented.length = length; ImageToBeSegmented.width = width;

	Segmentation_Parameters seg_params;
	seg_params.tau = 1.0; seg_params.h = 1.0; seg_params.omega_c = 1.4; seg_params.mod = 2;
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

	/*================= Free memory ==========================*/

	delete[] center_seg;
	free(ctDataCointainer);
	free(countingArray);
	for (k = 0; k < Height; k++) {
		delete[] maskThreshold[k];
		delete[] distanceMap[k];
		delete[] ball[k];
		delete[] segment[k];
		delete[] status[k];
	}
	delete[] maskThreshold;
	delete[] distanceMap;
	delete[] ball;
	delete[] segment;
	delete[] status;

	for (k = 0; k < height; k++) {
		delete[] imageDown[k];
		delete[] initial_seg[k];
	}
	delete[] imageDown;
	delete[] initial_seg;

	return EXIT_SUCCESS;
}
