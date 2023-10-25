// Source code for fast marching and path finding (11/09/2023)

#include <iostream>
#include <vector>
#include <string.h> 

#include<template_functions.h>
#include "../src/vtk_params.h"
#include "common_vtk.h"
#include "distanceForPathFinding.h"
#include "common_functions.h"
#include "../src/imageInterpolation.h"
#include "filtering.h"
#include "../src/segmentation3D_subsurf.h"
#include "../src/segmentation3d_gsubsurf.h"

using namespace std;

int main() {

	string inputPath = "C:/Users/Konan Allaly/Documents/Tests/input/vtk/";
	string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";

	OrientationMatrix orientation;
	orientation.v1 = { 1.0, 0.0, 0.0 }; orientation.v2 = { 0.0, 1.0, 0.0 }; orientation.v3 = { 0.0, 0.0, 1.0 };

	Vtk_File_Info* ctContainer = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	ctContainer->operation = copyFrom;

	std::string inputImagePath = inputPath + "Patient7_ct.vtk";
	//std::string inputImagePath = inputPath + "aortaP2_real_origin.vtk";
	readVtkFile(inputImagePath.c_str(), ctContainer);

	int Height = ctContainer->dimensions[2];
	int Length = ctContainer->dimensions[1];
	int Width = ctContainer->dimensions[0];
	int dim2D = Length * Width;
	cout << "Input image dim : " << ctContainer->dimensions[0] << " x " << ctContainer->dimensions[1] << " x " << ctContainer->dimensions[2] << "" << endl;

	cout << "Image origin : (" << ctContainer->origin[0] << ", " << ctContainer->origin[1] << ", " << ctContainer->origin[2] << ")" << endl;
	Point3D ctOrigin = { ctContainer->origin[0], ctContainer->origin[1], ctContainer->origin[2] };
	VoxelSpacing ctSpacing = { ctContainer->spacing[0], ctContainer->spacing[1], ctContainer->spacing[2] };

	int i = 0, j = 0, k = 0, x = 0, x2 = 0;
	
	/*================  Data container for fast marching and path finding ============*/

	dataType** imageFiltered = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageFiltered[k] = new dataType[dim2D]{0};
		if (imageFiltered[k] == NULL)
			return false;
	}
	if (imageFiltered == NULL)
		return false;

	inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/raw/filteredGMC_p7.raw";
	manageRAWFile3D<dataType>(imageFiltered, Length, Width, Height, inputImagePath.c_str(), LOAD_DATA, false);

	Image_Data CT; CT.imageDataPtr = imageFiltered;
	CT.height = Height; CT.length = Length; CT.width = Width;
	CT.orientation = orientation; CT.origin = ctOrigin; CT.spacing = ctSpacing;

	/*================== Image Filtering =============================================*/

	//Image_Data ImageData; ImageData.imageDataPtr = imageData; ImageData.height = Height;
	//ImageData.length = Length; ImageData.width = Width;
	//Filter_Parameters filter_parameters; filter_parameters.edge_detector_coefficient = 1;
	//filter_parameters.h = 1.0; filter_parameters.timeStepSize = 1.2; filter_parameters.eps2 = 1e-6;
	//filter_parameters.omega_c = 1.5; filter_parameters.tolerance = 1e-2; filter_parameters.maxNumberOfSolverIteration = 100;
	//filter_parameters.timeStepsNum = 1; filter_parameters.coef = 1e-6; 
	//filterImage(ImageData, filter_parameters, GEODESIC_MEAN_CURVATURE_FILTER);
	//outputImagePath = outputPath + "filteredGMC.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	//=================== Interpolation ====================================

	const size_t height_new = (size_t)((ctContainer->dimensions[2] * ctContainer->spacing[2]) / ctContainer->spacing[0]);
	cout << "Interpolated height : " << height_new << endl;
	cout << "spacing : " << ctContainer->spacing[0] << endl;

	dataType** imageInterpolated = new dataType * [height_new];
	dataType** distance = new dataType * [height_new];
	dataType** path = new dataType * [height_new];
	dataType** potential = new dataType * [height_new];
	dataType** image_mean = new dataType * [height_new];
	for (k = 0; k < height_new; k++) {
		imageInterpolated[k] = new dataType[dim2D]{0};
		distance[k] = new dataType[dim2D]{0};
		path[k] = new dataType[dim2D]{0};
		potential[k] = new dataType[dim2D]{0};
		image_mean[k] = new dataType[dim2D]{0};
		if (imageInterpolated[k] == NULL || distance[k] == NULL || path[k] == NULL || potential[k] == NULL || image_mean[k] == NULL)
			return false;
	}
	if (imageInterpolated == NULL || distance == NULL || path == NULL || potential == NULL || image_mean == NULL)
		return false;

	Image_Data Interpolated; Interpolated.imageDataPtr = imageInterpolated;
	Interpolated.height = height_new; Interpolated.length = Length; Interpolated.width = Width;
	Interpolated.orientation = orientation; Interpolated.origin = ctOrigin;
	Interpolated.spacing.sx = ctSpacing.sx; Interpolated.spacing.sy = ctSpacing.sy; Interpolated.spacing.sz = ctSpacing.sx;
	imageInterpolation3D(CT, Interpolated, NEAREST_NEIGHBOR);

	//inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/initial_seg.raw";
	//inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/path_p3.raw";
	//manageRAWFile3D<dataType>(distance, Length, Width, height_new, inputImagePath.c_str(), LOAD_DATA, false);

	////create initial segment
	//size_t l, m, n;
	//for (k = 0; k < height_new; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			x = x_new(i, j, Length);
	//			if (distance[k][x] == 1.0) {
	//				Point3D current_point = { (dataType)i, (dataType)j, (dataType)k };
	//				for (l = 0; l < height_new; l++) {
	//					for (m = 0; m < Length; m++) {
	//						for (n = 0; n < Width; n++) {
	//							Point3D temporary_point = { (dataType)m, (dataType)n, (dataType)l };
	//							double dist = getPoint3DDistance(current_point, temporary_point);
	//							if (dist <= 3.0) {
	//								path[l][x_new(m, n, Length)] = 1.0;
	//							}
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}
	//}

	//outputImagePath = outputPath + "initial_seg.raw";
	//manageRAWFile3D<dataType>(path, Length, Width, height_new, outputImagePath.c_str(), STORE_DATA, false);

	//inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/output/interpolated.raw";
	//manageRAWFile3D<dataType>(Interpolated.imageDataPtr, Length, Width, height_new, inputImagePath.c_str(), STORE_DATA, false);

	//================== Downsampling =====================================

	//const size_t height = (size_t)(height_new / 4);
	//cout << "heigth down : " << height << endl;
	//const size_t length = (size_t)(Length / 4);
	//const size_t width = (size_t)(Width / 4);

	//dataType** downsampled_image = new dataType * [height];
	//for (k = 0; k < height; k++) {
	//	downsampled_image[k] = new dataType[length * width];
	//}

	////Initialization
	//for (k = 0; k < height; k++) {
	//	for (i = 0; i < length * width; i++) {
	//		downsampled_image[k][i] = 0.0;
	//	}
	//}

	//Image_Data Downsampled; 
	//Downsampled.imageDataPtr = downsampled_image; 
	//Downsampled.length = length; Downsampled.width = width; Downsampled.height = height;
	//Downsampled.origin = ctOrigin; Downsampled.orientation = orientation;
	//Downsampled.spacing.sx = ctSpacing.sx * 4; Downsampled.spacing.sy = ctSpacing.sy * 4; Downsampled.spacing.sz = ctSpacing.sx * 4;
	////cout << "spacing down : " << Downsampled.spacing.sx << endl;
	//
	//imageInterpolation3D(Interpolated, Downsampled, TRILINEAR);

	//outputImagePath = outputPath + "path_down.raw";
	//manageRAWFile3D<dataType>(downsampled_image, length, width, height, outputImagePath.c_str(), STORE_DATA, false);

	//================== Fast Marching and Path finding ===================

	inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/interpolated/patient7/Im_mean_r3.raw";
	manageRAWFile3D<dataType>(image_mean, Length, Width, height_new, inputImagePath.c_str(), LOAD_DATA, false);

	////patient 2
	//Point3D final_point = { 262, 254, 250 }; // top 
	//Point3D initial_point = { 263, 257, 146 }; // bottom
	////Point3D middle_point = { 285, 284, 278 }; // middle

	//////Patient3
	////Point3D initial_point = { 259, 255, 244 };
	////Point3D final_point = { 260, 248, 348 };
	////Point3D middle_point = { 285, 301, 372 };

	////Patient4
	////Point3D middle_point = { 292, 277, 240 };
	//Point3D final_point = { 255, 224, 222 };
	//Point3D initial_point = { 269, 231, 111 };

	////Patient5
	////Point3D middle_point = { 260, 307, 247 };
	//Point3D final_point = { 273, 229, 230 };
	//Point3D initial_point = { 281, 229, 132 };

	////Patient6
	////Point3D middle_point = { 285, 266, 662 };
	////Point3D final_point = { 243, 221, 637 }; // test 1 
	//Point3D final_point = { 241, 214, 627 }; // test 2
	//Point3D initial_point = { 265, 244, 486 };// test 1

	//Patient7
	//Point3D middle_point = { 277, 346, 476 };
	Point3D final_point = { 252, 288, 442 };
	Point3D initial_point = { 254, 299, 261 };

	//Get corresponding point in interpolated image
	final_point = getRealCoordFromImageCoord3D(final_point, ctOrigin, ctSpacing, orientation);
	final_point = getImageCoordFromRealCoord3D(final_point, ctOrigin, Interpolated.spacing, orientation);
	final_point.x = (size_t)final_point.x; final_point.y = (size_t)final_point.y; final_point.z = (size_t)final_point.z;

	//Get corresponding point in interpolated image
	initial_point = getRealCoordFromImageCoord3D(initial_point, ctOrigin, ctSpacing, orientation);
	initial_point = getImageCoordFromRealCoord3D(initial_point, ctOrigin, Interpolated.spacing, orientation);
	initial_point.x = (size_t)initial_point.x; initial_point.y = (size_t)initial_point.y; initial_point.z = (size_t)initial_point.z;

	Point3D* seedPoints = new Point3D[2]; 
	seedPoints[0] = initial_point;
	seedPoints[1] = final_point;

	Potential_Parameters parameters;
	parameters.K = 0.005; parameters.epsilon = 0.01;
	parameters.c_ct = 0.0; parameters.c_pet = 0.0;
	parameters.c_max = 0.0; parameters.c_min = 0.0; parameters.c_mean = 1.0; parameters.c_sd = 0.0;

	findPathFromOneGivenPoint(Interpolated, image_mean, path, seedPoints, parameters);
	//findPathBetweenTwoGivenPoints(Interpolated, image_mean, path, seedPoints, parameters);
	//string outputImagePath = outputPath + "path.raw";
	//manageRAWFile3D<dataType>(path, Length, Width, height_new, outputImagePath.c_str(), STORE_DATA, false);

	//================== Automatic Cropping =========================================

	//imageMetaData croppedImageMetaData;
	//croppedImageMetaData = croppVolume(CT, 10);
	//cout << "Cropped image origin : (" << croppedImageMetaData.origin.x << ", " << croppedImageMetaData.origin.y << ", " << croppedImageMetaData.origin.z << ")" << endl;

	//croppedImageMetaData.origin = getImageCoordFromRealCoord3D(croppedImageMetaData.origin, ctOrigin, ctSpacing, orientation);
	//size_t i0 = (size_t)croppedImageMetaData.origin.x, j0 = (size_t)croppedImageMetaData.origin.y, k0 = (size_t)croppedImageMetaData.origin.z;

	//size_t height = croppedImageMetaData.dimension[2];
	//size_t length = croppedImageMetaData.dimension[0];
	//size_t width = croppedImageMetaData.dimension[1];
	//cout << "Cropped image dim : " << length << " x " << width << " x " << height << "" << endl;

	//dataType** croppImage = new dataType * [height];
	//for (k = 0; k < height; k++) {
	//	croppImage[k] = new dataType[length * width]{0};
	//}

	//size_t k2, i2, j2;
	//for (k = 0, k2 = k0; k < height; k++, k2++) {
	//	for (i = 0, i2 = i0; i < length; i++, i2++) {
	//		for (j = 0, j2 = j0; j < width; j++, j2++) {
	//			croppImage[k][x_new(i, j, length)] = ctContainer->dataPointer[k2][x_new(i2, j2, Length)];
	//		}
	//	}
	//}

	//string outputImagePath = outputPath + "croppedImage.raw";
	//manageRAWFile3D<dataType>(croppImage, length, width, height, outputImagePath.c_str(), STORE_DATA, false);

	//for (k = 0; k < height; k++) {
	//	delete[] croppImage[k];
	//}
	//delete[] croppImage;

	//================== Downsampling =====================================

	//const size_t height_down = (size_t)(height / 4);
	//const size_t length_down = (size_t)(length / 2);
	//const size_t width_down = (size_t)(width / 2);

	//cout << "down sampled height : " << height_down << endl;
	//cout << "down sampled lenght : " << length_down << endl;
	//cout << "down sampled width : " << width_down << endl;

	//dataType** downsampled_image = new dataType * [height_down];
	//dataType** initialSeg = new dataType * [height_down];
	//for (k = 0; k < height_down; k++) {
	//	downsampled_image[k] = new dataType[length_down * width_down]{0};
	//	initialSeg[k] = new dataType[length_down * width_down]{0};
	//}

	//Image_Data Downsampled; 
	//Downsampled.length = length_down; Downsampled.width = width_down; Downsampled.height = height_down;
	//Downsampled.origin = cropped.origin; Downsampled.orientation = orientation;
	//Downsampled.spacing.sx = cropped.spacing.sx * 2; Downsampled.spacing.sy = cropped.spacing.sy * 2; Downsampled.spacing.sz = cropped.spacing.sz * 4;
	//cout << "new spacing sx = sy =  " << Downsampled.spacing.sx << ", sz = " << Downsampled.spacing.sz << endl;

	//Downsampled.imageDataPtr = downsampled_image;
	//cropped.imageDataPtr = croppedImage;
	//imageInterpolation3D(cropped, Downsampled, NEAREST_NEIGHBOR);

	//outputImagePath = outputPath + "new_down_image.raw";
	////manageRAWFile3D<dataType>(downsampled_image, length_down, width_down, height_down, outputImagePath.c_str(), STORE_DATA, false);

	//Downsampled.imageDataPtr = initialSeg;
	//cropped.imageDataPtr = croppedPath;
	//imageInterpolation3D(cropped, Downsampled, NEAREST_NEIGHBOR);
	//
	//outputImagePath = outputPath + "down_path.raw";
	////manageRAWFile3D<dataType>(initialSeg, length_down, width_down, height_down, outputImagePath.c_str(), STORE_DATA, false);

	//outputImagePath = outputPath + "path_down.raw";
	//manageRAWFile3D<dataType>(downsampled_image, length, width, height, outputImagePath.c_str(), STORE_DATA, false);

	//================== Segmentation =====================================

	//Point3D* center_seg = new Point3D[1];
	//center_seg->x = 0; center_seg->y = 0; center_seg->z = 0; // not used in the current version of the code

	//dataType** initial_seg = new dataType * [height];
	//for (k = 0; k < height; k++) {
	//	initial_seg[k] = new dataType[length * width]{0};
	//}

	//inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/path_down.raw";
	//manageRAWFile3D<dataType>(initial_seg, length, width, height, inputImagePath.c_str(), LOAD_DATA, false);

	//outputImagePath = outputPath + "/segmentation/_seg_func_00.raw";
	//manageRAWFile3D<dataType>(initialSeg, length_down, width_down, height_down, outputImagePath.c_str(), STORE_DATA, false);

	//Image_Data ImageToBeSegmented;
	//ImageToBeSegmented.imageDataPtr = downsampled_image;
	//ImageToBeSegmented.height = height_down; ImageToBeSegmented.length = length_down; ImageToBeSegmented.width = width_down;

	//Segmentation_Parameters seg_params;
	//seg_params.tau = 1.0; seg_params.h = 1.0; seg_params.omega_c = 1.4; seg_params.mod = 10;
	//seg_params.maxNoGSIteration = 100; seg_params.coef = 100000; seg_params.gauss_seidelTolerance = 1e-3;
	//seg_params.maxNoOfTimeSteps = 500; seg_params.segTolerance = 1e-6; seg_params.eps2 = 1e-6;
	//seg_params.coef_conv = 1.0; seg_params.coef_dif = 1.0;

	////Smoothing by heat equation so we don't need all the parameters
	//Filter_Parameters parameters; parameters.coef = 1e-4;
	//parameters.h = seg_params.h; parameters.timeStepSize = parameters.h / 16.0; parameters.timeStepsNum = 1; //10
	////Filter_Parameters parameters; parameters.coef = 1e-4; parameters.edge_detector_coefficient = 1000; parameters.eps2 = 1e-4;
	////parameters.h = 6.0; parameters.maxNumberOfSolverIteration = 100; parameters.omega_c = 1.4; parameters.p = 1;
	////parameters.sigma = 0.1; parameters.timeStepSize = 1.0; parameters.timeStepsNum = 1; parameters.tolerance = 1e-6;

	//double firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
	//unsigned char segmentPath[] = "C:/Users/Konan Allaly/Documents/Tests/output/segmentation/";
	//generalizedSubsurfSegmentation(ImageToBeSegmented, initialSeg, seg_params, parameters, center_seg, 1, segmentPath);
	//double secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);

	//cout << "Time : " << (secondCpuTime - firstCpuTime) / 60 << " minutes" << endl;

	//========================= Free memory ==================================

	delete[] seedPoints;

	free(ctContainer);

	//for (k = 0; k < height_down; k++) {
	//	delete[] downsampled_image[k];
	//	delete[] initialSeg[k];
	//}
	//delete[] downsampled_image;
	//delete[] initialSeg;

	//for (k = 0; k < height; k++) {
	//	delete[] croppedImage[k];
	//	delete[] croppedPath[k];
	//}
	//delete[] croppedImage;
	//delete[] croppedPath;

	//delete[] center_seg;
	
	for (k = 0; k < Height; k++) {
		delete[] imageFiltered[k];
	}
	delete[] imageFiltered;

	for (k = 0; k < height_new; k++) {
		delete[] imageInterpolated[k];
		delete[] distance[k];
		delete[] path[k];
		delete[] potential[k];
		delete[] image_mean[k];
	}
	delete[] imageInterpolated;
	delete[] distance;
	delete[] path;
	delete[] potential;
	delete[] image_mean;

	return EXIT_SUCCESS;
}