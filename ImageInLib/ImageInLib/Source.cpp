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

	//std::string inputImagePath = inputPath + "Patient2_ct.vtk";
	std::string inputImagePath = inputPath + "aortaP2_real_origin.vtk";
	readVtkFile(inputImagePath.c_str(), ctContainer);

	int Height = ctContainer->dimensions[2];
	int Length = ctContainer->dimensions[1];
	int Width = ctContainer->dimensions[0];
	cout << "Input image dim : " << ctContainer->dimensions[0] << " x " << ctContainer->dimensions[1] << " x " << ctContainer->dimensions[2] << "" << endl;

	cout << "Image origin : (" << ctContainer->origin[0] << ", " << ctContainer->origin[1] << ", " << ctContainer->origin[2] << ")" << endl;
	Point3D ctOrigin = { ctContainer->origin[0], ctContainer->origin[1], ctContainer->origin[2] };
	VoxelSpacing ctSpacing = { ctContainer->spacing[0], ctContainer->spacing[1], ctContainer->spacing[2] };

	int i = 0, j = 0, k = 0, x = 0, x2 = 0, dim2D = Length * Width;

	//Vtk_File_Info* aortaContainer = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	//aortaContainer->operation = copyFrom;
	//inputImagePath = inputPath + "aortaP2.vtk";
	//readVtkFile(inputImagePath.c_str(), aortaContainer);
	//aortaContainer->origin[0] = ctContainer->origin[0];
	//aortaContainer->origin[1] = ctContainer->origin[1];
	//aortaContainer->origin[2] = ctContainer->origin[2];

	//aortaContainer->operation = copyTo;
	//vtkDataForm dataForm = dta_binary;
	//string outputImagePath = outputPath + "aortaP2_real_origin.vtk";
	//storeVtkFile(outputImagePath.c_str(), aortaContainer, dataForm);
	//free(aortaContainer);

	////Find min max data in order to rescalle data range positive
	//dataType min_ct = 100000000.0, max_ct = -100000000.0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length * Width; i++) {
	//		if (ctContainer->dataPointer[k][i] > max_ct) {
	//			max_ct = ctContainer->dataPointer[k][i];
	//		}
	//		if (ctContainer->dataPointer[k][i] < min_ct) {
	//			min_ct = ctContainer->dataPointer[k][i];
	//		}
	//	}
	//}
	//cout << "min data = " << min_ct << ", max data = " << max_ct << endl;

	//if (min_ct < 0) {
	//	//Rescale to positive range
	//	min_ct = (-1) * min_ct;
	//	for (k = 0; k < Height; k++) {
	//		for (i = 0; i < dim2D; i++) {
	//			ctContainer->dataPointer[k][i] = ctContainer->dataPointer[k][i] + min_ct;
	//		}
	//	}
	//}
	//rescaleNewRange(ctContainer->dataPointer, Length, Width, Height, 0.0, 1.0, max_ct, min_ct);

	string outputImagePath = outputPath + "aorta_p2.raw";
	manageRAWFile3D<dataType>(ctContainer->dataPointer, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);
	
	/*================  Data container for fast marching and path finding ============*/

	//dataType ** imageData = new dataType * [Height];
	//dataType** imageFiltered = new dataType * [Height];
	//for (k = 0; k < Height; k++) {
	//	imageData[k] = new dataType[dim2D]{0};
	//	imageFiltered[k] = new dataType[dim2D]{0};
	//	if (imageData[k] == NULL || imageFiltered[k] == NULL)
	//		return false;
	//}
	//if (imageData == NULL || imageFiltered == NULL)
	//	return false;

	//inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/raw/filteredGMC_p2.raw";
	//manageRAWFile3D<dataType>(imageFiltered, Length, Width, Height, inputImagePath.c_str(), LOAD_DATA, false);

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

	//Vtk_File_Info* petContainer = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	//petContainer->operation = copyFrom;

	//inputImagePath = inputPath + "Patient3_pet.vtk";
	//readVtkFile(inputImagePath.c_str(), petContainer);

	//int height = petContainer->dimensions[2];
	//int length = petContainer->dimensions[1];
	//int width = petContainer->dimensions[0];

	////Find min max data in order to rescalle data range positive
	//dataType min_pet = 100000000.0, max_pet = -100000000.0;
	//for (k = 0; k < height; k++) {
	//	for (i = 0; i < length * width; i++) {
	//		if (petContainer->dataPointer[k][i] > max_pet) {
	//			max_pet = petContainer->dataPointer[k][i];
	//		}
	//		if (petContainer->dataPointer[k][i] < min_pet) {
	//			min_pet = petContainer->dataPointer[k][i];
	//		}
	//	}
	//}

	//if (min_pet < 0) {
	//	//Rescale to positive range
	//	min_pet = (-1) * min_pet;
	//	for (k = 0; k < height; k++) {
	//		for (i = 0; i < length * width; i++) {
	//			petContainer->dataPointer[k][i] = petContainer->dataPointer[k][i] + min_pet;
	//		}
	//	}
	//}
	//rescaleToIntervalZeroOne(petContainer->dataPointer, length, width, height);

	//Image_Data CT; CT.imageDataPtr = imageFiltered;
	//CT.height = Height; CT.length = Length; CT.width = Width;
	//CT.orientation = orientation; CT.origin = ctOrigin; CT.spacing = ctSpacing;

	//Point3D petOrigin = { petContainer->origin[0], petContainer->origin[1], petContainer->origin[2] };
	//VoxelSpacing petSpacing = { petContainer->spacing[0], petContainer->spacing[1], petContainer->spacing[2] };

	//Image_Data PET; PET.imageDataPtr = petContainer->dataPointer; 
	//PET.height = height; PET.length = length; PET.width = width;
	//PET.orientation = orientation; PET.origin = petOrigin; PET.spacing = petSpacing;

	//const size_t height_new = (size_t)((ctContainer->dimensions[2] * ctContainer->spacing[2]) / ctContainer->spacing[0]);
	//cout << "Interpolated height : " << height_new << endl;
	//cout << "spacing : " << ctContainer->spacing[0] << endl;

	//dataType** imageInterpolated = new dataType * [height_new];
	//dataType** distance = new dataType * [height_new];
	//dataType** path = new dataType * [height_new];
	//dataType** potential = new dataType * [height_new];
	//dataType** image_mean = new dataType * [height_new];
	//for (k = 0; k < height_new; k++) {
	//	imageInterpolated[k] = new dataType[dim2D]{0};
	//	distance[k] = new dataType[dim2D]{0};
	//	path[k] = new dataType[dim2D]{0};
	//	potential[k] = new dataType[dim2D]{0};
	//	image_mean[k] = new dataType[dim2D]{0};
	//	if (imageInterpolated[k] == NULL || distance[k] == NULL || path[k] == NULL || potential[k] == NULL || image_mean[k] == NULL)
	//		return false;
	//}
	//if (imageInterpolated == NULL || distance == NULL || path == NULL || potential == NULL || image_mean == NULL)
	//	return false;

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

	//Image_Data Interpolated; Interpolated.imageDataPtr = imageInterpolated;
	//Interpolated.height = height_new; Interpolated.length = Length; Interpolated.width = Width;
	//Interpolated.orientation = orientation; Interpolated.origin = ctOrigin; 
	//Interpolated.spacing.sx = ctSpacing.sx; Interpolated.spacing.sy = ctSpacing.sy; Interpolated.spacing.sz = ctSpacing.sx;
	//imageInterpolation3D(CT, Interpolated, NEAREST_NEIGHBOR);

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

	//inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/interpolated/patient2/Im_mean_r5.raw";
	//manageRAWFile3D<dataType>(image_mean, Length, Width, height_new, inputImagePath.c_str(), LOAD_DATA, false);

	////patient 2
	//Point3D final_point = { 262, 254, 250 }; // top 
	//Point3D initial_point = { 263, 257, 146 }; // bottom
	////Point3D initial_point = { 274, 289, 452 }; // problematic point
	////Point3D middle_point = { 285, 284, 278 }; // middle

	//////Patient3
	////Point3D initial_point = { 259, 255, 244 };
	//////Point3D initial_point = { 262, 251, 261 }; // new test
	////Point3D final_point = { 260, 248, 348 };
	////Point3D middle_point = { 285, 301, 372 };

	//////Patient4
	//////Point3D middle_point = { 292, 277, 240 };
	////Point3D final_point = { 255, 224, 222 };
	////Point3D initial_point = { 269, 231, 111 };

	//////Patient5
	//////Point3D middle_point = { 260, 307, 247 };
	////Point3D final_point = { 273, 229, 230 };
	////Point3D initial_point = { 281, 229, 132 };

	//////Patient6
	//////Point3D middle_point = { 285, 266, 662 };
	//////Point3D final_point = { 243, 221, 637 }; // test 1 
	////Point3D final_point = { 241, 214, 627 }; // test 2
	////Point3D initial_point = { 265, 244, 486 };// test 1

	//////Patient7
	//////Point3D middle_point = { 277, 346, 476 };
	////Point3D final_point = { 252, 288, 442 };
	////Point3D initial_point = { 254, 299, 261 };

	////cout << "initial point : (" << initial_point.x << ", " << initial_point.y << ", " << initial_point.z << ")" << endl;

	////Get corresponding point in interpolated image
	//final_point = getRealCoordFromImageCoord3D(final_point, ctOrigin, ctSpacing, orientation);
	//final_point = getImageCoordFromRealCoord3D(final_point, ctOrigin, Interpolated.spacing, orientation);
	//final_point.x = (size_t)final_point.x; final_point.y = (size_t)final_point.y; final_point.z = (size_t)final_point.z;

	////Get corresponding point in interpolated image
	//initial_point = getRealCoordFromImageCoord3D(initial_point, ctOrigin, ctSpacing, orientation);
	//initial_point = getImageCoordFromRealCoord3D(initial_point, ctOrigin, Interpolated.spacing, orientation);
	//initial_point.x = (size_t)initial_point.x; initial_point.y = (size_t)initial_point.y; initial_point.z = (size_t)initial_point.z;

	////cout << "initial point : (" << initial_point.x << ", " << initial_point.y << ", " << initial_point.z << ")" << endl;
	//
	//Point3D* seedPoints = new Point3D[2]; 
	//seedPoints[0] = initial_point;
	//seedPoints[1] = final_point;

	//Potential_Parameters parameters;
	//parameters.K = 0.005; parameters.epsilon = 0.01;
	//parameters.c_ct = 0.0; parameters.c_pet = 0.0;
	//parameters.c_max = 0.0; parameters.c_min = 0.0; parameters.c_mean = 1.0; parameters.c_sd = 0.0;

	////cout << "initial point : (" << seedPoints[0].x << ", " << seedPoints[0].y << ", " << seedPoints[0].z << ")" << endl;
	//findPathBetweenTwoGivenPoints(Interpolated, image_mean, path, seedPoints, parameters);
	////findPathFromOneGivenPoint(Interpolated, image_mean, path, seedPoints, parameters);
	//string outputImagePath = outputPath + "path.raw";
	//manageRAWFile3D<dataType>(path, Length, Width, height_new, outputImagePath.c_str(), STORE_DATA, false);

	//computePotentialNew(Interpolated, image_mean, potential, seedPoints, 3.0, parameters);
	//computePotential_N(Interpolated, PET, statsImage, potential, seedPoints, 3.0, parameters);
	//fastMarching3D_N(imageInterpolated, distance, potential, Length, Width, height_new, seedPoints);
	//outputImagePath = outputPath + "distance_from_initial_point.raw";
	//manageRAWFile3D<dataType>(distance, Length, Width, height_new, outputImagePath.c_str(), STORE_DATA, false);

	//========================= Manual cropping ===========================

	//Image_Data cropped;

	////croppVolume(Interpolated, cropped, seedPoints, 10);
	////string outputImagePath = outputPath + "cropped.raw";
	////manageRAWFile3D<dataType>(cropped.imageDataPtr, cropped.length, cropped.width, cropped.height, outputImagePath.c_str(), STORE_DATA, false);
	////for (k = 0; k < cropped.height; k++) {
	////	delete[] cropped.imageDataPtr[k];
	////}
	////delete[] cropped.imageDataPtr;

	////const size_t offset = 10;
	////size_t i_min = (size_t)(seedPoints[1].y - 6 * offset);
	////size_t j_min = (size_t)(seedPoints[1].x - 6 * offset);
	////size_t k_min = (size_t)(seedPoints[0].z - 10 * offset);

	////size_t i_max = (size_t)(seedPoints[0].y + 6 * offset);
	////size_t j_max = (size_t)(seedPoints[0].x + 6 * offset);
	////size_t k_max = (size_t)(seedPoints[1].z + 2 * offset);

	//size_t i_min = 200, i_max = 350, j_min = 200, j_max = 350, k_min = 190, k_max = 630;
	//
	//size_t height = (size_t)(k_max - k_min);
	//size_t length = (size_t)(i_max - i_min);
	//size_t width = (size_t)(j_max - j_min);

	//cropped.height = height; //(size_t)(k_max - k_min + 1);
	//cropped.length = length; //(size_t)(i_max - i_min + 1);
	//cropped.width = width; //(size_t)(j_max - j_min + 1);

	//cout << "cropped height : " << height << endl;
	//cout << "cropped lenght : " << length << endl;
	//cout << "cropped width : " << width << endl;

	//cropped.spacing = Interpolated.spacing;
	//cropped.orientation = Interpolated.orientation;

	//cropped.origin.x = i_min; cropped.origin.y = j_min; cropped.origin.z = k_min;
	//cropped.origin = getRealCoordFromImageCoord3D(cropped.origin, Interpolated.origin, Interpolated.spacing, Interpolated.orientation);
	//cout << "origin cropped image : (" << cropped.origin.x << ", " << cropped.origin.y << ", " << cropped.origin.z << ")" << endl;

	//dataType** croppedImage = new dataType * [height];
	//dataType** croppedPath = new dataType * [height];
	//for (k = 0; k < height; k++) {
	//	croppedImage[k] = new dataType[length * width]{0};
	//	croppedPath[k] = new dataType[length * width]{0};
	//}

	//size_t k2, i2, j2;
	//for (k = 0, k2 = k_min; k < height; k++, k2++) {
	//	for (i = 0, i2 = i_min; i < length; i++, i2++) {
	//		for (j = 0, j2 = j_min; j < width; j++, j2++) {
	//			x = x_new(i, j, length);
	//			x2 = x_new(i2, j2, Length);
	//			croppedImage[k][x] = Interpolated.imageDataPtr[k2][x2]; // Interpolated.imageDataPtr[k2][x2];
	//		}
	//	}
	//}
	
	//string outputImagePath = outputPath + "cropped.raw";
	//manageRAWFile3D<dataType>(croppedImage, cropped.length, cropped.width, cropped.height, outputImagePath.c_str(), STORE_DATA, false);

	//outputImagePath = outputPath + "croppedPath.raw";
	//manageRAWFile3D<dataType>(croppedPath, cropped.length, cropped.width, cropped.height, outputImagePath.c_str(), STORE_DATA, false);

	////create initial segment
	//size_t l, m, n;
	//for (k = 0; k < height; k++) {
	//	for (i = 0; i < length; i++) {
	//		for (j = 0; j < width; j++) {
	//			x = x_new(i, j, length);
	//			if (croppedPath[k][x] == 1.0) {
	//				Point3D current_point = { i, j, k };
	//				for (l = 0; l < height; l++) {
	//					for (m = 0; m < length; m++) {
	//						for (n = 0; n < width; n++) {
	//							Point3D temporary_point = { m, n, l };
	//							double dist = getPoint3DDistance(current_point, temporary_point);
	//							if (dist <= 3.0) {
	//								initialSeg[l][x_new(m, n, length)] = 1.0;
	//							}
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}
	//}

	//outputImagePath = outputPath + "initial_segment.raw";
	//manageRAWFile3D<dataType>(initialSeg, cropped.length, cropped.width, cropped.height, outputImagePath.c_str(), STORE_DATA, false);
	//manageRAWFile3D<dataType>(initialSeg, length, width, height, outputImagePath.c_str(), LOAD_DATA, false);

	//string outputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/initial_segment.raw";
	//manageRAWFile3D<dataType>(croppedPath, length, width, height, outputImagePath.c_str(), LOAD_DATA, false);

	//outputImagePath = outputPath + "path.raw";
	////manageRAWFile3D<dataType>(croppedPath, length, width, height, outputImagePath.c_str(), STORE_DATA, false);

	//outputImagePath = outputPath + "croppedImage.raw";
	////manageRAWFile3D<dataType>(croppedImage, length, width, height, outputImagePath.c_str(), STORE_DATA, false);

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

	//delete[] seedPoints;

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
	
	//for (k = 0; k < Height; k++) {
	//	delete[] imageData[k];
	//	delete[] imageFiltered[k];
	//}
	//delete[] imageData;
	//delete[] imageFiltered;

	//for (k = 0; k < height_new; k++) {
	//	delete[] imageInterpolated[k];
	//	delete[] distance[k];
	//	delete[] path[k];
	//	delete[] potential[k];
	//	delete[] image_mean[k];
	//}
	//delete[] imageInterpolated;
	//delete[] distance;
	//delete[] path;
	//delete[] potential;
	//delete[] image_mean;

	return EXIT_SUCCESS;
}