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
	////std::string inputImagePath = inputPath + "LiverP3.vtk";
	//readVtkFile(inputImagePath.c_str(), ctContainer);

	//int Height = ctContainer->dimensions[2];
	//int Length = ctContainer->dimensions[1];
	//int Width = ctContainer->dimensions[0];
	//int dim2D = Length * Width;
	//cout << "Input image dim : " << ctContainer->dimensions[0] << " x " << ctContainer->dimensions[1] << " x " << ctContainer->dimensions[2] << "" << endl;

	//cout << "Image origin : (" << ctContainer->origin[0] << ", " << ctContainer->origin[1] << ", " << ctContainer->origin[2] << ")" << endl;
	//Point3D ctOrigin = { ctContainer->origin[0], ctContainer->origin[1], ctContainer->origin[2] };
	//VoxelSpacing ctSpacing = { ctContainer->spacing[0], ctContainer->spacing[1], ctContainer->spacing[2] };

	//int i = 0, j = 0, k = 0, x = 0, x2 = 0;

	//inputImagePath = inputPath + "Aorta_P3.vtk";
	//readVtkFile(inputImagePath.c_str(), ctContainer);
	//ctContainer->origin[0] = ctOrigin.x; 
	//ctContainer->origin[1] = ctOrigin.y;
	//ctContainer->origin[2] = ctOrigin.z;
	//vtkDataForm dataForm = dta_binary;
	//ctContainer->operation = copyTo;
	//inputImagePath = inputPath + "aorta_P3_real_origin.vtk";
	//storeVtkFile(inputImagePath.c_str(), ctContainer, dataForm);
	//exit(0);

	/*================  Data container for fast marching and path finding ============*/

	//dataType** imageFiltered = new dataType * [Height];
	//for (k = 0; k < Height; k++) {
	//	imageFiltered[k] = new dataType[dim2D]{0};
	//	if (imageFiltered[k] == NULL)
	//		return false;
	//}
	//if (imageFiltered == NULL)
	//	return false;

	//inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/raw/filteredGMC_p7.raw";
	////inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/_edge_detector_p3.raw";
	////inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/raw/patient2_ct_rescaled.raw";
	//manageRAWFile3D<dataType>(imageFiltered, Length, Width, Height, inputImagePath.c_str(), LOAD_DATA, false);

	//Image_Data CT; CT.imageDataPtr = imageFiltered;
	//CT.height = Height; CT.length = Length; CT.width = Width;
	//CT.orientation = orientation; CT.origin = ctOrigin; CT.spacing = ctSpacing;

	//FILE* file;
	//inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/points_path_plus_curvature.csv";
	//if (fopen_s(&file, inputImagePath.c_str(), "r") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}

	//const size_t datasize = 5 * sizeof(double);
	////const size_t number_of_element = datasize;
	//char first_line[datasize];
	//fgets(first_line, datasize, file);
	//char* token;
	////token = strtok_s(first_line, ",");
	//while (token != NULL) {
	//	printf("Token : %s", token);
	//	//token = strtok_s(NULL, ",");
	//}
	////printf("first line : %s", first_line);
	
	//while (feof(file) != true) {
	//	count++;
	//	fgets(first_line, datasize, file);
	//	if (count == 1) {
	//		printf("first line : %s", first_line);
	//	}
	//}

	//fclose(file);

	/*================ Find the slice with the largest piece of the liver ==============*/

	////size_t max_number = 0, count = 0, k_max = 0;
	////for (k = 0; k < Height; k++) {
	////	//count foreground pixels
	////	count = 0;
	////	for (i = 0; i < dim2D; i++) {
	////		if (ctContainer->dataPointer[k][i] == 1.0) {
	////			count++;
	////		}
	////	}
	////	//find the slice with maximum forground pixels
	////	if (count > max_number) {
	////		max_number = count;
	////		k_max = k;
	////	}
	////}
	////cout << "The slice with the large liver piece is : " << k_max << endl;

	////Store the largest slice
	//size_t k_max = 306;	
	//dataType* largestLiverSlice = new dataType[dim2D]{0};
	//for (i = 0; i < dim2D; i++) {
	//	largestLiverSlice[i] = imageFiltered[k_max][i];
	//}
	//string outputImagePath = outputPath + "largest_slice_edge_detector.raw";
	//store2dRawData<dataType>(largestLiverSlice, Length, Width, outputImagePath.c_str());
	//
	//delete[] largestLiverSlice;

	/*================== Image Filtering =============================================*/

	//Image_Data ImageData; ImageData.imageDataPtr = imageFiltered; ImageData.height = Height;
	//ImageData.length = Length; ImageData.width = Width;
	//Filter_Parameters filter_parameters; filter_parameters.edge_detector_coefficient = 1;
	//filter_parameters.h = 1.0; filter_parameters.timeStepSize = 1.2; filter_parameters.eps2 = 1e-6;
	//filter_parameters.omega_c = 1.5; filter_parameters.tolerance = 1e-2; filter_parameters.maxNumberOfSolverIteration = 100;
	//filter_parameters.timeStepsNum = 1; filter_parameters.coef = 1e-6; 
	//double firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
	//filterImage(ImageData, filter_parameters, GEODESIC_MEAN_CURVATURE_FILTER);
	//double secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
	//cout << "Time elapsed : " << (secondCpuTime - firstCpuTime) / 60 << " minutes" << endl;
	//string outputImagePath = outputPath + "filteredGMC_V3.raw";
	//manageRAWFile3D<dataType>(imageFiltered, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	//=================== Interpolation ====================================

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

	//Image_Data Interpolated; Interpolated.imageDataPtr = imageInterpolated;
	//Interpolated.height = height_new; Interpolated.length = Length; Interpolated.width = Width;
	//Interpolated.orientation = orientation; Interpolated.origin = ctOrigin;
	//Interpolated.spacing.sx = ctSpacing.sx; Interpolated.spacing.sy = ctSpacing.sy; Interpolated.spacing.sz = ctSpacing.sx;
	//imageInterpolation3D(CT, Interpolated, NEAREST_NEIGHBOR);

	////string outputImagePath = outputPath + "interpolated.raw";
	////manageRAWFile3D<dataType>(imageInterpolated, Length, Width, height_new, outputImagePath.c_str(), STORE_DATA, false);

	//////inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/initial_seg.raw";
	//inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/path_p3.raw";
	//manageRAWFile3D<dataType>(distance, Length, Width, height_new, inputImagePath.c_str(), LOAD_DATA, false);

	////create initial segment
	//size_t l = 0, m = 0, n = 0;
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

	//inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/mask_p5.raw";
	//manageRAWFile3D<dataType>(image_mean, Length, Width, height_new, inputImagePath.c_str(), LOAD_DATA, false);

	//for (k = 0; k < height_new; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		if (image_mean[k][i] == 10000) {
	//			image_mean[k][i] = 1.0;
	//		}
	//		else {
	//			image_mean[k][i] = 0.0;
	//		}
	//	}
	//}

	//================== Downsampling =====================================

	//const size_t height = (size_t)(height_new / 4);
	//cout << "heigth down : " << height << endl;
	//const size_t length = (size_t)(Length / 4);
	//const size_t width = (size_t)(Width / 4);

	//cout << "down sampled image dim : (" << length << ", " << width << ", " << height << ")" << endl;

	//dataType** downsampled_image = new dataType * [height];
	//dataType** downsampled_mask = new dataType * [height];
	//for (k = 0; k < height; k++) {
	//	downsampled_image[k] = new dataType[length * width]{0};
	//	downsampled_mask[k] = new dataType [length * width]{0};
	//}

	//Image_Data Downsampled; 
	//Downsampled.length = length; Downsampled.width = width; Downsampled.height = height;
	//Downsampled.origin = ctOrigin; Downsampled.orientation = orientation;
	//Downsampled.spacing.sx = ctSpacing.sx * 4; Downsampled.spacing.sy = ctSpacing.sy * 4; Downsampled.spacing.sz = ctSpacing.sx * 4;
	//cout << "spacing down : " << Downsampled.spacing.sx << endl;
	//cout << "slice thickness : " << Downsampled.spacing.sz << endl;

	//Downsampled.imageDataPtr = downsampled_image;
	//imageInterpolation3D(Interpolated, Downsampled, NEAREST_NEIGHBOR);

	////outputImagePath = outputPath + "down_image.raw";
	////manageRAWFile3D<dataType>(downsampled_image, length, width, height, outputImagePath.c_str(), STORE_DATA, false);
	//
	//Interpolated.imageDataPtr = path;
	//Downsampled.imageDataPtr = downsampled_mask;
	//imageInterpolation3D(Interpolated, Downsampled, NEAREST_NEIGHBOR);

	//outputImagePath = outputPath + "down_mask.raw";
	//manageRAWFile3D<dataType>(downsampled_mask, length, width, height, outputImagePath.c_str(), STORE_DATA, false);

	//string outputImagePath = outputPath + "mask_down.raw";
	////manageRAWFile3D<dataType>(downsampled_image, length, width, height, outputImagePath.c_str(), STORE_DATA, false);

	//================== Fast Marching and Path finding ===================

	//inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/interpolated/patient7/Im_mean_r3.raw";
	//manageRAWFile3D<dataType>(image_mean, Length, Width, height_new, inputImagePath.c_str(), LOAD_DATA, false);

	//inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/path_p7.raw";
	//manageRAWFile3D<dataType>(path, Length, Width, height_new, inputImagePath.c_str(), LOAD_DATA, false);

	//////patient 2
	////Point3D final_point = { 262, 254, 250 }; // top 
	////Point3D initial_point = { 263, 257, 146 }; // bottom
	//////Point3D middle_point = { 285, 284, 278 }; // middle

	//////Patient3
	////Point3D initial_point = { 259, 255, 244 };
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

	////Patient7
	////Point3D middle_point = { 277, 346, 476 };
	//Point3D final_point = { 252, 288, 442 };
	//Point3D initial_point = { 254, 299, 261 };

	////Get corresponding point in interpolated image
	//initial_point = getRealCoordFromImageCoord3D(initial_point, ctOrigin, ctSpacing, orientation);
	//initial_point = getImageCoordFromRealCoord3D(initial_point, ctOrigin, Interpolated.spacing, orientation);
	//initial_point.x = (size_t)initial_point.x; initial_point.y = (size_t)initial_point.y; initial_point.z = (size_t)initial_point.z;

	////Get corresponding point in interpolated image
	//final_point = getRealCoordFromImageCoord3D(final_point, ctOrigin, ctSpacing, orientation);
	//final_point = getImageCoordFromRealCoord3D(final_point, ctOrigin, Interpolated.spacing, orientation);
	//final_point.x = (size_t)final_point.x; final_point.y = (size_t)final_point.y; final_point.z = (size_t)final_point.z;

	//Point3D* seedPoints = new Point3D[2]; 
	//seedPoints[0] = initial_point;
	//seedPoints[1] = final_point;

	//Potential_Parameters parameters;
	//parameters.K = 0.005; parameters.epsilon = 0.01;
	//parameters.c_ct = 1.0; parameters.c_pet = 0.0;
	//parameters.c_max = 0.0; parameters.c_min = 0.0; parameters.c_mean = 1.0; parameters.c_sd = 0.0;

	////FILE* ct;
	////if (fopen_s(&ct, "C:/Users/Konan Allaly/Documents/Tests/output/image_value_path.csv", "w") != 0) {
	////	printf("Enable to open");
	////	return false;
	////}
	////FILE* pot;
	////if (fopen_s(&pot, "C:/Users/Konan Allaly/Documents/Tests/output/potential_value_path.csv", "w") != 0) {
	////	printf("Enable to open");
	////	return false;
	////}
	//computePotentialNew(Interpolated, image_mean, potential, seedPoints, 3.0, parameters);
	//string outputImagePath = outputPath + "potential_p7.raw";
	//manageRAWFile3D<dataType>(distance, Length, Width, height_new, outputImagePath.c_str(), STORE_DATA, false);
	//for (k = 0; k < height_new; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			x = x_new(j, i, Width);
	//			if (path[k][x] == 1.0 && potential[k][x] >= 0.01 && potential[k][x] <= 0.013) {
	//				distance[k][x] = 1.0;
	//				//fprintf(ct, "%d, %d, %d, %f \n", j, i, k, Interpolated.imageDataPtr[k][x]);
	//				//fprintf(pot, "%d, %d, %d, %f \n", j, i, k, potential[k][x]);
	//			}
	//		}
	//	}
	//}
	////fclose(ct);
	////fclose(pot);

	////Point3D max_potential = { 319, 272, 597 };
	////for (k = 0; k < height_new; k++) {
	////	for (i = 0; i < Length; i++) {
	////		for (j = 0; j < Width; j++) {
	////			x = x_new(j, i, Width);
	////			distance[k][x] = path[k][x];
	////			Point3D current_point = { j, i, k };
	////			double dist = getPoint3DDistance(max_potential, current_point);
	////			if (dist <= 3.0) {
	////				distance[k][x] = 1.0;
	////			}
	////		}
	////	}
	////}
	//outputImagePath = outputPath + "path_potential.raw";
	//manageRAWFile3D<dataType>(distance, Length, Width, height_new, outputImagePath.c_str(), STORE_DATA, false);

	//findPathFromOneGivenPoint(Interpolated, image_mean, path, seedPoints, parameters);
	
	//////findPathBetweenTwoGivenPoints(Interpolated, image_mean, path, seedPoints, parameters);
	//string outputImagePath = outputPath + "path.raw";
	//manageRAWFile3D<dataType>(path, Length, Width, height_new, outputImagePath.c_str(), STORE_DATA, false);

	////////Point3D point_max = {276.6567, 237.4935, 615.6494}; // max p5
	////////Point3D point_max = { 281.7054, 277.0704, 621.4479 }; // max p4
	//////Point3D point_max = { 244.5328, 242.9775, 692.7566 }; // max p3
	////Point3D point_max = { 266.6567, 249.7204, 718.3799 }; // max p2
	////Point3D point_max = { 260.947, 309.217, 926.9254 }; // max p6
	//Point3D point_max = { 241.9134, 303.439, 689.771 }; // max p7
	//manageRAWFile3D<dataType>(path, Length, Width, height_new, outputImagePath.c_str(), LOAD_DATA, false);

	//for (k = 0; k < height_new; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			x = x_new(j, i, Width);
	//			Point3D current_point = { j, i, k };
	//			double dist = getPoint3DDistance(point_max, current_point);
	//			if (dist <= 3.0) {
	//				path[k][x] = 1.0;
	//			}
	//		}
	//	}
	//}

	//outputImagePath = outputPath + "path_with_max.raw";
	//manageRAWFile3D<dataType>(path, Length, Width, height_new, outputImagePath.c_str(), STORE_DATA, false);

	//================== Automatic Cropping =========================================

	//imageMetaData croppedImageMetaData;
	//croppedImageMetaData = croppImage3D(CT, 10);
	//cout << "Cropped image origin RW CS : (" << croppedImageMetaData.origin.x << ", " << croppedImageMetaData.origin.y << ", " << croppedImageMetaData.origin.z << ")" << endl;

	//croppedImageMetaData.origin = getImageCoordFromRealCoord3D(croppedImageMetaData.origin, ctOrigin, ctSpacing, orientation);
	//size_t i0 = (size_t)croppedImageMetaData.origin.x, j0 = (size_t)croppedImageMetaData.origin.y, k0 = (size_t)croppedImageMetaData.origin.z;
	//cout << "Cropped image origin IMG CS : (" << i0 << ", " << j0 << ", " << k0 << ")" << endl;

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

	//================== Segmentation =====================================

	//Point3D* center_seg = new Point3D[1];
	//center_seg->x = 0; center_seg->y = 0; center_seg->z = 0; // not used in the current version of the code

	////////dataType** initial_seg = new dataType * [height];
	////////for (k = 0; k < height; k++) {
	////////	initial_seg[k] = new dataType[length * width]{0};
	////////}

	////////inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/path_down.raw";
	////////manageRAWFile3D<dataType>(initial_seg, length, width, height, inputImagePath.c_str(), LOAD_DATA, false);

	////string outputImagePath = outputPath + "/segmentation/_seg_func_00.raw";
	////manageRAWFile3D<dataType>(downsampled_mask, length, width, height, outputImagePath.c_str(), STORE_DATA, false);
	//////////manageRAWFile3D<dataType>(initialSeg, length_down, width_down, height_down, outputImagePath.c_str(), STORE_DATA, false);

	//Image_Data ImageToBeSegmented;
	//ImageToBeSegmented.imageDataPtr = imageFiltered; // to be updated
	//ImageToBeSegmented.height = Height; ImageToBeSegmented.length = Length; ImageToBeSegmented.width = Width;

	//Segmentation_Parameters seg_params;
	//seg_params.tau = 1.2; seg_params.h = 1.0; seg_params.omega_c = 1.5; seg_params.mod = 1;
	//seg_params.maxNoGSIteration = 100; seg_params.coef = 100000; seg_params.gauss_seidelTolerance = 1e-6;
	//seg_params.maxNoOfTimeSteps = 1; seg_params.segTolerance = 1e-6; seg_params.eps2 = 1e-6;
	//seg_params.coef_conv = 1.0; seg_params.coef_dif = 1.0;

	////Smoothing by heat equation so we don't need all the parameters
	//Filter_Parameters parameters; parameters.coef = 1e-6;
	//parameters.h = seg_params.h; parameters.timeStepSize = parameters.h / 4.0; parameters.timeStepsNum = 1; parameters.omega_c = 1.5; 
	//parameters.tolerance = 1e-6, parameters.maxNumberOfSolverIteration = 100;
	////Filter_Parameters parameters; parameters.coef = 1e-4; parameters.edge_detector_coefficient = 1000; parameters.eps2 = 1e-4;
	////parameters.h = 6.0; parameters.maxNumberOfSolverIteration = 100; parameters.omega_c = 1.4; parameters.p = 1;
	////parameters.sigma = 0.1; parameters.timeStepSize = 1.0; parameters.timeStepsNum = 1; parameters.tolerance = 1e-6;

	//double firstCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);
	//unsigned char segmentPath[] = "C:/Users/Konan Allaly/Documents/Tests/output/segmentation/";
	///*                                                 update required                     */
	//generalizedSubsurfSegmentation(ImageToBeSegmented, imageFiltered, seg_params, parameters, center_seg, 1, segmentPath);
	//double secondCpuTime = clock() / (dataType)(CLOCKS_PER_SEC);

	//cout << "Time : " << (secondCpuTime - firstCpuTime) / 60 << " minutes" << endl;

	//========================= Save CSV =====================================

	//inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/path_p2_ct_average.raw";
	//manageRAWFile3D<dataType>(path, Length, Width, height_new, inputImagePath.c_str(), LOAD_DATA, false);

	//FILE* file;
	//if (fopen_s(&file, "C:/Users/Konan Allaly/Documents/Tests/output/points_found.csv", "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}

	//for (k = 0; k < height_new; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			x = x_new(j, i, Width);
	//			if (path[k][x] == 1.0) {
	//				fprintf(file, "%d, %d, %d \n", j, i, k);
	//			}
	//		}
	//	}
	//}
	//fclose(file);

	//========================= Compute curvature for each path point ========

	//FILE* file;
	//if (fopen_s(&file, "C:/Users/Konan Allaly/Documents/Tests/output/curvature.csv", "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}

	//vector<Point3D> list_of_path_point;
	//Point3D current_point = { 0.0, 0.0, 0.0 }, previous_point = { 0.0, 0.0, 0.0 }, next_point = {0.0, 0.0, 0.0};

	////inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/path_p2_ct_average.raw";
	//inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/path_0.raw";
	//manageRAWFile3D<dataType>(path, Length, Width, height_new, inputImagePath.c_str(), LOAD_DATA, false);

	//for (k = 0; k < height_new; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			x = x_new(j, i, Width);
	//			if (path[k][x] == 1.0) {
	//				current_point.x = j; 
	//				current_point.y = i;
	//				current_point.z = k;
	//				list_of_path_point.push_back(current_point);
	//			}
	//		}
	//	}
	//}

	//int number_of_point = list_of_path_point.size();
	//double* distance_between_point = new double[number_of_point] {0};
	//for (i = 1; i < number_of_point; i++) {
	//	distance_between_point[i] = getPoint3DDistance(list_of_path_point[i], list_of_path_point[i - 1]);
	//}

	//double coef = 0.0, x_component = 0.0, y_component = 0.0, z_component = 0.0, curvature = 0.0;
	//for (i = 1; i < number_of_point - 1; i++) {
	//	previous_point = list_of_path_point[i - 1];
	//	current_point = list_of_path_point[i];
	//	next_point = list_of_path_point[i + 1];
	//	coef = 2 / (distance_between_point[i + 1] + distance_between_point[i]);
	//	x_component = coef * (((next_point.x - current_point.x) / distance_between_point[i + 1]) - ((current_point.x - previous_point.x) / distance_between_point[i]));
	//	y_component = coef * (((next_point.y - current_point.y) / distance_between_point[i + 1]) - ((current_point.y - previous_point.y) / distance_between_point[i]));
	//	z_component = coef * (((next_point.z - current_point.z) / distance_between_point[i + 1]) - ((current_point.z - previous_point.z) / distance_between_point[i]));
	//	curvature = sqrt(x_component * x_component + y_component * y_component + z_component * z_component);
	//	fprintf(file, "%d, %lf \n", i, curvature);
	//}
	//fclose(file);
	//delete[] distance_between_point;

	//========================= Free memory ==================================

	//delete[] seedPoints;

	free(ctContainer);

	//for (k = 0; k < height; k++) {
	//	delete[] downsampled_image[k];
	//	delete[] downsampled_mask[k];
	//	//delete[] initialSeg[k];
	//}
	//delete[] downsampled_image;
	//delete[] downsampled_mask;
	////delete[] initialSeg;

	//for (k = 0; k < height; k++) {
	//	delete[] croppedImage[k];
	//	delete[] croppedPath[k];
	//}
	//delete[] croppedImage;
	//delete[] croppedPath;

	//delete[] center_seg;
	
	//for (k = 0; k < Height; k++) {
	//	delete[] imageFiltered[k];
	//}
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