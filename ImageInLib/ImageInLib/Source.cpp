// Source code for fast marching and path finding (11/09/2023)

#include <iostream>
#include <sstream>  
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
#include "segmentation2d.h"
#include "../src/thresholding.h"
#include "../src/distance_function.h"
#include "morphological_change.h"
#include "Labelling.h"
#include"hough_transform.h"

using namespace std;

int main() {

	string inputPath = "C:/Users/Konan Allaly/Documents/Tests/input/";
	string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";
	
	string saving_name, extension;

	//===================== Load 3D patient data (.vtk) ================================

	OrientationMatrix orientation;
	orientation.v1 = { 1.0, 0.0, 0.0 }; orientation.v2 = { 0.0, 1.0, 0.0 }; orientation.v3 = { 0.0, 0.0, 1.0 };

	Vtk_File_Info* ctContainer = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	ctContainer->operation = copyFrom;

	std::string inputImagePath = inputPath + "vtk/Patient2_ct.vtk";
	//std::string inputImagePath = inputPath + "vtk/Aorta_p4.vtk";
	readVtkFile(inputImagePath.c_str(), ctContainer);

	int Height = ctContainer->dimensions[2];
	int Length = ctContainer->dimensions[1];
	int Width = ctContainer->dimensions[0];
	int dim2D = Length * Width;
	std::cout << "Input image dim : " << ctContainer->dimensions[0] << " x " << ctContainer->dimensions[1] << " x " << ctContainer->dimensions[2] << "" << std::endl;

	std::cout << "Image origin : (" << ctContainer->origin[0] << ", " << ctContainer->origin[1] << ", " << ctContainer->origin[2] << ")" << std::endl;
	Point3D ctOrigin = { ctContainer->origin[0], ctContainer->origin[1], ctContainer->origin[2] };
	VoxelSpacing ctSpacing = { ctContainer->spacing[0], ctContainer->spacing[1], ctContainer->spacing[2] };
	std::cout << "Image spacing : (" << ctContainer->spacing[0] << ", " << ctContainer->spacing[1] << ", " << ctContainer->spacing[2] << ")" << std::endl;

	//define 3D data structure
	int k, i, j, xd;
	dataType** imageData = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{0};
	}

	////find minimum and maximum
	//dataType min_data = 100000.0, max_data = 0.0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		if (ctContainer->dataPointer[k][i] > max_data) {
	//			max_data = ctContainer->dataPointer[k][i];
	//		}
	//		if (ctContainer->dataPointer[k][i] < min_data) {
	//			min_data = ctContainer->dataPointer[k][i];
	//		}
	//	}
	//}
	//std::cout << "min data = " << min_data << ", max data = " << max_data << std::endl;

	////rescal to positive range
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		imageData[k][i] = ctContainer->dataPointer[k][i] + (- 1) * min_data;
	//	}
	//}

	//===================== Interpolation =================================
	
	//int height = (int)((Height * ctSpacing.sz) / ctSpacing.sx);
	//dataType** image_interp = new dataType * [height];
	//dataType** path = new dataType * [height];
	//for (k = 0; k < height; k++) {
	//	image_interp[k] = new dataType[dim2D]{ 0 };
	//	path[k] = new dataType[dim2D]{ 0 };
	//}
	//Image_Data CT; CT.imageDataPtr = imageData; CT.height = Height; CT.length = Length; CT.width = Width;
	//CT.origin = ctOrigin; CT.spacing = ctSpacing; CT.orientation = orientation;
	//
	//Image_Data interpolated; interpolated.imageDataPtr = image_interp;
	//interpolated.height = height; interpolated.length = Length; interpolated.width = Width;
	//interpolated.orientation = orientation; interpolated.origin = ctOrigin;
	//interpolated.spacing.sx = ctSpacing.sx; interpolated.spacing.sy = ctSpacing.sy; interpolated.spacing.sz = ctSpacing.sx;
	//cout << "interpolated image dim : " << ctContainer->dimensions[0] << " x " << ctContainer->dimensions[1] << " x " << height << "" << endl;
	//imageInterpolation3D(CT, interpolated, NEAREST_NEIGHBOR);

	//===================== 2D Hough Transform to detect circle on path image ========================

	/*
	inputImagePath = inputPath + "path_p2.raw";
	manageRAWFile3D<dataType>(path, Length, Width, height, inputImagePath.c_str(), LOAD_DATA, false);

	vector<Point3D> point_path;
	Point3D current_point = { 0.0, 0.0, 0.0 };
	for (k = 0; k < height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				xd = x_new(i, j, Length);
				if (path[k][xd] == 1.0) {
					current_point.x = i;
					current_point.y = j;
					current_point.z = k;
					point_path.push_back(current_point);
				}
			}
		}
	}
	
	dataType* imageSlice = new dataType[dim2D]{0};
	dataType* houghSpace = new dataType[dim2D]{0};
	dataType* voteArray = new dataType[dim2D]{0};

	Point2D seed2D = { 0.0, 0.0 };

	//filtering parameters
	Filter_Parameters filter_parameters;
	filter_parameters.h = 1.0; filter_parameters.timeStepSize = 0.25; filter_parameters.eps2 = 1e-6;
	filter_parameters.omega_c = 1.5; filter_parameters.tolerance = 1e-3; filter_parameters.maxNumberOfSolverIteration = 100;
	filter_parameters.timeStepsNum = 5; filter_parameters.coef = 1e-6;
	Image_Data2D IMAGE; IMAGE.imageDataPtr = imageSlice;
	IMAGE.height = Length; IMAGE.width = Width;

	HoughParameters parameters = { 6.0, 15.0, 0.5, 5.0, 1000, 1.0, 0.15, 0.5 };

	int l = point_path.size();
	for (i = 0; i < l; i++) {
		
		extension = to_string(i);
		k = (size_t)point_path[i].z;
		seed2D.x = point_path[i].x;
		seed2D.y = point_path[i].y;

		copyDataToAnother2dArray(image_interp[k], imageSlice, Length, Width);

		//Slice rescaling
		rescaleNewRange2D(imageSlice, Length, Width, 0.0, 1.0);
		//Slice filtering
		heatImplicit2dScheme(IMAGE, filter_parameters);

		localHoughTransform(seed2D, imageSlice, voteArray, houghSpace, Length, Width, parameters);
		if (i < 10) {
			saving_name = outputPath + "slice_00" + extension + ".raw";
			store2dRawData(imageSlice, Length, Width, saving_name.c_str());
			saving_name = outputPath + "found_circle_00" + extension + ".raw";
			store2dRawData(houghSpace, Length, Width, saving_name.c_str());
		}
		else {
			if (i < 100) {
				saving_name = outputPath + "slice_0" + extension + ".raw";
				store2dRawData(imageSlice, Length, Width, saving_name.c_str());
				saving_name = outputPath + "found_circle_0" + extension + ".raw";
				store2dRawData(houghSpace, Length, Width, saving_name.c_str());
			}
			else {
				saving_name = outputPath + "slice_" + extension + ".raw";
				store2dRawData(imageSlice, Length, Width, saving_name.c_str());
				saving_name = outputPath + "found_circle_" + extension + ".raw";
				store2dRawData(houghSpace, Length, Width, saving_name.c_str());
			}
		}
	}

	delete[] imageSlice;
	delete[] houghSpace;
	delete[] voteArray;
	*/

	//===================== Fast marching and 2d Hough transform == ===========
	
	int height = (int)((Height * ctSpacing.sz) / ctSpacing.sx);
	dataType** image_interp = new dataType * [height];
	dataType** path = new dataType * [height];
	dataType** image_mean = new dataType * [height];
	for (k = 0; k < height; k++) {
		image_interp[k] = new dataType[dim2D]{ 0 };
		path[k] = new dataType[dim2D]{ 0 };
		image_mean[k] = new dataType[dim2D]{ 0 };
	}

	Image_Data CT; CT.imageDataPtr = imageData; CT.height = Height; CT.length = Length; CT.width = Width;
	CT.origin = ctOrigin; CT.spacing = ctSpacing; CT.orientation = orientation;

	inputImagePath = inputPath + "/raw/filteredGMC_p2.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, inputImagePath.c_str(), LOAD_DATA, false);

	Image_Data interpolated; interpolated.imageDataPtr = image_interp;
	interpolated.height = height; interpolated.length = Length; interpolated.width = Width;
	interpolated.orientation = orientation; interpolated.origin = ctOrigin;
	interpolated.spacing.sx = ctSpacing.sx; interpolated.spacing.sy = ctSpacing.sy; interpolated.spacing.sz = ctSpacing.sx;
	cout << "interpolated image dim : " << ctContainer->dimensions[0] << " x " << ctContainer->dimensions[1] << " x " << height << "" << endl;
	imageInterpolation3D(CT, interpolated, NEAREST_NEIGHBOR);
	
	//patient 2
	Point3D final_point = { 262, 254, 250 }; // top 
	//Point3D initial_point = { 263, 257, 146 }; // bottom

	Point3D initial_point = { 263, 257, 146 }; // bottom

	////Patient3
	//Point3D initial_point = { 259, 255, 244 };
	//Point3D final_point = { 260, 248, 348 };

	////Patient4
	//Point3D final_point = { 255, 224, 222 };
	//Point3D initial_point = { 269, 231, 111 };

	inputImagePath = inputPath + "/interpolated/patient2/Im_mean_r5.raw";
	manageRAWFile3D<dataType>(image_mean, Length, Width, height, inputImagePath.c_str(), LOAD_DATA, false);

	////Get corresponding point in interpolated image
	//initial_point = getRealCoordFromImageCoord3D(initial_point, ctOrigin, ctSpacing, orientation);
	//initial_point = getImageCoordFromRealCoord3D(initial_point, ctOrigin, interpolated.spacing, orientation);
	//initial_point.x = (size_t)initial_point.x; initial_point.y = (size_t)initial_point.y; initial_point.z = (size_t)initial_point.z;

	//Get corresponding point in interpolated image
	final_point = getRealCoordFromImageCoord3D(final_point, ctOrigin, ctSpacing, orientation);
	final_point = getImageCoordFromRealCoord3D(final_point, ctOrigin, interpolated.spacing, orientation);
	final_point.x = (size_t)final_point.x; final_point.y = (size_t)final_point.y; final_point.z = (size_t)final_point.z;

	Point3D* seedPoints = new Point3D[2]; 
	seedPoints[0] = initial_point;
	seedPoints[1] = final_point;

	Potential_Parameters parameters;
	parameters.K = 0.005; parameters.epsilon = 0.01;
	parameters.c_ct = 0.0; parameters.c_pet = 0.0;
	parameters.c_max = 0.0; parameters.c_min = 0.0; parameters.c_mean = 1.0; parameters.c_sd = 0.0;

	//findPathFromOneGivenPointWithStopingCriterium(interpolated, image_mean, path, seedPoints, parameters);
	//findPathFromOneGivenPointWithCircleDetection(interpolated, image_mean, path, seedPoints, parameters);

	findPathFromOneGivenPoint(interpolated, image_mean, path, seedPoints, parameters);

	//===================== 2D Hough transform from path points ==========

	/*
	//const size_t Length = 512, Width = 512, dim2D = Length * Width;
	dataType* imageSlice = new dataType[dim2D]{ 0 };
	dataType* houghSpace = new dataType[dim2D]{ 0 };
	dataType* voteArray = new dataType[dim2D]{ 0 };

	Point2D seed2D = { 0.0, 0.0 };
	size_t i_min, i_max, j_min, j_max, k_min, k_max;

	//filtering parameters
	Filter_Parameters filter_parameters;
	filter_parameters.h = 1.0; filter_parameters.timeStepSize = 0.25; filter_parameters.eps2 = 1e-6;
	filter_parameters.omega_c = 1.5; filter_parameters.tolerance = 1e-3; filter_parameters.maxNumberOfSolverIteration = 100;
	filter_parameters.timeStepsNum = 5; filter_parameters.coef = 1e-6;
	Image_Data2D IMAGE; IMAGE.imageDataPtr = imageSlice;
	IMAGE.height = Length; IMAGE.width = Width;

	////Load path point from file .csv
	FILE* file_load;
	inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/path_points.csv";
	if (fopen_s(&file_load, inputImagePath.c_str(), "r") != 0) {
		printf("Enable to open");
		return false;
	}
	int length_array = 708;

	FILE* file_store;
	inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/output/path_points_ratio.csv";
	if (fopen_s(&file_store, inputImagePath.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}

	dataType max_ratio = 0.0;
	dataType coordX, coordY, coordZ;
	HoughParameters parameters = { 4.0, 15.0, 0.5, 5.0, 1000, 1.0, 0.15, 0.5 };
	for (i = 0; i < length_array; i++) {
		
		fscanf_s(file_load, "%f", &(coordX));
		fscanf_s(file_load, ",");
		fscanf_s(file_load, "%f", &(coordY));
		fscanf_s(file_load, ",");
		fscanf_s(file_load, "%f", &(coordZ));
		fscanf_s(file_load, "\n");
		//path[(size_t)coordZ][x_new((size_t)coordX, (size_t)coordY, Length)] = 1.0;
		
		//set seed
		seed2D.x = coordX;
		seed2D.y = coordY;
		k = (size_t)coordZ;

		image_interp[(size_t)coordZ][x_new((size_t)coordX, (size_t)coordY, Length)] = 1.0;

		if (i >= 330 && i <= 535) {
			Point3D current = { coordX, coordY, coordZ };
			i_min = (size_t)(coordX - 10); i_max = (size_t)(coordX + 10);
			j_min = (size_t)(coordY - 10); j_max = (size_t)(coordY + 10);
			k_min = (size_t)(coordZ - 10); k_max = (size_t)(coordZ + 10);
			image_interp[(size_t)coordZ][x_new((size_t)coordX, (size_t)coordY, Length)] = 1.0;
			for (size_t t = k_min; t <= k_max; t++) {
				for (size_t u = i_min; u <= i_max; u++) {
					for (size_t v = j_min; v <= j_max; v++) {
						size_t tvu = x_new(u, v, Length);
						Point3D point = { (dataType)u, (dataType)v, (dataType)t };
						double distance = getPoint3DDistance(current, point);
						if (distance <= 3) {
							image_interp[t][tvu] = 1.0;
						}
					}
				}
			}
		}

		////extract slice
		//copyDataToAnother2dArray(interpolated.imageDataPtr[k], imageSlice, Length, Width);

		////Slice rescaling
		//rescaleNewRange2D(imageSlice, Length, Width, 0.0, 1.0);

		////Slice filtering
		//heatImplicit2dScheme(IMAGE, filter_parameters);
		
		//extension = to_string(i);
		//if (i < 10) {
		//	saving_name = outputPath + "threshold_000" + extension + ".raw";
		//}
		//else {
		//	if (i < 100) {
		//		saving_name = outputPath + "threshold_00" + extension + ".raw";
		//	}
		//	else {
		//		if (i < 1000) {
		//			saving_name = outputPath + "threshold_0" + extension + ".raw";
		//		}
		//		else {
		//			saving_name = outputPath + "threshold_" + extension + ".raw";
		//		}
		//	}
		//}

		//localHoughTransform(seed2D, imageSlice, voteArray, houghSpace, Length, Width, parameters, saving_name);

		//if (i < 10) {
		//	saving_name = outputPath + "slice_000" + extension + ".raw";
		//	store2dRawData(imageSlice, Length, Width, saving_name.c_str());
		//	saving_name = outputPath + "found_circle_000" + extension + ".raw";
		//	store2dRawData(houghSpace, Length, Width, saving_name.c_str());
		//	saving_name = outputPath + "ratio_000" + extension + ".raw";
		//	store2dRawData(voteArray, Length, Width, saving_name.c_str());
		//}
		//else {
		//	if (i < 100) {
		//		saving_name = outputPath + "slice_00" + extension + ".raw";
		//		store2dRawData(imageSlice, Length, Width, saving_name.c_str());
		//		saving_name = outputPath + "found_circle_00" + extension + ".raw";
		//		store2dRawData(houghSpace, Length, Width, saving_name.c_str());
		//		saving_name = outputPath + "ratio_00" + extension + ".raw";
		//		store2dRawData(voteArray, Length, Width, saving_name.c_str());
		//	}
		//	else {
		//		if (i < 1000) {
		//			saving_name = outputPath + "slice_0" + extension + ".raw";
		//			store2dRawData(imageSlice, Length, Width, saving_name.c_str());
		//			saving_name = outputPath + "found_circle_0" + extension + ".raw";
		//			store2dRawData(houghSpace, Length, Width, saving_name.c_str());
		//			saving_name = outputPath + "ratio_0" + extension + ".raw";
		//			store2dRawData(voteArray, Length, Width, saving_name.c_str());
		//		}
		//		else {
		//			saving_name = outputPath + "slice_" + extension + ".raw";
		//			store2dRawData(imageSlice, Length, Width, saving_name.c_str());
		//			saving_name = outputPath + "found_circle_" + extension + ".raw";
		//			store2dRawData(houghSpace, Length, Width, saving_name.c_str());
		//			saving_name = outputPath + "ratio_" + extension + ".raw";
		//			store2dRawData(voteArray, Length, Width, saving_name.c_str());
		//		}
		//	}
		//}

		//max_ratio = getTheMaxValue(voteArray, Length, Width);
		//fprintf(file_store, "%f,%f,%f,%f\n", coordX, coordY, coordZ, max_ratio);
	}
	fclose(file_load);
	fclose(file_store);

	saving_name = outputPath + "updated_path.raw";
	store3dRawData<dataType>(image_interp, Length, Width, height, saving_name.c_str());

	delete[] imageSlice;
	delete[] houghSpace;
	delete[] voteArray;
	*/

	//===================== Hough transform on 2D Slice ==================

	/*
	const size_t Length = 512, Width = 512, dim2D = Length * Width;

	dataType* imageSlice = new dataType[dim2D]{ 0 };
	dataType* houghSpace = new dataType[dim2D]{ 0 };
	dataType* voteArray = new dataType[dim2D]{ 0 };

	//slice 7 --> (263, 256, 317)
	//slice 38 --> (262, 255, 341)
	//Point2D seed2D = { 263.0, 256.0 };
	//Point2D seed2D = { 262.0, 255.0 }; // slice_0038
	//Point2D seed2D = { 264.0, 266.0 }; // slice_0115
	//Point2D seed2D = { 264.0, 267.0 }; // slice_0120
	Point2D seed2D = { 276.0, 253.0 }; // slice_0560
	//Point2D seed2D = { 283.0, 252.0 }; // slice_0570

	//load Slice
	//string inputImagePath = inputPath + "slice_0038.raw";
	//string inputImagePath = inputPath + "slice_0115.raw";
	//string inputImagePath = inputPath + "slice_0038.raw";
	string inputImagePath = inputPath + "slice_0560.raw";
	load2dArrayRAW<dataType>(imageSlice, Length, Width, inputImagePath.c_str(), false);

	HoughParameters parameters = { 7.0, 15.0, 0.5, 5.0, 1000, 1.0, 0.15, 0.5 };

	saving_name = outputPath + "threshold.raw";
	localHoughTransform(seed2D, imageSlice, voteArray, houghSpace, Length, Width, parameters, saving_name);
	
	saving_name = outputPath + "found_circle_without.raw";
	store2dRawData(houghSpace, Length, Width, saving_name.c_str());
	saving_name = outputPath + "ratio_without.raw";
	store2dRawData(voteArray, Length, Width, saving_name.c_str());

	delete[] imageSlice;
	delete[] houghSpace;
	delete[] voteArray;
	*/

	//=================== Test 2D erosion =============================================
	
	/*
	const size_t Length = 512, Width = 512, dim2D = Length * Width;
	size_t i;

	dataType* imageSlice = new dataType[dim2D]{ 0 };
	dataType* edgeDetector = new dataType[dim2D]{ 0 };
	dataType* gradientX = new dataType[dim2D]{ 0 };
	dataType* gradientY = new dataType[dim2D]{ 0 };

	//load Slice
	string inputImagePath = inputPath + "slice_0038.raw";
	load2dArrayRAW<dataType>(imageSlice, Length, Width, inputImagePath.c_str(), false);

	saving_name = outputPath + "loaded.raw";
	store2dRawData(imageSlice, Length, Width, saving_name.c_str());

	//Compute edge detector
	dataType norm_of_gradient = 0.0, K = 1000.0;
	computeImageGradient(imageSlice, gradientX, gradientY, Length, Width, 1.0);
	for (i = 0; i < dim2D; i++) {
		norm_of_gradient = sqrt(gradientX[i] * gradientX[i] + gradientY[i] * gradientY[i]);
		edgeDetector[i] = gradientFunction(norm_of_gradient, K);
	}

	dataType thres = 0.15;
	thresholding2DFunction(edgeDetector, Length, Width, thres, thres);

	saving_name = outputPath + "threshold.raw";
	store2dRawData(edgeDetector, Length, Width, saving_name.c_str());

	erosion2D(edgeDetector, Length, Width, 1.0, 0.0);
	erosion2D(edgeDetector, Length, Width, 1.0, 0.0);
	saving_name = outputPath + "eroded_two_times.raw";
	store2dRawData(edgeDetector, Length, Width, saving_name.c_str());
	*/

	//=================== Free memory =================================================

	delete[] seedPoints;
	for (k = 0; k < height; k++) {
		delete[] image_interp[k];
		delete[] path[k];
		delete[] image_mean[k];
	}
	delete[] image_interp;
	delete[] path;
	delete[] image_mean;

	free(ctContainer);
	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
	}
	delete[] imageData;

	return EXIT_SUCCESS;
}