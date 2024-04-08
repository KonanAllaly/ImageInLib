// Source code for localization of the heart (06/02/2024)

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
#include"../src/edgedetection.h"

#define thres_min 995
#define thres_max 1213
#define thres -500

int main() {

	string inputPath = "C:/Users/Konan Allaly/Documents/Tests/input/";
	string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";
	
	string loading_path, storing_path, extension;

	size_t i, j, k, xd;

	//===================== Load 3D patient data (.vtk) ================================
	
	OrientationMatrix orientation;
	orientation.v1 = { 1.0, 0.0, 0.0 }; 
	orientation.v2 = { 0.0, 1.0, 0.0 }; 
	orientation.v3 = { 0.0, 0.0, 1.0 };

	
	//aorta
	//loading_path = inputPath + "vtk/aorta/aortaP2_real_origin.vtk";
	//loading_path = inputPath + "vtk/aorta/aorta_P3_real_origin.vtk";
	//loading_path = inputPath + "vtk/lungs/cropped_ct_lungs_p6.vtk";
	//loading_path = inputPath + "vtk/lungs/box_lungs_ct_p6.vtk";
	//loading_path = inputPath + "vtk/ct/Patient6_ct.vtk";
	//loading_path = inputPath + "vtk/pet/Patient3_pet.vtk";
	//loading_path = inputPath + "patient/after/patient_ct.vtk";

	Vtk_File_Info* ctContainer = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	ctContainer->operation = copyFrom;

	loading_path = inputPath + "vtk/ct/Patient2_ct.vtk"; // patient before treatment
	//loading_path = inputPath + "vtk/ct/Patient3_ct.vtk"; // patient after treatment
	readVtkFile(loading_path.c_str(), ctContainer);

	int Height = ctContainer->dimensions[2];
	int Length = ctContainer->dimensions[0];
	int Width = ctContainer->dimensions[1];
	int dim2D = Length * Width;
	std::cout << "CT image dim : " << ctContainer->dimensions[0] << " x " << ctContainer->dimensions[1] << " x " << ctContainer->dimensions[2] << "" << std::endl;

	std::cout << "CT origin : (" << ctContainer->origin[0] << ", " << ctContainer->origin[1] << ", " << ctContainer->origin[2] << ")" << std::endl;
	Point3D ctOrigin = { ctContainer->origin[0], ctContainer->origin[1], ctContainer->origin[2] };
	VoxelSpacing ctSpacing = { ctContainer->spacing[0], ctContainer->spacing[1], ctContainer->spacing[2] };
	std::cout << "CT spacing : (" << ctContainer->spacing[0] << ", " << ctContainer->spacing[1] << ", " << ctContainer->spacing[2] << ")" << std::endl;

	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		ctContainer->dataPointer[k][i] = 0.0;
	//	}
	//}

	//define new array to load input images
	dataType** imageData = new dataType*[Height];
	dataType** maskThreshold = new dataType*[Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{0};
		maskThreshold[k] = new dataType[dim2D]{0};
		if (imageData[k] == NULL || maskThreshold[k] == NULL)
			return false;
	}
	if (imageData == NULL || maskThreshold == NULL)
		return false;

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
	////cout << "min data = " << min_data << ", max data = " << max_data << endl;

	////rescale to positive range
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		imageData[k][i] = ctContainer->dataPointer[k][i] + (-1) * min_data;
	//		//maskThreshold[k][i] = ctContainer->dataPointer[k][i] + (-1) * min_data;
	//		//maskThreshold[k][i] = ctContainer->dataPointer[k][i];
	//	}
	//}
	//max_data = max_data + (-1) * min_data;
	//min_data = (-1) * min_data + min_data;
	
	//cout << "new min data = " << min_data << ", new max data = " << max_data << endl;
	//rescaleNewRange(imageData, Length, Width, Height, 0.0, 1.0, max_data, min_data);

	Image_Data CT; 
	CT.imageDataPtr = imageData;
	CT.height = Height; 
	CT.length = Length; 
	CT.width = Width;
	CT.orientation = orientation; 
	CT.origin = ctOrigin; 
	CT.spacing = ctSpacing;

	size_t height = (size_t)((ctSpacing.sz / ctSpacing.sx) * Height);
	cout << "interpol height : " << height << endl;

	dataType** interpolated = new dataType * [height];
	for (k = 0; k < height; k++) {
		interpolated[k] = new dataType[dim2D]{ 0 };
	}

	//data structure for interpolated
	Image_Data INTERPOLATE;
	INTERPOLATE.imageDataPtr = interpolated;
	INTERPOLATE.height = height;
	INTERPOLATE.length = Length;
	INTERPOLATE.width = Width;
	INTERPOLATE.orientation = orientation;
	INTERPOLATE.origin = ctOrigin;
	INTERPOLATE.spacing.sx = ctSpacing.sx;
	INTERPOLATE.spacing.sy = ctSpacing.sy;
	INTERPOLATE.spacing.sz = ctSpacing.sx;

	loading_path = inputPath + "raw/interpolated/patient2/filtered_p2.raw";
	load3dArrayRAW<dataType>(interpolated, Length, Width, height, loading_path.c_str(), false);
	

	////======================= Update dimension ==========================

	/*
	Vtk_File_Info* newDimContainer = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	newDimContainer->operation = copyFrom;
	//loading_path = inputPath + "vasculitis/patient_before_treatment/manual/cropped/liver_before_m.vtk";
	//loading_path = inputPath + "vasculitis/patient_before_treatment/manual/cropped/seg5_before.vtk";
	//loading_path = inputPath + "vasculitis/patient_before_treatment/manual/cropped/ball_liver_before.vtk";
	loading_path = inputPath + "vasculitis/patient_before_treatment/manual/cropped/aorta_before.vtk";
	readVtkFile(loading_path.c_str(), newDimContainer);
	Point3D newOrigin = { newDimContainer->origin[0], newDimContainer->origin[1], newDimContainer->origin[2] };
	int height = newDimContainer->dimensions[2];
	int length = newDimContainer->dimensions[0];
	int width = newDimContainer->dimensions[1];
	std::cout << "New dim : " << newDimContainer->dimensions[0] << " x " << newDimContainer->dimensions[1] << " x " << newDimContainer->dimensions[2] << "" << std::endl;

	size_t l = 0, w = 0, h = 0, xdo = 0;
	for (k = 0; k < height; k++) {
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {
				
				Point3D point_cropped = { i, j, k };
				point_cropped = getRealCoordFromImageCoord3D(point_cropped, newOrigin, ctSpacing, orientation);
				Point3D point = getImageCoordFromRealCoord3D(point_cropped, ctOrigin, ctSpacing, orientation);
				
				//if (point.x < 0) {
				//	l = 0;
				//}
				//else {
				//	if (point.x > length - 1) {
				//		l = length - 1;
				//	}
				//	else {
				//		l = (size_t)point.x;
				//	}
				//}
				//if (point.y < 0) {
				//	w = 0;
				//}
				//else {
				//	if (point.y > width - 1) {
				//		w = width - 1;
				//	}
				//	else {
				//		w = (size_t)point.y;
				//	}
				//}
				//if (point.z < 0) {
				//	h = 0;
				//}
				//else {
				//	if (point.z > height - 1) {
				//		h = height - 1;
				//	}
				//	else {
				//		h = (size_t)point.z;
				//	}
				//}

				l = (size_t)point.x;
				w = (size_t)point.y;
				h = (size_t)point.z;

				xdo = x_new(i, j, length);
				xd = x_new(l, w, Length);
				if (newDimContainer->dataPointer[k][xdo] != 0.0) {
					ctContainer->dataPointer[h][xd] = 1.0;
				}
			}
		}
	}
	
	storing_path = outputPath + "aorta_before_m.vtk";
	vtkDataForm dataForm = dta_binary;
	newDimContainer->operation = copyTo;
	storeVtkFile(storing_path.c_str(), ctContainer, dataForm);
	free(newDimContainer);
	*/

	//======================== Distance map computed on segmented aorta =========

	//fastSweepingFunction_3D(maskThreshold, ctContainer->dataPointer, Length, Width, Height, 1.0, 10000000.0, 0.0);
	//storing_path = outputPath + "distance_map.raw";
	//manageRAWFile3D<dataType>(maskThreshold, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	//======================== Load path .csv ===================================

	//FILE* path_file;
	//loading_path = inputPath + "vasculitis/patient_after_treatment/path_automatic.csv";
	//if (fopen_s(&path_file, loading_path.c_str(), "r") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//dataType x = 0, y = 0, z = 0;

	//while (feof(path_file) == 0) {
	//	fscanf_s(path_file, "%f", &x);
	//	fscanf_s(path_file, ",");
	//	fscanf_s(path_file, "%f", &y);
	//	fscanf_s(path_file, ",");
	//	fscanf_s(path_file, "%f", &z);
	//	fscanf_s(path_file, "\n");
	//	std::cout << "Point = (" << x << "," << y << "," << z << ")" << std::endl;
	//}
	//fclose(path_file);

	//======================== find Liver centroid ==============================

	/*
	thresholding3dFunctionN(imageData, Length, Width, Height, thres_min, thres_max, 0.0, 1.0);
	fastSweepingFunction_3D(maskThreshold, imageData, Length, Width, Height, 1.0, 10000000.0, 0.0);
	
	Point3D point_max = { 0, 0, 0 };
	dataType max_dist = 0.0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				xd = x_new(i, j, Length);
				if (maskThreshold[k][xd] > max_dist) {
					max_dist = maskThreshold[k][xd];
					point_max.x = i; 
					point_max.y = j;
					point_max.z = k;
				}
			}
		}
	}

	point_max = getRealCoordFromImageCoord3D(point_max, ctOrigin, ctSpacing, orientation);
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				xd = x_new(i, j, Length);
				Point3D current_point = { i, j, k };
				current_point = getRealCoordFromImageCoord3D(current_point, ctOrigin, ctSpacing, orientation);
				double dist = getPoint3DDistance(point_max, current_point);
				if (dist <= 10.0) {
					imageData[k][xd] = 1.0;
				}
				else {
					imageData[k][xd] = 0.0;
				}
			}
		}
	}

	storing_path = outputPath + "distance_map_p2.raw";
	manageRAWFile3D<dataType>(maskThreshold, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	storing_path = outputPath + "ball_p2.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);
	*/

	//======================== Refinement by Dilation ===========================

	/*
	dataType** interpol = new dataType * [height_int];
	for (k = 0; k < height_int; k++) {
		interpol[k] = new dataType[dim2D]{ 0 };
	}

	Image_Data interpolated;
	interpolated.imageDataPtr = interpol;
	interpolated.height = height_int; interpolated.length = Length; interpolated.width = Width;
	interpolated.orientation = orientation; interpolated.origin = ctOrigin;
	interpolated.spacing.sx = ctSpacing.sx;
	interpolated.spacing.sy = ctSpacing.sy;
	interpolated.spacing.sz = ctSpacing.sx;
	imageInterpolation3D(CT, interpolated, NEAREST_NEIGHBOR);

	//size_t k_min = 143, k_max = 236, kn = 0; //original p2
	//size_t k_min = 317, k_max = 500, kn = 0; //interpolated p2
	size_t k_min = 536, k_max = 712, kn = 0; //interpolated p3
	const size_t height = k_max - k_min + 1;
	std::cout << "cropped height : " << height << std::endl;
	dataType** imageCropData = new dataType * [height];
	dataType** ball_liver = new dataType * [height];
	for (k = 0; k < height; k++) {
		imageCropData[k] = new dataType[dim2D]{ 0 };
		ball_liver[k] = new dataType[dim2D]{ 0 };
	}

	loading_path = inputPath + "liver_p3.raw";
	load3dArrayRAW<dataType>(imageCropData, Length, Width, height, loading_path.c_str(), false);

	VoxelSpacing ctInt = { ctContainer->spacing[0], ctContainer->spacing[1], ctContainer->spacing[0] };
	Point3D newOrigin = { 0, 0, k_min };
	newOrigin = getRealCoordFromImageCoord3D(newOrigin, ctOrigin, ctInt, orientation);
	std::cout << "New origin : (" << newOrigin.x << "," << newOrigin.y << "," << newOrigin.z << ")" << std::endl;
	
	dilatation3dHeighteenNeigbours(imageCropData, Length, Width, height, 1.0, 0.0);
	dilatation3dHeighteenNeigbours(imageCropData, Length, Width, height, 1.0, 0.0);
	dilatation3dHeighteenNeigbours(imageCropData, Length, Width, height, 1.0, 0.0);
	dilatation3dHeighteenNeigbours(imageCropData, Length, Width, height, 1.0, 0.0);
	dilatation3dHeighteenNeigbours(imageCropData, Length, Width, height, 1.0, 0.0);
	dilatation3dHeighteenNeigbours(imageCropData, Length, Width, height, 1.0, 0.0);
	//dilatation3dHeighteenNeigbours(imageCropData, Length, Width, height, 1.0, 0.0);

	//crop refinement
	for (kn = 0, k = k_min; kn < height; kn++, k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				xd = x_new(i, j, Length);
				if (imageCropData[kn][xd] == 1.0 && (interpol[k][xd] < thres_min || interpol[k][xd] > thres_max)) {
					imageCropData[kn][xd] = 0.0;
				}
			}
		}
	}

	dataType* centroid = (dataType*)malloc(3 * sizeof(dataType));
	centroidImage(imageCropData, centroid, height, Length, Width, 0.0);
	std::cout << "Liver centroid : (" << (size_t)centroid[0] << "," << (size_t)centroid[1] << "," << (size_t)centroid[2] << ")" << std::endl;
	Point3D centroid_liver = { (size_t)centroid[0], (size_t)centroid[1], (size_t)centroid[2] };
	centroid_liver = getRealCoordFromImageCoord3D(centroid_liver, newOrigin, ctInt, orientation);

	for (k = 0; k < height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				xd = x_new(i, j, Length);
				Point3D current_point = { i, j, k };
				current_point = getRealCoordFromImageCoord3D(current_point, newOrigin, ctInt, orientation);
				double dist_point = getPoint3DDistance(centroid_liver, current_point);
				if (dist_point <= 10) {
					ball_liver[k][xd] = 1.0;
				}
			}
		}
	}

	storing_path = outputPath + "refined_p3.raw";
	store3dRawData<dataType>(imageCropData, Length, Width, height, storing_path.c_str());

	storing_path = outputPath + "ball_p3.raw";
	store3dRawData<dataType>(ball_liver, Length, Width, height, storing_path.c_str());

	free(centroid);
	for (k = 0; k < height; k++) {
		delete[] imageCropData[k];
		delete[] ball_liver[k];
	}
	delete[] imageCropData;
	delete[] ball_liver;

	for (k = 0; k < height_int; k++) {
		delete[] interpol[k];
	}
	delete[] interpol;
	*/

	//======================== Collect PET values ===============================
	
	/*
	string beforePath = "C:/Users/Konan Allaly/Documents/Tests/input/vasculitis/patient_before_treatment/";
	string afterPath = "C:/Users/Konan Allaly/Documents/Tests/input/vasculitis/patient_after_treatment/";
	
	// Load PET
	Vtk_File_Info* petContainer = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	petContainer->operation = copyFrom;
	//loading_path = inputPath + "vtk/pet/Patient2_pet.vtk";
	loading_path = inputPath + "vtk/pet/Patient3_pet.vtk";
	readVtkFile(loading_path.c_str(), petContainer);

	int height = petContainer->dimensions[2];
	int length = petContainer->dimensions[0];
	int width = petContainer->dimensions[1];
	int dim2d = length * width;
	std::cout << "#######################" << std::endl;
	std::cout << "PET image dim : " << petContainer->dimensions[0] << " x " << petContainer->dimensions[1] << " x " << petContainer->dimensions[2] << "" << std::endl;

	std::cout << "PET origin : (" << petContainer->origin[0] << ", " << petContainer->origin[1] << ", " << petContainer->origin[2] << ")" << std::endl;
	Point3D petOrigin = { petContainer->origin[0], petContainer->origin[1], petContainer->origin[2] };
	VoxelSpacing petSpacing = { petContainer->spacing[0], petContainer->spacing[1], petContainer->spacing[2] };
	std::cout << "PET spacing : (" << petContainer->spacing[0] << ", " << petContainer->spacing[1] << ", " << petContainer->spacing[2] << ")" << std::endl;
	std::cout << "#######################" << std::endl;

	Image_Data PET; 
	PET.imageDataPtr = petContainer->dataPointer;
	PET.height = height;
	PET.length = length;
	PET.width = width;
	PET.orientation = orientation;
	PET.origin = petOrigin;
	PET.spacing = petSpacing;

	//Save uptake ratio
	//string saving_csv = outputPath + "ratio_before_auto_r20.csv";
	//string saving_csv = outputPath + "ratio_before_manual.csv";
	//string saving_csv = outputPath + "ratio_after_auto_r20.csv";
	string saving_csv = outputPath + "ratio_after_manual.csv";

	FILE* file;
	if (fopen_s(&file, saving_csv.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	fprintf(file, "x,y\n");

	//Load segments
	Vtk_File_Info* ballLiver = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	ballLiver->operation = copyFrom;
	//loading_path = beforePath + "ball_before_auto.vtk";
	//loading_path = beforePath + "ball_liver_before_m.vtk";
	//loading_path = afterPath + "ball_after_auto.vtk";
	loading_path = afterPath + "ball_liver_after_m.vtk";
	readVtkFile(loading_path.c_str(), ballLiver);
	//std::cout << "Origin ball : (" << ballLiver->origin[0] << ", " << ballLiver->origin[1] << ", " << ballLiver->origin[2] << ")" << std::endl;
	//std::cout << "dim ball liver : " << ballLiver->dimensions[0] << "x" << ballLiver->dimensions[1] << "x" << ballLiver->dimensions[2] << std::endl;

	Vtk_File_Info* segment1 = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	segment1->operation = copyFrom;
	//loading_path = beforePath + "seg1_before.vtk";
	loading_path = afterPath + "seg1_after.vtk";
	readVtkFile(loading_path.c_str(), segment1);
	//std::cout << "Origin seg 1 : (" << segment1->origin[0] << ", " << segment1->origin[1] << ", " << segment1->origin[2] << ")" << std::endl;
	//std::cout << "dim segment 1 : " << segment1->dimensions[0] << "x" << segment1->dimensions[1] << "x" << segment1->dimensions[2] << std::endl;

	Vtk_File_Info* segment2 = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	segment2->operation = copyFrom;
	//loading_path = beforePath + "seg2_before.vtk";
	loading_path = afterPath + "seg2_after.vtk";
	readVtkFile(loading_path.c_str(), segment2);
	//std::cout << "Origin seg 2 : (" << segment2->origin[0] << ", " << segment2->origin[1] << ", " << segment2->origin[2] << ")" << std::endl;
	//std::cout << "dim segment 2 : " << segment2->dimensions[0] << "x" << segment2->dimensions[1] << "x" << segment2->dimensions[2] << std::endl;

	Vtk_File_Info* segment3 = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	segment3->operation = copyFrom;
	//loading_path = beforePath + "seg3_before.vtk";
	loading_path = afterPath + "seg3_after.vtk";
	readVtkFile(loading_path.c_str(), segment3);
	//std::cout << "Origin seg 3 : (" << segment3->origin[0] << ", " << segment3->origin[1] << ", " << segment3->origin[2] << ")" << std::endl;
	//std::cout << "dim segment 3 : " << segment3->dimensions[0] << "x" << segment3->dimensions[1] << "x" << segment3->dimensions[2] << std::endl;

	Vtk_File_Info* segment4 = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	segment4->operation = copyFrom;
	//loading_path = beforePath + "seg4_before.vtk";
	loading_path = afterPath + "seg4_after.vtk";
	readVtkFile(loading_path.c_str(), segment4);
	//std::cout << "Origin seg 4 : (" << segment4->origin[0] << ", " << segment4->origin[1] << ", " << segment4->origin[2] << ")" << std::endl;
	//std::cout << "dim segment 4 : " << segment4->dimensions[0] << "x" << segment4->dimensions[1] << "x" << segment4->dimensions[2] << std::endl;

	Vtk_File_Info* segment5 = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	segment5->operation = copyFrom;
	//loading_path = beforePath + "seg5_before.vtk";
	loading_path = afterPath + "seg5_after.vtk";
	readVtkFile(loading_path.c_str(), segment5);
	//std::cout << "Origin seg 5 : (" << segment5->origin[0] << ", " << segment5->origin[1] << ", " << segment5->origin[2] << ")" << std::endl;
	//std::cout << "dim segment 5 : " << segment5->dimensions[0] << "x" << segment5->dimensions[1] << "x" << segment5->dimensions[2] << std::endl;

	size_t l = 0, w = 0, h = 0;

	//find max uptake in Liver
	dataType max_uptake_liver = 0.0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				xd = x_new(i, j, Length);
				if (ballLiver->dataPointer[k][xd] != 0.0) {
					Point3D point_ball = { i, j, k };
					point_ball = getRealCoordFromImageCoord3D(point_ball, ctOrigin, ctSpacing, orientation);
					point_ball = getImageCoordFromRealCoord3D(point_ball, petOrigin, petSpacing, orientation);
					l = (size_t)point_ball.x;
					w = (size_t)point_ball.y;
					h = (size_t)point_ball.z;
					dataType value = petContainer->dataPointer[h][x_new(l, w, length)];
					if (value > max_uptake_liver) {
						max_uptake_liver = value;
					}
				}
			}
		}
	}
	std::cout << "Max uptake liver = " << max_uptake_liver << std::endl;
	if (max_uptake_liver == 0) {
		max_uptake_liver = 1.0; //avoid division by 0 in next step
	}

	////Load path points : for automatic segment
	//FILE* path_file;
	////loading_path = beforePath + "path_before.csv";
	//loading_path = afterPath + "path_after.csv";
	//if (fopen_s(&path_file, loading_path.c_str(), "r") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//vector<Point3D> points_on_path;
	//Point3D current_point = { 0.0, 0.0, 0.0 };
	//while (feof(path_file) == 0) {
	//	fscanf_s(path_file, "%f", &current_point.x);
	//	fscanf_s(path_file, ",");
	//	fscanf_s(path_file, "%f", &current_point.y);
	//	fscanf_s(path_file, ",");
	//	fscanf_s(path_file, "%f", &current_point.z);
	//	fscanf_s(path_file, "\n");
	//	points_on_path.push_back(current_point);
	//	//std::cout << "Point = (" << current_point.x << "," << current_point.y << "," << current_point.z << ")" << std::endl;
	//}
	//fclose(path_file);

	//Profil ratio segment manual experiment
	dataType max1 = 0.0, max2 = 0.0, max3 = 0.0, max4 = 0.0, max5 = 0.0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				xd = x_new(i, j, Length);
				Point3D point_pet = { i, j, k };
				point_pet = getRealCoordFromImageCoord3D(point_pet, ctOrigin, ctSpacing, orientation);
				point_pet = getImageCoordFromRealCoord3D(point_pet, petOrigin, petSpacing, orientation);
				if (point_pet.x < 0) {
					l = 0;
				}
				else {
					if (point_pet.x > length - 1) {
						l = length - 1;
					}
					else {
						l = (size_t)point_pet.x;
					}
				}
				if (point_pet.y < 0) {
					w = 0;
				}
				else {
					if (point_pet.y > width - 1) {
						w = width - 1;
					}
					else {
						w = (size_t)point_pet.y;
					}
				}
				if (point_pet.z < 0) {
					h = 0;
				}
				else {
					if (point_pet.z > height - 1) {
						h = height - 1;
					}
					else {
						h = (size_t)point_pet.z;
					}
				}
				dataType value = petContainer->dataPointer[h][x_new(l, w, length)];
				//find the max in segment 1 if the point lies in segment 1
				if (segment1->dataPointer[k][xd] != 0) {
					if (value > max1) {
						max1 = value;
					}
				}
				//find the max in segment 2 if the point lies in segment 2
				if (segment2->dataPointer[k][xd] != 0) {
					if (value > max2) {
						max2 = value;
					}
				}
				//find the max in segment 3 if the point lies in segment 3
				if (segment3->dataPointer[k][xd] != 0) {
					if (value > max3) {
						max3 = value;
					}
				}
				//find the max in segment 4 if the point lies in segment 4
				if (segment4->dataPointer[k][xd] != 0) {
					if (value > max4) {
						max4 = value;
					}
				}
				//find the max in segment 5 if the point lies in segment 5
				if (segment5->dataPointer[k][xd] != 0) {
					if (value > max5) {
						max5 = value;
					}
				}
			}
		}
	}

	//int n = 0;
	//double radius = 20.0;
	//Statistics stats = { 0.0, 0.0, 0.0, 0.0 };
	//size_t count1 = 0, count2 = 0, count3 = 0, count4 = 0, count5 = 0;

	//for (n = 0; n < points_on_path.size(); n++) {
	//	//std::cout << "Point realworld = (" << points_on_path[n].x << "," << points_on_path[n].y << "," << points_on_path[n].z << ")" << std::endl;
	//	Point3D point_pet = getImageCoordFromRealCoord3D(points_on_path[n], petOrigin, petSpacing, orientation);
	//	//std::cout << "Point image PET = (" << point_pet.x << "," << point_pet.y << "," << point_pet.z << ")" << std::endl;
	//	Point3D point_ct = getImageCoordFromRealCoord3D(points_on_path[n], ctOrigin, ctSpacing, orientation);
	//	//std::cout << "Point image CT = (" << point_ct.x << "," << point_ct.y << "," << point_ct.z << ")" << std::endl;
	//	i = (size_t)point_ct.x;
	//	j = (size_t)point_ct.y;
	//	k = (size_t)point_ct.z;
	//	xd = x_new(i, j, Length);
	//	//std::cout << "Point image CT = (" << i << "," << j << "," << k << ")" << std::endl;
	//	stats = getPointNeighborhoodStats(PET, point_pet, radius);
	//	//std::cout << "min = " << stats.min_data << std::endl;
	//	//std::cout << "mean = " << stats.mean_data << std::endl;
	//	//std::cout << "max = " << stats.max_data << std::endl;
	//	//find the max in segment 1 if the point lies in segment 1
	//	if (segment1->dataPointer[k][xd] != 0) {
	//		count1++;
	//		if (stats.max_data > max1) {
	//			max1 = stats.max_data;
	//		}
	//	}
	//	//find the max in segment 2 if the point lies in segment 2
	//	if (segment2->dataPointer[k][xd] != 0) {
	//		count2++;
	//		if (stats.max_data > max2) {
	//			max2 = stats.max_data;
	//		}
	//	}	
	//	//find the max in segment 3 if the point lies in segment 3
	//	if (segment3->dataPointer[k][xd] != 0) {
	//		count3++;;
	//		if (stats.max_data > max3) {
	//			max3 = stats.max_data;
	//		}
	//	}
	//	//find the max in segment 4 if the point lies in segment 4
	//	if (segment4->dataPointer[k][xd] != 0) {
	//		count4++;
	//		if (stats.max_data > max4) {
	//			max4 = stats.max_data;
	//		}
	//	}
	//	//find the max in segment 5 if the point lies in segment 5
	//	if (segment5->dataPointer[k][xd] != 0) {
	//		count5++;
	//		if (stats.max_data > max5) {
	//			max5 = stats.max_data;
	//		}
	//	}
	//}
	
	//std::cout << count1 << " elements in segment 1" << std::endl;
	//std::cout << count2 << " elements in segment 2" << std::endl;
	//std::cout << count3 << " elements in segment 3" << std::endl;
	//std::cout << count4 << " elements in segment 4" << std::endl;
	//std::cout << count5 << " elements in segment 5" << std::endl;
	
	std::cout << "max1 = " << max1 << ", max2 = " << max2 << ", max3 = " << max3 << ", max4 = " << max4 << ", max5 = " << max5 << std::endl;

	size_t count = 1;
	fprintf(file, "%d,%f\n", count, max1 / max_uptake_liver);
	count++;
	fprintf(file, "%d,%f\n", count, max2 / max_uptake_liver);
	count++;
	fprintf(file, "%d,%f\n", count, max3 / max_uptake_liver);
	count++;
	fprintf(file, "%d,%f\n", count, max4 / max_uptake_liver);
	count++;
	fprintf(file, "%d,%f\n", count, max5 / max_uptake_liver);
	fclose(file);

	free(ballLiver);
	free(segment1);
	free(segment2);
	free(segment3);
	free(segment4);
	free(segment5);
	free(petContainer);

	//while (points_on_path.size() > 0) {
	//	points_on_path.pop_back();
	//}
	*/
	
	//======================== 2D path finding =================
	
	/*
	const size_t length = 256, width = 256;
	dataType* image = new dataType[length * width] {0};
	
	//Point2D seed = { 128, 200 };
	//for (i = 0; i < length; i++) {
	//	for (j = 0; j < 175; j++) {
	//		Point2D current_point = { i, j };
	//		double distance = getPoint2DDistance(seed, current_point);
	//		if (distance > 75 && distance < 80) {
	//			image[x_new(i - 15, j, length)] = 1.0;
	//		}
	//	}
	//}
	//i = 185;
	//for (j = 10; j < 174; j++) {
	//	image[x_new(i - 1, j, length)] = 1.0;
	//	image[x_new(i - 2, j, length)] = 1.0;
	//	image[x_new(i, j, length)] = 1.0;
	//	image[x_new(i + 1, j, length)] = 1.0;
	//	image[x_new(i + 2, j, length)] = 1.0;
	//}

	storing_path = outputPath + "imageSin.raw";
	store2dRawData<dataType>(image, length, width, storing_path.c_str());
	
	//loading_path = inputPath + "raw/slice/retine2.raw";
	//load2dArrayRAW<short>(imageSlice, length, width, loading_path.c_str(), false);

	//xd = x_new(100, 100, length);
	//cout << "before : " << imageSlice[xd] << endl;

	//for (i = 0; i < length * width; i++) {
	//	image[i] = (dataType)imageSlice[i];
	//}

	//cout << "after : " << image[xd] << endl;

	//rescaleNewRange2D(image, length, width, 0.0, 1.0);

	//storing_path = outputPath + "retine.raw";
	//store2dRawData<dataType>(image, length, width, storing_path.c_str());

	//delete[] imageSlice;
	//delete[] image;

	delete[] image;
	*/

	//======================== Interplation ====================

	/*
	dataType** interpolated = new dataType * [height];
	for (k = 0; k < height; k++) {
		interpolated[k] = new dataType[dim2D];
	}

	Image_Data INTERPOL;
	INTERPOL.imageDataPtr = interpolated;
	INTERPOL.height = height; INTERPOL.length = Length; INTERPOL.width = Width;
	INTERPOL.orientation = orientation; INTERPOL.origin = ctOrigin;
	INTERPOL.spacing.sx = ctSpacing.sx;
	INTERPOL.spacing.sy = ctSpacing.sy;
	INTERPOL.spacing.sz = ctSpacing.sx;
	imageInterpolation3D(CT, INTERPOL, NEAREST_NEIGHBOR);
	*/

	//======================== Filtering =======================

	/*
	Filter_Parameters filter_parameters; filter_parameters.edge_detector_coefficient = 1;
	filter_parameters.h = 1.0; filter_parameters.timeStepSize = 1.2; filter_parameters.eps2 = 1e-6;
	filter_parameters.omega_c = 1.5; filter_parameters.tolerance = 1e-2; filter_parameters.maxNumberOfSolverIteration = 100;
	filter_parameters.timeStepsNum = 1; filter_parameters.coef = 1e-6; 
	filterImage(INTERPOL, filter_parameters, GEODESIC_MEAN_CURVATURE_FILTER);
	
	storing_path = outputPath + "filteredGMC.raw";
	manageRAWFile3D<dataType>(interpolated, Length, Width, height, storing_path.c_str(), STORE_DATA, false);

	for (k = 0; k < height; k++) {
		delete[] interpolated[k];
	}
	delete[] interpolated;
	*/

	//======================== Segment the Liver ===============

	/*
	dataType thres1 = 995, thres2 = 1213;
	thresholding3dFunctionN(maskThreshold, Length, Width, Height, thres1, thres2, 0.0, 1.0);

	erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);

	int** labelArray = new int* [Height];
	bool** status = new bool* [Height];
	for (k = 0; k < Height; k++) {
		labelArray[k] = new int[dim2D] {0};
		status[k] = new bool[dim2D];
	}
	if (labelArray == NULL || status == NULL)
		return false;

	//Initialization
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			status[k][i] = false;
		}
	}

	int numberOfRegionCells = countNumberOfRegionsCells(maskThreshold, Length, Width, Height, 1.0);

	labelling3D(maskThreshold, labelArray, status, Length, Width, Height, 1.0);

	//Counting
	int* countingArray = new int[numberOfRegionCells] {0};
	if (countingArray == NULL)
		return false;

	//Get regions sizes
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (labelArray[k][i] > 0) {
				countingArray[labelArray[k][i]]++;
			}
		}
	}

	//Number of regions
	int numberOfRegions = 0;
	for (i = 0; i < numberOfRegionCells; i++) {
		if (countingArray[i] > 0) {
			numberOfRegions++;
		}
	}
	cout << numberOfRegions << " regions found" << endl;

	//Find the biggest region
	size_t regionSize = 0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (countingArray[labelArray[k][i]] > regionSize) {
				regionSize = countingArray[labelArray[k][i]];
			}
		}
	}

	//Keep and save the biggest region
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (countingArray[labelArray[k][i]] == regionSize) {
				maskThreshold[k][i] = 1.0;
			}
			else {
				maskThreshold[k][i] = 0.0;
			}
		}
	}

	dilatation3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	dilatation3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);

	dataType* centroid = (dataType*)malloc(3 * sizeof(dataType));
	centroidImage(maskThreshold, centroid, Height, Length, Width, 0.0);
	cout << "Liver centroid : (" << centroid[0] << "," << centroid[1] << "," << centroid[2] << ")" << endl;

	storing_path = outputPath + "liver.raw";
	store3dRawData<dataType>(maskThreshold, Length, Width, Height, storing_path.c_str());

	free(centroid);
	for (k = 0; k < Height; k++) {
		delete[] labelArray[k];
		delete[] status[k];
	}
	delete[] labelArray;
	delete[] status;
	*/

	//======================== Segment Liver on cropped data ======================

	/*
	dataType** interpol = new dataType * [height_int];
	for (k = 0; k < height_int; k++) {
		interpol[k] = new dataType[dim2D]{0};
	}

	Image_Data interpolated;
	interpolated.imageDataPtr = interpol;
	interpolated.height = height_int; interpolated.length = Length; interpolated.width = Width;
	interpolated.orientation = orientation; interpolated.origin = ctOrigin;
	interpolated.spacing.sx = ctSpacing.sx;
	interpolated.spacing.sy = ctSpacing.sy;
	interpolated.spacing.sz = ctSpacing.sx;
	imageInterpolation3D(CT, interpolated, NEAREST_NEIGHBOR);

	//size_t k_min = 143, k_max = 236, kn = 0; //original p2
	//size_t k_min = 317, k_max = 500, kn = 0; //interpolated p2
	size_t k_min = 536, k_max = 712, kn = 0; //interpolated p3
	const size_t height = k_max - k_min + 1;
	std::cout << "cropped height : " << height << std::endl;

	dataType** imageCropData = new dataType * [height];
	dataType** maskCropData = new dataType * [height];
	for (k = 0; k < height; k++) {
		imageCropData[k] = new dataType[dim2D]{0};
		maskCropData[k] = new dataType[dim2D]{0};
	}

	//cropping
	for (kn = 0, k = k_min; kn < height; kn++, k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				xd = x_new(i, j, Length);
				//imageCropData[kn][xd] = imageData[k][xd];
				//maskCropData[kn][xd] = imageData[k][xd];
				imageCropData[kn][xd] = interpol[k][xd];
				maskCropData[kn][xd] = interpol[k][xd];
			}
		}
	}

	VoxelSpacing ctInt = { ctContainer->spacing[0], ctContainer->spacing[1], ctContainer->spacing[0] };
	Point3D newOrigin = { 0, 0, k_min };
	newOrigin = getRealCoordFromImageCoord3D(newOrigin, ctOrigin, ctInt, orientation);
	std::cout << "New origin : (" << newOrigin.x << "," << newOrigin.y << "," << newOrigin.z << ")" << std::endl;

	thresholding3dFunctionN(maskCropData, Length, Width, height, thres_min, thres_max, 0.0, 1.0);

	erosion3dHeighteenNeigbours(maskCropData, Length, Width, height, 1.0, 0.0);
	erosion3dHeighteenNeigbours(maskCropData, Length, Width, height, 1.0, 0.0);
	erosion3dHeighteenNeigbours(maskCropData, Length, Width, height, 1.0, 0.0);
	erosion3dHeighteenNeigbours(maskCropData, Length, Width, height, 1.0, 0.0);
	erosion3dHeighteenNeigbours(maskCropData, Length, Width, height, 1.0, 0.0);
	erosion3dHeighteenNeigbours(maskCropData, Length, Width, height, 1.0, 0.0);
	erosion3dHeighteenNeigbours(maskCropData, Length, Width, height, 1.0, 0.0);
	erosion3dHeighteenNeigbours(maskCropData, Length, Width, height, 1.0, 0.0);

	int** labelArray = new int* [height];
	bool** status = new bool* [height];
	for (k = 0; k < height; k++) {
		labelArray[k] = new int[dim2D] {0};
		status[k] = new bool[dim2D];
	}
	if (labelArray == NULL || status == NULL)
		return false;

	//Initialization
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			status[k][i] = false;
		}
	}

	int numberOfRegionCells = countNumberOfRegionsCells(maskCropData, Length, Width, height, 1.0);

	labelling3D(maskCropData, labelArray, status, Length, Width, height, 1.0);

	//Counting
	int* countingArray = new int[numberOfRegionCells] {0};
	if (countingArray == NULL)
		return false;

	//Get regions sizes
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (labelArray[k][i] > 0) {
				countingArray[labelArray[k][i]]++;
			}
		}
	}

	//Number of regions
	int numberOfRegions = 0;
	for (i = 0; i < numberOfRegionCells; i++) {
		if (countingArray[i] > 0) {
			numberOfRegions++;
		}
	}
	cout << numberOfRegions << " regions found" << endl;

	//Find the biggest region
	size_t regionSize = 0;
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (countingArray[labelArray[k][i]] > regionSize) {
				regionSize = countingArray[labelArray[k][i]];
			}
		}
	}

	//Keep and save the biggest region
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (countingArray[labelArray[k][i]] == regionSize) {
				maskCropData[k][i] = 1.0;
			}
			else {
				maskCropData[k][i] = 0.0;
			}
		}
	}

	//dilatation3dHeighteenNeigbours(maskCropData, Length, Width, height, 1.0, 0.0);
	//dilatation3dHeighteenNeigbours(maskCropData, Length, Width, height, 1.0, 0.0);

	dataType* centroid = (dataType*)malloc(3 * sizeof(dataType));
	centroidImage(maskCropData, centroid, height, Length, Width, 0.0);
	std::cout << "Liver centroid : (" << centroid[0] << "," << centroid[1] << "," << centroid[2] << ")" << std::endl;

	storing_path = outputPath + "liver_p3.raw";
	store3dRawData<dataType>(maskCropData, Length, Width, height, storing_path.c_str());

	free(centroid);
	for (k = 0; k < height; k++) {
		delete[] labelArray[k];
		delete[] status[k];
		delete[] imageCropData[k];
		delete[] maskCropData[k];
	}
	delete[] labelArray;
	delete[] status;
	delete[] imageCropData;
	delete[] maskCropData;

	for (k = 0; k < height_int; k++) {
		delete[] interpol[k];
	}
	delete[] interpol;
	*/

	//======================== Refinement by generalized subjective surface =======
	
	/*
	dataType** filtered = new dataType * [height_int];
	for (k = 0; k < height_int; k++) {
		filtered[k] = new dataType[dim2D];
	}

	loading_path = inputPath + "patient/before/filteredGMC.raw";
	load3dArrayRAW<dataType>(filtered, Length, Width, height_int, loading_path.c_str(), false);

	//size_t k_min = 143, k_max = 236, kn = 0; //original p2
	size_t k_min = 317, k_max = 500, kn = 0; //interpolated p2
	//size_t k_min = 536, k_max = 712, kn = 0; //interpolated p3
	const size_t height = k_max - k_min + 1;
	std::cout << "cropped height : " << height << std::endl;

	dataType** imageCropData = new dataType * [height];
	dataType** initial_seg = new dataType * [height];
	for (k = 0; k < height; k++) {
		imageCropData[k] = new dataType[dim2D]{ 0 };
		initial_seg[k] = new dataType[dim2D]{ 0 };
	}

	loading_path = inputPath + "liver_p2.raw";
	load3dArrayRAW<dataType>(initial_seg, Length, Width, height, loading_path.c_str(), false);

	//cropping
	for (kn = 0, k = k_min; kn < height; kn++, k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				xd = x_new(i, j, Length);
				imageCropData[kn][xd] = filtered[k][xd];
			}
		}
	}

	VoxelSpacing ctInt = { ctContainer->spacing[0], ctContainer->spacing[1], ctContainer->spacing[0] };
	Point3D newOrigin = { 0, 0, k_min };
	newOrigin = getRealCoordFromImageCoord3D(newOrigin, ctOrigin, ctInt, orientation);
	std::cout << "New origin : (" << newOrigin.x << "," << newOrigin.y << "," << newOrigin.z << ")" << std::endl;

	Image_Data imageToBeSegmented;
	imageToBeSegmented.imageDataPtr = imageCropData;
	imageToBeSegmented.height = height; imageToBeSegmented.length = Length; imageToBeSegmented.width = Width;
	imageToBeSegmented.orientation = orientation; 
	imageToBeSegmented.origin = newOrigin;
	imageToBeSegmented.spacing = ctInt;
	
	Point3D* center_seg = new Point3D[1];
	center_seg->x = 0.0, center_seg->y = 0.0; center_seg->z = 0.0;

	Segmentation_Parameters seg_params;
	seg_params.tau = 1.2; seg_params.h = 1.0; seg_params.omega_c = 1.5; seg_params.mod = 5;
	seg_params.maxNoGSIteration = 40; seg_params.coef = 100000; seg_params.gauss_seidelTolerance = 1e-6;
	seg_params.maxNoOfTimeSteps = 100; seg_params.segTolerance = 1e-6; seg_params.eps2 = 1e-6;
	seg_params.coef_conv = 1.0; seg_params.coef_dif = 1.0;

	//Smoothing by heat equation so we don't need all the parameters
	Filter_Parameters smooth_params;
	smooth_params.h = seg_params.h; 
	smooth_params.timeStepSize = smooth_params.h; 
	smooth_params.timeStepsNum = 1; 
	smooth_params.omega_c = 1.5;
	smooth_params.tolerance = 1e-3, smooth_params.maxNumberOfSolverIteration = 100;

	unsigned char segmentPath[] = "C:/Users/Konan Allaly/Documents/Tests/output/segmentation/p2/";
	subsurfSegmentation(imageToBeSegmented, initial_seg, seg_params, smooth_params, center_seg, 1, segmentPath);
	//generalizedSubsurfSegmentation(imageToBeSegmented, initial_seg, seg_params, smooth_params, center_seg, 1, segmentPath);

	delete[] center_seg;
	for (k = 0; k < height; k++) {
		delete[] imageCropData[k];
		delete[] initial_seg[k];
	}
	delete[] imageCropData;
	delete[] initial_seg;

	for (k = 0; k < height_int; k++) {
		delete[] filtered[k];
	}
	delete[] filtered;
	*/

	//======================== Path finfing =======================================

	
	dataType** potential = new dataType * [height];
	dataType** action_field = new dataType * [height];
	dataType** pathPtr = new dataType * [height];
	for (k = 0; k < height; k++) {
		potential[k] = new dataType[dim2D]{ 0 };
		action_field[k] = new dataType[dim2D]{ 0 };
		pathPtr[k] = new dataType[dim2D]{ 0 };
	}
	if (potential == NULL || action_field == NULL || pathPtr == NULL)
		return false;

	////Patient 2
	//Point3D seed1 = { 263, 257, 146 };
	//Point3D seed2 = { 294, 315, 261 };
	//Point3D seed3 = { 252, 254, 261 };

	Point3D seed1 = { 263, 257, 146 }; // p2
	//Point3D seed1 = { 262, 243, 469 }; //p6
	//Point3D seed1 = { 269, 231, 111 }; //p4
	//Point3D seed1 = { 264, 263, 176 }; //p2 liver
	seed1 = getRealCoordFromImageCoord3D(seed1, ctOrigin, ctSpacing, orientation);
	seed1 = getImageCoordFromRealCoord3D(seed1, ctOrigin, INTERPOLATE.spacing, orientation);

	Point3D seed2 = { 293, 315, 257 }; // p2 : automatic detection
	//Point3D seed2 = { 294, 315, 261 };//p2
	//Point3D seed2 = { 281, 284, 644 };//p6
	//Point3D seed2 = { 287, 288, 229 };//p4
	seed2 = getRealCoordFromImageCoord3D(seed2, ctOrigin, ctSpacing, orientation);
	seed2 = getImageCoordFromRealCoord3D(seed2, ctOrigin, INTERPOLATE.spacing, orientation);

	Point3D seed3 = { 254, 248, 257 };// p2 : automatic detection
	//Point3D seed3 = { 252, 254, 261 };//p2
	//Point3D seed3 = { 245, 215, 644 };//p6
	//Point3D seed3 = { 254, 220, 229 };//p4
	seed3 = getRealCoordFromImageCoord3D(seed3, ctOrigin, ctSpacing, orientation);
	seed3 = getImageCoordFromRealCoord3D(seed3, ctOrigin, INTERPOLATE.spacing, orientation);

	////Patient before
	//Point3D seed1 = { 260, 254, 318 };
	//Point3D seed2 = { 295, 316, 565 };
	//Point3D seed3 = { 256, 256, 565 };

	////Patient after
	//Point3D seed1 = { 258, 255, 527 };
	//Point3D seed2 = { 288, 315, 774 };
	//Point3D seed3 = { 248, 251, 774 };

	Point3D* seedPoints = new Point3D[3];
	seedPoints[0] = seed1;
	seedPoints[1] = seed2;
	seedPoints[2] = seed3;

	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			xd = x_new(i, j, Length);
	//			Point3D currentPoint = { i, j, k };
	//			currentPoint = getRealCoordFromImageCoord3D(currentPoint, ctOrigin, ctSpacing, orientation);
	//			double d1 = getPoint3DDistance(seed1, currentPoint);
	//			double d2 = getPoint3DDistance(seed2, currentPoint);
	//			double d3 = getPoint3DDistance(seed3, currentPoint);
	//			if (d1 <= 4) {
	//				maskThreshold[k][xd] = 1.0;
	//			}
	//			if (d2 <= 4) {
	//				maskThreshold[k][xd] = 1.0;
	//			}
	//			if (d3 <= 4) {
	//				maskThreshold[k][xd] = 1.0;
	//			}
	//		}
	//	}
	//}

	//storing_path = outputPath + "keys_points_p2_liver.raw";
	//store3dRawData<dataType>(maskThreshold, Length, Width, Height, storing_path.c_str());

	Potential_Parameters parameters;
	parameters.K = 0.005; parameters.epsilon = 0.01;
	parameters.c_ct = 0.0; parameters.c_pet = 0.0;
	parameters.c_max = 0.0; parameters.c_min = 0.0; 
	parameters.c_mean = 1.0; parameters.c_sd = 0.0;
	double radius = 3.0;
	
	//computePotentialNew(interpolated, mean_image, potential, seed1, radius, parameters);
	//fastMarching3D_N(action_field, potential, Length, Width, height, seed1);
	//partialFrontPropagation(action_field, potential, Length, Width, height, seedPoints);

	//storing_path = outputPath + "partial_action_field.raw";
	//store3dRawData<dataType>(action_field, Length, Width, height, storing_path.c_str());

	//Point3D* seed = new Point3D[2];
	//seed[0] = seed1; seed[1] = seed3;
	//vector<Point3D> path_points;
	//shortestPath3D(action_field, pathPtr, Length, Width, height, 1.0, seed, path_points);

	//string saving_csv = outputPath + "min_per_slice.csv";
	//FILE* file;
	//if (fopen_s(&file, saving_csv.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}

	//Point3D point_of_max;
	//for (k = 0; k < height; k++) {
	//	dataType min_per_slice = 100000.0;
	//	point_of_max = { 0.0, 0.0, 0.0 };
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			xd = x_new(i, j, Length);
	//			if (action_field[k][xd] < min_per_slice) {
	//				min_per_slice = action_field[k][xd];
	//				point_of_max.x = i;
	//				point_of_max.y = j;
	//				point_of_max.z = k;
	//			}
	//		}
	//	}
	//	point_of_max = getRealCoordFromImageCoord3D(point_of_max, ctOrigin, interpolated.spacing, orientation);
	//	fprintf(file, "%f,%f,%f\n", point_of_max.x, point_of_max.y, point_of_max.z);
	//}
	//fclose(file);

	//for (i = 0; i < path_points.size(); i++) {
	//	path_points[i] = getRealCoordFromImageCoord3D(path_points[i], ctOrigin, interpolated.spacing, orientation);
	//	fprintf(file, "%f,%f,%f\n", path_points[i].x, path_points[i].y, path_points[i].z);
	//}
	//fclose(file);

	findPathTwoSteps(INTERPOLATE, pathPtr, seedPoints, parameters);

	//storing_path = outputPath + "path_before.raw";
	//store3dRawData<dataType>(pathPtr, Length, Width, height, storing_path.c_str());

	//delete[] seed;
	delete[] seedPoints;
	for (k = 0; k < height; k++) {
		delete[] potential[k];
		delete[] action_field[k];
		delete[] pathPtr[k];
	}
	delete[] potential;
	delete[] action_field;
	delete[] pathPtr;
	

	//======================== 2D Hough transform on Slice ===========================
	
	/*
	//loading_path = inputPath + "raw/filtered/filteredGMC_p2.raw";
	//load3dArrayRAW(imageData, Length, Width, Height, loading_path.c_str(), false);
	
	dataType* imageSlice = new dataType[dim2D]{ 0 };
	dataType* foundCircle = new dataType[dim2D]{ 0 };
	dataType* houghSpace = new dataType[dim2D]{ 0 };

	//size_t k_trachea = 257; //--->P2 : automatic detection
	//size_t k_abdo_aorta = 146; //--->P2 : abdominal aorta
	//size_t k_trachea = 261; //--->P2 : automatic detection
	//size_t k_trachea = 230; //--->P4
	size_t k_liver = 176; //--->P2
	for (i = 0; i < dim2D; i++) {
		imageSlice[i] = imageData[k_liver][i];
		foundCircle[i] = 0.0;
		houghSpace[i] = 0.0;
	}

	rescaleNewRange2D(imageSlice, Length, Width, 0.0, 1.0);

	//storing_path = outputPath + "filtered.raw";
	//store2dRawData<dataType>(imageSlice, Length, Width, storing_path.c_str());

	Image_Data2D IMAGE;
	IMAGE.imageDataPtr = imageSlice;
	IMAGE.height = Length; IMAGE.width = Width;
	
	PixelSpacing spacing = { ctSpacing.sx, ctSpacing.sy };
	HoughParameters parameters = { 7.5, 19.0, 0.5, 30.0, 1000.0, 1.0, 0.07, 0.5, spacing };
	
	//filtering parameters
	Filter_Parameters filter_parameters;
	filter_parameters.h = 1.0; filter_parameters.timeStepSize = 0.3; filter_parameters.eps2 = 1e-6;
	filter_parameters.omega_c = 1.5; filter_parameters.tolerance = 1e-3; filter_parameters.maxNumberOfSolverIteration = 1000;
	filter_parameters.timeStepsNum = 5; filter_parameters.coef = 1e-6;

	//Slice filtering
	//heatImplicit2dScheme(IMAGE, filter_parameters);

	//storing_path = outputPath + "smoothed.raw";
	//store2dRawData<dataType>(imageSlice, Length, Width, storing_path.c_str());
	
	//Point2D seed = { 248, 298 }; //---> P2 : automatic detection
	//Point2D seed = { 263, 257 }; //---> P2 : abdominal aorta
	//Point2D seed = { 263, 296 }; //---> P2
	//Point2D seed = { 263, 296 }; //---> P2
	//Point2D seed = { 262, 269 }; //---> P4
	Point2D seed = { 257, 244 }; //---> P2, liver
	Point2D found_point = { 0.0, 0.0 };

	Point2D found1 = { 292,314 }, found2 = { 253,255 }; // p2
	//imageSlice[x_new(292, 314, Length)] = 0.0;
	//imageSlice[x_new(253,255, Length)] = 0.0;

	////loading_path = inputPath + "threshold.raw";
	//loading_path = inputPath + "loaded.raw";
	//load2dArrayRAW<dataType>(imageSlice, Length, Width, loading_path.c_str(), false);
	
	//houghSpace[x_new(292, 314, Length)] = 1.0;
	//houghSpace[x_new(253,255, Length)] = 1.0;
	//houghSpace[x_new(264, 263, Length)] = 1.0;//p2,liver
	//imageSlice[x_new(264, 263, Length)] = 0.0;//p2,liver
	//imageSlice[x_new(292, 314, Length)] = 0.0;
	//imageSlice[x_new(253,255, Length)] = 0.0;
	//storing_path = outputPath + "input_plus_centers.raw";
	//store2dRawData<dataType>(imageSlice, Length, Width, storing_path.c_str());

	//string saving_csv = outputPath + "file_.csv";
	//FILE* file;
	//if (fopen_s(&file, saving_csv.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	
	found_point = localHoughTransform(seed, imageSlice, houghSpace, foundCircle, Length, Width, parameters, storing_path.c_str());
	cout << "found center : (" << found_point.x << "," << found_point.y << ")" << endl;

	//fclose(file);

	delete[] imageSlice;
	delete[] foundCircle;
	delete[] houghSpace;
	*/
	
	//======================== Segment lungs conneted to trachea ===================================================================
	
	/*
	////bounding box of soft tissue
	//size_t i_min = Length, i_max = 0, j_min = Width, j_max = 0, k_min = Height, k_max = 0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			xd = x_new(i, j, Length);
	//			
	//			if (ctContainer->dataPointer[k][xd] >= thres) {
	//				//find i_min
	//				if (i < i_min) {
	//					i_min = i;
	//				}
	//				//find j_min
	//				if (j < j_min) {
	//					j_min = j;
	//				}
	//				//find k_min
	//				if (k < k_min) {
	//					k_min = k;
	//				}
	//				//find i_max
	//				if (i > i_max) {
	//					i_max = i;
	//				}
	//				//find j_max
	//				if (j > j_max) {
	//					j_max = j;
	//				}
	//				//find k_max
	//				if (k > k_max) {
	//					k_max = k;
	//				}
	//			}
	//		}
	//	}
	//}

	//thresholding3dFunctionN(maskThreshold, Length, Width, Height, thres, thres, 0.0, 1.0);
	////Dilatation to close nose's wholes
	//dilatation3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 0.0, 1.0);

	////Dilatation to take back removed pixels during first dilatation
	//dilatation3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	//dilatation3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	//erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	
	////storing_path = outputPath + "threshold.raw";
	////store3dRawData<dataType>(maskThreshold, Length, Width, Height, storing_path.c_str());

	int** labelArray = new int* [Height];
	bool** status = new bool* [Height];
	for (k = 0; k < Height; k++) {
		labelArray[k] = new int[dim2D] {0};
		status[k] = new bool[dim2D];
	}
	if (labelArray == NULL || status == NULL)
		return false;

	////Initialization
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		status[k][i] = false;
	//	}
	//}

	////Erosion to remove trachea
	//erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	//labelling3D(maskThreshold, labelArray, status, Length, Width, Height, 1.0);

	////number of region cells
	//int numberOfRegionsCells = countNumberOfRegionsCells(maskThreshold, Length, Width, Height, 1.0);
	////cout << "Objects voxels number : " << numberOfRegionsCells << endl;

	////Counting
	//int* countingArray = new int[numberOfRegionsCells] {0};
	//if (countingArray == NULL)
	//	return false;

	////Get regions sizes
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		if (labelArray[k][i] > 0) {
	//			countingArray[labelArray[k][i]]++;
	//		}
	//	}
	//}

	////remove regions lying out of soft tissue bounding box
	//int label = 1;
	//size_t in, jn, kn;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (labelArray[k][x_new(i, j, Length)] == label) {
	//				if (i < i_min || i > i_max || j < j_min || j > j_max || k < k_min || k > k_max) {
	//					for (kn = 0; kn < Height; kn++) {
	//						for (in = 0; in < Length; in++) {
	//							for (jn = 0; jn < Width; jn++) {
	//								if (labelArray[kn][x_new(in, jn, Length)] == label) {
	//									labelArray[kn][x_new(in, jn, Length)] = 0;
	//								}
	//							}
	//						}
	//					}
	//					label++;
	//				}
	//			}
	//		}
	//	}
	//}

	////Find the biggest region
	//int max1 = 0, max2 = 0, max3 = 0, max4 = 0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		if (countingArray[labelArray[k][i]] > max1) {
	//			max1 = countingArray[labelArray[k][i]];
	//		}
	//		if (countingArray[labelArray[k][i]] > max2 && countingArray[labelArray[k][i]] < max1) {
	//			max2 = countingArray[labelArray[k][i]];
	//		}
	//		if (countingArray[labelArray[k][i]] > max3 && countingArray[labelArray[k][i]] < max2) {
	//			max3 = countingArray[labelArray[k][i]];
	//		}
	//		if (countingArray[labelArray[k][i]] > max4 && countingArray[labelArray[k][i]] < max3) {
	//			max4 = countingArray[labelArray[k][i]];
	//		}
	//	}
	//}
	////cout << "max1 = " << max1 << ", max2 = " << max2 << ", max3 = " << max3 << ", max4 = " << max4 << endl;

	////Keep biggest region
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		if (countingArray[labelArray[k][i]] == max1 || countingArray[labelArray[k][i]] == max2) {
	//			maskThreshold[k][i] = 1.0;
	//		}
	//		else {
	//			maskThreshold[k][i] = 0.0;
	//		}
	//	}
	//}

	////save labeled image
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		imageData[k][i] = (dataType)labelArray[k][i];
	//		//if (labelArray[k][i] != 0) {
	//		//	maskThreshold[k][i] = 1.0;
	//		//}
	//		//else {
	//		//	maskThreshold[k][i] = 0.0;
	//		//}
	//	}
	//}
	//storing_path = outputPath + "segment.raw";
	//store3dRawData<dataType>(imageData, Length, Width, Height, storing_path.c_str());

	////Dilatation to take lungs pixels back
	//dilatation3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);

	//storing_path = outputPath + "lungs.raw";
	//store3dRawData<dataType>(maskThreshold, Length, Width, Height, storing_path.c_str());

	//delete[] countingArray;
	for (k = 0; k < Height; k++) {
		delete[] labelArray[k];
		delete[] status[k];
	}
	delete[] labelArray;
	delete[] status;
	*/

	//======================== Segment Trachea =================
	
	/*
	dataType** lungsAndTrachea = new dataType * [Height];
	dataType** trachea = new dataType * [Height];
	int** labelArray = new int* [Height];
	bool** status = new bool* [Height];
	for (k = 0; k < Height; k++) {
		lungsAndTrachea[k] = new dataType[dim2D]{ 0 };
		trachea[k] = new dataType[dim2D]{ 0 };
		labelArray[k] = new int[dim2D] {0};
		status[k] = new bool[dim2D];
	}

	loading_path = inputPath + "lungs_and_trachea.raw";
	load3dArrayRAW<dataType>(lungsAndTrachea, Length, Width, Height, loading_path.c_str(), false);

	loading_path = inputPath + "lungs.raw";
	load3dArrayRAW<dataType>(imageData, Length, Width, Height, loading_path.c_str(), false);

	//Difference
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (imageData[k][i] == 1.0) {
				trachea[k][i] = 0.0;
			}
			else {
				trachea[k][i] = lungsAndTrachea[k][i];
			}
		}
	}

	erosion3dHeighteenNeigbours(trachea, Length, Width, Height, 1.0, 0.0);

	//storing_path = outputPath + "difference.raw";
	//store3dRawData<dataType>(trachea, Length, Width, Height, storing_path.c_str());

	//number of region cells
	int numberOfRegionsCells = countNumberOfRegionsCells(trachea, Length, Width, Height, 1.0);

	labelling3D(trachea, labelArray, status, Length, Width, Height, 1.0);

	//Counting
	int* countingArray = new int[numberOfRegionsCells] {0};
	if (countingArray == NULL)
		return false;

	//Get region size
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (labelArray[k][i] > 0) {
				countingArray[labelArray[k][i]]++;
			}
		}
	}

	//Find the biggest region
	int max1 = 0, max2 = 0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (countingArray[labelArray[k][i]] > max1) {
				max1 = countingArray[labelArray[k][i]];
			}
			if (countingArray[labelArray[k][i]] > max2 && countingArray[labelArray[k][i]] < max1) {
				max2 = countingArray[labelArray[k][i]];
			}
		}
	}

	//Keep biggest region
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (countingArray[labelArray[k][i]] == max1) {
				trachea[k][i] = 1.0;
			}
			else {
				trachea[k][i] = 0.0;
			}
		}
	}

	//dilatation to fill trachea
	dilatation3dHeighteenNeigbours(trachea, Length, Width, Height, 1.0, 0.0);

	storing_path = outputPath + "difference.raw";
	store3dRawData<dataType>(trachea, Length, Width, Height, storing_path.c_str());

	for (k = 0; k < Height; k++) {
		delete[] lungsAndTrachea[k];
		delete[] trachea[k];
		delete[] labelArray[k];
		delete[] status[k];
	}
	delete[] lungsAndTrachea;
	delete[] trachea;
	delete[] labelArray;
	delete[] status;
	*/
	
	//======================== Detect carina =======================
	
	/*
	loading_path = inputPath + "trachea.raw";
	load3dArrayRAW<dataType>(imageData, Length, Width, Height, loading_path.c_str(), false);
	
	int* labelArray = new int[dim2D];
	bool* status = new bool[dim2D];
	dataType* mask = new dataType[dim2D];
	if (labelArray == NULL || status == NULL || mask == NULL)
		return false;

	size_t n, count_slice = 0, slice_of_interest = 0;
	bool previous_slice = false;

	storing_path = outputPath + "info_regions.csv";
	FILE* file;
	if (fopen_s(&file, storing_path.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	fprintf(file, "slice number,number of regions\n");

	//size_t k_previous = 0, k_current = 0, k_next = 0;

	for (k = 0; k < Height; k++) {

		//if (k == 0) {
		//	k_previous = 0;
		//	k_current = k_previous;
		//	k_next = 1;
		//}
		//else {
		//	if (k == Height - 1) {
		//		k_previous = Height - 1;
		//		k_current = Height - 1;
		//		k_next = k_current;
		//	}
		//	else {
		//		k_previous = k - 1;
		//		k_current = k;
		//		k_next = k + 1;
		//	}
		//}

		//cout << "Slice " << k << ":" << endl;
		fprintf(file, "%d,", k);
		
		//Initialization
		for (i = 0; i < dim2D; i++) {
			mask[i] = imageData[k][i];
			labelArray[i] = 0;
			status[i] = false;
		}

		labeling(mask, labelArray, status, Length, Width, 1.0);

		//number of region cells
		int numberOfRegionsCells = 0;
		for (i = 0; i < dim2D; i++) {
			if (mask[i] == 1.0) {
				numberOfRegionsCells++;
			}
		}

		if (numberOfRegionsCells != 0) {


			int* countingArray = new int[numberOfRegionsCells] {0};

			//Get regions sizes
			for (i = 0; i < dim2D; i++) {
				if (labelArray[i] > 0) {
					countingArray[labelArray[i]]++;
				}
			}

			//Get number of regions
			int numb_region = 0, label = 1;
			for (i = 0; i < dim2D; i++) {
				if (labelArray[i] == label) {
					numb_region++;
					label++;
				}
			}
			
			////Print number of elements per region
			//for (i = 1; i <= numb_region; i++) {
			//	std::cout << "region " << i << " has " << countingArray[i] << " elements" << std::endl;
			//}
			//std::cout << "############" << std::endl;

			//std::cout << "Slice " << k  << " has " << numb_region << " region(s)" << std::endl;
			fprintf(file, "%d\n", numb_region);

			if (numb_region > 1) {
				count_slice++;
				previous_slice = true;
			}
			else {
				count_slice = 0;
				previous_slice = false;
			}

			if (count_slice == 5 && previous_slice == true) {
				slice_of_interest = k - 5;
				count_slice = 0;
				std::cout << "The slice of interest is : " << slice_of_interest << std::endl;
			}

			delete[] countingArray;

		}
		else {
			fprintf(file, "%d\n", numberOfRegionsCells);
			previous_slice = false;
			count_slice = 0;
		}
		//std::cout << "count slice " << count_slice << std::endl;

	}

	fclose(file);
	delete[] labelArray;
	delete[] status;
	*/

	//======================== Detect seed from carina ==============
	
	/*
	loading_path = inputPath + "trachea.raw";
	load3dArrayRAW<dataType>(imageData, Length, Width, Height, loading_path.c_str(), false);

	size_t carina_slice = 257; // patient2

	dataType** distanceMap = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		distanceMap[k] = new dataType[dim2D]{ 0 };
	}

	fastSweepingFunction_3D(distanceMap, imageData, Length, Width, Height, 1.0, 1000000.0, 0.0);

	dataType max_distance = 0.0;
	Point2D point_max = { 0.0, 0.0 };
	for (i = 0; i < Length; i++) {
		for (j = 0; j < Width; j++) {
			if (distanceMap[carina_slice][x_new(i, j, Length)] > max_distance) {
				point_max.x = (dataType)i;
				point_max.y = (dataType)j;
			}
		}
	}

	std::cout << "seed : (" << point_max.x << "," << point_max.y << ")" << std::endl;

	BoundingBox2D box = findBoundingBox2D(point_max, Length, Width, 8.5, 70);

	dataType* imageSlice = new dataType[dim2D]{ 0 };
	for (i = 0; i < dim2D; i++) {
		imageSlice[i] = maskThreshold[carina_slice][i];
	}

	//draw bounding box around seed in carina
	for (i = box.i_min; i <= box.i_max; i++) {
		for (j = box.j_min; j <= box.j_max; j++) {
			if (i == box.i_min || i == box.i_max || j == box.j_min || j == box.j_max) {
				imageSlice[x_new(i, j, Length)] = min_data;
			}
		}
	}

	storing_path = outputPath + "slice_of_interest.raw";
	store2dRawData<dataType>(imageSlice, Length, Width, storing_path.c_str());

	delete[] imageSlice;
	for (k = 0; k < Height; k++) {
		delete[] distanceMap[k];
	}
	delete[] distanceMap;
	*/
	
	//======================== Slice Cropping ================================================================
	
	/*
	imageMetaData croppedLungs;
	const size_t offset = 5;
	croppedLungs =  croppImage3D(ctImageData, offset);

	size_t length = croppedLungs.length, width = croppedLungs.width, height = croppedLungs.height;
	std::cout << "Cropped lungs image dim : " << length << " x " << width << " x " << height << "" << std::endl;
	Point3D originLungs = croppedLungs.origin;
	std::cout << "Cropped Lungs origin : (" << originLungs.x << ", " << originLungs.y << ", " << originLungs.z << ")" << std::endl;
	VoxelSpacing lungsSpacing = croppedLungs.spacing;

	dataType** cropped_ct = new dataType*[height];
	dataType** cropped_lungs = new dataType * [height];
	for (k = 0; k < height; k++) {
		cropped_ct[k] = new dataType[dim2D]{0};
		cropped_lungs[k] = new dataType[dim2D]{0};
		if (cropped_ct[k] == NULL || cropped_lungs[k] == NULL)
			return false;
	}
	if (cropped_ct == NULL || cropped_lungs == NULL)
		return false;

	//===
	//Image_Data Cropped;
	//Cropped.imageDataPtr = croppedCT;
	//Cropped.height = height; Cropped.length = length; Cropped.width = width;
	//Cropped.orientation = orientation; Cropped.origin = croppedOrigin; Cropped.spacing = croppedSpacing;
	//===

	Point3D origin_cropped = getImageCoordFromRealCoord3D(originLungs, ctOrigin, ctSpacing, orientation);
	size_t i0 = (size_t)origin_cropped.x, j0 = (size_t)origin_cropped.y, k0 = (size_t)origin_cropped.z;
	
	Point3D point_test = { 0.0, 0.0, k0 };
	point_test = getRealCoordFromImageCoord3D(point_test, ctOrigin, ctSpacing, orientation);
	cout << "Origin cropped image : (" << point_test.x << "," << point_test.y << "," << point_test.z << ")" << endl;
	
	size_t i1, j1, k1;
	for (k = 0, k1 = k0; k < height; k++, k1++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				cropped_ct[k][x_new(i, j, Length)] = imageData[k1][x_new(i, j, Length)];
				cropped_lungs[k][x_new(i, j, Length)] = lungsData[k1][x_new(i, j, Length)];
			}
		}
	}

	saving_name = outputPath + "cropped_ct_p6.raw";
	//store3dRawData<dataType>(cropped_ct, Length, Width, height, saving_name.c_str());

	saving_name = outputPath + "cropped_lungs_p6.raw";
	//store3dRawData<dataType>(cropped_lungs, Length, Width, height, saving_name.c_str());

	Vtk_File_Info* Container = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	Container->dataPointer = cropped_lungs; //cropped_ct;
	Container->dimensions[0] = Length; Container->dimensions[1] = Width; Container->dimensions[2] = height;
	Container->origin[0] = point_test.x; Container->origin[1] = point_test.y; Container->origin[2] = point_test.z;
	Container->spacing[0] = ctSpacing.sx; Container->spacing[1] = ctSpacing.sy; Container->spacing[2] = ctSpacing.sz;
	Container->operation = copyTo; Container->vDataType = dta_Flt;

	saving_name = outputPath + "cropped_lungs_p6.vtk";
	vtkDataForm data_form = dta_binary;
	storeVtkFile(saving_name.c_str(), Container, data_form);

	free(Container);
	for (k = 0; k < height; k++) {
		delete[] cropped_ct[k];
		delete[] cropped_lungs[k];
	}
	delete[] cropped_ct;
	delete[] cropped_lungs;
	*/

	//======================== Artificial 4D fast marching and path finding ================
	
	/*
	//create 3D artificial image
	const size_t length = 50, width = 50, height = 50;
	const size_t dim2D = length * width;

	dataType** imageData = new dataType * [height];
	dataType** path3D = new dataType * [height];
	dataType** potential3D = new dataType * [height];
	dataType** action3D = new dataType * [height];
	for (k = 0; k < height; k++) {
		imageData[k] = new dataType[dim2D]{ 0 };
		path3D[k] = new dataType[dim2D]{ 0 };
		potential3D[k] = new dataType[dim2D]{ 0 };
		action3D[k] = new dataType[dim2D]{ 0 };
	}

	//size_t step = 5;
	//double distance = 0.0;

	////create the image version 1
	//for (k = 0; k < height; k++) {
	//	Point3D seed = { 25, 25, k };
	//	if (k <= 5) {
	//		double coef = 1.0;
	//		double radius = step * coef;
	//		for (i = 0; i < length; i++) {
	//			for (j = 0; j < width; j++) {
	//				Point3D current = { i, j, k };
	//				distance = getPoint3DDistance(current, seed);
	//				if (distance <= radius) {
	//					imageData[k][x_new(i, j, length)] = 1.0;
	//				}
	//			}
	//		}
	//	}
	//	if (k > 5 && k <= 10) {
	//		double coef = 1.1;
	//		double radius = step * coef;
	//		for (i = 0; i < length; i++) {
	//			for (j = 0; j < width; j++) {
	//				Point3D current = { i, j, k };
	//				distance = getPoint3DDistance(current, seed);
	//				if (distance <= radius) {
	//					imageData[k][x_new(i, j, length)] = 1.0;
	//				}
	//			}
	//		}
	//	}
	//	if (k > 10 && k <= 15) {
	//		double coef = 1.2;
	//		double radius = step * coef;
	//		for (i = 0; i < length; i++) {
	//			for (j = 0; j < width; j++) {
	//				Point3D current = { i, j, k };
	//				distance = getPoint3DDistance(current, seed);
	//				if (distance <= radius) {
	//					imageData[k][x_new(i, j, length)] = 1.0;
	//				}
	//			}
	//		}
	//	}
	//	if (k > 15 && k <= 20) {
	//		double coef = 1.3;
	//		double radius = step * coef;
	//		for (i = 0; i < length; i++) {
	//			for (j = 0; j < width; j++) {
	//				Point3D current = { i, j, k };
	//				distance = getPoint3DDistance(current, seed);
	//				if (distance <= radius) {
	//					imageData[k][x_new(i, j, length)] = 1.0;
	//				}
	//			}
	//		}
	//	}
	//	if (k > 20 && k <= 25) {
	//		double coef = 1.4;
	//		double radius = step * coef;
	//		for (i = 0; i < length; i++) {
	//			for (j = 0; j < width; j++) {
	//				Point3D current = { i, j, k };
	//				distance = getPoint3DDistance(current, seed);
	//				if (distance <= radius) {
	//					imageData[k][x_new(i, j, length)] = 1.0;
	//				}
	//			}
	//		}
	//	}
	//	if (k > 25 && k <= 30) {
	//		double coef = 1.5;
	//		double radius = step * coef;
	//		for (i = 0; i < length; i++) {
	//			for (j = 0; j < width; j++) {
	//				Point3D current = { i, j, k };
	//				distance = getPoint3DDistance(current, seed);
	//				if (distance <= radius) {
	//					imageData[k][x_new(i, j, length)] = 1.0;
	//				}
	//			}
	//		}
	//	}
	//	if (k > 30 && k <= 35) {
	//		double coef = 1.6;
	//		double radius = step * coef;
	//		for (i = 0; i < length; i++) {
	//			for (j = 0; j < width; j++) {
	//				Point3D current = { i, j, k };
	//				distance = getPoint3DDistance(current, seed);
	//				if (distance <= radius) {
	//					imageData[k][x_new(i, j, length)] = 1.0;
	//				}
	//			}
	//		}
	//	}
	//	if (k > 35 && k <= 40) {
	//		double coef = 1.7;
	//		double radius = step * coef;
	//		for (i = 0; i < length; i++) {
	//			for (j = 0; j < width; j++) {
	//				Point3D current = { i, j, k };
	//				distance = getPoint3DDistance(current, seed);
	//				if (distance <= radius) {
	//					imageData[k][x_new(i, j, length)] = 1.0;
	//				}
	//			}
	//		}
	//	}
	//	if (k > 40 && k <= 45) {
	//		double coef = 1.8;
	//		double radius = step * coef;
	//		for (i = 0; i < length; i++) {
	//			for (j = 0; j < width; j++) {
	//				Point3D current = { i, j, k };
	//				distance = getPoint3DDistance(current, seed);
	//				if (distance <= radius) {
	//					imageData[k][x_new(i, j, length)] = 1.0;
	//				}
	//			}
	//		}
	//	}
	//	if (k > 45) {
	//		double coef = 1.9;
	//		double radius = step * coef;
	//		for (i = 0; i < length; i++) {
	//			for (j = 0; j < width; j++) {
	//				Point3D current = { i, j, k };
	//				distance = getPoint3DDistance(current, seed);
	//				if (distance <= radius) {
	//					imageData[k][x_new(i, j, length)] = 1.0;
	//				}
	//			}
	//		}
	//	}
	//}

	//create the input image (version 2)
	double radius_image = 5.0;
	for (k = 0; k < height; k++) {
		Point3D seed = { 25, 25, k };
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {
				Point3D current = { i, j, k };
				double distance = getPoint3DDistance(seed, current);
				if (distance <= radius_image) {
					imageData[k][x_new(i, j, length)] = 1.0;
				}
			}
		}
	}

	storing_path = outputPath + "artificial.raw";
	store3dRawData<dataType>(imageData, length, width, height, storing_path.c_str());

	Image_Data ctImageData;
	ctImageData.imageDataPtr = imageData;
	ctImageData.height = height; ctImageData.width = width; ctImageData.length = length;
	ctImageData.origin.x = 0.0; ctImageData.origin.y = 0.0; ctImageData.origin.z = 0.0;
	ctImageData.spacing.sx = 1.0; ctImageData.spacing.sy = 1.0; ctImageData.spacing.sz = 0.0;
	ctImageData.orientation = orientation;

	//Point3D seed1 = { 25, 25, 0 }, seed2 = { 25, 25, 49 };
	//Point3D seed1 = { 28, 27, 0 }, seed2 = { 25, 25, 49 };
	//Point3D seed1 = { 28, 27, 0 }, seed2 = { 28, 17, 49 };
	//Point3D seed1 = { 28, 28, 0 }, seed2 = { 30, 18, 49 };
	//Point3D seed1 = { 27, 27, 0 }, seed2 = { 23, 22, 49 };
	Point3D seed1 = { 25, 25, 4 }, seed2 = { 23, 22, 49 };
	Point3D* seedPoints = new Point3D[2];
	seedPoints[0] = seed1; 
	seedPoints[1] = seed2;

	Potential_Parameters parameters;
	parameters.K = 0.005; parameters.epsilon = 0.01;
	parameters.c_ct = 0.0; parameters.c_pet = 0.0;
	parameters.c_max = 0.0; parameters.c_min = 0.0;
	parameters.c_mean = 1.0; parameters.c_sd = 0.0;
	double radius = 3.0;

	dataType raduisScale = 5.0, radiusInitial = 5.0, radiusFinal = 15.0;
	size_t radLength = 4;
	size_t newHeight = height * radLength;

	dataType** potential4D = new dataType * [newHeight];
	dataType** action4D = new dataType * [newHeight];
	for (k = 0; k < newHeight; k++) {
		potential4D[k] = new dataType[dim2D]{ 0 };
		action4D[k] = new dataType[dim2D]{ 0 };
	}

	//constant potential 4D
	for (k = 0; k < newHeight; k++) {
		for (i = 0; i < dim2D; i++) {
			potential4D[k][i] = 1.0;
		}
	}

	//computePotentialNew(ctImageData, imageData, potential3D, seed1, radius, parameters);
	//computePotential4D(ctImageData, potentialPtr, seed1, radLength, rad1, rad2);
	//storing_path = outputPath + "potential.raw";
	//store3dRawData<dataType>(potentialPtr, length, width, height, storing_path.c_str());
	//constant potential 3D
	//for (k = 0; k < height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		potential3D[k][i] = 1.0;
	//	}
	//}

	//fastMarching3D_N(action3D, potential3D, length, width, height, seed1);
	//computeDistanceToOnePoint(action3D, length, width, height, seed1);
	fastMarching4D(action4D, potential4D, length, width, height, radLength, seed1, raduisScale, radiusInitial);
	//storing_path = outputPath + "action3D_euclid.raw";
	//store3dRawData<dataType>(action3D, length, width, height, storing_path.c_str());

	//storing_path = outputPath + "dist_eucl.raw";
	//store3dRawData<dataType>(actionPtr, length, width, newHeight, storing_path.c_str());

	dataType h = 1.0;
	//vector<Point3D> path_points;
	//shortestPath3D(action3D, path3D, length, width, height, h, seedPoints, path_points);
	//storing_path = outputPath + "path.raw";
	//store3dRawData<dataType>(pathPtr, length, width, height, storing_path.c_str());
	//shortestPath3DPlus(actionPtr, length, width, height, h, seedPoints, radLength, raduisScale, radiusInitial, radiusFinal);
	shortestPath4D(action4D, length, width, height, h, seedPoints, radLength, raduisScale, radiusInitial, radiusFinal);
	
	//string saving_csv = outputPath + "path4D.csv";
	//FILE* file;
	//if (fopen_s(&file, saving_csv.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//for (i = 0; i < path_points.size(); i++) {
	//	fprintf(file, "%f,%f,%f\n", path_points[i].x, path_points[i].y, path_points[i].z);
	//}
	//fclose(file);

	delete[] seedPoints;
	for (k = 0; k < height; k++) {
		delete[] imageData[k];
		delete[] path3D[k];
		delete[] potential3D[k];
		delete[] action3D[k];
	}
	delete[] imageData;
	delete[] path3D;
	delete[] potential3D;
	delete[] action3D;

	for (k = 0; k < newHeight; k++) {
		delete[] potential4D[k];
		delete[] action4D[k];
	}
	delete[] potential4D;
	delete[] action4D;
	*/
	
	//======================== 4D fast marching and path finding ================
	
	/*
	const size_t min_radius = 7, max_radius = 8, initial_radius = 7, final_radius = 7; //  max_radius = 20
	size_t lenght_radius = max_radius - min_radius + 1;
	std::cout << "number of radius : " << lenght_radius << std::endl;

	const size_t height = lenght_radius * height_int;
	dataType** potential = new dataType * [height];
	dataType** action_map = new dataType * [height];
	for (k = 0; k < height; k++) {
		potential[k] = new dataType[dim2D]{ 0 };
		action_map[k] = new dataType[dim2D]{ 0 };
	}

	//for (k = 0; k < height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		potential[k][i] = 1.0;
	//	}
	//}

	Point3D seed1 = { 263, 257, 146 }; //p2
	seed1 = getRealCoordFromImageCoord3D(seed1, ctOrigin, ctSpacing, orientation);
	seed1 = getImageCoordFromRealCoord3D(seed1, ctOrigin, INTERPOLATE.spacing, orientation);
	
	Point3D seed2 = { 294, 315, 261 };
	seed2 = getRealCoordFromImageCoord3D(seed2, ctOrigin, ctSpacing, orientation);
	seed2 = getImageCoordFromRealCoord3D(seed2, ctOrigin, INTERPOLATE.spacing, orientation);

	Point3D* seedPoints = new Point3D[2];
	seedPoints[0] = seed1;
	seedPoints[1] = seed2;

	//double start_time = clock();
	//computePotential4D(INTERPOLATE, potential, seed1, initial_radius, min_radius, max_radius);
	//double end_time = clock();
	//std::cout << "The computation lasts : " << (end_time - start_time) / CLOCKS_PER_SEC << std::endl;
	//storing_path = outputPath + "potential.raw";
	//store3dRawData<dataType>(potential, Length, Width, height, storing_path.c_str());

	loading_path = inputPath + "potential.raw";
	load3dArrayRAW<dataType>(potential, Length, Width, height, loading_path.c_str(), false);

	fastMarching4D(action_map, potential, Length, Width, height_int, seed1, min_radius, max_radius, initial_radius);
	storing_path = outputPath + "action_map.raw";
	store3dRawData<dataType>(action_map, Length, Width, height, storing_path.c_str());


	shortestPath4D(action_map, Length, Width, height_int, 1.0, seedPoints, min_radius, max_radius, initial_radius, final_radius);

	delete[] seedPoints;
	for (k = 0; k < height; k++) {
		delete[] action_map[k];
		delete[] potential[k];
	}
	delete[] action_map;
	delete[] potential;
	*/

	//======================== Free Memory ======================================

	
	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
		delete[] maskThreshold[k];
	}
	delete[] imageData;
	delete[] maskThreshold;
	
	for (k = 0; k < height; k++) {
		delete[] interpolated[k];
	}
	delete[] interpolated;
	free(ctContainer);
	

	return EXIT_SUCCESS;
}