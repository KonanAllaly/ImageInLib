#include <iostream>
#include <sstream>  
#include <vector>
#include <string.h> 
#include "common_math.h"

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

#define thres_min -30 // 995
#define thres_max 200// 199 //1213
//#define thres -500

int main() {

	string inputPath = "C:/Users/Konan Allaly/Documents/Tests/input/";
	string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";
	
	string loading_path, storing_path, extension;

	size_t i = 0, j = 0, k = 0, xd = 0;

	//===================== Load 3D patient data (.vtk) ================================
	
	
	/*
	We load the .vtk files containing the original pixels value
	and informations related the dimensions, the spacings and
	the origins. We also define the orientation matrix in 2D/3D
	needed when we need to perform interpolation.
	*/
	
	/*
	OrientationMatrix orientation;
	orientation.v1 = { 1.0, 0.0, 0.0 }; 
	orientation.v2 = { 0.0, 1.0, 0.0 }; 
	orientation.v3 = { 0.0, 0.0, 1.0 };

	//OrientationMatrix2D orientation2D;
	//orientation2D.v1 = { 1.0, 0.0 };
	//orientation2D.v2 = { 0.0, 1.0 };

	Vtk_File_Info* ctContainer = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	ctContainer->operation = copyFrom;
	
	loading_path = inputPath + "vtk/ct/Patient2_ct.vtk";
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
	*/

	//========================= Path finding in Spiral ==========================

	/*
	const size_t height = 100, length = 100, width = 100;
	dataType** imageData = new dataType * [height];
	dataType** actionMap = new dataType * [height];
	dataType** potential = new dataType * [height];
	for (k = 0; k < height; k++) {
		imageData[k] = new dataType[length * width]{ 0 };
		actionMap[k] = new dataType[length * width]{ 0 };
		potential[k] = new dataType[length * width]{ 0 };
	}

	//loading_path = inputPath + "shape/diagonal_tube.raw";
	loading_path = inputPath + "shape/u_shape_radius15.raw";
	//loading_path = inputPath + "shape/spiral_radius15.raw";
	load3dArrayRAW<dataType>(imageData, length, width, height, loading_path.c_str(), false);

	//Compute potential
	for (k = 0; k < height; k++) {
		for (i = 0; i < length * width; i++) {
			if (imageData[k][i] == 1.0) {
				potential[k][i] = 0.001;
			}
			else {
				potential[k][i] = 1.0;
			}
		}
	}

	Point3D* seedPoints = new Point3D[2];
	//Point3D pt1 = { 7, 9, 5 }, pt2 = { 94, 90, 18 }; //tube
	Point3D pt1 = { 30, 18, 46 }, pt2 = { 72, 21, 51 }; //u shape
	//Point3D pt1 = { 81, 54, 7 }, pt2 = { 80, 48, 87 }; //spiral

	seedPoints[0] = pt1;
	seedPoints[1] = pt2;

	//fastMarching3D_N(actionMap, potential, length, width, height, pt1);
	partialFrontPropagation(actionMap, potential, length, width, height, seedPoints);

	vector<Point3D> path_points;
	VoxelSpacing spacing = { 1.0, 1.0, 1.0 };
	shortestPath3D(actionMap, length, width, height, spacing, seedPoints, path_points);

	//string path_save = outputPath + "path_spiral_radius15.csv";
	string path_save = outputPath + "path_ushape_radius15.csv";
	FILE* pFile;
	if (fopen_s(&pFile, path_save.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	int n = 0;

	fprintf(pFile, "%f,%f,%f\n", pt1.x, pt1.y, pt1.z);
	for (n = path_points.size() - 1; n > -1; n--) {
		fprintf(pFile, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	}
	fprintf(pFile, "%f,%f,%f\n", pt2.x, pt2.y, pt2.z);
	fclose(pFile);

	delete[] seedPoints;
	for (k = 0; k < height; k++) {
		delete[] imageData[k];
		delete[] actionMap[k];
		delete[] potential[k];
	}
	delete[] imageData;
	delete[] actionMap;
	delete[] potential;
	*/

	//========================= Interpolated the aorta ==========================

	/*
	dataType** imageDataPtr = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageDataPtr[k] = new dataType[dim2D]{ 0 };
	}

	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (ctContainer->dataPointer[k][i] != 0) {
				imageDataPtr[k][i] = 1.0;
			}
			else {
				imageDataPtr[k][i] = 0.0;
			}
		}
	}

	Image_Data original_aorta = {
		Height,
		Length,
		Width,
		imageDataPtr,
		ctOrigin,
		ctSpacing,
		orientation
	};

	//loading_path = inputPath + "raw/aorta/aorta_p2.raw";
	//load3dArrayRAW<dataType>(imageDataPtr, Length, Width, Height, loading_path.c_str(), false);

	const size_t height = 1083;//844;//1351;//866;

	dataType** interpolImageDataPtr = new dataType * [height];
	for (k = 0; k < height; k++) {
		interpolImageDataPtr[k] = new dataType[dim2D]{ 0 };
	}

	VoxelSpacing intSpacing = { ctSpacing.sx, ctSpacing.sy, ctSpacing.sx };
	Image_Data interpolated_aorta = {
		height,
		Length,
		Width,
		interpolImageDataPtr,
		ctOrigin,
		intSpacing,
		orientation
	};

	imageInterpolation3D(original_aorta, interpolated_aorta, NEAREST_NEIGHBOR);

	storing_path = inputPath + "raw/interpolated/aorta_p3.raw";
	store3dRawData<dataType>(interpolImageDataPtr, Length, Width, height, storing_path.c_str());

	for (k = 0; k < Height; k++) {
		delete[] imageDataPtr[k];
	}
	delete[] imageDataPtr;

	for (k = 0; k < height; k++) {
		delete[] interpolImageDataPtr[k];
	}
	delete[] interpolImageDataPtr;
	*/

	//========================= Generate velocity field =========================

	/*
	const size_t height = 866;
	dataType** imageDataPtr = new dataType * [height];
	dataType** pathImageData = new dataType * [height];
	dataType** edgeDetector = new dataType * [height];
	for (k = 0; k < height; k++) {
		imageDataPtr[k] = new dataType[dim2D]{ 0 };
		pathImageData[k] = new dataType[dim2D]{ 0 };
		edgeDetector[k] = new dataType[dim2D]{ 0 };
	}

	VoxelSpacing intSpacing = {ctSpacing.sx, ctSpacing.sy, ctSpacing.sx};
	Image_Data CT = {
		height,
		Length,
		Width,
		imageDataPtr,
		ctOrigin,
		intSpacing,
		orientation
	};

	loading_path = inputPath + "raw/interpolated/filtered_p2.raw";
	load3dArrayRAW<dataType>(imageDataPtr, Length, Width, height, loading_path.c_str(), false);

	FILE* path_file;
	loading_path = inputPath + "csv/Path/path_ordered_p2.csv";
	if (fopen_s(&path_file, loading_path.c_str(), "r") != 0) {
		printf("Enable to open");
		return false;
	}
	dataType x = 0, y = 0, z = 0;
	vector<Point3D> path_points;

	while (feof(path_file) == 0) {
		fscanf_s(path_file, "%f", &x);
		fscanf_s(path_file, ",");
		fscanf_s(path_file, "%f", &y);
		fscanf_s(path_file, ",");
		fscanf_s(path_file, "%f", &z);
		fscanf_s(path_file, "\n");
		Point3D current_point = { x, y, z };
		current_point = getImageCoordFromRealCoord3D(current_point, CT.origin, CT.spacing, CT.orientation);
		path_points.push_back(current_point);
	}
	fclose(path_file);

	BoundingBox3D box;
	double radius = 25;

	for (size_t n = 0; n < path_points.size(); n++) {
		box = findBoundingBox3D(path_points[n], Length, Width, height, radius, 0.0);
		for (k = box.k_min; k <= box.k_max; k++) {
			for (i = box.i_min; i <= box.i_max; i++) {
				for (j = box.j_min; j <= box.j_max; j++) {
					xd = x_new(i, j, Length);
					pathImageData[k][xd] = imageDataPtr[k][xd];
				}
			}
		}
	}

	////compute3dImageGradient(imageDataPtr, gradientVectorX, gradientVectorY, gradientVectorZ, Length, Width, Height, ctSpacing);
	////fastSweepingFunction_3D(distance, maskThreshold, length, width, height, 1.0, 10000000.0, 0.0);
	//Point3D grad;
	//dataType norm_grad = 0.0, K = 1000.0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			//getGradient3D(imageDataPtr, Width, Length, Height, j, i, k, ctSpacing, &grad);
	//			getGradient3D(pathImageData, Width, Length, Height, j, i, k, ctSpacing, &grad);
	//			norm_grad = sqrt(grad.x * grad.x + grad.y * grad.y + grad.z * grad.z);
	//			edgeDetector[k][x_new(j, i, Width)] = gradientFunction(norm_grad, K);
	//		}
	//	}
	//}

	//thresholding3dFunctionN(edgeDetector, Length, Width, Height, 0.15, 0.15, 0.0, 1.0);
	storing_path = outputPath + "pathImage_p2.raw";
	store3dRawData<dataType>(pathImageData, Length, Width, height, storing_path.c_str());

	for (k = 0; k < Height; k++) {
		delete[] imageDataPtr[k];
		delete[] pathImageData[k];
		delete[] edgeDetector[k];
	}
	delete[] imageDataPtr;
	delete[] pathImageData;
	delete[] edgeDetector;
	*/

	//========================= Test translation ================================

	/*
	Vtk_File_Info* afterTreatment = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	afterTreatment->operation = copyFrom;

	loading_path = inputPath + "vtk/ct/Patient3_ct.vtk";
	readVtkFile(loading_path.c_str(), afterTreatment);
	Point3D ctOrigin = { afterTreatment->origin[0], afterTreatment->origin[1], afterTreatment->origin[2] };
	VoxelSpacing ctSpacing = { afterTreatment->spacing[0], afterTreatment->spacing[1], afterTreatment->spacing[0] };

	Point3D seedP2 = { 5.169647, 71.255859, -656.131958 };
	Point3D seedP3 = { 257, 255, 522 };
	seedP3 = getRealCoordFromImageCoord3D(seedP3, ctOrigin, ctSpacing, orientation);

	Point3D translate = { seedP3.x - seedP2.x, seedP3.y - seedP2.y, seedP3.z - seedP2.z };

	//Load path point
	FILE* path_file;
	loading_path = inputPath + "inputsForPhDStudentConference/path_ordered_p2.csv";
	if (fopen_s(&path_file, loading_path.c_str(), "r") != 0) {
		printf("Enable to open");
		return false;
	}

	FILE* translated_path;
	loading_path = inputPath + "inputsForPhDStudentConference/translated.csv";
	if (fopen_s(&translated_path, loading_path.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}

	dataType x = 0, y = 0, z = 0;
	dataType xt = 0, yt = 0, zt = 0;
	size_t count_point = 0;

	while (feof(path_file) == 0) {
		fscanf_s(path_file, "%f", &x);
		fscanf_s(path_file, ",");
		fscanf_s(path_file, "%f", &y);
		fscanf_s(path_file, ",");
		fscanf_s(path_file, "%f", &z);
		fscanf_s(path_file, "\n");
		count_point++;

		xt = x + translate.x;
		yt = y + translate.y;
		zt = z + translate.z;

		fprintf(translated_path, "%f,%f,%f\n", xt, yt, zt);
		//if (count_point == 1) {
		//	std::cout << "The first point is : (" << x << "," << y << "," << z << ")" << endl;
		//}

	}
	fclose(path_file);
	fclose(translated_path);
	free(afterTreatment);
	*/

	//========================= Test function histogram =========================
	
	/*
	const size_t Length = 512, Width = 512, Height = 406;
	size_t dim2D = Length * Width;

	dataType** imageDataPtr = new dataType * [Height];
	dataType** histogramPtr = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageDataPtr[k] = new dataType[dim2D]{ 0 };
		histogramPtr[k] = new dataType[dim2D]{ 0 };
	}

	loading_path = inputPath + "raw/ct/patient2_ct.raw";
	load3dArrayRAW<dataType>(imageDataPtr, Length, Width, Height, loading_path.c_str(), false);

	//const size_t binsCount = 6000;
	//computeHistogram(imageDataPtr, histogramPtr, Length, Width, Height, binsCount);
	//thresholdingOTSU(imageDataPtr, Length, Width, Height, 0.0, 1.0);
	//thresholdingOTSUNew(imageDataPtr, Length, Width, Height, binsCount, 0.0, 1.0);
	//storing_path = outputPath + "thresholdOtsu.raw";
	//store3dRawData<dataType>(imageDataPtr, Length, Width, Height, storing_path.c_str());

	//storing_path = outputPath + "histogram.raw";
	//store3dRawData<dataType>(histogramPtr, Length, Width, Height, storing_path.c_str());

	//2D OTSU method
	//const size_t Length = 260, Width = 260;
	//const size_t dim2D = Length * Width;

	dataType* image2D = new dataType[dim2D]{ 0 };
	dataType* threshold = new dataType[dim2D]{ 0 };

	//2D image extraction
	const size_t slice_number = 210;
	for (i = 0; i < dim2D; i++) {
		image2D[i] = imageDataPtr[slice_number][i];
	}

	//loading_path = inputPath + "raw/slice/memboa.raw";
	//loading_path = inputPath + "raw/slice/eye_1024_1024.raw";
	//load2dArrayRAW<dataType>(image2D, Length, Width, loading_path.c_str(), false);

	//storing_path = outputPath + "loaded.raw";
	//store2dRawData<dataType>(image2D, Length, Width, storing_path.c_str());

	//Point2D seed = {382, 454}; // eye image
	//Point2D seed = { 154, 77 }; // memboa
	//Point2D seed = { 188, 282 };
	Point2D seed1 = { 187, 283 };
	Point2D seed2 = { 249, 216 };
	double radius = 25;

	////thresholdingOTSU2D(image2D, 260, 260, 1.0, 0.0);
	//dataType thres1 = localOTSU2D(image2D, Length, Width, seed1, radius, 1.0, 0.0);
	//dataType thres2 = localOTSU2D(image2D, Length, Width, seed2, radius, 1.0, 0.0);

	//dataType thres = 0.0;
	//if (thres1 > thres2) {
	//	thres = thres1;
	//}
	//else {
	//	thres = thres2;
	//}
	//
	//for (i = 0; i < dim2D; i++) {
	//	if (image2D[i] <= thres) {
	//		threshold[i] = 0.0;
	//	}
	//	else {
	//		threshold[i] = 1.0;
	//	}
	//}

	double offset = 3;
	BoundingBox2D box = findBoundingBox2D(seed1, Length, Width, radius, offset);
	size_t i_min = box.i_min, i_max = box.i_max;
	size_t j_min = box.j_min, j_max = box.j_max;

	//Find the minimal and the maximal pixel value
	dataType min_data = 1000000, max_data = -1000000;
	for (i = i_min; i <= i_max; i++) {
		for (j = j_min; j <= j_max; j++) {
			xd = x_new(i, j, Length);
			if (image2D[xd] < min_data) {
				min_data = image2D[xd];
			}
			if (image2D[xd] > max_data) {
				max_data = image2D[xd];
			}
		}
	}
	cout << "min data = " << min_data << ", max data = " << max_data << endl;

	for (i = 0; i < dim2D; i++) {
		if (image2D[i] >= min_data && image2D[i] <= max_data) {
			threshold[i] = 1.0;
		}
		else {
			threshold[i] = 0.0;
		}
	}

	storing_path = outputPath + "threshold_box.raw";
	store2dRawData<dataType>(threshold, Length, Width, storing_path.c_str());

	delete[] image2D;
	delete[] threshold;

	for (k = 0; k < Height; k++) {
		delete[] imageDataPtr[k];
		delete[] histogramPtr[k];
	}
	delete[] imageDataPtr;
	delete[] histogramPtr;
	*/

	//========================= Load segmented liver and ajust dimension ========

	/*
	////Patient 2 : interpolated
	//size_t i_min = 120, i_max = 375;
	//size_t j_min = 160, j_max = 415;
	//size_t k_min = 335, k_max = 494;

	//Patient 3 : interpolated
	size_t i_min = 90, i_max = 345;
	size_t j_min = 130, j_max = 385;
	size_t k_min = 532, k_max = 711;

	Point3D LiverOrigin = { i_min, j_min, k_min };
	dataType sx = 4.6875;
	
	VoxelSpacing interpolSpacing = { ctSpacing.sx, ctSpacing.sy, ctSpacing.sx };
	LiverOrigin = getRealCoordFromImageCoord3D(LiverOrigin, ctOrigin, interpolSpacing, orientation);
	//const size_t length = 64, width = 64, height = 40; //P2
	const size_t length = 64, width = 64, height = 45; //P3

	dataType** liverShape = new dataType * [height];
	for (k = 0; k < height; k++) {
		liverShape[k] = new dataType[length * width]{ 0 };
	}

	loading_path = inputPath + "inputsForPhDStudentConference/_seg_func_5000_p3.raw";
	load3dArrayRAW<dataType>(liverShape, length, width, height, loading_path.c_str(), false);

	VoxelSpacing LiverSpacing = { sx, sx, sx };
	Image_Data LiverContainer = { height, length, width, liverShape, LiverOrigin, LiverSpacing, orientation };

	// Container for saving
	dataType** imageData = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{ 0 };
	}
	Image_Data savingContainer = { Height, Length, Width, imageData, ctOrigin, ctSpacing, orientation };

	imageInterpolation3D(LiverContainer, savingContainer, NEAREST_NEIGHBOR);

	// Ajust the contour
	// We need to ajust the segmentation because it the direct result
	// GSUBSURF so, we have different pixel value

	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (imageData[k][i] < 0.2) {
				imageData[k][i] = 1.0;
			}
			else {
				imageData[k][i] = 0.0;
			}
		}
	}

	storing_path = outputPath + "segment_liver_p3.raw";
	store3dRawData<dataType>(imageData, Length, Width, Height, storing_path.c_str());

	for (k = 0; k < height; k++) {
		delete[] liverShape[k];
	}
	delete[] liverShape;

	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
	}
	delete[] imageData;
	*/

	//======================== Load path .csv ===================================
	
	/*
	size_t height = 866;
	dataType** imageDataPtr = new dataType * [height];
	for (k = 0; k < height; k++) {
		imageDataPtr[k] = new dataType[dim2D]{ 0 };
	}

	VoxelSpacing intSpacing = { ctSpacing.sx, ctSpacing.sy, ctSpacing.sx };
	Image_Data CT = {
		Length,
		Width,
		height,
		imageDataPtr,
		ctOrigin,
		intSpacing,
		orientation
	};
	
	FILE* path_file;
	loading_path = inputPath + "csv/Path/path_ordered_p2.csv";
	if (fopen_s(&path_file, loading_path.c_str(), "r") != 0) {
		printf("Enable to open");
		return false;
	}
	dataType x = 0.0, y = 0.0, z = 0.0;

	size_t count_point = 0;
	while (feof(path_file) == 0) {
		fscanf_s(path_file, "%f", &x);
		fscanf_s(path_file, ",");
		fscanf_s(path_file, "%f", &y);
		fscanf_s(path_file, ",");
		fscanf_s(path_file, "%f", &z);
		fscanf_s(path_file, "\n");
		if (count_point % 1 == 0) {
			Point3D current_point = { x, y, z };
			Point3D point_box = getImageCoordFromRealCoord3D(current_point, ctOrigin, intSpacing, orientation);
			BoundingBox3D box = findBoundingBox3D(point_box, Length, Width, height, 25.0, 0.0);
			for (k = box.k_min; k <= box.k_max; k++) {
				for (i = box.i_min; i <= box.i_max; i++) {
					for (j = box.j_min; j <= box.j_max; j++) {
						Point3D point_sphere = { i, j, k };
						point_sphere = getRealCoordFromImageCoord3D(point_sphere, ctOrigin, intSpacing, orientation);
						double pDistance = getPoint3DDistance(current_point, point_sphere);
						if (pDistance <= 24.0) {
							imageDataPtr[k][x_new(i, j, Length)] = 1.0;
						}
					}
				}
			}
		}
		count_point++;
	}
	fclose(path_file);

	storing_path = outputPath + "path_shape.raw";
	store3dRawData<dataType>(imageDataPtr, Length, Width, height, storing_path.c_str());

	//Point3D point_before = { 6.378235, 71.263824, -647.862732};
	//std::cout << "Input point P = {" << point_before.x << ", " << point_before.y << ", " << point_before.z << ")" << std::endl;
	//point_before = getImageCoordFromRealCoord3D(point_before, ctOrigin, intSpacing, orientation);
	//std::cout << "Image Coord P = {" << point_before.x << ", " << point_before.y << ", " << point_before.z << ")" << std::endl;
	//Point3D point_after = getRealCoordFromImageCoord3D(point_before, ctOrigin, intSpacing, orientation);
	//std::cout << "Input point P = {" << point_after.x << ", " << point_after.y << ", " << point_after.z << ")" << std::endl;


	for (k = 0; k < height; k++) {
		delete[] imageDataPtr[k];
	}
	delete[] imageDataPtr;
	*/

	//======================== find Liver Approximative centroid ==============================

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

	//======================== Automatic quantitative analysis ==================

	/*
	// Load PET
	Vtk_File_Info* petContainer = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	petContainer->operation = copyFrom;

	loading_path = inputPath + "vtk/pet/Patient3_pet.vtk";
	readVtkFile(loading_path.c_str(), petContainer);

	int height = petContainer->dimensions[2];
	int length = petContainer->dimensions[0];
	int width = petContainer->dimensions[1];
	int dim2d = length * width;
	std::cout << "PET image dim : " << petContainer->dimensions[0] << " x " << petContainer->dimensions[1] << " x " << petContainer->dimensions[2] << "" << std::endl;

	std::cout << "PET origin : (" << petContainer->origin[0] << ", " << petContainer->origin[1] << ", " << petContainer->origin[2] << ")" << std::endl;
	Point3D petOrigin = { petContainer->origin[0], petContainer->origin[1], petContainer->origin[2] };
	VoxelSpacing petSpacing = { petContainer->spacing[0], petContainer->spacing[1], petContainer->spacing[2] };
	std::cout << "PET spacing : (" << petContainer->spacing[0] << ", " << petContainer->spacing[1] << ", " << petContainer->spacing[2] << ")" << std::endl;

	dataType** liverShape = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		liverShape[k] = new dataType[dim2D]{ 0 };
	}

	// Load liver segment
	loading_path = inputPath + "inputsForPhDStudentConference/segment_liver_p3.raw";
	load3dArrayRAW<dataType>(liverShape, Length, Width, Height, loading_path.c_str(), false);

	//Find Liver centroid
	dataType* centroid = new dataType[3];
	centroidImage(liverShape, centroid, Height, Length, Width, 1.0);

	Point3D center_liver_img_coord = { centroid[0], centroid[1], centroid[2] };
	Point3D center_liver_real_coord = getRealCoordFromImageCoord3D(center_liver_img_coord, ctOrigin, ctSpacing, orientation);
	double radius_cetroid_liver = 20;
	
	//Draw ball inside the liver
	dataType** ball = new dataType*[Height];
	for (k = 0; k < Height; k++) {
		ball[k] = new dataType[dim2D]{0};
	}

	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			Point3D current_point = { i, j, k };
	//			current_point = getRealCoordFromImageCoord3D(current_point, ctOrigin, ctSpacing, orientation);
	//			double dist = getPoint3DDistance(center_liver_real_coord, current_point);
	//			if (dist <= radius_cetroid_liver) {
	//				ball[k][x_new(i, j, Length)] = 1.0;
	//			}
	//		}
	//	}
	//}
	//storing_path = outputPath + "ball_liver_p2.raw";
	//store3dRawData<dataType>(ball, Length, Width, Height, storing_path.c_str());

	//Get the SUV mean around the liver centroid
	dataType mean_liver = 1.0, sum_val = 0;
	int count_point = 0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				Point3D current_point = { i, j, k };
				current_point = getRealCoordFromImageCoord3D(current_point, ctOrigin, ctSpacing, orientation);
				double d = getPoint3DDistance(center_liver_real_coord, current_point);
				if (d <= radius_cetroid_liver) {
					Point3D point_pet = getImageCoordFromRealCoord3D(current_point, petOrigin, petSpacing, orientation);
					dataType val_pet = petContainer->dataPointer[(size_t)point_pet.z][x_new((size_t)point_pet.x, (size_t)point_pet.y, length)];
					count_point++;
					sum_val += val_pet;
				}
			}
		}
	}
	if (count_point != 0) {
		mean_liver = sum_val / (dataType)count_point;
	}
	else {
		mean_liver = 1.0;
	}

	cout << "Mean inside liver is : " << mean_liver << endl;
	
	// Load the path points
	FILE* path_file;
	//loading_path = inputPath + "inputsForPhDStudentConference/path_ordered_p2.csv";
	loading_path = inputPath + "inputsForPhDStudentConference/translated.csv";
	if (fopen_s(&path_file, loading_path.c_str(), "r") != 0) {
		printf("Enable to open");
		return false;
	}

	dataType x = 0, y = 0, z = 0;
	vector<Point3D> path_points;
	while (feof(path_file) == 0) {
		fscanf_s(path_file, "%f", &x);
		fscanf_s(path_file, ",");
		fscanf_s(path_file, "%f", &y);
		fscanf_s(path_file, ",");
		fscanf_s(path_file, "%f", &z);
		fscanf_s(path_file, "\n");
		Point3D current_point = { x, y, z };
		path_points.push_back(current_point);
	}
	fclose(path_file);

	//File for the ratio max
	string saving_csv = outputPath + "ratio_max_35_after.csv";
	FILE* file;
	if (fopen_s(&file, saving_csv.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	fprintf(file, "x,y\n");

	//File for the ratio mean
	saving_csv = outputPath + "ratio_mean_35_after.csv";
	FILE* rmean;
	if (fopen_s(&rmean, saving_csv.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	fprintf(rmean, "x,y\n");

	dataType** shape_aorta_ct = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		shape_aorta_ct[k] = new dataType[dim2D]{ 0 };
	}
	
	double search_radius = 35;
	BoundingBox3D box_path_point;
	size_t n = 0;
	dataType maxSuv = 0.0, sumSuv = 0.0, meanSuv = 0.0;
	dataType ratio_max_aorta_mean_liver = 0.0, ratio_mean_aorta_mean_liver = 0.0;
	vector<dataType> ratio_list;
	
	//vector<Point3D> points_max_list;

	for (n = 0; n < path_points.size(); n++) {
		Point3D point_ct = getImageCoordFromRealCoord3D(path_points[n], ctOrigin, ctSpacing, orientation);
		box_path_point = findBoundingBox3D(point_ct, Length, Width, Height, search_radius, 2 * ctSpacing.sx);
		maxSuv = 0.0;
		sumSuv = 0.0;
		ratio_max_aorta_mean_liver = 0.0;
		count_point = 0;
		Point3D point_max = { 0.0, 0.0, 0.0 };
		for (k = box_path_point.k_min; k <= box_path_point.k_max; k++) {
			for (i = box_path_point.i_min; i <= box_path_point.i_max; i++) {
				for (j = box_path_point.j_min; j <= box_path_point.j_max; j++) {
					Point3D current_point = { i, j, k };
					current_point = getRealCoordFromImageCoord3D(current_point, ctOrigin, ctSpacing, orientation);
					double dist = getPoint3DDistance(current_point, path_points[n]);
					if (dist <= search_radius) {
						shape_aorta_ct[k][x_new(i, j, Length)] = 1.0;
						Point3D point_pet = getImageCoordFromRealCoord3D(current_point, petOrigin, petSpacing, orientation);
						dataType val_pet = petContainer->dataPointer[(size_t)point_pet.z][x_new((size_t)point_pet.x, (size_t)point_pet.y, length)];
						sumSuv += val_pet;
						count_point++;
						if (maxSuv < val_pet) {
							maxSuv = val_pet;
						}
					}
				}
			}
		}
		meanSuv = sumSuv / (dataType)count_point;
		if (mean_liver != 0) {
			ratio_max_aorta_mean_liver = maxSuv / mean_liver;
			ratio_mean_aorta_mean_liver = meanSuv / mean_liver;
		}
		else {
			ratio_max_aorta_mean_liver = maxSuv;
			ratio_mean_aorta_mean_liver = meanSuv;
		}

		//ratio_list.push_back(ratio_aorta_liver);
		fprintf(file, "%d,%f\n", n, ratio_max_aorta_mean_liver);
		fprintf(rmean, "%d,%f\n", n, ratio_mean_aorta_mean_liver);

		
		//fprintf(rmean, "%d,%f\n", n, meanSuv);

		//dataType vis = 1.0;
		//fprintf(visual, "%d,%f\n", n, vis);
	}

	storing_path = outputPath + "shape_aorta_35_p3.raw";
	store3dRawData<dataType>(shape_aorta_ct, Length, Width, Height, storing_path.c_str());

	fclose(file);
	fclose(rmean);
	//fclose(rpeak);

	////Find the max in the list
	//dataType max_ratio_list = 0;
	//for (n = 0; n < ratio_list.size(); n++) {
	//	if (max_ratio_list < ratio_list[n]) {
	//		max_ratio_list = ratio_list[n];
	//	}
	//}
	////Find points with max ratio and save them
	//saving_csv = outputPath + "ratio_max.csv";
	//FILE* file_ratio;
	//if (fopen_s(&file_ratio, saving_csv.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(file_ratio, "x,y,z\n");
	//for (n = 0; n < ratio_list.size(); n++) {
	//	if (ratio_list[n] >= max_ratio_list - 0.15 && ratio_list[n] <= max_ratio_list + 0.15) {
	//		fprintf(file_ratio, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	//	}
	//}
	//fclose(file_ratio);

	//storing_path = outputPath + "shape_aorta_p2.raw";
	//store3dRawData<dataType>(shape_aorta_pet, length, width, height, storing_path.c_str());
	
	//for (k = 0; k < height; k++) {
	//	delete[] shape_aorta_pet[k];
	//}
	//delete[] shape_aorta_pet;

	delete[] centroid;
	for (k = 0; k < Height; k++) {
		delete[] liverShape[k];
		delete[] ball[k];
		delete[] shape_aorta_ct[k];
	}
	delete[] liverShape;
	delete[] ball;
	delete[] shape_aorta_ct;

	free(petContainer);
	*/

	//======================== Collect PET values ===============================
	
	/*
	string beforePath = "C:/Users/Konan Allaly/Documents/Tests/input/vasculitis/patient_before_treatment/";
	string afterPath = "C:/Users/Konan Allaly/Documents/Tests/input/vasculitis/patient_after_treatment/";
	
	// Load PET
	Vtk_File_Info* petContainer = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	petContainer->operation = copyFrom;
	
	loading_path = inputPath + "vtk/pet/Patient2_pet.vtk";
	//loading_path = inputPath + "vtk/pet/Patient3_pet.vtk";
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

	////Load segments
	//Vtk_File_Info* ballLiver = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	//ballLiver->operation = copyFrom;
	////loading_path = beforePath + "ball_before_auto.vtk";
	////loading_path = beforePath + "ball_liver_before_m.vtk";
	////loading_path = afterPath + "ball_after_auto.vtk";
	////loading_path = afterPath + "ball_liver_after_m.vtk";
	//readVtkFile(loading_path.c_str(), ballLiver);
	////std::cout << "Origin ball : (" << ballLiver->origin[0] << ", " << ballLiver->origin[1] << ", " << ballLiver->origin[2] << ")" << std::endl;
	////std::cout << "dim ball liver : " << ballLiver->dimensions[0] << "x" << ballLiver->dimensions[1] << "x" << ballLiver->dimensions[2] << std::endl;

	//Vtk_File_Info* segment1 = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	//segment1->operation = copyFrom;
	////loading_path = beforePath + "seg1_before.vtk";
	//loading_path = afterPath + "seg1_after.vtk";
	//readVtkFile(loading_path.c_str(), segment1);
	////std::cout << "Origin seg 1 : (" << segment1->origin[0] << ", " << segment1->origin[1] << ", " << segment1->origin[2] << ")" << std::endl;
	////std::cout << "dim segment 1 : " << segment1->dimensions[0] << "x" << segment1->dimensions[1] << "x" << segment1->dimensions[2] << std::endl;

	//Vtk_File_Info* segment2 = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	//segment2->operation = copyFrom;
	////loading_path = beforePath + "seg2_before.vtk";
	//loading_path = afterPath + "seg2_after.vtk";
	//readVtkFile(loading_path.c_str(), segment2);
	////std::cout << "Origin seg 2 : (" << segment2->origin[0] << ", " << segment2->origin[1] << ", " << segment2->origin[2] << ")" << std::endl;
	////std::cout << "dim segment 2 : " << segment2->dimensions[0] << "x" << segment2->dimensions[1] << "x" << segment2->dimensions[2] << std::endl;

	//Vtk_File_Info* segment3 = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	//segment3->operation = copyFrom;
	////loading_path = beforePath + "seg3_before.vtk";
	//loading_path = afterPath + "seg3_after.vtk";
	//readVtkFile(loading_path.c_str(), segment3);
	////std::cout << "Origin seg 3 : (" << segment3->origin[0] << ", " << segment3->origin[1] << ", " << segment3->origin[2] << ")" << std::endl;
	////std::cout << "dim segment 3 : " << segment3->dimensions[0] << "x" << segment3->dimensions[1] << "x" << segment3->dimensions[2] << std::endl;

	//Vtk_File_Info* segment4 = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	//segment4->operation = copyFrom;
	////loading_path = beforePath + "seg4_before.vtk";
	//loading_path = afterPath + "seg4_after.vtk";
	//readVtkFile(loading_path.c_str(), segment4);
	////std::cout << "Origin seg 4 : (" << segment4->origin[0] << ", " << segment4->origin[1] << ", " << segment4->origin[2] << ")" << std::endl;
	////std::cout << "dim segment 4 : " << segment4->dimensions[0] << "x" << segment4->dimensions[1] << "x" << segment4->dimensions[2] << std::endl;

	//Vtk_File_Info* segment5 = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	//segment5->operation = copyFrom;
	////loading_path = beforePath + "seg5_before.vtk";
	//loading_path = afterPath + "seg5_after.vtk";
	//readVtkFile(loading_path.c_str(), segment5);
	////std::cout << "Origin seg 5 : (" << segment5->origin[0] << ", " << segment5->origin[1] << ", " << segment5->origin[2] << ")" << std::endl;
	////std::cout << "dim segment 5 : " << segment5->dimensions[0] << "x" << segment5->dimensions[1] << "x" << segment5->dimensions[2] << std::endl;

	size_t l = 0, w = 0, h = 0;

	////find max uptake in Liver
	//dataType max_uptake_liver = 0.0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			xd = x_new(i, j, Length);
	//			if (ballLiver->dataPointer[k][xd] != 0.0) {
	//				Point3D point_ball = { i, j, k };
	//				point_ball = getRealCoordFromImageCoord3D(point_ball, ctOrigin, ctSpacing, orientation);
	//				point_ball = getImageCoordFromRealCoord3D(point_ball, petOrigin, petSpacing, orientation);
	//				l = (size_t)point_ball.x;
	//				w = (size_t)point_ball.y;
	//				h = (size_t)point_ball.z;
	//				dataType value = petContainer->dataPointer[h][x_new(l, w, length)];
	//				if (value > max_uptake_liver) {
	//					max_uptake_liver = value;
	//				}
	//			}
	//		}
	//	}
	//}
	//std::cout << "Max uptake liver = " << max_uptake_liver << std::endl;
	//if (max_uptake_liver == 0) {
	//	max_uptake_liver = 1.0; //avoid division by 0 in next step
	//}

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

	////Profil ratio segment manual experiment
	//dataType max1 = 0.0, max2 = 0.0, max3 = 0.0, max4 = 0.0, max5 = 0.0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			xd = x_new(i, j, Length);
	//			Point3D point_pet = { i, j, k };
	//			point_pet = getRealCoordFromImageCoord3D(point_pet, ctOrigin, ctSpacing, orientation);
	//			point_pet = getImageCoordFromRealCoord3D(point_pet, petOrigin, petSpacing, orientation);
	//			if (point_pet.x < 0) {
	//				l = 0;
	//			}
	//			else {
	//				if (point_pet.x > length - 1) {
	//					l = length - 1;
	//				}
	//				else {
	//					l = (size_t)point_pet.x;
	//				}
	//			}
	//			if (point_pet.y < 0) {
	//				w = 0;
	//			}
	//			else {
	//				if (point_pet.y > width - 1) {
	//					w = width - 1;
	//				}
	//				else {
	//					w = (size_t)point_pet.y;
	//				}
	//			}
	//			if (point_pet.z < 0) {
	//				h = 0;
	//			}
	//			else {
	//				if (point_pet.z > height - 1) {
	//					h = height - 1;
	//				}
	//				else {
	//					h = (size_t)point_pet.z;
	//				}
	//			}
	//			dataType value = petContainer->dataPointer[h][x_new(l, w, length)];
	//			
	//			//find the max in segment 1 if the point lies in segment 1
	//			if (segment1->dataPointer[k][xd] != 0) {
	//				if (value > max1) {
	//					max1 = value;
	//				}
	//			}
	//			//find the max in segment 2 if the point lies in segment 2
	//			if (segment2->dataPointer[k][xd] != 0) {
	//				if (value > max2) {
	//					max2 = value;
	//				}
	//			}
	//			//find the max in segment 3 if the point lies in segment 3
	//			if (segment3->dataPointer[k][xd] != 0) {
	//				if (value > max3) {
	//					max3 = value;
	//				}
	//			}
	//			//find the max in segment 4 if the point lies in segment 4
	//			if (segment4->dataPointer[k][xd] != 0) {
	//				if (value > max4) {
	//					max4 = value;
	//				}
	//			}
	//			//find the max in segment 5 if the point lies in segment 5
	//			if (segment5->dataPointer[k][xd] != 0) {
	//				if (value > max5) {
	//					max5 = value;
	//				}
	//			}
	//		}
	//	}
	//}

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
	
	//std::cout << "max1 = " << max1 << ", max2 = " << max2 << ", max3 = " << max3 << ", max4 = " << max4 << ", max5 = " << max5 << std::endl;

	size_t count = 1;
	//fprintf(file, "%d,%f\n", count, max1 / max_uptake_liver);
	//count++;
	//fprintf(file, "%d,%f\n", count, max2 / max_uptake_liver);
	//count++;
	//fprintf(file, "%d,%f\n", count, max3 / max_uptake_liver);
	//count++;
	//fprintf(file, "%d,%f\n", count, max4 / max_uptake_liver);
	//count++;
	//fprintf(file, "%d,%f\n", count, max5 / max_uptake_liver);
	fclose(file);

	//free(ballLiver);
	//free(segment1);
	//free(segment2);
	//free(segment3);
	//free(segment4);
	//free(segment5);
	free(petContainer);

	//while (points_on_path.size() > 0) {
	//	points_on_path.pop_back();
	//}
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
	dataType** maskThreshold = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		maskThreshold[k] = new dataType[dim2D] {0};
	}
	if (maskThreshold == NULL)
		return false;

	// Copy
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			maskThreshold[k][i] = ctContainer->dataPointer[k][i];
		}
	}
	
	thresholding3dFunctionN(maskThreshold, Length, Width, Height, thres_min, thres_max, 0.0, 1.0);

	size_t n = 6;
	for (i = 0; i < n; i++) {
		//erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
		erosion3D(maskThreshold, Length, Width, Height, 1.0, 0.0);
	}
	//storing_path = outputPath + "threshold_p3.raw";
	//manageRAWFile3D<dataType>(maskThreshold, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

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

	//dilatation3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	//dilatation3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);

	//dataType* centroid = (dataType*)malloc(3 * sizeof(dataType));
	//centroidImage(maskThreshold, centroid, Height, Length, Width, 0.0);
	//cout << "Liver centroid : (" << centroid[0] << "," << centroid[1] << "," << centroid[2] << ")" << endl;
	//free(centroid);

	storing_path = outputPath + "liver_thres_p3.raw";
	store3dRawData<dataType>(maskThreshold, Length, Width, Height, storing_path.c_str());

	for (k = 0; k < Height; k++) {
		delete[] maskThreshold[k];
		delete[] labelArray[k];
		delete[] status[k];
	}
	delete[] maskThreshold;
	delete[] labelArray;
	delete[] status;
	*/

	//======================== Interpolate input and segment ======================
		
	/*
	//We load the rough segmentation of the liver, obtained by thresholding
	//connected component labeling and erosion. We then interpolate the original
	//image and the segment in order to have the same spacing in x, y and z.
	//Information related to the dimension, sapacing and origin are avalaible
	//by loading the .vtk file.
	
	dataType** imageData = new dataType * [Height];
	dataType** segment = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{ 0 };
		segment[k] = new dataType[dim2D]{ 0 };
	}

	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			imageData[k][i] = ctContainer->dataPointer[k][i];
		}
	}

	loading_path = inputPath + "liver/liver_p6.raw";
	load3dArrayRAW<dataType>(segment, Length, Width, Height, loading_path.c_str(), false);
	
	size_t height = (size_t)((ctSpacing.sz / ctSpacing.sx) * Height);
	cout << "interpol height : " << height << endl;

	dataType** interpolate = new dataType * [height];
	dataType** interpolateSegment = new dataType * [height];
	for (k = 0; k < height; k++) {
		interpolate[k] = new dataType[dim2D]{ 0 };
		interpolateSegment[k] = new dataType[dim2D]{ 0 };
	}

	Image_Data CT; 
	CT.imageDataPtr = imageData;
	CT.height = Height; 
	CT.length = Length; 
	CT.width = Width;
	CT.orientation = orientation; 
	CT.origin = ctOrigin; 
	CT.spacing.sx = ctSpacing.sx;
	CT.spacing.sy = ctSpacing.sy;
	CT.spacing.sz = ctSpacing.sz;

	Image_Data interpolated;
	interpolated.imageDataPtr = interpolate;
	interpolated.height = height; 
	interpolated.length = Length; 
	interpolated.width = Width;
	interpolated.orientation = orientation; 
	interpolated.origin = ctOrigin;
	interpolated.spacing.sx = ctSpacing.sx;
	interpolated.spacing.sy = ctSpacing.sy;
	interpolated.spacing.sz = ctSpacing.sx;
	
	//Interpolate the input image
	imageInterpolation3D(CT, interpolated, NEAREST_NEIGHBOR);

	storing_path = outputPath + "interpolate_p6.raw";
	store3dRawData<dataType>(interpolate, Length, Width, height, storing_path.c_str());

	// Interpolate the segment
	CT.imageDataPtr = segment;
	interpolated.imageDataPtr = interpolateSegment;

	imageInterpolation3D(CT, interpolated, NEAREST_NEIGHBOR);

	storing_path = outputPath + "interpolateSegment_p6.raw";
	store3dRawData<dataType>(interpolateSegment, Length, Width, height, storing_path.c_str());

	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
		delete[] segment[k];
	}
	delete[] imageData;
	delete[] segment;

	for (k = 0; k < height; k++) {
		delete[] interpolate[k];
		delete[] interpolateSegment[k];
	}
	delete[] interpolate;
	delete[] interpolateSegment;
	*/
	
	//======================= Crop interpolated image data =======================

	/*
	////Patient 1 : interpolated
	//size_t i_min = 90, i_max = 359;
	//size_t j_min = 90, j_max = 359;
	//size_t k_min = 307, k_max = 496;
	
	////Patient 2 : interpolated
	//size_t i_min = 120, i_max = 375;
	//size_t j_min = 160, j_max = 415;
	//size_t k_min = 335, k_max = 494;

	//Patient 3 : interpolated
	size_t i_min = 90, i_max = 345;
	size_t j_min = 130, j_max = 385;
	size_t k_min = 532, k_max = 711;

	////Patient 4 : interpoolated
	//size_t i_min = 120, i_max = 375;
	//size_t j_min = 150, j_max = 405;
	//size_t k_min = 311, k_max = 518;

	////Patient 6 : interpolated
	//size_t i_min = 120, i_max = 375;
	//size_t j_min = 180, j_max = 435;
	//size_t k_min = 703, k_max = 902;

	////Patient 3 : input
	//size_t i_min = 90, i_max = 345;
	//size_t j_min = 130, j_max = 385;
	//size_t k_min = 250, k_max = 333;

	const size_t length_crop = i_max - i_min + 1;
	const size_t width_crop = j_max - j_min + 1;
	const size_t height_crop = k_max - k_min + 1;
	const size_t dim2D_crop = length_crop * width_crop;
	std::cout << "cropped dim : " << length_crop << "x" << width_crop << "x" << height_crop << std::endl;

	dataType** imageCrop = new dataType * [height_crop];
	dataType** segmentCrop = new dataType * [height_crop];
	for (k = 0; k < height_crop; k++) {
		imageCrop[k] = new dataType[dim2D_crop]{ 0 };
		segmentCrop[k] = new dataType[dim2D_crop]{ 0 };
	}

	//const dataType height_interp = 866; // patient2
	const dataType height_interp = 1083; // patient3
	dataType** imageDataPtr = new dataType * [height_interp];
	dataType** segmentInterpPtr = new dataType * [height_interp];
	for (k = 0; k < height_interp; k++) {
		imageDataPtr[k] = new dataType[dim2D]{ 0 };
		segmentInterpPtr[k] = new dataType[dim2D]{ 0 };
	}

	dataType** segmentPtr = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		segmentPtr[k] = new dataType[dim2D]{ 0 };
	}

	VoxelSpacing interpolate_spacing = { ctSpacing.sx, ctSpacing.sy, ctSpacing.sx };
	Image_Data Interpolation = { 
		height_interp,
		Length,
		Width,
		segmentInterpPtr,
		ctOrigin,
		interpolate_spacing,
		orientation 
	};

	Image_Data inputStructure = {
	Height,
	Length,
	Width,
	segmentPtr,
	ctOrigin,
	ctSpacing,
	orientation
	};

	loading_path = inputPath + "raw/interpolated/patient3/filtered_p3.raw";
	load3dArrayRAW<dataType>(imageDataPtr, Length, Width, height_interp, loading_path.c_str(), false);

	loading_path = inputPath + "inputsForPhDStudentConference/liver_seg_p3.raw";
	load3dArrayRAW<dataType>(segmentPtr, Length, Width, Height, loading_path.c_str(), false);

	imageInterpolation3D(inputStructure, Interpolation, NEAREST_NEIGHBOR);

	//cropping
	size_t in = 0, jn = 0, kn = 0, xdn = 0;
	for (kn = 0, k = k_min; kn < height_crop; kn++, k++) {
		for (in = 0, i = i_min; in < length_crop; in++, i++) {
			for (jn = 0, j = j_min; jn < width_crop; jn++, j++) {
				xd = x_new(i, j, Length);
				xdn = x_new(in, jn, length_crop);
				imageCrop[kn][xdn] = imageDataPtr[k][xd];
				segmentCrop[kn][xdn] = segmentInterpPtr[k][xd];
			}
		}
	}

	Point3D newOrigin = { i_min, j_min, k_min };
	newOrigin = getRealCoordFromImageCoord3D(newOrigin, ctOrigin, interpolate_spacing, orientation);
	std::cout << "New origin : (" << newOrigin.x << "," << newOrigin.y << "," << newOrigin.z << ")" << std::endl;

	storing_path = outputPath + "cropImage_p3.raw";
	store3dRawData<dataType>(imageCrop, length_crop, width_crop, height_crop, storing_path.c_str());

	storing_path = outputPath + "cropSegment_p3.raw";
	store3dRawData<dataType>(segmentCrop, length_crop, width_crop, height_crop, storing_path.c_str());

	for (k = 0; k < height_crop; k++) {
		delete[] imageCrop[k];
		delete[] segmentCrop[k];
	}
	delete[] imageCrop;
	delete[] segmentCrop;

	for (k = 0; k < Height; k++) {
		delete[] segmentPtr[k];
	}
	delete[] segmentPtr;

	for (k = 0; k < height_interp; k++) {
		delete[] imageDataPtr[k];
		delete[] segmentInterpPtr[k];
	}
	delete[] imageDataPtr;
	delete[] segmentInterpPtr;
	*/

	//======================== Refinement by generalized subjective surface =======
	
	/*
	// We work with the cropped data to speed up the computation
	
	////Patient 1 : interpolated
	//size_t i_min = 90, i_max = 359;
	//size_t j_min = 90, j_max = 359;
	//size_t k_min = 307, k_max = 496;
	
	//Patient 2 : interpolated
	size_t i_min = 120, i_max = 375;
	size_t j_min = 160, j_max = 415;
	size_t k_min = 335, k_max = 494;

	////Patient 3 : interpolated
	//size_t i_min = 90, i_max = 345;
	//size_t j_min = 130, j_max = 385;
	//size_t k_min = 532, k_max = 711;

	////Patient 4 : interpolated
	//size_t i_min = 120, i_max = 375;
	//size_t j_min = 150, j_max = 405;
	//size_t k_min = 311, k_max = 518;

	////Patient 6: interpolated
	//size_t i_min = 120, i_max = 375;
	//size_t j_min = 180, j_max = 435;
	//size_t k_min = 703, k_max = 902;

	const size_t length = i_max - i_min + 1;
	const size_t width = j_max - j_min + 1;
	const size_t height = k_max - k_min + 1;
	const size_t dim2d = length * width;

	dataType** imageData = new dataType * [height];
	dataType** initialSegment = new dataType * [height];
	for (k = 0; k < height; k++) {
		imageData[k] = new dataType[dim2d]{0};
		initialSegment[k] = new dataType[dim2d]{0};
	}

	loading_path = inputPath + "inputsForPhDStudentConference/cropImage_p2.raw";
	load3dArrayRAW<dataType>(imageData, length, width, height, loading_path.c_str(), false);

	loading_path = inputPath + "inputsForPhDStudentConference/cropSegment_p2.raw";
	load3dArrayRAW<dataType>(initialSegment, length, width, height, loading_path.c_str(), false);

	VoxelSpacing ctInt = { ctContainer->spacing[0], ctContainer->spacing[1], ctContainer->spacing[0] };
	Point3D newOrigin = { i_min, j_min, k_min };
	newOrigin = getRealCoordFromImageCoord3D(newOrigin, ctOrigin, ctInt, orientation);
	std::cout << "New origin : (" << newOrigin.x << "," << newOrigin.y << "," << newOrigin.z << ")" << std::endl;

	Image_Data inputImage{
		height,
		length,
		width,
		imageData,
		newOrigin,
		ctInt,
		orientation
	};
	Image_Data inputSegment{
		height,
		length,
		width,
		initialSegment,
		newOrigin,
		ctInt,
		orientation
	};
	
	Point3D* center_seg = new Point3D[1];
	center_seg->x = 0.0, center_seg->y = 0.0; center_seg->z = 0.0;

	//Down sampling
	dataType down_factor = 4.0;
	const size_t h_down = (size_t)(height / down_factor);
	const size_t l_down = (size_t)(length / down_factor);
	const size_t w_down = (size_t)(width / down_factor);
	std::cout << "Downsampled image dim : " << l_down << " x " << w_down << " x " << h_down << "" << std::endl;
	dataType** imageDown = new dataType * [h_down];
	dataType** initialDown = new dataType * [h_down];
	for (k = 0; k < h_down; k++) {
		imageDown[k] = new dataType[l_down * w_down]{ 0 };
		initialDown[k] = new dataType[l_down * w_down]{ 0 };
	}

	VoxelSpacing downSampleSapcing = { ctInt.sx * down_factor, ctInt.sy * down_factor, ctInt.sz * down_factor };
	
	Image_Data downSampledImage{
	h_down,
	l_down,
	w_down,
	imageDown,
	newOrigin,
	downSampleSapcing,
	orientation
	};
	Image_Data downSampledSegment{
	h_down,
	l_down,
	w_down,
	initialDown,
	newOrigin,
	downSampleSapcing,
	orientation
	};
	std::cout << "Downsampled image spacing : (" << downSampleSapcing.sx << ", " << downSampleSapcing.sy << ", " << downSampleSapcing.sz << ")" << std::endl;

	imageInterpolation3D(inputImage, downSampledImage, NEAREST_NEIGHBOR);
	imageInterpolation3D(inputSegment, downSampledSegment, NEAREST_NEIGHBOR);

	//storing_path = outputPath + "downSeg.raw";
	//store3dRawData<dataType>(initialDown, l_down, w_down, h_down, storing_path.c_str());
	//storing_path = outputPath + "downImg.raw";
	//store3dRawData<dataType>(imageDown, l_down, w_down, h_down, storing_path.c_str());

	//Smoothing by heat equation, we don't need all the parameters
	Filter_Parameters smoothing_parameters{
		0.2, //tau
		1.0, //h
		0.0, //sigma
		0.0, //edge detector coeficient
		1.5, //omega
		1e-3,//tolerance
		0.0, // epsilon Evan-Spruck regularization
		0.0, // coef
		(size_t)1, //p
		(size_t)1, //time step number
		(size_t)100 //max number of solver iterations
	};

	Segmentation_Parameters segmentation_parameters{
		(size_t)30,  // Maximum number of Gauss - Seidel iterations
		2000,         // K : edge detector coef
		1e-6,         //epsilon : Evan-Spruck regularization
		(size_t)5000,   //number of current time step
		(size_t)5000,   //max number of time step
		(size_t)100,    //frequence au saving
		1e-4,         //segmentation tolerance
		down_factor,  //tau
		1.0,          //h
		1.5,          //omega
		1e-3,         //gauss_seidel tolerance
		0.01,         //convection coeff
		1.0           //diffusion coeff
	};
	
	unsigned char segmentPath[] = "C:/Users/Konan Allaly/Documents/Tests/output/segmentation/patient2/";
	//subsurfSegmentation(imageToBeSegmented, initialSegment, seg_params, smooth_params, center_seg, 1, segmentPath);
	generalizedSubsurfSegmentation(downSampledImage, initialDown, segmentation_parameters, smoothing_parameters, center_seg, 1, segmentPath);

	for (k = 0; k < h_down; k++) {
		delete[] imageDown[k];
		delete[] initialDown[k];
	}
	delete[] imageDown;
	delete[] initialDown;

	delete[] center_seg;
	for (k = 0; k < height; k++) {
		delete[] imageData[k];
		delete[] initialSegment[k];
	}
	delete[] imageData;
	delete[] initialSegment;
	*/

	//======================== Path finding from one seed in the ascending aorta ====
	
	/*
	//const size_t height = 866; // P2
	//const size_t height = 1351; // P6
	//const size_t height = 880; // P6
	//const size_t height = 844; // P4
	//const size_t height = 1083; // P3
	const size_t height = 866; // P2

	dataType** potential = new dataType * [height];
	dataType** action_field = new dataType * [height];
	dataType** pathPtr = new dataType * [height];
	dataType** imageData = new dataType * [height];
	for (k = 0; k < height; k++) {
		potential[k] = new dataType[dim2D]{ 0 };
		action_field[k] = new dataType[dim2D]{ 0 };
		pathPtr[k] = new dataType[dim2D]{ 0 };
		imageData[k] = new dataType[dim2D]{ 0 };
	}
	if (potential == NULL || action_field == NULL || pathPtr == NULL || imageData == NULL)
		return false;

	loading_path = inputPath + "raw/interpolated/patient2/filtered_p2.raw";
	//loading_path = inputPath + "raw/interpolated/patient2/filtered_p2.raw";
	//loading_path = inputPath + "raw/interpolated/patient4/filtered_p4.raw";
	//loading_path = inputPath + "raw/filtered/filteredGMC_p6.raw";
	load3dArrayRAW<dataType>(imageData, Length, Width, height, loading_path.c_str(), false);
	
	//VoxelSpacing intSpacing = { ctSpacing.sx, ctSpacing.sy, ctSpacing.sx };
	//Image_Data CTinterpolated = { height, Length, Width, imageData, ctOrigin, intSpacing, orientation };

	////Patient 2
	//Point3D seed1 = { 260, 257, 312 };
	//Point3D seed2 = { 295, 317, 565 };

	////Patient 3
	//Point3D seed1 = { 259, 254, 536 }; //old
	//Point3D seed2 = { 288, 319, 769 };//old

	////Patient 3 PhD student conference
	//Point3D seed1 = { 257, 255, 522 };
	//Point3D seed2 = { 288, 313, 779 };
	//Point3D seed3 = { 253, 251, 779 };

	//Patient 2 PhD student conference
	Point3D seed1 = { 260, 257, 306 };//*
	Point3D seed2 = { 266, 268, 393 };
	Point3D seed3 = { 288, 304, 480 };
	Point3D seed4 = { 295, 313, 566 };//*
	Point3D seed5 = { 279, 282, 586 };
	Point3D seed6 = { 256, 258, 566 };//*
	Point3D seed7 = { 253, 254, 539 };

	//Point3D seed1 = { 257, 251, 530 };
	//Point3D seed1 = { 257, 251, 526 };
	//Point3D seed1 = { 260, 257, 312 };
	//Point3D seed2 = { 295, 317, 565 };
	//Point3D seed2 = { 0, 0, 0 };

	////Patient 6
	//Point3D seed1 = { 266, 242, 734 };
	
	////Test 1
	//Point3D seed1 = { 263, 244, 717 };
	//Point3D seed2 = { 283, 283, 993 };

	//Test 2
	//Point3D seed1 = { 261, 246, 717 };
	//Point3D seed2 = { 283, 283, 993 };

	////Test 3
	//Point3D seed1 = { 261, 246, 717 };
	//Point3D seed2 = { 283, 283, 993 };
	//Point3D seed3 = { 243, 219, 993 };

	////Patient 4
	//Point3D seed1 = { 268, 232, 291 };
	//Point3D seed2 = { 289, 297, 591 };
	////Point3D seed3 = { 256, 224, 591 };

	//std::cout << "P1 : (" << seed1.x << ", " << seed1.y << ", " << seed1.z << ")" << std::endl;
	//std::cout << "P2 : (" << seed2.x << ", " << seed2.y << ", " << seed2.z << ")" << std::endl;
	//seed1 = getRealCoordFromImageCoord3D(seed1, ctOrigin, intSpacing, orientation);
	//seed2 = getRealCoordFromImageCoord3D(seed2, ctOrigin, intSpacing, orientation);
	//seed1 = getImageCoordFromRealCoord3D(seed1, ctOrigin, spacingCT, orientation);
	//seed2 = getImageCoordFromRealCoord3D(seed2, ctOrigin, spacingCT, orientation);
	//std::cout << "P1 updated : (" << seed1.x << ", " << seed1.y << ", " << seed1.z << ")" << std::endl;
	//std::cout << "P2 updated : (" << seed2.x << ", " << seed2.y << ", " << seed2.z << ")" << std::endl;

	VoxelSpacing imageSpacing = { ctSpacing.sx, ctSpacing.sy, ctSpacing.sx };
	Image_Data ctImageData{
		height,
		Length,
		Width,
		imageData,
		ctOrigin,
		imageSpacing,
		orientation
	};
	
	double radius = 3.0;
	dataType h = 1.0;
	Potential_Parameters potential_parameters{
		0.005, //K
		0.01,  //epsilon
		0.0,   //coef min
		0.0,   // coef max
		0.0,   // coef mean
		0.0,   //ceof sd
		0.0,   //coef ct
		0.0   //coef pet
	};

	compute3DPotential(ctImageData, potential, seed1, radius, potential_parameters);
	
	vector<Point3D> path_points;
	Point3D* seedPoints = new Point3D[2];
	int n = 0;

	string saving_csv = outputPath + "pathFromSevenPoints_p2.csv";
	FILE* file_real;
	if (fopen_s(&file_real, saving_csv.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}

	
	///////////////////// First segment ////////////////////////////
	seedPoints[0] = seed1;
	seedPoints[1] = seed2;
	partialFrontPropagation(action_field, potential, Length, Width, height, seedPoints);
	shortestPath3D(action_field, Length, Width, height, imageSpacing, seedPoints, path_points);

	for (n = path_points.size() - 1; n > -1; n--) {
		path_points[n] = getRealCoordFromImageCoord3D(path_points[n], ctOrigin, imageSpacing, orientation);
		fprintf(file_real, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	}
	while (path_points.size() > 0) {
		path_points.pop_back();
	}

	///////////////// Second segment ///////////////////////////
	seedPoints[0] = seed2;
	seedPoints[1] = seed3;
	partialFrontPropagation(action_field, potential, Length, Width, height, seedPoints);
	shortestPath3D(action_field, Length, Width, height, imageSpacing, seedPoints, path_points);

	for (n = path_points.size() - 1; n > -1; n--) {
		path_points[n] = getRealCoordFromImageCoord3D(path_points[n], ctOrigin, imageSpacing, orientation);
		fprintf(file_real, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	}
	while (path_points.size() > 0) {
		path_points.pop_back();
	}

	///////////////// Third segment ///////////////////////////
	seedPoints[0] = seed3;
	seedPoints[1] = seed4;
	partialFrontPropagation(action_field, potential, Length, Width, height, seedPoints);
	shortestPath3D(action_field, Length, Width, height, imageSpacing, seedPoints, path_points);

	for (n = path_points.size() - 1; n > -1; n--) {
		path_points[n] = getRealCoordFromImageCoord3D(path_points[n], ctOrigin, imageSpacing, orientation);
		fprintf(file_real, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	}
	while (path_points.size() > 0) {
		path_points.pop_back();
	}

	///////////////// Fourth segment ///////////////////////////
	seedPoints[0] = seed4;
	seedPoints[1] = seed5;
	partialFrontPropagation(action_field, potential, Length, Width, height, seedPoints);
	shortestPath3D(action_field, Length, Width, height, imageSpacing, seedPoints, path_points);

	for (n = path_points.size() - 1; n > -1; n--) {
		path_points[n] = getRealCoordFromImageCoord3D(path_points[n], ctOrigin, imageSpacing, orientation);
		fprintf(file_real, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	}
	while (path_points.size() > 0) {
		path_points.pop_back();
	}

	///////////////// Fifth segment ///////////////////////////
	seedPoints[0] = seed5;
	seedPoints[1] = seed6;
	partialFrontPropagation(action_field, potential, Length, Width, height, seedPoints);
	shortestPath3D(action_field, Length, Width, height, imageSpacing, seedPoints, path_points);

	for (n = path_points.size() - 1; n > -1; n--) {
		path_points[n] = getRealCoordFromImageCoord3D(path_points[n], ctOrigin, imageSpacing, orientation);
		fprintf(file_real, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	}
	while (path_points.size() > 0) {
		path_points.pop_back();
	}
	

	///////////////// Sixth segment ///////////////////////////
	seedPoints[0] = seed6;
	seedPoints[1] = seed7;
	partialFrontPropagation(action_field, potential, Length, Width, height, seedPoints);
	shortestPath3D(action_field, Length, Width, height, imageSpacing, seedPoints, path_points);

	for (n = path_points.size() - 1; n > -1; n--) {
		path_points[n] = getRealCoordFromImageCoord3D(path_points[n], ctOrigin, imageSpacing, orientation);
		fprintf(file_real, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	}
	while (path_points.size() > 0) {
		path_points.pop_back();
	}
	fclose(file_real);
	delete[] seedPoints;

	for (k = 0; k < height; k++) {
		delete[] potential[k];
		delete[] action_field[k];
		delete[] pathPtr[k];
		delete[] imageData[k];
	}
	delete[] potential;
	delete[] action_field;
	delete[] pathPtr;
	delete[] imageData;
	*/

	//======================== Path finding Multiple seeds =======================================
	
	/**
	//const size_t height = 866; // P2
	//const size_t height = 844; // P4
	//const size_t height = 1351; // P6
	const size_t height = 1083; // P3

	dataType** potential = new dataType * [height];
	dataType** action_field = new dataType * [height];
	dataType** pathPtr = new dataType * [height];
	dataType** imageData = new dataType * [height];
	for (k = 0; k < height; k++) {
		potential[k] = new dataType[dim2D]{ 0 };
		action_field[k] = new dataType[dim2D]{ 0 };
		pathPtr[k] = new dataType[dim2D]{ 0 };
		imageData[k] = new dataType[dim2D]{ 0 };
	}
	if (potential == NULL || action_field == NULL || pathPtr == NULL || imageData == NULL)
		return false;
	
	loading_path = inputPath + "raw/interpolated/filtered_p3.raw";
	//loading_path = inputPath + "raw/interpolated/patient2/Im_mean_r5.raw";
	//loading_path = inputPath + "raw/interpolated/patient4/Im_mean_r3.raw";
	//loading_path = inputPath + "raw/interpolated/patient6/Im_mean_r3.raw";
	load3dArrayRAW<dataType>(imageData, Length, Width, height, loading_path.c_str(), false);

	double radius = 3.0;
	dataType h = 1.0;
	Potential_Parameters parameters{
		0.005,//K
		0.01,//eps
		h,
		radius
	};

	//Point3D* seeds = new Point3D[3];
	////Patient 2
	//Point3D seed1 = { 260, 257, 312 };
	//Point3D seed2 = { 295, 317, 565 };
	//Point3D seed3 = { 255, 258, 565 };

	////Patient 4
	//Point3D seed1 = { 268, 232, 291 };
	//Point3D seed2 = { 289, 297, 591 };
	//Point3D seed3 = { 256, 224, 591 };

	////Patient 6
	//Point3D seed1 = { 266, 242, 734 };
	//Point3D seed2 = { 283, 283, 993 };
	//Point3D seed3 = { 243, 219, 993 };

	//seeds[0] = seed1;
	//seeds[1] = seed2;
	//seeds[2] = seed3;

	VoxelSpacing intSpacing = { ctSpacing.sx, ctSpacing.sy, ctSpacing.sx };

	Image_Data CT = { height, Length, Width, imageData, ctOrigin, intSpacing, orientation };

	////Patient2
	//Point3D seed1 = { 263, 260, 316 };
	//Point3D seed2 = { 264, 255, 344 };
	//Point3D seed3 = { 263, 269, 399 };
	//Point3D seed4 = { 291, 306, 487 };
	//Point3D seed5 = { 296, 316, 567 };
	//Point3D seed6 = { 282, 287, 591 };
	//Point3D seed7 = { 259, 260, 569 };
	//Point3D seed8 = { 257, 252, 538 };

	////Patient6
	//Point3D seed1 = { 267, 243, 733 };
	//Point3D seed2 = { 266, 249, 768 };
	//Point3D seed3 = { 262, 259, 809 };
	//Point3D seed4 = { 270, 282, 849 };
	//Point3D seed5 = { 266, 302, 932 };
	//Point3D seed6 = { 276, 245, 1029 };
	//Point3D seed7 = { 245, 217, 990 };
	//Point3D seed8 = { 247, 224, 959 };

	////Patient6
	//Point3D seed1 = { 269, 230, 290 };
	//Point3D seed2 = { 267, 241, 333 };
	//Point3D seed3 = { 262, 262, 390 };
	//Point3D seed4 = { 275, 286, 464 };
	//Point3D seed5 = { 288, 295, 540 };
	//Point3D seed6 = { 282, 259, 614 };
	//Point3D seed7 = { 262, 222, 577 };
	//Point3D seed8 = { 274, 225, 551 };

	//Patient3
	Point3D seed1 = { 258, 255, 529 };
	Point3D seed2 = { 261, 263, 596 };
	Point3D seed3 = { 266, 278, 645 };
	Point3D seed4 = { 288, 321, 743 };
	Point3D seed5 = { 278, 282, 802 };
	Point3D seed6 = { 247, 250, 766 };
	Point3D seed7 = { 253, 248, 750 };
	//Point3D seed8 = { 274, 225, 551 };

	compute3DPotential(CT, potential, seed1, radius, parameters);

	Point3D* path1 = new Point3D[2];
	Point3D* path2 = new Point3D[2];
	Point3D* path3 = new Point3D[2];
	Point3D* path4 = new Point3D[2];
	Point3D* path5 = new Point3D[2];
	Point3D* path6 = new Point3D[2];
	//Point3D* path7 = new Point3D[2];

	path1[0] = seed1; path1[1] = seed2;
	path2[0] = seed2; path2[1] = seed3;
	path3[0] = seed3; path3[1] = seed4;
	path4[0] = seed4; path4[1] = seed5;
	path5[0] = seed5; path5[1] = seed6;
	path6[0] = seed6; path6[1] = seed7;
	//path7[0] = seed7; path7[1] = seed8;
	
	vector<Point3D> path_points1;
	vector<Point3D> path_points2;
	vector<Point3D> path_points3;
	vector<Point3D> path_points4;
	vector<Point3D> path_points5;
	vector<Point3D> path_points6;
	//vector<Point3D> path_points7;

	//Segment 1
	partialFrontPropagation(action_field, potential, Length, Width, height, path1);
	shortestPath3D(action_field, Length, Width, height, intSpacing, path1, path_points1);

	//Segment 2
	partialFrontPropagation(action_field, potential, Length, Width, height, path2);
	shortestPath3D(action_field, Length, Width, height, intSpacing, path2, path_points2);

	//Segment 3
	partialFrontPropagation(action_field, potential, Length, Width, height, path3);
	shortestPath3D(action_field, Length, Width, height, intSpacing, path3, path_points3);

	//Segment 4
	partialFrontPropagation(action_field, potential, Length, Width, height, path4);
	shortestPath3D(action_field, Length, Width, height, intSpacing, path4, path_points4);

	//Segment 5
	partialFrontPropagation(action_field, potential, Length, Width, height, path5);
	shortestPath3D(action_field, Length, Width, height, intSpacing, path5, path_points5);

	//Segment 6
	partialFrontPropagation(action_field, potential, Length, Width, height, path6);
	shortestPath3D(action_field, Length, Width, height, intSpacing, path6, path_points6);

	////Segment 7
	//partialFrontPropagation(action_field, potential, Length, Width, height, path7);
	//shortestPath3D(action_field, Length, Width, height, intSpacing, path7, path_points7);

	string saving_real = outputPath + "path_point_real_coord_p3.csv";
	FILE* file_real;
	if (fopen_s(&file_real, saving_real.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}

	string saving_img = outputPath + "path_point_img_coord_p3.csv";
	FILE* file_img;
	if (fopen_s(&file_img, saving_img.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}

	int n = 0;
	
	for (n = path_points1.size() - 1; n > -1; n--) {
		fprintf(file_img, "%f,%f,%f\n", path_points1[n].x, path_points1[n].y, path_points1[n].z);
		path_points1[n] = getRealCoordFromImageCoord3D(path_points1[n], CT.origin, CT.spacing, orientation);
		fprintf(file_real, "%f,%f,%f\n", path_points1[n].x, path_points1[n].y, path_points1[n].z);
	}

	for (n = path_points2.size() - 1; n > -1; n--) {
		fprintf(file_img, "%f,%f,%f\n", path_points2[n].x, path_points2[n].y, path_points2[n].z);
		path_points2[n] = getRealCoordFromImageCoord3D(path_points2[n], CT.origin, CT.spacing, orientation);
		fprintf(file_real, "%f,%f,%f\n", path_points2[n].x, path_points2[n].y, path_points2[n].z);
	}

	for (n = path_points3.size() - 1; n > -1; n--) {
		fprintf(file_img, "%f,%f,%f\n", path_points3[n].x, path_points3[n].y, path_points3[n].z);
		path_points3[n] = getRealCoordFromImageCoord3D(path_points3[n], CT.origin, CT.spacing, orientation);
		fprintf(file_real, "%f,%f,%f\n", path_points3[n].x, path_points3[n].y, path_points3[n].z);
	}

	for (n = path_points4.size() - 1; n > -1; n--) {
		fprintf(file_img, "%f,%f,%f\n", path_points4[n].x, path_points4[n].y, path_points4[n].z);
		path_points4[n] = getRealCoordFromImageCoord3D(path_points4[n], CT.origin, CT.spacing, orientation);
		fprintf(file_real, "%f,%f,%f\n", path_points4[n].x, path_points4[n].y, path_points4[n].z);
	}

	for (n = path_points5.size() - 1; n > -1; n--) {
		fprintf(file_img, "%f,%f,%f\n", path_points5[n].x, path_points5[n].y, path_points5[n].z);
		path_points5[n] = getRealCoordFromImageCoord3D(path_points5[n], CT.origin, CT.spacing, orientation);
		fprintf(file_real, "%f,%f,%f\n", path_points5[n].x, path_points5[n].y, path_points5[n].z);
	}

	for (n = path_points6.size() - 1; n > -1; n--) {
		fprintf(file_img, "%f,%f,%f\n", path_points6[n].x, path_points6[n].y, path_points6[n].z);
		path_points6[n] = getRealCoordFromImageCoord3D(path_points6[n], CT.origin, CT.spacing, orientation);
		fprintf(file_real, "%f,%f,%f\n", path_points6[n].x, path_points6[n].y, path_points6[n].z);
	}

	//for (n = path_points7.size() - 1; n > -1; n--) {
	//	fprintf(file_img, "%f,%f,%f\n", path_points7[n].x, path_points7[n].y, path_points7[n].z);
	//	path_points7[n] = getRealCoordFromImageCoord3D(path_points7[n], CT.origin, CT.spacing, orientation);
	//	fprintf(file_real, "%f,%f,%f\n", path_points7[n].x, path_points7[n].y, path_points7[n].z);
	//}

	fclose(file_img);
	fclose(file_real);

	delete[] path1;
	delete[] path2;
	delete[] path3;
	delete[] path4;
	delete[] path5;
	delete[] path6;
	//delete[] path7;

	for (k = 0; k < height; k++) {
		delete[] potential[k];
		delete[] action_field[k];
		delete[] pathPtr[k];
		delete[] imageData[k];
	}
	delete[] potential;
	delete[] action_field;
	delete[] pathPtr;
	delete[] imageData;
	*/
	
	//======================== 2D Hough transform on Slice ===========================
	
	const size_t Height = 406, Length = 512, Width = 512;
	const size_t dim2D = Length * Width;

	dataType** imageData = new dataType*[Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{ 0 };
	}
	
	//loading_path = inputPath + "raw/filtered/filteredGMC_p2.raw";
	//loading_path = inputPath + "raw/ct/patient2_ct.raw";
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
	}

	//rescaleNewRange2D(imageSlice, Length, Width, 0.0, 1.0);
	//storing_path = outputPath + "filtered.raw";
	//store2dRawData<dataType>(imageSlice, Length, Width, storing_path.c_str());

	PixelSpacing spacing = { 1.171875, 1.171875 };
	Image_Data2D IMAGE;
	IMAGE.imageDataPtr = imageSlice;
	IMAGE.height = Length; IMAGE.width = Width;
	
	HoughParameters parameters = 
	{ 
		7.5, //minimal radius
		19.0, //maximal radius
		0.5, //epsilon to search points belonging the circle
		30.0, //offset to find the bounding box
		1000.0, //edge detector coefficient K 
		1.0, //h, space discretization for gradient computation
		0.07, //threshold edge detector
		0.5, //radius step
		spacing //pixel size
	};

	double radius = 12.0;
	double perimeter = 2 * M_PI * radius;

	//storing_path = outputPath + "filtered.raw";
	//store2dRawData<dataType>(imageSlice, Length, Width, storing_path.c_str());
	
	////filtering parameters
	//Filter_Parameters filter_parameters;
	//filter_parameters.h = 1.0; filter_parameters.timeStepSize = 0.3; filter_parameters.eps2 = 1e-6;
	//filter_parameters.omega_c = 1.5; filter_parameters.tolerance = 1e-3; filter_parameters.maxNumberOfSolverIteration = 1000;
	//filter_parameters.timeStepsNum = 5; filter_parameters.coef = 1e-6;

	////Slice filtering
	//heatImplicit2dScheme(IMAGE, filter_parameters);

	//storing_path = outputPath + "smoothed.raw";
	//store2dRawData<dataType>(imageSlice, Length, Width, storing_path.c_str());
	
	//Point2D seed = { 248, 298 }; //---> P2 : automatic detection
	//Point2D seed = { 263, 257 }; //---> P2 : abdominal aorta
	//Point2D seed = { 263, 296 }; //---> P2
	//Point2D seed = { 263, 296 }; //---> P2
	//Point2D seed = { 262, 269 }; //---> P4
	//Point2D seed = { 257, 244 }; //---> P2, liver
	//Point2D found_point = { 0.0, 0.0 };
	//Point2D found1 = { 292,314 }, found2 = { 253,255 }; // p2
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
	
	//found_point = localHoughTransform(seed, imageSlice, houghSpace, foundCircle, Length, Width, parameters, storing_path.c_str());
	//cout << "found center : (" << found_point.x << "," << found_point.y << ")" << endl;
	//fclose(file);

	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
	}
	delete[] imageData;

	delete[] imageSlice;
	delete[] foundCircle;
	delete[] houghSpace;
	

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

	//======================== Segment Trachea ===========================================
	
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
	
	//======================== Detect carina =============================================
	
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

	//======================== Detect seed from carina ===================================
	
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
	
	//======================== 2D path finding ==================================
	
	/*
	//const size_t height = 256, width = 256, dim2D = height * width;
	const size_t height = 1024, width = 1024, dim2D = height * width;

	dataType* imageData = new dataType[dim2D]{ 0 };
	dataType* actionMap = new dataType[dim2D]{ 0 };
	dataType* potential = new dataType[dim2D]{ 0 };
	dataType* path = new dataType[dim2D]{ 0 };

	//loading_path = inputPath + "artificial.raw";
	//loading_path = inputPath + "raw/slice/eye_1024_1024.raw";
	//loading_path = inputPath + "raw/slice/eye_rescaled.raw";
	loading_path = inputPath + "raw/slice/eye_filtered.raw";
	//loading_path = inputPath + "raw/slice/threshold_eye.raw";
	//loading_path = inputPath + "raw/slice/distanceMap.raw";
	load2dArrayRAW<dataType>(imageData, height, width, loading_path.c_str(), false);

	//storing_path = outputPath + "loaded.raw";
	//store2dRawData<dataType>(imageData, height, width, storing_path.c_str());

	Point2D* seedPoints = new Point2D[2];
	////Point2D seed1 = { 8, 128 }, seed2 = { 131, 66 }; // artificial
	////Point2D seed2 = { 284, 540 }, seed1 = { 22, 470 };
	////Point2D seed2 = { 356, 304 }, seed1 = { 558, 102 }; //small vessel
	
	Point2D seed1 = { 314, 537 }, seed2 = { 509, 824 }; //bigger vessel
	seedPoints[0] = seed1; 
	seedPoints[1] = seed2;

	storing_path = outputPath + "seed.csv";
	FILE* file_seed;
	if (fopen_s(&file_seed, storing_path.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	fprintf(file_seed, "x,y\n");
	fprintf(file_seed, "%f,%f\n", seed1.x, seed1.y);
	fprintf(file_seed, "%f,%f\n", seed1.x, seed1.y + 1);
	fprintf(file_seed, "%f,%f\n", seed1.x, seed1.y - 1);
	fprintf(file_seed, "%f,%f\n", seed1.x - 1, seed1.y - 1);
	fprintf(file_seed, "%f,%f\n", seed1.x - 1, seed1.y + 1);
	fprintf(file_seed, "%f,%f\n", seed1.x - 1, seed1.y);
	fprintf(file_seed, "%f,%f\n", seed1.x + 1, seed1.y);
	fprintf(file_seed, "%f,%f\n", seed1.x + 1, seed1.y - 1);
	fprintf(file_seed, "%f,%f\n", seed1.x + 1, seed1.y + 1);
	fprintf(file_seed, "%f,%f\n", seed2.x, seed2.y);
	fprintf(file_seed, "%f,%f\n", seed2.x, seed2.y + 1);
	fprintf(file_seed, "%f,%f\n", seed2.x, seed2.y - 1);
	fprintf(file_seed, "%f,%f\n", seed2.x - 1, seed2.y - 1);
	fprintf(file_seed, "%f,%f\n", seed2.x - 1, seed2.y + 1);
	fprintf(file_seed, "%f,%f\n", seed2.x - 1, seed2.y);
	fprintf(file_seed, "%f,%f\n", seed2.x + 1, seed2.y - 1);
	fprintf(file_seed, "%f,%f\n", seed2.x + 1, seed2.y + 1);
	fprintf(file_seed, "%f,%f\n", seed2.x + 1, seed2.y);
	fclose(file_seed);

	////generate stat potential
	//dataType radiusInitial = 2.0, radiusMax = 4.0, radiusStep = 1.0;
	//Statistics stats = { 0.0, 0.0, 0.0, 0.0, 0.0 };
	//Point2D origin = { 0.0, 0.0 };
	//PixelSpacing spacing = { 1.0, 1.0 };
	//Image_Data2D ctImageData = { height, width, imageData, origin, spacing, orientation2D };
	
	//stats = get2DPointNeighborhoodStats(ctImageData, seed1, radI);
	//dataType meanI = stats.mean_data / radI, varianceI = stats.variance / radI;
	//for (i = 0; i < height; i++) {
	//	for (j = 0; j < width; j++) {
	//		Point2D current = { i, j };
	//		stats = get2DPointNeighborhoodStats(ctImageData, current, radI);
	//		stats.mean_data = stats.mean_data / radI;
	//		stats.variance = stats.variance / radI;
	//		potential[x_new(i, j, height)] = 0.00001 + 1.0 * pow(meanI - stats.mean_data, 2) + 1.0 * pow(varianceI - stats.variance, 2);
	//	}
	//}

	//double start_time = clock();
	//computePotentialMeanVariance(ctImageData, potential, seed1, radiusInitial, radiusMax, radiusStep);
	//double end_time = clock();
	//std::cout << "The computation lasts : " << (end_time - start_time) / CLOCKS_PER_SEC << " seconds." << std::endl;

	partialFrontPropagation2D(imageData, actionMap, potential, height, width, seedPoints);

	storing_path = outputPath + "action.raw";
	store2dRawData<dataType>(actionMap, height, width, storing_path.c_str());
	
	shortestPath2d(actionMap, path, height, width, 1.0, seedPoints);

	//storing_path = outputPath + "loaded.raw";
	//store2dRawData<dataType>(imageData, height, width, storing_path.c_str());

	//fastSweepingFunction_2D(actionMap, imageData, height, width, 1.0, 1000000.0, 0.0);
	//dataType* gradientVectorX = new dataType[dim2D]{ 0 };
	//dataType* gradientVectorY = new dataType[dim2D]{ 0 };
	//computeImageGradient(imageData, gradientVectorX, gradientVectorY, height, width, 1.0);
	//dataType norm_of_gradient = 0.0, ux = 0.0, uy = 0.0;
	//for (k = 0; k < dim2D; k++) {
	//	ux = gradientVectorX[k];
	//	uy = gradientVectorY[k];
	//	norm_of_gradient = sqrt(ux * ux + uy * uy);
	//	actionMap[k] = gradientFunction(norm_of_gradient, 1000);
	//}
	//for (k = 0; k < dim2D; k++) {
	//	if (actionMap[k] < 0.15) {
	//		actionMap[k] = 1.0;
	//	}
	//	else {
	//		actionMap[k] = 0.0;
	//	}
	//}
	//delete[] gradientVectorX;
	//delete[] gradientVectorY;

	//rescaleNewRange2D(imageData, height, width, 0, 1);
	////filtering parameters
	//Filter_Parameters filter_parameters;
	//filter_parameters.h = 1.0; filter_parameters.timeStepSize = 0.25; filter_parameters.eps2 = 1e-6;
	//filter_parameters.omega_c = 1.5; filter_parameters.tolerance = 1e-3; filter_parameters.maxNumberOfSolverIteration = 100;
	//filter_parameters.timeStepsNum = 10; filter_parameters.coef = 1e-6;
	//Image_Data2D IMAGE;
	//IMAGE.imageDataPtr = imageData;
	//IMAGE.height = height; IMAGE.width = width;
	//heatImplicit2dScheme(IMAGE, filter_parameters);
	//storing_path = outputPath + "filtered.raw";
	//store2dRawData<dataType>(imageData, height, width, storing_path.c_str());

	//The potential is computed inside the fast marching function 
	//fastMarching2D(imageData, actionMap, potential, height, width, seedPoints);
	//double dPoint1 = 0.0, dPoint2 = 0.0;
	//for (i = 0; i < height; i++) {
	//	for (j = 0; j < width; j++) {
	//		Point2D current = { i, j };
	//		dPoint1 = getPoint2DDistance(seed1, current);
	//		dPoint2 = getPoint2DDistance(seed2, current);
	//		if (dPoint1 <= 10) {
	//			potential[x_new(i, j, height)] = 1.0;
	//		}
	//		if (dPoint2 <= 10) {
	//			potential[x_new(i, j, height)] = 1.0;
	//		}
	//	}
	//}

	//const size_t radiusLength = 10;
	//dataType radiusInitial = 3.0, radiusFinal = 8.0, radiusScale = 1.0;

	//dataType** potentialPtr = new dataType * [radiusLength];
	//dataType** actionPtr = new dataType * [radiusLength];
	//dataType** pathPtr = new dataType * [radiusLength];
	//for (k = 0; k < radiusLength; k++) {
	//	potentialPtr[k] = new dataType[dim2D] {0};
	//	actionPtr[k] = new dataType[dim2D]{0};
	//	pathPtr[k] = new dataType[dim2D]{ 0 };
	//}

	//Statistics stats = { 0.0, 0.0, 0.0, 0.0, 0.0 };
	//double radI = 5.0;
	//Point2D origin = { 0.0, 0.0 };
	//PixelSpacing spacing = { 1.0, 1.0 };
	//Image_Data2D ctImageData = { height, width, imageData, origin, spacing, orientation2D };
	//stats = get2DPointNeighborhoodStats(ctImageData, seed1, radI);

	//computePotentialMeanVariance(ctImageData, potentialPtr, radiusLength, seed1, radiusInitial, radiusScale);
	//storing_path = outputPath + "potential3D.raw";
	//store3dRawData<dataType>(potentialPtr, height, width, radiusLength, storing_path.c_str());
	//load3dArrayRAW<dataType>(potentialPtr, height, width, radiusLength, storing_path.c_str(), false);

	//computePotentialMeanVariance(ctImageData, potentialPtr, radiusLength, seed1, radiusInitial, radiusScale);
	//Point3D p1 = { seed1.x, seed1.y, radiusInitial }, p2 = { seed2.x, seed2.y, radiusFinal };
	//Point3D* seeds = new Point3D[2];
	//seeds[0] = p1; seeds[1] = p2;
	//partialFrontPropagation(actionPtr, potentialPtr, height, width, radiusLength, seeds);

	//storing_path = outputPath + "action3D.raw";
	//store3dRawData<dataType>(actionPtr, height, width, radiusLength, storing_path.c_str());
	//load3dArrayRAW<dataType>(actionPtr, height, width, radiusLength, storing_path.c_str(), false);

	//vector<Point3D> path_points;
	//shortestPath3D(actionPtr, pathPtr, height, width, radiusLength, 1.0, seeds, path_points);
	//string saving_name = outputPath + "path2DPlus.csv";
	//FILE* file;
	//if (fopen_s(&file, saving_name.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(file, "x,y,r\n");
	//for (k = 0; k < path_points.size(); k++) {
	//	fprintf(file, "%f,%f,%f\n", path_points[k].x, path_points[k].y, path_points[k].z);
	//}
	//fclose(file);
	
	//for (i = 0; i < 192; i++) {
	//	for (j = 0; j < width; j++) {
	//		if (j == 128) {
	//			imageData[x_new(i, j - 2, height)] = 1.0;
	//			imageData[x_new(i, j - 1, height)] = 1.0;
	//			imageData[x_new(i, j, height)] = 1.0;
	//			imageData[x_new(i, j + 1, height)] = 1.0;
	//			imageData[x_new(i, j + 2, height)] = 1.0;
	//		}
	//	}
	//}
	//size_t i_min = 128, i_max = 192;
	//size_t j_min = 64, j_max = 128;
	//size_t in = 0, jn = 0;
	//for (i = i_min, in = 0; i < i_max; i++, in++) {
	//	for (j = j_min, jn = 0; j < j_max; j++, jn++) {
	//		if (in == jn) {
	//			imageData[x_new(i, j - 2, height)] = 1.0;
	//			imageData[x_new(i, j - 1, height)] = 1.0;
	//			imageData[x_new(i, j, height)] = 1.0;
	//			imageData[x_new(i, j + 1, height)] = 1.0;
	//			imageData[x_new(i, j + 2, height)] = 1.0;
	//		}
	//	}
	//}
	//i_min = 128, i_max = 192;
	//j_min = 128, j_max = 192;
	//for (i = i_max, in = 0; i > i_min; i--, in++) {
	//	for (j = j_min, jn = 0; j < j_max; j++, jn++) {
	//		if (in == jn) {
	//			imageData[x_new(i, j - 2, height)] = 1.0;
	//			imageData[x_new(i, j - 1, height)] = 1.0;
	//			imageData[x_new(i, j, height)] = 1.0;
	//			imageData[x_new(i, j + 1, height)] = 1.0;
	//			imageData[x_new(i, j + 2, height)] = 1.0;
	//		}
	//	}
	//}
	//storing_path = outputPath + "artificial.raw";
	//store2dRawData<dataType>(imageData, height, width, storing_path.c_str());
	
	//for (k = 0; k < radiusLength; k++) {
	//	delete[] potentialPtr[k];
	//	delete[] actionPtr[k];
	//	delete[] pathPtr[k];
	//}
	//delete[] potentialPtr;
	//delete[] actionPtr;
	//delete[] pathPtr;

	//delete[] seeds;
	//delete[] seedPoints;

	delete[] seedPoints;
	delete[] imageData;
	delete[] actionMap;
	delete[] potential;
	delete[] path;
	*/

	//======================== Detection of break points ==============
	
	/*
	//The coordinates of found points are :
	Point3D f1 = { 263.251465, 258.381989, 313.890533 };
	Point3D f2 = { 268.761658, 269.798309, 396.832947 };
	Point3D f3 = { 280.994934, 297.958435, 476.050629 };
	Point3D f4 = { 294.67923, 310.155334, 554.720276 };

	string saving_csv = outputPath + "break_points.csv";
	FILE* file;
	if (fopen_s(&file, saving_csv.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	fprintf(file, "%f,%f,%f\n", f1.x, f1.y, f1.z);
	fprintf(file, "%f,%f,%f\n", f2.x, f2.y, f2.z);
	fprintf(file, "%f,%f,%f\n", f3.x, f3.y, f3.z);
	fprintf(file, "%f,%f,%f\n", f4.x, f4.y, f4.z);
	fclose(file);
	*/
	
	//======================== Free Memory ======================================

	//free(ctContainer);
	
	return EXIT_SUCCESS;
}