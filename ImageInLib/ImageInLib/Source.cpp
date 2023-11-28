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

#define slice_num 146

//#define thres_min 2000.0 //p1, p10
//#define thres_max 2200.0 //p1, p10
//#define thres_min 995.0 //p2, p3, p6, p7
//#define thres_max 1213.0 //p2, p3, p6, p7
////#define thres_min 3000.0 //p4, p5, p9
////#define thres_max 3200.0 //p4, p5, p9

#define level_0 900.0
#define level_1 1300.0

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

	std::string inputImagePath = inputPath + "vtk/Patient5_ct.vtk";
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
	dataType** maskThresh = new dataType * [Height];
	dataType** distanceMap = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{0};
		maskThresh[k] = new dataType[dim2D]{0};
		distanceMap[k] = new dataType[dim2D]{0};
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
	//		maskThresh[k][i] = imageData[k][i];
	//	}
	//}

	//max_data = (-1) * min_data + max_data;
	//std::cout << "new min = 0, new max = " << max_data << std::endl;
	//rescaleNewRange(imageData, Length, Width, Height, 0.0, 1.0, 4000.0, 0.0);

	//inputImagePath = inputPath + "raw/filteredGMC_p2.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, inputImagePath.c_str(), LOAD_DATA, false);

	Image_Data CT; CT.imageDataPtr = imageData; CT.height = Height; CT.length = Length; CT.width = Width;
	CT.origin = ctOrigin; CT.spacing = ctSpacing; CT.orientation = orientation;

	//==================== Window lovel for visualization =============================

	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		if (imageData[k][i] <= level_0) {
	//			imageData[k][i] = level_0;
	//		}
	//		if (imageData[k][i] >= level_1) {
	//			imageData[k][i] = level_1;
	//		}
	//	}
	//}

	//std::string outputImagePath = outputPath + "window_p2.raw";
	//store3dRawData<dataType>(imageData, Length, Width, Height, outputImagePath.c_str());

	//outputImagePath = outputPath + "loaded_p2.raw";
	//store3dRawData<dataType>(maskThresh, Length, Width, Height, outputImagePath.c_str());

	//==================== Fast Sweeping to locate the liver ==========================
	
	////Threshold
	//thresholding3dFunctionN(maskThresh, Length, Width, Height, thres_min, thres_max, 0.0, 1.0);
	//
	//string outputImagePath = outputPath + "threshold.raw";
	//store3dRawData<dataType>(maskThresh, Length, Width, Height, outputImagePath.c_str());

	////distance map by fast sweeping
	//dataType h = 1.0, large_value = 100000000.0;
	//fastSweepingFunction_3D(distanceMap, maskThresh, Length, Width, Height, h, large_value, 0.0);
	//outputImagePath = outputPath + "distance_map.raw";
	////store3dRawData<dataType>(distanceMap, Length, Width, Height, outputImagePath.c_str());

	////Point inside the liver
	//Point3D point_liver = { 0.0, 0.0, 0.0 }, current_point = { 0.0, 0.0, 0.0 };
	//point_liver = getPointWithTheHighestValue(distanceMap, Length, Width, Height);
	//std::cout << "Point find inside the liver : (" << point_liver.x << ", " << point_liver.y << ", " << point_liver.z << ")" << std::endl;
	//point_liver = getRealCoordFromImageCoord3D(point_liver, ctOrigin, ctSpacing, orientation);
	//double dist_point = 0.0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			xd = x_new(i, j, Length);
	//			current_point.x = i;
	//			current_point.y = j;
	//			current_point.z = k;
	//			current_point = getRealCoordFromImageCoord3D(current_point, ctOrigin, ctSpacing, orientation);
	//			dist_point = getPoint3DDistance(point_liver, current_point);
	//			if (dist_point <= 10.0) {
	//				maskThresh[k][xd] = 1.0;
	//			}
	//			else {
	//				maskThresh[k][xd] = 0.0;
	//			}
	//		}
	//	}
	//}
	//outputImagePath = outputPath + "ball_inside_liver.raw";
	//store3dRawData<dataType>(maskThresh, Length, Width, Height, outputImagePath.c_str());

	//==================== Segment the liver by labeling ================================

	////Erosion in order to isolate the liver
	//erosion3dHeighteenNeigbours(maskThresh, Length, Width, Height, 1.0, 0.0);
	//erosion3dHeighteenNeigbours(maskThresh, Length, Width, Height, 1.0, 0.0);
	//erosion3dHeighteenNeigbours(maskThresh, Length, Width, Height, 1.0, 0.0);
	//erosion3dHeighteenNeigbours(maskThresh, Length, Width, Height, 1.0, 0.0);
	//erosion3dHeighteenNeigbours(maskThresh, Length, Width, Height, 1.0, 0.0);
	//erosion3dHeighteenNeigbours(maskThresh, Length, Width, Height, 1.0, 0.0);
	//erosion3dHeighteenNeigbours(maskThresh, Length, Width, Height, 1.0, 0.0);
	//erosion3dHeighteenNeigbours(maskThresh, Length, Width, Height, 1.0, 0.0);
	//erosion3dHeighteenNeigbours(maskThresh, Length, Width, Height, 1.0, 0.0);
	//
	////outputImagePath = outputPath + "erode_six_times.raw";
	////store3dRawData<dataType>(maskThresh, Length, Width, Height, outputImagePath.c_str());

	////Labeling to find the biggest region
	//int** segment = new int* [Height];
	//bool** status = new bool* [Height];
	//for (k = 0; k < Height; k++) {
	//	segment[k] = new int[dim2D]{0};
	//	status[k] = new bool[dim2D]{0};
	//}
	//labelling3D(maskThresh, segment, status, Length, Width, Height, 1.0);

	////Number of regions
	//int numberOfRegionsCells = 0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		if (maskThresh[k][i] == 1) {
	//			numberOfRegionsCells++;
	//		}
	//	}
	//}
	//
	////Counting regions
	//int* counting = new int[numberOfRegionsCells] {0};
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		if (segment[k][i] > 0) {
	//			counting[segment[k][i]]++;
	//		}
	//	}
	//}

	////Find the largest region 
	//int maxVolume = 0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		if(counting[segment[k][i]] > maxVolume){
	//			maxVolume = counting[segment[k][i]];
	//		}
	//	}
	//}

	////Save the largest region
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		if (counting[segment[k][i]] != maxVolume) {
	//			maskThresh[k][i] = 0;
	//		}
	//		else {
	//			maskThresh[k][i] = 1;
	//		}
	//	}
	//}
	//
	//string outputImagePath = outputPath + "largest_region.raw";
	////store3dRawData<dataType>(maskThresh, Length, Width, Height, outputImagePath.c_str());

	//================================= Interpolation ======================================
	
	//int interp_Height = (int)((Height * ctSpacing.sz) / ctSpacing.sx);
	//dataType** interp_image = new dataType * [interp_Height];
	//for (k = 0; k < interp_Height; k++) {
	//	interp_image[k] = new dataType[dim2D]{0};
	//}
	//Image_Data interpolated; interpolated.imageDataPtr = interp_image;
	//interpolated.height = interp_Height; interpolated.length = Length; interpolated.width = Width;
	//interpolated.orientation = orientation; interpolated.origin = ctOrigin;
	//interpolated.spacing.sx = ctSpacing.sx; interpolated.spacing.sy = ctSpacing.sy; interpolated.spacing.sz = ctSpacing.sx;
	//cout << "interpolated image dim : " << ctContainer->dimensions[0] << " x " << ctContainer->dimensions[1] << " x " << interp_Height << "" << endl;

	//imageInterpolation3D(CT, interpolated, NEAREST_NEIGHBOR);

	//================================= Segmentation ========================================

	////Down sampling for segmentation
	//int length = (int)(Length / 4), width = (int)(Width / 4), height = (int)(interp_Height / 4);
	//dataType** down_image = new dataType * [height];
	//dataType** initial_seg = new dataType * [height];
	//for (k = 0; k < height; k++) {
	//	down_image[k] = new dataType[length * width]{0};
	//	initial_seg[k] = new dataType[length * width]{0};
	//}
	//cout << "down sampled image dim : " << length << " x " << width << " x " << height << "" << endl;

	//Image_Data imageToBeSegmeneted;
	//imageToBeSegmeneted.height = height; imageToBeSegmeneted.length = length; imageToBeSegmeneted.width = width;
	//imageToBeSegmeneted.orientation = orientation; imageToBeSegmeneted.origin = ctOrigin;
	//imageToBeSegmeneted.spacing.sx = ctSpacing.sx * 4; imageToBeSegmeneted.spacing.sy = ctSpacing.sy * 4; imageToBeSegmeneted.spacing.sz = interpolated.spacing.sz * 4;
	//cout << "down sampled image spacing : (" << imageToBeSegmeneted.spacing.sx << ", " << imageToBeSegmeneted.spacing.sy << ", " << imageToBeSegmeneted.spacing.sz << ")" << endl;

	////fill the initial segment
	//CT.imageDataPtr = maskThresh; // contain the segmented liver
	//interpolated.imageDataPtr = interp_image;
	//// interpolate to have same spacing in each direction
	//imageInterpolation3D(CT, interpolated, NEAREST_NEIGHBOR);
	//
	//// down sample the segmented liver
	//imageToBeSegmeneted.imageDataPtr = initial_seg;
	//imageInterpolation3D(interpolated, imageToBeSegmeneted, NEAREST_NEIGHBOR);

	////rescal to 0-1 in order to use the current parameters
	//rescaleNewRange(imageData, Length, Width, Height, 0.0, 1.0, max_data, min_data);
	//
	////interpolate the input image in order to have same spacing in each direction
	//CT.imageDataPtr = imageData; // contain the input image rescaled
	//interpolated.imageDataPtr = interp_image;
	//imageInterpolation3D(CT, interpolated, NEAREST_NEIGHBOR);

	//// down sample the input image
	//imageToBeSegmeneted.imageDataPtr = down_image;
	//imageInterpolation3D(interpolated, imageToBeSegmeneted, NEAREST_NEIGHBOR);

	//Point3D* center_seg = new Point3D[1];
	//center_seg->x = 0.0, center_seg->y = 0.0; center_seg->z = 0.0;

	//Segmentation_Parameters seg_params;
	//seg_params.tau = 1.2; seg_params.h = 1.0; seg_params.omega_c = 1.5; seg_params.mod = 10;
	//seg_params.maxNoGSIteration = 40; seg_params.coef = 100000; seg_params.gauss_seidelTolerance = 1e-6;
	//seg_params.maxNoOfTimeSteps = 1000; seg_params.segTolerance = 1e-6; seg_params.eps2 = 1e-6;
	//seg_params.coef_conv = 1.0; seg_params.coef_dif = 1.0;

	////Smoothing by heat equation so we don't need all the parameters
	//Filter_Parameters smooth_params;
	//smooth_params.h = seg_params.h; smooth_params.timeStepSize = smooth_params.h / 4.0; smooth_params.timeStepsNum = 1; smooth_params.omega_c = 1.5;
	//smooth_params.tolerance = 1e-6, smooth_params.maxNumberOfSolverIteration = 100;

	//unsigned char segmentPath[] = "C:/Users/Konan Allaly/Documents/Tests/output/segmentation/";
	//generalizedSubsurfSegmentation(imageToBeSegmeneted, initial_seg, seg_params, smooth_params, center_seg, 1, segmentPath);

	//====================== Save the information needed in the tests ======================
	
	//FILE* info_test;
	//saving_name = outputPath + "_info_test_done.txt";
	//if (fopen_s(&info_test, saving_name.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(info_test, "patient 10\n");
	//fprintf(info_test, "Image origin : ( %f, %f, %f )\n", ctOrigin.x, ctOrigin.y, ctOrigin.z);
	//fprintf(info_test, "Length = %d, Width = %d, Height = %d\n", Length, Width, Height);
	//fprintf(info_test, "spacing : sx = %f, sy = %f, sz = %f\n", ctSpacing.sx, ctSpacing.sy, ctSpacing.sz);
	//fprintf(info_test, "Original data range : min = %f, max = %f\n", min_data, max_data + min_data);
	//fprintf(info_test, "New data range : min = %f, max = %f\n", min_data - min_data, max_data);
	//fprintf(info_test, "Threshold to find the liver : thres min = %f , thres max = %f\n", thres_min, thres_max);
	////fprintf(info_test, "8 Erosions were needed to isolate the liver\n");
	////fprintf(info_test, "Smooting parameters : tau = %f, h = %f, omega = %f, tolerance = %f, max gauss seidel iteration = %d\n", smooth_params.timeStepSize, smooth_params.h, smooth_params.omega_c, smooth_params.tolerance, smooth_params.maxNumberOfSolverIteration);
	////fprintf(info_test, "Segmentation parameters : tau = %f, h = %f, omega = %f, gauss seidel tolerance = %f, max gauss seidel iteration = %d, number of time steps = %d, edge detector coef = %f, coef conv = %f, coef dif = %f\n", 
	//	//seg_params.tau, seg_params.h, seg_params.omega_c, seg_params.gauss_seidelTolerance, seg_params.maxNoGSIteration, seg_params.maxNoOfTimeSteps, seg_params.coef, seg_params.coef_conv, seg_params.coef_dif);
	//fclose(info_test);

	//===================== Data containers for local Hough transform ==================
	
	int height = (int)((Height * ctSpacing.sz) / ctSpacing.sx);
	dataType** image_interp = new dataType * [height];
	dataType** path = new dataType * [height];
	for (k = 0; k < height; k++) {
		image_interp[k] = new dataType[dim2D]{0};
		path[k] = new dataType[dim2D]{0};
	}

	inputImagePath = inputPath + "raw/filteredGMC_p5.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, inputImagePath.c_str(), LOAD_DATA, false);

	Image_Data interpolated; interpolated.imageDataPtr = image_interp;
	interpolated.height = height; interpolated.length = Length; interpolated.width = Width;
	interpolated.orientation = orientation; interpolated.origin = ctOrigin;
	interpolated.spacing.sx = ctSpacing.sx; interpolated.spacing.sy = ctSpacing.sy; interpolated.spacing.sz = ctSpacing.sx;
	cout << "interpolated image dim : " << ctContainer->dimensions[0] << " x " << ctContainer->dimensions[1] << " x " << height << "" << endl;
	imageInterpolation3D(CT, interpolated, NEAREST_NEIGHBOR);

	inputImagePath = inputPath + "raw/path_p5.raw";
	manageRAWFile3D<dataType>(path, Length, Width, height, inputImagePath.c_str(), LOAD_DATA, false);

	Point3D seed = { 0,0,0 };
	vector<Point3D> list_of_path_point;

	for (k = 0; k < height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				xd = x_new(j, i, Width);
				if (path[k][xd] == 1.0) {
					seed.x = j;
					seed.y = i;
					seed.z = k;
					list_of_path_point.push_back(seed);
				}
			}
		}
	}

	int number_of_point = list_of_path_point.size();

	//===================== 2D Hough Transform to detect circle ========================
	
	size_t count = 0;

	dataType* imageSlice = new dataType[dim2D]{0};
	dataType* houghSpace = new dataType[dim2D]{0};
	dataType* voteArray = new dataType[dim2D]{0};
	dataType* gradientX = new dataType[dim2D]{0};
	dataType* gradientY = new dataType[dim2D]{0};
	dataType* edgeDetector = new dataType[dim2D]{0};

	size_t in, kn;
	dataType norm_gradient = 0.0, K = 1000.0;
	size_t step = 1, step_of_max = 0;
	dataType max_ratio = 0, current_ratio = 0;
	double radius = 7.0, radius_of_max = 0.0;
	Point2D seed2D = { 0.0, 0.0 };
	dataType thres = 0.15;
	for (in = 0; in < number_of_point; in++) {
		
		kn = (size_t)list_of_path_point[in].z;
		
		//Slice extraction
		for (i = 0; i < dim2D; i++) {
			imageSlice[i] = image_interp[kn][i];
		}

		//Edge detector
		computeImageGradient(imageSlice, gradientX, gradientY, Length, Width, 1.0);
		for (i = 0; i < dim2D; i++) {
			norm_gradient = sqrt(gradientX[i] * gradientX[i] + gradientY[i] * gradientY[i]);
			edgeDetector[i] = gradientFunction(norm_gradient, K);
		}

		//manual threshold
		for (i = 0; i < dim2D; i++) {
			if (edgeDetector[i] < thres) {
				edgeDetector[i] = 1.0;
			}
			else {
				edgeDetector[i] = 0.0;
			}
		}

		//saving for paraview
		extension = to_string(in);
		if (in < 10) {
			saving_name = outputPath + "slice_000" + extension + ".raw";
			store2dRawData(imageSlice, Length, Width, saving_name.c_str());
			saving_name = outputPath + "_threshold_slice_000" + extension + ".raw";
			store2dRawData(edgeDetector, Length, Width, saving_name.c_str());
		}
		else {
			if (in < 100) {
				saving_name = outputPath + "slice_00" + extension + ".raw";
				store2dRawData(imageSlice, Length, Width, saving_name.c_str());
				saving_name = outputPath + "_threshold_slice_00" + extension + ".raw";
				store2dRawData(edgeDetector, Length, Width, saving_name.c_str());
			}
			else {
				if (in < 1000) {
					saving_name = outputPath + "slice_0" + extension + ".raw";
					store2dRawData(imageSlice, Length, Width, saving_name.c_str());
					saving_name = outputPath + "_threshold_slice_0" + extension + ".raw";
					store2dRawData(edgeDetector, Length, Width, saving_name.c_str());
				}
				else {
					saving_name = outputPath + "slice_" + extension + ".raw";
					store2dRawData(imageSlice, Length, Width, saving_name.c_str());
					saving_name = outputPath + "_threshold_slice_" + extension + ".raw";
					store2dRawData(edgeDetector, Length, Width, saving_name.c_str());
				}
			}
		}

		//local circular Hough transform
		step = 1, step_of_max = 0;
		max_ratio = 0.0, current_ratio = 0.0;
		radius = 7.0, radius_of_max = 0.0;
		seed2D.x = list_of_path_point[in].x;
		seed2D.y = list_of_path_point[in].y;

		while (radius <= 15.0) {
			localCircularHoughTransform(seed2D, edgeDetector, houghSpace, voteArray, Length, Width, radius);
			current_ratio = getTheMaxValue(voteArray, Length, Width);
			if (current_ratio > max_ratio) {
				max_ratio = current_ratio;
				radius_of_max = radius;
				step_of_max = step;
				
				//saving for paraview
				if (in < 10) {
					saving_name = outputPath + "found_circle_local_000" + extension + ".raw";
					store2dRawData(houghSpace, Length, Width, saving_name.c_str());
				}
				else {
					if (in < 100) {
						saving_name = outputPath + "found_circle_local_00" + extension + ".raw";
						store2dRawData(houghSpace, Length, Width, saving_name.c_str());
					}
					else {
						if (in < 1000) {
							saving_name = outputPath + "found_circle_local_0" + extension + ".raw";
							store2dRawData(houghSpace, Length, Width, saving_name.c_str());
						}
						else {
							saving_name = outputPath + "found_circle_local_" + extension + ".raw";
							store2dRawData(houghSpace, Length, Width, saving_name.c_str());
						}
					}
				}
			}
			radius += 0.5;
			step++;
		}
		std::cout << "The maximal ratio is on step " << step_of_max << " and the radius is : " << radius_of_max << std::endl;

	}

	////Slice extraction
	//for (i = 0; i < dim2D; i++) {
	//	imageSlice[i] = imageData[slice_num - 1][i];
	//}

	//extension = to_string(slice_num);
	//saving_name = outputPath + "_loaded_slice_" + extension + "_p2.raw";
	//store2dRawData(imageSlice, Length, Width, saving_name.c_str());

	//rescaleNewRange2D(imageSlice, Length, Width, 0.0, 1.0);

	////smoothing parameters
	//Filter_Parameters filter_parameters;
	//filter_parameters.h = 1.0; filter_parameters.timeStepSize = 0.25; filter_parameters.eps2 = 1e-6;
	//filter_parameters.omega_c = 1.5; filter_parameters.tolerance = 1e-3; filter_parameters.maxNumberOfSolverIteration = 100;
	//filter_parameters.timeStepsNum = 5; filter_parameters.coef = 1e-6; 

	//Image_Data2D IMAGE; IMAGE.imageDataPtr = imageSlice;
	//IMAGE.height = Length; IMAGE.width = Width;
	//heatImplicit2dScheme(IMAGE, filter_parameters);

	//saving_name = outputPath + "_filtered_slice" + extension + "_p2.raw";
	//store2dRawData(imageSlice, Length, Width, saving_name.c_str());

	////Edge detector
	//computeImageGradient(imageSlice, gradientX, gradientY, Length, Width, 1.0);
	//dataType norm_gradient = 0.0, K = 1000.0;
	//for (i = 0; i < dim2D; i++) {
	//	norm_gradient = sqrt(gradientX[i] * gradientX[i] + gradientY[i] * gradientY[i]);
	//	edgeDetector[i] = gradientFunction(norm_gradient, K);
	//}

	//saving_name = outputPath + "_edge_detector_slice" + extension + "_p2.raw";
	//store2dRawData(edgeDetector, Length, Width, saving_name.c_str());

	////manual threshold
	//dataType thres = 0.15;
	//for (i = 0; i < dim2D; i++) {
	//	if (edgeDetector[i] < thres) {
	//		edgeDetector[i] = 1.0;
	//	}
	//	else {
	//		edgeDetector[i] = 0.0;
	//	}
	//}

	//saving_name = outputPath + "_threshold_slice" + extension + "_p2.raw";
	//store2dRawData(edgeDetector, Length, Width, saving_name.c_str());

	//size_t step = 1, step_of_max = 0;
	//dataType max_ratio = 0, current_ratio = 0;
	
	//extension = to_string(step);
	//while (radius <= 15.0) {
	//	//circularHoughTransform(edgeDetector, houghSpace, voteArray, Length, Width, radius);
	//	localCircularHoughTransform(seed, edgeDetector, houghSpace, voteArray, Length, Width, radius);
	//	current_ratio = getTheMaxValue(voteArray, Length, Width);
	//	if (current_ratio > max_ratio) {
	//		max_ratio = current_ratio;
	//		radius_of_max = radius;
	//		step_of_max = step;
	//		saving_name = outputPath + "max_cumulative_ratio_local.raw";
	//		store2dRawData(voteArray, Length, Width, saving_name.c_str());
	//		saving_name = outputPath + "found_circle_local.raw";
	//		store2dRawData(houghSpace, Length, Width, saving_name.c_str());
	//	}
	//	radius += 0.5;
	//	step++;
	//}
	//std::cout << "The maximal ratio is on step " << step_of_max << " and the radius is : " << radius_of_max << std::endl;

	//void localCircularHoughTransform(Point2D seed, dataType* imageDataPtr, dataType* houghSpacePtr, dataType* votingArray, const size_t length, const size_t width, double radius);
	
	////FILE* slice_info;
	////extension = to_string(slice_num);
	////saving_name = outputPath + "_info_slice_" + extension + ".txt";
	////if (fopen_s(&slice_info, saving_name.c_str(), "w") != 0) {
	////	printf("Enable to open");
	////	return false;
	////}
	////fprintf(slice_info, "patient 9\n");
	////fprintf(slice_info, "slice  %s\n", extension);
	//////fprintf(slice_info, "Test done without smoothing\n");
	////fprintf(slice_info, "Smoothed by implicit heat\n");
	////fprintf(slice_info, "tau = %f, h = %f, omega = %f\n", filter_parameters.timeStepSize, filter_parameters.h, filter_parameters.omega_c);
	////fprintf(slice_info, "edge detector coef = %f\n", K);
	////fprintf(slice_info, "threshold = %f\n", thres);
	////fprintf(slice_info, "The maximal ratio is obtain for radius : %f\n", radius_of_max);
	////fprintf(slice_info, "Radius range  [7,15] with step 0.5 were considered.\n");
	////fprintf(slice_info, "Neighbors ratio : < 0.1\n");
	////fclose(slice_info);

	delete[] imageSlice;
	delete[] houghSpace;
	delete[] voteArray;
	delete[] gradientX;
	delete[] gradientY;
	delete[] edgeDetector;

	for (k = 0; k < height; k++) {
		delete[] image_interp[k];
		delete[] path[k];
	}
	delete[] image_interp;
	delete[] path;

	//=================== Fast marching and path finding ==============================

	////patient 2
	//Point3D final_point = { 262, 254, 250 }; // top 
	//Point3D initial_point = { 263, 257, 146 }; // bottom

	//dataType** image_mean = new dataType * [interp_Height];
	//dataType** path = new dataType * [interp_Height];
	//for (k = 0; k < interp_Height; k++) {
	//	image_mean[k] = new dataType[dim2D]{0};
	//	path[k] = new dataType[dim2D]{0};
	//}

	//inputImagePath = inputPath + "/interpolated/patient2/Im_mean_r5.raw";
	//manageRAWFile3D<dataType>(image_mean, Length, Width, interp_Height, inputImagePath.c_str(), LOAD_DATA, false);

	////Get corresponding point in interpolated image
	//initial_point = getRealCoordFromImageCoord3D(initial_point, ctOrigin, ctSpacing, orientation);
	//initial_point = getImageCoordFromRealCoord3D(initial_point, ctOrigin, interpolated.spacing, orientation);
	//initial_point.x = (size_t)initial_point.x; initial_point.y = (size_t)initial_point.y; initial_point.z = (size_t)initial_point.z;

	////Get corresponding point in interpolated image
	//final_point = getRealCoordFromImageCoord3D(final_point, ctOrigin, ctSpacing, orientation);
	//final_point = getImageCoordFromRealCoord3D(final_point, ctOrigin, interpolated.spacing, orientation);
	//final_point.x = (size_t)final_point.x; final_point.y = (size_t)final_point.y; final_point.z = (size_t)final_point.z;

	//Point3D* seedPoints = new Point3D[2]; 
	//seedPoints[0] = initial_point;
	//seedPoints[1] = final_point;

	//Potential_Parameters parameters;
	//parameters.K = 0.005; parameters.epsilon = 0.01;
	//parameters.c_ct = 1.0; parameters.c_pet = 0.0;
	//parameters.c_max = 0.0; parameters.c_min = 0.0; parameters.c_mean = 1.0; parameters.c_sd = 0.0;

	//findPathFromOneGivenPoint(interpolated, image_mean, path, seedPoints, parameters);
	//std::string outputImagePath = outputPath + "path_p2.raw";
	//store3dRawData<dataType>(path, Length, Width, interp_Height, outputImagePath.c_str());

	//=================== Free memory =================================================

	//delete[] center_seg;
	//for (k = 0; k < height; k++) {
	//	delete[] down_image[k];
	//	delete[] initial_seg[k];
	//}
	//delete[] down_image;
	//delete[] initial_seg;

	//delete[] counting;
	//for (k = 0; k < Height; k++) {
	//	delete[] segment[k];
	//	delete[] status[k];
	//}
	//delete[] segment;
	//delete[] status;

	free(ctContainer);
	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
		delete[] maskThresh[k];
		delete[] distanceMap[k];
	}
	delete[] imageData;
	delete[] maskThresh;
	delete[] distanceMap;

	//delete[] seedPoints;
	//for (k = 0; k < interp_Height; k++) {
	//	delete[] interp_image[k];
	//	delete[] image_mean[k];
	//	delete[] path[k];
	//}
	//delete[] interp_image;
	//delete[] image_mean;
	//delete[] path;

	return EXIT_SUCCESS;
}