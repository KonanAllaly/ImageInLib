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

using namespace std;

int main() {

	string inputPath = "C:/Users/Konan Allaly/Documents/Tests/input/vtk/";
	string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";

	Vtk_File_Info* dataCointainer = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	dataCointainer->operation = copyFrom;

	std::string inputImagePath = inputPath + "Patient6_ct.vtk";
	//readVtkFile(inputImagePath.c_str(), dataCointainer);

	int Height = 406; //dataCointainer->dimensions[2];
	int Length = 512; //dataCointainer->dimensions[1];
	int Width = 512; //dataCointainer->dimensions[0];

	int i, j, k, x, dim2D = Length * Width;

	////Find min max data in order to rescalle data range positive
	//dataType min_data = 100000000.0, max_data = -100000000.0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length * Width; i++) {
	//		if (dataCointainer->dataPointer[k][i] > max_data) {
	//			max_data = dataCointainer->dataPointer[k][i];
	//		}
	//		if (dataCointainer->dataPointer[k][i] < min_data) {
	//			min_data = dataCointainer->dataPointer[k][i];
	//		}
	//	}
	//}
	//cout << "min data = " << min_data << ", max data = " << max_data << endl;

	////Rescalle to positive range
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		dataCointainer->dataPointer[k][i] = dataCointainer->dataPointer[k][i] + (-1) * min_data;
	//	}
	//}

	/*================  Data container for fast marching and path finding ============*/

	dataType ** imageData = new dataType * [Height];
	dataType** distance = new dataType * [Height];
	dataType** path = new dataType * [Height];
	dataType** potential = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D];
		distance[k] = new dataType[dim2D];
		path[k] = new dataType[dim2D];
		potential[k] = new dataType[dim2D];

		if (imageData[k] == NULL || distance[k] == NULL || path[k] == NULL || potential[k] == NULL)
			return false;
	}
	if (imageData == NULL || distance == NULL || path == NULL || potential == NULL)
		return false;

	//Initialization
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			imageData[k][i] = 0.0; //dataCointainer->dataPointer[k][i];
			distance[k][i] = 0.0;
			path[k][i] = 0.0;
			potential[k][i] = 0.0;
		}
	}

	string outputImagePath = outputPath + "loaded_p6.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	inputImagePath = outputPath + "filteredGMC.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, inputImagePath.c_str(), LOAD_DATA, false);

	//rescaleNewRange(imageData, Length, Width, Height, 0.0, 1.0, 4095.0, 0.0);

	/*================== Image Filtering =============================================*/

	Image_Data ImageData; ImageData.imageDataPtr = imageData; ImageData.height = Height;
	ImageData.length = Length; ImageData.width = Width;
	Filter_Parameters filter_parameters; filter_parameters.edge_detector_coefficient = 1;
	filter_parameters.h = 1.0; filter_parameters.timeStepSize = 1.2; filter_parameters.eps2 = 1e-6;
	filter_parameters.omega_c = 1.5; filter_parameters.tolerance = 1e-2; filter_parameters.maxNumberOfSolverIteration = 100;
	filter_parameters.timeStepsNum = 1; filter_parameters.coef = 1e-6; 
	//filterImage(ImageData, filter_parameters, GEODESIC_MEAN_CURVATURE_FILTER);
	//outputImagePath = outputPath + "filteredGMC.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	/*================== Fast marching and path finding ==============================*/

	point3d* seeds = new point3d[2];

	////Patient 2
	seeds[1].x = 293; seeds[1].y = 312; seeds[1].z = 265;
	seeds[0].x = 264; seeds[0].y = 261; seeds[0].z = 170;

	////Patient 3
	//seeds[0].x = 285; seeds[0].y = 306; seeds[0].z = 369;
	//seeds[1].x = 259; seeds[1].y = 258; seeds[1].z = 274;

	////Patient 4
	//seeds[1].x = 292; seeds[1].y = 286; seeds[1].z = 235;
	//seeds[0].x = 265; seeds[0].y = 265; seeds[0].z = 152;

	////Patient 4
	//seeds[1].x = 292; seeds[1].y = 286; seeds[1].z = 235;
	//seeds[0].x = 265; seeds[0].y = 265; seeds[0].z = 152;

	////Patient 6
	//seeds[1].x = 268; seeds[1].y = 246; seeds[1].z = 486;
	//seeds[0].x = 284; seeds[0].y = 276; seeds[0].z = 652;

	//=================== Potential ====================================

	//compute3dPotential(imageData, potential, Length, Width, Height, seeds);

	//====================================================================

	//Point3D imageOrigin = { dataCointainer->origin[0], dataCointainer->origin[1], dataCointainer->origin[2] };
	//VoxelSpacing imageSpacing = { dataCointainer->spacing[0], dataCointainer->spacing[1], dataCointainer->spacing[2] };
	//OrientationMatrix orientation;
	//orientation.v1 = { 1.0, 0.0, 0.0 }; orientation.v2 = { 0.0, 1.0, 0.0 }; orientation.v3 = { 0.0, 0.0, 1.0 };

	//Point3D startPoint = { 258, 254, 251 }, endPoint = { 262, 260, 147 }, currentPoint = {0.0, 0.0, 0.0};
	//startPoint = getRealCoordFromImageCoord3D(startPoint, imageOrigin, imageSpacing, orientation);
	//endPoint = getRealCoordFromImageCoord3D(endPoint, imageOrigin, imageSpacing, orientation);
	//double dist1 = 0.0, dist2 = 0.0, radius = 5.0;

	////Verify if selected point are inside the aorta
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			x = x_new(i, j, Length);
	//			currentPoint.x = i; currentPoint.y = j; currentPoint.z = k;
	//			currentPoint = getRealCoordFromImageCoord3D(currentPoint, imageOrigin, imageSpacing, orientation);
	//			dist1 = getPoint3DDistance(startPoint, currentPoint);
	//			dist2 = getPoint3DDistance(endPoint, currentPoint);
	//			if (dist1 <= radius) {
	//				path[k][x] = 1.0;
	//			}
	//			if (dist2 <= radius) {
	//				path[k][x] = 1.0;
	//			}
	//		}
	//	}
	//}
	//outputImagePath = outputPath + "startAndEndPoints.raw";
	//manageRAWFile3D<dataType>(path, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	fastMarching3D_N(imageData, distance, potential, Length, Width, Height, seeds);

	outputImagePath = outputPath + "distance.raw";
	//manageRAWFile3D<dataType>(distance, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	outputImagePath = outputPath + "potential.raw";
	manageRAWFile3D<dataType>(potential, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	shortestPath3d(distance, path, Length, Width, Height, 1.0, seeds);

	outputImagePath = outputPath + "ResultPath.raw";
	manageRAWFile3D<dataType>(path, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	delete[] seeds;
	free(dataCointainer);
	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
		delete[] distance[k];
		delete[] path[k];
		delete[] potential[k];
	}
	delete[] imageData;
	delete[] distance;
	delete[] path;
	delete[] potential;

	return EXIT_SUCCESS;
}