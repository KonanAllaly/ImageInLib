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

using namespace std;

int main() {

	string inputPath = "C:/Users/Konan Allaly/Documents/Tests/input/vtk/";
	string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";

	OrientationMatrix orientation;
	orientation.v1 = { 1.0, 0.0, 0.0 }; orientation.v2 = { 0.0, 1.0, 0.0 }; orientation.v3 = { 0.0, 0.0, 1.0 };

	Vtk_File_Info* ctContainer = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	ctContainer->operation = copyFrom;

	std::string inputImagePath = inputPath + "Patient2_ct.vtk";
	readVtkFile(inputImagePath.c_str(), ctContainer);

	int Height = ctContainer->dimensions[2];
	int Length = ctContainer->dimensions[1];
	int Width = ctContainer->dimensions[0];

	int i, j, k, x, dim2D = Length * Width;

	//Find min max data in order to rescalle data range positive
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
	////cout << "min data = " << min_data << ", max data = " << max_data << endl;

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
	
	/*================  Data container for fast marching and path finding ============*/

	dataType ** imageData = new dataType * [Height];
	dataType** imageFiltered = new dataType * [Height];
	dataType** distance = new dataType * [Height];
	dataType** path = new dataType * [Height];
	dataType** potential = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D];
		distance[k] = new dataType[dim2D];
		path[k] = new dataType[dim2D];
		potential[k] = new dataType[dim2D];
		imageFiltered[k] = new dataType[dim2D];

		if (imageData[k] == NULL || distance[k] == NULL || path[k] == NULL || potential[k] == NULL || imageFiltered[k] == NULL)
			return false;
	}
	if (imageData == NULL || distance == NULL || path == NULL || potential == NULL || imageFiltered == NULL)
		return false;

	//Initialization
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			imageData[k][i] = 0;//ctContainer->dataPointer[k][i];
			distance[k][i] = 0.0;
			path[k][i] = 0.0;
			potential[k][i] = 0.0;
			imageFiltered[k][i] = 0.0;
		}
	}

	//inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/raw/patient2_ct_rescaled.raw";
	//inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/Im_sd_r4.raw";
	inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/filteredGMC.raw";
	manageRAWFile3D<dataType>(imageFiltered, Length, Width, Height, inputImagePath.c_str(), LOAD_DATA, false);

	inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/raw/patient2_ct_rescaled.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, inputImagePath.c_str(), LOAD_DATA, false);

	/*================== Image Filtering =============================================*/

	//Image_Data ImageData; ImageData.imageDataPtr = imageData; ImageData.height = Height;
	//ImageData.length = Length; ImageData.width = Width;
	//Filter_Parameters filter_parameters; filter_parameters.edge_detector_coefficient = 1;
	//filter_parameters.h = 1.0; filter_parameters.timeStepSize = 1.2; filter_parameters.eps2 = 1e-6;
	//filter_parameters.omega_c = 1.5; filter_parameters.tolerance = 1e-2; filter_parameters.maxNumberOfSolverIteration = 100;
	//filter_parameters.timeStepsNum = 1; filter_parameters.coef = 1e-6; 
	//filterImage(ImageData, filter_parameters, GEODESIC_MEAN_CURVATURE_FILTER);
	//string outputImagePath = outputPath + "filteredGMC.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	/*================== Fast marching and path finding ==============================*/

	//point3d* seeds = new point3d[2];
	//Patient 2
	//seeds[1].x = 262; seeds[1].y = 254; seeds[1].z = 250; // top
	
	//seeds[0].x = 256; seeds[0].y = 251; seeds[0].z = 251;
	//seeds[0].x = 256; seeds[0].y = 249; seeds[0].z = 254;
	//seeds[1].x = 293; seeds[1].y = 312; seeds[1].z = 265;
	//down
	//seeds[1].x = 264; seeds[1].y = 261; seeds[1].z = 170;
	//seeds[0].x = 295; seeds[0].y = 314; seeds[0].z = 265;

	//seeds[1].x = 284; seeds[1].y = 285; seeds[1].z = 278; // middle
	//seeds[0].x = 265; seeds[0].y = 259; seeds[0].z = 143; // down

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

	//=================== Interpolation ====================================

	Vtk_File_Info* petContainer = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	petContainer->operation = copyFrom;

	inputImagePath = inputPath + "Patient2_pet.vtk";
	//readVtkFile(inputImagePath.c_str(), petContainer);

	int height = petContainer->dimensions[2];
	int length = petContainer->dimensions[1];
	int width = petContainer->dimensions[0];

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
	////rescaleNewRange(petContainer->dataPointer, length, width, height, 0.0, 1.0, max_pet, min_pet);
	
	//rescaleToIntervalZeroOne(petContainer->dataPointer, length, width, height);

	Point3D ctOrigin = { ctContainer->origin[0], ctContainer->origin[1], ctContainer->origin[2] };
	VoxelSpacing ctSpacing = { ctContainer->spacing[0], ctContainer->spacing[1], ctContainer->spacing[2] };

	Image_Data CT; CT.imageDataPtr = imageData;
	CT.height = Height; CT.length = Length; CT.width = Width;
	CT.orientation = orientation; CT.origin = ctOrigin; CT.spacing = ctSpacing;

	//Point3D petOrigin = { petContainer->origin[0], petContainer->origin[1], petContainer->origin[2] };
	//VoxelSpacing petSpacing = { petContainer->spacing[0], petContainer->spacing[1], petContainer->spacing[2] };

	//Image_Data PET; PET.imageDataPtr = petContainer->dataPointer; 
	//PET.height = height; PET.length = length; PET.width = width;
	//PET.orientation = orientation; PET.origin = petOrigin; PET.spacing = petSpacing;

	//================== New Potential ====================================
	dataType max_pot = 0.0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			potential[k][i] = abs(imageData[k][i] - imageFiltered[k][i]);
			if (potential[k][i] > max_pot) {
				max_pot = potential[k][i];
			}
		}
	}

	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			potential[k][i] = 1.0 / (0.01 + abs(potential[k][i] - max_pot));
		}
	}

	string outputImagePath = outputPath + "new_potential.raw";
	manageRAWFile3D<dataType>(potential, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	//================== Fast Marching and Path finding ===================

	dataType** image_max = new dataType * [Height];
	dataType** image_min = new dataType * [Height];
	dataType** image_mean = new dataType * [Height];
	dataType** image_sd = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		image_max[k] = new dataType[dim2D];
		image_min[k] = new dataType[dim2D];
		image_mean[k] = new dataType[dim2D];
		image_sd[k] = new dataType[dim2D];
	}

	//Initialization
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			image_max[k][i] = 0.0;
			image_min[k][i] = 0.0;
			image_mean[k][i] = 0.0;
			image_sd[k][i] = 0.0;
		}
	}

	statictics_Pointers statsImage; statsImage.maximum = image_max; statsImage.minimum = image_min;
	statsImage.mean = image_mean; statsImage.sd = image_sd;

	//inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/Im_mean_r5.raw";
	//manageRAWFile3D<dataType>(image_mean, Length, Width, Height, inputImagePath.c_str(), LOAD_DATA, false);

	//inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/Im_max_r5.raw";
	//manageRAWFile3D<dataType>(image_max, Length, Width, Height, inputImagePath.c_str(), LOAD_DATA, false);

	//inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/Im_min_r5.raw";
	//manageRAWFile3D<dataType>(image_min, Length, Width, Height, inputImagePath.c_str(), LOAD_DATA, false);

	//inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/Im_sd_r5.raw";
	//manageRAWFile3D<dataType>(image_sd, Length, Width, Height, inputImagePath.c_str(), LOAD_DATA, false);

	Point3D point1 = { 262, 254, 250 }; // top 
	//Point3D point_new1 = { 263, 259, 170 }; // point new 1
	//Point3D point_new2 = { 267, 269, 192 }; // point new 2
	//Point3D point_new3 = { 272, 277, 205 }; // point new 3
	//Point3D point_new4 = { 296, 317, 250 }; // point new 4
	Point3D point2 = { 263, 257, 146 }; // bottom
	//Point3D point2 = { 285, 284, 278 }; // middle
	Point3D* seedPoints = new Point3D[2];
	//seedPoints[0] = point2; 
	seedPoints[0] = point2;
	seedPoints[1] = point1;

	Potential_Parameters parameters;
	parameters.K = 0.005; parameters.epsilon = 0.01;
	parameters.c_ct = 0.0; parameters.c_pet = 0.0;
	parameters.c_max = 0.0; parameters.c_min = 0.0; parameters.c_mean = 1.0; parameters.c_sd = 1.0;

	//compute3dPotential(imageData, potential, Length, Width, Height, seedPoints);

	//computePotential_N(CT, PET, statsImage, potential, seedPoints, 5.0, parameters);

	//newPotential(CT, potential, seedPoints, 3.0);

	//potentialOnEdgeImage(imageData, potential, seedPoints, Length, Width, Height);

	fastMarching3D_N(imageData, distance, potential, Length, Width, Height, seedPoints);

	//string outputImagePath = outputPath + "potential_eucl.raw";
	//manageRAWFile3D<dataType>(potential, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	outputImagePath = outputPath + "distance.raw";
	manageRAWFile3D<dataType>(distance, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	//inputImagePath = "C:/Users/Konan Allaly/Documents/Tests/input/distance_ct_mean.raw";
	//manageRAWFile3D<dataType>(distance, Length, Width, Height, inputImagePath.c_str(), LOAD_DATA, false);

	shortestPath3d(distance, path, Length, Width, Height, 1.0, seedPoints);

	outputImagePath = outputPath + "path.raw";
	manageRAWFile3D<dataType>(path, Length, Width, Height, outputImagePath.c_str(), STORE_DATA, false);

	//ctContainer->operation = copyTo;
	//ctContainer->dataPointer = path;
	//outputImagePath = outputPath + "path.vtk";
	//storeVtkFile(outputImagePath.c_str(), ctContainer, dta_binary);

	//========================= Free memory ============================

	for (k = 0; k < Height; k++) {
		delete[] image_max[k];
		delete[] image_min[k];
		delete[] image_mean[k];
		delete[] image_sd[k];
	}
	delete[] image_max;
	delete[] image_min;
	delete[] image_mean;
	delete[] image_sd;

	delete[] seedPoints;
	free(ctContainer);
	free(petContainer);
	
	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
		delete[] distance[k];
		delete[] path[k];
		delete[] potential[k];
		delete[] imageFiltered[k];
	}
	delete[] imageData;
	delete[] distance;
	delete[] path;
	delete[] potential;
	delete[] imageFiltered;

	return EXIT_SUCCESS;
}