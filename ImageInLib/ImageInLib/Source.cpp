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

using namespace std;

int main() {

	string inputPath = "C:/Users/Konan Allaly/Documents/Tests/input/";
	string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";
	
	string loading_name, saving_name, extension;

	size_t i, j, k;

	//===================== Load 3D patient data (.vtk) ================================

	OrientationMatrix orientation;
	orientation.v1 = { 1.0, 0.0, 0.0 }; orientation.v2 = { 0.0, 1.0, 0.0 }; orientation.v3 = { 0.0, 0.0, 1.0 };

	Vtk_File_Info* ctContainer = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	ctContainer->operation = copyFrom;

	std::string inputImagePath = inputPath + "vtk/ct/Patient6_ct.vtk";
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

	//define new array to load mean and filtered image
	dataType** imageData = new dataType*[Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{0};
	}

	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			imageData[k][i] = ctContainer->dataPointer[k][i];
		}
	}


	//======================== Cropping ================================================================

	//dataType** croppedCT = new dataType*[height];
	//for (k = 0; k < height; k++) {
	//	croppedCT[k] = new dataType[dim2D]{0};
	//	if (croppedCT[k] == NULL)
	//		return false;
	//}
	//if (croppedCT == NULL)
	//	return false;

	//===
	//Image_Data Cropped;
	//Cropped.imageDataPtr = croppedCT;
	//Cropped.height = height; Cropped.length = length; Cropped.width = width;
	//Cropped.orientation = orientation; Cropped.origin = croppedOrigin; Cropped.spacing = croppedSpacing;
	//===

	//Point3D origin_cropped = getImageCoordFromRealCoord3D(croppedOrigin, ctOrigin, ctSpacing, orientation);
	//size_t i0 = (size_t)origin_cropped.x, j0 = (size_t)origin_cropped.y, k0 = (size_t)origin_cropped.z;
	//Point3D point_test = { 0.0, 0.0, k0 };
	//point_test = getRealCoordFromImageCoord3D(point_test, ctOrigin, ctSpacing, orientation);
	//cout << "Origin cropped image : (" << point_test.x << "," << point_test.y << "," << point_test.z << ")" << endl;
	//
	//size_t i1, j1, k1;
	//for (k = 0, k1 = k0; k < height; k++, k1++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			croppedCT[k][x_new(i, j, Length)] = ctContainer->dataPointer[k1][x_new(i, j, Length)];
	//			//croppedCT[k][x_new(i, j, Length)] = imageData[k1][x_new(i, j, Length)];
	//			//croppedCT[k][x_new(j, i, width)] = ctContainer->dataPointer[k1][x_new(j1, i1, Width)];
	//		}
	//	}
	//}
	//saving_name = outputPath + "mean_cropped_p6.raw";
	//store3dRawData<dataType>(croppedCT, Length, Width, height, saving_name.c_str());

	//Vtk_File_Info* Container = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	//Container->dataPointer = croppedCT;
	//Container->dimensions[0] = Length; Container->dimensions[1] = Width; Container->dimensions[2] = height;
	//Container->origin[0] = point_test.x; Container->origin[1] = point_test.y; Container->origin[2] = point_test.z;
	//Container->spacing[0] = croppedContainer->spacing[0]; Container->spacing[1] = croppedContainer->spacing[1]; Container->spacing[2] = croppedContainer->spacing[2];
	//Container->operation = copyTo; Container->vDataType = dta_Flt;

	//saving_name = outputPath + "cropped_ct_p6.vtk";
	//vtkDataForm data_form = dta_binary;
	//storeVtkFile(saving_name.c_str(), Container, data_form);

	//for (k = 0; k < height; k++) {
	//	delete[] croppedCT[k];
	//}
	//delete[] croppedCT;


	//======================== Free Memory ======================================

	free(ctContainer);
	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
	}
	delete[] imageData;

	return EXIT_SUCCESS;
}