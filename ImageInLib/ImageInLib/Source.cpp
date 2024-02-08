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

	//std::string inputImagePath = inputPath + "vtk/ct/Patient6_ct.vtk";
	//std::string inputImagePath = inputPath + "vtk/lungs/cropped_ct_lungs_p6.vtk";
	std::string inputImagePath = inputPath + "vtk/lungs/box_lungs_ct_p6.vtk";
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

	//define new array to load input images
	dataType** imageData = new dataType*[Height];
	dataType** lungsData = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{0};
		lungsData[k] = new dataType[dim2D]{0};
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

	////rescale to positive range
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		//ctContainer->dataPointer[k][i] = ctContainer->dataPointer[k][i] + (-1) * min_data;
	//		imageData[k][i] = ctContainer->dataPointer[k][i] + (-1) * min_data;
	//	}
	//}

	//Copy to array
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			imageData[k][i] = ctContainer->dataPointer[k][i];
			lungsData[k][i] = ctContainer->dataPointer[k][i];
		}
	}

	//thresholdingOTSU(lungsData, Length, Width, Height, 0.0, 1.0);

	//inputImagePath = inputPath + "segment.raw";
	//load3dArrayRAW<dataType>(lungsData, Length, Width, Height, inputImagePath.c_str(), false);

	//inputImagePath = inputPath + "raw/lungs/lungs_p6.raw";
	//load3dArrayRAW<dataType>(lungsData, Length, Width, Height, inputImagePath.c_str(), false);

	//Image_Data ctImageData;
	//ctImageData.imageDataPtr = lungsData;
	//ctImageData.height = Height; ctImageData.length = Length; ctImageData.width = Width;
	//ctImageData.orientation = orientation; ctImageData.origin = ctOrigin; ctImageData.spacing = ctSpacing;

	//======================== Slice Cropping ================================================================

	//imageMetaData croppedLungs;
	//const size_t offset = 5;
	//croppedLungs =  croppImage3D(ctImageData, offset);

	//size_t length = croppedLungs.length, width = croppedLungs.width, height = croppedLungs.height;
	//std::cout << "Cropped lungs image dim : " << length << " x " << width << " x " << height << "" << std::endl;
	//Point3D originLungs = croppedLungs.origin;
	//std::cout << "Cropped Lungs origin : (" << originLungs.x << ", " << originLungs.y << ", " << originLungs.z << ")" << std::endl;
	//VoxelSpacing lungsSpacing = croppedLungs.spacing;

	//dataType** cropped_ct = new dataType*[height];
	//dataType** cropped_lungs = new dataType * [height];
	//for (k = 0; k < height; k++) {
	//	cropped_ct[k] = new dataType[dim2D]{0};
	//	cropped_lungs[k] = new dataType[dim2D]{0};
	//	if (cropped_ct[k] == NULL || cropped_lungs[k] == NULL)
	//		return false;
	//}
	//if (cropped_ct == NULL || cropped_lungs == NULL)
	//	return false;

	////===
	////Image_Data Cropped;
	////Cropped.imageDataPtr = croppedCT;
	////Cropped.height = height; Cropped.length = length; Cropped.width = width;
	////Cropped.orientation = orientation; Cropped.origin = croppedOrigin; Cropped.spacing = croppedSpacing;
	////===

	//Point3D origin_cropped = getImageCoordFromRealCoord3D(originLungs, ctOrigin, ctSpacing, orientation);
	//size_t i0 = (size_t)origin_cropped.x, j0 = (size_t)origin_cropped.y, k0 = (size_t)origin_cropped.z;
	//
	//Point3D point_test = { 0.0, 0.0, k0 };
	//point_test = getRealCoordFromImageCoord3D(point_test, ctOrigin, ctSpacing, orientation);
	//cout << "Origin cropped image : (" << point_test.x << "," << point_test.y << "," << point_test.z << ")" << endl;
	//
	//size_t i1, j1, k1;
	//for (k = 0, k1 = k0; k < height; k++, k1++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			cropped_ct[k][x_new(i, j, Length)] = imageData[k1][x_new(i, j, Length)];
	//			cropped_lungs[k][x_new(i, j, Length)] = lungsData[k1][x_new(i, j, Length)];
	//		}
	//	}
	//}

	//saving_name = outputPath + "cropped_ct_p6.raw";
	////store3dRawData<dataType>(cropped_ct, Length, Width, height, saving_name.c_str());

	//saving_name = outputPath + "cropped_lungs_p6.raw";
	////store3dRawData<dataType>(cropped_lungs, Length, Width, height, saving_name.c_str());

	//Vtk_File_Info* Container = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	//Container->dataPointer = cropped_lungs; //cropped_ct;
	//Container->dimensions[0] = Length; Container->dimensions[1] = Width; Container->dimensions[2] = height;
	//Container->origin[0] = point_test.x; Container->origin[1] = point_test.y; Container->origin[2] = point_test.z;
	//Container->spacing[0] = ctSpacing.sx; Container->spacing[1] = ctSpacing.sy; Container->spacing[2] = ctSpacing.sz;
	//Container->operation = copyTo; Container->vDataType = dta_Flt;

	//saving_name = outputPath + "cropped_lungs_p6.vtk";
	//vtkDataForm data_form = dta_binary;
	//storeVtkFile(saving_name.c_str(), Container, data_form);

	//free(Container);
	//for (k = 0; k < height; k++) {
	//	delete[] cropped_ct[k];
	//	delete[] cropped_lungs[k];
	//}
	//delete[] cropped_ct;
	//delete[] cropped_lungs;

	//======================== Cropp box ====================================

	//imageMetaData croppedLungs;
	//const size_t offset = 5;
	//croppedLungs =  croppImage3D(ctImageData, offset);

	//size_t length = croppedLungs.length, width = croppedLungs.width, height = croppedLungs.height;
	//std::cout << "Cropped lungs image dim : " << length << " x " << width << " x " << height << "" << std::endl;
	//Point3D originLungs = croppedLungs.origin;
	//std::cout << "Cropped Lungs origin : (" << originLungs.x << ", " << originLungs.y << ", " << originLungs.z << ")" << std::endl;
	//VoxelSpacing lungsSpacing = croppedLungs.spacing;

	//dataType** cropped_ct = new dataType*[height];
	//dataType** cropped_lungs = new dataType * [height];
	//for (k = 0; k < height; k++) {
	//	cropped_ct[k] = new dataType[dim2D]{0};
	//	cropped_lungs[k] = new dataType[dim2D]{0};
	//	if (cropped_ct[k] == NULL || cropped_lungs[k] == NULL)
	//		return false;
	//}
	//if (cropped_ct == NULL || cropped_lungs == NULL)
	//	return false;

	//Point3D origin_cropped = getImageCoordFromRealCoord3D(originLungs, ctOrigin, ctSpacing, orientation);
	//size_t i0 = (size_t)origin_cropped.x, j0 = (size_t)origin_cropped.y, k0 = (size_t)origin_cropped.z;
	//
	//Point3D point_test = { 0.0, 0.0, 0.0 };
	//point_test = getRealCoordFromImageCoord3D(point_test, ctOrigin, ctSpacing, orientation);
	//cout << "Origin cropped image : (" << point_test.x << "," << point_test.y << "," << point_test.z << ")" << endl;
	//
	//size_t i1, j1, k1;
	//for (k = 0, k1 = k0; k < height; k++, k1++) {
	//	for (i = 0, i1 = i0; i < length; i++, i1++) {
	//		for (j = 0, j1 = j0; j < width; j++, j1++) {
	//			cropped_ct[k][x_new(i, j, length)] = imageData[k1][x_new(i1, j1, Length)];
	//			cropped_lungs[k][x_new(i, j, length)] = lungsData[k1][x_new(i1, j1, Length)];
	//		}
	//	}
	//}

	//Vtk_File_Info* Container = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	//Container->dataPointer = cropped_lungs; // cropped_ct; //
	//Container->dimensions[0] = length; Container->dimensions[1] = width; Container->dimensions[2] = height;
	//Container->origin[0] = originLungs.x; Container->origin[1] = originLungs.y; Container->origin[2] = originLungs.z;
	//Container->spacing[0] = ctSpacing.sx; Container->spacing[1] = ctSpacing.sy; Container->spacing[2] = ctSpacing.sz;
	//Container->operation = copyTo; Container->vDataType = dta_Flt;

	//saving_name = outputPath + "box_lungs_p6.vtk";
	//vtkDataForm data_form = dta_binary;
	//storeVtkFile(saving_name.c_str(), Container, data_form);

	//free(Container);
	//for (k = 0; k < height; k++) {
	//	delete[] cropped_ct[k];
	//	delete[] cropped_lungs[k];
	//}
	//delete[] cropped_ct;
	//delete[] cropped_lungs;

	//======================== Extraction of Heart from CT data =================

	dataType** heartData = new dataType * [Height];
	dataType** maskThreshold = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		heartData[k] = new dataType[dim2D]{0};
		maskThreshold[k] = new dataType[dim2D]{0};
	}

	////Manual segmentation
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		if (lungsData[k][i] == 1.0) {
	//			maskThreshold[k][i] = imageData[k][i];
	//		}
	//		else {
	//			maskThreshold[k][i] = 0.0;
	//		}
	//	}
	//}

	thresholding3dFunctionN(imageData, Length, Width, Height, 100, 700, 0.0, 1.0);
	//erosion3dHeighteenNeigbours(imageData, Length, Width, Height, 1.0, 0.0);
	//dilatation3dHeighteenNeigbours(imageData, Length, Width, Height, 1.0, 0.0);
	 
	saving_name = outputPath + "threshold_lungs.raw";
	store3dRawData<dataType>(imageData, Length, Width, Height, saving_name.c_str());
	//store3dRawData<dataType>(maskThreshold, Length, Width, Height, saving_name.c_str());
	//store3dRawData<dataType>(lungsData, Length, Width, Height, saving_name.c_str());
	//free(histo);

	for (k = 0; k < Height; k++) {
		delete[] heartData[k];
		delete[] maskThreshold[k];
	}
	delete[] heartData;
	delete[] maskThreshold;

	//======================== Free Memory ======================================

	free(ctContainer);
	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
		delete[] lungsData[k];
	}
	delete[] imageData;
	delete[] lungsData;

	return EXIT_SUCCESS;
}