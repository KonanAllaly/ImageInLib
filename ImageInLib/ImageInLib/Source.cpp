#pragma warning(disable : 6011)
#pragma warning(disable : 4244)
#pragma warning(disable : 4477)
#pragma warning(disable : 4267)
#pragma warning(disable : 4700)

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string.h> 
#include <math.h>
#include <time.h>

#include "file.h"
#include "common_functions.h"
#include "Labelling.h"
#include "template_functions.h"
#include "../src/vtk_params.h"
#include "common_vtk.h"
#include "src/imageInterpolation.h"

#define _sx 0.976562
#define _sy 0.976562
#define _sz 2.5


int main() {

	size_t i, j, k;

	//image Dimensions
	size_t Width = 512;
	size_t Length = 512;
	const size_t dim2D = Width * Length;
	const size_t Height = 330;

	std::string inputPath = "C:/Users/Konan Allaly/Documents/Tests/input/";
	std::string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";

	//std::string inputImagePath = inputPath + "patient2_pet.raw";
	std::string inputImagePath = inputPath + "patient4_pet.vtk";

	//Operation operation = LOAD_DATA;
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, inputImagePath.c_str(), operation, false);
	//std::string  outputImagePath = outputPath + "load_PET.raw";
	//operation = STORE_DATA;
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, outputImagePath.c_str(), operation, false);

	Vtk_File_Info* vtkInfo = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));

	//double spacing_pet[3] = {4.0, 4.0, 4.0};
	//double origin_pet[3] = {286.586, 216.586, -1338.80};
	//int dim[3] = { Length, Width, Height };
	//VtkDataType DType = dta_Dbl;
	//vtkOperation ope = copyFrom;

	//typedef struct {
	//	double spacing[3];
	//	double origin[3];
	//	int dimensions[3]; // length, width, height
	//	VtkDataType vDataType;
	//	dataType** dataPointer;
	//	vtkOperation operation;
	//} Vtk_File_Info;
	//vtkInfo = { spacing_pet, origin_pet, dim, DType, imageData, ope };

	//vtkInfo->dataPointer = imageData;
	//vtkInfo->dimensions[0] = Length; vtkInfo->dimensions[1] = Width; vtkInfo->dimensions[2] = Height;
	//vtkInfo->spacing[0] = sx_pet; vtkInfo->spacing[1] = sy_pet; vtkInfo->spacing[2] = sz_pet;
	//vtkInfo->origin[0] = 250.0; vtkInfo->origin[1] = 250.0; vtkInfo->origin[2] = -878.5;
	//vtkInfo->vDataType = dta_Dbl;
	//vtkInfo->operation = copyFrom;

	readVtkFile(inputImagePath.c_str(), vtkInfo);

	//std::string  outputImagePath = outputPath + "load_PET.vtk";
	//Operation operation = STORE_DATA;
	//manageRAWFile3D<dataType>(vtkInfo->dataPointer, Length, Width, Height, outputImagePath.c_str(), operation, false);
	//storeVtkFile(outputImagePath.c_str(), vtkInfo, dta_binary);

	Image_Data image1;
	image1.imageDataPtr = vtkInfo->dataPointer;
	image1.origin.x = vtkInfo->origin[0]; image1.origin.y = vtkInfo->origin[1]; image1.origin.z = vtkInfo->origin[2];
	image1.height = vtkInfo->dimensions[2]; image1.length = vtkInfo->dimensions[1]; image1.width = vtkInfo->dimensions[0];
	image1.spacing.sx = (float)vtkInfo->spacing[0]; image1.spacing.sy = (float)vtkInfo->spacing[1]; image1.spacing.sz = (float)vtkInfo->spacing[2];
	image1.orientation.v1 = { 1, 0, 0 }; image1.orientation.v2 = { 0, 1, 0 }; image1.orientation.v3 = { 0, 0, -1 };

	dataType** imageData = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D];
	}
	if (imageData == NULL)
		return false;

	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			imageData[k][i] = 0.0;
		}
	}

	Image_Data image2;
	image2.imageDataPtr = imageData;
	image2.origin = {250.0, 250.0, -876.50 };
	image2.spacing = { _sx, _sy, _sz };
	image2.orientation.v1 = { 1, 0, 0 }; image2.orientation.v2 = {0, 1, 0}; image2.orientation.v3 = { 0, 0, -1 };
	image2.height = Height; image2.width = Width; image2.length = Length;

	interpolationMethod method = TRILINEAR;

	imageInterpolation3D(image1, image2, method);

	Operation operation = STORE_DATA;
	std::string  outputImagePath = outputPath + "interpolated.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, outputImagePath.c_str(), operation, false);
	
	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
	}
	delete[] imageData;

	free(vtkInfo);

	return EXIT_SUCCESS;
}
