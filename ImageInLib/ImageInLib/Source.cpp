#include <iostream>
#include <sstream>  
#include <vector>
#include <string.h> 
#include <time.h>
#include <cmath>

#include "common_math.h"
#include <template_functions.h>
#include "../src/vtk_params.h"
#include "common_vtk.h"
#include "distanceForPathFinding.h"
#include "common_functions.h"
#include "../src/distance_function.h"

#include "enhancement.h"
#include "eigen_systems.h"
#include "percentile.h"
#include "../src/thresholding.h"

#include "../src/heat_equation.h"
#include "segmentation2d.h"
#include "../src/segmentation3d_gsubsurf.h" 
#include "../src/non_linear_heat_equation.h"

#include "gaussian_distribution.h"

#include "labelling.h"
#include "morphological_change.h"

#define MAX_LINE_LENGTH 1024
#define epsilon 1e-6 

int main() {

	string root = "C:/Users/Konan Allaly/Documents/Tests/";
	
	string loading_path, storing_path, extension;

	size_t i = 0, j = 0, k = 0, xd = 0;

	//===================== Load 3D patient data (.vtk) ==============================================
	
	/*
	We load the .vtk files containing the original pixels value
	and informations related the dimensions, the spacings and
	the origins. We also define the orientation matrix in 2D/3D
	needed when we need to perform interpolation.
	*/
	
	OrientationMatrix orientation = { { 1.0, 0.0, 0.0 } , { 0.0, 1.0, 0.0 } , { 0.0, 0.0, 1.0 } };
	
	Vtk_File_Info* ctContainer = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	ctContainer->operation = copyFrom;

	const Filter_Parameters smoothParameters
	{
		1.0,// timeStepSize;
		1,// h;
		1.0,// sigma;
		100,// edge detector coefficient;
		1.4,// omega_c;
		1e-3,// tolerance;
		1e-6,// eps2;
		1e-6,// coef;
		1,// p;
		1,// timeStepsNum;
		1000// maxNumberOfSolverIteration;
	};

	Segmentation_Parameters segParameters =
	{
		30,//Maximum number of Gauss-Seidel iterations
		100000,//edge detector coef
		1e-6,//epsilon is the regularization factor (Evans-Spruck)
		5001,//Number of current time step
		5001,//Maximum number of time step
		10,//saving frequency
		1e-6,//segmentation tolerance
		1.0,//tau
		1,//h
		1.4,//omega_c
		1e-6,//tolerance
		1.0,//convection coef
		0.05,//diffusion coef
	};

	size_t Length = 0, Width = 0, Height = 0, dim2D;
	size_t length = 0, width = 0, height = 0, dim2d;
	size_t half_length = 0, half_width = 0, half_height = 0;
	size_t k_min = 0, k_max = 0, i_min = 0, i_max = 0, j_min = 0, j_max = 0;
	size_t k_f = 0, i_f = 0, j_f = 0, xd_f = 0;
	dataType max_d = 0;
	Point3D imageOrigin;
	VoxelSpacing imageSpacing, downSpacing;
	Image_Data imageToSegment, imageToSegmentDown, imageDataDistance;

	
	//================== Patient 6 ============================================
	
	std::cout << "Patient 1" << std::endl;
	loading_path = root + "input/vtk/petct/ct/Patient6_ct.vtk";
	readVtkFile(loading_path.c_str(), ctContainer);

	Height = (size_t)ctContainer->dimensions[2];
	Length = (size_t)ctContainer->dimensions[0];
	Width = (size_t)ctContainer->dimensions[1];
	dim2D = Length * Width;
	std::cout << "CT image dim : " << ctContainer->dimensions[0] << " x " << ctContainer->dimensions[1] << " x " << ctContainer->dimensions[2] << "" << std::endl;

	std::cout << "CT origin : (" << ctContainer->origin[0] << ", " << ctContainer->origin[1] << ", " << ctContainer->origin[2] << ")" << std::endl;
	imageOrigin = { (dataType)ctContainer->origin[0], (dataType)ctContainer->origin[1], (dataType)ctContainer->origin[2] };
	imageSpacing = { (dataType)ctContainer->spacing[0], (dataType)ctContainer->spacing[1], (dataType)ctContainer->spacing[2] };
	std::cout << "CT spacing : (" << ctContainer->spacing[0] << ", " << ctContainer->spacing[1] << ", " << ctContainer->spacing[2] << ")" << std::endl; 

	dataType** imageDataFull = new dataType * [Height];
	dataType** initialsegmentFull = new dataType * [Height];
	for (k = 0; k < Height; k++)
	{
		imageDataFull[k] = new dataType[dim2D]{ 0 };
		initialsegmentFull[k] = new dataType[dim2D]{ 0 };
	}
	//loading_path = root + "input/raw/filtered/filtered_p6.raw";
	//manageRAWFile3D<dataType>(imageDataFull, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	loading_path = root + "input/raw/liver/liver_p6.raw";
	manageRAWFile3D<dataType>(initialsegmentFull, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	////Cropping
	////P1
	//k_min = 150, k_max = 240;
	//i_min = 70, i_max = 370;
	//j_min = 125, j_max = 425;

	////P2
	//k_min = 250; k_max = 340;
	//i_min = 90; i_max = 390;
	//j_min = 130; j_max = 430;

	////P3
	//k_min = 115, k_max = 205;
	//i_min = 100, i_max = 400;
	//j_min = 130, j_max = 430;

	////P4
	//k_min = 120, k_max = 200;
	//i_min = 120, i_max = 400;
	//j_min = 130, j_max = 410;

	////P5
	//k_min = 455, k_max = 595;
	//i_min = 130, i_max = 380;
	//j_min = 140, j_max = 390;

	//P6
	k_min = 280, k_max = 410;
	i_min = 90, i_max = 390;
	j_min = 200, j_max = 500;

	length = i_max - i_min;
	width = j_max - j_min;
	height = k_max - k_min;
	dim2d = length * width;

	std::cout << "Cropped Height : " << height << ", Length : " << length << ", Width : " << width << std::endl;

	dataType** imageData = new dataType * [height];
	dataType** initialsegment = new dataType * [height];
	dataType** distanceMap = new dataType * [height];
	for (k = 0; k < height; k++)
	{
		imageData[k] = new dataType[dim2d]{ 0 };
		initialsegment[k] = new dataType[dim2d]{ 0 };
		distanceMap[k] = new dataType[dim2d]{ 0 };
	}

	//k_f = 0, i_f = 0, j_f = 0;
	//for(k = 0, k_f = k_min; k < height; k++, k_f++)
	//{
	//	for (i = 0, i_f = i_min; i < length; i++, i_f++)
	//	{
	//		for(j = 0, j_f = j_min; j < width; j++, j_f++)
	//		{
	//			xd = x_new(i, j, length);
	//			size_t xd_f = x_new(i_f, j_f, Length);
	//			imageData[k][xd] = imageDataFull[k_f][xd_f];
	//			initialsegment[k][xd] = initialsegmentFull[k_f][xd_f];
	//		}
	//	}
	//}

	//dilatation3dHeighteenNeigbours(initialsegment, length, width, height, 1, 0);
	//dilatation3dHeighteenNeigbours(initialsegment, length, width, height, 1, 0);
	//dilatation3dHeighteenNeigbours(initialsegment, length, width, height, 1, 0);

	//storing_path = root + "output/cropped_liver_p6.raw";
	//manageRAWFile3D<dataType>(initialsegment, length, width, height, storing_path.c_str(), STORE_DATA, false);

	//storing_path = root + "output/cropped_image_p6.raw";
	//manageRAWFile3D<dataType>(imageData, length, width, height, storing_path.c_str(), STORE_DATA, false);

	Point3D croppedOrigin = { i_min, j_min, k_min };
	croppedOrigin = getRealCoordFromImageCoord3D(croppedOrigin, imageOrigin, imageSpacing, orientation);
	std::cout << "Cropped image origin : (" << croppedOrigin.x << ", " << croppedOrigin.y << ", " << croppedOrigin.z << ")" << std::endl;

	imageDataDistance = { height, length, width, initialsegment, croppedOrigin, imageSpacing, orientation };
	//fastMarchingDistanceMap(imageDataDistance, distanceMap, 0.0);
	//storing_path = root + "output/distance_map_liver_p6.raw";
	//manageRAWFile3D<dataType>(distanceMap, length, width, height, storing_path.c_str(), STORE_DATA, false);

	////Copy
	//max_d = 0.0;
	//for (k = 0; k < height; k++)
	//{
	//	for (i = 0; i < dim2d; i++)
	//	{
	//		initialsegment[k][i] = distanceMap[k][i];
	//		if (max_d < distanceMap[k][i])
	//		{
	//			max_d = distanceMap[k][i];
	//		}
	//	}
	//}
	//rescaleNewRange(initialsegment, length, width, height, 0.0, max_d, 0.0, 1.0);

	// DownSampling
	half_length = length / 2;
	half_width = width / 2;
	half_height = height / 2;
	downSpacing = { (dataType)(ctContainer->spacing[0] * 2.0), (dataType)(ctContainer->spacing[1] * 2.0), (dataType)(ctContainer->spacing[2] * 2.0) };
	std::cout << "CT spacing : (" << downSpacing.sx << ", " << downSpacing.sy << ", " << downSpacing.sz << ")" << std::endl;
	std::cout << "Down Sample Height : " << half_height << ", Length : " << half_length << ", Width : " << half_width << std::endl;

	dataType** imageDataDown = new dataType * [half_height];
	dataType** initialsegmentDown = new dataType * [half_height];
	for (k = 0; k < half_height; k++)
	{
		imageDataDown[k] = new dataType[half_length * half_width]{ 0 };
		initialsegmentDown[k] = new dataType[half_length * half_width]{ 0 };
	}

	//imageToSegment = { height, length, width, initialsegment, croppedOrigin, imageSpacing, orientation };
	//imageToSegmentDown = { half_height, half_length, half_width, initialsegmentDown, croppedOrigin, downSpacing, orientation };
	//imageInterpolation3D(imageToSegment, imageToSegmentDown, NEAREST_NEIGHBOR);

	//imageToSegment.imageDataPtr = imageData;
	//imageToSegmentDown.imageDataPtr = imageDataDown;
	//imageInterpolation3D(imageToSegment, imageToSegmentDown, NEAREST_NEIGHBOR);

	//storing_path = root + "output/Segmentation/Liver/P6/";
	//GSUBSURF(imageToSegmentDown, initialsegmentDown, storing_path.c_str(), smoothParameters, segParameters);

	storing_path = root + "output/Segmentation/Liver/P6/_seg_func_2000.raw";
	manageRAWFile3D<dataType>(initialsegment, half_length, half_width, half_height, storing_path.c_str(), LOAD_DATA, false);

	for(k = 0; k < half_height; k++)
	{
		for (i = 0; i < half_length * half_width; i++)
		{
			if(initialsegment[k][i] > 0.015)
			{
				initialsegment[k][i] = 1.0;
			}
			else
			{
				initialsegment[k][i] = 0.0;
			}
		}
	}

	storing_path = root + "output/segment_test_p6.raw";
	manageRAWFile3D<dataType>(initialsegment, half_length, half_width, half_height, storing_path.c_str(), STORE_DATA, false);

	dataType* centroid_rough = new dataType[3];
	dataType* centroid_refined = new dataType[3];

	centroidImage(initialsegmentFull, centroid_rough, Height, Length, Width, 0.0);
	centroidImage(initialsegment, centroid_refined, half_height, half_length, half_width, 0.0);

	Point3D centroid_rough_p = { centroid_rough[0], centroid_rough[1], centroid_rough[2] };
	Point3D centroid_refined_p = { centroid_refined[0], centroid_refined[1], centroid_refined[2] };

	centroid_rough_p = getRealCoordFromImageCoord3D(centroid_rough_p, imageOrigin, imageSpacing, orientation);
	centroid_refined_p = getRealCoordFromImageCoord3D(centroid_refined_p, croppedOrigin, downSpacing, orientation);

	double distance_centroid = getPoint3DDistance(centroid_rough_p, centroid_refined_p);
	std::cout << "Distance between rough and refined centroid : " << distance_centroid << std::endl;

	delete[] centroid_rough;
	delete[] centroid_refined;

	for(k = 0; k < half_height; k++)
	{
		delete[] imageDataDown[k];
		delete[] initialsegmentDown[k];
	}
	delete[] imageDataDown;
	delete[] initialsegmentDown;

	for(k = 0; k < height; k++)
	{
		delete[] imageData[k];
		delete[] initialsegment[k];
	}
	delete[] imageData;	
	delete[] initialsegment;

	for (k = 0; k < Height; k++)
	{
		delete[] initialsegmentFull[k];
		delete[] imageDataFull[k];
	}
	delete[] imageDataFull;
	delete[] initialsegmentFull;
	
	
	/*
	//====================== Patient 2 ============================================

	std::cout << "Patient 2" << std::endl;
	loading_path = root + "input/vtk/petct/ct/Patient2_ct.vtk";
	readVtkFile(loading_path.c_str(), ctContainer);

	Height = (size_t)ctContainer->dimensions[2];
	Length = (size_t)ctContainer->dimensions[0];
	Width = (size_t)ctContainer->dimensions[1];
	dim2D = Length * Width;
	std::cout << "CT image dim : " << ctContainer->dimensions[0] << " x " << ctContainer->dimensions[1] << " x " << ctContainer->dimensions[2] << "" << std::endl;

	std::cout << "CT origin : (" << ctContainer->origin[0] << ", " << ctContainer->origin[1] << ", " << ctContainer->origin[2] << ")" << std::endl;
	imageOrigin = { (dataType)ctContainer->origin[0], (dataType)ctContainer->origin[1], (dataType)ctContainer->origin[2] };
	imageSpacing = { (dataType)ctContainer->spacing[0], (dataType)ctContainer->spacing[1], (dataType)ctContainer->spacing[2] };
	std::cout << "CT spacing : (" << ctContainer->spacing[0] << ", " << ctContainer->spacing[1] << ", " << ctContainer->spacing[2] << ")" << std::endl;

	dataType** imageDataFull = new dataType * [Height];
	dataType** initialsegmentFull = new dataType * [Height];
	for (k = 0; k < Height; k++)
	{
		imageDataFull[k] = new dataType[dim2D]{ 0 };
		initialsegmentFull[k] = new dataType[dim2D]{ 0 };
	}
	loading_path = root + "input/raw/filtered/filtered_p2.raw";
	manageRAWFile3D<dataType>(imageDataFull, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	loading_path = root + "input/raw/liver/liver_p2.raw";
	manageRAWFile3D<dataType>(initialsegmentFull, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	////Cropping

	//P2
	k_min = 250, k_max = 340;
	i_min = 90, i_max = 390;
	j_min = 130, j_max = 430;

	length = i_max - i_min;
	width = j_max - j_min;
	height = k_max - k_min;
	dim2d = length * width;

	std::cout << "Cropped Height : " << height << ", Length : " << length << ", Width : " << width << std::endl;

	dataType** imageData = new dataType * [height];
	dataType** initialsegment = new dataType * [height];
	dataType** distanceMap = new dataType * [height];
	for (k = 0; k < height; k++)
	{
		imageData[k] = new dataType[dim2d]{ 0 };
		initialsegment[k] = new dataType[dim2d]{ 0 };
		distanceMap[k] = new dataType[dim2d]{ 0 };
	}

	k_f = 0, i_f = 0, j_f = 0;
	for (k = 0, k_f = k_min; k < height; k++, k_f++)
	{
		for (i = 0, i_f = i_min; i < length; i++, i_f++)
		{
			for (j = 0, j_f = j_min; j < width; j++, j_f++)
			{
				xd = x_new(i, j, length);
				size_t xd_f = x_new(i_f, j_f, Length);
				imageData[k][xd] = imageDataFull[k_f][xd_f];
				initialsegment[k][xd] = initialsegmentFull[k_f][xd_f];
			}
		}
	}

	dilatation3dHeighteenNeigbours(initialsegment, length, width, height, 1, 0);

	storing_path = root + "output/cropped_liver_p2.raw";
	manageRAWFile3D<dataType>(initialsegment, length, width, height, storing_path.c_str(), STORE_DATA, false);

	storing_path = root + "output/cropped_image_p2.raw";
	manageRAWFile3D<dataType>(imageData, length, width, height, storing_path.c_str(), STORE_DATA, false);

	Point3D croppedOrigin_p2 = { i_min, j_min, k_min };
	croppedOrigin_p2 = getRealCoordFromImageCoord3D(croppedOrigin_p2, imageOrigin, imageSpacing, orientation);
	std::cout << "Cropped image origin : (" << croppedOrigin_p2.x << ", " << croppedOrigin_p2.y << ", " << croppedOrigin_p2.z << ")" << std::endl;

	imageDataDistance = { height, length, width, initialsegment, croppedOrigin_p2, imageSpacing, orientation };
	fastMarchingDistanceMap(imageDataDistance, distanceMap, 0.0);
	storing_path = root + "output/distance_map_liver_p2.raw";
	manageRAWFile3D<dataType>(distanceMap, length, width, height, storing_path.c_str(), STORE_DATA, false);

	//Copy
	max_d = 0.0;
	for (k = 0; k < height; k++)
	{
		for (i = 0; i < dim2d; i++)
		{
			initialsegment[k][i] = distanceMap[k][i];
			if (max_d < distanceMap[k][i])
			{
				max_d = distanceMap[k][i];
			}
		}
	}
	rescaleNewRange(initialsegment, length, width, height, 0.0, max_d, 0.0, 1.0);

	// DownSampling
	half_length = length / 2;
	half_width = width / 2;
	half_height = height / 2;
	downSpacing = { (dataType)(ctContainer->spacing[0] * 2.0), (dataType)(ctContainer->spacing[1] * 2.0), (dataType)(ctContainer->spacing[2] * 2.0) };

	dataType** imageDataDown = new dataType * [half_height];
	dataType** initialsegmentDown = new dataType * [half_height];
	for (k = 0; k < half_height; k++)
	{
		imageDataDown[k] = new dataType[half_length * half_width]{ 0 };
		initialsegmentDown[k] = new dataType[half_length * half_width]{ 0 };
	}

	imageToSegment = { height, length, width, initialsegment, croppedOrigin_p2, imageSpacing, orientation };
	imageToSegmentDown = { half_height, half_length, half_width, initialsegmentDown, croppedOrigin_p2, downSpacing, orientation };

	imageInterpolation3D(imageToSegment, imageToSegmentDown, NEAREST_NEIGHBOR);

	imageToSegment.imageDataPtr = imageData;
	imageToSegmentDown.imageDataPtr = imageDataDown;

	imageInterpolation3D(imageToSegment, imageToSegmentDown, NEAREST_NEIGHBOR);

	storing_path = root + "output/Segmentation/Liver/P2/";
	GSUBSURF(imageToSegmentDown, initialsegmentDown, storing_path.c_str(), smoothParameters, segParameters);

	for (k = 0; k < half_height; k++)
	{
		delete[] imageDataDown[k];
		delete[] initialsegmentDown[k];
	}
	delete[] imageDataDown;
	delete[] initialsegmentDown;

	for (k = 0; k < height; k++)
	{
		delete[] imageData[k];
		delete[] initialsegment[k];
	}
	delete[] imageData;
	delete[] initialsegment;

	for (k = 0; k < Height; k++)
	{
		delete[] initialsegmentFull[k];
		delete[] imageDataFull[k];
	}
	delete[] imageDataFull;
	delete[] initialsegmentFull;

	//====================== Patient 3 ============================================

	std::cout << "Patient 3" << std::endl;
	loading_path = root + "input/vtk/petct/ct/Patient3_ct.vtk";
	readVtkFile(loading_path.c_str(), ctContainer);

	Height = (size_t)ctContainer->dimensions[2];
	Length = (size_t)ctContainer->dimensions[0];
	Width = (size_t)ctContainer->dimensions[1];
	dim2D = Length * Width;
	std::cout << "CT image dim : " << ctContainer->dimensions[0] << " x " << ctContainer->dimensions[1] << " x " << ctContainer->dimensions[2] << "" << std::endl;

	std::cout << "CT origin : (" << ctContainer->origin[0] << ", " << ctContainer->origin[1] << ", " << ctContainer->origin[2] << ")" << std::endl;
	imageOrigin = { (dataType)ctContainer->origin[0], (dataType)ctContainer->origin[1], (dataType)ctContainer->origin[2] };
	imageSpacing = { (dataType)ctContainer->spacing[0], (dataType)ctContainer->spacing[1], (dataType)ctContainer->spacing[2] };
	std::cout << "CT spacing : (" << ctContainer->spacing[0] << ", " << ctContainer->spacing[1] << ", " << ctContainer->spacing[2] << ")" << std::endl;

	dataType** imageDataFull = new dataType * [Height];
	dataType** initialsegmentFull = new dataType * [Height];
	for (k = 0; k < Height; k++)
	{
		imageDataFull[k] = new dataType[dim2D]{ 0 };
		initialsegmentFull[k] = new dataType[dim2D]{ 0 };
	}
	loading_path = root + "input/raw/filtered/filtered_p3.raw";
	manageRAWFile3D<dataType>(imageDataFull, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	loading_path = root + "input/raw/liver/liver_p3.raw";
	manageRAWFile3D<dataType>(initialsegmentFull, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	////Cropping

	//P3
	dataType k_min = 115, k_max = 205;
	dataType i_min = 100, i_max = 400;
	dataType j_min = 130, j_max = 430;

	length = i_max - i_min;
	width = j_max - j_min;
	height = k_max - k_min;
	dim2d = length * width;

	std::cout << "Cropped Height : " << height << ", Length : " << length << ", Width : " << width << std::endl;

	dataType** imageData = new dataType * [height];
	dataType** initialsegment = new dataType * [height];
	dataType** distanceMap = new dataType * [height];
	for (k = 0; k < height; k++)
	{
		imageData[k] = new dataType[dim2d]{ 0 };
		initialsegment[k] = new dataType[dim2d]{ 0 };
		distanceMap[k] = new dataType[dim2d]{ 0 };
	}

	k_f = 0, i_f = 0, j_f = 0;
	for (k = 0, k_f = k_min; k < height; k++, k_f++)
	{
		for (i = 0, i_f = i_min; i < length; i++, i_f++)
		{
			for (j = 0, j_f = j_min; j < width; j++, j_f++)
			{
				xd = x_new(i, j, length);
				size_t xd_f = x_new(i_f, j_f, Length);
				imageData[k][xd] = imageDataFull[k_f][xd_f];
				initialsegment[k][xd] = initialsegmentFull[k_f][xd_f];
			}
		}
	}

	dilatation3dHeighteenNeigbours(initialsegment, length, width, height, 1, 0);

	storing_path = root + "output/cropped_liver_p3.raw";
	manageRAWFile3D<dataType>(initialsegment, length, width, height, storing_path.c_str(), STORE_DATA, false);

	storing_path = root + "output/cropped_image_p3.raw";
	manageRAWFile3D<dataType>(imageData, length, width, height, storing_path.c_str(), STORE_DATA, false);

	Point3D croppedOrigin_p3 = { i_min, j_min, k_min };
	croppedOrigin_p3 = getRealCoordFromImageCoord3D(croppedOrigin_p3, imageOrigin, imageSpacing, orientation);
	std::cout << "Cropped image origin : (" << croppedOrigin_p3.x << ", " << croppedOrigin_p3.y << ", " << croppedOrigin_p3.z << ")" << std::endl;

	imageDataDistance = { height, length, width, initialsegment, croppedOrigin_p3, imageSpacing, orientation };
	fastMarchingDistanceMap(imageDataDistance, distanceMap, 0.0);
	storing_path = root + "output/distance_map_liver_p3.raw";
	manageRAWFile3D<dataType>(distanceMap, length, width, height, storing_path.c_str(), STORE_DATA, false);

	//Copy
	max_d = 0.0;
	for (k = 0; k < height; k++)
	{
		for (i = 0; i < dim2d; i++)
		{
			initialsegment[k][i] = distanceMap[k][i];
			if (max_d < distanceMap[k][i])
			{
				max_d = distanceMap[k][i];
			}
		}
	}
	rescaleNewRange(initialsegment, length, width, height, 0.0, max_d, 0.0, 1.0);

	// DownSampling
	half_length = length / 2;
	half_width = width / 2;
	half_height = height / 2;
	downSpacing = { (dataType)(ctContainer->spacing[0] * 2.0), (dataType)(ctContainer->spacing[1] * 2.0), (dataType)(ctContainer->spacing[2] * 2.0) };

	dataType** imageDataDown = new dataType * [half_height];
	dataType** initialsegmentDown = new dataType * [half_height];
	for (k = 0; k < half_height; k++)
	{
		imageDataDown[k] = new dataType[half_length * half_width]{ 0 };
		initialsegmentDown[k] = new dataType[half_length * half_width]{ 0 };
	}

	imageToSegment = { height, length, width, initialsegment, croppedOrigin_p3, imageSpacing, orientation };
	imageToSegmentDown = { half_height, half_length, half_width, initialsegmentDown, croppedOrigin_p3, downSpacing, orientation };

	imageInterpolation3D(imageToSegment, imageToSegmentDown, NEAREST_NEIGHBOR);

	imageToSegment.imageDataPtr = imageData;
	imageToSegmentDown.imageDataPtr = imageDataDown;

	imageInterpolation3D(imageToSegment, imageToSegmentDown, NEAREST_NEIGHBOR);

	storing_path = root + "output/Segmentation/Liver/P3/";
	GSUBSURF(imageToSegmentDown, initialsegmentDown, storing_path.c_str(), smoothParameters, segParameters);

	for (k = 0; k < half_height; k++)
	{
		delete[] imageDataDown[k];
		delete[] initialsegmentDown[k];
	}
	delete[] imageDataDown;
	delete[] initialsegmentDown;

	for (k = 0; k < height; k++)
	{
		delete[] imageData[k];
		delete[] initialsegment[k];
	}
	delete[] imageData;
	delete[] initialsegment;

	for (k = 0; k < Height; k++)
	{
		delete[] initialsegmentFull[k];
		delete[] imageDataFull[k];
	}
	delete[] imageDataFull;
	delete[] initialsegmentFull;

	//====================== Patient 4 ============================================

	std::cout << "Patient 4" << std::endl;
	loading_path = root + "input/vtk/petct/ct/Patient4_ct.vtk";
	readVtkFile(loading_path.c_str(), ctContainer);

	Height = (size_t)ctContainer->dimensions[2];
	Length = (size_t)ctContainer->dimensions[0];
	Width = (size_t)ctContainer->dimensions[1];
	dim2D = Length * Width;
	std::cout << "CT image dim : " << ctContainer->dimensions[0] << " x " << ctContainer->dimensions[1] << " x " << ctContainer->dimensions[2] << "" << std::endl;

	std::cout << "CT origin : (" << ctContainer->origin[0] << ", " << ctContainer->origin[1] << ", " << ctContainer->origin[2] << ")" << std::endl;
	imageOrigin = { (dataType)ctContainer->origin[0], (dataType)ctContainer->origin[1], (dataType)ctContainer->origin[2] };
	imageSpacing = { (dataType)ctContainer->spacing[0], (dataType)ctContainer->spacing[1], (dataType)ctContainer->spacing[2] };
	std::cout << "CT spacing : (" << ctContainer->spacing[0] << ", " << ctContainer->spacing[1] << ", " << ctContainer->spacing[2] << ")" << std::endl;

	dataType** imageDataFull = new dataType * [Height];
	dataType** initialsegmentFull = new dataType * [Height];
	for (k = 0; k < Height; k++)
	{
		imageDataFull[k] = new dataType[dim2D]{ 0 };
		initialsegmentFull[k] = new dataType[dim2D]{ 0 };
	}
	loading_path = root + "input/raw/filtered/filtered_p4.raw";
	manageRAWFile3D<dataType>(imageDataFull, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	loading_path = root + "input/raw/liver/liver_p4.raw";
	manageRAWFile3D<dataType>(initialsegmentFull, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	////Cropping

	//P4
	k_min = 120, k_max = 200;
	i_min = 120, i_max = 400;
	j_min = 130, j_max = 410;

	length = i_max - i_min;
	width = j_max - j_min;
	height = k_max - k_min;
	dim2d = length * width;

	std::cout << "Cropped Height : " << height << ", Length : " << length << ", Width : " << width << std::endl;

	dataType** imageData = new dataType * [height];
	dataType** initialsegment = new dataType * [height];
	dataType** distanceMap = new dataType * [height];
	for (k = 0; k < height; k++)
	{
		imageData[k] = new dataType[dim2d]{ 0 };
		initialsegment[k] = new dataType[dim2d]{ 0 };
		distanceMap[k] = new dataType[dim2d]{ 0 };
	}

	k_f = 0, i_f = 0, j_f = 0;
	for (k = 0, k_f = k_min; k < height; k++, k_f++)
	{
		for (i = 0, i_f = i_min; i < length; i++, i_f++)
		{
			for (j = 0, j_f = j_min; j < width; j++, j_f++)
			{
				xd = x_new(i, j, length);
				size_t xd_f = x_new(i_f, j_f, Length);
				imageData[k][xd] = imageDataFull[k_f][xd_f];
				initialsegment[k][xd] = initialsegmentFull[k_f][xd_f];
			}
		}
	}

	dilatation3dHeighteenNeigbours(initialsegment, length, width, height, 1, 0);

	storing_path = root + "output/cropped_liver_p4.raw";
	manageRAWFile3D<dataType>(initialsegment, length, width, height, storing_path.c_str(), STORE_DATA, false);

	storing_path = root + "output/cropped_image_p4.raw";
	manageRAWFile3D<dataType>(imageData, length, width, height, storing_path.c_str(), STORE_DATA, false);

	Point3D croppedOrigin_p4 = { i_min, j_min, k_min };
	croppedOrigin_p4 = getRealCoordFromImageCoord3D(croppedOrigin_p4, imageOrigin, imageSpacing, orientation);
	std::cout << "Cropped image origin : (" << croppedOrigin_p4.x << ", " << croppedOrigin_p4.y << ", " << croppedOrigin_p4.z << ")" << std::endl;

	imageDataDistance = { height, length, width, initialsegment, croppedOrigin_p4, imageSpacing, orientation };
	fastMarchingDistanceMap(imageDataDistance, distanceMap, 0.0);
	storing_path = root + "output/distance_map_liver_p4.raw";
	manageRAWFile3D<dataType>(distanceMap, length, width, height, storing_path.c_str(), STORE_DATA, false);

	//Copy
	max_d = 0.0;
	for (k = 0; k < height; k++)
	{
		for (i = 0; i < dim2d; i++)
		{
			initialsegment[k][i] = distanceMap[k][i];
			if (max_d < distanceMap[k][i])
			{
				max_d = distanceMap[k][i];
			}
		}
	}
	rescaleNewRange(initialsegment, length, width, height, 0.0, max_d, 0.0, 1.0);

	// DownSampling
	half_length = length / 2;
	half_width = width / 2;
	half_height = height / 2;
	downSpacing = { (dataType)(ctContainer->spacing[0] * 2.0), (dataType)(ctContainer->spacing[1] * 2.0), (dataType)(ctContainer->spacing[2] * 2.0) };

	dataType** imageDataDown = new dataType * [half_height];
	dataType** initialsegmentDown = new dataType * [half_height];
	for (k = 0; k < half_height; k++)
	{
		imageDataDown[k] = new dataType[half_length * half_width]{ 0 };
		initialsegmentDown[k] = new dataType[half_length * half_width]{ 0 };
	}

	imageToSegment = { height, length, width, initialsegment, croppedOrigin_p4, imageSpacing, orientation };
	imageToSegmentDown = { half_height, half_length, half_width, initialsegmentDown, croppedOrigin_p4, downSpacing, orientation };

	imageInterpolation3D(imageToSegment, imageToSegmentDown, NEAREST_NEIGHBOR);

	imageToSegment.imageDataPtr = imageData;
	imageToSegmentDown.imageDataPtr = imageDataDown;

	imageInterpolation3D(imageToSegment, imageToSegmentDown, NEAREST_NEIGHBOR);

	storing_path = root + "output/Segmentation/Liver/P4/";
	GSUBSURF(imageToSegmentDown, initialsegmentDown, storing_path.c_str(), smoothParameters, segParameters);

	for (k = 0; k < half_height; k++)
	{
		delete[] imageDataDown[k];
		delete[] initialsegmentDown[k];
	}
	delete[] imageDataDown;
	delete[] initialsegmentDown;

	for (k = 0; k < height; k++)
	{
		delete[] imageData[k];
		delete[] initialsegment[k];
	}
	delete[] imageData;
	delete[] initialsegment;

	for (k = 0; k < Height; k++)
	{
		delete[] initialsegmentFull[k];
		delete[] imageDataFull[k];
	}
	delete[] imageDataFull;
	delete[] initialsegmentFull;

	//====================== Patient 5 ============================================

	std::cout << "Patient 5" << std::endl;
	loading_path = root + "input/vtk/petct/ct/Patient5_ct.vtk";
	readVtkFile(loading_path.c_str(), ctContainer);

	Height = (size_t)ctContainer->dimensions[2];
	Length = (size_t)ctContainer->dimensions[0];
	Width = (size_t)ctContainer->dimensions[1];
	dim2D = Length * Width;
	std::cout << "CT image dim : " << ctContainer->dimensions[0] << " x " << ctContainer->dimensions[1] << " x " << ctContainer->dimensions[2] << "" << std::endl;

	std::cout << "CT origin : (" << ctContainer->origin[0] << ", " << ctContainer->origin[1] << ", " << ctContainer->origin[2] << ")" << std::endl;
	imageOrigin = { (dataType)ctContainer->origin[0], (dataType)ctContainer->origin[1], (dataType)ctContainer->origin[2] };
	imageSpacing = { (dataType)ctContainer->spacing[0], (dataType)ctContainer->spacing[1], (dataType)ctContainer->spacing[2] };
	std::cout << "CT spacing : (" << ctContainer->spacing[0] << ", " << ctContainer->spacing[1] << ", " << ctContainer->spacing[2] << ")" << std::endl;

	dataType** imageDataFull = new dataType * [Height];
	dataType** initialsegmentFull = new dataType * [Height];
	for (k = 0; k < Height; k++)
	{
		imageDataFull[k] = new dataType[dim2D]{ 0 };
		initialsegmentFull[k] = new dataType[dim2D]{ 0 };
	}
	loading_path = root + "input/raw/filtered/filtered_p5.raw";
	manageRAWFile3D<dataType>(imageDataFull, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	loading_path = root + "input/raw/liver/liver_p5.raw";
	manageRAWFile3D<dataType>(initialsegmentFull, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	////Cropping

	//P5
	k_min = 460, k_max = 590;
	i_min = 130, i_max = 380;
	j_min = 140, j_max = 390;

	length = i_max - i_min;
	width = j_max - j_min;
	height = k_max - k_min;
	dim2d = length * width;

	std::cout << "Cropped Height : " << height << ", Length : " << length << ", Width : " << width << std::endl;

	dataType** imageData = new dataType * [height];
	dataType** initialsegment = new dataType * [height];
	dataType** distanceMap = new dataType * [height];
	for (k = 0; k < height; k++)
	{
		imageData[k] = new dataType[dim2d]{ 0 };
		initialsegment[k] = new dataType[dim2d]{ 0 };
		distanceMap[k] = new dataType[dim2d]{ 0 };
	}

	k_f = 0, i_f = 0, j_f = 0;
	for (k = 0, k_f = k_min; k < height; k++, k_f++)
	{
		for (i = 0, i_f = i_min; i < length; i++, i_f++)
		{
			for (j = 0, j_f = j_min; j < width; j++, j_f++)
			{
				xd = x_new(i, j, length);
				size_t xd_f = x_new(i_f, j_f, Length);
				imageData[k][xd] = imageDataFull[k_f][xd_f];
				initialsegment[k][xd] = initialsegmentFull[k_f][xd_f];
			}
		}
	}

	dilatation3dHeighteenNeigbours(initialsegment, length, width, height, 1, 0);

	storing_path = root + "output/cropped_liver_p5.raw";
	manageRAWFile3D<dataType>(initialsegment, length, width, height, storing_path.c_str(), STORE_DATA, false);

	storing_path = root + "output/cropped_image_p5.raw";
	manageRAWFile3D<dataType>(imageData, length, width, height, storing_path.c_str(), STORE_DATA, false);

	Point3D croppedOrigin_p5 = { i_min, j_min, k_min };
	croppedOrigin_p5 = getRealCoordFromImageCoord3D(croppedOrigin_p5, imageOrigin, imageSpacing, orientation);
	std::cout << "Cropped image origin : (" << croppedOrigin_p5.x << ", " << croppedOrigin_p5.y << ", " << croppedOrigin_p5.z << ")" << std::endl;

	imageDataDistance = { height, length, width, initialsegment, croppedOrigin_p5, imageSpacing, orientation };
	fastMarchingDistanceMap(imageDataDistance, distanceMap, 0.0);
	storing_path = root + "output/distance_map_liver_p5.raw";
	manageRAWFile3D<dataType>(distanceMap, length, width, height, storing_path.c_str(), STORE_DATA, false);

	//Copy
	max_d = 0.0;
	for (k = 0; k < height; k++)
	{
		for (i = 0; i < dim2d; i++)
		{
			initialsegment[k][i] = distanceMap[k][i];
			if (max_d < distanceMap[k][i])
			{
				max_d = distanceMap[k][i];
			}
		}
	}
	rescaleNewRange(initialsegment, length, width, height, 0.0, max_d, 0.0, 1.0);

	// DownSampling
	half_length = length / 2;
	half_width = width / 2;
	half_height = height / 2;
	downSpacing = { (dataType)(ctContainer->spacing[0] * 2.0), (dataType)(ctContainer->spacing[1] * 2.0), (dataType)(ctContainer->spacing[2] * 2.0) };

	dataType** imageDataDown = new dataType * [half_height];
	dataType** initialsegmentDown = new dataType * [half_height];
	for (k = 0; k < half_height; k++)
	{
		imageDataDown[k] = new dataType[half_length * half_width]{ 0 };
		initialsegmentDown[k] = new dataType[half_length * half_width]{ 0 };
	}

	imageToSegment = { height, length, width, initialsegment, croppedOrigin_p5, imageSpacing, orientation };
	imageToSegmentDown = { half_height, half_length, half_width, initialsegmentDown, croppedOrigin_p5, downSpacing, orientation };

	imageInterpolation3D(imageToSegment, imageToSegmentDown, NEAREST_NEIGHBOR);

	imageToSegment.imageDataPtr = imageData;
	imageToSegmentDown.imageDataPtr = imageDataDown;

	imageInterpolation3D(imageToSegment, imageToSegmentDown, NEAREST_NEIGHBOR);

	storing_path = root + "output/Segmentation/Liver/P5/";
	GSUBSURF(imageToSegmentDown, initialsegmentDown, storing_path.c_str(), smoothParameters, segParameters);

	for (k = 0; k < half_height; k++)
	{
		delete[] imageDataDown[k];
		delete[] initialsegmentDown[k];
	}
	delete[] imageDataDown;
	delete[] initialsegmentDown;

	for (k = 0; k < height; k++)
	{
		delete[] imageData[k];
		delete[] initialsegment[k];
	}
	delete[] imageData;
	delete[] initialsegment;

	for (k = 0; k < Height; k++)
	{
		delete[] initialsegmentFull[k];
		delete[] imageDataFull[k];
	}
	delete[] imageDataFull;
	delete[] initialsegmentFull;

	//====================== Patient 6 ============================================

	std::cout << "Patient 6" << std::endl;
	loading_path = root + "input/vtk/petct/ct/Patient6_ct.vtk";
	readVtkFile(loading_path.c_str(), ctContainer);

	Height = (size_t)ctContainer->dimensions[2];
	Length = (size_t)ctContainer->dimensions[0];
	Width = (size_t)ctContainer->dimensions[1];
	dim2D = Length * Width;
	std::cout << "CT image dim : " << ctContainer->dimensions[0] << " x " << ctContainer->dimensions[1] << " x " << ctContainer->dimensions[2] << "" << std::endl;

	std::cout << "CT origin : (" << ctContainer->origin[0] << ", " << ctContainer->origin[1] << ", " << ctContainer->origin[2] << ")" << std::endl;
	imageOrigin = { (dataType)ctContainer->origin[0], (dataType)ctContainer->origin[1], (dataType)ctContainer->origin[2] };
	imageSpacing = { (dataType)ctContainer->spacing[0], (dataType)ctContainer->spacing[1], (dataType)ctContainer->spacing[2] };
	std::cout << "CT spacing : (" << ctContainer->spacing[0] << ", " << ctContainer->spacing[1] << ", " << ctContainer->spacing[2] << ")" << std::endl;

	dataType** imageDataFull = new dataType * [Height];
	dataType** initialsegmentFull = new dataType * [Height];
	for (k = 0; k < Height; k++)
	{
		imageDataFull[k] = new dataType[dim2D]{ 0 };
		initialsegmentFull[k] = new dataType[dim2D]{ 0 };
	}
	loading_path = root + "input/raw/filtered/filtered_p6.raw";
	manageRAWFile3D<dataType>(imageDataFull, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	loading_path = root + "input/raw/liver/liver_p6.raw";
	manageRAWFile3D<dataType>(initialsegmentFull, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	////Cropping

	//P6
	k_min = 280, k_max = 410;
	i_min = 90, i_max = 390;
	j_min = 200, j_max = 500;

	length = i_max - i_min;
	width = j_max - j_min;
	height = k_max - k_min;
	dim2d = length * width;

	std::cout << "Cropped Height : " << height << ", Length : " << length << ", Width : " << width << std::endl;

	dataType** imageData = new dataType * [height];
	dataType** initialsegment = new dataType * [height];
	dataType** distanceMap = new dataType * [height];
	for (k = 0; k < height; k++)
	{
		imageData[k] = new dataType[dim2d]{ 0 };
		initialsegment[k] = new dataType[dim2d]{ 0 };
		distanceMap[k] = new dataType[dim2d]{ 0 };
	}

	k_f = 0, i_f = 0, j_f = 0;
	for (k = 0, k_f = k_min; k < height; k++, k_f++)
	{
		for (i = 0, i_f = i_min; i < length; i++, i_f++)
		{
			for (j = 0, j_f = j_min; j < width; j++, j_f++)
			{
				xd = x_new(i, j, length);
				size_t xd_f = x_new(i_f, j_f, Length);
				imageData[k][xd] = imageDataFull[k_f][xd_f];
				initialsegment[k][xd] = initialsegmentFull[k_f][xd_f];
			}
		}
	}

	dilatation3dHeighteenNeigbours(initialsegment, length, width, height, 1, 0);

	storing_path = root + "output/cropped_liver_p6.raw";
	manageRAWFile3D<dataType>(initialsegment, length, width, height, storing_path.c_str(), STORE_DATA, false);

	storing_path = root + "output/cropped_image_p6.raw";
	manageRAWFile3D<dataType>(imageData, length, width, height, storing_path.c_str(), STORE_DATA, false);

	Point3D croppedOrigin_p6 = { i_min, j_min, k_min };
	croppedOrigin_p6 = getRealCoordFromImageCoord3D(croppedOrigin_p6, imageOrigin, imageSpacing, orientation);
	std::cout << "Cropped image origin : (" << croppedOrigin_p6.x << ", " << croppedOrigin_p6.y << ", " << croppedOrigin_p6.z << ")" << std::endl;

	imageDataDistance = { height, length, width, initialsegment, croppedOrigin_p6, imageSpacing, orientation };
	fastMarchingDistanceMap(imageDataDistance, distanceMap, 0.0);
	storing_path = root + "output/distance_map_liver_p6.raw";
	manageRAWFile3D<dataType>(distanceMap, length, width, height, storing_path.c_str(), STORE_DATA, false);

	//Copy
	dataType max_d = 0.0;
	for (k = 0; k < height; k++)
	{
		for (i = 0; i < dim2d; i++)
		{
			initialsegment[k][i] = distanceMap[k][i];
			if (max_d < distanceMap[k][i])
			{
				max_d = distanceMap[k][i];
			}
		}
	}
	rescaleNewRange(initialsegment, length, width, height, 0.0, max_d, 0.0, 1.0);

	// DownSampling
	half_length = length / 2;
	half_width = width / 2;
	half_height = height / 2;
	downSpacing = { (dataType)(ctContainer->spacing[0] * 2.0), (dataType)(ctContainer->spacing[1] * 2.0), (dataType)(ctContainer->spacing[2] * 2.0) };

	dataType** imageDataDown = new dataType * [half_height];
	dataType** initialsegmentDown = new dataType * [half_height];
	for (k = 0; k < half_height; k++)
	{
		imageDataDown[k] = new dataType[half_length * half_width]{ 0 };
		initialsegmentDown[k] = new dataType[half_length * half_width]{ 0 };
	}

	imageToSegment = { height, length, width, initialsegment, croppedOrigin_p6, imageSpacing, orientation };
	imageToSegmentDown = { half_height, half_length, half_width, initialsegmentDown, croppedOrigin_p6, downSpacing, orientation };

	imageInterpolation3D(imageToSegment, imageToSegmentDown, NEAREST_NEIGHBOR);

	imageToSegment.imageDataPtr = imageData;
	imageToSegmentDown.imageDataPtr = imageDataDown;

	imageInterpolation3D(imageToSegment, imageToSegmentDown, NEAREST_NEIGHBOR);

	storing_path = root + "output/Segmentation/Liver/P6/";
	GSUBSURF(imageToSegmentDown, initialsegmentDown, storing_path.c_str(), smoothParameters, segParameters);

	for (k = 0; k < half_height; k++)
	{
		delete[] imageDataDown[k];
		delete[] initialsegmentDown[k];
	}
	delete[] imageDataDown;
	delete[] initialsegmentDown;

	for (k = 0; k < height; k++)
	{
		delete[] imageData[k];
		delete[] initialsegment[k];
	}
	delete[] imageData;
	delete[] initialsegment;

	for (k = 0; k < Height; k++)
	{
		delete[] initialsegmentFull[k];
		delete[] imageDataFull[k];
	}
	delete[] imageDataFull;
	delete[] initialsegmentFull;

	//========================================================
	*/

	//==================== Liver centroid estimation =========================

	free(ctContainer);
	return EXIT_SUCCESS;
}