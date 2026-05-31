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

	
	//================== Liver  ============================================
	
	std::cout<<"Original data"<<std::endl;
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
	std::cout << "================================" << std::endl;

	dataType** imageDataFull = new dataType * [Height];
	dataType** initialsegmentFull = new dataType * [Height];
	for (k = 0; k < Height; k++)
	{
		imageDataFull[k] = new dataType[dim2D]{ 0 };
		initialsegmentFull[k] = new dataType[dim2D]{ 0 };
	}
	//loading_path = root + "input/raw/filtered/filtered_p6.raw";
	//manageRAWFile3D<dataType>(imageDataFull, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	loading_path = root + "input/raw/liver/liver_p3.raw";
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

	//P3
	k_min = 115, k_max = 205;
	i_min = 100, i_max = 400;
	j_min = 130, j_max = 430;

	////P4
	//k_min = 120, k_max = 200;
	//i_min = 120, i_max = 400;
	//j_min = 130, j_max = 410;

	////P5
	//k_min = 455, k_max = 595;
	//i_min = 130, i_max = 380;
	//j_min = 140, j_max = 390;

	////P6
	//k_min = 280, k_max = 410;
	//i_min = 90, i_max = 390;
	//j_min = 200, j_max = 500;

	length = i_max - i_min;
	width = j_max - j_min;
	height = k_max - k_min;
	dim2d = length * width;

	std::cout<< "Cropped data" << std::endl;
	std::cout<< "Cropped Height : " << height << ", Length : " << length << ", Width : " << width << std::endl;
	std::cout << "================================" << std::endl;

	dataType** imageData = new dataType * [height];
	dataType** initialsegment = new dataType * [height];
	dataType** distanceMap = new dataType * [height];
	for (k = 0; k < height; k++)
	{
		imageData[k] = new dataType[dim2d]{ 0 };
		initialsegment[k] = new dataType[dim2d]{ 0 };
		distanceMap[k] = new dataType[dim2d]{ 0 };
	}

	/*
	k_f = 0, i_f = 0, j_f = 0;
	for(k = 0, k_f = k_min; k < height; k++, k_f++)
	{
		for (i = 0, i_f = i_min; i < length; i++, i_f++)
		{
			for(j = 0, j_f = j_min; j < width; j++, j_f++)
			{
				xd = x_new(i, j, length);
				size_t xd_f = x_new(i_f, j_f, Length);
				//imageData[k][xd] = imageDataFull[k_f][xd_f];
				initialsegment[k][xd] = initialsegmentFull[k_f][xd_f];
			}
		}
	}

	//dilatation3dHeighteenNeigbours(initialsegment, length, width, height, 1, 0);
	//dilatation3dHeighteenNeigbours(initialsegment, length, width, height, 1, 0);
	//dilatation3dHeighteenNeigbours(initialsegment, length, width, height, 1, 0);

	storing_path = root + "output/rough_segment_p5.raw";
	manageRAWFile3D<dataType>(initialsegment, length, width, height, storing_path.c_str(), STORE_DATA, false);

	//storing_path = root + "output/cropped_image_p6.raw";
	//manageRAWFile3D<dataType>(imageData, length, width, height, storing_path.c_str(), STORE_DATA, false);

	Point3D croppedOrigin = { i_min, j_min, k_min };
	croppedOrigin = getRealCoordFromImageCoord3D(croppedOrigin, imageOrigin, imageSpacing, orientation);
	std::cout << "Cropped image origin : (" << croppedOrigin.x << ", " << croppedOrigin.y << ", " << croppedOrigin.z << ")" << std::endl;

	//imageDataDistance = { height, length, width, initialsegment, croppedOrigin, imageSpacing, orientation };
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


	std::cout<<"Down sampled data" << std::endl;
	// DownSampling
	half_length = length / 2;
	half_width = width / 2;
	half_height = height / 2;
	downSpacing = { (dataType)(ctContainer->spacing[0] * 2.0), (dataType)(ctContainer->spacing[1] * 2.0), (dataType)(ctContainer->spacing[2] * 2.0) };
	std::cout << "Down sampled spacing : (" << downSpacing.sx << ", " << downSpacing.sy << ", " << downSpacing.sz << ")" << std::endl;
	std::cout << "Down sampled Height : " << half_height << ", Length : " << half_length << ", Width : " << half_width << std::endl;

	dataType** imageDataDown = new dataType * [half_height];
	dataType** initialsegmentDown = new dataType * [half_height];
	dataType** distanceMapDown = new dataType * [height];
	for (k = 0; k < half_height; k++)
	{
		imageDataDown[k] = new dataType[half_length * half_width]{ 0 };
		initialsegmentDown[k] = new dataType[half_length * half_width]{ 0 };
		distanceMapDown[k] = new dataType[half_length * half_width]{ 0 };
	}

	//imageToSegment = { height, length, width, initialsegment, croppedOrigin, imageSpacing, orientation };
	//imageToSegmentDown = { half_height, half_length, half_width, initialsegmentDown, croppedOrigin, downSpacing, orientation };
	//imageInterpolation3D(imageToSegment, imageToSegmentDown, NEAREST_NEIGHBOR);

	//imageToSegment.imageDataPtr = imageData;
	//imageToSegmentDown.imageDataPtr = imageDataDown;
	//imageInterpolation3D(imageToSegment, imageToSegmentDown, NEAREST_NEIGHBOR);

	//storing_path = root + "output/Segmentation/Liver/P6/";
	//GSUBSURF(imageToSegmentDown, initialsegmentDown, storing_path.c_str(), smoothParameters, segParameters);

	storing_path = root + "output/Segmentation/Liver/P5/_seg_func_2000.raw";
	manageRAWFile3D<dataType>(initialsegmentDown, half_length, half_width, half_height, storing_path.c_str(), LOAD_DATA, false);

	for(k = 0; k < half_height; k++)
	{
		for (i = 0; i < half_length * half_width; i++)
		{
			if(initialsegmentDown[k][i] > 0.015)
			{
				initialsegmentDown[k][i] = 1.0;
			}
			else
			{
				initialsegmentDown[k][i] = 0.0;
			}
		}
	}

	storing_path = root + "output/refined_segment_p5.raw";
	manageRAWFile3D<dataType>(initialsegmentDown, half_length, half_width, half_height, storing_path.c_str(), STORE_DATA, false);

	dataType* centroid_rough = new dataType[3];
	dataType* centroid_refined = new dataType[3];

	centroidImage(initialsegment, centroid_rough, height, length, width, 0.0);
	//std::cout << "Rough Centroid in image coordinate : (" << centroid_rough[0] << ", " << centroid_rough[1] << ", " << centroid_rough[2] << ")" << std::endl;
	centroidImage(initialsegmentDown, centroid_refined, half_height, half_length, half_width, 0.0);
	//std::cout << "Refined Centroid in image coordinate : (" << centroid_refined[0] << ", " << centroid_refined[1] << ", " << centroid_refined[2] << ")" << std::endl;

	Point3D centroid_rough_p = { centroid_rough[0], centroid_rough[1], centroid_rough[2] };
	Point3D centroid_refined_p = { centroid_refined[0], centroid_refined[1], centroid_refined[2] };

	centroid_rough_p = getRealCoordFromImageCoord3D(centroid_rough_p, croppedOrigin, imageSpacing, orientation);
	centroid_refined_p = getRealCoordFromImageCoord3D(centroid_refined_p, croppedOrigin, downSpacing, orientation);

	double distance_centroid = getPoint3DDistance(centroid_rough_p, centroid_refined_p);
	std::cout << "Distance between rough and refined centroid : " << distance_centroid << std::endl;

	Image_Data distanceSegmentRough = { height, length, width, initialsegment, croppedOrigin, imageSpacing, orientation };
	Image_Data distanceSegmentRefined = { half_height, half_length, half_width, initialsegmentDown, croppedOrigin, downSpacing, orientation };

	fastMarchingDistanceMap(distanceSegmentRough, distanceMap, 0.0);
	storing_path = root + "output/dmap_rough_segment_p5.raw";
	manageRAWFile3D<dataType>(distanceMap, length, width, height, storing_path.c_str(), STORE_DATA, false);
	fastMarchingDistanceMap(distanceSegmentRefined, distanceMapDown, 0.0);
	storing_path = root + "output/dmap_refined_segment_p5.raw";
	manageRAWFile3D<dataType>(distanceMapDown, half_length, half_width, half_height, storing_path.c_str(), STORE_DATA, false);

	Point3D roughCentroidInDistanceMap = getPointWithTheHighestValue(distanceMap, length, width, height);
	//std::cout << "Rough Centroid in distance map coordinate : (" << roughCentroidInDistanceMap.x << ", " << roughCentroidInDistanceMap.y << ", " << roughCentroidInDistanceMap.z << ")" << std::endl;
	Point3D refinedCentroidInDistanceMap = getPointWithTheHighestValue(distanceMapDown, half_length, half_width, half_height);
	//std::cout << "Refined Centroid in distance map coordinate : (" << refinedCentroidInDistanceMap.x << ", " << refinedCentroidInDistanceMap.y << ", " << refinedCentroidInDistanceMap.z << ")" << std::endl;

	roughCentroidInDistanceMap = getRealCoordFromImageCoord3D(roughCentroidInDistanceMap, croppedOrigin, imageSpacing, orientation);
	refinedCentroidInDistanceMap = getRealCoordFromImageCoord3D(refinedCentroidInDistanceMap, croppedOrigin, downSpacing, orientation);

	double distance_centroid_distanceMap = getPoint3DDistance(roughCentroidInDistanceMap, refinedCentroidInDistanceMap);
	std::cout << "Distance between rough and refined centroid in distance map : " << distance_centroid_distanceMap << std::endl;

	double d_rough_centroid_distanceMap = getPoint3DDistance(centroid_rough_p, roughCentroidInDistanceMap);
	std::cout << "Distance between rough centroid and distance map centroid : " << d_rough_centroid_distanceMap << std::endl;

	double d_refined_centroid_distanceMap = getPoint3DDistance(centroid_refined_p, refinedCentroidInDistanceMap);
	std::cout << "Distance between refined centroid and distance map centroid : " << d_refined_centroid_distanceMap << std::endl;

	
	//Save ball initial
	for(k = 0; k < height; k++)
	{
		for (i = 0; i < length; i++) 
		{
			for (j = 0; j < width; j++) 
			{
				xd = x_new(i, j, length);
				initialsegment[k][xd] = 0.0;
				Point3D currentPoint = { (dataType)i, (dataType)j, (dataType)k };
				currentPoint = getRealCoordFromImageCoord3D(currentPoint, croppedOrigin, imageSpacing, orientation);
				double distanceToRoughCentroid = getPoint3DDistance(currentPoint, centroid_rough_p);
				if (distanceToRoughCentroid <= 10.0)
				{
					initialsegment[k][xd] = 1.0;
				}

				distanceMap[k][xd] = 0.0;
				double distanceToRoughCentroidInDistanceMap = getPoint3DDistance(currentPoint, roughCentroidInDistanceMap);
				if (distanceToRoughCentroidInDistanceMap <= 10.0)
				{
					distanceMap[k][xd] = 1.0;
				}
			}
		}
	}
	storing_path = root + "output/ball_initial_rough_centroid_p5.raw";
	manageRAWFile3D<dataType>(initialsegment, length, width, height, storing_path.c_str(), STORE_DATA, false);
	storing_path = root + "output/ball_initial_rough_distance_map_p5.raw";
	manageRAWFile3D<dataType>(distanceMap, length, width, height, storing_path.c_str(), STORE_DATA, false);

	//Ball down
	for (k = 0; k < half_height; k++)
	{
		for (i = 0; i < half_length; i++)
		{
			for (j = 0; j < half_width; j++)
			{
				xd = x_new(i, j, half_length);
				initialsegmentDown[k][xd] = 0.0;
				Point3D currentPoint = { (dataType)i, (dataType)j, (dataType)k };
				currentPoint = getRealCoordFromImageCoord3D(currentPoint, croppedOrigin, downSpacing, orientation);
				double distanceToRefinedCentroid = getPoint3DDistance(currentPoint, centroid_refined_p);
				if (distanceToRefinedCentroid <= 10.0)
				{
					initialsegmentDown[k][xd] = 1.0;
				}
				distanceMapDown[k][xd] = 0.0;
				double distanceToRefinedCentroidInDistanceMap = getPoint3DDistance(currentPoint, refinedCentroidInDistanceMap);
				if (distanceToRefinedCentroidInDistanceMap <= 10.0)
				{
					distanceMapDown[k][xd] = 1.0;
				}
			}
		}
	}
	storing_path = root + "output/ball_down_refined_centroid_p5.raw";
	manageRAWFile3D<dataType>(initialsegmentDown, half_length, half_width, half_height, storing_path.c_str(), STORE_DATA, false);
	storing_path = root + "output/ball_down_refined_distance_map_p5.raw";
	manageRAWFile3D<dataType>(distanceMapDown, half_length, half_width, half_height, storing_path.c_str(), STORE_DATA, false);
	*/

	/*
	//======================== PET data ============================================
	Vtk_File_Info* petContainer = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	petContainer->operation = copyFrom;
	loading_path = root + "input/vtk/petct/pet/Patient5_pet.vtk";
	readVtkFile(loading_path.c_str(), petContainer);
	std::cout << "================================" << std::endl;
	std::cout << "PET data" << std::endl;
	std::cout << "PET image dim : " << petContainer->dimensions[0] << " x " << petContainer->dimensions[1] << " x " << petContainer->dimensions[2] << "" << std::endl;
	std::cout << "PET origin : (" << petContainer->origin[0] << ", " << petContainer->origin[1] << ", " << petContainer->origin[2] << ")" << std::endl;
	std::cout << "PET spacing : (" << petContainer->spacing[0] << ", " << petContainer->spacing[1] << ", " << petContainer->spacing[2] << ")" << std::endl;

	size_t petLength = (size_t)petContainer->dimensions[0];
	Point3D petOrigin = { (dataType)petContainer->origin[0], (dataType)petContainer->origin[1], (dataType)petContainer->origin[2] };
	VoxelSpacing petSpacing = { (dataType)petContainer->spacing[0], (dataType)petContainer->spacing[1], (dataType)petContainer->spacing[2] };

	//Compute mean suv in liver
	double radius_liver_roi = 20 + 1.5 * fmax(petSpacing.sx, fmax(petSpacing.sy, petSpacing.sz));// approx 15mm

	double suv_centroid_rough = 0.0;
	size_t count_centroid_rough = 0;

	double suv_centroid_refined = 0.0;
	size_t count_centroid_refined = 0;

	double suv_rough_distanceMap = 0.0;
	size_t count_rough_distanceMap = 0;

	double suv_refined_distanceMap = 0.0;
	size_t count_refined_distanceMap = 0;

	double suv_rough_global = 0.0;
	size_t count_rough_global = 0;

	double suv_refined_global = 0.0;
	size_t count_refined_global = 0;

	for (k = 0; k < height; k++) 
	{
		for(i = 0; i < length; i++)
		{
			for (j = 0; j < width; j++) 
			{
				xd = x_new(i, j, length);
				Point3D currentPoint = { (dataType)i, (dataType)j, (dataType)k };
				currentPoint = getRealCoordFromImageCoord3D(currentPoint, croppedOrigin, imageSpacing, orientation);
				
				Point3D petCoord = getImageCoordFromRealCoord3D(currentPoint, petOrigin, petSpacing, orientation);
				double suvValue = petContainer->dataPointer[(size_t)petCoord.z][x_new((size_t)petCoord.x, (size_t)petCoord.y, petLength)];

				double d_current_centroid_rough = getPoint3DDistance(currentPoint, centroid_rough_p);
				if (d_current_centroid_rough <= radius_liver_roi)
				{
					suv_centroid_rough += suvValue;
					count_centroid_rough++;
				}

				double d_current_centroid_refined = getPoint3DDistance(currentPoint, centroid_refined_p);
				if (d_current_centroid_refined <= radius_liver_roi)
				{
					suv_centroid_refined += suvValue;
					count_centroid_refined++;
				}

				double d_current_rough_distanceMap = getPoint3DDistance(currentPoint, roughCentroidInDistanceMap);
				if (d_current_rough_distanceMap <= radius_liver_roi)
				{
					suv_rough_distanceMap += suvValue;
					count_rough_distanceMap++;
				}

				double d_current_refined_distanceMap = getPoint3DDistance(currentPoint, refinedCentroidInDistanceMap);
				if (d_current_refined_distanceMap <= radius_liver_roi)
				{
					suv_refined_distanceMap += suvValue;
					count_refined_distanceMap++;
				}

				if(initialsegment[k][xd] > 0.0)
				{
					suv_rough_global += suvValue;
					count_rough_global++;
				}
			}
		}
	}

	suv_centroid_rough /= (double)count_centroid_rough;
	std::cout << "Mean SUV in rough centroid ROI : " << suv_centroid_rough << std::endl;
	suv_centroid_refined /= (double)count_centroid_refined;
	std::cout << "Mean SUV in refined centroid ROI : " << suv_centroid_refined << std::endl;
	suv_rough_distanceMap /= (double)count_rough_distanceMap;
	std::cout << "Mean SUV in rough distance map centroid ROI : " << suv_rough_distanceMap << std::endl;
	suv_refined_distanceMap /= (double)count_refined_distanceMap;
	std::cout << "Mean SUV in refined distance map centroid ROI : " << suv_refined_distanceMap << std::endl;
	suv_rough_global /= (double)count_rough_global;
	std::cout << "Mean SUV in rough segmentation : " << suv_rough_global << std::endl;
	
	free(petContainer);
	*/

	//==================== Estimate the diaphragm =========================
	
	dataType** lungs = new dataType * [Height];
	for (k = 0; k < Height; k++)
	{
		lungs[k] = new dataType[dim2D]{ 0 };
	}
	loading_path = root + "input/raw/lungs/-500HU/lungs_p3.raw";
	manageRAWFile3D<dataType>(lungs, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	size_t z_torax_abdomen = Height;
	for(k = 0; k < Height; k++)
	{
		for(i = 0; i < dim2D; i++)
		{
			if(lungs[k][i] > 0.0)
			{
				if(k < z_torax_abdomen)
				{
					z_torax_abdomen = k;
				}
			}
		}
	}
	std::cout << "Estimated diaphragm z coordinate in image coordinate : " << z_torax_abdomen << std::endl;

	Point3D origin_torax_abdomen = { 0.0, 0.0, (dataType)z_torax_abdomen };
	origin_torax_abdomen = getRealCoordFromImageCoord3D(origin_torax_abdomen, imageOrigin, imageSpacing, orientation);
	std::cout << "Estimated diaphragm z coordinate in real coordinate : (" << origin_torax_abdomen.x << ", " << origin_torax_abdomen.y << ", " << origin_torax_abdomen.z << ")" << std::endl;

	for(k = 0; k < Height; k++)
	{
		delete[] lungs[k];
	}
	delete[] lungs;

	//delete[] centroid_rough;
	//delete[] centroid_refined;

	//for(k = 0; k < half_height; k++)
	//{
	//	delete[] imageDataDown[k];
	//	delete[] initialsegmentDown[k];
	//	delete[] distanceMapDown[k];
	//}
	//delete[] imageDataDown;
	//delete[] initialsegmentDown;
	//delete[] distanceMapDown;

	//======================= Split abdomen and torax ============================================

	string root_path = root + "output/Segmentation/Aorta/centered paths/";
	// Input centered path
	FILE* path_file;
	loading_path = root_path + "finalCurve_p3.csv";
	if (fopen_s(&path_file, loading_path.c_str(), "r") != 0)
	{
		printf("Enable to open");
		return false;
	}

	//File to save thorax centerline points 
	string saving_thorax_csv = root_path + "path in thorax/thorax_p3.csv";
	FILE* f_thorax;
	if (fopen_s(&f_thorax, saving_thorax_csv.c_str(), "w") != 0)
	{
		printf("Enable to open");
		return false;
	}

	//File to save abdminal centerline points 
	string saving_abdomen_csv = root_path + "path in abdomen/abdomen_p3.csv";
	FILE* f_abdomen;
	if (fopen_s(&f_abdomen, saving_abdomen_csv.c_str(), "w") != 0)
	{
		printf("Enable to open");
		return false;
	}

	dataType x, y, z;
	while (feof(path_file) == 0)
	{
		fscanf_s(path_file, "%f", &x);
		fscanf_s(path_file, ",");
		fscanf_s(path_file, "%f", &y);
		fscanf_s(path_file, ",");
		fscanf_s(path_file, "%f", &z);
		fscanf_s(path_file, "\n");
		Point3D current_point = { x, y, z };

		Point3D current_point_ct = getImageCoordFromRealCoord3D(current_point, imageOrigin, imageSpacing, orientation);
		size_t z_current = (size_t)current_point_ct.z;
		if(z_current > z_torax_abdomen)
		{
			fprintf(f_thorax, "%f,%f,%f\n", x, y, z);
		}
		else
		{
			fprintf(f_abdomen, "%f,%f,%f\n", x, y, z);
		}
	}

	fclose(path_file);
	fclose(f_thorax);
	fclose(f_abdomen);

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

	free(ctContainer);
	return EXIT_SUCCESS;
}