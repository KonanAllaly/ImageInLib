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
	
	std::cout<<"Original data"<<std::endl;
	loading_path = root + "input/vtk/petct/ct/Patient1_ct.vtk";
	readVtkFile(loading_path.c_str(), ctContainer);

	size_t Height = (size_t)ctContainer->dimensions[2];
	size_t Length = (size_t)ctContainer->dimensions[0];
	size_t Width = (size_t)ctContainer->dimensions[1];
	size_t dim2D = Length * Width;
	std::cout << "CT image dim : " << ctContainer->dimensions[0] << " x " << ctContainer->dimensions[1] << " x " << ctContainer->dimensions[2] << "" << std::endl;

	std::cout << "CT origin : (" << ctContainer->origin[0] << ", " << ctContainer->origin[1] << ", " << ctContainer->origin[2] << ")" << std::endl;
	Point3D imageOrigin = { (dataType)ctContainer->origin[0], (dataType)ctContainer->origin[1], (dataType)ctContainer->origin[2] };
	VoxelSpacing imageSpacing = { (dataType)ctContainer->spacing[0], (dataType)ctContainer->spacing[1], (dataType)ctContainer->spacing[2] };
	std::cout << "CT spacing : (" << ctContainer->spacing[0] << ", " << ctContainer->spacing[1] << ", " << ctContainer->spacing[2] << ")" << std::endl; 
	std::cout << "================================" << std::endl;

	/*
	//================== Liver  ============================================
	VoxelSpacing downSpacing;
	Image_Data imageToSegment, imageToSegmentDown, imageDataDistance;
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
	size_t length = 0, width = 0, height = 0, dim2d;
	size_t half_length = 0, half_width = 0, half_height = 0;
	size_t k_min = 0, k_max = 0, i_min = 0, i_max = 0, j_min = 0, j_max = 0;
	size_t k_f = 0, i_f = 0, j_f = 0, xd_f = 0;
	dataType max_d = 0;

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

	*/

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

	/*
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
	*/

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

	/*
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
	*/

	/*
	//=================== Loadt PET    ================================
	
	Vtk_File_Info* petContainer = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	petContainer->operation = copyFrom;
	std::cout << "Original data" << std::endl;
	loading_path = root + "input/vtk/petct/pet/Patient6_pet.vtk";
	readVtkFile(loading_path.c_str(), petContainer);

	Point3D petOrigin = { (dataType)petContainer->origin[0], (dataType)petContainer->origin[1], (dataType)petContainer->origin[2] };
	std::cout << "PET origin : (" << petContainer->origin[0] << ", " << petContainer->origin[1] << ", " << petContainer->origin[2] << ")" << std::endl;
	VoxelSpacing petSpacing = { (dataType)petContainer->spacing[0], (dataType)petContainer->spacing[1], (dataType)petContainer->spacing[2] };
	std::cout << "PET spacing : (" << petContainer->spacing[0] << ", " << petContainer->spacing[1] << ", " << petContainer->spacing[2] << ")" << std::endl;
	size_t petLength = (size_t)petContainer->dimensions[0];
	size_t petWidth = (size_t)petContainer->dimensions[1];
	size_t petHeight = (size_t)petContainer->dimensions[2];
	std::cout << "PET image dim : " << petContainer->dimensions[0] << " x " << petContainer->dimensions[1] << " x " << petContainer->dimensions[2] << "" << std::endl;
	
	//=================== SUV in Liver ===============================

	dataType** liver = new dataType * [Height];
	for (k = 0; k < Height; k++)
	{
		liver[k] = new dataType[dim2D]{ 0 };
	}
	loading_path = root + "input/raw/liver/liver_p6.raw";
	manageRAWFile3D<dataType>(liver, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	bool** status_liver = new bool* [petHeight];
	for (k = 0; k < petHeight; k++)
	{
		status_liver[k] = new bool[petLength * petWidth]{ false };
	}

	//Compute mean suv in liver
	dataType* centroid_liver = new dataType[3];
	centroidImage(liver, centroid_liver, Height, Length, Width, 0.0);
	Point3D centroid_liver_p = { centroid_liver[0], centroid_liver[1], centroid_liver[2] };
	centroid_liver_p = getRealCoordFromImageCoord3D(centroid_liver_p, imageOrigin, imageSpacing, orientation);

	double radius_liver_roi = 20 + 1.5 * fmax(petSpacing.sx, fmax(petSpacing.sy, petSpacing.sz));// approx 15mm

	double mean_suv_liver = 0.0;
	size_t count_suv_liver = 0;
	for (k = 0; k < petHeight; k++)
	{
		for (i = 0; i < petLength; i++)
		{
			for (j = 0; j < petWidth; j++)
			{
				xd = x_new(i, j, petLength);
				Point3D currentPoint = { (dataType)i, (dataType)j, (dataType)k };
				currentPoint = getRealCoordFromImageCoord3D(currentPoint, petOrigin, petSpacing, orientation);
				double suvValue = petContainer->dataPointer[k][x_new(i, j, petLength)];
				double d_current_centroid_liver = getPoint3DDistance(currentPoint, centroid_liver_p);
				if (d_current_centroid_liver <= radius_liver_roi && !status_liver[k][x_new(i, j, petLength)])
				{
					mean_suv_liver += suvValue;
					count_suv_liver++;
					status_liver[k][x_new(i, j, petLength)] = true;
				}
			}
		}
	}
	mean_suv_liver /= (double)count_suv_liver;
	std::cout << "Mean SUV in liver : " << mean_suv_liver << std::endl;

	for(k = 0; k < Height; k++)
	{
		delete[] liver[k];
		if(k < petHeight)
		{
			delete[] status_liver[k];
		}
	}
	delete[] liver;
	delete[] status_liver;

	//================== Load aorta segment and compute distance map =

	dataType** aorta = new dataType * [Height];
	dataType** ascending = new dataType * [Height];
	dataType** arch = new dataType * [Height];
	dataType** descending = new dataType * [Height];
	dataType** abdominal = new dataType * [Height];
	dataType** distanceMapAorta = new dataType * [Height];
	for (k = 0; k < Height; k++)
	{
		aorta[k] = new dataType[dim2D]{ 0 };
		ascending[k] = new dataType[dim2D]{ 0 };
		arch[k] = new dataType[dim2D]{ 0 };
		descending[k] = new dataType[dim2D]{ 0 };
		abdominal[k] = new dataType[dim2D]{ 0 };
		distanceMapAorta[k] = new dataType[dim2D]{ 0 };
	}
	loading_path = root + "output/Segmentation/Aorta/p6/segment_aorta_full_dim_p6.raw";
	manageRAWFile3D<dataType>(aorta, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	Image_Data aortaSegment = { Height, Length, Width, aorta, imageOrigin, imageSpacing, orientation };
	fastMarchingDistanceMap(aortaSegment, distanceMapAorta, 0.0);
	storing_path = root + "output/distance_map_p6.raw";
	manageRAWFile3D<dataType>(distanceMapAorta, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	//=================== SUV in segments ============================

	dataType x, y, z, index = 0;
	std::vector<Point3D> path_points_ascending, path_points_descending, path_points_arch, path_points_abdominal;
	dataType pet_suv = 0.0, max_suv_roi = 0.0;

	FILE* ratio_file;
	storing_path = root + "output/ratio_p6.csv";
	if (fopen_s(&ratio_file, storing_path.c_str(), "w") != 0)
	{
		printf("Unable to open");
		return false;
	}
	fprintf(ratio_file, "index,reference,ratio\n");

	string root_path = root + "output/Segmentation/Aorta/centered paths/";
	
	//Treat ascending aorta
	FILE* path_ascending;
	loading_path = root_path + "ascending/ascending_p6.csv";
	if (fopen_s(&path_ascending, loading_path.c_str(), "r") != 0)
	{
		printf("Unable to open");
		return false;
	}
	size_t k_max_ascending = 0;
	while (fscanf_s(path_ascending, "%f,%f,%f", &x, &y, &z) == 3)
	{
		Point3D current_point = { x, y, z };
		path_points_ascending.push_back(current_point);
		current_point = getImageCoordFromRealCoord3D(current_point, imageOrigin, imageSpacing, orientation);
		if ((size_t)current_point.z > k_max_ascending)
		{
			k_max_ascending = (size_t)current_point.z;
		}
	}
	fclose(path_ascending);

	dataType max_suv_ascending = 0.0;
	for (size_t n = 0; n < path_points_ascending.size(); n++) 
	{
		index++;
		Point3D current_point = path_points_ascending[n];
		Point3D petCoord = getImageCoordFromRealCoord3D(current_point, petOrigin, petSpacing, orientation);
		Point3D ctCoord = getImageCoordFromRealCoord3D(current_point, imageOrigin, imageSpacing, orientation);

		double radius = (double)distanceMapAorta[(size_t)ctCoord.z][x_new((size_t)ctCoord.x, (size_t)ctCoord.y, Length)];
		double radius_path_roi = radius + 2.0 * fmax(petSpacing.sx, fmax(petSpacing.sy, petSpacing.sz));
		
		BoundingBox3D box_roi = findBoundingBox3D(ctCoord, Length, Width, Height, radius_path_roi, 2.0);
		max_suv_roi = 0.0;
		for (k = box_roi.k_min; k < box_roi.k_max; k++)
		{
			for (i = box_roi.i_min; i < box_roi.i_max; i++)
			{
				for (j = box_roi.j_min; j < box_roi.j_max; j++)
				{
					xd = x_new(i, j, Length);
					Point3D currentPoint = { (dataType)i, (dataType)j, (dataType)k };
					currentPoint = getRealCoordFromImageCoord3D(currentPoint, imageOrigin, imageSpacing, orientation);
					Point3D currentPointpet = getImageCoordFromRealCoord3D(currentPoint, petOrigin, petSpacing, orientation);
					size_t k_pet = (size_t)currentPointpet.z;
					size_t i_pet = (size_t)currentPointpet.x;
					size_t j_pet = (size_t)currentPointpet.y;
					double d_current_path_point = getPoint3DDistance(currentPoint, current_point);
					pet_suv = petContainer->dataPointer[k_pet][x_new(i_pet, j_pet, petLength)];
					if (d_current_path_point <= radius_path_roi && aorta[k][xd] > 0.0 && k < k_max_ascending)
					{
						ascending[k][x_new(i, j, Length)] = 1.0;
						if(max_suv_ascending < pet_suv)
						{
							max_suv_ascending = pet_suv;
						}
						if(max_suv_roi < pet_suv)
						{
							max_suv_roi = pet_suv;
						}
					}
				}
			}
		}
		
		dataType reference_suv = 1.0;
		max_suv_roi /= mean_suv_liver;
		fprintf(ratio_file, "%d,%f,%f\n", index, reference_suv, max_suv_roi);
	}
	std::cout << "Max SUV in ascending aorta : " << max_suv_ascending / mean_suv_liver << std::endl;

	storing_path = root + "output/ascending_p6.raw";
	manageRAWFile3D<dataType>(ascending, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	//Treat aortic arch
	FILE* path_arch;
	loading_path = root_path + "arch/arch_p6.csv";
	if (fopen_s(&path_arch, loading_path.c_str(), "r") != 0)
	{
		printf("Unable to open");
		return false;
	}
	while (fscanf_s(path_arch, "%f,%f,%f", &x, &y, &z) == 3)
	{
		Point3D current_point = { x, y, z };
		path_points_arch.push_back(current_point);
	}
	fclose(path_arch);

	dataType max_suv_arch = 0.0;
	for (size_t n = 0; n < path_points_arch.size(); n++)
	{
		index++;
		Point3D current_point = path_points_arch[n];
		Point3D petCoord = getImageCoordFromRealCoord3D(current_point, petOrigin, petSpacing, orientation);
		Point3D ctCoord = getImageCoordFromRealCoord3D(current_point, imageOrigin, imageSpacing, orientation);

		double radius = (double)distanceMapAorta[(size_t)ctCoord.z][x_new((size_t)ctCoord.x, (size_t)ctCoord.y, Length)];
		double radius_path_roi = radius + 2.0 * fmax(petSpacing.sx, fmax(petSpacing.sy, petSpacing.sz));

		BoundingBox3D box_roi = findBoundingBox3D(ctCoord, Length, Width, Height, radius_path_roi, 2.0);
		max_suv_roi = 0.0;
		for (k = box_roi.k_min; k < box_roi.k_max; k++)
		{
			for (i = box_roi.i_min; i < box_roi.i_max; i++)
			{
				for (j = box_roi.j_min; j < box_roi.j_max; j++)
				{
					xd = x_new(i, j, Length);
					Point3D currentPoint = { (dataType)i, (dataType)j, (dataType)k };
					currentPoint = getRealCoordFromImageCoord3D(currentPoint, imageOrigin, imageSpacing, orientation);
					Point3D currentPointpet = getImageCoordFromRealCoord3D(currentPoint, petOrigin, petSpacing, orientation);
					size_t k_pet = (size_t)currentPointpet.z;
					size_t i_pet = (size_t)currentPointpet.x;
					size_t j_pet = (size_t)currentPointpet.y;
					double d_current_path_point = getPoint3DDistance(currentPoint, current_point);
					if ((d_current_path_point <= radius_path_roi) && (aorta[k][xd] > 0.0))
					{
						pet_suv = petContainer->dataPointer[k_pet][x_new(i_pet, j_pet, petLength)];
						arch[k][xd] = 1.0;
						if (max_suv_arch < pet_suv)
						{
							max_suv_arch = pet_suv;
						}
						if (max_suv_roi < pet_suv)
						{
							max_suv_roi = pet_suv;
						}
					}
				}
			}
		}

		dataType reference_suv = 1.0;
		max_suv_roi /= mean_suv_liver;
		fprintf(ratio_file, "%d,%f,%f\n", index, reference_suv, max_suv_roi);
	}
	std::cout << "Max SUV in aortic arch : " << max_suv_arch / mean_suv_liver << std::endl;

	storing_path = root + "output/arch_p6.raw";
	manageRAWFile3D<dataType>(arch, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	//Treat descending aorta
	FILE* path_descending;
	loading_path = root_path + "descending/descending_p6.csv";
	if (fopen_s(&path_descending, loading_path.c_str(), "r") != 0)
	{
		printf("Enable to open");
		return false;
	}
	size_t k_max_descending = 0;
	while (fscanf_s(path_descending, "%f,%f,%f", &x, &y, &z) == 3)
	{
		Point3D current_point = { x, y, z };
		path_points_descending.push_back(current_point);

		current_point = getImageCoordFromRealCoord3D(current_point, imageOrigin, imageSpacing, orientation);
		if ((size_t)current_point.z > k_max_descending)
		{
			k_max_descending = (size_t)current_point.z;
		}
	}
	fclose(path_descending);

	dataType max_suv_descending = 0.0;
	for (size_t n = 0; n < path_points_descending.size(); n++)
	{
		index++;
		Point3D current_point = path_points_descending[n];
		Point3D petCoord = getImageCoordFromRealCoord3D(current_point, petOrigin, petSpacing, orientation);
		Point3D ctCoord = getImageCoordFromRealCoord3D(current_point, imageOrigin, imageSpacing, orientation);

		double radius = (double)distanceMapAorta[(size_t)ctCoord.z][x_new((size_t)ctCoord.x, (size_t)ctCoord.y, Length)];
		double radius_path_roi = radius + 2.0 * fmax(petSpacing.sx, fmax(petSpacing.sy, petSpacing.sz));

		BoundingBox3D box_roi = findBoundingBox3D(ctCoord, Length, Width, Height, radius_path_roi, 2.0);
		max_suv_roi = 0.0;
		for (k = box_roi.k_min; k < box_roi.k_max; k++)
		{
			for (i = box_roi.i_min; i < box_roi.i_max; i++)
			{
				for (j = box_roi.j_min; j < box_roi.j_max; j++)
				{
					xd = x_new(i, j, Length);
					Point3D currentPoint = { (dataType)i, (dataType)j, (dataType)k };
					currentPoint = getRealCoordFromImageCoord3D(currentPoint, imageOrigin, imageSpacing, orientation);
					Point3D currentPointpet = getImageCoordFromRealCoord3D(currentPoint, petOrigin, petSpacing, orientation);
					size_t k_pet = (size_t)currentPointpet.z;
					size_t i_pet = (size_t)currentPointpet.x;
					size_t j_pet = (size_t)currentPointpet.y;
					double d_current_path_point = getPoint3DDistance(currentPoint, current_point);
					if (d_current_path_point <= radius_path_roi && aorta[k][x_new(i, j, Length)] > 0.0 && k < k_max_descending)
					{
						pet_suv = petContainer->dataPointer[k_pet][x_new(i_pet, j_pet, petLength)];
						descending[k][xd] = 1.0;
						if (max_suv_descending < pet_suv)
						{
							max_suv_descending = pet_suv;
						}
						if (max_suv_roi < pet_suv)
						{
							max_suv_roi = pet_suv;
						}
					}
				}
			}
		}

		dataType reference_suv = 1.0;
		max_suv_roi /= mean_suv_liver;
		fprintf(ratio_file, "%d,%f,%f\n", index, reference_suv, max_suv_roi);
	}
	std::cout << "Max SUV in descending aorta : " << max_suv_descending / mean_suv_liver << std::endl;

	storing_path = root + "output/descending_p6.raw";
	manageRAWFile3D<dataType>(descending, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	//Treat abdominal aorta
	FILE* path_abdominal;
	loading_path = root_path + "abdominal/abdomen_p6.csv";
	if (fopen_s(&path_abdominal, loading_path.c_str(), "r") != 0)
	{
		printf("Unable to open");
		return false;
	}
	size_t k_max_abdominal = 0;
	while (fscanf_s(path_abdominal, "%f,%f,%f", &x, &y, &z) == 3)
	{
		Point3D current_point = { x, y, z };
		path_points_abdominal.push_back(current_point);

		current_point = getImageCoordFromRealCoord3D(current_point, imageOrigin, imageSpacing, orientation);
		if ((size_t)current_point.z > k_max_abdominal)
		{
			k_max_abdominal = (size_t)current_point.z;
		}
	}
	fclose(path_abdominal);

	dataType max_suv_abdominal = 0.0;
	for (size_t n = 0; n < path_points_abdominal.size(); n++)
	{
		index++;
		Point3D current_point = path_points_abdominal[n];
		Point3D petCoord = getImageCoordFromRealCoord3D(current_point, petOrigin, petSpacing, orientation);
		Point3D ctCoord = getImageCoordFromRealCoord3D(current_point, imageOrigin, imageSpacing, orientation);

		double radius = (double)distanceMapAorta[(size_t)ctCoord.z][x_new((size_t)ctCoord.x, (size_t)ctCoord.y, Length)];
		double radius_path_roi = radius + 2.0 * fmax(petSpacing.sx, fmax(petSpacing.sy, petSpacing.sz));

		BoundingBox3D box_roi = findBoundingBox3D(ctCoord, Length, Width, Height, radius_path_roi, 2.0);
		max_suv_roi = 0.0;
		for (k = box_roi.k_min; k < box_roi.k_max; k++)
		{
			for (i = box_roi.i_min; i < box_roi.i_max; i++)
			{
				for (j = box_roi.j_min; j < box_roi.j_max; j++)
				{
					xd = x_new(i, j, Length);
					Point3D currentPoint = { (dataType)i, (dataType)j, (dataType)k };
					currentPoint = getRealCoordFromImageCoord3D(currentPoint, imageOrigin, imageSpacing, orientation);
					Point3D currentPointpet = getImageCoordFromRealCoord3D(currentPoint, petOrigin, petSpacing, orientation);
					size_t k_pet = (size_t)currentPointpet.z;
					size_t i_pet = (size_t)currentPointpet.x;
					size_t j_pet = (size_t)currentPointpet.y;
					double d_current_path_point = getPoint3DDistance(currentPoint, current_point);
					if (d_current_path_point <= radius_path_roi && aorta[k][x_new(i, j, Length)] > 0.0 && k < k_max_abdominal)
					{
						pet_suv = petContainer->dataPointer[k_pet][x_new(i_pet, j_pet, petLength)];
						abdominal[k][xd] = 1.0;
						if (max_suv_abdominal < pet_suv)
						{
							max_suv_abdominal = pet_suv;
						}
						if (max_suv_roi < pet_suv)
						{
							max_suv_roi = pet_suv;
						}
					}
				}
			}
		}

		dataType reference_suv = 1.0;
		max_suv_roi /= mean_suv_liver;
		fprintf(ratio_file, "%d,%f,%f\n", index, reference_suv, max_suv_roi);
	}
	std::cout << "Max SUV in abdominal aorta : " << max_suv_abdominal / mean_suv_liver << std::endl;

	storing_path = root + "output/abdominal_p6.raw";
	manageRAWFile3D<dataType>(abdominal, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);
	
	fclose(ratio_file);
	//============================

	path_points_ascending.clear();
	path_points_arch.clear();
	path_points_descending.clear();
	path_points_abdominal.clear();

	for(k = 0; k < Height; k++)
	{
		delete[] aorta[k];
		delete[] ascending[k];
		delete[] arch[k];
		delete[] descending[k];
		delete[] abdominal[k];
		delete[] distanceMapAorta[k];
	}
	delete[] aorta;
	delete[] distanceMapAorta;
	delete[] ascending;
	delete[] arch;
	delete[] descending;
	delete[] abdominal;
	*/

	/*
	//Load ascending aorta path
	FILE* path_arch;
	loading_path = root_path + "arch/arch_p1.csv";
	if (fopen_s(&path_arch, loading_path.c_str(), "r") != 0)
	{
		printf("Enable to open");
		return false;
	}
	size_t k_max_arch = 0;
	while (feof(path_arch) == 0)
	{
		fscanf_s(path_arch, "%f", &x);
		fscanf_s(path_arch, ",");
		fscanf_s(path_arch, "%f", &y);
		fscanf_s(path_arch, ",");
		fscanf_s(path_arch, "%f", &z);
		fscanf_s(path_arch, "\n");
		Point3D current_point = { x, y, z };
		path_points_arch.push_back(current_point);

		current_point = getImageCoordFromRealCoord3D(current_point, imageOrigin, imageSpacing, orientation);
		if ((size_t)current_point.z > k_max_arch)
		{
			k_max_arch = (size_t)current_point.z;
		}
	}
	fclose(path_arch);

	//Load descending aorta path
	FILE* path_descending;
	loading_path = root_path + "descending/descending_p1.csv";
	if (fopen_s(&path_descending, loading_path.c_str(), "r") != 0)
	{
		printf("Enable to open");
		return false;
	}
	size_t k_max_descending = 0;
	while (feof(path_descending) == 0)
	{
		fscanf_s(path_descending, "%f", &x);
		fscanf_s(path_descending, ",");
		fscanf_s(path_descending, "%f", &y);
		fscanf_s(path_descending, ",");
		fscanf_s(path_descending, "%f", &z);
		fscanf_s(path_descending, "\n");
		Point3D current_point = { x, y, z };
		path_points_descending.push_back(current_point);

		current_point = getImageCoordFromRealCoord3D(current_point, imageOrigin, imageSpacing, orientation);
		if ((size_t)current_point.z > k_max_descending)
		{
			k_max_descending = (size_t)current_point.z;
		}
	}
	fclose(path_descending);

	//Load abdominal aorta path
	FILE* path_abdominal;
	loading_path = root_path + "abdominal/abdominal_p1.csv";
	if (fopen_s(&path_abdominal, loading_path.c_str(), "r") != 0)
	{
		printf("Enable to open");
		return false;
	}
	size_t k_max_abdominal = 0;
	while (feof(path_abdominal) == 0)
	{
		fscanf_s(path_abdominal, "%f", &x);
		fscanf_s(path_abdominal, ",");
		fscanf_s(path_abdominal, "%f", &y);
		fscanf_s(path_abdominal, ",");
		fscanf_s(path_abdominal, "%f", &z);
		fscanf_s(path_abdominal, "\n");
		Point3D current_point = { x, y, z };
		path_points_abdominal.push_back(current_point);

		current_point = getImageCoordFromRealCoord3D(current_point, imageOrigin, imageSpacing, orientation);
		if ((size_t)current_point.z > k_max_abdominal)
		{
			k_max_abdominal = (size_t)current_point.z;
		}
	}
	fclose(path_abdominal);
	*/

	//free(petContainer);

	//for(k = 0; k < height; k++)
	//{
	//	delete[] imageData[k];
	//	delete[] initialsegment[k];
	//}
	//delete[] imageData;	
	//delete[] initialsegment;

	//for (k = 0; k < Height; k++)
	//{
	//	delete[] initialsegmentFull[k];
	//	delete[] imageDataFull[k];
	//}
	//delete[] imageDataFull;
	//delete[] initialsegmentFull;

	//=========== Potential function for path extraction =============

	dataType** potential = new dataType * [Height];
	dataType** aorta = new dataType * [Height];
	for (k = 0; k < Height; k++)
	{
		potential[k] = new dataType[dim2D]{ 0 };
		aorta[k] = new dataType[dim2D]{ 0 };
	}
	string root_potential = root + "output/Data journal paper submission/";
	//loading_path = root_potential + "P1/potential_p1.raw";
	//manageRAWFile3D<dataType>(potential, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	//storing_path = root + "output/slice_thorax_edge_image.raw";
	//storing_path = root + "output/slice_abdomen_edge_image.raw";
	//storing_path = root + "output/slice_heart_edge_image.raw";
	//storing_path = root + "output/slice_thorax_potential.raw";
	//storing_path = root + "output/slice_abdomen_potential.raw";
	//storing_path = root + "output/slice_heart_potential.raw";
	//manageRAWFile2D<dataType>(potential[229], Length, Width, storing_path.c_str(), STORE_DATA, false);

	string root_aorta = root + "output/Segmentation/Aorta/";
	//loading_path = root_aorta + "p1/segment_aorta_full_dim_p1.raw";
	//manageRAWFile3D<dataType>(aorta, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	//storing_path = root + "output/slice_thorax_aorta_segment.raw";
	//storing_path = root + "output/slice_abdomen_aorta_segment.raw";
	//storing_path = root + "output/slice_heart_aorta_segment.raw";
	//manageRAWFile2D<dataType>(aorta[255], Length, Width, storing_path.c_str(), STORE_DATA, false);

	//storing_path = root + "output/slice_thorax_ct_image.raw";
	//storing_path = root + "output/slice_abdomen_ct_image.raw";
	//storing_path = root + "output/slice_heart_ct_image.raw";
	//manageRAWFile2D<dataType>(ctContainer->dataPointer[255], Length, Width, storing_path.c_str(), STORE_DATA, false);

	//dataType z_view = -450.0;
	//for(k = 0; k < Height; k++)
	//{
	//	Point3D currentPoint = { 0.0, 0.0, (dataType)k };
	//	currentPoint = getRealCoordFromImageCoord3D(currentPoint, imageOrigin, imageSpacing, orientation);
	//	dataType z_min = currentPoint.z - 0.5 * imageSpacing.sz;
	//	dataType z_max = currentPoint.z + 0.5 * imageSpacing.sz;
	//	if (z_min <= z_view && z_max >= z_view)
	//	{
	//		std::cout << "The thorax slice is: " << k << std::endl;
	//	}
	//}

	/*
	//restore original image dimension for segment image
	size_t length = 150, width = 150, height = 350;
	dataType** imageData = new dataType * [height];
	for(k = 0; k < height; k++)
	{
		imageData[k] = new dataType[length * width]{ 0 };
	}
	loading_path = root_aorta + "p1/_seg_func_05000.raw";
	manageRAWFile3D<dataType>(imageData, length, width, height, loading_path.c_str(), LOAD_DATA, false);
	Point3D segmentOrigin = { -65.625, 4.375, -700.234 };
	VoxelSpacing segmentSpacing = { 1.171875, 1.171875, 1.171875 };
	Image_Data segmentImageStr = { height, length, width, imageData, segmentOrigin, segmentSpacing, orientation };

	for (k = 0; k < height; k++) 
	{
		for(i = 0; i < length; i++)
		{
			for(j = 0; j < width; j++)
			{
				xd = x_new(i, j, length);
				Point3D currentPoint = { (dataType)i, (dataType)j, (dataType)k };
				currentPoint = getRealCoordFromImageCoord3D(currentPoint, segmentOrigin, segmentSpacing, orientation);
				Point3D currentPointCt = getImageCoordFromRealCoord3D(currentPoint, imageOrigin, imageSpacing, orientation);
				size_t x = (size_t)currentPointCt.x, y = (size_t)currentPointCt.y, z = (size_t)currentPointCt.z;
				aorta[z][x_new(x, y, Length)] = imageData[k][xd];
			}
		}
	}

	storing_path = root_potential + "P1/level_sets_full_dim_p1.raw";
	manageRAWFile3D<dataType>(aorta, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	for(k = 0; k < height; k++)
	{
		delete[] imageData[k];
	}
	delete[] imageData;
	*/
	
	for (k = 0; k < Height; k++)
	{
		delete[] potential[k];
		delete[] aorta[k];
	}
	delete[] potential;
	delete[] aorta;

	free(ctContainer);
	return EXIT_SUCCESS;
}