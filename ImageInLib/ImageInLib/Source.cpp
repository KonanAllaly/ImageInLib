#include <iostream>
#include <sstream>  
#include <vector>
#include <string.h> 
#include <time.h> //measure time
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

#include <chrono>

#define MAX_LINE_LENGTH 1024
#define epsilon 1e-6 

int main() 
{

	using Clock = std::chrono::steady_clock;

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

	//================= Estimate the filtering time for each patient ==========================================

	dataType** imageData = new dataType*[Height];
	//dataType** testImageData = new dataType * [Height];
	dataType** potential = new dataType * [Height];
	dataType** action = new dataType * [Height];
	for (k = 0; k < Height; k++) 
	{
		imageData[k] = new dataType[dim2D]{0};
		//testImageData[k] = new dataType[dim2D]{0};
		potential[k] = new dataType[dim2D]{0};
		action[k] = new dataType[dim2D]{0};
	}

	//storing_path = root + "output/Data journal paper submission/P1/filtered_p1.raw";
	////storing_path = root + "output/filtered_p5.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, storing_path.c_str(), LOAD_DATA, false);
	
	/*
	//copy
	for(k = 0; k < Height; k++)
	{
		memcpy(imageData[k], ctContainer->dataPointer[k], dim2D * sizeof(dataType));
	}

	//Find min and max values in the image
	dataType minValue = imageData[0][0];
	dataType maxValue = imageData[0][0];
	for(k = 0; k < Height; k++)
	{
		for(i = 0; i < dim2D; i++)
		{
			if(imageData[k][i] < minValue) minValue = imageData[k][i];
			if(imageData[k][i] > maxValue) maxValue = imageData[k][i];
		}
	}
	std::cout << "Min value : " << minValue << std::endl;
	std::cout << "Max value : " << maxValue << std::endl;

	rescaleNewRange(imageData, Length, Width, Height, 0.0, 1.0, maxValue, minValue);

	storing_path = root + "output/rescaled_p3.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	Filter_Parameters smoothParameters =
	{
		1.0,// tau;
		1.0,// h not used here
		1.0,// sigma
		1000,// edge_detector_coefficient
		1.4,// omega_c;
		1e-3,// tolerance;
		1e-6,// eps2
		1,// coef
		1,// linked to sigma
		1,//number of time step;
		100// max number of iteration;
	};
	Image_Data inputImage = { Height, Length, Width, imageData, imageOrigin, imageSpacing, orientation };

	//clock_t start_p = clock();
	geodesicMeanCurvature(inputImage, smoothParameters);
	//clock_t end_p = clock();
	//std::cout << "Filtering time: " << double(end_p - start_p) / CLOCKS_PER_SEC << " seconds" << std::endl;

	std::vector<double> ex_times;
	for (size_t n = 0; n < 10; ++n)
	{
		manageRAWFile3D<dataType>(imageData, Length, Width, Height, storing_path.c_str(), LOAD_DATA, false);
		
		auto start = std::chrono::steady_clock::now();
		geodesicMeanCurvature(inputImage, smoothParameters);
		auto end = std::chrono::steady_clock::now();

		ex_times.push_back(std::chrono::duration<double>(end - start).count());
	}

	dataType mean_time = 0.0, max_time = 0.0, min_time = 10000.0;
	for(size_t n = 0; n < ex_times.size(); n++)
	{
		mean_time += ex_times[n];
		if(ex_times[n] > max_time) max_time = ex_times[n];
		if(ex_times[n] < min_time) min_time = ex_times[n];
	}
	//std::cout << "Filtering time for 10 iterations : " << mean_time << " seconds" << std::endl;
	mean_time /= (dataType)ex_times.size();
	std::cout << "Mean filtering time: " << mean_time << " seconds" << std::endl;
	std::cout << "Max filtering time: " << max_time << " seconds" << std::endl;
	std::cout << "Min filtering time: " << min_time << " seconds" << std::endl;
	
	storing_path = root + "output/filtered_p3.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);
	
	storing_path = root + "output/Data journal paper submission/P3/filtered_p3.raw";
	manageRAWFile3D<dataType>(testImageData, Length, Width, Height, storing_path.c_str(), LOAD_DATA, false);

	//Compute similarity between the two images
	size_t nb_voxel_different = 0;
	dataType eps = 1e-6, max_diff = 0;
	for (k = 0; k < Height; k++) 
	{
		for(i = 0; i < dim2D; i++)
		{
			max_diff = fabs(imageData[k][i] - testImageData[k][i]);
			if(max_diff > eps)
			{
				nb_voxel_different++;
			}
		}
	}
	std::cout << "Number of different voxels between the two images: " << nb_voxel_different << std::endl;
	std::cout<< "Maximum difference between the two images: " << max_diff << std::endl;
	*/

	//================= Estimate the potential computation duration for each patient ===============================

	/*
	storing_path = root + "output/Data journal paper submission/P5/filtered_p5.raw";
	////storing_path = root + "output/filtered_p5.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, storing_path.c_str(), LOAD_DATA, false);

	////Patient 1
	//Point3D seed1 = { 261, 257, 145 };
	//Point3D seed2 = { 259, 250, 246 };

	////Patient 2
	//Point3D seed1 = { 260, 254, 246 };
	//Point3D seed2 = { 255, 245, 350 };

	////Patient 3
	//Point3D seed1 = { 268, 231, 112 };
	//Point3D seed2 = { 254, 224, 221 };

	////Patient 4
	//Point3D seed1 = { 279.0, 229.0, 134.0 };
	//Point3D seed2 = { 280.0, 235.0, 223.0 };

	//Patient 5
	Point3D seed1 = { 265, 244, 470 };
	Point3D seed2 = { 235, 219, 626 };

	//////Patient 6
	//Point3D seed1 = { 249.0, 299.0, 256.0 };
	//Point3D seed2 = { 258.0, 286.0, 443.0 };

	Point3D* endPoints = new Point3D[2];
	endPoints[0] = seed1;
	endPoints[1] = seed2;
	//endPoints[2] = { 0.0, 0.0, 0.0 };

	//seed1 = getRealCoordFromImageCoord3D(seed1, imageOrigin, imageSpacing, orientation);
	//seed2 = getRealCoordFromImageCoord3D(seed2, imageOrigin, imageSpacing, orientation);

	//std::cout << "Seed 1 in real world coord : " << seed1.x << "," << seed1.y << "," << seed1.z << std::endl;
	//std::cout << "Seed 2 in real world coord : " << seed2.x << "," << seed2.y << "," << seed2.z << std::endl;
	
	double radius = 3.0;
	Potential_Parameters parameters
	{
		1000, //edge detector coefficient
		0.25, //threshold
		0.001,//epsilon
		radius
	};
	Image_Data inputImageStr = { Height, Length, Width, imageData, imageOrigin, imageSpacing, orientation };
	
	//clock_t start_p = clock();
	//compute3DPotential(inputImageStr, potential, endPoints, parameters);
	//clock_t end_p = clock();
	//std::cout << "Potential computation time: " << double(end_p - start_p) / CLOCKS_PER_SEC << " seconds" << std::endl;

	
	std::vector<double> ex_times;
	for (size_t n = 0; n < 10; ++n)
	{
		auto start = std::chrono::steady_clock::now();
		compute3DPotential(inputImageStr, potential, endPoints, parameters);
		auto end = std::chrono::steady_clock::now();
		ex_times.push_back(std::chrono::duration<double>(end - start).count());
	}

	dataType mean_time = 0.0, max_time = 0.0, min_time = 10000.0;
	for (size_t n = 0; n < ex_times.size(); n++)
	{
		mean_time += ex_times[n];
		if (ex_times[n] > max_time) max_time = ex_times[n];
		if (ex_times[n] < min_time) min_time = ex_times[n];
	}
	//std::cout << "Filtering time for 10 iterations : " << mean_time << " seconds" << std::endl;
	mean_time /= (dataType)ex_times.size();
	std::cout << "Mean time: " << mean_time << " seconds" << std::endl;
	std::cout << "Max time: " << max_time << " seconds" << std::endl;
	std::cout << "Min time: " << min_time << " seconds" << std::endl;
	string saving_csv = root + "output/execution_time_potential_p4.csv";
	FILE* f_exc_potential;
	if (fopen_s(&f_exc_potential, saving_csv.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	fprintf(f_exc_potential, "Mean time, Max time, Min time\n");
	fprintf(f_exc_potential, "%f, %f, %f\n", mean_time, max_time, min_time);
	fclose(f_exc_potential);	
	

	storing_path = root + "output/Data journal paper submission/P5/potential_p5.raw";
	////storing_path = root + "output/potential_new_p5.raw";
	manageRAWFile3D<dataType>(potential, Length, Width, Height, storing_path.c_str(), LOAD_DATA, false);

	
	////storing_path = root + "output/Data journal paper submission/P1/filtered_p1.raw";
	//storing_path = root + "output/Data journal paper submission/P4/potential_p4.raw";
	//manageRAWFile3D<dataType>(testImageData, Length, Width, Height, storing_path.c_str(), LOAD_DATA, false);

	////Compute similarity between the two images
	//size_t nb_voxel_different = 0;
	//dataType eps = 1e-6, max_diff = 0;
	//for (k = 0; k < Height; k++)
	//{
	//	for (i = 0; i < dim2D; i++)
	//	{
	//		max_diff = fabs(potential[k][i] - testImageData[k][i]);
	//		if (max_diff > eps)
	//		{
	//			nb_voxel_different++;
	//		}
	//	}
	//}
	//std::cout << "Number of different voxels between the two images: " << nb_voxel_different << std::endl;
	//std::cout << "Maximum difference between the two images: " << max_diff << std::endl;
	//
	*/

	//================= Estimate time for front propagation and path finding for each patient ========================

	/*
	//loading_path = root + "output/potential_p6.raw";
	//manageRAWFile3D<dataType>(potential, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	
	Image_Data actionPtr = { Height, Length, Width, action, imageOrigin, imageSpacing, orientation };
	//partialFrontPropagation(actionPtr, potential, endPoints);

	vector<Point3D> key_points;
	frontPropagationWithKeyPointDetection(actionPtr, potential, endPoints, 50.0, key_points);

	Path_Parameters parameters_path
	{
		0.8, // tau
		1000,// max number of iterations
		0.8 // tolerance
	};

	//vector<Point3D> path_points;
	//Image_Data toPathExtraction = { Height, Length, Width, action, imageOrigin, imageSpacing, orientation };
	//shortestPath3D(toPathExtraction, endPoints, path_points, parameters_path);

	
	//FILE* path_points_file;
	//string save_path_file = root + "output/path_points_test_p5.csv";
	//if (fopen_s(&path_points_file, save_path_file.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(path_points_file, "x,y,z\n");

	//for(size_t it = 0; it < path_points.size(); it++) 
	//{
	//	//convert to real world coordinates
	//	path_points[it] = getRealCoordFromImageCoord3D(path_points[it], imageOrigin, imageSpacing, orientation);
	//	fprintf(path_points_file, "%f,%f,%f\n", path_points[it].x, path_points[it].y, path_points[it].z);
	//}
	//fclose(path_points_file);
	
	//vector<Point3D> path_points;
	//Image_Data toPathExtraction = { Height, Length, Width, action, imageOrigin, imageSpacing, orientation };
	//shortestPath3D(toPathExtraction, endPoints, path_points, parameters_path);

	
	//FILE* key_points_file;
	//string save_path_file = root + "output/key_points_p5.csv";
	//if (fopen_s(&key_points_file, save_path_file.c_str(), "w") != 0) 
	//{
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(key_points_file, "x,y,z\n");

	//for (size_t it = 0; it < key_points.size(); it++)
	//{
	//	//convert to real world coordinates
	//	key_points[it] = getRealCoordFromImageCoord3D(key_points[it], imageOrigin, imageSpacing, orientation);
	//	fprintf(key_points_file, "%f,%f,%f\n", key_points[it].x, key_points[it].y, key_points[it].z);
	//}
	//fclose(key_points_file);
	

	//string new_storing_path;
	vector<Point3D> path_key_points, path_temporary;
	Image_Data toPathExtraction = { Height, Length, Width, action, imageOrigin, imageSpacing, orientation };
	//FILE* path_points_file;
	//////string save_path_file = outputPath + "Data journal paper submission/end_points_p1.csv";
	//string save_path_file = root + "output/path_points_kp_p5.csv";
	//if (fopen_s(&path_points_file, save_path_file.c_str(), "w") != 0) 
	//{
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(path_points_file, "x,y,z\n");
	
	//storing_path = root + "P5/partial/front_";
	for(int i_n = key_points.size() - 1; i_n > 0; i_n--)
	{
		endPoints[0] = key_points[i_n];
		endPoints[1] = key_points[i_n - 1];
		partialFrontPropagation(actionPtr, potential, endPoints);
		shortestPath3D(toPathExtraction, endPoints, path_temporary, parameters_path);
		for(int it = path_temporary.size() - 1; it > -1; it--)
		{
			path_key_points.push_back(path_temporary[it]);
			////convert to real world coordinates
			//path_temporary[it] = getRealCoordFromImageCoord3D(path_temporary[it], imageOrigin, imageSpacing, orientation);
			//fprintf(path_points_file, "%f,%f,%f\n", path_temporary[it].x, path_temporary[it].y, path_temporary[it].z);
		}
		//path_key_points.clear();
		path_temporary.clear();
	}
	//fclose(path_points_file);

	////reverse the path points vector to have the path from the first key point to the last key point
	//vector<Point3D> path_ordered;
	//for (size_t in = 0; in < path_key_points.size(); in++)
	//{
	//	path_ordered.push_back(path_key_points[in]);
	//}

	std::cout << "Number of path points: " << path_key_points.size() << std::endl;

	
	size_t count_iter = 0, id_path = 0, current_id = 0;
	for(size_t i_n = 0; i_n < path_key_points.size(); i_n++)
	{
		count_iter++;
		if (count_iter % 10 == 0)
		{
			id_path++;
			string save_points = root + "output/Animation 21-06/P5/Succes/path csv/path_key_points_" + std::to_string(id_path) + ".csv";
			FILE* path_points_save;
			if (fopen_s(&path_points_save, save_points.c_str(), "w") != 0)
			{
				printf("Enable to open");
				return false;
			}
			fprintf(path_points_save, "x,y,z\n");
			for (size_t p = 0; p < i_n; p++)
			{
				Point3D processed = getRealCoordFromImageCoord3D(path_key_points[p], imageOrigin, imageSpacing, orientation);
				fprintf(path_points_save, "%f,%f,%f\n", processed.x, processed.y, processed.z);
			}
			fclose(path_points_save);
			//current_id = i_n;
		}
	}

	id_path++;
	for (size_t i_n = 0; i_n < path_key_points.size(); i_n++)
	{
		string save_points = root + "output/Animation 21-06/P5/Succes/path csv/path_key_points_" + std::to_string(id_path) + ".csv";
		FILE* path_points_save;
		if (fopen_s(&path_points_save, save_points.c_str(), "w") != 0)
		{
			printf("Enable to open");
			return false;
		}
		fprintf(path_points_save, "x,y,z\n");
		for (size_t p = 0; p < i_n; p++)
		{
			Point3D processed = getRealCoordFromImageCoord3D(path_key_points[p], imageOrigin, imageSpacing, orientation);
			fprintf(path_points_save, "%f,%f,%f\n", processed.x, processed.y, processed.z);
		}
		fclose(path_points_save);
	}
	

	path_key_points.clear();
	//path_ordered.clear();
	key_points.clear();
	*/

	//================= Aorta Segmentation ========================================

	/*
	////Croping : p1
	size_t i_min = 200;
	size_t j_min = 200;
	size_t k_min = 120;
	size_t length = 150;
	size_t width = 150;
	size_t height = 180;

	dataType** segImageData = new dataType * [height];
	dataType** initialSegment = new dataType * [height] { 0 };
	for (k = 0; k < height; k++) {
		segImageData[k] = new dataType[length * width]{ 0 };
		initialSegment[k] = new dataType[length * width]{ 0 };
	}

	Point3D newOrigin = { i_min, j_min, k_min };
	newOrigin = getRealCoordFromImageCoord3D(newOrigin, imageOrigin, imageSpacing, orientation);
	std::cout << "New origin : " << newOrigin.x << ", " << newOrigin.y << ", " << newOrigin.z << std::endl;
	Image_Data inputImage = { height, length, width, segImageData, newOrigin, imageSpacing, orientation };

	size_t i_ext, j_ext, k_ext;
	for(k = 0, k_ext = k_min; k < height; k++, k_ext++)
	{
		for (i = 0, i_ext = i_min; i < length; i++, i_ext++) 
		{
			for (j = 0, j_ext = j_min; j < width; j++, j_ext++) 
			{
				segImageData[k][x_new(i, j, length)] = imageData[k_ext][x_new(i_ext, j_ext, Length)];
			}
		}
	}

	
	//Generate the initial segmentation
	FILE* path_file;
	loading_path = root + "output/Segmentation/Aorta/centered paths/finalCurve_p1.csv";
	if (fopen_s(&path_file, loading_path.c_str(), "r") != 0) 
	{
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
	size_t n;
	Point3D pPoints, pCurrent;
	double dist, min_dist;
	for (k = 0; k < height; k++)
	{
		for (i = 0; i < length; i++)
		{
			for (j = 0; j < width; j++)
			{
				pCurrent = { (dataType)i, (dataType)j, (dataType)k };
				pCurrent = getRealCoordFromImageCoord3D(pCurrent, newOrigin, imageSpacing, orientation);
				min_dist = (double)(length * width);
				for (n = 0; n < path_points.size(); n++)
				{
					pPoints = path_points[n];
					dist = getPoint3DDistance(pPoints, pCurrent);
					if (min_dist > dist)
					{
						min_dist = dist;
					}
				}
				initialSegment[k][x_new(i, j, length)] = 1.0 / (1.0 + min_dist);
			}
		}
	}
	

	storing_path = root + "output/Segmentation 21-08/initial_segment_p1.raw";
	manageRAWFile3D<dataType>(initialSegment, length, width, height, storing_path.c_str(), LOAD_DATA, false);

	Filter_Parameters smoothParameters =
	{
		1.0,// tau;
		1.0,// h not used here
		1.0,// sigma
		1000,// edge_detector_coefficient
		1.4,// omega_c;
		1e-3,// tolerance;
		1e-6,// eps2
		1,// coef
		1,// linked to sigma
		1,//number of time step;
		100// max number of iteration;
	};

	Segmentation_Parameters segParameters =
	{
		30,//Maximum number of Gauss-Seidel iterations
		100000,//edge detector coef
		1e-6,//epsilon is the regularization factor (Evans-Spruck)
		2000,//Number of current time step
		2000,//Maximum number of time step
		10,//saving frequency
		1e-6,//segmentation tolerance
		0.5,//tau
		1,//h
		1.4,//omega_c
		0.001,//tolerance
		1.0,//convection coef
		0.05,//diffusion coef
	};

	storing_path = root + "output/Segmentation 21-08/test/";
	Image_Data imageToSegment = { height, length, width, segImageData, newOrigin, imageSpacing, orientation };
	
	clock_t start_p = clock();
	GSUBSURF(imageToSegment, initialSegment, storing_path.c_str(), smoothParameters, segParameters);
	clock_t end_p = clock();
	std::cout << "GSUBSURF time: " << double(end_p - start_p) / CLOCKS_PER_SEC << " seconds" << std::endl;

	for (k = 0; k < height; k++) {
		delete[] segImageData[k];
		delete[] initialSegment[k];
	}
	delete[] segImageData;
	delete[] initialSegment;
	*/
	
	//================== Old versus new potential ==================================
	
	/*
	//Patient 1
	Point3D seed1 = { 261, 257, 145 };
	Point3D seed2 = { 259, 250, 246 };
	Point3D* endPoints = new Point3D[2];
	endPoints[0] = seed1;
	endPoints[1] = seed2;

	double radius = 3.0;
	Potential_Parameters parameters
	{
		1000, //edge detector coefficient
		0.2, //threshold
		0.001,//epsilon
		radius
	};
	Image_Data inputImageStr = { Height, Length, Width, imageData, imageOrigin, imageSpacing, orientation };

	compute3DPotential(inputImageStr, potential, endPoints, parameters);

	//storing_path = root + "output/potential_new.raw";
	//manageRAWFile3D<dataType>(potential, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	Image_Data actionPtr = { Height, Length, Width, action, imageOrigin, imageSpacing, orientation };
	partialFrontPropagation(actionPtr, potential, endPoints);

	Path_Parameters parameters_path
	{
		0.8, // tau
		1000,// max number of iterations
		0.8 // tolerance
	};

	vector<Point3D> path_points;
	Image_Data toPathExtraction = { Height, Length, Width, action, imageOrigin, imageSpacing, orientation };
	shortestPath3D(toPathExtraction, endPoints, path_points, parameters_path);

	FILE* path_points_file;
	string save_path_file = root + "output/path_points_old_3.csv";
	if (fopen_s(&path_points_file, save_path_file.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	fprintf(path_points_file, "x,y,z\n");

	for(size_t it = 0; it < path_points.size(); it++) 
	{
		//convert to real world coordinates
		path_points[it] = getRealCoordFromImageCoord3D(path_points[it], imageOrigin, imageSpacing, orientation);
		fprintf(path_points_file, "%f,%f,%f\n", path_points[it].x, path_points[it].y, path_points[it].z);
	}
	fclose(path_points_file);

	path_points.clear();
	*/

	//================== Quantitative analysis ====================================
	
	/*
	storing_path = root + "output/Data journal paper submission/P1/segment_aorta_full_dim_p1.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, storing_path.c_str(), LOAD_DATA, false);

	storing_path = root + "output/Data journal paper submission/P1/distance_map_aorta_p1.raw";
	manageRAWFile3D<dataType>(potential, Length, Width, Height, storing_path.c_str(), LOAD_DATA, false);

	Point3D center = { 42.760132, 140.8125, -387.464905 };
	Point3D centerImageCoord = getImageCoordFromRealCoord3D(center, imageOrigin, imageSpacing, orientation);
	size_t ic = (size_t)round(centerImageCoord.x);
	size_t jc = (size_t)round(centerImageCoord.y);
	size_t kc = (size_t)round(centerImageCoord.z);
	double radius = (double)(potential[kc][x_new(ic, jc, Length)] + 2 * 4.0);

	for(k = 0; k < Height; k++)
	{
		for(i = 0; i < Length; i++)
		{
			for(j = 0; j < Width; j++)
			{
				Point3D currentPoint = { (dataType)i, (dataType)j, (dataType)k };
				currentPoint = getRealCoordFromImageCoord3D(currentPoint, imageOrigin, imageSpacing, orientation);
				double dist = getPoint3DDistance(center, currentPoint);
				if (dist <= radius) 
				{
					action[k][x_new(i, j, Length)] = 1.0;
					if(imageData[k][x_new(i, j, Length)] > 0.5)
					{
						testImageData[k][x_new(i, j, Length)] = 1.0;
					}
				}

			}
		}
	}

	storing_path = root + "output/ball.raw";
	manageRAWFile3D<dataType>(action, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	storing_path = root + "output/piece_aorta.raw";
	manageRAWFile3D<dataType>(testImageData, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);
	*/ 
	//================== Segment the liver ============================

	int** labelArray = new int* [Height];
	bool** status = new bool* [Height];
	for(k = 0; k < Height; k++)
	{
		labelArray[k] = new int[dim2D]{ 0 };
		status[k] = new bool[dim2D] { false };
	}

	//dataType thres_min = -30, thres_max = 200;
	//copy the data
	//for(k = 0; k < Height; k++)
	//{
	//	for (i = 0; i < dim2D; i++)
	//	{
	//		imageData[k][i] = (dataType)ctContainer->dataPointer[k][i];
	//	}
	//}
	//thresholding3dFunctionN(imageData, Length, Width, Height, thres_min, thres_max, 0.0, 1.0);
	//storing_path = root + "output/threshold.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);
	//dataType min_radius = 1.5 * fmax(imageSpacing.sx, fmax(imageSpacing.sy, imageSpacing.sz));
	//for (size_t n = 0; n < 5; n++) 
	//{
	//	Image_Data inputImage = { Height, Length, Width, imageData, imageOrigin, imageSpacing, orientation };
	//	fastMarchingDistanceMap(inputImage, action, 0.0);
	//	for (k = 0; k < Height; k++)
	//	{
	//		for (i = 0; i < dim2D; i++)
	//		{
	//			if (action[k][i] <= min_radius)
	//			{
	//				imageData[k][i] = 0.0;
	//			}
	//			else
	//			{
	//				imageData[k][i] = 1.0;
	//			}
	//		}
	//	}
	//}
	//storing_path = root + "output/liver_segmented.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	
	//copying data
	for (k = 0; k < Height; k++)
	{
		for (i = 0; i < dim2D; i++)
		{
			imageData[k][i] = (dataType)ctContainer->dataPointer[k][i];
		}
	}

	Image_Data toDistanceMap = { Height, Length, Width, imageData, imageOrigin, imageSpacing, orientation };
	dataType min_distance = 1.5 * fmax(imageSpacing.sz, fmax(imageSpacing.sy, imageSpacing.sx));

	dataType thres_min = -30, thres_max = 200, background = 0, foreground = 1;
	thresholding3dFunctionN(imageData, Length, Width, Height, thres_min, thres_max, background, foreground);

	//storing_path = root + "output/threshold.raw";
	//manageRAWFile3D<dataType>(segmentedLiver, Length, Width, newHeight, storing_path.c_str(), STORE_DATA, false);

	int numberOfRegionCells = countNumberOfRegionsCells(imageData, Length, Width, Height, foreground);
	std::cout << "Initial volume size = " << numberOfRegionCells << std::endl;

	labelling3D(imageData, labelArray, status, Length, Width, Height, foreground);

	//Counting
	int* countingArray = new int[numberOfRegionCells] {0};
	if (countingArray == NULL)
		return false;

	//Get regions sizes
	for (k = 0; k < Height; k++)
	{
		for (i = 0; i < dim2D; i++)
		{
			if (labelArray[k][i] > 0)
			{
				countingArray[labelArray[k][i]]++;
			}
		}
	}

	//Find the largest region
	int largestRegionSize = 0;
	for (k = 0; k < Height; k++)
	{
		for (i = 0; i < dim2D; i++)
		{
			if (countingArray[labelArray[k][i]] > largestRegionSize)
			{
				largestRegionSize = countingArray[labelArray[k][i]];
			}
		}
	}

	int initialLargestRegionSize = largestRegionSize;
	//std::cout << "Initial largest volume size = " << initialLargestRegionSize << std::endl;

	//Keep the only the largest region
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++)
		{
			if (countingArray[labelArray[k][i]] == largestRegionSize)
			{
				imageData[k][i] = foreground;
			}
			else {
				imageData[k][i] = background;
			}

			//Initialize labelArray and status for the next iteration
			labelArray[k][i] = 0;
			status[k][i] = false;
		}
	}

	delete[] countingArray;

	size_t num_iterations = 5;
	dataType max_distance = 0.0;
	size_t iteration = 0;

	while (iteration < num_iterations)
	{

		iteration++;
		fastMarchingDistanceMap(toDistanceMap, action, background);
		for (k = 0; k < Height; k++)
		{
			for (i = 0; i < dim2D; i++)
			{
				if (action[k][i] <= min_distance)
				{
					imageData[k][i] = background;
				}
			}
		}

		numberOfRegionCells = countNumberOfRegionsCells(imageData, Length, Width, Height, foreground);
		labelling3D(imageData, labelArray, status, Length, Width, Height, foreground);

		//Counting
		int* countingArray = new int[numberOfRegionCells] {0};
		if (countingArray == NULL)
			return false;

		//Get regions sizes
		for (k = 0; k < Height; k++)
		{
			for (i = 0; i < dim2D; i++)
			{
				if (labelArray[k][i] > 0)
				{
					countingArray[labelArray[k][i]]++;
				}
			}
		}

		//Find the largest region
		largestRegionSize = 0;
		for (k = 0; k < Height; k++)
		{
			for (i = 0; i < dim2D; i++)
			{
				if (countingArray[labelArray[k][i]] > largestRegionSize)
				{
					largestRegionSize = countingArray[labelArray[k][i]];
				}
			}
		}

		//Keep the only the largest region
		for (k = 0; k < Height; k++)
		{
			for (i = 0; i < dim2D; i++)
			{
				if (countingArray[labelArray[k][i]] == largestRegionSize)
				{
					imageData[k][i] = foreground;
				}
				else
				{
					imageData[k][i] = background;
				}
			}
		}

		//Reinitialization
		for (k = 0; k < Height; k++)
		{
			for (i = 0; i < dim2D; i++)
			{
				labelArray[k][i] = 0;
				status[k][i] = false;
			}
		}

		//dataType volume_ratio = (dataType)largestRegionSize / (dataType)initialLargestRegionSize;
		//initialLargestRegionSize = largestRegionSize;
		//fprintf(ratio_file, "%zu,%f\n", iteration, volume_ratio);
		//std::cout << "Largest volume size = " << largestRegionSize << " at step " << iteration << std::endl;
		//std::cout << "Volume ratio = " << volume_ratio << " at step " << iteration << std::endl;
		//std::cout << "================================" << std::endl;

		delete[] countingArray;
	}

	storing_path = root + "output/liver_segmented.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	/*
	//invert
	for(k = 0; k < Height; k++)
	{
		for (i = 0; i < dim2D; i++)
		{
			if (imageData[k][i] == foreground)
			{
				imageData[k][i] = background;
			}
			else
			{
				imageData[k][i] = foreground;
			}
		}
	}

	//Refinement of the segmented liver
	for(int n = 0; n < 5; n++)
	{
		fastMarchingDistanceMap(toDistanceMap, action, background);
		for (k = 0; k < Height; k++)
		{
			for (i = 0; i < dim2D; i++)
			{
				if (action[k][i] > 0 && action[k][i] <= min_distance)
				{
					imageData[k][i] = foreground;
					action[k][i] = 0.0;
				}
			}
		}
	}

	storing_path = root + "output/liver_segmented_refined.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);
	*/

	for(k = 0; k < Height; k++)
	{
		delete[] labelArray[k];
		delete[] status[k];
	}
	delete[] labelArray;
	delete[] status;	

	//======================== Refinement ================================================

	//dataType** segmentedLiver = new dataType * [Height];
	//for (k = 0; k < Height; k++)
	//{
	//	segmentedLiver[k] = new dataType[dim2D]{ 0 };
	//}
	//loading_path = root + "output/segment_dmap_p5.raw";
	//manageRAWFile3D<dataType>(segmentedLiver, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	//Get bounding box of the segmented liver
	size_t k_min = Height, k_max = 0, i_min = Length, i_max = 0, j_min = Width, j_max = 0;
	for(k = 0; k < Height; k++)
	{
		for(i = 0; i < Length; i++)
		{
			for(j = 0; j < Width; j++)
			{
				xd = x_new(i, j, Length);
				if(imageData[k][xd] > 0)
				{
					if(k < k_min) k_min = k;
					if(k > k_max) k_max = k;
					if(i < i_min) i_min = i;
					if(i > i_max) i_max = i;
					if(j < j_min) j_min = j;
					if(j > j_max) j_max = j;
				}
			}
		}
	}

	i_min -= 25;
	i_max += 25;
	j_min -= 25;
	j_max += 25;
	k_min -= 25;
	k_max += 25;

	Point3D cropOrigin = { (dataType)i_min, (dataType)j_min, (dataType)k_min };
	cropOrigin = getRealCoordFromImageCoord3D(cropOrigin, imageOrigin, imageSpacing, orientation);
	std::cout << "Cropped origin : (" << cropOrigin.x << ", " << cropOrigin.y << ", " << cropOrigin.z << ")" << std::endl;

	size_t length = i_max - i_min;
	size_t width = j_max - j_min;
	size_t height = k_max - k_min;
	std::cout << "Cropped dimensions : " << length << " x " << width << " x " << height << std::endl;

	dataType** croppedLiver = new dataType * [height];
	dataType** distanceMap = new dataType * [height];
	dataType** cropCT = new dataType * [height];
	for(k = 0; k < height; k++)
	{
		croppedLiver[k] = new dataType[length * width]{ 0 };
		distanceMap[k] = new dataType[length * width]{ 0 };
		cropCT[k] = new dataType[length * width]{ 0 };
	}

	//Cropping the segmented liver
	size_t i_l, j_l, k_l, xd_l;
	for(k = 0, k_l = k_min; k < height; k++, k_l++)
	{
		for(i = 0, i_l = i_min; i < length; i++, i_l++)
		{
			for(j = 0, j_l = j_min; j < width; j++, j_l++)
			{
				xd = x_new(i, j, length);
				xd_l = x_new(i_l, j_l, Length);
				if(k_l < Height && i_l < Length && j_l < Width)
				{
					croppedLiver[k][xd] = imageData[k_l][xd_l];
					cropCT[k][xd] = ctContainer->dataPointer[k_l][xd_l];
					
					//if (imageData[k_l][xd_l] == foreground) 
					//{
					//	croppedLiver[k][xd] = background;
					//}else
					//{
					//	croppedLiver[k][xd] = foreground;
					//}
				}
			}
		}
	}

	string save_slice_root = root + "output/slice ct/";
	for (k = 0; k < height; k++)
	{
		storing_path = save_slice_root + "slice_ct_" + std::to_string(k) + ".raw";
		manageRAWFile2D<dataType>(cropCT[k], length, width, storing_path.c_str(), STORE_DATA, false);

		storing_path = save_slice_root + "slice_seg_" + std::to_string(k) + ".raw";
		manageRAWFile2D<dataType>(croppedLiver[k], length, width, storing_path.c_str(), STORE_DATA, false);

		//storing_path = save_slice_root + "slice_seg_ref_" + std::to_string(k) + ".raw";
		//manageRAWFile2D<dataType>(distanceMap[k], length, width, storing_path.c_str(), STORE_DATA, false);
	}

	//storing_path = root + "output/crop_liver_ct_p1.raw";
	//manageRAWFile3D<dataType>(cropCT, length, width, height, storing_path.c_str(), STORE_DATA, false);

	//storing_path = root + "output/refined_segment_p1.raw";
	//manageRAWFile3D<dataType>(distanceMap, length, width, height, storing_path.c_str(), LOAD_DATA, false);

	//storing_path = root + "output/crop_segment_p1.raw";
	//manageRAWFile3D<dataType>(croppedLiver, length, width, height, storing_path.c_str(), STORE_DATA, false);

	Image_Data cDistanceMap = { height, length, width, croppedLiver, cropOrigin, imageSpacing, orientation };
	min_distance = 1.0 * fmax(imageSpacing.sx, fmax(imageSpacing.sy, imageSpacing.sz));

	//size_t num_iterations = 6, iteration = 0;
	//size_t foreground = 1, background = 0;

	iteration = 0;
	while(iteration < 6)
	{
		fastMarchingDistanceMap(cDistanceMap, distanceMap, foreground);
		for (k = 0; k < height; k++)
		{
			for (i = 0; i < length * width; i++)
			{
				if (distanceMap[k][i] > 0 && distanceMap[k][i] <= min_distance)
				{
					croppedLiver[k][i] = foreground;
				}
				distanceMap[k][i] = 0.0;
			}
		}
		iteration++;
	}

	save_slice_root = root + "output/slice ct/";
	for (k = 0; k < height; k++)
	{
		//storing_path = save_slice_root + "slice_ct_" + std::to_string(k) + ".raw";
		//manageRAWFile2D<dataType>(cropCT[k], length, width, storing_path.c_str(), STORE_DATA, false);

		//storing_path = save_slice_root + "slice_seg_" + std::to_string(k) + ".raw";
		//manageRAWFile2D<dataType>(croppedLiver[k], length, width, storing_path.c_str(), STORE_DATA, false);

		storing_path = save_slice_root + "slice_seg_ref_" + std::to_string(k) + ".raw";
		manageRAWFile2D<dataType>(croppedLiver[k], length, width, storing_path.c_str(), STORE_DATA, false);
	}

	for(k = 0; k < height; k++)
	{
		delete[] croppedLiver[k];
		delete[] distanceMap[k];
		delete[] cropCT[k];
	}
	delete[] croppedLiver;
	delete[] distanceMap;
	delete[] cropCT;

	//for(k = 0; k < Height; k++)
	//{
	//	delete[] segmentedLiver[k];
	//
	//}
	//delete[] segmentedLiver;



	//================== Clean up memory ==========================================

	//delete[] endPoints;
	for(k = 0; k < Height; k++)
	{
		delete[] imageData[k];
		//delete[] testImageData[k];
		delete[] potential[k];
		delete[] action[k];
	}
	delete[] imageData;
	//delete[] testImageData;
	delete[] potential;
	delete[] action;

	free(ctContainer);

	return EXIT_SUCCESS;
}