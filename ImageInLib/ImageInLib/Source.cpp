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
	loading_path = root + "input/vtk/petct/ct/Patient5_ct.vtk";
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
	dataType** testImageData = new dataType * [Height];
	dataType** potential = new dataType * [Height];
	dataType** action = new dataType * [Height];
	for (k = 0; k < Height; k++) 
	{
		imageData[k] = new dataType[dim2D]{0};
		testImageData[k] = new dataType[dim2D]{0};
		potential[k] = new dataType[dim2D]{0};
		action[k] = new dataType[dim2D]{0};
	}
	
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

	/*
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
	*/

	storing_path = root + "output/Data journal paper submission/P5/potential_p5.raw";
	////storing_path = root + "output/potential_new_p5.raw";
	manageRAWFile3D<dataType>(potential, Length, Width, Height, storing_path.c_str(), LOAD_DATA, false);

	/*
	//storing_path = root + "output/Data journal paper submission/P1/filtered_p1.raw";
	storing_path = root + "output/Data journal paper submission/P4/potential_p4.raw";
	manageRAWFile3D<dataType>(testImageData, Length, Width, Height, storing_path.c_str(), LOAD_DATA, false);

	//Compute similarity between the two images
	size_t nb_voxel_different = 0;
	dataType eps = 1e-6, max_diff = 0;
	for (k = 0; k < Height; k++)
	{
		for (i = 0; i < dim2D; i++)
		{
			max_diff = fabs(potential[k][i] - testImageData[k][i]);
			if (max_diff > eps)
			{
				nb_voxel_different++;
			}
		}
	}
	std::cout << "Number of different voxels between the two images: " << nb_voxel_different << std::endl;
	std::cout << "Maximum difference between the two images: " << max_diff << std::endl;
	*/

	//================= Estimate time for front propagation and path finding for each patient ========================

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

	/*
	FILE* path_points_file;
	string save_path_file = root + "output/path_points_test_p5.csv";
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
	*/

	//vector<Point3D> path_points;
	//Image_Data toPathExtraction = { Height, Length, Width, action, imageOrigin, imageSpacing, orientation };
	//shortestPath3D(toPathExtraction, endPoints, path_points, parameters_path);

	/*
	FILE* key_points_file;
	string save_path_file = root + "output/key_points_p5.csv";
	if (fopen_s(&key_points_file, save_path_file.c_str(), "w") != 0) 
	{
		printf("Enable to open");
		return false;
	}
	fprintf(key_points_file, "x,y,z\n");

	for (size_t it = 0; it < key_points.size(); it++)
	{
		//convert to real world coordinates
		key_points[it] = getRealCoordFromImageCoord3D(key_points[it], imageOrigin, imageSpacing, orientation);
		fprintf(key_points_file, "%f,%f,%f\n", key_points[it].x, key_points[it].y, key_points[it].z);
	}
	fclose(key_points_file);
	*/

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
	

	//================== Clean up memory ==========================================
	delete[] endPoints;
	for(k = 0; k < Height; k++)
	{
		delete[] imageData[k];
		delete[] testImageData[k];
		delete[] potential[k];
		delete[] action[k];
	}
	delete[] imageData;
	delete[] testImageData;
	delete[] potential;
	delete[] action;

	free(ctContainer);

	return EXIT_SUCCESS;
}