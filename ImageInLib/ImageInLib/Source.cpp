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

#include "../src/heat_equation.h"
#include "segmentation2d.h"
#include "../src/segmentation3d_gsubsurf.h"
#include "../src/non_linear_heat_equation.h"

#include "gaussian_distribution.h"

int main() {

	string inputPath = "C:/Users/Konan Allaly/Documents/Tests/input/";
	string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";
	
	string loading_path, storing_path, extension;

	size_t i = 0, j = 0, k = 0, xd = 0;

	//===================== Load 3D patient data (.vtk) ================================
	
	/*
	We load the .vtk files containing the original pixels value
	and informations related the dimensions, the spacings and
	the origins. We also define the orientation matrix in 2D/3D
	needed when we need to perform interpolation.
	*/
	
	OrientationMatrix orientation = { { 1.0, 0.0, 0.0 } , { 0.0, 1.0, 0.0 } , { 0.0, 0.0, 1.0 } };
	
	Vtk_File_Info* ctContainer = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	ctContainer->operation = copyFrom;
	loading_path = inputPath + "vtk/petct/ct/Patient5_ct.vtk";
	readVtkFile(loading_path.c_str(), ctContainer);

	std::cout << "============ Input ================ " << std::endl;

	int Height = ctContainer->dimensions[2];
	int Length = ctContainer->dimensions[0];
	int Width = ctContainer->dimensions[1];
	int dim2D = Length * Width;
	std::cout << "CT image dim : " << ctContainer->dimensions[0] << " x " << ctContainer->dimensions[1] << " x " << ctContainer->dimensions[2] << "" << std::endl;

	std::cout << "CT origin : (" << ctContainer->origin[0] << ", " << ctContainer->origin[1] << ", " << ctContainer->origin[2] << ")" << std::endl;
	Point3D ctOrigin = { ctContainer->origin[0], ctContainer->origin[1], ctContainer->origin[2] };
	VoxelSpacing ctSpacing = { ctContainer->spacing[0], ctContainer->spacing[1], ctContainer->spacing[2] };
	std::cout << "CT spacing : (" << ctContainer->spacing[0] << ", " << ctContainer->spacing[1] << ", " << ctContainer->spacing[2] << ")" << std::endl; 
	

	//==================== Compute Hausdoff distance and Ratio =======================================
	
	/*
	dataType img_f, n0, n1, n2, nmg, x, y, z, ptmg, ptid, scal;
	char header[MAX_LINE_LENGTH];
	string inputPointCloud = inputPath + "vtk/petct/aorta/Hausdoff 31-03/Isolines/";

	FILE* file_hausdoff;
	storing_path = inputPointCloud + "h distance/hausdoff_distance_patient2.csv";
	if (fopen_s(&file_hausdoff, storing_path.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	fprintf(file_hausdoff, "manual,gsubsurf,HD,ratio,mean_manual,mean_gsubsurf,MHD,ratio_mean\n");

	
	FILE* file_manual;
	loading_path = inputPointCloud + "manual/_patient_2.csv";
	if (fopen_s(&file_manual, loading_path.c_str(), "r") != 0) {
		printf("Enable to open");
		return false;
	}
	fgets(header, MAX_LINE_LENGTH, file_manual);
	fscanf_s(file_manual, ",");

	
	vector<Point3D> points_manual;
	while (feof(file_manual) == 0) {

		fscanf_s(file_manual, "%f", &n0);
		fscanf_s(file_manual, ",");
		fscanf_s(file_manual, "%f", &n1);
		fscanf_s(file_manual, ",");
		fscanf_s(file_manual, "%f", &n2);
		fscanf_s(file_manual, ",");
		fscanf_s(file_manual, "%f", &nmg);
		fscanf_s(file_manual, ",");

		fscanf_s(file_manual, "%f", &x);
		fscanf_s(file_manual, ",");
		fscanf_s(file_manual, "%f", &y);
		fscanf_s(file_manual, ",");
		fscanf_s(file_manual, "%f", &z);
		fscanf_s(file_manual, ",");

		fscanf_s(file_manual, "%f", &ptmg);
		fscanf_s(file_manual, ",");
		fscanf_s(file_manual, "%f", &scal);
		fscanf_s(file_manual, ",");
		fscanf_s(file_manual, "%f", &ptid);

		fscanf_s(file_manual, "\n");

		Point3D current_point = { x, y, z };
		points_manual.push_back(current_point);
	}
	fclose(file_manual);



	vector<Point3D>points_gsubsurf;

	FILE* file_gsubsurf;
	//std::string path_root = inputPointCloud + "gsubsurf/_02_patient_1.csv";
	std::string path_root = inputPointCloud + "gsubsurf/_02_patient_2.csv";
	//std::string path_root = inputPointCloud + "gsubsurf/_03_patient_3.csv";
	//std::string path_root = inputPointCloud + "gsubsurf/_03_patient_4.csv";
	//std::string path_root = inputPointCloud + "gsubsurf/_017_patient_5.csv";
	//std::string path_root = inputPointCloud + "gsubsurf/_03_patient_6.csv";
	if (fopen_s(&file_gsubsurf, path_root.c_str(), "r") != 0) {
		printf("Enable to open the file");
		return false;
	}
	fgets(header, MAX_LINE_LENGTH, file_gsubsurf);
	fscanf_s(file_gsubsurf, ",");

	while (feof(file_gsubsurf) == 0) {

		fscanf_s(file_gsubsurf, "%f", &img_f);
		fscanf_s(file_gsubsurf, ",");
		fscanf_s(file_gsubsurf, "%f", &n0);
		fscanf_s(file_gsubsurf, ",");
		fscanf_s(file_gsubsurf, "%f", &n1);
		fscanf_s(file_gsubsurf, ",");
		fscanf_s(file_gsubsurf, "%f", &n2);
		fscanf_s(file_gsubsurf, ",");
		fscanf_s(file_gsubsurf, "%f", &nmg);
		fscanf_s(file_gsubsurf, ",");

		fscanf_s(file_gsubsurf, "%f", &x);
		fscanf_s(file_gsubsurf, ",");
		fscanf_s(file_gsubsurf, "%f", &y);
		fscanf_s(file_gsubsurf, ",");
		fscanf_s(file_gsubsurf, "%f", &z);
		fscanf_s(file_gsubsurf, ",");

		fscanf_s(file_gsubsurf, "%f", &ptmg);
		fscanf_s(file_gsubsurf, ",");
		fscanf_s(file_gsubsurf, "%f", &ptid);

		fscanf_s(file_gsubsurf, "\n");

		Point3D current_point = { x, y, z };
		points_gsubsurf.push_back(current_point);
	}
	fclose(file_gsubsurf);

	size_t count_similar = 0;

	double max_manual = 0.0;
	double max_gsubsurf = 0.0;
	double mean_manual = 0.0;
	double mean_gsubsurf = 0.0;
	for (int mn = 0; mn < points_manual.size(); mn++) {
		double min_distance = 1000000000000000000000.0;
		for (int gs = 0; gs < points_gsubsurf.size(); gs++) {
			double pDistance = getPoint3DDistance(points_manual[mn], points_gsubsurf[gs]);
			if (min_distance > pDistance) {
				min_distance = pDistance;
			}
		}
		if (max_manual < min_distance) {
			max_manual = min_distance;
		}
		if (min_distance <= 1.171875) {
			count_similar++;
		}
		mean_manual += min_distance;
	}

	for (int gs = 0; gs < points_gsubsurf.size(); gs++) {
		double min_distance = 1000000000000000000000.0;
		for (int mn = 0; mn < points_manual.size(); mn++) {
			double pDistance = getPoint3DDistance(points_manual[mn], points_gsubsurf[gs]);
			if (min_distance > pDistance) {
				min_distance = pDistance;
			}
		}
		if (max_gsubsurf < min_distance) {
			max_gsubsurf = min_distance;
		}
		mean_gsubsurf += min_distance;
	}

	size_t count_manual = points_manual.size();
	size_t count_gsubsurf = points_gsubsurf.size();
	mean_manual /= count_manual;
	mean_gsubsurf /= count_gsubsurf;
	dataType ratio = 0.0;
	dataType ratio_mean = 0.0;
	ratio = max_gsubsurf / max_manual;
	ratio_mean = mean_gsubsurf / mean_manual;

	dataType HD = 0, MHD = 0;
	if (max_manual >= max_gsubsurf) {
		HD = max_manual;
	}
	else {
		HD = max_gsubsurf;
	}
	if (mean_manual >= mean_gsubsurf) {
		MHD = mean_manual;
	}
	else {
		MHD = mean_gsubsurf;
	}

	std::cout << "We have : " << count_manual << " point in the manual segmentation" << std::endl;
	std::cout << "We have : " << count_gsubsurf << " point in the segmentation" << std::endl;
	std::cout << "We have : " << count_similar << " similar_points" << std::endl;

	fprintf(file_hausdoff, "%f,%f,%f,%f,%f,%f,%f,%f\n", max_manual, max_gsubsurf, HD, ratio, mean_manual, mean_gsubsurf, MHD, ratio_mean);
	fclose(file_hausdoff);
	*/

	/*
	for (size_t n = 1; n <= 40; n++) {

		FILE* file_gsubsurf;
		loading_path = path_root + to_string(n) + ".csv";
		if (fopen_s(&file_gsubsurf, loading_path.c_str(), "r") != 0) {
			printf("Enable to open the file");
			return false;
		}
		fgets(header, MAX_LINE_LENGTH, file_gsubsurf);
		fscanf_s(file_gsubsurf, ",");

		while (feof(file_gsubsurf) == 0) {

			fscanf_s(file_gsubsurf, "%f", &img_f);
			fscanf_s(file_gsubsurf, ",");
			fscanf_s(file_gsubsurf, "%f", &n0);
			fscanf_s(file_gsubsurf, ",");
			fscanf_s(file_gsubsurf, "%f", &n1);
			fscanf_s(file_gsubsurf, ",");
			fscanf_s(file_gsubsurf, "%f", &n2);
			fscanf_s(file_gsubsurf, ",");
			fscanf_s(file_gsubsurf, "%f", &nmg);
			fscanf_s(file_gsubsurf, ",");

			fscanf_s(file_gsubsurf, "%f", &x);
			fscanf_s(file_gsubsurf, ",");
			fscanf_s(file_gsubsurf, "%f", &y);
			fscanf_s(file_gsubsurf, ",");
			fscanf_s(file_gsubsurf, "%f", &z);
			fscanf_s(file_gsubsurf, ",");

			fscanf_s(file_gsubsurf, "%f", &ptmg);
			fscanf_s(file_gsubsurf, ",");
			fscanf_s(file_gsubsurf, "%f", &ptid);

			fscanf_s(file_gsubsurf, "\n");

			Point3D current_point = { x, y, z };
			points_gsubsurf.push_back(current_point);
		}
		fclose(file_gsubsurf);

		double max_manual = 0.0;
		double max_gsubsurf = 0.0;
		double mean_manual = 0.0;
		double mean_gsubsurf = 0.0;
		for (int mn = 0; mn < points_manual.size(); mn++) {
			double min_distance = 1000000000000000000000.0;
			for (int gs = 0; gs < points_gsubsurf.size(); gs++) {
				double pDistance = getPoint3DDistance(points_manual[mn], points_gsubsurf[gs]);
				if (min_distance > pDistance) {
					min_distance = pDistance;
				}
			}
			if (max_manual < min_distance) {
				max_manual = min_distance;
			}
			mean_manual += min_distance;
		}

		for (int gs = 0; gs < points_gsubsurf.size(); gs++) {
			double min_distance = 1000000000000000000000.0;
			for (int mn = 0; mn < points_manual.size(); mn++) {
				double pDistance = getPoint3DDistance(points_manual[mn], points_gsubsurf[gs]);
				if (min_distance > pDistance) {
					min_distance = pDistance;
				}
			}
			if (max_gsubsurf < min_distance) {
				max_gsubsurf = min_distance;
			}
			mean_gsubsurf += min_distance;
		}
		
		size_t count_manual = points_manual.size();
		size_t count_gsubsurf = points_gsubsurf.size();
		mean_manual /= count_manual;
		mean_gsubsurf /= count_gsubsurf;
		dataType ratio = 0.0;
		dataType ratio_mean = 0.0;
		//if (max_manual < max_gsubsurf) {
		//	ratio = max_manual / max_gsubsurf;
		//}
		//else {
		//	ratio = max_gsubsurf / max_manual;
		//}
		ratio = max_gsubsurf / max_manual;
		ratio_mean = mean_gsubsurf / mean_manual;
		
		dataType HD = 0, MHD = 0;
		if (max_manual >= max_gsubsurf) {
			HD = max_manual;
		}
		else {
			HD = max_gsubsurf;
		}
		if (mean_manual >= mean_gsubsurf) {
			MHD = mean_manual;
		}
		else {
			MHD = mean_gsubsurf;
		}
		fprintf(file_hausdoff, "%f,%f,%f,%f,%f,%f,%f,%f\n", max_manual, max_gsubsurf, HD, ratio, mean_manual, mean_gsubsurf, MHD, ratio_mean);

		while (points_gsubsurf.size() > 0)
		{
			points_gsubsurf.pop_back();
		}

	}
	fclose(file_hausdoff);
	*/

	//==================== Test Potential function ====================================================
	
	/*
	dataType** imageData = new dataType * [Height];
	dataType** toSmoothImageData = new dataType * [Height];
	dataType** action = new dataType * [Height];
	dataType** potential = new dataType * [Height];
	dataType** actionFirstFront = new dataType * [Height];
	dataType** actionSecondFront = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{ 0 };
		toSmoothImageData[k] = new dataType[dim2D]{ 0 };
		action[k] = new dataType[dim2D]{ 0 };
		potential[k] = new dataType[dim2D]{ 0 };
		actionFirstFront[k] = new dataType[dim2D]{ 0 };
		actionSecondFront[k] = new dataType[dim2D]{ 0 };
	}
	////loading_path = inputPath + "raw/filtered/New/filtered_p1.raw";
	//loading_path = inputPath + "raw/filtered/filteredGMC_p3.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	//storing_path = outputPath + "input_p3.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, storing_path.c_str(), LOAD_DATA, false);

	//const dataType sigma = 2.0;
	//copyDataToAnotherArray(ctContainer->dataPointer, toSmoothImageData, Height, Length, Width);
	//rescaleNewRange(toSmoothImageData, Length, Width, Height, 0.0, 1.0, maxValue, minValue);
	//gaussianFiltering(toSmoothImageData, imageData, Length, Width, Height,sigma);

	////Patient 1
	//Point3D seed1 = { 262, 258, 146 };
	//Point3D seed2 = { 266, 256, 245 };

	////Patient 2
	//Point3D seed1 = { 257, 254, 249 };
	//Point3D seed2 = { 257, 243, 350 };

	//Patient 3
	Point3D seed1 = { 268, 230, 116 };
	Point3D seed2 = { 266, 221, 218 };
	
	////Patient 4
	//Point3D seed1 = { 280, 229, 135 };
	//Point3D seed2 = { 285, 234, 220 };

	////Patient 5
	//Point3D seed1 = { 265, 243, 471 };
	//Point3D seed2 = { 239, 225, 625 };

	////Patient 6
	//Point3D seed1 = { 250, 298, 258 };
	//Point3D seed2 = { 268, 288, 443 };

	Point3D* endPoints = new Point3D[3];
	endPoints[0] = seed1;
	endPoints[1] = seed2;
	endPoints[2] = { 0.0, 0.0, 0.0 };
	double radius = 5.0;
	Potential_Parameters parameters{
		1000, //edge detector coefficient
		0.25, //threshold
		0.001,//epsilon
		radius
	};
	Image_Data inputImageStr = { Height, Length, Width, imageData, ctOrigin, ctSpacing, orientation };
	//compute3DPotential(inputImageStr, potential, endPoints, parameters);
	
	storing_path = outputPath + "potential_p3.raw";
	manageRAWFile3D<dataType>(potential, Length, Width, Height, storing_path.c_str(), LOAD_DATA, false);

	////Draw initial seed ball
	//string saving_ball = outputPath + "initial_ball_p1.csv";
	//FILE* initial_ball;
	//if (fopen_s(&initial_ball, saving_ball.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//BoundingBox3D box = findBoundingBox3D(seed1, Length, Width, Height, radius, 2);
	//Point3D centerBall = getRealCoordFromImageCoord3D(seed1, ctOrigin, ctSpacing, orientation);
	//for (size_t ik = box.k_min; ik <= box.k_max; ik++) {
	//	for (size_t ii = box.i_min; ii <= box.i_max; ii++) {
	//		for (size_t ij = box.j_min; ij <= box.j_max; ij++) {
	//			Point3D pointBall = { ii, ij, ik };
	//			pointBall = getRealCoordFromImageCoord3D(pointBall, ctOrigin, ctSpacing, orientation);
	//			double pDistance = getPoint3DDistance(pointBall, centerBall);
	//			if (pDistance <= radius) {
	//				fprintf(initial_ball, "%f,%f,%f\n", pointBall.x, pointBall.y, pointBall.z);
	//			}
	//		}
	//	}
	//}
	//fclose(initial_ball);

	Image_Data actionMapStr = { Height, Length, Width, action, ctOrigin, ctSpacing, orientation };
	storing_path = outputPath + "";
	partialFrontPropagation(actionMapStr, potential, endPoints);

	//////Save the end points in files
	//string saving_csv = outputPath + "endPoints_p3.csv";
	//FILE* f_end_points;
	//if (fopen_s(&f_end_points, saving_csv.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//Point3D seed_save = getRealCoordFromImageCoord3D(seed1, ctOrigin, ctSpacing, orientation);
	//fprintf(f_end_points, "%f,%f,%f\n", seed_save.x, seed_save.y, seed_save.z);
	//seed_save = getRealCoordFromImageCoord3D(seed2, ctOrigin, ctSpacing, orientation);
	//fprintf(f_end_points, "%f,%f,%f\n", seed_save.x, seed_save.y, seed_save.z);
	//fclose(f_end_points);

	//storing_path = outputPath + "action_map_p3.raw";
	//manageRAWFile3D<dataType>(action, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	////Extract and save the path points
	//dataType tau = 0.8, tolerance = 1.0;
	//Path_Parameters pathParameters = { tau, 1000, tolerance };
	//vector<Point3D> path_points;
	//shortestPath3D(actionMapStr, endPoints, path_points, pathParameters);
	//FILE* path_file;
	//storing_path = outputPath + "path_points_p3.csv";
	//if (fopen_s(&path_file, storing_path.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(path_file, "x,y,z\n");
	//for(int it = 0; it < path_points.size(); it++) {
	//	Point3D current_point = path_points[it];
	//	current_point = getRealCoordFromImageCoord3D(current_point, ctOrigin, ctSpacing, orientation);
	//	fprintf(path_file, "%f,%f,%f\n", current_point.x, current_point.y, current_point.z);
	//}
	//fclose(path_file);

	//loading_path = inputPath + "edge_image.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);
	////crop
	//size_t imin = 180, jmin = 170, kmin = 90;
	//const size_t length = 150, width = 150, height = 170;
	//dataType** imageCrop = new dataType * [height];
	//for(k = 0; k < height; k++) {
	//	imageCrop[k] = new dataType[length * width]{ 0 };
	//}

	//Point3D originCrop = { imin, jmin, kmin };
	//originCrop = getRealCoordFromImageCoord3D(originCrop, ctOrigin, ctSpacing, orientation);
	//std::cout << "CT origin : (" << originCrop.x << ", " << originCrop.y << ", " << originCrop.z << ")" << std::endl;

	//size_t i_new = 0, j_new = 0, k_new = 0;
	//for (k = 0, k_new = kmin; k < height; k++, k_new++) {
	//	for (i = 0, i_new = imin; i < length; i++, i_new++) {
	//		for (j = 0, j_new = jmin; j < width; j++, j_new++) {
	//			imageCrop[k][x_new(i, j, length)] = imageData[k_new][x_new(i_new, j_new, Length)];
	//		}
	//	}
	//}
	//storing_path = outputPath + "edge_image_crop_p3.raw";
	//manageRAWFile3D<dataType>(imageCrop, length, width, height, storing_path.c_str(), STORE_DATA, false);
	//for (k = 0; k < height; k++) {
	//	delete[] imageCrop[k];
	//}
	//delete[] imageCrop[k];

	delete[] endPoints;
	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
		delete[] action[k];
		delete[] potential[k];
		delete[] actionFirstFront[k];
		delete[] actionSecondFront[k];
	}
	delete[] imageData;
	delete[] action;
	delete[] potential;
	delete[] actionFirstFront;
	delete[] actionSecondFront;
	free(ctContainer);
	*/
	
	/*
    //Cropped
	const size_t Height = 170, Width = 150, Length = 150;
	const size_t dim2D = Length * Width;
	dataType** imageData = new dataType * [Height];
	dataType** action = new dataType * [Height];
	dataType** potential = new dataType * [Height];
	dataType** difference = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{ 0 };
		action[k] = new dataType[dim2D]{ 0 };
		potential[k] = new dataType[dim2D]{ 0 };
		difference[k] = new dataType[dim2D]{0};
	}
	loading_path = inputPath + "crop_p3.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);
	
	Point3D ctOrigin = { -74.2188, -83.9845, -651.5 };
	std::cout << "CT origin : (" << ctOrigin.x << ", " << ctOrigin.y << ", " << ctOrigin.z << ")" << std::endl;
	VoxelSpacing ctSpacing = { 0.976562, 0.976562, 2.5 };
	std::cout << "CT origin : (" << ctSpacing.sx << ", " << ctSpacing.sy << ", " << ctSpacing.sz << ")" << std::endl;

	Point3D seed1 = { 86, 60, 23 };
	Point3D seed2 = { 93, 57, 125 };

	Point3D* endPoints = new Point3D[2];
	endPoints[0] = seed1;
	endPoints[1] = seed2;
	double radius = 3.0;
	Potential_Parameters parameters{
		1000, //edge detector coefficient
		0.15, //threshold
		0.001,//epsilon
		radius
	};
	//OrientationMatrix orientation = { { 1.0, 0.0, 0.0 } , { 0.0, 1.0, 0.0 } , { 0.0, 0.0, 1.0 } };
	Image_Data inputImageStr = { Height, Length, Width, imageData, ctOrigin, ctSpacing, orientation };
	compute3DPotential(inputImageStr, potential, endPoints, parameters);

	//storing_path = outputPath + "crop/potential_crop_p3.raw";
	//manageRAWFile3D<dataType>(potential, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	const double LengthKeyPoints = 50;
	vector<Point3D> key_points;

	Image_Data actionMapStr = { Height, Length, Width, action, ctOrigin, ctSpacing, orientation };
	storing_path = outputPath + "Action Spiral/csv/action_";
	//partialFrontPropagation(actionMapStr, potential, endPoints, storing_path);
	////fastMarching3dWithSpacing(inputImageStr, action, potential, seed1);
	frontPropagationWithKeyPointDetection(actionMapStr, potential, endPoints, LengthKeyPoints, key_points, storing_path);

	////Extract and save the path points
	//dataType tau = 0.8, tolerance = 1.0;
	//Path_Parameters pathParameters = { tau, 1000, tolerance };
	//vector<Point3D> path_points;
	//shortestPath3D(actionMapStr, endPoints, path_points, pathParameters);
	//FILE* path_file;
	////storing_path = outputPath + "Action Spiral/path_points_partial.csv";
	//storing_path = outputPath + "Action Spiral/path_points_full.csv";
	//if (fopen_s(&path_file, storing_path.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(path_file, "x,y,z\n");
	//for(int it = 0; it < path_points.size(); it++) {
	//	Point3D current_point = path_points[it];
	//	current_point = getRealCoordFromImageCoord3D(current_point, ctOrigin, ctSpacing, orientation);
	//	fprintf(path_file, "%f,%f,%f\n", current_point.x, current_point.y, current_point.z);
	//}
	//fclose(path_file);

	//loading_path = outputPath + "distance_fm.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);
	//loading_path = outputPath + "distance_fs.raw";
	//manageRAWFile3D<dataType>(action, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);
	//loading_path = outputPath + "distance_rt.raw";
	//manageRAWFile3D<dataType>(potential, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);
	////Compute the difference between the two distance maps
	//dataType fs_norm = 0.0, rt_norm = 0.0;
	//for(k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		fs_norm += pow(imageData[k][i] - action[k][i], 2);
	//		rt_norm += pow(imageData[k][i] - potential[k][i], 2);
	//		difference[k][i] = imageData[k][i] - potential[k][i];
	//	}
	//}
	//fs_norm /= (dataType)(Height * dim2D);
	//rt_norm /= (dataType)(Height * dim2D);
	//std::cout << "Means Square FS " << fs_norm << std::endl;
	//std::cout << "Means Square RT " << rt_norm << std::endl;
	//storing_path = outputPath + "difference_fm_rt.raw";
	//manageRAWFile3D<dataType>(difference, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	delete[] endPoints;
	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
		delete[] action[k];
		delete[] potential[k];
		delete[] difference[k];
	}
	delete[] imageData;
	delete[] action;
	delete[] potential;
	delete[] difference;
	*/
	
	/*
	const size_t Length = 100, Width = 100;
	const size_t dim2D = Length * Width;
	dataType* imageData = new dataType[dim2D] { 0 };
	dataType* action = new dataType[dim2D] { 0 };
	dataType* potential = new dataType[dim2D] { 0 };
	dataType* dMapBruteForce = new dataType[dim2D]{ 0 };
	dataType* dMapRouyTourin = new dataType[dim2D]{ 0 };
	dataType* dMapFastSweeping = new dataType[dim2D]{ 0 };
	dataType* dMapFastMarching = new dataType[dim2D]{ 0 };
	dataType* differencedMap = new dataType[dim2D]{ 0 };

	Point2D center = { 50.0, 50.0 };
	dataType foregroundValue = 1.0;
	int scale = 0;
	//imageData[x_new(25, 25, Length)] = foregroundValue;

	////Square
	//for (i = 10; i <= 90; i++) {
	//	for (j = 10; j <= 90; j++) {
	//		if (i == 10 || i == 90 || j == 10 || j == 90) {
	//			imageData[x_new(i, j, Length)] = foregroundValue;
	//		}
	//	}
	//}

	//Circle
	double radius = 40.0, smallRadius = 20.0;
	double step = 2 * M_PI / 1000;
	double phi, coord_x, coord_y;
	size_t ind_x, ind_y;
	for(size_t it = 0; it < 1000; it++) {
		
		phi = it * step;
		coord_x = center.x + radius * cos(phi);
		coord_y = center.y + radius * sin(phi);
		if (coord_x < 0) {
			ind_x = 0;
		}
		else if (coord_x >= Length) {
			ind_x = Length - 1;
		}
		else {
			ind_x = (size_t)coord_x;
		}
		
		if (coord_y < 0) {
			ind_y = 0;
		}
		else if (coord_y >= Width) {
			ind_y = Width - 1;
		}
		else {
			ind_y = (size_t)(coord_y);
		}
		
		imageData[x_new(ind_x, ind_y, Length)] = 1.0;
		////Create holes in the circle
		//if (phi >= (M_PI / 6.0) && phi <= (M_PI / 3.0)) {
		//	imageData[x_new(ind_x, ind_y, Length)] = 0.0;
		//}
		//else {
		//	imageData[x_new(ind_x, ind_y, Length)] = 1.0;
		//}
	}
	//storing_path = outputPath + "input_circle.raw";
	//manageRAWFile2D<dataType>(imageData, Length, Width, storing_path.c_str(), STORE_DATA, false);

	////2D Spiral
	//dataType a = 0.8;
	//size_t nb_turns = 10;
	//for (dataType it = 0.001; it <= 2 * nb_turns * M_PI; it += 0.001) {
	//	size_t indx = (size_t)(a * it * cos(it) + 50);
	//	size_t indy = (size_t)(a * it * sin(it) + 50);
	//	if(indx < Length && indy < Width) {
	//		imageData[x_new(indx, indy, Length)] = foregroundValue;
	//	}
	//}
	//storing_path = outputPath + "input_spiral.raw";
	//manageRAWFile2D<dataType>(imageData, Length, Width, storing_path.c_str(), STORE_DATA, false);

	////Tube
	//for (i = 5; i < Length - 5; i++) 
	//{
	//	for (j = 5; j < Width - 5; j++) 
	//	{
	//		Point2D pPoint = { i, j };
	//		if (i == j) 
	//		{
	//			imageData[x_new(i - 5, j, Length)] = 1.0;
	//			imageData[x_new(i + 5, j, Length)] = 1.0;
	//			imageData[x_new(i, j - 5, Length)] = 1.0;
	//			imageData[x_new(i, j + 5, Length)] = 1.0;
	//			//for (size_t in = 0; in < Length; in++) 
	//			//{
	//			//	for (size_t jn = 0; jn < Width; jn++) 
	//			//	{
	//			//		Point2D current_point = { in, jn };
	//			//		double pDistance = getPoint2DDistance(current_point, pPoint);
	//			//		if (pDistance <= 10) 
	//			//		{
	//			//			imageData[x_new(in, jn, Length)] = 1.0;
	//			//		}
	//			//	}
	//			//}
	//		}
	//	}
	//}
	////storing_path = outputPath + "tube.raw";
	////manageRAWFile2D<dataType>(imageData, Length, Width, storing_path.c_str(), STORE_DATA, false);

	Point2D pOrigin = { 0.0, 0.0 };
	PixelSpacing pSpacing = { 1.0, 1.0 };
	OrientationMatrix2D pOrientation = { {1.0, 0.0}, {0.0, 1.0} };
	Image_Data2D toDistanceMap = { Length, Width, imageData, pOrigin, pSpacing, pOrientation };
	////fastMarchingForDistanceMap(toDistanceMap, dMapFastMarching, foregroundValue);
	bruteForceDistanceMap2D(toDistanceMap, dMapFastMarching, foregroundValue);

	rescaleNewRange2D(differencedMap, Length, Width, 0.0, 1.0);
	
	//for (i = 0; i < dim2D; i++) 
	//{
	//	dMapFastMarching[i] = 1.0 / (1.0 + 1.0 * dMapFastMarching[i]);
	//}

	//storing_path = outputPath + "f_distance_map.raw";
	//manageRAWFile2D<dataType>(dMapFastMarching, Length, Width, storing_path.c_str(), STORE_DATA, false);

	storing_path = outputPath + "f_distance_map.vtk";
	save2DImageAs3Dvtk(dMapFastMarching, Length, Width, (char*)storing_path.c_str(), scale);

	////storing_path = outputPath + "object_with_holes_four_pixels.raw";
	//storing_path = outputPath + "semi_circle_with_four_pixels_holes.raw";
	//manageRAWFile2D<dataType>(imageData, Length, Width, storing_path.c_str(), STORE_DATA, false);

	////Save the end points in files
	//FILE* p_end;
	//string end_points = outputPath + "end_points_tube.csv";
	////string end_points = outputPath + "end_u_tube.csv";
	//if (fopen_s(&p_end, end_points.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(p_end, "x,y\n");
	//fprintf(p_end, "%f,%f\n", endPoints[0].x, endPoints[0].y);
	//fprintf(p_end, "%f,%f\n", endPoints[1].x, endPoints[1].y);
	//fclose(p_end);

	Point2D seed1 = { 50.0, 50.0 };
	Point2D seed2 = { 80.0, 50.0 };
	Point2D* endPoints = new Point2D[2];
	endPoints[0] = seed1;
	endPoints[1] = seed2;

	double radius_potential = 3.0;
	Potential_Parameters parameters{
		1,//1000, //edge detector coefficient
		1,//0.15, //threshold
		0.01,//epsilon
		radius_potential
	};
	Image_Data2D imageDataStr = { Length, Width, imageData, pOrigin, pSpacing, pOrientation };
	//computePotential(imageDataStr, potential, endPoints, parameters);
	//storing_path = outputPath + "combined_potential_circle.raw";
	//manageRAWFile2D<dataType>(potential, Length, Width, storing_path.c_str(), STORE_DATA, false);

	
	//storing_path = outputPath + "combined_potential_circle.vtk";
	//save2DImageAs3Dvtk(potential, Length, Width, (char*)storing_path.c_str(), scale);
	
	//////storing_path = outputPath + "action empty tube/pot1_2/action_";
	//////storing_path = outputPath + "action empty tube with holes/pot2/action_";
	////storing_path = outputPath + "u shape/pot2/action_";
	//storing_path = outputPath + "action tube/pot1/full/action_";
	////partialFrontPropagation2D(imageDataStr, action, potential, endPoints, storing_path);

	//Path_Parameters parameters = { 0.8, 1000, 1.0 };
	//vector<Point2D> path_points;
	//Image_Data2D toActionStr = { Length, Width, action, origin, pSpacing, NULL };
	//shortestPath2d(toActionStr, endPoints, path_points, parameters);

	////Save the end points in files
	//FILE* p_points;
	//string file_points = outputPath + "path_points_empty_tube_pot1.csv";
	//string file_points = outputPath + "path_points_empty_tube.csv";
	//string file_points = outputPath + "u_tube_points_p2.csv";
	//if (fopen_s(&p_points, file_points.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(p_points, "x,y\n");
	//for(size_t it = 0; it < path_points.size(); it++) {
	//	fprintf(p_points, "%f,%f\n", path_points[it].x, path_points[it].y);
	//}
	//fclose(p_points);

	//delete[] endPoints;
	delete[] imageData;
	delete[] action;
	delete[] potential;

	delete[] dMapBruteForce;
	delete[] dMapRouyTourin;
	delete[] dMapFastSweeping;
	delete[] dMapFastMarching;
	delete[] differencedMap;
	*/
	
	/*
	//2D real images
	const size_t Length = 512;
	const size_t Width = 512;
	const size_t dim2D = Length * Width;
	dataType* imageData = new dataType[dim2D]{ 0 };
	dataType* smoothedImage = new dataType[dim2D]{ 0 };
	dataType* potential = new dataType[dim2D]{ 0 };
	dataType* action = new dataType[dim2D]{ 0 };
	*/

	//loading_path = inputPath + "raw/slice/slice_204_p1.raw";
	////loading_path = inputPath + "raw/slice/eye_image_512_512.raw";
	//manageRAWFile2D<dataType>(imageData, Length, Width, loading_path.c_str(), LOAD_DATA, false);
	//rescaleNewRange2D(imageData, Length, Width, 0.0, 1.0);
	
	/*
	storing_path = outputPath + "input_image_2D.raw";
	manageRAWFile2D<dataType>(imageData, Length, Width, storing_path.c_str(), LOAD_DATA, false);
	Image_Data2D imageDataStr = { Length, Width, imageData, {0.0, 0.0}, {1.0, 1.0}, {{1.0, 0.0}, {0.0, 1.0}} };

	//dataType min_val = imageData[0], max_val = imageData[0];
	//for(i = 0; i < dim2D; i++) {
	//	if(imageData[i] < min_val) {
	//		min_val = imageData[i];
	//	}
	//	if(imageData[i] > max_val) {
	//		max_val = imageData[i];
	//	}
	//}
	//std::cout << "Before Filtering Min: " << min_val << " Max: " << max_val << std::endl;
	
	Filter_Parameters filtering_parameters =
	{
		0.5,// timeStepSize;
		1.171875,// h not used here
		1.0,// sigma not used here
		0,// edge_detector_coefficient not used here
		1.4,// omega_c;
		0.001,// tolerance;
		0,// eps2 not used here
		0,// coef not used here
		1,// linked to sigma and not used here
		5,//number of time step;
		100// max number of iteration;
	};
	heatImplicit2dScheme(imageDataStr, filtering_parameters);

	//const dataType sigma = 3.0;
	//gaussianSmoothing2D(imageData, smoothedImage, Length, Width, sigma);
	//rescaleNewRange2D(smoothedImage, Length, Width, 0.0, 1.0);
	*/
	
	/*
	storing_path = outputPath + "filtered.raw";
	manageRAWFile2D<dataType>(imageData, Length, Width, storing_path.c_str(), LOAD_DATA, false);

	//min_val = smoothedImage[0], max_val = smoothedImage[0];
	//for (i = 0; i < dim2D; i++) {
	//	if (smoothedImage[i] < min_val) {
	//		min_val = smoothedImage[i];
	//	}
	//	if (smoothedImage[i] > max_val) {
	//		max_val = smoothedImage[i];
	//	}
	//}
	//std::cout << "After Filtering Min: " << min_val << " Max: " << max_val << std::endl;

	//int scale = 0;
	//storing_path = outputPath + "input_3D_view.vtk";
	//save2DImageAs3Dvtk(smoothedImage, Length, Width, (char*)storing_path.c_str(), scale);

	//End Points sllice 
	Point2D* endPoints = new Point2D[2];
	endPoints[0] = { 175.0, 310.0 };
	endPoints[1] = { 255.0, 208.0 };// One point is used for this test

	//endPoints[0] = { 161.0, 301.0 };
	//endPoints[1] = { 159.0, 237.0 };
	
	////Aorta End Points
	//endPoints[0] = { 240.0, 209.0 };
	//endPoints[1] = { 182.0, 340.0 };

	////End Points Distance Map
	//endPoints[0] = { 271.0, 186.0 };
	//endPoints[1] = { 181.0, 348.0 };

	////Segment liver slice
	//endPoints[0] = { 261.0, 238.0 };
	//endPoints[1] = { 194.0, 197.0 };

	////Retinal image End Points
	//endPoints[0] = { 188.0, 64.0 };
	//endPoints[1] = { 370.0, 320.0 };

	////Retinal image End Points second Tests
	//endPoints[1] = { 181.0, 50.0 };
	//endPoints[0] = { 355.0, 461.0 };
	////endPoints[1] = { 476.0, 420.0 };
	////endPoints[1] = { 470.0, 343.0 };
	//endPoints[2] = { 0.0, 0.0 };

	//endPoints[0] = { 170.0, 12.0 };
	//endPoints[1] = { 500.0, 433.0 };
	////endPoints[2] = { 0.0, 0.0 };

	////Save the end points in files
	//FILE* end_points_file;
	//string file_points = outputPath + "end_points.csv";
	//if (fopen_s(&end_points_file, file_points.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(end_points_file, "x,y\n");
	//for(size_t it = 0; it < 2; it++) {
	//	fprintf(end_points_file, "%f,%f\n", endPoints[it].x, endPoints[it].y);
	//}
	//fclose(end_points_file);

	OrientationMatrix2D orientation2D = { {1.0, 0.0}, {0.0, 1.0} };
	Point2D iOrigin = { 0.0, 0.0 };
	PixelSpacing spacing = { 1.0, 1.0 };
	double radius = 2.0;
	Potential_Parameters parameters{
		100, //edge detector coefficient
		0.4, //threshold
		0.1,//epsilon
		radius
	};
	Image_Data2D toPotentialStr = { Length, Width, imageData, iOrigin, spacing, orientation2D };
	computePotential(toPotentialStr, potential, endPoints, parameters);

	storing_path = outputPath + "potential.raw";
	//manageRAWFile2D<dataType>(potential, Length, Width, storing_path.c_str(), LOAD_DATA, false);

	Image_Data2D toActionStr = { Length, Width, imageData, iOrigin, spacing, orientation2D };
	//fastMarching2D(toActionStr, action, potential, endPoints);
	////partialFrontPropagation2D(toActionStr, action, potential, endPoints);
	////rouyTourinFrontPropagation2D(toActionStr, action, potential, 0.001, 5000);
	
	//clock_t start, end;
	////start = clock();
	////fastMarching2D(toActionStr, action, potential, endPoints);
	//partialFrontPropagation2D(toActionStr, action, potential, endPoints, storing_path);
	//end = clock();
	//double laps = (end - start);// / CLOCKS_PER_SEC;
	//std::cout << "execution time : " << CLOCKS_PER_SEC << std::endl;
	//start = clock();
	//rouyTourinFrontPropagation2D(toActionStr, action, potential, 0.01, 5000);
	//end = clock();
	//double laps = (end - start) / CLOCKS_PER_SEC;
	//std::cout << "execution time : " << laps << std::endl;

	////storing_path = outputPath + "actionFM.raw";
	//storing_path = outputPath + "actionRT.raw";
	//manageRAWFile2D<dataType>(action, Length, Width, storing_path.c_str(), STORE_DATA, false);

	//Path_Parameters parameters_path = { 0.8, 1000, 1.0 };
	//vector<Point2D> path_points;
	//Image_Data2D toExtractPath = { Length, Width, action, iOrigin, spacing, orientation2D };
	//shortestPath2d(toExtractPath, endPoints, path_points, parameters_path);

	//FILE* path_points_file;
	////string save_path_file = outputPath + "path_points_RT.csv";
	//string save_path_file = outputPath + "path_points_FM.csv";
	//if (fopen_s(&path_points_file, save_path_file.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(path_points_file, "x,y\n");
	//for(size_t it = 0; it < path_points.size(); it++) {
	//	fprintf(path_points_file, "%f,%f\n", path_points[it].x, path_points[it].y);
	//}
	//fclose(path_points_file);

	////Save a sequence of path points
	//size_t id_save = 0;
	//for (size_t it = 0; it < path_points.size(); it++) {
	//	if ((it > 0) && (it % 5 == 0)) {
	//		id_save++;
	//		string save_path_file = outputPath + "partial action/p_points/path_points_eye_" + to_string(id_save) + ".csv";
	//		FILE* path_points_file;
	//		if (fopen_s(&path_points_file, save_path_file.c_str(), "w") != 0) {
	//			printf("Enable to open");
	//			return false;
	//		}
	//		fprintf(path_points_file, "x,y\n");
	//		for (size_t ik = 0; ik < it; ik++) {
	//			fprintf(path_points_file, "%f,%f\n", path_points[ik].x, path_points[ik].y);
	//		}
	//		fclose(path_points_file);
	//	}
	//}

	////Threshold
	//dataType thres_min = 0.245, thres_max = 0.29;
	//for (i = 0; i < dim2D; i++) {
	//	if( imageData[i] >= thres_min && imageData[i] <= thres_max ) {
	//		imageData[i] = 1.0;
	//	}
	//	else {
	//		imageData[i] = 0.0;
	//	}
	//}
	//storing_path = outputPath + "threshold.raw";
	//manageRAWFile2D<dataType>(imageData, Length, Width, storing_path.c_str(), STORE_DATA, false);
	//Image_Data2D toDistanceMapStr = { Length, Width, imageData, iOrigin, spacing, orientation2D };
	//fastMarchingForDistanceMap(toDistanceMapStr, action, 0.0);
	////storing_path = outputPath + "distance_map.raw";
	////manageRAWFile2D<dataType>(action, Length, Width, storing_path.c_str(), STORE_DATA, false);
	//for(i = 0; i < dim2D; i++) {
	//	action[i] = 1.0 / (1.0 + 100 * action[i]);
	//}
	//storing_path = outputPath + "slice_aorta.vtk";
	//save2DImageAs3Dvtk(action, Length, Width, (char*)storing_path.c_str(), 2);

	//storing_path = outputPath + "potential_3D_view.vtk";
	//save2DImageAs3Dvtk(potential, Length, Width, (char*)storing_path.c_str(), scale);

	
	//Read the file of points with negative discriminant
	vector<Point2D> nd;
	dataType x, y;
	string path_discriminant = "C:/Users/Konan Allaly/Documents/Tests/output/negative_discrinant_simple_file.csv";
	FILE* pFile;
	if (fopen_s(&pFile, path_discriminant.c_str(), "r") != 0) {
		printf("Enable to open");
		return false;
	}
	while (feof(pFile) == 0) {
		fscanf_s(pFile, "%f", &x);
		fscanf_s(pFile, ",");
		fscanf_s(pFile, "%f", &y);
		fscanf_s(pFile, "\n");
		Point2D pt = { x, y };
		nd.push_back(pt);
	}
	//std::cout << "Number of points with negative discriminant: " << nd.size() << std::endl;
	fclose(pFile);

	dataType* actionFM = new dataType[dim2D]{ 0 };
	dataType* actionRT = new dataType[dim2D]{ 0 };
	storing_path = outputPath + "actionFM.raw";
	manageRAWFile2D<dataType>(actionFM, Length, Width, storing_path.c_str(), LOAD_DATA, false);
	storing_path = outputPath + "actionRT.raw";
	manageRAWFile2D<dataType>(actionRT, Length, Width, storing_path.c_str(), LOAD_DATA, false);
	dataType mean_diff = 0.0;
	for(i = 0; i < dim2D; i++)
	{
		mean_diff += pow(actionFM[i] - actionRT[i], 2);
	}
	mean_diff /= (dataType)dim2D;
	std::cout << "The mean diff is : " << mean_diff << std::endl;
	dataType diff_action = 0.0;
	dataType max_diff = 0.0, min_diff = 1000000.0;
	for (k = 0; k < nd.size(); k++) 
	{
		xd = x_new((size_t)nd[k].x, (size_t)nd[k].y, Length);
		diff_action = fabs(actionFM[xd] - actionRT[xd]);
		if (diff_action > max_diff) {
			max_diff = diff_action;
		}
		if (diff_action < min_diff) {
			min_diff = diff_action;
		}
		if(k < 10) {
			std::cout << "Point " << k << " : (" << nd[k].x << "," << nd[k].y << ") FM: " << actionFM[xd] << " RT: " << actionRT[xd] << " Diff: " << diff_action << std::endl;
		}
	}
	std::cout << "Max diff : " << max_diff << std::endl;
	std::cout << "Min diff : " << min_diff << std::endl;

	delete[] actionFM;
	delete[] actionRT;
	

	//vector<Point2D> fm, rt;
	//dataType x, y;
	//FILE* path_file_FM;
	//loading_path = outputPath + "fm_path_points.csv";
	//if (fopen_s(&path_file_FM, loading_path.c_str(), "r") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//while (feof(path_file_FM) == 0) {
	//	fscanf_s(path_file_FM, "%f", &x);
	//	fscanf_s(path_file_FM, ",");
	//	fscanf_s(path_file_FM, "%f", &y);
	//	fscanf_s(path_file_FM, "\n");
	//	Point2D pt = { x, y };
	//	fm.push_back(pt);
	//}
	//fclose(path_file_FM);
	//FILE* path_file_RT;
	//loading_path = outputPath + "rt_path_points.csv";
	//if (fopen_s(&path_file_RT, loading_path.c_str(), "r") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//while (feof(path_file_RT) == 0) {
	//	fscanf_s(path_file_RT, "%f", &x);
	//	fscanf_s(path_file_RT, ",");
	//	fscanf_s(path_file_RT, "%f", &y);
	//	fscanf_s(path_file_RT, "\n");
	//	Point2D pt = { x, y };
	//	rt.push_back(pt);
	//}
	//fclose(path_file_RT);
	//
	//dataType H_fm = 0.0;
	//for(i = 0; i < fm.size(); i++)
	//{
	//	dataType min_fm = 1000000;
	//	for(j = 0; j < rt.size(); j++)
	//	{
	//		dataType d = getPoint2DDistance(fm[i], rt[j]);
	//		if(min_fm > d)
	//		{
	//			min_fm = d;
	//		}
	//	}
	//	if (min_fm > H_fm) {
	//		H_fm = min_fm;
	//	}
	//}
	//
	//dataType H_rt = 0.0;
	//for (i = 0; i < rt.size(); i++)
	//{
	//	dataType min_rt = 1000000;
	//	for (j = 0; j < fm.size(); j++)
	//	{
	//		dataType d = getPoint2DDistance(rt[i], fm[j]);
	//		if (min_rt > d)
	//		{
	//			min_rt = d;
	//		}
	//	}
	//	if (min_rt > H_rt) {
	//		H_rt = min_rt;
	//	}
	//}
	//std::cout << "FM : " << H_fm / (dataType)fm.size() << std::endl;
	//std::cout << "RT : " << H_rt / (dataType)rt.size() << std::endl;

	delete[] endPoints;
	delete[] imageData;
	delete[] smoothedImage;
	delete[] potential;
	delete[] action;
	*/

	//==================== Test segmentation 2D =======================================================

	/*
	//const size_t Length = 512, Width = 512;
	//const size_t dim2D = Length * Width;
	dataType* imageData = new dataType[dim2D] {0};
	dataType* initialSegment = new dataType[dim2D]{ 0 };
	copyDataToAnother2dArray(ctContainer->dataPointer[186], imageData, Length, Width);

	storing_path = outputPath + "input.raw";
	manageRAWFile2D<dataType>(imageData, Length, Width, storing_path.c_str(), STORE_DATA, false);
	rescaleNewRange2D(imageData, Length, Width, 0.0, 1.0);

	Image_Data2D imageDataStr = { Length, Width, imageData, {0.0, 0.0}, {1.0, 1.0}, {{1.0, 0.0},{0.0, 1.0}} };
	dataType h = 1.0;
	const Filter_Parameters implicitParameters
	{
		0.2,// timeStepSize;
		h,// h;
		1.0,// sigma;
		1000,// edge detector coefficient;
		1.4,// omega_c;
		1e-3,// tolerance;
		1e-6,// eps2;
		1e-6,// coef;
		1,// p;
		3,// timeStepsNum;
		1000// maxNumberOfSolverIteration;
	};
	
	//heatImplicit2dScheme(imageDataStr, implicitParameters);
	//storing_path = outputPath + "filtered.raw";
	//manageRAWFile2D<dataType>(imageData, Length, Width, storing_path.c_str(), STORE_DATA, false);

	Point2D* center = new Point2D[1];
	center[0] = {165.0, 279.0};
	dataType v = 0.5, R = 30.0;
	generateInitialSegmentationFunction(initialSegment, Length, Width, center, v, R);
	//storing_path = outputPath + "seg00.raw";
	//manageRAWFile2D<dataType>(initialSegment, Length, Width, storing_path.c_str(), STORE_DATA, false);

	Segmentation_Parameters segmentation_parms
	{
		50, // Maximum number of Gauss-Seidel iterations
		10000, // constant K in the Perona-Malik function G for the image
		1e-6, // epsilon is the regularization factor (Evans-Spruck)
		2000,// Number of current time step
		2000,// Maximum number of time step
		10, // Kind of writing density
		1e-6, // Tolerance for stopping of the segmentation process
		5.0 * h, //tau
		h, //h
		1.4, //omega_c
		1e-6, // gauss seidelTolerance;
		1.0,// coef_conv;// controle advection
		0.5//coef_dif; // controle curvature
	};
	//string segmentPath = outputPath + "seg/segment_";
	//subsurf(imageDataStr, initialSegment, segmentPath.c_str(), implicitParameters, segmentation_parms);
	
	string segmentPath = outputPath + "seg/IIOE/segment_";
	gsubsurf_iioe(imageDataStr, initialSegment, segmentPath.c_str(), implicitParameters, segmentation_parms);
	//string segmentPath = outputPath + "seg/SUBSURF/segment_";
	//subsurf(imageDataStr, initialSegment, (const char*)segmentPath.c_str(), implicitParameters, segmentation_parms);
	//string segmentPath = outputPath + "seg/GSUBSURF/segment_";
	//gsubsurf(imageDataStr, initialSegment, (const char*)segmentPath.c_str(), implicitParameters, segmentation_parms);

	delete[] center;
	delete[] imageData;
	delete[] initialSegment;
	
	free(ctContainer);
	*/

	//==================== Test segmentation regtangular grid ==================================

	/*
	dataType** inputImageData = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		inputImageData[k] = new dataType[dim2D]{ 0 };
	}

	//copyDataToAnotherArray(ctContainer->dataPointer, inputImageData, Height, Length, Width);

	//Croping : p1
	size_t i_min = 200;
	size_t j_min = 200;
	size_t k_min = 120;
	size_t length = 150;
	size_t width = 150;
	size_t height = 180;
	dataType** imageData = new dataType * [height];
	dataType** initialSegment = new dataType * [height]{ 0 };
	for (k = 0; k < height; k++) {
		imageData[k] = new dataType[length * width]{ 0 };
		initialSegment[k] = new dataType[length * width]{0};
	}
	
	//size_t i_ext, j_ext, k_ext;
	//for(k = 0, k_ext = k_min; k < height; k++, k_ext++)
	//{
	//	for (i = 0, i_ext = i_min; i < length; i++, i_ext++) 
	//	{
	//		for (j = 0, j_ext = j_min; j < width; j++, j_ext++) 
	//		{
	//			imageData[k][x_new(i, j, length)] = inputImageData[k_ext][x_new(i_ext, j_ext, Length)];
	//		}
	//	}
	//}
	////Find min and max
	//dataType minValue = 1000000.0, maxValue = -1000000.0;
	//for (k = 0; k < height; k++) {
	//	for (i = 0; i < length; i++) {
	//		for (j = 0; j < width; j++) {
	//			dataType value = imageData[k][x_new(i, j, length)];
	//			if (value < minValue) minValue = value;
	//			if (value > maxValue) maxValue = value;
	//		}
	//	}
	//}
	//std::cout << "Min new : " << minValue << ", Max new: " << maxValue << std::endl;
	//rescaleNewRange(imageData, length, width, height, 0.0, 1.0, maxValue, minValue);
	
	storing_path = outputPath + "cropped_p1.raw";
	//manageRAWFile3D<dataType>(imageData, length, width, height, storing_path.c_str(), STORE_DATA, false);

	Point3D newOrigin = {i_min, j_min, k_min};
	newOrigin = getRealCoordFromImageCoord3D(newOrigin, ctOrigin, ctSpacing, orientation);
	std::cout << "New origin : " << newOrigin.x << ", " << newOrigin.y << ", " << newOrigin.z << std::endl;
	Image_Data inputImage = { height, length, width, imageData, newOrigin, ctSpacing, orientation };

	////Generate the initial segmentation
	//FILE* path_file;
	//loading_path = inputPath + "paths/path segmentation/centered_path_p1.csv";
	//if (fopen_s(&path_file, loading_path.c_str(), "r") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//dataType x = 0, y = 0, z = 0;
	//vector<Point3D> path_points;
	//while (feof(path_file) == 0) {
	//	fscanf_s(path_file, "%f", &x);
	//	fscanf_s(path_file, ",");
	//	fscanf_s(path_file, "%f", &y);
	//	fscanf_s(path_file, ",");
	//	fscanf_s(path_file, "%f", &z);
	//	fscanf_s(path_file, "\n");
	//	Point3D current_point = { x, y, z };
	//	//current_point = getImageCoordFromRealCoord3D(current_point, newOrigin, ctSpacing, orientation);
	//	path_points.push_back(current_point);
	//	//maskSegment[(size_t)current_point.z][x_new((size_t)current_point.x, (size_t)current_point.y, length)] = 1.0;
	//}
	//fclose(path_file);
	////storing_path = outputPath + "mask_segment_p1.raw";
	////manageRAWFile3D<dataType>(maskSegment, length, width, height, storing_path.c_str(), STORE_DATA, false);
	//size_t n;
	//Point3D pPoints, pCurrent;
	//double dist, min_dist;
	//for (k = 0; k < height; k++)
	//{
	//	for (i = 0; i < length; i++)
	//	{
	//		for (j = 0; j < width; j++)
	//		{
	//			pCurrent = { (dataType)i, (dataType)j, (dataType)k };
	//			pCurrent = getRealCoordFromImageCoord3D(pCurrent, newOrigin, ctSpacing, orientation);
	//			min_dist = (double)(length * width);
	//			for (n = 0; n < path_points.size(); n++)
	//			{
	//				pPoints = path_points[n];
	//				dist = getPoint3DDistance(pPoints, pCurrent);
	//				if (min_dist > dist)
	//				{
	//					min_dist = dist;
	//				}
	//			}
	//			initialSegment[k][x_new(i, j, length)] = 1.0 / (1.0 + min_dist);
	//		}
	//	}
	//}
	
	storing_path = outputPath + "initial_segment_p1.raw";
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

	////geodesicMeanCurvature(inputImage, smoothParameters);
	storing_path = outputPath + "filtered_GMCF_p1.raw";
	manageRAWFile3D<dataType>(imageData, length, width, height, storing_path.c_str(), LOAD_DATA, false);

	Segmentation_Parameters segParameters = 
	{
		30,//Maximum number of Gauss-Seidel iterations
		100000,//edge detector coef
		1e-6,//epsilon is the regularization factor (Evans-Spruck)
		200,//Number of current time step
		200,//Maximum number of time step
		10,//saving frequency
		1e-6,//segmentation tolerance
		1.0,//tau
		1,//h
		1.4,//omega_c
		0.001,//tolerance
		1.0,//convection coef
		0.05,//diffusion coef
	};

	//storing_path = outputPath + "segmentation rectangular/gsubsurf/";
	////GSUBSURF_IIOE(inputImage, initialSegment, storing_path.c_str(), smoothParameters, segParameters);
	//GSUBSURF(inputImage, initialSegment, storing_path.c_str(), smoothParameters, segParameters);
	
	//storing_path = outputPath + "segmentation rectangular/gsubsurf_iioe/";
	//GSUBSURF_IIOE(inputImage, initialSegment, storing_path.c_str(), smoothParameters, segParameters);

	storing_path = outputPath + "segmentation rectangular/gsubsurf_s_one_iioe/";
	GSUBSURF_S_ONE_IIOE(inputImage, initialSegment, storing_path.c_str(), smoothParameters, segParameters);

	for(k = 0; k < Height; k++)
	{
		if(k < height)
		{
			delete[] imageData[k];
			delete[] initialSegment[k];
		}
		delete[] inputImageData[k];
	}
	delete[] imageData;
	delete[] initialSegment;
	delete[] inputImageData;

	free(ctContainer);
	*/

	//==================== Path Extraction 3D image ==================================

	//3D real image
	
	dataType** imageData = new dataType * [Height];
	dataType** potential = new dataType * [Height];
	dataType** action = new dataType * [Height];
	for(k = 0; k < Height; k++)
	{
		imageData[k] = new dataType[dim2D]{ 0 };
		potential[k] = new dataType[dim2D]{ 0 };
		action[k] = new dataType[dim2D]{ 0 };
	}

	//dataType minValue = 1000000.0, maxValue = -1000000.0;
	//for (k = 0; k < Height; k++) 
	//{
	//	for (i = 0; i < dim2D; i++) 
	//	{
	//		imageData[k][i] = ctContainer->dataPointer[k][i];
	//		if (imageData[k][i] < minValue) minValue = imageData[k][i];
	//		if (imageData[k][i] > maxValue) maxValue = imageData[k][i];
	//	}
	//}
	//std::cout << "Min data : " << minValue << ", Max data: " << maxValue << std::endl;
	//rescaleNewRange(imageData, Length, Width, Height, 0.0, 1.0, maxValue, minValue);
	//storing_path = outputPath + "rescaled_p5.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	//Filter_Parameters smoothParameters =
	//{
	//	0.2,// tau;
	//	1.0,// h not used here
	//	1.0,// sigma
	//	1000,// edge_detector_coefficient
	//	1.4,// omega_c;
	//	1e-3,// tolerance;
	//	1e-6,// eps2
	//	1,// coef
	//	1,// linked to sigma
	//	1,//number of time step;
	//	100// max number of iteration;
	//};
	Image_Data inputImage = { Height, Length, Width, imageData, ctOrigin, ctSpacing, orientation };
	//geodesicMeanCurvature(inputImage, smoothParameters);
	//storing_path = outputPath + "P3/filtered_p3.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, storing_path.c_str(), LOAD_DATA, false);

	////epsilon setting
	//for(k = 0; k < Height; k++)
	//{
	//	for(i = 0; i < dim2D; i++)
	//	{
	//		if(imageData[k][i] < 1.0e-6)
	//		{
	//			imageData[k][i] = 0.0;
	//		}
	//	}
	//}

	//Point3D* endPoints = new Point3D[2];
	//endPoints[0] = { 261.0, 257.0, 145.0 };
	//endPoints[1] = { 259.0, 250.0, 246.0 };

	//Point3D* endPoints = new Point3D[2];//p2
	//endPoints[0] = { 260.0, 254.0, 246.0 };
	//endPoints[1] = { 255.0, 245.0, 350.0 };

	//Point3D* endPoints = new Point3D[2];//p3
	//endPoints[0] = { 268.0, 231.0, 112.0 };
	//endPoints[1] = { 254.0, 224.0, 221.0 };

	//Point3D* endPoints = new Point3D[2];//p4
	//endPoints[0] = { 279.0, 229.0, 134.0 };
	//endPoints[1] = { 280.0, 235.0, 223.0 };

	Point3D* endPoints = new Point3D[2];//p5
	endPoints[0] = { 265.0, 244.0, 470.0 };
	endPoints[1] = { 235.0, 219.0, 626.0 };

	//Point3D* endPoints = new Point3D[2];//p6
	//endPoints[0] = { 249.0, 299.0, 256.0 };
	//endPoints[1] = { 258.0, 286.0, 443.0 };

	double radius = 3.0;
	Potential_Parameters parameters{
		1000, //edge detector coefficient
		0.25, //threshold, (0.15 --> p1, p2)
		0.001,//epsilon
		radius
	};
	//compute3DPotential(inputImage, potential, endPoints, parameters);

	storing_path = outputPath + "P5/potential_p5.raw";
	manageRAWFile3D<dataType>(potential, Length, Width, Height, storing_path.c_str(), LOAD_DATA, false);

	Image_Data toAction = { Height, Length, Width, action, ctOrigin, ctSpacing, orientation };
	//////frontPropagation(inputImage, action, potential, endPoints[0]);
	//partialFrontPropagation(toAction, potential, endPoints);

	storing_path = outputPath + "action_map_p5.raw";
	//manageRAWFile3D<dataType>(action, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	vector<Point3D> key_points;
	const double LengthKeyPoints = 50.0;
	frontPropagationWithKeyPointDetection(toAction, potential, endPoints, LengthKeyPoints, key_points);

	storing_path = outputPath + "keyp_action_map_p5.raw";
	//manageRAWFile3D<dataType>(action, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	//FILE* key_points_file;
	//string save_key_file = outputPath + "key_points_p4.csv";
	//if (fopen_s(&key_points_file, save_key_file.c_str(), "w") != 0) 
	//{
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(key_points_file, "x,y,z\n");
	//for (size_t it = 0; it < key_points.size(); it++)
	//{
	//	//convert to real world coordinates
	//	Point3D kp = getRealCoordFromImageCoord3D(key_points[it], ctOrigin, ctSpacing, orientation);
	//	fprintf(key_points_file, "%f,%f,%f\n", kp.x, kp.y, kp.z);
	//}
	//fclose(key_points_file);

	Path_Parameters parameters_path
	{
		0.8, // tau
		1000,// max number of iterations
		0.8 // tolerance
	};
	vector<Point3D> path_points;
	//path_points.push_back(endPoints[1]);
	//path_points.push_back(endPoints[0]);
	Image_Data toPathExtraction = { Height, Length, Width, action, ctOrigin, ctSpacing, orientation };
	//shortestPath3D(toPathExtraction, endPoints, path_points, parameters_path);

	FILE* path_points_file;
	//string save_path_file = outputPath + "end_points_p5.csv";
	string save_path_file = outputPath + "path_points_new_p5.csv";
	//string save_path_file = outputPath + "path_points_p5.csv";
	if (fopen_s(&path_points_file, save_path_file.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	fprintf(path_points_file, "x,y,z\n");

	//for(size_t it = 0; it < path_points.size(); it++) 
	//{
	//	//convert to real world coordinates
	//	path_points[it] = getRealCoordFromImageCoord3D(path_points[it], ctOrigin, ctSpacing, orientation);
	//	fprintf(path_points_file, "%f,%f,%f\n", path_points[it].x, path_points[it].y, path_points[it].z);
	//}
	//fclose(path_points_file);

	storing_path = outputPath + "action_p5_";
	string new_storing_path;
	for(int i_n = key_points.size() - 1; i_n > 0; i_n--)
	{
		endPoints[0] = key_points[i_n];
		endPoints[1] = key_points[i_n - 1];
		partialFrontPropagation(toAction, potential, endPoints);
		shortestPath3D(toPathExtraction, endPoints, path_points, parameters_path);
		extension = to_string(i_n);
		new_storing_path = storing_path + extension + ".raw";
		manageRAWFile3D<dataType>(action, Length, Width, Height, new_storing_path.c_str(), STORE_DATA, false);
		for(size_t it = 0; it < path_points.size(); it++) 
		{
			//convert to real world coordinates
			path_points[it] = getRealCoordFromImageCoord3D(path_points[it], ctOrigin, ctSpacing, orientation);
			fprintf(path_points_file, "%f,%f,%f\n", path_points[it].x, path_points[it].y, path_points[it].z);
		}
		path_points.clear();
	}
	fclose(path_points_file);
	
	delete[] endPoints;	
	for (k = 0; k < Height; k++) 
	{
		delete[] imageData[k];
		delete[] potential[k];
		delete[] action[k];
	}
	delete[] imageData;
	delete[] potential;
	delete[] action;

	free(ctContainer);
	

	/*
	//3D artificial image
	
	////Artificial image
	//const size_t Height = 100, Width = 60, Length = 60;//short spiral
	const size_t Height = 100, Width = 100, Length = 100;
	const size_t dim2D = Length * Width;
	dataType** imageData = new dataType * [Height];
	dataType** newImageData = new dataType * [Height];
	dataType** action = new dataType * [Height];
	dataType** potential = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{ 0 };
		newImageData[k] = new dataType[dim2D]{ 0 };
		action[k] = new dataType[dim2D]{ 0 };
		potential[k] = new dataType[dim2D]{ 0 };
	}
	
	////loading_path = inputPath + "spiral.raw";
	////loading_path = inputPath + "low_contrast_spiral.raw";
	//loading_path = inputPath + "shape/spiral/spiral.raw";
	////loading_path = inputPath + "shape/Action/empty_spiral.raw";
	////loading_path = inputPath + "shape/Action/empty_spiral_with_hole.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	////create holes
	//Point3D hole_center = { 27, 31, 50 };
	//double hole_radius = 10.0;
	//for(k = 0; k < Height; k++)
	//{
	//	for(i = 0; i < Length; i++)
	//	{
	//		for(j = 0; j < Width; j++)
	//		{
	//			Point3D p = { (dataType)i, (dataType)j, (dataType)k };
	//			double dist = getPoint3DDistance(p, hole_center);
	//			if(dist <= hole_radius)
	//			{
	//				imageData[k][x_new(i, j, Length)] = 0.0;
	//			}
	//		}
	//	}
	//}
	
	//storing_path = outputPath + "with_hole.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, storing_path.c_str(), LOAD_DATA, false);
	//////decrease the contrast
	
	//dataType new_value;
	//for (k = 0; k < Height; k++) 
	//{
	//	for(i = 0; i < dim2D; i++) 
	//	{
	//		if(imageData[k][i] > 0.0) 
	//		{
	//			new_value = 0.8;//generateRandNormal(0.8, 0.01);
	//		}
	//		else 
	//		{
	//			new_value = 0.7;//generateRandNormal(0.7, 0.01);
	//		}
	//		imageData[k][i] = new_value;
	//		
	//		//if (imageData[k][i] == 1.0) 
	//		//{
	//		//	imageData[k][i] = 0.8;
	//		//}
	//		//else {
	//		//	imageData[k][i] = 0.7;
	//		//}
	//	}
	//}
	//storing_path = outputPath + "noisy_image.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);
	
	Point3D iOrigin = { 0.0, 0.0, 0.0 };
	VoxelSpacing iSpacing = { 1.0, 1.0, 1.0 };
	OrientationMatrix orientation = { { 1.0, 0.0, 0.0 } , { 0.0, 1.0, 0.0 } , { 0.0, 0.0, 1.0 } };
	Image_Data inputImageStr = { Height, Length, Width, imageData, iOrigin, iSpacing, orientation };

	//Filter_Parameters smoothParameters =
	//{
	//	1.0,// tau;
	//	1.0,// h not used here
	//	1.0,// sigma
	//	1000,// edge_detector_coefficient
	//	1.4,// omega_c;
	//	1e-3,// tolerance;
	//	1e-6,// eps2
	//	1,// coef
	//	1,// linked to sigma
	//	1,//number of time step;
	//	100// max number of iteration;
	//};
	//heatImplicitRectangularScheme(inputImageStr, smoothParameters);
	storing_path = outputPath + "edge_image.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, storing_path.c_str(), LOAD_DATA, false);

	//End Points for the spiral
	Point3D seed1 = { 81, 52, 14 };
	Point3D seed2 = { 79, 47, 84 };
	Point3D* endPoints = new Point3D[2];
	endPoints[0] = seed1;
	endPoints[1] = seed2;
	vector<Point3D> key_points;
	const double LengthKeyPoints = 15;

	////create hole close to point 2
	Point3D hole_center1 = { 79, 47, 70 };
	//Point3D hole_center1 = { 49, 20, 60 };
	Point3D hole_center2 = { 29, 47, 55 };
	for(k = 0; k < Height; k++)
	{
		for(i = 0; i < Length; i++)
		{
			for(j = 0; j < Width; j++)
			{
				Point3D p = { (dataType)i, (dataType)j, (dataType)k };
				double dist = getPoint3DDistance(p, hole_center1);
				if(dist <= 8.0)
				{
					imageData[k][x_new(i, j, Length)] = 0.0;
				}
				dist = getPoint3DDistance(p, hole_center2);
				if (dist <= 8.0)
				{
					imageData[k][x_new(i, j, Length)] = 0.0;
				}
			}
		}
	}
	storing_path = outputPath + "with_hole.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	////End Points for the spiral: short spiral
	//Point3D seed1 = { 40, 30, 85 };
	//Point3D seed2 = { 42, 30, 15 };
	//Point3D* endPoints = new Point3D[2];
	//endPoints[0] = seed1;
	//endPoints[1] = seed2;
	//vector<Point3D> key_points;
	//const double LengthKeyPoints = 15;

	////End Points for the spiral
	//Point3D seed1 = { 90, 53, 15 };
	//Point3D seed2 = { 90, 45, 85 };
	//Point3D* endPoints = new Point3D[2];
	//endPoints[0] = seed1;
	//endPoints[1] = seed2;
	//vector<Point3D> key_points;
	//const double LengthKeyPoints = 15;

	////Empty the long spiral
	//Point3D p_grad;
	//for (k = 0; k < Height; k++)
	//{
	//	for (i = 0; i < Length; i++)
	//	{
	//		for (j = 0; j < Width; j++)
	//		{
	//			if (imageData[k][x_new(i, j, Length)] > 0.0)
	//			{
	//				p_grad = { (dataType)i, (dataType)j, (dataType)k };
	//				getGradient3D(inputImageStr, i, j, k, &p_grad);
	//				double norm = sqrt(p_grad.x * p_grad.x + p_grad.y * p_grad.y + p_grad.z * p_grad.z);
	//				if (norm > 0.0) 
	//				{
	//					newImageData[k][x_new(i, j, Length)] = 1.0;
	//				}
	//			}
	//		}
	//	}
	//}

	//Image_Data ImageStr = { Height, Length, Width, newImageData, iOrigin, iSpacing, orientation };
	//storing_path = "C:/Users/Konan Allaly/Documents/Tests/output/empty_spiral.raw";
	//manageRAWFile3D<dataType>(newImageData, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	////End points spiral thin
	//Point3D seed1 = { 90, 51, 15 };
	//Point3D seed2 = { 90, 47, 84 };
	//Point3D* endPoints = new Point3D[2];
	//endPoints[0] = seed1;
	//endPoints[1] = seed2;
	//vector<Point3D> key_points;
	//const double LengthKeyPoints = 15;

	double radius_potential = 3.0;
	Potential_Parameters parameters{
		1,//1000, //edge detector coefficient
		1,//0.15, //threshold
		0.01,//epsilon
		radius_potential
	};
	Image_Data imageStr = { Height, Length, Width, imageData, iOrigin, iSpacing, orientation };
	//compute3DPotential(imageStr, potential, endPoints, parameters);

	//storing_path = outputPath + "potential_v1.raw";
	//manageRAWFile3D<dataType>(potential, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	//Image_Data toDistanceMap = { Height, Length, Width, imageData, iOrigin, iSpacing, orientation };
	//fastMarching3dForDistanceMap(toDistanceMap, action, 1.0);
	//storing_path = outputPath + "distance_map.raw";
	//manageRAWFile3D<dataType>(action, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	//storing_path = outputPath + "potential.raw";
	//manageRAWFile3D<dataType>(potential, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	Image_Data actionMapStr = { Height, Length, Width, action, iOrigin, iSpacing, orientation };
	//partialFrontPropagation(actionMapStr, potential, endPoints);

	storing_path = outputPath + "action.raw";
	//manageRAWFile3D<dataType>(action, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	////key_points.clear();
	//key_points.push_back(seed1);
	//key_points.push_back(seed2);
	////Save the keys points in files
	//string saving_csv = outputPath + "endpoints_spiral.csv";
	//FILE* f_key_point;
	//if (fopen_s(&f_key_point, saving_csv.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(f_key_point, "x,y,z\n");
	//for (int n = 0; n < key_points.size(); n++) {
	//	fprintf(f_key_point, "%f,%f,%f\n", key_points[n].x, key_points[n].y, key_points[n].z);
	//}
	//fclose(f_key_point);
	//key_points.clear();

	//Extract and save the path points
	dataType tau = 0.95, tolerance = 1.0;
	Path_Parameters pathParameters = { tau, 1000, tolerance };
	vector<Point3D> path_points;
	//path_points.push_back(endPoints[1]);
	//path_points.push_back(endPoints[0]);
	//shortestPath3D(actionMapStr, endPoints, path_points, pathParameters);
	FILE* path_file;
	//storing_path = outputPath + "endPoints.csv";
	//storing_path = outputPath + "path_points_classic.csv";
	//storing_path = outputPath + "path_points_distance_map.csv";
	//storing_path = outputPath + "path_points_shortcut.csv";
	//storing_path = outputPath + "path_points_shortcut_V2.csv";
	//if (fopen_s(&path_file, storing_path.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(path_file, "x,y,z\n");
	//for(int it = 0; it < path_points.size(); it++) {
	//	Point3D current_point = path_points[it];
	//	fprintf(path_file, "%f,%f,%f\n", current_point.x, current_point.y, current_point.z);
	//}
	//fclose(path_file);
	

	////Compare the distance for different methods
	//dataType** distanceMapBruteForce = new dataType * [Height];
	//dataType** distanceMapFastMarching = new dataType * [Height];
	//dataType** distanceMapFastSweeping = new dataType * [Height];
	//dataType** distanceMapRouyTourin = new dataType * [Height];
	//dataType** compare = new dataType * [Height];
	//for (k = 0; k < Height; k++) {
	//	distanceMapBruteForce[k] = new dataType[dim2D]{ 0 };
	//	distanceMapFastMarching[k] = new dataType[dim2D]{ 0 };
	//	distanceMapFastSweeping[k] = new dataType[dim2D]{ 0 };
	//	distanceMapRouyTourin[k] = new dataType[dim2D]{ 0 };
	//	compare[k] = new dataType[dim2D]{ 0 };
	//}
	//imageData[50][x_new(30, 30, Length)] = 1.0;

	//////bruteForceFunction_3D(distanceMapBruteForce, imageData, Length, Width, Height, 1000000000000.0, 0.0);
	////bruteForceDistanceMap(inputImageStr, distanceMapBruteForce, 1.0);
	//loading_path = outputPath + "distance_map_brtforce.raw";
	//manageRAWFile3D<dataType>(distanceMapBruteForce, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	//////fastMarching(distanceMapFastMarching, imageData, Height, Length, Width, 1.0);
	////fastMarching3dForDistanceMap(inputImageStr, distanceMapFastMarching, 1.0);
	//loading_path = outputPath + "distance_map_FM.raw";
	//manageRAWFile3D<dataType>(distanceMapFastMarching, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	//fastSweepingFunction_3D(distanceMapFastSweeping, imageData, Length, Width, Height, 1.0, 10000000000.0, 0.0);
	//fastSweepingDistanceMap(inputImageStr, distanceMapFastSweeping, 1.0);
	//loading_path = outputPath + "distance_map_fstswp.raw";
	//manageRAWFile3D<dataType>(compare, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	//rouyTourinFunction_3D(distanceMapRouyTourin, imageData, 0.5, Length, Width, Height, 0.4, 1.0);
	//rouyTourinDistanceMap(inputImageStr, distanceMapRouyTourin, 1.0, 0.4, 0.5);
	//loading_path = outputPath + "distance_map_RT.raw";
	//manageRAWFile3D<dataType>(distanceMapRouyTourin, Length, Width, Height, loading_path.c_str(), STORE_DATA, false);

	//dataType norm_fm = 0.0, norm_fstswp = 0.0, norm_rt = 0, norm_compare = 0.0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		norm_fm += pow(distanceMapBruteForce[k][i] - distanceMapFastMarching[k][i], 2);
	//		//norm_fstswp += pow(distanceMapBruteForce[k][i] - distanceMapFastSweeping[k][i], 2);
	//		//norm_rt += pow(distanceMapBruteForce[k][i] - distanceMapRouyTourin[k][i], 2);
	//		//norm_compare += pow(distanceMapRouyTourin[k][i] - compare[k][i], 2);

	//		compare[k][i] = distanceMapFastMarching[k][i] - distanceMapBruteForce[k][i];
	//	}
	//}
	////std::cout << "Residual Fast Marching: " << sqrt(norm_fm) << std::endl;
	////std::cout << "Residual Fast Sweeping: " << sqrt(norm_fstswp) << std::endl;
	////std::cout << "Residual Rouy Tourin: " << sqrt(norm_rt) << std::endl;
	//////std::cout << "Residual Compare: " << sqrt(norm_compare) << std::endl;

	//std::cout << "Diagonal lenght: " << sqrt( Length * Length + Width * Width ) << std::endl;
	//loading_path = outputPath + "difference.raw";
	//manageRAWFile3D<dataType>(compare, Length, Width, Height, loading_path.c_str(), STORE_DATA, false);

	//for (k = 0; k < Height; k++) {
	//	delete[] distanceMapBruteForce[k];
	//	delete[] distanceMapFastMarching[k];
	//	delete[] distanceMapFastSweeping[k];
	//	delete[] distanceMapRouyTourin[k];
	//}
	//delete[] distanceMapBruteForce;
	//delete[] distanceMapFastMarching;
	//delete[] distanceMapFastSweeping;
	//delete[] distanceMapRouyTourin;

	delete[] endPoints;
	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
		delete[] action[k];
		delete[] potential[k];
		delete[] newImageData[k];
	}
	delete[] imageData;
	delete[] action;
	delete[] potential;
	delete[] newImageData;
	*/

	//==================== Compare distance map 3D ==================================

	/*
	dataType** distanceMapFM = new dataType * [Height];
	dataType** distanceMapFS = new dataType * [Height];
	for(k = 0; k < Height; k++)
	{
		distanceMapFM[k] = new dataType[dim2D]{ 0 };
		distanceMapFS[k] = new dataType[dim2D]{ 0 };
	}

	storing_path = outputPath + "distance_map_fm.raw";
	manageRAWFile3D<dataType>(distanceMapFM, Length, Width, Height, storing_path.c_str(), LOAD_DATA, false);

	storing_path = outputPath + "distance_map_fs.raw";
	manageRAWFile3D<dataType>(distanceMapFS, Length, Width, Height, storing_path.c_str(), LOAD_DATA, false);

	dataType norm = 0.0, diff = 0.0;
	for(k = 0; k < Height; k++)
	{
		for(i = 0; i < dim2D; i++)
		{
			//norm += pow(distanceMapFM[k][i] - distanceMapFS[k][i], 2);
			diff = fabs(distanceMapFM[k][i] - distanceMapFS[k][i]);
			if(norm < diff)
			{
				norm = diff;
			}
		}
	}
	//std::cout << "The mean square diff is : " << norm / (dataType)(Height * dim2D) << std::endl;
	std::cout << "Max difference " << norm << std::endl;

	for(k = 0; k < Height; k++)
	{
		delete[] distanceMapFM[k];
		delete[] distanceMapFS[k];
	}
	delete[] distanceMapFM;
	delete[] distanceMapFS;
	free(ctContainer);
	*/

	return EXIT_SUCCESS;
}