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
	loading_path = root + "input/vtk/petct/ct/Patient1_ct.vtk";
	readVtkFile(loading_path.c_str(), ctContainer);

	std::cout << "============ Input CT ================ " << std::endl;

	int Height = ctContainer->dimensions[2];
	int Length = ctContainer->dimensions[0];
	int Width = ctContainer->dimensions[1];
	int dim2D = Length * Width;
	std::cout << "CT image dim : " << ctContainer->dimensions[0] << " x " << ctContainer->dimensions[1] << " x " << ctContainer->dimensions[2] << "" << std::endl;

	std::cout << "CT origin : (" << ctContainer->origin[0] << ", " << ctContainer->origin[1] << ", " << ctContainer->origin[2] << ")" << std::endl;
	Point3D imageOrigin = { ctContainer->origin[0], ctContainer->origin[1], ctContainer->origin[2] };
	VoxelSpacing imageSpacing = { ctContainer->spacing[0], ctContainer->spacing[1], ctContainer->spacing[2] };
	std::cout << "CT spacing : (" << ctContainer->spacing[0] << ", " << ctContainer->spacing[1] << ", " << ctContainer->spacing[2] << ")" << std::endl; 
	std::cout << "=========================================" << std::endl;

	/*
	float** imageDataF = new float* [Height];
	double** imageDataD = new double* [Height];
	for (k = 0; k < Height; k++) {
		imageDataF[k] = new float[dim2D] {0};
		imageDataD[k] = new double[dim2D] {0};
	}

	const char * load_ptr = "C:/Users/Konan Allaly/Documents/Tests/input/raw/filtered/filtered_p1.raw";
	manageRAWFile3D<float>(imageDataF, Length, Width, Height, load_ptr, LOAD_DATA, false);

	for (k = 0; k < Height; k++) 
	{
		for (i = 0; i < dim2D; i++) 
		{
			imageDataD[k][i] = (double)imageDataF[k][i];
		}
	}
	const char* store_ptr = "C:/Users/Konan Allaly/Documents/Tests/input/raw/filtered/filtered_p1_double.raw";
	manageRAWFile3D<double>(imageDataD, Length, Width, Height, store_ptr, STORE_DATA, false);

	for(k = 0; k < Height; k++)
	{
		delete[] imageDataF[k];
		delete[] imageDataD[k];
	}
	delete[] imageDataF;
	delete[] imageDataD;
	free(ctContainer);
	*/
	
	//==================== Translate for registration ================================================

	/*
	Point3D origin_before = { -300, -230, -1022.5 };
	Point3D origin_after = { -300, -230, -1339.5 };
	VoxelSpacing ctSpacing = { 1.171875, 1.171875, 2.5 };

	Point3D t1 = { 260, 258, 142 };
	t1 = getRealCoordFromImageCoord3D(t1, origin_before, ctSpacing, orientation);
	
	Point3D t2 = { 260, 255, 243 };
	t2 = getRealCoordFromImageCoord3D(t2, origin_after, ctSpacing, orientation);

	Point3D translation = { t1.x - t2.x, t1.y - t2.y, t1.z - t2.z };

	t2.x = t2.x + translation.x;
	t2.y = t2.y + translation.y;
	t2.z = t2.z + translation.z;

	//t1 = getRealCoordFromImageCoord3D(t1, origin_before, ctSpacing, orientation);
	//t2 = getRealCoordFromImageCoord3D(t2, origin_after, ctSpacing, orientation);

	//Point3D newOrigin = { origin_after.x + translation.x, origin_after.y + translation.y , origin_after.z + translation.z };
	//newOrigin = getRealCoordFromImageCoord3D(newOrigin, origin_after, ctSpacing, orientation);
	//std::cout << "New origin : (" << newOrigin.x << ", " << newOrigin.y << ", " << newOrigin.z << ")" << std::endl;

	Point3D origin_after_pet = { -286.586, -216.586, -1338.5 };
	Point3D newOriginPET = { origin_after_pet.x + translation.x, origin_after_pet.y + translation.y , origin_after_pet.z + translation.z };
	std::cout << "New origin : (" << newOriginPET.x << ", " << newOriginPET.y << ", " << newOriginPET.z << ")" << std::endl;

	////Save points
	//string saving_csv = outputPath + "points.csv";
	//FILE* f_point;
	//if (fopen_s(&f_point, saving_csv.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(f_point, "x,y,z\n");
	//fprintf(f_point, "%f,%f,%f\n", t1.x, t1.y, t1.z);
	//fprintf(f_point, "%f,%f,%f\n", t2.x, t2.y, t2.z);
	//fclose(f_point);

	free(ctContainer);
	*/
	
	//==================== Compute Hausdoff distance and Ratio =======================================
	
	/*
	dataType img_f, n0, n1, n2, nmg, x, y, z, ptmg, ptid, scal;
	char header[MAX_LINE_LENGTH];
	string inputPointCloud = inputPath + "vtk/petct/aorta/Hausdoff 31-03/Isolines/";

	//FILE* file_hausdoff;
	//storing_path = inputPointCloud + "h distance/hausdoff_distance_patient2.csv";
	//if (fopen_s(&file_hausdoff, storing_path.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(file_hausdoff, "manual,gsubsurf,HD,ratio,mean_manual,mean_gsubsurf,MHD,ratio_mean\n");

	
	FILE* file_manual;
	loading_path = inputPointCloud + "manual/_patient_6.csv";
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
	std::string path_root = inputPointCloud + "gsubsurf/_03_patient_6.csv";
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

	//fprintf(file_hausdoff, "%f,%f,%f,%f,%f,%f,%f,%f\n", max_manual, max_gsubsurf, HD, ratio, mean_manual, mean_gsubsurf, MHD, ratio_mean);
	//fclose(file_hausdoff);

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
		//fprintf(file_hausdoff, "%f,%f,%f,%f,%f,%f,%f,%f\n", max_manual, max_gsubsurf, HD, ratio, mean_manual, mean_gsubsurf, MHD, ratio_mean);

		while (points_gsubsurf.size() > 0)
		{
			points_gsubsurf.pop_back();
		}

	}
	//fclose(file_hausdoff);
	*/

	//==================== Histogram on Segmentation Mask ============================================

	/*
	//const size_t length = 150, width = 150, height = 350;//p1,p2
	//const size_t length = 150, width = 150, height = 380;//p3
	//const size_t length = 150, width = 150, height = 350;//p4
	const size_t length = 160, width = 160, height = 360;//p5
	//const size_t length = 150, width = 150, height = 400;//p6
	loading_path = outputPath + "Segmentation/Aorta/p5/_seg_func_05000.raw";

	//Point3D segmentOrigin = { -63.9648, -236.996, 1354.21 };

	dataType** segFunc = new dataType * [height];
	for (k = 0; k < height; k++) 
	{
		segFunc[k] = new dataType[width * length]{ 0 };
	}
	manageRAWFile3D<dataType>(segFunc, length, width, height, loading_path.c_str(), LOAD_DATA, false);

	//find max and min value
	dataType maxValue = segFunc[0][0];
	dataType minValue = segFunc[0][0];
	for (k = 0; k < height; k++)
	{
		for (i = 0; i < length; i++)
		{
			for (j = 0; j < width; j++)
			{
				xd = x_new(i, j, length);
				if (segFunc[k][xd] > maxValue) {
					maxValue = segFunc[k][xd];
				}
				if (segFunc[k][xd] < minValue) {
					minValue = segFunc[k][xd];
				}
			}
		}
	}
	std::cout << "Max value = " << maxValue << std::endl;
	std::cout << "Min value = " << minValue << std::endl;
	//rescaleNewRange(segFunc, length, width, height, 0.0, 256.0, maxValue, minValue);

	const size_t binCount = 100;
	size_t* histogram = new size_t[binCount]{ 0 };

	dataType sizeClass = (maxValue - minValue) / (dataType)binCount;

	computeHistogram(segFunc, histogram, length, width, height, binCount);

	dataType peak1 = 0, peak2 = 0;
	size_t index_peak1 = 0, index_peak2 = 0;
	for(i = 0; i < binCount; i++)
	{
		if(histogram[i] > peak1)
		{
			peak1 = histogram[i];
			index_peak1 = i;
		}
		else if(histogram[i] > peak2 && histogram[i] != histogram[index_peak1])
		{
			peak2 = histogram[i];
			index_peak2 = i;
		}
	}
	std::cout << "Peak 1 = " << peak1 << ", index = " << index_peak1 << std::endl;
	std::cout << "Peak 2 = " << peak2 << ", index = " << index_peak2 << std::endl;

	FILE* file_histogram;
	storing_path = outputPath + "histogram_seg_patient5.csv";
	if (fopen_s(&file_histogram, storing_path.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}

	dataType isovalue = 0.0;
	//fprintf(file_histogram, "%f,%d\n", isovalue, 0);

	fprintf(file_histogram, "isosurface,count_voxels\n");
	for(i = 0; i < binCount; i++)
	{
		isovalue = (i + 1) * sizeClass;
		//fprintf(file_histogram, "%f,%d\n", isovalue, histogram[i]);
		if(isovalue > 0.09)
		{
			fprintf(file_histogram, "%f,%d\n", isovalue, histogram[i]);
		}
	}
	fclose(file_histogram);

	delete[] histogram;
	for(k = 0; k < height; k++)
	{
		delete[] segFunc[k];
	}
	delete[] segFunc;
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
	
	////2D real images
	//const size_t Length = 512;
	//const size_t Width = 512;
	//const size_t dim2D = Length * Width;
	//dataType* imageData = new dataType[dim2D]{ 0 };
	//dataType* smoothedImage = new dataType[dim2D]{ 0 };
	//dataType* potential = new dataType[dim2D]{ 0 };
	//dataType* action = new dataType[dim2D]{ 0 };

	//loading_path = inputPath + "raw/slice/slice_204_p1.raw";
	////loading_path = inputPath + "raw/slice/eye_image_512_512.raw";
	//manageRAWFile2D<dataType>(imageData, Length, Width, loading_path.c_str(), LOAD_DATA, false);
	//rescaleNewRange2D(imageData, Length, Width, 0.0, 1.0);
	
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
	const size_t Length = 512, Width = 512;
	const size_t dim2D = Length * Width;
	dataType* imageData = new dataType[dim2D] {0};
	dataType* initialSegment = new dataType[dim2D]{ 0 };
	//copyDataToAnother2dArray(ctContainer->dataPointer[206], imageData, Length, Width);
	//rescaleNewRange2D(imageData, Length, Width, 0.0, 1.0);

	storing_path = outputPath + "input_slice.raw";
	//manageRAWFile2D<dataType>(imageData, Length, Width, storing_path.c_str(), LOAD_DATA, false);
	
	Point2D sliceOrigin = { 0.0, 0.0 };
	PixelSpacing sliceSpacing = { 1.171875, 1.171875 };
	Image_Data2D imageDataStr = { Length, Width, imageData, sliceOrigin, sliceSpacing, {{1.0, 0.0},{0.0, 1.0}} };
	
	const Filter_Parameters filteringParameters
	{
		1.5,// timeStepSize;
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
	
	////heatImplicit2dScheme(imageDataStr, filteringParameters);
	//geodesicMeanCurvature2D(imageDataStr, filteringParameters);
	storing_path = outputPath + "filtered.raw";
	manageRAWFile2D<dataType>(imageData, Length, Width, storing_path.c_str(), LOAD_DATA, false);

	Point2D* center = new Point2D[1];
	center[0] = { 191, 275 };// { 165.0, 279.0 };
	dataType v = 0.5, R = 200.0;
	generateInitialSegmentationFunction(initialSegment, Length, Width, center, v, R);
	storing_path = outputPath + "Segmentation 2D/seg00.raw";
	manageRAWFile2D<dataType>(initialSegment, Length, Width, storing_path.c_str(), STORE_DATA, false);

	Segmentation_Parameters segmentation_parms
	{
		50, // Maximum number of Gauss-Seidel iterations
		20000, // constant K in the Perona-Malik function G for the image
		1e-6, // epsilon is the regularization factor (Evans-Spruck)
		2000,// Number of current time step
		2000,// Maximum number of time step
		10, // Kind of writing density
		1e-6, // Tolerance for stopping of the segmentation process
		sliceSpacing.sx, //tau
		1.0, //h
		1.4, //omega_c
		1e-6, // gauss seidelTolerance;
		5.0 * sliceSpacing.sx,// coef_conv;// controle advection
		2.0 * sliceSpacing.sx//coef_dif; // controle curvature
	};
	
	//string segmentPath = outputPath + "seg/SUBSURF/_";
	//subsurf(imageDataStr, initialSegment, segmentPath.c_str(), filteringParameters, segmentation_parms);
	
	//string segmentPath = outputPath + "Segmentation 2D/implicit/_";
	//gsubsurf_implicit(imageDataStr, initialSegment, (const char*)segmentPath.c_str(), filteringParameters, segmentation_parms);
	
	//gsubsurf_iioe(imageDataStr, initialSegment, (const char*)segmentPath.c_str(), filteringParameters, segmentation_parms);
	
	string segmentPath = outputPath + "Segmentation 2D/s1iioe/_";
	gsubsurf_s_one_iioe(imageDataStr, initialSegment, (const char*)segmentPath.c_str(), filteringParameters, segmentation_parms);
	
	//string segmentPath = outputPath + "Segmentation 2D/s2iioe/_";
	//gsubsurf_s_two_iioe(imageDataStr, initialSegment, (const char*)segmentPath.c_str(), filteringParameters, segmentation_parms);
	
	//string segmentPath = outputPath + "seg/New folder/_";
	//gsubsurf(imageDataStr, initialSegment, (const char*)segmentPath.c_str(), filteringParameters, segmentation_parms);
	

	delete[] center;
	delete[] imageData;
	delete[] initialSegment;
	
	//free(ctContainer);
	*/

	//==================== Test segmentation 3D =======================================================

	/*
	dataType** inputImageData = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		inputImageData[k] = new dataType[dim2D]{ 0 };
	}

	//copyDataToAnotherArray(ctContainer->dataPointer, inputImageData, Height, Length, Width);

	////Croping : p1
	//size_t i_min = 200;
	//size_t j_min = 200;
	//size_t k_min = 120;
	//size_t length = 150;
	//size_t width = 150;
	//size_t height = 180;

	////Croping : p3
	//size_t i_min = 200;
	//size_t j_min = 185;
	//size_t k_min = 90;
	//size_t length = 150;
	//size_t width = 150;
	//size_t height = 170;

	//Croping : p5
	size_t i_min = 170;
	size_t j_min = 170;
	size_t k_min = 435;
	size_t length = 170;
	size_t width = 170;
	size_t height = 260;

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
	
	//storing_path = outputPath + "Segmentation/cropped_p5.raw";
	//manageRAWFile3D<dataType>(imageData, length, width, height, storing_path.c_str(), STORE_DATA, false);

	Point3D newOrigin = {i_min, j_min, k_min};
	newOrigin = getRealCoordFromImageCoord3D(newOrigin, imageOrigin, imageSpacing, orientation);
	std::cout << "New origin : " << newOrigin.x << ", " << newOrigin.y << ", " << newOrigin.z << std::endl;
	Image_Data inputImage = { height, length, width, imageData, newOrigin, imageSpacing, orientation };

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
	//geodesicMeanCurvature(inputImage, smoothParameters);

	storing_path = outputPath + "Segmentation/filtered_GMCF_p5.raw";
	manageRAWFile3D<dataType>(imageData, length, width, height, storing_path.c_str(), LOAD_DATA, false);

	////Generate the initial segmentation
	//FILE* path_file;
	////loading_path = inputPath + "paths/path segmentation/centered_path_p3.csv";
	//loading_path = outputPath + "Segmentation/Aorta/centered paths/finalCurve_p5.csv";
	//if (fopen_s(&path_file, loading_path.c_str(), "r") != 0) 
	//{
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
	//	path_points.push_back(current_point);
	//}
	//fclose(path_file);

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
	//			pCurrent = getRealCoordFromImageCoord3D(pCurrent, newOrigin, imageSpacing, orientation);
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

	storing_path = outputPath + "Segmentation/initial_segment_p5.raw";
	manageRAWFile3D<dataType>(initialSegment, length, width, height, storing_path.c_str(), LOAD_DATA, false);

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
		0.001,//diffusion coef
	};

	storing_path = outputPath + "Segmentation/P5/";
	Image_Data imageToSegment = { height, length, width, imageData, newOrigin, imageSpacing, orientation };
	//generalizedSubsurfSegmentation(imageToSegment, initialSegment, segParameters, smoothParameters, (unsigned char*)storing_path.c_str());
	GSUBSURF(imageToSegment, initialSegment, storing_path.c_str(), smoothParameters, segParameters);
	//GSUBSURF_S_ONE_IIOE(inputImage, initialSegment, storing_path.c_str(), smoothParameters, segParameters);
	////GSUBSURF(inputImage, initialSegment, storing_path.c_str(), smoothParameters, segParameters);

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

	//==================== Path Extraction 3D image ===================================================

	/*
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

	//Point3D* endPoints = new Point3D[2];//p1
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
	//geodesicMeanCurvature(inputImage, smoothParameters);
	storing_path = outputPath + "Data journal paper submission/P5/filtered_p5.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, storing_path.c_str(), LOAD_DATA, false);

	double radius = 3.0;
	Potential_Parameters parameters{
		1000, //edge detector coefficient
		0.3,  //threshold, (0.15 --> p1, p2), threshold (0.2 --> p3, p4, p5, p6)
		0.001,//epsilon
		radius
	};
	//compute3DPotential(inputImage, potential, endPoints, parameters);

	//storing_path = outputPath + "Data journal paper submission/potential_p6.raw";
	storing_path = outputPath + "potential_p5.raw";
	manageRAWFile3D<dataType>(potential, Length, Width, Height, storing_path.c_str(), LOAD_DATA, false);

	Image_Data toAction = { Height, Length, Width, action, imageOrigin, imageSpacing, orientation };
	partialFrontPropagation(toAction, potential, endPoints);
	//storing_path = outputPath + "Data journal paper submission/action_map_partial_p6.raw";
	//manageRAWFile3D<dataType>(action, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	vector<Point3D> key_points;
	const double LengthKeyPoints = 50.0;
	frontPropagationWithKeyPointDetection(toAction, potential, endPoints, LengthKeyPoints, key_points);
	FILE* key_points_file;
	string save_key_file = outputPath + "Data journal paper submission/key_points_p5.csv";
	if (fopen_s(&key_points_file, save_key_file.c_str(), "w") != 0) 
	{
		printf("Enable to open");
		return false;
	}
	fprintf(key_points_file, "x,y,z\n");
	for (size_t it = 0; it < key_points.size(); it++)
	{
		//convert to real world coordinates
		Point3D kp = getRealCoordFromImageCoord3D(key_points[it], imageOrigin, imageSpacing, orientation);
		fprintf(key_points_file, "%f,%f,%f\n", kp.x, kp.y, kp.z);
	}
	fclose(key_points_file);
	
	Path_Parameters parameters_path
	{
		0.8, // tau
		1000,// max number of iterations
		0.8 // tolerance
	};

	vector<Point3D> path_points;
	//path_points.push_back(endPoints[1]);
	//path_points.push_back(endPoints[0]);
	Image_Data toPathExtraction = { Height, Length, Width, action, imageOrigin, imageSpacing, orientation };
	//shortestPath3D(toPathExtraction, endPoints, path_points, parameters_path);

	//FILE* path_points_file;
	////string save_path_file = outputPath + "Data journal paper submission/end_points_p6.csv";
	//string save_path_file = outputPath + "path_points_p6.csv";
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
	
	//string new_storing_path;
	//vector<Point3D> path_points;
	//Image_Data toPathExtraction = { Height, Length, Width, action, imageOrigin, imageSpacing, orientation };
	//FILE* path_points_file;
	////string save_path_file = outputPath + "Data journal paper submission/end_points_p1.csv";
	//string save_path_file = outputPath + "Data journal paper submission/path_points_kp_p6.csv";
	//if (fopen_s(&path_points_file, save_path_file.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(path_points_file, "x,y,z\n");
	//storing_path = outputPath + "Data journal paper submission/key_p_action_";
	//for(int i_n = key_points.size() - 1; i_n > 0; i_n--)
	//{
	//	endPoints[0] = key_points[i_n];
	//	endPoints[1] = key_points[i_n - 1];
	//	partialFrontPropagation(toAction, potential, endPoints);
	//	shortestPath3D(toPathExtraction, endPoints, path_points, parameters_path);
	//	extension = to_string(i_n);
	//	new_storing_path = storing_path + extension + ".raw";
	//	manageRAWFile3D<dataType>(action, Length, Width, Height, new_storing_path.c_str(), STORE_DATA, false);
	//	for(int it = path_points.size() - 1; it > -1; it--) 
	//	{
	//		//convert to real world coordinates
	//		path_points[it] = getRealCoordFromImageCoord3D(path_points[it], imageOrigin, imageSpacing, orientation);
	//		fprintf(path_points_file, "%f,%f,%f\n", path_points[it].x, path_points[it].y, path_points[it].z);
	//	}
	//	path_points.clear();
	//}
	//fclose(path_points_file);
	
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
	*/
	
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

	//==================== Quantitative analysis study 3D =============================================

	/*
    //Translated origin
	//CT = {-300, -226.484, -1275}
	Point3D imageOrigin = { -300, -226.484, -1275 }; //for analysis after treatment
	//PET = {-286.586, -213.07, -1274}
	
    
	//==================== PET ================
	Vtk_File_Info* petContainer = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	petContainer->operation = copyFrom;
	loading_path = inputPath + "vtk/petct/pet/Patient2_pet.vtk";
	readVtkFile(loading_path.c_str(), petContainer);
	//Point3D PETimageOrigin = { petContainer->origin[0], petContainer->origin[1], petContainer->origin[2] };
	Point3D PETimageOrigin = { -286.586, -213.07, -1274 };// for analysis after treatment
	std::cout << "PET origin : (" << petContainer->origin[0] << ", " << petContainer->origin[1] << ", " << petContainer->origin[2] << ")" << std::endl;
	VoxelSpacing PETimageSpacing = { petContainer->spacing[0], petContainer->spacing[1], petContainer->spacing[2] };
	std::cout << "PET spacing : (" << petContainer->spacing[0] << ", " << petContainer->spacing[1] << ", " << petContainer->spacing[2] << ")" << std::endl;
	size_t height_pet = petContainer->dimensions[2];
	size_t width_pet = petContainer->dimensions[1];
	size_t length_pet = petContainer->dimensions[0];
	std::cout << "PET image dim : " << petContainer->dimensions[0] << " x " << petContainer->dimensions[1] << " x " << petContainer->dimensions[2] << "" << std::endl;
	std::cout << "=========================================" << std::endl;
	//==========================================
	
	// Load liver CT
	dataType** maskLiver = new dataType * [Height];
	dataType** imageDataFull = new dataType * [Height];
	dataType** maskAorta = new dataType * [Height];
	dataType** distanceMap = new dataType * [Height]{ 0 };
	for (k = 0; k < Height; k++)
	{
		maskLiver[k] = new dataType[dim2D]{ 0 };
		imageDataFull[k] = new dataType[dim2D]{ 0 };
		maskAorta[k] = new dataType[dim2D] { 0 };
		distanceMap[k] = new dataType[dim2D]{ 0 };
	}
	Image_Data inputImageStr = { Height, Length, Width, imageDataFull, imageOrigin, imageSpacing, orientation };
	loading_path = inputPath + "raw/liver/liver_p2.raw";
	manageRAWFile3D<dataType>(maskLiver, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	//Ball Liver PET
	dataType** ballLiverPet = new dataType * [height_pet];
	bool** statusPetBall = new bool * [height_pet];//ensure that each voxel is used once to compute the mean value
	for(k = 0; k < height_pet; k++)
	{
		ballLiverPet[k] = new dataType[length_pet * width_pet]{ 0 };
		statusPetBall[k] = new bool[length_pet * width_pet]{ false };
	}

	//Find Liver centroid CT
	dataType* centroid = new dataType[3];
	centroidImage(maskLiver, centroid, Height, Length, Width, 0.0);
	Point3D center_liver_img_coord_ct = { centroid[0], centroid[1], centroid[2] };
	Point3D center_liver_real_coord = getRealCoordFromImageCoord3D(center_liver_img_coord_ct, imageOrigin, imageSpacing, orientation);
	Point3D center_liver_img_coord_pet = getImageCoordFromRealCoord3D(center_liver_real_coord, PETimageOrigin, PETimageSpacing, orientation);
	std::cout << "Liver centroid : (" << center_liver_real_coord.x << ", " << center_liver_real_coord.y << ", " << center_liver_real_coord.z << ")" << std::endl;
	//double radius_liver_roi = 20;// approx 15mm
	double radius_liver_roi = 20 + 1.5 * fmax(PETimageSpacing.sx, fmax(PETimageSpacing.sy, PETimageSpacing.sz));// approx 15mm

	////get mean SUV liver
	//dataType mean_SUV_Liver = 2.39835;//p1
	//dataType mean_SUV_Liver = 2.37516;
	//dataType mean_SUV_Liver = 2.57092;//p2

	dataType mean_SUV_Liver = 0.0;
	
	size_t count_voxels = 0;
	BoundingBox3D box = findBoundingBox3D(center_liver_img_coord_pet, length_pet, width_pet, height_pet, radius_liver_roi, 2.0);
	for (size_t kk = box.k_min; kk <= box.k_max; kk++)
	{
		for (size_t ii = box.i_min; ii <= box.i_max; ii++)
		{
			for (size_t jj = box.j_min; jj <= box.j_max; jj++)
			{
				
				Point3D p_pet = { (dataType)ii, (dataType)jj, (dataType)kk };
				Point3D p = getRealCoordFromImageCoord3D(p_pet, PETimageOrigin, PETimageSpacing, orientation);
				double dist = getPoint3DDistance(p, center_liver_real_coord);
				if (dist <= radius_liver_roi)
				{
					size_t ii_pet = (size_t)p_pet.x;
					size_t jj_pet = (size_t)p_pet.y;
					size_t kk_pet = (size_t)p_pet.z;
					if(statusPetBall[kk_pet][x_new(ii_pet, jj_pet, length_pet)] == false)
					{
						ballLiverPet[kk_pet][x_new(ii_pet, jj_pet, length_pet)] = 1.0;
						mean_SUV_Liver += petContainer->dataPointer[kk_pet][x_new(ii_pet, jj_pet, length_pet)];
						count_voxels++;
						statusPetBall[kk_pet][x_new(ii_pet, jj_pet, length_pet)] == true;
					}
				}
			}
		}
	}
	
	mean_SUV_Liver /= (dataType)count_voxels;
	std::cout << "Liver mean SUV : " << mean_SUV_Liver << std::endl;
	//storing_path = outputPath + "P2 translated/ball_liver_p1.raw";
	//manageRAWFile3D<dataType>(ballLiverPet, length_pet, width_pet, height_pet, storing_path.c_str(), STORE_DATA, false);

	for(k = 0; k < height_pet; k++)
	{
		delete[] ballLiverPet[k];
		delete[] statusPetBall[k];
	}
	delete[] ballLiverPet;
	delete[] statusPetBall;

	////==================== Segmented aorta distance map ==================================
	//const size_t length = 150, width = 150, height = 350;//p1,p2
	////const size_t length = 150, width = 150, height = 380;//p3
	////const size_t length = 150, width = 150, height = 350;//p4
	////const size_t length = 160, width = 160, height = 360;//p5
	////const size_t length = 150, width = 150, height = 400;//p6
	//dataType** imageData = new dataType * [height];
	//dataType** maskData = new dataType * [height];
	//bool** statusPoints = new bool * [height];
	//for (k = 0; k < height; k++)
	//{
	//	imageData[k] = new dataType[length * width]{ 0 };
	//	maskData[k] = new dataType[length * width]{ 0 };
	//	statusPoints[k] = new bool[length * width]{ 0 };
	//}

	//////Patient 1
	//Point3D segmentOrigin = { -65.625, 4.375, -700.234 };
	//VoxelSpacing segmentSpacing = { 1.171875, 1.171875, 1.171875 };
	//Image_Data inputImageStr = { height, length, width, imageData, segmentOrigin, segmentSpacing, orientation };

	////Patient 2
	//Point3D segmentOrigin = { -77.3438, 16.0938, -780.316 };
	//VoxelSpacing segmentSpacing = { 1.171875, 1.171875, 1.171875 };
	//Image_Data inputImageStr = { height, length, width, imageData, segmentOrigin, segmentSpacing, orientation };

	////Patient 3
	//Point3D segmentOrigin = { -64.4532, -74.2188, -622.594 };
	//VoxelSpacing segmentSpacing = { 0.976562, 0.976562, 0.976562 };
	//Image_Data inputImageStr = { height, length, width, imageData, segmentOrigin, segmentSpacing, orientation };

	////Patient 4
	//Point3D segmentOrigin = { -64.4532, -64.4532, -537.266 };
	//VoxelSpacing segmentSpacing = { 0.976562, 0.976562, 0.976562 };
	//Image_Data inputImageStr = { height, length, width, imageData, segmentOrigin, segmentSpacing, orientation };

	////Patient 5
	//Point3D segmentOrigin = { -63.9648, -236.996, 1354.21 };
	//VoxelSpacing segmentSpacing = { 0.976562, 0.976562, 0.976562 };
	//Image_Data inputImageStr = { height, length, width, imageData, segmentOrigin, segmentSpacing, orientation };

	////Patient 6
	//Point3D segmentOrigin = { -73.7305, -206.754, -767.055 };
	//VoxelSpacing segmentSpacing = { 0.976562, 0.976562, 0.976562 };
	//Image_Data inputImageStr = { height, length, width, imageData, segmentOrigin, segmentSpacing, orientation };

	//loading_path = outputPath + "Segmentation/3D image/p1/_seg_func_05000.raw";
	//manageRAWFile3D<dataType>(imageData, length, width, height, loading_path.c_str(), LOAD_DATA, false);
	//for (k = 0; k < height; k++) 
	//{
	//	for(i = 0; i < length * width; i++) 
	//	{
	//		if(imageData[k][i] >= 0.2) 
	//		{
	//			imageData[k][i] = 1.0;
	//		}
	//		else 
	//		{
	//			imageData[k][i] = 0.0;
	//		}
	//	}
	//}

	//storing_path = outputPath + "P1/seg_aorta_isoline_02.raw";
	//manageRAWFile3D<dataType>(imageData, length, width, height, storing_path.c_str(), LOAD_DATA, false);
	////Observe mask at aorta edge in PET
	//size_t min_i, max_i, min_j, max_j, min_k, max_k;
	//size_t l_mask = 0;
	//for(k = 0; k < height; k++)
	//{
	//	for(i = 0; i < length; i++)
	//	{
	//		for(j = 0; j < width; j++)
	//		{
	//			if(imageData[k][x_new(i, j, length)] == 1.0)
	//			{
	//				if (isEdgeVoxel(imageData, length, width, height, i, j, k, 0.0) == true) 
	//				{
	//					Point3D p = { (dataType)i, (dataType)j, (dataType)k };
	//					Point3D p_real = getRealCoordFromImageCoord3D(p, segmentOrigin, segmentSpacing, orientation);
	//					Point3D p_pet = getImageCoordFromRealCoord3D(p_real, PETimageOrigin, PETimageSpacing, orientation);
	//					size_t i_pet = (size_t)p_pet.x;
	//					size_t j_pet = (size_t)p_pet.y;
	//					size_t k_pet = (size_t)p_pet.z;
	//					min_i = (i_pet - l_mask > 0) ? i_pet - l_mask : 0;
	//					max_i = (i_pet + l_mask < length_pet - 1) ? i_pet + l_mask : length_pet - 1;
	//					min_j = (j_pet - l_mask > 0) ? j_pet - l_mask : 0;
	//					max_j = (j_pet + l_mask < width_pet - 1) ? j_pet + l_mask : width_pet - 1;
	//					min_k = (k_pet - l_mask > 0) ? k_pet - l_mask : 0;
	//					max_k = (k_pet + l_mask < height_pet - 1) ? k_pet + l_mask : height_pet - 1;
	//					for (size_t kk = min_k; kk <= max_k; kk++) 
	//					{
	//						for (size_t ii = min_i; ii <= max_i; ii++)
	//						{
	//							for (size_t jj = min_j; jj <= max_j; jj++)
	//							{
	//								ballLiverPet[kk][x_new(ii, jj, length_pet)] = 1.0;
	//							}
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}
	//}
	//storing_path = outputPath + "P1/mask_edgeOnevoxel.raw";
	//manageRAWFile3D<dataType>(ballLiverPet, length_pet, width_pet, height_pet, storing_path.c_str(), STORE_DATA, false);
	
	////Get statistiques at aorta edge
	//size_t count_voxels = 0;
	//dataType min_value = 1e6, max_value = -1e6, mean_value = 0.0;
	//for(k = 0; k < height; k++)
	//{
	//	for(i = 0; i < length; i++)
	//	{
	//		for(j = 0; j < width; j++)
	//		{
	//			if(imageData[k][x_new(i, j, length)] == 1.0 && isEdgeVoxel(imageData, length, width, height, i, j, k, 0.0) == true)
	//			{
	//				//maskData[k][x_new(i, j, length)] = 1.0;
	//				Point3D p = { (dataType)i, (dataType)j, (dataType)k };
	//				Point3D p_real = getRealCoordFromImageCoord3D(p, segmentOrigin, segmentSpacing, orientation);
	//				Point3D p_pet = getImageCoordFromRealCoord3D(p_real, PETimageOrigin, PETimageSpacing, orientation);
	//				size_t i_pet = (size_t)p_pet.x;
	//				size_t j_pet = (size_t)p_pet.y;
	//				size_t k_pet = (size_t)p_pet.z;
	//				
	//				//dataType pet_value = petContainer->dataPointer[k_pet][x_new(i_pet, j_pet, length_pet)];
	//				//if (pet_value < min_value)
	//				//{
	//				//	min_value = pet_value;
	//				//}
	//				//if (pet_value > max_value)
	//				//{
	//				//	max_value = pet_value;
	//				//}
	//				//if (statusPetBall[k_pet][x_new(i_pet, j_pet, length_pet)] == false)
	//				//{
	//				//	mean_value += pet_value;
	//				//	count_voxels++;
	//				//	statusPetBall[k_pet][x_new(i_pet, j_pet, length_pet)] == true;
	//				//}
	//				Statistics statistics = getStatisticsInNeighborhood3D(petContainer->dataPointer, length_pet, width_pet, height_pet, i_pet, j_pet, k_pet, 4);
	//				dataType pet_value_min = statistics.min_data;
	//				dataType pet_value_max = statistics.max_data;
	//				dataType pet_value_mean = statistics.mean_data;
	//				
	//				if(pet_value_min < min_value)
	//				{
	//					min_value = pet_value_min;
	//				}
	//				if (pet_value_max > max_value)
	//				{
	//					max_value = pet_value_max;
	//				}
	//				if (statusPetBall[k_pet][x_new(i_pet, j_pet, length_pet)] == false) 
	//				{
	//					mean_value += pet_value_mean;
	//					count_voxels++;
	//					statusPetBall[k_pet][x_new(i_pet, j_pet, length_pet)] == true;
	//				}
	//			}
	//		}
	//	}
	//}
	//std::cout << "Nb aorta edge voxels: " << count_voxels << std::endl;
	//mean_value /= (dataType)count_voxels;
	//std::cout << "Aorta edge PET values: min = " << min_value << ", max = " << max_value << ", mean = " << mean_value << std::endl;

	loading_path = outputPath + "Segmentation/Aorta/p2/segment_aorta_full_dim_p2.raw";
	manageRAWFile3D<dataType>(imageDataFull, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	//fastMarchingDistanceMap(inputImageStr, distanceMap, 0.0);
	storing_path = outputPath + "distance_map_aorta_p2.raw";
	manageRAWFile3D<dataType>(distanceMap, Length, Width, Height, storing_path.c_str(), LOAD_DATA, false);

	////File to save ratios
	//string saving_ratios_csv = outputPath + "ratios_after.csv";
	//FILE* f_ratios;
	//if (fopen_s(&f_ratios, saving_ratios_csv.c_str(), "w") != 0)
	//{
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(f_ratios, "index,radius,ratio_max\n");

	// Input centered path
	FILE* path_file;
	loading_path = outputPath + "Segmentation/Aorta/centered paths/finalCurve_p1.csv";
	if (fopen_s(&path_file, loading_path.c_str(), "r") != 0) 
	{
		printf("Enable to open");
		return false;
	}
	
	size_t index = 0;
	dataType x = 0.0, y = 0.0, z = 0.0;
	dataType aorta_max_SUV = 0, aorta_pet_value = 0.0;
	
	while (feof(path_file) == 0)
	{
		fscanf_s(path_file, "%f", &x);
		fscanf_s(path_file, ",");
		fscanf_s(path_file, "%f", &y);
		fscanf_s(path_file, ",");
		fscanf_s(path_file, "%f", &z);
		fscanf_s(path_file, "\n");
		Point3D current_point = { x, y, z };
		index++;

		Point3D current_point_ct = getImageCoordFromRealCoord3D(current_point, imageOrigin, imageSpacing, orientation);
		
		size_t i_seg = (size_t)(current_point_ct.x);
		size_t j_seg = (size_t)(current_point_ct.y);
		size_t k_seg = (size_t)(current_point_ct.z);
		double dist_seg = distanceMap[k_seg][x_new(i_seg, j_seg, Length)];
		dataType offset_distance = 2.0 * fmax(PETimageSpacing.sx, fmax(PETimageSpacing.sy, PETimageSpacing.sz));
		double radius_aorta_roi = dist_seg + offset_distance;

		//index <= 25 ----> ascending aorta
		//index > 35 && index <= 150 ----> aortic arch
		//index > 175 && index <= 275 ----> descending aorta
		//index > 290 ----> abdominal aorta
		if(index <= 25)
		{
			//aorta_max_SUV = 0.0;
			//count_voxels = 0;
			BoundingBox3D box = findBoundingBox3D(current_point_ct, Length, Width, Height, radius_aorta_roi, offset_distance);
			for (size_t kk = box.k_min; kk <= box.k_max; kk++)
			{
				for (size_t ii = box.i_min; ii <= box.i_max; ii++)
				{
					for (size_t jj = box.j_min; jj <= box.j_max; jj++)
					{
						Point3D p_ct = { (dataType)ii, (dataType)jj, (dataType)kk };
						Point3D p_real = getRealCoordFromImageCoord3D(p_ct, imageOrigin, imageSpacing, orientation);
						double dist = getPoint3DDistance(p_real, current_point);

						if (dist <= radius_aorta_roi)
						{
							if (imageDataFull[kk][x_new(ii, jj, Length)] == 1.0)
							{
								maskAorta[kk][x_new(ii, jj, Length)] = 1.0;
								Point3D p_pet = getImageCoordFromRealCoord3D(p_real, PETimageOrigin, PETimageSpacing, orientation);
								size_t ii_pet = (size_t)p_pet.x;
								size_t jj_pet = (size_t)p_pet.y;
								size_t kk_pet = (size_t)p_pet.z;
								aorta_pet_value = petContainer->dataPointer[kk_pet][x_new(ii_pet, jj_pet, length_pet)];
								count_voxels++;
								if (aorta_pet_value > aorta_max_SUV)
								{
									aorta_max_SUV = aorta_pet_value;
								}
							}
						}
					}
				}
			}
			//dataType ratio_max_aorta_mean_liver = aorta_max_SUV / mean_SUV_Liver;
			//fprintf(f_ratios, "%d,%lf,%lf,%lf\n", index, radius_aorta_roi, ratio_max_aorta_mean_liver);
		}
		
		////aorta_max_SUV = 0.0;
		////mean_SUV_aorta = 0.0;
		////count_voxels = 0;
		//BoundingBox3D box = findBoundingBox3D(current_point_ct, length, width, height, radius_aorta_roi, offset_distance);
		//for (size_t kk = box.k_min; kk <= box.k_max; kk++)
		//{
		//	for (size_t ii = box.i_min; ii <= box.i_max; ii++)
		//	{
		//		for (size_t jj = box.j_min; jj <= box.j_max; jj++)
		//		{
		//			Point3D p_ct = { (dataType)ii, (dataType)jj, (dataType)kk };
		//			Point3D p_real = getRealCoordFromImageCoord3D(p_ct, segmentOrigin, segmentSpacing, orientation);
		//			double dist = getPoint3DDistance(p_real, current_point);
		//			
		//			if (dist <= radius_aorta_roi)
		//			{
		//				//maskData[kk][x_new(ii, jj, length)] = 1.0;
		//				if(imageData[kk][x_new(ii, jj, length)] == 1.0)
		//				{
		//					//maskAorta[kk][x_new(ii, jj, length)] = 1.0;
		//					Point3D p_pet = getImageCoordFromRealCoord3D(p_real, PETimageOrigin, PETimageSpacing, orientation);
		//					size_t ii_pet = (size_t)p_pet.x;
		//					size_t jj_pet = (size_t)p_pet.y;
		//					size_t kk_pet = (size_t)p_pet.z;
		//					aorta_pet_value = petContainer->dataPointer[kk_pet][x_new(ii_pet, jj_pet, length_pet)];
		//					if (aorta_pet_value > aorta_max_SUV)
		//					{
		//						aorta_max_SUV = aorta_pet_value;
		//					}
		//				}
		//			}
		//			
		//			////First test: full sphere
		//			//if (dist <= radius_aorta_roi)
		//			//{
		//			//	Point3D p_pet = getImageCoordFromRealCoord3D(p_real, PETimageOrigin, PETimageSpacing, orientation);
		//			//	size_t ii_pet = (size_t)p_pet.x;
		//			//	size_t jj_pet = (size_t)p_pet.y;
		//			//	size_t kk_pet = (size_t)p_pet.z;
		//			//	maskData[kk][x_new(ii, jj, length)] = 1.0;
		//			//	aorta_pet_value = petContainer->dataPointer[kk_pet][x_new(ii_pet, jj_pet, length_pet)];
		//			//	count_voxels++;
		//			//	mean_SUV_aorta += aorta_pet_value;
		//			//	if (aorta_pet_value > aorta_max_SUV)
		//			//	{
		//			//		aorta_max_SUV = aorta_pet_value;
		//			//	}
		//			//}
		//			////Second test: intersection sphere-segmented aorta
		//			//if (dist <= radius_aorta_roi)
		//			//{
		//			//	if (imageData[kk][x_new(ii, jj, length)] == 1.0)
		//			//	{
		//			//		Point3D p_pet = getImageCoordFromRealCoord3D(p_real, PETimageOrigin, PETimageSpacing, orientation);
		//			//		size_t ii_pet = (size_t)p_pet.x;
		//			//		size_t jj_pet = (size_t)p_pet.y;
		//			//		size_t kk_pet = (size_t)p_pet.z;
		//			//		maskData[kk][x_new(ii, jj, length)] = 1.0;
		//			//		aorta_pet_value = petContainer->dataPointer[kk_pet][x_new(ii_pet, jj_pet, length_pet)];
		//			//		count_voxels++;
		//			//		mean_SUV_aorta += aorta_pet_value;
		//			//		if (aorta_pet_value > aorta_max_SUV)
		//			//		{
		//			//			aorta_max_SUV = aorta_pet_value;
		//			//		}
		//			//	}
		//			//}
		//			////Third test: shell
		//			//if (dist >= half_radius && dist <= radius_aorta_roi)
		//			//{
		//			//	Point3D p_pet = getImageCoordFromRealCoord3D(p_real, PETimageOrigin, PETimageSpacing, orientation);
		//			//	size_t ii_pet = (size_t)p_pet.x;
		//			//	size_t jj_pet = (size_t)p_pet.y;
		//			//	size_t kk_pet = (size_t)p_pet.z;
		//			//	maskData[kk][x_new(ii, jj, length)] = 1.0;
		//			//	aorta_pet_value = petContainer->dataPointer[kk_pet][x_new(ii_pet, jj_pet, length_pet)];
		//			//	count_voxels++;
		//			//	mean_SUV_aorta += aorta_pet_value;
		//			//	if (aorta_pet_value > aorta_max_SUV)
		//			//	{
		//			//		aorta_max_SUV = aorta_pet_value;
		//			//	}
		//			//}
		//			////Fouth test: shell min max aorta radiuses
		//			//size_t xdp = x_new(ii, jj, length);
		//			////dist >= min_dist_aorta &&
		//			//if (dist >= half_radius && dist <= radius_aorta_roi)
		//			//{
		//			//	Point3D p_pet = getImageCoordFromRealCoord3D(p_real, PETimageOrigin, PETimageSpacing, orientation);
		//			//	size_t ii_pet = (size_t)p_pet.x;
		//			//	size_t jj_pet = (size_t)p_pet.y;
		//			//	size_t kk_pet = (size_t)p_pet.z;
		//			//	maskData[kk][xdp] = 1.0;
		//			//	aorta_pet_value = petContainer->dataPointer[kk_pet][x_new(ii_pet, jj_pet, length_pet)];
		//			//	count_voxels++;
		//			//	mean_SUV_aorta += aorta_pet_value;
		//			//	if (aorta_pet_value > aorta_max_SUV)
		//			//	{
		//			//		aorta_max_SUV = aorta_pet_value;
		//			//	}
		//			//	//statusPoints[kk][xdp] = true;
		//			//}
		//			////Get Statistics along the aorta
		//			//if(dist <= radius_aorta_roi && imageData[kk][x_new(ii, jj, length)] == 1.0)
		//			//{
		//			//	Point3D p_pet = getImageCoordFromRealCoord3D(p_real, PETimageOrigin, PETimageSpacing, orientation);
		//			//	size_t ii_pet = (size_t)p_pet.x;
		//			//	size_t jj_pet = (size_t)p_pet.y;
		//			//	size_t kk_pet = (size_t)p_pet.z;
		//			//	aorta_pet_value = petContainer->dataPointer[kk_pet][x_new(ii_pet, jj_pet, length_pet)];
		//			//	if(aorta_pet_value > aorta_max_SUV)
		//			//	{
		//			//		aorta_max_SUV = aorta_pet_value;
		//			//	}
		//			//	if (aorta_pet_value < aorta_min_SUV)
		//			//	{
		//			//		aorta_min_SUV = aorta_pet_value;
		//			//	}
		//			//	//if(statusPetBall[kk_pet][x_new(ii_pet, jj_pet, length_pet)] == false)
		//			//	//{
		//			//	//	mean_SUV_aorta += aorta_pet_value;
		//			//	//	count_voxels++;
		//			//	//	statusPetBall[kk_pet][x_new(ii_pet, jj_pet, length_pet)] = true;
		//			//	//}
		//			//	mean_SUV_aorta += aorta_pet_value;
		//			//	count_voxels++;
		//			//}
		//		}
		//	}
		//}

		//mean_SUV_aorta /= (dataType)count_voxels;
		//std::cout << "Statistics along the aorta: min = " << aorta_min_SUV << ", max = " << aorta_max_SUV << ", mean = " << mean_SUV_aorta << std::endl;
		//dataType ratio_mean_aorta_mean_liver = mean_SUV_aorta / mean_SUV_Liver;
		//dataType ratio_max_aorta_mean_liver = aorta_max_SUV / mean_SUV_Liver;
		//fprintf(f_ratios, "%d,%lf,%lf,%lf\n", index, radius_aorta_roi, ratio_max_aorta_mean_liver, ratio_mean_aorta_mean_liver);	
		
	}
	//fclose(f_ratios);
	fclose(path_file);

	std::cout << "Max aorta ratio : " << aorta_max_SUV / mean_SUV_Liver << std::endl;

	//std::cout << "Abdominal max ratio : " << aorta_max_SUV / mean_SUV_Liver << std::endl;
	//std::cout << "Descending max ratio : " << aorta_max_SUV / mean_SUV_Liver << std::endl;
	//std::cout << "Aortic arch max ratio : " << aorta_max_SUV / mean_SUV_Liver << std::endl;
	
	//storing_path = outputPath + "P1/ascending_aorta.raw";
	//storing_path = outputPath + "P1/aortic_arch.raw";
	//storing_path = outputPath + "P2 translated/abdominal_aorta.raw";
	//manageRAWFile3D<dataType>(maskAorta, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	////std::cout << "Aorta min distance to centerline: " << min_dist_aorta << " mm" << std::endl;
	////std::cout << "Aorta max distance to centerline: " << max_dist_aorta << " mm" << std::endl;
	////storing_path = outputPath + "P1/mask_aorta_fourth.raw";
	////manageRAWFile3D<dataType>(maskData, length, width, height, storing_path.c_str(), STORE_DATA, false);
	//////Analysis at the point of interest
	////std::cout << "Point of interest at index 265: (" << point_of_interest.x << ", " << point_of_interest.y << ", " << point_of_interest.z << ")" << std::endl;
	////Point3D point_ct = getImageCoordFromRealCoord3D(point_of_interest, segmentOrigin, segmentSpacing, orientation);
	////size_t i_poi = (size_t)(point_ct.x);
	////size_t j_poi = (size_t)(point_ct.y);
	////size_t k_poi = (size_t)(point_ct.z);
	////double dist_poi = distanceMap[k_poi][x_new(i_poi, j_poi, length)] + offset_distance;
	////std::cout << "Distance to aorta at point of interest: " << dist_poi << " mm" << std::endl;
	//////Draw spheres at point of interest
	////string saving_sphere_csv = outputPath + "P1/spheres_tests.csv";
	////FILE* f_sphere;
	////if (fopen_s(&f_sphere, saving_sphere_csv.c_str(), "w") != 0)
	////{
	////	printf("Enable to open");
	////	return false;
	////}
	////fprintf(f_ratios, "x,y,z\n");
	////fprintf(f_sphere, "%f,%f,%f\n", point_of_interest.x, point_of_interest.y, point_of_interest.z);
	////box = findBoundingBox3D(point_ct, length, width, height, max_dist_aorta, 2.0);
	//
	////for(size_t kk = box.k_min; kk <= box.k_max; kk++)
	////{
	////	for (size_t ii = box.i_min; ii <= box.i_max; ii++)
	////	{
	////		for (size_t jj = box.j_min; jj <= box.j_max; jj++)
	////		{
	////			Point3D p_ct = { (dataType)ii, (dataType)jj, (dataType)kk };
	////			Point3D p_real = getRealCoordFromImageCoord3D(p_ct, segmentOrigin, segmentSpacing, orientation);
	////			double dist = getPoint3DDistance(p_real, point_of_interest);
	////
	////			if (dist >= min_dist_aorta && dist <= max_dist_aorta)
	////			{
	////				fprintf(f_sphere, "%f,%f,%f\n", p_real.x, p_real.y, p_real.z);
	////			}
	////		}
	////	}
	////}
	////fclose(f_sphere);

	//for (k = 0; k < Height; k++)
	//{
	//	if (k < height_pet)
	//	{
	//		delete[] ballLiverPet[k];
	//		delete[] statusPetBall[k];
	//	}
	//	delete[] maskLiver[k];
	//}
	//delete[] ballLiverPet;
	//delete[] statusPetBall;
	//delete[] maskLiver;

	for (k = 0; k < Height; k++)
	{
		delete[] distanceMap[k];
		delete[] maskAorta[k];
		delete[] imageDataFull[k];
		delete[] maskLiver[k];
	}
	delete[] imageDataFull;
	delete[] distanceMap;
	delete[] maskAorta;
	delete[] maskLiver;

	free(petContainer);
	free(ctContainer);
	*/

	//==================== Adjust segment for quantitative analysis ===================================
	
	/*
	const size_t length = 150, width = 150, height = 350;

	dataType** imageData = new dataType * [height];
	for (k = 0; k < height; k++) 
	{
		imageData[k] = new dataType[length * width]{ 0 };
	}

	loading_path = outputPath + "Segmentation/3D image/p1/_seg_func_05000.raw";
	manageRAWFile3D<dataType>(imageData, length, width, height, loading_path.c_str(), LOAD_DATA, false);

	dataType** imageDataFull = new dataType * [Height];
	for (k = 0; k < Height; k++)
	{
		imageDataFull[k] = new dataType[Length * Width]{ 0 };
	}

	////Patient 1
	Point3D segmentOrigin = { -65.625, 4.375, -700.234 };
	VoxelSpacing segmentSpacing = { 1.171875, 1.171875, 1.171875 };

	////Patient 2
	//Point3D segmentOrigin = { -77.3438, 16.0938, -780.316 };
	//VoxelSpacing segmentSpacing = { 1.171875, 1.171875, 1.171875 };

	for (k = 0; k < height; k++) 
	{
		for(i = 0; i < length; i++) 
		{
			for (j = 0; j < width; j++) 
			{
				Point3D p_ct = { (dataType)i, (dataType)j, (dataType)k };
				Point3D p_real = getRealCoordFromImageCoord3D(p_ct, segmentOrigin, segmentSpacing, orientation);
				Point3D p_full = getImageCoordFromRealCoord3D(p_real, imageOrigin, imageSpacing, orientation);
				size_t i_full = (size_t)(p_full.x);
				size_t j_full = (size_t)(p_full.y);
				size_t k_full = (size_t)(p_full.z);
				size_t xd_full = x_new(i_full, j_full, Length);

				xd = x_new(i, j, length);
				if (imageData[k][xd] >= 0.2)
				{
					imageDataFull[k_full][xd_full] = 1.0;
				}
				else
				{
					imageDataFull[k_full][xd_full] = 0.0;
				}
			}
		}
	}

	loading_path = outputPath + "Segmentation/Aorta/p1/segment_aorta_full_dim_p1.raw";
	manageRAWFile3D<dataType>(imageDataFull, Length, Width, Height, loading_path.c_str(), STORE_DATA, false);

	for(k = 0; k < Height; k++) 
	{
		if (k < height)
		{
			delete[] imageData[k];
		}
		delete[] imageDataFull[k];
	}
	delete[] imageData;
	delete[] imageDataFull;

	free(ctContainer);
	*/

	//==================== Analyse extracted path curvature ===========================================
	
	/*
	FILE* path_file;
	//loading_path = outputPath + "Segmentation/Aorta/centered paths/finalCurve_p1.csv";
	loading_path = root + "Curves/Output/smooth_v2.csv";
	if (fopen_s(&path_file, loading_path.c_str(), "r") != 0)
	{
		printf("Enable to open");
		return false;
	}

	vector<Point3D> pPoints;
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
		pPoints.push_back(current_point);
	}
	fclose(path_file);

	const size_t nb_path_points = pPoints.size();

	////Multiple
	//string saving_csv = outputPath + "P5/curvature.csv";
	//FILE* f_curvature;
	//if (fopen_s(&f_curvature, saving_csv.c_str(), "w") != 0)
	//{
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(f_curvature, "Index,x,y,z,Curvature,Tx,Ty,Tz\n");
	//
	//dataType* norm_save = new dataType[nb_path_points]{0};
	//if(norm_save == NULL)
	//{
	//	return false;
	//}
	//
	//dataType** tangent = new dataType*[3];
	//for(i = 0; i < 3; i++)
	//{
	//	tangent[i] = new dataType[nb_path_points]{0};
	//}
	//
	//string save_curvature_max;
	//dataType h_i, h_i_plus, coef, norm_r;
	//dataType r_x, r_y, r_z, Tx, Ty, Tz;
	//size_t index = 0;
	//size_t iminus, iplus, icurrent, offset = 0;
	//size_t res, qt;
	//Point3D pointMax = {0, 0, 0};
	//dataType maxNorm = 0.0;
	//while(offset < 10)
	//{
	//	offset++;
	//	index = 0;
	//	maxNorm = 0.0;
	//	pointMax = { 0, 0, 0 };
	//	res = (pPoints.size() - 1) % offset;
	//	qt = (pPoints.size() - 1) / offset;
	//	for (i = 1; i < qt; i++)
	//	{
	//		index++;
	//		icurrent = i * offset;
	//		iminus = (i - 1) * offset;
	//		iplus = (i + 1) * offset;
	//		h_i = getPoint3DDistance(pPoints[iminus], pPoints[icurrent]);
	//		h_i_plus = getPoint3DDistance(pPoints[iplus], pPoints[icurrent]);
	//		coef = 2.0 / (h_i + h_i_plus);
	//		r_x = coef * (((pPoints[iplus].x - pPoints[icurrent].x) / h_i_plus) - ((pPoints[icurrent].x - pPoints[iminus].x) / h_i));
	//		r_y = coef * (((pPoints[iplus].y - pPoints[icurrent].y) / h_i_plus) - ((pPoints[icurrent].y - pPoints[iminus].y) / h_i));
	//		r_z = coef * (((pPoints[iplus].z - pPoints[icurrent].z) / h_i_plus) - ((pPoints[icurrent].z - pPoints[iminus].z) / h_i));
	//		norm_r = sqrt(r_x * r_x + r_y * r_y + r_z * r_z);
	//		if (maxNorm < norm_r)
	//		{
	//			maxNorm = norm_r;
	//			pointMax = pPoints[icurrent];
	//		}	
	//
	//		dataType coef2 = 1.0 / (h_i + h_i_plus);
	//		tangent[0][i] = coef2 * (pPoints[iplus].x - pPoints[iminus].x);
	//		tangent[1][i] = coef2 * (pPoints[iplus].y - pPoints[iminus].y);
	//		tangent[2][i] = coef2 * (pPoints[iplus].z - pPoints[iminus].z);
	//		norm_save[i] = norm_r;
	//	}
	//
	//	save_curvature_max = outputPath + "P5/New/curvature_max_" + to_string(offset) + ".csv";
	//	FILE* f_curvature_max;
	//	if (fopen_s(&f_curvature_max, save_curvature_max.c_str(), "w") != 0)
	//	{
	//		printf("Enable to open");
	//		return false;
	//	}
	//	fprintf(f_curvature_max, "x,y,z\n");
	//	fprintf(f_curvature_max, "%lf,%lf,%lf\n", pointMax.x, pointMax.y, pointMax.z);
	//	fclose(f_curvature_max);
	//}
	//
	//for(i = 1; i < (nb_path_points - 1); i++)
	//{
	//	fprintf(f_curvature, "%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", i, pPoints[icurrent].x, pPoints[icurrent].y, pPoints[icurrent].z, norm_save[i], tangent[0][i], tangent[1][i], tangent[2][i]);
	//}
	//fclose(f_curvature);
	//	for(i = 0; i < 3; i++)
	//{
	//	delete[] tangent[i];
	//}
	//delete[] tangent;
	//delete[] norm_save;
	
	string saving_csv = root + "curvature_smooth.csv";
	//string saving_csv = root + "output/curvature_float.csv";
	FILE* f_curvature;
	if (fopen_s(&f_curvature, saving_csv.c_str(), "w") != 0)
	{
		printf("Enable to open");
		return false;
	}

	dataType h_i, h_i_plus, coef, norm_r, max_curv = 0.0;
	dataType r_x, r_y, r_z;
	size_t index = 0, icurrent, iminus, iplus;
	Point3D pointMax = { 0.0, 0.0, 0.0 };
	for(i = 1; i < (nb_path_points - 1); i++)
	{
		index++;
		icurrent = i;
		iminus = (i - 1);
		iplus = (i + 1);
		h_i = getPoint3DDistance(pPoints[iminus], pPoints[icurrent]);
		h_i_plus = getPoint3DDistance(pPoints[iplus], pPoints[icurrent]);
		coef = 2.0 / (h_i + h_i_plus);
		r_x = coef * (((pPoints[iplus].x - pPoints[icurrent].x) / h_i_plus) - ((pPoints[icurrent].x - pPoints[iminus].x) / h_i));
		r_y = coef * (((pPoints[iplus].y - pPoints[icurrent].y) / h_i_plus) - ((pPoints[icurrent].y - pPoints[iminus].y) / h_i));
		r_z = coef * (((pPoints[iplus].z - pPoints[icurrent].z) / h_i_plus) - ((pPoints[icurrent].z - pPoints[iminus].z) / h_i));
		norm_r = sqrt(r_x * r_x + r_y * r_y + r_z * r_z);
		fprintf(f_curvature, "%d,%lf\n", index, norm_r);
		//if(norm_r > max_curv)
		//{
		//	max_curv = norm_r;
		//	pointMax = pPoints[icurrent];
		//}
	}
	fclose(f_curvature);

	//string saving_max_curv_pt_csv = outputPath + "max_curvature_pm.csv";
	//FILE* f_curv_max;
	//if (fopen_s(&f_curv_max, saving_max_curv_pt_csv.c_str(), "w") != 0)
	//{
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(f_curv_max, "x,y,z\n");
	//fprintf(f_curv_max, "%lf,%lf,%lf\n", pointMax.x, pointMax.y, pointMax.z);
	//fclose(f_curv_max);
	*/
	
	
	srand(time(NULL));
	const char* test_circle = "C:/Users/Konan Allaly/Documents/Tests/output/test_curve.csv";
	FILE* file_test;
	if (fopen_s(&file_test, test_circle, "w") != 0) {
		printf("Enable to open");
		return false;
	}

	vector<Point3D> circlePoints;

	/*
	size_t N = 100;
	double angle_step = M_PI / (double)N;
	double angle_rad = 0.0;
	//half cricle in the xy plane
	for (i = 0; i <= N; i++) {
		angle_rad = i * angle_step;
		double x = 50.0 + radius * cos(angle_rad);
		double y = 50.0 + radius * sin(angle_rad);
		double z = 50.0;
		Point3D point = { x, y, z };
		
		//if (i > 0 && i < N)
		//{
		//	double r_x = (double)rand() / RAND_MAX;
		//	x += 0.4 * r_x;
		//	double r_y = (double)rand() / RAND_MAX;
		//	y += 0.4 * r_y;
		//}
		
		fprintf(file_test, "%lf,%lf,%lf\n", x, y, z);
		circlePoints.push_back(point);
	}
	
	//straight line
	size_t N2 = 50;
	for(i = 1; i < N2 - 2; i++)
	{
		double x = 50.0 + radius * cos(angle_rad);
		double y = 50.0 + radius * sin(angle_rad) - i * angle_step;
		double z = 50.0;
		Point3D point = { x, y, z };
		fprintf(file_test, "%f,%f,%f\n", x, y, z);
		circlePoints.push_back(point);
	}

	//straight line
	for (i = 1; i < (N2 - 2); i++)
	{
		double x = 50.0 + radius * cos(angle_rad) - i * angle_step;
		double y = 50.0 + radius * sin(angle_rad) - (N2 - 3) * angle_step;
		double z = 50.0;
		Point3D point = { x, y, z };
		fprintf(file_test, "%f,%f,%f\n", x, y, z);
		circlePoints.push_back(point);
	}

	for (i = 1; i < N; i++) {
		angle_rad = i * angle_step;
		double x = 50.0 - 2.5 * radius;
		double y = 50.0 - radius * (1.5 - sin(angle_rad));
		double z = 50.0 + radius * (-1.0 + cos(angle_rad));
		Point3D point = { x, y, z };

		fprintf(file_test, "%lf,%lf,%lf\n", x, y, z);
		circlePoints.push_back(point);
	}
	*/

	Point3D p1 = { 0.0, 0.0, 0.0 };
	Point3D p2 = { 2.0, 0.0, 0.0 };
	double radius = 0.5 * getPoint3DDistance(p1, p2);

	Point3D center = { 0.5 * (p1.x + p2.x), 0.5 * (p1.y + p2.y), 0.5 * (p1.z + p2.z) };
	double angle_step = M_PI / 100.0;
	for (i = 0; i <= 100; i++) 
	{
		double angle_rad = (double)i * angle_step;
		double x = center.x + radius * cos(angle_rad);
		double y = center.y + radius * sin(angle_rad);
		double z = center.z;
		//add noise
		if (i > 0 && i < 100 && (i % 5 == 0))
		{
			double r_x = (double)rand() / RAND_MAX;
			x += 0.02 * r_x;
			double r_y = (double)rand() / RAND_MAX;
			y += 0.02 * r_y;
		}
		Point3D point = { x, y, z };
		circlePoints.push_back(point);
	}

	for (i = 1; i <= 50; i++) 
	{
		double x = p1.x;
		double y = p1.y - (double)i * angle_step;
		double z = p1.z;
		//add noise
		if (i > 0 && i < 50 && (i % 5 == 0))
		{
			double r_x = (double)rand() / RAND_MAX;
			x += 0.02 * r_x;
			double r_y = (double)rand() / RAND_MAX;
			y += 0.02 * r_y;
		}
		Point3D point = { x, y, z };
		circlePoints.push_back(point);
	}

	Point3D p4 = { 2.0, -0.5 * M_PI * radius, 0.0 };
	Point3D p3 = { 0.0, -0.5 * M_PI * radius, 0.0 };

	size_t N2 = (size_t)(getPoint3DDistance(p3, p4) / angle_step);

	for (i = 1; i <= N2; i++)
	{
		double x = p4.x - (double)i * angle_step;
		double y = p4.y;
		double z = p4.z;
		//add noise
		if (i > 0 && i < N2 && (i % 5 == 0))
		{
			double r_x = (double)rand() / RAND_MAX;
			x += 0.02 * r_x;
			double r_y = (double)rand() / RAND_MAX;
			y += 0.02 * r_y;
		}
		Point3D point = { x, y, z };
		circlePoints.push_back(point);
	}

	Point3D p5 = { 2.0, -0.5 * M_PI * radius, -2.0 };
	Point3D center2 = { 0.5 * (p4.x + p5.x), 0.5 * (p4.y + p5.y), 0.5 * (p4.z + p5.z) };
	
	for (i = 0; i < 100; i++)
	{
		double angle_rad = i * angle_step;
		double x = center2.x;
		double y = center2.y - radius * sin(angle_rad);
		double z = center2.z + radius * cos(angle_rad);
		//add noise
		if (i > 0 && i < 100 && (i % 5 == 0))
		{
			double r_x = (double)rand() / RAND_MAX;
			x += 0.02 * r_x;
			double r_y = (double)rand() / RAND_MAX;
			y += 0.02 * r_y;
		}
		Point3D point = { x, y, z };
		circlePoints.push_back(point);
	}
	
	for(i = 0; i < circlePoints.size(); i++)
	{
		fprintf(file_test, "%lf,%lf,%lf\n", circlePoints[i].x, circlePoints[i].y, circlePoints[i].z);
	}
	
	fclose(file_test);

	//Compute and save curvature
	const char* test_curvature = "C:/Users/Konan Allaly/Documents/Tests/output/curvature_double.csv";
	FILE* file_curvature;
	if (fopen_s(&file_curvature, test_curvature, "w") != 0) 
	{
		printf("Enable to open");
		return false;
	}

	dataType ref_curv = 1.0 / radius;
	dataType max_error = 0.0, min_error = 1.0;
	dataType mean_error = 0.0;
	size_t size_pts = circlePoints.size();
	dataType* curvature = new dataType[size_pts]{ 0 };
	size_t counter = 0;
	for(i = 1; i < size_pts - 1; i++)
	{
		curvature[i] = computeCurvatureThreePoints(circlePoints[i - 1], circlePoints[i], circlePoints[i + 1]);
		dataType error = fabs(curvature[i] - ref_curv);
		if(error > max_error)
		{
			max_error = error;
		}
		if(error < min_error)
		{
			min_error = error;
		}
		mean_error += error;
		fprintf(file_curvature, "%d,%f\n", i, curvature[i]);
		counter++;
	}
	mean_error /= (dataType)counter;
	std::cout << "Curvature error: mean = " << mean_error << std::endl;
	std::cout << "Curvature error: min = " << min_error << ", max = " << max_error << std::endl;
	fclose(file_curvature);
	
	delete[] curvature;

	//Compute and save curvature
	const char* test_spacing = "C:/Users/Konan Allaly/Documents/Tests/output/sapcing_double.csv";
	FILE* file_spacing;
	if (fopen_s(&file_spacing, test_spacing, "w") != 0)
	{
		printf("Enable to open");
		return false;
	}
	for(i = 0; i < circlePoints.size() - 1; i++)
	{
		dataType spacing = getPoint3DDistance(circlePoints[i], circlePoints[i + 1]);
		fprintf(file_spacing, "%d,%lf\n", i, spacing);
	}
	fclose(file_spacing);

	circlePoints.clear();

	//==================== Liver Cropping Test ========================================================

	/*
	dataType** imageData = new dataType * [Height];
	for (k = 0; k < Height; k++) 
	{
		imageData[k] = new dataType[Length * Width]{ 0 };
	}
	Image_Data inputImageStr = { Height, Length, Width, imageData, imageOrigin, imageSpacing, orientation };
	loading_path = inputPath + "raw/liver/liver_p2.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	////P1
	//size_t i_min = 120, i_max = 300;
	//size_t j_min = 180, j_max = 340;
	//size_t k_min = 160, k_max = 230;

	//P2
	size_t i_min = 110, i_max = 290;
	size_t j_min = 195, j_max = 355;
	size_t k_min = 260, k_max = 330;

	const size_t length = i_max - i_min;
	const size_t width = j_max - j_min;
	const size_t height = k_max - k_min;

	dataType** croppedImageData = new dataType * [height];
	for (k = 0; k < height; k++) 
	{
		croppedImageData[k] = new dataType[length * width]{ 0 };
	}

	size_t ic, jc, kc;
	for(kc = 0, k = k_min; kc < height; kc++, k++)
	{
		for (ic = 0, i = i_min; ic < length; ic++, i++) 
		{
			for (jc = 0, j = j_min; jc < width; jc++, j++) 
			{
				croppedImageData[kc][x_new(ic, jc, length)] = imageData[k][x_new(i, j, Length)];
			}
		}
	}

	storing_path = outputPath + "liver_cropped_p2.raw";
	manageRAWFile3D<dataType>(croppedImageData, length, width, height, storing_path.c_str(), STORE_DATA, false);

	Point3D croppedOrigin = { i_min, j_min, k_min };
	croppedOrigin = getRealCoordFromImageCoord3D(croppedOrigin, imageOrigin, imageSpacing, orientation);
	std::cout << "Cropped Origin: (" << croppedOrigin.x << ", " << croppedOrigin.y << ", " << croppedOrigin.z << ")" << std::endl;
	std::cout << "Cropped dimensions: (" << length << ", " << width << ", " << height << ")" << std::endl;


	for (k = 0; k < Height; k++) 
	{
		if(k < height)
		{
			delete[] croppedImageData[k];
		}
		delete[] imageData[k];
	}
	delete[] imageData;
	delete[] croppedImageData;
	*/

	//==================== Test Filtering ============================================================

	/*
	dataType** imageData = new dataType * [Height];
	dataType** potential = new dataType * [Height];
	dataType** action = new dataType * [Height];
	for (k = 0; k < Height; k++) 
	{
		imageData[k] = new dataType[dim2D]{ 0 };
		potential[k] = new dataType[dim2D]{ 0 };
		action[k] = new dataType[dim2D]{ 0 };
	}
	Image_Data inputImage = { Height, Length, Width, imageData, imageOrigin, imageSpacing, orientation };

	////Rescaling the image to the range [0, 1] for filtering
	//for (k = 0; k < Height; k++) 
	//{
	//	for (i = 0; i < dim2D; i++) 
	//	{
	//		imageData[k][i] = ctContainer->dataPointer[k][i];
	//	}
	//}
	//rescaleNewRange(imageData, Length, Width, Height, 0, 1, 2076, -1024);

	////Filtering
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
	////heatImplicitRectangularScheme(inputImage, smoothParameters);
	//geodesicMeanCurvature(inputImage, smoothParameters);

	//storing_path = root + "output/gmcf_float.raw";
	storing_path = root + "/output/gmcf_double.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, storing_path.c_str(), LOAD_DATA, false);

	Point3D* endPoints = new Point3D[2];//p1
	endPoints[0] = { 261.0, 257.0, 145.0 };
	endPoints[1] = { 259.0, 250.0, 246.0 };

	double radius = 3.0;
	Potential_Parameters parameters{
		1000, //edge detector coefficient
		0.15,  //threshold, (0.15 --> p1, p2), threshold (0.2 --> p3, p4, p5, p6)
		0.001,//epsilon
		radius
	};
	compute3DPotential(inputImage, potential, endPoints, parameters);

	//storing_path = root + "output/potential_float.raw";
	storing_path = root + "/output/potential_double.raw";
	manageRAWFile3D<dataType>(potential, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	Image_Data toAction = { Height, Length, Width, action, imageOrigin, imageSpacing, orientation };
	partialFrontPropagation(toAction, potential, endPoints);
	//storing_path = root + "output/action_map_float.raw";
	storing_path = root + "output/action_map_double.raw";
	manageRAWFile3D<dataType>(action, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	Path_Parameters parameters_path
	{
		0.8, // tau
		1000,// max number of iterations
		0.8 // tolerance
	};
	Image_Data toPathExtraction = { Height, Length, Width, action, imageOrigin, imageSpacing, orientation };
	vector<Point3D> path_points;
	shortestPath3D(toPathExtraction, endPoints, path_points, parameters_path);
	FILE* path_points_file;
	string save_path_file = root + "output/path_points_double.csv";
	if (fopen_s(&path_points_file, save_path_file.c_str(), "w") != 0) 
	{
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

	delete[] endPoints;
	for (k = 0; k < Height; k++)
	{
		delete[] imageData[k];
		delete[] potential[k];
	}
	delete[] imageData;
	delete[] potential;
	*/
	
	/*
	float** imageDataF = new float * [Height];
	double** imageDataD = new double* [Height];
	for (k = 0; k < Height; k++) 
	{
		imageDataF[k] = new float[dim2D]{ 0 };
		imageDataD[k] = new double[dim2D] { 0 };
	}
	
	storing_path = root + "/output/distance_float.raw";
	manageRAWFile3D<float>(imageDataF, Length, Width, Height, storing_path.c_str(), LOAD_DATA, false);

	storing_path = root + "/output/distance_double.raw";
	manageRAWFile3D<double>(imageDataD, Length, Width, Height, storing_path.c_str(), LOAD_DATA, false);

	double min_diff = 0.0, max_diff = 0.0, mean_diff = 0.0;
	size_t count = 0;
	for(k = 0; k < Height; k++) 
	{
		for (i = 0; i < dim2D; i++) 
		{
			double diff = fabs(imageDataF[k][i] - imageDataD[k][i]);
			if (diff > max_diff)
			{
				max_diff = diff;
			}
			if (diff < min_diff)
			{
				min_diff = diff;
			}
			mean_diff += diff;
			count++;
		}
	}
	mean_diff /= (double)count;
	std::cout << "Difference between float and double images: mean = " << mean_diff << ", min = " << min_diff << ", max = " << max_diff << std::endl;
		for(k = 0; k < Height; k++)
	{
		delete[] imageDataF[k];
		delete[] imageDataD[k];
	}
	delete[] imageDataF;
	delete[] imageDataD;
	*/

	return EXIT_SUCCESS;
}