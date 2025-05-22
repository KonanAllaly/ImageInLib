#include <iostream>
#include <sstream>  
#include <vector>
#include <string.h> 
#include <time.h>

#include "common_math.h"
#include <template_functions.h>
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
#include "hough_transform.h"
#include "../src/edgedetection.h"
#include "enhancement.h"
#include "../src/fast_marching.h"
#include "../src/non_linear_heat_equation.h"
#include "eigen_systems.h"

#define MAX_LINE_LENGTH 1024

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

	/*
	Vtk_File_Info* ctContainer = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	ctContainer->operation = copyFrom;
	loading_path = inputPath + "vtk/petct/ct/Patient1_ct.vtk";
	//loading_path = inputPath + "vtk/petct/aorta/AortaPatient3.vtk";
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
	*/

	//========================= Detect Heart region ==================================
	
	/*
	//Prepare inputs
	dataType** imageData = new dataType * [Height];
	dataType** maskLungs = new dataType * [Height];
	dataType** maskLiver = new dataType * [Height];
	dataType** maskHeart = new dataType * [Height];
	dataType** ball = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{ 0 };
		maskLungs[k] = new dataType[dim2D]{ 0 };
		maskLiver[k] = new dataType[dim2D]{ 0 };
		maskHeart[k] = new dataType[dim2D]{ 0 };
		ball[k] = new dataType[Length * Width]{ 0 };
	}

	//Load filtered input image
	loading_path = inputPath + "/raw/filtered/filteredGMC_p1.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	//Load segmented lungs
	loading_path = inputPath + "/raw/lungs/lungs_p1.raw";
	//loading_path = inputPath + "/raw/lungs/lungs_plus_trachea_p1.raw";
	manageRAWFile3D<dataType>(maskLungs, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	//Load segmented liver
	loading_path = inputPath + "/raw/liver/liver_p1.raw";
	manageRAWFile3D<dataType>(maskLiver, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);
	dilatation3D(maskLiver, Length, Width, Height, 1.0, 0.0);
	dilatation3D(maskLiver, Length, Width, Height, 1.0, 0.0);
	dilatation3D(maskLiver, Length, Width, Height, 1.0, 0.0);
	dilatation3D(maskLiver, Length, Width, Height, 1.0, 0.0);
	dilatation3D(maskLiver, Length, Width, Height, 1.0, 0.0);
	dilatation3D(maskLiver, Length, Width, Height, 1.0, 0.0);

	storing_path = outputPath + "p1_1/liver_test.raw";
	manageRAWFile3D<dataType>(maskLiver, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);
	
	//Find lungs centroid
	dataType* centroid_lungs = new dataType[3]{ 0 };
	centroidImage(maskLungs, centroid_lungs, Height, Length, Width, 0.0);
	Point3D pCentroid = { centroid_lungs[0], centroid_lungs[1] , centroid_lungs[2] };
	double radius = 10;

	////remove slices containing the trachea
	//const size_t trachea_k1 = 266, trachea_k2 = 354;//p1
	//const size_t trachea_k1 = 365, trachea_k2 = 449;//p2
	//const size_t trachea_k1 = 238, trachea_k2 = 285;//p4
	//const size_t trachea_k1 = 649, trachea_k2 = 737;//p5
	//const size_t trachea_k1 = 232, trachea_k2 = 292;//p3
	//const size_t trachea_k1 = 461, trachea_k2 = 519;//p6

	const size_t k_centroid = (size_t)pCentroid.z;

	//Find the farest point to the centroid
	Point3D prCentroid = getRealCoordFromImageCoord3D(pCentroid, ctOrigin, ctSpacing, orientation);
	double farest_distance = 0.0;
	Point3D farest_point;
	for (k = 0; k <= k_centroid; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				if (maskLungs[k][x_new(i, j, Length)] == 1.0) {
					Point3D current_point = { i, j, k };
					current_point = getRealCoordFromImageCoord3D(current_point, ctOrigin, ctSpacing, orientation);
					double fDistance = getPoint3DDistance(current_point, prCentroid);
					if (farest_distance < fDistance) {
						farest_distance = fDistance;
						farest_point.x = i;
						farest_point.y = j;
						farest_point.z = k;
					}
				}
			}
		}
	}
	std::cout << "max distance = " << farest_distance << std::endl;

	//Find lungs bounding box
	size_t k_min = Height, k_max = 0;
	size_t i_min = Length, i_max = 0;
	size_t j_min = Width, j_max = 0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				if (maskLungs[k][x_new(i, j, Length)] == 1.0) {
					if (i_min > i) {
						i_min = i;
					}
					if (j_min > j) {
						j_min = j;
					}
					if (k_min > k) {
						k_min = k;
					}

					if (i_max < i) {
						i_max = i;
					}
					if (j_max < j) {
						j_max = j;
					}
					if (k_max < k) {
						k_max = k;
					}
				}
			}
		}
	}

	//cropped dimension
	const size_t height = k_max - k_min + 1;
	const size_t length = i_max - i_min + 1;
	const size_t width = j_max - j_min + 1;
	std::cout << "New dim : " << length << " x " << width << " x " << height << std::endl;

	Point3D croppedOrigin = { i_min, j_min, k_min };
	croppedOrigin = getRealCoordFromImageCoord3D(croppedOrigin, ctOrigin, ctSpacing, orientation);
	std::cout << "Cropped image origin : (" << croppedOrigin.x << ", " << croppedOrigin.y << ", " << croppedOrigin.z << ")" << std::endl;

	dataType** croppedImage = new dataType * [height];
	dataType** croppedLungs = new dataType * [height];
	dataType** segmentCrop = new dataType * [height];
	dataType** maskThreshold = new dataType * [height];
	for (k = 0; k < height; k++) {
		croppedImage[k] = new dataType[length * width]{ 0 };
		croppedLungs[k] = new dataType[length * width]{ 0 };
		segmentCrop[k] = new dataType[length * width]{ 0 };
		maskThreshold[k] = new dataType[length * width]{ 0 };
	}

	size_t k_ext, i_ext, j_ext;
	for (k = 0, k_ext = k_min; k < height; k++, k_ext++) {
		for (i = 0, i_ext = i_min; i < length; i++, i_ext++) {
			for (j = 0, j_ext = j_min; j < width; j++, j_ext++) {

				croppedImage[k][x_new(i, j, length)] = imageData[k_ext][x_new(i_ext, j_ext, Length)];

				//Cropping according lungs farest point to the centroid
				Point3D current_point = { i, j, k };
				current_point = getRealCoordFromImageCoord3D(current_point, croppedOrigin, ctSpacing, orientation);
				double fDistance = getPoint3DDistance(current_point, prCentroid);
				if (fDistance <= farest_distance) {
					if (maskLungs[k_ext][x_new(i_ext, j_ext, Length)] == 1.0) {
						maskThreshold[k][x_new(i, j, length)] = 0.0;
					}
					else {
						maskThreshold[k][x_new(i, j, length)] = 1.0;
					}
				}

				////remove trachea slices when having trachea segmentation
				//if (k_ext >= trachea_k1 && k_ext <= trachea_k2) {
				//	maskThreshold[k][x_new(i, j, length)] = 0.0;
				//}
				
				//remove trachea slices when having just lungs centroid
				if (k_ext >= k_centroid) {
					maskThreshold[k][x_new(i, j, length)] = 0.0;
				}
				//remove liver pixels
				if (maskLiver[k_ext][x_new(i_ext, j_ext, Length)] == 1.0) {
					maskThreshold[k][x_new(i, j, length)] = 0.0;
				}
			}
		}
	}
	storing_path = outputPath + "p1_1/threshold.raw";
	manageRAWFile3D<dataType>(maskThreshold, length, width, height, storing_path.c_str(), STORE_DATA, false);
	
	
	//Save data information
	string store_information = outputPath + "p1_1/info_cropping_p1.txt";
	FILE* info_test;
	if (fopen_s(&info_test, store_information.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	////Original image
	fprintf(info_test, "Input : patient p1\n");
	fprintf(info_test, "Dimension : Length = %d, Width = %d, Height = %d\n", Length, Width, Height);
	fprintf(info_test, "Origin (%f, %f, %f)\n", ctOrigin.x, ctOrigin.y, ctOrigin.z);
	fprintf(info_test, "Spacing (%f, %f, %f)\n", ctSpacing.sx, ctSpacing.sy, ctSpacing.sz);
	////Description
	fprintf(info_test, "\nThe lungs segmentation is used to find lungs bounding box\n");
	fprintf(info_test, "Cropped image dimension : Length = %d, Width = %d, Height = %d\n", length, width, height);
	fprintf(info_test, "Cropped image origin : (%f, %f, %f)\n", croppedOrigin.x, croppedOrigin.y, croppedOrigin.z);
	fclose(info_test);
	

	//Interpolation for distance map
	size_t new_height = (size_t)((ctSpacing.sz / ctSpacing.sx) * height);
	std::cout << "New height = " << new_height << std::endl;
	VoxelSpacing interpolatedSpacing = { ctSpacing.sx, ctSpacing.sy, ctSpacing.sx };

	dataType** interpolImage = new dataType * [new_height];
	dataType** interpolMask = new dataType * [new_height];
	dataType** distanceMap = new dataType * [new_height];
	dataType** ballMaxDist = new dataType * [new_height];
	for (k = 0; k < new_height; k++) {
		interpolImage[k] = new dataType[length * width]{ 0 };
		interpolMask[k] = new dataType[length * width]{ 0 };
		distanceMap[k] = new dataType[length * width]{ 0 };
		ballMaxDist[k] = new dataType[length * width]{ 0 };
	}

	Image_Data imageDataCRP
	{
		height,
		length,
		width,
		maskThreshold,
		croppedOrigin,
		ctSpacing,
		orientation
	};

	Image_Data imageDataINT
	{
		new_height,
		length,
		width,
		interpolMask,
		croppedOrigin,
		interpolatedSpacing,
		orientation
	};

	imageInterpolation3D(imageDataCRP, imageDataINT, NEAREST_NEIGHBOR);
	storing_path = outputPath + "p1_1/lungs_crop.raw";
	manageRAWFile3D<dataType>(interpolMask, length, width, new_height, storing_path.c_str(), STORE_DATA, false);

	rouyTourinFunction_3D(distanceMap, interpolMask, 0.5, length, width, new_height, 0.4, ctSpacing.sx);

	storing_path = outputPath + "p1_1/distance_map.raw";
	manageRAWFile3D<dataType>(distanceMap, length, width, new_height, storing_path.c_str(), STORE_DATA, false);

	dataType max_distance = getTheMaxValue3D(distanceMap, length, width, new_height);

	Point3D point_max = getPointWithTheHighestValue(distanceMap, length, width, new_height);
	Point3D point_max_real_coord = getRealCoordFromImageCoord3D(point_max, croppedOrigin, interpolatedSpacing, orientation);

	//Draw ball for visualization
	for (k = 0; k < new_height; k++) {
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {
				Point3D current_point = { i, j, k };
				current_point = getRealCoordFromImageCoord3D(current_point, croppedOrigin, interpolatedSpacing, orientation);
				double pDistance = getPoint3DDistance(current_point, point_max_real_coord);
				if (pDistance <= max_distance) {
					ballMaxDist[k][x_new(i, j, length)] = 1.0;
				}
			}
		}
	}
	storing_path = outputPath + "p1_1/ball_max_distance.raw";
	manageRAWFile3D<dataType>(ballMaxDist, length, width, new_height, storing_path.c_str(), STORE_DATA, false);

	////Segment the heart by region growing
	//Point3D seedRG = getImageCoordFromRealCoord3D(point_max_real_coord, croppedOrigin, ctSpacing, orientation);
	//
	//Image_Data imageDataRG
	//{
	//	height,
	//	length,
	//	width,
	//	croppedImage,
	//	croppedOrigin,
	//	ctSpacing,
	//	orientation
	//};
	//
	//regionGrowing(imageDataRG, segmentCrop, length, width, height, seedRG, radius);
	//storing_path = outputPath + "p1/segment.raw";
	//manageRAWFile3D<dataType>(segmentCrop, length, width, height, storing_path.c_str(), STORE_DATA, false);

	for (k = 0; k < new_height; k++) {
		delete[] interpolImage[k];
		delete[] interpolMask[k];
		delete[] distanceMap[k];
		delete[] ballMaxDist[k];
	}
	delete[] interpolImage;
	delete[] interpolMask;
	
	delete[] distanceMap;
	delete[] ballMaxDist;
	
	delete[] centroid_lungs;

	for (k = 0; k < height; k++) {
		delete[] croppedImage[k];
		delete[] croppedLungs[k];
		delete[] segmentCrop[k];
		delete[] maskThreshold[k];
	}
	delete[] croppedImage;
	delete[] croppedLungs;
	delete[] segmentCrop;
	delete[] maskThreshold;

	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
		delete[] maskLungs[k];
		delete[] maskLiver[k];
		delete[] maskHeart[k];
		delete[] ball[k];
	}
	delete[] imageData;
	delete[] maskLungs;
	delete[] maskLiver;
	delete[] maskHeart;
	delete[] ball;

	free(ctContainer);
	*/
	
	//========================= Translate as path ===============================
	
	/*
	Image_Data beforeTreatment = {
		ctContainer->dimensions[2],
		ctContainer->dimensions[0],
		ctContainer->dimensions[1],
		ctContainer->dataPointer,
		beforeOrigin,
		beforeSpacing,
		orientation
	};
	
	Vtk_File_Info* afterTreatmentContainer = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	afterTreatmentContainer->operation = copyFrom;

	loading_path = inputPath + "vtk/ct/Patient3_ct.vtk";
	readVtkFile(loading_path.c_str(), afterTreatmentContainer);
	Point3D afterOrigin = { afterTreatmentContainer->origin[0], afterTreatmentContainer->origin[1], afterTreatmentContainer->origin[2] };
	VoxelSpacing afterSpacing = { afterTreatmentContainer->spacing[0], afterTreatmentContainer->spacing[1], afterTreatmentContainer->spacing[0] };

	Image_Data afterTreatment = {
		afterTreatmentContainer->dimensions[2],
		afterTreatmentContainer->dimensions[0],
		afterTreatmentContainer->dimensions[1],
		afterTreatmentContainer->dataPointer,
		afterOrigin,
		afterSpacing,
		orientation
	};

	//Point3D seedP2 = { 5.169647, 71.255859, -656.131958 };//path_ordered_p2.csv
	//Point3D seedP3 = { 257, 255, 522 };

	//Point3D seedP2 = { 5.233856, 71.174133, -660.895447 };//first point of path_points_p2.csv, conference Jasna
	Point3D seedP2 = { 260.466217, 257.001923, 308.569214 };//first point of final_points_p2.csv, conference Jasna
	//seedP2 = getRealCoordFromImageCoord3D(seedP2, beforeOrigin, beforeSpacing, orientation);
	Point3D seedP3 = { 257, 255, 522 }; //chosen manually inside patient 3 aorta 
	//seedP3 = getRealCoordFromImageCoord3D(seedP3, afterOrigin, afterSpacing, orientation); // get correponding point for patient 3

	Point3D translate = { seedP3.x - seedP2.x, seedP3.y - seedP2.y, seedP3.z - seedP2.z };

	//Load path point
	FILE* path_file;
	//loading_path = inputPath + "inputsForPhDStudentConference/path_ordered_p2.csv";
	loading_path = "C:/Users/Konan Allaly/Documents/Presentation Jasna/final_path_p2.csv";
	if (fopen_s(&path_file, loading_path.c_str(), "r") != 0) {
		printf("Enable to open");
		return false;
	}

	FILE* translated_path;
	loading_path = "C:/Users/Konan Allaly/Documents/Presentation Jasna/final_translated_p2.csv";
	if (fopen_s(&translated_path, loading_path.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}

	dataType x = 0, y = 0, z = 0;
	dataType xt = 0, yt = 0, zt = 0;
	size_t count_point = 0;

	while (feof(path_file) == 0) {
		fscanf_s(path_file, "%f", &x);
		fscanf_s(path_file, ",");
		fscanf_s(path_file, "%f", &y);
		fscanf_s(path_file, ",");
		fscanf_s(path_file, "%f", &z);
		fscanf_s(path_file, "\n");
		count_point++;

		xt = x + translate.x;
		yt = y + translate.y;
		zt = z + translate.z;

		fprintf(translated_path, "%f,%f,%f\n", xt, yt, zt);
		//if (count_point == 1) {
		//	std::cout << "The first point is : (" << x << "," << y << "," << z << ")" << endl;
		//}

	}
	fclose(path_file);
	fclose(translated_path);
	free(afterTreatmentContainer);
	*/

	//======================== Automatic quantitative analysis ==================

	/*
	// Load PET
	Vtk_File_Info* petContainer = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	petContainer->operation = copyFrom;

	loading_path = inputPath + "vtk/pet/Patient3_pet.vtk";
	readVtkFile(loading_path.c_str(), petContainer);

	int height = petContainer->dimensions[2];
	int length = petContainer->dimensions[0];
	int width = petContainer->dimensions[1];
	int dim2d = length * width;
	std::cout << "PET image dim : " << petContainer->dimensions[0] << " x " << petContainer->dimensions[1] << " x " << petContainer->dimensions[2] << "" << std::endl;

	std::cout << "PET origin : (" << petContainer->origin[0] << ", " << petContainer->origin[1] << ", " << petContainer->origin[2] << ")" << std::endl;
	Point3D petOrigin = { petContainer->origin[0], petContainer->origin[1], petContainer->origin[2] };
	VoxelSpacing petSpacing = { petContainer->spacing[0], petContainer->spacing[1], petContainer->spacing[2] };
	std::cout << "PET spacing : (" << petContainer->spacing[0] << ", " << petContainer->spacing[1] << ", " << petContainer->spacing[2] << ")" << std::endl;

	dataType** liverShape = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		liverShape[k] = new dataType[dim2D]{ 0 };
	}

	// Load liver segment
	loading_path = inputPath + "inputsForPhDStudentConference/segment_liver_p3.raw";
	load3dArrayRAW<dataType>(liverShape, Length, Width, Height, loading_path.c_str(), false);

	//Find Liver centroid
	dataType* centroid = new dataType[3];
	centroidImage(liverShape, centroid, Height, Length, Width, 1.0);

	Point3D center_liver_img_coord = { centroid[0], centroid[1], centroid[2] };
	Point3D center_liver_real_coord = getRealCoordFromImageCoord3D(center_liver_img_coord, ctOrigin, ctSpacing, orientation);
	double radius_cetroid_liver = 20;
	
	////Draw ball inside the liver
	//dataType** ball = new dataType*[Height];
	//for (k = 0; k < Height; k++) {
	//	ball[k] = new dataType[dim2D]{0};
	//}

	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			Point3D current_point = { i, j, k };
	//			current_point = getRealCoordFromImageCoord3D(current_point, ctOrigin, ctSpacing, orientation);
	//			double dist = getPoint3DDistance(center_liver_real_coord, current_point);
	//			if (dist <= radius_cetroid_liver) {
	//				ball[k][x_new(i, j, Length)] = 1.0;
	//			}
	//		}
	//	}
	//}
	//storing_path = outputPath + "ball_liver_p2.raw";
	//store3dRawData<dataType>(ball, Length, Width, Height, storing_path.c_str());

	//Get the SUV mean around the liver centroid
	BoundingBox3D box_liver = findBoundingBox3D(center_liver_img_coord, Length, Width, Height, radius_cetroid_liver, 2 * radius_cetroid_liver);
	dataType mean_liver = 1.0, sum_val = 0;
	int count_point = 0;
	for (k = box_liver.k_min; k <= box_liver.k_max; k++) {
		for (i = box_liver.i_min; i <= box_liver.i_max; i++) {
			for (j = box_liver.j_min; j <= box_liver.j_max; j++) {
				Point3D current_point = { i, j, k };
				current_point = getRealCoordFromImageCoord3D(current_point, ctOrigin, ctSpacing, orientation);
				double d = getPoint3DDistance(center_liver_real_coord, current_point);
				if (d <= radius_cetroid_liver) {
					Point3D point_pet = getImageCoordFromRealCoord3D(current_point, petOrigin, petSpacing, orientation);
					dataType val_pet = petContainer->dataPointer[(size_t)point_pet.z][x_new((size_t)point_pet.x, (size_t)point_pet.y, length)];
					count_point++;
					sum_val += val_pet;
				}
			}
		}
	}
	if (count_point != 0) {
		mean_liver = sum_val / (dataType)count_point;
	}
	else {
		mean_liver = 1.0;
	}

	cout << "Mean inside liver is : " << mean_liver << endl;
	
	// Load the path points
	FILE* path_file;
	////loading_path = inputPath + "inputsForPhDStudentConference/path_ordered_p2.csv";
	////loading_path = inputPath + "inputsForPhDStudentConference/translated.csv";
	//loading_path = "C:/Users/Konan Allaly/Documents/Presentation Jasna/path_points_p2.csv";
	//loading_path = "C:/Users/Konan Allaly/Documents/Presentation Jasna/final_path_p2.csv";
	loading_path = "C:/Users/Konan Allaly/Documents/Presentation Jasna/final_translated_p2.csv";
	if (fopen_s(&path_file, loading_path.c_str(), "r") != 0) {
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

		//Saved without the spacing
		//TODO: Export the path considering the spacing or save in real world coordinates 
		//x *= ctSpacing.sx;
		//y *= ctSpacing.sy;
		//z *= ctSpacing.sz;
		Point3D current_point = { x, y, z };

		//update point dimension
		current_point = getRealCoordFromImageCoord3D(current_point, ctOrigin, intSpacing, orientation);
		
		path_points.push_back(current_point);
	}
	fclose(path_file);

	//File for the ratio max
	double search_radius = 25;
	//string saving_csv = outputPath + "ratio_max_5_before.csv";
	string saving_csv = outputPath + "ratio_max_25_after.csv";
	FILE* file;
	if (fopen_s(&file, saving_csv.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	fprintf(file, "x,y\n");

	//File for the ratio mean
	//saving_csv = outputPath + "ratio_mean_5_before.csv";
	saving_csv = outputPath + "ratio_mean_25_after.csv";
	FILE* rmean;
	if (fopen_s(&rmean, saving_csv.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	fprintf(rmean, "x,y\n");

	dataType** shape_aorta_ct = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		shape_aorta_ct[k] = new dataType[dim2D]{ 0 };
	}

	//Aorta Shape pet
	dataType** shape_aorta_pet = new dataType * [height];
	for (k = 0; k < height; k++) {
		shape_aorta_pet[k] = new dataType[dim2d]{ 0 };
	}

	BoundingBox3D box_path_point;
	size_t n = 0;
	dataType maxSuv = 0.0, sumSuv = 0.0, meanSuv = 0.0;
	dataType ratio_max_aorta_mean_liver = 0.0, ratio_mean_aorta_mean_liver = 0.0;
	vector<dataType> ratio_list;
	
	//vector<Point3D> points_max_list;

	Point3D point_ct;
	size_t i_pet, j_pet, k_pet, xd_pet;
	for (n = 0; n < path_points.size(); n++) {
		point_ct = getImageCoordFromRealCoord3D(path_points[n], ctOrigin, ctSpacing, orientation);
		box_path_point = findBoundingBox3D(point_ct, Length, Width, Height, search_radius, 2 * ctSpacing.sx);
		point_ct = getRealCoordFromImageCoord3D(point_ct, ctOrigin, ctSpacing, orientation);
		maxSuv = 0.0;
		sumSuv = 0.0;
		ratio_max_aorta_mean_liver = 0.0;
		count_point = 0;
		Point3D point_max = { 0.0, 0.0, 0.0 };
		for (k = box_path_point.k_min; k <= box_path_point.k_max; k++) {
			for (i = box_path_point.i_min; i <= box_path_point.i_max; i++) {
				for (j = box_path_point.j_min; j <= box_path_point.j_max; j++) {
					Point3D current_point = { i, j, k };
					current_point = getRealCoordFromImageCoord3D(current_point, ctOrigin, ctSpacing, orientation);
					double dist = getPoint3DDistance(current_point, point_ct);
					if (dist <= search_radius) {
						shape_aorta_ct[k][x_new(i, j, Length)] = 1.0; // save an visualize the shape in CT
						Point3D point_pet = getImageCoordFromRealCoord3D(current_point, petOrigin, petSpacing, orientation);
						i_pet = (size_t)point_pet.x;
						j_pet = (size_t)point_pet.y;
						k_pet = (size_t)point_pet.z;
						xd_pet = x_new(i_pet, j_pet, length);
						dataType val_pet = petContainer->dataPointer[k_pet][xd_pet];
						shape_aorta_pet[k_pet][xd_pet] = 1.0; // save an visualize the shape in PET
						sumSuv += val_pet;
						count_point++;
						if (maxSuv < val_pet) {
							maxSuv = val_pet;
						}
					}
				}
			}
		}
		meanSuv = sumSuv / (dataType)count_point;
		if (mean_liver != 0) {
			ratio_max_aorta_mean_liver = maxSuv / mean_liver;
			ratio_mean_aorta_mean_liver = meanSuv / mean_liver;
		}
		else {
			ratio_max_aorta_mean_liver = maxSuv;
			ratio_mean_aorta_mean_liver = meanSuv;
		}

		//ratio_list.push_back(ratio_aorta_liver);
		fprintf(file, "%d,%f\n", n, ratio_max_aorta_mean_liver);
		fprintf(rmean, "%d,%f\n", n, ratio_mean_aorta_mean_liver);

		//fprintf(rmean, "%d,%f\n", n, meanSuv);
		//dataType vis = 1.0;
		//fprintf(visual, "%d,%f\n", n, vis);
	}

	storing_path = outputPath + "shape_aorta_25_ct_p3.raw";
	store3dRawData<dataType>(shape_aorta_ct, Length, Width, Height, storing_path.c_str());

	storing_path = outputPath + "shape_aorta_25_pet_p3.raw";
	store3dRawData<dataType>(shape_aorta_pet, length, width, height, storing_path.c_str());

	fclose(file);
	fclose(rmean);
	//fclose(rpeak);

	////Find the max in the list
	//dataType max_ratio_list = 0;
	//for (n = 0; n < ratio_list.size(); n++) {
	//	if (max_ratio_list < ratio_list[n]) {
	//		max_ratio_list = ratio_list[n];
	//	}
	//}
	////Find points with max ratio and save them
	//saving_csv = outputPath + "ratio_max.csv";
	//FILE* file_ratio;
	//if (fopen_s(&file_ratio, saving_csv.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(file_ratio, "x,y,z\n");
	//for (n = 0; n < ratio_list.size(); n++) {
	//	if (ratio_list[n] >= max_ratio_list - 0.15 && ratio_list[n] <= max_ratio_list + 0.15) {
	//		fprintf(file_ratio, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	//	}
	//}
	//fclose(file_ratio);
	
	for (k = 0; k < height; k++) {
		delete[] shape_aorta_pet[k];
	}
	delete[] shape_aorta_pet;

	delete[] centroid;
	for (k = 0; k < Height; k++) {
		delete[] liverShape[k];
		//delete[] ball[k];
		delete[] shape_aorta_ct[k];
	}
	delete[] liverShape;
	//delete[] ball;
	delete[] shape_aorta_ct;

	free(petContainer);
	*/

	//======================== Filtering ========================================
	
	/*
	dataType** imageData = new dataType * [Height];
	dataType** filtered = new dataType * [Height];
	dataType** edgeDetector = new dataType * [Height];
	dataType** gradX = new dataType * [Height];
	dataType** gradY = new dataType * [Height];
	dataType** gradZ = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{ 0 };
		filtered[k] = new dataType[dim2D]{ 0 };
		edgeDetector[k] = new dataType[dim2D]{ 0 };
		gradX[k] = new dataType[dim2D]{ 0 };
		gradY[k] = new dataType[dim2D]{ 0 };
		gradZ[k] = new dataType[dim2D]{ 0 };
	}

	//Shift + copy
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			imageData[k][i] = ctContainer->dataPointer[k][i];
		}
	}
	rescaleNewRange(imageData, Length, Width, Height, 0.0, 1.0, minData, maxData);

	//storing_path = outputPath + "rescal.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	Image_Data inputImageData = { Height, Length, Width, imageData, ctOrigin, ctSpacing, orientation };

	size_t hauteur = (size_t)((ctSpacing.sz / ctSpacing.sx) * Height);
	VoxelSpacing intSpacing = { ctSpacing.sx, ctSpacing.sy, ctSpacing.sx };

	dataType** interpolate = new dataType * [hauteur];
	for (k = 0; k < hauteur; k++) {
		interpolate[k] = new dataType[dim2D]{ 0 };
	}

	Image_Data interpolateImageData = { hauteur, Length, Width, interpolate, ctOrigin, intSpacing, orientation };

	imageInterpolation3D(inputImageData, interpolateImageData, NEAREST_NEIGHBOR);

	storing_path = outputPath + "interpolate.raw";
	manageRAWFile3D<dataType>(interpolate, Length, Width, hauteur, storing_path.c_str(), STORE_DATA, false);

	dataType K = 1000;
	Filter_Parameters filter_parameters{
		1.5 * ctSpacing.sx,//tau
		ctSpacing.sx,//h
		0.0,//sigma
		K,//edge detector coefficient
		1.5,//omega_c
		1e-3,//tolerance
		1e-6,//epsilon2
		1e-6,//coef
		1,//p
		1,//number of time step
		100,// max solver iteration
	};

	geodesicMeanCurvatureTimeStep(interpolateImageData, filter_parameters);

	Image_Data filteredImageData = { Height, Length, Width, filtered, ctOrigin, ctSpacing, orientation };
	imageInterpolation3D(interpolateImageData, filteredImageData, NEAREST_NEIGHBOR);
	
	storing_path = outputPath + "filtered.raw";
	manageRAWFile3D<dataType>(filtered, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	////Compute edge detector
	//compute3dImageGradient(filtered, gradX, gradY, gradZ, Length, Width, Height, ctSpacing);
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			xd = x_new(i, j, Length);
	//			dataType norm_grad = gradX[k][xd] * gradX[k][xd] + gradY[k][xd] * gradY[k][xd] + gradZ[k][xd] * gradZ[k][xd];
	//			edgeDetector[k][xd] = gradientFunction(norm_grad, 20000);
	//		}
	//	}
	//}

	//storing_path = outputPath + "edge_detector.raw";
	//manageRAWFile3D<dataType>(edgeDetector, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
		delete[] filtered[k];
		delete[] edgeDetector[k];
		delete[] gradX[k];
		delete[] gradY[k];
		delete[] gradZ[k];
	}
	delete[] imageData;
	delete[] filtered;
	delete[] edgeDetector;
	delete[] gradX;
	delete[] gradY;
	delete[] gradZ;

	for (k = 0; k < hauteur; k++) {
		delete[] interpolate[k];
	}
	delete[] interpolate;

	free(ctContainer);
	*/
	
	//======================== Segment the Liver ================================================
	
	/*
	dataType** maskThreshold = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		maskThreshold[k] = new dataType[dim2D] {0};
	}
	if (maskThreshold == NULL)
		return false;

	// Copy
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			maskThreshold[k][i] = ctContainer->dataPointer[k][i];
			////remove blood
			//if (ctContainer->dataPointer[k][i] >= 30 && ctContainer->dataPointer[k][i] <= 45) {
			//	maskThreshold[k][i] = minData;
			//}
			//remove water
			if (ctContainer->dataPointer[k][i] == 0) {
				maskThreshold[k][i] = minData;
			}
			//remove fat
			if (ctContainer->dataPointer[k][i] >= -100 && ctContainer->dataPointer[k][i] <= -60) {
				maskThreshold[k][i] = minData;
			}
		}
	}
	
	//dataType thres_min = -30, thres_max = 190;
	dataType thres_min = -50, thres_max = 199;
	
	thresholding3dFunctionN(maskThreshold, Length, Width, Height, thres_min, thres_max, 0.0, 1.0);

	size_t n = 14;
	for (i = 0; i < n; i++) {
		erosion3D(maskThreshold, Length, Width, Height, 1.0, 0.0);
	}
	
	//storing_path = outputPath + "threshold_p4.raw";
	//manageRAWFile3D<dataType>(maskThreshold, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	int** labelArray = new int* [Height];
	bool** status = new bool* [Height];
	for (k = 0; k < Height; k++) {
		labelArray[k] = new int[dim2D] { 0 };
		status[k] = new bool[dim2D];
	}
	if (labelArray == NULL || status == NULL)
		return false;

	//Initialization
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			status[k][i] = false;
		}
	}

	int numberOfRegionCells = countNumberOfRegionsCells(maskThreshold, Length, Width, Height, 1.0);

	labelling3D(maskThreshold, labelArray, status, Length, Width, Height, 1.0);

	//Counting
	int* countingArray = new int[numberOfRegionCells] {0};
	if (countingArray == NULL)
		return false;

	//Get regions sizes
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (labelArray[k][i] > 0) {
				countingArray[labelArray[k][i]]++;
			}
		}
	}

	//Number of regions
	int numberOfRegions = 0;
	for (i = 0; i < numberOfRegionCells; i++) {
		if (countingArray[i] > 0) {
			numberOfRegions++;
		}
	}

	//Find the largest region
	size_t regionSize = 0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (countingArray[labelArray[k][i]] > regionSize) {
				regionSize = countingArray[labelArray[k][i]];
			}
		}
	}

	//Keep and save the largest region
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (countingArray[labelArray[k][i]] == regionSize) {
				maskThreshold[k][i] = 1.0;
			}
			else {
				maskThreshold[k][i] = 0.0;
			}
		}
	}

	storing_path = outputPath + "liver_p6.raw";
	store3dRawData<dataType>(maskThreshold, Length, Width, Height, storing_path.c_str());

	for (k = 0; k < Height; k++) {
		delete[] maskThreshold[k];
		delete[] labelArray[k];
		delete[] status[k];
	}
	delete[] maskThreshold;
	delete[] labelArray;
	delete[] status;
	free(ctContainer);
	*
	
	//======================== Refinement by generalized subjective surface ======================
	
	/*
	// We work with the cropped data to speed up the computation
	
	////Patient 1 : interpolated
	//size_t i_min = 90, i_max = 359;
	//size_t j_min = 90, j_max = 359;
	//size_t k_min = 307, k_max = 496;
	
	//Patient 2 : interpolated
	size_t i_min = 120, i_max = 375;
	size_t j_min = 160, j_max = 415;
	size_t k_min = 335, k_max = 494;

	////Patient 3 : interpolated
	//size_t i_min = 90, i_max = 345;
	//size_t j_min = 130, j_max = 385;
	//size_t k_min = 532, k_max = 711;

	////Patient 4 : interpolated
	//size_t i_min = 120, i_max = 375;
	//size_t j_min = 150, j_max = 405;
	//size_t k_min = 311, k_max = 518;

	////Patient 6: interpolated
	//size_t i_min = 120, i_max = 375;
	//size_t j_min = 180, j_max = 435;
	//size_t k_min = 703, k_max = 902;

	const size_t length = i_max - i_min + 1;
	const size_t width = j_max - j_min + 1;
	const size_t height = k_max - k_min + 1;
	const size_t dim2d = length * width;

	dataType** imageData = new dataType * [height];
	dataType** initialSegment = new dataType * [height];
	for (k = 0; k < height; k++) {
		imageData[k] = new dataType[dim2d]{0};
		initialSegment[k] = new dataType[dim2d]{0};
	}

	loading_path = inputPath + "inputsForPhDStudentConference/cropImage_p2.raw";
	load3dArrayRAW<dataType>(imageData, length, width, height, loading_path.c_str(), false);

	loading_path = inputPath + "inputsForPhDStudentConference/cropSegment_p2.raw";
	load3dArrayRAW<dataType>(initialSegment, length, width, height, loading_path.c_str(), false);

	VoxelSpacing ctInt = { ctContainer->spacing[0], ctContainer->spacing[1], ctContainer->spacing[0] };
	Point3D newOrigin = { i_min, j_min, k_min };
	newOrigin = getRealCoordFromImageCoord3D(newOrigin, ctOrigin, ctInt, orientation);
	std::cout << "New origin : (" << newOrigin.x << "," << newOrigin.y << "," << newOrigin.z << ")" << std::endl;

	Image_Data inputImage{
		height,
		length,
		width,
		imageData,
		newOrigin,
		ctInt,
		orientation
	};
	Image_Data inputSegment{
		height,
		length,
		width,
		initialSegment,
		newOrigin,
		ctInt,
		orientation
	};
	
	Point3D* center_seg = new Point3D[1];
	center_seg->x = 0.0, center_seg->y = 0.0; center_seg->z = 0.0;

	//Down sampling
	dataType down_factor = 4.0;
	const size_t h_down = (size_t)(height / down_factor);
	const size_t l_down = (size_t)(length / down_factor);
	const size_t w_down = (size_t)(width / down_factor);
	std::cout << "Downsampled image dim : " << l_down << " x " << w_down << " x " << h_down << "" << std::endl;
	dataType** imageDown = new dataType * [h_down];
	dataType** initialDown = new dataType * [h_down];
	for (k = 0; k < h_down; k++) {
		imageDown[k] = new dataType[l_down * w_down]{ 0 };
		initialDown[k] = new dataType[l_down * w_down]{ 0 };
	}

	VoxelSpacing downSampleSapcing = { ctInt.sx * down_factor, ctInt.sy * down_factor, ctInt.sz * down_factor };
	
	Image_Data downSampledImage{
	h_down,
	l_down,
	w_down,
	imageDown,
	newOrigin,
	downSampleSapcing,
	orientation
	};
	Image_Data downSampledSegment{
	h_down,
	l_down,
	w_down,
	initialDown,
	newOrigin,
	downSampleSapcing,
	orientation
	};
	std::cout << "Downsampled image spacing : (" << downSampleSapcing.sx << ", " << downSampleSapcing.sy << ", " << downSampleSapcing.sz << ")" << std::endl;

	imageInterpolation3D(inputImage, downSampledImage, NEAREST_NEIGHBOR);
	imageInterpolation3D(inputSegment, downSampledSegment, NEAREST_NEIGHBOR);

	//storing_path = outputPath + "downSeg.raw";
	//store3dRawData<dataType>(initialDown, l_down, w_down, h_down, storing_path.c_str());
	//storing_path = outputPath + "downImg.raw";
	//store3dRawData<dataType>(imageDown, l_down, w_down, h_down, storing_path.c_str());

	//Smoothing by heat equation, we don't need all the parameters
	Filter_Parameters smoothing_parameters{
		0.2, //tau
		1.0, //h
		0.0, //sigma
		0.0, //edge detector coeficient
		1.5, //omega
		1e-3,//tolerance
		0.0, // epsilon Evan-Spruck regularization
		0.0, // coef
		(size_t)1, //p
		(size_t)1, //time step number
		(size_t)100 //max number of solver iterations
	};

	Segmentation_Parameters segmentation_parameters{
		(size_t)30,  // Maximum number of Gauss - Seidel iterations
		2000,         // K : edge detector coef
		1e-6,         //epsilon : Evan-Spruck regularization
		(size_t)5000,   //number of current time step
		(size_t)5000,   //max number of time step
		(size_t)100,    //frequence au saving
		1e-4,         //segmentation tolerance
		down_factor,  //tau
		1.0,          //h
		1.5,          //omega
		1e-3,         //gauss_seidel tolerance
		0.01,         //convection coeff
		1.0           //diffusion coeff
	};
	
	unsigned char segmentPath[] = "C:/Users/Konan Allaly/Documents/Tests/output/segmentation/patient2/";
	//subsurfSegmentation(imageToBeSegmented, initialSegment, seg_params, smooth_params, center_seg, 1, segmentPath);
	generalizedSubsurfSegmentation(downSampledImage, initialDown, segmentation_parameters, smoothing_parameters, center_seg, 1, segmentPath);

	for (k = 0; k < h_down; k++) {
		delete[] imageDown[k];
		delete[] initialDown[k];
	}
	delete[] imageDown;
	delete[] initialDown;

	delete[] center_seg;
	for (k = 0; k < height; k++) {
		delete[] imageData[k];
		delete[] initialSegment[k];
	}
	delete[] imageData;
	delete[] initialSegment;
	*/

	//======================== Path extraction from multiple seeds ===============================
	
	/*
	dataType** imageData = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{ 0 };
	}
	if (imageData == NULL)
		return false;
	
	////Load filtered input image
	loading_path = inputPath + "raw/filtered/New/filtered_p2.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	////copy
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		imageData[k][i] = ctContainer->dataPointer[k][i];
	//	}
	//}
	//rescaleNewRange(imageData, Length, Width, Height, 0.0, 1.0, maxData, minData);
	//storing_path = outputPath + "recalled_p2.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	//string test_meta_data = outputPath + "metaDataTest_p6.txt";
	//FILE* test_description;
	//if (fopen_s(&test_description, test_meta_data.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(test_description, "This test is performed on Thursday 16th, January 2025\n");
	//fprintf(test_description, "The objective of the test is to extract paths inside the aorta\n");
	//fprintf(test_description, "The input data is patient number 6\n");
	//fprintf(test_description, "Due to the sensitivity of the approach to noise\n");
	//fprintf(test_description, "filtering by GMCF is performed on the input image\n");
	//fprintf(test_description, "\nPatient 6: origin = (%f, %f, %f)\n", ctOrigin.x, ctOrigin.y, ctOrigin.z);
	//fprintf(test_description, "Voxel sizes = (%f, %f, %f)\n", ctSpacing.sx, ctSpacing.sy, ctSpacing.sz);
	//fprintf(test_description, "Dimension : Length = %d , Width = %d and Height = %d\n", Length, Width, Height);
	//fprintf(test_description, "\nThe path extraction is performed into 2 steps using 3 input points\n");
	//fprintf(test_description, "The current implementation of fast marching use isotropic voxel size\n");
	//const size_t height = (size_t)((ctSpacing.sz / ctSpacing.sx) * Height);
	//std::cout << "Interpolated Height = " << height << std::endl;
	//fprintf(test_description, "Interpolated image dimension: Length = %d, Width = %d, Height = %d\n", Length, Width, height);
	
	const size_t height = (size_t)((ctSpacing.sz / ctSpacing.sx) * Height);
	dataType** potential = new dataType * [height];
	dataType** actionField = new dataType * [height];
	dataType** imageInterpolated = new dataType * [height];
	for (k = 0; k < height; k++) {
		potential[k] = new dataType[dim2D]{ 0 };
		actionField[k] = new dataType[dim2D]{ 0 };
		imageInterpolated[k] = new dataType[dim2D]{ 0 };
	}
	if (potential == NULL || actionField == NULL || imageInterpolated == NULL)
		return false;

	Image_Data ctImageData
	{
		Height,
		Length,
		Width,
		imageData,
		ctOrigin,
		ctSpacing,
		orientation
	};
	
	VoxelSpacing intSpacing = { ctSpacing.sx, ctSpacing.sy, ctSpacing.sx };
	Image_Data IntCtImageData
	{
		height,
		Length,
		Width,
		imageInterpolated,
		ctOrigin,
		intSpacing,
		orientation
	};
	
	//fprintf(test_description, "We perform interpolation on the zth direction to have cubic voxels\n");
	//fprintf(test_description, "Interpolated image Dimension : Length = %d , Width = %d and Height = %d\n", Length, Width, height);
	
	//imageInterpolation3D(ctImageData, IntCtImageData, NEAREST_NEIGHBOR);

	double radius = 3.0;
	Potential_Parameters parameters{
		1000, //edge detector coefficient
		0.15, //threshold
		0.005,//epsilon
		radius
	};
	
	//fprintf(test_description, "\nThe potential function is computed only on the initial points\n");
	//fprintf(test_description, "To avoid noise or point location effects, we take the average value in spherical (3 mm of radius) volume around the initial point\n");
	//fprintf(test_description, "Potential function parameters : \n");
	//fprintf(test_description, "radius = %f define the volume where to compute the seed value\n", radius);
	//fprintf(test_description, "epsilon = %f is the path smoothing parameter\n", parameters.eps);
	//fprintf(test_description, "thres = %f is the threshold for distance map\n", parameters.eps);
	
	//Points in original image coordinates
	////Patient 6
	//Point3D seed1 = {250, 300, 255};
	//Point3D seed2 = {277, 351, 459};
	//Point3D seed3 = {250, 295, 459};

	////Patient 5
	//Point3D seed1 = { 268, 241, 476 };
	//Point3D seed2 = { 283, 282, 646 };
	//Point3D seed3 = { 243, 216, 646 };

	////Patient 4
	//Point3D seed1 = { 284, 230, 132 };
	//Point3D seed2 = { 259, 313, 236 };
	//Point3D seed3 = { 279, 242, 236 };

	////Patient 3
	//Point3D seed1 = { 267, 232, 111 };
	//Point3D seed2 = { 290, 288, 232 };
	//Point3D seed3 = { 262, 227, 232 };

	////Patient 2
	//Point3D seed1 = { 260, 256, 245 };
	//Point3D seed2 = { 290, 315, 363 };
	//Point3D seed3 = { 249, 253, 363 };

	////Patient 1
	//Point3D seed1 = { 262, 255, 147 };
	//Point3D seed2 = { 298, 315, 264 };
	//Point3D seed3 = { 255, 258, 264 };

	//fprintf(test_description, "\nThe points used to extract the current path are set manually\n");
	//fprintf(test_description, "The points are defined in original image then interpolated\n");
	//fprintf(test_description, "Point 1 : (%f, %f, %f)\n", seed1.x, seed1.y, seed1.z);
	//fprintf(test_description, "Point 2 : (%f, %f, %f)\n", seed2.x, seed2.y, seed2.z);
	//fprintf(test_description, "Point 3 : (%f, %f, %f)\n", seed3.x, seed3.y, seed3.z);
	//fclose(test_description);

	////To real world coordinate
	//seed1 = getRealCoordFromImageCoord3D(seed1, ctOrigin, ctSpacing, orientation);
	//seed2 = getRealCoordFromImageCoord3D(seed2, ctOrigin, ctSpacing, orientation);
	//seed3 = getRealCoordFromImageCoord3D(seed3, ctOrigin, ctSpacing, orientation);

	//seed1 = getImageCoordFromRealCoord3D(seed1, ctOrigin, intSpacing, orientation);
	//seed2 = getImageCoordFromRealCoord3D(seed2, ctOrigin, intSpacing, orientation);
	//seed3 = getImageCoordFromRealCoord3D(seed3, ctOrigin, intSpacing, orientation);

	
	//To interpolated image coordinates
	Point3D* seedPoints = new Point3D[2];
	
	//Point3D seed1 = { 260, 256, 143 };//p1
	Point3D seed1 = { 259, 255, 244 };//p2, 03/03
	//Point3D seed1 = { 258, 254, 247 };//p2
	//Point3D seed1 = { 263, 232, 112 };//p3
	//Point3D seed1 = { 282, 231, 133 };//p4
	//Point3D seed1 = { 263, 243, 472 };//p5
	//Point3D seed1 = { 250, 299, 256 };//p6
	seed1 = getRealCoordFromImageCoord3D(seed1, ctOrigin, ctSpacing, orientation);
	seed1 = getImageCoordFromRealCoord3D(seed1, ctOrigin, intSpacing, orientation);

	//Point3D seed2 = { 292, 316, 245 };//p1
	Point3D seed2 = { 289, 317, 343 };//p2, 03/03
	//Point3D seed2 = { 287, 317, 352 };//p2
	//Point3D seed2 = { 286, 294, 195 };//p3
	//Point3D seed2 = { 254, 313, 224 };//p4
	//Point3D seed2 = { 263, 300, 602 };//p5
	//Point3D seed2 = { 281, 349, 394 };//p6
	seed2 = getRealCoordFromImageCoord3D(seed2, ctOrigin, ctSpacing, orientation);
	seed2 = getImageCoordFromRealCoord3D(seed2, ctOrigin, intSpacing, orientation);

	//Point3D seed3 = { 279, 275, 276 };//p1
	Point3D seed3 = { 276, 281, 371 };//p2, 03/03
	//Point3D seed3 = { 270, 278, 373 };//p2
	//Point3D seed3 = { 285, 254, 240 };//p3
	//Point3D seed3 = { 271, 283, 247 };//p4
	//Point3D seed3 = { 278, 249, 655 };//p5
	//Point3D seed3 = { 272, 318, 480 };//p6
	seed3 = getRealCoordFromImageCoord3D(seed3, ctOrigin, ctSpacing, orientation);
	seed3 = getImageCoordFromRealCoord3D(seed3, ctOrigin, intSpacing, orientation);

	//Point3D seed4 = { 253, 252, 262 };//p1
	Point3D seed4 = { 249, 252, 362 };//p2, 03/03
	//Point3D seed4 = { 248, 254, 360 };//p2
	//Point3D seed4 = { 255, 220, 224 };//p3
	//Point3D seed4 = { 273, 237, 234 };//p4
	//Point3D seed4 = { 246, 214, 646 };//p5
	//Point3D seed4 = { 250, 295, 452 };//p6
	seed4 = getRealCoordFromImageCoord3D(seed4, ctOrigin, ctSpacing, orientation);
	seed4 = getImageCoordFromRealCoord3D(seed4, ctOrigin, intSpacing, orientation);

	//Point3D seed5 = { 265, 254, 244 };//p1
	Point3D seed5 = { 258, 252, 344 };//p2, 03/03
	//Point3D seed5 = { 285, 251, 335 };//p2
	//Point3D seed5 = { 283, 215, 201 };//p3
	//Point3D seed5 = { 309, 233, 209 };//p4
	//Point3D seed5 = { 238, 227, 611 };//p5
	//Point3D seed5 = { 294, 294, 416 };//p6
	seed5 = getRealCoordFromImageCoord3D(seed5, ctOrigin, ctSpacing, orientation);
	seed5 = getImageCoordFromRealCoord3D(seed5, ctOrigin, intSpacing, orientation);

	////Branches ---> tested just for patient 1 and 2
	//Point3D seed6 = { 208, 246, 88 };//p1
	Point3D seed6 = { 235, 254, 225 };//p2
	seed6 = getRealCoordFromImageCoord3D(seed6, ctOrigin, ctSpacing, orientation);
	seed6 = getImageCoordFromRealCoord3D(seed6, ctOrigin, intSpacing, orientation);

	//Point3D seed7 = { 302, 244, 88 };//p1
	Point3D seed7 = { 279, 266, 225 };//p2
	seed7 = getRealCoordFromImageCoord3D(seed7, ctOrigin, ctSpacing, orientation);
	seed7 = getImageCoordFromRealCoord3D(seed7, ctOrigin, intSpacing, orientation);

	Image_Data inputImageData = { height, Length, Width, imageInterpolated, ctOrigin, intSpacing, orientation };

	//compute3DPotential(inputImageData, potential, seed1, parameters);
	storing_path = outputPath + "potential_p2.raw";
	manageRAWFile3D<dataType>(potential, Length, Width, height, storing_path.c_str(), LOAD_DATA, false);

	//fastMarching3D_N(actionField, potential, Length, Width, height, seed3);
	//fastMarching3dWithSpacing(inputImageData, actionField, potential, seed1, ctSpacing);
	//storing_path = outputPath + "action.raw";
	//manageRAWFile3D<dataType>(actionField, Length, Width, height, storing_path.c_str(), STORE_DATA, false);

	
	//============ first segment =========
	vector<Point3D> path_points;
	Image_Data actionDataStr = { height, Length, Width, actionField, ctOrigin, intSpacing, orientation };
	dataType tau = 0.8;
	dataType tolerance = ctSpacing.sx;

	//copy points to file
	string saving_csv = outputPath + "path_point_with_branches.csv";
	FILE* file_multiple;
	if (fopen_s(&file_multiple, saving_csv.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}

	seedPoints[0] = seed1;
	seedPoints[1] = seed2;
	partialFrontPropagation(actionField, potential, Length, Width, height, seedPoints);
	
	shortestPath3D(actionDataStr, seedPoints, tau, tolerance, path_points);
	for (int n = path_points.size() - 1; n > -1; n--) {
		path_points[n] = getRealCoordFromImageCoord3D(path_points[n], ctOrigin, intSpacing, orientation);
		fprintf(file_multiple, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	}
	while (path_points.size() > 0) {
		path_points.pop_back();
	}

	
	//======== branches aorta =========
	
	seedPoints[1] = seed6;
	shortestPath3D(actionDataStr, seedPoints, tau, tolerance, path_points);
	for (int n = path_points.size() - 1; n > -1; n--) {
		path_points[n] = getRealCoordFromImageCoord3D(path_points[n], ctOrigin, intSpacing, orientation);
		fprintf(file_multiple, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	}
	while (path_points.size() > 0) {
		path_points.pop_back();
	}

	seedPoints[1] = seed7;
	shortestPath3D(actionDataStr, seedPoints, tau, tolerance, path_points);
	for (int n = path_points.size() - 1; n > -1; n--) {
		path_points[n] = getRealCoordFromImageCoord3D(path_points[n], ctOrigin, intSpacing, orientation);
		fprintf(file_multiple, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	}
	while (path_points.size() > 0) {
		path_points.pop_back();
	}
	
	////============ second segment ==========
	

	seedPoints[0] = seed2;
	seedPoints[1] = seed3;
	partialFrontPropagation(actionField, potential, Length, Width, height, seedPoints);
	shortestPath3D(actionDataStr, seedPoints, tau, tolerance, path_points);
	for (int n = path_points.size() - 1; n > -1; n--) {
		path_points[n] = getRealCoordFromImageCoord3D(path_points[n], ctOrigin, intSpacing, orientation);
		fprintf(file_multiple, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	}
	while (path_points.size() > 0) {
		path_points.pop_back();
	}

	////=========== Third segment ==============

	seedPoints[0] = seed3;
	seedPoints[1] = seed4;
	partialFrontPropagation(actionField, potential, Length, Width, height, seedPoints);
	shortestPath3D(actionDataStr, seedPoints, tau, tolerance, path_points);

	for (int n = path_points.size() - 1; n > -1; n--) {
		path_points[n] = getRealCoordFromImageCoord3D(path_points[n], ctOrigin, intSpacing, orientation);
		fprintf(file_multiple, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	}
	while (path_points.size() > 0) {
		path_points.pop_back();
	}

	////=========== Fouth segment ==============

	seedPoints[0] = seed4;
	seedPoints[1] = seed5;
	partialFrontPropagation(actionField, potential, Length, Width, height, seedPoints);
	shortestPath3D(actionDataStr, seedPoints, tau, tolerance, path_points);

	for (int n = path_points.size() - 1; n > -1; n--) {
		path_points[n] = getRealCoordFromImageCoord3D(path_points[n], ctOrigin, intSpacing, orientation);
		fprintf(file_multiple, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	}
	while (path_points.size() > 0) {
		path_points.pop_back();
	}

	fclose(file_multiple);
	
	//findPathTwoSteps(IntCtImageData, seedPoints, parameters);
	
	
	delete[] seedPoints;
	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
	}
	delete[] imageData;

	for (k = 0; k < height; k++) {
		delete[] potential[k];
		delete[] actionField[k];
		delete[] imageInterpolated[k];
	}
	delete[] potential;
	delete[] actionField;
	delete[] imageInterpolated;
	
	free(ctContainer);
	*/
	
	//======================== Segment lungs connected to trachea ==============

	/*
	dataType** maskThreshold = new dataType * [Height];
	dataType** imageData = new dataType * [Height];
	dataType** maskLungs = new dataType * [Height];
	int** segment = new int * [Height];
	bool** status = new bool * [Height];
	for (k = 0; k < Height; k++) {
		maskThreshold[k] = new dataType[dim2D]{ 0 };
		imageData[k] = new dataType[dim2D]{ 0 };
		maskLungs[k] = new dataType[dim2D]{ 0 };
		segment[k] = new int[dim2D]{ 0 };
		status[k] = new bool[dim2D]{ false };
	}

	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			//if (minData <= -1024) {
			//	maskThreshold[k][i] = -1024;
			//}
			//else {
			//	maskThreshold[k][i] = ctContainer->dataPointer[k][i];
			//}
			maskThreshold[k][i] = ctContainer->dataPointer[k][i];
			//imageData[k][i] = ctContainer->dataPointer[k][i];
		}
	}

	dataType thres = -150;
	thresholding3dFunctionN(maskThreshold, Length, Width, Height, thres, thres, 0.0, 1.0);
	//iterativeThreshold(maskThreshold, Length, Width, Height, 0.0, 1.0);
	
	//dilatation3D(maskThreshold, Length, Width, Height, 1.0, 0.0);
	erosion3D(maskThreshold, Length, Width, Height, 1.0, 0.0);
	//erosion3D(maskThreshold, Length, Width, Height, 1.0, 0.0);

	labelling3D(maskThreshold, segment, status, Length, Width, Height, 1.0);

	int numberOfRegionCells = countNumberOfRegionsCells(maskThreshold, Length, Width, Height, 1.0);

	//Counting
	int* countingArray = new int[numberOfRegionCells] {0};
	if (countingArray == NULL)
		return false;

	//Get regions sizes
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (segment[k][i] > 0) {
				countingArray[segment[k][i]]++;
			}
		}
	}

	//Number of regions
	int numberOfRegions = 0;
	for (i = 0; i < numberOfRegionCells; i++) {
		if (countingArray[i] > 0) {
			numberOfRegions++;
		}
	}
	//std::cout << "We have " << numberOfRegions << " regions" << std::endl;

	//Find the largest region
	size_t max1 = 0;
	size_t max2 = 0;
	size_t max3 = 0;
	size_t max4 = 0;
	size_t max5 = 0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (countingArray[segment[k][i]] > max1) {
				max1 = countingArray[segment[k][i]];
			}
			if (countingArray[segment[k][i]] < max1 && countingArray[segment[k][i]] > max2){
				max2 = countingArray[segment[k][i]];
			}
			if (countingArray[segment[k][i]] < max2 && countingArray[segment[k][i]] > max3) {
				max3 = countingArray[segment[k][i]];
			}
			if (countingArray[segment[k][i]] < max3 && countingArray[segment[k][i]] > max4) {
				max4 = countingArray[segment[k][i]];
			}
			if (countingArray[segment[k][i]] < max4 && countingArray[segment[k][i]] > max5) {
				max5 = countingArray[segment[k][i]];
			}
		}
	}

	//size_t elmt = (size_t)(max2 / max3);
	//std::cout << "Ratio regions : " << elmt << std::endl;

	//Keep region
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (countingArray[segment[k][i]] == max3) {
				maskLungs[k][i] = 1.0;
			}
			else {
				maskLungs[k][i] = 0.0;
			}
			//if (elmt >= 2) {
			//	if (countingArray[segment[k][i]] == max2) {
			//		maskLungs[k][i] = 1.0;
			//	}
			//	else {
			//		maskLungs[k][i] = 0.0;
			//	}
			//}
			//else {
			//	if (countingArray[segment[k][i]] == max2 || countingArray[segment[k][i]] == max3) {
			//		maskLungs[k][i] = 1.0;
			//	}
			//	else {
			//		maskLungs[k][i] = 0.0;
			//	}
			//}
		}
	}
	
	//Keep the same number of erosion versus dilatation
	dilatation3D(maskLungs, Length, Width, Height, 1.0, 0.0);

	delete[] countingArray;
	
	storing_path = outputPath + "lungs_p3.raw";
	manageRAWFile3D<dataType>(maskLungs, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	//Find Lungs centroid
	dataType* centroid_lungs = new dataType[3]{ 0 };
	centroidImage(maskLungs, centroid_lungs, Height, Length, Width, 0.0);
	Point3D pCentroid = { centroid_lungs[0], centroid_lungs[1] , centroid_lungs[2] };
	std::cout << "Lungs centroid : (" << pCentroid.x << "," << pCentroid.y << "," << pCentroid.z << ")" << std::endl;
	//Point3D prCentroid = getRealCoordFromImageCoord3D(pCentroid, ctOrigin, ctSpacing, orientation);
	//double radius = 10.0;

	//Find the aorta in the lungs centroid slice
	size_t k_slice = (size_t)pCentroid.z;
	Point2D centroid_seed = {pCentroid.x, pCentroid.y};

	dataType* imageSlice = new dataType[dim2D]{ 0 };
	dataType* gradientX = new dataType[dim2D]{ 0 };
	dataType* gradientY = new dataType[dim2D]{ 0 };
	dataType* threshold = new dataType[dim2D]{ 0 };
	dataType* foundCirclePtr = new dataType[dim2D]{ 0 };

	//get the slice
	copyDataToAnother2dArray(ctContainer->dataPointer[k_slice], imageSlice, Length, Width);

	//Rescale to range 0-1
	rescaleNewRange2D(imageSlice, Length, Width, 0.0, 1.0);

	dataType norm_of_gradient = 0.0, edge_value;

	PixelSpacing spacing = { ctSpacing.sx, ctSpacing.sy };
	Point2D origin = { 0.0, 0.0 };
	OrientationMatrix2D orientation2D;
	orientation2D.v1 = { 1.0, 0.0 };
	orientation2D.v2 = { 0.0, 1.0 };

	HoughParameters h_parameters =
	{
		3.0,   //minimal radius
		15.0,  //maximal radius
		0.5,   //radius step
		0.5,   //epsilon
		100,    //offset 
		1000,  //K edge detector coefficient
		0.15,  //threshold edge detector
	};
	Image_Data2D IMAGE_HOUGH
	{
		Length,
		Width,
		threshold,
		origin,
		spacing,
		orientation2D
	};

	//Slice filtering
	Filter_Parameters filter_parameters{
		0.25 * spacing.sx * spacing.sx, //tau
		spacing.sx, //h
		0.0, //sigma
		0.0, //edge detector coef (not used)
		1.5, //omega c
		1e-3, //tolerance
		1e-6, //epsilon2
		1e-6, //coef
		1, //p Neumann boundary condition
		5, //time step
		1000, //max solver iteration
	};
	Image_Data2D IMAGE{
		Length,
		Width,
		imageSlice,
		origin,
		spacing,
		orientation2D
	};
	heatImplicit2dScheme(IMAGE, filter_parameters);

	storing_path = outputPath + "image_slice.raw";
	manageRAWFile2D<dataType>(imageSlice, Length, Width, storing_path.c_str(), STORE_DATA, false);

	//Compute gradient
	computeImageGradient(IMAGE, gradientX, gradientY);

	//Compute edge detector and threshold
	norm_of_gradient = 0.0;
	edge_value = 0.0;
	for (i = 0; i < dim2D; i++) {
		norm_of_gradient = sqrt(gradientX[i] * gradientX[i] + gradientY[i] * gradientY[i]);
		edge_value = gradientFunction(norm_of_gradient, h_parameters.K);
		if (edge_value <= h_parameters.thres) {
			threshold[i] = 1.0;
		}
	}

	Point2D point_found = localCircleDetection(IMAGE_HOUGH, foundCirclePtr, centroid_seed, h_parameters);
	std::cout << "Found center : (" << point_found.x << "," << point_found.y << ")" << std::endl;

	storing_path = outputPath + "threshold.raw";
	manageRAWFile2D<dataType>(threshold, Length, Width, storing_path.c_str(), STORE_DATA, false);
	storing_path = outputPath + "found_circle.raw";
	manageRAWFile2D<dataType>(foundCirclePtr, Length, Width, storing_path.c_str(), STORE_DATA, false);

	delete[] imageSlice;
	delete[] gradientX;
	delete[] gradientY;
	delete[] threshold;
	delete[] foundCirclePtr;

	////Draw ball for visualization
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			Point3D current_point = { i, j, k };
	//			current_point = getRealCoordFromImageCoord3D(current_point, ctOrigin, ctSpacing, orientation);
	//			double pDistance = getPoint3DDistance(current_point, prCentroid);
	//			if (pDistance <= radius) {
	//				maskTrachea[k][x_new(i, j, Length)] = 1.0;
	//			}
	//		}
	//	}
	//}
	//storing_path = outputPath + "ball_centroid_lungs_p1.raw";
	//manageRAWFile3D<dataType>(maskTrachea, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	//dataType* imageSlice = new dataType[dim2D]{ 0 };
	//size_t k_slice = (size_t)pCentroid.z;

	//copyDataToAnother2dArray(imageData[k_slice], imageSlice, Length, Width);
	//storing_path = outputPath + "slice_trachea_p1.raw";
	//manageRAWFile2D(imageSlice, Length, Width, storing_path.c_str(), STORE_DATA, false);
	
	vector<Point3D> path_max_distance;
	Image_Data ctImageData = {
		Height,
		Length,
		Width,
		maskLungs,
		ctOrigin,
		ctSpacing,
		orientation
	};
	//rouyTourinDistanceMap(ctImageData, difference, 1.0, 0.5, 0.4);

	storing_path = outputPath + "distance_map_lungs_p1.raw";
	//manageRAWFile3D<dataType>(difference, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	dataType max_per_slice = 0;
	Point3D point_max = { 0, 0, 0 };
	for (k = 0; k < Height; k++) {
		max_per_slice = 0.0;
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				xd = x_new(i, j, Length);
				if (difference[k][xd] > max_per_slice) {
					max_per_slice = difference[k][xd];
					point_max.x = i;
					point_max.y = j;
					point_max.z = k;
				}
			}
		}
		if (max_per_slice != 0) {
			path_max_distance.push_back(point_max);
		}
	}
	
	string saving_csv = outputPath + "point_max_p1.csv";
	FILE* file_path;
	if (fopen_s(&file_path, saving_csv.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}

	//copy points to file
	int n = 0;
	for (n = path_max_distance.size() - 1; n > -1; n--) {
		path_max_distance[n] = getRealCoordFromImageCoord3D(path_max_distance[n], ctOrigin, ctSpacing, orientation);
		fprintf(file_path, "%f,%f,%f\n", path_max_distance[n].x, path_max_distance[n].y, path_max_distance[n].z);
	}
	fclose(file_path);

	//======================= get the lungs without the trachea ==============

	//copy data
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			//copy the mask threshsold for next processings
			maskThreshold[k][i] = maskLungs[k][i];
			//maskTrachea[k][i] = maskLungs[k][i];
			//initialize for next processings
			maskLungs[k][i] = 0.0;
		}
	}

	//erosion to disconnet the lungs to the trachea
	erosion3D(maskThreshold, Length, Width, Height, 1.0, 0.0);
	erosion3D(maskThreshold, Length, Width, Height, 1.0, 0.0);

	//Labeling to find the lungs without the trachea
	//initialize the arrays
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			segment[k][i] = 0;
			status[k][i] = false;
		}
	}
	
	labelling3D(maskThreshold, segment, status, Length, Width, Height, 1.0);
	numberOfRegionCells = countNumberOfRegionsCells(maskThreshold, Length, Width, Height, 1.0);

	//Counting
	int* countingArraySecond = new int[numberOfRegionCells] {0};
	if (countingArraySecond == NULL)
		return false;

	//Get regions sizes
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (segment[k][i] > 0) {
				countingArraySecond[segment[k][i]]++;
			}
		}
	}

	//Number of regions
	numberOfRegions = 0;
	for (i = 0; i < numberOfRegionCells; i++) {
		if (countingArraySecond[i] > 0) {
			numberOfRegions++;
		}
	}

	//Find the largest region
	max1 = 0;
	max2 = 0;
	max3 = 0;
	max4 = 0;
	max5 = 0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (countingArraySecond[segment[k][i]] > max1) {
				max1 = countingArraySecond[segment[k][i]];
			}
			if (countingArraySecond[segment[k][i]] < max1 && countingArraySecond[segment[k][i]] > max2) {
				max2 = countingArraySecond[segment[k][i]];
			}
			if (countingArraySecond[segment[k][i]] < max2 && countingArraySecond[segment[k][i]] > max3) {
				max3 = countingArraySecond[segment[k][i]];
			}
			if (countingArraySecond[segment[k][i]] < max3 && countingArraySecond[segment[k][i]] > max4) {
				max4 = countingArray[segment[k][i]];
			}
			if (countingArraySecond[segment[k][i]] < max4 && countingArraySecond[segment[k][i]] > max5) {
				max5 = countingArraySecond[segment[k][i]];
			}
		}
	}

	//Keep region
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (countingArraySecond[segment[k][i]] == max1 || countingArraySecond[segment[k][i]] == max2) {
				maskLungs[k][i] = 1.0;
			}
			else {
				maskLungs[k][i] = 0.0;
			}
		}
	}

	//Keep the same number of erosions
	dilatation3D(maskLungs, Length, Width, Height, 1.0, 0.0);
	dilatation3D(maskLungs, Length, Width, Height, 1.0, 0.0);

	delete[] countingArraySecond;

	storing_path = outputPath + "p3/lungs_without_trachea_p3.raw";
	manageRAWFile3D<dataType>(maskLungs, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	//===== Apply labeling on the difference to extract the trachea ===============

	//Difference
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			//difference[k][i] = maskTrachea[k][i] - maskLungs[k][i];
			//maskThreshold[k][i] = maskTrachea[k][i] - maskLungs[k][i];
		}
	}
	
	storing_path = outputPath + "p3/difference_p3.raw";
	//manageRAWFile3D<dataType>(difference, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	//initialize the arrays
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			//maskTrachea[k][i] = 0.0;
			segment[k][i] = 0;
			status[k][i] = false;
		}
	}

	labelling3D(maskThreshold, segment, status, Length, Width, Height, 1.0);
	numberOfRegionCells = countNumberOfRegionsCells(maskThreshold, Length, Width, Height, 1.0);

	//Counting
	int* countingArrayThird = new int[numberOfRegionCells] {0};
	if (countingArrayThird == NULL)
		return false;

	//Get regions sizes
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (segment[k][i] > 0) {
				countingArrayThird[segment[k][i]]++;
			}
		}
	}

	//Number of regions
	numberOfRegions = 0;
	for (i = 0; i < numberOfRegionCells; i++) {
		if (countingArrayThird[i] > 0) {
			numberOfRegions++;
		}
	}

	//Find the largest region
	max1 = 0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (countingArrayThird[segment[k][i]] > max1) {
				max1 = countingArrayThird[segment[k][i]];
			}
		}
	}

	//Keep region
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (countingArrayThird[segment[k][i]] == max1) {
				//maskTrachea[k][i] = 1.0;
			}
			else {
				//maskTrachea[k][i] = 0.0;
			}
		}
	}
	
	//dilatation3D(maskTrachea, Length, Width, Height, 1.0, 0.0);

	storing_path = outputPath + "p3/trachea_p3.raw";
	//manageRAWFile3D<dataType>(maskTrachea, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	delete[] countingArrayThird;

	for (k = 0; k < Height; k++) {
		delete[] segment[k];
		delete[] status[k];
	}
	delete[] segment;
	delete[] status;

	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
		delete[] maskThreshold[k];
		delete[] maskLungs[k];
		delete[] segment[k];
		delete[] status[k];
	}
	delete[] imageData;
	delete[] maskThreshold;
	delete[] maskLungs;
	delete[] segment;
	delete[] status;
	free(ctContainer);
	*/
	
	//======================== 2D Hough transform on Slice ===========================
	
	/*
	//const size_t Length = 512, Width = 512;
	//const size_t dim2D = Length * Width;

	dataType** imageData = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{ 0 };
	}

	loading_path = inputPath + "raw/liver/liver_p1.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	dataType* imageSlice = new dataType[dim2D]{ 0 };
	bool* statusPixel = new bool[dim2D] { false };

	////Lungs Centroid
	//Point3D centroid_p1 = { 264, 289, 249 };
	//Point3D centroid_p2 = { 257, 289, 348 };
	//Point3D centroid_p3 = { 252, 280, 220 };
	//Point3D centroid_p4 = { 237, 272, 225 };
	//Point3D centroid_p5 = { 255, 266, 623 };
	//Point3D centroid_p6 = { 242, 328, 441 };
	//Point3D slice_liver_p1 = { 214, 273, 202 };

	//size_t k_slice = (size_t)centroid_p1.z;
	////size_t k_slice = (size_t)slice_liver_p1.z;

	//Point2D centroid_seed = {centroid_p1.x, centroid_p1.y};
	////Point2D centroid_seed = { slice_liver_p1.x, slice_liver_p1.y };

	PixelSpacing spacing = { ctSpacing.sx, ctSpacing.sy };
	Point2D origin = { 0.0, 0.0 };
	OrientationMatrix2D orientation2D;
	orientation2D.v1 = { 1.0, 0.0 };
	orientation2D.v2 = { 0.0, 1.0 };

	Image_Data2D IMAGE{
		Length,
		Width,
		imageSlice,
		origin,
		spacing,
		orientation2D
	};

	//filtering parameters
	Filter_Parameters filter_parameters{
		0.25 * spacing.sx * spacing.sx, //tau
		spacing.sx, //h
		0.0, //sigma
		0.0, //edge detector coef (not used)
		1.5, //omega c
		1e-3, //tolerance
		1e-6, //epsilon2
		1e-6, //coef
		1, //p Neumann boundary condition
		5, //time step
		1000, //max solver iteration
	};

	HoughParameters h_parameters =
	{
		5.0,   //minimal radius
		15.0,  //maximal radius
		0.5,   //radius step
		0.5,   //epsilon
		100,    //offset 
		1000,  //K edge detector coefficient
		0.15,  //threshold edge detector
	};

	dataType* gradientX = new dataType[dim2D]{ 0 };
	dataType* gradientY = new dataType[dim2D]{ 0 };
	dataType* threshold = new dataType[dim2D]{ 0 };
	dataType norm_of_gradient = 0.0, edge_value;

	Image_Data2D IMAGE_HOUGH
	{
		Length,
		Width,
		threshold,
		origin,
		spacing,
		orientation2D
	};

	dataType* houghSpacePtr = new dataType[dim2D]{ 0 };
	dataType* foundCirclePtr = new dataType[dim2D]{ 0 };
	size_t count_slice = 0;

	for (k = 0; k < Height; k++) {

		Point2D centroid_seed = get2dImagecentroid(imageData[k], Length, Width, 0.0);
		std::cout << "Found Centroid = (" << centroid_seed.x << "," << centroid_seed.y << ")" << std::endl;
		initialize2dArray(foundCirclePtr, Length, Width);
		initialize2dArray(threshold, Length, Width);
		initialize2dArray(gradientX, Length, Width);
		initialize2dArray(gradientY, Length, Width);

		if (centroid_seed.x != 0 && centroid_seed.y != 0) {
			
			count_slice++;

			//get the liver slice
			copyDataToAnother2dArray(ctContainer->dataPointer[k], imageSlice, Length, Width);
			
			//Rescale to range 0-1
			rescaleNewRange2D(imageSlice, Length, Width, 0.0, 1.0);

			//Slice filtering
			heatImplicit2dScheme(IMAGE, filter_parameters);

			//Compute gradient
			computeImageGradient(IMAGE, gradientX, gradientY);

			//Compute edge detector and threshold
			norm_of_gradient = 0.0;
			edge_value = 0.0;
			for (i = 0; i < dim2D; i++) {
				norm_of_gradient = sqrt(gradientX[i] * gradientX[i] + gradientY[i] * gradientY[i]);
				edge_value = gradientFunction(norm_of_gradient, h_parameters.K);
				if (edge_value <= h_parameters.thres) {
					threshold[i] = 1.0;
				}
			}

			Point2D point_found = localCircleDetection(IMAGE_HOUGH, foundCirclePtr, centroid_seed, h_parameters);

			storing_path = outputPath + "hough image/input/slice_" + to_string(count_slice) + ".raw";
			manageRAWFile2D<dataType>(imageSlice, Length, Width, storing_path.c_str(), STORE_DATA, false);

			storing_path = outputPath + "hough image/threshold/threshold_" + to_string(count_slice) + ".raw";
			manageRAWFile2D<dataType>(threshold, Length, Width, storing_path.c_str(), STORE_DATA, false);

			storing_path = outputPath + "hough image/found_circle/circle_" + to_string(count_slice) + ".raw";
			manageRAWFile2D<dataType>(foundCirclePtr, Length, Width, storing_path.c_str(), STORE_DATA, false);

		}

		std::cout << "Slice " << k << std::endl;
	}

	//copyDataToAnother2dArray(ctContainer->dataPointer[k_slice], imageSlice, Length, Width);
	
	//loading_path = inputPath + "raw/slice lungs centroid/slice_p6.raw";
	//manageRAWFile2D(imageSlice, Length, Width, loading_path.c_str(), LOAD_DATA, false);
	rescaleNewRange2D(imageSlice, Length, Width, 0.0, 1.0);

	storing_path = outputPath + "input.raw";
	manageRAWFile2D<dataType>(imageSlice, Length, Width, storing_path.c_str(), STORE_DATA, false);

	PixelSpacing spacing = { ctSpacing.sx, ctSpacing.sy };
	Point2D origin = { 0.0, 0.0 };
	OrientationMatrix2D orientation2D;
	orientation2D.v1 = { 1.0, 0.0};
	orientation2D.v2 = { 0.0, 1.0};

	Image_Data2D IMAGE {
		Length,
		Width,
		imageSlice,
		origin,
		spacing,
		orientation2D
	};

	//filtering parameters
	Filter_Parameters filter_parameters{
		0.25 * spacing.sx * spacing.sx, //tau
		spacing.sx, //h
		0.0, //sigma
		0.0, //edge detector coef (not used)
		1.5, //omega c
		1e-3, //tolerance
		1e-6, //epsilon2
		1e-6, //coef
		1, //p Neumann boundary condition
		5, //time step
		1000, //max solver iteration
	};

	//Slice filtering
	heatImplicit2dScheme(IMAGE, filter_parameters);

	//storing_path = outputPath + "filtered.raw";
	//manageRAWFile2D<dataType>(imageSlice, Length, Width, storing_path.c_str(), STORE_DATA, false);

	HoughParameters h_parameters =
	{
		5.0,   //minimal radius
		15.0,  //maximal radius
		0.5,   //radius step
		0.5,   //epsilon
		100,    //offset 
		1000,  //K edge detector coefficient
		0.15,  //threshold edge detector
	};

	dataType* gradientX = new dataType[dim2D]{ 0 };
	dataType* gradientY = new dataType[dim2D]{ 0 };
	dataType* threshold = new dataType[dim2D]{ 0 };
	computeImageGradient(IMAGE, gradientX, gradientY);
	dataType norm_of_gradient = 0.0, edge_value;
	for (i = 0; i < dim2D; i++) {
		norm_of_gradient = sqrt(gradientX[i] * gradientX[i] + gradientY[i] * gradientY[i]);
		edge_value = gradientFunction(norm_of_gradient, h_parameters.K);
		if (edge_value <= h_parameters.thres) {
			threshold[i] = 1.0;
		}
	}
	
	storing_path = outputPath + "threshold.raw";
	manageRAWFile2D<dataType>(threshold, Length, Width, storing_path.c_str(), STORE_DATA, false);

	Image_Data2D IMAGE_HOUGH
	{
		Length,
		Width,
		threshold,
		origin,
		spacing,
		orientation2D
	};

	dataType* houghSpacePtr = new dataType[dim2D]{ 0 };
	dataType* foundCirclePtr = new dataType[dim2D]{ 0 };

	//Point2D point_found = localCircleDetection(IMAGE_HOUGH, foundCirclePtr, centroid_seed, h_parameters);

	//Point2D point_found = localHoughTransform(IMAGE_HOUGH, centroid_seed, houghSpacePtr, foundCirclePtr, h_parameters);
	//Point2D point_found = localHoughWithCanny(IMAGE_HOUGH, centroid_seed, houghSpacePtr, foundCirclePtr, h_parameters);

	storing_path = outputPath + "found_circle.raw";
	manageRAWFile2D<dataType>(foundCirclePtr, Length, Width, storing_path.c_str(), STORE_DATA, false);
	
	//localCircleDetection(IMAGE_HOUGH, centroid_seed, h_parameters);
	//circleDetection(IMAGE_HOUGH, h_parameters);
	//storing_path = outputPath + "smoothed.raw";
	//store2dRawData<dataType>(imageSlice, Length, Width, storing_path.c_str());
	
	//Point2D seed = { 248, 298 }; //---> P2 : automatic detection
	//Point2D seed = { 263, 257 }; //---> P2 : abdominal aorta
	//Point2D seed = { 263, 296 }; //---> P2
	//Point2D seed = { 263, 296 }; //---> P2
	//Point2D seed = { 262, 269 }; //---> P4
	//Point2D seed = { 257, 244 }; //---> P2, liver
	//Point2D found_point = { 0.0, 0.0 };
	//Point2D found1 = { 292,314 }, found2 = { 253,255 }; // p2
	//imageSlice[x_new(292, 314, Length)] = 0.0;
	//imageSlice[x_new(253,255, Length)] = 0.0;

	////loading_path = inputPath + "threshold.raw";
	//loading_path = inputPath + "loaded.raw";
	//load2dArrayRAW<dataType>(imageSlice, Length, Width, loading_path.c_str(), false);
	
	//houghSpace[x_new(292, 314, Length)] = 1.0;
	//houghSpace[x_new(253,255, Length)] = 1.0;
	//houghSpace[x_new(264, 263, Length)] = 1.0;//p2,liver
	//imageSlice[x_new(264, 263, Length)] = 0.0;//p2,liver
	//imageSlice[x_new(292, 314, Length)] = 0.0;
	//imageSlice[x_new(253,255, Length)] = 0.0;
	//storing_path = outputPath + "input_plus_centers.raw";
	//store2dRawData<dataType>(imageSlice, Length, Width, storing_path.c_str());

	//string saving_csv = outputPath + "file_.csv";
	//FILE* file;
	//if (fopen_s(&file, saving_csv.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	
	//found_point = localHoughTransform(seed, imageSlice, houghSpace, foundCircle, Length, Width, parameters, storing_path.c_str());
	//cout << "found center : (" << found_point.x << "," << found_point.y << ")" << endl;
	//fclose(file);

	delete[] imageSlice;
	delete[] statusPixel;
	delete[] gradientX;
	delete[] gradientY;
	delete[] threshold;	

	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
	}
	delete[] imageData;
	free(ctContainer);
	*/
	
	//======================== Detect carina from segment =============================
	
	/*
	dataType** imageData = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{ 0 };
	}

	loading_path = inputPath + "trachea.raw";
	load3dArrayRAW<dataType>(imageData, Length, Width, Height, loading_path.c_str(), false);
	
	int* labelArray = new int[dim2D];
	bool* status = new bool[dim2D];
	dataType* mask = new dataType[dim2D];
	if (labelArray == NULL || status == NULL || mask == NULL)
		return false;

	size_t n, count_slice = 0, slice_of_interest = 0;
	bool previous_slice = false;

	storing_path = outputPath + "info_regions.csv";
	FILE* file;
	if (fopen_s(&file, storing_path.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	fprintf(file, "slice number,number of regions\n");

	//size_t k_previous = 0, k_current = 0, k_next = 0;

	for (k = 0; k < Height; k++) {

		//if (k == 0) {
		//	k_previous = 0;
		//	k_current = k_previous;
		//	k_next = 1;
		//}
		//else {
		//	if (k == Height - 1) {
		//		k_previous = Height - 1;
		//		k_current = Height - 1;
		//		k_next = k_current;
		//	}
		//	else {
		//		k_previous = k - 1;
		//		k_current = k;
		//		k_next = k + 1;
		//	}
		//}

		//cout << "Slice " << k << ":" << endl;
		fprintf(file, "%d,", k);
		
		//Initialization
		for (i = 0; i < dim2D; i++) {
			mask[i] = imageData[k][i];
			labelArray[i] = 0;
			status[i] = false;
		}

		labeling(mask, labelArray, status, Length, Width, 1.0);

		//number of region cells
		int numberOfRegionsCells = 0;
		for (i = 0; i < dim2D; i++) {
			if (mask[i] == 1.0) {
				numberOfRegionsCells++;
			}
		}

		if (numberOfRegionsCells != 0) {

			int* countingArray = new int[numberOfRegionsCells] {0};

			//Get regions sizes
			for (i = 0; i < dim2D; i++) {
				if (labelArray[i] > 0) {
					countingArray[labelArray[i]]++;
				}
			}

			//Get number of regions
			int numb_region = 0, label = 1;
			for (i = 0; i < dim2D; i++) {
				if (labelArray[i] == label) {
					numb_region++;
					label++;
				}
			}
			
			////Print number of elements per region
			//for (i = 1; i <= numb_region; i++) {
			//	std::cout << "region " << i << " has " << countingArray[i] << " elements" << std::endl;
			//}
			//std::cout << "############" << std::endl;

			//std::cout << "Slice " << k  << " has " << numb_region << " region(s)" << std::endl;
			fprintf(file, "%d\n", numb_region);

			if (numb_region > 1) {
				count_slice++;
				previous_slice = true;
			}
			else {
				count_slice = 0;
				previous_slice = false;
			}

			if (count_slice == 5 && previous_slice == true) {
				slice_of_interest = k - 5;
				count_slice = 0;
				std::cout << "The slice of interest is : " << slice_of_interest << std::endl;
			}

			delete[] countingArray;

		}
		else {
			fprintf(file, "%d\n", numberOfRegionsCells);
			previous_slice = false;
			count_slice = 0;
		}
		//std::cout << "count slice " << count_slice << std::endl;

	}

	fclose(file);
	delete[] labelArray;
	delete[] status;
	*/

	//======================== Detect carina using path extraction ====================

	/*
	dataType** imageData = new dataType * [Height];
	dataType** maskLungs = new dataType * [Height];
	dataType** ball_centroid = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{ 0 };
		maskLungs[k] = new dataType[dim2D]{ 0 };
		ball_centroid[k] = new dataType[dim2D]{ 0 };
	}

	//loading_path = inputPath + "raw/filtered/filteredGMC_p2.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	loading_path = inputPath + "raw/lungs/-150HU/lungs_p6.raw";
	//loading_path = inputPath + "raw/lungs/-500HU/lungs_p3.raw";
	//loading_path = inputPath + "raw/lungs/iterative/lungs_p3.raw";
	manageRAWFile3D<dataType>(maskLungs, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	for (size_t it = 1; it <= 5; it++) {
		dilatation3D(maskLungs, Length, Width, Height, 1.0, 0.0);
		erosion3D(maskLungs, Length, Width, Height, 1.0, 0.0);
		//storing_path = outputPath + "dilatation_" + to_string(it) + ".raw";
		//manageRAWFile3D<dataType>(maskLungs, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);
	}

	dataType* centroid_lungs = new dataType[3]{0};
	centroidImage(maskLungs, centroid_lungs, Height, Length, Width, 0.0);

	Point3D pCentroid = { centroid_lungs[0], centroid_lungs[1], centroid_lungs[2] };
	Point3D pCentroidReal = getRealCoordFromImageCoord3D(pCentroid, ctOrigin, ctSpacing, orientation);

	string test_meta_data = outputPath + "metaDataTest_p6.txt";
	FILE* test_description;
	if (fopen_s(&test_description, test_meta_data.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	fprintf(test_description, "This test is performed on Monday 10th, February 2025\n");
	fprintf(test_description, "The objective is to find points inside the right and left lungs to detect the carina\n");
	fprintf(test_description, "\nPatient 6: origin = (%f, %f, %f)\n", ctOrigin.x, ctOrigin.y, ctOrigin.z);
	fprintf(test_description, "Voxel sizes = (%f, %f, %f)\n", ctSpacing.sx, ctSpacing.sy, ctSpacing.sz);
	fprintf(test_description, "Dimension : Length = %d , Width = %d and Height = %d\n", Length, Width, Height);
	fprintf(test_description, "Segmented lungs centroid : (%f,%f,%f)\n", pCentroidReal.x, pCentroidReal.y, pCentroidReal.z);
	fclose(test_description);
	
	//Draw ball centroid
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				Point3D current_point = { i, j, k };
				current_point = getRealCoordFromImageCoord3D(current_point, ctOrigin, ctSpacing, orientation);
				double distance_point = getPoint3DDistance(pCentroidReal, current_point);
				if (distance_point <= 10) {
					maskLungs[k][x_new(i, j, Length)] = 0;
					ball_centroid[k][x_new(i, j, Length)] = 1.0;
				}
			}
		}
	}
	storing_path = outputPath + "p6/ball_centroid.raw";
	manageRAWFile3D<dataType>(ball_centroid, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	//================== Approx left and right lung =============
	dataType** leftLung = new dataType * [Height];
	dataType** rightLung = new dataType * [Height];
	dataType** ball_left = new dataType * [Height];
	dataType** ball_right = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		leftLung[k] = new dataType[dim2D]{ 0 };
		ball_left[k] = new dataType[dim2D]{ 0 };
		rightLung[k] = new dataType[dim2D]{ 0 };
		ball_right[k] = new dataType[dim2D]{0};
	}

	//Separate left and right lungs
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				xd = x_new(i, j, Length);
				if (i < pCentroid.x) {
					leftLung[k][xd] = maskLungs[k][xd];
					rightLung[k][xd] = 0.0;
				}
				else {
					leftLung[k][xd] = 0.0;
					rightLung[k][xd] = maskLungs[k][xd];
				}
			}
		}
	}
	storing_path = outputPath + "p6/left.raw";
	manageRAWFile3D<dataType>(leftLung, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);
	storing_path = outputPath + "p6/right.raw";
	manageRAWFile3D<dataType>(rightLung, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	
	dataType** distanceMap = new dataType * [Height];
	dataType** actionMap = new dataType * [Height];
	dataType** potential = new dataType * [Height];
	dataType** ball_maximum = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		distanceMap[k] = new dataType[dim2D]{ 0 };
		actionMap[k] = new dataType[dim2D]{ 0 };
		potential[k] = new dataType[dim2D]{ 0 };
		ball_maximum[k] = new dataType[dim2D]{ 0 };
	}

	////========= Interpolation ====
	//Image_Data inputImageStr { Height, Length, Width, imageData, ctOrigin, ctSpacing, orientation};
	//Image_Data intImageStr{ height, Length, Width, intImageData, ctOrigin, intSpacing, orientation };
	//imageInterpolation3D(inputImageStr, intImageStr, NEAREST_NEIGHBOR);
	//Image_Data inputLungsStr { Height, Length, Width, maskLungs, ctOrigin, ctSpacing, orientation};
	//Image_Data intLungsStr{ height, Length, Width, intMaskLungs, ctOrigin, intSpacing, orientation };
	//imageInterpolation3D(inputLungsStr, intLungsStr, NEAREST_NEIGHBOR);
	//============================

	Image_Data distanceMapStr
	{
		Height,
		Length,
		Width,
		maskLungs,
		ctOrigin,
		ctSpacing,
		orientation
	};

	dataType tau = 0.4;
	dataType tolerance = 0.5;
	rouyTourinDistanceMap(distanceMapStr, distanceMap, 1.0, tolerance, tau);
	//rouyTourinFunction_3D(distanceMap, intMaskLungs, 0.5, Length, Width, height, 0.4, intSpacing.sx);
	//fastSweepingFunction_3D(distanceMap, intMaskLungs, Length, Width, height, intSpacing.sx, 10000000.0, 0.0);
	storing_path = outputPath + "distance_map_RT.raw";
	manageRAWFile3D<dataType>(distanceMap, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	Point3D pMaxDist = { 0.0, 0.0, 0.0 };
	
	dataType max_distance = 0.0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				xd = x_new(i, j, Length);
				if (distanceMap[k][xd] > max_distance) {
					max_distance = distanceMap[k][xd];
					pMaxDist.x = i;
					pMaxDist.y = j;
					pMaxDist.z = k;
				}
			}
		}
	}

	//Point3D first_point = { 267, 243, 315 };//p1
	//Point3D first_point = { 251, 243, 434 };//p2
	Point3D first_point = { 271, 231, 274 };//p3
	//Point3D first_point = { 265, 231, 270 };//p4
	//Point3D first_point = { 241, 205, 704 };//p5
	//Point3D first_point = { 248, 299, 520 };//p6
	//first_point = getRealCoordFromImageCoord3D(first_point, ctOrigin, ctSpacing, orientation);
	//first_point = getImageCoordFromRealCoord3D(first_point, ctOrigin, intSpacing, orientation);
	Point3D* seedPoints = new Point3D[2];
	seedPoints[0] = first_point;
	seedPoints[1] = pMaxDist;

	Image_Data inputImageData
	{
		Height,
		Length,
		Width,
		imageData,
		ctOrigin,
		ctSpacing,
		orientation
	};
	
	double radius = 3.0;
	Potential_Parameters parameters{
		1000, //edge detector coefficient
		0.15, //threshold
		0.005,//epsilon
		radius
	};
	compute3DPotential(inputImageData, potential, first_point, parameters);
	//storing_path = outputPath + "potential.raw";
	//manageRAWFile3D<dataType>(potential, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);
	fastMarching3dWithSpacing(inputImageData, actionMap, potential, first_point, ctSpacing);
	//fastMarching3D_N(actionMap, potential, Length, Width, height, first_point);
	storing_path = outputPath + "action.raw";
	manageRAWFile3D<dataType>(actionMap, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	vector<Point3D> path_points;
	//shortestPath3D(actionMap, Length, Width, Height, ctSpacing, seedPoints, path_points);

	//copy points to file
	string saving_csv = outputPath + "path_points_L.csv";
	FILE* file_left;
	if (fopen_s(&file_left, saving_csv.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	int n = 0;
	for (n = path_points.size() - 1; n > -1; n--) {
		path_points[n] = getRealCoordFromImageCoord3D(path_points[n], ctOrigin, ctSpacing, orientation);
		fprintf(file_left, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	}
	fclose(file_left);
	while (path_points.size() > 0) {
		path_points.pop_back();
	}

	//Remove maximal region
	Point3D pMaxDistReal = getRealCoordFromImageCoord3D(pMaxDist, ctOrigin, ctSpacing, orientation);
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				distanceMap[k][x_new(i, j, Length)] = 0;
				Point3D current_point = { i, j, k };
				current_point = getRealCoordFromImageCoord3D(current_point, ctOrigin, ctSpacing, orientation);
				double distance_point = getPoint3DDistance(pMaxDistReal, current_point);
				if (distance_point <= max_distance) {
					maskLungs[k][x_new(i, j, Length)] = 0;
					ball_maximum[k][x_new(i, j, Length)] = 1.0;
				}
			}
		}
	}
	storing_path = outputPath + "ball_max.raw";
	manageRAWFile3D<dataType>(ball_maximum, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);
	
	rouyTourinDistanceMap(distanceMapStr, distanceMap, 1.0, tolerance, tau);
	//rouyTourinFunction_3D(distanceMap, intMaskLungs, 0.5, Length, Width, height, 0.4, intSpacing.sx);
	//fastSweepingFunction_3D(distanceMap, intMaskLungs, Length, Width, height, intSpacing.sx, 10000000.0, 0.0);
	pMaxDist = getPointWithTheHighestValue(distanceMap, Length, Width, Height);
	seedPoints[1] = pMaxDist;
	//shortestPath3D(actionMap, Length, Width, Height, ctSpacing, seedPoints, path_points);

	saving_csv = outputPath + "path_points_R.csv";
	FILE* file_right;
	if (fopen_s(&file_right, saving_csv.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	n = 0;
	for (n = path_points.size() - 1; n > -1; n--) {
		path_points[n] = getRealCoordFromImageCoord3D(path_points[n], ctOrigin, ctSpacing, orientation);
		fprintf(file_right, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	}
	fclose(file_right);
	while (path_points.size() > 0) {
		path_points.pop_back();
	}

	delete[] seedPoints;

	for (k = 0; k < Height; k++) {
		delete[] distanceMap[k];
		delete[] actionMap[k];
		delete[] potential[k];
		delete[] ball_maximum[k];
	}
	delete[] distanceMap;
	delete[] actionMap;
	delete[] potential;
	delete[] ball_maximum;

	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
		delete[] maskLungs[k];
		delete[] ball_centroid[k];
		//delete[] leftLung[k];
		//delete[] rightLung[k];
		//delete[] ball_left[k];
		//delete[] ball_right[k];
	}
	delete[] imageData;
	delete[] maskLungs;
	delete[] ball_centroid;
	//delete[] leftLung;
	//delete[] rightLung;
	//delete[] ball_left;
	//delete[] ball_right;

	free(ctContainer);
	*/
	
	//===================== double distance map lungs =====================
	
	/*
	dataType** imageData = new dataType * [Height];
	dataType** maskLungs = new dataType * [Height];
	dataType** distanceMap = new dataType * [Height];
	dataType** action = new dataType * [Height];
	dataType** potential = new dataType * [Height];
	dataType** ball_max = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{ 0 };
		maskLungs[k] = new dataType[dim2D]{ 0 };
		distanceMap[k] = new dataType[dim2D]{ 0 };
		action[k] = new dataType[dim2D]{ 0 };
		potential[k] = new dataType[dim2D]{ 0 };
		ball_max[k] = new dataType[dim2D]{ 0 };
	}

	//loading_path = inputPath + "Left Right Lungs/p3/left.raw";
	//loading_path = inputPath + "raw/lungs/-150HU/lungs_p3.raw";
	//loading_path = inputPath + "raw/lungs/-500HU/lungs_p3.raw";
	loading_path = inputPath + "raw/lungs/iterative/lungs_p3.raw";
	manageRAWFile3D<dataType>(maskLungs, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);
	
	dilatation3D(maskLungs, Length, Width, Height, 1.0, 0.0);
	dilatation3D(maskLungs, Length, Width, Height, 1.0, 0.0);
	dilatation3D(maskLungs, Length, Width, Height, 1.0, 0.0);
	erosion3D(maskLungs, Length, Width, Height, 1.0, 0.0);
	erosion3D(maskLungs, Length, Width, Height, 1.0, 0.0);
	erosion3D(maskLungs, Length, Width, Height, 1.0, 0.0);
	dilatation3D(maskLungs, Length, Width, Height, 1.0, 0.0);
	dilatation3D(maskLungs, Length, Width, Height, 1.0, 0.0);
	dilatation3D(maskLungs, Length, Width, Height, 1.0, 0.0);
	erosion3D(maskLungs, Length, Width, Height, 1.0, 0.0);
	erosion3D(maskLungs, Length, Width, Height, 1.0, 0.0);
	erosion3D(maskLungs, Length, Width, Height, 1.0, 0.0);
	
	//storing_path = outputPath + "update_input.raw";
	//manageRAWFile3D<dataType>(maskLungs, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	dataType tau = 0.4;
	dataType tolerance = 0.5;
	Image_Data distanceStr = { Height, Length, Width, maskLungs, ctOrigin, ctSpacing, orientation };
	rouyTourinDistanceMap(distanceStr, distanceMap, 1.0, tolerance, tau);

	//storing_path = outputPath + "distance.raw";
	//manageRAWFile3D<dataType>(distanceMap, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	Point3D pMaxDistLeft = { 0.0, 0.0, 0.0 };
	
	//Find maximum distance
	dataType max_distance = 0.0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				xd = x_new(i, j, Length);
				if (distanceMap[k][xd] > max_distance) {
					max_distance = distanceMap[k][xd];
					pMaxDistLeft.x = i;
					pMaxDistLeft.y = j;
					pMaxDistLeft.z = k;
				}
			}
		}
	}

	Point3D prMaxDist = getRealCoordFromImageCoord3D(pMaxDistLeft, ctOrigin, ctSpacing, orientation);

	//remove distance point
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				xd = x_new(i, j, Length);
				Point3D current_point = { i, j, k };
				//distanceMap[k][xd] = 0.0;
				current_point = getRealCoordFromImageCoord3D(current_point, ctOrigin, ctSpacing, orientation);
				double pDistance = getPoint3DDistance(current_point, prMaxDist);
				if (pDistance <= max_distance) {
					distanceMap[k][xd] = 0.0;
					ball_max[k][xd] = 1.0;
					//maskLungs[k][i] = 0.0;
				}
			}
		}
	}
	//storing_path = outputPath + "visualize_max.raw";
	//manageRAWFile3D<dataType>(ball_max, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	//storing_path = outputPath + "updated_lungs.raw";
	//manageRAWFile3D<dataType>(maskLungs, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	//rouyTourinDistanceMap(distanceStr, distanceMap, 1.0, tolerance, tau);
	Point3D pMaxDistRight = getPointWithTheHighestValue(distanceMap, Length, Width, Height);

	//storing_path = outputPath + "updated_distance.raw";
	//manageRAWFile3D<dataType>(distanceMap, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	//Point3D pMaxDistRight = getPointWithTheHighestValue(distanceMap, Length, Width, Height);
	//Point3D first_point = { 263, 260, 294 };//trachea, p1
	//Point3D first_point = { 251, 243, 434 };//trachea, p2
	Point3D first_point = { 271, 231, 274 };//trachea, p3
	//Point3D first_point = { 265, 231, 270 };//Trachea, p4
	//Point3D first_point = { 241, 205, 704 };//Trachea, p5
	//Point3D first_point = { 248, 298, 520 };//Trachea, p6
	double radius = 3.0;
	Potential_Parameters parameters{
		1000, //edge detector coefficient
		0.15, //threshold
		0.005,//epsilon
		radius
	};
	Image_Data inputImage = { Height, Length, Width, imageData, ctOrigin, ctSpacing, orientation };

	//loading_path = inputPath + "raw/filtered/filteredGMC_p3.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	////Copy
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		imageData[k][i] = ctContainer->dataPointer[k][i];
	//	}
	//}
	//rescaleNewRange(imageData, Length, Width, Height, 0.0, 1.0, maxData, minData);

	//compute3DPotential(inputImage, potential, first_point, parameters);
	//storing_path = outputPath + "potential.raw";
	//manageRAWFile3D<dataType>(potential, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	//fastMarching3dWithSpacing(inputImage, action, potential, first_point, ctSpacing);
	storing_path = outputPath + "action.raw";
	manageRAWFile3D<dataType>(action, Length, Width, Height, storing_path.c_str(), LOAD_DATA, false);

	vector<Point3D> path_points;

	string saving_csv = outputPath + "path_left.csv";
	FILE* file_segment_left;
	if (fopen_s(&file_segment_left, saving_csv.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}

	Point3D* seedPoints = new Point3D[2];
	seedPoints[0] = first_point;
	seedPoints[1] = pMaxDistLeft;
	std::cout << "Max left : (" << pMaxDistLeft.x << "," << pMaxDistLeft.y << "," << pMaxDistLeft.z << ")" << std::endl;
	Image_Data actionMapStr = { Height, Length, Width, action, ctOrigin, ctSpacing, orientation };
	shortestPath3D(actionMapStr, seedPoints, 0.8, 1.0, path_points);
	for (int n = path_points.size() - 1; n > -1; n--) {
		path_points[n] = getRealCoordFromImageCoord3D(path_points[n], ctOrigin, ctSpacing, orientation);
		fprintf(file_segment_left, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	}
	while (path_points.size() != 0) {
		path_points.pop_back();
	}
	fclose(file_segment_left);

	saving_csv = outputPath + "path_right.csv";
	FILE* file_segment_right;
	if (fopen_s(&file_segment_right, saving_csv.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	seedPoints[1] = pMaxDistRight;
	std::cout << "Max Right : (" << pMaxDistRight.x << "," << pMaxDistRight.y << "," << pMaxDistRight.z << ")" << std::endl;
	shortestPath3D(actionMapStr, seedPoints, 0.8, 1.0, path_points);
	for (int n = path_points.size() - 1; n > -1; n--) {
		path_points[n] = getRealCoordFromImageCoord3D(path_points[n], ctOrigin, ctSpacing, orientation);
		fprintf(file_segment_right, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	}
	while (path_points.size() != 0) {
		path_points.pop_back();
	}
	fclose(file_segment_right);

	delete[] seedPoints;

	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
		delete[] maskLungs[k];
		delete[] distanceMap[k];
		delete[] action[k];
		delete[] potential[k];
		delete[] ball_max[k];
	}
	delete[] imageData;
	delete[] maskLungs;
	delete[] distanceMap;
	delete[] action;
	delete[] potential;
	delete[] ball_max;
	*/

	//==================== Analyze the lung connected to the trachea slice by slice =============
	
	/*
	dataType** imageData = new dataType * [Height];
	dataType** maskLungs = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{ 0 };
		maskLungs[k] = new dataType[dim2D]{ 0 };
	}

	loading_path = inputPath + "raw/lungs/-150HU/lungs_p4.raw";
	//loading_path = inputPath + "raw/lungs/-500HU/lungs_p4.raw";
	//loading_path = inputPath + "raw/lungs/iterative/lungs_p4.raw";
	manageRAWFile3D<dataType>(maskLungs, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	//storing_path = outputPath + "loaded_P4.raw";
	//manageRAWFile3D<dataType>(maskLungs, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);
	for (size_t iter = 0; iter < 1; iter++) {
		dilatation3D(maskLungs, Length, Width, Height, 1.0, 0.0);
		erosion3D(maskLungs, Length, Width, Height, 1.0, 0.0);
	}
	//dilatation3D(maskLungs, Length, Width, Height, 1.0, 0.0);
	//dilatation3D(maskLungs, Length, Width, Height, 1.0, 0.0);
	//erosion3D(maskLungs, Length, Width, Height, 1.0, 0.0);
	//erosion3D(maskLungs, Length, Width, Height, 1.0, 0.0);
	//storing_path = outputPath + "updated_P6.raw";
	//manageRAWFile3D<dataType>(maskLungs, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	dataType* lungSlice = new dataType[dim2D]{ 0 };
	int* segment = new int[dim2D]{ 0 };
	bool* status = new bool[dim2D] { false };

	//string saving_info_regions = outputPath + "Lungs Trachea Labeling/number_of_regions_p6.csv";
	//FILE* file_regions;
	//if (fopen_s(&file_regions, saving_info_regions.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(file_regions, "slice,nbr\n");

	size_t consecutive_grp1 = 0, count_grp = 0;
	size_t consecutive_grp2 = 0;
	size_t consecutive_grp3 = 0;
	size_t previous_slice = 0, current_slice = 0, next_slice = 0;
	size_t total_slice_1_elt = 0;

	vector<size_t> slice_id;

	for (k = 0; k < Height; k++) {

		//copy
		copyDataToAnother2dArray(maskLungs[k], lungSlice, Length, Width);

		//get number of forepixel
		int numberOfForegroundPixels = 0;
		for (i = 0; i < dim2D; i++) {
			if (lungSlice[i] == 1.0) {
				numberOfForegroundPixels++;
			}
		}

		//Initialize array
		for (i = 0; i < dim2D; i++) {
			segment[i] = 0;
			status[i] = false;
		}

		if (numberOfForegroundPixels > 0 && numberOfForegroundPixels < dim2D) {
			
			int* countingArray = new int[numberOfForegroundPixels] {0};

			labeling(lungSlice, segment, status, Length, Width, 1.0);

			//Get regions size
			for (i = 0; i < dim2D; i++) {
				if (segment[i] > 0) {
					countingArray[segment[i]]++;
				}
			}
			
			////Remove too small regions
			//size_t regionSize = 15;
			//for (i = 0; i < dim2D; i++) {
			//	if (countingArray[segment[i]] < regionSize) {
			//		countingArray[segment[i]] = 0;
			//	}
			//}

			//get number of regions
			size_t numberOfRegions = 0;
			for (i = 0; i < numberOfForegroundPixels; i++) {
				if (countingArray[i] > 0) {
					numberOfRegions++;
				}
			}
			////std::cout << numberOfRegions << " regions for slice " << k << std::endl;
			//fprintf(file_regions, "%d,%d\n", k, numberOfRegions);

			//Deal with consecutive slices
			if (numberOfRegions == 1) {
				
				slice_id.push_back(k);

				//total_slice_1_elt++;
				//consecutive_grp1++;
				//current_slice = k;
				//if (consecutive_grp1 == 1) {
				//	count_grp++;
				//	next_slice = k + 1;
				//	slice_id.push_back(k);
				//}
				//else {
				//	if (next_slice == current_slice) {
				//		next_slice = current_slice + 1;
				//		slice_id.push_back(k);
				//	}
				//	else {
				//		consecutive_grp1 = 0;
				//		std::cout << "Slices group " << count_grp << " : " << std::endl;
				//		for (size_t itn = 0; itn < slice_id.size(); itn++) {
				//			std::cout << slice_id[itn] << ", ";
				//		}
				//		std::cout << "" << std::endl;
				//		while (slice_id.size() > 0) {
				//			slice_id.pop_back();
				//		}
				//		//slice_id.push_back(k);
				//	}
				//}
				
				//storing_path = outputPath + "one region/slice_" + to_string(k) + ".raw";
				//manageRAWFile2D<dataType>(lungSlice, Length, Width, storing_path.c_str(), STORE_DATA, false);

			}
			
			//if (numberOfRegions == 2) {
			//	storing_path = outputPath + "two regions/slice_" + to_string(k) + ".raw";
			//	manageRAWFile2D<dataType>(lungSlice, Length, Width, storing_path.c_str(), STORE_DATA, false);
			//}
			//if (numberOfRegions == 3) {
			//	storing_path = outputPath + "three regions/slice_" + to_string(k) + ".raw";
			//	manageRAWFile2D<dataType>(lungSlice, Length, Width, storing_path.c_str(), STORE_DATA, false);
			//}

			delete[] countingArray;
		}
		//else {
		//	//std::cout << "0 region in slice " << k << std::endl;
		//	size_t numberOfRegions = 0;
		//	fprintf(file_regions, "%d,%d\n", k, numberOfRegions);
		//}
	}
	//fclose(file_regions);

	size_t count_grp_elt = 0;
	count_grp = 1;
	std::cout << "" << std::endl;
	for (size_t itn = 0; itn < slice_id.size() ; itn++) {

		if (itn == 0) {
			previous_slice = slice_id[itn];
			current_slice = slice_id[itn];
			next_slice = slice_id[itn + 1];
			count_grp_elt++;
		}else{
			if (itn == slice_id.size() - 1) {
				previous_slice = slice_id[itn - 1];
				current_slice = slice_id[itn];
				next_slice = slice_id[itn];
				count_grp_elt++;
			}
			else {
				previous_slice = slice_id[itn - 1];
				current_slice = slice_id[itn];
				next_slice = slice_id[itn + 1];
				count_grp_elt++;
			}
		}

		if ((current_slice - previous_slice) > 1) { 
			std::cout << "Group " << count_grp << " has " << count_grp_elt << " elements" << std::endl;
			count_grp++;
			count_grp_elt = 0;
		}
	}
	std::cout << "Group " << count_grp << " has " << count_grp_elt << " elements" << std::endl;
	std::cout << "\nWe found " << count_grp << " consecutive slices" << std::endl;
	std::cout << "\We found " << slice_id.size() << " slices with 1 segment" << std::endl;

	delete[] lungSlice;
	delete[] segment;
	delete[] status;

	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
		delete[] maskLungs[k];
	}
	delete[] imageData;
	delete[] maskLungs;
	free(ctContainer);
	*/

	//==================== Path extraction in Lungs ==================

	/*
	dataType** imageData = new dataType * [Height];
	dataType** potential = new dataType * [Height];
	dataType** action = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{ 0 };
		potential[k] = new dataType[dim2D]{ 0 };
		action[k] = new dataType[dim2D]{ 0 };
	}

	//Copy
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			imageData[k][i] = ctContainer->dataPointer[k][i];
		}
	}
	rescaleNewRange(imageData, Length, Width, Height, 0.0, 1.0, minData, maxData);

	Image_Data inputImageData = { Height, Length, Width, imageData, ctOrigin, ctSpacing, orientation };

	//Patient 1 : drawback in aorta
	Point3D first_point = { 261, 258, 147 };
	Point3D second_point = { 249, 251, 259 };

	////Patient 1
	//Point3D first_point = { 262, 254, 297 };
	//Point3D left = { 200, 294, 231 };
	//Point3D right = { 345, 304, 231 };

	////Patient 2
	//Point3D first_point = { 257, 251, 397 };
	//Point3D left = { 192, 281, 331 };
	//Point3D right = { 338, 288, 331 };

	////Patient 3
	//Point3D first_point = { 272, 231, 273 };
	//Point3D left = { 174, 278, 205 };
	////Point3D right = { 361, 288, 205 };
	//Point3D right = { 313, 286, 251 };//drawback

	////Patient 4
	//Point3D first_point = { 265, 231, 271 };
	//Point3D left = { 177, 295, 211 };
	//Point3D right = { 299, 336, 211 };

	////Patient 5
	//Point3D first_point = { 241, 200, 714 };
	//Point3D left = { 179, 278, 610 };
	//Point3D right = { 317, 267, 610 };

	////Patient 6
	//Point3D first_point = { 249, 291, 531 };
	//Point3D left = { 183, 318, 424 };
	//Point3D right = { 350, 341, 424 };
	
	Point3D* seedPoints = new Point3D[2];
	seedPoints[0] = first_point;

	double radius = 3.0;
	Potential_Parameters parameters{
		1000, //edge detector coefficient
		0.15, //threshold
		0.005,//epsilon
		radius
	};
	compute3DPotential(inputImageData, potential, first_point, parameters);

	fastMarching3dWithSpacing(inputImageData, action, potential, first_point, ctSpacing);

	Image_Data actionStr = { Height, Length, Width, action, ctOrigin, ctSpacing, orientation };
	vector<Point3D> path_points;

	//seedPoints[1] = left;
	seedPoints[1] = second_point;
	shortestPath3D(actionStr, seedPoints, 0.8, 1.0, path_points);

	//copy points to file
	string saving_csv = outputPath + "drawback_aorta.csv";
	FILE* file_left;
	if (fopen_s(&file_left, saving_csv.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	for (int n = path_points.size() - 1; n > -1; n--) {
		path_points[n] = getRealCoordFromImageCoord3D(path_points[n], ctOrigin, ctSpacing, orientation);
		fprintf(file_left, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	}
	fclose(file_left);
	while (path_points.size() > 0) {
		path_points.pop_back();
	}

	//seedPoints[1] = right;
	//shortestPath3D(actionStr, seedPoints, 0.8, 1.0, path_points);

	////copy points to file
	//saving_csv = outputPath + "path_points_right_drawback.csv";
	//FILE* file_right;
	//if (fopen_s(&file_right, saving_csv.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//for (int n = path_points.size() - 1; n > -1; n--) {
	//	path_points[n] = getRealCoordFromImageCoord3D(path_points[n], ctOrigin, ctSpacing, orientation);
	//	fprintf(file_right, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	//}
	//fclose(file_right);
	//while (path_points.size() > 0) {
	//	path_points.pop_back();
	//}

	////Save inputs test
	//string test_meta_data = outputPath + "info_test_carina_p3_drawback.txt";
	//FILE* test_description;
	//if (fopen_s(&test_description, test_meta_data.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(test_description, "The test is performed on Thursday 13th February\n");
	//fprintf(test_description, "The objective is to show that when points are chosen properly in the trachea, the left and rigth lungs, we can detect the carina in robust way using path extraction\n");
	//fprintf(test_description, "For this prototyping the points are set manually. One in the trachea, one in the left lung and one in the right lung\n");
	//fprintf(test_description, "The points should be carefully chosen under the carina for optimal results\n");
	//fprintf(test_description, "The original image is used whithout filtering and interpolation\n");
	//fprintf(test_description, "\Dimension : (Length = %d, Width = %d, Height = %d)\n", Length, Width, Height);
	//fprintf(test_description, "Origin : (%.3f, %.3f, %.3f)\n", ctOrigin.x, ctOrigin.y, ctOrigin.z);
	//fprintf(test_description, "Spacing : (%.3f, %.3f, %.3f)\n", ctSpacing.sx, ctSpacing.sy, ctSpacing.sz);
	//fprintf(test_description, "\nPoint in trachea : (%.1f, %.1f, %.1f)\n", first_point.x, first_point.y, first_point.z);
	//fprintf(test_description, "Point in left lung : (%.1f, %.1f, %.1f)\n", left.x, left.y, left.z);
	//fprintf(test_description, "Point in right lung : (%.1f, %.1f, %.1f)\n", right.x, right.y, right.z);
	//fclose(test_description);

	delete[] seedPoints;

	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
		delete[] potential[k];
		delete[] action[k];
	}
	delete[] imageData;
	delete[] potential;
	delete[] action;
	free(ctContainer);
	*/

	//==================== Segment Trachea ===========================

	/*
	dataType** maskLungs = new dataType * [Height];
	dataType** maskTrachea = new dataType * [Height];
	dataType** difference = new dataType * [Height];
	int** segment = new int* [Height];
	bool** status = new bool* [Height];
	for (k = 0; k < Height; k++) {
		maskLungs[k] = new dataType[dim2D]{ 0 };
		maskTrachea[k] = new dataType[dim2D]{ 0 };
		difference[k] = new dataType[dim2D]{ 0 };
		segment[k] = new int[dim2D]{ 0 };
		status[k] = new bool[dim2D]{ 0 };
	}

	//loading_path = inputPath + "raw/lungs/iterative/lungs_p4.raw";
	loading_path = inputPath + "raw/lungs/-150HU/lungs_p5.raw";
	//loading_path = inputPath + "raw/lungs/-500HU/lungs_p5.raw";
	manageRAWFile3D<dataType>(maskLungs, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	copyDataToAnotherArray(maskLungs, maskTrachea, Height, Length, Width);

	erosion3D(maskTrachea, Length, Width, Height, 1.0, 0.0);
	dilatation3D(maskTrachea, Length, Width, Height, 1.0, 0.0);

	//Difference
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			difference[k][i] = maskLungs[k][i] - maskTrachea[k][i];
		}
	}

	//dilatation3D(maskTrachea, Length, Width, Height, 1.0, 0.0);

	storing_path = outputPath + "difference.raw";
	manageRAWFile3D<dataType>(difference, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	storing_path = outputPath + "remain.raw";
	manageRAWFile3D<dataType>(maskTrachea, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	for (k = 0; k < Height; k++) {
		delete[] maskLungs[k];
		delete[] maskTrachea[k];
		delete[] difference[k];
		delete[] segment[k];
		delete[] status[k];
	}
	delete[] maskLungs;
	delete[] maskTrachea;
	delete[] difference;
	delete[] segment;
	delete[] status;
	free(ctContainer);
	*/

	//==================== GSUBSURF with Peak function ==========

	/*
	//=================== 2D image ==============================
	
	
	//const size_t Length = 512, Width = 512;
	//size_t dim2D = Length * Width;
	
	dataType* imageData = new dataType[dim2D]{ 0 };
	dataType* initialSegment = new dataType[dim2D]{ 0 };

	//loading_path = inputPath + "raw/slice/image2.raw";
	//manageRAWFile2D<dataType>(imageData, Length, Width, loading_path.c_str(), LOAD_DATA, false);
	//rescaleNewRange2D(imageData, Length, Width, 0.0, 1.0);

	copyDataToAnother2dArray(ctContainer->dataPointer[263], imageData, Length, Width);
	rescaleNewRange2D(imageData, Length, Width, 0.0, 1.0);

	dataType v = 1.0;
	dataType radius = 20;
	Point2D sSeed = { 294, 317 };//new
	double pDistance = 0.0;
	dataType value = 0.0;

	//copy points to file
	string saving_csv = outputPath + "peak_function_v2.csv";
	FILE* file_peak;
	if (fopen_s(&file_peak, saving_csv.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	fprintf(file_peak, "x,y\n");
	
	for (i = 0; i < Length; i++) {
		for (j = 0; j < Width; j++) {
			xd = x_new(i, j, Length);
			Point2D current_point = { i, j };
			pDistance = getPoint2DDistance(sSeed, current_point);
			if (pDistance <= radius) {
				value = 1.0 / (pDistance + v);
			}
			else {
				value = 0;//1.0 / (radius + v);
			}
			fprintf(file_peak, "%d,%f\n", xd, value);
			initialSegment[xd] = value;
		}
	}
	fclose(file_peak);

	storing_path = outputPath + "initialSegment_v2.raw";
	manageRAWFile2D<dataType>(initialSegment, Length, Width, storing_path.c_str(), STORE_DATA, false);

	storing_path = outputPath + "loaded_new.raw";
	manageRAWFile2D<dataType>(imageData, Length, Width, storing_path.c_str(), STORE_DATA, false);

	const Filter_Parameters smoothing_parameters = {
		0.25,//tau
		1.171875,//h
		0.0,//sigma
		0,//K--> we use heat implicit
		1.4,//omega
		0.001,//tolerance
		0.001,//eps 2
		0.001,//coef
		1,//p
		5,//time step number
		100,//max solver iteration
	};

	Segmentation_Parameters segmentation_parameters = {
		200,//max iteration number
		10000,//edge detector coef
		0.0001,//eps 2
		100,//number of current time step
		1000,//number of time step
		10,//saving frequency
		0.000001,//segmentation tolerance
		0.1,//tau
		1.171875,//h
		1.4,//omega_c
		0.01,//tolerance
		1.0,//convection coef
		0.1//diffusion coef
	};

	Image_Data2D inputImageData = {Length, Width, imageData};
	
	//storing_path = outputPath + "segmentation/2d slice/";
	//subsurf(inputImageData, initialSegment, storing_path.c_str(), smoothing_parameters, segmentation_parameters);
	storing_path = outputPath + "segmentation/2d slice/gsubsurf/";
	gsubsurf(inputImageData, initialSegment, storing_path.c_str(), smoothing_parameters, segmentation_parameters);

	delete[] imageData;
	delete[] initialSegment;
	*/
	
	//==================== 3D image ====================
	
	/*
	dataType** imageData = new dataType * [Height];
	dataType** maskAorta = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{ 0 };
		maskAorta[k] = new dataType[dim2D]{ 0 };
	}
	 
	loading_path = inputPath + "raw/filtered/New/filtered_p6.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);
	
	Image_Data inputImageData = { Height, Length, Width, imageData, ctOrigin, ctSpacing, orientation };
	
	std::cout << "============ Interpolated ================ " << std::endl;
	size_t hauteur = (size_t)((ctSpacing.sz / ctSpacing.sx) * Height);
	std::cout << "Interpolated Height = " << hauteur << std::endl;
	VoxelSpacing intSpacing = { ctSpacing.sx, ctSpacing.sy, ctSpacing.sx };
	std::cout << "Interpolated Spacing : (" << intSpacing.sx << ", " << intSpacing.sy << ", " << intSpacing.sz << ")" << std::endl;
	
	
	dataType** interpolatedImage = new dataType * [hauteur];
	dataType** maskInterpolation = new dataType * [hauteur];
	for (k = 0; k < hauteur; k++) {
		interpolatedImage[k] = new dataType[dim2D]{ 0 };
		maskInterpolation[k] = new dataType[dim2D]{ 0 };
	}
	Image_Data interpolatedImageData = { hauteur, Length, Width, interpolatedImage, ctOrigin, intSpacing, orientation };
	//imageInterpolation3D(inputImageData, interpolatedImageData, NEAREST_NEIGHBOR);
	//==========================================================
	*/
	
	////p1
	////Croping
	//size_t i_min = 200;
	//size_t j_min = 200;
	//size_t k_min = 275; // 240 for neww 
	//size_t length = 150;
	//size_t width = 150;
	//size_t height = 350; // 380 for new
	
	////p2
	////Croping
	//size_t i_min = 190;
	//size_t j_min = 210;
	//size_t k_min = 477;
	//size_t length = 150;
	//size_t width = 150;
	//size_t height = 350;

	////p3
	////Croping
	//size_t i_min = 190;
	//size_t j_min = 180;
	//size_t k_min = 260;
	//size_t length = 150;
	//size_t width = 150;
	//size_t height = 380;

	////p4
	////Croping
	//size_t i_min = 190;
	//size_t j_min = 190;
	//size_t k_min = 310;
	//size_t length = 150;
	//size_t width = 150;
	//size_t height = 350;

	////p5
	//size_t i_min = 190;
	//size_t j_min = 170;
	//size_t k_min = 695;
	//size_t length = 160;
	//size_t width = 160;
	//size_t height = 360;
	
	////p6
	//size_t i_min = 180;
	//size_t j_min = 245;
	//size_t k_min = 365;
	//size_t length = 150;
	//size_t width = 150;
	//size_t height = 400;
	
	/*
	std::cout << "============ Cropped ================ " << std::endl;
	
	std::cout << "New dimensions : Length = " << length << ", Width = " << width << ", Height = " << height << std::endl;
	std::cout << "Cropped Spacing : (" << intSpacing.sx << ", " << intSpacing.sy << ", " << intSpacing.sz << ")" << std::endl;

	Point3D newOrigin = { i_min, j_min, k_min };
	newOrigin = getRealCoordFromImageCoord3D(newOrigin, ctOrigin, intSpacing, orientation);
	std::cout << "Cropped origin : (" << newOrigin.x << ", " << newOrigin.y << ", " << newOrigin.z << ")" << std::endl;
	*/
	
	/*
	dataType** imageCrop = new dataType * [height];
	dataType** initialSegment = new dataType * [height];
	dataType** maskSegment = new dataType * [height];
	for (k = 0; k < height; k++) {
		imageCrop[k] = new dataType[length * width]{ 0 };
		initialSegment[k] = new dataType[length * width]{ 0 };
		maskSegment[k] = new dataType[length * width]{ 0 };
	}
	Image_Data croppedImageData = { height, length, width, imageCrop, newOrigin, intSpacing, orientation };
	
	size_t k_n = 0, i_n = 0, j_n = 0;
	for (k = 0, k_n = k_min; k < height; k++, k_n++) {
		for (i = 0, i_n = i_min; i < length; i++, i_n++) {
			for (j = 0, j_n = j_min; j < width; j++, j_n++) {
				imageCrop[k][x_new(i, j, length)] = interpolatedImage[k_n][x_new(i_n, j_n, Length)];
			}
		}
	}
	*/

	////========================================================

	/*
	//========== Downsampling =============================
	std::cout << "============ Downsampling ================ " << std::endl;
	size_t h = (size_t)(height / 2);
	size_t w = (size_t)(width / 2);
	size_t l = (size_t)(length / 2);
	std::cout << "New dimensions : Lenght = " << l << ", Width = " << w << ", Height = " << h << std::endl;
	dataType** imageDown = new dataType * [h];
	dataType** initial_segment = new dataType * [h];
	dataType** maskPath = new dataType * [h];
	for (k = 0; k < h; k++) {
		imageDown[k] = new dataType[l * w]{ 0 };
		initial_segment[k] = new dataType[l * w]{ 0 };
		maskPath[k] = new dataType[l * w]{ 0 };
	}
	VoxelSpacing downSpacing = { 2 * ctSpacing.sx, 2 * ctSpacing.sy, 2 * ctSpacing.sx };
	std::cout << "DownSampled Spacing : (" << downSpacing.sx << ", " << downSpacing.sy << ", " << downSpacing.sz << ")" << std::endl;
	Image_Data downImageData = { h, l, w, imageDown, newOrigin, downSpacing, orientation };
	imageInterpolation3D(croppedImageData, downImageData, NEAREST_NEIGHBOR);
	//==========================================================
	*/
	
	/*
	FILE* path_file;
	//loading_path = inputPath + "paths/path segmentation/centered_path_p1.csv";
	loading_path = inputPath + "paths/automatic key points/centered/path_p2.csv";
	if (fopen_s(&path_file, loading_path.c_str(), "r") != 0) {
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
		current_point = getImageCoordFromRealCoord3D(current_point, newOrigin, intSpacing, orientation);
		path_points.push_back(current_point);
		maskSegment[(size_t)current_point.z][x_new((size_t)current_point.x, (size_t)current_point.y, length)] = 1.0;
	}
	fclose(path_file);

	Image_Data toDistanceMap = { height, length, width, maskSegment, newOrigin, intSpacing, orientation };
	rouyTourinDistanceMap(toDistanceMap, initialSegment, 0, 0.5, 0.4);

	storing_path = outputPath + "distance_map_p2.raw";
	manageRAWFile3D<dataType>(initialSegment, length, width, height, storing_path.c_str(), STORE_DATA, false);

	for (k = 0; k < height; k++) {
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {
				initialSegment[k][x_new(i, j, length)] = 1.0 / (initialSegment[k][x_new(i, j, length)] + 1.0);
			}
		}
	}

	storing_path = outputPath + "initial_segment_p2.raw";
	manageRAWFile3D<dataType>(initialSegment, length, width, height, storing_path.c_str(), STORE_DATA, false);
	*/
	
	//////generate the initial segment
	//double radius = 25;
	//double offset = 0.5;
	//dataType v = 1.0;
	//vector<double> distance_path_points;
	//for (k = 0; k < height; k++) {
	//	for (i = 0; i < length; i++) {
	//		for (j = 0; j < width; j++) {
	//			Point3D current_point = { i, j, k };
	//			current_point = getRealCoordFromImageCoord3D(current_point, newOrigin, intSpacing, orientation);
	//			dataType minDistance = 100000000;
	//			for (int n = 0; n < path_points.size(); n++) {
	//				Point3D pPoint = getRealCoordFromImageCoord3D(path_points[n], newOrigin, intSpacing, orientation);
	//				double pDistance = getPoint3DDistance(current_point, pPoint);
	//				if (minDistance > pDistance) {
	//					minDistance = pDistance;
	//				}
	//			}
	//			if (minDistance == 0) {
	//				std::cout << "Zero distance at point (" << i << ", " << j << ", " << k << ")" << std::endl;
	//			}
	//			initialSegment[k][x_new(i, j, length)] = 1.0 / (minDistance + 1.0);
	//		}
	//	}
	//}

	//////generate the initial segment
	//double radius = 25;
	//double offset = 0.5;
	//dataType v = 1.0;
	//for (k = 0; k < h; k++) {
	//	for (i = 0; i < l; i++) {
	//		for (j = 0; j < w; j++) {
	//			Point3D current_point = { i, j, k };
	//			current_point = getRealCoordFromImageCoord3D(current_point, newOrigin, downSpacing, orientation);
	//			dataType minDistance = 100000000;
	//			for (int n = 0; n < path_points.size(); n++) {
	//				double pDistance = getPoint3DDistance(current_point, path_points[n]);
	//				if (minDistance > pDistance) {
	//					minDistance = pDistance;
	//				}
	//			}
	//			initial_segment[k][x_new(i, j, l)] = 1.0 / (minDistance + v);
	//		}
	//	}
	//}

	//FILE* path_value;
	//storing_path = outputPath + "initial_value_path_pts_2.csv";
	//if (fopen_s(&path_value, storing_path.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//for (int n = 0; n < path_points.size(); n++) {
	//	//Point3D pPoint = getImageCoordFromRealCoord3D(path_points[n], newOrigin, intSpacing, orientation);
	//	fprintf(path_value, "%d,%f\n", n, initialSegment[(size_t)path_points[n].z][x_new((size_t)path_points[n].x, (size_t)path_points[n].y, length)]);
	//}
	//fclose(path_value);

	//==================== GSUBSURF ===========================

	/*
	Segmentation_Parameters segParameters = {
		30,//Maximum number of Gauss-Seidel iterations
		100000,//edge detector coef
		0.000001,//epsilon is the regularization factor (Evans-Spruck)
		5000,//Number of current time step
		5000,//Maximum number of time step
		10,//saving frequency
		0.00000001,//segmentation tolerance
		2 * intSpacing.sx,//tau
		intSpacing.sx,//h
		1.5,//omega_c
		0.001,//tolerance
		1.0,//convection coef
		0.005,//diffusion coef
	};

	Filter_Parameters smoothParameters = {
		0.25 * intSpacing.sx * intSpacing.sx,//tau
		intSpacing.sx,//h
		0.0,//sigma
		0,//K--> we use heat implicit
		1.5,//omega
		0.001,//tolerance
		0.0001,//eps 2
		0.001,//coef
		1,//p
		5,//time step number
		100,//max solver iteration
	};


	Image_Data segmentInputData = { height, length, width, imageCrop, newOrigin, intSpacing, orientation };
	//Image_Data segmentInputData = { h, l, w, imageDown, newOrigin, downSpacing, orientation };
	
	storing_path = outputPath + "segmentation/3D image/p2/31-03-2025/";
	double cpu_start = clock();
	generalizedSubsurfSegmentation(segmentInputData, initialSegment, segParameters, smoothParameters, (unsigned char*)storing_path.c_str());
	//generalizedSubsurfSegmentation(segmentInputData, initial_segment, segParameters, smoothParameters, segCenters, no_of_centers, (unsigned char*)storing_path.c_str());
	double cpu_end = clock();
	double needed_cpu_time = ((double)(cpu_end - cpu_start)) / CLOCKS_PER_SEC;
	printf("Hundred iterations need : %.3f s", needed_cpu_time);
	
	
	//for (k = 0; k < h; k++) {
	//	delete[] imageDown[k];
	//	delete[] initial_segment[k];
	//	delete[] maskPath[k];
	//}
	//delete[] imageDown;
	//delete[] initial_segment;
	//delete[] maskPath;
	
	
	for (k = 0; k < height; k++) {
		delete[] imageCrop[k];
		delete[] initialSegment[k];
		delete[] maskSegment[k];
	}
	delete[] imageCrop;
	delete[] initialSegment;
	delete[] maskSegment;
	*/
	
	//for (k = 0; k < hauteur; k++) {
	//	delete[] interpolatedImage[k];
	//	delete[] maskInterpolation[k];
	//}
	//delete[] interpolatedImage;
	//delete[] maskInterpolation;
	//
	//for (k = 0; k < Height; k++) {
	//	delete[] imageData[k];
	//	delete[] maskAorta[k];
	//}
	//delete[] imageData;
	//delete[] maskAorta;
	//free(ctContainer);

	//==================== Compute Hausdoff distance and Ratio ===========================
	
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

	//==================== Test Potential function ========================================
	
	/*
	dataType** inputImageData = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		inputImageData[k] = new dataType[dim2D]{ 0 };
	}
	loading_path = inputPath + "raw/filtered/New/filtered_p1.raw";
	//manageRAWFile3D<dataType>(inputImageData, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	const size_t height = (size_t)((ctSpacing.sz / ctSpacing.sx) * Height);
	dataType** imageData = new dataType * [height];
	dataType** action = new dataType * [height];
	dataType** potential = new dataType * [height];
	for (k = 0; k < height; k++) {
		imageData[k] = new dataType[dim2D]{ 0 };
		action[k] = new dataType[dim2D]{ 0 };
		potential[k] = new dataType[dim2D]{ 0 };
	}

	Image_Data inputImageStr = { Height, Length, Width, inputImageData, ctOrigin, ctSpacing, orientation };
	VoxelSpacing intSpacing = { ctSpacing.sx, ctSpacing.sy, ctSpacing.sx };
	Image_Data imageStr = { height, Length, Width, imageData, ctOrigin, intSpacing, orientation };
	//imageInterpolation3D(inputImageStr, imageStr, TRILINEAR);
	////imageInterpolation3D(inputImageStr, imageStr, NEAREST_NEIGHBOR);
	storing_path = outputPath + "interpolated_p1.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, height, storing_path.c_str(), STORE_DATA, false);

	////real image p1
	//Point3D seed1 = { 261, 257, 311 };
	//Point3D seed2 = { 257, 250, 535 };

	//Patient 1
	Point3D seed1 = { 262, 258, 146 };
	seed1 = getRealCoordFromImageCoord3D(seed1, ctOrigin, ctSpacing, orientation);
	seed1 = getImageCoordFromRealCoord3D(seed1, ctOrigin, intSpacing, orientation);
	Point3D seed2 = { 266, 256, 245 };
	seed2 = getRealCoordFromImageCoord3D(seed2, ctOrigin, ctSpacing, orientation);
	seed2 = getImageCoordFromRealCoord3D(seed2, ctOrigin, intSpacing, orientation);

	////Patient 2
	//Point3D seed1 = { 257, 254, 249 };
	//seed1 = getRealCoordFromImageCoord3D(seed1, ctOrigin, ctSpacing, orientation);
	//seed1 = getImageCoordFromRealCoord3D(seed1, ctOrigin, intSpacing, orientation);
	//Point3D seed2 = { 257, 243, 350 };
	//seed2 = getRealCoordFromImageCoord3D(seed2, ctOrigin, ctSpacing, orientation);
	//seed2 = getImageCoordFromRealCoord3D(seed2, ctOrigin, intSpacing, orientation);

	////Patient 3
	//Point3D seed1 = { 268, 230, 116 };
	//seed1 = getRealCoordFromImageCoord3D(seed1, ctOrigin, ctSpacing, orientation);
	//seed1 = getImageCoordFromRealCoord3D(seed1, ctOrigin, intSpacing, orientation);
	//Point3D seed2 = { 266, 221, 218 };
	//seed2 = getRealCoordFromImageCoord3D(seed2, ctOrigin, ctSpacing, orientation);
	//seed2 = getImageCoordFromRealCoord3D(seed2, ctOrigin, intSpacing, orientation);
	
	////Patient 4
	//Point3D seed1 = { 280, 229, 135 };
	//seed1 = getRealCoordFromImageCoord3D(seed1, ctOrigin, ctSpacing, orientation);
	//seed1 = getImageCoordFromRealCoord3D(seed1, ctOrigin, intSpacing, orientation);
	//Point3D seed2 = { 285, 234, 220 };
	//seed2 = getRealCoordFromImageCoord3D(seed2, ctOrigin, ctSpacing, orientation);
	//seed2 = getImageCoordFromRealCoord3D(seed2, ctOrigin, intSpacing, orientation);

	////Patient 5
	//Point3D seed1 = { 265, 243, 471 };
	//seed1 = getRealCoordFromImageCoord3D(seed1, ctOrigin, ctSpacing, orientation);
	//seed1 = getImageCoordFromRealCoord3D(seed1, ctOrigin, intSpacing, orientation);
	//Point3D seed2 = { 239, 225, 612 };
	//seed2 = getRealCoordFromImageCoord3D(seed2, ctOrigin, ctSpacing, orientation);
	//seed2 = getImageCoordFromRealCoord3D(seed2, ctOrigin, intSpacing, orientation);

	////Patient 6
	//Point3D seed1 = { 250, 298, 258 };
	//seed1 = getRealCoordFromImageCoord3D(seed1, ctOrigin, ctSpacing, orientation);
	//seed1 = getImageCoordFromRealCoord3D(seed1, ctOrigin, intSpacing, orientation);
	//Point3D seed2 = { 268, 288, 433 };
	//seed2 = getRealCoordFromImageCoord3D(seed2, ctOrigin, ctSpacing, orientation);
	//seed2 = getImageCoordFromRealCoord3D(seed2, ctOrigin, intSpacing, orientation);

	Point3D* endPoints = new Point3D[2];
	endPoints[0] = seed1;
	endPoints[1] = seed2;

	double radius = 3.0;
	Potential_Parameters parameters{
		1000, //edge detector coefficient
		0.15, //threshold
		0.000001,//epsilon
		radius
	};
	//compute3DPotential(imageStr, potential, endPoints, parameters);

	storing_path = outputPath + "potential_p1.raw";
	manageRAWFile3D<dataType>(potential, Length, Width, height, storing_path.c_str(), LOAD_DATA, false);

	vector<Point3D> key_points;
	const double LengthKeyPoints = 50;
	
	Image_Data actionMapStr = { height, Length, Width, action, ctOrigin, intSpacing, orientation };
	storing_path = outputPath + "partial front/p1/action_";
	partialFrontPropagation(actionMapStr, potential, endPoints, storing_path);
	//storing_path = outputPath + "key points/p1/action_";
	//frontPropagationWithKeyPointDetection(actionMapStr, potential, endPoints, LengthKeyPoints, key_points, storing_path);

	////Save the keys points in files
	string saving_csv = outputPath + "seed_p6.csv";
	key_points.push_back(seed1);
	key_points.push_back(seed2);
	//string saving_csv = outputPath + "key_points_p1.csv";
	FILE* f_key_point;
	if (fopen_s(&f_key_point, saving_csv.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	for (int n = 0; n < key_points.size(); n++) {
		key_points[n] = getRealCoordFromImageCoord3D(key_points[n], ctOrigin, intSpacing, orientation);
		fprintf(f_key_point, "%f,%f,%f\n", key_points[n].x, key_points[n].y, key_points[n].z);
	}
	fclose(f_key_point);
	key_points.clear();
	
	delete[] endPoints;
	for (k = 0; k < height; k++) {
		if (k < Height) {
			delete[] inputImageData[k];
		}
		delete[] imageData[k];
		delete[] action[k];
		delete[] potential[k];
	}
	delete[] inputImageData;
	delete[] imageData;
	delete[] action;
	delete[] potential;
	free(ctContainer);
	*/
	
	
	//Artificial image
	const size_t Height = 40, Width = 100, Length = 100;
	const size_t dim2D = Length * Width;
	dataType** imageData = new dataType * [Height];
	dataType** action = new dataType * [Height];
	dataType** potential = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{ 0 };
		action[k] = new dataType[dim2D]{ 0 };
		potential[k] = new dataType[dim2D]{ 0 };
	}
	//loading_path = inputPath + "shape/spiral/spiral_radius5.raw";
	//loading_path = inputPath + "shape/spiral/spiral.raw";
	loading_path = inputPath + "shape/tube3D.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	////create 3D tubular tube
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			Point3D pCenter = { i, j, k };
	//			if (i == j) {
	//				for(size_t in = 0; in < Length; in++) {
	//					for (size_t jn = 0; jn < Width; jn++) {
	//						Point3D current_point = { in, jn, k };
	//						double pDistance = getPoint3DDistance(current_point, pCenter);
	//						if (pDistance < 8) {
	//							imageData[k][x_new(in, jn, Length)] = 1.0;
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}
	//}
	//storing_path = outputPath + "tube.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);
	
	////input 1
	//Point3D seed1 = { 82, 54, 8 };
	//Point3D seed2 = { 81, 45, 85 };

	////input 2
	//Point3D seed1 = { 81, 51, 15 };
	//Point3D seed2 = { 81, 47, 88 };

	//diagonal tube
	Point3D seed1 = { 89, 90, 24 };
	Point3D seed2 = { 9, 9, 16 };

	Point3D* endPoints = new Point3D[2];
	endPoints[0] = seed1;
	endPoints[1] = seed2;

	double radius = 3.0;
	Potential_Parameters parameters{
		1,//1000, //edge detector coefficient
		1,//0.15, //threshold
		0.01,//epsilon
		radius
	};
	Point3D iOrigin = { 0.0, 0.0, 0.0 };
	VoxelSpacing iSpacing = { 1.0, 1.0, 1.0 };
	Image_Data imageStr = { Height, Length, Width, imageData, iOrigin, iSpacing, orientation };
	compute3DPotential(imageStr, potential, endPoints, parameters);

	//storing_path = outputPath + "potential_artificial_img2.raw";
	////manageRAWFile3D<dataType>(potential, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	vector<Point3D> key_points;
	key_points.push_back(seed1);
	key_points.push_back(seed2);
	const double LengthKeyPoints = 10;

	storing_path = outputPath + "partial action/N/action_";
	Image_Data actionMapStr = { Height, Length, Width, action, iOrigin, iSpacing, orientation };
	////frontPropagationWithKeyPointDetection(actionMapStr, potential, endPoints, LengthKeyPoints, key_points, storing_path);
	partialFrontPropagation(actionMapStr, potential, endPoints, storing_path);
	
	//Save the keys points in files
	string saving_csv = outputPath + "endpoints.csv";
	FILE* f_key_point;
	if (fopen_s(&f_key_point, saving_csv.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	for (int n = 0; n < key_points.size(); n++) {
		fprintf(f_key_point, "%f,%f,%f\n", key_points[n].x, key_points[n].y, key_points[n].z);
	}
	fclose(f_key_point);
	key_points.clear();

	delete[] endPoints;
	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
		delete[] action[k];
		delete[] potential[k];
	}
	delete[] imageData;
	delete[] action;
	delete[] potential;
	

    /*
	//2D experiments
	const size_t Length = 100, Width = 100;
	dataType* imageData = new dataType[Length * Width] { 0 };
	dataType* action = new dataType[Length * Width] { 0 };
	dataType* potential = new dataType[Length * Width] { 0 };

	Point2D center = { 50.0, 50.0 };

	////Tube
	//for (i = 0; i < Length; i++) {
	//	for (j = 0; j < Width; j++) {
	//		Point2D pPoint = { i, j };
	//		if (i == j) 
	//		{
	//			for (size_t in = 0; in < Length; in++) {
	//				for (size_t jn = 0; jn < Width; jn++) {
	//					Point2D current_point = { in, jn };
	//					double pDistance = getPoint2DDistance(current_point, pPoint);
	//					if (pDistance < 8) {
	//						imageData[x_new(in, jn, Length)] = 1.0;
	//					}
	//				}
	//			}
	//		}
	//	}
	//}

	////create one pixel missing edges in tube
	//for (i = 0; i < Length; i++) {
	//	for (j = 0; j < Width; j++) {
	//		if (j >= 50 && j <= 50) {
	//			imageData[x_new(i, j, Length)] = 0.0;
	//		}
	//	}
	//}
	
	for (i = 0; i < Length; i++) {
		for (j = (Width / 2); j < Width; j++) {
			Point2D current_point = { i, j };
			double pDistance = getPoint2DDistance(current_point, center);
			if (pDistance > 25.0 && pDistance < 35) {
				imageData[x_new(i, j - 15, Length)] = 1.0;
			}
		}
	}

	for (i = 0; i < Length; i++) {
		for (j = 0; j < Width; j++) {
			if (i >= 50 && i <= 53) {
				imageData[x_new(i, j, Length)] = 0.0;
			}
		}
	}

	////storing_path = outputPath + "object_with_holes_four_pixels.raw";
	//storing_path = outputPath + "semi_circle_with_four_pixels_holes.raw";
	//manageRAWFile2D<dataType>(imageData, Length, Width, storing_path.c_str(), STORE_DATA, false);

	Point2D* endPoints = new Point2D[2];
	endPoints[0] = { 20.0, 40.0 };
	endPoints[1] = { 82.0, 41.0 };
	PixelSpacing pSpacing = { 1.0, 1.0 };
	Point2D origin = { 0.0, 0.0 };
	Image_Data2D imageDataStr = { Length, Width, imageData, origin, pSpacing, NULL };

	////Tube
	//Point2D* endPoints = new Point2D[2];
	//endPoints[0] = { 5.0, 7.0 };
	//endPoints[1] = {97.0, 90.0};// { 99.0, 94.0 };
	//PixelSpacing pSpacing = { 1.0, 1.0 };
	//Point2D origin = { 0.0, 0.0 };
	//Image_Data2D imageDataStr = { Length, Width, imageData, origin, pSpacing, NULL };

	////Save the end points in files
	//FILE* p_end;
	//string end_points = outputPath + "end_points_semi_circle.csv";
	//if (fopen_s(&p_end, end_points.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(p_end, "%f,%f\n", endPoints[0].x, endPoints[0].y);
	//fprintf(p_end, "%f,%f\n", endPoints[1].x, endPoints[1].y);
	//fclose(p_end);

	computePotential(imageDataStr, potential, endPoints);
	partialFrontPropagation2D(imageData, action, potential, Length, Width, endPoints);
	
	////storing_path = outputPath + "action_tube_v1.raw";
	////manageRAWFile2D<dataType>(action, Length, Width, storing_path.c_str(), STORE_DATA, false);

	//delete[] endPoints;
	delete[] imageData;
	delete[] action;
	delete[] potential;
	*/

	//==================== Aorta bifurcation detection ===================================
	
	/*
	const size_t hauteur = (size_t)((ctSpacing.sz / ctSpacing.sx)* Height);
	
	dataType** imageData = new dataType * [Height];
	dataType** inputImageData = new dataType * [hauteur];
	for (k = 0; k < hauteur; k++) {
		if (k < Height) {
			imageData[k] = new dataType[dim2D]{ 0 };
		}
		inputImageData[k] = new dataType[dim2D]{ 0 };
	}

	loading_path = inputPath + "raw/filtered/New/filtered_p1.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	VoxelSpacing intSpacing = { ctSpacing.sx, ctSpacing.sy, ctSpacing.sx };
	Image_Data interpolate = { hauteur, Length, Width, inputImageData, ctOrigin, intSpacing, orientation };
	Image_Data inputD = { Height, Length, Width, imageData, ctOrigin, ctSpacing, orientation };
	
	////imageInterpolation3D(inputD, interpolate, NEAREST_NEIGHBOR);
	//imageInterpolation3D(inputD, interpolate, TRILINEAR);
	storing_path = outputPath + "interpolated_p1.raw";
	manageRAWFile3D<dataType>(inputImageData, Length, Width, hauteur, storing_path.c_str(), LOAD_DATA, false);

	dataType* imageSlice = new dataType[dim2D]{ 0 };
	dataType* foundCirclePtr = new dataType[dim2D]{ 0 };
	Point2D pOrigin = { 0.0, 0.0 };
	PixelSpacing pSpacing = {ctSpacing.sx, ctSpacing.sy};
	OrientationMatrix2D pOrientation = { {1.0, 0.0}, { 0.0, 1.0 } };
	Image_Data2D sliceData = { Length, Width, imageSlice, pOrigin, pSpacing, pOrientation};
	dataType imageBackground = 0.0;
	dataType foreGroundValue = 1.0;
	
	HoughParameters hParameters =
	{
		7.0,   //minimal radius
		15.0,  //maximal radius
		0.5,   //radius step
		0.5,   //epsilon
		40,   //offset 
		1000,  //K edge detector coefficient
		0.2,  //threshold edge detector
	};

	//filtering parameters
	Filter_Parameters filter_parameters
	{
		0.25,//tau
		ctSpacing.sx,//h
		0.0,//sigma
		0,//K--> we use heat implicit
		1.5,//omega
		0.001,//tolerance
		0.000001,//eps 2
		0.000001,//coef
		1,//p
		10,//time step number, p1-->10
		100,//max solver iteration
	};

	PixelSpacing spacing = { ctSpacing.sx, ctSpacing.sy };
	Image_Data2D IMAGE
	{
		Length,
		Width,
		imageSlice,
		pOrigin,
		spacing
	};

	Point2D found_point = { 0.0, 0.0 };
	Point2D liver_centroid = { 256.0, 256.0 };

	std::string path_found_circles;
	std::string path_hough = outputPath + "Hough Transform 02-05/";

	////Save the points in clusters wise
	//FILE* p_centers;
	//string accum_centers = path_hough + "accumulated_ceneters.csv";
	//if (fopen_s(&p_centers, accum_centers.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}

	////Single slice
	//size_t k_slice = 450;
	//copyDataToAnother2dArray(inputImageData[k_slice], imageSlice, Length, Width);
	//storing_path = outputPath + "Hough Transform 02-05/threshold.raw";
	//path_found_circles = outputPath + "Hough Transform 02-05/found_centers.csv";
	//found_point = localCircleDetection(sliceData, foundCirclePtr, liver_centroid, hParameters, storing_path, path_found_circles);
	//storing_path = outputPath + "Hough Transform 02-05/found_circles.raw";
	//manageRAWFile2D<dataType>(foundCirclePtr, Length, Width, storing_path.c_str(), STORE_DATA, false);
	
	//p1 : 127 to 282 , interpolated 279 to 611
	//p5 " 453 to 690
	//p1, interpolated, 339 to 491 (liver)
	//k = 425, 633

	size_t k_slice = 393;
	copyDataToAnother2dArray(inputImageData[k_slice], imageSlice, Length, Width);
	storing_path = outputPath + "Hough Transform 02-05/found_s_390.csv";
	found_point = localCircleDetectionWithOptimization(sliceData, foundCirclePtr, liver_centroid, hParameters, storing_path);
	
	
	//dataType x, y;
	//for (k = 390; k <= 400; k++) {
	//	
	//	Point2D liver_centroid = { 256, 256 };

	//	copyDataToAnother2dArray(inputImageData[k], imageSlice, Length, Width);

	//	//rescaleNewRange2D(imageSlice, Length, Width, 0.0, 1.0);
	//	//heatImplicit2dScheme(IMAGE, filter_parameters);
	//	////copyDataToAnother2dArray(liverData[k], liverSlice, Length, Width);

	//	extension = to_string(k);
	//	path_found_circles = outputPath + "Hough Transform 02-05/centers/found_centers" + extension + ".csv";
	//	storing_path = outputPath + "Hough Transform 02-05/threshold/threshold_" + extension + ".raw";
	//	found_point = localCircleDetection(sliceData, foundCirclePtr, liver_centroid, hParameters, storing_path, path_found_circles);

	//	FILE* current_p_centers;
	//	if (fopen_s(&current_p_centers, path_found_circles.c_str(), "r") != 0) {
	//		printf("Enable to open");
	//		return false;
	//	}
	//	while (feof(current_p_centers) == 0) {
	//		fscanf_s(current_p_centers, "%f", &x);
	//		fscanf_s(current_p_centers, ",");
	//		fscanf_s(current_p_centers, "%f", &y);
	//		fscanf_s(current_p_centers, "\n");
	//		
	//		fprintf(p_centers, "%f,%f\n", x, y);
	//	}
	//	fclose(current_p_centers);

	//	////save slice and bounding box
	//	//BoundingBox2D box = findBoundingBox2D(liver_centroid, Length, Width, hParameters.radius_max, hParameters.offset);
	//	//for (i = box.i_min; i <= box.i_max; i++) {
	//	//	for (j = box.j_min; j <= box.j_max; j++) {
	//	//		if (i == box.i_min || i == box.i_max || j == box.j_min || j == box.j_max) {
	//	//			imageSlice[x_new(i, j, Length)] = 1.0;
	//	//		}
	//	//	}
	//	//}

	//	//storing_path = path_hough + "input/slice_" + extension + ".raw";
	//	//manageRAWFile2D<dataType>(imageSlice, Length, Width, storing_path.c_str(), STORE_DATA, false);

	//}
	//fclose(p_centers);
	
	
	delete[] imageSlice;
	delete[] foundCirclePtr;
	
	for (k = 0; k < hauteur; k++) {
		if (k < Height) {
			delete[] imageData[k];
		}
		delete inputImageData[k];
	}
	delete[] imageData;
	delete[] inputImageData;
	free(ctContainer);
	*/

	//==================== Filter out points with no neighbors in previous and next slice ==============

	/*
	string current_slice;// = outputPath + "test_slice/centers_slice_378.csv";
	string previous_slice;// = outputPath + "test_slice/centers_slice_377.csv";
	string next_slice;// = outputPath + "test_slice/centers_slice_379.csv";

	std::vector<Point2D> points_current, points_previous, points_next;
	dataType x = 0, y = 0;

	//One slice
	for (k = 340; k <= 490; k++) {
		
		extension = to_string(k);
		current_slice = outputPath + "test_slice/centers_slice_" + extension + ".csv";
		previous_slice = outputPath + "test_slice/centers_slice_" + to_string(k - 1) + ".csv";
		next_slice = outputPath + "test_slice/centers_slice_" + to_string(k + 1) + ".csv";

		FILE* current_p_centers;
		if (fopen_s(&current_p_centers, current_slice.c_str(), "r") != 0) {
			printf("Enable to open");
			return false;
		}
		while (feof(current_p_centers) == 0) {
			fscanf_s(current_p_centers, "%f", &x);
			fscanf_s(current_p_centers, ",");
			fscanf_s(current_p_centers, "%f", &y);
			fscanf_s(current_p_centers, "\n");
			Point2D p = { x, y };
			points_current.push_back(p);
		}
		fclose(current_p_centers);

		FILE* previous_p_centers;
		if (fopen_s(&previous_p_centers, previous_slice.c_str(), "r") != 0) {
			printf("Enable to open");
			return false;
		}
		while (feof(previous_p_centers) == 0) {
			fscanf_s(previous_p_centers, "%f", &x);
			fscanf_s(previous_p_centers, ",");
			fscanf_s(previous_p_centers, "%f", &y);
			fscanf_s(previous_p_centers, "\n");
			Point2D p = { x, y };
			points_previous.push_back(p);
		}
		fclose(previous_p_centers);

		FILE* next_p_centers;
		if (fopen_s(&next_p_centers, next_slice.c_str(), "r") != 0) {
			printf("Enable to open");
			return false;
		}
		while (feof(next_p_centers) == 0) {
			fscanf_s(next_p_centers, "%f", &x);
			fscanf_s(next_p_centers, ",");
			fscanf_s(next_p_centers, "%f", &y);
			fscanf_s(next_p_centers, "\n");
			Point2D p = { x, y };
			points_next.push_back(p);
		}
		fclose(next_p_centers);

		size_t nb_previous = 0, nb_next = 0;
		size_t ind_previous = 0, ind_next = 0;
		Point2D pNorth = { 0.0, 0.0 }, pSouth = { 0.0, 0.0 }, pEast = { 0.0, 0.0 }, pWest = { 0.0, 0.0 };
		Point2D pNorthEast = { 0.0, 0.0 }, pNorthWest = { 0.0, 0.0 }, pSouthEast = { 0.0, 0.0 }, pSouthWest = { 0.0, 0.0 };
		Point2D pCurrent = { 0.0, 0.0 };

		//Save the remaining points in current slice
		FILE* pFiltered_centers;
		string filtered_centers = outputPath + "filtered_one_slice/filtered_centers_" + extension + ".csv";
		if (fopen_s(&pFiltered_centers, filtered_centers.c_str(), "w") != 0) {
			printf("Enable to open");
			return false;
		}

		for (size_t it = 0; it < points_current.size(); it++) {

			i = (size_t)points_current[it].x;
			j = (size_t)points_current[it].y;
			pCurrent.x = i; 
			pCurrent.y = j;
			nb_next = 0; 
			nb_previous = 0;

			if (i > 0 && i < Length - 1 && j > 0 && j < Width - 1) {

				pNorth.x = i; 
				pNorth.y = j + 1;
				
				pSouth.x = i; 
				pSouth.y = j - 1;
				
				pEast.x = i + 1; 
				pEast.y = j;
				
				pWest.x = i - 1; 
				pWest.y = j;
				
				pNorthEast.x = i + 1; 
				pNorthEast.y = j + 1;
				
				pNorthWest.x = i - 1; 
				pNorthWest.y = j + 1;
				
				pSouthEast.x = i + 1; 
				pSouthEast.y = j - 1;
				
				pSouthWest.x = i - 1; 
				pSouthWest.y = j - 1;

			}

			//Check if the current point has neighboors in next slice
			for (ind_next = 0; ind_next < points_next.size(); ind_next++) {
				if (points_next[ind_next].x == pCurrent.x && points_next[ind_next].y == pCurrent.y) {
					nb_next++;
				}
				if (points_next[ind_next].x == pNorth.x && points_next[ind_next].y == pNorth.y) {
					nb_next++;
				}
				if (points_next[ind_next].x == pSouth.x && points_next[ind_next].y == pSouth.y) {
					nb_next++;
				}
				if (points_next[ind_next].x == pEast.x && points_next[ind_next].y == pEast.y) {
					nb_next++;
				}
				if (points_next[ind_next].x == pWest.x && points_next[ind_next].y == pWest.y) {
					nb_next++;
				}
				if (points_next[ind_next].x == pNorthEast.x && points_next[ind_next].y == pNorthEast.y) {
					nb_next++;
				}
				if (points_next[ind_next].x == pNorthWest.x && points_next[ind_next].y == pNorthWest.y) {
					nb_next++;
				}
				if (points_next[ind_next].x == pSouthEast.x && points_next[ind_next].y == pSouthEast.y) {
					nb_next++;
				}
				if (points_next[ind_next].x == pSouthWest.x && points_next[ind_next].y == pSouthWest.y) {
					nb_next++;
				}
			}

			//Chech if the current point has neighboors in previous slice
			for (ind_previous = 0; ind_previous < points_previous.size(); ind_previous++) {
				if (points_previous[ind_previous].x == pCurrent.x && points_previous[ind_previous].y == pCurrent.y) {
					nb_previous++;
				}
				if (points_previous[ind_previous].x == pNorth.x && points_previous[ind_previous].y == pNorth.y) {
					nb_previous++;
				}
				if (points_previous[ind_previous].x == pSouth.x && points_previous[ind_previous].y == pSouth.y) {
					nb_previous++;
				}
				if (points_previous[ind_previous].x == pEast.x && points_previous[ind_previous].y == pEast.y) {
					nb_previous++;
				}
				if (points_previous[ind_previous].x == pWest.x && points_previous[ind_previous].y == pWest.y) {
					nb_previous++;
				}
				if (points_previous[ind_previous].x == pNorthEast.x && points_previous[ind_previous].y == pNorthEast.y) {
					nb_previous++;
				}
				if (points_previous[ind_previous].x == pNorthWest.x && points_previous[ind_previous].y == pNorthWest.y) {
					nb_previous++;
				}
				if (points_previous[ind_previous].x == pSouthEast.x && points_previous[ind_previous].y == pSouthEast.y) {
					nb_previous++;
				}
				if (points_previous[ind_previous].x == pSouthWest.x && points_previous[ind_previous].y == pSouthWest.y) {
					nb_previous++;
				}
			}

			if (nb_previous != 0 && nb_next != 0) {
				fprintf(pFiltered_centers, "%f,%f\n", points_current[it].x, points_current[it].y);
			}

		}

		fclose(pFiltered_centers);

		while (points_current.size() > 0) {
			points_current.pop_back();
		}
		while (points_previous.size() > 0) {
			points_previous.pop_back();
		}
		while (points_next.size() > 0) {
			points_next.pop_back();
		}

	}
	*/
	
	/*
	const size_t Length = 512, Width = 512;
	dataType x, y, ratio, radius;

	//five slices
	std::vector<Point2D> points_current;
	std::vector<Point2D> points_previous_1, points_previous_2, points_previous_3, points_previous_4, points_previous_5;
	std::vector<Point2D> points_next_1, points_next_2, points_next_3, points_next_4, points_next_5;
	std::vector<dataType> ratios, radiuses;

	dataType max_ratio = 0, max_radius = 0, initial_max = 0, initial_max_radius = 0;
	Point2D max_point = { 0.0, 0.0 }, initial_point_max = {0.0, 0.0};
	//string to_slice = outputPath + "filter_hough_points/f_centers/center_ratio_radius_";
	string to_slice = outputPath + "filter_hough_29_04/";
	for (k = 425; k <= 633; k++) {

		extension = to_string(k);
		string current_slice = to_slice + "found_centers/center_ratio_radius_" + extension + ".csv";
		
		
		string previous_slice_5 = to_slice + to_string(k - 5) + ".csv";
		string previous_slice_4 = to_slice + to_string(k - 4) + ".csv";
		string previous_slice_3 = to_slice + to_string(k - 3) + ".csv";
		string previous_slice_2 = to_slice + to_string(k - 2) + ".csv";
		string previous_slice_1 = to_slice + to_string(k - 1) + ".csv";
		string next_slice_1 = to_slice + to_string(k + 1) + ".csv";
		string next_slice_2 = to_slice + to_string(k + 2) + ".csv";
		string next_slice_3 = to_slice + to_string(k + 3) + ".csv";
		string next_slice_4 = to_slice + to_string(k + 4) + ".csv";
		string next_slice_5 = to_slice + to_string(k + 5) + ".csv";
		

		FILE* current_p_centers;
		if (fopen_s(&current_p_centers, current_slice.c_str(), "r") != 0) {
			printf("Enable to open");
			return false;
		}
		while (feof(current_p_centers) == 0) {
			fscanf_s(current_p_centers, "%f", &x);
			fscanf_s(current_p_centers, ",");
			fscanf_s(current_p_centers, "%f", &y);
			fscanf_s(current_p_centers, ",");
			fscanf_s(current_p_centers, "%f", &ratio);
			fscanf_s(current_p_centers, ",");
			fscanf_s(current_p_centers, "%f", &radius);
			fscanf_s(current_p_centers, "\n");
			ratios.push_back(ratio);
			radiuses.push_back(radius);
			Point2D p = { x, y };
			points_current.push_back(p);
		}
		fclose(current_p_centers);
		
		
		FILE* previous_p_centers_5;
		if (fopen_s(&previous_p_centers_5, previous_slice_5.c_str(), "r") != 0) {
			printf("Enable to open");
			return false;
		}
		while (feof(previous_p_centers_5) == 0) {
			fscanf_s(previous_p_centers_5, "%f", &x);
			fscanf_s(previous_p_centers_5, ",");
			fscanf_s(previous_p_centers_5, "%f", &y);
			fscanf_s(previous_p_centers_5, ",");
			fscanf_s(previous_p_centers_5, "%f", &ratio);
			fscanf_s(previous_p_centers_5, ",");
			fscanf_s(previous_p_centers_5, "%f", &radius);
			fscanf_s(previous_p_centers_5, "\n");
			Point2D p = { x, y };
			points_previous_5.push_back(p);
		}
		fclose(previous_p_centers_5);

		FILE* previous_p_centers_4;
		if (fopen_s(&previous_p_centers_4, previous_slice_4.c_str(), "r") != 0) {
			printf("Enable to open");
			return false;
		}
		while (feof(previous_p_centers_4) == 0) {
			fscanf_s(previous_p_centers_4, "%f", &x);
			fscanf_s(previous_p_centers_4, ",");
			fscanf_s(previous_p_centers_4, "%f", &y);
			fscanf_s(previous_p_centers_4, ",");
			fscanf_s(previous_p_centers_4, "%f", &ratio);
			fscanf_s(previous_p_centers_4, ",");
			fscanf_s(previous_p_centers_4, "%f", &radius);
			fscanf_s(previous_p_centers_4, "\n");
			Point2D p = { x, y };
			points_previous_4.push_back(p);
		}
		fclose(previous_p_centers_4);

		FILE* previous_p_centers_3;
		if (fopen_s(&previous_p_centers_3, previous_slice_3.c_str(), "r") != 0) {
			printf("Enable to open");
			return false;
		}
		while (feof(previous_p_centers_3) == 0) {
			fscanf_s(previous_p_centers_3, "%f", &x);
			fscanf_s(previous_p_centers_3, ",");
			fscanf_s(previous_p_centers_3, "%f", &y);
			fscanf_s(previous_p_centers_3, ",");
			fscanf_s(previous_p_centers_3, "%f", &ratio);
			fscanf_s(previous_p_centers_3, ",");
			fscanf_s(previous_p_centers_3, "%f", &radius);
			fscanf_s(previous_p_centers_3, "\n");
			Point2D p = { x, y };
			points_previous_3.push_back(p);
		}
		fclose(previous_p_centers_3);

		FILE* previous_p_centers_2;
		if (fopen_s(&previous_p_centers_2, previous_slice_2.c_str(), "r") != 0) {
			printf("Enable to open");
			return false;
		}
		while (feof(previous_p_centers_2) == 0) {
			fscanf_s(previous_p_centers_2, "%f", &x);
			fscanf_s(previous_p_centers_2, ",");
			fscanf_s(previous_p_centers_2, "%f", &y);
			fscanf_s(previous_p_centers_2, ",");
			fscanf_s(previous_p_centers_2, "%f", &ratio);
			fscanf_s(previous_p_centers_2, ",");
			fscanf_s(previous_p_centers_2, "%f", &radius);
			fscanf_s(previous_p_centers_2, "\n");
			Point2D p = { x, y };
			points_previous_2.push_back(p);
		}
		fclose(previous_p_centers_2);

		FILE* previous_p_centers_1;
		if (fopen_s(&previous_p_centers_1, previous_slice_1.c_str(), "r") != 0) {
			printf("Enable to open");
			return false;
		}
		while (feof(previous_p_centers_1) == 0) {
			fscanf_s(previous_p_centers_1, "%f", &x);
			fscanf_s(previous_p_centers_1, ",");
			fscanf_s(previous_p_centers_1, "%f", &y);
			fscanf_s(previous_p_centers_1, ",");
			fscanf_s(previous_p_centers_1, "%f", &ratio);
			fscanf_s(previous_p_centers_1, ",");
			fscanf_s(previous_p_centers_1, "%f", &radius);
			fscanf_s(previous_p_centers_1, "\n");
			Point2D p = { x, y };
			points_previous_1.push_back(p);
		}
		fclose(previous_p_centers_1);

		FILE* next_p_centers_1;
		if (fopen_s(&next_p_centers_1, next_slice_1.c_str(), "r") != 0) {
			printf("Enable to open");
			return false;
		}
		while (feof(next_p_centers_1) == 0) {
			fscanf_s(next_p_centers_1, "%f", &x);
			fscanf_s(next_p_centers_1, ",");
			fscanf_s(next_p_centers_1, "%f", &y);
			fscanf_s(next_p_centers_1, ",");
			fscanf_s(next_p_centers_1, "%f", &ratio);
			fscanf_s(next_p_centers_1, ",");
			fscanf_s(next_p_centers_1, "%f", &radius);
			fscanf_s(next_p_centers_1, "\n");
			Point2D p = { x, y };
			points_next_1.push_back(p);
		}
		fclose(next_p_centers_1);

		FILE* next_p_centers_2;
		if (fopen_s(&next_p_centers_2, next_slice_2.c_str(), "r") != 0) {
			printf("Enable to open");
			return false;
		}
		while (feof(next_p_centers_2) == 0) {
			fscanf_s(next_p_centers_2, "%f", &x);
			fscanf_s(next_p_centers_2, ",");
			fscanf_s(next_p_centers_2, "%f", &y);
			fscanf_s(next_p_centers_2, ",");
			fscanf_s(next_p_centers_2, "%f", &ratio);
			fscanf_s(next_p_centers_2, ",");
			fscanf_s(next_p_centers_2, "%f", &radius);
			fscanf_s(next_p_centers_2, "\n");
			Point2D p = { x, y };
			points_next_2.push_back(p);
		}
		fclose(next_p_centers_2);

		FILE* next_p_centers_3;
		if (fopen_s(&next_p_centers_3, next_slice_3.c_str(), "r") != 0) {
			printf("Enable to open");
			return false;
		}
		while (feof(next_p_centers_3) == 0) {
			fscanf_s(next_p_centers_3, "%f", &x);
			fscanf_s(next_p_centers_3, ",");
			fscanf_s(next_p_centers_3, "%f", &y);
			fscanf_s(next_p_centers_3, ",");
			fscanf_s(next_p_centers_3, "%f", &ratio);
			fscanf_s(next_p_centers_3, ",");
			fscanf_s(next_p_centers_3, "%f", &radius);
			fscanf_s(next_p_centers_3, "\n");
			Point2D p = { x, y };
			points_next_3.push_back(p);
		}
		fclose(next_p_centers_3);

		FILE* next_p_centers_4;
		if (fopen_s(&next_p_centers_4, next_slice_4.c_str(), "r") != 0) {
			printf("Enable to open");
			return false;
		}
		while (feof(next_p_centers_4) == 0) {
			fscanf_s(next_p_centers_4, "%f", &x);
			fscanf_s(next_p_centers_4, ",");
			fscanf_s(next_p_centers_4, "%f", &y);
			fscanf_s(next_p_centers_4, ",");
			fscanf_s(next_p_centers_4, "%f", &ratio);
			fscanf_s(next_p_centers_4, ",");
			fscanf_s(next_p_centers_4, "%f", &radius);
			fscanf_s(next_p_centers_4, "\n");
			Point2D p = { x, y };
			points_next_4.push_back(p);
		}
		fclose(next_p_centers_4);

		FILE* next_p_centers_5;
		if (fopen_s(&next_p_centers_5, next_slice_5.c_str(), "r") != 0) {
			printf("Enable to open");
			return false;
		}
		while (feof(next_p_centers_5) == 0) {
			fscanf_s(next_p_centers_5, "%f", &x);
			fscanf_s(next_p_centers_5, ",");
			fscanf_s(next_p_centers_5, "%f", &y);
			fscanf_s(next_p_centers_5, ",");
			fscanf_s(next_p_centers_5, "%f", &ratio);
			fscanf_s(next_p_centers_5, ",");
			fscanf_s(next_p_centers_5, "%f", &radius);
			fscanf_s(next_p_centers_5, "\n");
			Point2D p = { x, y };
			points_next_5.push_back(p);
		}
		fclose(next_p_centers_5);
		

		Point2D pNorth = { 0.0, 0.0 }, pSouth = { 0.0, 0.0 }, pEast = { 0.0, 0.0 }, pWest = { 0.0, 0.0 };
		Point2D pNorthEast = { 0.0, 0.0 }, pNorthWest = { 0.0, 0.0 }, pSouthEast = { 0.0, 0.0 }, pSouthWest = { 0.0, 0.0 };
		Point2D pCurrent = { 0.0, 0.0 };
		
		
		//Save the remaining points in current slice
		FILE* pFiltered_centers;
		string filtered_centers = outputPath + "filter_hough_points/filtered/one slice/f_points_" + extension + ".csv";
		if (fopen_s(&pFiltered_centers, filtered_centers.c_str(), "w") != 0) {
			printf("Enable to open");
			return false;
		}
		
		size_t nb_previous_1 = 0, nb_previous_2 = 0, nb_previous_3 = 0, nb_previous_4 = 0, nb_previous_5 = 0;
		size_t nb_next_1 = 0, nb_next_2 = 0, nb_next_3 = 0, nb_next_4 = 0, nb_next_5 = 0;
		size_t ind = 0;

		max_ratio = 0; max_radius = 0;
		initial_max = 0; initial_max_radius = 0;
		for (size_t it = 0; it < points_current.size(); it++) {

			if (initial_max < ratios[it]) {
				initial_max = ratios[it];
				initial_max_radius = radiuses[it];
				initial_point_max = points_current[it];
			}

			i = (size_t)points_current[it].x;
			j = (size_t)points_current[it].y;
			pCurrent.x = i;
			pCurrent.y = j;

			nb_previous_1 = 0; nb_previous_2 = 0; nb_previous_3 = 0; nb_previous_4 = 0; nb_previous_5 = 0;
			nb_next_1 = 0; nb_next_2 = 0; nb_next_3 = 0; nb_next_4 = 0; nb_next_5 = 0;

			//Compute the neigbors in given slice
			if (i > 0 && i < Length - 1 && j > 0 && j < Width - 1) {

				pNorth.x = i;
				pNorth.y = j + 1;

				pSouth.x = i;
				pSouth.y = j - 1;

				pEast.x = i + 1;
				pEast.y = j;

				pWest.x = i - 1;
				pWest.y = j;

				pNorthEast.x = i + 1;
				pNorthEast.y = j + 1;

				pNorthWest.x = i - 1;
				pNorthWest.y = j + 1;

				pSouthEast.x = i + 1;
				pSouthEast.y = j - 1;

				pSouthWest.x = i - 1;
				pSouthWest.y = j - 1;

			}

			//Check neighbors in slices above
			for (ind = 0; ind < points_next_1.size(); ind++) {
				if (points_next_1[ind].x == pCurrent.x && points_next_1[ind].y == pCurrent.y) {
					nb_next_1++;
				}
				if (points_next_1[ind].x == pNorth.x && points_next_1[ind].y == pNorth.y) {
					nb_next_1++;
				}
				if (points_next_1[ind].x == pSouth.x && points_next_1[ind].y == pSouth.y) {
					nb_next_1++;
				}
				if (points_next_1[ind].x == pEast.x && points_next_1[ind].y == pEast.y) {
					nb_next_1++;
				}
				if (points_next_1[ind].x == pWest.x && points_next_1[ind].y == pWest.y) {
					nb_next_1++;
				}
				if (points_next_1[ind].x == pNorthEast.x && points_next_1[ind].y == pNorthEast.y) {
					nb_next_1++;
				}
				if (points_next_1[ind].x == pNorthWest.x && points_next_1[ind].y == pNorthWest.y) {
					nb_next_1++;
				}
				if (points_next_1[ind].x == pSouthEast.x && points_next_1[ind].y == pSouthEast.y) {
					nb_next_1++;
				}
				if (points_next_1[ind].x == pSouthWest.x && points_next_1[ind].y == pSouthWest.y) {
					nb_next_1++;
				}
			}

			for (ind = 0; ind < points_next_2.size(); ind++) {
				if (points_next_2[ind].x == pCurrent.x && points_next_2[ind].y == pCurrent.y) {
					nb_next_2++;
				}
				if (points_next_2[ind].x == pNorth.x && points_next_2[ind].y == pNorth.y) {
					nb_next_2++;
				}
				if (points_next_2[ind].x == pSouth.x && points_next_2[ind].y == pSouth.y) {
					nb_next_2++;
				}
				if (points_next_2[ind].x == pEast.x && points_next_2[ind].y == pEast.y) {
					nb_next_2++;
				}
				if (points_next_2[ind].x == pWest.x && points_next_2[ind].y == pWest.y) {
					nb_next_2++;
				}
				if (points_next_2[ind].x == pNorthEast.x && points_next_2[ind].y == pNorthEast.y) {
					nb_next_2++;
				}
				if (points_next_2[ind].x == pNorthWest.x && points_next_2[ind].y == pNorthWest.y) {
					nb_next_2++;
				}
				if (points_next_2[ind].x == pSouthEast.x && points_next_2[ind].y == pSouthEast.y) {
					nb_next_2++;
				}
				if (points_next_2[ind].x == pSouthWest.x && points_next_2[ind].y == pSouthWest.y) {
					nb_next_2++;
				}
			}

			for (ind = 0; ind < points_next_3.size(); ind++) {
				if (points_next_3[ind].x == pCurrent.x && points_next_3[ind].y == pCurrent.y) {
					nb_next_3++;
				}
				if (points_next_3[ind].x == pNorth.x && points_next_3[ind].y == pNorth.y) {
					nb_next_3++;
				}
				if (points_next_3[ind].x == pSouth.x && points_next_3[ind].y == pSouth.y) {
					nb_next_3++;
				}
				if (points_next_3[ind].x == pEast.x && points_next_3[ind].y == pEast.y) {
					nb_next_3++;
				}
				if (points_next_3[ind].x == pWest.x && points_next_3[ind].y == pWest.y) {
					nb_next_3++;
				}
				if (points_next_3[ind].x == pNorthEast.x && points_next_3[ind].y == pNorthEast.y) {
					nb_next_3++;
				}
				if (points_next_3[ind].x == pNorthWest.x && points_next_3[ind].y == pNorthWest.y) {
					nb_next_3++;
				}
				if (points_next_3[ind].x == pSouthEast.x && points_next_3[ind].y == pSouthEast.y) {
					nb_next_3++;
				}
				if (points_next_3[ind].x == pSouthWest.x && points_next_3[ind].y == pSouthWest.y) {
					nb_next_3++;
				}
			}

			for (ind = 0; ind < points_next_4.size(); ind++) {
				if (points_next_4[ind].x == pCurrent.x && points_next_4[ind].y == pCurrent.y) {
					nb_next_4++;
				}
				if (points_next_4[ind].x == pNorth.x && points_next_4[ind].y == pNorth.y) {
					nb_next_4++;
				}
				if (points_next_4[ind].x == pSouth.x && points_next_4[ind].y == pSouth.y) {
					nb_next_4++;
				}
				if (points_next_4[ind].x == pEast.x && points_next_4[ind].y == pEast.y) {
					nb_next_4++;
				}
				if (points_next_4[ind].x == pWest.x && points_next_4[ind].y == pWest.y) {
					nb_next_4++;
				}
				if (points_next_4[ind].x == pNorthEast.x && points_next_4[ind].y == pNorthEast.y) {
					nb_next_4++;
				}
				if (points_next_4[ind].x == pNorthWest.x && points_next_4[ind].y == pNorthWest.y) {
					nb_next_4++;
				}
				if (points_next_4[ind].x == pSouthEast.x && points_next_4[ind].y == pSouthEast.y) {
					nb_next_4++;
				}
				if (points_next_4[ind].x == pSouthWest.x && points_next_4[ind].y == pSouthWest.y) {
					nb_next_4++;
				}
			}

			for (ind = 0; ind < points_next_5.size(); ind++) {
				if (points_next_5[ind].x == pCurrent.x && points_next_5[ind].y == pCurrent.y) {
					nb_next_5++;
				}
				if (points_next_5[ind].x == pNorth.x && points_next_5[ind].y == pNorth.y) {
					nb_next_5++;
				}
				if (points_next_5[ind].x == pSouth.x && points_next_5[ind].y == pSouth.y) {
					nb_next_5++;
				}
				if (points_next_5[ind].x == pEast.x && points_next_5[ind].y == pEast.y) {
					nb_next_5++;
				}
				if (points_next_5[ind].x == pWest.x && points_next_5[ind].y == pWest.y) {
					nb_next_5++;
				}
				if (points_next_5[ind].x == pNorthEast.x && points_next_5[ind].y == pNorthEast.y) {
					nb_next_5++;
				}
				if (points_next_5[ind].x == pNorthWest.x && points_next_5[ind].y == pNorthWest.y) {
					nb_next_5++;
				}
				if (points_next_5[ind].x == pSouthEast.x && points_next_5[ind].y == pSouthEast.y) {
					nb_next_5++;
				}
				if (points_next_5[ind].x == pSouthWest.x && points_next_5[ind].y == pSouthWest.y) {
					nb_next_5++;
				}
			}

			//Chech neighbors in slices below 
			for (ind = 0; ind < points_previous_1.size(); ind++) {
				if (points_previous_1[ind].x == pCurrent.x && points_previous_1[ind].y == pCurrent.y) {
					nb_previous_1++;
				}
				if (points_previous_1[ind].x == pNorth.x && points_previous_1[ind].y == pNorth.y) {
					nb_previous_1++;
				}
				if (points_previous_1[ind].x == pSouth.x && points_previous_1[ind].y == pSouth.y) {
					nb_previous_1++;
				}
				if (points_previous_1[ind].x == pEast.x && points_previous_1[ind].y == pEast.y) {
					nb_previous_1++;
				}
				if (points_previous_1[ind].x == pWest.x && points_previous_1[ind].y == pWest.y) {
					nb_previous_1++;
				}
				if (points_previous_1[ind].x == pNorthEast.x && points_previous_1[ind].y == pNorthEast.y) {
					nb_previous_1++;
				}
				if (points_previous_1[ind].x == pNorthWest.x && points_previous_1[ind].y == pNorthWest.y) {
					nb_previous_1++;
				}
				if (points_previous_1[ind].x == pSouthEast.x && points_previous_1[ind].y == pSouthEast.y) {
					nb_previous_1++;
				}
				if (points_previous_1[ind].x == pSouthWest.x && points_previous_1[ind].y == pSouthWest.y) {
					nb_previous_1++;
				}
			}

			for (ind = 0; ind < points_previous_2.size(); ind++) {
				if (points_previous_2[ind].x == pCurrent.x && points_previous_2[ind].y == pCurrent.y) {
					nb_previous_2++;
				}
				if (points_previous_2[ind].x == pNorth.x && points_previous_2[ind].y == pNorth.y) {
					nb_previous_2++;
				}
				if (points_previous_2[ind].x == pSouth.x && points_previous_2[ind].y == pSouth.y) {
					nb_previous_2++;
				}
				if (points_previous_2[ind].x == pEast.x && points_previous_2[ind].y == pEast.y) {
					nb_previous_2++;
				}
				if (points_previous_2[ind].x == pWest.x && points_previous_2[ind].y == pWest.y) {
					nb_previous_2++;
				}
				if (points_previous_2[ind].x == pNorthEast.x && points_previous_2[ind].y == pNorthEast.y) {
					nb_previous_2++;
				}
				if (points_previous_2[ind].x == pNorthWest.x && points_previous_2[ind].y == pNorthWest.y) {
					nb_previous_2++;
				}
				if (points_previous_2[ind].x == pSouthEast.x && points_previous_2[ind].y == pSouthEast.y) {
					nb_previous_2++;
				}
				if (points_previous_2[ind].x == pSouthWest.x && points_previous_2[ind].y == pSouthWest.y) {
					nb_previous_2++;
				}
			}

			for (ind = 0; ind < points_previous_3.size(); ind++) {
				if (points_previous_3[ind].x == pCurrent.x && points_previous_3[ind].y == pCurrent.y) {
					nb_previous_3++;
				}
				if (points_previous_3[ind].x == pNorth.x && points_previous_3[ind].y == pNorth.y) {
					nb_previous_3++;
				}
				if (points_previous_3[ind].x == pSouth.x && points_previous_3[ind].y == pSouth.y) {
					nb_previous_3++;
				}
				if (points_previous_3[ind].x == pEast.x && points_previous_3[ind].y == pEast.y) {
					nb_previous_3++;
				}
				if (points_previous_3[ind].x == pWest.x && points_previous_3[ind].y == pWest.y) {
					nb_previous_3++;
				}
				if (points_previous_3[ind].x == pNorthEast.x && points_previous_3[ind].y == pNorthEast.y) {
					nb_previous_3++;
				}
				if (points_previous_3[ind].x == pNorthWest.x && points_previous_3[ind].y == pNorthWest.y) {
					nb_previous_3++;
				}
				if (points_previous_3[ind].x == pSouthEast.x && points_previous_3[ind].y == pSouthEast.y) {
					nb_previous_3++;
				}
				if (points_previous_3[ind].x == pSouthWest.x && points_previous_3[ind].y == pSouthWest.y) {
					nb_previous_3++;
				}
			}

			for (ind = 0; ind < points_previous_4.size(); ind++) {
				if (points_previous_4[ind].x == pCurrent.x && points_previous_4[ind].y == pCurrent.y) {
					nb_previous_4++;
				}
				if (points_previous_4[ind].x == pNorth.x && points_previous_4[ind].y == pNorth.y) {
					nb_previous_4++;
				}
				if (points_previous_4[ind].x == pSouth.x && points_previous_4[ind].y == pSouth.y) {
					nb_previous_4++;
				}
				if (points_previous_4[ind].x == pEast.x && points_previous_4[ind].y == pEast.y) {
					nb_previous_4++;
				}
				if (points_previous_4[ind].x == pWest.x && points_previous_4[ind].y == pWest.y) {
					nb_previous_4++;
				}
				if (points_previous_4[ind].x == pNorthEast.x && points_previous_4[ind].y == pNorthEast.y) {
					nb_previous_4++;
				}
				if (points_previous_4[ind].x == pNorthWest.x && points_previous_4[ind].y == pNorthWest.y) {
					nb_previous_4++;
				}
				if (points_previous_4[ind].x == pSouthEast.x && points_previous_4[ind].y == pSouthEast.y) {
					nb_previous_4++;
				}
				if (points_previous_4[ind].x == pSouthWest.x && points_previous_4[ind].y == pSouthWest.y) {
					nb_previous_4++;
				}
			}

			for (ind = 0; ind < points_previous_5.size(); ind++) {
				if (points_previous_5[ind].x == pCurrent.x && points_previous_5[ind].y == pCurrent.y) {
					nb_previous_5++;
				}
				if (points_previous_5[ind].x == pNorth.x && points_previous_5[ind].y == pNorth.y) {
					nb_previous_5++;
				}
				if (points_previous_5[ind].x == pSouth.x && points_previous_5[ind].y == pSouth.y) {
					nb_previous_5++;
				}
				if (points_previous_5[ind].x == pEast.x && points_previous_5[ind].y == pEast.y) {
					nb_previous_5++;
				}
				if (points_previous_5[ind].x == pWest.x && points_previous_5[ind].y == pWest.y) {
					nb_previous_5++;
				}
				if (points_previous_5[ind].x == pNorthEast.x && points_previous_5[ind].y == pNorthEast.y) {
					nb_previous_5++;
				}
				if (points_previous_5[ind].x == pNorthWest.x && points_previous_5[ind].y == pNorthWest.y) {
					nb_previous_5++;
				}
				if (points_previous_5[ind].x == pSouthEast.x && points_previous_5[ind].y == pSouthEast.y) {
					nb_previous_5++;
				}
				if (points_previous_5[ind].x == pSouthWest.x && points_previous_5[ind].y == pSouthWest.y) {
					nb_previous_5++;
				}
			}

			////Keep the center if it has neighboors in all slices
			//if (nb_next_1 != 0 && nb_next_2 != 0 && nb_next_3 != 0 && nb_next_4 != 0 && nb_next_5 != 0 &&
			//	nb_previous_1 != 0 && nb_previous_2 != 0 && nb_previous_3 != 0 && nb_previous_4 != 0 && nb_previous_5 != 0) {
			//	fprintf(pFiltered_centers, "%f,%f\n", points_current[it].x, points_current[it].y);
			//	if (max_ratio < ratios[it]) {
			//		max_ratio = ratios[it];
			//		max_radius = radiuses[it];
			//		max_point = points_current[it];
			//	}
			//}

			////Keep the center if it has neighboors in all slices
			//if (nb_next_1 != 0 && nb_next_2 != 0 && nb_next_3 != 0 && nb_next_4 != 0 &&
			//	nb_previous_1 != 0 && nb_previous_2 != 0 && nb_previous_3 != 0 && nb_previous_4 != 0) {
			//	fprintf(pFiltered_centers, "%f,%f\n", points_current[it].x, points_current[it].y);
			//	if (max_ratio < ratios[it]) {
			//		max_ratio = ratios[it];
			//		max_radius = radiuses[it];
			//		max_point = points_current[it];
			//	}
			//}

			//Keep the center if it has neighboors in all slices
			if (nb_next_1 != 0 && nb_previous_1 != 0) {
				fprintf(pFiltered_centers, "%f,%f\n", points_current[it].x, points_current[it].y);
				if (max_ratio < ratios[it]) {
					max_ratio = ratios[it];
					max_radius = radiuses[it];
					max_point = points_current[it];
				}
			}

		}
		
		//empty the vectors
		while (points_current.size() > 0) {
			points_current.pop_back();
		}
		while (points_previous_1.size() > 0) {
			points_previous_1.pop_back();
		}
		while (points_previous_2.size() > 0) {
			points_previous_2.pop_back();
		}
		while (points_previous_3.size() > 0) {
			points_previous_3.pop_back();
		}
		while (points_previous_4.size() > 0) {
			points_previous_4.pop_back();
		}
		while (points_previous_5.size() > 0) {
			points_previous_5.pop_back();
		}
		while (points_next_1.size() > 0) {
			points_next_1.pop_back();
		}
		while (points_next_2.size() > 0) {
			points_next_2.pop_back();
		}
		while (points_next_3.size() > 0) {
			points_next_3.pop_back();
		}
		while (points_next_4.size() > 0) {
			points_next_4.pop_back();
		}
		while (points_next_5.size() > 0) {
			points_next_5.pop_back();
		}

		initial_max = 0; initial_max_radius = 0;
		for (size_t it = 0; it < points_current.size(); it++) {
			if (initial_max < ratios[it]) {
				initial_max = ratios[it];
				initial_max_radius = radiuses[it];
				initial_point_max = points_current[it];
			}
		}

		//draw initial circle for current slice
		FILE* pInitial_centers;
		//string found_centers_initial = outputPath + "filter_hough_points/filtered/initial/initial_circle_" + extension + ".csv";
		string found_centers_initial = to_slice + "initial/initial_circle_" + extension + ".csv";
		if (fopen_s(&pInitial_centers, found_centers_initial.c_str(), "w") != 0) {
			printf("Enable to open");
			return false;
		}
		double perimeter = 2 * M_PI * initial_max_radius;
		size_t number_of_circle_points = (size_t)(perimeter / 1.0 + 0.5);
		double step = 2 * M_PI / (double)number_of_circle_points;
		BoundingBox2D box = findBoundingBox2D(initial_point_max, Length, Width, initial_max_radius, 1.0);
		for (size_t np = 0; np < number_of_circle_points; np++)
		{
			double phi = (double)np * step;
			double coordx = 0.5 + initial_point_max.x + initial_max_radius * cos(phi);
			double coordy = 0.5 + initial_point_max.y + initial_max_radius * sin(phi);
			size_t indx, indy;

			if (coordx < 0) {
				indx = 0;
			}
			else {
				if (coordx > Length - 1) {
					indx = Length - 1;
				}
				else {
					indx = (size_t)coordx;
				}
			}

			if (coordy < 0) {
				indy = 0;
			}
			else {
				if (coordy > Width - 1) {
					indy = Width - 1;
				}
				else {
					indy = (size_t)coordy;
				}
			}

			fprintf(pInitial_centers, "%d,%d\n", indx, indy);
		}
		fclose(pInitial_centers);
		
		//draw found circle for current slice
		FILE* pFound_centers;
		string found_centers = outputPath + "filter_hough_points/filtered/one slice/circle_" + extension + ".csv";
		if (fopen_s(&pFound_centers, found_centers.c_str(), "w") != 0) {
			printf("Enable to open");
			return false;
		}
		perimeter = 2 * M_PI * max_radius;
		number_of_circle_points = (size_t)(perimeter / 1.0 + 0.5);
		step = 2 * M_PI / (double)number_of_circle_points;
		box = findBoundingBox2D(max_point, Length, Width, max_radius, 1.0);
		for (size_t np = 0; np < number_of_circle_points; np++)
		{
			double phi = (double)np * step;
			//adding +0.5 make it looks more like circle
			double coordx = 0.5 + max_point.x + max_radius * cos(phi);
			double coordy = 0.5 + max_point.y + max_radius * sin(phi);
			size_t indx, indy;

			if (coordx < 0) {
				indx = 0;
			}
			else {
				if (coordx > Length - 1) {
					indx = Length - 1;
				}
				else {
					indx = (size_t)coordx;
				}
			}

			if (coordy < 0) {
				indy = 0;
			}
			else {
				if (coordy > Width - 1) {
					indy = Width - 1;
				}
				else {
					indy = (size_t)coordy;
				}
			}

			fprintf(pFound_centers, "%d,%d\n", indx, indy);
		}
		
		while (ratios.size() > 0) {
			ratios.pop_back();
		}
		while (radiuses.size() > 0) {
			radiuses.pop_back();
		}	

		fclose(pFound_centers);
		fclose(pFiltered_centers);
		

	}
	*/

	//==================== Hough Transform in Lung region =========================

	/*
	dataType** imageData = new dataType * [Height];
	dataType** lungsData = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{ 0 };
		lungsData[k] = new dataType[dim2D]{ 0 };
	}

	loading_path = inputPath + "raw/filtered/New/filtered_p6.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);
	
	loading_path = inputPath + "raw/lungs/iterative/lungs_p6.raw";
	manageRAWFile3D<dataType>(lungsData, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	//get lungs centroid
	dataType* centroid_lungs = new dataType[3]{ 0 };
	centroidImage(lungsData, centroid_lungs, Height, Length, Width, 0.0);
	Point3D centroid_lungs_3D = { centroid_lungs[0], centroid_lungs[1], centroid_lungs[2] };
	Point2D slice_centroid = { centroid_lungs_3D.x, centroid_lungs_3D.y };
	size_t index_slice = (size_t)centroid_lungs_3D.z;
	centroid_lungs_3D = getRealCoordFromImageCoord3D(centroid_lungs_3D, ctOrigin, ctSpacing, orientation);

	FILE* cLungs;
	string lungs_centroid = outputPath + "lungs_centroid_p6.csv";
	if (fopen_s(&cLungs, lungs_centroid.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	fprintf(cLungs, "%f,%f,%f\n", centroid_lungs_3D.x, centroid_lungs_3D.y, centroid_lungs_3D.z);
	fclose(cLungs);

	dataType* imageSlice = new dataType[dim2D]{ 0 };
	dataType* foundCirclePtr = new dataType[dim2D]{ 0 };

	copyDataToAnother2dArray(imageData[index_slice], imageSlice, Length, Width);

	Point2D sOrigin = { 0.0, 0.0 };
	PixelSpacing sSpacing = { ctSpacing.sx, ctSpacing.sy };
	OrientationMatrix2D orientation2D = { {1.0, 0.0}, {0.0, 1.0 } };
	Image_Data2D imageDataToHoughTransform = { Length, Width, imageSlice, sOrigin, sSpacing, orientation2D};

	HoughParameters hParameters =
	{
		10,   //minimal radius
		20.0,  //maximal radius
		0.25,   //radius step
		0.5,   //epsilon
		60,   //offset 
		1000,  //K edge detector coefficient
		0.15,  //threshold edge detector
	};
	
	storing_path = outputPath + "centers_p6.csv";
	//storing_path = outputPath + "threshold_p6.raw";
	Point2D found_point = localCircleDetection(imageDataToHoughTransform, foundCirclePtr, slice_centroid, hParameters, storing_path);

	//draw search area
	BoundingBox2D box = findBoundingBox2D(slice_centroid, Length, Width, hParameters.radius_max, hParameters.offset);
	for (size_t in = box.i_min; in <= box.i_max; in++) {
		for (size_t jn = box.j_min; jn <= box.j_max; jn++) {
			if (in == box.i_min || in == box.i_max || jn == box.j_min || jn == box.j_max) {
				imageSlice[x_new(in, jn, Length)] = 0.0;
			}
		}
	}

	storing_path = outputPath + "slice_centroid_p6.raw";
	manageRAWFile2D<dataType>(imageSlice, Length, Width, storing_path.c_str(), STORE_DATA, false);

	storing_path = outputPath + "found_circle_p6.raw";
	manageRAWFile2D<dataType>(foundCirclePtr, Length, Width, storing_path.c_str(), STORE_DATA, false);

	////save all the found circles
	//FILE* found_points;
	//storing_path = outputPath + "point.csv";
	//if (fopen_s(&found_points, storing_path.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//}
	//FILE* found_pixels;
	//storing_path = outputPath + "pixels.csv";
	//if (fopen_s(&found_pixels, storing_path.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//}
	
	//Point2D center = { 256.0, 256.0 };
	//dataType radius = 15.5;
	//dataType perimetre = 2 * M_PI * radius;
	//dataType step = 2 * M_PI / perimetre;
	//size_t nb_point = (size_t)(perimetre + 1.0);
	//
	//size_t indx, indy;
	//for (size_t np = 0; np < nb_point; np++)
	//{
	//	dataType phi = (dataType)np * step;
	//	dataType coordx = 0.5 + center.x + radius * cos(phi);
	//	dataType coordy = 0.5 + center.y + radius * sin(phi);
	//
	//	if (coordx < 0) {
	//		indx = 0;
	//	}
	//	else {
	//		if (coordx > Length - 1) {
	//			indx = Length - 1;
	//		}
	//		else {
	//			indx = (size_t)round(coordx);
	//		}
	//	}
	//
	//	if (coordy < 0) {
	//		indy = 0;
	//	}
	//	else {
	//		if (coordy > Width - 1) {
	//			indy = Width - 1;
	//		}
	//		else {
	//			indy = (size_t)round(coordy);
	//		}
	//	}
	//
	//	fprintf(found_points, "%f,%f\n", coordx, coordy);
	//	fprintf(found_pixels, "%d,%d\n", indx, indy);
	//}
	//fclose(found_pixels);
	//fclose(found_points);
	
	delete[] imageSlice;
	delete[] foundCirclePtr;

	for (k = 0; k < Height; k++) {
		//if (k < height) {
		//	delete[] crop_lungs[k];
		//	delete[] distance_map[k];
		//	delete[] image_for_hough[k];
		//}
		delete[] imageData[k];
		delete[] lungsData[k];
	}
	delete[] imageData;
	delete[] lungsData;
	free(ctContainer);
	*/

	//==================== Hough Transform with optimization =====================

	/*
	dataType** imageData = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{ 0 };
	}

	//copyDataToAnotherArray(ctContainer->dataPointer, imageData, Height, Length, Width);
	loading_path = inputPath + "raw/filtered/New/filtered_p1.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	PixelSpacing spacing = { ctSpacing.sx, ctSpacing.sy };
	Point2D pOrigin = { 0.0, 0.0 };

	dataType* imageSlice = new dataType[dim2D]{ 0 };
	dataType* foundCircles = new dataType[dim2D]{ 0 };
	dataType* houghSpacePtr = new dataType[dim2D]{ 0 };

	HoughParameters hParameters =
	{
		8.0,   //minimal radius
		15.0,  //maximal radius
		0.5,   //radius step
		0.5,   //epsilon
		60,   //offset 
		1000,  //K edge detector coefficient
		0.2,  //threshold edge detector
	};

	Image_Data2D sliceData = { Length, Width, imageSlice, pOrigin, spacing };

	//size_t kslice = 242;
	//copyDataToAnother2dArray(imageData[kslice], imageSlice, Length, Width);

	Point2D seed = { 256.0, 256.0 };

	//std::string p_threshold = outputPath + "edge_image_025.raw";
	//std::string p_points = outputPath + "potential_centers_025.csv";
	//Point2D found_point = localCircleDetection(sliceData, foundCircles, seed, hParameters, p_threshold, p_points);
	//Point2D found_point = localHoughTransform(sliceData, seed, houghSpacePtr, foundCircles, hParameters);
	//Point2D found_point = localHoughWithCanny(sliceData, seed, houghSpacePtr, foundCircles, hParameters);

	////Draw bounding box
	BoundingBox2D box = findBoundingBox2D(seed, Length, Width, hParameters.radius_max, hParameters.offset);
	//for (i = box.i_min; i <= box.i_max; i++) {
	//	for (j = box.j_min; j <= box.j_max; j++) {
	//		if (i == box.i_min || i == box.i_max || j == box.j_min || j == box.j_max) {
	//			foundCircles[x_new(i, j, Length)] = 1.0;
	//		}
	//	}
	//}

	//save all the found circles
	FILE* found_points_hough;
	std::string points_test = outputPath + "Test Hough 06-05/centers_p1.csv";
	if (fopen_s(&found_points_hough, points_test.c_str(), "w") != 0) {
		printf("Enable to open");
	}

	for (k = 201; k <= 231; k++) {
		
		extension = to_string(k);
		std::string p_threshold = outputPath + "Test Hough 06-05/threshold/edge_image_" + extension + ".raw";
		std::string p_points = outputPath + "Test Hough 06-05/center/potential_centers_" + extension + ".csv";

		copyDataToAnother2dArray(imageData[k], imageSlice, Length, Width);
		Point2D found_point = localCircleDetection(sliceData, foundCircles, seed, hParameters, p_threshold, p_points);
		
		fprintf(found_points_hough, "%f,%f\n", found_point.x, found_point.y);
		for (i = box.i_min; i <= box.i_max; i++) {
			for (j = box.j_min; j <= box.j_max; j++) {
				if (i == box.i_min || i == box.i_max || j == box.j_min || j == box.j_max) {
					imageSlice[x_new(i, j, Length)] = 1.0;
				}
			}
		}
		storing_path = outputPath + "Test Hough 06-05/slice/slices_" + extension + ".raw";
		manageRAWFile2D<dataType>(imageSlice, Length, Width, storing_path.c_str(), STORE_DATA, false);

	}
	fclose(found_points_hough);

	delete[] imageSlice;
	delete[] foundCircles;

	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
	}
	delete[] imageData;

	free(ctContainer);
	*/

	return EXIT_SUCCESS;
}