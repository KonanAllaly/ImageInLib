#include <iostream>
#include <sstream>  
#include <vector>
#include <string.h> 
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
	
	OrientationMatrix orientation;
	orientation.v1 = { 1.0, 0.0, 0.0 }; 
	orientation.v2 = { 0.0, 1.0, 0.0 }; 
	orientation.v3 = { 0.0, 0.0, 1.0 };

	Vtk_File_Info* ctContainer = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	ctContainer->operation = copyFrom;
	loading_path = inputPath + "vtk/petct/ct/Patient2_ct.vtk";
	readVtkFile(loading_path.c_str(), ctContainer);

	int Height = ctContainer->dimensions[2];
	int Length = ctContainer->dimensions[0];
	int Width = ctContainer->dimensions[1];
	int dim2D = Length * Width;
	std::cout << "CT image dim : " << ctContainer->dimensions[0] << " x " << ctContainer->dimensions[1] << " x " << ctContainer->dimensions[2] << "" << std::endl;

	std::cout << "CT origin : (" << ctContainer->origin[0] << ", " << ctContainer->origin[1] << ", " << ctContainer->origin[2] << ")" << std::endl;
	Point3D ctOrigin = { ctContainer->origin[0], ctContainer->origin[1], ctContainer->origin[2] };
	VoxelSpacing ctSpacing = { ctContainer->spacing[0], ctContainer->spacing[1], ctContainer->spacing[2] };
	std::cout << "CT spacing : (" << ctContainer->spacing[0] << ", " << ctContainer->spacing[1] << ", " << ctContainer->spacing[2] << ")" << std::endl; 

	//Find the minimum and maximum values to perform shiftting
	dataType minData = ctContainer->dataPointer[0][0];
	dataType maxData = ctContainer->dataPointer[0][0];
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length * Width; i++) {
			if (ctContainer->dataPointer[k][i] > maxData) {
				maxData = ctContainer->dataPointer[k][i];
			}
			if (ctContainer->dataPointer[k][i] < minData) {
				minData = ctContainer->dataPointer[k][i];
			}
		}
	}
	std::cout << "Min = " << minData << ", Max = " << maxData << std::endl;

	exit(0);
	
	/*
	Vtk_File_Info* aortaContainer = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	aortaContainer->operation = copyFrom;
	loading_path = inputPath + "vtk/aorta/Aorta_p5.vtk";
	readVtkFile(loading_path.c_str(), aortaContainer);
	aortaContainer->origin[0] = ctContainer->origin[0];
	aortaContainer->origin[1] = ctContainer->origin[1];
	aortaContainer->origin[2] = ctContainer->origin[2];
	aortaContainer->operation = copyTo;
	vtkDataForm dataForm = dta_binary;
	loading_path = inputPath + "vtk/aorta/aorta_p5_real_origin.vtk";
	storeVtkFile(loading_path.c_str(), aortaContainer, dataForm);
	free(aortaContainer);
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
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{ 0 };
	}

	//Shift + copy
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			imageData[k][i] = ctContainer->dataPointer[k][i] - minData;
		}
	}

	maxData -= minData;
	minData -= minData;

	rescaleNewRange(imageData, Length, Width, Height, 0.0, 1.0, minData, maxData);

	Filter_Parameters filter_parameters{
		0.01,//tau
		1.0,//h
		0.0,//sigma
		1000.0,//K : edge detector coefficient
		1.5,//omega_c
		1e-2,//tolerance
		1e-2,//epsilon2
		1e-6,//coef
		1,//p
		1,//number of time step
		100,// max solver iteration
	};

	//Filtering_Parameters filter_parameters {
	//	1.0 / 6.0,//tau
	//	ctSpacing.sx,//hx
	//	ctSpacing.sy,//hy
	//	ctSpacing.sz,//hz
	//	0.0,//sigma
	//	1000.0,//K : edge detector coefficient
	//	1.5,//omega_c
	//	1e-2,//tolerance
	//	1e-2,//epsilon2
	//	1e-6,//coef
	//	1,//number of time step
	//	100,// max solver iteration
	//}; 

	Image_Data ctImageData{
		Height,
		Length,
		Width,
		imageData,
		ctOrigin,
		ctSpacing,
		orientation
	};

	//heatExplicitScheme(ctImageData, filter_parameters);
	//heatExplicitRectangularScheme(ctImageData, filter_parameters);
	//geodesicMeanCurvatureRectangularTimeStep(ctImageData, filter_parameters);
	geodesicMeanCurvatureTimeStep(ctImageData, filter_parameters);
	
	storing_path = outputPath + "filtered3.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	////Smooting parameters : tau = 0.25, h = 1.0, omega = 1.5, tolerance = 0.000001, max gauss seidel iteration = 100
	//Segmentation parameters : tau = 1.2, h = 1.0, omega = 1.5, gauss seidel tolerance = 0.000001, 
	//max gauss seidel iteration = 100, number of time steps = 500, edge detector coef = 100000.0, 
	//coef conv = 1.0, coef dif = 1.0

	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
	}
	delete[] imageData;

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
			maskThreshold[k][i] = ctContainer->dataPointer[k][i];// -minData;
		}
	}
	
	dataType thres_min = -30, thres_max = 190;
	//dataType thres_min = -50, thres_max = 199;
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

	storing_path = outputPath + "liver_p7.raw";
	store3dRawData<dataType>(maskThreshold, Length, Width, Height, storing_path.c_str());

	for (k = 0; k < Height; k++) {
		delete[] maskThreshold[k];
		delete[] labelArray[k];
		delete[] status[k];
	}
	delete[] maskThreshold;
	delete[] labelArray;
	delete[] status;
	*/
	
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

	//======================== Path extraction from three seeds ===============================
	
	/*
	dataType** imageData = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{ 0 };
	}
	if (imageData == NULL)
		return false;
	
	//Load filtered input image
	loading_path = inputPath + "raw/filtered/filteredGMC_p6.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	string test_meta_data = outputPath + "metaDataTest_p6.txt";
	FILE* test_description;
	if (fopen_s(&test_description, test_meta_data.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	fprintf(test_description, "This test is performed on Thursday 16th, January 2025\n");
	fprintf(test_description, "The objective of the test is to extract paths inside the aorta\n");
	fprintf(test_description, "The input data is patient number 6\n");
	fprintf(test_description, "Due to the sensitivity of the approach to noise\n");
	fprintf(test_description, "filtering by GMCF is performed on the input image\n");
	fprintf(test_description, "\nPatient 6: origin = (%f, %f, %f)\n", ctOrigin.x, ctOrigin.y, ctOrigin.z);
	fprintf(test_description, "Voxel sizes = (%f, %f, %f)\n", ctSpacing.sx, ctSpacing.sy, ctSpacing.sz);
	fprintf(test_description, "Dimension : Length = %d , Width = %d and Height = %d\n", Length, Width, Height);
	fprintf(test_description, "\nThe path extraction is performed into 2 steps using 3 input points\n");
	fprintf(test_description, "The current implementation of fast marching use isotropic voxel size\n");

	const size_t height = (size_t)((ctSpacing.sz / ctSpacing.sx) * Height);
	std::cout << "Interpolated Height = " << height << std::endl;
	fprintf(test_description, "Interpolated image dimension: Length = %d, Width = %d, Height = %d\n", Length, Width, height);
	
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
	fprintf(test_description, "We perform interpolation on the zth direction to have cubic voxels\n");
	fprintf(test_description, "Interpolated image Dimension : Length = %d , Width = %d and Height = %d\n", Length, Width, height);
	imageInterpolation3D(ctImageData, IntCtImageData, NEAREST_NEIGHBOR);

	double radius = 3.0;
	Potential_Parameters parameters{
		1000, //edge detector coefficient
		0.15, //threshold
		0.005,//epsilon
		radius
	};
	fprintf(test_description, "\nThe potential function is computed only on the initial points\n");
	fprintf(test_description, "To avoid noise or point location effects, we take the average value in spherical (3 mm of radius) volume around the initial point\n");
	fprintf(test_description, "Potential function parameters : \n");
	fprintf(test_description, "radius = %f define the volume where to compute the seed value\n", radius);
	fprintf(test_description, "epsilon = %f is the path smoothing parameter\n", parameters.eps);
	fprintf(test_description, "thres = %f is the threshold for distance map\n", parameters.eps);
	
	//Points in original image coordinates
	
	//Patient 6
	Point3D seed1 = {250, 300, 255};
	Point3D seed2 = {277, 351, 459};
	Point3D seed3 = {250, 295, 459};

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

	fprintf(test_description, "\nThe points used to extract the current path are set manually\n");
	fprintf(test_description, "The points are defined in original image then interpolated\n");
	fprintf(test_description, "Point 1 : (%f, %f, %f)\n", seed1.x, seed1.y, seed1.z);
	fprintf(test_description, "Point 2 : (%f, %f, %f)\n", seed2.x, seed2.y, seed2.z);
	fprintf(test_description, "Point 3 : (%f, %f, %f)\n", seed3.x, seed3.y, seed3.z);

	//To real world coordinate
	seed1 = getRealCoordFromImageCoord3D(seed1, ctOrigin, ctSpacing, orientation);
	seed2 = getRealCoordFromImageCoord3D(seed2, ctOrigin, ctSpacing, orientation);
	seed3 = getRealCoordFromImageCoord3D(seed3, ctOrigin, ctSpacing, orientation);

	seed1 = getImageCoordFromRealCoord3D(seed1, ctOrigin, intSpacing, orientation);
	seed2 = getImageCoordFromRealCoord3D(seed2, ctOrigin, intSpacing, orientation);
	seed3 = getImageCoordFromRealCoord3D(seed3, ctOrigin, intSpacing, orientation);

	//To interpolated image coordinates

	Point3D* seedPoints = new Point3D[3];
	seedPoints[0] = seed1;
	seedPoints[1] = seed2;
	seedPoints[2] = seed3;

	fclose(test_description);

	findPathTwoSteps(IntCtImageData, seedPoints, parameters);
	
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
	*/

	/*
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
	*/

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
	
	/*
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
	rouyTourinDistanceMap(ctImageData, difference, 1.0, 0.5, 0.4);

	storing_path = outputPath + "distance_map_lungs_p1.raw";
	manageRAWFile3D<dataType>(difference, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

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
	*/

	/*
	//======================= get the lungs without the trachea ==============

	//copy data
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			//copy the mask threshsold for next processings
			maskThreshold[k][i] = maskLungs[k][i];
			maskTrachea[k][i] = maskLungs[k][i];
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
			maskThreshold[k][i] = maskTrachea[k][i] - maskLungs[k][i];
		}
	}
	
	storing_path = outputPath + "p3/difference_p3.raw";
	manageRAWFile3D<dataType>(difference, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	//initialize the arrays
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			maskTrachea[k][i] = 0.0;
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
				maskTrachea[k][i] = 1.0;
			}
			else {
				maskTrachea[k][i] = 0.0;
			}
		}
	}
	
	dilatation3D(maskTrachea, Length, Width, Height, 1.0, 0.0);

	storing_path = outputPath + "p3/trachea_p3.raw";
	manageRAWFile3D<dataType>(maskTrachea, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	delete[] countingArrayThird;

	//=====================================

	for (k = 0; k < Height; k++) {
		delete[] segment[k];
		delete[] status[k];
	}
	delete[] segment;
	delete[] status;


	*/

	//for (k = 0; k < Height; k++) {
	//	delete[] imageData[k];
	//	delete[] maskThreshold[k];
	//	delete[] maskLungs[k];
	//	delete[] segment[k];
	//	delete[] status[k];
	//}
	//delete[] imageData;
	//delete[] maskThreshold;
	//delete[] maskLungs;
	//delete[] segment;
	//delete[] status;
	//free(ctContainer);
	
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
	*/

	/*
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
	*/

	/*
	copyDataToAnother2dArray(ctContainer->dataPointer[k_slice], imageSlice, Length, Width);
	
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

	Point2D point_found = localCircleDetection(IMAGE_HOUGH, foundCirclePtr, centroid_seed, h_parameters);

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
	*/

	/*
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
	
	//======================== Detect carina =============================================
	
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

	//======================== Test new fast marching ==========================
	
	/*
	dataType** imageData = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{ 0 };
	}

	loading_path = inputPath + "raw/filtered/filteredGMC_p1.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	Image_Data imageDataStr
	{
		Height,
		Length,
		Width,
		imageData,
		ctOrigin,
		ctSpacing,
		orientation
	};

	////dataType tau = 0.4;
	////dataType tolerance = 0.5;
	////rouyTourinDistanceMap(imageDataStr, actionMap, 1.0, tolerance, tau);
	////storing_path = outputPath + "distance_map.raw";
	////manageRAWFile3D<dataType>(actionMap, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);
	
	//Interpolation
	const size_t height = 406;//(size_t)((ctSpacing.sz / ctSpacing.sx) * Height);
	dataType** intImageData = new dataType * [height];
	dataType** actionMap = new dataType * [height];
	dataType** potential = new dataType * [height];
	for (k = 0; k < height; k++) {
		intImageData[k] = new dataType[dim2D]{ 0 };
		actionMap[k] = new dataType[dim2D]{ 0 };
		potential[k] = new dataType[dim2D]{ 0 };
	}

	VoxelSpacing intSpacing = { ctSpacing.sx, ctSpacing.sy, ctSpacing.sz };

	//////Patient 1
	//Point3D seed1 = { 262, 255, 147 }; //manual setting
	//Point3D seed2 = { 298, 315, 264 }; //obtained by Hough Transform using the lungs centroid location
	////Point3D seed2 = { 298, 315, 264 }; //manual setting
	////Point3D seed3 = { 255, 258, 264 };

	//In Lungs
	Point3D seed1 = { 267, 243, 315 };//trachea
	//Point3D seed2 = { 339, 326, 215 };//left lungs
	Point3D seed2 = { 232, 335, 215 };//right lungs

	////To real world coordinate
	//seed1 = getRealCoordFromImageCoord3D(seed1, ctOrigin, ctSpacing, orientation);
	//seed2 = getRealCoordFromImageCoord3D(seed2, ctOrigin, ctSpacing, orientation);
	////seed3 = getRealCoordFromImageCoord3D(seed3, ctOrigin, ctSpacing, orientation);

	//seed1 = getImageCoordFromRealCoord3D(seed1, ctOrigin, intSpacing, orientation);
	//seed2 = getImageCoordFromRealCoord3D(seed2, ctOrigin, intSpacing, orientation);
	////seed3 = getImageCoordFromRealCoord3D(seed3, ctOrigin, intSpacing, orientation);

	Point3D* seedPoints = new Point3D[2];
	seedPoints[0] = seed1;
	seedPoints[1] = seed2;
	//seedPoints[2] = seed3;

	Image_Data intImageDataStr
	{
		height,
		Length,
		Width,
		imageData,//intImageData,
		ctOrigin,
		intSpacing,
		orientation
	};

	//imageInterpolation3D(imageDataStr, intImageDataStr, NEAREST_NEIGHBOR);
	//////imageInterpolation3D(imageDataStr, intImageDataStr, TRILINEAR);

	double radius = 3.0;
	Potential_Parameters parameters{
		1000, //edge detector coefficient
		0.15, //threshold
		0.005,//epsilon
		radius
	};
	compute3DPotential(intImageDataStr, potential, seed1, parameters);

	storing_path = outputPath + "potential_lungs.raw";
	manageRAWFile3D<dataType>(potential, Length, Width, height, storing_path.c_str(), STORE_DATA, false);

	fastMarching3dWithSpacing(intImageDataStr, actionMap, potential, seed1, intSpacing);
	//fastMarching3D_N(actionMap, potential, Length, Width, height, seed1);
	
	//storing_path = outputPath + "action.raw";
	//manageRAWFile3D<dataType>(actionMap, Length, Width, height, storing_path.c_str(), STORE_DATA, false);

	vector<Point3D> path_points;
	shortestPath3D(actionMap, Length, Width, height, intSpacing, seedPoints, path_points);

	string saving_csv = outputPath + "path_point_lungs2.csv";
	FILE* file_path;
	if (fopen_s(&file_path, saving_csv.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}

	//copy points to file
	int n = 0;
	for (n = path_points.size() - 1; n > -1; n--) {
		path_points[n] = getRealCoordFromImageCoord3D(path_points[n], ctOrigin, intSpacing, orientation);
		fprintf(file_path, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	}
	while (path_points.size() != 0) {
		path_points.pop_back();
	}

	delete[] seedPoints;

	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
	}
	delete[] imageData;

	for (k = 0; k < height; k++) {
		delete[] intImageData[k];
		delete[] actionMap[k];
		delete[] potential[k];
	}
	delete[] intImageData;
	delete[] actionMap;
	delete[] potential;

	free(ctContainer);
	*/

	//======================== Threshold edge detector for trachea ============
	
	/*
	dataType** imageData = new dataType * [Height];
	dataType** maskThreshold = new dataType * [Height];
	dataType** gradientX = new dataType * [Height];
	dataType** gradientY = new dataType * [Height];
	dataType** gradientZ = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{ 0 };
		maskThreshold[k] = new dataType[dim2D]{ 0 };
		gradientX[k] = new dataType[dim2D]{ 0 };
		gradientY[k] = new dataType[dim2D]{ 0 };
		gradientZ[k] = new dataType[dim2D]{ 0 };
	}

	loading_path = inputPath + "raw/filtered/filteredGMC_p1.raw";
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

	compute3dImageGradient(imageData, gradientX, gradientY, gradientZ, Length, Width, Height, ctSpacing);

	//Compute edge detector and threshold
	dataType thres = 0.05;
	dataType K = 1000;
	dataType norm_of_gradient = 0.0;
	dataType edge_value = 0.0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				xd = x_new(i, j, Length);
				norm_of_gradient = sqrt(gradientX[k][xd] * gradientX[k][xd] + gradientY[k][xd] * gradientY[k][xd] + gradientZ[k][xd] * gradientZ[k][xd]);
				edge_value = gradientFunction(norm_of_gradient, K);
				if (edge_value <= thres) {
					maskThreshold[k][xd] = 1.0;
				}
			}
		}
	}
	storing_path = outputPath + "threshold_edge_p1.raw";
	manageRAWFile3D<dataType>(maskThreshold, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
		delete[] maskThreshold[k];
		delete[] gradientX[k];
		delete[] gradientY[k];
		delete[] gradientZ[k];
	}
	delete[] imageData;
	delete[] maskThreshold;
	delete[] gradientX;
	delete[] gradientY;
	delete[] gradientZ;
	*/

	//========== Find lungs segmentation first point ===============
	
	//Point3D first_point = { 0.0, 0.0, 0.0 };
	//size_t count_point = 0;
	//int k_n;
	//for (k_n = height - 1; k_n > -1; k_n--) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			xd = x_new(i, j, Length);
	//			if (interpolateImageData[k_n][xd] == 1) {
	//				count_point++;
	//				if (count_point == 1) {
	//					first_point.x = i;
	//					first_point.y = j;
	//					first_point.z = k_n;
	//				}
	//			}
	//		}
	//	}
	//}
	//std::cout << "First point : (" << first_point.x << "," << first_point.y << "," << first_point.z << ")" << std::endl;

	//============ Path planning from point with max distance and lungs segmentation first point ===

	
	//================
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

	storing_path = outputPath + "p6/mask_lungs.raw";
	manageRAWFile3D<dataType>(maskLungs, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	//for (size_t it = 1; it <= 5; it++) {
	//	erosion3D(maskLungs, Length, Width, Height, 1.0, 0.0);
	//	//storing_path = outputPath + "erosion_" + to_string(it) + ".raw";
	//	//manageRAWFile3D<dataType>(maskLungs, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);
	//}

	////////Erosion and dilatation to remove all the internal background pixels
	//dilatation3D(maskLungs, Length, Width, Height, 1.0, 0.0);
	//dilatation3D(maskLungs, Length, Width, Height, 1.0, 0.0);
	//dilatation3D(maskLungs, Length, Width, Height, 1.0, 0.0);
	//erosion3D(maskLungs, Length, Width, Height, 1.0, 0.0);
	//erosion3D(maskLungs, Length, Width, Height, 1.0, 0.0);
	//erosion3D(maskLungs, Length, Width, Height, 1.0, 0.0);

	//storing_path = outputPath + "lungs_morphology2.raw";
	//manageRAWFile3D<dataType>(maskLungs, Length, Width, Height, storing_path.c_str(), STORE_DATA, false);

	dataType* centroid_lungs = new dataType[3]{0};
	centroidImage(maskLungs, centroid_lungs, Height, Length, Width, 0.0);

	Point3D pCentroid = { centroid_lungs[0], centroid_lungs[1], centroid_lungs[2] };

	//Draw ball centroid
	Point3D pCentroidReal = getRealCoordFromImageCoord3D(pCentroid, ctOrigin, ctSpacing, orientation);
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


	//===========================================================

	/*
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
	shortestPath3D(actionMap, Length, Width, Height, ctSpacing, seedPoints, path_points);

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
	shortestPath3D(actionMap, Length, Width, Height, ctSpacing, seedPoints, path_points);

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
	*/

	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
		delete[] maskLungs[k];
		delete[] ball_centroid[k];
		delete[] leftLung[k];
		delete[] rightLung[k];
		delete[] ball_left[k];
		delete[] ball_right[k];
	}
	delete[] imageData;
	delete[] maskLungs;
	delete[] ball_centroid;
	delete[] leftLung;
	delete[] rightLung;
	delete[] ball_left;
	delete[] ball_right;

	free(ctContainer);
	
	//=============== Compare fast marching implementation =============

	/*
	dataType** imageData = new dataType * [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{ 0 };
	}

	//loading_path = inputPath + "raw/filtered/filteredGMC_p4.raw";
	//manageRAWFile3D<dataType>(imageData, Length, Width, Height, loading_path.c_str(), LOAD_DATA, false);

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

	size_t height = (size_t)((ctSpacing.sz / ctSpacing.sx) * Height);
	VoxelSpacing intSpacing = { ctSpacing.sx, ctSpacing.sy, ctSpacing.sx };

	dataType** potential = new dataType * [height];
	dataType** actionMap = new dataType * [height];
	dataType** interpolate = new dataType * [height];
	for (k = 0; k < height; k++) {
		potential[k] = new dataType[dim2D]{ 0 };
		actionMap[k] = new dataType[dim2D]{ 0 };
		interpolate[k] = new dataType[dim2D]{ 0 };
	}

	//Image_Data IntImageData
	//{
	//	height,
	//	Length,
	//	Width,
	//	interpolate,
	//	ctOrigin,
	//	intSpacing,
	//	orientation
	//};

	//imageInterpolation3D(inputImageData, IntImageData, NEAREST_NEIGHBOR);
	
	Point3D first_point = { 265, 231, 270 };//Trachea
	first_point = getRealCoordFromImageCoord3D(first_point, ctOrigin, ctSpacing, orientation);
	first_point = getImageCoordFromRealCoord3D(first_point, ctOrigin, intSpacing, orientation);

	Point3D second_point = { 179, 302, 239 };//Left lung
	second_point = getRealCoordFromImageCoord3D(second_point, ctOrigin, ctSpacing, orientation);
	second_point = getImageCoordFromRealCoord3D(second_point, ctOrigin, intSpacing, orientation);

	//double radius = 3.0;
	//Potential_Parameters parameters{
	//	1000, //edge detector coefficient
	//	0.15, //threshold
	//	0.005,//epsilon
	//	radius
	//};
	////compute3DPotential(IntImageData, potential, first_point, parameters);
	////storing_path = outputPath + "potential2.raw";
	////manageRAWFile3D<dataType>(potential, Length, Width, height, storing_path.c_str(), LOAD_DATA, false);

	//fastMarching3dWithSpacing(IntImageData, actionMap, potential, first_point, intSpacing);
	//fastMarching3D_N(actionMap, potential, Length, Width, height, first_point);
	storing_path = outputPath + "action4.raw";
	manageRAWFile3D<dataType>(actionMap, Length, Width, height, storing_path.c_str(), LOAD_DATA, false);

	Point3D* seedPoints = new Point3D[2];
	seedPoints[0] = first_point;
	seedPoints[1] = second_point;

	vector<Point3D> path_points;
	shortestPath3D(actionMap, Length, Width, height, intSpacing, seedPoints, path_points);

	//copy points to file
	string saving_csv = outputPath + "path_points_4.csv";
	FILE* file_left;
	if (fopen_s(&file_left, saving_csv.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	int n = 0;
	for (n = path_points.size() - 1; n > -1; n--) {
		path_points[n] = getRealCoordFromImageCoord3D(path_points[n], ctOrigin, intSpacing, orientation);
		fprintf(file_left, "%f,%f,%f\n", path_points[n].x, path_points[n].y, path_points[n].z);
	}
	fclose(file_left);
	while (path_points.size() > 0) {
		path_points.pop_back();
	}

	delete[] seedPoints;

	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
	}
	delete[] imageData;

	for (k = 0; k < height; k++) {
		delete[] interpolate[k];
		delete[] actionMap[k];
		delete[] potential[k];

	}
	delete[] interpolate;
	delete[] actionMap;
	delete[] potential;

	free(ctContainer);
	*/
	
	return EXIT_SUCCESS;
}