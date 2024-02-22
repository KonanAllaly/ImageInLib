// Source code for localization of the heart (06/02/2024)

#include <iostream>
#include <sstream>  
#include <vector>
#include <string.h> 

#include<template_functions.h>
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
#include"hough_transform.h"
#include"../src/edgedetection.h"

#define thres_min 40 //-700
#define thres_max 60 //225
#define thres -500

using namespace std;

int main() {

	string inputPath = "C:/Users/Konan Allaly/Documents/Tests/input/";
	string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";
	
	string loading_path, storing_path, extension;

	size_t i, j, k, xd;

	//===================== Load 3D patient data (.vtk) ================================

	OrientationMatrix orientation;
	orientation.v1 = { 1.0, 0.0, 0.0 }; orientation.v2 = { 0.0, 1.0, 0.0 }; orientation.v3 = { 0.0, 0.0, 1.0 };

	Vtk_File_Info* ctContainer = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	ctContainer->operation = copyFrom;

	//loading_path = inputPath + "vtk/ct/Patient6_ct.vtk";
	//loading_path = inputPath + "vtk/lungs/cropped_ct_lungs_p6.vtk";
	//loading_path = inputPath + "vtk/lungs/box_lungs_ct_p6.vtk";
	loading_path = inputPath + "vtk/ct/Patient2_ct.vtk";
	readVtkFile(loading_path.c_str(), ctContainer);

	int Height = ctContainer->dimensions[2];
	int Length = ctContainer->dimensions[0];
	int Width = ctContainer->dimensions[1];
	int dim2D = Length * Width;
	cout << "Input image dim : " << ctContainer->dimensions[0] << " x " << ctContainer->dimensions[1] << " x " << ctContainer->dimensions[2] << "" << endl;

	cout << "Image origin : (" << ctContainer->origin[0] << ", " << ctContainer->origin[1] << ", " << ctContainer->origin[2] << ")" << endl;
	Point3D ctOrigin = { ctContainer->origin[0], ctContainer->origin[1], ctContainer->origin[2] };
	VoxelSpacing ctSpacing = { ctContainer->spacing[0], ctContainer->spacing[1], ctContainer->spacing[2] };
	cout << "Image spacing : (" << ctContainer->spacing[0] << ", " << ctContainer->spacing[1] << ", " << ctContainer->spacing[2] << ")" << endl;

	//define new array to load input images
	dataType** imageData = new dataType*[Height];
	dataType** maskThreshold = new dataType*[Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D]{0};
		maskThreshold[k] = new dataType[dim2D]{0};
		if (imageData[k] == NULL || maskThreshold[k] == NULL)
			return false;
	}
	if (imageData == NULL || maskThreshold == NULL)
		return false;

	//find minimum and maximum
	dataType min_data = 100000.0, max_data = 0.0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (ctContainer->dataPointer[k][i] > max_data) {
				max_data = ctContainer->dataPointer[k][i];
			}
			if (ctContainer->dataPointer[k][i] < min_data) {
				min_data = ctContainer->dataPointer[k][i];
			}
		}
	}
	cout << "min data = " << min_data << ", max data = " << max_data << endl;

	//thresholding3dFunctionN(ctContainer->dataPointer, Length, Width, Height, thres, thres, 0.0, 1.0);
	//storing_path = outputPath + "threshold.raw";
	//store3dRawData<dataType>(ctContainer->dataPointer, Length, Width, Height, storing_path.c_str());

	////rescale to positive range
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		//imageData[k][i] = ctContainer->dataPointer[k][i] + (-1) * min_data;
	//		maskThreshold[k][i] = ctContainer->dataPointer[k][i] + (-1) * min_data;
	//		//maskThreshold[k][i] = ctContainer->dataPointer[k][i];
	//	}
	//}

	//Copy to array
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			imageData[k][i] = ctContainer->dataPointer[k][i];
		}
	}

	//======================== Path finding ===========================
	
	//======================== Bounding box soft tissue ===============
	/*
	size_t i_min = Length, i_max = 0, j_min = Width, j_max = 0, k_min = Height, k_max = 0;
	
	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				xd = x_new(i, j, Length);
				if (ctContainer->dataPointer[k][xd] >= thres) {
					//find i_min
					if (i < i_min) {
						i_min = i;
					}
					//find j_min
					if (j < j_min) {
						j_min = j;
					}
					//find k_min
					if (k < k_min) {
						k_min = k;
					}
					//find i_max
					if (i > i_max) {
						i_max = i;
					}
					//find j_max
					if (j > j_max) {
						j_max = j;
					}
					//find k_max
					if (k > k_max) {
						k_max = k;
					}
				}
			}
		}
	}

	cout << "i_min = " << i_min << ", j_min = " << j_min << " , k_min = " << k_min << endl;
	cout << "i_max = " << i_max << ", j_max = " << j_max << " , k_max = " << k_max << endl;

	//size_t length = i_max - i_min + 1;
	//size_t width = j_max - j_min + 1;
	//size_t height = k_max - k_min + 1;

	//Point3D cropped_origin = { i_min, j_min, k_min };
	//cropped_origin = getRealCoordFromImageCoord3D(cropped_origin, ctOrigin, ctSpacing, orientation);
	//std::cout << "cropped image origin : (" << cropped_origin.x << ", " << cropped_origin.x << ", " << cropped_origin.z << ")" << endl;

	//cout << "cropped image dimension : " << length << "x" << width << "x" << height << endl;

	//dataType** imageData = new dataType * [height];
	//for (k = 0; k < height; k++) {
	//	imageData[k] = new dataType[length * width]{ 0 };
	//}

	////cropping
	//size_t kn, in, jn;
	//for (kn = 0, k = k_min; kn < height; kn++, k++) {
	//	for (in = 0, i = i_min; in < length; in++, i++) {
	//		for (jn = 0, j = j_min; jn < width; jn++, j++) {
	//			imageData[kn][x_new(in, jn, length)] = ctContainer->dataPointer[k][x_new(i, j, Length)];
	//		}
	//	}
	//}
	*/
	//======================== Region growing ==========================
	/*
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (maskThreshold[k][i] >= -700 && maskThreshold[k][i] <= -600) {
				maskThreshold[k][i] = 1.0;
			}
			else {
				maskThreshold[k][i] = 0.0;
			}
		}
	}
	storing_path = outputPath + "threshold.raw";
	store3dRawData<dataType>(maskThreshold, Length, Width, Height, storing_path.c_str());
	*/
	//======================== Labelling 3D ===================================================================
	/*
	////Labelling
	//thresholdingOTSU(maskThreshold, Length, Width, Height, 0.0, 1.0);
	//iterativeThreshold(maskThreshold, Length, Width, Height, 0.0, 1.0);
	//storing_path = outputPath + "threshold.raw";
	//store3dRawData<dataType>(imageData, Length, Width, Height, storing_path.c_str());
	thresholding3dFunctionN(imageData, Length, Width, Height, thres, thres, 0.0, 1.0);

	////Remove foreground pixel in soft tissue bounding box
	//copyDataToAnotherArray(imageData, maskThreshold, Height, Length, Width);
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (i < i_min || i > i_max || j < j_min || j > j_max || k < k_min || k > k_max) {
	//				maskThreshold[k][x_new(i, j, Length)] = 0.0;
	//			}
	//		}
	//	}
	//}
	
	//storing_path = outputPath + "threshold_and_box_500.raw";
	//store3dRawData<dataType>(maskThreshold, Length, Width, Height, storing_path.c_str());

	//for (i = 0; i < 1; i++) {
	//	dilatation3dHeighteenNeigbours(imageData, Length, Width, Height, 0.0, 1.0);
	//	//dilatation3dHeighteenNeigbours(imageData, Length, Width, Height, 1.0, 0.0);
	//	//dilatation3dHeighteenNeigbours(imageData, Length, Width, Height, 1.0, 0.0);
	//	//erosion3dHeighteenNeigbours(imageData, Length, Width, Height, 1.0, 0.0);
	//}

	//storing_path = outputPath + "threshold_dilate2.raw";
	//store3dRawData<dataType>(imageData, Length, Width, Height, storing_path.c_str());

	int** labelArray = new int* [Height];
	bool** status = new bool* [Height];
	for (k = 0; k < Height; k++) {
		labelArray[k] = new int[dim2D] {0};
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

	labelling3D(imageData, labelArray, status, Length, Width, Height, 1.0);

	//number of region cells
	int numberOfRegionsCells = 0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (imageData[k][i] == 1.0) {
				numberOfRegionsCells++;
			}
		}
	}
	//cout << "Objects voxels number : " << numberOfRegionsCells << endl;

	//Counting
	int* countingArray = new int[numberOfRegionsCells] {0};
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

	////Number of regions
	//int numberOfRegions = 0;
	//for (i = 0; i < numberOfRegionsCells; i++) {
	//	if (countingArray[i] > 0) {
	//		numberOfRegions++;
	//	}
	//}
	//cout << numberOfRegions << " regions found" << endl;

	////remove regions lying out of soft tissue bounding box
	//int label = 1;
	//size_t in, jn, kn;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (labelArray[k][x_new(i, j, Length)] == label) {
	//				if (i < i_min || i > i_max || j < j_min || j > j_max || k < k_min || k > k_max) {
	//					for (kn = 0; kn < Height; kn++) {
	//						for (in = 0; in < Length; in++) {
	//							for (jn = 0; jn < Width; jn++) {
	//								if (labelArray[kn][x_new(in, jn, Length)] == label) {
	//									labelArray[kn][x_new(in, jn, Length)] = 0;
	//								}
	//							}
	//						}
	//					}
	//					label++;
	//				}
	//			}
	//		}
	//	}
	//}

	////Save regions info
	//storing_path = outputPath + "info_regions.csv";
	//FILE* file;
	//if (fopen_s(&file, storing_path.c_str(), "w") != 0) {
	//	printf("Enable to open");
	//	return false;
	//}
	//fprintf(file, "Region,Size\n");
	//for (i = 0; i < numberOfRegionsCells; i++) {
	//	if (countingArray[i] > 1000) {
	//		fprintf(file, "%d,%d\n", i, countingArray[i]);
	//	}
	//}
	//fclose(file);

	//Find the biggest region
	int max1 = 0, max2 = 0, max3 = 0, max4 = 0;
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (countingArray[labelArray[k][i]] > max1) {
				max1 = countingArray[labelArray[k][i]];
			}
			if (countingArray[labelArray[k][i]] > max2 && countingArray[labelArray[k][i]] < max1) {
				max2 = countingArray[labelArray[k][i]];
			}
			if (countingArray[labelArray[k][i]] > max3 && countingArray[labelArray[k][i]] < max2) {
				max3 = countingArray[labelArray[k][i]];
			}
			if (countingArray[labelArray[k][i]] > max4 && countingArray[labelArray[k][i]] < max3) {
				max4 = countingArray[labelArray[k][i]];
			}
		}
	}
	//cout << "max1 = " << max1 << ", max2 = " << max2 << ", max3 = " << max3 << ", max4 = " << max4 << endl;

	//Keep biggest region
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			if (countingArray[labelArray[k][i]] == max3 || countingArray[labelArray[k][i]] == max2) {
				maskThreshold[k][i] = 1.0;
			}
			else {
				maskThreshold[k][i] = 0.0;
			}
		}
	}

	////delete small regions
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < Length; i++) {
	//		for (j = 0; j < Width; j++) {
	//			if (countingArray[labelArray[k][x_new(i, j, Length)]] < 1000) {
	//				labelArray[k][x_new(i, j, Length)] = 0;
	//			}
	//		}
	//	}
	//}

	//save labeled image
	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			imageData[k][i] = (dataType)labelArray[k][i];
			//if (labelArray[k][i] != 0) {
			//	maskThreshold[k][i] = 1.0;
			//}
			//else {
			//	maskThreshold[k][i] = 0.0;
			//}
		}
	}

	storing_path = outputPath + "segment.raw";
	store3dRawData<dataType>(imageData, Length, Width, Height, storing_path.c_str());

	storing_path = outputPath + "biggest_region.raw";
	store3dRawData<dataType>(maskThreshold, Length, Width, Height, storing_path.c_str());

	////Re-Initialization
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		status[k][i] = false;
	//		labelArray[k][i] = 0;
	//	}
	//}
	//
	//for (i = 0; i < 2; i++) {
	//	dilatation3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	//}
	//for (i = 0; i < 2; i++) {
	//	erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	//}
	//labelling3D(maskThreshold, labelArray, status, Length, Width, Height, 1.0);

	////save labeled image
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		maskThreshold[k][i] = (dataType)labelArray[k][i];
	//	}
	//}

	////Re-initialization
	//for (i = 0; i < numberOfRegionsCells; i++) {
	//	countingArray[i] = 0;
	//}

	////Get region size
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		if (labelArray[k][i] > 0) {
	//			countingArray[labelArray[k][i]]++;
	//		}
	//	}
	//}

	////Find the biggest region
	//int max1 = 0, max2 = 0, max3 = 0, max4 = 0;
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		if (countingArray[labelArray[k][i]] > max1) {
	//			max1 = countingArray[labelArray[k][i]];
	//		}
	//		if (countingArray[labelArray[k][i]] > max2 && countingArray[labelArray[k][i]] < max1) {
	//			max2 = countingArray[labelArray[k][i]];
	//		}
	//		//if (countingArray[labelArray[k][i]] > max3 && countingArray[labelArray[k][i]] < max2) {
	//		//	max3 = countingArray[labelArray[k][i]];
	//		//}
	//		//if (countingArray[labelArray[k][i]] > max4 && countingArray[labelArray[k][i]] < max3) {
	//		//	max4 = countingArray[labelArray[k][i]];
	//		//}
	//	}
	//}

	////Keep biggest region
	//for (k = 0; k < Height; k++) {
	//	for (i = 0; i < dim2D; i++) {
	//		if (countingArray[labelArray[k][i]] == max1) {
	//			maskThreshold[k][i] = 1.0;
	//		}
	//		else {
	//			maskThreshold[k][i] = 0.0;
	//		}
	//	}
	//}

	//for (i = 0; i < 1; i++) {
	//	erosion3dHeighteenNeigbours(maskThreshold, Length, Width, Height, 1.0, 0.0);
	//}

	//storing_path = outputPath + "segment.raw";
	//store3dRawData<dataType>(maskThreshold, Length, Width, Height, storing_path.c_str());

	//delete[] countingArray;
	for (k = 0; k < Height; k++) {
		delete[] labelArray[k];
		delete[] status[k];
	}
	delete[] labelArray;
	delete[] status;
	*/
	//======================== Tests on Trachea ================
	/*
	loading_path = inputPath + "trachea.raw";
	load3dArrayRAW<dataType>(tracheaData, Length, Width, Height, loading_path.c_str(), false);
	
	size_t n_op = 5;
	for (i = 0; i < n_op; i++) {
		dilatation3dHeighteenNeigbours(tracheaData, Length, Width, Height, 1.0, 0.0);
	}
	for (i = 0; i < n_op; i++) {
		erosion3dHeighteenNeigbours(tracheaData, Length, Width, Height, 1.0, 0.0);
	}
	storing_path = outputPath + "trachea_dilated5_eroded5.raw";
	store3dRawData<dataType>(tracheaData, Length, Width, Height, storing_path.c_str());
	
	fastSweepingFunction_3D(distanceMap, tracheaData, Length, Width, Height, 1.0, 10000000.0, 0.0);
	storing_path = outputPath + "distance.raw";
	store3dRawData<dataType>(distanceMap, Length, Width, Height, storing_path.c_str());
	
	storing_path = outputPath + "max_dist_per_slice.csv";
	FILE* file;
	if (fopen_s(&file, storing_path.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}

	dataType max_dist = 0.0, max_dist_slice = 0.0;
	dataType mx = 0.0, my = 0.0, mz = 0.0;
	Point3D point_max = { 0.0, 0.0, 0.0 };
	for (k = 0; k < Height; k++) {
		
		max_dist_slice = 0;
		mx = 0; my = 0; mz = k;
		
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				
				xd = x_new(i, j, Length);
				if (distanceMap[k][xd] > max_dist_slice) {
					max_dist_slice = distanceMap[k][xd];
					mx = i;
					my = j;
				}
			}
		}
		
		if (mx != 0 && my != 0) {
			Point3D point_save = { mx, my, mz };
			point_save = getRealCoordFromImageCoord3D(point_save, ctOrigin, ctSpacing, orientation);
			fprintf(file, "%f,%f,%f,%f\n", point_save.x, point_save.y, point_save.z, max_dist_slice);
		}
	}
	fclose(file);

	point_max = getPointWithTheHighestValue(distanceMap, Length, Width, Height);
	point_max = getRealCoordFromImageCoord3D(point_max, ctOrigin, ctSpacing, orientation);

	for (k = 0; k < Height; k++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				Point3D current_point = { i, j, k };
				current_point = getRealCoordFromImageCoord3D(current_point, ctOrigin, ctSpacing, orientation);
				double dist = getPoint3DDistance(point_max, current_point);
				if (dist <= 3) {
					imageData[k][x_new(i, j, Length)] = 1.0;
				}
				else {
					imageData[k][x_new(i, j, Length)] = 0.0;
				}
			}
		}
	}

	storing_path = outputPath + "ball_max.raw";
	store3dRawData<dataType>(imageData, Length, Width, Height, storing_path.c_str());
	*/
	//======================== Labeling 2D =========================
	/*
	loading_path = inputPath + "trachea_dilated5_eroded5.raw";
	load3dArrayRAW<dataType>(maskThreshold, Length, Width, Height, loading_path.c_str(), false);

	int* labelArray = new int[dim2D];
	bool* status = new bool[dim2D];
	dataType* mask = new dataType[dim2D];
	if (labelArray == NULL || status == NULL || mask == NULL)
		return false;

	size_t n;

	storing_path = outputPath + "info_regions.csv";
	FILE* file;
	if (fopen_s(&file, storing_path.c_str(), "w") != 0) {
		printf("Enable to open");
		return false;
	}
	fprintf(file, "slice number,number of regions\n");

	for (k = 0; k < Height; k++) {

		cout << "Slice " << k << ":" << endl;
		fprintf(file, "%d,", k);
		
		//Initialization
		for (i = 0; i < dim2D; i++) {
			mask[i] = maskThreshold[k][i];
			labelArray[i] = 0;
			status[i] = false;
		}

		labeling(mask, labelArray, status, Length, Width, 1.0);

		//number of region cells
		int numberOfRegionsCells = 0;
		for (i = 0; i < dim2D; i++) {
			if (maskThreshold[k][i] == 1.0) {
				numberOfRegionsCells++;
			}
		}
		cout << "Number of object pixels : " << numberOfRegionsCells << endl;

		if (numberOfRegionsCells != 0) {

			//Get regions sizes
			int* countingArray = new int[numberOfRegionsCells] {0};
			if (countingArray == NULL)
				return false;
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

			cout << numb_region << " region found" << endl;
			fprintf(file, "%d\n", numb_region);

			//Get size of each region
			size_t region_volume;
			label = 1;
			for (n = 0; n < numb_region; n++) {
				region_volume = 0;
				for (i = 0; i < dim2D; i++) {
					if (labelArray[i] == label) {
						region_volume++;
					}
				}
				cout << "Region " << label << " contains " << region_volume << " elements" << endl;
				label++;
			}

			delete[] countingArray;

		}
		else {
			fprintf(file, "%d\n", numberOfRegionsCells);
		}

		cout << "#########################"<< endl;

	}
	fclose(file);
	delete[] labelArray;
	delete[] status;
	*/
	//======================== Slice Cropping ================================================================
	/*
	imageMetaData croppedLungs;
	const size_t offset = 5;
	croppedLungs =  croppImage3D(ctImageData, offset);

	size_t length = croppedLungs.length, width = croppedLungs.width, height = croppedLungs.height;
	std::cout << "Cropped lungs image dim : " << length << " x " << width << " x " << height << "" << std::endl;
	Point3D originLungs = croppedLungs.origin;
	std::cout << "Cropped Lungs origin : (" << originLungs.x << ", " << originLungs.y << ", " << originLungs.z << ")" << std::endl;
	VoxelSpacing lungsSpacing = croppedLungs.spacing;

	dataType** cropped_ct = new dataType*[height];
	dataType** cropped_lungs = new dataType * [height];
	for (k = 0; k < height; k++) {
		cropped_ct[k] = new dataType[dim2D]{0};
		cropped_lungs[k] = new dataType[dim2D]{0};
		if (cropped_ct[k] == NULL || cropped_lungs[k] == NULL)
			return false;
	}
	if (cropped_ct == NULL || cropped_lungs == NULL)
		return false;

	//===
	//Image_Data Cropped;
	//Cropped.imageDataPtr = croppedCT;
	//Cropped.height = height; Cropped.length = length; Cropped.width = width;
	//Cropped.orientation = orientation; Cropped.origin = croppedOrigin; Cropped.spacing = croppedSpacing;
	//===

	Point3D origin_cropped = getImageCoordFromRealCoord3D(originLungs, ctOrigin, ctSpacing, orientation);
	size_t i0 = (size_t)origin_cropped.x, j0 = (size_t)origin_cropped.y, k0 = (size_t)origin_cropped.z;
	
	Point3D point_test = { 0.0, 0.0, k0 };
	point_test = getRealCoordFromImageCoord3D(point_test, ctOrigin, ctSpacing, orientation);
	cout << "Origin cropped image : (" << point_test.x << "," << point_test.y << "," << point_test.z << ")" << endl;
	
	size_t i1, j1, k1;
	for (k = 0, k1 = k0; k < height; k++, k1++) {
		for (i = 0; i < Length; i++) {
			for (j = 0; j < Width; j++) {
				cropped_ct[k][x_new(i, j, Length)] = imageData[k1][x_new(i, j, Length)];
				cropped_lungs[k][x_new(i, j, Length)] = lungsData[k1][x_new(i, j, Length)];
			}
		}
	}

	saving_name = outputPath + "cropped_ct_p6.raw";
	//store3dRawData<dataType>(cropped_ct, Length, Width, height, saving_name.c_str());

	saving_name = outputPath + "cropped_lungs_p6.raw";
	//store3dRawData<dataType>(cropped_lungs, Length, Width, height, saving_name.c_str());

	Vtk_File_Info* Container = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	Container->dataPointer = cropped_lungs; //cropped_ct;
	Container->dimensions[0] = Length; Container->dimensions[1] = Width; Container->dimensions[2] = height;
	Container->origin[0] = point_test.x; Container->origin[1] = point_test.y; Container->origin[2] = point_test.z;
	Container->spacing[0] = ctSpacing.sx; Container->spacing[1] = ctSpacing.sy; Container->spacing[2] = ctSpacing.sz;
	Container->operation = copyTo; Container->vDataType = dta_Flt;

	saving_name = outputPath + "cropped_lungs_p6.vtk";
	vtkDataForm data_form = dta_binary;
	storeVtkFile(saving_name.c_str(), Container, data_form);

	free(Container);
	for (k = 0; k < height; k++) {
		delete[] cropped_ct[k];
		delete[] cropped_lungs[k];
	}
	delete[] cropped_ct;
	delete[] cropped_lungs;
	*/
	//======================== Cropp box ====================================
	/*
	imageMetaData croppedLungs;
	const size_t offset = 5;
	croppedLungs =  croppImage3D(ctImageData, offset);

	size_t length = croppedLungs.length, width = croppedLungs.width, height = croppedLungs.height;
	std::cout << "Cropped lungs image dim : " << length << " x " << width << " x " << height << "" << std::endl;
	Point3D originLungs = croppedLungs.origin;
	std::cout << "Cropped Lungs origin : (" << originLungs.x << ", " << originLungs.y << ", " << originLungs.z << ")" << std::endl;
	VoxelSpacing lungsSpacing = croppedLungs.spacing;

	dataType** cropped_ct = new dataType*[height];
	dataType** cropped_lungs = new dataType * [height];
	for (k = 0; k < height; k++) {
		cropped_ct[k] = new dataType[dim2D]{0};
		cropped_lungs[k] = new dataType[dim2D]{0};
		if (cropped_ct[k] == NULL || cropped_lungs[k] == NULL)
			return false;
	}
	if (cropped_ct == NULL || cropped_lungs == NULL)
		return false;

	Point3D origin_cropped = getImageCoordFromRealCoord3D(originLungs, ctOrigin, ctSpacing, orientation);
	size_t i0 = (size_t)origin_cropped.x, j0 = (size_t)origin_cropped.y, k0 = (size_t)origin_cropped.z;
	
	Point3D point_test = { 0.0, 0.0, 0.0 };
	point_test = getRealCoordFromImageCoord3D(point_test, ctOrigin, ctSpacing, orientation);
	cout << "Origin cropped image : (" << point_test.x << "," << point_test.y << "," << point_test.z << ")" << endl;
	
	size_t i1, j1, k1;
	for (k = 0, k1 = k0; k < height; k++, k1++) {
		for (i = 0, i1 = i0; i < length; i++, i1++) {
			for (j = 0, j1 = j0; j < width; j++, j1++) {
				cropped_ct[k][x_new(i, j, length)] = imageData[k1][x_new(i1, j1, Length)];
				cropped_lungs[k][x_new(i, j, length)] = lungsData[k1][x_new(i1, j1, Length)];
			}
		}
	}

	Vtk_File_Info* Container = (Vtk_File_Info*)malloc(sizeof(Vtk_File_Info));
	Container->dataPointer = cropped_lungs; // cropped_ct; //
	Container->dimensions[0] = length; Container->dimensions[1] = width; Container->dimensions[2] = height;
	Container->origin[0] = originLungs.x; Container->origin[1] = originLungs.y; Container->origin[2] = originLungs.z;
	Container->spacing[0] = ctSpacing.sx; Container->spacing[1] = ctSpacing.sy; Container->spacing[2] = ctSpacing.sz;
	Container->operation = copyTo; Container->vDataType = dta_Flt;

	saving_name = outputPath + "box_lungs_p6.vtk";
	vtkDataForm data_form = dta_binary;
	storeVtkFile(saving_name.c_str(), Container, data_form);

	free(Container);
	for (k = 0; k < height; k++) {
		delete[] cropped_ct[k];
		delete[] cropped_lungs[k];
	}
	delete[] cropped_ct;
	delete[] cropped_lungs;
	*/
	//======================= Hough tranform on whole volume =======================
	/*
	//const size_t Length = 512, Width = 512, dim2D = Length * Width;
	dataType* imageSlice = new dataType[dim2D]{ 0 };
	dataType* foundCircle = new dataType[dim2D]{ 0 };
	dataType* houghSpace = new dataType[dim2D]{ 0 };
	
	//PixelSpacing spacing = { 1.171875, 1.171875 };
	PixelSpacing spacing = { ctSpacing.sx, ctSpacing.sy };
	HoughParameters parameters = { 4.0, 15.0, 0.5, 2.0, 1000, 1.0, 0.15, 0.5, spacing };

	//filtering parameters
	Filter_Parameters filter_parameters;
	filter_parameters.h = 1.0; filter_parameters.timeStepSize = 0.25; filter_parameters.eps2 = 1e-6;
	filter_parameters.omega_c = 1.5; filter_parameters.tolerance = 1e-3; filter_parameters.maxNumberOfSolverIteration = 100;
	filter_parameters.timeStepsNum = 5; filter_parameters.coef = 1e-6;

	Image_Data2D IMAGE;
	IMAGE.imageDataPtr = imageSlice;
	IMAGE.height = Length; IMAGE.width = Width;

	////Slice loading
	//loading_path = inputPath + "raw/slice/Image2.raw";
	//load2dArrayRAW<dataType>(imageSlice, Length, Width, loading_path.c_str(), false);

	for (k = 0; k < Height; k++) {

		extension = to_string(k);

		cout << "Slice " << k << endl;

		for (i = 0; i < dim2D; i++) {
			imageSlice[i] = imageData[k][i];
		}

		//storing_path = outputPath + "loaded.raw";
		//store2dRawData<dataType>(imageSlice, Length, Width, storing_path.c_str());

		rescaleNewRange2D(imageSlice, Length, Width, 0.0, 1.0);

		//storing_path = outputPath + "rescaled.raw";
		//store2dRawData<dataType>(imageSlice, Length, Width, storing_path.c_str());

		heatImplicit2dScheme(IMAGE, filter_parameters);

		storing_path = outputPath + "slice_" + extension + ".raw";
		store2dRawData<dataType>(imageSlice, Length, Width, storing_path.c_str());

		//storing_path = outputPath + "filtered.raw";
		//store2dRawData<dataType>(imageSlice, Length, Width, storing_path.c_str());

		//loading_path = inputPath + "slice_145.raw";
		//load2dArrayRAW<dataType>(imageSlice, Length, Width, loading_path.c_str(), false);
		
		storing_path = outputPath + "threshold_" + extension + ".raw";
		houghTransform(imageSlice, foundCircle, Length, Width, parameters, storing_path);

		storing_path = outputPath + "found_circle_" + extension + ".raw";
		store2dRawData<dataType>(foundCircle, Length, Width, storing_path.c_str());

	}

	delete[] imageSlice;
	delete[] foundCircle;
	delete[] houghSpace;
	*/
	//======================= Hough transform to detect bronchial tree ==========
	/*
	//const size_t Length = 512, Width = 512, dim2D = Length * Width;
	dataType* imageSlice = new dataType[dim2D]{ 0 };
	dataType* maskThreshold = new dataType[dim2D]{ 0 };
	dataType* foundCircle = new dataType[dim2D]{ 0 };
	dataType* gradientX = new dataType[dim2D]{ 0 };
	dataType* gradientY = new dataType[dim2D]{ 0 };
	dataType* gradientAngle = new dataType[dim2D]{ 0 };
	dataType* normGradient = new dataType[dim2D]{ 0 };
	size_t* statusPixel = new size_t[dim2D]{ 0 };

	//PixelSpacing spacing = { 1.171875, 1.171875 };
	PixelSpacing spacing = { ctSpacing.sx, ctSpacing.sy };
	HoughParameters parameters = { 5.0, 15.0, 0.5, 5.0, 1000, 1.0, 0.15, 0.5, spacing };

	//filtering parameters
	Filter_Parameters filter_parameters;
	filter_parameters.h = 1.0; filter_parameters.timeStepSize = 0.25; filter_parameters.eps2 = 1e-6;
	filter_parameters.omega_c = 1.5; filter_parameters.tolerance = 1e-3; filter_parameters.maxNumberOfSolverIteration = 100;
	filter_parameters.timeStepsNum = 5; filter_parameters.coef = 1e-6;

	Image_Data2D IMAGE;
	IMAGE.imageDataPtr = imageSlice;
	IMAGE.height = Length; IMAGE.width = Width;

	////Slice loading
	//loading_path = inputPath + "raw/slice/Image2.raw";
	//load2dArrayRAW<dataType>(imageSlice, Length, Width, loading_path.c_str(), false);

	for (k = 0; k < Height; k++) {
		
		extension = to_string(k);

		cout << "Slice " << k << endl;
		
		//Extract slice
		for (i = 0; i < dim2D; i++) {
			imageSlice[i] = imageData[k][i];
		}
		//storing_path = outputPath + "loaded.raw";
		//store2dRawData<dataType>(imageSlice, Length, Width, storing_path.c_str());

		//Slice rescaling
		rescaleNewRange2D(imageSlice, Length, Width, 0.0, 1.0);
		//storing_path = outputPath + "rescaled.raw";
		//store2dRawData<dataType>(imageSlice, Length, Width, storing_path.c_str());

		//Slice filtering
		heatImplicit2dScheme(IMAGE, filter_parameters);
		//storing_path = outputPath + "filtered.raw";
		//store2dRawData<dataType>(imageSlice, Length, Width, storing_path.c_str());

		copyDataToAnother2dArray(imageSlice, maskThreshold, Length, Width);

		computeImageGradient(maskThreshold, gradientX, gradientY, Length, Width, parameters.h);
		computeAngleFromGradient(gradientAngle, gradientX, gradientY, Length, Width);

		//norm of gradient
		for (i = 0; i < dim2D; i++) {
			normGradient[i] = sqrt(gradientX[i] * gradientX[i] + gradientY[i] * gradientY[i]);
		}

		nonMaximumSuppression(maskThreshold, normGradient, gradientAngle, Length, Width);
		thresholdByHyteresis(maskThreshold, normGradient, statusPixel, Length, Width, 0.005, 0.009);

		//Erosion and Dilation
		erosion2D(maskThreshold, Length, Width, 1.0, 0.0);
		dilatation2D(maskThreshold, Length, Width, 1.0, 0.0);

		storing_path = outputPath + "threshold_" + extension + ".raw";
		store2dRawData<dataType>(maskThreshold, Length, Width, storing_path.c_str());

		houghTransform(maskThreshold, foundCircle, Length, Width, parameters);

		storing_path = outputPath + "found_circle_" + extension + ".raw";
		store2dRawData<dataType>(foundCircle, Length, Width, storing_path.c_str());

		initialize2dArray(foundCircle, Length, Width);
		initialize2dArray(maskThreshold, Length, Width);
		initialize2dArray(normGradient, Length, Width);
		initialize2dArray(gradientAngle, Length, Width);
		initialize2dArray(gradientX, Length, Width);
		initialize2dArray(gradientY, Length, Width);
		for (i = 0; i < dim2D; i++) {
			statusPixel[i] = false;
		}

	}

	//delete[] imageSlice;
	//delete[] maskThreshold;
	//delete[] foundCircle;
	//delete[] gradientX;
	//delete[] gradientY;
	//delete[] gradientAngle;
	//delete[] normGradient;
	//delete[] statusPixel;
	*/
	//======================== Free Memory ======================================

	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
		delete[] maskThreshold[k];
	}
	delete[] imageData;
	delete[] maskThreshold;
	free(ctContainer);

	return EXIT_SUCCESS;
}