#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h> 
#include <memory.h>
#include <math.h>
#include <string>

#include"hough_transform.h"
#include"morphological_change.h"
#include"../src/imageInterpolation.h"

#define M_PI 3.14159265358979323846
#define foreground 1.0
#define background 0.0

void initialize2dArray(dataType* array2D, const size_t length, const size_t width) {
	size_t i, dim2D = length * width;
	for (i = 0; i < dim2D; i++) {
		array2D[i] = 0.0;
	}
}

Point2D getPointWithMaximalValue2D(dataType* arrayPtr2D, const size_t length, const size_t width) {
	size_t i, j, xd;
	Point2D max_point = { 0.0, 0.0};
	dataType max_value = 0.0;
	for (i = 0; i < length; i++) {
		for (j = 0; j < width; j++) {
			xd = x_new(i, j, length);
			if (arrayPtr2D[xd] > max_value) {
				max_value = arrayPtr2D[xd];
				max_point.x = i;
				max_point.y = j;
			}
		}
	}
	return max_point;
}

void circularHoughTransform(dataType* imageDataPtr, dataType* houghSpacePtr, dataType* foundCirclePtr, const size_t length, const size_t width, HoughParameters params) {
	
	size_t i, j, xd, dim2D = length * width;
	size_t i_new, j_new, xd_new, count_vote = 0, count_neighbor = 0;
	BoundingBox2D box = { 0, 0, 0, 0 };
	Point2D center_circle = { 0.0, 0.0 }, current_point = { 0.0, 0.0 }, found_center = {0.0, 0.0};

	double dist_point = 0.0;
	double hough_ratio = 0.0, neighbor_ratio = 0.0, radius = 0.0;
	double small_radius_surface = 0.0, band_surface = 0.0;
	double max_ratio = 0.0, found_ratio = 0.0, found_radius = 0.0;
	bool test_max = false;

	std::string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";
	std::string saving_path;

	dataType* edgeDetector = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* maskThreshold = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* houghSpaceMax = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* gradientX = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* gradientY = (dataType*)malloc(sizeof(dataType) * dim2D);

	initialize2dArray(edgeDetector, length, width);
	initialize2dArray(maskThreshold, length, width);
	initialize2dArray(houghSpaceMax, length, width);
	initialize2dArray(gradientX, length, width);
	initialize2dArray(gradientY, length, width);

	//compute the edge detector
	computeImageGradient(imageDataPtr, gradientX, gradientY, length, width, params.h);
	dataType norm_of_gradient = 0.0;
	for (i = 0; i < dim2D; i++) {
		for (i = 0; i < dim2D; i++) {
			norm_of_gradient = sqrt(gradientX[i] * gradientX[i] + gradientY[i] * gradientY[i]);
			edgeDetector[i] = gradientFunction(norm_of_gradient, params.K);
		}
	}

	saving_path = outputPath + "edge_detector.raw";
	store2dRawData(edgeDetector, length, width, saving_path.c_str());

	//copy edge detector to mask threshold
	copyDataToAnother2dArray(edgeDetector, maskThreshold, length, width);

	//thresholding of the edge detector
	thresholding2DFunction(maskThreshold, length, width, params.thres, params.thres);

	saving_path = outputPath + "threshold.raw";
	store2dRawData(maskThreshold, length, width, saving_path.c_str());

	//erosion2D(maskThreshold, length, width, 1.0, 0.0);

	radius = params.radius_min;
	while (radius <= params.radius_max) {
		//compute the surface of the band of interest
		band_surface = 4 * M_PI * radius * params.epsilon;
		//compute surface of smaller circle
		small_radius_surface = M_PI * (radius - params.epsilon) * (radius - params.epsilon);
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {
				xd = x_new(i, j, length);
				center_circle.x = i;
				center_circle.y = j;
				//find the bounding box for the center of interest
				box = findBoundingBox2D(center_circle, length, width, radius, params.offset);
				//vote
				count_vote = 0;
				count_neighbor = 0;
				if (maskThreshold[xd] == background) {
					for (i_new = box.i_min; i_new <= box.i_max; i_new++) {
						for (j_new = box.j_min; j_new <= box.j_max; j_new++) {
							xd_new = x_new(i_new, j_new, length);
							current_point.x = i_new;
							current_point.y = j_new;
							dist_point = getPoint2DDistance(center_circle, current_point);
							//count foreground pixels in the band
							if (maskThreshold[xd_new] == foreground && dist_point >= (radius - params.epsilon) && dist_point <= (radius + params.epsilon)) {
								count_vote++;
							}
							//count foreground pixels close to the current center
							if (dist_point < (radius - params.epsilon) && maskThreshold[xd_new] == foreground) {
								count_neighbor++;
							}
						}
					}
				}
				//compute ratios
				hough_ratio = (double)count_vote / band_surface;
				if (count_neighbor == 0) {
					neighbor_ratio = (double)count_neighbor / small_radius_surface;
					houghSpacePtr[xd] = hough_ratio;
				}
			}
		}
		//find the maximal ratio for current radius
		max_ratio = getTheMaxValue(houghSpacePtr, length, width);
		//find the maximal ratio for whole radius range
		if (max_ratio > found_ratio) {
			found_ratio = max_ratio;
			found_radius = radius;
			test_max = true;
			copyDataToAnother2dArray(houghSpacePtr, houghSpaceMax, length, width);
		}
		//update the found_center if the value of found_ratio was updated
		if (test_max == true) {
			found_center = getPointWithMaximalValue2D(houghSpaceMax, length, width);
			test_max = false;
		}
		radius += params.radius_step;
	}

	printf("Found radius = %f\n", found_radius);
	printf("Found center : (%f, %f)\n", found_center.x, found_center.y);

	//draw the found circle
	box = findBoundingBox2D(found_center, length, width, found_radius, params.offset);
	for (i = box.i_min; i <= box.i_max; i++) {
		for (j = box.j_min; j <= box.j_max; j++) {
			xd = x_new(i, j, length);
			current_point.x = i;
			current_point.y = j;
			dist_point = getPoint2DDistance(found_center, current_point);
			if (dist_point >= (found_radius - params.epsilon) && dist_point <= (found_radius + params.epsilon) && dist_point == 0.0) {
				foundCirclePtr[xd] = 1.0;
			}
		}
	}

	saving_path = outputPath + "found_circle.raw";
	store2dRawData(foundCirclePtr, length, width, saving_path.c_str());

	//copy back
	copyDataToAnother2dArray(houghSpaceMax, houghSpacePtr, length, width);

	free(edgeDetector);
	free(maskThreshold);
	free(houghSpaceMax);
	free(gradientX);
	free(gradientY);

}

void localHoughTransform(Point2D seed, dataType* imageDataPtr, dataType* houghSpacePtr, dataType* foundCirclePtr, const size_t length, const size_t width, HoughParameters params, std::string savingPath, FILE* saveInfo) {

	size_t i, j, xd, dim2D = length * width;
	size_t i_new, j_new, xd_new;
	BoundingBox2D box = { 0, 0, 0, 0 };
	Point2D center_circle = { 0.0, 0.0 }, current_point = { 0.0, 0.0 }, found_center = { 0.0, 0.0 };

	double dist_point = 0.0, found_radius, radius = 0.0;
	bool test_max = false;
	
	dataType* edgeDetector = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* maskThreshold = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* houghSpaceMax = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* gradientX = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* gradientY = (dataType*)malloc(sizeof(dataType) * dim2D);

	initialize2dArray(gradientX, length, width);
	initialize2dArray(gradientY, length, width);

	//compute the edge detector
	initialize2dArray(edgeDetector, length, width);
	computeImageGradient(imageDataPtr, gradientX, gradientY, length, width, params.h);
	dataType norm_of_gradient = 0.0;
	for (i = 0; i < dim2D; i++) {
		norm_of_gradient = sqrt(gradientX[i] * gradientX[i] + gradientY[i] * gradientY[i]);
		edgeDetector[i] = gradientFunction(norm_of_gradient, params.K);
	}

	//copy edge detector to mask threshold
	initialize2dArray(maskThreshold, length, width);
	copyDataToAnother2dArray(edgeDetector, maskThreshold, length, width);

	//thresholding of the edge detector
	thresholding2DFunction(maskThreshold, length, width, params.thres, params.thres);

	////Erosion
	erosion2D(maskThreshold, length, width, foreground, background);
	dilatation2D(maskThreshold, length, width, foreground, background);
	//erosion2D(maskThreshold, length, width, foreground, background);
	//dilatation2D(maskThreshold, length, width, foreground, background);
	//erosion2D(maskThreshold, length, width, foreground, background);
	//dilatation2D(maskThreshold, length, width, foreground, background);

	//saving_path = outputPath + "threshold.raw";
	store2dRawData(maskThreshold, length, width, savingPath.c_str());

	initialize2dArray(houghSpaceMax, length, width);
	initialize2dArray(foundCirclePtr, length, width);

	size_t count_vote = 0, count_neighbor = 0, total_neighbor = 0;
	dataType hough_ratio = 0.0, found_ratio = 0.0, max_ratio = 0.0;
	size_t i_min = 0, i_max = 0, j_min = 0, j_max = 0;
	radius = params.radius_min;
	found_radius = 0.0;
	while (radius <= params.radius_max) {
		initialize2dArray(houghSpacePtr, length, width);
		//find the bounding box of the seed point
		box = findBoundingBox2D(seed, length, width, radius, params.offset);
		//small bounding box
		i_min = (size_t)(box.i_min + radius + params.epsilon);
		if (box.i_max >= radius) {
			i_max = (size_t)(box.i_max - radius - params.epsilon);
		}
		j_min = (size_t)(box.j_min + radius + params.epsilon);
		if (box.j_max >= radius) {
			j_max = (size_t)(box.j_max - radius - params.epsilon);
		}
		
		//visit all the pixels in the bounding box
		for (i = i_min; i <= i_max; i++) {
			for (j = j_min; j <= j_max; j++) {
				xd = x_new(i, j, length);
				center_circle.x = (dataType)i;
				center_circle.y = (dataType)j;
				
				//vote
				count_vote = 0;
				count_neighbor = 0;
				total_neighbor = 0;
				if (maskThreshold[xd] == background) {
					for (i_new = box.i_min; i_new <= box.i_max; i_new++) {
						for (j_new = box.j_min; j_new <= box.j_max; j_new++) {
							xd_new = x_new(i_new, j_new, length);
							current_point.x = (dataType)i_new;
							current_point.y = (dataType)j_new;
							dist_point = getPoint2DDistance(center_circle, current_point);
							//count foreground pixels in the band
							if (dist_point >= (radius - params.epsilon) && dist_point <= (radius + params.epsilon)) {
								total_neighbor++;
								if (maskThreshold[xd_new] == foreground) {
									count_vote++;
								}
							}
							//count foreground pixels close to the current center
							if (dist_point < (radius - params.epsilon) && maskThreshold[xd_new] == foreground) {
								count_neighbor++;
							}

						}
					}
				}
				//compute ratios
				hough_ratio = (dataType)count_vote / (dataType)total_neighbor;
				//houghSpacePtr[xd] = hough_ratio;
				if (count_neighbor == 0) {
					houghSpacePtr[xd] = hough_ratio;
				}
			}
		}

		//keep the maximal ratio for each pixel
		for (i = box.i_min; i <= box.i_max; i++) {
			for (j = box.j_min; j <= box.j_max; j++) {
				xd = x_new(i, j, length);
				if (houghSpacePtr[xd] > houghSpaceMax[xd]) {
					houghSpaceMax[xd] = houghSpacePtr[xd];
				}
			}
		}

		//find the maximal ratio for current radius
		max_ratio = getTheMaxValue(houghSpacePtr, length, width);
		//find the maximal ratio for whole radius range
		if (max_ratio > found_ratio) {
			found_ratio = max_ratio;
			found_radius = radius;
			test_max = true;
			//copyDataToAnother2dArray(houghSpacePtr, houghSpaceMax, length, width);
		}
		//update the found_center if the value of found_ratio was updated
		if (test_max == true) {
			found_center = getPointWithMaximalValue2D(houghSpaceMax, length, width);
			test_max = false;
		}
		radius += params.radius_step;
	}

	printf("Found radius = %lf\n", found_radius);
	//printf("Found center : (%f, %f)\n", found_center.x, found_center.y);

	//add the seed
	foundCirclePtr[x_new((size_t)seed.x, (size_t)seed.y, length)] = 1.0;

	//draw maximal bounding box
	double radius_max_box = params.radius_max + params.epsilon;
	box = findBoundingBox2D(seed, length, width, radius_max_box, params.offset);
	for (i = box.i_min; i <= box.i_max; i++) {
		for (j = box.j_min; j <= box.j_max; j++) {
			xd = x_new(i, j, length);
			if (i == box.i_min || i == box.i_max || j == box.j_min || j == box.j_max) {
				foundCirclePtr[xd] = 1.0;
			}
		}
	}

	if (found_center.x != 0 && found_center.y != 0) {
		//add found center
		foundCirclePtr[x_new((size_t)found_center.x, (size_t)found_center.y, length)] = 1.0;
		//draw the found circle
		box = findBoundingBox2D(found_center, length, width, found_radius, params.offset);
		for (i = box.i_min; i <= box.i_max; i++) {
			for (j = box.j_min; j <= box.j_max; j++) {
				xd = x_new(i, j, length);
				current_point.x = (dataType)i;
				current_point.y = (dataType)j;
				dist_point = getPoint2DDistance(found_center, current_point);
				if (dist_point >= (found_radius - params.epsilon) && dist_point <= (found_radius + params.epsilon)) {
					foundCirclePtr[xd] = 1.0;
				}
			}
		}
	}

	//copy back
	copyDataToAnother2dArray(houghSpaceMax, houghSpacePtr, length, width);

	max_ratio = getTheMaxValue(houghSpaceMax, length, width);
	fprintf(saveInfo, "%f,%f\n", max_ratio, found_radius);

	free(edgeDetector);
	free(maskThreshold);
	free(houghSpaceMax);
	free(gradientX);
	free(gradientY);
}

Point2D localHoughWithCanny(Point2D seed, dataType* imageDataPtr, dataType* houghSpacePtr, dataType* foundCirclePtr, const size_t length, const size_t width, HoughParameters params, std::string savingPath, FILE* saveInfo) {

	size_t i, j, xd, dim2D = length * width;
	size_t i_new, j_new, xd_new;
	BoundingBox2D box = { 0, 0, 0, 0 };
	Point2D center_circle = { 0.0, 0.0 }, current_point = { 0.0, 0.0 }, found_center = { 0.0, 0.0 };

	double dist_point = 0.0, found_radius, radius = 0.0;
	bool test_max = false;

	Point2D originImage = { 0.0, 0.0 };
	OrientationMatrix2D orientation;
	orientation.v1 = { 1.0, 0.0 }; orientation.v2 = { 0.0, 1.0 };

	dataType* maskThreshold = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* houghSpaceMax = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* gradientX = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* gradientY = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* gradientAngle = (dataType*)malloc(sizeof(dataType) * dim2D);
	dataType* normGradient = (dataType*)malloc(sizeof(dataType) * dim2D);
	size_t* statusPixel = (size_t*)malloc(sizeof(size_t) * dim2D);

	initialize2dArray(maskThreshold, length, width);
	initialize2dArray(houghSpaceMax, length, width);
	initialize2dArray(gradientX, length, width);
	initialize2dArray(gradientY, length, width);
	initialize2dArray(gradientAngle, length, width);
	initialize2dArray(normGradient, length, width);
	initialize2dArray(foundCirclePtr, length, width);

	computeImageGradient(imageDataPtr, gradientX, gradientY, length, width, params.h);
	computeAngleFromGradient(gradientAngle, gradientX, gradientY, length, width);
	
	//norm of gradient
	for (i = 0; i < dim2D; i++) {
		normGradient[i] = sqrt(gradientX[i] * gradientX[i] + gradientY[i] * gradientY[i]);
	}
	
	nonMaximumSuppression(maskThreshold, normGradient, gradientAngle, length, width);
	thresholdByHyteresis(maskThreshold, normGradient, statusPixel, length, width, 0.005, 0.009);

	////Erosion
	erosion2D(maskThreshold, length, width, foreground, background);
	dilatation2D(maskThreshold, length, width, foreground, background);

	//saving_path = outputPath + "threshold.raw";
	store2dRawData(maskThreshold, length, width, savingPath.c_str());

	size_t count_vote = 0, count_neighbor = 0, total_neighbor = 0;
	dataType hough_ratio = 0.0, found_ratio = 0.0, max_ratio = 0.0;
	size_t i_min = 0, i_max = 0, j_min = 0, j_max = 0;
	radius = params.radius_min;
	found_radius = 0.0;

	while (radius <= params.radius_max) {
		//initialize the 2D array
		initialize2dArray(houghSpacePtr, length, width);
		//find the bounding box of the seed point
		box = findBoundingBox2D(seed, length, width, radius, params.offset);
		//small bounding box
		i_min = (size_t)(box.i_min + radius + params.epsilon);
		if (box.i_max >= radius) {
			i_max = (size_t)(box.i_max - radius - params.epsilon);
		}
		j_min = (size_t)(box.j_min + radius + params.epsilon);
		if (box.j_max >= radius) {
			j_max = (size_t)(box.j_max - radius - params.epsilon);
		}

		//visit all the pixels in the bounding box
		for (i = i_min; i <= i_max; i++) {
			for (j = j_min; j <= j_max; j++) {
				xd = x_new(i, j, length);

				center_circle.x = (dataType)i;
				center_circle.y = (dataType)j;
				getRealCoordFromImageCoord2D(center_circle, originImage, params.spacing, orientation);

				//vote
				count_vote = 0;
				count_neighbor = 0;
				total_neighbor = 0;
				if (maskThreshold[xd] == background) {
					for (i_new = box.i_min; i_new <= box.i_max; i_new++) {
						for (j_new = box.j_min; j_new <= box.j_max; j_new++) {
							xd_new = x_new(i_new, j_new, length);
							
							current_point.x = (dataType)i_new;
							current_point.y = (dataType)j_new;
							getRealCoordFromImageCoord2D(current_point, originImage, params.spacing, orientation);
							dist_point = getPoint2DDistance(center_circle, current_point);

							//count foreground pixels in the band
							if (dist_point >= (radius - params.epsilon) && dist_point <= (radius + params.epsilon)) {
								total_neighbor++;
								if (maskThreshold[xd_new] == foreground) {
									count_vote++;
								}
							}
							//count foreground pixels close to the current center
							if (dist_point < (radius - params.epsilon) && maskThreshold[xd_new] == foreground) {
								count_neighbor++;
							}

						}
					}
				}

				//compute ratios
				hough_ratio = (dataType)count_vote / (dataType)total_neighbor;
				if (count_neighbor == 0) {
					houghSpacePtr[xd] = hough_ratio;
				}
			}
		}

		//keep the maximal ratio for each pixel
		for (i = box.i_min; i <= box.i_max; i++) {
			for (j = box.j_min; j <= box.j_max; j++) {
				xd = x_new(i, j, length);
				if (houghSpacePtr[xd] > houghSpaceMax[xd]) {
					houghSpaceMax[xd] = houghSpacePtr[xd];
				}
			}
		}

		//find the maximal ratio for current radius
		max_ratio = getTheMaxValue(houghSpacePtr, length, width);
		//find the maximal ratio for whole radius range
		if (max_ratio > found_ratio) {
			found_ratio = max_ratio;
			found_radius = radius;
			test_max = true;
			//copyDataToAnother2dArray(houghSpacePtr, houghSpaceMax, length, width);
		}
		//update the found_center if the value of found_ratio was updated
		if (test_max == true) {
			found_center = getPointWithMaximalValue2D(houghSpaceMax, length, width);
			test_max = false;
		}
		radius += params.radius_step;
	}

	printf("Found radius = %lf\n", found_radius);
	//printf("Found center : (%f, %f)\n", found_center.x, found_center.y);

	//add the seed
	foundCirclePtr[x_new((size_t)seed.x, (size_t)seed.y, length)] = 1.0;

	//draw maximal bounding box
	double radius_max_box = params.radius_max + params.epsilon;
	box = findBoundingBox2D(seed, length, width, radius_max_box, params.offset);
	for (i = box.i_min; i <= box.i_max; i++) {
		for (j = box.j_min; j <= box.j_max; j++) {
			xd = x_new(i, j, length);
			if (i == box.i_min || i == box.i_max || j == box.j_min || j == box.j_max) {
				foundCirclePtr[xd] = 1.0;
			}
		}
	}

	if (found_center.x != 0 && found_center.y != 0) {
		
		//add found center
		foundCirclePtr[x_new((size_t)found_center.x, (size_t)found_center.y, length)] = 1.0;
		
		//draw the found circle
		box = findBoundingBox2D(found_center, length, width, found_radius, params.offset);
		found_center = getRealCoordFromImageCoord2D(found_center, originImage, params.spacing, orientation);

		for (i = box.i_min; i <= box.i_max; i++) {
			for (j = box.j_min; j <= box.j_max; j++) {
				xd = x_new(i, j, length);
				
				current_point.x = (dataType)i;
				current_point.y = (dataType)j;
				current_point = getRealCoordFromImageCoord2D(current_point, originImage, params.spacing, orientation);

				dist_point = getPoint2DDistance(found_center, current_point);

				if (dist_point >= (found_radius - params.epsilon) && dist_point <= (found_radius + params.epsilon)) {
					foundCirclePtr[xd] = 1.0;
				}
			}
		}
	}

	//copy back
	copyDataToAnother2dArray(houghSpaceMax, houghSpacePtr, length, width);

	max_ratio = getTheMaxValue(houghSpaceMax, length, width);
	fprintf(saveInfo, "%f,%f\n", max_ratio, found_radius);

	free(maskThreshold);
	free(houghSpaceMax);
	free(gradientX);
	free(gradientY);
	free(normGradient);
	free(gradientAngle);
	free(statusPixel);

	return found_center;
}