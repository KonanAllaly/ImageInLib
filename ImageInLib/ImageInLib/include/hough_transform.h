#pragma once

#include"common_functions.h"
#include"../src/thresholding.h"
#include"template_functions.h"

typedef struct {
	const double radius_min, radius_max; // define the range of radius
	double epsilon, offset;
	dataType K; // coefficient of the edge detector
	dataType h; // space discretisation for gradien computation
	dataType thres; //threshold of the edge detector
	dataType radius_step; //step for changing the radius
} HoughParameters;

void initialize2dArray(dataType* array2D, const size_t length, const size_t width);

Point2D getPointWithMaximalValue2D(dataType* arrayPtr2D, const size_t length, const size_t width);

void circularHoughTransform(dataType* imageDataPtr, dataType* houghSpacePtr, dataType* foundCirclePtr, const size_t length, const size_t width, HoughParameters params);

void localHoughTransform(Point2D seed, dataType* imageDataPtr, dataType* houghSpacePtr, dataType* foundCirclePtr, const size_t length, const size_t width, HoughParameters params, std::string savingPath);
