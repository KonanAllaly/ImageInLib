#pragma once

#include"common_functions.h"
#include"../src/thresholding.h"
#include"template_functions.h"
#include"../src/edgedetection.h"

typedef struct {
	double radius_min, radius_max; // define the range of radius
	double radius_step; //step for changing the radius 
	double epsilon, offset;
	dataType K; // coefficient of the edge detector
	dataType thres; //threshold of the edge detector
} HoughParameters;

void initialize2dArrayBool(bool* array2D, const size_t length, const size_t width);

/// <summary>
/// Initialize 2D array with 0.0
/// </summary>
/// <param name="array2D">2D array</param>
/// <param name="length">array length</param>
/// <param name="width">array width</param>
void initialize2dArray(dataType* array2D, const size_t length, const size_t width);

/// <summary>
/// Get the coordinates of the pixel with the maximal value
/// </summary>
/// <param name="arrayPtr2D">2D array</param>
/// <param name="length">array length</param>
/// <param name="width">array width</param>
/// <returns></returns>
Point2D getPointWithMaximalValue2D(dataType* arrayPtr2D, const size_t length, const size_t width);

/// <summary>
/// Detect circles in given image
/// </summary>
/// <param name="imageDataPtr">input image data</param>
/// <param name="houghSpacePtr">voting array for hough transform</param>
/// <param name="foundCirclePtr">found circles</param>
/// <param name="length">image length</param>
/// <param name="width">image width</param>
/// <param name="params">hough transform parameters</param>
//void circularHoughTransform(dataType* imageDataPtr, dataType* houghSpacePtr, dataType* foundCirclePtr, const size_t length, const size_t width, HoughParameters params);

/// <summary>
/// Detect cicles around given point
/// </summary>
/// <param name="seed">given point</param>
/// <param name="imageDataPtr">input image data</param>
/// <param name="houghSpacePtr">voting array for hough transform</param>
/// <param name="foundCirclePtr">found circles</param>
/// <param name="length">image length</param>
/// <param name="width">image width</param>
/// <param name="params">Hough transform parameters</param>
/// <param name="savingPath">saving path</param>
Point2D localHoughTransform(Image_Data2D imageDataStr, Point2D seed, dataType* houghSpacePtr, dataType* foundCirclePtr, HoughParameters params);

Point2D localHoughWithCanny(Image_Data2D imageDataStr, Point2D seed, dataType* houghSpacePtr, dataType* foundCirclePtr, HoughParameters params);

void houghTransform(Image_Data2D imageDataStr, dataType* foundCirclePtr, HoughParameters parameters);

bool circleDetection(Image_Data2D imageDataPtr, const HoughParameters hParameters);

Point2D localCircleDetection(Image_Data2D imageDataPtr, dataType* foundCirclePtr, Point2D seed, const HoughParameters hParameters);

Point2D get2dImagecentroid(dataType* imageDataPtr, size_t length, size_t width, dataType imageBackground);
