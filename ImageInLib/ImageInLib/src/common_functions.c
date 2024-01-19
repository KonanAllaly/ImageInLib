#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h> 
#include <memory.h>
#include <math.h>
#include "common_functions.h"

#define M_PI 3.14159265358979323846
#define foreground 1.0
#define background 0.0

//==============================================================================
// Local Function Prototype
//==============================================================================
void reflection3DB(dataType** toReflectImage, size_t imageHeight, size_t imageLength, size_t imageWidth, size_t p)
{
	size_t k, i, j;
	size_t height = imageHeight, length = imageLength, width = imageWidth; // Actual Z, X, Y dimensions
	size_t rowLength = length + 2 * p;
	// Modified Dimension less 1
	size_t mheight = height - 1, mlength = length - 1, mwidth = width - 1;
	// Z Reflection
	for (i = p; i <= (mlength)+p; i++)
	{
		for (j = p; j <= (mwidth)+p; j++)
		{
			// If
			toReflectImage[p - 1][x_new(i, j, rowLength)] = toReflectImage[p + 1][x_new(i, j, rowLength)];
			toReflectImage[(mheight)+p + 1][x_new(i, j, rowLength)] = toReflectImage[(mheight)+p - 1][x_new(i, j, rowLength)];
		}
	}
	// Y reflection
	for (j = p; j <= (mwidth)+p; j++)
	{
		for (k = 0; k <= (mheight)+2 * p; k++)
		{
			toReflectImage[k][x_new(p - 1, j, rowLength)] = toReflectImage[k][x_new(p + 1, j, rowLength)];
			toReflectImage[k][x_new((mlength)+p + 1, j, rowLength)] = toReflectImage[k][x_new((mlength)+p - 1, j, rowLength)];
		}
	}
	// X Direction
	for (i = 0; i <= (mlength)+2 * p; i++)
	{
		for (k = 0; k <= (mheight)+2 * p; k++)
		{
			toReflectImage[k][x_new(i, p - 1, rowLength)] = toReflectImage[k][x_new(i, p + 1, rowLength)];
			toReflectImage[k][x_new(i, (mwidth)+p + 1, rowLength)] = toReflectImage[k][x_new(i, (mwidth)+p - 1, rowLength)];
		}
	}
}
//==============================================================================
//imageHeight and imageLength is considered as data extension together with boundary
void reflection3D(dataType** toReflectImage, size_t imageHeight, size_t imageLength, size_t imageWidth)
{
	size_t k, i, j;
	size_t length = imageLength, width = imageWidth; // Actual X, Y dimensions

	const size_t heightMin = imageHeight - 1;
	const size_t widthMin = width - 1;
	const size_t lengthMin = length - 1;
	const size_t sliceSize = imageLength * imageWidth;

	size_t x;

	// Y reflection
	for (k = 1; k <= heightMin; k++)
	{
		for (i = 0; i < length; i++)
		{
			toReflectImage[k][x_new(i, 0, length)] = toReflectImage[k][x_new(i, 1, length)];
			toReflectImage[k][x_new(i, widthMin, length)] = toReflectImage[k][x_new(i, widthMin - 1, length)];
		}
	}

	// X Direction
	for (k = 1; k <= heightMin; k++)
	{
		for (j = 0; j < width; j++)
		{
			x = x_new(0, j, length);
			toReflectImage[k][x] = toReflectImage[k][x + 1];

			x = x_new(lengthMin, j, length);
			toReflectImage[k][x] = toReflectImage[k][x - 1];
		}
	}

	// Z Reflection
	for (i = 0; i < sliceSize; i++)
	{
		toReflectImage[0][i] = toReflectImage[1][i];
		toReflectImage[heightMin][i] = toReflectImage[heightMin - 1][i];
	}
}
//==============================================================================
/*
* Evaluates the diffusivity value
*/
dataType gradientFunction(dataType value, dataType coef)
{
	return (dataType)(1.0 / (1 + coef * value));
}
//==============================================================================
size_t x_new(const size_t rowIndex, const size_t columnIndex, const size_t rowLength)
{
	return rowIndex + columnIndex * rowLength; // x + y*DimX
}
//==============================================================================
void copyDataToExtendedArea(dataType** originalDataPtr, dataType** extendedDataPtr, const size_t originalHeight, const size_t originalLength, const size_t originalWidth)
{
	const size_t height_ext = originalHeight + 2;
	const size_t length_ext = originalLength + 2;
	const size_t width_ext = originalWidth + 2;

	size_t sliceBound = (length_ext - 1) * width_ext;
	size_t i, k, k_ext, i_d = 0;

	for (k = 0, k_ext = 1; k < originalHeight; k++, k_ext++)
	{
		i_d = 0;
		for (i = length_ext + 1; i < sliceBound; i += length_ext)
		{
			memcpy(&(extendedDataPtr[k_ext][i]), &(originalDataPtr[k][i_d]), originalLength * sizeof(dataType));
			i_d += originalLength;
		}
	}
}
//==============================================================================
void copyDataToReducedArea(dataType** originalDataPtr, const dataType** extendedDataPtr, const size_t originalHeight, const size_t originalLength, const size_t originalWidth)
{
	const size_t height_ext = originalHeight + 2;
	const size_t length_ext = originalLength + 2;
	const size_t width_ext = originalWidth + 2;

	size_t sliceBound = (length_ext - 1) * width_ext;
	size_t i, k, k_ext, i_d = 0;

	for (k = 0, k_ext = 1; k < originalHeight; k++, k_ext++)
	{
		i_d = 0;
		for (i = length_ext + 1; i < sliceBound; i += length_ext)
		{
			memcpy(&(originalDataPtr[k][i_d]), &(extendedDataPtr[k_ext][i]), originalLength * sizeof(dataType));
			i_d += originalLength;
		}
	}
}
//==============================================================================
void copyDataToAnotherArray(dataType** source, dataType** destination, size_t height, size_t length, size_t width)
{
	size_t k, i, j, xd;
	for (k = 0; k < height; k++)
	{
		for (i = 0; i < length; i++)
		{
			for (j = 0; j < width; j++)
			{
				// 2D flattening
				xd = x_new(i, j, length);
				destination[k][xd] = source[k][xd];
			}
		}
	}
}
//==============================================================================
size_t x_flat(const size_t rowIndex, const size_t columnIndex, const size_t heightIndex, const size_t rowLength, const size_t columnLength)
{
	return rowIndex + rowLength * (columnIndex + columnLength * heightIndex); // x + xlen * (y + ylen * z)
}
//==============================================================================
void rescaleNewRange(dataType** imageDataPtr, size_t imageLength, size_t imageWidth, size_t imageHeight, dataType minNew, dataType maxNew, dataType max_dta, dataType min_dta) {
	size_t k, i, j, xd;

	// Find the Min and Max Intensity
	for (k = 0; k < imageHeight; k++) {
		for (i = 0; i < imageLength; i++) {
			for (j = 0; j < imageWidth; j++) {

				// 1D Conversion of row and column
				xd = x_new(i, j, imageLength);
				if (imageDataPtr[k][xd] > max_dta) {
					max_dta = imageDataPtr[k][xd];
				}
				if (imageDataPtr[k][xd] < min_dta) {
					min_dta = imageDataPtr[k][xd];
				}
			}
		}
	}
	// Rescale from min_new to max_new
	dataType diffOld = max_dta - min_dta;
	dataType diffNew = maxNew - minNew;
	dataType scale_factor = (diffNew) / (diffOld);
	for (k = 0; k < imageHeight; k++) {
		for (i = 0; i < imageLength; i++) {
			for (j = 0; j < imageWidth; j++) {
				// 1D Conversion of row and column
				xd = x_new(i, j, imageLength);
				imageDataPtr[k][xd] = scale_factor * (imageDataPtr[k][xd] - max_dta) + maxNew;
				// Alternatively
				//dta[k][xd] = scale_factor * (dta[k][xd] - min_dta) + min_new; }
			}
		}
	}
}
//==============================================================================
void copyDataToAnother2dArray(dataType* source, dataType* destination, size_t imageHeight, size_t imageWidth) {
	size_t i, dim2D = imageHeight * imageWidth;
	for (i = 0; i < dim2D; i++) {
		destination[i] = source[i];
	}
}
//==============================================================================
void copyDataTo2dExtendedArea(dataType* originalDataPtr, dataType* extendedDataPtr, const size_t originalHeight, const size_t originalWidth)
{
	const size_t height_ext = originalHeight + 2;
	const size_t width_ext = originalWidth + 2;

	size_t sliceBound = (height_ext - 1) * width_ext;
	size_t i, j, i_ext, j_ext, i_d = 0;

	for (i = 0, i_ext = 1; i < originalHeight; i++, i_ext++) {
		for (j = 0, j_ext = 1; j < originalWidth; j++, j_ext++) {
			extendedDataPtr[x_new(i_ext, j_ext, height_ext)] = originalDataPtr[x_new(i, j, originalHeight)];
		}
	}

}
//==============================================================================
void copyDataTo2dReducedArea(dataType* originalDataPtr, const dataType* extendedDataPtr, const size_t originalHeight, const size_t originalWidth)
{
	const size_t height_ext = originalHeight + 2;
	const size_t width_ext = originalWidth + 2;

	size_t sliceBound = (height_ext - 1) * width_ext;
	size_t i, j, i_ext, j_ext, i_d = 0;

	for (i = 0, i_ext = 1; i < originalHeight; i++, i_ext++) {
		for (j = 0, j_ext = 1; j < originalWidth; j++, j_ext++) {
			originalDataPtr[x_new(i, j, originalHeight)] = extendedDataPtr[x_new(i_ext, j_ext, height_ext)];
		}
	}

}
//==============================================================================
void reflection2D(dataType* toReflectImage, size_t imageHeight, size_t imageWidth)
{
	size_t i, j;
	size_t height = imageHeight, width = imageWidth;

	const size_t heightMin = height - 1;
	const size_t widthMin = width - 1;

	size_t x;

	// Y reflection
	for (i = 0; i < height; i++)
	{
		toReflectImage[x_new(i, 0, height)] = toReflectImage[x_new(i, 1, height)];
		toReflectImage[x_new(i, widthMin, height)] = toReflectImage[x_new(i, widthMin - 1, height)];
	}

	// X Direction
	for (j = 0; j < width; j++)
	{
		x = x_new(0, j, height);
		toReflectImage[x] = toReflectImage[x + 1];

		x = x_new(heightMin, j, height);
		toReflectImage[x] = toReflectImage[x - 1];
	}

}
//==============================================================================
double getPoint2DDistance(const Point2D a, const Point2D b)
{
	return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}
//==============================================================================
double getPoint3DDistance(Point3D a, Point3D b)
{
	return sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}
//==============================================================================
Point3D getPointWithTheHighestValue(dataType** distanceMapPtr, const size_t length, const size_t width, const size_t height) {
	
	Point3D result = { 0.0, 0.0, 0.0 };
	dataType max_value = 0.0;

	for (size_t k = 0; k < height; k++) {
		for (size_t i = 0; i < length; i++) {
			for (size_t j = 0; j < width; j++) {
				size_t x = x_new(i, j, length);
				if (distanceMapPtr[k][x] > max_value) {
					max_value = distanceMapPtr[k][x];
					result.x = i; result.y = j; result.z = k;
				}
			}
		}
	}
	return result;
}
//==============================================================================
/*
void circularHoughTransform(dataType* imageDataPtr, dataType* houghSpacePtr, dataType* votingArray, const size_t length, const size_t width, double radius) {
	
	size_t i, j, xd, dim2D = length * width;

	double epsilon = 0.0, offset = 5.0;
	size_t i_min, i_max, j_min, j_max;
	Point2D center_circle = { 0.0, 0.0 }, current_point = { 0.0, 0.0 };

	double ind_x = 0.0, ind_y = 0.0, dist_point = 0.0;
	size_t n, count_vote = 0, max_count = 0, N = 6;
	size_t i_new, j_new, xd_new;
	double perimeter = 0.0; //(dataType)(2 * 3.14 * radius);

	//initialization
	for (i = 0; i < dim2D; i++) {
		houghSpacePtr[i] = 0.0;
		votingArray[i] = 0.0;
	}

	//Find the number of foreground pixel lying on circle with given radius
	size_t count_neighbor = 0;
	double ratio = 0.0, neighbor_ratio = 0.0, surface = 0.0;

	for (n = 1; n < N; n++) {
		
		epsilon = 0.1 * n;
		perimeter = (double)(((2 * M_PI * (radius - epsilon)) + (2 * 3.14 * (radius + epsilon))) / 2);
		surface = (double)(M_PI * (radius - epsilon) * (radius - epsilon));
		
		//vote for all the pixels
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {
				xd = x_new(i, j, length);
				center_circle.x = i; 
				center_circle.y = j;
				if (imageDataPtr[xd] == 0.0) {

					//find the bounding box
					ind_x = (double)i - (radius + offset);
					if (ind_x <= 0) {
						i_min = 0;
					}
					else {
						i_min = (size_t)ind_x;
					}
					ind_y = (double)j - (radius + offset);
					if (ind_y <= 0) {
						j_min = 0;
					}
					else {
						j_min = (size_t)ind_y;
					}
					ind_x = (double)i + (radius + offset);
					if (ind_x >= (length - 1)) {
						i_max = length - 1;
					}
					else {
						i_max = (size_t)ind_x;
					}
					ind_y = (double)j + (radius + offset);
					if (ind_y >= (width - 1)) {
						j_max = width - 1;
					}
					else {
						j_max = (size_t)ind_y;
					}
					
					//vote
					count_vote = 0;
					count_neighbor = 0;
					for (i_new = i_min; i_new <= i_max; i_new++) {
						for (j_new = j_min; j_new <= j_max; j_new++) {
							
							xd_new = x_new(i_new, j_new, length);
							current_point.x = i_new;
							current_point.y = j_new;
							dist_point = getPoint2DDistance(center_circle, current_point);

							//count foreground pixels in the band
							if (imageDataPtr[xd_new] == 1.0 && dist_point >= (radius - epsilon) && dist_point <= (radius + epsilon)) {
								count_vote++;
							}

							//count foreground pixels close to the current center
							if (dist_point < (radius - epsilon) && imageDataPtr[xd_new] == 1.0) {
								count_neighbor++;
							}

						}
					}

					//printf("count : %d\n", count_vote);
					//printf("count neighboring foreground pixels : %d\n", count_neighbor);
					//printf("######################################################\n");

					ratio = (double)count_vote / perimeter;
					neighbor_ratio = (double)count_neighbor / surface;

					if (neighbor_ratio == 0.0) {
						votingArray[xd] += ratio;
					}
					else {
						votingArray[xd] += 0.0;
					}

				}
				else {
					votingArray[xd] += 0.0;
				}
			}
		}

	}

	//Normalization of the ratio
	for (i = 0; i < dim2D; i++) {
		votingArray[i] = votingArray[i] / (N - 1);
	}

	//find the max ratio
	double max_ratio = 0.0;
	for (i = 0; i < dim2D; i++) {
		if (votingArray[i] > max_ratio) {
			max_ratio = votingArray[i];
		}
	}
	
	//Draw the circle with maximum ratio
	for (i = 0; i < length; i++) {
		for (j = 0; j < width; j++) {
			
			xd = x_new(i, j, length);
			if (votingArray[xd] == max_ratio) {
				center_circle.x = i; 
				center_circle.y = j;
				ind_x = (double)i - (radius + offset);
				if (ind_x <= 0) {
					i_min = 0;
				}
				else {
					i_min = (size_t)ind_x;
				}
				ind_y = (double)j - (radius + offset);
				if (ind_y <= 0) {
					j_min = 0;
				}
				else {
					j_min = (size_t)ind_y;
				}

				ind_x = (double)i + (radius + offset);
				if (ind_x >= (length - 1)) {
					i_max = length - 1;
				}
				else {
					i_max = (size_t)ind_x;
				}
				ind_y = (double)j + (radius + offset);
				if (ind_y >= (width - 1)) {
					j_max = width - 1;
				}
				else {
					j_max = (size_t)ind_y;
				}
				
				for (i_new = i_min; i_new <= i_max; i_new++) {
					for (j_new = j_min; j_new <= j_max; j_new++) {
						xd_new = x_new(i_new, j_new, length);
						current_point.x = i_new; 
						current_point.y = j_new;
						dist_point = getPoint2DDistance(center_circle, current_point);
						if (dist_point <= radius) {
							houghSpacePtr[xd_new] = 1.0;
						}
					}
				}

			}
		}
	}
	
}
*/
//==============================================================================
void localCircularHoughTransform(Point2D seed, dataType* imageDataPtr, dataType* houghSpacePtr, dataType* votingArray, const size_t length, const size_t width, double radius, double offset) {

	size_t i, j, xd, dim2D = length * width;
	size_t i_new, j_new, xd_new;

	//initialization
	for (i = 0; i < dim2D; i++) {
		houghSpacePtr[i] = 0.0;
		votingArray[i] = 0.0;
	}

	double epsilon = 0.5, dist_point = 0.0;
	//double perimeter = 4 * M_PI * radius * epsilon; // surface of the band
	//double surface = pow(M_PI * (radius - epsilon), 2); // surface of the smaller circle

	Point2D center_circle = { 0.0, 0.0 }, current_point = { 0.0, 0.0 };
	BoundingBox2D box = { 0, 0, 0, 0 };

	size_t total_point = 0, count_vote = 0, count_neighbor = 0;
	dataType ratio = 0.0, neighbor_ratio = 0.0;
	
	//Find the bounding box for the seed point
	box = findBoundingBox2D(seed, length, width, radius, offset);

	//voting process : count number of forground pixels matching with criterium
	for (i = box.i_min; i <= box.i_max; i++) {
		for (j = box.j_min; j <= box.j_max; j++) {
			xd = x_new(i, j, length);
			center_circle.x = (dataType)i; 
			center_circle.y = (dataType)j;
			//vote for just background pixels
			if (imageDataPtr[xd] == background) {
				//vote
				count_vote = 0;
				count_neighbor = 0;
				for (i_new = box.i_min; i_new <= box.i_max; i_new++) {
					for (j_new = box.j_min; j_new <= box.j_max; j_new++) {
						xd_new = x_new(i_new, j_new, length);
						current_point.x = (dataType)i_new;
						current_point.y = (dataType)j_new;
						dist_point = getPoint2DDistance(center_circle, current_point);
						//count foreground pixels in the band
						if (dist_point >= (radius - epsilon) && dist_point <= (radius + epsilon)) {
							count_vote++;
							if (imageDataPtr[xd_new] == foreground) {
								total_point++;
							}
						}
						//count foreground pixels close to the current center
						if (dist_point < (radius - epsilon) && imageDataPtr[xd_new] == foreground) {
							count_neighbor++;
						}
					}
				}
				//compute ratios
				ratio = (dataType)count_vote / (dataType)total_point;
				neighbor_ratio = (dataType)count_neighbor;
				//votingArray[xd] = ratio;
				//filter out false centers
				if (neighbor_ratio == 0.0) {
					votingArray[xd] = ratio;
				}
			}
		}
	}

	//find the max ratio
	dataType max_ratio = 0.0;
	for (i = box.i_min; i <= box.i_max; i++) {
		for (j = box.j_min; j <= box.j_max; j++) {
			xd = x_new(i, j, length);
			if (votingArray[xd] > max_ratio) {
				max_ratio = votingArray[xd];
				center_circle.x = (dataType)i;
				center_circle.y = (dataType)j;
			}
		}
	}

	if (max_ratio > 0.0) {

		//==mark the center of circle
		houghSpacePtr[x_new((size_t)center_circle.x, (size_t)center_circle.y, length)] = 1.0;
		//=====================
		box = findBoundingBox2D(center_circle, length, width, radius, offset);

		//Draw the circle with maximum ratio
		for (i = box.i_min; i <= box.i_max; i++) {
			for (j = box.j_min; j <= box.j_max; j++) {
				xd = x_new(i, j, length);
				current_point.x = (dataType)i;
				current_point.y = (dataType)j;
				dist_point = getPoint2DDistance(center_circle, current_point);
				if (dist_point >= (radius - epsilon) && dist_point <= (radius + epsilon)) {
					houghSpacePtr[xd] = 1.0;
				}
			}
		}
	}

}
//==============================================================================
dataType getTheMaxValue(dataType* imageDataPtr, const size_t length, const size_t width) {
	dataType max_value = 0.0;
	for (size_t i = 0; i < length * width; i++) {
		if (imageDataPtr[i] > max_value) {
			max_value = imageDataPtr[i];
		}
	}
	return max_value;
}
//==============================================================================
void rescaleNewRange2D(dataType* imageDataPtr, size_t imageLength, size_t imageWidth, dataType minNew, dataType maxNew) {
	
	size_t i, dim2D = imageLength * imageWidth;
	dataType min_data = 100000.0, max_data = 0.0;
	// Find the Min and Max Intensity
	for (i = 0; i < dim2D; i++) {
		if (imageDataPtr[i] > max_data) {
			max_data = imageDataPtr[i];
		}
		if (imageDataPtr[i] < min_data) {
			min_data = imageDataPtr[i];
		}
	}
	// Rescale from min_new to max_new
	dataType diffOld = max_data - min_data;
	dataType diffNew = maxNew - minNew;
	dataType scale_factor = (diffNew) / (diffOld);
	for (i = 0; i < dim2D; i++) {
		imageDataPtr[i] = scale_factor * (imageDataPtr[i] - max_data) + maxNew;
	}
}
//==============================================================================
void houghTransform(dataType* imageDataPtr, dataType* houghSpace, const size_t length, const size_t width) {

	size_t i, j, n, xd, dim2D = length * width;
	dataType radius = 0.0, epsilon = 7.0, step = 0.5;
	const size_t N = 21;

	dataType*** accumulator = (dataType***)malloc(sizeof(dataType**) * length);
	for (i = 0; i < length; i++) {
		accumulator[i] = (dataType**)malloc(sizeof(dataType*) * width);
		for (j = 0; j < width; j++) {
			accumulator[i][j] = (dataType*)malloc(sizeof(dataType) * N);
		}
	}

	//initilize accumulator
	for (i = 0; i < length; i++) {
		for (j = 0; j < width; j++) {
			for (n = 0; n < N; n++) {
				accumulator[i][j][n] = 0;
			}
		}
	}

	int theta = 0;
	dataType x = 0.0, y = 0.0, xCord = 0.0, yCord = 0.0;
	size_t k, l;
	for (i = 0; i < length; i++) {
		for (j = 0; j < width; j++) {
			xCord = (dataType)i;
			yCord = (dataType)j;
			xd = x_new(i, j, length);
			for (n = 0; n < N; n++) {
				radius = epsilon + (dataType)n * step;
				//printf("radius = %f", radius);

				//for (k = 0; k < length; k++) {
				//	for (l = 0; l < width; l++) {
				//
				//	}
				//}

				for (theta = 0; theta <= 360; theta++) {
					x = (dataType)i - radius * cos(theta);
					y = (dataType)j - radius * sin(theta);
					
					//TO DO !!!
					if (imageDataPtr[xd] == 1.0 && x == xCord && y == yCord) {
						accumulator[i][j][n] += 1;
						printf("This condition happens");
					}
				}

			}
		}
	}

	//find the maximum of accumulator
	Point2D found_center = { 0.0, 0.0 };
	dataType found_radius = 0.0, max_accum = 0.0;
	for (i = 0; i < length; i++) {
		for (j = 0; j < width; j++) {
			for (n = 0; n < N; n++) {
				radius = epsilon + (dataType)n * step;
				if (accumulator[i][j][n] > max_accum) {
					max_accum = accumulator[i][j][n];
					found_radius = radius;
					found_center.x = i;
					found_center.y = j;
				}
			}
		}
	}

	printf("max accumulator = %f", max_accum);

	//draw the circle with found radius
	Point2D current = { 0.0, 0.0 };
	double dist_point = 0.0;
	for (i = 0; i < length; i++) {
		for (j = 0; j < width; j++) {
			xd = x_new(i, j, length);
			current.x = i; 
			current.y = j;
			dist_point = getPoint2DDistance(found_center, current);
			if (dist_point < found_radius) {
				houghSpace[xd] = 1.0;
			}
			else {
				houghSpace[xd] = 0.0;
			}
		}
	}

	for (i = 0; i < length; i++) {
		for (j = 0; j < width; j++) {
			free(accumulator[i][j]);
		}
		free(accumulator[i]);
	}
	free(accumulator);
}
//==============================================================================
BoundingBox2D findBoundingBox2D(Point2D point, const size_t length, const size_t width, double radius, double offset) {
	
	BoundingBox2D box = {0, 0, 0, 0};

	//find i_min
	double ind_x = (double)point.x - (radius + offset);
	if (ind_x < 0) {
		box.i_min = 0;
	}
	else {
		if (ind_x < length - 1) {
			box.i_min = (size_t)ind_x;
		}
		else {
			box.i_min = length - 1;
		}
	}

	//find i_max
	ind_x = (double)point.x + (radius + offset);
	if (ind_x >= length - 1) {
		box.i_max = length - 1;
	}
	else {
		box.i_max = (size_t)ind_x;
	}

	//find j_min
	double ind_y = (double)point.y - (radius + offset);
	if (ind_y < 0) {
		box.j_min = 0;
	}
	else {
		if (ind_y < width - 1) {
			box.j_min = (size_t)ind_y;
		}
		else {
			box.j_min = width - 1;
		}
	}

	//find j_max
	ind_y = (double)point.y + (radius + offset);
	if (ind_y >= width - 1) {
		box.j_max = width - 1;
	}
	else {
		box.j_max = (size_t)ind_y;
	}

	return box;
}
//==============================================================================
void computeImageGradient(dataType* imageDataPtr, dataType* gradientVectorX, dataType* gradientVectorY, const size_t length, const size_t width, dataType h) {
	size_t i, j, xd;
	dataType ux = 0.0, uy = 0.0;
	for (i = 0; i < length; i++) {
		for (j = 0; j < width; j++) {
			xd = x_new(i, j, length);
			//x direction
			if (i == 0) {
				ux = (imageDataPtr[x_new(i + 1, j, length)] - imageDataPtr[x_new(i, j, length)]) / h;
			}
			else {
				if (i == length - 1) {
					ux = (imageDataPtr[x_new(i, j, length)] - imageDataPtr[x_new(i - 1, j, length)]) / h;
				}
				else {
					ux = (imageDataPtr[x_new(i + 1, j, length)] - imageDataPtr[x_new(i - 1, j, length)]) / (2 * h);
				}
			}
			//x direction
			if (j == 0) {
				uy = (imageDataPtr[x_new(i, j + 1, length)] - imageDataPtr[x_new(i, j, length)]) / h;
			}
			else {
				if (j == width - 1) {
					uy = (imageDataPtr[x_new(i, j, length)] - imageDataPtr[x_new(i, j - 1, length)]) / h;
				}
				else {
					uy = (imageDataPtr[x_new(i, j + 1, length)] - imageDataPtr[x_new(i, j - 1, length)]) / (2 * h);
				}
			}

			gradientVectorX[xd] = ux;
			gradientVectorY[xd] = uy;
		}
	}
}