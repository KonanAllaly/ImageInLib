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
void copyDataToReducedArea(dataType** originalDataPtr, dataType** extendedDataPtr, const size_t originalHeight, const size_t originalLength, const size_t originalWidth)
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
					result.x = i; 
					result.y = j; 
					result.z = k;
				}
			}
		}
	}
	return result;
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
	//dataType min_data = 100000.0, max_data = 0.0;
	dataType min_data = imageDataPtr[0], max_data = imageDataPtr[0];
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
BoundingBox3D findBoundingBox3D(Point3D point, const size_t length, const size_t width, const size_t height, double radius, double offset) {

	if (point.x < 0 || point.y < 0 || point.z < 0) {
		printf("The coodinates of current point should be in image coordinates\n");
		return (BoundingBox3D) { 0, 0, 0, 0, 0, 0 };
	}

	size_t x_min = 0, x_max = 0, y_min = 0, y_max = 0, z_min = 0, z_max = 0;

	//find i_min
	double ind_min_x = (double)point.x - (radius + 0.5);
	if (ind_min_x < 0) {
		x_min = 0;
	}
	else {
		x_min = (size_t)ind_min_x;
	}

	//find i_max
	double ind_max_x = (double)point.x + (radius + 0.5);
	if (ind_max_x >= length - 1) {
		x_max = length - 1;
	}
	else {
		x_max = (size_t)ind_max_x;
	}

	//find j_min
	double ind_min_y = (double)point.y - (radius + 0.5);
	if (ind_min_y < 0) {
		y_min = 0;
	}
	else {
		y_min = (size_t)ind_min_y;
	}

	//find j_max
	double ind_max_y = (double)point.y + (radius + 0.5);
	if (ind_max_y >= width - 1) {
		y_max = width - 1;
	}
	else {
		y_max = (size_t)ind_max_y;
	}

	//find k_min
	double ind_min_z = (double)point.z - (radius + 0.5);
	if (ind_min_z < 0) {
		z_min = 0;
	}
	else {
		z_min = (size_t)ind_min_z;
	}

	//find k_max
	double ind_max_z = (double)point.z + (radius + 0.5);
	if (ind_max_z >= height - 1) {
		z_max = height - 1;
	}
	else {
		z_max = (size_t)ind_max_z;
	}

	return (BoundingBox3D) {x_min, x_max, y_min, y_max, z_min, z_max};
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
//==============================================================================
bool getGradient2D(dataType* imageDataPtr, const size_t width, const size_t length, const size_t ind_x, const size_t ind_y, const PixelSpacing fVolume, Point2D* grad)
{
	
	if (imageDataPtr == NULL || width < 2 || length < 2 || fVolume.sx == 0 || fVolume.sy == 0)
	{
		return false;
	}
	
	if (imageDataPtr == NULL || width < 2 || length < 2)
	{
		return false;
	}

	size_t x = ind_x, y = ind_y;
	if (x >= width)
	{
		x = width - 1;
	}

	if (y >= length)
	{
		y = length - 1;
	}

	dataType dx = 0.0, dy = 0.0;
	dataType hx_c = 2 * fVolume.sx;
	dataType hy_c = 2 * fVolume.sy;

	if (x == 0)
	{
		dx = (imageDataPtr[x_new(x + 1, y, width)] - imageDataPtr[x_new(x, y, width)]) / fVolume.sx;
	}
	else if (x == width - 1)
	{
		dx = (imageDataPtr[x_new(x, y, width)] - imageDataPtr[x_new(x - 1, y, width)]) / fVolume.sx;
	}
	else
	{
		dx = (imageDataPtr[x_new(x + 1, y, width)] - imageDataPtr[x_new(x - 1, y, width)]) / hx_c;
	}

	if (y == 0)
	{
		dy = (imageDataPtr[x_new(x, y + 1, width)] - imageDataPtr[x_new(x, y, width)]) / fVolume.sy;
	}
	else if (y == length - 1)
	{
		dy = (imageDataPtr[x_new(x, y, width)] - imageDataPtr[x_new(x, y - 1, width)]) / fVolume.sy;
	}
	else
	{
		dy = (imageDataPtr[x_new(x, y + 1, width)] - imageDataPtr[x_new(x, y - 1, width)]) / hy_c;
	}

	grad->x = dx;
	grad->y = dy;
	
	return true;
}
//==============================================================================
bool getGradient3D(dataType** imageDataPtr, const size_t width, const size_t length, const size_t height, const size_t ind_x, const size_t ind_y, const size_t ind_z, const VoxelSpacing fVolume, Point3D* grad)
{

	if (imageDataPtr == NULL || width < 2 || length < 2 || height < 2 ||
		fVolume.sx == 0 || fVolume.sy == 0 || fVolume.sz == 0)
	{
		return false;
	}


	if (imageDataPtr == NULL || width < 2 || length < 2 || height < 2)
	{
		return false;
	}

	size_t x = ind_x, y = ind_y, z = ind_z;
	if (x >= width)
	{
		x = width - 1;
	}

	if (y >= length)
	{
		y = length - 1;
	}

	if (z >= height)
	{
		z = height - 1;
	}

	dataType dx = 0.0, dy = 0.0, dz = 0.0;
	dataType hx_c = 2 * fVolume.sx;
	dataType hy_c = 2 * fVolume.sy;
	dataType hz_c = 2 * fVolume.sz;

	if (x == 0)
	{
		dx = (imageDataPtr[z][x_new(x + 1, y, width)] - imageDataPtr[z][x_new(x, y, width)]) / fVolume.sx;
	}
	else if (x == width - 1)
	{
		dx = (imageDataPtr[z][x_new(x, y, width)] - imageDataPtr[z][x_new(x - 1, y, width)]) / fVolume.sx;
	}
	else
	{
		dx = (imageDataPtr[z][x_new(x + 1, y, width)] - imageDataPtr[z][x_new(x - 1, y, width)]) / hx_c;
	}

	if (y == 0)
	{
		dy = (imageDataPtr[z][x_new(x, y + 1, width)] - imageDataPtr[z][x_new(x, y, width)]) / fVolume.sy;
	}
	else if (y == length - 1)
	{
		dy = (imageDataPtr[z][x_new(x, y, width)] - imageDataPtr[z][x_new(x, y - 1, width)]) / fVolume.sy;
	}
	else
	{
		dy = (imageDataPtr[z][x_new(x, y + 1, width)] - imageDataPtr[z][x_new(x, y - 1, width)]) / hy_c;
	}

	size_t xd = x_new(x, y, width);
	if (z == 0)
	{
		dz = (imageDataPtr[z + 1][xd] - imageDataPtr[z][xd]) / fVolume.sz;
	}
	else if (z == height - 1)
	{
		dz = (imageDataPtr[z][xd] - imageDataPtr[z - 1][xd]) / fVolume.sz;
	}
	else
	{
		dz = (imageDataPtr[z + 1][xd] - imageDataPtr[z - 1][xd]) / hz_c;
	}

	grad->x = dx;
	grad->y = dy;
	grad->z = dz;

	return true;
}


