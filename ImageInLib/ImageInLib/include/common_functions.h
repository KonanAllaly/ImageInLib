
#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef COMMON_FUNCTIONS
#define COMMON_FUNCTIONS
	//==============================================================================
	//debug constants
#ifndef MEASURE_TIME
#define MEASURE_TIME
#endif
// Show/Hide Output in Console
#ifndef CONSOLE_OUTPUT
#define CONSOLE_OUTPUT
#endif
#undef CONSOLE_OUTPUT
//==============================================================================
// Typedefs
// Data Type

#ifndef USE_DOUBLE_PRECISION
	typedef float dataType;
#else
	typedef double dataType;
#endif // USE_DOUBLE_PRECISION

	//==============================================================================
	// Includes
#include <stddef.h>
#include <omp.h>
#include <stdbool.h>
//==============================================================================
// MACROs
//==============================================================================
// STRUCTs
// Common 2D Points - {x,y}
	typedef struct {
		dataType x, y;
	} Point2D;

	// Common 3D Points - {x,y,z}
	typedef struct {
		dataType x, y, z;
	} Point3D;

	//Structure to handle image spacing
	typedef struct {
		dataType sx, sy, sz;
	} VoxelSpacing;

	//Matrix for rotation
	typedef struct {
		Point3D v1, v2, v3;
	}OrientationMatrix;

	// Image Container and Properties
	typedef struct {
		// Image Dimensions
		size_t height, length, width; // Absolute Dimension
		dataType** imageDataPtr; // Image Data Containers
		Point3D origin; // image origin
		VoxelSpacing spacing; // pixel size and distance between slice
		OrientationMatrix orientation;
	} Image_Data;

	// Generate Random Points
	typedef struct {
		size_t k, xd, p;
	}Random3dPoints;

	typedef struct {
		size_t xd, p;
	}Random2dPoints;

	typedef struct {
		Point2D v1, v2;
	}OrientationMatrix2D;

	typedef struct {
		dataType sx, sy;
	} PixelSpacing;

	typedef struct {
		// Image Dimensions
		size_t height, width; // Absolute Dimension
		dataType* imageDataPtr; // Image Data Containers
		Point2D origin; // image origin
		PixelSpacing spacing; // pixel size
		OrientationMatrix2D orientation;
	} Image_Data2D;

	typedef struct {
		dataType min_data, max_data, mean_data, sd_data, variance;
	} Statistics;

	//==============================================================================
	// Shapes Container
	typedef struct {
		// Shapes
		dataType** shp;
		size_t num; // shape number
	} Shapes;
	//==============================================================================
	// Enumeration
	enum FiniteDifference
	{
		FINITE_FORWARD = 1, FINITE_BACKWARD, FINITE_CENTRAL
	};
	//==============================================================================
	// FUNCTION PROTOTYPES
	/*
	* Calculate to 1D representation from 2D Array
	* x from column and row indices
	*/
	size_t x_new(const size_t rowIndex, const size_t columnIndex, const size_t rowLength);
	//==============================================================================
	/*
	Flattens 3D array to 1D array
	*/
	size_t x_flat(const size_t rowIndex, const size_t columnIndex, const size_t heightIndex, const size_t rowLength, const size_t columnLength);
	//==============================================================================
	/*
	* Reflection Function Zuska
	* toReflectImage is the data to perform reflection on
	* imageHeight is the actual Z dimension
	* imageLength is the actual X dimension
	* imageWidth is the actual Y dimension
	* P is the boundary constant/distance
	* Reflection does not include the current point but only its adjacent neighbor
	* Preserves Neumann Boundary condition - Target is for those methods
	*/
	void reflection3DB(dataType** toReflectImage, size_t imageHeight, size_t imageLength, size_t imageWidth, size_t p);
	//==============================================================================
	/*
	* Reflection Function Jozef
	* Reflection includes the current point together with its adjacent neighbor
	* Targets Distance Function types
	*/
	void reflection3D(dataType** toReflectImage, size_t imageHeight, size_t imageLength, size_t imageWidth);
	//==============================================================================
	/*
	* Gradient Calculation function
	*/
	dataType gradientFunction(dataType value, dataType coef);
	//==============================================================================
	/*
	void copyDataToExtendedArea(const dataType ** originalDataPtr, dataType ** extendedDataPtr, const size_t originalHeight, const size_t originalLength, const size_t originalWidth);

	The function copies data from source given by originalDataPtr
	representing data of dimensions given by originalHeight, originalLength and originalWidth
	to an array extendedDataPtr representing enlarged data (by 1 vx in each direction)
	of dimensions originalHeight + 2, originalLength + 2 and originalWidth + 2
	*/
	void copyDataToExtendedArea(dataType** originalDataPtr, dataType** extendedDataPtr, const size_t originalHeight, const size_t originalLength, const size_t originalWidth);
	//==============================================================================
	void copyDataToReducedArea(dataType** originalDataPtr, dataType** extendedDataPtr, const size_t originalHeight, const size_t originalLength, const size_t originalWidth);
	//==============================================================================
	// Copy from one pointer to another pointer
	void copyDataToAnotherArray(dataType** source, dataType** destination, size_t height, size_t length, size_t width);
	//==============================================================================
	void rescaleNewRange(dataType** imageDataPtr, size_t imageLength, size_t imageWidth, size_t imageHeight, dataType minNew, dataType maxNew, dataType max_dta, dataType min_dta);
	//==============================================================================
	typedef struct {
		size_t k_min, i_min, j_min, k_max, i_max, j_max;
	} ClipBox;
	//==============================================================================
	// Calc. centroid of image data
	void centroidImage(dataType** imageDataPtr, dataType* centroid, size_t imageHeight, size_t imageLength, size_t imageWidth, dataType imageBackground);
	void centroidClipBox(dataType* centroid, ClipBox coord, dataType** imageDataPtr, size_t imageLength, dataType imageBackground);
	//==============================================================================
	void copyDataToAnother2dArray(dataType* source, dataType* destination, size_t imageHeight, size_t imageWidth);
	//==============================================================================
	void copyDataTo2dExtendedArea(dataType* originalDataPtr, dataType* extendedDataPtr, const size_t originalHeight, const size_t originalWidth);
	//==============================================================================
	void copyDataTo2dReducedArea(dataType* originalDataPtr, const dataType* extendedDataPtr, const size_t originalHeight, const size_t originalWidth);
	//==============================================================================
	void reflection2D(dataType* toReflectImage, size_t imageHeight, size_t imageWidth);
	//==============================================================================
	double getPoint2DDistance(const Point2D a, const Point2D b);
	//==============================================================================
	double getPoint3DDistance(Point3D a, Point3D b);
	//==============================================================================
	Point3D getPointWithTheHighestValue(dataType** distanceMapPtr, const size_t length, const size_t width, const size_t height);
	//==============================================================================
	dataType getTheMaxValue(dataType* imageDataPtr, const size_t length, const size_t width);
	//==============================================================================
	void rescaleNewRange2D(dataType* imageDataPtr, size_t imageLength, size_t imageWidth, dataType minNew, dataType maxNew);
	//==============================================================================
	void computeImageGradient(dataType* imageDataPtr, dataType* gradientVectorX, dataType* gradientVectorY, const size_t length, const size_t width, dataType h);
	//==============================================================================
	typedef struct {
		size_t i_min, i_max, j_min, j_max;
	} BoundingBox2D;

	typedef struct {
		size_t i_min, i_max, j_min, j_max, k_min, k_max;
	} BoundingBox3D;

	/// <summary>
	/// Find the bounding box for given 2D point
	/// </summary>
	/// <param name="point">point of interest</param>
	/// <param name="length">length of the image</param>
	/// <param name="width">width of the image</param>
	/// <param name="radius">radius to be considered for the bounding box</param>
	/// <param name="offset">offset</param>
	/// <returns></returns>
	BoundingBox2D findBoundingBox2D(Point2D point, const size_t length, const size_t width, double radius, double offset);

	/// <summary>
	/// Find the bounding box of 3D given point
	/// </summary>
	/// <param name="point">point of interest</param>
	/// <param name="length">image length</param>
	/// <param name="width">image width</param>
	/// <param name="height">image height</param>
	/// <param name="radius">radius of the box</param>
	/// <param name="offset">extend</param>
	/// <returns></returns>
	BoundingBox3D findBoundingBox3D(Point3D point, const size_t length, const size_t width, const size_t height, double radius, double offset);

	//================================================================================
	//void computeHistogram(dataType** imageDataPtr, const size_t length, const size_t width, const size_t height, const size_t bins);
	//================================================================================
	//void maximumIntensityProjection(dataType** imageData, dataType** resultImage, const size_t length, const size_t width, const size_t height);
	//===============

	bool getGradient2D(dataType* imageDataPtr, const size_t width, const size_t length, const size_t ind_x, const size_t ind_y, const PixelSpacing fVolume, Point2D* grad);

	bool getGradient3D(dataType** imageDataPtr, const size_t width, const size_t length, const size_t height, const size_t ind_x, const size_t ind_y, const size_t ind_z, const VoxelSpacing fVolume, Point3D* grad);

#endif // !COMMON_FUNCTIONS

#ifdef __cplusplus
}
#endif