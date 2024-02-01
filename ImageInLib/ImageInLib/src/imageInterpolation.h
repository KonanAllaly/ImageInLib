#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include <stdio.h>
#include <stdbool.h>
#include "common_functions.h"
#include "transformation.h"
#include "interpolations.h"

	typedef enum
	{
		NEAREST_NEIGHBOR = 1,
		BILINEAR, // for 2D image
		TRILINEAR // for 3D image
	} interpolationMethod;

	typedef struct
	{
		dataType** minimum;
		dataType** maximum;
		dataType** mean;
		dataType** sd;
	} statictics_Pointers;

	//====================
	//3D Functions

	/*
	* Point3D getImageCoordFromRealCoord3D(Point3D src_point, Point3D real_origin, VoxelSpacing real_spacing, OrientationMatrix orientation)
	* src_point : source 3D points in image coordinates system
	* real_origin : origin in real world
	* real_spacing : - distance between voxels in x and y direction
	*                 - distance between slices
	* orientation : matrix for orientation/rotation IJK((1, 0, 0),(0,1,0),(0,0,1)), RAS((-1,0,0),(0,-1,0),(0,0,1))
	* or LPS((0,0,1),(0,-1,0),(0,0,-1))
	* For given point in real world cordinates system the function provide
	* the corresponding point in image coordinates system
	*/
	Point3D getRealCoordFromImageCoord3D(Point3D image_point, Point3D real_origin, VoxelSpacing real_spacing, OrientationMatrix orientation);

	/*
	* Point3D getRealCoordFomImageCoord3D(Point3D src_point, Point3D real_origin, VoxelSpacing real_spacing, OrientationMatrix orientation)
	* src_point : source 3D points in image coordinates system
	* real_origin : origin in real world
	* real_spacing : - distance between voxels in x and y direction
	*                 - distance between slices
	* orientation : matrix for orientation/rotation IJK((1, 0, 0),(0,1,0),(0,0,1)), RAS((-1,0,0),(0,-1,0),(0,0,1))
	* or LPS((0,0,1),(0,-1,0),(0,0,-1))
	* For given point in image cordinates system the function provide
	* the corresponding point in real world coordinates system
	*/
	Point3D getImageCoordFromRealCoord3D(Point3D real_point, Point3D real_origin, VoxelSpacing real_spacing, OrientationMatrix orientation);

	/*
	* This function return the interpolated value by nearest neighbor approach
	* dataType getInterpolatedValueNearestNeighbor3D(Image_Data src_image, Point3D point)
	* Image_Data src_image : source image
	* Point3D point : the current point
	*/
	dataType getInterpolatedValueNearestNeighbor3D(Image_Data src_image, Point3D point);

	/*
	* This function return the interpolated value by bilinear interpolation approach
	* dataType getInterpolatedValueTrilinear3D(Image_Data src_image, Point3D point)
	* Image_Data src_image : source image
	* Point3D point : the current point
	*/
	dataType getInterpolatedValueTrilinear3D(Image_Data src_image, Point3D point);

	/*
	* This function perform image interpolation
	* bool imageInterpolation3D(Image_Data src_image, Image_Data dest_image, dataType scale_factor, interpolationMethod method)
	* Image_Data src_image : source image structure
	*                          - height, lenght, width : image dimension
	*                          - imageDataPtr : pointer for pixels value
	*                          - origin : image origin
	*                          - spacing : pixel size
	* Image_Data dest_image : interpolated image structure
	*                          - height, Lenght, width : interpolated image dimension
	*                          - imageDataPtr : pointer for interpolated pixels value
	*                          - origin : interpolated image origin
	*                          - spacing : interpolated pixel size
	* interpolationMethod method : nearest neighbor, bilinear (for 2D image), trilinear 
	*/
	bool imageInterpolation3D(Image_Data src_image, Image_Data dest_image, interpolationMethod method);

	//=====================================================
	//2D Functions

	/*
	* Point2D getImageCoordFromRealCoord2D(Point2D src_point, Point2D real_origin, PixelSpacing real_spacing, OrientationMatrix2D orientation)
	* Point2D src_point : source 2D points in image coordinates system
	* Point2D real_origin : origin in real world coordinates system
	* PixelSpacing real_spacing : distance between voxels in x and y direction
	* orientation : matrix for orientation
	* For given point in image cordinates system, the function provides
	* the corresponding point in image real world coordinates system
	*/
	Point2D getRealCoordFromImageCoord2D(Point2D src_point, Point2D real_origin, PixelSpacing real_spacing, OrientationMatrix2D orientation);

	/*
	* Point2D getImageCoordFromRealCoord2D(Point2D src_point, Point2D real_origin, PixelSpacing real_spacing, OrientationMatrix2D orientation)
	* Point2D src_point : source 2D points in image coordinates system
	* Point2D real_origin : origin in real world coordinates system
	* PixelSpacing real_spacing : distance between voxels in x and y direction
	* OrientationMatrix2D orientation : matrix for orientation
	* For given point in real world cordinates system, the function provides
	* the corresponding point in image coordinates system
	*/
	Point2D getImageCoordFromRealCoord2D(Point2D src_point, Point2D image_origin, PixelSpacing real_spacing, OrientationMatrix2D orientation);

	/*
	* This function return the interpolated value by nearest neighbor approach
	* dataType getInterpolatedValueNearestNeighbor2D(Image_Data2D src_image, Point2D point)
	* Image_Data2D src_image : source image
	* Point2D point : the current point
	*/
	dataType getInterpolatedValueNearestNeighbor2D(Image_Data2D src_image, Point2D point);

	/*
	* This function return the interpolated value by bilinear interpolation approach
	* dataType getInterpolatedValueNearestNeighbor2D(Image_Data2D src_image, Point2D point)
	* Image_Data2D src_image : source image
	* Point2D point : the current point
	*/
	dataType getInterpolatedValueBilinear2D(Image_Data2D src_image, Point2D point);

	/*
	* This function perform image interpolation
	* bool imageInterpolation2D(Image_Data2D src_image, Image_Data2D dest_image, dataType scale_factor, interpolationMethod method)
	* Image_Data2D src_image : source image structure
	*                          - height, width : image dimension
	*                          - imageDataPtr : pointer for pixels value
	*                          - origin : image origin
	*                          - spacing : pixel size
	* Image_Data2D dest_image : interpolated image structure
	*                          - height, width : interpolated image dimension
	*                          - imageDataPtr : pointer for interpolated pixels value
	*                          - origin : interpolated image origin
	*                          - spacing : interpolated pixel size
	* interpolationMethod method : nearest neighbor, bilinear, trilinear(for 3D images) 
	*/
	bool imageInterpolation2D(Image_Data2D src_image, Image_Data2D dest_image, interpolationMethod method);

	/*
	* This function return statistics (min, max, mean, standart deviation) in image
	* around given point in real world coordinates system
	* getStatistics(Image_Data imageData, Point3D point_ct, dataType radius);
	* imageData : structure containing the information of the image (lenght, width, heigth, origin, voxels value...)
	* point : the given point in real world coordinates system
	* radius : radium of the region (ball) where to compute the statistics
	*/
	Statistics getStatistics(Image_Data imageData, Point3D point_ct, dataType radius);

	Statistics getStats(Image_Data imageData, Point3D point_of_interest, double radius);

	Statistics getPointNeighborhoodStats(Image_Data imageData, Point3D point_of_interest, double radius);

	/// <summary>
	/// : This function create statistics images from statistics computed in each voxel
	/// in a small ball defined by the radius
	/// </summary>
	/// <param name="imageData"> : Structure to hold the input image</param>
	/// <param name="statsImage"> : Structure to hold the statistics images</param>
	/// <param name="radius"> : radius to be considered for the statistics</param>
	/// <returns></returns>
	bool generateStatisticsImages(Image_Data imageData, statictics_Pointers statsImage, double radius);

#ifdef __cplusplus
}

#endif