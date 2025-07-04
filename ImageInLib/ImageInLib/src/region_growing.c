#include "region_growing.h"

bool regionGrowing2D(Image_Data2D inputImageData, dataType* segment, Point2D seed, dataType radius, dataType foreGroundValue)
{
	if (inputImageData.imageDataPtr == NULL || segment == NULL)
		return false;

	const size_t length = inputImageData.height;
	const size_t width = inputImageData.width;
	if (seed.x < 0 || seed.x > length || seed.y < 0 || seed.y > width)
		return false;

	//Find the threshold in region around the seed point
	BoundingBox2D box = findBoundingBox2D(seed, length, width, radius);

	//get the threshold value in region around seed point
	return true;
}

bool regionGrowing3D(Image_Data inputImageData, dataType** segment, Point3D seed, dataType radius, dataType foreGroundValue)
{
	return true;
}

void regionGrowingSegmentation(void* pInputImageData, void* segmentationResult, void* seed, const dataType radius, const dataType foreGroundValue, const dimRG dimension)
{
	switch (dimension)
	{
	case RG_2D:
	{
		Image_Data2D inputImageData = *(Image_Data2D*)pInputImageData;
		dataType* segment = (dataType*)segmentationResult;
		Point2D pStart = *(Point2D*)seed;
		regionGrowing2D(inputImageData, segment, pStart, radius, foreGroundValue);
		break;
	}
	case RG_3D:
	{
		Image_Data inputImageData = *(Image_Data*)pInputImageData;
		dataType** segment = (dataType**)segmentationResult;
		Point3D pStart = *(Point3D*)seed;
		regionGrowing3D(inputImageData, segment, pStart, radius, foreGroundValue);
		break;
	}
	default:
		break;
	}
}