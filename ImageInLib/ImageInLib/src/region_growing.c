#include "region_growing.h"

bool regionGrowing2D(Image_Data2D inputImageData, dataType* segment, dataType foreGroundValue)
{
	return true;
}

bool regionGrowing3D(Image_Data inputImageData, dataType** segment, dataType foreGroundValue)
{
	return true;
}

void regionGrowingSegmentation(void* pInputImageData, void* segmentationResult, void* seed, dataType foreGroundValue, const dimRG dimension)
{
	switch (dimension)
	{
	case RG_2D:
	{
		Image_Data2D inputImageData = *(Image_Data2D*)pInputImageData;
		dataType* segment = (dataType*)segmentationResult;
		Point2D* seed = (Point2D*)seed;
		regionGrowing2D(inputImageData, segment, seed, foreGroundValue);
		break;
	}
	case RG_3D:
	{
		Image_Data inputImageData = *(Image_Data*)pInputImageData;
		dataType** segment = (dataType**)segmentationResult;
		Point3D* seed = (Point3D*)seed;
		regionGrowing3D(inputImageData, segment, foreGroundValue);
		break;
	}
	default:
		break;
	}
}