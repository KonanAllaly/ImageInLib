#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef LOCAL_STATISTICS
#define LOCAL_STATISTICS

#include "labeling.h"

	typedef enum
	{
		SQUARE = 1,
		SPHERE,
	} seachDomainShape;

	bool statistics2D(Image_Data2D inputImageData, Point2D pointOfInterest, const dataType radiusOfSearchDomain);

	bool statistics3D(Image_Data inputImageData, Point3D pointOfInterest, dataType radiusOfSearchDomain);

	void getStatisticsAroundPoint(void* pInputImageData, void* pointOfInterest, const pDimension dim, const seachDomainShape shape);

#endif // !LOCAL_STATISTICS

#ifdef __cplusplus
}
#endif
