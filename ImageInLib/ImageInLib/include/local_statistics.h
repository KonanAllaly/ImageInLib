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

	Statistics statistics2D(Image_Data2D inputImageData, Point2D pointOfInterest, const dataType radiusOfSearchDomain, const seachDomainShape shape);

	Statistics statistics3D(Image_Data inputImageData, Point3D pointOfInterest, const dataType radiusOfSearchDomain, const seachDomainShape shape);

	void getStatisticsAroundPoint(void* pInputImageData, void* pointOfInterest, Statistics stats, const dataType radiusOfSearchDomain, const seachDomainShape shape, const pDimension dim);

#endif // !LOCAL_STATISTICS

#ifdef __cplusplus
}
#endif
