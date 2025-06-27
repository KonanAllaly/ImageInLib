
#pragma once
#include "front_propagation.h"

	typedef enum
	{
		PARTIAL_FRONT_PROPAGATION = 1,
		FRONT_PROPAGATION = 2,
		DOUBLE_FRONT_PROPAGATION = 3,
		KEY_POINT_DETECTION = 4,
		DISTANCE_MAP = 5
	} PropagationType;

	void fastMarchingfrontPropagation2D(Image_Data2D inputImageData, dataType* potential, Point2D* endPoints, const PropagationType pType);

	void fastMarchingfrontPropagation(Image_Data inputImageData, dataType** potential, Point3D* endPoints, const PropagationType pType);