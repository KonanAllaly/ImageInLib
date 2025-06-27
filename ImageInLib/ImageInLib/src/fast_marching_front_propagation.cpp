#include "fast_marching_front_propagation.h"

void fastMarchingfrontPropagation2D(Image_Data2D inputImageData, dataType* actionMapPtr, dataType* potentialPtr, Point2D* endPoints, const PropagationType pType)
{
	switch (pType)
	{
    case PARTIAL_FRONT_PROPAGATION:
		partialFrontPropagation2D(inputImageData, actionMapPtr, potentialPtr, endPoints);
		break;
	/*
	case FRONT_PROPAGATION:
		frontPropagation(inputImageData, potential, endPoints);
		break;
	case DOUBLE_FRONT_PROPAGATION:
		doubleFrontPropagation(inputImageData, potential, endPoints);
		break;
	case KEY_POINT_DETECTION:
		frontPropagationsWithKeyPointDetection(inputImageData, potential, endPoints);
		break;
	case DISTANCE_MAP:
		distanceMapFastMarching(inputImageData, potential, endPoints);
		break*/;
	default:
		break;
	}
}

void fastMarchingfrontPropagation(Image_Data inputImageData, dataType** actionMapPtr, dataType** potentialPtr, Point3D* endPoints, const PropagationType pType)
{
	switch (pType)
	{
		case PARTIAL_FRONT_PROPAGATION:
			partialFrontPropagation(inputImageData, actionMapPtr, potentialPtr, endPoints);
			break;
			/*
		case FRONT_PROPAGATION:
			frontPropagation(inputImageData, potential, endPoints);
			break;
		case DOUBLE_FRONT_PROPAGATION:
			doubleFrontPropagation(inputImageData, potential, endPoints);
			break;
		case KEY_POINT_DETECTION:
			frontPropagationsWithKeyPointDetection(inputImageData, potential, endPoints);
			break;
		case DISTANCE_MAP:
			distanceMapFastMarching(inputImageData, potential, endPoints);
			break*/;
	default:
		break;
	}
}