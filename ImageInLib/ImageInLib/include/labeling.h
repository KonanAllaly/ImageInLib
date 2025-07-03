#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef LABELING
#define LABELING

#include <stdbool.h>
#include "common_functions.h"

	typedef enum
	{
		LABEL_2D_ARRAY = 1,
		LABEL_3D_ARRAY,
	} labelingApproach;

	typedef struct labelingPoint
	{
		dataType x;
		dataType y;
		dataType z;
		struct labelingPoint* next;
		struct labelingPoint* previous;
		int id;
	} labelingPoint;

	typedef struct labelingList
	{
		size_t number_of_points;
		labelingPoint* first_point;
	} labelingList;

	typedef struct LinkedList
	{
		size_t number_of_points;
		labelingPoint* first_point;
	} LinkedList;

	LinkedList createLinkedList();

	void popFirstElement(labelingList* linked_list);

	bool labeling2D(Image_Data2D inputImageData, dataType* segment, dataType foreGroundValue);
	
	bool labeling3D(Image_Data inputImageData, dataType** segment, dataType foreGroundValue);
	
	void connectedComponentsLabeling(void* pInputImageData, void* segmentationResult, dataType foreGroundValue, const labelingApproach dimension);

#endif // !LABELING

#ifdef __cplusplus
}
#endif
