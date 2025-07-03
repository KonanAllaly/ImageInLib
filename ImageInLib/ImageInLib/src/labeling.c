#include <stdlib.h>
#include "labeling.h"

void popFirstEltList(LinkedCurve* linked_curve)
{
	if (linked_curve == NULL || linked_curve->first_point == NULL)
	{
		return;
	}
	LinkedPoint* current_point = linked_curve->first_point;
	LinkedPoint* next_point = current_point->next;
	if (next_point != NULL)
	{
		next_point->previous = NULL;
	}
	free(current_point);
	linked_curve->first_point = next_point;
	linked_curve->number_of_points--;
}

LinkedPoint* pushPointToList(LinkedCurve* linked_curve, LinkedPoint* linked_point, const double x, const double y)
{
	if (linked_curve == NULL || linked_point == NULL)
	{
		return NULL; // Safety check
	}

	LinkedPoint* new_linked_point = (LinkedPoint*)malloc(sizeof(LinkedPoint));
	if (new_linked_point == NULL)
	{
		return NULL; // malloc failed
	}
	new_linked_point->x = x;
	new_linked_point->y = y;
	new_linked_point->next = linked_point->next;
	new_linked_point->previous = linked_point;
	new_linked_point->distance_to_next = 0;
	new_linked_point->id = getNextID();

	if (linked_point->next != NULL)
	{
		linked_point->next->previous = new_linked_point;
	}

	linked_point->next = new_linked_point;
	linked_curve->number_of_points++;

	return new_linked_point;
}

bool labeling2D(Image_Data2D inputImageData, dataType* segment, dataType foreGroundValue)
{
	if (inputImageData.imageDataPtr == NULL || segment == NULL)
		return false;

	const size_t length = inputImageData.height;
	const size_t width = inputImageData.width;
	const size_t dim2D = length * width;

	short* statusArray = (short*)malloc(sizeof(short) * dim2D);
	if (statusArray == NULL) {
		return false;
	}
	// Status meaning:
	// 3: unvisited
	// 2: discovered (in list)
	// 1: processed (done)

	//Initialize
	for (size_t i = 0; i < dim2D; i++) {
		statusArray[i] = 3;
	}

	int label = 0;
	for (size_t i = 0; i < length; i++) 
	{
		for (size_t j = 0; j < width; j++) 
		{
			size_t xd = x_new(i, j, length);
			if(statusArray[xd] == 3)
			{
				if(inputImageData.imageDataPtr[xd] == foreGroundValue)
				{	
					//new segment
					label++;
					LinkedCurve neighbours_list = createLinkedCurve();
					LinkedPoint* current_point = (LinkedPoint*)malloc(sizeof(LinkedPoint));
					current_point->x = i;
					current_point->y = j;
					current_point->previous = NULL;
					current_point->next = NULL;
					neighbours_list.first_point = current_point;

					//North
					if (i > 0) {
						size_t xn = x_new(i - 1, j, length);
						if (statusArray[xn] == 3)
						{
							if (inputImageData.imageDataPtr[xn] == foreGroundValue)
							{
								current_point = pushPointToList(&neighbours_list, current_point, i - 1, j);
								statusArray[xn] = 2;
							}
							else {
								statusArray[xn] = 1;
							}
						}
					}
					//South
					if (i < (length - 1)) {
						size_t xs = x_new(i + 1, j, length);
						if (statusArray[xs] == 3)
						{
							if (inputImageData.imageDataPtr[xs] == foreGroundValue)
							{
								current_point = pushPointToList(&neighbours_list, current_point, i + 1, j);
								statusArray[xs] = 2;
							}
							else {
								statusArray[xs] = 1;
							}
						}
					}
					//West
					if (j > 0) {
						size_t xw = x_new(i, j - 1, length);
						if (statusArray[xw] == 3)
						{
							if (inputImageData.imageDataPtr[xw] == foreGroundValue)
							{
								current_point = pushPointToList(&neighbours_list, current_point, i, j - 1);
								statusArray[xw] = 2;
							}
							else {
								statusArray[xw] = 1;
							}
						}
					}
					//East
					if (j < (width - 1)) {
						size_t xe = x_new(i, j + 1, length);
						if (statusArray[xe] == 3)
						{
							if (inputImageData.imageDataPtr[xe] == foreGroundValue)
							{
								current_point = pushPointToList(&neighbours_list, current_point, i, j + 1);
								statusArray[xe] = 2;
							}
							else {
								statusArray[xe] = 1;
							}
						}
					}

					//fix the current point
					segment[xd] = label; //set the label
					statusArray[xd] = 1; //set the status

					while (neighbours_list.number_of_points > 0) {
						
						LinkedPoint* pcurrent_point = neighbours_list.first_point;
						size_t it = (size_t)pcurrent_point->x;
						size_t jt = (size_t)pcurrent_point->y;
						size_t xdt = x_new(it, jt, length);
						popFirstEltList(&neighbours_list);
						
						//North
						if(it > 0)
						{
							size_t xn = x_new(it - 1, jt, length);
							if(statusArray[xn] == 3)
							{
								if(inputImageData.imageDataPtr[xn] == foreGroundValue)
								{
									current_point = pushPointToList(&neighbours_list, current_point, it - 1, jt);
									statusArray[xn] = 2;
								}else
								{
									statusArray[xn] = 1;
								}
							}
						}
						//South
						if (it < (length - 1))
						{
							size_t xs = x_new(it + 1, jt, length);
							if (statusArray[xs] == 3)
							{
								if (inputImageData.imageDataPtr[xs] == foreGroundValue)
								{
									current_point = pushPointToList(&neighbours_list, current_point, it + 1, jt);
									statusArray[xs] = 2;
								}
								else
								{
									statusArray[xs] = 1;
								}
							}
						}
						//West
						if (jt > 0)
						{
							size_t xw = x_new(it, jt - 1, length);
							if (statusArray[xw] == 3)
							{
								if (inputImageData.imageDataPtr[xw] == foreGroundValue)
								{
									current_point = pushPointToList(&neighbours_list, current_point, it, jt - 1);
									statusArray[xw] = 2;
								}
								else
								{
									statusArray[xw] = 1;
								}
							}
						}
						//East
						if (jt < (width - 1))
						{
							size_t xe = x_new(it, jt + 1, length);
							if (statusArray[xe] == 3)
							{
								if (inputImageData.imageDataPtr[xe] == foreGroundValue)
								{
									current_point = pushPointToList(&neighbours_list, current_point, it, jt + 1);
									statusArray[xe] = 2;
								}
								else
								{
									statusArray[xe] = 1;
								}
							}
						}
	
						//fix the current point
						segment[xdt] = label; //set the label
						statusArray[xdt] = 1; //set the status
					}
				}
				else{
					statusArray[xd] == 1;
				}
			}
		}
	}

	free(statusArray);
	return true;
}

bool labeling3D(Image_Data inputImageData, dataType** segment, dataType foreGroundValue)
{
	if (inputImageData.imageDataPtr == NULL || segment == NULL)
		return false;

	const size_t length = inputImageData.length;
	const size_t width = inputImageData.width;
	const size_t height = inputImageData.height;
	return true;
}

void connectedComponentsLabeling(void* pInputImageData, void* segmentationResult, dataType foreGroundValue, const labelingApproach dimension) {
	switch(dimension)
	{
	case LABEL_2D_ARRAY:
	{
		Image_Data2D inputImageData = *(Image_Data2D*)pInputImageData;
		dataType* segment = (dataType*)segmentationResult;
		labeling2D(inputImageData, segment, foreGroundValue);
		break;
	}
	case LABEL_3D_ARRAY:
	{
		Image_Data inputImageData = *(Image_Data*)pInputImageData;
		dataType** segment = (dataType**)segmentationResult;
		labeling3D(inputImageData, segment, foreGroundValue);
		break;
	}
	default:
		break;
	}
}