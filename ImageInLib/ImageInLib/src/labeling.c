#include <stdlib.h>
#include "labeling.h"

void popFirstEltList(LinkedCurve* linked_curve)
{
	if (linked_curve == NULL)
	{
		return;
	}
	LinkedPoint* current_point = linked_curve->first_point;
	LinkedPoint* next_point = current_point->next;
	free(current_point);
	current_point = NULL;
	linked_curve->number_of_points--;
}

LinkedPoint* pushPointToList(LinkedCurve* linked_curve, LinkedPoint* linked_point, const size_t x, const size_t y)
{
	LinkedPoint* new_linked_point = (LinkedPoint*)malloc(sizeof(LinkedPoint));
	new_linked_point->x = x;
	new_linked_point->y = y;
	new_linked_point->next = NULL;
	new_linked_point->previous = NULL;
	new_linked_point->distance_to_next = 0;
	new_linked_point->id = 0;//getNextID();

	if (linked_point->next != NULL)
	{
		new_linked_point->next = linked_point->next;
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

	bool* statusArray = (bool*)malloc(sizeof(bool) * dim2D);
	if (statusArray == NULL) {
		return false;
	}

	//Initialize
	for (size_t i = 0; i < dim2D; i++) {
		statusArray[i] = false;
	}

	int label = 0;
	for (size_t i = 0; i < length; i++) 
	{
		for (size_t j = 0; j < width; j++) 
		{
			size_t xd = x_new(i, j, length);
			if(statusArray[xd] == false)
			{
				if(inputImageData.imageDataPtr[xd] == foreGroundValue)
				{	
					//new segment
					label++;
					LinkedCurve neighbours_list = createLinkedCurve();
					LinkedPoint* current_point = (LinkedPoint*)malloc(sizeof(LinkedPoint));
					current_point->x = i;
					current_point->y = j;
					current_point->next = NULL;
					neighbours_list.first_point = current_point;

					//North
					if (i > 0 && i < length && j >= 0 && j < width) {
						if (statusArray[x_new(i - 1, j, length)] == false)
						{
							if (inputImageData.imageDataPtr[x_new(i - 1, j, length)] == foreGroundValue)
							{
								current_point = pushPointToList(&neighbours_list, current_point, i - 1, j);
							}
							else {
								statusArray[x_new(i - 1, j, length)] == true;
							}
						}
					}
					//South
					if (i >= 0 && i < (length - 1) && j >= 0 && j < width) {
						if (statusArray[x_new(i + 1, j, length)] == false)
						{
							if (inputImageData.imageDataPtr[x_new(i + 1, j, length)] == foreGroundValue)
							{
								current_point = pushPointToList(&neighbours_list, current_point, i + 1, j);
							}
							else {
								statusArray[x_new(i + 1, j, length)] == true;
							}
						}
					}
					//West
					if (i >= 0 && i < length && j > 0 && j < width) {
						if (statusArray[x_new(i, j - 1, length)] == false)
						{
							if (inputImageData.imageDataPtr[x_new(i, j - 1, length)] == foreGroundValue)
							{
								current_point = pushPointToList(&neighbours_list, current_point, i, j - 1);
							}
							else {
								statusArray[x_new(i, j - 1, length)] == true;
							}
						}
					}
					//East
					if (i >= 0 && i < length && j >= 0 && j < (width - 1)) {
						if (statusArray[x_new(i, j + 1, length)] == false)
						{
							if (inputImageData.imageDataPtr[x_new(i, j + 1, length)] == foreGroundValue)
							{
								current_point = pushPointToList(&neighbours_list, current_point, i, j + 1);
							}
							else {
								statusArray[x_new(i, j + 1, length)] == true;
							}
						}
					}

					//fix the current point
					segment[xd] = label; //set the label
					statusArray[xd] = true; //set the status

					LinkedPoint* pcurrent_point = neighbours_list.first_point;
					LinkedPoint* next_point = NULL;
					size_t it, jt, xdt;
					while (neighbours_list.number_of_points > 0) {
						next_point = pcurrent_point->next;
						it = (size_t)pcurrent_point->x;
						jt = (size_t)pcurrent_point->y;
						xdt = x_new(it, jt, length);
						popFirstEltList(&neighbours_list);
						//North
						if(it > 0 && it < length && jt >= 0 && jt < width)
						{
							if(statusArray[x_new(it - 1, jt, length)] == false)
							{
								if(inputImageData.imageDataPtr[x_new(it - 1, jt, length)] == foreGroundValue)
								{
									current_point = pushPointToList(&neighbours_list, current_point, it - 1, jt);
								}else
								{
									statusArray[x_new(it - 1, jt, length)] = true;
								}
							}
						}
						//South
						if (it >= 0 && it < (length - 1) && jt >= 0 && jt < width)
						{
							if (statusArray[x_new(it + 1, jt, length)] == false)
							{
								if (inputImageData.imageDataPtr[x_new(it + 1, jt, length)] == foreGroundValue)
								{
									current_point = pushPointToList(&neighbours_list, current_point, it + 1, jt);
								}
								else
								{
									statusArray[x_new(it + 1, jt, length)] = true;
								}
							}
						}
						//West
						if (it >= 0 && it < length && jt > 0 && jt < width)
						{
							if (statusArray[x_new(it, jt - 1, length)] == false)
							{
								if (inputImageData.imageDataPtr[x_new(it, jt - 1, length)] == foreGroundValue)
								{
									current_point = pushPointToList(&neighbours_list, current_point, it, jt - 1);
								}
								else
								{
									statusArray[x_new(it, jt - 1, length)] = true;
								}
							}
						}
						//East
						if (it >= 0 && it < length && jt >= 0 && jt < (width - 1))
						{
							if (statusArray[x_new(it, jt + 1, length)] == false)
							{
								if (inputImageData.imageDataPtr[x_new(it, jt + 1, length)] == foreGroundValue)
								{
									current_point = pushPointToList(&neighbours_list, current_point, it, jt + 1);
								}
								else
								{
									statusArray[x_new(it, jt + 1, length)] = true;
								}
							}
						}
						if (next_point != NULL) {
							pcurrent_point = next_point;
						}
						//fix the current point
						segment[xdt] = label; //set the label
						statusArray[xdt] = true; //set the status
					}
				}
				else{
					statusArray[xd] == true;
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
