#include <stdlib.h>
#include "labeling.h"

LinkedList createLinkedList()
{
	LinkedList linkedList;
	linkedList.number_of_points = 0;
	return linkedList;
}

void popFirstElement(labelingList* linked_list)
{
	if (linked_list == NULL || linked_list->first_point == NULL)
	{
		return;
	}
	labelingPoint* current_point = linked_list->first_point;
	labelingPoint* next_point = current_point->next;
	if (next_point != NULL)
	{
		next_point->previous = NULL;
	}
	free(current_point);
	linked_list->first_point = next_point;
	linked_list->number_of_points--;
}

labelingPoint* pushElementToList(labelingList* linked_list, labelingPoint* linked_point, const dataType x, const dataType y, const dataType z)
{
	if (linked_list == NULL || linked_point == NULL)
	{
		return NULL; // Safety check
	}

	labelingPoint* new_linked_point = (labelingPoint*)malloc(sizeof(labelingPoint));
	if (new_linked_point == NULL)
	{
		return NULL; // malloc failed
	}
	new_linked_point->x = x;
	new_linked_point->y = y;
	new_linked_point->z = z;
	new_linked_point->next = linked_point->next;
	new_linked_point->previous = linked_point;
	new_linked_point->id = getNextID();
	if (linked_point->next != NULL)
	{
		linked_point->next->previous = new_linked_point;
	}
	linked_point->next = new_linked_point;
	linked_list->number_of_points++;
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

	//default value for z coordinate
	dataType z = -1;

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

					label++;
					LinkedList neighbours_list = createLinkedList();
					labelingPoint* current_point = (labelingPoint*)malloc(sizeof(labelingPoint));
					current_point->x = i;
					current_point->y = j;
					current_point->z = z;//In 2D : set the third coordinate to -1
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
								current_point = pushElementToList(&neighbours_list, current_point, i - 1, j, z);
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
								current_point = pushElementToList(&neighbours_list, current_point, i + 1, j, z);
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
								current_point = pushElementToList(&neighbours_list, current_point, i, j - 1, z);
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
								current_point = pushElementToList(&neighbours_list, current_point, i, j + 1, z);
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

						labelingPoint* pcurrent_point = neighbours_list.first_point;
						size_t it = (size_t)pcurrent_point->x;
						size_t jt = (size_t)pcurrent_point->y;
						size_t xdt = x_new(it, jt, length);
						popFirstElement(&neighbours_list);
						
						//North
						if(it > 0)
						{
							size_t xn = x_new(it - 1, jt, length);
							if(statusArray[xn] == 3)
							{
								if(inputImageData.imageDataPtr[xn] == foreGroundValue)
								{
									current_point = pushElementToList(&neighbours_list, current_point, it - 1, jt, z);
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
									current_point = pushElementToList(&neighbours_list, current_point, it + 1, jt, z);
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
									current_point = pushElementToList(&neighbours_list, current_point, it, jt - 1, z);
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
									current_point = pushElementToList(&neighbours_list, current_point, it, jt + 1, z);
									statusArray[xe] = 2;
								}
								else
								{
									statusArray[xe] = 1;
								}
							}
						}
	
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

	const size_t height = inputImageData.height;
	const size_t length = inputImageData.length;
	const size_t width = inputImageData.width;
	const size_t dim2D = length * width;

	short** statusArray = (short**)malloc(sizeof(short*) * height);
	for(size_t k = 0; k < height; k++)
	{
		statusArray[k] = (short*)malloc(sizeof(short) * dim2D);
	}
	if (statusArray == NULL) {
		return false;
	}
	// Status meaning:
	// 3: unvisited
	// 2: discovered (in list)
	// 1: processed (done)

	//Initialize
	for(size_t k = 0; k < height; k++)
	{
		for (size_t i = 0; i < dim2D; i++) 
		{
			statusArray[k][i] = 3;
		}
	}
	
	int label = 0;
	for(size_t k = 0; k < height; k++)
	{
		for (size_t i = 0; i < length; i++)
		{
			for (size_t j = 0; j < width; j++)
			{
				size_t xd = x_new(i, j, length);
				if (statusArray[k][xd] == 3)
				{
					if (inputImageData.imageDataPtr[k][xd] == foreGroundValue)
					{
						label++;
						LinkedList neighbours_list = createLinkedList();
						labelingPoint* current_point = (labelingPoint*)malloc(sizeof(labelingPoint));
						current_point->x = i;
						current_point->y = j;
						current_point->z = k;
						current_point->previous = NULL;
						current_point->next = NULL;
						neighbours_list.first_point = current_point;

						//Top
						if (k > 0) {
							if (statusArray[k - 1][xd] == 3)
							{
								if (inputImageData.imageDataPtr[k - 1][xd] == foreGroundValue)
								{
									current_point = pushElementToList(&neighbours_list, current_point, i, j, k - 1);
									statusArray[k - 1][xd] = 2;
								}
								else {
									statusArray[k - 1][xd] = 1;
								}
							}
						}
						//Bottom
						if (k < (height - 1)) {
							if (statusArray[k + 1][xd] == 3)
							{
								if (inputImageData.imageDataPtr[k + 1][xd] == foreGroundValue)
								{
									current_point = pushElementToList(&neighbours_list, current_point, i, j, k + 1);
									statusArray[k + 1][xd] = 2;
								}
								else {
									statusArray[k + 1][xd] = 1;
								}
							}
						}
						//North
						if (i > 0) {
							size_t xn = x_new(i - 1, j, length);
							if (statusArray[k][xn] == 3)
							{
								if (inputImageData.imageDataPtr[k][xn] == foreGroundValue)
								{
									current_point = pushElementToList(&neighbours_list, current_point, i - 1, j, k);
									statusArray[k][xn] = 2;
								}
								else {
									statusArray[k][xn] = 1;
								}
							}
						}
						//South
						if (i < (length - 1)) {
							size_t xs = x_new(i + 1, j, length);
							if (statusArray[xs] == 3)
							{
								if (inputImageData.imageDataPtr[k][xs] == foreGroundValue)
								{
									current_point = pushElementToList(&neighbours_list, current_point, i + 1, j, k);
									statusArray[k][xs] = 2;
								}
								else {
									statusArray[k][xs] = 1;
								}
							}
						}
						//West
						if (j > 0) {
							size_t xw = x_new(i, j - 1, length);
							if (statusArray[xw] == 3)
							{
								if (inputImageData.imageDataPtr[k][xw] == foreGroundValue)
								{
									current_point = pushElementToList(&neighbours_list, current_point, i, j - 1, k);
									statusArray[k][xw] = 2;
								}
								else {
									statusArray[k][xw] = 1;
								}
							}
						}
						//East
						if (j < (width - 1)) {
							size_t xe = x_new(i, j + 1, length);
							if (statusArray[k][xe] == 3)
							{
								if (inputImageData.imageDataPtr[k][xe] == foreGroundValue)
								{
									current_point = pushElementToList(&neighbours_list, current_point, i, j + 1, k);
									statusArray[k][xe] = 2;
								}
								else {
									statusArray[k][xe] = 1;
								}
							}
						}

						//fix the current point
						segment[k][xd] = label; //set the label
						statusArray[k][xd] = 1; //set the status
						while (neighbours_list.number_of_points > 0) {

							labelingPoint* pcurrent_point = neighbours_list.first_point;
							size_t it = (size_t)pcurrent_point->x;
							size_t jt = (size_t)pcurrent_point->y;
							size_t kt = (size_t)pcurrent_point->z;
							size_t xdt = x_new(it, jt, length);
							popFirstElement(&neighbours_list);

							//Top
							if (kt > 0)
							{
								if (statusArray[kt - 1][xdt] == 3)
								{
									if (inputImageData.imageDataPtr[kt - 1][xdt] == foreGroundValue)
									{
										current_point = pushElementToList(&neighbours_list, current_point, it, jt, kt - 1);
										statusArray[kt - 1][xdt] = 2;
									}
									else
									{
										statusArray[kt - 1][xdt] = 1;
									}
								}
							}
							//Bottom
							if (kt < (height - 1))
							{
								if (statusArray[kt + 1][xdt] == 3)
								{
									if (inputImageData.imageDataPtr[kt + 1][xdt] == foreGroundValue)
									{
										current_point = pushElementToList(&neighbours_list, current_point, it, jt, kt + 1);
										statusArray[kt + 1][xdt] = 2;
									}
									else
									{
										statusArray[kt + 1][xdt] = 1;
									}
								}
							}
							//North
							if (it > 0)
							{
								size_t xn = x_new(it - 1, jt, length);
								if (statusArray[kt][xn] == 3)
								{
									if (inputImageData.imageDataPtr[kt][xn] == foreGroundValue)
									{
										current_point = pushElementToList(&neighbours_list, current_point, it - 1, jt, kt);
										statusArray[kt][xn] = 2;
									}
									else
									{
										statusArray[kt][xn] = 1;
									}
								}
							}
							//South
							if (it < (length - 1))
							{
								size_t xs = x_new(it + 1, jt, length);
								if (statusArray[kt][xs] == 3)
								{
									if (inputImageData.imageDataPtr[kt][xs] == foreGroundValue)
									{
										current_point = pushElementToList(&neighbours_list, current_point, it + 1, jt, kt);
										statusArray[kt][xs] = 2;
									}
									else
									{
										statusArray[kt][xs] = 1;
									}
								}
							}
							//West
							if (jt > 0)
							{
								size_t xw = x_new(it, jt - 1, length);
								if (statusArray[kt][xw] == 3)
								{
									if (inputImageData.imageDataPtr[kt][xw] == foreGroundValue)
									{
										current_point = pushElementToList(&neighbours_list, current_point, it, jt - 1, kt);
										statusArray[kt][xw] = 2;
									}
									else
									{
										statusArray[kt][xw] = 1;
									}
								}
							}
							//East
							if (jt < (width - 1))
							{
								size_t xe = x_new(it, jt + 1, length);
								if (statusArray[kt][xe] == 3)
								{
									if (inputImageData.imageDataPtr[kt][xe] == foreGroundValue)
									{
										current_point = pushElementToList(&neighbours_list, current_point, it, jt + 1, kt);
										statusArray[kt][xe] = 2;
									}
									else
									{
										statusArray[kt][xe] = 1;
									}
								}
							}

							segment[kt][xdt] = label; //set the label
							statusArray[kt][xdt] = 1; //set the status
						}
					}
					else {
						statusArray[k][xd] == 1;
					}
				}
			}
		}
	}

	for (size_t k = 0; k < height; k++) {
		free(statusArray[k]);
	}
	free(statusArray);

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