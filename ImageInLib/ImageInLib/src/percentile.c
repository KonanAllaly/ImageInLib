#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "percentile.h"

void bubbleSort(dataType* pArray, const size_t height, const size_t width)
{
	if(pArray == NULL)
		return;
	size_t i, j;
	size_t dim2D = height * width;
	for (i = dim2D - 1; i > 0; i--) 
	{
		for(j = 0; j < i; j++) 
		{
			if (pArray[j + 1] < pArray[j]) 
			{
				swap_elts(&pArray[j], &pArray[j + 1], sizeof(dataType));
			}
		}
	}
}

int partition(dataType* pArray, size_t low, size_t high)
{
	////choose the pivot
	//dataType pivot = pArray[(size_t)((high + low) / 2)];
	dataType pivot = pArray[high];

	int i = (low - 1); //index of smaller element
	for (int j = low; j < high; j++)
	{
		//if current element is smaller than or equal to pivot
		if (pArray[j] < pivot)
		{
			i++;
			swap_elts(&pArray[i], &pArray[j], sizeof(dataType));
		}
	}

	//Move pivot after the last smaller element
	swap_elts(&pArray[i + 1], &pArray[high], sizeof(dataType));
	//return the index of the pivot
	return (i + 1);
}

void quickSort(dataType* pArray, size_t low, size_t high)
{
	if(pArray == NULL)
		return;

	if (low < high) 
	{
		//pi is partitioning index, pArray[pi] is now at right place
		int pi = partition(pArray, low, high);
		//Recursively sort elements before partition and after partition
		quickSort(pArray, low, pi - 1);
		quickSort(pArray, pi + 1, high);
	}
}

dataType computePercentile(dataType* pArray, const size_t height, const size_t width, const double perc)
{
	if (pArray == NULL || perc < 0.0 || perc > 100.0)
		return;
	size_t dim2D = height * width;
	size_t index = (size_t)round((perc / 100.0) * (dim2D - 1));
	bubbleSort(pArray, height, width);
	return pArray[index];
}