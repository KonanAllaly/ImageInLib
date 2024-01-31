/*
* Author: Markjoe Olunna UBA
* Purpose: ImageInLife project - 4D Image Segmentation Methods
* Language:  C
*/
#include <stdio.h>
#include <stdlib.h>
#include "edgedetection.h"
#include "data_initialization.h"
#include "thresholding.h"
#include "common_functions.h"
#include "math.h"

#define M_PI 3.14159265358979323846

void edgeDetection3dFunctionD(dataType** image3DPtr, dataType** edge3DPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType bgroundvalue, dataType fgroundvalue, dataType insideValue)
{
	size_t i, k;//loop counter for z dimension
	const size_t dim2D = xDim * yDim;
	//==============================================================================
	//thresholding of voxel values
	thresholding3dFunctionD(image3DPtr, xDim, yDim, zDim, (bgroundvalue + fgroundvalue) / 2, bgroundvalue, fgroundvalue);
	//==============================================================================
	// initialize the edge array to background value
	initialize3dArrayD(edge3DPtr, xDim, yDim, zDim, bgroundvalue);
	//==============================================================================
	// filling the 2d array with number=value
	for (k = 0; k < zDim; k++)
	{
		for (i = 0; i < dim2D; i++) {
			if (image3DPtr[k][i] == fgroundvalue)//while on the object
			{
				if (k == 0)
				{
					if (i == 0)//if at the top left corner,
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
						{
							edge3DPtr[k][i] = fgroundvalue; //check if pixel neighbors has background value
						}
						else if ((image3DPtr[k][i + 1] == fgroundvalue) || (image3DPtr[k][i + xDim] == fgroundvalue) || (image3DPtr[k + 1][i] == fgroundvalue))
						{
							edge3DPtr[k][i] = insideValue; //check if pixel neighbors has foreground value
						}
					}
					else if (i == (xDim - 1)) //if at the top right coner,
					{
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
						{
							edge3DPtr[k][i] = fgroundvalue; //check if pixel neighbors has background value
						}
						else if ((image3DPtr[k][i - 1] == fgroundvalue) || (image3DPtr[k][i + xDim] == fgroundvalue) || (image3DPtr[k + 1][i] == fgroundvalue))
						{
							edge3DPtr[k][i] = insideValue; //check if pixel neighbors has foreground value
						}
					}
					else if (i == ((yDim - 1) * xDim)) //if at the bottom left corner,
					{
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
						{
							edge3DPtr[k][i] = fgroundvalue; //check if pixel neighbors has background value
						}
						else if ((image3DPtr[k][i + 1] == fgroundvalue) || (image3DPtr[k][i - xDim] == fgroundvalue) || (image3DPtr[k + 1][i] == fgroundvalue))
						{
							edge3DPtr[k][i] = insideValue; //check if pixel neighbors has foreground value
						}
					}
					else if (i == (dim2D - 1)) //if at the bottom right corner,
					{
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
						{
							edge3DPtr[k][i] = fgroundvalue; // check if pixel neighbors has background value
						}
						else if ((image3DPtr[k][i - 1] == fgroundvalue) || (image3DPtr[k][i - xDim] == fgroundvalue) || (image3DPtr[k + 1][i] == fgroundvalue))
						{
							edge3DPtr[k][i] = insideValue; //check if pixel neighbors has foreground value
						}
					}
					else if ((i > 0) && (i < (xDim - 1)))//if at first row (excluding left and right corner),
					{
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
						{
							edge3DPtr[k][i] = fgroundvalue; // //check if pixel neighbors has background value
						}
						else if ((image3DPtr[k][i - 1] == fgroundvalue) || (image3DPtr[k][i + 1] == fgroundvalue) || (image3DPtr[k][i + xDim] == fgroundvalue) || (image3DPtr[k + 1][i] == fgroundvalue))
						{
							edge3DPtr[k][i] = insideValue; //check if pixel neighbors has foreground value
						}
					}
					else if ((i > ((yDim - 1) * xDim)) && (i < (dim2D - 1))) //if at last row
					{ //(excluding left and right corner)
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
						{
							edge3DPtr[k][i] = fgroundvalue; // check if pixel neighbors has background value
						}
						else if ((image3DPtr[k][i - 1] == fgroundvalue) || (image3DPtr[k][i + 1] == fgroundvalue) || (image3DPtr[k][i - xDim] == fgroundvalue) || (image3DPtr[k + 1][i] == fgroundvalue))
						{
							edge3DPtr[k][i] = insideValue; //check if pixel neighbors has foreground value
						}
					}
					else if ((i != 0) && (i != ((yDim - 1) * xDim)) && ((i % xDim) == 0))//if at first column
					{ //(excluding top and bottom left corner)
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) || (image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
						{
							edge3DPtr[k][i] = fgroundvalue; // check if pixel neighbors has background value
						}
						else if ((image3DPtr[k][i + 1] == fgroundvalue) || (image3DPtr[k][i + xDim] == fgroundvalue) || (image3DPtr[k][i - xDim] == fgroundvalue) || (image3DPtr[k + 1][i] == fgroundvalue))
						{
							edge3DPtr[k][i] = insideValue; //check if pixel neighbors has foreground value
						}
					}
					else if ((i != (xDim - 1)) && (i != (dim2D - 1)) && ((i % xDim) == (xDim - 1)))//if at last column
					{//(excluding top and bottom right corner)
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) || (image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
						{
							edge3DPtr[k][i] = fgroundvalue; // check if pixel neighbors has background value
						}
						else if ((image3DPtr[k][i - 1] == fgroundvalue) || (image3DPtr[k][i + xDim] == fgroundvalue) || (image3DPtr[k][i - xDim] == fgroundvalue) || (image3DPtr[k + 1][i] == fgroundvalue))
						{
							edge3DPtr[k][i] = insideValue; //check if pixel neighbors has foreground value
						}
					}
					else {
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
						{
							edge3DPtr[k][i] = fgroundvalue; // check if pixel neighbors has background value
						}
						else if ((image3DPtr[k][i - 1] == fgroundvalue) || (image3DPtr[k][i + 1] == fgroundvalue) || (image3DPtr[k][i - xDim] == fgroundvalue) || (image3DPtr[k][i + xDim] == fgroundvalue) || (image3DPtr[k + 1][i] == fgroundvalue))
						{
							edge3DPtr[k][i] = insideValue; //check if pixel neighbors has foreground value
						}
					}
				}

				else if (k == (zDim - 1))
				{
					if (i == 0)//if at the top left corner,
					{
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue))
						{
							edge3DPtr[k][i] = fgroundvalue; // check if pixel neighbors has background value
						}
						else if ((image3DPtr[k][i + 1] == fgroundvalue) || (image3DPtr[k][i + xDim] == fgroundvalue) || (image3DPtr[k - 1][i] == fgroundvalue))
						{
							edge3DPtr[k][i] = insideValue; //check if pixel neighbors has foreground value
						}
					}
					else if (i == (xDim - 1)) //if at the top right coner,
					{
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue))
						{
							edge3DPtr[k][i] = fgroundvalue; // check if pixel neighbors has background value
						}
						else if ((image3DPtr[k][i - 1] == fgroundvalue) || (image3DPtr[k][i + xDim] == fgroundvalue) || (image3DPtr[k - 1][i] == fgroundvalue))
						{
							edge3DPtr[k][i] = insideValue; //check if pixel neighbors has foreground value
						}
					}
					else if (i == ((yDim - 1) * xDim)) //if at the bottom left corner,
					{
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue))
						{
							edge3DPtr[k][i] = fgroundvalue; // check if pixel neighbors has background value
						}
						else if ((image3DPtr[k][i + 1] == fgroundvalue) || (image3DPtr[k][i - xDim] == fgroundvalue) || (image3DPtr[k - 1][i] == fgroundvalue))
						{
							edge3DPtr[k][i] = insideValue; //check if pixel neighbors has foreground value
						}
					}
					else if (i == (dim2D - 1)) //if at the bottom right corner,
					{
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue))
						{
							edge3DPtr[k][i] = fgroundvalue; // check if pixel neighbors has background value
						}
						else if ((image3DPtr[k][i - 1] == fgroundvalue) || (image3DPtr[k][i - xDim] == fgroundvalue) || (image3DPtr[k - 1][i] == fgroundvalue))
						{
							edge3DPtr[k][i] = insideValue; //check if pixel neighbors has foreground value
						}
					}
					else if ((i > 0) && (i < (xDim - 1)))//if at first row (excluding left and right corner),
					{
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue))
						{
							edge3DPtr[k][i] = fgroundvalue; //check if pixel neighbors has background value
						}
						else if ((image3DPtr[k][i - 1] == fgroundvalue) || (image3DPtr[k][i + 1] == fgroundvalue) || (image3DPtr[k][i + xDim] == fgroundvalue) || (image3DPtr[k - 1][i] == fgroundvalue))
						{
							edge3DPtr[k][i] = insideValue; //check if pixel neighbors has foreground value
						}
					}
					else if ((i > ((yDim - 1) * xDim)) && (i < (dim2D - 1))) //if at last row
					{//(excluding left and right corner),
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue))
						{
							edge3DPtr[k][i] = fgroundvalue; // check if pixel neighbors has background value
						}
						else if ((image3DPtr[k][i - 1] == fgroundvalue) || (image3DPtr[k][i + 1] == fgroundvalue) || (image3DPtr[k][i - xDim] == fgroundvalue) || (image3DPtr[k - 1][i] == fgroundvalue))
						{
							edge3DPtr[k][i] = insideValue; //check if pixel neighbors has foreground value
						}
					}
					else if ((i != 0) && (i != ((yDim - 1) * xDim)) && ((i % xDim) == 0))//if at first column
					{//(excluding top and bottom left corner),
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) || (image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue))
						{
							edge3DPtr[k][i] = fgroundvalue; // check if pixel neighbors has background value
						}
						else if ((image3DPtr[k][i + 1] == fgroundvalue) || (image3DPtr[k][i + xDim] == fgroundvalue) || (image3DPtr[k][i - xDim] == fgroundvalue) || (image3DPtr[k - 1][i] == fgroundvalue))
						{
							edge3DPtr[k][i] = insideValue; //check if pixel neighbors has foreground value
						}
					}
					else if ((i != (xDim - 1)) && (i != (dim2D - 1)) && ((i % xDim) == (xDim - 1)))//if at last column
					{//(excluding top and bottom right corner),
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) || (image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue))
						{
							edge3DPtr[k][i] = fgroundvalue; // check if pixel neighbors has background value
						}
						else if ((image3DPtr[k][i - 1] == fgroundvalue) || (image3DPtr[k][i + xDim] == fgroundvalue) || (image3DPtr[k][i - xDim] == fgroundvalue) || (image3DPtr[k - 1][i] == fgroundvalue))
						{
							edge3DPtr[k][i] = insideValue; //check if pixel neighbors has foreground value
						}
					}
					else {//otherwise,
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue))
						{
							edge3DPtr[k][i] = fgroundvalue; // check if pixel neighbors has background value
						}
						else if ((image3DPtr[k][i - 1] == fgroundvalue) || (image3DPtr[k][i + 1] == fgroundvalue) || (image3DPtr[k][i - xDim] == fgroundvalue) || (image3DPtr[k][i + xDim] == fgroundvalue) || (image3DPtr[k - 1][i] == fgroundvalue))
						{
							edge3DPtr[k][i] = insideValue; //check if pixel neighbors has foreground value
						}
					}
				}

				else
				{
					if (i == 0)//if at the top left corner,
					{
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
						{
							edge3DPtr[k][i] = fgroundvalue; //check if pixel neighbors has background value
						}
						else if ((image3DPtr[k][i + 1] == fgroundvalue) || (image3DPtr[k][i + xDim] == fgroundvalue) || (image3DPtr[k - 1][i] == fgroundvalue) || (image3DPtr[k + 1][i] == fgroundvalue))
						{
							edge3DPtr[k][i] = insideValue; //check if pixel neighbors has foreground value
						}
					}
					else if (i == (xDim - 1)) //if at the top right coner,
					{
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
						{
							edge3DPtr[k][i] = fgroundvalue; //check if pixel neighbors has background value
						}
						else if ((image3DPtr[k][i - 1] == fgroundvalue) || (image3DPtr[k][i + xDim] == fgroundvalue) || (image3DPtr[k - 1][i] == fgroundvalue) || (image3DPtr[k + 1][i] == fgroundvalue))
						{
							edge3DPtr[k][i] = insideValue; //check if pixel neighbors has foreground value
						}
					}
					else if (i == ((yDim - 1) * xDim)) //if at the bottom left corner,
					{
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
						{
							edge3DPtr[k][i] = fgroundvalue; // //check if pixel neighbors has background value
						}
						else if ((image3DPtr[k][i + 1] == fgroundvalue) || (image3DPtr[k][i - xDim] == fgroundvalue) || (image3DPtr[k - 1][i] == fgroundvalue) || (image3DPtr[k + 1][i] == fgroundvalue))
						{
							edge3DPtr[k][i] = insideValue; //check if pixel neighbors has foreground value
						}
					}
					else if (i == (dim2D - 1)) //if at the bottom right corner,
					{
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
						{
							edge3DPtr[k][i] = fgroundvalue; // check if pixel neighbors has background value
						}
						else if ((image3DPtr[k][i - 1] == fgroundvalue) || (image3DPtr[k][i - xDim] == fgroundvalue) || (image3DPtr[k - 1][i] == fgroundvalue) || (image3DPtr[k + 1][i] == fgroundvalue))
						{
							edge3DPtr[k][i] = insideValue; //check if pixel neighbors has foreground value
						}
					}
					else if ((i > 0) && (i < (xDim - 1)))//if at first row (excluding left and right corner),
					{
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue))
						{
							edge3DPtr[k][i] = fgroundvalue; //check if pixel neighbors has background value
						}
						else if ((image3DPtr[k][i - 1] == fgroundvalue) || (image3DPtr[k][i + 1] == fgroundvalue) || (image3DPtr[k][i + xDim] == fgroundvalue) || (image3DPtr[k + 1][i] == fgroundvalue) || (image3DPtr[k - 1][i] == fgroundvalue))
						{
							edge3DPtr[k][i] = insideValue; //check if pixel neighbors has foreground value
						}
					}
					else if ((i > ((yDim - 1) * xDim)) && (i < (dim2D - 1))) //if at last row
					{//(excluding left and right corner),
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue))
						{
							edge3DPtr[k][i] = fgroundvalue; // check if pixel neighbors has background value
						}
						else if ((image3DPtr[k][i - 1] == fgroundvalue) || (image3DPtr[k][i + 1] == fgroundvalue) || (image3DPtr[k][i - xDim] == fgroundvalue) || (image3DPtr[k + 1][i] == fgroundvalue) || (image3DPtr[k - 1][i] == fgroundvalue))
						{
							edge3DPtr[k][i] = insideValue; //check if pixel neighbors has foreground value
						}
					}
					else if ((i != 0) && (i != ((yDim - 1) * xDim)) && ((i % xDim) == 0))//if at first column
					{//(excluding top and bottom left corner),
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) || (image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue))
						{
							edge3DPtr[k][i] = fgroundvalue; // check if pixel neighbors has background value
						}
						else if ((image3DPtr[k][i + 1] == fgroundvalue) || (image3DPtr[k][i + xDim] == fgroundvalue) || (image3DPtr[k][i - xDim] == fgroundvalue) || (image3DPtr[k + 1][i] == fgroundvalue) || (image3DPtr[k - 1][i] == fgroundvalue))
						{
							edge3DPtr[k][i] = insideValue; //check if pixel neighbors has foreground value
						}
					}
					else if ((i != (xDim - 1)) && (i != (dim2D - 1)) && ((i % xDim) == (xDim - 1)))//if at last column
					{//(excluding top and bottom right corner),
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) || (image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue))
						{
							edge3DPtr[k][i] = fgroundvalue; // check if pixel neighbors has background value
						}
						else if ((image3DPtr[k][i - 1] == fgroundvalue) || (image3DPtr[k][i + xDim] == fgroundvalue) || (image3DPtr[k][i - xDim] == fgroundvalue) || (image3DPtr[k + 1][i] == fgroundvalue) || (image3DPtr[k - 1][i] == fgroundvalue))
						{
							edge3DPtr[k][i] = insideValue; //check if pixel neighbors has foreground value
						}
					}
					else {//otherwise
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i - xDim] == bgroundvalue) || (image3DPtr[k][i + xDim] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
						{
							edge3DPtr[k][i] = fgroundvalue; // check if pixel neighbors has background value
						}
						else if ((image3DPtr[k][i - 1] == fgroundvalue) || (image3DPtr[k][i + 1] == fgroundvalue) || (image3DPtr[k][i - xDim] == fgroundvalue) || (image3DPtr[k][i + xDim] == fgroundvalue) || (image3DPtr[k - 1][i] == fgroundvalue) || (image3DPtr[k + 1][i] == fgroundvalue))
						{
							edge3DPtr[k][i] = insideValue; //check if pixel neighbors has foreground value
						}
					}
				}
			}
		}
	}
}


bool edgeDetection3dFunctionUC(unsigned char ** image3DPtr, unsigned char ** edge3DPtr, const size_t xDim, 
	const size_t yDim, 	const size_t zDim, const unsigned char bgroundvalue, const unsigned char fgroundvalue)
{
	size_t i, k;//loop counter for z dimension
	const size_t dim2D = xDim * yDim;

	//checks if the memory was allocated
	if (image3DPtr == NULL)
		return false;
	if (edge3DPtr == NULL)
		return false;

	// initialize the edge array to background value
	initialize3dArrayUC(edge3DPtr, xDim, yDim, zDim, bgroundvalue);

	// filling the 2d array with number=value
	for (k = 0; k < zDim; k++)
	{
		for (i = 0; i < dim2D; i++) {
			if (image3DPtr[k][i] == fgroundvalue)//while on the object
			{ 
				if (k == 0)
				{
					if (i == 0)//if at the top left corner, 
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue)
							|| (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if (i == (yDim - 1)) //if at the top right coner, 
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue)
							|| (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if (i == ((xDim - 1) * yDim)) //if at the bottom left corner, 
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i - yDim] == bgroundvalue)
							|| (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if (i == (dim2D - 1)) //if at the bottom right corner,
					{// check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i - yDim] == bgroundvalue)
							|| (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i > 0) && (i < (yDim - 1)))//if at first row (excluding left and right corner),
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) ||
							(image3DPtr[k][i + yDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i >((xDim - 1) * yDim)) && (i < (dim2D - 1))) //if at last row 
					{//(excluding left and right corner), check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) ||
							(image3DPtr[k][i - yDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i != 0) && (i != ((xDim - 1) * yDim)) && ((i % yDim) == 0))//if at first column 
					{//(excluding top and bottom left corner), check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue) ||
							(image3DPtr[k][i - yDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i != (yDim - 1)) && (i != (dim2D - 1)) && ((i % yDim) == (yDim - 1)))//if at last column 
					{//(excluding top and bottom right corner), check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue) ||
							(image3DPtr[k][i - yDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else {//otherwise, check if pixel neighbors has background value

						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) ||
							(image3DPtr[k][i - yDim] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue) ||
							(image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
				}
				
				else if (k == (zDim - 1))
				{
					if (i == 0)//if at the top left corner, 
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue)
							|| (image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if (i == (yDim - 1)) //if at the top right coner, 
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue)
							|| (image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if (i == ((xDim - 1) * yDim)) //if at the bottom left corner, 
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i - yDim] == bgroundvalue)
							|| (image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if (i == (dim2D - 1)) //if at the bottom right corner,
					{// check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i - yDim] == bgroundvalue)
							|| (image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i > 0) && (i < (yDim - 1)))//if at first row (excluding left and right corner),
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) ||
							(image3DPtr[k][i + yDim] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i >((xDim - 1) * yDim)) && (i < (dim2D - 1))) //if at last row 
					{//(excluding left and right corner), check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) ||
							(image3DPtr[k][i - yDim] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i != 0) && (i != ((xDim - 1) * yDim)) && ((i % yDim) == 0))//if at first column 
					{//(excluding top and bottom left corner), check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue) ||
							(image3DPtr[k][i - yDim] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i != (yDim - 1)) && (i != (dim2D - 1)) && ((i % yDim) == (yDim - 1)))//if at last column 
					{//(excluding top and bottom right corner), check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue) ||
							(image3DPtr[k][i - yDim] == bgroundvalue) || (image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else {//otherwise, check if pixel neighbors has background value

						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) ||
							(image3DPtr[k][i - yDim] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
				}

				else
				{
					if (i == 0)//if at the top left corner, 
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if (i == (yDim - 1)) //if at the top right coner, 
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if (i == ((xDim - 1) * yDim)) //if at the bottom left corner, 
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i - yDim] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if (i == (dim2D - 1)) //if at the bottom right corner,
					{// check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i - yDim] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i > 0) && (i < (yDim - 1)))//if at first row (excluding left and right corner),
					{//check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) ||
							(image3DPtr[k][i + yDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i >((xDim - 1) * yDim)) && (i < (dim2D - 1))) //if at last row 
					{//(excluding left and right corner), check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) ||
							(image3DPtr[k][i - yDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i != 0) && (i != ((xDim - 1) * yDim)) && ((i % yDim) == 0))//if at first column 
					{//(excluding top and bottom left corner), check if pixel neighbors has background value
						if ((image3DPtr[k][i + 1] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue) ||
							(image3DPtr[k][i - yDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else if ((i != (yDim - 1)) && (i != (dim2D - 1)) && ((i % yDim) == (yDim - 1)))//if at last column 
					{//(excluding top and bottom right corner), check if pixel neighbors has background value
						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue) ||
							(image3DPtr[k][i - yDim] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
					else {//otherwise, check if pixel neighbors has background value

						if ((image3DPtr[k][i - 1] == bgroundvalue) || (image3DPtr[k][i + 1] == bgroundvalue) ||
							(image3DPtr[k][i - yDim] == bgroundvalue) || (image3DPtr[k][i + yDim] == bgroundvalue) ||
							(image3DPtr[k - 1][i] == bgroundvalue) || (image3DPtr[k + 1][i] == bgroundvalue))
							edge3DPtr[k][i] = fgroundvalue;
						else
							edge3DPtr[k][i] = bgroundvalue;
					}
				}
			}

		}

	}

	return true;
}

//function for edgedetection in 2D.
bool edgeDetection2dFunctionUC(unsigned char * image2DPtr, unsigned char * edge2DPtr, const size_t xDim, 
	const size_t yDim, 	const unsigned char bgroundvalue, const unsigned char fgroundvalue)
{
	const size_t dim2D = xDim * yDim;
	size_t i;

	if (image2DPtr == NULL)
		return false;
	if (edge2DPtr == NULL)
		return false;

	// Edge detection in 2D
	for (i = 0; i < dim2D; i++) {
		if (image2DPtr[i] == fgroundvalue)//while on the object
		{
			if (i == 0)//if at the top left corner, 
			{//check if pixel neighbors has background value
				if ((image2DPtr[i + 1] == bgroundvalue) || (image2DPtr[i + yDim] == bgroundvalue))
					edge2DPtr[i] = fgroundvalue;
				else
					edge2DPtr[i] = bgroundvalue;
			}
			else if (i == (yDim - 1)) //if at the top right coner, 
			{//check if pixel neighbors has background value
				if ((image2DPtr[i - 1] == bgroundvalue) || (image2DPtr[i + yDim] == bgroundvalue))
					edge2DPtr[i] = fgroundvalue;
				else
					edge2DPtr[i] = bgroundvalue;
			}
			else if (i == ((xDim - 1) * yDim)) //if at the bottom left corner, 
			{//check if pixel neighbors has background value
				if ((image2DPtr[i + 1] == bgroundvalue) || (image2DPtr[i - yDim] == bgroundvalue))
					edge2DPtr[i] = fgroundvalue;
				else
					edge2DPtr[i] = bgroundvalue;
			}
			else if (i == (dim2D - 1)) //if at the bottom right corner,
			{// check if pixel neighbors has background value
				if ((image2DPtr[i - 1] == bgroundvalue) || (image2DPtr[i - yDim] == bgroundvalue))
					edge2DPtr[i] = fgroundvalue;
				else
					edge2DPtr[i] = bgroundvalue;
			}
			else if ((i > 0) && (i < (yDim - 1)))//if at first row (excluding left and right corner),
			{//check if pixel neighbors has background value
				if ((image2DPtr[i - 1] == bgroundvalue) || (image2DPtr[i + 1] == bgroundvalue) ||
					(image2DPtr[i + yDim] == bgroundvalue))
					edge2DPtr[i] = fgroundvalue;
				else
					edge2DPtr[i] = bgroundvalue;
			}
			else if ((i >((xDim - 1) * yDim)) && (i < (dim2D - 1))) //if at last row 
			{//(excluding left and right corner), check if pixel neighbors has background value
				if ((image2DPtr[i - 1] == bgroundvalue) || (image2DPtr[i + 1] == bgroundvalue) ||
					(image2DPtr[i - yDim] == bgroundvalue))
					edge2DPtr[i] = fgroundvalue;
				else
					edge2DPtr[i] = bgroundvalue;
			}
			else if ((i != 0) && (i != ((xDim - 1) * yDim)) && ((i % yDim) == 0))//if at first column 
			{//(excluding top and bottom left corner), check if pixel neighbors has background value
				if ((image2DPtr[i + 1] == bgroundvalue) || (image2DPtr[i + yDim] == bgroundvalue) ||
					(image2DPtr[i - yDim] == bgroundvalue))
					edge2DPtr[i] = fgroundvalue;
				else
					edge2DPtr[i] = bgroundvalue;
			}
			else if ((i != (yDim - 1)) && (i != (dim2D - 1)) && ((i % yDim) == (yDim - 1)))//if at last column 
			{//(excluding top and bottom right corner), check if pixel neighbors has background value
				if ((image2DPtr[i - 1] == bgroundvalue) || (image2DPtr[i + yDim] == bgroundvalue) ||
					(image2DPtr[i - yDim] == bgroundvalue))
					edge2DPtr[i] = fgroundvalue;
				else
					edge2DPtr[i] = bgroundvalue;
			}
			else {//otherwise, check if pixel neighbors has background value

				if ((image2DPtr[i - 1] == bgroundvalue) || (image2DPtr[i + 1] == bgroundvalue) ||
					(image2DPtr[i - yDim] == bgroundvalue) || (image2DPtr[i + yDim] == bgroundvalue))
					edge2DPtr[i] = fgroundvalue;
				else
					edge2DPtr[i] = bgroundvalue;
			}

		}

	}


	return true;
}

//Functions for Canny edge detector

void computeAngleFromGradient(dataType* anglePtr, dataType* gradientX, dataType* gradientY, const size_t length, const size_t width) {
	size_t i, j, xd;

	for (i = 0; i < length; i++) {
		for (j = 0; j < width; j++) {
			xd = x_new(i, j, length);
			anglePtr[xd] = atan(gradientY[xd] / (gradientX[xd] + 0.000001));
		}
	}

	for (i = 0; i < length * width; i++) {
		anglePtr[i] = anglePtr[i] * 180 / M_PI;
		if (anglePtr[i] < 0) {
			anglePtr[i] += 180;
		}
	}
}

void nonMaximumSuppression(dataType* resultPtr, dataType* normOfGradient, dataType* anglePtr, const size_t length, const size_t width) {
	
	dataType q, r;
	size_t i, j, xd;

	for (i = 1; i < length - 1; i++) {
		for (j = 1; j < width - 1; j++) {
			
			q = 1.0;
			r = 1.0;
			xd = x_new(i, j, length);

			if ((anglePtr[xd] >= 0 && anglePtr[xd] < 22.5) || (anglePtr[xd] >= 157.5 && anglePtr[xd] <= 180)) {
				r = normOfGradient[x_new(i, j - 1, length)];
				q = normOfGradient[x_new(i, j + 1, length)];
			}
			else {
				if (anglePtr[xd] >= 22.5 && anglePtr[xd] < 67.5) {
					r = normOfGradient[x_new(i - 1, j + 1, length)];
					q = normOfGradient[x_new(i + 1, j - 1, length)];
				}
				else {
					if (anglePtr[xd] >= 67.5 && anglePtr[xd] < 112.5) {
						r = normOfGradient[x_new(i - 1, j, length)];
						q = normOfGradient[x_new(i + 1, j - 1, length)];
					}
					else {
						if (anglePtr[xd] >= 112.5 && anglePtr[xd] < 157.5) {
							r = normOfGradient[x_new(i + 1, j + 1, length)];
							q = normOfGradient[x_new(i - 1, j - 1, length)];
						}
					}
				}
			}

			if (normOfGradient[xd] >= q && normOfGradient[xd] >= r) {
				resultPtr[xd] = normOfGradient[xd];
			}
			else {
				resultPtr[xd] = 0.0;
			}
		}
	}
}

void thresholdByHyteresis(dataType* resultPtr, dataType* normOfGradient, size_t* status, const size_t length, const size_t width, dataType thres_min, dataType thres_max) {
	size_t i, j, xd, count_neighbor = 0;

	for (i = 0; i < length * width; i++) {
		if (normOfGradient[i] > thres_max) {
			status[i] = 3;
			resultPtr[i] = 1.0;
		}
		else {
			if (normOfGradient[i] < thres_min) {
				status[i] = 1;
				resultPtr[i] = 0.0;
			}
			else {
				status[i] = 2;
			}
		}
	}

	for (i = 1; i < length - 1; i++) {
		for (j = 1; j < width - 1; j++) {
			
			xd = x_new(i, j, length);
			count_neighbor = 0;

			if (status[xd] == 2) {

				if (status[x_new(i - 1, j - 1, length)] == 3) {
					count_neighbor++;
				}
				if (status[x_new(i, j - 1, length)] == 3) {
					count_neighbor++;
				}
				if (status[x_new(i + 1, j - 1, length)] == 3) {
					count_neighbor++;
				}

				if (status[x_new(i - 1, j, length)] == 3) {
					count_neighbor++;
				}
				if (status[x_new(i + 1, j, length)] == 3) {
					count_neighbor++;
				}

				if (status[x_new(i - 1, j + 1, length)] == 3) {
					count_neighbor++;
				}
				if (status[x_new(i, j + 1, length)] == 3) {
					count_neighbor++;
				}
				if (status[x_new(i + 1, j + 1, length)] == 3) {
					count_neighbor++;
				}

				if (count_neighbor > 0) {
					resultPtr[xd] = 1.0;
				}
				else {
					resultPtr[xd] = 0.0;
				}

			}
		}
	}
}