#include <stdlib.h>
#include <math.h>
#include "fast_marching_front_propagation_2d.h"

dataType upwindFiniteDifference2dX(dataType* action, const size_t length, const size_t width, const size_t ind_x, const size_t ind_y) {
	dataType x_minus, x_plus;
	if (ind_x == 0) {
		x_minus = INFINITY;
	}
	else {
		x_minus = action[x_new(ind_x - 1, ind_y, length)];
	}
	if (ind_x == length - 1) {
		x_plus = INFINITY;
	}
	else {
		x_plus = action[x_new(ind_x + 1, ind_y, length)];
	}
	return fmin(x_minus, x_plus);
}

dataType upwindFiniteDifference2dY(dataType* action, const size_t length, const size_t width, const size_t ind_x, const size_t ind_y) {
	dataType y_minus, y_plus;
	if (ind_y == 0) {
		y_minus = INFINITY;
	}
	else {
		y_minus = action[x_new(ind_x, ind_y - 1, length)];
	}
	if (ind_y == width - 1) {
		y_plus = INFINITY;
	}
	else {
		y_plus = action[x_new(ind_x, ind_y + 1, length)];
	}
	return fmin(y_minus, y_plus);
}

dataType solve2dQuadratic(dataType ux, dataType uy, dataType p, PixelSpacing h) {

	dataType solution = 0.0, a = 0.0, b = 0.0, c = 0.0, delta = 0.0;
	dataType p_2 = p * p;

	dataType hx = h.sx;
	dataType hx_2 = h.sx * h.sx;

	dataType hy = h.sy;
	dataType hy_2 = h.sy * h.sy;

	if (p <= 0.0) {
		printf("Error: P must be positive.");
		return INFINITY; // Return 0 if P is not positive
	}

	if (ux == INFINITY && uy != INFINITY)
	{
		a = 1.0;
		b = -2 * uy;
		c = (dataType)(uy * uy - hy_2 * p_2);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			solution = (dataType)((-b + sqrt(delta)) / (2 * a));
			if (solution >= uy) {
				return solution;
			}
			else {
				return (dataType)(uy + hy * p);
			}
		}
		else {
			return (dataType)(uy + hy * p);
		}
	}

	if (ux != INFINITY && uy == INFINITY)
	{
		a = 1.0;
		b = -2 * ux;
		c = (dataType)(ux * ux - hx_2 * p_2);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0)
		{
			solution = (dataType)((-b + sqrt(delta)) / (2 * a));
			if (solution >= ux)
			{
				return solution;
			}
			else {
				return (dataType)(ux + hx * p);
			}
		}
		else {
			return (dataType)(ux + hx * p);
		}
	}

	if (ux != INFINITY && uy != INFINITY)
	{
		a = hx_2 + hy_2;
		b = -2 * (hy_2 * ux + hx_2 * uy);
		c = (dataType)(hy_2 * ux * ux + hx_2 * uy * uy - hx_2 * hy_2 * p_2);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0)
		{
			solution = (dataType)((-b + sqrt(delta)) / (2 * a));
			if (solution >= fmax(ux, uy)) {
				return solution;
			}
			else {
				return (dataType)(fmin(ux + hx * p, uy + hy * p));
			}
		}
		else {
			return (dataType)(fmin(ux + hx * p, uy + hy * p));
		}
	}

}

void fastMarchingFrontPropagation2D(void* inputImageData, void* actionPtr, void* potentialPtr, void* endPoints, const PropagationType pType)
{
	switch (pType)
	{
	case PARTIAL_FRONT_PROPAGATION:
	{
		Image_Data2D inputImage = *(Image_Data2D*)inputImageData;
		dataType* action = (dataType*)actionPtr;
		dataType* potential = (dataType*)potentialPtr;
		Point2D* endPoint = (Point2D*)endPoints;
		partialFrontPropagation2D(inputImage, action, potential, endPoint);
		break;
	}
	default:
		break;
	}
}

bool partialFrontPropagation2D(Image_Data2D inputImage, dataType* action, dataType* potential, Point2D* endPoint)
{
	if(inputImage.imageDataPtr == NULL || action == NULL || potential == NULL || endPoint == NULL)
	{
		return false;
	}

	const size_t length = inputImage.height;
	const size_t width = inputImage.width;
	PixelSpacing spacing = inputImage.spacing;

	size_t dim2D = length * width;
	size_t length_minus = length - 1;
	size_t width_minus = width - 1;

	short* labelArray = (short*)malloc(sizeof(short) * dim2D);
	if(labelArray == NULL)
	{
		return false;
	}

	heapStructure* narrowBand = createHeap(dim2D);

	size_t i = (size_t)endPoint[0].x;
	size_t j = (size_t)endPoint[0].y;
	size_t currentIndx = x_new(i, j, length);
	if(i < 0 || i > length || j < 0 || j > width)
	{
		false;
	}

	//STEP 1
	//In labelAray we have : 1 ---> already processed, 2 ---> narrow band and 3 ---> not processed
	for (size_t k = 0; k < dim2D; k++) {
		action[k] = INFINITY;
		labelArray[k] = 3;
	}

	action[currentIndx] = 0.0;
	labelArray[currentIndx] = 1;

	//North
	if (i > 0) {
		size_t iminus = i - 1;
		dataType ux = upwindFiniteDifference2dX(action, length, width, iminus, j);
		dataType uy = upwindFiniteDifference2dY(action, length, width, iminus, j);
		size_t indxNorth = x_new(iminus, j, length);
		dataType coefSpeed = potential[indxNorth];
		dataType dNorth = solve2dQuadratic(ux, uy, coefSpeed, spacing);
		pointFastMarching NorthNeighbor = { iminus, j, -1, dNorth, indxNorth };
		pushToHeap(&narrowBand, NorthNeighbor);
		action[indxNorth] = dNorth;
		labelArray[indxNorth] = 2;
	}

	//South
	if (i < length_minus) {
		size_t iplus = i + 1;
		dataType ux = upwindFiniteDifference2dX(action, length, width, iplus, j);
		dataType uy = upwindFiniteDifference2dY(action, length, width, iplus, j);
		size_t indxSouth = x_new(iplus, j, length);
		dataType coefSpeed = potential[indxSouth];
		dataType dSouth = solve2dQuadratic(ux, uy, coefSpeed, spacing);
		pointFastMarching SouthNeighbor = { iplus, j, -1, dSouth, indxSouth };
		pushToHeap(&narrowBand, SouthNeighbor);
		action[indxSouth] = dSouth;
		labelArray[indxSouth] = 2;
	}

	//East
	if (j < width_minus) {
		size_t jplus = j + 1;
		dataType ux = upwindFiniteDifference2dX(action, length, width, i, jplus);
		dataType uy = upwindFiniteDifference2dY(action, length, width, i, jplus);
		size_t indxEast = x_new(i, jplus, length);
		dataType coefSpeed = potential[indxEast];
		dataType dEast = solve2dQuadratic(ux, uy, coefSpeed, spacing);
		pointFastMarching EastNeighbor = {i, jplus, -1, dEast, indxEast};
		pushToHeap(&narrowBand, EastNeighbor);
		action[indxEast] = dEast;
		labelArray[indxEast] = 2;
	}

	//West
	if (j > 0) {
		size_t jminus = j - 1;
		dataType ux = upwindFiniteDifference2dX(action, length, width, i, jminus);
		dataType uy = upwindFiniteDifference2dY(action, length, width, i, jminus);
		size_t indxWest = x_new(i, jminus, length);
		dataType coefSpeed = potential[indxWest];
		dataType dWest = solve2dQuadratic(ux, uy, coefSpeed, spacing);
		pointFastMarching WestNeighbor = { i, jminus, -1, dWest, indxWest };
		pushToHeap(&narrowBand, WestNeighbor);
		action[indxWest] = dWest;
		labelArray[indxWest] = 2;
	}

	dataType max_save_action = 0.0;

	size_t x_final_point = endPoint[1].x;
	size_t y_final_point = endPoint[1].y;
	if (x_final_point < 0 || x_final_point > length || y_final_point < 0 || y_final_point > width)
	{
		false;
	}

	while (narrowBand->size > 0) {

		pointFastMarching current = getPointWithMinimalArrival(&narrowBand);
		size_t i = current.x;
		size_t j = current.y;
		size_t currentIndx = x_new(i, j, length);
		labelArray[currentIndx] = 1;

		// Exit when the final point is reached by the front
		if (i == x_final_point && j == y_final_point) {
			labelArray[currentIndx] = 1;
			break;
		}

		if (action[currentIndx] > max_save_action) {
			max_save_action = action[currentIndx];
		}

		//West
		if (i > 0) {
			size_t iminus = i - 1;
			size_t indxWest = x_new(iminus, j, length);
			short label = labelArray[indxWest];
			if (label != 1) {
				dataType ux = upwindFiniteDifference2dX(action, length, width, iminus, j);
				dataType uy = upwindFiniteDifference2dY(action, length, width, iminus, j);
				dataType coefSpeed = potential[indxWest];
				dataType dWest = solve2dQuadratic(ux, uy, coefSpeed, spacing);
				pointFastMarching WestNeighbor = { iminus, j, -1, dWest, dWest };
				if (label == 3) {
					pushToHeap(&narrowBand, WestNeighbor);
					action[indxWest] = dWest;
					labelArray[indxWest] = 2;
				}
				else {
					if (dWest < action[indxWest]) {
						action[indxWest] = dWest;
						int pIndex = getPointPosition(&narrowBand, indxWest);
						if (pIndex != -1) {
							heapifyUp(&narrowBand, pIndex);
						}
					}
				}
			}
		}

		//East
		if (i < length_minus) {
			size_t iplus = i + 1;
			size_t indxEast = x_new(iplus, j, length);
			short label = labelArray[indxEast];
			if (label != 1) {
				dataType ux = upwindFiniteDifference2dX(action, length, width, iplus, j);
				dataType uy = upwindFiniteDifference2dY(action, length, width, iplus, j);
				dataType coefSpeed = potential[indxEast];
				dataType dEast = solve2dQuadratic(ux, uy, coefSpeed, spacing);
				pointFastMarching EastNeighbor = { iplus, j, -1, dEast, indxEast };
				if (label == 3) {
					pushToHeap(&narrowBand, EastNeighbor);
					action[indxEast] = dEast;
					labelArray[indxEast] = 2;
				}
				else {
					if (dEast < action[indxEast]) {
						action[indxEast] = dEast;
						int pIndex = getPointPosition(&narrowBand, indxEast);
						if (pIndex != -1) {
							heapifyUp(&narrowBand, pIndex);
						}
					}
				}
			}
		}

		//North
		if (j > 0) {
			size_t jminus = j - 1;
			size_t indxNorth = x_new(i, jminus, length);
			short label = labelArray[indxNorth];
			if (label != 1) {
				dataType ux = upwindFiniteDifference2dX(action, length, width, i, jminus);
				dataType uy = upwindFiniteDifference2dY(action, length, width, i, jminus);
				dataType coefSpeed = potential[indxNorth];
				dataType dNorth = solve2dQuadratic(ux, uy, coefSpeed, spacing);
				pointFastMarching NorthNeighbor = { i, jminus, -1, dNorth, indxNorth };
				if (label == 3) {
					pushToHeap(&narrowBand, NorthNeighbor);
					action[indxNorth] = dNorth;
					labelArray[indxNorth] = 2;
				}
				else {
					if (dNorth < action[indxNorth]) {
						action[indxNorth] = dNorth;
						int pIndex = getPointPosition(&narrowBand, indxNorth);
						if (pIndex != -1) {
							heapifyUp(&narrowBand, pIndex);
						}
					}
				}
			}
		}

		//South
		if (j < width_minus) {
			size_t jplus = j + 1;
			size_t indxSouth = x_new(i, jplus, length);
			short label = labelArray[indxSouth];
			if (label != 1) {
				dataType ux = upwindFiniteDifference2dX(action, length, width, i, jplus);
				dataType uy = upwindFiniteDifference2dY(action, length, width, i, jplus);
				dataType coefSpeed = potential[indxSouth];
				dataType dSouth = solve2dQuadratic(ux, uy, coefSpeed, spacing);
				pointFastMarching SouthNeighbor = { i, jplus, -1, dSouth, indxSouth };
				if (label == 3) {
					pushToHeap(&narrowBand, SouthNeighbor);
					action[indxSouth] = dSouth;
					labelArray[indxSouth] = 2;
				}
				else {
					if (dSouth < action[indxSouth]) {
						action[indxSouth] = dSouth;
						int pIndex = getPointPosition(&narrowBand, indxSouth);
						if (pIndex != -1) {
							heapifyUp(&narrowBand, pIndex);
						}
					}
				}
			}
		}
	}

	//Set the action of the non-processed points to maximum action value + 1
	for (size_t k = 0; k < dim2D; k++) {
		if (labelArray[k] == 3) {
			action[k] = max_save_action + 1;
		}
	}

	//release memory
	free(narrowBand->data);
	free(narrowBand);
	free(labelArray);

	return true;
}