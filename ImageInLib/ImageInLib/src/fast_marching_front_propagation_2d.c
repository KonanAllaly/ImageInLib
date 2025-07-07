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

dataType solve2dQuadratic(dataType dx, dataType dy, dataType p, PixelSpacing h) {

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

	if (dx == INFINITY && dy != INFINITY)
	{
		a = 1.0;
		b = -2 * dy;
		c = (dataType)(dy * dy - hy_2 * p_2);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0) {
			solution = (dataType)((-b + sqrt(delta)) / (2 * a));
			if (solution >= dy) {
				return solution;
			}
			else {
				return (dataType)(dy + hy * p);
			}
		}
		else {
			return (dataType)(dy + hy * p);
		}
	}

	if (dx != INFINITY && dy == INFINITY)
	{
		a = 1.0;
		b = -2 * dx;
		c = (dataType)(dx * dx - hx_2 * p_2);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0)
		{
			solution = (dataType)((-b + sqrt(delta)) / (2 * a));
			if (solution >= dx)
			{
				return solution;
			}
			else {
				return (dataType)(dx + hx * p);
			}
		}
		else {
			return (dataType)(dx + hx * p);
		}
	}

	if (dx != INFINITY && dy != INFINITY)
	{
		a = hx_2 + hy_2;
		b = -2 * (hy_2 * dx + hx_2 * dy);
		c = (dataType)(hy_2 * dx * dx + hx_2 * dy * dy - hx_2 * hy_2 * p_2);
		delta = (dataType)(b * b - 4 * a * c);
		if (delta >= 0)
		{
			solution = (dataType)((-b + sqrt(delta)) / (2 * a));
			if (solution >= fmax(dx, dy)) {
				return solution;
			}
			else {
				return (dataType)(fmin(dx + hx * p, dy + hy * p));
			}
		}
		else {
			return (dataType)(fmin(dx + hx * p, dy + hy * p));
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

	//std::vector<pointFastMarching2D> narrowBand;
	labelingList narrowBand = createlabelingList();
	labelingPoint* current_point = (labelingPoint*)malloc(sizeof(labelingPoint));
	if(current_point == NULL)
	{
		return;
	}

	size_t i = (size_t)endPoint[0].x;
	size_t j = (size_t)endPoint[0].y;
	size_t currentIndx = x_new(i, j, length);

	current_point->x = endPoint[0].x;
	current_point->y = endPoint[1].y;
	current_point->z = -1;//In 2D : set the third coordinate to -1
	current_point->arrival = 0.0;
	current_point->previous = NULL;
	current_point->next = NULL;
	narrowBand.first_point = current_point;

	//STEP 1
	//In labelAray we have : 1 ---> already processed, 2 ---> narrow band and 3 ---> not processed
	for (size_t k = 0; k < dim2D; k++) {
		action[k] = INFINITY;
		labelArray[k] = 3;
	}

	action[currentIndx] = 0.0;
	labelArray[currentIndx] = 1;

	//East
	if (j < width - 1 && i >= 0 && i < length) {
		size_t jplus = j + 1;
		//dataType x = selectX(action, length, width, i, jplus);
		//dataType y = selectY(action, length, width, i, jplus);
		size_t indxEast = x_new(i, jplus, length);
		dataType coefSpeed = potential[indxEast];
		//dataType dEast = solve2dQuadratic(x, y, coefSpeed, spacing);
		//pointFastMarching2D EastNeighbor = { i, jplus, dEast };
		//action[indxEast] = dEast;
		//narrowBand.push_back(EastNeighbor);
		labelArray[indxEast] = 2;
	}

	/*
	//West
	if (j > 0 && i >= 0 && i < length) {
		size_t jminus = j - 1;
		dataType x = selectX(actionMapPtr, length, width, i, jminus);
		dataType y = selectY(actionMapPtr, length, width, i, jminus);
		size_t indxWest = x_new(i, jminus, length);
		dataType coefSpeed = potentialPtr[indxWest];
		dataType dWest = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D WestNeighbor = { i, jminus, dWest };
		actionMapPtr[indxWest] = dWest;
		narrowBand.push_back(WestNeighbor);
		labelArray[indxWest] = 2;
	}

	//North
	if (j >= 0 && j < width && i > 0) {
		size_t iminus = i - 1;
		dataType x = selectX(actionMapPtr, length, width, iminus, j);
		dataType y = selectY(actionMapPtr, length, width, iminus, j);
		size_t indxNorth = x_new(iminus, j, length);
		dataType coefSpeed = potentialPtr[indxNorth];
		dataType dNorth = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D NorthNeighbor = { iminus, j, dNorth };
		actionMapPtr[indxNorth] = dNorth;
		narrowBand.push_back(NorthNeighbor);
		labelArray[indxNorth] = 2;
	}

	//South
	if (j >= 0 && j < width && i < length - 1) {
		size_t iplus = i + 1;
		dataType x = selectX(actionMapPtr, length, width, iplus, j);
		dataType y = selectY(actionMapPtr, length, width, iplus, j);
		size_t indxSouth = x_new(iplus, j, length);
		dataType coefSpeed = potentialPtr[indxSouth];
		dataType dSouth = solve2dQuadratic(x, y, coefSpeed, spacing);
		pointFastMarching2D SouthNeighbor = { iplus, j, dSouth };
		actionMapPtr[indxSouth] = dSouth;
		narrowBand.push_back(SouthNeighbor);
		labelArray[indxSouth] = 2;
	}

	//heapify 2D vector
	heapifyVector2D(narrowBand);

	//Save points for visualization
	size_t x_final_point = (size_t)endPoints[1].x;
	size_t y_final_point = (size_t)endPoints[1].y;

	dataType max_save_action = 0.0;

	while (narrowBand.size() > 0) {

		pointFastMarching2D current = narrowBand[0];
		size_t i = current.x;
		size_t j = current.y;
		size_t currentIndx = x_new(i, j, length);
		labelArray[currentIndx] = 1;

		// Exit when the final point is reached by the front
		if (i == x_final_point && j == y_final_point) {
			labelArray[currentIndx] = 1; // Mark the seed point as processed
			break;
		}

		if (actionMapPtr[currentIndx] > max_save_action) {
			max_save_action = actionMapPtr[currentIndx];
		}

		deleteRootHeap2D(narrowBand);

		//West
		if (i > 0 && i < length && j >= 0 && j < width) {
			size_t iminus = i - 1;
			size_t indxWest = x_new(iminus, j, length);
			short label = labelArray[indxWest];
			if (label != 1) {
				dataType x = selectX(actionMapPtr, length, width, iminus, j);
				dataType y = selectY(actionMapPtr, length, width, iminus, j);
				dataType coefSpeed = potentialPtr[indxWest];
				dataType dWest = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D WestNeighbor = { iminus, j, dWest };
				if (label == 3) {
					actionMapPtr[indxWest] = dWest;
					labelArray[indxWest] = 2;
					addPointHeap2D(narrowBand, WestNeighbor);
				}
				else {
					if (dWest < actionMapPtr[indxWest]) {
						actionMapPtr[indxWest] = dWest;
						size_t pIndex = getIndexFromHeap2D(narrowBand, iminus, j);
						if (pIndex != -1) {
							heapifyUp2D(narrowBand, pIndex);
						}
					}
				}
			}
		}

		//East
		if (i >= 0 && i < length_minus && j >= 0 && j < width) {
			size_t iplus = i + 1;
			size_t indxEast = x_new(iplus, j, length);
			short label = labelArray[indxEast];
			if (label != 1) {
				dataType x = selectX(actionMapPtr, length, width, iplus, j);
				dataType y = selectY(actionMapPtr, length, width, iplus, j);
				dataType coefSpeed = potentialPtr[indxEast];
				dataType dEast = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D EastNeighbor = { iplus, j, dEast };
				if (label == 3) {
					actionMapPtr[indxEast] = dEast;
					labelArray[indxEast] = 2;
					addPointHeap2D(narrowBand, EastNeighbor);
				}
				else {
					if (dEast < actionMapPtr[indxEast]) {
						actionMapPtr[indxEast] = dEast;
						size_t pIndex = getIndexFromHeap2D(narrowBand, iplus, j);
						if (pIndex != -1) {
							heapifyUp2D(narrowBand, pIndex);
						}
					}
				}
			}
		}

		//North
		if (j > 0 && j < width && i >= 0 && i < length) {
			size_t jminus = j - 1;
			size_t indxNorth = x_new(i, jminus, length);
			short label = labelArray[indxNorth];
			if (label != 1) {
				dataType x = selectX(actionMapPtr, length, width, i, jminus);
				dataType y = selectY(actionMapPtr, length, width, i, jminus);
				dataType coefSpeed = potentialPtr[indxNorth];
				dataType dNorth = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D NorthNeighbor = { i, jminus, dNorth };
				if (label == 3) {
					actionMapPtr[indxNorth] = dNorth;
					labelArray[indxNorth] = 2;
					addPointHeap2D(narrowBand, NorthNeighbor);
				}
				else {
					if (dNorth < actionMapPtr[indxNorth]) {
						actionMapPtr[indxNorth] = dNorth;
						size_t pIndex = getIndexFromHeap2D(narrowBand, i, jminus);
						if (pIndex != -1) {
							heapifyUp2D(narrowBand, pIndex);
						}
					}
				}
			}
		}

		//South
		if (j >= 0 && j < width_minus && i >= 0 && i < length) {
			size_t jplus = j + 1;
			size_t indxSouth = x_new(i, jplus, length);
			short label = labelArray[indxSouth];
			if (label != 1) {
				dataType x = selectX(actionMapPtr, length, width, i, jplus);
				dataType y = selectY(actionMapPtr, length, width, i, jplus);
				dataType coefSpeed = potentialPtr[indxSouth];
				dataType dSouth = solve2dQuadratic(x, y, coefSpeed, spacing);
				pointFastMarching2D SouthNeighbor = { i, jplus, dSouth };
				if (label == 3) {
					actionMapPtr[indxSouth] = dSouth;
					labelArray[indxSouth] = 2;
					addPointHeap2D(narrowBand, SouthNeighbor);
				}
				else {
					if (dSouth < actionMapPtr[indxSouth]) {
						actionMapPtr[indxSouth] = dSouth;
						size_t pIndex = getIndexFromHeap2D(narrowBand, i, jplus);
						if (pIndex != -1) {
							heapifyUp2D(narrowBand, pIndex);
						}
					}
				}
			}
		}
	}

	//Set the action of the non-processed points to maximum action value + 1
	for (size_t k = 0; k < dim2D; k++) {
		if (labelArray[k] == 3) {
			actionMapPtr[k] = max_save_action + 1;
		}
	}
	narrowBand.clear();
	*/

	free(current_point);
	free(labelArray);
	return true;
}