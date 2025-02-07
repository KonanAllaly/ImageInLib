
#include <iostream>
#include <climits>
#include <crtdbg.h>
#include <corecrt_malloc.h>
#include <omp.h>
#include<vector>
#include <common_vtk.h>
#include "Labelling.h"
#include "morphological_change.h"
#include "imageInterpolation.h"
#include "../distanceForPathFinding.h"
#include "../src/distance_function.h"
#include "template_functions.h"
#include "thresholding.h"

using namespace std;

//for 2D image using 2D arrays
bool labelling2D(int** imageDataPtr, int** segmentedImage, bool** statusArray, const size_t length, const size_t width, int object) {

	if (imageDataPtr == NULL || segmentedImage == NULL || statusArray == NULL)
		return false;

	vector<int> iStack, jStack, iTmpStack, jTmpStack;
	int i = 0, j = 0, iNew = 0, jNew = 0, k = 0, xd = 0, label = 1;

	//statusArray has false everywhere at the beginning
	//false---> non-processed and true---> processed
	//segmentedImage has 0 everywhere at the begining
	for (i = 0; i < length; i++) {
		for (j = 0; j < width; j++) {

			if (statusArray[i][j] == false) {
				if (imageDataPtr[i][j] == object) {
					//top neighbor
					if (i > 0 && statusArray[i - 1][j] == false) {
						if (imageDataPtr[i - 1][j] == object) {
							//element is in region add its coordinates in stacks
							iStack.push_back(i - 1);
							jStack.push_back(j);
						}
						else {
							//it's not object, so, update its status
							statusArray[i - 1][j] = true;
						}
					}
					//right neighbor
					if (j < width - 1 && statusArray[i][j + 1] == false) {
						if (imageDataPtr[i][j + 1] == object) {
							//element is in region add its coordinates in stacks
							iStack.push_back(i);
							jStack.push_back(j + 1);
						}
						else {
							statusArray[i][j + 1] = true;
						}
					}
					//bottom neighbor
					if (i < length - 1 && statusArray[i + 1][j] == false) {
						if (imageDataPtr[i + 1][j] == object) {
							//if element is in region add its coodinates in stacks
							iStack.push_back(i + 1);
							jStack.push_back(j);
						}
						else {
							statusArray[i + 1][j] = true;
						}
					}
					//left neighbor
					if (j > 0 && statusArray[i][j - 1] == false) {
						if (imageDataPtr[i][j - 1] == object) {
							//element is in region add its coodinates in stacks
							iStack.push_back(i);
							jStack.push_back(j - 1);
						}
						else {
							statusArray[i][j - 1] = true;
						}
					}
					//after checking all neighbors of current element
					//we give its label, and update its status
					segmentedImage[i][j] = label;
					statusArray[i][j] = true;
					//Now we check neighbors of its neighbors saved in the stacks
					//we start by the last added element in stacks
					//If there is no neighbor iStack.size() = jStack.size() = 0, and the while loop will no be ran
					while (iStack.size() > 0 && jStack.size() > 0) { //One is enought because they have same size
						iNew = iStack.size() - 1;
						jNew = jStack.size() - 1;
						//top neighbor
						if (iStack[iNew] > 0 && statusArray[ iStack[iNew] - 1 ][ jStack[jNew] ] == false) {
							if (imageDataPtr[iStack[iNew] - 1][jStack[jNew]] == object) {
								//If the element is in the current region, then save its coordinates in temporary stacks
								iTmpStack.push_back(iStack[iNew] - 1);
								jTmpStack.push_back(jStack[jNew]);
							}
							else {
								//Not in region, update status and go to the next neighbor
								statusArray[iStack[iNew] - 1][jStack[jNew]] = true;
							}
						}
						//right neighbor
						if (jStack[jNew] < width - 1 && statusArray[ iStack[iNew] ][ jStack[jNew] + 1 ] == false) {
							if (imageDataPtr[ iStack[iNew] ][ jStack[jNew] + 1 ] == object) {
								//If the element is in the current region, then save its coordinates in temporary stacks
								iTmpStack.push_back(iStack[iNew]);
								jTmpStack.push_back(jStack[jNew] + 1);
							}
							else {
								statusArray[ iStack[iNew] ][ jStack[jNew] + 1 ] = true;
							}
						}
						//bottom neighbor
						if (iStack[iNew] < length - 1 && statusArray[ iStack[iNew] + 1 ][ jStack[jNew] ] == false) {
							if (imageDataPtr[ iStack[iNew] + 1 ][ jStack[jNew] ] == object) {
								//If the element is in the current region, then save its coordinates in temporary stacks
								iTmpStack.push_back(iStack[iNew] + 1);
								jTmpStack.push_back(jStack[jNew]);
							}
							else {
								statusArray[ iStack[iNew] + 1 ][ jStack[jNew] ] = true;
							}
						}
						//left neighbor
						if (jStack[jNew] > 0 && statusArray[ iStack[iNew] ][ jStack[jNew] - 1 ] == false) {
							if (imageDataPtr[ iStack[iNew] ][ jStack[jNew] - 1 ] == object) {
								//If the element is in the current region, then save its coordinates in temporary stacks
								iTmpStack.push_back(iStack[iNew]);
								jTmpStack.push_back(jStack[jNew] - 1);
							}
							else {
								statusArray[ iStack[iNew] ][ jStack[jNew] - 1 ] = true;
							}
						}
						//updating of processed element befor removal
						segmentedImage[ iStack[iNew] ][ jStack[jNew] ] = label;
						statusArray[ iStack[iNew] ][ jStack[jNew] ] = true;
						//Remove the processed element of the initial stacks
						iStack.pop_back();
						jStack.pop_back();
						//Add new found neighbors in initial stacks
						//we can do it once because they have the same size
						for (k = 0; k < iTmpStack.size(); k++) {
							//if any neighbors have been found iTmpStack.size() = 0
							//and nothing will happen
							iStack.push_back(iTmpStack[k]);
						}
						for (k = 0; k < jTmpStack.size(); k++) {
							//if any neighbors have been found jTmpStack.size() = 0
							//and nothing will happen
							jStack.push_back(jTmpStack[k]);
						}
						//empty the temporary stacks
						//we can do it once because they have the same size
						while (iTmpStack.size() > 0) {
							iTmpStack.pop_back();
						}
						while (jTmpStack.size() > 0) {
							jTmpStack.pop_back();
						}
					}
					label++;
				}
				else {
					statusArray[i][j] = true;
				}
			}
		}
	}
	return true;
}

bool labeling(dataType* imageDataPtr, int* labelArray, bool* statusArray, const size_t length, const size_t width, dataType object) {

	if (imageDataPtr == NULL || labelArray == NULL || statusArray == NULL)
		return false;

	int label = 1;
	size_t i, j, k, xd;

	vector< point2dLabelling> initialStack, tempStack;

	//statusArray has false everywhere at the beginning
	//false---> non-processed and true---> processed
	//segmentedImage has 0 everywhere at the begining
	for (i = 0; i < length; i++) {
		for (j = 0; j < width; j++) {

			xd = x_new(i, j, length);

			if (statusArray[xd] == false) {
				if (imageDataPtr[xd] == object) {

					//East
					if (i > 0 && statusArray[x_new(i - 1, j, length)] == false) {
						if (imageDataPtr[x_new(i - 1, j, length)] == object) {
							//object pixel, so added to initial stack
							point2dLabelling current_point = { i - 1, j, false };
							initialStack.push_back(current_point);
						}
						else {
							//not object, update its status
							statusArray[x_new(i - 1, j, length)] = true;
						}
					}
					//South
					if (j < width - 1 && statusArray[x_new(i, j + 1, length)] == false) {
						if (imageDataPtr[x_new(i, j + 1, length)] == object) {
							//object pixel, so added to initial stack
							point2dLabelling current_point = { i, j + 1, false };
							initialStack.push_back(current_point);
						}
						else {
							//not object, update its status
							statusArray[x_new(i, j + 1, length)] = true;
						}
					}
					//West
					if (i < length - 1 && statusArray[x_new(i + 1, j, length)] == false) {
						if (imageDataPtr[x_new(i + 1, j, length)] == object) {
							//object pixel, so added to initial stack
							point2dLabelling current_point = { i + 1, j, false };
							initialStack.push_back(current_point);
						}
						else {
							//not object, update its status
							statusArray[x_new(i + 1, j, length)] = true;
						}
					}
					//North
					if (j > 0 && statusArray[x_new(i, j - 1, length)] == false) {
						if (imageDataPtr[x_new(i, j - 1, length)] == object) {
							//object pixel, so added to initial stack
							point2dLabelling current_point = { i, j - 1, false };
							initialStack.push_back(current_point);
						}
						else {
							//not object, update its status
							statusArray[x_new(i, j - 1, length)] = true;
						}
					}

					//after checking all neighbors of current element
					//we set its label, and update its status
					labelArray[xd] = label;
					statusArray[xd] = true;

					//Now we check neighbors of elemts in initial stack
					//we start by the last added element in stacks
					while (initialStack.size() > 0) {

						//Processed the last element in the initial stack
						size_t l = initialStack.size() - 1; // len
						size_t i_s = initialStack[l].x;
						size_t j_s = initialStack[l].y;
						size_t xd_s = x_new(i_s, j_s, length);

						//East
						if (i_s > 0 && statusArray[x_new(i_s - 1, j_s, length)] == false) {
							if (imageDataPtr[x_new(i_s - 1, j_s, length)] == object) {
								//object pixel, so added to temporary stack
								point2dLabelling current_point = { i_s - 1, j_s, false };
								tempStack.push_back(current_point);
							}
							else {
								//Not in region, update status and go to the next neighbor
								statusArray[x_new(i_s - 1, j_s, length)] = true;
							}
						}

						//South
						if (j_s < width - 1 && statusArray[x_new(i_s, j_s + 1, length)] == false) {
							if (imageDataPtr[x_new(i_s, j_s + 1, length)] == object) {
								//object pixel, so added to temporary stack
								point2dLabelling current_point = { i_s, j_s + 1, false };
								tempStack.push_back(current_point);
							}
							else {
								statusArray[x_new(i_s, j_s + 1, length)] = true;
							}
						}

						//East
						if (i_s < length - 1 && statusArray[x_new(i_s + 1, j_s, length)] == false) {
							if (imageDataPtr[x_new(i_s + 1, j_s, length)] == object) {
								//object pixel, so added to temporary stack
								point2dLabelling current_point = { i_s + 1, j_s, false };
								tempStack.push_back(current_point);
							}
							else {
								statusArray[x_new(i_s + 1, j_s, length)] = true;
							}
						}

						//North
						if (j_s > 0 && statusArray[x_new(i_s, j_s - 1, length)] == false) {
							if (imageDataPtr[x_new(i_s, j_s - 1, length)] == object) {
								//object pixel, so added to temporary stack
								point2dLabelling current_point = { i_s, j_s - 1, false };
								tempStack.push_back(current_point);
							}
							else {
								statusArray[x_new(i_s, j_s - 1, length)] = true;
							}
						}

						//updating of processed element before removal
						labelArray[xd_s] = label;
						statusArray[xd_s] = true;

						//Remove the processed element from the initial stack
						initialStack.pop_back();

						//Add new found neighbors in initial stack
						if (tempStack.size() != 0) {
							//if no new element is found tempStack is empty 
							// and nothing will happend
							for (k = 0; k < tempStack.size(); k++) {
								initialStack.push_back(tempStack[k]);
							}
						}

						//empty the temporary stack
						while (tempStack.size() > 0) {
							tempStack.pop_back();
						}
					}

					label++;
				}
				else {
					statusArray[xd] = true;
				}
			}
		}
	}
	return true;
}

//for 3D images
bool labelling3D(dataType** imageDataPtr, int** segmentedImage, bool** statusArray, const size_t length, const size_t width, const size_t height, dataType object) {

	if (imageDataPtr == NULL || segmentedImage == NULL || statusArray == NULL)
		return false;

	vector<size_t> iStack, jStack, kStack, iTmpStack, jTmpStack, kTmpStack;
	size_t i = 0, j = 0, k = 0, iNew = 0, jNew = 0, kNew = 0, n = 0;
	int label = 1;

	//statusArray has false everywhere in the beginning
	//false---> non-processed and true---> processed
	//segmentedImage has 0 everywhere
	for (k = 0; k < height; k++) {
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {

				if (statusArray[k][x_new(i, j, length)] == false) {
					if (imageDataPtr[k][x_new(i, j, length)] == object) {

						//top neighbor
						if (k > 0 && statusArray[k - 1][x_new(i, j, length)] == false) {
							if (imageDataPtr[k - 1][x_new(i, j, length)] == object) {
								//the element is in current region, then add its coordinates in stacks
								kStack.push_back(k - 1);
								iStack.push_back(i);
								jStack.push_back(j);
							}
							else {
								//it's not object, so, update its status
								statusArray[k - 1][x_new(i, j, length)] = true;
							}
						}
						//down neighbor
						if (k < height - 1 && statusArray[k + 1][x_new(i, j, length)] == false) {
							if (imageDataPtr[k + 1][x_new(i, j, length)] == object) {
								//the element is in current region, then add its coordinates in stacks
								kStack.push_back(k + 1);
								iStack.push_back(i);
								jStack.push_back(j);
							}
							else {
								//it's not object, so, update its status
								statusArray[k + 1][x_new(i, j, length)] = true;
							}
						}
						//left neighbor
						if (i > 0 && statusArray[k][x_new(i - 1, j, length)] == false) {
							if (imageDataPtr[k][x_new(i - 1, j, length)] == object) {
								//the element is in region, then add its coordinates in stacks
								kStack.push_back(k);
								iStack.push_back(i - 1);
								jStack.push_back(j);
							}
							else {
								statusArray[k][x_new(i - 1, j, length)] = true;
							}
						}
						//right neighbor
						if (i < length - 1 && statusArray[k][x_new(i + 1, j, length)] == false) {
							if (imageDataPtr[k][x_new(i + 1, j, length)] == object) {
								//if element is in region add its coordinates in stacks
								kStack.push_back(k);
								iStack.push_back(i + 1);
								jStack.push_back(j);
							}
							else {
								statusArray[k][x_new(i + 1, j, length)] = true;
							}
						}
						//front neighbor
						if (j > 0 && statusArray[k][x_new(i, j - 1, length)] == false) {
							if (imageDataPtr[k][x_new(i, j - 1, length)] == object) {
								//if element is in region add its coodinates in stacks
								kStack.push_back(k);
								iStack.push_back(i);
								jStack.push_back(j - 1);
							}
							else {
								statusArray[k][x_new(i, j - 1, length)] = true;
							}
						}
						//behind neighbor
						if (j < width - 1 && statusArray[k][x_new(i, j + 1, length)] == false) {
							if (imageDataPtr[k][x_new(i, j + 1, length)] == object) {
								//the element is in region, then add its coodinates in stacks
								kStack.push_back(k);
								iStack.push_back(i);
								jStack.push_back(j + 1);
							}
							else {
								statusArray[k][x_new(i, j + 1, length)] = true;
							}
						}
						
						//after checking all neighbors of current element
						//we give its label, and update its status
						segmentedImage[k][x_new(i, j, length)] = label;
						statusArray[k][x_new(i, j, length)] = true;

						//Now we check neighbors of its neighbors saved in the stacks
						//I start by the last added element in stacks
						//If there is no neighbor iStack.size() = jStack.size() = kStack.size() = 0 , and the while loop will no be ran
						while (iStack.size() > 0 && jStack.size() > 0 && kStack.size() > 0) {
							//One is enought because they have same size

							//We work with the last element in the initial stacks
							kNew = kStack.size() - 1;
							iNew = iStack.size() - 1;
							jNew = jStack.size() - 1;

							//top neighbor
							if (kStack[kNew] > 0 && statusArray[kStack[kNew] - 1][x_new(iStack[iNew], jStack[jNew], length)] == false) {
								if (imageDataPtr[kStack[kNew] - 1][x_new(iStack[iNew], jStack[jNew], length)] == object) {
									//If the element is in the current region, then save its coordinates in temporary stacks
									iTmpStack.push_back(iStack[iNew]);
									jTmpStack.push_back(jStack[jNew]);
									kTmpStack.push_back(kStack[kNew] - 1);
								}
								else {
									//Not in region, update status and go to the next neighbor
									statusArray[kStack[kNew] - 1][x_new(iStack[iNew], jStack[jNew], length)] = true;
								}
							}
							//down neighbor
							if (kStack[kNew] < height - 1 && statusArray[kStack[kNew] + 1][x_new(iStack[iNew], jStack[jNew], length)] == false) {
								if (imageDataPtr[kStack[kNew] + 1][x_new(iStack[iNew], jStack[jNew], length)] == object) {
									//the element is in the current region, then save its coordinates in temporary stacks
									iTmpStack.push_back(iStack[iNew]);
									jTmpStack.push_back(jStack[jNew]);
									kTmpStack.push_back(kStack[kNew] + 1);
								}
								else {
									//Not in region, update status and go to the next neighbor
									statusArray[kStack[kNew] + 1][x_new(iStack[iNew], jStack[jNew], length)] = true;
								}
							}
							//right neighbor
							if (iStack[iNew] > 0 && statusArray[kStack[kNew]][x_new(iStack[iNew] - 1, jStack[jNew], length)] == false) {
								if (imageDataPtr[kStack[kNew]][x_new(iStack[iNew] - 1, jStack[jNew], length)] == object) {
									//the element is in the current region, then save its coordinates in temporary stacks
									iTmpStack.push_back(iStack[iNew] - 1);
									jTmpStack.push_back(jStack[jNew]);
									kTmpStack.push_back(kStack[kNew]);
								}
								else {
									//Not in region, update status and go to the next neighbor
									statusArray[kStack[kNew]][x_new(iStack[iNew] - 1, jStack[jNew], length)] = true;
								}
							}
							//left neighbor
							if (iStack[iNew] < length - 1 && statusArray[kStack[kNew]][x_new(iStack[iNew] + 1, jStack[jNew], length)] == false) {
								if (imageDataPtr[kStack[kNew]][x_new(iStack[iNew] + 1, jStack[jNew], length)] == object) {
									//If the element is in the current region, then save its coordinates in temporary stacks
									iTmpStack.push_back(iStack[iNew] + 1);
									jTmpStack.push_back(jStack[jNew]);
									kTmpStack.push_back(kStack[kNew]);
								}
								else {
									//Not in region, update status and go to the next neighbor
									statusArray[kStack[kNew]][x_new(iStack[iNew] + 1, jStack[jNew], length)] = true;
								}
							}
							//front neighbor
							if (jStack[jNew] > 0 && statusArray[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] - 1, length)] == false) {
								if (imageDataPtr[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] - 1, length)] == object) {
									//If the element is in the current region, then save its coordinates in temporary stacks
									iTmpStack.push_back(iStack[iNew]);
									jTmpStack.push_back(jStack[jNew] - 1);
									kTmpStack.push_back(kStack[kNew]);
								}
								else {
									//Not in region, update status and go to the next neighbor
									statusArray[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] - 1, length)] = true;
								}
							}
							//behind neighbor
							if (jStack[jNew] < width - 1 && statusArray[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] + 1, length)] == false) {
								if (imageDataPtr[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] + 1, length)] == object) {
									//the element is in the current region, then save its coordinates in temporary stacks
									iTmpStack.push_back(iStack[iNew]);
									jTmpStack.push_back(jStack[jNew] + 1);
									kTmpStack.push_back(kStack[kNew]);
								}
								else {
									//Not in region, update status and go to the next neighbor
									statusArray[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] + 1, length)] = true;
								}
							}


							//updating of processed element befor removal
							segmentedImage[kStack[kNew]][x_new(iStack[iNew], jStack[jNew], length)] = label;
							statusArray[kStack[kNew]][x_new(iStack[iNew], jStack[jNew], length)] = true;

							//Remove the processed element of the initial stacks
							kStack.pop_back();
							iStack.pop_back();
							jStack.pop_back();

							//Add new found neighbors in initial stacks
							//we can do it once because they have the same size
							for (n = 0; n < iTmpStack.size(); n++) {
								//if any neighbors have been found iTmpStack.size() = 0
								//and nothing will happen
								iStack.push_back(iTmpStack[n]);
							}
							for (n = 0; n < jTmpStack.size(); n++) {
								//if any neighbors have been found jTmpStack.size() = 0
								//and nothing will happen
								jStack.push_back(jTmpStack[n]);
							}
							for (n = 0; n < kTmpStack.size(); n++) {
								//if any neighbors have been found jTmpStack.size() = 0
								//and nothing will happen
								kStack.push_back(kTmpStack[n]);
							}

							//empty the temporary stacks
							//we can do it once because they have the same size
							while (iTmpStack.size() > 0) {
								iTmpStack.pop_back();
							}
							while (jTmpStack.size() > 0) {
								jTmpStack.pop_back();
							}
							while (kTmpStack.size() > 0) {
								kTmpStack.pop_back();
							}

							//End of big while loop
						}
						label++;
					}
					else {
						statusArray[k][x_new(i, j, length)] = true;
					}
				}
			}
		}
	}

	return true;
}

int countNumberOfRegionsCells(dataType** binaryImagePtr, const size_t length, const size_t width, const size_t height, dataType object) {
	if (binaryImagePtr == NULL)
		return 0;
	int numberOfRegionCells = 0;
	for (size_t k = 0; k < height; k++) {
		for (size_t i = 0; i < length * width; i++) {
			if (binaryImagePtr[k][i] == object) {
				numberOfRegionCells++;
			}
		}
	}
	return numberOfRegionCells;
}

//===============================================================================================================

bool regionGrowing(dataType** maskThreshold, dataType** segmentedImage, const size_t length, const size_t width, const size_t height, Point3D pSeed) {
	
	if (maskThreshold == NULL || segmentedImage == NULL) {
		return false;
	}

	size_t dim2D = length * width;

	//array to verify given pixel status, wether processed or not yet processed
	bool** pStatus = new bool* [height];
	for (size_t s = 0; s < height; s++) 
	{
		pStatus[s] = new bool[dim2D] {false};
		if (pStatus[s] == NULL)
			return false;
	}
	if (pStatus == NULL)
		return false;

	vector<point3dLabelling> pSegment, tempSegment;

	size_t height_minus = height - 1;
	size_t length_minus = length - 1;
	size_t width_minus = width - 1;

	size_t iSeed, jSeed, kSeed;
	if ((pSeed.x >= 0 && pSeed.x <= length_minus) && (pSeed.y >= 0 && pSeed.y <= width_minus) && (pSeed.z >= 0 && pSeed.z <= height_minus)) {
		iSeed = (size_t)pSeed.x;
		jSeed = (size_t)pSeed.y;
		kSeed = (size_t)pSeed.z;
	}
	else {
		cout << "Seed out of the bounds" << endl;
		return false;
	}

	//West neighbor
	if (iSeed > 0) {
		size_t i_minus = iSeed - 1;
		if (maskThreshold[kSeed][x_new(i_minus, jSeed, length)] == 1.0)
		{
			point3dLabelling pWest = { i_minus, jSeed, kSeed, 1 };
			pSegment.push_back(pWest);
		}
		else {
			pStatus[kSeed][x_new(i_minus, jSeed, length)] = true;
		}
	}
	//East neighnor
	if (iSeed < length_minus) {
		size_t i_plus = iSeed + 1;
		if (maskThreshold[kSeed][x_new(i_plus, jSeed, length)] == 1.0)
		{
			point3dLabelling pEast = { i_plus, jSeed, kSeed, 1 };
			pSegment.push_back(pEast);
		}
		else {
			pStatus[kSeed][x_new(i_plus, jSeed, length)] = true;
		}
	}
	//North
	if (jSeed > 0) {
		size_t j_minus = jSeed - 1;
		if (maskThreshold[kSeed][x_new(iSeed, j_minus, length)] == 1.0)
		{
			point3dLabelling pNorth = { iSeed, j_minus, kSeed, 1 };
			pSegment.push_back(pNorth);
		}
		else {
			pStatus[kSeed][x_new(iSeed, j_minus, length)] = true;
		}
	}
	//South
	if (jSeed < width_minus) {
		size_t j_plus = jSeed + 1;
		if (maskThreshold[kSeed][x_new(iSeed, j_plus, length)] == 1.0)
		{
			point3dLabelling pSouth = { iSeed, j_plus, kSeed, 1 };
			pSegment.push_back(pSouth);
		}
		else {
			pStatus[kSeed][x_new(iSeed, j_plus, length)] = true;
		}
	}
	//Top
	if (kSeed > 0) {
		size_t k_minus = kSeed - 1;
		if (maskThreshold[kSeed][x_new(iSeed, jSeed, length)] == 1.0)
		{
			point3dLabelling pTop = { iSeed, jSeed, k_minus, 1 };
			pSegment.push_back(pTop);
		}
		else {
			pStatus[k_minus][x_new(iSeed, jSeed, length)] = true;
		}
	}
	//Bottom
	if (kSeed < height_minus) {
		size_t k_plus = kSeed + 1;
		if (maskThreshold[kSeed][x_new(iSeed, jSeed, length)] == 1.0)
		{
			point3dLabelling pBottom = { iSeed, jSeed, k_plus, 1 };
			pSegment.push_back(pBottom);
		}
		else {
			pStatus[k_plus][x_new(iSeed, jSeed, length)] = true;
		}
	}

	size_t xd = x_new(iSeed, jSeed, length);
	pStatus[kSeed][xd] = true;
	if (maskThreshold[kSeed][xd] == 1) {
		segmentedImage[kSeed][xd] = 1.0;
	}

	if (pSegment.size() == 0) {
		cout << "Seed wrongly chosen" << endl;
		return false;
	}
	else {
		while(pSegment.size() > 0)
		{
			int l = pSegment.size();
			iSeed = (size_t)pSegment[l - 1].x;
			jSeed = (size_t)pSegment[l - 1].y;
			kSeed = (size_t)pSegment[l - 1].z;
			xd = x_new(iSeed, jSeed, length);

			if (pStatus[kSeed][xd] == false) {
				
				//West neighbor
				if (iSeed > 0) {
					size_t i_minus = iSeed - 1;
					if (maskThreshold[kSeed][x_new(i_minus, jSeed, length)] == 1.0)
					{
						if (pStatus[kSeed][x_new(i_minus, jSeed, length)] == false) {
							point3dLabelling pWest = { i_minus, jSeed, kSeed, 1 };
							tempSegment.push_back(pWest);
						}
					}
					else
					{
						pStatus[kSeed][x_new(i_minus, jSeed, length)] = true;
					}
				}
				//East neighbor
				if (iSeed < length_minus) {
					size_t i_plus = iSeed + 1;
					if (maskThreshold[kSeed][x_new(i_plus, jSeed, length)] == 1.0)
					{
						if (pStatus[kSeed][x_new(i_plus, jSeed, length)] == false) {
							point3dLabelling pEast = { i_plus, jSeed, kSeed, 1 };
							tempSegment.push_back(pEast);
						}
					}
					else
					{
						pStatus[kSeed][x_new(i_plus, jSeed, length)] = true;
					}
				}
				//North neighbor
				if (jSeed > 0) {
					size_t j_minus = jSeed - 1;
					if (maskThreshold[kSeed][x_new(iSeed, j_minus, length)] == 1.0)
					{
						if (pStatus[kSeed][x_new(iSeed, j_minus, length)] == false) {
							point3dLabelling pNorth = { iSeed, j_minus, kSeed, 1 };
							tempSegment.push_back(pNorth);
						}
					}
					else
					{
						pStatus[kSeed][x_new(iSeed, j_minus, length)] = true;
					}
				}
				//South neighbor
				if (jSeed < width_minus) {
					size_t j_plus = jSeed + 1;
					if (maskThreshold[kSeed][x_new(iSeed, j_plus, length)] == 1.0)
					{
						if (pStatus[kSeed][x_new(iSeed, j_plus, length)] == false) {
							point3dLabelling pSouth = { iSeed, j_plus, kSeed, 1 };
							tempSegment.push_back(pSouth);
						}
					}
					else
					{
						pStatus[kSeed][x_new(iSeed, j_plus, length)] = true;
					}
				}
				//Top neighbor
				if (kSeed > 0) {
					size_t k_minus = kSeed - 1;
					if (maskThreshold[k_minus][x_new(iSeed, jSeed, length)] == 1.0)
					{
						if (pStatus[k_minus][x_new(iSeed, jSeed, length)] == false) {
							point3dLabelling pTop = { iSeed, jSeed, k_minus, 1 };
							tempSegment.push_back(pTop);
						}
					}
					else
					{
						pStatus[k_minus][x_new(iSeed, jSeed, length)] = true;
					}
				}
				//Bottom neighbor
				if (kSeed < height_minus) {
					size_t k_plus = kSeed + 1;
					if (maskThreshold[k_plus][x_new(iSeed, jSeed, length)] == 1.0)
					{
						if (pStatus[k_plus][x_new(iSeed, jSeed, length)] == false) {
							point3dLabelling pBottom = { iSeed, jSeed, k_plus, 1 };
							tempSegment.push_back(pBottom);
						}
					}
					else
					{
						pStatus[k_plus][x_new(iSeed, jSeed, length)] = true;
					}
				}

				pStatus[kSeed][xd] = true;
				segmentedImage[kSeed][xd] = 1.0;
				pSegment.pop_back();

				//Move the elements from the temporary segment to the initial segment
				for (size_t s = 0; s < tempSegment.size(); s++)
				{
					pSegment.push_back(tempSegment[s]);
				}

				//Empty the temporary segment
				while (tempSegment.size() > 0)
				{
					tempSegment.pop_back();
				}

			}
			else{
				pSegment.pop_back();
			}
		}
	}

	for (size_t s = 0; s < height; s++) {
		delete[] pStatus[s];
	}
	delete[] pStatus;

	return true;
}

//===============================================================================================================

/*
bool lungsSegmentation(Image_Data imageData, dataType** segmentedImage, Point3D pSeed, dataType tMin, dataType tMax) {
	
	if (imageData.imageDataPtr == NULL || segmentedImage == NULL) {
		return false;
	}

	const size_t length = imageData.length;
	const size_t width = imageData.width;
	const size_t height = imageData.height;
	size_t dim2D = length * width;

	bool** pStatus = new bool* [height];
	dataType** maskThreshold = new dataType * [height];
	for (size_t s = 0; s < height; s++)
	{
		pStatus[s] = new bool[dim2D] {false};
		maskThreshold[s] = new dataType[dim2D]{ 0 };
		if (pStatus[s] == NULL || maskThreshold[s] == NULL)
			return false;
	}
	if (pStatus == NULL || maskThreshold == NULL)
		return false;

	copyDataToAnotherArray(imageData.imageDataPtr, maskThreshold, height, length, width);
	thresholding3dFunctionN(maskThreshold, length, width, height, tMin, tMax, 0.0, 1.0);
	dilatation3D(maskThreshold, length, width, height, 0.0, 1.0);

	regionGrowing(maskThreshold, segmentedImage, length, width, height, pSeed);

	//delete largest region
	for (size_t k = 0; k < height; k++) {
		for (size_t i = 0; i < length; i++) {
			for (size_t j = 0; j < width; j++) {
				size_t xd = x_new(i, j, length);
				if (segmentedImage[k][xd] == 1.0) {
					maskThreshold[k][xd] = 0.0;

					//initial the array for next step
					segmentedImage[k][xd] = 0.0;
				}
			}
		}
	}

	//find new seed
	size_t count_forground = 0;
	Point3D nSeed = { 0, 0, 0 };
	for (size_t k = 0; k < height; k++) {
		for (size_t i = 0; i < length; i++) {
			for (size_t j = 0; j < width; j++) {
				size_t xd = x_new(i, j, length);
				if (maskThreshold[k][xd] == 1.0) {
					count_forground++;
					if (count_forground == 1) {
						nSeed.x = i;
						nSeed.y = j;
						nSeed.z = k;
					}
				}
			}
		}
	}
	cout << count_forground << " foreground points have been tested = " << endl;
	cout << "Exit with k = " << nSeed.z << ", i = " << nSeed.x << ", j = " << nSeed.y << endl;

	regionGrowing(maskThreshold, segmentedImage, length, width, height, nSeed);

	for (size_t s = 0; s < height; s++) {
		delete[] pStatus[s];
		delete[] maskThreshold[s];
	}
	delete[] pStatus;
	delete[] maskThreshold;

	return true;
}
*/

