/*
* Author: Konan ALLALY
* Purpose: INFLANET project - Image Processing in Nuclear Medicine (2D/3D)
* Language:  C/C++
*/
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

using namespace std;

//for 2D image using 2D arrays
bool labelling2D(int** imageDataPtr, int** segmentedImage, bool** statusArray, const size_t length, const size_t width, int object) {

	if (imageDataPtr == NULL)
		return false;
	if (segmentedImage == NULL)
		return false;
	if (statusArray == NULL)
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

//for 3D images
bool labelling3D(dataType** imageDataPtr, int** segmentedImage, bool** statusArray, const size_t length, const size_t width, const size_t height, dataType object) {

	if (imageDataPtr == NULL)
		return false;
	if (segmentedImage == NULL)
		return false;
	if (statusArray == NULL)
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

bool regionGrowing(dataType** imageDataPtr, dataType** segmentedImage, bool** statusArray, const size_t length, const size_t width, const size_t height, dataType thres_min, dataType thres_max, Point3D * seedPoint) {

	if (imageDataPtr == NULL || segmentedImage == NULL || statusArray == NULL)
		return false;

	vector<size_t> iStack, jStack, kStack, iTmpStack, jTmpStack, kTmpStack;
	size_t i, j, k, iNew = 0, jNew = 0, kNew = 0, n = 0, dim2D = length * width;

	i = (size_t)seedPoint->x; j = (size_t)seedPoint->y; k = (size_t)seedPoint->z;

	//Processed the seed point
	//top neighbor
	if (k > 0 && statusArray[k - 1][x_new(i, j, length)] == false) {
		if (imageDataPtr[k - 1][x_new(i, j, length)] >= thres_min && imageDataPtr[k - 1][x_new(i, j, length)] <= thres_max) {
			//the element is in region, then add its coordinates in stacks
			kStack.push_back(k - 1);
			iStack.push_back(i);
			jStack.push_back(j);
		}
		else {
			statusArray[k - 1][x_new(i, j, length)] = true;
			segmentedImage[k - 1][x_new(i, j, length)] = 0;
		}
	}
	//bottom neighbor
	if (k < height - 1 && statusArray[k + 1][x_new(i, j, length)] == false) {
		if (imageDataPtr[k + 1][x_new(i, j, length)] >= thres_min && imageDataPtr[k + 1][x_new(i, j, length)] <= thres_max) {
			//the element is in region, then add its coordinates in stacks
			kStack.push_back(k + 1);
			iStack.push_back(i);
			jStack.push_back(j);
		}
		else {
			statusArray[k + 1][x_new(i, j, length)] = true;
			segmentedImage[k + 1][x_new(i, j, length)] = 0;
		}
	}
	//left neighbor
	if (i > 0 && statusArray[k][x_new(i - 1, j, length)] == false) {
		if (imageDataPtr[k][x_new(i - 1, j, length)] >= thres_min && imageDataPtr[k][x_new(i - 1, j, length)] <= thres_max) {
			//the element is in region, then add its coordinates in stacks
			kStack.push_back(k);
			iStack.push_back(i - 1);
			jStack.push_back(j);
		}
		else {
			statusArray[k][x_new(i - 1, j, length)] = true;
			segmentedImage[k][x_new(i - 1, j, length)] = 0;
		}
	}
	//right neighbor
	if (i < length - 1 && statusArray[k][x_new(i + 1, j, length)] == false) {
		if (imageDataPtr[k][x_new(i + 1, j, length)] >= thres_min && imageDataPtr[k][x_new(i + 1, j, length)] <= thres_max) {
			//the element is in region, then add its coordinates in stacks
			kStack.push_back(k);
			iStack.push_back(i + 1);
			jStack.push_back(j);
		}
		else {
			statusArray[k][x_new(i + 1, j, length)] = true;
			segmentedImage[k][x_new(i + 1, j, length)] = 0;
		}
	}
	//front neighbor
	if (j > 0 && statusArray[k][x_new(i, j - 1, length)] == false) {
		if (imageDataPtr[k][x_new(i, j - 1, length)] >= thres_min && imageDataPtr[k][x_new(i, j - 1, length)] <= thres_max) {
			//the element is in region, then add its coordinates in stacks
			kStack.push_back(k);
			iStack.push_back(i);
			jStack.push_back(j - 1);
		}
		else {
			statusArray[k][x_new(i, j - 1, length)] = true;
			segmentedImage[k][x_new(i, j - 1, length)] = 0;
		}
	}
	//behind neighbor
	if (j < width - 1 && statusArray[k][x_new(i, j + 1, length)] == false) {
		if (imageDataPtr[k][x_new(i, j + 1, length)] >= thres_min && imageDataPtr[k][x_new(i, j + 1, length)] <= thres_max) {
			//the element is in region, then add its coordinates in stacks
			kStack.push_back(k);
			iStack.push_back(i);
			jStack.push_back(j + 1);
		}
		else {
			statusArray[k][x_new(i, j + 1, length)] = true;
			segmentedImage[k][x_new(i, j + 1, length)] = 0;
		}
	}
	//Update the seed point status
	statusArray[k][x_new(i, j, length)] = true;
	//set  the seed point's label
	segmentedImage[k][x_new(i, j, length)] = imageDataPtr[k][x_new(i, j, length)];

	//Processed the other points
	while (iStack.size() > 0 && jStack.size() > 0 && kStack.size() > 0) {
		//One is enought because they have same size

		//We work with the last element in the initial stacks
		kNew = kStack.size() - 1;
		iNew = iStack.size() - 1;
		jNew = jStack.size() - 1;

		//top neighbor
		if (kStack[kNew] > 0 && statusArray[kStack[kNew] - 1][x_new(iStack[iNew], jStack[jNew], length)] == false) {
			if (imageDataPtr[kStack[kNew] - 1][x_new(iStack[iNew], jStack[jNew], length)] >= thres_min && imageDataPtr[kStack[kNew] - 1][x_new(iStack[iNew], jStack[jNew], length)] <= thres_max) {
				//If the element is in the current region, then save its coordinates in temporary stacks
				iTmpStack.push_back(iStack[iNew]);
				jTmpStack.push_back(jStack[jNew]);
				kTmpStack.push_back(kStack[kNew] - 1);
			}
			else {
				//Not in region, update status and go to the next neighbor
				statusArray[kStack[kNew] - 1][x_new(iStack[iNew], jStack[jNew], length)] = true;
				segmentedImage[kStack[kNew] - 1][x_new(iStack[iNew], jStack[jNew], length)] = 0;
			}
		}
		//down neighbor
		if (kStack[kNew] < height - 1 && statusArray[kStack[kNew] + 1][x_new(iStack[iNew], jStack[jNew], length)] == false) {
			if (imageDataPtr[kStack[kNew] + 1][x_new(iStack[iNew], jStack[jNew], length)] >= thres_min && imageDataPtr[kStack[kNew] + 1][x_new(iStack[iNew], jStack[jNew], length)] <= thres_max) {
				//the element is in the current region, then save its coordinates in temporary stacks
				iTmpStack.push_back(iStack[iNew]);
				jTmpStack.push_back(jStack[jNew]);
				kTmpStack.push_back(kStack[kNew] + 1);
			}
			else {
				//Not in region, update status and go to the next neighbor
				statusArray[kStack[kNew] + 1][x_new(iStack[iNew], jStack[jNew], length)] = true;
				segmentedImage[kStack[kNew] + 1][x_new(iStack[iNew], jStack[jNew], length)] = 0;
			}
		}
		//right neighbor
		if (iStack[iNew] > 0 && statusArray[kStack[kNew]][x_new(iStack[iNew] - 1, jStack[jNew], length)] == false) {
			if (imageDataPtr[kStack[kNew]][x_new(iStack[iNew] - 1, jStack[jNew], length)] >= thres_min && imageDataPtr[kStack[kNew]][x_new(iStack[iNew] - 1, jStack[jNew], length)] <= thres_max) {
				//the element is in the current region, then save its coordinates in temporary stacks
				iTmpStack.push_back(iStack[iNew] - 1);
				jTmpStack.push_back(jStack[jNew]);
				kTmpStack.push_back(kStack[kNew]);
			}
			else {
				//Not in region, update status and go to the next neighbor
				statusArray[kStack[kNew]][x_new(iStack[iNew] - 1, jStack[jNew], length)] = true;
				segmentedImage[kStack[kNew]][x_new(iStack[iNew] - 1, jStack[jNew], length)] = 0;
			}
		}
		//left neighbor
		if (iStack[iNew] < length - 1 && statusArray[kStack[kNew]][x_new(iStack[iNew] + 1, jStack[jNew], length)] == false) {
			if (imageDataPtr[kStack[kNew]][x_new(iStack[iNew] + 1, jStack[jNew], length)] >= thres_min && imageDataPtr[kStack[kNew]][x_new(iStack[iNew] + 1, jStack[jNew], length)] <= thres_max) {
				//If the element is in the current region, then save its coordinates in temporary stacks
				iTmpStack.push_back(iStack[iNew] + 1);
				jTmpStack.push_back(jStack[jNew]);
				kTmpStack.push_back(kStack[kNew]);
			}
			else {
				//Not in region, update status and go to the next neighbor
				statusArray[kStack[kNew]][x_new(iStack[iNew] + 1, jStack[jNew], length)] = true;
				segmentedImage[kStack[kNew]][x_new(iStack[iNew] + 1, jStack[jNew], length)] = 0;
			}
		}
		//front neighbor
		if (jStack[jNew] > 0 && statusArray[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] - 1, length)] == false) {
			if (imageDataPtr[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] - 1, length)] >= thres_min && imageDataPtr[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] - 1, length)] <= thres_max) {
				//If the element is in the current region, then save its coordinates in temporary stacks
				iTmpStack.push_back(iStack[iNew]);
				jTmpStack.push_back(jStack[jNew] - 1);
				kTmpStack.push_back(kStack[kNew]);
			}
			else {
				//Not in region, update status and give it label
				statusArray[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] - 1, length)] = true;
				segmentedImage[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] - 1, length)] = 0;
			}
		}
		//behind neighbor
		if (jStack[jNew] < width - 1 && statusArray[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] + 1, length)] == false) {
			if (imageDataPtr[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] + 1, length)] >= thres_min && imageDataPtr[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] + 1, length)] <= thres_max) {
				//the element is in the current region, then save its coordinates in temporary stacks
				iTmpStack.push_back(iStack[iNew]);
				jTmpStack.push_back(jStack[jNew] + 1);
				kTmpStack.push_back(kStack[kNew]);
			}
			else {
				//Not in region, update status and go to the next neighbor
				statusArray[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] + 1, length)] = true;
				segmentedImage[kStack[kNew]][x_new(iStack[iNew], jStack[jNew] + 1, length)] = 0;
			}
		}

		//updating of processed element before removal
		statusArray[kStack[kNew]][x_new(iStack[iNew], jStack[jNew], length)] = true;
		segmentedImage[kStack[kNew]][x_new(iStack[iNew], jStack[jNew], length)] = imageDataPtr[kStack[kNew]][x_new(iStack[iNew], jStack[jNew], length)];

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

	return true;
}

bool regionGrowing3D_N(Image_Data ctImageData, dataType** segmentedImage, double radius, Point3D seedPoint) {
	
	if (ctImageData.imageDataPtr == NULL || segmentedImage == NULL)
		return false;

	size_t length = ctImageData.length, width = ctImageData.width, height = ctImageData.height;

	Statistics seedStats = { 0.0, 0.0, 0.0, 0.0 };

	vector<point3dLabelling> initialStack, tempStack;

	size_t i_s = (size_t)seedPoint.x, j_s = (size_t)seedPoint.y, k_s = (size_t)seedPoint.z;
	size_t indCurrent = x_new(i_s, j_s, length);
	size_t i, j, k, xd, n, l, dim2D = length * width;

	dataType** gradientX = new dataType * [height];
	dataType** gradientY = new dataType * [height];
	dataType** gradientZ = new dataType * [height];
	dataType** edgeDetector = new dataType * [height];
	dataType** distanceMap = new dataType * [height];
	bool** statusArray = new bool * [height];
	for (k = 0; k < height; k++) {
		gradientX[k] = new dataType[dim2D]{ 0 };
		gradientY[k] = new dataType[dim2D]{ 0 };
		gradientZ[k] = new dataType[dim2D]{ 0 };
		edgeDetector[k] = new dataType[dim2D]{ 0 };
		distanceMap[k] = new dataType[dim2D]{ 0 };
		statusArray[k] = new bool [dim2D]{ false };
	}
	if (gradientX == NULL || gradientY == NULL || gradientZ == NULL || edgeDetector == NULL || distanceMap == NULL || statusArray == NULL)
		return false;

	//compute input image gradient
	compute3dImageGradient(ctImageData.imageDataPtr, gradientX, gradientY, gradientZ, length, width, height, 1.0);

	dataType norm_of_gradient = 0.0, edge_coef = 1000.0;
	dataType ux = 0.0, uy = 0.0, uz = 0.0;

	//compute edge detector
	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			ux = gradientX[k][i];
			uy = gradientY[k][i];
			uz = gradientZ[k][i];
			norm_of_gradient = sqrt(ux * ux + uy * uy + uz * uz);
			edgeDetector[k][i] = gradientFunction(norm_of_gradient, edge_coef);
			//threshold
			if (edgeDetector[k][i] < 0.15) {
				edgeDetector[k][i] = 0.0;
			}
			else {
				edgeDetector[k][i] = 1.0;
			}
		}
	}

	//string path_name = "C:/Users/Konan Allaly/Documents/Tests/output/";
	//string saving_name = path_name + "edgeDetector.raw";
	//store3dRawData<dataType>(edgeDetector, length, width, height, saving_name.c_str());

	fastSweepingFunction_3D(distanceMap, edgeDetector, length, width, height, 1.0, 10000000.0, 0.0);

	//saving_name = path_name + "distanceMap.raw";
	//store3dRawData<dataType>(distanceMap, length, width, height, saving_name.c_str());

	for (k = 0; k < height; k++) {
		for (i = 0; i < dim2D; i++) {
			//keep points with distance lower than 15 (Aorta maximal radius)
			if (distanceMap[k][i] > 15) {
				distanceMap[k][i] = 0.0;
			}
			//normalization
			distanceMap[k][i] = distanceMap[k][i] / 15;
			//inversion
			distanceMap[k][i] = 1.0 / (1.0 + distanceMap[k][i]);
		}
	}

	////Find threshold range according to the seed point
	//ctImageData.imageDataPtr = distanceMap;
	//Point3D initial_point = getRealCoordFromImageCoord3D(seedPoint, ctImageData.origin, ctImageData.spacing, ctImageData.orientation);
	//seedStats = getStats(ctImageData, initial_point, radius);
	
	dataType thres_min = distanceMap[k_s][indCurrent] - 0.025, thres_max = distanceMap[k_s][indCurrent] + 0.025;
	cout << "thres min = " << thres_min << ", thres max = " << thres_max << endl;

	//find seed point neighbors
	//nord
	if (j_s > 0) {
		size_t indNorth = x_new(i_s, j_s - 1, length);
		if (distanceMap[k_s][indNorth] >= thres_min && distanceMap[k_s][indNorth] <= thres_max) {
			point3dLabelling currentPoint = {i_s, j_s - 1, k_s, false};
			initialStack.push_back(currentPoint);
		}
		else {
			statusArray[k_s][indNorth] = true;
			segmentedImage[k_s][indNorth] = 0.0;
		}
	}

	//south
	if (j_s < width - 1) {
		size_t indSouth = x_new(i_s, j_s + 1, length);
		if (distanceMap[k_s][indSouth] >= thres_min && distanceMap[k_s][indSouth] <= thres_max) {
			point3dLabelling currentPoint = { i_s, j_s + 1, k_s, false };
			initialStack.push_back(currentPoint);
		}
		else {
			statusArray[k_s][indSouth] = true;
			segmentedImage[k_s][indSouth] = 0.0;
		}
	}

	//west
	if (i_s > 0) {
		size_t indWest = x_new(i_s - 1, j_s, length);
		if (distanceMap[k_s][indWest] >= thres_min && distanceMap[k_s][indWest] <= thres_max) {
			point3dLabelling currentPoint = { i_s - 1, j_s, k_s, false };
			initialStack.push_back(currentPoint);
		}
		else {
			statusArray[k_s][indWest] = true;
			segmentedImage[k_s][indWest] = 0.0;
		}
	}

	//east
	if (i_s < length - 1) {
		size_t indEast = x_new(i_s + 1, j_s, length);
		if (distanceMap[k_s][indEast] >= thres_min && distanceMap[k_s][indEast] <= thres_max) {
			point3dLabelling currentPoint = { i_s + 1, j_s, k_s, false };
			initialStack.push_back(currentPoint);
		}
		else {
			statusArray[k_s][indEast] = true;
			segmentedImage[k_s][indEast] = 0.0;
		}
	}

	//top
	if (k_s > 0) {
		if (distanceMap[k_s - 1][indCurrent] >= thres_min && distanceMap[k_s - 1][indCurrent] <= thres_max) {
			point3dLabelling currentPoint = { i_s, j_s, k_s - 1, false };
			initialStack.push_back(currentPoint);
		}
		else {
			statusArray[k_s - 1][indCurrent] = true;
			segmentedImage[k_s - 1][indCurrent] = 0.0;
		}
	}

	//bottom
	if (k_s < height - 1) {
		if (distanceMap[k_s + 1][indCurrent] >= thres_min && distanceMap[k_s + 1][indCurrent] <= thres_max) {
			point3dLabelling currentPoint = { i_s, j_s, k_s + 1, false };
			initialStack.push_back(currentPoint);
		}
		else {
			statusArray[k_s + 1][indCurrent] = true;
			segmentedImage[k_s + 1][indCurrent] = 0.0;
		}
	}

	//we don't set the seed point, so if it doesn't match the theshold condition it will be skipped

	while (initialStack.size() > 0) {
		
		//Processed the last element in the initial stack
		l = initialStack.size() - 1;
		i_s = initialStack[l].x;
		j_s = initialStack[l].y;
		k_s = initialStack[l].z;
		indCurrent = x_new(i_s, j_s, length);

		//north
		if (j_s > 0) {
			size_t indNorth = x_new(i_s, j_s - 1, length);
			if(statusArray[k_s][indNorth] == false){
				if (distanceMap[k_s][indNorth] >= thres_min && distanceMap[k_s][indNorth] <= thres_max) {
					point3dLabelling currentPoint = { i_s, j_s - 1, k_s, false };
					tempStack.push_back(currentPoint);
				}
				else {
					statusArray[k_s][indNorth] = true;
					segmentedImage[k_s][indNorth] = 0.0;
				}
			}
		}

		//south
		if (j_s < width - 1) {
			size_t indSouth = x_new(i_s, j_s + 1, length);
			if (statusArray[k_s][indSouth] == false) {
				if (distanceMap[k_s][indSouth] >= thres_min && distanceMap[k_s][indSouth] <= thres_max) {
					point3dLabelling currentPoint = { i_s, j_s + 1, k_s, false };
					tempStack.push_back(currentPoint);
				}
				else {
					statusArray[k_s][indSouth] = true;
					segmentedImage[k_s][indSouth] = 0.0;
				}
			}
		}

		//east
		if (i_s > 0) {
			size_t indEast = x_new(i_s - 1, j_s, length);
			if (statusArray[k_s][indEast] == false) {
				if (distanceMap[k_s][indEast] >= thres_min && distanceMap[k_s][indEast] <= thres_max) {
					point3dLabelling currentPoint = { i_s - 1, j_s, k_s, false };
					tempStack.push_back(currentPoint);
				}
				else {
					statusArray[k_s][indEast] = true;
					segmentedImage[k_s][indEast] = 0.0;
				}
			}
		}

		//west
		if (i_s < length - 1) {
			size_t indWest = x_new(i_s + 1, j_s, length);
			if (statusArray[k_s][indWest] == false) {
				if (distanceMap[k_s][indWest] >= thres_min && distanceMap[k_s][indWest] <= thres_max) {
					point3dLabelling currentPoint = { i_s + 1, j_s, k_s, false };
					tempStack.push_back(currentPoint);
				}
				else {
					statusArray[k_s][indWest] = true;
					segmentedImage[k_s][indWest] = 0.0;
				}
			}
		}

		//top
		if (k_s > 0) {
			if (statusArray[k_s - 1][indCurrent] == false) {
				if (distanceMap[k_s - 1][indCurrent] >= thres_min && distanceMap[k_s - 1][indCurrent] <= thres_max) {
					point3dLabelling currentPoint = { i_s, j_s, k_s - 1, false };
					tempStack.push_back(currentPoint);
				}
				else {
					statusArray[k_s - 1][indCurrent] = true;
					segmentedImage[k_s - 1][indCurrent] = 0.0;
				}
			}
		}

		//bottom
		if (k_s < height - 1) {
			if (statusArray[k_s + 1][indCurrent] == false) {
				if (distanceMap[k_s + 1][indCurrent] >= thres_min && distanceMap[k_s + 1][indCurrent] <= thres_max) {
					point3dLabelling currentPoint = { i_s, j_s, k_s + 1, false };
					tempStack.push_back(currentPoint);
				}
				else {
					statusArray[k_s + 1][indCurrent] = true;
					segmentedImage[k_s + 1][indCurrent] = 0.0;
				}
			}
		}

		//change the status of the processed point
		statusArray[k_s][indCurrent] = true;
		segmentedImage[k_s][indCurrent] = 1.0;

		//remove the processed point from initial stack
		initialStack.pop_back();

		//move points from temporary stack to initial stack
		for (n = 0; n < tempStack.size(); n++) {
			initialStack.push_back(tempStack[n]);
		}

		//empty the temporary stack
		while (tempStack.size() > 0) {
			tempStack.pop_back();
		}

	}

	for (k = 0; k < height; k++) {
		delete[] gradientX[k];
		delete[] gradientY[k];
		delete[] gradientZ[k];
		delete[] edgeDetector[k];
		delete[] distanceMap[k];
		delete[] statusArray[k];
	}
	delete[] gradientX;
	delete[] gradientY;
	delete[] gradientZ;
	delete[] edgeDetector;
	delete[] distanceMap;
	delete[] statusArray;

	return true;
}
