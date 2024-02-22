/*
* Author: Konan ALLALY
* Purpose: INFLANET project - Image Processing in Nuclear Medicine (2D/3D)
* Language:  C/C++
*/
#pragma warning(disable : 6011)
#pragma warning(disable : 4244)
#pragma warning(disable : 4477)
#pragma warning(disable : 4267)

#include<stdio.h>
#include<stdlib.h>
#include <string.h>
#include <common_vtk.h>
#include "heat_equation.h"
#include "../include/morphological_change.h"

//3D functions

//Erosion removes floating pixels and thin lines so that only substantive objects remain
bool erosion3D(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType object, dataType background) {
	if (imageDataPtr == NULL)
		return false;

	size_t i, j, k, cpt, dim2d = xDim * yDim;
	size_t xDim_ext = xDim + 2, yDim_ext = yDim + 2, zDim_ext = zDim + 2, dim2d_ext = xDim_ext * yDim_ext;

	dataType** extendedImage = (dataType**)malloc(zDim_ext * sizeof(dataType*));
	dataType** maskImage = (dataType**)malloc(zDim_ext * sizeof(dataType*));
	for (k = 0; k < zDim_ext; k++) {
		extendedImage[k] = (dataType*)malloc(dim2d_ext * sizeof(dataType));
		maskImage[k] = (dataType*)malloc(dim2d_ext * sizeof(dataType));
	}
	if (extendedImage == NULL || maskImage == NULL) return false;

	//Initialization
	for (k = 0; k < zDim_ext; k++) {
		for (i = 0; i < xDim_ext; i++) {
			for (j = 0; j < yDim_ext; j++) {
				extendedImage[k][x_new(i, j, xDim_ext)] = background;
				maskImage[k][x_new(i, j, xDim_ext)] = background;
			}
		}
	}
	copyDataToExtendedArea(imageDataPtr, extendedImage, zDim, xDim, yDim);
	copyDataToExtendedArea(imageDataPtr, maskImage, zDim, xDim, yDim);
	reflection3D(extendedImage, zDim_ext, xDim_ext, yDim_ext);

	//neigbors check
	for (k = 1; k < zDim_ext - 1; k++) {
		for (i = 1; i < xDim_ext - 1; i++) {
			for (j = 1; j < yDim_ext - 1; j++) {

				cpt = 0;

				if (extendedImage[k - 1][x_new(i - 1, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i - 1, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i - 1, j + 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i, j + 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i + 1, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i + 1, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i + 1, j + 1, xDim_ext)] == object) {
					cpt++;
				}

				if (extendedImage[k][x_new(i - 1, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i - 1, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i - 1, j + 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i, j + 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i + 1, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i + 1, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i + 1, j + 1, xDim_ext)] == object) {
					cpt++;
				}

				if (extendedImage[k + 1][x_new(i - 1, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i - 1, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i - 1, j + 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i, j + 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i + 1, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i + 1, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i + 1, j + 1, xDim_ext)] == object) {
					cpt++;
				}

				if (cpt < 26) {
					maskImage[k][x_new(i, j, xDim_ext)] = background;
				}
			}
		}
	}
	copyDataToReducedArea(imageDataPtr, maskImage, zDim, xDim, yDim);

	//clean memory
	for (k = 0; k < zDim_ext; k++) {
		free(extendedImage[k]); 
		free(maskImage[k]);
	}
	free(extendedImage); free(maskImage);

	return true;
}

bool erosion3dHeighteenNeigbours(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType object, dataType background) {
	if (imageDataPtr == NULL)
		return false;

	size_t i, j, k, cpt, dim2d = xDim * yDim;
	size_t xDim_ext = xDim + 2, yDim_ext = yDim + 2, zDim_ext = zDim + 2, dim2d_ext = xDim_ext * yDim_ext;

	dataType** extendedImage = (dataType**)malloc(zDim_ext * sizeof(dataType*));
	dataType** maskImage = (dataType**)malloc(zDim_ext * sizeof(dataType*));
	for (k = 0; k < zDim_ext; k++) {
		extendedImage[k] = (dataType*)malloc(dim2d_ext * sizeof(dataType));
		maskImage[k] = (dataType*)malloc(dim2d_ext * sizeof(dataType));
	}
	if (extendedImage == NULL || maskImage == NULL) return false;
	//Initialization
	for (k = 0; k < zDim_ext; k++) {
		for (i = 0; i < xDim_ext; i++) {
			for (j = 0; j < yDim_ext; j++) {
				extendedImage[k][x_new(i, j, xDim_ext)] = background;
				maskImage[k][x_new(i, j, xDim_ext)] = background;
			}
		}
	}
	copyDataToExtendedArea(imageDataPtr, extendedImage, zDim, xDim, yDim);
	copyDataToExtendedArea(imageDataPtr, maskImage, zDim, xDim, yDim);
	reflection3D(extendedImage, zDim_ext, xDim_ext, yDim_ext);

	//neigbors check ---> 18 neigbors
	for (k = 1; k < zDim_ext - 1; k++) {
		for (i = 1; i < xDim_ext - 1; i++) {
			for (j = 1; j < yDim_ext - 1; j++) {

				cpt = 0;

				if (extendedImage[k][x_new(i, j, xDim_ext)] != background) {

					if (extendedImage[k - 1][x_new(i - 1, j + 1, xDim_ext)] != background) {
						cpt++;
					}
					if (extendedImage[k - 1][x_new(i, j, xDim_ext)] != background) {
						cpt++;
					}
					if (extendedImage[k - 1][x_new(i, j + 1, xDim_ext)] != background) {
						cpt++;
					}
					if (extendedImage[k - 1][x_new(i + 1, j - 1, xDim_ext)] != background) {
						cpt++;
					}
					if (extendedImage[k - 1][x_new(i + 1, j, xDim_ext)] != background) {
						cpt++;
					}

					if (extendedImage[k][x_new(i - 1, j - 1, xDim_ext)] != background) {
						cpt++;
					}
					if (extendedImage[k][x_new(i - 1, j, xDim_ext)] != background) {
						cpt++;
					}
					if (extendedImage[k][x_new(i - 1, j + 1, xDim_ext)] != background) {
						cpt++;
					}
					if (extendedImage[k][x_new(i, j - 1, xDim_ext)] != background) {
						cpt++;
					}
					if (extendedImage[k][x_new(i, j + 1, xDim_ext)] != background) {
						cpt++;
					}
					if (extendedImage[k][x_new(i + 1, j - 1, xDim_ext)] != background) {
						cpt++;
					}
					if (extendedImage[k][x_new(i + 1, j, xDim_ext)] != background) {
						cpt++;
					}
					if (extendedImage[k][x_new(i + 1, j + 1, xDim_ext)] != background) {
						cpt++;
					}

					if (extendedImage[k + 1][x_new(i - 1, j + 1, xDim_ext)] != background) {
						cpt++;
					}
					if (extendedImage[k + 1][x_new(i, j, xDim_ext)] != background) {
						cpt++;
					}
					if (extendedImage[k + 1][x_new(i, j + 1, xDim_ext)] != background) {
						cpt++;
					}
					if (extendedImage[k + 1][x_new(i + 1, j - 1, xDim_ext)] != background) {
						cpt++;
					}
					if (extendedImage[k + 1][x_new(i + 1, j, xDim_ext)] != background) {
						cpt++;
					}

					if (cpt < 18) {
						maskImage[k][x_new(i, j, xDim_ext)] = background;
					}
				}
			}
		}
	}
	copyDataToReducedArea(imageDataPtr, maskImage, zDim, xDim, yDim);

	for (k = 0; k < zDim_ext; k++) {
		free(extendedImage[k]);
		free(maskImage[k]);
	}
	free(extendedImage);
	free(maskImage);

	return true;
}

//Dilation makes objects more visible and fills in small holes in objects.
bool dilatation3D(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType object, dataType background) {
	if (imageDataPtr == NULL)
		return false;

	size_t i, j, k, cpt, dim2d = xDim * yDim;
	size_t xDim_ext = xDim + 2, yDim_ext = yDim + 2, zDim_ext = zDim + 2, dim2d_ext = xDim_ext * yDim_ext;

	dataType** extendedImage = (dataType**)malloc(zDim_ext * sizeof(dataType*));
	dataType** maskImage = (dataType**)malloc(zDim_ext * sizeof(dataType*));
	for (k = 0; k < zDim_ext; k++) {
		extendedImage[k] = (dataType*)malloc(dim2d_ext * sizeof(dataType));
		maskImage[k] = (dataType*)malloc(dim2d_ext * sizeof(dataType));
	}
	if (extendedImage == NULL || maskImage == NULL) return false;

	//Initialization
	for (k = 0; k < zDim_ext; k++) {
		for (i = 0; i < xDim_ext; i++) {
			for (j = 0; j < yDim_ext; j++) {
				extendedImage[k][x_new(i, j, xDim_ext)] = background;
				maskImage[k][x_new(i, j, xDim_ext)] = background;
			}
		}
	}
	copyDataToExtendedArea(imageDataPtr, extendedImage, zDim, xDim, yDim);
	reflection3D(extendedImage, zDim_ext, xDim_ext, yDim_ext);
	copyDataToAnotherArray(extendedImage, maskImage, zDim_ext, xDim_ext, yDim_ext);

	//neigbors check ---> 26 neighbors
	for (k = 1; k < zDim_ext - 1; k++) {
		for (i = 1; i < xDim_ext - 1; i++) {
			for (j = 1; j < yDim_ext - 1; j++) {

				cpt = 0;

				if (extendedImage[k - 1][x_new(i - 1, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i - 1, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i - 1, j + 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i, j + 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i + 1, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i + 1, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i + 1, j + 1, xDim_ext)] == object) {
					cpt++;
				}

				if (extendedImage[k][x_new(i - 1, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i - 1, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i - 1, j + 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i, j + 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i + 1, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i + 1, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k][x_new(i + 1, j + 1, xDim_ext)] == object) {
					cpt++;
				}

				if (extendedImage[k + 1][x_new(i - 1, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i - 1, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i - 1, j + 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i, j + 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i + 1, j - 1, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i + 1, j, xDim_ext)] == object) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i + 1, j + 1, xDim_ext)] == object) {
					cpt++;
				}

				if (cpt > 0) {
					maskImage[k][x_new(i, j, xDim_ext)] = object;
				}
			}
		}
	}
	copyDataToReducedArea(imageDataPtr, maskImage, zDim, xDim, yDim);

	//clean memory
	for (k = 0; k < zDim_ext; k++) {
		free(extendedImage[k]);
		free(maskImage[k]);
	}
	free(extendedImage); free(maskImage);

	return true;
}

bool dilatation3dHeighteenNeigbours(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType object, dataType background) {
	if (imageDataPtr == NULL)
		return false;

	size_t i, j, k, cpt, dim2d = xDim * yDim;
	size_t xDim_ext = xDim + 2, yDim_ext = yDim + 2, zDim_ext = zDim + 2, dim2d_ext = xDim_ext * yDim_ext;

	dataType** extendedImage = (dataType**)malloc(zDim_ext * sizeof(dataType*));
	dataType** maskImage = (dataType**)malloc(zDim_ext * sizeof(dataType*));
	for (k = 0; k < zDim_ext; k++) {
		extendedImage[k] = (dataType*)malloc(dim2d_ext * sizeof(dataType));
		maskImage[k] = (dataType*)malloc(dim2d_ext * sizeof(dataType));
	}
	if (extendedImage == NULL || maskImage == NULL) return false;
	
	//Initialization
	for (k = 0; k < zDim_ext; k++) {
		for (i = 0; i < dim2d_ext; i++) {
			extendedImage[k][i] = background;
			maskImage[k][i] = background;
		}
	}

	////copy to extended area
	//size_t kn, in, jn;
	//for (k = 0, kn = 1; k < zDim; k++, kn++) {
	//	for (i = 0, in = 1; i < xDim; i++, in++) {
	//		for (j = 0, jn = 1; j < yDim; j++, jn++) {
	//			extendedImage[kn][x_new(in, jn, xDim_ext)] = imageDataPtr[k][x_new(i, j, xDim)];
	//		}
	//	}
	//}
	//reflection3D(extendedImage, zDim_ext, xDim_ext, yDim_ext);
	
	copyDataToExtendedArea(imageDataPtr, extendedImage, zDim, xDim, yDim);
	copyDataToExtendedArea(imageDataPtr, maskImage, zDim, xDim, yDim);
	reflection3D(extendedImage, zDim_ext, xDim_ext, yDim_ext);
	reflection3D(maskImage, zDim_ext, xDim_ext, yDim_ext);

	//neigbors check ---> 18 neigbors
	for (k = 1; k < zDim_ext - 1; k++) {
		for (i = 1; i < xDim_ext - 1; i++) {
			for (j = 1; j < yDim_ext - 1; j++) {

				cpt = 0;

				if (extendedImage[k - 1][x_new(i - 1, j + 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i, j, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i, j + 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i + 1, j - 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i + 1, j, xDim_ext)] != background) {
					cpt++;
				}

				if (extendedImage[k][x_new(i - 1, j - 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k][x_new(i - 1, j, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k][x_new(i - 1, j + 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k][x_new(i, j - 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k][x_new(i, j + 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k][x_new(i + 1, j - 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k][x_new(i + 1, j, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k][x_new(i + 1, j + 1, xDim_ext)] != background) {
					cpt++;
				}

				if (extendedImage[k + 1][x_new(i - 1, j + 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i, j, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i, j + 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i + 1, j - 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i + 1, j, xDim_ext)] != background) {
					cpt++;
				}

				if (cpt > 0) {
					maskImage[k][x_new(i, j, xDim_ext)] = object;
				}

			}
		}
	}
	
	copyDataToReducedArea(imageDataPtr, maskImage, zDim, xDim, yDim);

	////copy to reduced area
	//for (k = 0, kn = 1; k < zDim; k++, kn++) {
	//	for (i = 0, in = 1; i < xDim; i++, in++) {
	//		for (j = 0, jn = 1; j < yDim; j++, jn++) {
	//			imageDataPtr[k][x_new(i, j, xDim)] = maskImage[kn][x_new(in, jn, xDim_ext)];
	//		}
	//	}
	//}
	
	for (k = 0; k < zDim_ext; k++) {
		free(extendedImage[k]);
		free(maskImage[k]);
	}
	free(extendedImage);
	free(maskImage);

	return true;
}

//2D functions
bool erosion2D(dataType* imageDataPtr, const size_t length, const size_t width, dataType object, dataType background) {
	
	if (imageDataPtr == NULL)
		return false;

	const size_t length_ext = length + 2, width_ext = width + 2, dim2D_ext = length_ext * width_ext;
	size_t i, j, count_neighbor = 0;

	dataType* extendedArray = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
	dataType* mask = (dataType*)malloc(sizeof(dataType) * dim2D_ext);

	//initialization
	for (i = 0; i < dim2D_ext; i++) {
		extendedArray[i] = 0.0;
		mask[i] = 0.0;
	}

	copyDataTo2dExtendedArea(imageDataPtr, extendedArray,length, width);
	reflection2D(extendedArray, length_ext, width_ext);
	copyDataTo2dExtendedArea(imageDataPtr, mask, length, width);
	reflection2D(mask, length_ext, width_ext);

	for (i = 1; i < length_ext - 1; i++) {
		for (j = 1; j < width_ext - 1; j++) {
			
			count_neighbor = 0;

			if (extendedArray[x_new(i - 1, j - 1, length_ext)] == object) {
				count_neighbor++;
			}
			if (extendedArray[x_new(i, j - 1, length_ext)] == object) {
				count_neighbor++;
			}
			if (extendedArray[x_new(i + 1, j - 1, length_ext)] == object) {
				count_neighbor++;
			}

			if (extendedArray[x_new(i - 1, j, length_ext)] == object) {
				count_neighbor++;
			}
			if (extendedArray[x_new(i + 1, j, length_ext)] == object) {
				count_neighbor++;
			}

			if (extendedArray[x_new(i - 1, j + 1, length_ext)] == object) {
				count_neighbor++;
			}
			if (extendedArray[x_new(i, j + 1, length_ext)] == object) {
				count_neighbor++;
			}
			if (extendedArray[x_new(i + 1, j + 1, length_ext)] == object) {
				count_neighbor++;
			}

			if (count_neighbor < 8) {
				mask[x_new(i, j, length_ext)] = background;
			}
		}
	}

	copyDataTo2dReducedArea(imageDataPtr, mask, length, width);

	free(extendedArray);
	free(mask);

	return true;
}

bool dilatation2D(dataType* imageDataPtr, const size_t length, const size_t width, dataType object, dataType background) {

	if (imageDataPtr == NULL)
		return false;

	const size_t length_ext = length + 2, width_ext = width + 2, dim2D_ext = length_ext * width_ext;
	size_t i, j, count_neighbor = 0;

	dataType* extendedArray = (dataType*)malloc(sizeof(dataType) * dim2D_ext);
	dataType* mask = (dataType*)malloc(sizeof(dataType) * dim2D_ext);

	//initialization
	for (i = 0; i < dim2D_ext; i++) {
		extendedArray[i] = 0.0;
		mask[i] = 0.0;
	}

	copyDataTo2dExtendedArea(imageDataPtr, extendedArray, length, width);
	reflection2D(extendedArray, length_ext, width_ext);
	copyDataTo2dExtendedArea(imageDataPtr, mask, length, width);
	reflection2D(mask, length_ext, width_ext);

	for (i = 1; i < length_ext - 1; i++) {
		for (j = 1; j < width_ext - 1; j++) {

			count_neighbor = 0;

			if (extendedArray[x_new(i - 1, j - 1, length_ext)] == object) {
				count_neighbor++;
			}
			if (extendedArray[x_new(i, j - 1, length_ext)] == object) {
				count_neighbor++;
			}
			if (extendedArray[x_new(i + 1, j - 1, length_ext)] == object) {
				count_neighbor++;
			}

			if (extendedArray[x_new(i - 1, j, length_ext)] == object) {
				count_neighbor++;
			}
			if (extendedArray[x_new(i + 1, j, length_ext)] == object) {
				count_neighbor++;
			}

			if (extendedArray[x_new(i - 1, j + 1, length_ext)] == object) {
				count_neighbor++;
			}
			if (extendedArray[x_new(i, j + 1, length_ext)] == object) {
				count_neighbor++;
			}
			if (extendedArray[x_new(i + 1, j + 1, length_ext)] == object) {
				count_neighbor++;
			}

			if (count_neighbor >= 1) {
				mask[x_new(i, j, length_ext)] = object;
			}
		}
	}

	copyDataTo2dReducedArea(imageDataPtr, mask, length, width);

	free(extendedArray);
	free(mask);

	return true;
}

//========================================================

bool dilatationTwentySixNeigbours(dataType** imageDataPtr, const size_t xDim, const size_t yDim, const size_t zDim, dataType object, dataType background) {
	if (imageDataPtr == NULL)
		return false;

	size_t i, j, k, cpt, dim2d = xDim * yDim;
	size_t xDim_ext = xDim + 2, yDim_ext = yDim + 2, zDim_ext = zDim + 2, dim2d_ext = xDim_ext * yDim_ext;

	dataType** extendedImage = (dataType**)malloc(zDim_ext * sizeof(dataType*));
	dataType** maskImage = (dataType**)malloc(zDim_ext * sizeof(dataType*));
	for (k = 0; k < zDim_ext; k++) {
		extendedImage[k] = (dataType*)malloc(dim2d_ext * sizeof(dataType));
		maskImage[k] = (dataType*)malloc(dim2d_ext * sizeof(dataType));
	}
	if (extendedImage == NULL || maskImage == NULL) return false;

	//Initialization
	for (k = 0; k < zDim_ext; k++) {
		for (i = 0; i < dim2d_ext; i++) {
			extendedImage[k][i] = background;
			maskImage[k][i] = background;
		}
	}

	copyDataToExtendedArea(imageDataPtr, extendedImage, zDim, xDim, yDim);
	//copyDataToExtendedArea(imageDataPtr, maskImage, zDim, xDim, yDim);
	reflection3D(extendedImage, zDim_ext, xDim_ext, yDim_ext);
	//reflection3D(maskImage, zDim_ext, xDim_ext, yDim_ext);
	//copyDataToAnotherArray(extendedImage, maskImage, zDim_ext, xDim_ext, yDim_ext);

	//neigbors check ---> 18 neigbors
	for (k = 1; k < zDim_ext - 1; k++) {
		for (i = 1; i < xDim_ext - 1; i++) {
			for (j = 1; j < yDim_ext - 1; j++) {

				cpt = 0;

				if (extendedImage[k - 1][x_new(i - 1, j - 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i - 1, j, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i - 1, j + 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i, j - 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i, j, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i, j + 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i + 1, j - 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i + 1, j, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k - 1][x_new(i + 1, j + 1, xDim_ext)] != background) {
					cpt++;
				}

				if (extendedImage[k][x_new(i - 1, j - 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k][x_new(i - 1, j, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k][x_new(i - 1, j + 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k][x_new(i, j - 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k][x_new(i, j + 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k][x_new(i + 1, j - 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k][x_new(i + 1, j, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k][x_new(i + 1, j + 1, xDim_ext)] != background) {
					cpt++;
				}

				if (extendedImage[k + 1][x_new(i - 1, j - 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i - 1, j, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i - 1, j + 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i, j - 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i, j, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i, j + 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i + 1, j - 1, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i + 1, j, xDim_ext)] != background) {
					cpt++;
				}
				if (extendedImage[k + 1][x_new(i + 1, j + 1, xDim_ext)] != background) {
					cpt++;
				}

				if (cpt > 0) {
					maskImage[k][x_new(i, j, xDim_ext)] = object;
				}

			}
		}
	}
	copyDataToReducedArea(imageDataPtr, maskImage, zDim, xDim, yDim);

	for (k = 0; k < zDim_ext; k++) {
		free(extendedImage[k]);
		free(maskImage[k]);
	}
	free(extendedImage);
	free(maskImage);

	return true;
}