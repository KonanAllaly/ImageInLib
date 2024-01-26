#include "cropping.h"
#include "imageInterpolation.h"

imageMetaData croppImage3D(Image_Data ctImageData, const size_t offset) {

	imageMetaData croppedImageMetaData;

	croppedImageMetaData.spacing.sx = ctImageData.spacing.sx;
	croppedImageMetaData.spacing.sy = ctImageData.spacing.sy;
	croppedImageMetaData.spacing.sz = ctImageData.spacing.sz;

	size_t height = ctImageData.height;
	size_t length = ctImageData.length;
	size_t width = ctImageData.width;

	size_t k, i, j, x;
	size_t i_min = length, j_min = width, k_min = height, i_max = 0, j_max = 0, k_max = 0;

	for (k = 0; k < height; k++) {
		for (i = 0; i < length; i++) {
			for (j = 0; j < width; j++) {
				x = x_new(i, j, length);
				if (ctImageData.imageDataPtr[k][x] == 1.0) {

					//find i_min
					if (i < i_min) {
						i_min = i;
					}
					//find i_max
					if (i > i_max) {
						i_max = i;
					}

					//find j_min
					if (j < j_min) {
						j_min = j;
					}
					//find j_max
					if (j > j_max) {
						j_max = j;
					}

					//find k_min
					if (k < k_min) {
						k_min = k;
					}
					//find k_max
					if (k > k_max) {
						k_max = k;
					}
				}
			}
		}
	}

	croppedImageMetaData.origin.x = i_min;
	croppedImageMetaData.origin.y = j_min;
	croppedImageMetaData.origin.z = k_min;

	croppedImageMetaData.origin = getRealCoordFromImageCoord3D(croppedImageMetaData.origin, ctImageData.origin, ctImageData.spacing, ctImageData.orientation);

	croppedImageMetaData.dimension[0] = (size_t)(i_max - i_min + 2 * offset);
	croppedImageMetaData.dimension[1] = (size_t)(j_max - j_min + 2 * offset);
	croppedImageMetaData.dimension[2] = (size_t)(k_max - k_min + 2 * offset);

	return croppedImageMetaData;
}