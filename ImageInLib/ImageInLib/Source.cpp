#pragma warning(disable : 6011)
#pragma warning(disable : 4244)
#pragma warning(disable : 4477)
#pragma warning(disable : 4267)
#pragma warning(disable : 4700)

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <string.h> 
#include <math.h>
#include <time.h>

#include "file.h"
#include "common_functions.h"
#include "Labelling.h"
#include "template_functions.h"


int main() {

	size_t i, j, k;

	//image Dimensions
	size_t Width = 144;
	size_t Length = 144;
	const size_t dim2D = Width * Length;
	const size_t Height = 255;

	std::string inputPath = "C:/Users/Konan Allaly/Documents/Tests/input/";
	std::string outputPath = "C:/Users/Konan Allaly/Documents/Tests/output/";

	dataType** imageData = new dataType* [Height];
	for (k = 0; k < Height; k++) {
		imageData[k] = new dataType[dim2D];
	}
	if (imageData == NULL)
		return false;

	for (k = 0; k < Height; k++) {
		for (i = 0; i < dim2D; i++) {
			imageData[k][i] = 0.0;
		}
	}

	std::string inputImagePath = inputPath + "patient2_pet.raw";

	Operation operation = LOAD_DATA;
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, inputImagePath.c_str(), operation, false);

	std::string  outputImagePath = outputPath + "load_PET.raw";

	operation = STORE_DATA;
	manageRAWFile3D<dataType>(imageData, Length, Width, Height, outputImagePath.c_str(), operation, false);
	
	for (k = 0; k < Height; k++) {
		delete[] imageData[k];
	}
	delete[] imageData;

	return EXIT_SUCCESS;
}
