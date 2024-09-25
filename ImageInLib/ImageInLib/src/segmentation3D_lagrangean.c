#include "segmentation3D_lagrangean.h"
#include "common_functions.h"
#include <stdlib.h>
#include <math.h>
#include "file.h"

bool lagrangeanExplicit3DCurveSegmentation(Image_Data inputImage3D, const Lagrangean3DSegmentationParameters* pSegmentationParams,
    unsigned char* pOutputPathPtr, Curve3D* pResultSegmentation)
{

    if (pSegmentationParams == NULL || pResultSegmentation == NULL) {
        return false;
    }

    CurvePoint3D* pOldCurve = (CurvePoint3D*)malloc(sizeof(CurvePoint3D) * pResultSegmentation->numPoints);
    Curve3D oldSegmentation = { pOldCurve, pResultSegmentation->numPoints };

    const size_t sizeHeight = sizeof(dataType*) * inputImage3D.height;
    const size_t sizeLength = sizeof(dataType*) * inputImage3D.length;
    const size_t sizeWidth = sizeof(dataType*) * inputImage3D.width;

    const size_t sliceDimension = sizeLength * sizeWidth;
    const size_t sliceSize = sliceDimension * sizeof(dataType);

    dataType** pgrad_x = (dataType**)malloc(sizeHeight);         // velocity component x
    dataType** pgrad_y = (dataType**)malloc(sizeHeight);        // velocity component y
    dataType** pgrad_z = (dataType**)malloc(sizeHeight);       // velocity component z
    dataType** abs_val_grad = (dataType**)malloc(sizeHeight); // absolute value of gradient
    for (size_t k = 0; k < inputImage3D.height; k++) {
        pgrad_x[k] = (dataType*)malloc(sliceSize);
        pgrad_y[k] = (dataType*)malloc(sliceSize);
        pgrad_z[k] = (dataType*)malloc(sliceSize);
        abs_val_grad[k] = (dataType*)malloc(sliceSize);
        if (pgrad_x[k] == NULL || pgrad_y[k] == NULL || pgrad_z[k] == NULL || abs_val_grad[k] == NULL)
            return false;
    }
    if (pgrad_x == NULL || pgrad_y == NULL || pgrad_z == NULL || abs_val_grad == NULL)
        return false;

    dataType* edge_detector = (dataType*)malloc(sizeHeight); // edge detector
    const dataType edge_detector_coef = 10;
    dataType* similar_intensity_detector = NULL;
    //if (!pSegmentationParams->open_curve)
    //{
    //    similar_intensity_detector = (dataType*)malloc(dataSize); // similar intensity detector
    //}

    const dataType similar_intensity_detector_coef = 10;
    const dataType hx = 1.0, hy = 1.0, hz = 1.0;      //spatial discretization step
    const dataType hx_c = 1.0, hy_c = 1.0, hz_c = 1.0;  //h for central differences
    Point3D current_grad;
    FiniteVolumeSize3D finite_volume = { 1.0, 1.0 , 1.0 };

    dataType lambda = 1.0;

    if (!pSegmentationParams->open_curve) {
        lambda = pSegmentationParams->lambda;
    }

    //get absolute value of gradient
    size_t xd = 0;
    for (size_t k = 0; k < inputImage3D.height; k++)
    {
        for (size_t i = 0; i < inputImage3D.length; i++)
        {
            for (size_t j = 0; j < inputImage3D.width; j++)
            {
                getGradient3D(inputImage3D.imageDataPtr, inputImage3D.length, inputImage3D.width, inputImage3D.height, i, j, k, finite_volume, &current_grad);
                xd = x_new(i, j, inputImage3D.length);
                abs_val_grad[k][xd] = norm3D(current_grad);
            }
        }
    }

    return true;
}