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
    const size_t sliceDimension = inputImage3D.length * inputImage3D.width;
    const size_t sliceSize = sliceDimension * sizeof(dataType);

    dataType** pgrad_x = (dataType**)malloc(sizeHeight);       // velocity component x
    dataType** pgrad_y = (dataType**)malloc(sizeHeight);       // velocity component y
    dataType** pgrad_z = (dataType**)malloc(sizeHeight);       // velocity component z
    dataType** abs_val_grad = (dataType**)malloc(sizeHeight);  // absolute value of gradient
    dataType** edge_detector = (dataType**)malloc(sizeHeight); // edge detector
    for (size_t k = 0; k < inputImage3D.height; k++) {
        pgrad_x[k] = (dataType*)malloc(sliceSize);
        pgrad_y[k] = (dataType*)malloc(sliceSize);
        pgrad_z[k] = (dataType*)malloc(sliceSize);
        abs_val_grad[k] = (dataType*)malloc(sliceSize);
        edge_detector[k] = (dataType*)malloc(sliceSize);
        if (pgrad_x[k] == NULL || pgrad_y[k] == NULL || pgrad_z[k] == NULL || 
            abs_val_grad[k] == NULL || edge_detector[k] == NULL)
            return false;
    }
    if (pgrad_x == NULL || pgrad_y == NULL || pgrad_z == NULL || 
        abs_val_grad == NULL || edge_detector == NULL)
        return false;

    dataType** similar_intensity_detector = NULL;
    if (!pSegmentationParams->open_curve)
    {
        similar_intensity_detector = (dataType**)malloc(sizeHeight); // similar intensity detector
        for (size_t k = 0; k < inputImage3D.height; k++) {
            similar_intensity_detector[k] = (dataType*)malloc(sliceSize);
        }
    }

    const dataType hx = 1.0, hy = 1.0, hz = 1.0;                 //spatial discretization step
    const dataType hx_c = 1.0, hy_c = 1.0, hz_c = 1.0;           //h for central differences
    Point3D current_grad;
    FiniteVolumeSize3D finite_volume = { 1.0, 1.0, 1.0 };

    dataType lambda = 1.0;

    if (!pSegmentationParams->open_curve) {
        lambda = pSegmentationParams->lambda;
    }

    for (size_t k = 0; k < inputImage3D.height; k++)
    {
        for (size_t i = 0; i < inputImage3D.length; i++)
        {
            for (size_t j = 0; j < inputImage3D.width; j++)
            {
                getGradient3D(inputImage3D.imageDataPtr, inputImage3D.width, inputImage3D.length, inputImage3D.height, j, i, k, finite_volume, &current_grad);
                size_t xd = x_new(j, i, inputImage3D.width);

                //Absolute gradient
                abs_val_grad[k][xd] = norm3D(current_grad);

                //Edge detector
                edge_detector[k][xd] = edgeDetector(abs_val_grad[k][xd], pSegmentationParams->edge_detector_coef);

                //Similar intensity detector
                if (!pSegmentationParams->open_curve) {
                    similar_intensity_detector[k][xd] = similarIntensityDetector(inputImage3D.imageDataPtr[k][xd], pSegmentationParams->refence_intensity, pSegmentationParams->intensityCoef);
                }
            }
        }
    }
    //unsigned char _path_grad[] = "C:/Users/Konan Allaly/Documents/Tests/Curves/Output/gradient.raw";
    //manageFile(abs_val_grad, inputImage3D.length, inputImage3D.width, inputImage3D.height, _path_grad , STORE_DATA_RAW, BINARY_DATA, (Storage_Flags){ false, false });

    //get the velocity field
    for (size_t k = 0; k < inputImage3D.height; k++)
    {
        for (size_t i = 0; i < inputImage3D.length; i++)
        {
            for (size_t j = 0; j < inputImage3D.width; j++)
            {
                getGradient3D(edge_detector, inputImage3D.width, inputImage3D.length, inputImage3D.height, j, i, k, finite_volume, &current_grad);
                size_t xd = x_new(j, i, inputImage3D.width);
                pgrad_x[k][xd] = current_grad.x;
                pgrad_y[k][xd] = current_grad.y;
                pgrad_z[k][xd] = current_grad.z;
            }
        }
    }

    size_t current_i = 0;
    size_t current_j = 0;
    size_t current_k = 0;

    for (size_t i = 0; i < (pResultSegmentation->numPoints); i++)
    {
        pResultSegmentation->pPoints[i].x = pSegmentationParams->pinitial_condition->pPoints[i].x;
        pResultSegmentation->pPoints[i].y = pSegmentationParams->pinitial_condition->pPoints[i].y;
        pResultSegmentation->pPoints[i].z = pSegmentationParams->pinitial_condition->pPoints[i].z;
    }

    size_t iterPt = 0;
    size_t maxIterPt = oldSegmentation.numPoints;

    if (pSegmentationParams->open_curve) {
        iterPt = 1;
        maxIterPt = oldSegmentation.numPoints - 1;
    }

    //_l - lower, _g - greater, _c - current
    dataType vx, vy, vz;
    dataType rx_l, rx_g;
    dataType ry_l, ry_g;
    dataType rz_l, rz_g;
    dataType rx_c, ry_c, rz_c;
    dataType nx, ny, nz, dot, norm, tx, ty, tz;
    dataType h_g, h_c;
    dataType curv_x, curv_y, curv_z, mean_dist;

    dataType tau = pSegmentationParams->time_step_size;
    dataType mu = pSegmentationParams->mu;
    dataType eps = pSegmentationParams->eps;

    if (!pSegmentationParams->open_curve) 
    {
        for (size_t t = 0; t < pSegmentationParams->num_time_steps; t++)
        {
            for (size_t i = 0; i < (oldSegmentation.numPoints); i++)
            {
                oldSegmentation.pPoints[i].x = pResultSegmentation->pPoints[i].x;
                oldSegmentation.pPoints[i].y = pResultSegmentation->pPoints[i].y;
                oldSegmentation.pPoints[i].z = pResultSegmentation->pPoints[i].z;
            }

            for (size_t i = iterPt; i < maxIterPt; i++)
            {
                //let us move by curve points just in vector field
                current_i = (size_t)(oldSegmentation.pPoints[i].x + 0.5);
                current_j = (size_t)(oldSegmentation.pPoints[i].y + 0.5);
                current_k = (size_t)(oldSegmentation.pPoints[i].z + 0.5);

                //let us keep points inside the image
                if (current_i < 0) {
                    current_i = 0;
                }
                else if (current_i >= inputImage3D.length) {
                    current_i = inputImage3D.length - 1;
                }

                if (current_j < 0) {
                    current_j = 0;
                }
                else if (current_j >= inputImage3D.width) {
                    current_j = inputImage3D.width - 1;
                }

                if (current_k < 0) {
                    current_k = 0;
                }
                else if (current_k >= inputImage3D.height) {
                    current_k = inputImage3D.height - 1;
                }

                size_t xd = x_new(current_j, current_i, inputImage3D.width);

                rx_c = oldSegmentation.pPoints[i].x;
                ry_c = oldSegmentation.pPoints[i].y;
                rz_c = oldSegmentation.pPoints[i].z;

                if (i == 0) {
                    rx_l = oldSegmentation.pPoints[oldSegmentation.numPoints - 1].x;
                    ry_l = oldSegmentation.pPoints[oldSegmentation.numPoints - 1].y;
                    rz_l = oldSegmentation.pPoints[oldSegmentation.numPoints - 1].z;

                    rx_g = oldSegmentation.pPoints[1].x;
                    ry_g = oldSegmentation.pPoints[1].y;
                    rz_g = oldSegmentation.pPoints[1].z;
                }
                else if (i == oldSegmentation.numPoints - 1) {
                    rx_l = oldSegmentation.pPoints[i - 1].x;
                    ry_l = oldSegmentation.pPoints[i - 1].y;
                    rz_l = oldSegmentation.pPoints[i - 1].z;

                    rx_g = oldSegmentation.pPoints[0].x;
                    ry_g = oldSegmentation.pPoints[0].y;
                    rz_g = oldSegmentation.pPoints[0].z;
                }
                else {
                    rx_l = oldSegmentation.pPoints[i - 1].x;
                    ry_l = oldSegmentation.pPoints[i - 1].y;
                    rz_l = oldSegmentation.pPoints[i - 1].z;

                    rx_g = oldSegmentation.pPoints[i + 1].x;
                    ry_g = oldSegmentation.pPoints[i + 1].y;
                    rz_g = oldSegmentation.pPoints[i + 1].z;
                }

                h_c = (dataType)sqrt(pow(rx_c - rx_l, 2) + pow(ry_c - ry_l, 2) + pow(rz_c - rz_l, 2));
                h_g = (dataType)sqrt(pow(rx_g - rx_c, 2) + pow(ry_g - ry_c, 2) + pow(rz_g - ry_c, 2));

                norm = (dataType)sqrt(pow(rx_g - rx_l, 2) + pow(ry_g - ry_l, 2) + pow(rz_g - rz_l, 2));
                tx = (rx_g - rx_l) / norm;
                ty = (ry_g - ry_l) / norm;
                tz = (rz_g - rz_l) / norm;

                nx = ty * tz;
                ny = -2 * tx * tz;
                nz = tx * ty;

                dot = pgrad_x[current_k][xd] * nx + pgrad_y[current_k][xd] * ny + pgrad_z[current_k][xd] * nz;

                mean_dist = 2.0 / (h_g + h_c);
                curv_x = mean_dist * ((rx_g - rx_c) / h_g - ((rx_c - rx_l) / h_c));
                curv_y = mean_dist * ((ry_g - ry_c) / h_g - ((ry_c - ry_l) / h_c));
                curv_z = mean_dist * ((rz_g - rz_c) / h_g - ((rz_c - rz_l) / h_c));

                vx = mu * ((1 - lambda) * similar_intensity_detector[current_k][xd]
                    - lambda * dot) + eps * curv_x;

                vy = mu * ((1 - lambda) * similar_intensity_detector[current_k][xd]
                    - lambda * dot) + eps * curv_y;

                vz = mu * ((1 - lambda) * similar_intensity_detector[current_k][xd]
                    - lambda * dot) + eps * curv_z;

                pResultSegmentation->pPoints[i].x = oldSegmentation.pPoints[i].x + tau * vx * nx;
                pResultSegmentation->pPoints[i].y = oldSegmentation.pPoints[i].y + tau * vy * ny;
                pResultSegmentation->pPoints[i].z = oldSegmentation.pPoints[i].z + tau * vz * nz;

            }

        }
    }

    //Free memory
    free(pOldCurve);
    for (size_t k = 0; k < inputImage3D.height; k++) {
        free(pgrad_x[k]);
        free(pgrad_y[k]);
        free(pgrad_z[k]);
        free(abs_val_grad[k]);
        free(edge_detector[k]);

    }
    free(pgrad_x);
    free(pgrad_y);
    free(pgrad_z);
    free(abs_val_grad);
    free(edge_detector);

    if (!pSegmentationParams->open_curve)
    {
        for (size_t k = 0; k < inputImage3D.height; k++) {
            free(similar_intensity_detector[k]);
        }
        free(similar_intensity_detector);
    }

    return true;
}