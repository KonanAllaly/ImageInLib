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

    const dataType edge_detector_coef = 10;

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

    const dataType similar_intensity_detector_coef = 10;
    const dataType hx = 1.0, hy = 1.0, hz = 1.0;        //spatial discretization step
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

    //get edge detector
    for (size_t k = 0; k < inputImage3D.height; k++)
    {
        for (size_t i = 0; i < sliceDimension; i++)
        {
            edge_detector[k][i] = edgeDetector(abs_val_grad[k][i], edge_detector_coef);
        }
    }

    Point3D centroid = get3dCurveCentroid(pSegmentationParams->pinitial_condition);

    size_t centroid_i = (size_t)(centroid.y + 0.5);
    size_t centroid_j = (size_t)(centroid.x + 0.5);
    size_t centroid_k = (size_t)(centroid.z + 0.5);
    dataType ref_intensity = inputImage3D.imageDataPtr[centroid_k][x_new(centroid_i, centroid_j, inputImage3D.length)];

    if (!pSegmentationParams->open_curve)
    {
        for (size_t k = 0; k < inputImage3D.height; k++)
        {
            for (size_t i = 0; i < sliceDimension; i++) 
            {
                similar_intensity_detector[k][i] = similarIntensityDetector(inputImage3D.imageDataPtr[k][i], ref_intensity, similar_intensity_detector_coef);
            }
        }
    }

    Image_Data edgeDetector = { inputImage3D.height, inputImage3D.length, inputImage3D.width, edge_detector };

    //get velocity
    for (size_t k = 0; k < inputImage3D.height; k++)
    {
        for (size_t i = 0; i < inputImage3D.length; i++)
        {
            for (size_t j = 0; j < inputImage3D.width; j++)
            {
                getGradient3D(edgeDetector.imageDataPtr, edgeDetector.length, edgeDetector.width, edgeDetector.height, i, j, k, finite_volume, &current_grad);
                xd = x_new(i, j, inputImage3D.length);
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
                current_j = inputImage3D.height - 1;
            }

            xd = x_new(current_i, current_j, inputImage3D.length);

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

            norm = (dataType)sqrt(pow(rx_l - rx_g, 2) + pow(ry_l - ry_g, 2) + pow(rz_l - ry_g, 2));
            tx = (rx_g - rx_l) / norm;
            ty = (ry_g - ry_l) / norm;
            tz = (rz_g - rz_l) / norm;

            // First choice of othogonal vector
            nx = ty * tz;
            ny = -2 * ty * tz;
            nz = tx * ty;

            //// Second choice of orthogonal vector
            //nx = ty;
            //ny = -tx;
            //nz = 0.0;

            dot = pgrad_x[current_k][xd] * nx + pgrad_y[current_k][xd] * ny + pgrad_y[current_k][xd] * nz;

            //basic velocity - gradients of edge detector
            vx = -pgrad_x[current_k][xd];
            vy = -pgrad_y[current_k][xd];
            vz = -pgrad_z[current_k][xd];

            //projection of velocity field in direction of normal vector
            dot = vx * nx + vy * ny + vz * nz;

            h_c = (dataType)sqrt(pow(rx_c - rx_l, 2) + pow(ry_c - ry_l, 2) + pow(rz_c - rz_l, 2));
            h_g = (dataType)sqrt(pow(rx_g - rx_c, 2) + pow(ry_g - ry_c, 2) + pow(rz_g - ry_c, 2));

            //mu * normal vector * projection + eps * normal * curvature
            vx = pSegmentationParams->mu * (lambda * dot * nx +
                ((dataType)1.0 - lambda) * similar_intensity_detector[current_k][xd] * nx) +
                pSegmentationParams->eps * (dataType)(2.0 / (h_g + h_c)) * ((rx_g - rx_c) / h_g - ((rx_c - rx_l) / h_c));
            vy = pSegmentationParams->mu * dot * ny + (lambda * dot * ny +
                ((dataType)1.0 - lambda) * similar_intensity_detector[current_k][xd] * ny) +
                pSegmentationParams->eps * (dataType)(2.0 / (h_g + h_c)) * ((ry_g - ry_c) / h_g - ((ry_c - ry_l) / h_c));
            vz = pSegmentationParams->mu * dot * nz + (lambda * dot * nz +
                ((dataType)1.0 - lambda) * similar_intensity_detector[current_k][xd] * nz) +
                pSegmentationParams->eps * (dataType)(2.0 / (h_g + h_c)) * ((rz_g - rz_c) / h_g - ((rz_c - rz_l) / h_c));

            //it is just simple motion in vector field

            pResultSegmentation->pPoints[i].x = oldSegmentation.pPoints[i].x +
                pSegmentationParams->time_step_size * vx;

            pResultSegmentation->pPoints[i].y = oldSegmentation.pPoints[i].y +
                pSegmentationParams->time_step_size * vy;

            pResultSegmentation->pPoints[i].z = oldSegmentation.pPoints[i].z +
                pSegmentationParams->time_step_size * vz;
        }
    }

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