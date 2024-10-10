#include "common_functions.h"
#include <stdlib.h>
#include <math.h>
#include "file.h"
#include "segmentation2D_lagrangean.h"
#include "solvers.h"

bool evolveBySingleStep3D(Image_Data* pimage, Image_Data* pedge, LinkedCurve3D* plinked_curve, SchemeData3D* pscheme_data, const Lagrangean3DSegmentationParameters* pparams);

void calculateCurvature3D(LinkedCurve3D* plinked_curve, SchemeData3D* pscheme_data);

void normal_velocity3D(Image_Data* pimage, Image_Data* pedge, LinkedCurve3D* plinked_curve, SchemeData3D* pscheme_data,
    void(*pget_velocity)(Image_Data*, double, double, double, double*, double*, double*),
    void(*pget_g2)(Image_Data*, double, double, double, double, double, double*),
    const double ref_intensity, const double g2_coef, const double eps, const double lambda);

void tang_velocity3D(LinkedCurve3D* plinked_curve, SchemeData* pscheme_data, const double omega);

bool semiCoefficients3D(LinkedCurve3D* plinked_curve, SchemeData* pscheme_data, const double eps, const double dt);

//===========================================================

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

    dataType lambda;
    if (!pSegmentationParams->open_curve) {
        lambda = pSegmentationParams->lambda;
    }
    else {
        lambda = 1.0;
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
    unsigned char _path_grad[] = "C:/Users/Konan Allaly/Documents/Tests/Curves/Output/edge_detector.raw";
    manageFile(edge_detector, inputImage3D.length, inputImage3D.width, inputImage3D.height, _path_grad , STORE_DATA_RAW, BINARY_DATA, (Storage_Flags){ false, false });

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

bool lagrangeanSemiImplicit3DCurveSegmentation(Image_Data inputImage3D, const Lagrangean3DSegmentationParameters* pSegmentationParams,
    unsigned char* pOutputPathPtr, Curve3D* pResultSegmentation)
{
    if (pSegmentationParams == NULL || pSegmentationParams->open_curve || pResultSegmentation == NULL) {
        return false;
    }

    if (pSegmentationParams->num_points < 3) {
        return false;
    }

    const dataType hx = 1.0, hy = 1.0, hz = 1.0;           //spatial discretization step
    const dataType hx_c = 1.0, hy_c = 1.0, hz_c = 1.0;    //h for central differences
    Point3D current_grad;
    FiniteVolumeSize3D finite_volume = { 1.0, 1.0, 1.0 };

    const size_t dataDimension = inputImage3D.length * inputImage3D.width;
    const size_t sliceSize = dataDimension * sizeof(dataType);
    const size_t sizeHeight = sizeof(dataType*) * inputImage3D.height;
    
    dataType** edge_detector = (dataType**)malloc(sizeHeight);
    for (size_t k = 0; k < inputImage3D.height; k++) {
        edge_detector[k] = (dataType*)malloc(sliceSize);
        if (edge_detector[k] == NULL)
            return false;
    }
    if (edge_detector == NULL)
        return false;

    //get edge detector
    for (size_t k = 0; k < inputImage3D.height; k++)
    {
        for (size_t i = 0; i < inputImage3D.length; i++)
        {
            for (size_t j = 0; j < inputImage3D.width; j++)
            {
                getGradient3D(inputImage3D.imageDataPtr, inputImage3D.width, inputImage3D.length, inputImage3D.height, j, i, k, finite_volume, &current_grad);
                size_t xd = x_new(j, i, inputImage3D.width);
                dataType norm_of_grad = norm3D(current_grad);
                edge_detector[k][xd] = edgeDetector(norm_of_grad, pSegmentationParams->edge_detector_coef);
            }
        }
    }
    Image_Data edge = { inputImage3D.height, inputImage3D.length, inputImage3D.width, edge_detector };

    resetIDGenerator();
    //let us consider single curve without topological changes

    bool isOrientedPositively = true; //does not metter for open curves

    if (!pSegmentationParams->open_curve)
    {
        isOrientedPositively = isCurveOrientedPositively(pSegmentationParams->pinitial_condition);
    }

    LinkedCurve linked_curve = createLinkedCurve();
    initializeLinkedCurve(pSegmentationParams->pinitial_condition, &linked_curve, !isOrientedPositively, !pSegmentationParams->open_curve);

    if (!pSegmentationParams->open_curve)
    {
        //it is still necessary to think which part of these data will be in the revised list
        size_t length_of_data = linked_curve.number_of_points + 2;
        SchemeData3D* pscheme_data = (SchemeData3D*)calloc(length_of_data, sizeof(SchemeData3D));

        for (size_t it = 1, res_it = 0; it <= pSegmentationParams->num_time_steps; it++)
        {
            if (length_of_data < linked_curve.number_of_points + 2)
            {
                free(pscheme_data);
                length_of_data = linked_curve.number_of_points + 2;
                pscheme_data = (SchemeData3D*)calloc(length_of_data, sizeof(SchemeData3D));
            }
            //evolve curve
            evolveBySingleStep3D(&inputImage3D, &edge, &linked_curve, pscheme_data, pSegmentationParams);
        }

        free(pscheme_data);
        for (size_t k = 0; k < inputImage3D.height; k++) {
            free(edge_detector[k]);
        }
        free(edge_detector);

        //LinkedPoint* pt = linked_curve.first_point;
        //for (size_t i = 0; i < linked_curve.number_of_points; i++)
        //{
        //    pResultSegmentation->pPoints[i].x = (dataType)pt->x;
        //    pResultSegmentation->pPoints[i].y = (dataType)pt->y;
        //    pResultSegmentation->pPoints[i].z = 1.0;//(dataType)pt->z;
        //    pt = pt->next;
        //}
        //releaseLinkedCurve(&linked_curve);

        return true;

    }

    return false;
}

bool evolveBySingleStep3D(Image_Data* pimage, Image_Data* pedge, LinkedCurve3D* plinked_curve, SchemeData3D* pscheme_data, const Lagrangean3DSegmentationParameters* pparams)
{
    if (plinked_curve == NULL || pscheme_data == NULL || pparams == NULL ||
        pimage->imageDataPtr == NULL || pedge->imageDataPtr == NULL) 
    {
        return false;
    }

    const double lambda = pparams->lambda;
    const double eps = pparams->eps;
    const double omega = pparams->omega;
    const double dt = pparams->time_step_size;

    calculateCurvature3D(plinked_curve, pscheme_data);
    normal_velocity3D(pimage, pedge, plinked_curve, pscheme_data, pparams->get_velocity, pparams->get_g2, pparams->refence_intensity, pparams->intensityCoef, eps, lambda);
    tang_velocity3D(plinked_curve, pscheme_data, omega);

    if (!semiCoefficients3D(plinked_curve, pscheme_data, eps, dt))
    {
        return false;
    }

    //TODO

    //////////////////////    X component ///////////////////////////////////////////////////////////
    //LinkedPoint3D* current_point = plinked_curve->first_point;
    //for (size_t i = 1; i <= plinked_curve->number_of_points; i++)
    //{
    //    pscheme_data[i].ps = pscheme_data[i].m * current_point->x +
    //        0.5 * (pscheme_data[i].beta_ps_expl * (current_point->next->y - current_point->previous->y)) +
    //        0.25 * (fmin(-pscheme_data[i].alfa, 0.0) * (current_point->previous->x - current_point->next->x) +
    //            fmin(pscheme_data[i].alfa, 0.0) * (current_point->next->x - current_point->previous->x));
    //    current_point = current_point->next;
    //}
    //sherman_morris(pscheme_data, plinked_curve->number_of_points);

    /////////////////////    Y component   ///////////////////////////////////////////////////////////
    //current_point = plinked_curve->first_point;
    //for (size_t i = 1; i <= plinked_curve->number_of_points; i++)
    //{
    //    pscheme_data[i].ps = pscheme_data[i].m * current_point->y -
    //        0.5 * (pscheme_data[i].beta_ps_expl * (current_point->next->x - current_point->previous->x)) +
    //        0.25 * (fmin(-pscheme_data[i].alfa, 0.0) * (current_point->previous->y - current_point->next->y) +
    //            fmin(pscheme_data[i].alfa, 0.0) * (current_point->next->y - current_point->previous->y));
    //    current_point = current_point->next;
    //}
    //current_point = plinked_curve->first_point;
    //for (size_t i = 1; i <= plinked_curve->number_of_points; i++)
    //{
    //    updatePoint(plinked_curve, current_point, pscheme_data[i].sol, current_point->y);
    //    current_point = current_point->next;
    //}
    //sherman_morris(pscheme_data, plinked_curve->number_of_points);

    /////////////////////    Z component   ///////////////////////////////////////////////////////////
    //current_point = plinked_curve->first_point;
    //for (size_t i = 1; i <= plinked_curve->number_of_points; i++)
    //{
    //    pscheme_data[i].ps = pscheme_data[i].m * current_point->z -
    //        0.5 * (pscheme_data[i].beta_ps_expl * (current_point->next->x - current_point->previous->x)) +
    //        0.25 * (fmin(-pscheme_data[i].alfa, 0.0) * (current_point->previous->y - current_point->next->y) +
    //            fmin(pscheme_data[i].alfa, 0.0) * (current_point->next->y - current_point->previous->y));
    //    current_point = current_point->next;
    //}
    //current_point = plinked_curve->first_point;
    //for (size_t i = 1; i <= plinked_curve->number_of_points; i++)
    //{
    //    updatePoint(plinked_curve, current_point, pscheme_data[i].sol, current_point->y);
    //    current_point = current_point->next;
    //}
    //sherman_morris(pscheme_data, plinked_curve->number_of_points);

    //current_point = plinked_curve->first_point;
    //for (size_t i = 1; i <= plinked_curve->number_of_points; i++)
    //{
    //    updatePoint(plinked_curve, current_point, current_point->x, pscheme_data[i].sol);
    //    current_point = current_point->next;
    //}

    return true;
}

//the mean curvature of three adjacent elements is calculated - 4 points
void calculateCurvature3D(LinkedCurve3D* plinked_curve, SchemeData3D* pscheme_data)
{
    if (plinked_curve == NULL || pscheme_data == NULL) 
    {
        return;
    }

    const size_t curve_length = plinked_curve->number_of_points;
    double phi;
    LinkedPoint3D* current_point = plinked_curve->first_point;

    double h_i_minus = 0.0;
    double h_i_plus = 0.0;
    double h_i = 0.0;

    double x_i_minus_2, x_i_minus_1, x_i, x_i_plus_1;
    double y_i_minus_2, y_i_minus_1, y_i, y_i_plus_1;
    double z_i_minus_2, z_i_minus_1, z_i, z_i_plus_1;

    //TODO : use 3D formulation, the current is 2D

    for (size_t i = 1; i <= curve_length; i++)
    {
        h_i_minus = current_point->previous->previous->distance_to_next;
        h_i_plus = current_point->distance_to_next;
        h_i = current_point->previous->distance_to_next;

        x_i_minus_2 = current_point->previous->previous->x;
        x_i_minus_1 = current_point->previous->x;
        x_i = current_point->x;
        x_i_plus_1 = current_point->next->x;

        y_i_minus_2 = current_point->previous->previous->y;
        y_i_minus_1 = current_point->previous->y;
        y_i = current_point->y;
        y_i_plus_1 = current_point->next->y;

        z_i_minus_2 = current_point->previous->previous->z;
        z_i_minus_1 = current_point->previous->z;
        z_i = current_point->z;
        z_i_plus_1 = current_point->next->z;

        //phi = (
        //    (x_i_minus_1 - x_i_minus_2) * (x_i_plus_1 - x_i) +
        //    (y_i_minus_1 - y_i_minus_2) * (y_i_plus_1 - y_i) +
        //    (z_i_minus_1 - z_i_minus_2) * (z_i_plus_1 - z_i)
        //    ) / (h_i_plus * h_i_minus);

        //if (phi > 1)
        //    phi = 1.;
        //else if (phi < -1)
        //    phi = -1.;

        //pscheme_data[i].curvature = 0.5 * signum(
        //    (x_i_minus_1 - x_i_minus_2) * (y_i_plus_1 - y_i) -
        //    (x_i_plus_1 - x_i) * (y_i_minus_1 - y_i_minus_2)
        //) * acos(phi) / h_i;

        //pscheme_data[i].curvature = 0.5 * signum(
        //    (y_i_minus_1 - y_i_minus_2) * (z_i_plus_1 - z_i) - (z_i_minus_1 - z_i_minus_2) * (y_i_plus_1 - y_i) + 
        //    (z_i_plus_1 - z_i) * (x_i_minus_1 - x_i_minus_2) - (x_i_minus_1 - x_i_minus_2) * (z_i_plus_1 - z_i) +
        //    (x_i_minus_1 - x_i_minus_2) * (y_i_plus_1 - y_i) - (y_i_minus_1 - y_i_minus_2) * (x_i_plus_1 - x_i)
        //) * acos(phi) / h_i;
        //
        //current_point = current_point->next;
    }

    //pscheme_data[0].curvature = pscheme_data[curve_length].curvature;
    //pscheme_data[curve_length + 1].curvature = pscheme_data[curve_length + 1].curvature;

}

//beta preparation
void normal_velocity3D(Image_Data* pimage, Image_Data* pedge, LinkedCurve3D* plinked_curve, SchemeData3D* pscheme_data,
    void(*pget_velocity)(Image_Data*, double, double, double, double*, double*, double*),
    void(*pget_g2)(Image_Data*, double, double, double, double, double, double*),
    const double ref_intensity, const double g2_coef, const double eps, const double lambda)
{
    
    if (plinked_curve == NULL || pscheme_data == NULL || pget_velocity == NULL ||
        pget_g2 == NULL || pimage == NULL || pedge == NULL)
    {
        return;
    }

    const size_t number_of_points = plinked_curve->number_of_points;

    LinkedPoint3D* current_point = plinked_curve->first_point;

    double h_i, h_i_plus;
    double tx, ty, tz;
    double vx, vy, vz;
    double nvx, nvy, nvz;
    double grad_x, grad_y, grad_z;
    double dot, norm;
    double n1_x, n1_y, n1_z;
    double n2_x, n2_y, n2_z;
    double curv_x, curv_y, curv_z;

    for (size_t i = 1; i <= number_of_points; i++)
    {
        h_i = current_point->previous->distance_to_next;
        h_i_plus = current_point->distance_to_next;

        tx = (current_point->next->x - current_point->previous->x) / (h_i_plus + h_i);
        ty = (current_point->next->y - current_point->previous->y) / (h_i_plus + h_i);
        tz = (current_point->next->z - current_point->previous->z) / (h_i_plus + h_i);

        (*pget_velocity)(pedge, current_point->x, current_point->y, current_point->z, &vx, &vy, &vz);

        dot = tx * vx + ty * vy + tz * vz;

        nvx = -vx - dot * vx;
        nvy = -vy - dot * vy;
        nvz = -vz - dot * vz;

        Point3D pnorm = { nvx, nvy, nvz };
        norm = norm3D(pnorm);

        n1_x = nvx / norm;
        n1_y = nvy / norm;
        n1_z = nvz / norm;

        n2_x = n1_y * tz - n1_z * ty;
        n2_y = n1_z * tx - n1_x * tz;
        n2_z = n1_x * ty - n1_y * tx;

        //curv_x = (2.0 / (h_i_plus + h_i) * ((current_point->next->x - current_point->x) / h_i_plus - (current_point->x - current_point->previous->x) / h_i));
        //curv_y = (2.0 / (h_i_plus + h_i) * ((current_point->next->y - current_point->y) / h_i_plus - (current_point->y - current_point->previous->y) / h_i));
        //curv_z = (2.0 / (h_i_plus + h_i) * ((current_point->next->z - current_point->z) / h_i_plus - (current_point->z - current_point->previous->z) / h_i));

        current_point = current_point->next;
    }

    //for (size_t i = 1; i <= number_of_points; i++)
    //{
    //    pscheme_data[i].beta_ps_expl = pscheme_data[i].f;
    //    pscheme_data[i].beta = pscheme_data[i].curvature * eps - pscheme_data[i].beta_ps_expl;
    //}

    //pscheme_data[0].f = pscheme_data[number_of_points].f;
    //pscheme_data[number_of_points + 1].f = pscheme_data[1].f;
    //pscheme_data[0].beta = pscheme_data[number_of_points].beta;
    //pscheme_data[number_of_points + 1].beta = pscheme_data[1].beta;
    //pscheme_data[0].beta_ps_expl = pscheme_data[number_of_points].beta_ps_expl;
    //pscheme_data[number_of_points + 1].beta_ps_expl = pscheme_data[1].beta_ps_expl;
}

void tang_velocity3D(LinkedCurve3D* plinked_curve, SchemeData* pscheme_data, const double omega)
{
    //TODO : adapt to 3D formulation.

    int i;
    double mean = 0;
    const size_t number_of_points = plinked_curve->number_of_points;
    const double curve_length = plinked_curve->length;
    LinkedPoint3D* current_point = plinked_curve->first_point;
    const double avg_length = curve_length / number_of_points;
    double h_i = -1;

    for (i = 1; i <= number_of_points; i++)
    {
        h_i = current_point->previous->distance_to_next;

        mean += pscheme_data[i].beta * pscheme_data[i].curvature * h_i;

        current_point = current_point->next;
    }
    mean /= curve_length;

    // the alpha of the first point in the sequence - it will therefore not move in the tangential direction
    pscheme_data[1].alfa = 0.0;

    double alpha_sum = 0;

    current_point = plinked_curve->first_point->next;
    for (i = 2; i <= number_of_points; i++)
    {
        h_i = current_point->previous->distance_to_next;

        pscheme_data[i].alfa = pscheme_data[i - 1].alfa +
            pscheme_data[i].beta * pscheme_data[i].curvature * h_i -
            mean * h_i + omega * (avg_length - h_i);

        alpha_sum += pscheme_data[i].alfa;

        current_point = current_point->next;
    }
    pscheme_data[0].alfa = pscheme_data[number_of_points].alfa;
    pscheme_data[number_of_points + 1].alfa = pscheme_data[1].alfa;
}

bool semiCoefficients3D(LinkedCurve3D* plinked_curve, SchemeData* pscheme_data,const double eps, const double dt)
{
    if (plinked_curve == NULL || pscheme_data == NULL) 
    {
        return false;
    }

    double h_i = -1;
    double h_i_plus = -1;
    LinkedPoint* current_point = plinked_curve->first_point;

    for (size_t i = 1; i <= plinked_curve->number_of_points; i++)
    {
        h_i = current_point->previous->distance_to_next;
        h_i_plus = current_point->distance_to_next;

        pscheme_data[i].b = -0.5 * fmax(-pscheme_data[i].alfa, 0.0) - 1.0 / h_i * eps;//lower diagonal
        pscheme_data[i].c = -0.5 * fmax(pscheme_data[i].alfa, 0.0) - 1.0 / h_i_plus * eps;//upper diagonal

        pscheme_data[i].m = (h_i_plus + h_i) / (2.0 * dt);
        pscheme_data[i].a = pscheme_data[i].m - (pscheme_data[i].b + pscheme_data[i].c);//stiffness matrix

        if (fabs(pscheme_data[i].b) + fabs(pscheme_data[i].c) > fabs(pscheme_data[i].a) || pscheme_data[i].a < 0) {
            //the matrix is not positive dominant
            return false;
        }

        current_point = current_point->next;
    }

    return true;
}