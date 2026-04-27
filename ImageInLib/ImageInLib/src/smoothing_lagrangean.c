#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "file.h"

#include "smoothing_lagrangean.h"
#include "../src/imageInterpolation.h"

bool smoothingByLagrangeanCurveEvolution(const LagrangeanSmoothingParameters* pSmoothingParameters,
    unsigned char* pOutputPathPtr, Curve3D* pResultCurve)
{
    //check if the pointers a well allocated
    if (pSmoothingParameters == NULL || pResultCurve == NULL) 
    {
        return false;
    }

    //check that we have at least three points at the biginning
    if (pSmoothingParameters->num_points < 3) 
    {
        return false;
    }

    //let us consider single curve without topological changes
    resetIDGenerator();

    bool isOrientedPositively = true;
	bool isCurveClosed = !pSmoothingParameters->open_curve;

    if (isCurveClosed)
    {
        isOrientedPositively = is3dCurveOrientedPositively(pSmoothingParameters->pinitial_condition);
    }

    //create initial linked curve
    LinkedCurve3D initial_curve = create3dLinkedCurve();

    //Initialize the evolving linked curve
    initialize3dLinkedCurve(pSmoothingParameters->pinitial_condition, &initial_curve, false, isCurveClosed);

    resetIDGenerator();

    //create evolving linked curve
    LinkedCurve3D evolving_curve = create3dLinkedCurve();

    //Initialize the evolving linked curve
    initialize3dLinkedCurve(pSmoothingParameters->pinitial_condition, &evolving_curve, false, isCurveClosed);

    //Data size for the scheme, 
    // we add 2 to be sure that we have enough memory
    // in case of topological changes (even if we do not expect any)
	size_t length_of_data = pSmoothingParameters->num_points + 2; 
    SchemeData3D* pscheme_data = (SchemeData3D*)calloc(length_of_data, sizeof(SchemeData3D));

	size_t it = 1; //count of time steps
	
    // Test if we are at the first time step, 
    // in this case we do not have to compute the closest point 
    // on the initial curve since it is the same as the evolving curve
    bool isFirstTimeStep = true; 

	double initilal_length = evolving_curve.length;
    double current_length = initilal_length, previous_length = initilal_length;
    double length_diff = 0.0;

    ////Save the total lenght
    //const char * savingPath = "C:/Users/Konan Allaly/Documents/Tests/Curves/Output/Torsion/curve_length_p3.csv";
    ////savingPath = outputPath + "centered_smoothed_p1_without_attr.csv";
    //FILE* file_save;
    //if (fopen_s(&file_save, savingPath, "w") != 0) {
    //    printf("Enable to open");
    //    return false;
    //}
    //fprintf(file_save, "%d,%lf\n", it, evolving_curve.length);

    do
    {
        if (it == 1)
        {
            isFirstTimeStep = true;
        }

        if (length_of_data < evolving_curve.number_of_points + 2)
        {
            free(pscheme_data);
            length_of_data = evolving_curve.number_of_points + 2;
            pscheme_data = (SchemeData3D*)calloc(length_of_data, sizeof(SchemeData3D));
        }
        
        //evolve curve
        evolveForSmoothingBySingleStep(&initial_curve, &evolving_curve, pscheme_data, pSmoothingParameters, isFirstTimeStep);
        it++;

        isFirstTimeStep = false;

		current_length = evolving_curve.length;
		length_diff = fabs(current_length - previous_length);

        //fprintf(file_save, "%d,%lf\n", it, evolving_curve.length);

		previous_length = current_length;

    } while (it <= pSmoothingParameters->num_time_steps);
    //it <= pSmoothingParameters->num_time_steps ||
    //length_diff > pSmoothingParameters->tolerance

    LinkedPoint3D* final_curve = evolving_curve.first_point;

    for (size_t i = 0; i < evolving_curve.number_of_points; i++)
    {
        pResultCurve->pPoints[i].x = (dataType)final_curve->x;
        pResultCurve->pPoints[i].y = (dataType)final_curve->y;
        pResultCurve->pPoints[i].z = (dataType)final_curve->z;
        final_curve = final_curve->next;
    }

    //fclose(file_save);

	final_curve = NULL;
    free(pscheme_data);
    release3dLinkedCurve(&evolving_curve);
    release3dLinkedCurve(&initial_curve);

    return true;
}

bool getClosestPointToCurve(LinkedCurve3D* initial_curve, LinkedPoint3D* evolving_point, bool isTheFirstTimeStep, double* px, double* py, double* pz)
{
    if(initial_curve == NULL || evolving_point == NULL || px == NULL || py == NULL || pz == NULL)
    {
        return false;
	}

    if(isTheFirstTimeStep)
    {
        *px = evolving_point->x;
        *py = evolving_point->y;
        *pz = evolving_point->z;
    }
    else
    {
        double dist_min = INFINITY, x, y, z;
        LinkedPoint3D* reference_point = initial_curve->first_point;
        for (size_t i = 1; i < initial_curve->number_of_points; i++)
        {
            //Step 1 : compute t

            //pt1 = r^{0}_{j} - r^{0}_{j+1}
            double x_pt1 = reference_point->next->x - reference_point->x;
            double y_pt1 = reference_point->next->y - reference_point->y;
            double z_pt1 = reference_point->next->z - reference_point->z;

            //pt2 = r^{n}_i - r^{0}_{j}
            double x_pt2 = evolving_point->x - reference_point->x;
            double y_pt2 = evolving_point->y - reference_point->y;
            double z_pt2 = evolving_point->z - reference_point->z;

            double numerator = x_pt1 * x_pt2 + y_pt1 * y_pt2 + z_pt1 * z_pt2;//scalar product : p1.p2
            double denominator = x_pt1 * x_pt1 + y_pt1 * y_pt1 + z_pt1 * z_pt1;//scalar product : p1.p1
            double t = numerator / denominator;

            if (denominator != 0.0)
            {
                double t_cl = fmax(0.0, fmin(1.0, t));
                x = reference_point->x + t_cl * x_pt1;
                y = reference_point->y + t_cl * y_pt1;
                z = reference_point->z + t_cl * z_pt1;
            }
            else
            {
                //For debugging, this should not happen since we should not have two coincident points in the initial curve
                printf("! Points merging !");
                return false;
            }

            double dist = sqrt((evolving_point->x - x) * (evolving_point->x - x)
                + (evolving_point->y - y) * (evolving_point->y - y)
                + (evolving_point->z - z) * (evolving_point->z - z));

            if (dist < dist_min)
            {
                dist_min = dist;
                *px = x;
                *py = y;
                *pz = z;
            }
            reference_point = reference_point->next;
        }
    }

    return true;
}

bool normalVelocitySmoothing(LinkedCurve3D* initial_curve, LinkedCurve3D* evolving_curve, SchemeData3D* pscheme_data,
    bool isTheFirstTimeStep, const double delta, const double lambda)
{

    //check if the pointers a well allocated
    if (initial_curve == NULL || evolving_curve == NULL || pscheme_data == NULL || evolving_curve->number_of_points < 2)
    {
        return false;
    }

	//No topological change, therefore the number of points is the same during the evolution
    const size_t number_of_points = evolving_curve->number_of_points;

    LinkedPoint3D* evolving_point = evolving_curve->first_point;
    LinkedPoint3D* p_ref_initial = NULL;

    double x = 0.0, y = 0.0, z = 0.0;
    
    double tx, ty, tz;
    double n1_x, n1_y, n1_z;
    double n2_x, n2_y, n2_z;
    double denominator;
    double r0_minus_r_x, r0_minus_r_y, r0_minus_r_z;
    double attract1, attract2;
    for (size_t i = 1; i <= number_of_points; i++)
    {
        if (i > 1 && i < number_of_points)
        {
            
            //Get the attraction vector
            x = evolving_point->x;
            y = evolving_point->y;
            z = evolving_point->z;
            if(!getClosestPointToCurve(initial_curve, evolving_point, isTheFirstTimeStep, &x, &y, &z))
            {
                return false;
            }
            r0_minus_r_x = x - evolving_point->x;
            r0_minus_r_y = y - evolving_point->y;
            r0_minus_r_z = z - evolving_point->z;

            double h_i = evolving_point->previous->distance_to_next;
            double h_i_plus = evolving_point->distance_to_next;

            if(h_i == 0.0|| h_i_plus == 0.0)
            {
                //TODO: think about why this happens and what to do...
                return false;
            }
            double som_dist = h_i_plus + h_i;

            //Compute the tangent vector components
            tx = (evolving_point->next->x - evolving_point->previous->x) / som_dist;
            ty = (evolving_point->next->y - evolving_point->previous->y) / som_dist;
            tz = (evolving_point->next->z - evolving_point->previous->z) / som_dist;
            denominator = 1.0 + tz;
            
            //Compute the normal plane vectors
            if(denominator == 0.0)
            {
                n1_x = 0.0;
                n1_y = -1.0;
                n1_z = 0.0;

                n2_x = -1.0;
                n2_y = 0.0;
                n2_z = 0.0;

            }
            else
            {
                //Compute N1 components
                n1_x = 1.0 - (tx * tx) / denominator;
                n1_y = -(tx * ty) / denominator;
                n1_z = -tx;

                //Compute N2 components
                n2_x = -(tx * ty) / denominator;
                n2_y = 1.0 - (ty * ty) / denominator;
                n2_z = -ty;
            }

            //Compute the curvature vector
            double coef = 2.0 / som_dist;
            double curv_x = coef * (((evolving_point->next->x - evolving_point->x) / h_i_plus) - ((evolving_point->x - evolving_point->previous->x) / h_i));
            double curv_y = coef * (((evolving_point->next->y - evolving_point->y) / h_i_plus) - ((evolving_point->y - evolving_point->previous->y) / h_i));
            double curv_z = coef * (((evolving_point->next->z - evolving_point->z) / h_i_plus) - ((evolving_point->z - evolving_point->previous->z) / h_i));
            
            //Compute the curvature component in the new normal basis
            pscheme_data[i].k1 = curv_x * n1_x + curv_y * n1_y + curv_z * n1_z;
            pscheme_data[i].k2 = curv_x * n2_x + curv_y * n2_y + curv_z * n2_z;

            attract1 = lambda * (r0_minus_r_x * n1_x + r0_minus_r_y * n1_y + r0_minus_r_z * n1_z);
            attract2 = lambda * (r0_minus_r_x * n2_x + r0_minus_r_y * n2_y + r0_minus_r_z * n2_z);

            pscheme_data[i].normal_x = attract1 * n1_x + attract2 * n2_x;
            pscheme_data[i].normal_y = attract1 * n1_y + attract2 * n2_y;
            pscheme_data[i].normal_z = attract1 * n1_z + attract2 * n2_z;

            //Compute the normal velocity component
            pscheme_data[i].u = delta * pscheme_data[i].k1 + r0_minus_r_x * n1_x + r0_minus_r_y * n1_y + r0_minus_r_z * n1_z;
            pscheme_data[i].v = delta * pscheme_data[i].k2 + r0_minus_r_x * n2_x + r0_minus_r_y * n2_y + r0_minus_r_z * n2_z;

        }
        else
        {
            //when this happens, we are at the end points of an open curve
            pscheme_data[i].k1 = 0.0;
            pscheme_data[i].k2 = 0.0;
			pscheme_data[i].u = 0.0;
			pscheme_data[i].v = 0.0;
            pscheme_data[i].normal_x = 0.0;
            pscheme_data[i].normal_y = 0.0;
            pscheme_data[i].normal_z = 0.0;
        }
        evolving_point = evolving_point->next;
    }

    return true;
}

bool tangentialVelocitySmoothing(LinkedCurve3D* evolving_curve, SchemeData3D* pscheme_data, const double omega, bool isCurveOpen)
{
    if(evolving_curve == NULL || pscheme_data == NULL || evolving_curve->number_of_points < 2)
    {
        return false;
	}

    dataType mean = 0.0;
    const size_t number_of_points = evolving_curve->number_of_points;
    const dataType curve_length = evolving_curve->length;
    dataType h_i = INFINITY;

    dataType avg_length = INFINITY;
    if (isCurveOpen)
    {
        avg_length = curve_length / (dataType)(number_of_points - 1);//The curve is open
    }
    else
    {
        avg_length = curve_length / (dataType)(number_of_points);//The curve is closed
    }

    LinkedPoint3D* current_point = evolving_curve->first_point;
    for (size_t i = 1; i <= number_of_points; i++)
    {
        if (i == 1)
        {
            //if it is not the first point
            h_i = current_point->distance_to_next;
        }
        else
        {
            //first point
            h_i = current_point->previous->distance_to_next;;
        }

        mean += h_i * (pscheme_data[i].u * pscheme_data[i].k1 + pscheme_data[i].v * pscheme_data[i].k2);

        current_point = current_point->next;
    }

    mean /= curve_length;
    
    pscheme_data[0].alfa = 0.0;
    pscheme_data[1].alfa = 0.0;

    current_point = evolving_curve->first_point;
    for (size_t i = 1; i <= number_of_points; i++)
    {
        if (i == 1)
        {
            h_i = current_point->distance_to_next;
        }
        else
        {
            h_i = current_point->previous->distance_to_next;
        }

        pscheme_data[i].alfa = pscheme_data[i - 1].alfa + h_i * (pscheme_data[i].u * pscheme_data[i].k1 + pscheme_data[i].v * pscheme_data[i].k2) - h_i * mean + omega * (avg_length - h_i);
        current_point = current_point->next;
    }

    pscheme_data[number_of_points].alfa = 0.0;
    pscheme_data[number_of_points + 1].alfa = 0.0;
    
    return true;
}

bool coefficientsSmoothing(LinkedCurve3D* evolving_curve, SchemeData3D* pscheme_data, const double delta, const double tau)
{
    if (evolving_curve == NULL || pscheme_data == NULL)
    {
        return false;
    }

    double h_i = -1;
    double h_i_plus = -1;
    LinkedPoint3D* current_point = evolving_curve->first_point;
    LinkedPoint3D* previous_point;

    bool is_curve_closed = evolving_curve->first_point->previous != NULL;

    for (size_t i = 1; i <= evolving_curve->number_of_points; i++)
    {

        if (i > 1 && i < evolving_curve->number_of_points)
        {
            previous_point = current_point->previous;
            h_i = previous_point->distance_to_next;

            if (is_curve_closed || i < evolving_curve->number_of_points)
            {
                h_i_plus = current_point->distance_to_next;
            }
            else
            {
                h_i_plus = h_i;
            }

            pscheme_data[i].a = -delta / h_i - 0.5 * fmax(-pscheme_data[i].alfa, 0);      //lower diagonal
            pscheme_data[i].c = -delta / h_i_plus - 0.5 * fmax(pscheme_data[i].alfa, 0); //upper diagonal
            pscheme_data[i].m = (h_i_plus + h_i) / (2.0 * tau);
            pscheme_data[i].b = pscheme_data[i].m - (pscheme_data[i].a + pscheme_data[i].c);//diagonal
        }
        else
        {
            if (i == 1)
            {
                h_i = current_point->distance_to_next;
                h_i_plus = h_i;
            }
            else
            {
                h_i = current_point->previous->distance_to_next;
                h_i_plus = h_i;
            }

            pscheme_data[i].a = 0.0;
            pscheme_data[i].c = 0.0;
            pscheme_data[i].m = (h_i_plus + h_i) / (2.0 * tau);
            pscheme_data[i].b = 1.0;
        }

        current_point = current_point->next;
    }

    return true;
}

bool evolveForSmoothingBySingleStep(LinkedCurve3D* initial_curve, LinkedCurve3D* evolving_curve, SchemeData3D* pscheme_data,
    const LagrangeanSmoothingParameters* pSmoothingParameters, bool isTheFirstTimeStep)
{
    //check if the pointers a well allocated
    if (initial_curve == NULL || evolving_curve == NULL || pscheme_data == NULL || pSmoothingParameters == NULL)
    {
        return false;
    }

    const double delta = pSmoothingParameters->delta;
	const double lambda = pSmoothingParameters->lambda;
    const double omega = pSmoothingParameters->omega;
    const double tau = pSmoothingParameters->time_step_size;

	bool is_curve_open = pSmoothingParameters->open_curve;

    if (!normalVelocitySmoothing(initial_curve, evolving_curve, pscheme_data, isTheFirstTimeStep, delta, lambda))
    {
        return false;
    }

    //function to compute the tangential velocity
    if(!tangentialVelocitySmoothing(evolving_curve, pscheme_data, omega, is_curve_open))
    {
        return false;
	}

	//Set the coefficients of the linear system for the implicit scheme
    if (!coefficientsSmoothing(evolving_curve, pscheme_data, delta, tau))
    {
        return false;
    }

    LinkedPoint3D* current_point = evolving_curve->first_point;

    //////////////////////    X component ///////////////////////////////////////////////////////////

    for (size_t i = 1; i <= evolving_curve->number_of_points; i++)
    {

        if (i == 1 || i == evolving_curve->number_of_points)
        {
            pscheme_data[i].ps = current_point->x;
        }
        else
        {
			//pscheme_data[i].ps = pscheme_data[i].m * current_point->x 
            //    - 0.5 * fmin(pscheme_data[i].alfa, 0) * (current_point->x - current_point->next->x)
            //    - 0.5 * fmin(-pscheme_data[i].alfa, 0) * (current_point->x - current_point->previous->x)
            //    + lambda * tau * pscheme_data[i].m * pscheme_data[i].w * pscheme_data[i].normal_x;
            
            pscheme_data[i].ps = pscheme_data[i].m * current_point->x 
                - 0.5 * fmin(pscheme_data[i].alfa, 0) * (current_point->x - current_point->next->x) 
                - 0.5 * fmin(-pscheme_data[i].alfa, 0) * (current_point->x - current_point->previous->x) 
                + pscheme_data[i].m * tau * pscheme_data[i].normal_x;

        }
        current_point = current_point->next;
    }

    //if (is_curve_closed)
    //{
    //    calculate_by_thomas3D(pscheme_data, evolving_curve->number_of_points);
    //}
    //else
    //{
    //    sherman_morris3D(pscheme_data, evolving_curve->number_of_points);
    //}

    //We are dealing only with open curves
    calculate_by_thomas3D(pscheme_data, evolving_curve->number_of_points);

    /////////////////////    Y component   ///////////////////////////////////////////////////////////

    current_point = evolving_curve->first_point;
    for (size_t i = 1; i <= evolving_curve->number_of_points; i++)
    {

        if (i == 1 || i == evolving_curve->number_of_points)
        {
            pscheme_data[i].ps = current_point->y;
        }
        else
        {
            //pscheme_data[i].ps = pscheme_data[i].m * current_point->y
            //    - 0.5 * fmin(pscheme_data[i].alfa, 0) * (current_point->y - current_point->next->y)
            //    - 0.5 * fmin(-pscheme_data[i].alfa, 0) * (current_point->y - current_point->previous->y)
            //    + lambda * tau * pscheme_data[i].m * pscheme_data[i].w * pscheme_data[i].normal_y;

            //pscheme_data[i].ps = pscheme_data[i].m * current_point->y
            //    - 0.5 * fmin(pscheme_data[i].alfa, 0) * (current_point->y - current_point->next->y)
            //    - 0.5 * fmin(-pscheme_data[i].alfa, 0) * (current_point->y - current_point->previous->y);
            //    + lambda * pscheme_data[i].w * 0.5 * (current_point->next->x - current_point->previous->x);

            pscheme_data[i].ps = pscheme_data[i].m * current_point->y 
                - 0.5 * fmin(pscheme_data[i].alfa, 0) * (current_point->y - current_point->next->y) 
                - 0.5 * fmin(-pscheme_data[i].alfa, 0) * (current_point->y - current_point->previous->y) 
                + pscheme_data[i].m * tau * pscheme_data[i].normal_y;

        }
        current_point = current_point->next;
    }

    current_point = evolving_curve->first_point;
    for (size_t i = 1; i <= evolving_curve->number_of_points; i++)
    {
        //the same scheme data variable is reused for all components
        //so we update the components before solving the next system
        update3dPoint(evolving_curve, current_point, pscheme_data[i].sol, current_point->y, current_point->z);
        current_point = current_point->next;
    }

    //if (is_curve_closed)
    //{
    //    calculate_by_thomas3D(pscheme_data, evolving_curve->number_of_points);
    //}
    //else
    //{
    //    sherman_morris3D(pscheme_data, evolving_curve->number_of_points);
    //}

    //We are dealing only with open curves
    calculate_by_thomas3D(pscheme_data, evolving_curve->number_of_points);

    /////////////////////    Z component   ///////////////////////////////////////////////////////////

    current_point = evolving_curve->first_point;
    for (size_t i = 1; i <= evolving_curve->number_of_points; i++)
    {

        if (i == 1 || i == evolving_curve->number_of_points)
        {
            pscheme_data[i].ps = current_point->z;
        }
        else
        {
            //pscheme_data[i].ps = pscheme_data[i].m * current_point->z
            //    - 0.5 * fmin(pscheme_data[i].alfa, 0) * (current_point->z - current_point->next->z)
            //    - 0.5 * fmin(-pscheme_data[i].alfa, 0) * (current_point->z - current_point->previous->z)
            //    + lambda * tau * pscheme_data[i].m * pscheme_data[i].w * pscheme_data[i].normal_z;

            //pscheme_data[i].ps = pscheme_data[i].m * current_point->z
            //    - 0.5 * fmin(pscheme_data[i].alfa, 0) * (current_point->z - current_point->next->z)
            //    - 0.5 * fmin(-pscheme_data[i].alfa, 0) * (current_point->z - current_point->previous->z);

            pscheme_data[i].ps = pscheme_data[i].m * current_point->z 
                - 0.5 * fmin(pscheme_data[i].alfa, 0) * (current_point->z - current_point->next->z) 
                - 0.5 * fmin(-pscheme_data[i].alfa, 0) * (current_point->z - current_point->previous->z) 
                + pscheme_data[i].m * tau * pscheme_data[i].normal_z;

        }
        current_point = current_point->next;
    }

    current_point = evolving_curve->first_point;
    for (size_t i = 1; i <= evolving_curve->number_of_points; i++)
    {
        //The same scheme_data variable is reused for all components
        //so we update the components before solving the next system
        update3dPoint(evolving_curve, current_point, current_point->x, pscheme_data[i].sol, current_point->z);
        current_point = current_point->next;
    }

    //if (is_curve_closed)
    //{
    //    calculate_by_thomas3D(pscheme_data, evolving_curve->number_of_points);
    //}
    //else
    //{
    //    sherman_morris3D(pscheme_data, evolving_curve->number_of_points);
    //}

    //We are dealing only with open curves
    calculate_by_thomas3D(pscheme_data, evolving_curve->number_of_points);

    current_point = evolving_curve->first_point;
    for (size_t i = 1; i <= evolving_curve->number_of_points; i++)
    {
        //The same scheme data variable is reused for all components
        //so we update the components before solving the next system
        update3dPoint(evolving_curve, current_point, current_point->x, current_point->y, pscheme_data[i].sol);
        current_point = current_point->next;
    }

    return true;
}

//Function to investiage curvature, torsion and tangent vector

bool computeCurvatureTorsionAndTangent(Curve3D* curve, const char * save_path)
{
    if (curve == NULL || curve->numPoints < 2 || save_path == NULL)
    {
        return false;
    }

    size_t new_path_length;

    //Save curvature
	const char* extension_curvature = "_curvature.csv";
    new_path_length = strlen(save_path) + strlen(extension_curvature) + 1;
    char* save_curvature = malloc(new_path_length);
    strcpy_s(save_curvature, new_path_length, save_path);
    strcat_s(save_curvature, new_path_length, extension_curvature);
    
    FILE* file_save_curvature;
    if (fopen_s(&file_save_curvature, save_curvature, "w") != 0) {
        printf("Enable to open");
        free(save_curvature);
        return false;
    }

    //Save tangent
    const char* extension_tangent = "_tangent.csv";
    new_path_length = strlen(save_path) + strlen(extension_tangent) + 1;
    char* save_tangent = malloc(new_path_length);
    strcpy_s(save_tangent, new_path_length, save_path);
    strcat_s(save_tangent, new_path_length, extension_tangent);

    FILE* file_save_tangent;
    if (fopen_s(&file_save_tangent, save_tangent, "w") != 0) {
        printf("Enable to open");
        free(save_curvature);
        free(save_tangent);
        return false;
    }

    //Save torsion
    const char* extension_torsion = "_torsion.csv";
    new_path_length = strlen(save_path) + strlen(extension_torsion) + 1;
    char* save_torsion = malloc(new_path_length);
    strcpy_s(save_torsion, new_path_length, save_path);
    strcat_s(save_torsion, new_path_length, extension_torsion);

    FILE* file_save_torsion;
    if (fopen_s(&file_save_torsion, save_torsion, "w") != 0) {
        printf("Enable to open");
        free(save_curvature);
        free(save_tangent);
        free(save_torsion);
        return false;
    }

    for (size_t i = 1; i < curve->numPoints - 1; i++)
    {
        Point3D currentPt = curve->pPoints[i].pt;
        Point3D prevPt = curve->pPoints[i - 1].pt;
        Point3D nextPt = curve->pPoints[i + 1].pt;
        double h1 = getPoint3DDistance(currentPt, prevPt);
        double h2 = getPoint3DDistance(currentPt, nextPt);
        double coef = 1.0 / (h1 + h2);

		//curvature vector
        double curv_x = 2.0 * coef * ((nextPt.x - currentPt.x) / h2 - (currentPt.x - prevPt.x) / h1);
        double curv_y = 2.0 * coef * ((nextPt.y - currentPt.y) / h2 - (currentPt.y - prevPt.y) / h1);
        double curv_z = 2.0 * coef * ((nextPt.z - currentPt.z) / h2 - (currentPt.z - prevPt.z) / h1);
        double curvature = sqrt(curv_x * curv_x + curv_y * curv_y + curv_z * curv_z);
        fprintf(file_save_curvature, "%d, %lf\n", i, curvature);

		//tangent vector
        double tan_x = coef * (nextPt.x - prevPt.x);
        double tan_y = coef * (nextPt.y - prevPt.y);
        double tan_z = coef * (nextPt.z - prevPt.z);
        double tan = sqrt(tan_x * tan_x + tan_y * tan_y + tan_z * tan_z);
        fprintf(file_save_tangent, "%d, %lf, %lf\n", i, tan_z, tan);
        
		//torsion vector
        if (i >= 2 && i < curve->numPoints - 2) 
        {
            Point3D i_c_minus_two = curve->pPoints[i - 2].pt;
            Point3D i_c_minus_one = curve->pPoints[i - 1].pt;
            Point3D i_c = curve->pPoints[i].pt;
            Point3D i_c_plus_one = curve->pPoints[i + 1].pt;
            Point3D i_c_plus_two = curve->pPoints[i + 2].pt;

            double e_i_x = (i_c_plus_one.x - i_c.x);
            double e_i_y = (i_c_plus_one.y - i_c.y);
            double e_i_z = (i_c_plus_one.z - i_c.z);

            double e_i_minus_1_x = (i_c.x - i_c_minus_one.x);
            double e_i_minus_1_y = (i_c.y - i_c_minus_one.y);
            double e_i_minus_1_z = (i_c.z - i_c_minus_one.z);

            double e_i_plus_x = (i_c_plus_two.x - i_c_plus_one.x);
            double e_i_plus_y = (i_c_plus_two.y - i_c_plus_one.y);
            double e_i_plus_z = (i_c_plus_two.z - i_c_plus_one.z);

            double torsion_numerator =
                e_i_minus_1_x * (e_i_y * e_i_plus_z - e_i_z * e_i_plus_y)
                + e_i_minus_1_y * (e_i_z * e_i_plus_x - e_i_x * e_i_plus_z)
                + e_i_minus_1_z * (e_i_x * e_i_plus_y - e_i_y * e_i_plus_x);

            double cros_x = e_i_minus_1_y * e_i_z - e_i_minus_1_z * e_i_y;
            double cros_y = e_i_minus_1_z * e_i_x - e_i_minus_1_x * e_i_z;
            double cros_z = e_i_minus_1_x * e_i_y - e_i_minus_1_y * e_i_x;

            double torsion_denominator = cros_x * cros_x + cros_y * cros_y + cros_z * cros_z;
            double torsion = torsion_numerator / torsion_denominator;

            fprintf(file_save_torsion, "%d, %lf\n", i, torsion);
        }
    }
    fclose(file_save_curvature);
    fclose(file_save_tangent);
    fclose(file_save_torsion);

	free(save_curvature);
	free(save_tangent);
	free(save_torsion);

    return true;
}