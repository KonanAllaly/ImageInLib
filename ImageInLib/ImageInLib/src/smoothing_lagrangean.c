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
    initialize3dLinkedCurve(pSmoothingParameters->pinitial_condition, &initial_curve, !isOrientedPositively, isCurveClosed);

    //create evolving linked curve
    LinkedCurve3D evolving_curve = create3dLinkedCurve();

    //Initialize the evolving linked curve
    initialize3dLinkedCurve(pSmoothingParameters->pinitial_condition, &evolving_curve, !isOrientedPositively, isCurveClosed);

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

    do
    {
        if (it == 1)
        {
            isFirstTimeStep = true;
        }
        else
        {
            isFirstTimeStep = false;
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

    } while (it < pSmoothingParameters->num_time_steps);

    LinkedPoint3D* final_curve = evolving_curve.first_point;

    for (size_t i = 0; i < evolving_curve.number_of_points; i++)
    {
        pResultCurve->pPoints[i].x = (dataType)final_curve->x;
        pResultCurve->pPoints[i].y = (dataType)final_curve->y;
        pResultCurve->pPoints[i].z = (dataType)final_curve->z;
        final_curve = final_curve->next;
    }

	final_curve = NULL;
    free(pscheme_data);
    release3dLinkedCurve(&evolving_curve);
    release3dLinkedCurve(&initial_curve);

    return true;
}

void normalVelocitySmoothing(LinkedCurve3D* initial_curve, LinkedCurve3D* evolving_curve, SchemeData3D* pscheme_data,
    const double delta, const double lambda, bool isFirstTimeStep)
{

    //check if the pointers a well allocated
    if (initial_curve == NULL || evolving_curve == NULL || pscheme_data == NULL)
    {
        return;
    }

	//No topological change, therefore the number of points is the same during the evolution
    const size_t number_of_points = evolving_curve->number_of_points;

    double h_i = -1, h_i_plus = -1;            // distance between two neighboring points
    double tx = 0, ty = 0, tz = 0;             //tangential vector components
    double nx = 0, ny = 0, nz = 0;             //tangential vector components
    double som_dist = 0;
    double curv_x = 0, curv_y = 0, curv_z = 0; //discrete curvature vector components

    double min_dist = 1e6;//evolving_curve->length;
    Point3D point_of_interest = { 0.0, 0.0, 0.0 };
    Point3D p1 = { 0.0, 0.0, 0.0 };
    Point3D p2 = { 0.0, 0.0, 0.0 };
    Point3D q = { 0.0, 0.0, 0.0 };
    double dist_min, dist, dist_new, t, t_cl, numerator, denominator, curvature, omega;

    LinkedPoint3D* evolving_point = evolving_curve->first_point;
    LinkedPoint3D* p_ref_initial = NULL;

    //Loop for the curve evolution
    for (size_t i = 1; i <= number_of_points; i++)
    { 
        if (i > 1 && i < number_of_points)
        {
            //Find the closest point on the initial curve  
            dist_min = 1e6;
            if (isFirstTimeStep == true)//This test should be moved up
            {
				//if true, we are at the first time step, 
                // therefore the point of interest is the same as the evolving point
                // which mean w = (x0 - x).N = 0
                point_of_interest.x = evolving_point->x;
                point_of_interest.y = evolving_point->y;
                point_of_interest.z = evolving_point->z;
            }
            else
            {   
                p_ref_initial = initial_curve->first_point->next;//r^{0}_{j}
                for (size_t j = 2; j <= number_of_points; j++)
                {
                    //Step 1 : compute t

                    //p1 = r^{0}_{j} - r^{0}_{j-1}
                    p1.x = p_ref_initial->x - p_ref_initial->previous->x;
                    p1.y = p_ref_initial->y - p_ref_initial->previous->y;
                    p1.z = p_ref_initial->z - p_ref_initial->previous->z;

                    //p2 = r^{n}_i - r^{0}_{j-1}
                    p2.x = evolving_point->x - p_ref_initial->previous->x;
                    p2.y = evolving_point->y - p_ref_initial->previous->y;
                    p2.z = evolving_point->z - p_ref_initial->previous->z;

                    numerator = p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;//scalar product : p1.p2
                    denominator = p1.x * p1.x + p1.y * p1.y + p1.z * p1.z;//scalar product : p1.p1
                    if (denominator != 0.0)
                    {
                        t = numerator / denominator;
                    }
                    else
                    {
						//For debugging, this should not happen 
                        // since we should not have two coincident 
                        // points in the initial curve
                        //printf("Warning: zero denominator\n");
                        t = 0.0;
                    }

					//If 0 <= t <= 1, the closest point is between r^{0}_{j-1} and r^{0}_{j}
					//If t < 0, the closest point is r^{0}_{j-1}
					//If t > 1, the closest point is r^{0}_{j}

                    //Step 2 : compute the closest point q
                    t_cl = fmax(0.0, fmin(1.0, t));
                    q.x = p_ref_initial->previous->x + t_cl * p1.x;
                    q.y = p_ref_initial->previous->y + t_cl * p1.y;
                    q.z = p_ref_initial->previous->z + t_cl * p1.z;

                    //Step 3 : compute the distance
                    dist = sqrt((evolving_point->x - q.x) * (evolving_point->x - q.x) +
                        (evolving_point->y - q.y) * (evolving_point->y - q.y) +
                        (evolving_point->z - q.z) * (evolving_point->z - q.z));

                    if (dist < dist_min)
                    {
                        dist_min = dist;
                        point_of_interest = q;
                    }

                    p_ref_initial = p_ref_initial->next;
                }
                
				p_ref_initial = NULL;
            }
            
            //Set the normal velocity towards the initial curve
            h_i = evolving_point->previous->distance_to_next;
            h_i_plus = evolving_point->distance_to_next;
            som_dist = h_i_plus + h_i;

            //Compute the discrete curvature vector components
            curv_x = (2.0 / som_dist) * (((evolving_point->next->x - evolving_point->x) / h_i_plus) - ((evolving_point->x - evolving_point->previous->x) / h_i));
            curv_y = (2.0 / som_dist) * (((evolving_point->next->y - evolving_point->y) / h_i_plus) - ((evolving_point->y - evolving_point->previous->y) / h_i));
            curv_z = (2.0 / som_dist) * (((evolving_point->next->z - evolving_point->z) / h_i_plus) - ((evolving_point->z - evolving_point->previous->z) / h_i));
            curvature = sqrt(curv_x * curv_x + curv_y * curv_y + curv_z * curv_z);

            //Set normal vector components
            if (curvature == 0)
            {
                evolving_point->nvx = 0.0;
                evolving_point->nvy = 0.0;
                evolving_point->nvz = 0.0;
            }
            else
            {
                evolving_point->nvx = curv_x / curvature;
                evolving_point->nvy = curv_y / curvature;
                evolving_point->nvz = curv_z / curvature;
            }

            pscheme_data[i].w = (point_of_interest.x - evolving_point->x) * evolving_point->nvx 
                + (point_of_interest.y - evolving_point->y) * evolving_point->nvy
                + (point_of_interest.z - evolving_point->z) * evolving_point->nvz;

            pscheme_data[i].beta = -delta * curvature + lambda * pscheme_data[i].w;
        }
        else
        {
            //when this happens, we are at the end points of an open curve
            evolving_point->nvx = 0.0;
            evolving_point->nvy = 0.0;
            evolving_point->nvz = 0.0;
            pscheme_data[i].w = 0.0;
            pscheme_data[i].beta = 0.0; 
        }
        evolving_point = evolving_point->next;
    }
    //fclose(file);
}

void tangentialVelocitySmoothing(LinkedCurve3D* evolving_curve, SchemeData3D* pscheme_data, const double omega, bool isCurveOpen)
{

    double mean = 0.0;
    const size_t number_of_points = evolving_curve->number_of_points;
    const double curve_length = evolving_curve->length;
    double h_i = -1;

	double avg_length = 1.0;
    if(isCurveClosed)
    {
        avg_length = curve_length / (double)(number_of_points - 1);//The curve is open
    }
    else
    {
        avg_length = curve_length / (double)(number_of_points);//The curve is closed
    }
    

    LinkedPoint3D* current_point = evolving_curve->first_point;
    for (size_t i = 1; i <= number_of_points; i++)
    {
        if (i > 1)
        {
            //if it is not the first point
            h_i = current_point->previous->distance_to_next;
        }
        else
        {
            //first point
            h_i = current_point->distance_to_next;
        }

        mean += h_i * pscheme_data[i].u * pscheme_data[i].k1;

        current_point = current_point->next;
    }

    mean /= curve_length;

    //it will therefore not move in the tangential direction
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

        pscheme_data[i].alfa = pscheme_data[i - 1].alfa + h_i * mean - h_i * pscheme_data[i].u * pscheme_data[i].k1 + omega * (avg_length - h_i);
        current_point = current_point->next;
    }

    pscheme_data[number_of_points].alfa = 0.0;
    pscheme_data[number_of_points + 1].alfa = 0.0;

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

    for (size_t iter = 1; iter <= evolving_curve->number_of_points; iter++)
    {

        if (is_curve_closed || (iter > 1 && iter < evolving_curve->number_of_points))
        {
            previous_point = current_point->previous;
            h_i = previous_point->distance_to_next;

            if (is_curve_closed || iter < evolving_curve->number_of_points)
            {
                h_i_plus = current_point->distance_to_next;
            }
            else
            {
                h_i_plus = h_i;
            }

            pscheme_data[iter].a = -delta / h_i - 0.5 * fmax(-pscheme_data[iter].alfa, 0);      //lower diagonal
            pscheme_data[iter].c = -delta / h_i_plus - 0.5 * fmax(pscheme_data[iter].alfa, 0); //upper diagonal
            pscheme_data[iter].m = (h_i_plus + h_i) / (2.0 * tau);
            pscheme_data[iter].b = pscheme_data[iter].m - (pscheme_data[iter].a + pscheme_data[iter].c);//diagonal
        }
        else
        {
            if (iter == 1)
            {
                h_i = current_point->distance_to_next;
                h_i_plus = h_i;
            }
            else
            {
                h_i = current_point->previous->distance_to_next;
                h_i_plus = h_i;
            }

            pscheme_data[iter].a = 0.0;
            pscheme_data[iter].c = 0.0;
            pscheme_data[iter].m = (h_i_plus + h_i) / (2.0 * tau);
            pscheme_data[iter].b = 1.0;
        }

        current_point = current_point->next;
    }

    return true;
}

bool evolveForSmoothingBySingleStep(LinkedCurve3D* initial_curve, LinkedCurve3D* evolving_curve, SchemeData3D* pscheme_data,
    const LagrangeanSmoothingParameters* pSmoothingParameters, bool isFirstTimeStep)
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

	bool is_curve_closed = pSmoothingParameters->open_curve;

    //function to compute the normal velocity
    normalVelocitySmoothing(initial_curve, evolving_curve, pscheme_data, delta, lambda, isFirstTimeStep);

    //function to compute the tangential velocity
    tangentialVelocitySmoothing(evolving_curve, pscheme_data, omega, is_curve_closed);

    if (!coefficientsSmoothing(evolving_curve, pscheme_data, delta, tau))
    {
        return false;
    }

    LinkedPoint3D* current_point = evolving_curve->first_point;

    //////////////////////    X component ///////////////////////////////////////////////////////////

    for (size_t i = 1; i <= evolving_curve->number_of_points; i++)
    {

        if (is_curve_closed && (i == 1 || i == evolving_curve->number_of_points))
        {
            pscheme_data[i].ps = current_point->x;
        }
        else
        {
            pscheme_data[i].ps = pscheme_data[i].m * current_point->x + lambda * pscheme_data[i].w * current_point->nvx
                - 0.5 * fmin(pscheme_data[i].alfa, 0) * (current_point->x - current_point->next->x)
                - 0.5 * fmin(-pscheme_data[i].alfa, 0) * (current_point->x - current_point->previous->x);
        }
        current_point = current_point->next;
    }

    if (is_curve_closed)
    {
        calculate_by_thomas3D(pscheme_data, evolving_curve->number_of_points);
    }
    else
    {
        sherman_morris3D(pscheme_data, evolving_curve->number_of_points);
    }

    /////////////////////    Y component   ///////////////////////////////////////////////////////////

    current_point = evolving_curve->first_point;
    for (size_t i = 1; i <= evolving_curve->number_of_points; i++)
    {

        if (is_curve_closed && (i == 1 || i == evolving_curve->number_of_points))
        {
            pscheme_data[i].ps = current_point->y;
        }
        else
        {
            pscheme_data[i].ps = pscheme_data[i].m * current_point->y + lambda * pscheme_data[i].w * current_point->nvy
                - 0.5 * fmin(pscheme_data[i].alfa, 0) * (current_point->y - current_point->next->y)
                - 0.5 * fmin(-pscheme_data[i].alfa, 0) * (current_point->y - current_point->previous->y);
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

    if (is_curve_closed)
    {
        calculate_by_thomas3D(pscheme_data, evolving_curve->number_of_points);
    }
    else
    {
        sherman_morris3D(pscheme_data, evolving_curve->number_of_points);
    }

    /////////////////////    Z component   ///////////////////////////////////////////////////////////

    current_point = evolving_curve->first_point;
    for (size_t i = 1; i <= evolving_curve->number_of_points; i++)
    {

        if (is_curve_closed && (i == 1 || i == evolving_curve->number_of_points))
        {
            pscheme_data[i].ps = current_point->z;
        }
        else
        {
            pscheme_data[i].ps = pscheme_data[i].m * current_point->z + lambda * pscheme_data[i].w * current_point->nvz
                - 0.5 * fmin(pscheme_data[i].alfa, 0) * (current_point->z - current_point->next->z)
                - 0.5 * fmin(-pscheme_data[i].alfa, 0) * (current_point->z - current_point->previous->z);
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

    if (is_curve_closed)
    {
        calculate_by_thomas3D(pscheme_data, evolving_curve->number_of_points);
    }
    else
    {
        sherman_morris3D(pscheme_data, evolving_curve->number_of_points);
    }

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