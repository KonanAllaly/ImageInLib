#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "common_functions.h"
//#include "../src/segmentation3D_subsurf.h"
//#include "segmentation2d.h"

	// Structure that holds the parameters used during 3D Lagrangean smoothing process.
	typedef struct
	{
		size_t num_time_steps;	        // Number of time steps
		Curve3D* pinitial_condition;	// Initial curve
		size_t num_points;              // Number of initial segmentation curve points
		dataType time_step_size;	    // discrete time step size
		dataType delta;                 // curvature weight
		dataType lambda;                // attraction weight
		dataType omega;                 // redistribution speed
		bool open_curve;
	} LagrangeanSmoothingParameters;

	//bool smoothingByLagrangeanCurveEvolution(const LagrangeanSmoothingParameters* pSmoothingParams, Curve3D* pResultSegmentation, unsigned char* pOutputPathPtr);

#ifdef __cplusplus
}
#endif
