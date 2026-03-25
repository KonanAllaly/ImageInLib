#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "common_functions.h"
#include "../src/solvers.h"

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

	bool getClosestPointToCurve(LinkedCurve3D* initial_curve, LinkedPoint3D* evolving_curve, bool isTheFirstTimeStep, double* px, double* py, double* pz);

	bool smoothingByLagrangeanCurveEvolution(const LagrangeanSmoothingParameters* pSmoothingParameters,
		unsigned char* pOutputPathPtr, Curve3D* pResultCurve);

	bool normalVelocitySmoothing(LinkedCurve3D* initial_curve, LinkedCurve3D* evolving_curve, SchemeData3D* pscheme_data,
		bool isTheFirstTimeStep, const double delta, const double lambda);

	bool tangentialVelocitySmoothing(LinkedCurve3D* evolving_curve, SchemeData3D* pscheme_data, const double omega, bool isCurveOpen);

	bool coefficientsSmoothing(LinkedCurve3D* plinked_curve, SchemeData3D* pscheme_data, const double delta, const double tau);

	bool evolveForSmoothingBySingleStep(LinkedCurve3D* pinitial_curve, LinkedCurve3D* plinked_curve, SchemeData3D* pscheme_data,
		const LagrangeanSmoothingParameters* pparams, bool isFirstTimeStep);

#ifdef __cplusplus
}
#endif
