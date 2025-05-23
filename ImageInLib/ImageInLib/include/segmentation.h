#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "common_functions.h"
#include "../src/segmentation3D_subsurf.h"
#include "segmentation2d.h"

	typedef enum
	{
		SUBSURF_MODEL = 1,
		GSUBSURF_MODEL,
		LABELING,
		GSUBSURF_ATLAS_MODEL,
		CURVE_2D_EXPLCIT,
		CURVE_2D_SEMI_IMPLICIT,
		CURVE_3D_EXPLICIT,
		CURVE_3D_SEMI_IMPLICIT
	} SegmentationMethod;

	// Structures and functions for 2D Langrangean segmentation

	// Structure that holds the parameters used during 2d Lagrangean segmentation process.
	typedef struct
	{
		size_t num_time_steps;	// Number of time steps
		Curve2D* pinitial_condition;	// Initial curve
		size_t num_points;// Number of initial segmentation curve points
		dataType time_step_size;	//discrete time step size
		dataType mu;	//influence of edge detector gradient field
		dataType lambda; //weight between projected gradient field and intensity similarity field
		dataType eps;	//influence of curvature
		dataType omega;	//redistribution speed
		dataType edgeCoef; // edge detector coefficient 
		dataType intensityCoef; // coeficient for the similar intensity detector
		dataType reference_intensity; //The reference intensity for G2
		void(*get_velocity)(Image_Data2D*, double, double, double*, double*);//pointer to the function returning the velocity for a given coordinate
		void(*get_g2)(Image_Data2D*, double, double, double, double, double*);//pointer to the function returning the g2 value for a given coordinate
		bool open_curve;
	} Lagrangean2DSegmentationParameters;

	// Structure that holds the parameters used during 3D Lagrangean segmentation process.
	typedef struct
	{
		size_t num_time_steps;	// Number of time steps
		Curve3D* pinitial_condition;	// Initial curve
		size_t num_points;// Number of initial segmentation curve points
		dataType time_step_size;	//discrete time step size
		dataType mu; //influence of edge detector gradient field
		dataType lambda; //weight between projected gradient field and intensity similarity field
		dataType eps; //influence of curvature
		dataType omega;	//redistribution speed
		dataType edge_detector_coef; //edge detector coefficient
		dataType intensity_coef; // coeficient for the similar intensity detector
		dataType reference_intensity; //The reference intensity for G2
		size_t n_save; //Frequence for saving the current curve
		dataType tolerance; //Stop the evolution when the total length is not increasing
		bool open_curve;
	} Lagrangean3DSegmentationParameters;

	/// <summary>
	/// The method segments the image according to given inoput paraneters
	/// </summary>
	/// <param name="pInputImageData">Pointer to Image_Data for 3D cases or pointer to Image_Data2D respectively representing an image to be segmented</param>
	/// <param name="pSegParameters">Pointer to Segmentation_Parameters representing segmentation parameters</param>
	/// <param name="pFilterParameters">Pointer to FilterParameters representing pre-filtering parameters</param>
	/// <param name="model">SegmentationMethod representing available segmentation methods</param>
	/// <param name="outputPathPtr">A pointer to unsigned char representing a destination path for resulting segmentation </param>
	/// <param name="pResultSegment">A pointer to Curve2D where results are returned. It is used by CURVE_2D_OPEN_EXPLCIT only</param>
	//TODO: utilize pResultSegment for another segmentation methods as well
	void segmentImage(void* pInputImageData, void* pSegParameters, void* pFilterParameters,
		const SegmentationMethod model, unsigned char* outputPathPtr, void* pResultSegment);

	void segment2dImage(Image_Data2D inputImageData, dataType* initialSegment, Segmentation_Parameters segParameters, FilterParameters filteringParameters,
		point2d* centers, const char* outputPathPtr, const SegmentationMethod model);

#ifdef __cplusplus
}
#endif