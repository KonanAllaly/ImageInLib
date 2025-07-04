
#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#ifndef COMMON_FUNCTIONS
#define COMMON_FUNCTIONS
	//==============================================================================
	//debug constants
#ifndef MEASURE_TIME
#define MEASURE_TIME
#endif
// Show/Hide Output in Console
#ifndef CONSOLE_OUTPUT
#define CONSOLE_OUTPUT
#endif
#undef CONSOLE_OUTPUT
//==============================================================================
// Typedefs
// Data Type

#ifndef USE_DOUBLE_PRECISION
	typedef float dataType;
#else
	typedef double dataType;
#endif // USE_DOUBLE_PRECISION

	//==============================================================================
	// Includes
#include <stddef.h>
#include <omp.h>
#include <stdbool.h>
//==============================================================================
// MACROs
//==============================================================================
// STRUCTs

	// Common 2D Points - {x,y}
	typedef struct ptstruct {
		dataType x;
		dataType y;
	} Point2D;

	// 2D point extended by flag indicating, if the point is end point (1st or last)
	typedef struct {
		union {
			struct ptstruct;
			Point2D pt;
		};
		bool isEndPoint;
		void* prev_point;
		void* next_point;
	} CurvePoint2D;

	// 2D Curve - list of CurvePoint2D
	typedef struct {
		CurvePoint2D* pPoints;
		size_t numPoints;
	} Curve2D;

	// LinkedCurve
	typedef struct LinkedPoint
	{
		struct LinkedPoint* next;
		struct LinkedPoint* previous;
		double x;
		double y;
		double distance_to_next;
		double distance_to_next_x;
		double distance_to_next_y;
		unsigned long long id;
	} LinkedPoint;

	typedef struct LinkedCurve
	{
		size_t number_of_points;
		LinkedPoint* first_point;
		double length;
	} LinkedCurve;

	// Common 3D Points - {x,y,z}
	typedef struct ptsstruct {
		dataType x;
		dataType y;
		dataType z;
	} Point3D;

	// 3D point extended by flag indicating, if the point is end point (1st or last)
	typedef struct {
		union {
			struct ptsstruct;
			Point3D pt;
		};
		bool isEndPoint;
		void* prev_point;
		void* next_point;
	} CurvePoint3D;

	// 3D Curve - list of CurvePoint3D
	typedef struct {
		CurvePoint3D* pPoints;
		size_t numPoints;
	} Curve3D;

	//Structure to handle image spacing
	typedef struct {
		dataType sx;
		dataType sy;
		dataType sz;
	} VoxelSpacing;

	//Matrix for rotation
	typedef struct {
		Point3D v1;
		Point3D v2;
		Point3D v3;
	}OrientationMatrix;

	typedef struct {
		dataType hx;
		dataType hy;
	} FiniteVolumeSize2D;

	typedef struct {
		dataType hx;
		dataType hy;
		dataType hz;
	} FiniteVolumeSize3D;

	// Image Container and Properties
	typedef struct {
		// Image Dimensions
		size_t height, length, width; // Absolute Dimension
		dataType** imageDataPtr; // Image Data Containers
		Point3D origin; // image origin
		VoxelSpacing spacing; // distance between pixels and distance between slice //--> voxel dimension
		OrientationMatrix orientation;
	} Image_Data;

	// Generate Random Points
	typedef struct {
		size_t k, xd, p;
	}Random3dPoints;

	typedef struct {
		size_t xd, p;
	}Random2dPoints;

	typedef struct {
		Point2D v1, v2;
	}OrientationMatrix2D;

	typedef struct {
		dataType sx, sy;
	} PixelSpacing;

	typedef struct {
		// Image Dimensions
		size_t height, width; // Absolute Dimension
		dataType* imageDataPtr; // Image Data Containers
		Point2D origin; // image origin
		PixelSpacing spacing; // pixel size
		OrientationMatrix2D orientation;
	} Image_Data2D;

	typedef struct {
		dataType min_data, max_data, mean_data, sd_data;
	} Statistics;

	typedef struct {
		size_t i_min, i_max;
		size_t j_min, j_max;
		size_t k_min, k_max;
	} BoundingBox;

	//==============================================================================
	// Shapes Container
	typedef struct {
		// Shapes
		dataType** shp;
		size_t num; // shape number
	} Shapes;
	//==============================================================================
	// Enumeration
	enum FiniteDifference
	{
		FINITE_FORWARD = 1, FINITE_BACKWARD, FINITE_CENTRAL
	};
	//==============================================================================
	// FUNCTION PROTOTYPES
	/*
	* Calculate to 1D representation from 2D Array
	* x from column and row indices
    #include <stddef.h> // Ensure this include is present at the top of the file  

    // Replace the problematic line with the following:  
    typedef size_t numPoints;
	*/
	size_t x_new(const size_t rowIndex, const size_t columnIndex, const size_t rowLength);
	//==============================================================================
	/* Flattens 3D array to 1D array
	*/
	size_t x_flat(const size_t rowIndex, const size_t columnIndex, const size_t heightIndex, const size_t rowLength, const size_t columnLength);
	//==============================================================================
	/*
	* Reflection Function Zuska
	* toReflectImage is the data to perform reflection on
	* imageHeight is the actual Z dimension
	* imageLength is the actual X dimension
	* imageWidth is the actual Y dimension
	* P is the boundary constant/distance
	* Reflection does not include the current point but only its adjacent neighbor
	* Preserves Neumann Boundary condition - Target is for those methods
	*/
	void reflection3DB(dataType** toReflectImage, size_t imageHeight, size_t imageLength, size_t imageWidth, size_t p);
	//==============================================================================
	/*
	* Reflection Function Jozef
	* Reflection includes the current point together with its adjacent neighbor
	* Targets Distance Function types
	*/
	void reflection3D(dataType** toReflectImage, size_t imageHeight, size_t imageLength, size_t imageWidth);
	//==============================================================================
	/*
	* Edge Detector Calculation function
	*/
	dataType edgeDetector(dataType value, dataType coef);
	//==============================================================================

	/// <summary>
	/// Evaluates the similar density detector (how much is the current value similar to reference value)
	/// </summary>
	/// <param name="currValue">current value</param>
	/// <param name="refValue">refference value</param>
	/// <param name="coef">coeficient</param>
	/// <returns></returns>
	dataType similarIntensityDetector(dataType currValue, dataType refValue, dataType coef);

	/*
	void copyDataToExtendedArea(const dataType ** originalDataPtr, dataType ** extendedDataPtr, const size_t originalHeight, const size_t originalLength, const size_t originalWidth);

	The function copies data from source given by originalDataPtr
	representing data of dimensions given by originalHeight, originalLength and originalWidth
	to an array extendedDataPtr representing enlarged data (by 1 vx in each direction)
	of dimensions originalHeight + 2, originalLength + 2 and originalWidth + 2
	*/
	void copyDataToExtendedArea(dataType** originalDataPtr, dataType** extendedDataPtr, const size_t originalHeight, const size_t originalLength, const size_t originalWidth);
	//==============================================================================
	void copyDataToReducedArea(dataType** originalDataPtr, const dataType** extendedDataPtr, const size_t originalHeight, const size_t originalLength, const size_t originalWidth);
	//==============================================================================
	// Copy from one pointer to another pointer
	void copyDataToAnotherArray(dataType** source, dataType** destination, size_t height, size_t length, size_t width);
	//==============================================================================
	void rescaleNewRange(dataType** imageDataPtr, size_t imageLength, size_t imageWidth, size_t imageHeight, dataType minNew, dataType maxNew, dataType max_dta, dataType min_dta);
	//==============================================================================
	void rescaleNewRange2D(dataType* imageDataPtr, size_t imageLength, size_t imageWidth, dataType minNew, dataType maxNew);
	//==============================================================================
	typedef struct {
		size_t k_min, i_min, j_min, k_max, i_max, j_max;
	} ClipBox;
	//==============================================================================
	// Calc. centroid of image data
	void centroidImage(dataType** imageDataPtr, dataType* centroid, size_t imageHeight, size_t imageLength, size_t imageWidth, dataType imageBackground);
	void centroidClipBox(dataType* centroid, ClipBox coord, dataType** imageDataPtr, size_t imageLength, dataType imageBackground);
	//==============================================================================
	void copyDataToAnother2dArray(dataType* source, dataType* destination, size_t imageHeight, size_t imageWidth);
	//==============================================================================
	void copyDataTo2dExtendedArea(dataType* originalDataPtr, dataType* extendedDataPtr, const size_t originalHeight, const size_t originalWidth);
	//==============================================================================
	void copyDataTo2dReducedArea(dataType* originalDataPtr, const dataType* extendedDataPtr, const size_t originalHeight, const size_t originalWidth);
	//==============================================================================
	void reflection2D(dataType* toReflectImage, size_t imageHeight, size_t imageWidth);
	//==============================================================================
	void reflection2DB(dataType* toReflectImage, size_t imageLength, size_t imageWidth, size_t p);
	//==============================================================================
	double getPoint2DDistance(const Point2D a, const Point2D b);
	//==============================================================================
	/// <summary>
	/// The points of pCurve are copied to pArray. The function expects same length of the particular objects (curve and array)
	/// </summary>
	/// <param name="pCurve">Pointer to the curve to be copied</param>
	/// <param name="array">Alocated array of expected size length (same as the curve) * 2</param>
	/// <returns></returns>
	bool copyCurve2DPointsToArray(const Curve2D* pCurve, dataType** pArray);

	/// <summary>
	/// The function check, if the given curve is closed or not
	/// </summary>
	/// <param name="pcurve">Input curve to be checked</param>
	/// <returns>Returns true, if the curve is closed, otherwise false.</returns>
	bool isCurveClosed(const Curve2D* pcurve);

	/// <summary>
	/// The function is implemented according to https://en.wikipedia.org/wiki/Polygon
	/// It returns indicates, if the given curve is oriented possitively.
	/// </summary>
	/// <param name="pcurve">pointer to given curve</param>
	/// <returns>true, if the curve given by pcurve is oriented positively (clock-wise) and false, otherwise (if the curve is unitialized as well)</returns>
	bool isCurveOrientedPositively(const Curve2D* pcurve);

	/// <summary>
	/// Returns approximation of curve center of gravity
	/// The function expects, that the curve points are distributed almost uniformly along the curve
	/// </summary>
	/// <param name="pcurve">a pointer to the input curve</param>
	/// <returns>calculted centroid</returns>
	Point2D getCurveCentroid(const Curve2D* pcurve);

	/// <summary>
	/// Calculates gradient in 2D point given by central difference on input data
	/// </summary>
	/// <param name="pbase_data">base (e.g. image) data</param>
	/// <param name="width">base data width</param>
	/// <param name="height">base data height</param>
	/// <param name="ind_x">x coordinate of the finite volume to calculate the gradient component</param>
	/// <param name="ind_y">y coordinate of the finite volume to calculate the gradient component</param>
	/// <param name="sz">size of finite volumes</param>
	/// <param name="grad">output - calculated gradient</param>
	/// <returns>True, if it was possible to estimate gradient</returns>
	bool getGradient2D(dataType* pbase_data, const size_t width, const size_t height, const size_t ind_x, const size_t ind_y, const FiniteVolumeSize2D sz, Point2D* grad);

	/// <summary>
	/// The function returns the distance to given point from orgin (0,0) - in other words, calculated a norm of the given vector 
	/// </summary>
	/// <param name="pt">Given input point</param>
	/// <returns>Returns the result of (sqrt(pt.x * pt.x + pt.y * pt.y))</returns>
	dataType norm(const Point2D pt);

	/// <summary>
	/// Compute the euclidian distance between two given 3D points
	/// </summary>
	/// <param name="a">first point</param>
	/// <param name="b">second point</param>
	/// <returns>distance between points</returns>
	double getPoint3DDistance(const Point3D a, const Point3D b);

	/// <summary>
	/// Find the point/voxel with the largest value in an array/image
	/// </summary>
	/// <param name="arrayPtr"></param>
	/// <param name="length">length of the array</param>
	/// <param name="width">width if the array</param>
	/// <param name="height">height of the array</param>
	/// <returns>return the corresponding 3D point</returns>
	Point3D getPointWithTheHighestValue(dataType** arrayPtr, const size_t length, const size_t width, const size_t height);

	/// <summary>
	/// The function returns the signum of the given value parameter represented
	/// by -1 for negative sign, 0 for zero value, 1 for possitive sign
	/// </summary>
	/// <param name="value">the value to be investigated</param>
	/// <returns></returns>
	double signum(const double value);

	// Functions for LinkedCurve management

	/// <summary>
	/// the function kjust creates empty LinkedCurve
	/// </summary>
	/// <returns>Returns the entity of empty Linked|Curve</returns>
	LinkedCurve createLinkedCurve();

	/// <summary>
	/// Creates a linked point given by coordinates, without any neighbours
	/// </summary>
	/// <param name="point_x">input x coordinate</param>
	/// <param name="point_y">input y coordinate</param>
	/// <returns>a pointer to created LinkedPoint</returns>
	LinkedPoint* createLinkedPoint(const double point_x, const double point_y);

	/// <summary>
	/// Pushes a point by the coordinates x and y immediately after the given linked point 
	/// </summary>
	/// <param name="linked_curve">given linked curve</param>
	/// <param name="linked_point">given linked point of the linked_curve (the function pushes a new point immediately after this point)</param>
	/// <param name="point_x">x coordinate of the point to be pushed after linked_point</param>
	/// <param name="point_y">y coordinate of the point to be pushed after linked_point</param>
	/// <returns>the pointer to newlu pushed LinkedPoint</returns>
	LinkedPoint* pushAfterPoint(LinkedCurve* linked_curve, LinkedPoint* linked_point, const double point_x, const double point_y);

	/// <summary>
	/// Releases memory reserved for the points of given linked_curve 
	/// </summary>
	/// <param name="linked_curve">pointer to given LinkedCurve to be released</param>
	void releaseLinkedCurve(LinkedCurve* linked_curve);		

	/// <summary>
	/// Initializes plinked_curve by data given by pcurve 
	/// </summary>
	/// <param name="pcurve">pointer to input Curve2D - serving as the base, plinked_curve is initialized based on</param>
	/// <param name="plinked_curve">pointer to LinkedCurve to be initialized</param>
	/// <param name="reverse">flag indicating, if the curve should be initialized in reverse order</param>
	/// <param name="close_curve">flag indicating, if the curve should be cloesed or open</param>
	/// <returns>Returns true, if the initialisation succeeds, false otherwide</returns>
	bool initializeLinkedCurve(Curve2D* pcurve, LinkedCurve* plinked_curve, const bool reverse, const bool close_curve);
	
	/// <summary>
	/// Updates the distance to the next linked point and the lenght of the whole curve as well
	/// </summary>
	/// <param name="linked_curve">pointer to given LionkedCurve</param>
	/// <param name="linked_point">pointer to given LinkedPoint</param>
	/// <returns>returns the new length </returns>
	double updateDistanceToNext(LinkedCurve* linked_curve, LinkedPoint* linked_point);

	/// <summary>
	/// update all local and global lengths of the given curve
	/// </summary>
	/// <param name="linked_curve">pointer to given curve</param>
	/// <returns>returns true, if the lenghts were upadted successfully, false otherwise</returns>
	bool updateLinkedCurveLengths(LinkedCurve* linked_curve);
  
	/// <summary>
	/// Updates thegiven LinkedPoint and corresponding lenths
	/// </summary>
	/// <param name="linked_curve">pointer to the input curve</param>
	/// <param name="linked_point">poiner to the point to be updated</param>
	/// <param name="x">input new x coordinate value</param>
	/// <param name="y">input new y coordinate value</param>
	/// <returns></returns>
	bool updatePoint(LinkedCurve* linked_curve, LinkedPoint* linked_point, const double x, const double y);

	//id generator	
	void resetIDGenerator();
	unsigned long long getNextID();

	/// <summary>
	/// Compute the reference intensity
	/// </summary>
	/// <param name="pimage">Image data structure</param>
	/// <param name="point1">Input point or initial circle center</param>
	/// <param name="radius">Radius of the circle to consider around the first point</>
	/// <returns>Return the mean pixel value in small circle around the first input point</returns>
	dataType getReferenceIntensity(Image_Data2D pimage, Point2D point1, double radius);

	/// <summary>
	/// Calculates gradient in 3D point given by central difference on input data
	/// </summary>
	/// <param name="pbase_data">base (e.g. image) data</param>
	/// <param name="width">base data width</param>
	/// <param name="length">base data length</param>
	/// <param name="height">base data height</param>
	/// <param name="ind_x">x coordinate of the finite volume to calculate the gradient component</param>
	/// <param name="ind_y">y coordinate of the finite volume to calculate the gradient component</param>
	/// <param name="ind_z">z coordinate of the finite volume to calculate the gradient component</param>
	/// <param name="fVolume">size of finite volume</param>
	/// <param name="grad">output - calculated gradient</param>
	/// <returns>True, if it was possible to estimate gradient</returns>
	bool getGradient3D(dataType ** pbase_data, const size_t width, const size_t length, const size_t height, const size_t ind_x, const size_t ind_y, const size_t ind_z, const FiniteVolumeSize3D fVolume, Point3D * grad);

	/// <summary>
	/// The function returns the distance to given point from orgin (0,0) - in other words, calculated a norm of the given vector 
	/// </summary>
	/// <param name="pt">Given input point</param>
	/// <returns>Returns the result of (sqrt(pt.x * pt.x + pt.y * pt.y + pt.z * pt.z))</returns>
	dataType norm3D(const Point3D pt);

	/// <summary>
	/// Returns approximation of curve center of gravity
	/// The function expects, that the curve points are distributed almost uniformly along the curve
	/// </summary>
	/// <param name="pcurve">a pointer to the input curve</param>
	/// <returns>calculted centroid</returns>
	Point3D get3dCurveCentroid(const Curve3D* pcurve);

	/// <summary>
	/// The points of pCurve are copied to pArray. The function expects same length of the particular objects (curve and array)
	/// </summary>
	/// <param name="pCurve">Pointer to the curve to be copied</param>
	/// <param name="array">Alocated array of expected size length (same as the curve) * 3</param>
	/// <returns>return true if everything went well, false otherwise</returns>
	bool copyCurve3DPointsToArray(const Curve3D* pCurve, dataType** pArray);

	/// <summary>
	/// Compute the reference intensity
	/// </summary>
	/// <param name="pimage">Image data structure</param>
	/// <param name="point1">Input point or initial circle center</param>
	/// <param name="radius">Radius of the sphere to consider around the first point</>
	/// <returns>Return the mean pixel value in small circle around the first input point</returns>
	dataType getReferenceIntensity3D(Image_Data pimage, Point3D point, double radius);

	// LinkedCurve3D
	typedef struct LinkedPoint3D
	{
		struct LinkedPoint3D* next;
		struct LinkedPoint3D* previous;
		double x;
		double y;
		double z;
		double distance_to_next;
		double distance_to_next_x;
		double distance_to_next_y;
		double distance_to_next_z;
		unsigned long long id;
		//extend the 3D points structure with the normal velocity vector
		double nvx;
		double nvy;
		double nvz;
		double average_distance_to_next; //use to observe the average distance after evolution
	} LinkedPoint3D;

	typedef struct LinkedCurve3D
	{
		size_t number_of_points;
		LinkedPoint3D* first_point;
		double length;
	} LinkedCurve3D;

	/// <summary>
	/// Updates the distance to the next linked point and the lenght of the whole curve as well
	/// </summary>
	/// <param name="linked_curve">pointer to given linked curve</param>
	/// <param name="linked_point">pointer to given linked curve</param>
	/// <returns>returns the new length </returns>
	double updateDistance3dToNext(LinkedCurve3D* linked_curve, LinkedPoint3D* linked_point);

	/// <summary>
	/// Updates the distance to the next linked point and the lenght of the whole curve as well
	/// </summary>
	/// <param name="linked_curve">pointer to given linked curve</param>
	/// <param name="linked_point">pointer to given linked curve</param>
	/// <param name="x">point x coordinate</param>
	/// <param name="y">point y coordinate</param>
	/// <param name="z">point z coordinate</param>
	/// <returns>return true if everything went well, false otherwise</returns>
	bool update3dPoint(LinkedCurve3D* linked_curve, LinkedPoint3D* linked_point, const double x, const double y, const double z);

	typedef enum
	{
		TWO_D = 1,
		THREE_D
	} pDimension;

	/// <summary>
	/// function to create 3D linked curve
	/// </summary>
	/// <returns>return empty 3D linked curve</returns>
	LinkedCurve3D create3dLinkedCurve();

	/// <summary>
	/// check if the 3D curve is oriented
	/// </summary>
	/// <param name="pcurve">input curve</param>
	/// <returns>return true if everything went well, false otherwise</returns>
	bool is3dCurveOrientedPositively(const Curve3D* pcurve);

	/// <summary>
	/// Create linked point
	/// </summary>
	/// <param name="point_x">point x coordinate</param>
	/// <param name="point_y">point y coordinate</param>
	/// <param name="point_z">point z coordinate</param>
	/// <returns>return linked point</returns>
	LinkedPoint3D* create3dLinkedPoint(const double point_x, const double point_y, const double point_z);

	/// <summary>
	/// Pushes a point by the coordinates x, y and z after the given linked point
	/// </summary>
	/// <param name="linked_curve">linked curve</param>
	/// <param name="linked_point">linked point of the linked_curve (the function pushes a new point after this point)</param>
	/// <param name="point_x">point x coordinate</param>
	/// <param name="point_y">point y coordinate</param>
	/// <param name="point_z">point z coordinate</param>
	/// <returns>return 3D linked points</returns>
	LinkedPoint3D* pushAfter3dPoint(LinkedCurve3D* linked_curve, LinkedPoint3D* linked_point, const double point_x, const double point_y, const double point_z);

	/// <summary>
	/// Initialize given 3d linked curve (opened or closed)
	/// </summary>
	/// <param name="pcurve">the input curve, contains all initial curve points</param>
	/// <param name="plinked_curve">linked curve</param>
	/// <param name="reverse">check if the linked curve should be reverted or not</param>
	/// <param name="close_curve">check whether the curve is opened or closed</param>
	/// <returns>return true if everything went well, false otherwise</returns>
	bool initialize3dLinkedCurve(Curve3D* pcurve, LinkedCurve3D* plinked_curve, const bool reverse, const bool close_curve);

	/// <summary>
	/// release a given linked curve
	/// </summary>
	/// <param name="linked_curve">linked</param>
	void release3dLinkedCurve(LinkedCurve3D* linked_curve);

	BoundingBox findPointBoundingBox(Image_Data imageDataPtr, Point3D point, double radius);

#endif // !COMMON_FUNCTIONS

#ifdef __cplusplus
}
#endif