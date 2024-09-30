#include <stdlib.h>
#include "common_math.h"
#include "common_functions.h"
#include "generate_2d_curves.h"
#include "generate_3d_curves.h"

size_t getNumberOfLatitudeDivision(const Point3D* pInitialPoints, const size_t initialPointsCount, const double pointsDistance) {
    if (initialPointsCount != 2)
    {
        return 0;
    }
    return (size_t)((abs(pInitialPoints[1].z - pInitialPoints[0].z) + 0.5) / pointsDistance);
}

size_t howManyPointsForSphereCurve(const Point3D* pInitialPoints, const size_t initialPointsCount, const double pointsDistance)
{
    if (initialPointsCount != 2)
    {
        return 0;
    }

    size_t numberOfPoints = 0;

    //Compute the number of latitude division
    size_t countLatitudeDivision = getNumberOfLatitudeDivision(pInitialPoints, initialPointsCount, pointsDistance);

    double radius = getPoint3DDistance(pInitialPoints[0], pInitialPoints[1]);

    //Compute the number of points per circle for each latitude
    double theta, radius_circle;
    for (size_t i = 0; i < countLatitudeDivision; i++)
    {
        theta = M_PI * i / countLatitudeDivision - M_PI / 2;
        radius_circle = radius * cos(theta);
        const double perimeter = getCirclePerimeter(radius_circle);
        numberOfPoints += (size_t)((perimeter / pointsDistance) + 0.5);
    }

    return numberOfPoints;
}

bool generateSphereCurve(Curve3D* pCurve, const Point3D* pInitialPoints, const size_t initialPointsCount, const double pointsDistance)
{
    if (pCurve == NULL || pInitialPoints == NULL)
    {
        return false;
    }

    //const size_t spherePointsCount = howManyPointsForSphereCurve(pInitialPoints, initialPointsCount, pointsDistance);
    //CurvePoint2D* pcircle_points = malloc(sizeof(CurvePoint2D) * circlePointsCount);
    //pcurve->pPoints = pcircle_points;
    //pcurve->numPoints = circlePointsCount;
    //const double radius = getPoint2DDistance(pinitial_points[0], pinitial_points[1]);

    //double anlgleStep = 2 * M_PI / (double)circlePointsCount;
    //Point2D center = pinitial_points[0];

    //for (size_t i = 0; i < circlePointsCount; i++)
    //{
    //    double phi = anlgleStep * (double)i;
    //    pcurve->pPoints[i].x = (dataType)(center.x + radius * cos(phi));
    //    pcurve->pPoints[i].y = (dataType)(center.y + radius * sin(phi));
    //    pcurve->pPoints[i].isEndPoint = false;
    //}

    /*
    //Compute the number of latitude division
    size_t countLatitudeDivision = getNumberOfLatitudeDivision(pInitialPoints, initialPointsCount, pointsDistance);
    //Compute number of longitude division
    size_t countLongitudeDivision = 10;

    double x, y, z, r_circle;
    double radius = getPoint3DDistance(pInitialPoints[0], pInitialPoints[1]);
    Point3D center = pInitialPoints[0];

    for (size_t i = 0; i < countLatitudeDivision; i++)
    {
        double theta = M_PI * i / countLatitudeDivision - M_PI / 2;
        //theta varies from - M_PI / 2 to M_PI / 2 (from South Pole to North Pole)

        z = radius * sin(theta);

        //Compute the radius of the circle at this latitude
        r_circle = radius * cos(theta);

        //Loop over longitude divisions
        for (size_t j = 0; j < countLongitudeDivision; j++)
        {
            double phi = 2 * M_PI * j / countLongitudeDivision;
            // phi varies from 0 to 2 * M_PI (full circle)

            x = r_circle * cos(phi);
            y = r_circle * sin(phi);

            size_t indexPoint = j + countLongitudeDivision * i;

            pCurve->pPoints[indexPoint].x = (dataType)(center.x + x);
            pCurve->pPoints[indexPoint].y = (dataType)(center.y + y);
            pCurve->pPoints[indexPoint].z = (dataType)(center.z + z);
            pCurve->pPoints[indexPoint].isEndPoint = false;
        }
    }
    */

    return true;
}

bool generate3DLineCurve(Curve3D* pCurve, const Point3D* pInitialPoints, const size_t initialPointsCount, const double pointsDistance)
{
	//TODO
	return true;
}

size_t howManyPointsFor3DLineCurve(const Point3D* pInitialPoints, const size_t initialPointsCount, const double pointsDistance)
{
    //TODO
    return 0;

    /*
    if (initialPointsCount != 2)
    {
        return 0;
    }

    const double length = getPoint3DDistance(pInitialPoints[0], pInitialPoints[1]);
    return (size_t)((length / pointsDistance) + 1.5);
    */
}