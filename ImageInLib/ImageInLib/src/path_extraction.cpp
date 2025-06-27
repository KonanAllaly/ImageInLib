#include "path_extraction"

bool shortestPath2d(Image_Data2D distanceFuncPtr, Point2D* seedPoints, vector<Point2D>& path_points, Path_Parameters parameters) {

	if (distanceFuncPtr.imageDataPtr == NULL || seedPoints == NULL)
		return false;

	const size_t height = distanceFuncPtr.height;
	const size_t width = distanceFuncPtr.width;
	dataType dist_min = 0.0;

	size_t i = (size_t)seedPoints[1].x;
	size_t j = (size_t)seedPoints[1].y;
	dataType x = seedPoints[1].x;
	dataType y = seedPoints[1].y;

	Point2D grad;
	dataType tau = parameters.tau;
	size_t count_iter = 0;

	do {

		getGradient2D(distanceFuncPtr.imageDataPtr, height, width, i, j, distanceFuncPtr.spacing, &grad);
		dataType gradNorm = sqrt(grad.x * grad.x + grad.y * grad.y);
		x -= tau * grad.x / gradNorm;
		y -= tau * grad.y / gradNorm;

		Point2D current = { x, y };
		dist_min = getPoint2DDistance(current, seedPoints[0]);
		path_points.push_back(current);

		i = (size_t)round(x);
		j = (size_t)round(y);

		count_iter++;
	} while (dist_min > parameters.tolerance && count_iter < parameters.max_iteration);

	return true;
}
