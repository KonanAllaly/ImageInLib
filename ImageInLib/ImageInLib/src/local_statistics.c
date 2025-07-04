#include "math.h"
#include "local_statistics.h"

Statistics statistics2D(Image_Data2D inputImageData, Point2D pointOfInterest, const dataType radiusOfSearchDomain, const seachDomainShape shape)
{
	if (inputImageData.imageDataPtr == NULL)
	{
		return (Statistics) { -1, -1, -1, -1 };
	}
		
	const size_t length = inputImageData.height;
	const size_t width = inputImageData.width;

	BoundingBox2D box = findBoundingBox2D(pointOfInterest, length, width, radiusOfSearchDomain);

	dataType min_data = 100000000;//large value
	dataType max_data = 0;
	dataType mean_data = 0.0, std = 0.0, sum = 0;
	size_t count = 0;

	if (shape == SQUARE) {
		//compute minimum, maximum and mean
		for (size_t i = box.x_min; i <= box.x_max; i++)
		{
			for (size_t j = box.y_min; j <= box.y_max; j++)
			{
				size_t xd = x_new(i, j, length);
				dataType value = inputImageData.imageDataPtr[xd];
				if(min_data > value)
				{
					min_data = value;
				}
				if (max_data < value)
				{
					max_data = value;
				}
				count++;
				mean_data += value;
			}
		}
		if (count != 0) {
			mean_data /= (dataType)count;
		}
		else {
			return (Statistics) { -1, -1, -1, -1 };
		}
		//compute standard deviation
		for (size_t i = box.x_min; i <= box.x_max; i++)
		{
			for (size_t j = box.y_min; j <= box.y_max; j++)
			{
				size_t xd = x_new(i, j, length);
				dataType value = inputImageData.imageDataPtr[xd];
				std += pow(value - mean_data, 2);
			}
		}
		std = sqrt(std / (dataType)count);
	}
	else {
		if(shape == SPHERE)
		{
			//compute minimum, maximum and mean
			for (size_t i = box.x_min; i <= box.x_max; i++)
			{
				for (size_t j = box.y_min; j <= box.y_max; j++)
				{
					Point2D current_point = { i, j };
					double pDistance = getPoint2DDistance(pointOfInterest, current_point);
					{
						if(pDistance <= radiusOfSearchDomain)
						{
							size_t xd = x_new(i, j, length);
							dataType value = inputImageData.imageDataPtr[xd];
							if (min_data > value)
							{
								min_data = value;
							}
							if (max_data < value)
							{
								max_data = value;
							}
							count++;
							mean_data += value;
						}
					}

				}
			}

			if (count != 0) {
				mean_data /= (dataType)count;
			}
			else {
				return (Statistics) { -1, -1, -1, -1 };
			}
			//compute standard deviation
			for (size_t i = box.x_min; i <= box.x_max; i++)
			{
				for (size_t j = box.y_min; j <= box.y_max; j++)
				{
					Point2D current_point = { i, j };
					double pDistance = getPoint2DDistance(pointOfInterest, current_point);
					if(pDistance <= radiusOfSearchDomain)
					{
						size_t xd = x_new(i, j, length);
						dataType value = inputImageData.imageDataPtr[xd];
						std += pow(value - mean_data, 2);
					}
				}
			}
			std = sqrt(std / (dataType)count);
		}
	}

	return (Statistics) { min_data, max_data, mean_data, std };
}

Statistics statistics3D(Image_Data inputImageData, Point3D pointOfInterest, const dataType radiusOfSearchDomain, const seachDomainShape shape)
{
	if (inputImageData.imageDataPtr == NULL)
	{
		return (Statistics) { -1, -1, -1, -1 };
	}

	return (Statistics) { -1, -1, -1, -1 };
}

void getStatisticsAroundPoint(void* pInputImageData, void* pointOfInterest, Statistics stats, const dataType radiusOfSearchDomain, const seachDomainShape shape, const pDimension dim)
{
	switch (dim)
	{
	case TWO_D:
	{
		Image_Data2D inputImageData = *(Image_Data2D*)pInputImageData;
		Point2D pStart = *(Point2D*)pointOfInterest;
		stats = statistics2D(inputImageData, pStart, radiusOfSearchDomain, shape);
		break;
	}
	case THREE_D:
	{
		Image_Data inputImageData = *(Image_Data*)pInputImageData;
		Point3D pStart = *(Point3D*)pointOfInterest;
		stats = statistics3D(inputImageData, pStart, radiusOfSearchDomain, shape);
		break;
	}
	default:
		break;
	}
}

