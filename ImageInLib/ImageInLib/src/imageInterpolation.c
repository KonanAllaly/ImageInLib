#include"../src/data_storage.h"
#include "interpolations.h"
#include "imageInterpolation.h"
#include<math.h>

bool nearestNeighborInterpolation(dataType** originalImage, dataType** newImage, size_t imageLength, size_t imageWidth, size_t imageHeight, dataType originalSpacing, dataType newSpacing)
{
    if (originalImage == NULL || newImage == NULL)
        return false;

    size_t i, j, k, kn, x;
    dataType k_int, k1, k2;

    k_int = 0; kn = 0;

    for (k = 0; k < imageHeight - 1; k++) {
        k1 = originalSpacing * k;
        k2 = originalSpacing * (k + 1);
        do {
            for (i = 0; i < imageLength; i++) {
                for (j = 0; j < imageWidth; j++) {
                    x = x_new(i, j, imageLength);
                    if (k_int == 0) {
                        newImage[kn][x] = originalImage[k][x];
                    }
                    else {
                        if ((k_int - k1) < (k2 - k_int)) {
                            newImage[kn][x] = originalImage[k][x];
                        }
                        else {
                            newImage[kn][x] = originalImage[k + 1][x];
                        }
                    }
                }
            }
            k_int = k_int + newSpacing;
            kn = kn + 1;
        } while (k_int < k2);
    }

    return true;
}

bool linear2dInterpolation(dataType** originalImage, dataType** newImage, size_t imageLength, size_t imageWidth, size_t imageHeight, dataType originalSpacing, dataType newSpacing)
{
    if (originalImage == NULL || newImage == NULL)
        return false;

    size_t i, j, k, x, kn;
    dataType k_int, k1, k2;
    dataType divisionByOriginalSpacing = 1.0 / originalSpacing;

    k_int = 0; kn = 0;
    for (k = 0; k < imageHeight - 1; k++) {
        k1 = k * originalSpacing;
        k2 = k1 + originalSpacing;
        do {
            for (i = 0; i < imageLength; i++) {
                for (j = 0; j < imageWidth; j++) {
                    x = x_new(i, j, imageLength);
                    newImage[kn][x] = (dataType)(originalImage[k][x] * ((k2 - k_int) * divisionByOriginalSpacing) +
                        originalImage[k + 1][x] * ((k_int - k1) * divisionByOriginalSpacing));
                }
            }
            k_int = k_int + newSpacing;
            kn = kn + 1;
        } while (k_int < k2);
    }

    return true;
}

bool downSampling(dataType** originalImage, dataType** newImage, size_t length, size_t width, size_t height) {

    if (originalImage == NULL || newImage == NULL)
        return false;

    size_t lengthNew = (length / 2), widthNew = (width / 2), heightNew = (height / 2);

    size_t i, j, k;

    for (k = 0; k < heightNew; k++) {
        for (i = 0; i < lengthNew; i++) {
            for (j = 0; j < widthNew; j++) {
                newImage[k][x_new(i, j, lengthNew)] = originalImage[k * 2][x_new(i * 2, j * 2, length)];
            }
        }
    }
    return true;
}

bool upSampling(dataType** originalImage, dataType** newImage, size_t length, size_t width, size_t height) {

    if (originalImage == NULL || newImage == NULL)
        return false;

    size_t i, j, k, in, jn, kn, indx;

    size_t lengthNew = length * 2, widthNew = width * 2, heightNew = height * 2;

    for (k = 0; k < height; k++) {
        for (i = 0; i < length; i++) {
            for (j = 0; j < width; j++) {

                indx = x_new(i, j, length);
                in = 2 * i; jn = j * 2; kn = k * 2;

                newImage[kn][x_new(in, jn, lengthNew)] = originalImage[k][indx];
                newImage[kn][x_new(in + 1, jn, lengthNew)] = originalImage[k][indx];
                newImage[kn][x_new(in, jn + 1, lengthNew)] = originalImage[k][indx];
                newImage[kn][x_new(in + 1, jn + 1, lengthNew)] = originalImage[k][indx];

                newImage[kn + 1][x_new(in, jn, lengthNew)] = originalImage[k][indx];
                newImage[kn + 1][x_new(in + 1, jn, lengthNew)] = originalImage[k][indx];
                newImage[kn + 1][x_new(in, jn + 1, lengthNew)] = originalImage[k][indx];
                newImage[kn + 1][x_new(in + 1, jn + 1, lengthNew)] = originalImage[k][indx];
            }
        }
    }

    return true;
}

//========================================================

//3D images

Point3D getImageCoordFromRealCoord3D(Point3D image_point, Point3D image_origin, VoxelSpacing image_spacing, OrientationMatrix orientation) {
    Point3D resultPoint;
    resultPoint.x = image_origin.x + image_spacing.sx * orientation.v1.x * image_point.x;
    resultPoint.y = image_origin.y + image_spacing.sy * orientation.v2.y * image_point.y;
    resultPoint.z = image_origin.z + image_spacing.sz * orientation.v3.z * image_point.z;
    return resultPoint;
}

Point3D getRealCoordFomImageCoord3D(Point3D real_point, Point3D real_origin, VoxelSpacing real_spacing, OrientationMatrix orientation) {
    Point3D resultPoint;
    resultPoint.x = (real_point.x - real_origin.x) / real_spacing.sx * orientation.v1.x;
    resultPoint.y = (real_point.y - real_origin.y) / real_spacing.sy * orientation.v2.y;
    resultPoint.z = (real_point.z - real_origin.z) / real_spacing.sz * orientation.v3.z;
    return resultPoint;
}

// ======================================================

bool resizeImage(Image_Data2D oldImage, Image_Data2D newImage) {

    if (oldImage.imageDataPtr == NULL || newImage.imageDataPtr == NULL)
        return false;

    const size_t height_old = oldImage.height, width_old = oldImage.width;
    const size_t height_new = newImage.height, width_new = newImage.width;

    int i, j;
    dataType i_new, j_new;

    int i_int, j_int;

    int i_floor, i_ceil;
    int j_floor, j_ceil;

    //Compute scale factor
    dataType sx = (dataType)height_old / height_new;
    dataType sy = (dataType)width_old / width_new;

    dataType val = 0.0;

    //Map back to original coordinates
    for (i = 0; i < height_new; i++) {
        for (j = 0; j < width_new; j++) {

            //Compute corresponding point to current pixel in old image
            i_new = i * sx;
            j_new = j * sy;

            //find neighbors
            i_floor = floor(i_new);
            if (ceil(i_new) <= height_old - 1) {
                i_ceil = ceil(i_new);
            }
            else {
                i_ceil = height_old - 1;
            }

            j_floor = floor(j_new);
            if (ceil(j_new) <= width_old - 1) {
                j_ceil = ceil(j_new);
            }
            else {
                j_ceil = width_old - 1;
            }

            i_int = (int)i_new;
            j_int = (int)j_new;

            //find the closest neighbor
            if ((i_floor == i_ceil) && (j_floor == j_ceil)) {
                val = oldImage.imageDataPtr[x_new(i_int, j_int, height_old)];
            }

            if ((i_floor == i_ceil) && (j_floor != j_ceil)) {

                if (abs(j - j_ceil) > abs(j - j_floor)) {
                    val = oldImage.imageDataPtr[x_new(i_int, j_floor, height_old)];
                }
                else {
                    val = oldImage.imageDataPtr[x_new(i_int, j_ceil, height_old)];
                }

            }

            if (i_floor != i_ceil && (j_floor == j_ceil)) {

                if (abs(i - i_ceil) > abs(i - i_floor)) {
                    val = oldImage.imageDataPtr[x_new(i_floor, j_int, height_old)];
                }
                else {
                    val = oldImage.imageDataPtr[x_new(i_ceil, j_int, height_old)];
                }

            }

            if (i_floor != i_ceil && (j_floor != j_ceil)) {

                dataType minDist = 100000000;

                //Top left neighbor
                dataType distTL = sqrt(pow(i - i_floor, 2) + pow(j - j_floor, 2));
                if (distTL < minDist) {
                    minDist = distTL;
                    val = oldImage.imageDataPtr[x_new(i_floor, j_floor, height_old)];
                }

                //Top right neighbor
                dataType distTR = sqrt(pow(i - i_ceil, 2) + pow(j - j_floor, 2));
                if (distTR < minDist) {
                    minDist = distTR;
                    val = oldImage.imageDataPtr[x_new(i_ceil, j_floor, height_old)];
                }

                //Bottom Left neighbor
                dataType distBL = sqrt(pow(i - i_floor, 2) + pow(j - j_ceil, 2));
                if (distBL < minDist) {
                    minDist = distBL;
                    val = oldImage.imageDataPtr[x_new(i_floor, j_ceil, height_old)];
                }

                //Bottom Right neighbor
                dataType distBR = sqrt(pow(i - i_ceil, 2) + pow(j - j_ceil, 2));
                if (distBR < minDist) {
                    minDist = distBR;
                    val = oldImage.imageDataPtr[x_new(i_ceil, j_ceil, height_old)];
                }
            }

            newImage.imageDataPtr[x_new(i, j, height_new)] = val;
        }
    }

    return true;
}

//=======================================================
//2D images

Point2D getRealCoordFromImageCoord2D(Point2D image_point, Point2D real_origin, PixelSpacing real_spacing, OrientationMatrix2D orientation) {
    Point2D result_point;
    result_point.x = real_origin.x + real_spacing.sx * orientation.v1.x * image_point.x;
    result_point.y = real_origin.y + real_spacing.sy * orientation.v2.y * image_point.y;
    return result_point;
}

Point2D getImageCoordFromRealCoord2D(Point2D real_point, Point2D real_origin, PixelSpacing real_spacing, OrientationMatrix2D orientation) {
    Point2D result_point;
    result_point.x = (real_point.x - real_origin.x) / real_spacing.sx * orientation.v1.x;
    result_point.y = (real_point.y - real_origin.y) / real_spacing.sy * orientation.v2.y;
    return result_point;
}

PointNeighbors2D getPointNeighbors2D(Point2D point, PixelSpacing spacing) {
    
    PointNeighbors2D neighbors;

    int i_floor = floor(point.x / spacing.sx), i_ceil = ceil(point.x/ spacing.sx);
    int j_floor = floor(point.y / spacing.sy), j_ceil = ceil(point.y/ spacing.sy);

    neighbors.top_left.x = i_floor * spacing.sx; neighbors.top_left.y = j_floor * spacing.sy;
    neighbors.bottom_left.x = i_floor * spacing.sx; neighbors.bottom_left.y = j_ceil * spacing.sy;
    neighbors.top_right.x = i_ceil * spacing.sx; neighbors.top_right.y = j_floor * spacing.sy;
    neighbors.bottom_right.x = i_ceil * spacing.sx; neighbors.bottom_left.y = j_ceil * spacing.sy;

    return neighbors;
}

dataType getDistance2D(Point2D point1, Point2D point2) {
    return sqrt(pow(point1.x - point2.x, 2) + pow(point1.y - point2.y, 2));
}

Point2D getNearestNeighbor2D(Point2D point, PixelSpacing spacing) {

    PointNeighbors2D neighbors;
    Point2D result_point;

    neighbors = getPointNeighbors2D(point, spacing);

    //Case 1 : 
    if (neighbors.top_left.x == neighbors.bottom_left.x && neighbors.top_left.y == neighbors.top_right.y) {
        result_point.x = point.x;
        result_point.y = point.y;
    }

    //Case 2 :
    if (neighbors.top_left.x == neighbors.top_right.x && neighbors.top_left.y != neighbors.bottom_left.y) {
        if (abs(point.y - neighbors.top_left.y) < abs(point.y - neighbors.bottom_left.y)) {
            result_point.x = point.x;
            result_point.y = neighbors.top_left.y;
        }
        else {
            result_point.x = point.x;
            result_point.y = neighbors.bottom_left.y;
        }
    }

    //Case 3 :
    if (neighbors.top_left.y == neighbors.bottom_left.y && neighbors.top_left.x != neighbors.top_right.x) {
        if (abs(point.x - neighbors.top_left.x) < abs(point.x - neighbors.top_right.x)) {
            result_point.x = neighbors.top_left.x;
            result_point.y = point.y;
        }
        else {
            result_point.x = neighbors.top_left.x;
            result_point.y = point.y;
        }
    }

    //Case 4 :
    if (neighbors.top_left.x != neighbors.top_right.x && neighbors.top_left.y != neighbors.bottom_left.y) {

        dataType distance[4];
        distance[0] = getDistance2D(point, neighbors.top_left);
        distance[1] = getDistance2D(point, neighbors.top_right);
        distance[2] = getDistance2D(point, neighbors.bottom_left);
        distance[3] = getDistance2D(point, neighbors.bottom_right);

        dataType min_dist = 1000000000000;
        int k, min_indice;

        for (k = 0; k < 4; k++) {
            if (distance[k] < min_dist) {
                min_dist = distance[k];
                min_indice = k;
            }
        }

        if (min_indice == 0) {
            result_point.x = neighbors.top_left.x;
            result_point.y = neighbors.top_left.y;
        }
        if (min_indice == 1) {
            result_point.x = neighbors.top_right.x;
            result_point.y = neighbors.top_right.y;
        }
        if (min_indice == 2) {
            result_point.x = neighbors.bottom_left.x;
            result_point.y = neighbors.bottom_left.y;
        }
        if (min_indice == 3) {
            result_point.x = neighbors.bottom_right.x;
            result_point.y = neighbors.bottom_right.y;
        }
    }

    return result_point;
}

/*
dataType getInterpolatedValueNearestNeighbor2D(Point2D point, PixelSpacing spacing) {

}
*/

bool interpolateImageFromImageCStoRealWorldCS(Image_Data2D src_image, Point2D real_origin, PixelSpacing real_spacing, OrientationMatrix2D real_orientation, Image_Data2D interp_image) {

    if (src_image.imageDataPtr == NULL || interp_image.imageDataPtr == NULL)
        return false;

    const size_t src_height = src_image.height, src_width = src_image.width;
    const size_t interp_height = interp_image.height, interp_width = interp_image.width;

    int i, j;
    Point2D src_image_point, interp_image_point, rw_point, interp_rw_point;

    PointNeighbors2D neighbor_points;
    OrientationMatrix2D orientation;

    //We need to recompute the scale factor from height and width because they were casted
    dataType scale_x = (dataType)src_height / interp_height;
    dataType scale_y = (dataType)src_width / interp_width;

    for (i = 0; i < interp_height; i++) {
        for (j = 0; j < interp_width; j++) {

            //Interpolated image point
            interp_image_point.x = i; 
            interp_image_point.y = j;

            //Get interpolated point in real world CS
            interp_rw_point = getRealCoordFromImageCoord2D(interp_image_point, real_origin, real_spacing, real_orientation);
           
            //Find the source point correponding to interpolated point in real world coordinates system
            rw_point.x = interp_rw_point.x * scale_x;
            rw_point.y = interp_rw_point.y * scale_y;

            //Find the nearest neighbor in real world coordinates system
            Point2D nearest_neighbor = getNearestNeighbor2D(rw_point, interp_image.spacing);

            //Find the corresponding point in image coordinates system
            interp_image_point = getImageCoordFromRealCoord2D(nearest_neighbor, real_origin, real_spacing, real_orientation);
            
            //TODO
            //getInterpolatedValueNearestNeighbor2D
            //interp_image.imageDataPtr[x_new(i, j, interp_height)] = src_image.imageDataPtr[x_new((int)int_point.x, (int)int_point.y, src_height)];

            //Set the interpolated value
            int i_int = (int)nearest_neighbor.x, j_int = (int)nearest_neighbor.y;

            if (i_int > src_height - 1) {
                i_int = src_height - 1;
            }
            if (j_int > src_width - 1) {
                j_int = src_width - 1;
            }

            interp_image.imageDataPtr[x_new(i, j, interp_height)] = src_image.imageDataPtr[x_new(i_int, j_int, src_height)];

        }
    }
    return true;
}

/*
bool resizeImageFromRealCoordToImageCoord(Image_Data2D src_image, Image_Data2D dest_image) {

    if (src_image.imageDataPtr == NULL || dest_image.imageDataPtr == NULL)
        return false;

    const size_t src_height = src_image.height, src_width = src_image.width;
    const size_t dest_height = dest_image.height, dest_width = dest_image.width;

    int i, j;
    Point2D src_point, dest_point, int_point;
    PointNeighbors2D neighbor_points;

    OrientationMatrix2D orientation;
    orientation.v1.x = src_image.orientation.v1.x * dest_image.orientation.v1.x + src_image.orientation.v2.x * dest_image.orientation.v1.y;
    orientation.v2.x = src_image.orientation.v1.x * dest_image.orientation.v2.x + src_image.orientation.v2.x * dest_image.orientation.v2.y;
    orientation.v1.y = src_image.orientation.v1.y * dest_image.orientation.v1.x + src_image.orientation.v2.y * dest_image.orientation.v1.y;
    orientation.v2.y = src_image.orientation.v1.y * dest_image.orientation.v2.x + src_image.orientation.v2.y * dest_image.orientation.v2.y;

    for (i = 0; i < dest_height; i++) {
        for (j = 0; j < dest_width; j++) {

            dest_point.x = i;
            dest_point.y = j;

            //find the corresponding point to the image CS points
            src_point = getRealCoordFromImageCoord2D(dest_point, src_image.origin, src_image.spacing, orientation);

            int_point = findNearestNeighbor2d(src_point, src_height, src_width);

            dest_image.imageDataPtr[x_new(i, j, dest_height)] = src_image.imageDataPtr[x_new((int)int_point.x, (int)int_point.y, src_height)];

        }
    }
    return true;
}
*/
