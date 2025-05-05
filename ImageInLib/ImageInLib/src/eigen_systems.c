#include <math.h>
#include "eigen_systems.h"
#include <common_math.h>

#define ROTATE(a,i,j,k,l) g=a[i][j]; h=a[k][l]; a[i][j]=g-s*(h+g*tau); a[k][l]=h+s*(g-h*tau);

//#define ROTATE(a, i, j, k, l, s, tau) \
//	do { \
//		dataType g_ = (a)[(i)][(j)]; \
//		dataType h_ = (a)[(k)][(l)]; \
//		(a)[(i)][(j)] = g_ - (s) * (h_ + g_ * (tau)); \
//		(a)[(k)][(l)] = h_ + (s) * (g_ - h_ * (tau)); \
//	} while (0)

//static inline void rotate(dataType** a, size_t i, size_t j, size_t k, size_t l, dataType s, dataType tau) {
//	dataType g = a[i][j];
//	dataType h = a[k][l];
//	a[i][j] = g - s * (h + g * tau);
//	a[k][l] = h + s * (g - h * tau);
//}

bool getHessianMatrix3D(dataType** imageDataPtr, dataType** hessianMatrix, const size_t length, const size_t width, const size_t height, const size_t ind_x, const size_t ind_y, const size_t ind_z, const VoxelSpacing fVolume)
{

	if (imageDataPtr == NULL || hessianMatrix == NULL || 
		fVolume.sx == 0 || fVolume.sy == 0 || fVolume.sz == 0 ||
		ind_x == 0 || ind_y == 0 || ind_z == 0 ||
		ind_x == length - 1 || ind_y == width - 1 || ind_z == height - 1)
	{
		return false;
	}

	size_t i = ind_x, j = ind_y, k = ind_z;

	size_t iplus = i + 1, iminus = i - 1;
	size_t jplus = j + 1, jminus = j - 1;
	size_t kplus = k + 1, kminus = k - 1;

	size_t xd = x_new(i, j, length);
	dataType u = imageDataPtr[k][xd];

	// Hessian Matrix
	// x direction
	hessianMatrix[0][0] = (imageDataPtr[k][x_new(iplus, j, length)] - 2 * u + imageDataPtr[k][x_new(iminus, j, length)]) 
		/ (fVolume.sx * fVolume.sx);

	hessianMatrix[0][1] = (imageDataPtr[k][x_new(iplus, jplus, length)] + imageDataPtr[k][x_new(iminus, jminus, length)]
		- imageDataPtr[k][x_new(iplus, jminus, length)] - imageDataPtr[k][x_new(iminus, jplus, length)]) 
		/ (4 * fVolume.sx * fVolume.sy);

	hessianMatrix[0][2] = (imageDataPtr[kplus][x_new(iplus, j, length)] + imageDataPtr[kminus][x_new(iminus, j, length)]
		- imageDataPtr[kminus][x_new(iplus, j, length)] - imageDataPtr[kplus][x_new(iminus, j, length)])
		/ (4 * fVolume.sx * fVolume.sz);
	
	// y direction
	hessianMatrix[1][1] = (imageDataPtr[k][x_new(i, jplus, length)] - 2 * u + imageDataPtr[k][x_new(i, jminus, length)]) 
		/ (fVolume.sy * fVolume.sy);

	hessianMatrix[1][0] = hessianMatrix[0][1];

	hessianMatrix[1][2] = (imageDataPtr[kplus][x_new(i, jplus, length)] + imageDataPtr[kminus][x_new(i, jminus, length)]
		- imageDataPtr[kplus][x_new(i, jminus, length)] - imageDataPtr[kminus][x_new(i, jplus, length)])
		/ (4 * fVolume.sy * fVolume.sz);
				
	// z direction						
	hessianMatrix[2][2] = (imageDataPtr[kplus][xd] - 2 * u + imageDataPtr[kminus][xd]) / (fVolume.sz * fVolume.sz);

	hessianMatrix[2][1] = hessianMatrix[1][2];

	hessianMatrix[2][0] = hessianMatrix[0][2];

	return true;
}

bool eigenSystems(dataType ** a, size_t n, dataType d[], dataType ** v, size_t *nrot) {
	
	if (a == NULL || v == NULL || d == NULL) {
		return false;
	}

	size_t j, iq, ip, i, it;
	dataType tresh, theta, tau, t, sm, s, h, g, c;
	
	dataType* b = (dataType*)malloc(n * sizeof(dataType));
	dataType* z = (dataType*)malloc(n * sizeof(dataType));
	if (b == NULL || z == NULL) {
		free(b);
		free(z);
		return false;
	}

	//Initialize identity matrix, corresponding to the eigenvectors
	for (ip = 0; ip < n; ip++) {
		for (iq = 0; iq < n; iq++) {
			v[ip][iq] = 0.0;
		}
		v[ip][ip] = 1.0;
	}

	//Initialize the diagonal elements of the matrix
	for (ip = 0; ip < n; ip++) {
		b[ip] = d[ip] = a[ip][ip];
		z[ip] = 0.0;
	}
	
	*nrot = 0;
	for (i = 0; i < 50; i++) {
		
		//Sum off-diagonal elements
		sm = 0.0;
		for (ip = 0; ip < n - 1; ip++) {
			for (iq = ip + 1; iq < n; iq++) {
				sm += fabs(a[ip][iq]);
			}
		}
		
		//The normal return, which relies on quadratic convergence to machine underflow 
		// is to check if the sum of off-diagonal elements is zero 
		if (sm == 0.0) {
			free(b);
			free(z);
			return true;
		}

		//On the first three sweeps
		// we use a threshold to skip the rotation if the off-diagonal element is small
		if (i < 3) {
			tresh = 1e-12;//0.2 * sm / (n * n);
		}
		else {
			tresh = 0.0; // ...thereafter  (there is no threshold)
		}

		//Sweep over the matrix, looking for nonzero off-diagonal elements
		for (ip = 0; ip < n - 1; ip++) {
			for (iq = ip + 1; iq < n; iq++) {
				
				g = 100.0 * fabs(a[ip][iq]);
				
				//After four sweeps, skip the rotation if the off-diagonal element is small
				if (i > 3 && (fabs(d[ip]) + g) == fabs(d[ip]) && (fabs(d[iq]) + g) == fabs(d[iq])) {
					a[ip][iq] = 0.0;
				}
				else if (fabs(a[ip][iq]) > tresh) {
					h = d[iq] - d[ip];
					if ((dataType)(fabs(h) + g) == (dataType)(fabs(h))) {
						t = a[ip][iq] / h;
					}
					else {
						theta = 0.5 * h / a[ip][iq];
						t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
					}
					c = 1.0 / sqrt(1 + t * t);
					s = t * c;
					tau = s / (1.0 + c);
					h = t * a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq] = 0.0;

					//case of rotation 1 <= j < p
					for (j = 0; j < ip - 1; j++) {
						ROTATE(a, j, ip, j, iq);
						//ROTATE(a, j, ip, j, iq, s, tau);
						//rotate(a, j, ip, j, iq, s, tau);
					}
					//case of rotation p < j < q
					for (j = ip + 1; j < iq - 1; j++) {
						ROTATE(a, ip, j, j, iq);
						//ROTATE(a, ip, j, j, iq, s, tau);
						//rotate(a, ip, j, j, iq, s, tau);
					}
					//case of rotation q < j <= n
					for (j = iq + 1; j < n; j++) {
						ROTATE(a, ip, j, iq, j);
						//ROTATE(a, ip, j, iq, j, s, tau);
						//rotate(a, ip, j, iq, j, s, tau);
					}
					for (j = 0; j < n; j++) {
						ROTATE(v, j, ip, j, iq);
						//ROTATE(v, j, ip, j, iq, s, tau);
						//rotate(v, j, ip, j, iq, s, tau);
					}
					++(*nrot);
				}
			}
		}
		
		//Update d with the new eigenvalues and reset the off-diagonal elements
		for (it = 0; it < n; it++) {
			b[it] += z[it];
			d[it] = b[it];
			z[it] = 0.0;
		}
		//nrerror("eigenSystems", "Eigenvalue calculation failed");
	}

	free(b);
	free(z);
	return true;
}

bool sortEigenValues(dataType* eigen_values, const size_t dim) {
	
	if (eigen_values == NULL || dim < 1) {
		return false;
	}

	size_t k, i, j;
	for(k = 0; k < dim; k++){
		for (i = 0; i < dim - 1; i++) {
			if (fabs(eigen_values[i]) > fabs(eigen_values[i + 1])) {
				dataType temp = eigen_values[i];
				eigen_values[i] = eigen_values[i + 1];
				eigen_values[i + 1] = temp;
			}
		}
	}

	return true;
}

bool generateGaussianKernel(dataType** kernel, const size_t kernelSize, const dataType sigma) {

	if (kernel == NULL || kernelSize < 1 || sigma <= 0) {
		return false;
	}

	size_t i, j, k, xd, dim = 3;
	dataType sum = 0.0;
	for (k = 0; k < kernelSize; k++) {
		for (i = 0; i < kernelSize; i++) {
			for (j = 0; j < kernelSize; j++) {
				xd = x_new(i, j, kernelSize);
				kernel[k][xd] = exp(-(pow(i, 2) + pow(j, 2) + pow(k, 2)) / (2.0 * sigma * sigma));
				sum += kernel[k][xd];
			}
		}
	}
	
	//The empirical normalization is considered to avoid accumulated errors
	// Analitical normalization ---> 1.0 / pow(sqrt(2.0 * M_PI * sigma * sigma), d);  d is the dimension
	for (k = 0; k < kernelSize; k++) {
		for (i = 0; i < kernelSize * kernelSize; i++) {
			kernel[k][i] /= sum;
		}
	}
	return true;
}

bool gaussianFiltering(dataType** imageDataPtr, dataType** filteredImageDataPtr, const size_t length, const size_t width, const size_t height, const dataType sigma) {
	
	if (imageDataPtr == NULL || filteredImageDataPtr == NULL || sigma <= 0) {
		return false;
	}
	size_t i, j, k, xd;
	const size_t kernelSize = 2 * (size_t)sigma + 1;
	const size_t p = (size_t)sigma;
	if (kernelSize < 3) {
		return false;
	}

	dataType** kernel = (dataType**)malloc(kernelSize * sizeof(dataType*));
	for (k = 0; k < kernelSize; k++) {
		kernel[k] = (dataType*)malloc(kernelSize * kernelSize * sizeof(dataType));
		if (kernel[k] == NULL) {
			for (size_t m = 0; m < k; m++) {
				free(kernel[m]);
			}
			free(kernel);
			return false;
		}
	}
	if (kernel == NULL) {
		return false;
	}

	generateGaussianKernel(kernel, kernelSize, sigma);

	const size_t pLength = length + 2 * p;
	const size_t pWidth = width + 2 * p;
	const size_t pHeight = height + 2 * p;
	dataType** imageDataExt = (dataType**)malloc(pHeight * sizeof(dataType*));
	for (k = 0; k < pHeight; k++) {
		imageDataExt[k] = (dataType*)malloc(pLength * pWidth * sizeof(dataType));
	}

	//Copy to extended area
	size_t i_ext, j_ext, k_ext, xd_ext;
	for (k = 0, k_ext = p; k < height; k++, k_ext++) {
		for (i = 0, i_ext = p; i < length; i++, i_ext++) {
			for (j = 0, j_ext = p; j < width; j++, j_ext++) {
				xd = x_new(i, j, length);
				xd_ext = x_new(i_ext, j_ext, pLength);
				imageDataExt[k_ext][xd_ext] = imageDataPtr[k][xd];
			}
		}
	}
	reflection3DB(imageDataExt, height, length, width, p);

	//Convolution
	size_t im, jm, km, xdm, xdn;
	for (k = 0, k_ext = p; k < height; k++, k_ext++) {
		for (i = 0, i_ext = p; i < length; i++, i_ext++) {
			for (j = 0, j_ext = p; j < width; j++, j_ext++) {
				xd = x_new(i, j, length);
				xd_ext = x_new(i_ext, j_ext, pLength);
				for (size_t kk = 0; kk < kernelSize; kk++) {
					for (size_t ki = 0; ki < kernelSize; ki++) {
						for (size_t kj = 0; kj < kernelSize; kj++) {
							km = k_ext + kk - p;
							im = i_ext + ki - p;
							jm = j_ext + kj - p;
							xdm = x_new(im, jm, pLength);
							xdn = x_new(ki, kj, kernelSize);
							filteredImageDataPtr[k][xd] += kernel[kk][xdn] * imageDataExt[km][xdm];
						}
					}
				}
			}
		}
	}
	
	for (k = 0; k < kernelSize; k++) {
		free(kernel[k]);
	}
	free(kernel);
	return true;
}

bool multiscaleFiltering(Image_Data inputImageData, dataType** vesselnessImage, dataType sigma) {

	if (inputImageData.imageDataPtr == NULL || vesselnessImage == NULL) {
		return false;
	}

	size_t k, i, j, xd, k_ext, i_ext, j_ext, xd_ext;
	size_t length = inputImageData.length;
	size_t width = inputImageData.width;
	size_t height = inputImageData.height;
	VoxelSpacing spacing = inputImageData.spacing;

	dataType** smootedhImage = (dataType**)malloc(height * sizeof(dataType*));
	for (k = 0; k < height; k++) {
		smootedhImage[k] = (dataType*)malloc(length * width * sizeof(dataType));
	}

	gaussianFiltering(inputImageData.imageDataPtr, smootedhImage, length, width, height, sigma);

	size_t dim = 3;
	dataType* eigen_values = (dataType*)malloc(dim * sizeof(dataType));
	dataType** eigen_vectors = (dataType**)malloc(dim * sizeof(dataType));
	for (k = 0; k < dim; k++) {
		eigen_vectors[k] = (dataType*)malloc(dim * sizeof(dataType));
	}
	size_t max_sweeps = 0;

	size_t length_ext = length + 2;
	size_t width_ext = width + 2;
	size_t height_ext = height + 2;

	if (length_ext < 3 || width_ext < 3 || height_ext < 3) {
		return false;
	}
	dataType** imageDataExt = (dataType**)malloc(height_ext * sizeof(dataType*));
	for (k = 0; k < height_ext; k++) {
		imageDataExt[k] = (dataType*)malloc(length_ext * width_ext * sizeof(dataType));
	}

	dataType** hessian_matrix = (dataType**)malloc(height * sizeof(dataType*));
	for (k = 0; k < height; k++) {
		hessian_matrix[k] = (dataType*)malloc(length * width * sizeof(dataType));
	}

	copyDataToExtendedArea(smootedhImage, imageDataExt, height, length, width);
	reflection3D(imageDataExt, height_ext, length_ext, width_ext);

	dataType rb, ra, s, c;
	dataType alpha2 = 0.25, beta2 = 0.25;

	size_t it, jt;

	for (k = 0, k_ext = 1; k < height; k++, k_ext++) {
		for (i = 0, i_ext = 1; i < length; i++, i_ext++) {
			for (j = 0, j_ext = 1; j < width; j++, j_ext++) {
				xd = x_new(i, j, length);

				//initialize arrays
				for (it = 0; it < dim; it++) {
					eigen_values[it] = 0.0;
					for (jt = 0; jt < dim; jt++) {
						hessian_matrix[it][jt] = 0.0;
						eigen_vectors[it][jt] = 0.0;
					}
				}

				getHessianMatrix3D(imageDataExt, hessian_matrix, length_ext, width_ext, height_ext, i_ext, j_ext, k_ext, spacing);
				eigenSystems(hessian_matrix, dim, eigen_values, eigen_vectors, &max_sweeps);
				sortEigenValues(eigen_values, dim);

				rb = fabs(eigen_values[0]) / sqrt(fabs(eigen_values[1]) * fabs(eigen_values[2]));
				ra = fabs(eigen_values[1]) / fabs(eigen_values[2]);
				s = pow(fabs(eigen_values[0]), 2) + pow(fabs(eigen_values[1]), 2) + pow(fabs(eigen_values[2]), 2);
				c = 0.5 * s;

				if (eigen_values[1] > 0 || eigen_values[2] > 0) {
					vesselnessImage[k][xd] = 0;
				}
				else {
					vesselnessImage[k][xd] = (1 - exp(-pow(ra, 2) / (2 * alpha2))) * exp(-pow(rb, 2) / (2 * beta2)) * (1 - exp(-pow(s, 2) / (2 * c * c)));
				}
			}
		}
	}

	free(eigen_values);
	for (k = 0; k < dim; k++) {
		free(hessian_matrix[k]);
		free(eigen_vectors[k]);
	}
	free(hessian_matrix);
	free(eigen_vectors);

	for (k = 0; k < height_ext; k++) {
		free(imageDataExt[k]);
	}
	free(imageDataExt);

	return true;
}

//=====================================================================================================

void reflection2DB(dataType* toReflectImage, size_t imageLength, size_t imageWidth, size_t p)
{
	size_t i, j;
	size_t length = imageLength, width = imageWidth; // Actual X, Y dimensions
	size_t rowLength = length + 2 * p;
	
	// Modified Dimension less 1
	size_t mlength = length - 1, mwidth = width - 1;
	
	// Y reflection
	for (j = p; j <= (mwidth) + p; j++)
	{
		toReflectImage[x_new(p - 1, j, rowLength)] = toReflectImage[x_new(p + 1, j, rowLength)];
		toReflectImage[x_new((mlength) + p + 1, j, rowLength)] = toReflectImage[x_new((mlength) + p - 1, j, rowLength)];
	}
	// X Direction
	for (i = 0; i <= (mlength) + 2 * p; i++)
	{
		toReflectImage[x_new(i, p - 1, rowLength)] = toReflectImage[x_new(i, p + 1, rowLength)];
		toReflectImage[x_new(i, (mwidth) + p + 1, rowLength)] = toReflectImage[x_new(i, (mwidth) + p - 1, rowLength)];
	}
}

bool generateGaussianKernel2D(dataType* kernel, const size_t kernelSize, const dataType sigma) {

	if (kernel == NULL || kernelSize < 1 || sigma <= 0) {
		return false;
	}

	dataType sum = 0.0;
	for (size_t i = 0; i < kernelSize; i++) {
		for (size_t j = 0; j < kernelSize; j++) {
			size_t xd = x_new(i, j, kernelSize);
			kernel[xd] = exp(-(pow(i, 2) + (j, 2)) / (2.0 * sigma * sigma));
			sum += kernel[xd];
		}
	}

	//The empirical normalization is considered to avoid accumulated errors
	// Analitical normalization ---> 1.0 / pow(sqrt(2.0 * M_PI * sigma * sigma), d);  d is the dimension
	for (size_t i = 0; i < kernelSize * kernelSize; i++) {
		kernel[i] /= sum;
	}
	return true;
}

bool gaussianSmoothing2D(dataType* imageDataPtr, dataType* smoothImageDataPtr, const size_t length, const size_t width, const dataType sigma) {

	if (imageDataPtr == NULL || smoothImageDataPtr == NULL || sigma <= 0) {
		return false;
	}

	size_t i, j, k, xd;
	const size_t kernelSize = 2 * (size_t)sigma + 1;
	const size_t p = (size_t)sigma;
	if (kernelSize < 3) {
		return false;
	}

	dataType* kernel = (dataType*)malloc(kernelSize * kernelSize * sizeof(dataType*));
	if (kernel == NULL) {
		return false;
	}

	//generate the Gaussian Kernel
	generateGaussianKernel2D(kernel, kernelSize, sigma);

	const size_t pLength = length + 2 * p;
	const size_t pWidth = width + 2 * p;
	dataType* imageDataExt = (dataType*)malloc(pLength * pWidth * sizeof(dataType));

	//Copy to extended area
	size_t i_ext, j_ext, k_ext, xd_ext;
	for (i = 0, i_ext = p; i < length; i++, i_ext++) {
		for (j = 0, j_ext = p; j < width; j++, j_ext++) {
			xd = x_new(i, j, length);
			xd_ext = x_new(i_ext, j_ext, pLength);
			imageDataExt[xd_ext] = imageDataPtr[xd];
		}
	}
	
	reflection2DB(imageDataExt, length, width, p);

	//Convolution
	size_t im, jm, km, xdm, xdn;
	for (i = 0, i_ext = p; i < length; i++, i_ext++) {
		for (j = 0, j_ext = p; j < width; j++, j_ext++) {
			xd = x_new(i, j, length);
			xd_ext = x_new(i_ext, j_ext, pLength);
			for (size_t ki = 0; ki < kernelSize; ki++) {
				for (size_t kj = 0; kj < kernelSize; kj++) {
					im = i_ext + ki - p;
					jm = j_ext + kj - p;
					xdm = x_new(im, jm, pLength);
					xdn = x_new(ki, kj, kernelSize);
					smoothImageDataPtr[xd] += kernel[xdn] * imageDataExt[xdm];
				}
			}
		}
	}

	free(kernel);
	
	return true;
}
