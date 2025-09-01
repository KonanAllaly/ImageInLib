#ifdef __cplusplus
extern "C" {
#endif

#pragma once
#include "common_functions.h"

void bubbleSort(dataType* pArray, const size_t height, const size_t width);

int partition(dataType* pArray, size_t low, size_t high);

void quickSort(dataType* pArray, size_t low, size_t high);

dataType computePercentile(dataType* pArray, const size_t height, const size_t width, const double perc);

#ifdef __cplusplus
}
#endif