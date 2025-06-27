#pragma once

#include<iostream>
#include<vector>
#include "common_functions.h"
#include "../src/data_load.h"
#include <stdio.h>
#include <string.h>


	typedef struct {
		size_t x, y;
		dataType arrival;
	}pointFastMarching2D;

	/// <summary>
	/// Swaps the values of two 2D points represented by pointFastMarching2D structures.
	/// </summary>
	/// <param name="a">Pointer to the first 2D point to swap.</param>
	/// <param name="b">Pointer to the second 2D point to swap.</param>
	void swap2dPoints(pointFastMarching2D* a, pointFastMarching2D* b);

	/// <summary>
	/// Restores the heap property by moving the element at the specified position down in a 2D fast marching heap.
	/// </summary>
	/// <param name="in_Process">A reference to the vector representing the heap of 2D fast marching points.</param>
	/// <param name="pos">The index of the element to heapify down.</param>
	void heapifyDown2D(std::vector<pointFastMarching2D>& in_Process, int pos);

	/// <summary>
	/// Restores the heap property by moving the element at the specified position up in a 2D fast marching heap.
	/// </summary>
	/// <param name="in_Process">A reference to the vector representing the heap of pointFastMarching2D elements.</param>
	/// <param name="pos">The index of the element to move up in the heap.</param>
	void heapifyUp2D(std::vector<pointFastMarching2D>& in_Process, int pos);

	/// <summary>
	/// Transforms a 2D vector of pointFastMarching2D objects into a heap structure.
	/// </summary>
	/// <param name="in_Process">A reference to a vector of pointFastMarching2D objects to be heapified.</param>
	void heapifyVector2D(std::vector<pointFastMarching2D>& in_Process);

	/// <summary>
	/// Deletes or processes the root element from a 2D fast marching heap represented by a vector of pointFastMarching2D objects.
	/// </summary>
	/// <param name="in_Process">A reference to a vector containing pointFastMarching2D objects that represent the heap.</param>
	void deleteRootHeap2D(std::vector<pointFastMarching2D>& in_Process);

	/// <summary>
	/// Adds a 2D point to the processing heap for fast marching algorithms.
	/// </summary>
	/// <param name="in_Process">A reference to the vector representing the current processing heap of 2D points.</param>
	/// <param name="point">The 2D point to be added to the processing heap.</param>
	void addPointHeap2D(std::vector<pointFastMarching2D>& in_Process, pointFastMarching2D point);
	
	int getIndexFromHeap2D(std::vector<pointFastMarching2D>& in_Process, size_t i, size_t j);