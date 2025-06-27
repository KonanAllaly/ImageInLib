#include "front_propagation.h"

void swap2dPoints(pointFastMarching2D* a, pointFastMarching2D* b) {
	pointFastMarching2D temp = *a;
	*a = *b;
	*b = temp;
}

//pos is the index of the element to be heapified
void heapifyDown2D(std::vector<pointFastMarching2D>& in_Process, int pos) {
	//we use type int for indexes because we do operations like pos--
	int length_array = in_Process.size();
	int current = pos;
	int left_child = 2 * pos + 1;
	int right_child = 2 * pos + 2;

	dataType val_current = 0.0, val_left = 0.0, val_right = 0.0;

	if (current >= 0 && current < length_array) {
		val_current = in_Process[current].arrival;
	}

	if (left_child < length_array) {
		val_left = in_Process[left_child].arrival;
		if (val_left < val_current) {
			current = left_child;
			val_current = in_Process[current].arrival;
		}
	}

	if (right_child < length_array) {
		val_right = in_Process[right_child].arrival;
		if (val_right < val_current) {
			current = right_child;
		}
	}

	if (current != pos) {
		swap2dPoints(&in_Process[pos], &in_Process[current]);
		heapifyDown2D(in_Process, current);
	}
}

void heapifyUp2D(std::vector<pointFastMarching2D>& in_Process, int i) {

	int current = i;

	if (i > 0) {
		int parent = (i - 1) / 2;
		dataType val_current = in_Process[current].arrival;
		dataType val_parent = in_Process[parent].arrival;
		if (val_current < val_parent) {
			current = parent;
		}
	}

	if (current != i) {
		swap2dPoints(&in_Process[current], &in_Process[i]);
		heapifyUp2D(in_Process, current);
	}

}

void heapifyVector2D(std::vector<pointFastMarching2D>& in_Process) {
	//we use type int for indexes because we do operations like pos--
	int length_array = in_Process.size();
	int indx, start = length_array / 2 - 1;
	for (indx = start; indx >= 0; indx--) {
		heapifyDown2D(in_Process, indx);
	}
}

void deleteRootHeap2D(std::vector<pointFastMarching2D>& in_Process) {
	//we use type int for indexes because we do operations like pos--
	int l = in_Process.size();
	swap2dPoints(&in_Process[0], &in_Process[l - 1]);
	in_Process.pop_back();
	heapifyDown2D(in_Process, 0);
}

void addPointHeap2D(std::vector<pointFastMarching2D>& in_Process, pointFastMarching2D point) {
	//we use type int for indexes because we do operations like pos--
	in_Process.push_back(point);
	int l = in_Process.size();
	heapifyUp2D(in_Process, l - 1);
}

int getIndexFromHeap2D(std::vector<pointFastMarching2D>& in_Process, size_t i, size_t j) {
	for (int ind = 0; ind < in_Process.size(); ind++) {
		if (in_Process[ind].x == i && in_Process[ind].y == j) {
			return ind;
		}
	}
	return -1; //not found
}