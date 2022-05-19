#ifndef __HEAPS__
#define __HEAPS__

#include <cmath>
#include <ostream>
#include <vector>
#include <utility>
#include <climits>

void swap(std::pair<int,double> *x, std::pair<int,double> *y);

class MHeap {
	std::vector<std::pair<int, double>> arr;
	int capacity;
	int heap_size;

public:
	MHeap(int capacity);

	void MinHeapify(int);

	int parent(int i) { return (i-1)/2; }
	int left(int i) { return (2*i + 1); }
	int right(int i) { return (2*i + 2); }

	int extractMin();
	void decreaseKey(int i, double new_weight);
	int getMin() { return arr[0].first; }
	void insertKey(int k, double weight);
	bool isEmpty() { return heap_size == 0; }
};
#endif
