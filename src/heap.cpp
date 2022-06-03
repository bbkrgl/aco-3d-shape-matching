#include "heap.h"

MHeap::MHeap(int cap)
{
	heap_size = 0;
	capacity = cap;
	arr = std::vector<std::pair<int, double>>(cap);
}

void MHeap::insertKey(int k, double weight)
{
	if (heap_size == capacity) {
		return;
	}

	heap_size++;
	int i = heap_size - 1;
	arr[i] = std::pair<int, double>(k,weight);

	while (i != 0 && arr[parent(i)].second > arr[i].second) {
			swap(&arr[i], &arr[parent(i)]);
			i = parent(i);
		}
}

void MHeap::decreaseKey(int i, double new_weight)
{
	arr[i] = std::pair<int, double>(i, new_weight);
	while (i != 0 && arr[parent(i)].second > arr[i].second) {
			swap(&arr[i], &arr[parent(i)]);
			i = parent(i);
		}
}
  
int MHeap::extractMin()
{
	if (heap_size <= 0)
		return INT_MAX;

	if (heap_size == 1) {
		heap_size--;
		return arr[0].first;
	}

	std::pair<int, double> root = arr[0];
	arr[0] = arr[heap_size-1];
	heap_size--;
	MinHeapify(0);

	return root.first;
}
  
void MHeap::MinHeapify(int i)
{
	int l = left(i);
	int r = right(i);
	int smallest = i;

	if (l < heap_size && arr[l].second < arr[i].second)
		smallest = l;

	if (r < heap_size && arr[r].second < arr[smallest].second)
		smallest = r;

	if (smallest != i) {
		swap(&arr[i], &arr[smallest]);
		MinHeapify(smallest);
	}
}
  
void swap(std::pair<int,double> *x, std::pair<int, double> *y)
{
	std::pair<int,double> temp = *x;
	*x = *y;
	*y = temp;
}
