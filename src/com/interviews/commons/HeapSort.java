package com.interviews.commons;

/*
 * Know the difference between HeapSort and QuickSort.
 */

public class HeapSort {

	private Node[] heapArray;
	private int maxSize;
	private int currentSize = -1;

	public class Node {
		public int key;
		public int index; // needed only to solve a specific problem.

		public Node(int key, int index) {
			this.index = index;
			this.key = key;
		}

		public void setIndex(int index) {
			this.index = index;
		}
	}

	public HeapSort(int maxSize) {
		heapArray = new Node[maxSize];
		this.maxSize = maxSize;
	}

	public HeapSort(Node[] heapArray) {
		buildMinHeap(heapArray);
	}

	public void insert(int key, int index) {
		if (maxSize == currentSize)
			return;
		currentSize++;
		Node node = new Node(Integer.MAX_VALUE, index);
		heapArray[currentSize] = node;
		decreaseKey(currentSize, key);
	}

	private void decreaseKey(int i, int val) {
		if (heapArray[i].key <= val)
			return;
		heapArray[i].key = val;
		while (i >= 0 && heapArray[i].key < heapArray[parent(i)].key) {
			int key = heapArray[parent(i)].key;
			int index = heapArray[parent(i)].index;
			heapArray[parent(i)].key = heapArray[i].key;
			heapArray[parent(i)].index = heapArray[i].index;
			heapArray[i].key = key;
			heapArray[i].index = index;
			i = parent(i);
		}
	}

	public void delete(int i) {
		int key = heapArray[i].key;
		int index = heapArray[parent(i)].index;
		heapArray[i].key = heapArray[currentSize].key;
		heapArray[currentSize].key = key;
		heapArray[i].index = index;
		minHeapify(heapArray, i);
		currentSize--;
	}

	public int left(int i) {
		return 2 * i;
	}

	public int right(int i) {
		return 2 * i + 1;
	}

	public int parent(int i) {
		return (int) Math.floor(i / 2);
	}

	public boolean isEmpty() {
		return currentSize < 0;
	}

	public void sort() {
		sort(heapArray);
	}

	private void sort(Node heapArray[]) {
		buildMinHeap(heapArray);
		for (int i = currentSize; i >= 1; i--) {
			int key = heapArray[0].key;
			int index = heapArray[0].index;
			heapArray[i].key = key;
			heapArray[i].index = index;
			heapArray[0].key = heapArray[i].key;
			heapArray[0].index = heapArray[i].index;
			minHeapify(heapArray, 0);
		}
	}

	public int extractMax(Node heapArray[]) {
		int max = heapArray[0].key;
		heapArray[0].key = heapArray[currentSize].key;
		heapArray[0].index = heapArray[currentSize].index;
		minHeapify(heapArray, 0);
		return max;
	}

	public Node extractMin() {
		int max = heapArray[0].key;
		int index = heapArray[0].index;
		heapArray[0].key = heapArray[currentSize].key;
		heapArray[0].index = heapArray[currentSize].index;
		minHeapify(heapArray, 0);
		return new Node(max, index);
	}

	private void minHeapify(Node heapArray[], int i) {
		int l = left(i);
		int r = right(i);
		int smallest = i;
		if (l < currentSize && heapArray[l].key <= heapArray[i].key) {
			smallest = l;
		}
		if (r < currentSize && heapArray[r].key <= heapArray[smallest].key) {
			smallest = r;
		}
		if (smallest != i) {
			int key = heapArray[i].key;
			heapArray[i].key = heapArray[smallest].key;
			heapArray[smallest].key = key;
			minHeapify(heapArray, smallest);
		}
	}

	private void buildMinHeap(Node heapArray[]) {
		for (int i = (int) Math.floor(heapArray.length / 2); i >= 0; i--) {
			minHeapify(heapArray, i);
		}
	}

}
