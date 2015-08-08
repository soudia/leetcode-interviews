package com.interviews.utils;

public class MatrixMergeSort {

	private int array[][];

	public MatrixMergeSort(int array[][]) {
		this.array = array;
	}

	public void sort() {
		if (array == null)
			return;
		mergeSort(0, array.length);

	}

	private void mergeSort(int lb, int ub) {
		if (ub <= lb)
			return;

		int pivot = (int) Math.floor((ub - lb) / 2);

		mergeSort(lb, pivot);
		mergeSort(pivot, ub);
		merge(lb, pivot, ub);
	}

	private void merge(int lb, int pivot, int ub) {
		int a1[][] = new int[pivot - lb][];
		int a2[][] = new int[ub - pivot][];

		for (int i = lb; i < pivot; i++) {
			a1[i] = array[i];
		}

		for (int i = pivot; i < ub; i++) {
			a2[i] = array[i];
		}

		int k = 0;
		int j = 0;

		for (int i = lb; i < ub; i++) {
			if (a1[k][a1[0].length - 1] <= a2[j][a2[0].length - 1])
				k++;
			else {
				array[i] = a2[j];
				j++;
			}
		}
	}
}
