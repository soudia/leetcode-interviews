package com.interviews.facebook;

public class MergeSort {

	private int array[];

	public static enum Order {
		ASC, DESC;
	}

	public MergeSort(int array[]) {
		this.array = array;
	}

	public void sort() {
		sort(Order.ASC);
	}

	public void sort(Order o) {
		if (array == null)
			return;
		mergeSort(0, array.length, Order.ASC);
	}

	private void mergeSort(int lb, int ub, Order order) {
		if (ub <= lb)
			return;

		int pivot = (int) Math.floor((ub - lb) / 2);

		mergeSort(lb, pivot, order);
		mergeSort(pivot, ub, order);
		merge(lb, pivot, ub, order);
	}

	private void merge(int lb, int pivot, int ub, Order o) {
		int a1[] = new int[pivot - lb];
		int a2[] = new int[ub - pivot];

		for (int i = lb; i < pivot; i++) {
			a1[i] = array[i];
		}

		for (int i = pivot; i < ub; i++) {
			a2[i] = array[i];
		}

		int k = 0;
		int j = 0;

		for (int i = lb; i < ub; i++) {
			if (o == Order.DESC) {
				if (a1[k] >= a2[j])
					k++;
				else {
					array[i] = a2[j];
					j++;
				}
			} else {
				if (a1[k] <= a2[j])
					k++;
				else {
					array[i] = a2[j];
					j++;
				}
			}
		}
	}
}
