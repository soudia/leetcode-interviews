package com.interviews.google;

public class MergeSort {

	private double array[];
	private double[] tempMergArr;

	public static enum Order {
		ASC, DESC;
	}

	public MergeSort(double array[]) {
		this.array = array;
	}

	public void sort() {
		sort(Order.ASC);
	}

	public void sort(Order o) {
		if (array == null)
			return;
		this.tempMergArr = new double[array.length];
		mergeSort(0, array.length - 1, o);
	}

	private void mergeSort(int lb, int ub, Order order) {
		if (ub <= lb)
			return;

		int pivot = lb + (ub - lb) / 2;

		mergeSort(lb, pivot, order);
		mergeSort(pivot + 1, ub, order);
		merge(lb, pivot, ub, order);
	}

	private void merge(int lb, int middle, int ub, Order o) {
		for (int i = lb; i <= ub; i++) {
			tempMergArr[i] = array[i];
		}
		int i = lb, k = lb;
		int j = middle + 1;
		while (i <= middle && j <= ub) {
			if (o == Order.DESC) {
				if (tempMergArr[i] >= tempMergArr[j]) {
					array[k] = tempMergArr[i];
					i++;
				} else {
					array[k] = tempMergArr[j];
					j++;
				}
				k++;
			} else {
				if (tempMergArr[i] <= tempMergArr[j]) {
					array[k] = tempMergArr[i];
					i++;
				} else {
					array[k] = tempMergArr[j];
					j++;
				}
				k++;
			}
		}
		while (i <= middle) {
			array[k] = tempMergArr[i];
			k++;
			i++;
		}

	}
}
