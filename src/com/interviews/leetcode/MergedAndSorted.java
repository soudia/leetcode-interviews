package com.interviews.leetcode;

public class MergedAndSorted {
	int[] a;
	int start1;
	int[] b;
	int start2;

	int[] insertionPoints;

	public MergedAndSorted(int[] a, int start1, int n1, int[] b, int start2,
			int n2) {
		this.a = a;
		this.start1 = start1;
		this.b = b;
		this.start2 = start2;

		insertionPoints = new int[n2];

		int searchStart = start1;
		for (int i = 0; i < n2; ++i) {
			int searchEnd = start1 + n1;

			// binary search for b[start2+i] inside array a,
			// starting at the previous insertion point
			while (searchStart < searchEnd) {
				int searchMid = (searchStart + searchEnd) / 2;
				if (b[start2 + i] < a[searchMid]) {
					searchEnd = searchMid;
				} else if (b[start2 + i] > a[searchMid]) {
					searchStart = searchMid + 1;
				} else {
					searchStart = searchMid;
					break;
				}
			}
			insertionPoints[i] = searchStart;
		}
	}

	/**
	 * Find the element at given index of the `merged' array
	 */
	public int get(int index) {
		index += start1;
		int i = 0;
		while (i < insertionPoints.length && index > insertionPoints[i] + i) {
			++i;
		}
		if (i < insertionPoints.length && index == insertionPoints[i] + i) {
			return b[start2 + i];
		} else {
			return a[index - i];
		}
	}
}
