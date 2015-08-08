package com.interviews.google;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.Stack;

import com.interviews.commons.Graph;
import com.interviews.commons.Node;
import com.interviews.commons.TreeTraversal;
import com.interviews.google.MergeSort.Order;

/**
 * Useful links
 * 
 * http://crackprogramming.blogspot.ca/?view=magazine
 * http://www.google.com/diversity/hiring.html
 * http://stackoverflow.com/questions
 * /21078445/find-connected-components-in-a-graph
 * 
 * http://stackoverflow.com/questions/20251962/fastest-way-to-find-n-closest-
 * nodes-to-a-given-node-in-a-directed-weighted-graph
 * 
 * @author Ousmane
 *
 */
public class GoogleInterviews {

	/**
	 * Write a function that would print all positive numbers smaller than n
	 * that can be expressed as the sum of two cubes in two different ways.
	 */
	public void posNumbersCubes(int n) {
		int size = (int) Math.pow((double) n, (double) 1 / 3);

		for (int i = 1; i <= size; i++) {
			for (int j = i; j <= size; j++) {
				int first = (int) Math.pow(i, 3);
				int second = (int) Math.pow(j, 3);
				if (first + second == n) {
					System.out.println("First " + first + " Second" + second);
				}
			}
		}
	}

	/**
	 * How would you sort 1 billion integers stored in an array
	 */
	public void sortBigArray(int[] array) {

	}

	/**
	 * You are given an unsorted sequence of integers 'a'. Find the longest
	 * subsequence 'b' such that elements of this subsequence are strictly
	 * increasing numbers. Elements in the subsequence 'b' must appear in the
	 * same relative order as in the subsequence 'a'. You may assume that 'a'
	 * can fit to the memory.
	 */
	public String getLongestSubs(int num[]) {
		int max = Integer.MIN_VALUE;
		int[] dp = new int[num.length];
		int track[] = new int[num.length];

		int start = 0;

		Arrays.fill(dp, 0);
		Arrays.fill(track, -1);

		for (int i = 1; i < num.length; i++) {
			for (int j = 0; j < i; j++) {
				if (num[i] > num[j]) {
					if (dp[j] + 1 > dp[i]) {
						dp[i] = dp[j] + 1;
						track[i] = j;
					}
				}
				if (dp[i] > max) {
					max = dp[i];
					track[i] = j;
				}
			}
		}
		StringBuilder st = new StringBuilder();
		int index = start;
		st.append(num[start]);
		while (track[index] != -1) {
			index = track[index];
			st.append(num[index]);
		}
		return st.reverse().toString();
	}

	/**
	 * Find popular item in a sorted array of natural numbers. An item is
	 * popular if its repeated n/4 times or more. Example: Input - 123444445678
	 * Popular item is 4.
	 * 
	 * Solution: If an item is repeated at least n/4 times, then
	 */
	public int mostPopular(int nums[]) {

		if (nums == null) {
			return Integer.MIN_VALUE;
		}

		for (int i = 0; i < nums.length / 4; i++) {
			int pos = (int) ((i + 1) * nums.length) / 4;
			pos = firstOccurrence(nums[pos], i, pos, nums);
			if (pos == i) {
				return pos;
			}
		}
		return -1;
	}

	private int firstOccurrence(int val, int lb, int ub, int nums[]) {
		if (lb <= ub) {
			int middle = lb + (ub - lb) / 2;
			if (nums[middle] == val) {
				return middle;
			}
			if (nums[middle] < val)
				return firstOccurrence(val, lb, middle - 1, nums);
			return firstOccurrence(val, middle + 1, ub, nums);
		}
		return -1;
	}

	/**
	 * Given an array A[1..N] of positive integers, you have to sort it in
	 * ascending order in the following manner: In every operation, select any
	 * two sub-arrays A[i...(i+k-1)] and A[j...(j+k-1)] such that (i+k-1) < j
	 * and swap A[i] with A[j], A[i+1] with A[j+1], ..., and A[i+k-1] with
	 * A[j+k-1]. How can we figure out minimum number of swaps in most effective
	 * way ?
	 */
	public int minSwaps(int nums[]) {
		if (nums == null) {
			return Integer.MIN_VALUE;
		}
		int k = 1, swap = 0;
		int length = nums.length - 1;
		for (int i = 0; i <= length; i++) {
			if (nums[i + 1] > nums[i]) {
				k++;
			} else {
				int j = (i + k) > length ? length : i + k;
				while ((j > k) && (nums[j] > nums[--j]))
					;
				if (j > i) {
					int temp = nums[i];
					nums[i] = nums[j];
					nums[j] = temp;
					swap++;
				}
			}
		}
		return swap;
	}

	/**
	 * On a table, there are four papers, one with the letter X written on top,
	 * one with the letter Y written on top, one with the number 1 written on
	 * top, and one with the number 2 written on top. Every paper has one of
	 * those letters on one side and one of those numbers on the other side. To
	 * prove that a paper with the letter X contains even number on the other
	 * side, what is the minimum number of conditions we have to check ?
	 * Solution: Check that the paper marked Y has 1 on its back.
	 */

	/**
	 * Design and implement a calculator that can calculate expressions like: +
	 * 2 4, * 8 (+ 7 12), (+ 7 (* 8 12) (* 2 (+ 9 4) 7) 3 ). All items are space
	 * delimited.
	 * 
	 * Solution: This expression is printed as if we were traversing a tree in
	 * pre-order traversal where the internal nodes of the tree are operators
	 * (*, +) and the leaves are numbers. We can reconstruct the tree and apply
	 * an in-order traversal on the tree to evaluate the expression.
	 * 
	 * To reconstruct the tree, operators are internal nodes, elements between
	 * parenthesis are right nodes and numbers before parenthesis are left nodes
	 */
	public int compute(String expr) {
		Node<Character> tree = parseTree(0, expr);
		if (tree == null)
			return Integer.MIN_VALUE;
		return compute(tree);
	}

	public int compute(Node<Character> tree) {
		int result = 0;
		char lVal = ' ', rVal = ' ', val = ' ';
		val = tree.val;
		Node<Character> left = tree.left;
		if (left != null) {
			lVal = left.val;
			if (val == '+') {
				rVal = '0';
			} else {
				rVal = '1';
			}
			compute(tree.left);
		}
		Node<Character> right = tree.right;
		if (right != null) {
			rVal = right.val;
			if (val == '+') {
				lVal = '0';
			} else {
				lVal = '1';
			}
			compute(tree.right);
		}
		result = calculate(result, val, lVal, rVal);
		if (right == null && left == null) {
			return result;
		}

		return result;
	}

	private int calculate(int result, char val, char lVal, char rVal) {
		switch (val) {
		case '+':
			result += (int) lVal + (int) rVal;
		case '*':
			if (result == 0)
				result = 1;
			result *= (int) lVal + (int) rVal;
		}
		return result;
	}

	private Node<Character> parseTree(int pos, String expr) {
		Node<Character> node = null;
		if (pos >= expr.length()) {
			return node;
		}

		char aChar = expr.charAt(pos);

		if (pos == expr.length() - 1 && aChar == ')') {
			return node;
		}

		if (aChar == '(' && pos == 0) {
			aChar = '+';
		}

		while (aChar == ')' && pos < expr.length()) {
			aChar = expr.charAt(++pos);
		}

		node = new Node<Character>(aChar);
		if (aChar != '(') {
			node.left = parseTree(pos + 2, expr);
		} else {
			node.right = parseTree(pos + 2, expr);
		}

		return node;
	}

	/**
	 * Write a function that receives a stream of pairs: id (some unique
	 * integer) + weight (positive real number). The stream ends with the
	 * special pair {0, 0}. The function should return one of the unique ids at
	 * random, with probability proportional to its weight (compared to the
	 * total weights sum when stream has ended). Stream size and weights sum are
	 * unknown in advance, and you are allowed O(1) memory consumption.
	 */
	public int getRandom(double[][] pair) {
		int result = 0;
		double totalWeight = pair[0][1];
		/*
		 * we assume the length to be big enough to simulate the stream.
		 */
		for (int i = 1; i < pair.length; i++) {
			if (pair[i][0] == 0 && pair[i][1] == 0) {
				return result;
			}
			totalWeight += pair[i][1];
			double chance = Math.random() * totalWeight;
			if (chance < totalWeight) {
				result = i;
			}
		}
		return result;
	}

	/**
	 * qaz is a value for a number that is less than the other next values which
	 * have indexes greater than the index of this number. Example: 33, 25, 26,
	 * 58, 41, 59. qaz of 33 is 3 where 33 is less than 3 numbers (58, 41, 59).
	 * qaz of (25) is 4. The question is to find the max qaz. You can solve it
	 * in O(n^2) but the requirement is O(nlog(n)).
	 * 
	 * Solution 1: Use queues,
	 * 
	 * Solution 2: Try merge sort
	 */

	public int maxQaz(int nums[]) {
		Queue<Integer> queue = new Queue<Integer>();
		for (int i = 0; i < nums.length; i++) {
			queue.enqueue(nums[i]);
		}

		Queue<Integer> larger = new Queue<Integer>();
		Queue<Integer> smaller = new Queue<Integer>();

		int i = 0, curr = queue.dequeue().val;
		int qaz[] = new int[nums.length];

		while (!queue.isEmpty()) {
			int val = queue.dequeue().val;
			if (val <= curr) {
				smaller.enqueue(val);
			} else {
				larger.enqueue(val);
			}
			if (queue.isEmpty()) {
				if (!larger.isEmpty()) {
					qaz[i++] = larger.size();
				} else {
					qaz[i++] = 0;
				}
				if (!smaller.isEmpty()) {
					while (!smaller.isEmpty()) {
						queue.enqueue(smaller.dequeue().val);
					}
				} else {
					while (!larger.isEmpty()) {
						queue.enqueue(larger.dequeue().val);
					}
				}
			}
			curr = val;
		}

		int max = qaz[0], j = 0;

		for (i = 1; i < qaz.length; i++) {
			if (qaz[i] > max) {
				max = qaz[i];
				j = i;
			}
		}
		return nums[j];
	}

	/*
	 * Right solution: from the Internet void UpdateQAZ(int low, int mid, int
	 * high) { int i,j;
	 * 
	 * added_qaz = 0;
	 * 
	 * i = mid; j = high;
	 * 
	 * while ( ( i >= 0 ) && ( j > mid ) ) { if ( A[ j ] > A[ i ] ) {
	 * added_qaz++; j--; } else { A[ i ].qaz += added_qaz; i--; } }
	 * 
	 * if ( j <= mid ) { for ( x = i ; x >=low ; x-- ) { A[ x ].qaz +=
	 * added_qaz; } } }
	 */

	private class Queue<T> {
		QNode<T> head = null, tail = null;

		public void enqueue(T val) {
			QNode<T> node = new QNode<T>(val);
			if (head == null) {
				head = node;
				tail = head;
			} else {
				tail.next = node;
				tail = node;
			}
		}

		public boolean isEmpty() {
			return head == null;
		}

		public QNode<T> dequeue() {
			QNode<T> node = head;
			head = node.next;
			return node;
		}

		public int size() {
			QNode<T> curr = head;
			int size = 0;
			while (curr != tail) {
				size++;
			}
			return size;
		}
	}

	private class QNode<T> {
		QNode<T> next;
		T val;

		public QNode(T val) {
			this.val = val;
		}
	}

	/**
	 * You are given a binary search tree and a positive integer K. Return the
	 * k-th element of the tree. No pre-processing or modifying of the tree is
	 * allowed.
	 */

	public Node<Character> kThElement(int k, Node<Character> root) {
		Node<Character> mininum = minimum(root);
		Node<Character> succ = mininum;
		for (int i = k - 1; i > 0; i--) {
			succ = successor(mininum);
		}
		return succ;
	}

	private Node<Character> successor(Node<Character> node) {
		if (node.right != null) {
			return minimum(node.right);
		}
		Node<Character> parent = node.parent;
		while (parent != null && parent.right == node) {
			node = parent;
			parent = parent.parent;
		}
		return node;
	}

	private Node<Character> minimum(Node<Character> right) {
		Node<Character> node = right;
		while (node.left != null) {
			node = node.left;
		}
		return node;
	}

	/**
	 * You are given an array of distinct numbers. You need to return an index
	 * to a "local minimum" element, which is defined as an element that is
	 * smaller than both its adjacent elements. In the case of the array edges,
	 * the condition is reduced to one adjacent element. Reach a solution with
	 * better complexity than the trivial solution of O(n).
	 */
	public int localMinTrivial(int nums[]) { // trivial solution
		for (int i = 0; i < nums.length; i++) {
			if (i == 0) {
				if (nums[i] <= nums[i + 1]) {
					return i;
				} else if (i == nums.length - 1) {
					if (nums[i] <= nums[i - 1]) {
						return i;
					}
				}
			} else if (nums[i] <= Math.min(nums[i - 1], nums[i + 1])) {
				return i;
			}
		}
		return Integer.MIN_VALUE;
	}

	public int localMin(int nums[]) {
		return localMin(0, nums.length, nums);
	}

	private int localMin(int lb, int ub, int nums[]) {

		if (lb > ub) {
			return Integer.MIN_VALUE;
		}

		int middle = (lb + ub) / 2;
		int left = nums[lb];
		int right = nums[ub];
		if (middle == lb + 1) {
			if (left <= nums[middle]) {
				return lb;
			}
		} else if (middle == ub - 1) {
			if (right <= nums[middle]) {
				return ub;
			}
		} else if (nums[middle] <= Math.min(nums[middle - 1], nums[middle + 1])) {
			return middle;
		} else {
			left = localMin(lb, middle - 1, nums);
			right = localMin(middle + 1, ub, nums);
		}
		return Integer.MIN_VALUE;
	}

	/**
	 * You are given a positive integer number N. Return the minimum number K,
	 * such that N can be represented as K integer squares. Example: 9 ---> 1 (9
	 * = 3 ^2), 8 ---> 2 (8 = 2^2 + 2^2), 15 ---> 4 (15 = 3^2 + 2^2 + 1^2 + 1^2)
	 */

	public int minSumSquareSeq(int num) {
		if (num == 0 || num == 1) {
			return num;
		}
		int sqrt = (int) Math.sqrt(num);
		int diff = num - (sqrt * sqrt);
		return 1 + minSumSquareSeq(diff);
	}

	/**
	 * You are given a double linked list and an array of pointers to elements
	 * in this list (no assumptions can be made on the array - number of
	 * pointers, order and duplicates allowed). Return the number of sequences
	 * of elements (groups of consecutive elements), pointed by the array.
	 * 
	 * For example, if this is the array (number relates to index in the list,
	 * not the actual pointer value): 9, 1, 3, 7, 8, 5, 2. Then, the output is
	 * 3, representing these sequences: [1, 2, 3], [5], [7, 8, 9].
	 */

	/**
	 * A book contains pages numbered from 1 - N. Imagine now that you
	 * concatenate all page numbers in the book such that you obtain a sequence
	 * of numbers which can be represented as a string. You can compute the
	 * number of occurrences 'k' of certain digit 'd' in this string. For
	 * example let N = 12, d = 1, hence, s = '123456789101112' => k = 5. Since,
	 * digit '1' occurs five times in that string.
	 * 
	 * Problem: write a method that, given a digit 'd' and number of its
	 * occurrences 'k', returns a number of pages N. More precisely, return a
	 * lower and an upper bound of this number.
	 * 
	 * 1234567891011121314151617181920212223242526, d = 2, k = 10
	 * 1234567891011121314151617181920212223242526272829303132333435;
	 * 123456789101112131415161718192021222324252627282930313233343536373839404142434445
	 * , d = 4, k = 10
	 */

	public int[] pageNumber(String seq, int digit, int occ) {
		int pages[] = new int[2]; // Try it again later.
		return pages;
	}

	/**
	 * You are given a string of number characters (S) like the following one
	 * for example: "132493820173849029382910382". Now, let us assume we tie the
	 * letters to numbers in order such that: A = "0", B = "1", C = "2", ..., M
	 * = "12", N = "13", ..., Y = "24", Z = "25". Write an algorithm to
	 * determine how many strings of letters we can make with S, by converting
	 * from numbers to letters.
	 */

	public int numWords(String s) {
		Map<String, Character> map = new HashMap<String, Character>();
		int value = 0;
		for (char c = 'A'; c <= 'Z'; c++) {
			map.put(String.valueOf(value), c);
			value++;
		}

		int dp[][] = new int[s.length()][s.length()];

		for (int i = 0; i < s.length(); i++) {
			dp[i][i] = 1;
		}

		for (int i = 0; i < s.length(); i++) {
			for (int j = i + 1; j < s.length(); j++) {
				if ((j - i) > 1) {
					dp[i][j] = 0;
				} else {
					if (map.containsKey(s.substring(i, j))) {
						dp[i][j] = 1;
					} else {
						dp[i][j] = 0;
					}
				}
			}
		}

		int count = 0;
		for (int i = 0; i < s.length(); i++) {
			if (dp[0][i] == 1) {
				count++;
			}
		}
		return count;
	}

	/**
	 * You are given an array of integers 'a' that can fit in a memory. Write a
	 * method that returns an array of same length such that each element 'i' of
	 * this array is a sum of 'a' except the element a[i]. You are not allowed
	 * to use '-' operator.
	 * 
	 * Complexity required O(n)
	 */
	public int[] sum(int[] a) {
		int[] sum = new int[a.length];
		int[] left = new int[a.length];
		int[] right = new int[a.length];
		int s = 0;
		for (int i = 0; i < a.length; i++) {
			left[i] = s;
			s += a[i];
		}
		s = 0;
		for (int i = a.length - 1; i >= 0; i--) {
			right[i] = s;
			s += a[i];
		}
		for (int i = 0; i < a.length; i++) {
			sum[i] = left[i] + right[i];
		}
		return sum;
	}

	/**
	 * Given two binary trees, A and B, w can define A as a subset of B if the
	 * root of A exists in B and if we can superimpose A on B at that location
	 * without changing B. That is the subtrees from the root of A and the node
	 * in B is the same as the root of A.
	 * 
	 * Write an algorithm to determine if one tree is a subset of another tree.
	 */

	public boolean isSubTree(Node<Character> r1, Node<Character> r2) {
		String str1 = treenode(r1);
		String str2 = treenode(r2);

		if (str1.length() < str2.length()) {
			return str2.indexOf(str1) >= 0;
		} else {
			return str1.indexOf(str2) >= 0;
		}
	}

	private String treenode(Node<Character> r) {
		String nodes = "";
		if (r == null) {
			return nodes;
		}
		nodes += r.val;
		treenode(r.left);
		treenode(r.right);
		return nodes;
	}

	/**
	 * Given a binary search tree (BST), write an algorithm that will convert
	 * this BST into a doubly linked list that is sorted (ascending or
	 * descending order) and return the first element in this list.
	 */

	public Node<Character> toBST(Node<Character> node) {
		Node<Character> head = null, curr = null;
		Stack<Node<Character>> stack = new Stack<Node<Character>>();
		while (!stack.isEmpty() || node != null) {
			if (node != null) {
				stack.push(node.left);
			} else {
				node = stack.pop();
				if (head == null) {
					head = node;
					curr = node;
				} else {
					curr.right = node;
					node.left = curr;
					curr = node;
				}
				node = node.right;
			}
		}
		return head;
	}

	/**
	 * Given a vector of integers, v[i] represents the stock price on day 1.
	 * Now, you may do at most K transactions. You must sell your stock before
	 * you buy it again and that means you cannot have two stocks at the same
	 * time. Write a program to find max profit you can get.
	 * 
	 * Type of Knapsack problem.
	 */
	public int maxProfit(int v[], int k) {

		return 0;
	}

	/**
	 * Assume we only take the least significant digit of each value in
	 * Fibonacci sequence, and form the sequence of digits into pairs. In those
	 * pairs, the first value of one pair is the same as second value of its
	 * predecessor.
	 */

	public List<int[]> outputPairs(int n) {
		List<int[]> pairs = new ArrayList<int[]>();
		fibonacci(n, pairs);
		return pairs;
	}

	private int fibonacci(int n, List<int[]> pairs) {
		int[] aPair = new int[2];
		if (n == 0) {
			return 0;
		} else if (n == 1) {
			aPair[0] = 1;
			aPair[1] = 1;
			pairs.add(aPair);
			return 1;
		} else {
			aPair[0] = fibonacci(n - 2, pairs);
			aPair[1] = fibonacci(n - 1, pairs);
			pairs.add(aPair);
			return aPair[0] + aPair[1];
		}
	}

	/**
	 * Given an unsorted array of natural numbers where numbers repeat in the
	 * array, output numbers in the order of frequency.
	 * 
	 * Example: [0, 0, 100, 3, 5, 4, 6, 4, 2, 100, 2, 100],
	 * 
	 * n == 2 => [100, 0] or [100, 2]
	 */
	public int[] freqOuput(int n, int array[]) {
		if (array == null || n == 0) {
			return null;
		}
		Map<Integer, Number> occ = new HashMap<Integer, Number>();
		for (int i = 0; i < array.length; i++) {
			if (!occ.containsKey(array[i])) {
				Number num = new Number();
				occ.put(array[i], num);
			} else {
				Number num = occ.get(array[i]);
				num.occ++;
				occ.put(array[i], num);
			}
		}
		int[] vals = new int[n];
		int maxIndex = array[0];
		Number num = occ.get(maxIndex);
		int max = num.occ, j = 0;

		while (true) {
			for (int i : occ.keySet()) {
				num = occ.get(i);
				if (max < num.occ && !num.visited) {
					max = num.occ;
					maxIndex = i;
				}
			}
			if (j == n) {
				return vals;
			}
			vals[j++] = maxIndex;
			max = Integer.MIN_VALUE;
			num = occ.get(maxIndex);
			num.visited = true;
			occ.put(maxIndex, num);
		}
	}

	class Number {
		int occ = 1;
		boolean visited = false;
	}

	/**
	 * Read a training set of strings. For each character in any of the strings,
	 * calculate the probability of seeing that character and store it for later
	 * use. Based on the training set, generate a random string of length N.
	 */
	public String markovString(List<String> training, int n) {
		if (training == null || n <= 0) {
			return null;
		}
		int size = 0;
		Map<Character, Integer> occ = new HashMap<Character, Integer>();
		for (String instance : training) {
			if (instance == null || instance == "") {
				continue;
			}
			for (int i = 0; i < instance.length(); i++) {
				char c = instance.charAt(i);
				if (!occ.containsKey(c)) {
					occ.put(c, 0);
				} else {
					int occur = occ.get(c) + 1;
					occ.put(c, occur);
				}
			}
			size += instance.length();
		}

		String str = "";

		while (true) {
			if (str.length() == n) {
				return str;
			}
			Random random = new Random();
			double rand = random.nextDouble();
			for (char c : occ.keySet()) {
				double prob = (1.0 * occ.get(c)) / size;
				if (rand < prob) {
					str += c;
				}
			}
		}
	}

	private static String print(int array[]) {
		if (array == null) {
			return null;
		}
		String str = "";
		for (int i = 0; i < array.length; i++) {
			if (i < array.length - 1) {
				str += array[i] + ",";
				continue;
			}
			str += array[i];
		}
		return str;
	}

	/**
	 * Print all different possible subsets summing a given number.
	 */
	public void allSubsets(int num) {
		allSubsets(4, "");
	}

	private void allSubsets(int num, String out) {
		if (num == 0) {
			System.out.println(out);
		} else if (num > 0) {
			for (int i = 1; i <= num; i++) {
				allSubsets(num - i, out + " " + Integer.toString(i));
			}
		}
	}

	/**
	 * You are given two arrays - A and B. The length of the arrays is the same
	 * (N). The values in the arrays are integers from 0 to N-1 in random order
	 * and no number is repeated. You need to list a sequence of two-element
	 * swaps required to reorder A into the order the array B. Additionally, at
	 * each step of the sequence you are allowed to swap two elements only if
	 * one of them is 0.
	 */

	public List<int[]> reorder(int A[], int B[]) {
		int i = 0;
		List<int[]> swapList = new ArrayList<int[]>();
		while (true) {
			if (i < 0 || i >= A.length) {
				return swapList;
			}
			if (A[i] != B[i]) {
				if (A[i] == 0) {
					int j = find(B[i], A);
					swapList.add(swap(A, i, j));
					i = j;
				} else {
					i++;
				}
			} else {
				if (A[i] == 0) {
					int j;
					if (i == 0) {
						j = i + 1;
					} else {
						j = i - 1;
					}
					swapList.add(swap(A, i, j));
					i = j;
				} else {
					i++;
				}
			}
		}
	}

	private int[] swap(int A[], int i, int j) {
		int[] swap = new int[2];
		int temp = A[i];
		A[i] = A[j];
		A[j] = temp;
		swap[0] = i;
		swap[1] = j;
		return swap;
	}

	private int find(int n, int A[]) {
		for (int i = 0; i < A.length; i++) {
			if (A[i] == n) {
				return i;
			}
		}
		throw new IllegalArgumentException("Element is not in the array");
	}

	/**
	 * Given an array of integers, build a new array of integers such that every
	 * 2nd element of the array is greater than its left and right element.
	 */

	public int[] newArray(int A[]) {
		int output[] = new int[Math.round(A.length / 2)];
		int i = 1;
		quickSort(A, 0, A.length);
		while (i < A.length - 1) {
			output[i - 1] = A[i - 1];
			output[i] = A[i + 1];
			output[i + 1] = A[i];
			i += 2;
		}
		return output;
	}

	private void quickSort(int A[], int lb, int ub) {
		int pivot = partition(A, lb, ub);
		if (lb < pivot - 1)
			quickSort(A, lb, pivot - 1);
		if (pivot < ub)
			quickSort(A, pivot, ub);
	}

	private int partition(int A[], int left, int right) {
		int i = left, j = right;
		int tmp;
		int pivot = A[(left + right) / 2];
		while (i <= j) {
			while (A[i] < pivot)
				i++;
			while (A[j] > pivot)
				j--;
			if (i <= j) {
				tmp = A[i];
				A[i] = A[j];
				A[j] = tmp;
				i++;
				j--;
			}
		}
		return i;
	}

	/**
	 * Given a number "n", find the least number of perfect square numbers sum
	 * needed to get "n". Correct implementation
	 */
	public int perfectSquares(int n) {
		int table[] = new int[n + 1];
		for (int i = 0; i <= n; i++) {
			table[i] = i;
		}
		int max = (int) Math.floor(Math.sqrt(n));
		for (int i = 2; i <= max; i++) {
			int squared = (int) Math.pow(i, 2);
			for (int j = 0; j <= n; j++) {
				if (squared <= j) {
					table[j] = Math.min(table[j], table[j - squared] + 1);
				}
			}
		}
		return table[n];
	}

	/**
	 * Find longest substring with "m" unique characters in a given string
	 */
	public String seqUniqChar(String str, int m) {
		if (str == null || str == "" || str.length() < m) {
			return null;
		}
		int beg = 0;
		String lStr = "";
		for (int i = m; i < str.length(); i++) {
			String subs = str.substring(beg, i);
			int uniq = uniqChars(subs);
			if (uniq == m) {
				if (lStr.length() < subs.length()) {
					lStr = subs;
				}
			} else if (uniq > m) {
				beg++;
			}
		}
		return lStr;
	}

	private int uniqChars(String str) {
		Map<Character, Boolean> uniq = new HashMap<Character, Boolean>();
		for (int i = 0; i < str.length(); i++) {
			if (!uniq.containsKey(str.charAt(i))) {
				uniq.put(str.charAt(i), true);
			}
		}
		return uniq.size();
	}

	/**
	 * Consider a social network graph like Facebook. You are throwing a party
	 * and want to invite some of your friends. Design an algorithm to select
	 * top "n" friends from "m" using the social graph.
	 * 
	 * Solution: Find the max edge weight between "me" and my friends and add
	 * the other node to the list of friends to invite. Proceed that way till
	 * you get all the m friends.
	 */

	public List<Graph<Long>.Vertex> topFriends(Graph<Long> graph, int m) {

		List<Graph<Long>.Edge> edges = new ArrayList<Graph<Long>.Edge>();
		List<Graph<Long>.Vertex> friends = new ArrayList<Graph<Long>.Vertex>();
		while (true) {
			int max = 0;
			Graph<Long>.Vertex vertex = null;
			for (Graph<Long>.Edge edge : edges) {
				if (edge.distance > max) {
					max = edge.distance;
					vertex = edge.end;
				}
			}
			edges.remove(vertex);
			if (friends.size() == m) {
				return friends;
			}
			friends.add(vertex);
		}
	}

	/**
	 * You are given printouts from an algorithm which ran over an unsorted
	 * binary tree. One printout is from an in-order run and another from a
	 * pre-order run. Can you reconstruct the tree ? If so, then write an
	 * algorithm.
	 * 
	 * Examples: pre-order - F, B, A, D, C, E, G, I, H in-order - A, B, C, D, E,
	 * F, G, H, H
	 */
	public Node<Character> build(String inOrder, String preOrder) {
		if (inOrder == null || preOrder == null || inOrder == ""
				|| preOrder == "") {
			return null;
		}
		char c = preOrder.charAt(0);
		Node<Character> root = new Node<Character>(c);
		Node<Character> curr = root;
		Node<Character> parent = null;
		for (int i = 1; i < preOrder.length(); i++) {
			c = preOrder.charAt(i);
			Node<Character> node = new Node<Character>(c);
			int k = position(c, inOrder);
			int j = position(preOrder.charAt(i - 1), inOrder);
			if (j < 0 || k < 0) {
				return null;
			}
			Node<Character> pointer = curr;
			while (pointer != null) {
				j = position(pointer.val, inOrder);
				if (!(k < j)) {
					if (pointer.left != null) {
						Node<Character> rpointer = pointer;
						while (rpointer.right != null) {
							rpointer = rpointer.right;
						}
						rpointer.right = node;
						break;
					}
					pointer = pointer.parent;
				} else {
					curr.left = node;
					parent = curr;
					curr = curr.left;
					curr.parent = parent;
					break;
				}
			}
		}
		return root;
	}

	private int position(char c, String str) {
		return str.indexOf(String.valueOf(c));
	}

	/**
	 * You have N points on a 2D surface. List K points at shortest distance to
	 * the point (0, 0).
	 */
	public List<TwoD> kNN(List<TwoD> datapoints, int k) {
		if (datapoints == null) {
			return null;
		}
		int i = 0;
		double distance[] = new double[datapoints.size()];
		for (TwoD dp : datapoints) {
			distance[i++] = Math.sqrt((Math.pow(dp.x, 2) + Math.pow(dp.y, 2)));
		}
		List<TwoD> knn = new ArrayList<GoogleInterviews.TwoD>();
		for (i = 0; i < k; i++) {
			int j = minIndex(distance, i);
			knn.add(datapoints.get(j));
		}
		return knn;
	}

	private int minIndex(double[] array, int pos) {
		if (array == null || pos < 0 || pos >= array.length) {
			throw new IllegalArgumentException();
		}
		double min = array[pos];
		int minIndex = pos;
		for (int i = pos + 1; i < array.length; i++) {
			if (array[i] < min) {
				min = array[i];
				minIndex = i;
			}
		}
		array[minIndex] = array[pos];
		array[pos] = min;
		return minIndex;
	}

	class TwoD {
		int x, y;
	}

	/**
	 * There is a company with 250 employees. Its records contain: EmployeeID,
	 * ManagerID (which is a reference to the EmployeeID of the manager).
	 * 
	 * 1. List direct reports for a given ID.
	 */
	public List<Long> directReports(List<Employee> employees, long managerId) {
		if (employees == null) {
			return null;
		}
		List<Long> reports;
		Map<Long, List<Long>> org = new HashMap<Long, List<Long>>();
		for (Employee emp : employees) {
			if (!org.containsKey(emp.managerId)) {
				reports = new ArrayList<Long>();
			} else {
				reports = org.get(emp.managerId);
			}
			reports.add(emp.id);
			org.put(emp.managerId, reports);
		}
		return org.get(managerId);
	}

	/**
	 * There is a company with 250 employees. Its records contain: EmployeeID,
	 * ManagerID (which is a reference to the EmployeeID of the manager).
	 * 
	 * List all (also indirectly) reporting employees to ID.
	 */
	public List<Long> allReports(List<Employee> employees, long managerId) {
		List<Long> dReports = directReports(employees, managerId);
		List<Long> allReports = new ArrayList<Long>(dReports);
		for (long id : dReports) {
			dReports = directReports(employees, id);
			if (dReports == null)
				continue;
			allReports.addAll(dReports);
		}
		return allReports;
	}

	class Employee {
		long id;
		long managerId;
	}

	/**
	 * You have a binary tree where each node knows the number of nodes in its
	 * sub-tree (including itself). Given a node "n" and an integer "k", write a
	 * function to return the kth node in an in-order traversal. Can you do this
	 * non-recursively ?
	 */

	class BTNode {
		BTNode left;
		BTNode right;
		int num;
	}

	public BTNode kNode(BTNode root, int k) {
		int visited = 0;
		BTNode node = null;
		inOrder(root, visited, k, node);
		return node;
	}

	public void inOrder(BTNode root, int visited, int k, BTNode node) {
		if (visited == k) {
			node = root;
			return;
		}
		if (root != null) {
			inOrder(root.left, visited, k, node);
			visited++;
			inOrder(root.right, visited, k, node);
		}
	}

	public BTNode nonRecurKNode(BTNode root, int k) {
		int visited = 0;
		java.util.Queue<BTNode> queue = new LinkedList<BTNode>();
		queue.add(root);
		while (!queue.isEmpty()) {
			BTNode node = queue.remove();
			if (node != null) {
				if (node.right != null) {
					queue.add(node.right);
				}
				node = node.left;
			} else {
				node = queue.remove();
				if (visited == node.num) {
					return node;
				}
				visited++;
				node = node.right;
			}
		}
		return null;
	}

	/**
	 * Write a function to get maximum and second maximum element of a stack.
	 * The function should be in O(1) complexity. Extend it for finding kth max
	 * in O(1).
	 */
	public int maxStack(Stack<Integer> s, int k) {
		return 0;
	}

	class MaxStack {

	}

	/**
	 * Given an arraylist of N integers, (1) find a non-empty subset whose sum
	 * is a multiple of N, (2) find a non-empty subset whose sum is a multiple
	 * of 2N. Complexity O(n^3)
	 */
	public List<Integer> subsetSum(int[] numbers, int M) {
		Map<Integer, LinkedList<Integer>> subsets;
		subsets = new HashMap<Integer, LinkedList<Integer>>();
		subsets.put(0, new LinkedList<Integer>());

		for (int n : numbers) {
			Set<Integer> keys = new HashSet<Integer>(subsets.keySet());
			for (int i : keys) {
				int sum = (i + n) % M;
				if (sum == 0) {
					subsets.get(i).add(n);
					return subsets.get(i);
				}
				if (subsets.containsKey(sum))
					continue;
				@SuppressWarnings("unchecked")
				LinkedList<Integer> list = (LinkedList<Integer>) subsets.get(i)
						.clone();
				list.add(n);
				subsets.put(sum, list);
			}
		}
		return null;
	}

	/**
	 * Given a BST and a number x, find two nodes in the BST whose sum is equal
	 * to x. You cannot use extra memory in converting BST into one array and
	 * then solve this like 2sum.
	 * 
	 * Logic: 1. Difference between root node value and sum, compare with left
	 * and right. 2. If no equality found, check if you should go left or right.
	 * 3. If left, find diff in left. If null, check successor of diff in left
	 * and do the same again. If null, return null. 4. If right, compare diff
	 * with root val. If less, find diff left
	 */
	public List<Node<Integer>> bst2NodeSum(Node<Integer> root, int sum) {
		if (root == null) {
			return null;
		}
		List<Node<Integer>> list = new ArrayList<Node<Integer>>();

		int val = root.val;
		int diff = sum - val;

		Node<Integer> left = root.left;
		Node<Integer> right = root.right;
		Node<Integer> succ, node;
		while (left != null) {
			if (left.val > diff) {
				if (left.left == null) {
					diff = sum - left.val;
					node = find(root, diff);
					if (node != null) {
						list.add(left);
						list.add(node);
						return list;
					} else
						return null;
				}
				left = left.left;
			} else if (left.val == diff) {
				list.add(root);
				list.add(left);
				return list;
			} else {
				left = left.parent;
				diff = sum - left.val;
				node = find(root, diff);
				if (node != null) {
					list.add(node);
					list.add(left);
					return list;
				} else
					succ = succ(node);
				diff = sum - succ.val;
				node = find(root, diff);
				if (node != null) {
					list.add(succ);
					list.add(node);
					return list;
				} else
					return null;
			}
		}
		Node<Integer> prev = null;
		while (right != null) {
			prev = right;
			if (right.val < diff) {
				right = right.right;
			} else if (right.val == diff) {
				list.add(root);
				list.add(right);
				return list;
			} else
				right = right.left;
			if (right == null) {
				diff = sum - prev.val;
				node = find(root, diff);
				if (node != null) {
					list.add(prev);
					list.add(node);
					return list;
				} else
					return null;
			}
		}
		return null;
	}

	public void insertBST(int val) {
		Node<Integer> node = new Node<Integer>(val);
		Node<Integer> curr = rTree, prev = null;
		while (curr != null) {
			prev = curr;
			if (curr.val <= val) {
				curr = curr.right;
			} else {
				curr = curr.left;
			}
		}
		node.parent = prev;
		if (prev == null) {
			rTree = node;
		} else {
			if (prev.val <= val) {
				prev.right = node;
			} else {
				prev.left = node;
			}
		}
	}

	public Node<Integer> find(Node<Integer> root, Integer val) {
		if (root == null)
			return null;
		if (root.val < val) {
			return find(root.right, val);
		} else if (root.val > val) {
			return find(root.left, val);
		} else {
			return root;
		}
	}

	private Node<Integer> succ(Node<Integer> node) {
		if (node.right != null) {
			return min(node.right);
		}
		Node<Integer> parent = node.parent;
		while (parent != null && parent.right == node) {
			node = parent;
			parent = parent.parent;
		}
		return null;
	}

	private Node<Integer> min(Node<Integer> right) {
		Node<Integer> node = right;
		while (node.left != null) {
			node = node.left;
		}
		return node;
	}

	/**
	 * Given a string S, you are allowed to convert it to a palindrome by adding
	 * 0 or more characters in front of it. Find the length of the shortest
	 * palindrome you can create from S by applying the above transformation.
	 */
	public int shortestPalindrome(String str) {
		int i = 0, j = str.length() - 1;
		if (str.charAt(i) != str.charAt(j)) {
			str = str.charAt(j) + str;
		}
		i += 1;
		j -= 1;
		while (i <= j) {
			if (str.charAt(i) != str.charAt(j)) {
				str = str.substring(0, i) + str.charAt(j) + str.substring(i);
			} else {
				i++;
				j--;
			}
		}
		return str.length();
	}

	/***
	 * Define a Queue (FIFO) structure using only stacks (LIFO)
	 */
	public class QueueStack<T> {
		Stack<T> s1 = new Stack<T>();
		Stack<T> s2 = new Stack<T>();

		public void push(T val) {
			s1.push(val);
		}

		public T pop() {
			if (!s2.isEmpty()) {
				return s2.pop();
			}
			while (!s1.isEmpty()) {
				s2.push(s1.pop());
			}
			return s2.pop();
		}

		public T peek() {
			if (!s2.isEmpty()) {
				return s2.peek();
			}
			while (!s1.isEmpty()) {
				s2.push(s1.pop());
			}
			return s2.peek();
		}

		public boolean isEmpty() {
			return s1.isEmpty() && s2.isEmpty();
		}
	}

	/**
	 * Given a NxN matrix which contains all distinct 1 to n^2 numbers, write
	 * code to print sequence of increasing adjacent sequential numbers.
	 */
	public String seqNumbers(int[][] matrix) {
		Set<Integer> set = new HashSet<Integer>();
		for (int i = 0; i < matrix.length - 1; i++) {
			for (int j = 0; j < matrix.length - 1; j++) {
				if (matrix[i + 1][j] == matrix[i][j] + 1) {
					set.add(matrix[i][j]);
					set.add(matrix[i + 1][j]);
				}
				if (matrix[i][j + 1] == matrix[i][j] + 1) {
					set.add(matrix[i][j]);
					set.add(matrix[i][j + 1]);
				}
			}
		}
		String str = "";
		for (int num : set) {
			str += num;
		}
		return str;
	}

	/**
	 * Given a complete binary tree, find a max element.
	 */
	public int maxElement(Node<Integer> root) {
		if (root == null)
			return Integer.MIN_VALUE;
		int max = Integer.MIN_VALUE;
		return maxElement(root, max);
	}

	private int maxElement(Node<Integer> root, int max) {
		if (root != null) {
			max = Math.max(max, root.val);
			int maxC = Math.max(maxElement(root.left, max),
					maxElement(root.right, max));
			max = Math.max(max, maxC);
		}
		return max;
	}

	/**
	 * Given an array of integers, sort the array into a wave like array, namely
	 * a1 >= a2 <= a3 >= a4 <= a5 ...
	 */
	public int[] waveSort(int array[]) {
		if (array == null)
			return null;
		quickSort(array, 0, array.length - 1);
		int mid = (int) Math.ceil(array.length / 2);
		int wave[] = new int[array.length];
		int i = mid, j = mid, k = 0;
		while (j < array.length || i > 0) {
			if (j < array.length) {
				wave[k++] = array[j++];
			}
			if (i > 0) {
				wave[k++] = array[--i];
			}
		}
		return wave;
	}

	/**
	 * Determine minimum sequence of adjacent values in the input parameter
	 * array that is equal to an input parameter sum.
	 */
	public int minSeq(int array[], int num) {
		if (array == null)
			throw new IllegalArgumentException();
		int i = 0, j = 0, sum = array[0], k = 0;
		int min = Integer.MAX_VALUE;
		while (i < array.length) {
			while (sum < num && k < i)
				sum += array[++k];
			if (k >= i) {
				i++;
			} else {
				if (sum == num)
					min = Math.min(min, k - j + 1);
				k = ++j;
				sum = array[k];
			}
		}
		return min;
	}

	/**
	 * Given two strings a and b, find whether any anagram of string a is a
	 * sub-string of string b. For e.g. if a = "xyz" and b = "afdgzyxksldfm"
	 * then the program should return true.
	 */
	public boolean isSubstr(String a, String b) {

		if (b.length() < a.length())
			return false;

		Set<String> anagrams = anagrams(a);

		for (String str : anagrams) {
			int i = 0, j = 0, count = 0;
			while (true) {
				if (str.charAt(i) != b.charAt(j)) {
					j++;
				} else {
					count = ++i;
					j++;
				}
				if (j >= b.length() || i >= str.length()) {
					if (count < str.length()) {
						return false;
					}
					break;
				}
			}
		}
		return true;
	}

	// TODO To master
	public boolean hasAnagramSubstring(String src, String target) {
		if (target.length() > src.length())
			return false;
		int srcLen = src.length(), targetLen = target.length();
		int targetCount[] = new int[128], count[] = new int[128], i, j;
		// initialize
		for (i = 0; i < target.length(); ++i) {
			++targetCount[target.charAt(i)];
			++count[src.charAt(i)];
		}
		// loop
		i = 0;
		while (true) {
			// check if substring is an anagram
			for (j = 0; j < targetLen; ++j) {
				if (count[target.charAt(j)] != targetCount[target.charAt(j)])
					break;
			}
			if (j == targetLen)
				return true;
			// slide window
			if (i + 1 + targetLen > srcLen)
				break;
			--count[src.charAt(i)];
			++count[src.charAt(i + targetLen)];
			++i;
		}

		return false;
	}

	private Set<String> anagrams(String a) {
		Set<String> results = new HashSet<String>();
		anagram(results, a, a.length());
		return results;
	}

	public void anagram(Set<String> list, String a, int length) {
		if (length == 1)
			return;
		for (int i = 0; i < length; i++) {
			anagram(list, a, length - 1);
			list.add(changeOrder(a, length));
		}
	}

	private String changeOrder(String str, int newsize) {
		char[] chars = str.toCharArray();
		int j;
		int pointAt = chars.length - newsize;
		char temp = chars[pointAt];

		for (j = pointAt + 1; j < chars.length; j++) {
			chars[j - 1] = chars[j];
		}
		chars[j - 1] = temp;
		String newStr = "";
		for (int i = 0; i < chars.length; i++) {
			newStr += chars[i];
		}
		return newStr;
	}

	/**
	 * Given a board with black (1) and white (0), blacks are all connected.
	 * Find the minimum rectangle that contains all blacks.
	 */
	public List<int[]> rectangle(int[][] board) {
		if (board == null)
			return null;
		int maxI = 0, maxJ = 0;
		int minI = board.length, minJ = minI;
		List<int[]> minRec = new ArrayList<int[]>();
		for (int i = 0; i < board.length; i++) {
			int j = 0;
			while (j < board[i].length && board[i][j++] != 1)
				;
			if (board[i][j - 1] == 1) {
				maxI = i;
				minI = Math.min(i, minI);
				minJ = Math.min(j - 1, minJ);
			}
			while (j < board[i].length && board[i][j++] == 1)
				;
			if (board[i][j - 1] == 0) {
				if (j < board.length)
					maxJ = Math.max(j - 1, maxJ);
			}
		}
		minRec.add(new int[] { minI, minJ });
		minRec.add(new int[] { maxI, maxJ });
		return minRec;
	}

	/**
	 * Given an array of integers, and a range (low, high), find all continuous
	 * subsequences in the array which have sum in the range. Is there a
	 * solution better than O(n^2).
	 */
	public List<List<Integer>> sumInRange(int[] array, int low, int high) {
		if ((array == null) || (high < low))
			return null;
		List<List<Integer>> rslt = new ArrayList<List<Integer>>();
		int i = 0, j = -1, sum = 0;
		while (i < array.length) {
			while (true) {
				if ((low <= sum) && (sum <= high)) {
					if ((j - i) >= 1) {
						record(rslt, array, i, j);
					}
				}
				if (sum > high) {
					sum -= array[i++];
				} else {
					if (j < array.length - 1) {
						sum += array[++j];
						continue;
					}
					break;
				}
			}
			sum -= array[i++];
		}
		return rslt;
	}

	private void record(List<List<Integer>> rslt, int[] array, int start,
			int end) {
		if (rslt == null) {
			rslt = new ArrayList<List<Integer>>();
		}
		List<Integer> temp = new ArrayList<Integer>();
		for (int i = start; i <= end; i++) {
			temp.add(array[i]);
		}
		rslt.add(temp);
	}

	/**
	 * Given a dictionary of words, and a set of characters, judge if all the
	 * characters can form the words from the dictionary, without any character
	 * left. For example, given the dictionary {hello, world, is, my first,
	 * program}, if the character set is "iiifrssst", you should return true
	 * because you can form {is, is, first} from the set. If the character set
	 * is "eiifrsst", you should return false because you cannot use all the
	 * characters from the set. PS: The dictionary may contain tens of thousands
	 * of words, and the char set could be up to hundreds.
	 */
	public boolean formWords(List<String> dict, char chars[]) {
		if (dict == null || chars == null)
			return false;
		Set<String> set;
		Map<Character, Set<String>> map = new HashMap<Character, Set<String>>();
		for (String word : dict) {
			for (int i = 0; i < word.length(); i++) {
				char aChar = word.charAt(i);
				if (!map.containsKey(aChar)) {
					set = new HashSet<String>();
					map.put(aChar, set);
				} else {
					set = map.get(aChar);
				}
				set.add(word);
				map.put(aChar, set);
			}
		}
		for (char aChar : chars) {
			if (!map.containsKey(aChar))
				return false;
		}

		return true;
	}

	/**
	 * Given a binary search tree of integer values, return the count of nodes
	 * where all the nodes under the sub-tree lies between a given range [x, y].
	 * Assume that, there are more than 100,000 nodes.
	 * 
	 * Follow-up question: do the same for a k-d (k dimensional tree).
	 */
	public int countNodesBetweenRange(Node<Integer> root, int lb, int ub) {
		if (root != null && lb <= ub) {
			if (root.val < lb)
				return countNodesBetweenRange(root.right, lb, ub);
			else if (root.val > ub)
				return countNodesBetweenRange(root.left, lb, ub);
			else {
				return 1 + countNodesBetweenRange(root.left, lb, ub)
						+ countNodesBetweenRange(root.right, lb, ub);
			}
		}
		return 0;
	}

	/**
	 * Assume you are playing a game where you need to guess a word from a
	 * dictionary. You are given a machine you can try to guess the word. The
	 * machine will return how many characters have been matched by your guess.
	 * Design a system to crack the word.
	 * 
	 * Solution: Start with all 26 alphabets in collection. Make 26 guesses for
	 * each character one by one. If machine replies 0, then remove that letter
	 * from the collection Otherwise, machine may return 1 or a larger value
	 * (e.g. for 'miss', when you query 's', it will return 2), maintain
	 * alphabet => frequency information.
	 * 
	 * At the end, you will have total length L of the word (which is summation
	 * of the frequency) and alphabet => frequency information
	 * 
	 * Now just go through each word from the dictionary of length L and match
	 * the frequency of characters. Repeat till all the words have been explored
	 * (we need to handle anagrams- 'cat' and 'act').
	 */

	/**
	 * You are given a string 's' and a dictionary of English words. Your goal
	 * is to write an algorithm that returns all words from the dictionary that
	 * can be formed by characters from that string 's'.
	 * 
	 * Example: s = "ogeg" following words can be formed from 's': egg, go, ego.
	 * 
	 * Create permutation of the string and lookup each permutation in the
	 * dictionary (represented here by a map)
	 */

	public void permutate(String[] words, int depth, String permutation) {
		if (depth == words.length) {
			System.out.println(permutation);
		} else {
			String w = words[depth];
			for (int i = 0; i < w.length(); i++) {
				permutate(words, depth + 1, permutation + w.charAt(i));
			}
		}
	}

	public void print(String[][] strings) {
		List<String> combinations = new ArrayList<String>();
		for (int i = 0; i < strings.length; i++) {
			print(strings, combinations, i, 0, strings[i].length - 1);
		}
		for (String str : combinations) {
			System.out.println(str);
		}
	}

	void print(String[][] strs, List<String> results, int index, int start,
			int end) {
		if (strs == null || start < 0 || start > end)
			return;
		List<String[]> adjs = adjacents(strs, index);
		if (adjs == null || adjs.size() == 0)
			return;
		for (int i = 0; i < strs[index].length; i++) {
			StringBuilder builder = new StringBuilder(strs[index][i])
					.append(" ");
			for (String[] array : adjs) {
				for (int j = 0; j < array.length; j++) {
					builder.append(array[j]).append(" ");
					if (results.indexOf(builder.toString()) >= 0)
						continue;
					results.add(builder.toString());
				}
				print(strs, results, index + 1, start, end);
			}
		}
	}

	List<String[]> adjacents(String[][] strs, int i) {
		List<String[]> adjs = new ArrayList<String[]>();
		if (strs == null || i < 0 || i >= strs.length)
			return null;
		for (int j = 0; j < strs.length; j++) {
			if (j == i)
				continue;
			adjs.add(strs[j]);
		}
		return adjs;
	}

	/**
	 * Given a sorted doubly linked list, write an algorithm to convert this
	 * list into a binary search tree (BST). The BST will be used to represent a
	 * set. You may expect that client will use it to search for presence of a
	 * key in the set. You may assume that you are given the above Node
	 * implementation. Your algorithm is not allowed to create any other
	 * instance of the Node.
	 */

	public static void main(String args[]) {
		GoogleInterviews gg = new GoogleInterviews();
		int n = 15;
		int val = gg.perfectSquares(n);
		System.out.println(val);

		System.out.println(gg.seqUniqChar("aabacbeabbed", 3));

		Node<Character> root = gg.build("ABCDEFGHI", "FBADCEGIH");
		TreeTraversal.preOrder(root);
		System.out.println();

		gg.allSubsets(4);

		int[] freq = gg.freqOuput(4, new int[] { 0, 0, 100, 3, 5, 4, 6, 4, 2,
				100, 2, 100 });
		System.out.println(print(freq));

		int numbers[] = { 2, 3, 5, 4, 1 };
		List<Integer> results = gg.subsetSum(numbers, 13);
		System.out.println(results);

		System.out.println(gg.shortestPalindrome("110"));

		int array[] = { 11, 8, 6, 7, 5, 10, 13, 12, 14 };
		for (int num : array) {
			gg.insertBST(num);
		}

		System.out.println(gg.maxElement(gg.rTree));
		int waive[] = gg.waveSort(array);

		System.out.print("[ ");
		for (int i = 0; i < array.length; i++) {
			if (i < array.length - 1) {
				System.out.print(waive[i] + ", ");
			} else {
				System.out.println(waive[i] + " ]");
			}
		}

		array = new int[] { 11, -6, 10, 14, 7, 8, 10 };
		System.out.println(gg.minSeq(array, 15));

		System.out.println(gg.anagrams("xyz"));

		System.out.println(gg.isSubstr("xyz", "abcdfyxzdflks"));

		System.out.println(gg.hasAnagramSubstring("abcdfyxzdflks", "yzx"));

		int[][] board = { { 0, 0, 0, 0, 0 }, { 0, 1, 1, 1, 0 },
				{ 0, 1, 1, 0, 0 }, { 0, 1, 0, 0, 0 }, { 0, 0, 0, 0, 0 } };
		List<int[]> minRec = gg.rectangle(board);
		for (int[] point : minRec) {
			System.out.print("[ ");
			for (int i = 0; i < point.length; i++) {
				if (i < point.length - 1) {
					System.out.print(point[i] + ", ");
				} else {
					System.out.println(point[i] + " ]");
				}
			}
		}

		List<List<Integer>> rslt = gg.sumInRange(
				new int[] { 5, 3, 9, 7, 4, 1 }, 5, 13);

		if (rslt != null) {
			for (List<Integer> list : rslt) {
				System.out.println(list);
			}
		}

		gg.permutate(new String[] { "red", "fox", "super" }, 0, "");

		List<String[]> tickets = new ArrayList<String[]>();
		tickets.add(new String[] { "MUC", "LHR" });
		tickets.add(new String[] { "CDG", "MUC" });
		tickets.add(new String[] { "SFO", "SJC" });
		tickets.add(new String[] { "LHR", "SFO" });

		String itinerary[] = gg.itinerary(tickets);
		for (String city : itinerary) {
			System.out.print(city + " ");
		}
		System.out.println();

		String words[] = new String[] { "tab", "cat", "bat" };
		System.out.println(gg.palindrome(words));

		String str = "CareerCup";
		gg.i18n(str, 1, str.length() - 1);

		String a = "8402", b = "9201";
		int l = 0;
		System.out.println(gg.streamSub(a, b, l, 0));

		TreeNode<String> t1 = new TreeNode<String>("A");
		t1.right = new TreeNode<String>("B");
		t1.right.right = new TreeNode<String>("C");

		TreeNode<String> t2 = new TreeNode<String>("B");
		t2.left = new TreeNode<String>("A");
		t2.right = new TreeNode<String>("C");

		System.out.println(gg.sameInOrder(t1, t2));

		int[][] pos = { { 5, 9, 17 }, { 4, 13, 18 }, { 16, 31, 32 } };
		System.out.println("Least Distance: " + gg.leastDistance(pos));
		System.out.println("Sum Num Less: " + gg.sumNumLessThan(10));

		array = new int[] { 1, 1, 2, 3, 3, 3, 4, 5 };
		System.out.println("Occurrence: " + gg.occurrence(array, 3));

		double[] ar = new double[] { -12.3, -12.9, -13.7, -14.5 };
		double[] range = gg.rangeWithMaxNum(ar);
		System.out.print("Range: [ ");
		for (int i = 0; i < range.length; i++) {
			System.out.print(range[i] + " ");
		}
		System.out.println("]");

		LLinkedList list = gg.new LLinkedList();
		list.add("(1)");
		list.add("(2)");
		list.add("(3(4))");
		list.add("(5)6");
		LLinkedList.Iterator iter = list.iterator();
		while (iter.hasNext()) {
			System.out.print(iter.next() + " ");
		}

		int matrix[] = new int[] { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
		int lengths[] = new int[] { 3, 3, 3 };
		matrix = gg.reverse(matrix, lengths);

		System.out.println();
		for (int i = 0; i < matrix.length; i++) {
			System.out.print(matrix[i] + " ");
		}

		System.out.println();
		int[] heap = new int[] { 16, 14, 10, 8, 7, 9, 3, 2, 4, 1 };
		System.out.println(gg.kThLargestHeapElement(heap, 4)); // TODO

		double sArray[] = new double[] { 1, 2, 3, 5, 9, 11 };
		System.out.println("Count: " + gg.countSum(sArray, 10));

		int[] nums = { 1, 3, 4, 5, 7, 8, 9 };
		System.out.println("Num nodes between range: "
				+ gg.countNodesBetweenRange(gg.buildBST(nums), 6, 8));

		File file = gg.serialize(gg.build(nums));
		gg.deserializeBFS(file);

		String strs[][] = { { "quick", "lazy" }, { "brown", "black", "grey" },
				{ "fox", "dog" } };
		gg.print(strs);
	}

	public Node<Integer> buildBST(int[] array) {
		Node<Integer> root = buildBST(array, 0, array.length - 1);
		return root;
	}

	public Node<Integer> buildBST(int[] array, int start, int end) {
		if (array == null || start > end) {
			return null;
		}
		int middle = start + (end - start) / 2;
		Node<Integer> node = new Node<Integer>(array[middle]);
		node.left = buildBST(array, start, middle - 1);
		node.right = buildBST(array, middle + 1, end);
		return node;
	}

	/**
	 * You are asked to transfer a structure of multiple trees from machine A to
	 * machine B. You could assume the structure is like this:
	 * 
	 * <pre>
	 * class TreeNode {
	 * 	String val;
	 * 	List&lt;TreeNode&gt; sons;
	 * }
	 * </pre>
	 * 
	 * In machine A, you will be given the pointer to the root node, and in
	 * machine B you should return the pointer to root node.
	 */
	class KdTreeNode {
		String val;
		List<KdTreeNode> sons;

		public KdTreeNode(String val) {
			this.val = val;
			this.sons = new ArrayList<KdTreeNode>();
		}
	}

	/**
	 * @see http://www.geeksforgeeks.org/serialize-deserialize-n-ary-tree/
	 */

	public File serialize(KdTreeNode root) {
		if (root == null)
			return null;
		String str = serializeBFS(root, new StringBuilder());
		File file = new File("/tmp/tree.csv");
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(file));
			writer.write(str);
			System.out.println(str);
		} catch (IOException e) {
			System.out.println(e.getLocalizedMessage());
		} finally {
			try {
				writer.close();
			} catch (IOException e) {
				System.out.println(e.getLocalizedMessage());
			}
		}
		return file;
	}

	private String serializeBFS(KdTreeNode node, StringBuilder str) {
		if (node == null)
			return null;
		Queue<KdTreeNode> q = new Queue<KdTreeNode>();
		q.enqueue(node);
		while (!q.isEmpty()) {
			node = q.dequeue().val;
			str.append(node.val).append(":");
			if (node.sons == null || node.sons.size() == 0) {
				str.append("\n");
				continue;
			}
			for (int i = 0; i < node.sons.size(); i++) {
				KdTreeNode son = node.sons.get(i);
				if (son == null)
					continue;
				if (i < node.sons.size() - 1)
					str.append(son.val).append(",");
				else
					str.append(son.val).append("\n");
				q.enqueue(son);
			}
		}
		return str.toString();
	}

	public KdTreeNode deserializeBFS(File file) {
		BufferedReader br = null;
		KdTreeNode root = null;
		boolean firstLine = true;
		Map<String, KdTreeNode> map = new HashMap<String, KdTreeNode>();
		try {
			br = new BufferedReader(new FileReader(file));
			String line = br.readLine();
			while (line != null) {
				if (line.isEmpty())
					continue;
				int pos = line.indexOf(":");
				String parentVal = line.substring(0, pos);
				KdTreeNode parent = map.get(parentVal);
				if (parent == null)
					parent = new KdTreeNode(parentVal);
				if (firstLine) {
					root = parent;
					firstLine = false;
				}
				String[] tokens = line.substring(pos + 1).split(",");
				for (String tok : tokens) {
					KdTreeNode son = new KdTreeNode(tok);
					parent.sons.add(son);
					map.put(tok, son);
				}
				line = br.readLine();
			}
		} catch (FileNotFoundException e) {
		} catch (IOException e) {
		} finally {
			try {
				br.close();
			} catch (IOException e) {
			}
		}
		return root;
	}

	public KdTreeNode build(int nums[]) {
		return build(nums, 0, nums.length - 1);
	}

	private KdTreeNode build(int nums[], int i, int j) {
		if (j < i)
			return null;
		int middle = i + (j - i) / 2;
		String val = String.valueOf(nums[middle]);
		KdTreeNode node = new KdTreeNode(val);
		node.sons.add(build(nums, i, middle - 1));
		node.sons.add(build(nums, middle + 1, j));
		return node;
	}

	/**
	 * We have words and their positions in a paragraph in sorted order. Write
	 * an algorithm to find the least distance for a given 3 words.
	 * 
	 * <pre>
	 * Example: job: 5, 9, 17, in: 4, 13, 18, Google: 16, 19, 20.
	 * Can you extend it to "n" words ?
	 * </pre>
	 * 
	 * Context: In Google search results, the search terms are highlighted in
	 * the sort paragraph that shows up. We need to find the shortest sentence
	 * that has all the words if we have the positions as mentioned above.
	 */
	public int leastDistance(int[][] pos) {
		int i = 0, min = Integer.MAX_VALUE;
		for (int j = 0; j < pos[i].length; j++) {
			min = Math.min(min, DFS(i, j, pos));
		}
		return min;
	}

	private int DFS(int i, int j, int[][] pos) {
		Queue<Integer> q = new Queue<Integer>();
		q.enqueue(pos[i][j]);
		int min = Integer.MAX_VALUE;
		while (!q.isEmpty()) {
			int val = q.dequeue().val;
			for (int k = j; k < pos[i + 1].length; k++) {
				q.enqueue(pos[i + 1][k]);
				int w = DFS_Visit(i + 1, k, val, pos);
				min = Math.min(min, w);
			}
			i++;
			if (i >= pos.length - 1) {
				break;
			}
		}
		return min;
	}

	private int DFS_Visit(int i, int j, int val, int[][] pos) {
		if (i >= pos.length) {
			return val;
		}
		val = Math.abs(val - pos[i][j]);
		return DFS_Visit(i + 1, j, val, pos);
	}

	/**
	 * Write an algorithm to find sum of numbers which are smaller than N and
	 * divisible by 3 or 5.
	 * 
	 * <pre>
	 * Example:
	 * N =  9 => 3 + 5 + 6 = 14
	 * N = 10 => 3 + 5 + 6 + 9 = 23
	 * </pre>
	 */
	public int sumNumLessThan(int n) {
		if (n == 0)
			return 0;
		int i = 0, sum = 0;
		while (i < n) {
			if (i % 3 == 0 || i % 5 == 0) {
				sum += i;
			}
			i++;
		}
		return sum;
	}

	/**
	 * Given an array that is rotated n times find the number where the peak
	 * happens. The array is sorted in increasing order. Follow up question: how
	 * will you rearrange the values in normal sorted order.
	 * 
	 * <pre>
	 * Example: Sorted: 2, 7, 32, 48, 55
	 *  2 pivot: 32, 48, 55, 2, 7
	 * 55 pivot: 2, 7, 55, 32, 48
	 *  7 pivot: 55, 48, 32, 7, 2
	 * </pre>
	 */
	public int peakRotated(int[] array) {
		return 0;
	}

	/**
	 * Given a max-heap represented as an array, return the k-th largest element
	 * without modifying the heap. Do it in linear or logarithmic time. Do-able
	 * in O(n) using quick-select (TODO
	 * http://en.wikipedia.org/wiki/Quickselect)
	 * 
	 * <pre>
	 * i = 1, l = 14, min = 14, s = 10
	 * i = 2, l =  8, min = 10, s = 7
	 * i = 3, l =  9, min =  9, s = 3
	 * i = 4, l =  4, min =  4, s = 2
	 * </pre>
	 * 
	 * TODO
	 */
	public int kThLargestHeapElement(int[] heap, int k) {
		if (k == 0)
			return heap[0];
		if (k == 1) {
			return Math.max(heap[left(k / 2)], heap[right(k / 2)]);
		}
		int i = 1, larger = Integer.MIN_VALUE;
		int smaller = Integer.MAX_VALUE, min = Integer.MAX_VALUE;
		while (i <= (heap.length - 1) / 2) {
			int maxL = Math.max(heap[left(i)], heap[right(i)]);
			larger = Math.max(maxL, larger);
			min = Math.min(larger, min);
			int minL = Math.min(heap[left(i)], heap[right(i)]);
			smaller = Math.min(smaller, minL);
			larger = Math.max(smaller, min);
			i++;
		}
		return min;
	}

	private int left(int i) {
		return 2 * i;
	}

	private int right(int i) {
		return 2 * i + 1;
	}

	/**
	 * There is a matrix. Each element in the matrix is a bit (0 or 1). Write a
	 * method to inverse the matrix. The matrix is stored in a one dimensional
	 * char array. The length of each row is given. How do you improve your
	 * solution when handling large amount of data ?
	 */
	public int[] reverse(int[] matrix, int[] lengths) {
		int j = 0, l = 0;
		for (int i = 0; i < matrix.length; i += l) {
			reverse(matrix, i, i + lengths[j]);
			l = lengths[j++];
		}
		return matrix;
	}

	private void reverse(int[] matrix, int i, int j) {
		if (matrix == null)
			return;
		int l = j - 1;
		for (int k = i; k < j; k++) {
			if (l == k)
				break;
			int temp = matrix[k];
			matrix[k] = matrix[l];
			matrix[l] = temp;
			l--;
		}
	}

	/**
	 * Program an iterator for a linked list which may include nodes which are
	 * nested within other nodes i.e. (1) -> (2) -> (3(4)) -> ((5)6). Iterator
	 * should return 1 -> 2 -> 3 -> 4 -> 5 -> 6.
	 */
	public class LLinkedList {
		QNode<String> current = null;
		Iterator iter = null;
		QNode<String> head = null, tail = null;

		public void add(String val) {
			QNode<String> node = new QNode<String>(val);
			if (head == null) {
				head = node;
				tail = node;
				head.next = tail;
			} else {
				tail.next = node;
				tail = node;
			}
			current = head;
		}

		public Iterator iterator() {
			if (iter == null) {
				iter = new Iterator();
			}
			return iter;
		}

		private class Iterator {

			Queue<String> queue = new Queue<String>();

			public Iterator() {
			}

			public boolean hasNext() {
				return current != null;
			}

			public String next() {
				if (queue.isEmpty()) {
					if (current == null)
						return null;
					String val = current.val;
					for (char c : val.toCharArray()) {
						if (c == '(' || c == ')')
							continue;
						queue.enqueue("" + c);
					}
					current = current.next;
				}
				return queue.dequeue().val;
			}
		}
	}

	/**
	 * You are given an array, divide it into two halves such that the sum of
	 * those 2 halves are equal.
	 */
	public int halves(int[] array) {
		int sum = 0;
		for (int i = 0; i < array.length; i++) {
			sum += array[i];
		}
		int curr = 0, i = 0;
		for (i = 0; i < array.length; i++) {
			curr += array[i];
			if (curr >= sum / 2) {
				return i;
			}
		}
		return -1;
	}

	/**
	 * Given a range of floats, find the range of size 1 with maximum number of
	 * elements. Example: [-13.7, -14.5, -12.3, -12.9], one possible answer
	 * [-12.3, -13.3]
	 * 
	 * <pre>
	 * Sorted array [-12.3, -12.9, -13.7, -14.5]
	 * </pre>
	 */
	public double[] rangeWithMaxNum(double[] range) {
		MergeSort sort = new MergeSort(range);
		sort.sort(Order.DESC);
		double max = Integer.MIN_VALUE;
		int l = range.length;
		double[] intv = new double[2];
		for (int i = 0; i < l; i++) {
			double n = range[i] - 1;
			int pos = binaryInsert(range, i, l - 1, n);
			if (pos >= 0 && (pos - i) > max) {
				max = pos - i;
				intv = new double[] { n, range[i] };
			}
		}
		return intv;
	}

	public int binarySearch(double[] array, int i, int j, double n) {
		if (i <= j) {
			int middle = i + (j - i) / 2;
			if (array[middle] == n)
				return middle;
			if (array[middle] < n)
				return binarySearch(array, i, middle - 1, n);
			return binarySearch(array, middle + 1, j, n);
		}
		return -1;
	}

	public int binaryInsert(double[] array, int i, int j, double n) {
		int curr = 0;
		while (true) {
			curr = i + (j - i) / 2;
			if (array[curr] == n) {
				return curr;
			}
			if (array[curr] < n) {
				i = curr + 1;
				if (i > j)
					return curr + 1;
			} else {
				j = curr - 1;
				if (i > j)
					return curr;
			}
		}
	}

	/**
	 * Given a sorted array and n, find the count of sum of 2 numbers greater
	 * than or equal to n.
	 */
	public int countSum(double[] array, double n) {
		if (array == null)
			throw new IllegalArgumentException();
		int length = array.length, count = 0;
		for (int i = 0; i < length - 1 && i <= n; i++) {
			double diff = n - array[i];
			int pos = binaryInsert(array, i + 1, length, diff);
			count += (length - pos);
		}
		return count;
	}

	/**
	 * Given a sorted array, find a way to find the number of occurrence of a
	 * number. Example: [1, 1, 2, 3, 3, 3, 4, 5]. Do it (let us say 3) in a time
	 * better than O(n).
	 */
	public int occurrence(int[] array, int n) {
		if (array == null)
			throw new IllegalArgumentException();
		int l = array.length - 1;
		return binary(array, 0, l, 0, n);
	}

	private int binary(int[] array, int i, int j, int val, int n) {
		if (i <= j) {
			int middle = (i + j) / 2;
			if (array[middle] == n)
				val++;
			val = binary(array, i, middle - 1, val, n);
			val = binary(array, middle + 1, j, val, n);
		}
		return val;
	}

	private static class TreeNode<T> {
		public T value;
		public TreeNode<T> left;
		public TreeNode<T> right;

		public TreeNode(T value) {
			this.value = value;
			left = null;
			right = null;
		}
	}

	private Node<Integer> rTree = null;

	/**
	 * You are given two streams 's_a' and 's_b' of digits. Each stream
	 * represents an integer 'a' and 'b' from its less significant digit to the
	 * most significant digit. For example, integer 2048 is represented in the
	 * stream as 8, 4, 0, 2.
	 * 
	 * Write a function that subtracts two integers 'a - b' and returns the
	 * result as string. You are not allowed to buffer entire 'a' or 'b' into a
	 * single string, i.e., you may access only a single digit per stream at
	 * time (imagine that 'a' and 'b' are stored in huge files). You may assume
	 * that 'a >= b'. Example: s_a: 8 4 0 2, s_b: 4 2 0 1, result 1024.
	 * 
	 */
	public String streamSub(String a, String b, int i, int r) {
		if (i >= a.length()) {
			return "";
		}
		int aval = Integer.parseInt("" + a.charAt(i));
		int bval = 0;
		if (i < b.length()) {
			bval = Integer.parseInt("" + b.charAt(i));
		}
		bval += r;
		r = aval < bval ? 1 : 0;
		aval = aval < bval ? 10 + aval : aval;
		int diff = aval - bval;
		return streamSub(a, b, i + 1, r) + diff;
	}

	/**
	 * i18n (where 18 stands for the number of letters between the first i and
	 * the last n in the word "internationalization"). Generate all such
	 * possible i18n strings for any given string. e.g. "CareerCup" = > "c7p",
	 * "ca6p", "c6up", "car5p", "ca5up", "care4p", "car4up", "care4p", "car4up",
	 * "caree3p", "care3up", career2p, car4cup, careerC1p, ca5rcup, careerCup,
	 * c6arcup... till the count is 0 which means its the complete string again.
	 * Recursion.
	 * 
	 * <pre>
	 * c7p, ca6p, car5p, care4p, caree3p, career2p, careerc1p, 
	 * c6up, ca5up, car4up, care3up, caree2up, career1up
	 * c5cup, ca4cup, car3cup, care2cup, caree1cup,
	 * c4rcup, ca3rcup, car2rcup, care1rcup,
	 * c3ercup, ca2ercup, car1ercup,
	 * c2eercup, ca1eercup
	 * c1reercup,
	 * </pre>
	 */
	public void i18n(String str, int i, int j) {
		if (i > j) {
			System.out.println(str);
			return;
		} else if (j > i) {
			String n = str.substring(0, i);
			String l = str.substring(j);
			String s = str.substring(i, j);
			String r = n + s.length() + l;
			System.out.println(r);
			i18n(str, i, j - 1);
		} else {
			j = str.length() - 1;
			i18n(str, i + 1, j);
		}
	}

	/**
	 * You are given a list of words. Find if two words can be joined together
	 * to form a palindrome e.g., consider a list {bat, tab, cat}, then bat and
	 * tab can be joined together to form a palindrome.
	 */
	public boolean palindrome(String[] words) {
		if (words == null)
			return false;
		Map<String, Integer> strs = new HashMap<String, Integer>();
		for (int i = 0; i < words.length; i++) {
			strs.put(words[i], i);
		}
		for (int i = 0; i < words.length; i++) {
			if (strs.containsKey(reverse(words[i]))) {
				return true;
			}
		}
		return false;
	}

	private String reverse(String str) {
		if (str == null || str == "") {
			return str;
		}
		int i = str.length() - 1;
		StringBuilder reverse = new StringBuilder();
		while (i >= 0) {
			reverse.append(str.charAt(i--));
		}
		return reverse.toString();
	}

	/**
	 * Given a bunch of airline of tickets with [from, to], for example [MUC,
	 * LHR], [CDG, MUC], [SFO, SJC], [LHR, SFO], please reconstruct the
	 * itinerary in order, for example [CDG, MUC, LHR, SFO, SJC];
	 * 
	 * 
	 * MUC:0, LHR:1, CDG:2, SFO:3, SJC:4
	 * 
	 * MUC:2, LHR:1, CDG:0, SFO:3, SJC:4
	 * 
	 * MUC:1, LHR:2, CDG:0, SFO:3, SJC:4
	 */

	public String[] itinerary(List<String[]> tickets) {
		if (tickets == null || tickets.size() == 0)
			return null;
		Map<String, Integer> map = new HashMap<String, Integer>();
		int index = 0;
		for (String ticks[] : tickets) {
			if (!map.containsKey(ticks[0]) && !map.containsKey(ticks[1])) {
				map.put(ticks[0], index++);
				map.put(ticks[1], index++);

			} else {
				if (!map.containsKey(ticks[0])) {
					int pos = map.get(ticks[1]);
					map.put(ticks[0], pos);
					map.put(ticks[1], index++);
				} else if (!map.containsKey(ticks[1])) {
					map.put(ticks[1], index++);
				} else {
					int lhs = map.get(ticks[0]);
					int rhs = map.get(ticks[1]);
					if (lhs > rhs) {
						map.put(ticks[0], rhs);
						map.put(ticks[1], lhs);
					}
				}
			}
		}

		for (String ticks[] : tickets) {
			int lhs = map.get(ticks[0]);
			int rhs = map.get(ticks[1]);
			if (lhs > rhs) {
				map.put(ticks[0], rhs);
				map.put(ticks[1], lhs);
			}
		}

		String[] itinerary = new String[tickets.size() + 1];

		for (String city : map.keySet()) {
			itinerary[map.get(city)] = city;
		}
		return itinerary;
	}

	/**
	 * If a canoe can hold 2 kids and a max weight of 150 LBS, write a function
	 * that returns the minimum number of canoe needed, given a list of kids and
	 * their weights.
	 * 
	 * 
	 * 2 5 3 7 9 [2, 5], [3, 7], [9]
	 */
	public int minCanoe(int[] weights, int maxWeight) {
		if (weights == null)
			throw new IllegalArgumentException();
		int i = 0, j = 0, sum = 0, num = 0;
		while (i < weights.length) {
			sum += weights[j];
			if (sum < maxWeight) {
				j++;
				if (j % 2 == 0) {
					i = j;
					sum = 0;
					num++;
				}
			} else {
				i++;
				sum = 0;
				num++;
			}
		}
		return num;
	}

	/**
	 * Given two binary trees (not BST), return true if both trees have the same
	 * in-order else return false. We can call in-order on both trees and
	 * compare their results, but we want to do this in parallel. // TODO Test
	 * this method
	 */
	public <T> boolean sameInOrder(TreeNode<T> t1, TreeNode<T> t2) {
		TreeStack<T> t1s = new TreeStack<T>(t1);
		TreeStack<T> t2s = new TreeStack<T>(t2);
		while (t1s.hasNext() && t2s.hasNext()) {
			TreeNode<T> node1 = t1s.next();
			TreeNode<T> node2 = t2s.next();
			if (node1.value != node2.value) {
				return false;
			}
		}
		return !(t1s.hasNext() || t2s.hasNext());
	}

	private class TreeStack<T> {

		private Stack<TreeNode<T>> stack;

		public TreeStack(TreeNode<T> node) {
			stack = new Stack<TreeNode<T>>();
			push(node);
		}

		private void push(TreeNode<T> node) {
			if (node == null) {
				return;
			}
			while (node != null) {
				stack.push(node);
				node = node.left;
			}
		}

		public TreeNode<T> next() {
			if (stack.isEmpty())
				return null;
			TreeNode<T> toReturn = stack.peek();
			TreeNode<T> node = stack.pop();
			node = node.right;
			while (node != null) {
				stack.push(node);
				node = node.left;
			}
			return toReturn;
		}

		public boolean hasNext() {
			return !stack.isEmpty();
		}
	}

	/**
	 * Given an integer, find the next highest and next lowest integers, with
	 * equal number of 1s in their binary representations as the original number
	 */

	/**
	 * Given a list of integers that fall within a known short but unknown range
	 * of values, how to find the median value ? TODO: Use selection sort (O(n))
	 * since the range is known to be small.
	 */

	/**
	 * How do you add two numbers that are larger/longer than integer type ?
	 */

	/**
	 * How would you get all of the unique IDs out of a single file ? What if
	 * the file was a very large file ?
	 */

	/**
	 * In GMAIL, while composing an email, upon adding a contact, related
	 * contacts are displayed. How would you implement that feature ? Write an
	 * algorithm for that. What data structure would you use to store the
	 * weights ? In what format, would you persist this data ?
	 */

	/**
	 * You have a sorted array containing the age of every person on earth. [0,
	 * 0, 0, 0, 1, 1, ..., 28, 28, ..., 110, ...]. Find out how many people have
	 * each age.
	 */

	/**
	 * What is the maximum number of edges you could add to n vertices to make
	 * an acyclic undirected graph ? Follow up: What is the maximum number of
	 * edges you could add to n vertices to make an acyclic direct graph ?
	 */

	/**
	 * Given a dictionary, and a list of letters (or consider as a string), find
	 * the longest word that only uses letters from the string.
	 */

	/**
	 * Given an array of numbers, find a pair whose sum is closest to zero.
	 * Solution: Sort array, do binary search to figure out position of zero,
	 * and return either adjacent values or preceding 2 values or following 2
	 * values. TODO disregard the hint
	 */

	// GlassDoor Interviews

	/**
	 * Given a binary tree where there may be duplicates, but all other logic of
	 * the BST is intact, determine the most frequently occurring element.
	 */

	/**
	 * You are given an array of integers and a number 'n'. List all sequences
	 * of that array that sum to 'n'. Example: {1, 2, 2, 3, 4, 5}, sum = 5;
	 * Results: {1, 2, 2}, {1, 4}, {2, 3}.
	 */

	/**
	 * Given a bunch of points and a convex polygon in R^2, determine if each
	 * point is inside the polygon.
	 */

	/**
	 * How to read last 100 lines in a billion lines file where each line has a
	 * different length ?
	 */

	/**
	 * Given a string 's', find the minimum cuts that partition 's' into
	 * substrings which are all palindrome.
	 */

	/**
	 * How to search a heap with recursion ?
	 */

	/**
	 * 1.) Find the k-th largest element of a n*n matrix, given every row or
	 * column of the matrix is pre-sorted in O(n*log n) time. 2.) Find the k-th
	 * largest element of two pre-sorted array in O(log n).
	 */
}
