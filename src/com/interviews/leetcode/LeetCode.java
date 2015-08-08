package com.interviews.leetcode;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;

import com.interviews.commons.Node;
import com.interviews.google.MergeSort;
import com.interviews.google.MergeSort.Order;

/**
 * 
 * http://n00tc0d3r.blogspot.ca/2013/05/recover-binary-search-tree.html
 *
 */

public class LeetCode {

	static int postIndex = 0, preIndex = 0;

	/**
	 * Word Search: Given a 2D board and a word, find if the word exists in the
	 * grid. The word can be constructed from letters of sequentially adjacent
	 * cells, where "adjacent" cells are those horizontally or vertically
	 * neighboring. The same cell may not be used more than once.
	 * 
	 * <pre>
	 * Example:
	 * 
	 * [ 
	 * 	  ["ABCE"],
	 * 	  ["SFCS"],
	 * 	  ["ADEE"]
	 * ]
	 * word = "ABCCDE" -> returns true,
	 * word = "SEE"    -> returns true,
	 * word = "ABCB"   -> returns false
	 * </pre>
	 */
	public boolean exists(char[][] grid, String word) {
		if (grid == null || grid.length == 0 || word == null
				|| word.length() == 0) {
			return false;
		}
		boolean visited[][] = new boolean[grid.length][grid.length];
		for (int i = 0; i < grid.length; i++) {
			for (int j = 0; j < grid[i].length; j++) {
				if (match(grid, word, visited, 0, i, j)) {
					return true;
				}
			}
		}
		return false;
	}

	private boolean match(char[][] grid, String word, boolean[][] visited,
			int i, int r, int c) {
		if (grid[r][c] != word.charAt(i) || visited[r][c]) {
			return false;
		}
		visited[r][c] = true;
		if (i == word.length() - 1) {
			return true;
		}
		if (r - 1 >= 0 && match(grid, word, visited, i + 1, r - 1, c)) {
			return true;
		}
		if (r + 1 < grid.length && match(grid, word, visited, i + 1, r + 1, c)) {
			return true;
		}
		if (c - 1 >= 0 && match(grid, word, visited, i + 1, r, c - 1)) {
			return true;
		}
		if (c + 1 < grid.length && match(grid, word, visited, i + 1, r, c + 1)) {
			return true;
		}
		visited[r][c] = false;
		return false;
	}

	/**
	 * Modified version of the previous question (holds if rows anterior to i
	 * are not to be considered
	 */
	public boolean isPresent(char[][] grid, String word) {
		int i = 0, j = 0, k = 0;
		boolean visited[] = new boolean[grid.length];
		while (i < grid.length) {
			if (k == word.length()) {
				return true;
			}
			if (j < grid.length && grid[i][j] == word.charAt(k) && !visited[j]) {
				k++;
				visited[j++] = true;
			} else {
				j = 0;
				i++;
				if (i < grid.length) {
					while (j < grid.length && grid[i][j] != word.charAt(k)) {
						j++;
					}
				}
				reset(visited);
			}
		}
		return false;
	}

	private void reset(boolean[] visited) {
		for (int i = 0; i < visited.length; i++) {
			visited[i] = false;
		}
	}

	/**
	 * Word Search II: Given a 2D board and a list of words from the dictionary,
	 * find all words in the board.
	 * 
	 * Each word must be constructed from letters of sequentially adjacent cell,
	 * where "adjacent" cells are those horizontally or vertically neighboring.
	 * The same letter cell may not be used more than once in a word.
	 * 
	 * <pre>
	 * For example,
	 * Given words = ["oath","pea","eat","rain"] and board =
	 * 
	 * [
	 *   ['o','a','a','n'],
	 *   ['e','t','a','e'],
	 *   ['i','h','k','r'],
	 *   ['i','f','l','v']
	 * ]
	 * 
	 * Return ["eat","oath"].
	 * </pre>
	 */
	public List<String> findWords(char[][] board, String[] words) {
		ArrayList<String> result = new ArrayList<String>();

		int m = board.length;
		int n = board[0].length;

		for (String word : words) {
			boolean flag = false;
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++) {
					char[][] newBoard = new char[m][n];
					for (int x = 0; x < m; x++)
						for (int y = 0; y < n; y++)
							newBoard[x][y] = board[x][y];

					if (dfs(newBoard, word, i, j, 0)) {
						flag = true;
					}
				}
			}
			if (flag) {
				result.add(word);
			}
		}

		return result;
	}

	public boolean dfs(char[][] board, String word, int i, int j, int k) {
		int m = board.length;
		int n = board[0].length;

		if (i < 0 || j < 0 || i >= m || j >= n || k > word.length() - 1) {
			return false;
		}

		if (board[i][j] == word.charAt(k)) {
			char temp = board[i][j];
			board[i][j] = '#';

			if (k == word.length() - 1) {
				return true;
			} else if (dfs(board, word, i - 1, j, k + 1)
					|| dfs(board, word, i + 1, j, k + 1)
					|| dfs(board, word, i, j - 1, k + 1)
					|| dfs(board, word, i, j + 1, k + 1)) {
				board[i][j] = temp;
				return true;
			}

		} else {
			return false;
		}

		return false;
	}

	/**
	 * Suppose a sorted array is rotated at some pivot unknown to you
	 * beforehand.
	 * 
	 * (i.e., 0 1 2 4 5 6 7 might become 4 5 6 7 0 1 2).
	 * 
	 * You are given a target value to search. If found in the array return its
	 * index, otherwise return -1. Assuming duplicates are allowed, write a
	 * function to determine if a given target is in the array..
	 */
	public boolean existTarget(int[] array, int target) {
		if (array == null)
			return false;
		int i = 0;
		while (i < array.length - 1) {
			if (array[i] > array[i + 1]) {
				break;
			}
			i++;
		}
		int j = i + 1;
		if (array[j] == target)
			return true;
		return bSearch(array, 0, j, target)
				|| bSearch(array, j, array.length - 1, target);
	}

	private boolean bSearch(int[] array, int lb, int ub, int elem) {
		if (lb <= ub) {
			int middle = lb + (ub - lb) / 2;
			if (array[middle] == elem) {
				return true;
			}
			if (array[middle] < elem) {
				return bSearch(array, middle + 1, ub, elem);
			}
			return bSearch(array, lb, middle - 1, elem);
		}
		return false;
	}

	/**
	 * Minimum Total: Given a triangle, find the minimum path sum from top to
	 * bottom. Each step, you may move to adjacent numbers on the row below. For
	 * example, given the following triangle
	 * 
	 * <pre>
	 * [
	 *      [2],
	 *     [3,4],
	 *    [6,5,7],
	 *   [4,1,8,3]
	 * ]
	 * </pre>
	 * 
	 * The minimum path sum from top to bottom is 11 (i.e., 2 + 3 + 5 + 1 = 11).
	 * 
	 * Note: Bonus point if you are able to do this using only O(n) extra space,
	 * where n is the total number of rows in the triangle.
	 */
	public int[] minimumTotal(int[][] triangle) {
		if (triangle == null)
			return null;
		int min = Integer.MAX_VALUE;
		int[] path = new int[triangle.length];
		for (int i = 0; i < triangle.length - 1; i++) {
			for (int j = 0; j < triangle[i + 1].length - 1; j++) {
				min = minVal(i, j, triangle, min);
				path[i + 1] = minIndex(i, j, triangle, min, path[i + 1]);
			}
			min = Integer.MAX_VALUE;
		}
		for (int i = 0; i < path.length; i++) {
			path[i] = triangle[i][path[i]];
		}
		return path;
	}

	private int minIndex(int i, int j, int[][] triangle, int min, int pos) {
		if (min == triangle[i + 1][j]) {
			return j;
		} else if (min == triangle[i + 1][j + 1]) {
			return j + 1;
		}
		return pos;
	}

	private int minVal(int i, int j, int[][] triangle, int min) {
		min = Math.min(min,
				Math.min(triangle[i + 1][j], triangle[i + 1][j + 1]));
		return min;
	}

	/**
	 * http://n00tc0d3r.blogspot.ca/2013/06/triangle.html
	 */
	public int minimumTotal(List<List<Integer>> triangle) {
		int[] preRow = new int[triangle.size()];
		// BFS
		for (List<Integer> row : triangle) {
			int last = row.size() - 1;
			if (last > 0)
				preRow[last] = preRow[last - 1] + row.get(last);
			int pre = preRow[0];
			preRow[0] += row.get(0);
			for (int i = 1; i < last; ++i) {
				// save the previous value so that we can use it for next number
				int temp = preRow[i];
				// either from i-1 or i from the previous row
				preRow[i] = Math.min(preRow[i] + row.get(i), pre + row.get(i));
				pre = temp;
			}
		}
		// find the min
		int res = preRow[0];
		for (int num : preRow) {
			if (num < res)
				res = num;
		}
		return res;
	}

	/**
	 * Search Insert Position: Given a sorted array and a target value, return
	 * the index if the target is found. If not, return the index where it would
	 * be if it were inserted in order.
	 */
	public int binaryInsert(int[] array, int val) {
		if (array == null)
			return -1;
		int length = array.length - 1;
		int pos = binaryInsert(array, 0, length, val);
		return pos;
	}

	private int binaryInsert(int[] array, int lb, int ub, int val) {
		while (lb <= ub) {
			int middle = lb + (ub - lb) / 2;
			if (array[middle] == val) {
				return middle;
			}
			if (array[middle] < val) {
				lb = middle + 1;
				if (array[lb] > val) {
					return lb;
				}
			} else {
				ub = middle - 1;
				if (array[ub] < val) {
					return middle;
				}
			}
		}
		return -1;
	}

	/**
	 * Search for a range: Given a sorted array of integers, find the starting
	 * and ending position of a given target value.
	 * 
	 * Your algorithm's runtime complexity must be in the order of O(log n).
	 * 
	 * If the target is not found in the array, return [-1, -1].
	 * 
	 * <pre>
	 * For example,
	 * Given [5, 7, 7, 8, 8, 10] and target value 8,
	 * return [3, 4].
	 * </pre>
	 */
	public int[] searchRange(int[] array, int val) { //
		if (array == null)
			return null;
		MinMax instance = new MinMax();
		int length = array.length - 1;
		binSearch(array, 0, length, val, instance);
		int min = instance.min == Integer.MAX_VALUE ? -1 : instance.min;
		int max = instance.max == Integer.MIN_VALUE ? -1 : instance.max;
		return new int[] { min, max };
	}

	private void binSearch(int[] array, int lb, int ub, int val, MinMax inst) {
		if (lb <= ub) {
			int middle = lb + (ub - lb) / 2;
			if (array[middle] == val) {
				inst.set(middle);
			}
			if (array[middle] < val) {
				binSearch(array, middle + 1, ub, val, inst);
			} else {
				binSearch(array, lb, middle - 1, val, inst);
			}
		}
	}

	private class MinMax {
		int min = Integer.MAX_VALUE;
		int max = Integer.MIN_VALUE;

		public void set(int val) {
			min = Math.min(min, val);
			max = Math.max(max, val);
		}
	}

	/**
	 * Write an efficient algorithm that searches for a value in an m x n
	 * matrix. This matrix has the following properties: 1. integers in each row
	 * are sorted from left to right. 2. the first integer of each row is
	 * greater than the last integer of the previous row.
	 * 
	 * <pre>
	 * Example:
	 * [
	 *   [1,   3,  5,  7],
	 *   [10, 11, 16, 20],
	 *   [23, 30, 34, 50]
	 * ]
	 * 
	 * </pre>
	 */
	public boolean find(int[][] array, int val) {
		if (array == null)
			return false;
		int i = 0, length = array.length - 1;
		for (; i <= length; i++) {
			if (array[i][length] == val) {
				return true;
			}
			if (array[i][length] < val) {
				i++;
			} else {
				break;
			}
		}
		for (int k = length; k >= 0; k--) {
			if (array[i][k] == val) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Given a sorted array, remove the duplicates in place such that each
	 * element appear only once and return the new length. Do not allocate extra
	 * space for another array, you must do this in place with constant memory.
	 * 
	 * For example, given input array [1, 1, 2], the function should return
	 * length = 2, with the first two elements of being 1 and 2 respectively. It
	 * doesn't matter what you leave beyond the new length.
	 */
	public int newLength(int[] array) {
		if (array == null)
			return -1;
		int length = array.length;
		int i = 0, j = 1, k = 0;
		while (true) {
			if (i >= length || j >= length) {
				break;
			}
			if (array[j] == array[i]) {
				j++;
			} else {
				k++;
				i = j++;
			}
		}
		return k + 1;
	}

	/**
	 * Unique Paths: A robot is located at the top-left corner of an m*n grid.
	 * The robot can only move either down or right at any point in time. The
	 * robot is trying to reach the bottom-right corner of the grid. How many
	 * possible unique paths are there.
	 */
	public void uniquePaths(int[][] grid) {
		if (grid == null)
			return;
		Queue<int[]> queue = new Queue<int[]>();
		queue.enqueue(new int[] { 0, 0 });
		Map<int[], Boolean> visited = new HashMap<int[], Boolean>();
		Map<int[], List<int[]>> adjs = new HashMap<int[], List<int[]>>();
		while (!queue.isEmpty()) {
			int val[] = queue.dequeue().val;
			visited.put(val, true);
			List<int[]> adj = adjacent(grid, val[0], val[1]);
			adjs.put(val, adj);
			if (adj != null) {
				for (int[] array : adj) {
					if (!visited.containsKey(array)) {
						queue.enqueue(array);
					}
				}
			}
		}
		List<int[]> path = new ArrayList<int[]>();
		Map<List<int[]>, Boolean> exist = new HashMap<List<int[]>, Boolean>();
		uniquePaths(new int[] { 0, 0 }, adjs, exist, path, new int[] {
				grid.length - 1, grid[0].length - 1 });
	}

	private List<int[]> adjacent(int[][] grid, int i, int j) {
		if (i >= grid.length - 1 || j >= grid[0].length - 1)
			return null;
		List<int[]> adjs = new ArrayList<int[]>();
		adjs.add(new int[] { i + 1, j });
		adjs.add(new int[] { i, j + 1 });
		return adjs;
	}

	private void uniquePaths(int[] key, Map<int[], List<int[]>> map,
			Map<List<int[]>, Boolean> exist, List<int[]> path, int[] end) {
		if (map == null || map.size() == 0)
			return;
		List<int[]> adjs = map.get(key);
		if (adjs == null) {
			if (!exist.containsKey(path) && key[0] == end[0]
					&& key[1] == end[1]) {
				exist.put(path, true);
				for (int[] k : path) {
					System.out.print("[" + k[0] + ", " + k[1] + "] ");
				}
				System.out.println();
			}
			path.remove(key);
		} else {
			for (int[] array : adjs) {
				path.add(array);
				uniquePaths(array, map, exist, path, end);
			}
		}
	}

	/**
	 * Given two words (start and end), and a dictionary, find all
	 * transformation sequence(s) from start to end, such that: 1. only one
	 * letter can be changed at a time, 2. each intermediate word must exist in
	 * the dictionary. Follow up: find all shortest transformation sequences
	 * (TODO Apply Dijkstra)
	 * 
	 * <pre>
	 * Example:
	 * 
	 * start = "hit", end = "cog"
	 * dict = ["hot", "dog", "lot", "log"]
	 * 
	 * Return 
	 * [ 
	 * 	 ["hit", "hot", "dot", "dog", "cog"]
	 * ]
	 * </pre>
	 */
	public List<List<String>> ladder(String[] dict, String start, String end) {
		if (dict == null || start == end) {
			return null;
		}
		Map<String, Boolean> visited = new HashMap<String, Boolean>();
		for (String word : dict) {
			visited.put(word, false);
		}
		Queue<String> queue = new Queue<String>();
		queue.enqueue(start);

		List<List<String>> paths = new ArrayList<List<String>>();
		while (!queue.isEmpty()) {
			String word = queue.dequeue().val;
			List<String> adj = adjacent(word, visited);
			for (int i = 0; i < adj.size(); i++) {
				List<String> path = new ArrayList<String>();
				path.add(start);
				queue.enqueue(adj.get(i));
				if (DFS_Visit(adj.get(i), end, path, visited)) {
					path.add(end);
					paths.add(path);
				}
			}
		}
		return paths;
	}

	private List<String> adjacent(String word, Map<String, Boolean> visited) {
		List<String> adj = new ArrayList<String>();
		for (Map.Entry<String, Boolean> entry : visited.entrySet()) {
			String key = entry.getKey();
			if (isOneEditWord(key, word) && !entry.getValue()) {
				adj.add(key);
			}
		}
		return adj;
	}

	private boolean DFS_Visit(String word, String end, List<String> ladder,
			Map<String, Boolean> visited) {
		if (word != end) {
			ladder.add(word);
			visited.put(word, true);
			List<String> adj = adjacent(word, visited);
			if (adj == null || adj.size() == 0) {
				return false;
			}
			for (int i = 0; i < adj.size(); i++) {
				DFS_Visit(adj.get(i), end, ladder, visited);
			}
		}
		return true;
	}

	private boolean isOneEditWord(String lhs, String rhs) {
		if (lhs == rhs) {
			return false;
		}
		if (lhs != null && rhs != null && lhs.length() != rhs.length()) {
			return false;
		}
		int count = 0;
		for (int i = 0; i < lhs.length(); i++) {
			if (lhs.charAt(i) != rhs.charAt(i)) {
				count++;
			}
		}
		return count == 1;
	}

	/**
	 * Basic Calculator: Implement a basic calculator to evaluate a simple
	 * expression string.
	 * 
	 * The expression string may contain open ( and closing parentheses ), the
	 * plus + or minus sign -, non-negative integers and empty spaces .
	 * 
	 * You may assume that the given expression is always valid.
	 * 
	 * <pre>
	 * Some examples:
	 * 
	 * "1 + 1" = 2
	 * " 2-1 + 2 " = 3
	 * "(1+(4+5+2)-3)+(6+8)" = 23
	 * </pre>
	 */
	public int calculate(String s) {
		return 0;
	}

	/**
	 * Shortest Palindrome: Given a string S, you are allowed to convert it to a
	 * palindrome by adding characters in front of it. Find and return the
	 * shortest palindrome you can find by performing this transformation.
	 * 
	 * For example:
	 * 
	 * Given "aacecaaa", return "aaacecaaa".
	 * 
	 * Given "abcd", return "dcbabcd".
	 */
	public String shortestPalindrome(String s) {
		return null;
	}

	/**
	 * Largest Number: Given a list of non negative integers, arrange them such
	 * that they form the largest number.
	 * 
	 * For example, given [3, 30, 34, 5, 9], the largest formed number is
	 * 9534330.
	 * 
	 * Note: The result may be very large, so you need to return a string
	 * instead of an integer.
	 */
	public String largestNumber(int[] nums) {
		return null;
	}

	public static void main(String args[]) {
		LeetCode lt = new LeetCode();
		char[][] grid = { { 'A', 'B', 'C', 'E' }, { 'S', 'F', 'C', 'S' },
				{ 'A', 'D', 'E', 'E' } };
		System.out.println("Exist: " + lt.exists(grid, "ABCCDE"));
		String dict[] = { "hot", "dog", "lot", "log" };
		List<List<String>> ladder = lt.ladder(dict, "hit", "cog");

		System.out.print("Word Ladder: ");
		for (List<String> path : ladder) {
			for (String str : path) {
				System.out.print(str + " ");
			}
			System.out.println();
		}

		Set<String> set = new HashSet<String>();
		for (String str : dict) {
			set.add(str);
		}

		System.out.println();
		int[][] triangle = { { 2 }, { 3, 4 }, { 6, 5, 7 }, { 4, 1, 8, 3 } };
		int path[] = lt.minimumTotal(triangle);

		System.out.print("Min Path Sum: ");
		for (int i = 0; i < path.length; i++) {
			System.out.print(path[i] + " ");
		}

		System.out.println();
		int array[] = { 4, 5, 6, 7, 0, 0, 1, 2, 2 };
		System.out.println("Search Rotated Array: " + lt.existTarget(array, 4));

		int matrix[][] = { { 0, 0, 0 }, { 0, 1, 0 }, { 0, 0, 0 } };
		lt.uniquePaths(matrix);

		array = new int[] { 1, 3, 7, 8, 10 };
		System.out.println("Position: " + lt.binaryInsert(array, 9));

		matrix = new int[][] { { 1, 3, 5, 7 }, { 10, 11, 16, 20 },
				{ 23, 30, 34, 50 } };
		System.out.println("Present: " + lt.find(matrix, 30));

		array = new int[] { 1, 3, 3, 3, 10, 11 };
		System.out.println("Length: " + lt.newLength(array));

		array = new int[] { 5, 7, 7, 8, 8, 10 };
		int[] rslt = lt.searchRange(array, 9);
		System.out.println("Range: [" + rslt[0] + ", " + rslt[1] + "]");

		array = new int[] { 1, 2, 3, 4, 5, 6, 7 };
		array = lt.rotate(array, 12);

		for (int i = 0; i < array.length; i++) {
			System.out.print(array[i] + " ");
		}

		System.out.println();
		array = new int[] { 5, 7, 7, 8, 8, 10 };
		System.out.println("New Length: " + lt.removeInstance(array, 8));

		array = new int[] { 100, 4, 200, 1, 3, 2 };
		System.out.println("Sequence Length: "
				+ lt.longestConsecutiveSequence(array));

		int lhs[] = new int[] { 1, 2, 3, 4 };
		int rhs[] = new int[] { 5, 6, 7, 9 };
		System.out.println("Median: " + lt.median(lhs, rhs));

		matrix = new int[][] { { 0, 0, 0, 0 }, { 0, 1, 1, 0 }, { 0, 1, 1, 0 } };
		List<int[]> list = lt.largestRectangle(matrix);
		for (int[] point : list) {
			System.out.print("[" + point[0] + ", " + point[1] + "] ");
		}
		System.out.println();

		System.out.print("Contiguous array sum: [ ");
		array = new int[] { -2, 1, -3, 4, -1, 2, 1, -5, 4 };
		int index[] = lt.contiguousArrayLargestSum(array);
		for (int i = index[0]; i <= index[1]; i++) {
			System.out.print(array[i] + " ");
		}
		System.out.println("]");

		System.out.print("Contiguous product array: [ ");
		array = new int[] { 2, -3, 0, 2, 4 };
		index = lt.contiguousArrayLargestProduct(array);
		for (int i = index[0]; i <= index[1]; i++) {
			System.out.print(array[i] + " ");
		}
		System.out.println("]");

		int[] preorder = new int[] { 6, 2, 1, 4, 3, 5, 7, 9, 8, 10 };
		int[] inorder = new int[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
		System.out.print("Built Tree I : ");
		Node<Integer> tree = lt.buildTree1(preorder, inorder);
		lt.printInOrder(tree);

		System.out.println();
		int[] postorder = new int[] { 1, 3, 5, 4, 2, 8, 9, 7, 6 };
		System.out.print("Built Tree II: ");
		LeetCode.postIndex = postorder.length - 1;
		lt.printInOrder(lt.buildTree2(inorder, postorder));

		System.out.println();
		int[] height = new int[] { 6, 2, 1, 4, 3, 5, 7, 9, 8, 10 };
		System.out.println("Max area: " + lt.maxArea(height));

		List<int[]> intvs = new ArrayList<int[]>();
		intvs.add(new int[] { 1, 2 });
		intvs.add(new int[] { 3, 5 });
		intvs.add(new int[] { 6, 7 });
		intvs.add(new int[] { 8, 10 });
		intvs.add(new int[] { 12, 16 });
		List<int[]> merge = lt.insert(intvs, new int[] { 4, 9 });
		System.out.print("After Insert: ");
		for (int[] intv : merge) {
			System.out.print("[" + intv[0] + ", " + intv[1] + "] ");
		}

		System.out.println();
		int nums[] = new int[] { 3, 4, -1, 1 };
		System.out.println("Missing number: " + lt.firstMissingPositive(nums));

		nums = new int[] { 1, 2, 3, 1 };
		System.out.println("Peak element: " + lt.getPeak(nums));

		System.out.print("Combination sum I: ");
		double[] candidates = new double[] { 2, 3, 6, 7 };
		List<List<Double>> rslts = lt.combinationSum1(candidates, 9);
		for (List<Double> comb : rslts) {
			System.out.print("[ ");
			for (double val : comb) {
				System.out.print(val + " ");
			}
			System.out.print("] ");
		}

		System.out.println();
		System.out.print("Combination sum III: ");
		rslts = lt.combinationSum3(3, 7);
		for (List<Double> comb : rslts) {
			System.out.print("[ ");
			for (double val : comb) {
				System.out.print(val + " ");
			}
			System.out.print("] ");
		}

		System.out.println();
		System.out.print("Combination sum II: ");
		candidates = new double[] { 10, 1, 2, 7, 6, 1, 5 };
		rslts = lt.combinationSum2(candidates, 8);
		for (List<Double> comb : rslts) {
			System.out.print("[ ");
			for (double val : comb) {
				System.out.print(val + " ");
			}
			System.out.print("] ");
		}

		System.out.println();
		int[] prices = new int[] { 5, 3, 10, 6, 8 };
		System.out.println("Max profit 1 Trans. : " + lt.maxProfit1(prices));
		System.out.println("Max profit 2 Trans. : " + lt.maxProfit3(prices));

		nums = new int[] { 4, 5, 6, 7, 0, 1, 2 };
		System.out.println("Minimum: " + lt.findMin(nums));
		System.out.println("Pow: " + (int) lt.pow(8, 5));
		System.out.println("Sqrt: " + (int) lt.mySqrt(139));
		System.out.println("Division: " + (int) lt.divide(15, 3));

		int[] nums1 = new int[] { 5, 7, 8, 11, 14, 15, 17 };
		int[] nums2 = new int[] { 1, 3, 9, 10, 12, 13, 14, 16 };
		System.out
				.println("Median: " + lt.findMedianSortedArrays(nums1, nums2));

		TreeNode newTree = lt.buildTree(nums2);
		System.out.println("Valid BST: " + lt.isValidBST(newTree));
		System.out.println("Num nodes complete BT: " + lt.countNodes(newTree));

		BSTIterator iter = lt.new BSTIterator(newTree);
		System.out.print("BST Iterator: [ ");
		while (iter.hasNext()) {
			System.out.print(iter.next() + " ");
		}
		System.out.println("]");

		preIndex = 0;
		preorder = new int[] { 1, 2, 3, 9, 11, 4, 5, 7, 5 };
		inorder = new int[] { 9, 3, 11, 2, 4, 1, 7, 5, 5 };
		tree = lt.buildTree1(preorder, inorder);

		System.out.println("Max path sum: " + lt.maxPathSum(tree));
		List<Integer> vals = lt.preorderTraversal(newTree);
		System.out.print("Iterative Preorder: [ ");
		for (int i = 0; i < vals.size(); i++) {
			System.out.print(vals.get(i) + " ");
		}
		System.out.println("]");
		newTree = lt.invertTree(newTree);
		System.out.print("Inverse Tree: ");
		lt.printInOrder(newTree);

		System.out.println();
		System.out.println("Min Depth: " + lt.minDepth(newTree));
		System.out.println("Max Depth: " + lt.maxDepth(newTree));
		System.out.println("Balanced BST: " + lt.isBalanced(newTree));

		nums2 = new int[] { 1, 9, 3, 4, 6, 7, 2 };
		newTree = lt.buildTree(nums2);
		lt.recoverTree2(newTree); //
		lt.printInOrder(newTree);

		nums = new int[] { 1, 2, 3, 4, 5, 6, 7, 8 };
		ListNode head = lt.create(nums);
		newTree = lt.sortedListToBST(head);
		lt.printInOrder(newTree);

		System.out.println();
		height = new int[] { 0, 1, 0, 2, 1, 0, 1, 3, 2, 1, 2, 1 };
		System.out.println("Water trapped after rain: " + lt.trap(height) + ":"
				+ lt.trap1(height));

		nums = new int[] { 0, 1, 2, 1, 1, 2, 0 };
		lt.sortColors2(nums);
		System.out.print("[ ");
		for (int i = 0; i < nums.length; i++) {
			System.out.print(nums[i] + " ");
		}
		System.out.println(" ]");

		nums = new int[] { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
		head = lt.create(nums);
		lt.print(lt.reverseKGroup(head, 3));

		nums = new int[] { 1, 2, 3, 4, 5 };
		head = lt.create(nums);
		lt.print(lt.reverseBetween(head, 2, 4));

		ListNode nodeA = lt.create(new int[] { 4, 5, 7, 8, 9 });
		ListNode nodeB = lt.create(new int[] { 6, 2, 5, 7, 8, 9 });
		lt.print(lt.getIntersectionNode(nodeA, nodeB));

		System.out.print("Remove Duplicates From List I: ");
		ListNode dups = lt.create(new int[] { 1, 2, 3, 3, 4, 4, 5 });
		lt.print(lt.deleteDuplicates1(dups));

		System.out.print("Remove Duplicates From List II: ");
		dups = lt.create(new int[] { 1, 2, 3, 3, 4, 4, 5 });
		lt.print(lt.deleteDuplicates(dups));

		System.out.println();
		int mat[][] = { { 1, 2, 3, 4 }, { 5, 6, 7, 8 }, { 9, 10, 11, 12 } };
		List<Integer> llist = lt.spiralOrder(mat);
		System.out.print("Spiral Order: [ ");
		for (int i = 0; i < llist.size(); i++) {
			System.out.print(llist.get(i) + " ");
		}
		System.out.println("]");

		System.out.print("Spiral Matrix: ");
		mat = lt.generateMatrix(3);
		for (int i = 0; i < mat.length; i++) {
			System.out.print("[ ");
			for (int j = 0; j < mat.length; j++) {
				System.out.print(mat[i][j] + " ");
			}
			System.out.print("] ");
		}

		System.out.println();
		nums1 = new int[] { 3, 2, 1, 0, 4 };
		nums2 = new int[] { 2, 3, 1, 1, 4 };
		System.out.println(lt.canJump(nums1) + " - " + lt.canJump(nums2));
		System.out.println("Min Jumps: " + lt.jump(nums2));
		System.out.println("Multiply Strings: " + lt.multiply("12345", "354"));

		String words[] = { "This", "is", "an", "example", "of", "text",
				"justification.", "Good", "luck", "with", "all", "your",
				"endeavors", "big", "or", "small", "crazy", "or", "realistic",
				"whatever", "really", "whatever." };
		List<String> lines = lt.fullJustify(words, 16);
		for (String ln : lines) {
			System.out.println(ln);
		}
		System.out.println("Scrambled: " + lt.isScramble("great", "rgtae"));
		System.out.println("Decoding Ways: " + lt.numDecodings("12"));

		nums = new int[] { 1, 2, 3, 4, 5 };
		head = lt.create(nums);
		lt.reorderList(head);
		System.out.print("Reorder List: ");
		lt.print(head);
		nums = new int[] { 3, 5, 9, 30, 34 };
		System.out.println("Maximum Gap: " + lt.maximumGap(nums));

		nums = new int[] { 5, 7, 2, 3, 4, 9, 6, 2 };
		System.out.println("House Robber: " + lt.rob(nums));

		System.out.print("Subsets: ");
		nums = new int[] { 1, 2, 3, 4 };
		List<List<Integer>> subsets = lt.subsets(nums);
		for (List<Integer> l : subsets) {
			if (l.size() == 0) {
				System.out.print("[ ]");
				continue;
			}
			System.out.print("[ ");
			for (int i : l) {
				System.out.print(i + " ");
			}
			System.out.print("] ");
		}

		System.out.println();
		Set<String> wordDict = new HashSet<String>();
		wordDict.add("leet");
		wordDict.add("code");
		System.out.println("Word Break: " + lt.wordBreak("setcode", wordDict));

		vals = lt.grayCode(2);
		System.out.print("Gray Code: [ ");
		for (int num : vals) {
			System.out.print(num + " ");
		}
		System.out.println("] ");

		lt.zigzagLevelOrder(lt.buildTree(new int[] { 8, 9, 5, 3, 15, 20, 7 }));
		List<TreeNode> trees = lt.generateTrees(4);
		for (TreeNode aTree : trees) {
			lt.printInOrder(aTree);
			System.out.println();
		}

		System.out.println("Min Candies: "
				+ lt.candy1(new int[] { 1, 3, 2, 5, 6, 8, 7 }));
	}

	/**
	 * Given two words word1 and word2, find the minimum number of steps
	 * required to convert word1 to word2. (each operation is counted as 1
	 * step.)
	 * 
	 * <pre>
	 * You have the following 3 operations permitted on a word:
	 * 
	 * a) Insert a character
	 * b) Delete a character
	 * c) Replace a character
	 * </pre>
	 */
	public int minDistance(String word1, String word2) {
		return 0;
	}

	/**
	 * Given a string s1, we may represent it as a binary tree by partitioning
	 * it to two non-empty substrings recursively.
	 * 
	 * Below is one possible representation of s1 = "great":
	 * 
	 * <pre>
	 * 
	 * 	   great               rgeat					rgtae
	 *    /    \ 			  /    \				   /     \
	 *   gr    eat			 rg    eat			      rg     tae
	 *  / \    /  \			/ \    /  \			     / \    /  \
	 * g   r  e   at 	   r   g  e    at	        r   g  ta   e
	 *            / \			  	  /  \			      / \
	 *           a   t				 a    t				 t   a
	 * 
	 * We say that "rgtae" is a scrambled string of "great".
	 * 
	 * Given two strings s1 and s2 of the same length, determine if s2 is a scrambled string of s1.
	 * </pre>
	 */
	public boolean isScramble(String s1, String s2) {
		if (s1 == null || s2 == null)
			return false;
		if (s1.length() != s2.length())
			return false;
		return isScramble(s1, s2, 0, s1.length() - 1);
	}

	private boolean isScramble(String s1, String s2, int i, int j) {
		if (i >= j)
			return false;
		int middle = i + (j - i) / 2;
		String sub = s1.substring(i, middle);
		String lstr = scramble(sub, i, sub.length() - 1);
		if (s2.equals(lstr + s1.substring(middle)))
			return true;
		sub = s1.substring(middle);
		String rstr = scramble(sub, i, sub.length() - 1);
		if (s2.equals(s1.substring(i, middle) + rstr))
			return true;
		if (s2.equals(lstr + rstr))
			return true;
		return isScramble(s1, s2, middle + 1, j)
				|| isScramble(s1, s2, i, middle);
	}

	private String scramble(String str, int i, int j) {
		if (str == null || str.length() == 0)
			return null;
		if (j < i)
			return null;
		else if (j > i) {
			int mid = i + (j - i) / 2;
			String scrambled = scramble(str, mid + 1, j)
					+ scramble(str, i, mid);
			return scrambled;
		}
		return str.charAt(i) + "";
	}

	/**
	 * Reorder list: Given a singly linked list L: L0→L1→…→Ln-1→Ln, reorder it
	 * to: L0→Ln→L1→Ln-1→L2→Ln→2→…
	 * 
	 * You must do this in-place without altering the nodes' values.
	 * 
	 * For example, Given {1,2,3,4}, reorder it to {1,4,2,3}
	 * 
	 * {1,2,3,4,5} => {1,5,2,4} {5,4,3,2,1} {1,5,2,4,3}
	 */
	public void reorderList(ListNode head) {
		if (head == null)
			return;
		int n = length(head);
		reorderList(head, 0, n - 1);
	}

	private void reorderList(ListNode head, int i, int j) {
		if (head == null || j < i)
			return;

		ListNode curr = head;
		ListNode node = head;
		ListNode prev = head;

		int middle = i + (j - i) / 2;

		while (middle > 0) {
			while (node.next != null) {
				prev = node;
				node = node.next;
			}
			node.next = curr.next;
			curr.next = node;
			prev.next = null;
			curr = node.next;
			node = curr;
			middle--;
		}
	}

	/**
	 * Unique Binary Search Tree: Given n, generate all structurally unique
	 * BST's (binary search trees) that store values 1...n.
	 * 
	 * <pre>
	 * For example,
	 * Given n = 3, your program should return all 5 unique BST's shown below.
	 * 
	 *    1         3     3      2      1
	 *     \       /     /      / \      \
	 *      3     2     1      1   3      2
	 *     /     /       \                 \
	 *    2     1         2                 3
	 * </pre>
	 */
	public List<TreeNode> generateTrees(int n) {
		return generateTrees(1, n);
	}

	private List<TreeNode> generateTrees(int start, int end) {
		List<TreeNode> nodes = new ArrayList<TreeNode>();
		if (start > end) {
			nodes.add(null);
			return nodes;
		}
		for (int i = start; i <= end; i++) {
			List<TreeNode> lefts = generateTrees(start, i - 1);
			List<TreeNode> rights = generateTrees(i + 1, end);
			for (TreeNode left : lefts) {
				for (TreeNode right : rights) {
					TreeNode node = new TreeNode(i);
					node.left = left;
					node.right = right;
					nodes.add(node);
				}
			}
		}
		return nodes;
	}

	public int numTrees(int n) {
		int[] trees = new int[n + 1];
		trees[0] = 1;
		for (int i = 1; i <= n; i++) {
			for (int j = 0; j < i; j++) {
				trees[i] += trees[j] * trees[i - j - 1];
			}
		}
		return trees[n];
	}

	public List<TreeNode> generateTrees(int nums[]) {
		com.interviews.facebook.MergeSort sort;
		sort = new com.interviews.facebook.MergeSort(nums);
		sort.sort(com.interviews.facebook.MergeSort.Order.ASC);
		return generateTrees(nums, 1, nums.length - 1);
	}

	public List<TreeNode> generateTrees(int nums[], int start, int end) {
		List<TreeNode> nodes = new ArrayList<TreeNode>();
		if (start > end) {
			nodes.add(null);
			return nodes;
		}
		for (int i = start; i <= end; i++) {
			List<TreeNode> lefts = generateTrees(nums, start, i - 1);
			List<TreeNode> rights = generateTrees(nums, i + 1, end);
			for (TreeNode left : lefts) {
				for (TreeNode right : rights) {
					TreeNode node = new TreeNode(nums[i]);
					node.left = left;
					node.right = right;
					nodes.add(node);
				}
			}
		}
		return nodes;
	}

	/**
	 * Given an unsorted array, find the maximum difference between the
	 * successive elements in its sorted form.
	 * 
	 * Try to solve it in linear time/space.
	 * 
	 * Return 0 if the array contains less than 2 elements.
	 * 
	 * You may assume all elements in the array are non-negative integers and
	 * fit in the 32-bit signed integer range.
	 */
	public int maximumGap(int[] nums) { // TODO
		if (nums == null || nums.length < 2)
			return 0;
		int n = nums.length;
		Bucket[] buckets = new Bucket[n + 1];

		for (int i = 0; i < buckets.length; i++) {
			buckets[i] = new Bucket();
		}
		int max = nums[0], min = nums[0];
		for (int i = 1; i < n; i++) {
			max = Math.max(max, nums[i]);
			min = Math.min(min, nums[i]);
		}
		double interval = (1.0 * n) / (max - min);
		for (int i = 0; i < n; i++) {
			int j = (int) ((nums[i] - min) * interval);
			if (buckets[j].low == -1) {
				buckets[j].low = nums[i];
				buckets[j].high = nums[i];
			} else {
				buckets[j].low = Math.min(buckets[j].low, nums[i]);
				buckets[j].high = Math.max(buckets[j].high, nums[i]);
			}
		}
		int result = 0, prev = buckets[0].high;
		for (int i = 1; i < n; i++) {
			if (buckets[i].low != -1) {
				result = Math.max(result, buckets[i].low - prev);
				prev = buckets[i].high;
			}
		}
		return result;
	}

	int getMax(int num[]) {
		int m = Integer.MIN_VALUE;
		for (int i = 0; i < num.length; i++) {
			m = Math.max(m, num[i]);
		}
		return m;
	}

	void countSort(int num[], int nd) {
		int n = num.length;
		int output[] = new int[n];
		int count[] = new int[10];

		for (int i = 0; i < n; i++) {
			count[(num[i] / nd) % 10]++;
		}

		for (int i = 1; i < 10; i++) {
			count[i] += count[i - 1];
		}

		for (int i = n - 1; i >= 0; i--) {
			output[count[(num[i] / nd) % 10] - 1] = num[i];
			count[(num[i] / nd) % 10]--;
		}

		for (int i = 0; i < n; i++) {
			num[i] = output[i];
		}
	}

	void radixsort(int num[]) {
		int max_n = getMax(num);
		for (int nd = 1; max_n / nd > 0; nd *= 10) {
			countSort(num, nd);
		}
	}

	public int maximumGap2(int num[]) {
		if (num.length < 2)
			return 0;
		radixsort(num);
		int res = Math.abs(num[1] - num[0]);
		for (int i = 2; i < num.length; i++) {
			if (num[i] - num[i - 1] > res) {
				res = Math.abs(num[i] - num[i - 1]);
			}
		}
		return res;
	}

	/**
	 * House Robber: You are a professional robber planning to rob houses along
	 * a street. Each house has a certain amount of money stashed, the only
	 * constraint stopping you from robbing each of them is that adjacent houses
	 * have security system connected and it will automatically contact the
	 * police if two adjacent houses were broken into on the same night.
	 * 
	 * Given a list of non-negative integers representing the amount of money of
	 * each house, determine the maximum amount of money you can rob tonight
	 * without alerting the police.
	 */
	public int rob(int[] nums) {
		if (nums == null || nums.length == 0)
			return 0;
		int d[] = new int[nums.length];
		for (int i = 0; i < nums.length; i++) {
			int pos = getMaxPosition(nums, d, i);
			if (!safe(d, pos))
				continue;
			d[pos] = nums[pos];
		}
		int amount = 0;
		for (int i = 0; i < d.length; i++) {
			amount = amount + d[i];
		}
		return amount;
	}

	public int rob1(int[] num) { // Dynamic Programming version
		if (num == null || num.length == 0)
			return 0;
		int n = num.length;
		int[] dp = new int[n + 1];
		dp[0] = 0;
		dp[1] = num[0];
		for (int i = 2; i < n + 1; i++) {
			dp[i] = Math.max(dp[i - 1], dp[i - 2] + num[i - 1]);
		}
		return dp[n];
	}

	private boolean safe(int d[], int i) {
		if (i < 0 || i > d.length || d[i] != 0)
			return false;
		if (i == 0)
			return d[i + 1] == 0;
		if (i == d.length - 1)
			return d[i - 1] == 0;
		return d[i - 1] == 0 && d[i + 1] == 0;
	}

	private int getMaxPosition(int[] nums, int[] d, int i) {
		if (nums == null || nums.length == 0 || i >= nums.length)
			return 0;
		int j = 0, k = 0;
		if (i == 0)
			j = 1;
		else if (i == nums.length - 1) {
			k = i - 1;
			j = i;
		} else {
			k = i - 1;
			j = i + 1;
		}
		int pos = -1;
		if (safe(d, k)) {
			pos = Math.max(nums[i], nums[k]) == nums[i] ? i : k;
		}
		if (safe(d, j)) {
			int tmp = Math.max(nums[i], nums[j]) == nums[i] ? i : j;
			pos = Math.max(pos, tmp);
		}
		return pos;
	}

	/**
	 * After robbing those houses on that street, the thief has found himself a
	 * new place for his thievery so that he will not get too much attention.
	 * This time, all houses at this place are arranged in a circle. That means
	 * the first house is the neighbor of the last one. Meanwhile, the security
	 * system for these houses remain the same as for those in the previous
	 * street.
	 * 
	 * Given a list of non-negative integers representing the amount of money of
	 * each house, determine the maximum amount of money you can rob tonight
	 * without alerting the police.
	 */
	public int rob2(int[] nums) { // TODO
		if (nums == null || nums.length == 0)
			return 0;

		int n = nums.length;

		if (n == 1) {
			return nums[0];
		}
		if (n == 2) {
			return Math.max(nums[1], nums[0]);
		}

		// include 1st element, and not last element
		int[] dp = new int[n + 1];
		dp[0] = 0;
		dp[1] = nums[0];

		for (int i = 2; i < n; i++) {
			dp[i] = Math.max(dp[i - 1], dp[i - 2] + nums[i - 1]);
		}

		// not include first element, and include last element
		int[] dr = new int[n + 1];
		dr[0] = 0;
		dr[1] = nums[1];

		for (int i = 2; i < n; i++) {
			dr[i] = Math.max(dr[i - 1], dr[i - 2] + nums[i]);
		}

		return Math.max(dp[n - 1], dr[n - 1]);
	}

	/**
	 * Binary Tree Zigzag Level Order Traversal: Given a binary tree, return the
	 * zigzag level order traversal of its nodes' values. (i.e, from left to
	 * right, then right to left for the next level and alternate between).
	 * 
	 * <pre>
	 * For example:
	 * Given binary tree {3,9,20,#,#,15,7},
	 * 
	 *     3
	 *    / \
	 *   9  20
	 *     /  \
	 *    15   7
	 * 
	 * return its zigzag level order traversal as:
	 * 
	 * [
	 *   [3],
	 *   [20,9],
	 *   [15,7]
	 * ]
	 * </pre>
	 */
	public List<List<Integer>> zigzagLevelOrder(TreeNode root) {

		Stack<TreeNode> s1 = new Stack<TreeNode>();
		Stack<TreeNode> s2 = new Stack<TreeNode>();
		s1.push(root);
		List<List<Integer>> zigzag = new ArrayList<List<Integer>>();
		while (!(s1.isEmpty() && s2.isEmpty())) {
			List<Integer> level = new ArrayList<Integer>();
			if (!s1.isEmpty()) {
				while (!s1.isEmpty()) {
					TreeNode node = s1.pop();
					level.add(node.val);
					if (node.left != null) {
						s2.push(node.left);
					}
					if (node.right != null) {
						s2.push(node.right);
					}
				}
			} else {
				while (!s2.isEmpty()) {
					TreeNode node = s2.pop();
					level.add(node.val);
					if (node.right != null) {
						s1.push(node.right);
					}
					if (node.left != null) {
						s1.push(node.left);
					}
				}
			}
			zigzag.add(level);
		}
		return zigzag;
	}

	/**
	 * Given a string s and a dictionary of words dict, determine if s can be
	 * segmented into a space-separated sequence of one or more dictionary
	 * words.
	 * 
	 * For example, given s = "leetcode", dict = ["leet", "code"].
	 * 
	 * Return true because "leetcode" can be segmented as "leet code".
	 */
	public boolean wordBreak(String s, Set<String> wordDict) {
		return wordBreak(s, wordDict, 0);
	}

	private boolean wordBreak(String s, Set<String> wordDict, int start) {
		if (wordDict == null || wordDict.size() == 0)
			return false;
		if (start == s.length())
			return true;
		for (String str : wordDict) {
			int len = str.length();
			int end = start + len;
			if (end > s.length())
				continue;
			if (s.substring(start, end).equals(str)) {
				if (wordBreak(s, wordDict, end))
					return true;
			}
		}
		return false;
	}

	/**
	 * Gray code: The gray code is a binary numeral system where two successive
	 * values differ in only one bit.
	 * 
	 * Given a non-negative integer n representing the total number of bits in
	 * the code, print the sequence of gray code. A gray code sequence must
	 * begin with 0.
	 * 
	 * For example, given n = 2, return [0,1,3,2]. Its gray code sequence is:
	 * 
	 * <pre>
	 * 00 - 0
	 * 01 - 1
	 * 11 - 3
	 * 10 - 2
	 * </pre>
	 * 
	 * Note: For a given n, a gray code sequence is not uniquely defined.
	 * 
	 * For example, [0,2,3,1] is also a valid gray code sequence according to
	 * the above definition.
	 * 
	 * For now, the judge is able to judge based on one instance of gray code
	 * sequence. Sorry about that.
	 * 
	 * We need to generate 2^n numbers, so the complexity is O(2^n).
	 */
	public List<Integer> grayCode(int n) { // TODO
		List<Integer> ret = new ArrayList<Integer>();
		ret.add(0);
		for (int i = 0; i < n; i++) {
			int size = ret.size();
			for (int j = size - 1; j >= 0; j--)
				ret.add(ret.get(j) + size);
		}
		return ret;
	}

	/**
	 * Subsets: Given a set of distinct integers, nums, return all possible
	 * subsets.
	 * 
	 * Note:
	 * 
	 * Elements in a subset must be in non-descending order. The solution set
	 * must not contain duplicate subsets.
	 * 
	 * For example, If nums = [1,2,3], a solution is:
	 * 
	 * <pre>
	 * [
	 *   [3],
	 *   [1],
	 *   [2],
	 *   [1,2,3],
	 *   [1,3],
	 *   [2,3],
	 *   [1,2],
	 *   []
	 * ]
	 * </pre>
	 */
	public List<List<Integer>> subsets(int[] nums) {
		if (nums == null)
			return null;
		if (nums.length == 0)
			return new ArrayList<List<Integer>>();
		Map<Integer, Boolean> visited;
		visited = new HashMap<Integer, Boolean>();
		List<List<Integer>> subsets = new ArrayList<List<Integer>>();
		for (int i = 0; i < nums.length; i++) {
			if (visited.containsKey(nums[i]))
				continue;
			List<Integer> list = new ArrayList<Integer>();
			subsets(nums, visited, subsets, list, i);
		}
		subsets.add(new ArrayList<Integer>());
		return subsets;
	}

	private void subsets(int[] nums, Map<Integer, Boolean> visited,
			List<List<Integer>> subsets, List<Integer> list, int i) {
		visited.put(nums[i], true);
		list.add(nums[i]);
		if (!subsets.contains(list))
			subsets.add(list);
		List<Integer> adjs = adjs(nums, visited, i);
		int ind = subsets.size() - 1;
		for (int j : adjs) {
			for (int k = 0; k <= ind; k++) {
				list = subsets.get(k);
				list = new ArrayList<Integer>(list);
				list.add(nums[j]);
				if (!subsets.contains(list)) {
					subsets.add(list);
				}
			}
			list = new ArrayList<Integer>();
			subsets(nums, visited, subsets, list, j);
		}
	}

	public List<Integer> adjs(int[] nums, Map<Integer, Boolean> visited, int i) {
		if (nums == null || i >= nums.length)
			return null;
		List<Integer> lst = new ArrayList<Integer>();
		for (int j = i + 1; j < nums.length; j++)
			if (!visited.containsKey(nums[j]))
				lst.add(j);
		return lst;
	}

	/**
	 * Subsets 2: Given a collection of integers that might contain duplicates,
	 * nums, return all possible subsets.
	 * 
	 * Note:
	 * 
	 * Elements in a subset must be in non-descending order. The solution set
	 * must not contain duplicate subsets.
	 * 
	 * For example, If nums = [1,2,2], a solution is:
	 * 
	 * <pre>
	 * [
	 *   [2],
	 *   [1],
	 *   [1,2,2],
	 *   [2,2],
	 *   [1,2],
	 *   []
	 * ]
	 * </pre>
	 * 
	 * TODO: Modify {@link LeetCode#adjs(int[], Map, int)} and store the indices
	 * instead of the values of the array.
	 */
	public List<List<Integer>> subsetsWithDup(int[] nums) {
		return null;
	}

	/**
	 * There are N children standing in a line. Each child is assigned a rating
	 * value.
	 * 
	 * You are giving candies to these children subjected to the following
	 * requirements:
	 * 
	 * Each child must have at least one candy. Children with a higher rating
	 * get more candies than their neighbors.
	 * 
	 * What is the minimum candies you must give?
	 */
	public int candy(int[] ratings) { // FIXME (See correct solution below)
		if (ratings == null || ratings.length == 0)
			return 0;
		int max = ratings[0], n = ratings.length;
		int d[] = new int[n + 1]; // candies
		for (int i = 0; i < n; i++) {
			d[i] = 1;
			if (ratings[i] > max)
				max = ratings[i];
		}
		for (int i = 1; i < n - 1; i++) {
			if (ratings[i] > ratings[i - 1]) {
				d[i] = d[i - 1] + 1;
			}
			if (ratings[i] > ratings[i + 1]) {
				if (d[i] < d[i + 1])
					d[i] = d[i + 1] + 1;
			}
		}
		int nums[] = new int[max + 1];
		for (int i = 0; i < max; i++) {
			nums[i] = -1;
		}
		for (int i = 0; i < n; i++) {
			nums[ratings[i]] = i;
		}
		int j = 0;
		for (int i = 1; i < max; i++) {
			if (nums[i] < 0 || nums[j] < 0) {
				if (nums[i] >= 0)
					j = i;
				continue;
			}
			if (d[nums[i]] < d[nums[j]])
				d[nums[i]] = d[nums[j]];
			j = i;
		}
		int candies = 0;
		for (int i = 0; i < n; i++) {
			candies += d[i];
		}
		return candies;
	}

	public int candy1(int[] ratings) {
		if (ratings == null || ratings.length == 0) {
			return 0;
		}

		int[] candies = new int[ratings.length];
		candies[0] = 1;

		// from let to right
		for (int i = 1; i < ratings.length; i++) {
			if (ratings[i] > ratings[i - 1]) {
				candies[i] = candies[i - 1] + 1;
			} else {
				// if not ascending, assign 1
				candies[i] = 1;
			}
		}

		int result = candies[ratings.length - 1];

		// from right to left
		for (int i = ratings.length - 2; i >= 0; i--) {
			int cur = 1;
			if (ratings[i] > ratings[i + 1]) {
				cur = candies[i + 1] + 1;
			}
			result += Math.max(cur, candies[i]);
			candies[i] = cur;
		}
		return result;
	}

	/**
	 * Distinct subsequences: Given a string S and a string T, count the number
	 * of distinct subsequences of T in S.
	 * 
	 * A subsequence of a string is a new string which is formed from the
	 * original string by deleting some (can be none) of the characters without
	 * disturbing the relative positions of the remaining characters. (ie, "ACE"
	 * is a subsequence of "ABCDE" while "AEC" is not).
	 * 
	 * Here is an example: S = "rabbbit", T = "rabbit"
	 * 
	 * Return 3.
	 */
	public int numDistinct(String s, String t) {
		return 0;
	}

	/**
	 * Decode ways: A message containing letters from A-Z is being encoded to
	 * numbers using the following mapping:
	 * 
	 * <pre>
	 * 'A' -> 1
	 * 'B' -> 2
	 * ...
	 * 'Z' -> 26
	 * </pre>
	 * 
	 * Given an encoded message containing digits, determine the total number of
	 * ways to decode it.
	 * 
	 * For example, Given encoded message "12", it could be decoded as "AB" (1
	 * 2) or "L" (12).
	 * 
	 * The number of ways decoding "12" is 2.
	 * 
	 * <pre>
	 *   123 => ABC  => AW  => LC
	 *  1234 => ABCD => AWD => LCD
	 *  2145 => BADE => BNE => UDE
	 * 21124 => BAABD => BAAX => BALD => BKBD => BKX => UABD => ULD => UAX
	 * </pre>
	 */
	public int numDecodings(String s) {
		// TODO: For this code to work on numbers like "012",
		// remove the leading 0s first.
		if (s == null || s.length() == 0)
			throw new IllegalArgumentException();
		int upperLimit = 26;
		int maxHeadSize = ("" + upperLimit).length();
		return numDecodings(s, upperLimit, maxHeadSize);
	}

	private int numDecodings(String str, int upperLimit, int max) {
		if (str.length() == 0)
			return 1;
		int sum = 0;
		for (int i = 1; i <= max && i <= str.length(); i++) {
			String head = str.substring(0, i);
			String tail = str.substring(i);
			if (Integer.parseInt(head) > upperLimit)
				break;
			sum += numDecodings(tail, upperLimit, max);
		}
		return sum;
	}

	/**
	 * Multiply Strings. Given two numbers represented as strings, return
	 * multiplication of the numbers as a string.
	 * 
	 * Note: The numbers can be arbitrarily large and are non-negative.
	 */
	public String multiply(String num1, String num2) {
		if (num1 == null || num1 == "")
			return num2;
		if (num2 == null || num2 == "")
			return num1;
		String sum = "0";
		int n = num2.length() - 1;
		for (int i = n; i >= 0; i--) {
			char c = num2.charAt(i);
			String mul = multiply(num1, c);
			sum = add(mul, sum);
			num1 += 0;
		}
		return sum;
	}

	private String multiply(String num, char dig) {
		if (num == null)
			return null;
		int n = num.length() - 1, rem = 0;
		String mult = "";
		int oprnd1 = dig - '0';
		for (int i = n; i >= 0; i--) {
			int oprnd2 = num.charAt(i) - '0';
			int sum = oprnd1 * oprnd2 + rem;
			rem = sum / 10;
			if (sum >= 10)
				sum = sum - (rem * 10);
			mult = sum + mult;
		}
		return mult;
	}

	private String add(String num1, String num2) {
		return add(num1, num2, num1.length() - 1, num2.length() - 1, 0);
	}

	private String add(String num1, String num2, int i, int j, int rem) {
		int lval = 0, rval = 0;
		if (i < 0 && j < 0)
			return "";
		if (i >= 0)
			lval = num1.charAt(i) - '0';
		if (j >= 0)
			rval = num2.charAt(j) - '0';
		int sum = lval + rval + rem;
		rem = sum >= 10 ? 1 : 0;
		sum = rem > 0 ? sum - 10 : sum;
		return add(num1, num2, i - 1, j - 1, rem) + sum;
	}

	/**
	 * Text Justification. Given an array of words and a length L, format the
	 * text such that each line has exactly L characters and is fully (left and
	 * right) justified.
	 * 
	 * You should pack your words in a greedy approach; that is, pack as many
	 * words as you can in each line. Pad extra spaces ' ' when necessary so
	 * that each line has exactly L characters.
	 * 
	 * Extra spaces between words should be distributed as evenly as possible.
	 * If the number of spaces on a line do not divide evenly between words, the
	 * empty slots on the left will be assigned more spaces than the slots on
	 * the right.
	 * 
	 * For the last line of text, it should be left justified and no extra space
	 * is inserted between words.
	 * 
	 * <pre>
	 * For example,
	 * words: ["This", "is", "an", "example", "of", "text", "justification."]
	 * L: 16.
	 * 
	 * Return the formatted lines as:
	 * 
	 * [
	 *    "This    is    an",
	 *    "example  of text",
	 *    "justification.  "
	 * ]
	 * </pre>
	 */
	public List<String> fullJustify(String[] words, int maxWidth) {
		if (words == null || maxWidth <= 0)
			return null;
		int nLength, j = 0, i, length = words[0].length();
		List<String> lines = new ArrayList<String>();
		for (i = 1; i < words.length; i++) {
			nLength = length + words[i].length() + (i - j);
			if (nLength >= maxWidth) {
				int diff = maxWidth - length;
				if (nLength == maxWidth) {
					diff -= words[i].length();
					lines.add(newLine(words, j, i, diff));
					j = ++i;
				} else {
					lines.add(newLine(words, j, i - 1, diff));
					j = i;
				}
				length = 0;
			}
			length += words[i].length();
		}
		if (length < maxWidth) {
			int diff = maxWidth - length;
			lines.add(newLine(words, j, i - 1, diff));
		}
		return lines;
	}

	private String newLine(String[] words, int j, int i, int diff) {
		if (words == null || i < j || diff < 0)
			return null;
		while (diff > 0) {
			if (i == j) {
				words[i] += " ";
				diff--;
			} else {
				for (int k = j; k < i && diff > 0; k++) {
					words[k] += " ";
					diff--;
				}
			}
		}
		StringBuilder line = new StringBuilder();
		for (int k = j; k <= i; k++) {
			line.append(words[k]);
		}
		return line.toString();
	}

	/**
	 * Given an array of non-negative integers, you are initially positioned at
	 * the first index of the array.
	 * 
	 * Each element in the array represents your maximum jump length at that
	 * position.
	 * 
	 * Determine if you are able to reach the last index.
	 * 
	 * <pre>
	 * For example:
	 * A = [2,3,1,1,4], return true.
	 * 
	 * A = [3,2,1,0,4], return false.
	 * </pre>
	 */
	public boolean canJump(int[] nums) {
		if (nums == null)
			return false;
		int j = nums.length - 1;
		Map<Integer, Boolean> visited;
		visited = new HashMap<Integer, Boolean>();
		for (int i = 0; i < j; i++) {
			if (visited.containsKey(i))
				continue;
			if (canReach(nums, visited, i, j))
				return true;
		}
		return false;
	}

	private boolean canReach(int[] nums, Map<Integer, Boolean> visited, int i,
			int j) {
		if (i == j)
			return true;
		visited.put(i, true);
		List<Integer> adj = adj(nums, visited, i);
		if (adj.isEmpty())
			return false;
		for (int k : adj) {
			if (canReach(nums, visited, k, j))
				return true;
		}
		return false;
	}

	public List<Integer> adj(int[] nums, Map<Integer, Boolean> visited, int i) {
		if (nums == null || i >= nums.length)
			return null;
		List<Integer> lst = new ArrayList<Integer>();
		for (int j = 1; j <= nums[i]; j++) {
			if (!visited.containsKey(i + j))
				lst.add(i + j);
		}
		return lst;
	}

	/**
	 * Given an array of non-negative integers, you are initially positioned at
	 * the first index of the array.
	 * 
	 * Each element in the array represents your maximum jump length at that
	 * position.
	 * 
	 * Your goal is to reach the last index in the minimum number of jumps.
	 * 
	 * For example: Given array A = [2,3,1,1,4]
	 * 
	 * The minimum number of jumps to reach the last index is 2. (Jump 1 step
	 * from index 0 to 1, then 3 steps to the last index.)
	 */
	public int jump(int[] nums) {
		if (nums == null)
			throw new IllegalArgumentException();
		int j = nums.length - 1, min = Integer.MAX_VALUE;
		Map<Integer, Boolean> visited;
		visited = new HashMap<Integer, Boolean>();
		List<Integer> jumps = new ArrayList<Integer>();
		for (int i = 0; i < j; i++) {
			if (visited.containsKey(i))
				continue;
			minJump(nums, visited, i, j, i, jumps);
		}
		if (jumps.isEmpty())
			return -1;
		for (int i : jumps) {
			if (i < min)
				min = i;
		}
		return min;
	}

	public void minJump(int[] nums, Map<Integer, Boolean> visited, int i,
			int j, int step, List<Integer> jumps) {
		if (i == j) {
			jumps.add(step);
			step--;
			return;
		}
		step++;
		visited.put(i, true);
		List<Integer> adj = adj(nums, visited, i);
		if (adj.isEmpty()) {
			return;
		}
		for (int k : adj) {
			minJump(nums, visited, k, j, step, jumps);
		}
	}

	/**
	 * Given a matrix of m x n elements (m rows, n columns), return all elements
	 * of the matrix in spiral order.
	 * 
	 * <pre>
	 * For example,
	 * Given the following matrix:
	 * 
	 * [
	 *  [ 1, 2, 3 ],
	 *  [ 4, 5, 6 ],
	 *  [ 7, 8, 9 ]
	 * ]
	 * 
	 * You should return [1, 2, 3, 6, 9, 8, 7, 4, 5]
	 * </pre>
	 */
	public List<Integer> spiralOrder(int[][] matrix) {
		if (matrix == null)
			return null;
		List<Integer> vals = new ArrayList<Integer>();
		int i = 0, j = matrix.length - 1;
		int k = matrix[i].length - 1, m, l;
		while (i <= j) {
			for (m = i; m < k; m++) {
				vals.add(matrix[i][m]);
			}
			for (l = i; l <= j; l++) {
				vals.add(matrix[l][m]);
			}
			for (l = m - 1; l > i; l--) {
				vals.add(matrix[j][l]);
			}
			for (m = j; m > i; m--) {
				vals.add(matrix[m][l]);
			}
			j--;
			i++;
			k--;
		}
		return vals;
	}

	/**
	 * Given an integer n, generate a square matrix filled with elements from 1
	 * to n2 in spiral order.
	 * 
	 * <pre>
	 * For example,
	 * Given n = 3,
	 * You should return the following matrix:
	 * 
	 * [
	 *  [ 1, 2, 3 ],
	 *  [ 8, 9, 4 ],
	 *  [ 7, 6, 5 ]
	 * ]
	 * </pre>
	 */
	public int[][] generateMatrix(int n) {
		if (n <= 0)
			return null;
		int nums[][] = new int[n][n];
		for (int i = 0; i < n * n; i++) {
			nums[i / n][i % n] = i + 1;
		}
		List<Integer> list = spiralOrder(nums);
		for (int i = 0; i < list.size(); i++) {
			nums[i / n][i % n] = list.get(i);
		}
		return nums;
	}

	private ListNode create(int nums[]) {
		ListNode head = null, tail = null;
		for (int i = 0; i < nums.length; i++) {
			ListNode node = new ListNode(nums[i]);
			if (head == null) {
				head = node;
				tail = head;
			} else {
				tail.next = node;
				tail = node;
			}
		}
		return head;
	}

	private void print(ListNode node) {
		if (node == null)
			return;
		System.out.print("[ ");
		while (node != null) {
			System.out.print(node.val + " ");
			node = node.next;
		}
		System.out.println("]");
	}

	/**
	 * 4Sum : Given an array S of n integers, are there elements a, b, c, and d
	 * in S such that a + b + c + d = target? Find all unique quadruplets in the
	 * array which gives the sum of target.
	 * 
	 * Note:
	 * 
	 * Elements in a quadruplet (a,b,c,d) must be in non-descending order. (i.e,
	 * a ≤ b ≤ c ≤ d) The solution set must not contain duplicate quadruplets.
	 * 
	 * <pre>
	 *     For example, given array S = {1 0 -1 0 -2 2}, and target = 0.
	 * 
	 *     A solution set is:
	 *     (-1,  0, 0, 1)
	 *     (-2, -1, 1, 2)
	 *     (-2,  0, 0, 2)
	 * </pre>
	 */
	public List<List<Integer>> fourSum(int[] nums, int target) {
		return null;
	}

	public List<List<Integer>> threeSum(int[] nums) {
		return null;
	}

	/**
	 * 3Sum Closest: Given an array S of n integers, find three integers in S
	 * such that the sum is closest to a given number, target. Return the sum of
	 * the three integers. You may assume that each input would have exactly one
	 * solution.
	 * 
	 * <pre>
	 *     For example, given array S = {-1 2 1 -4}, and target = 1.
	 *     The sum that is closest to the target is 2. (-1 + 2 + 1 = 2).
	 * </pre>
	 */
	public int threeSumClosest(int[] nums, int target) {
		return 0;
	}

	public TreeNode convert(Node<Integer> node) {
		TreeNode root = null;
		return convert(node, root);
	}

	private TreeNode convert(Node<Integer> nnode, TreeNode tnode) {
		if (nnode == null) {
			return null;
		}
		tnode = new TreeNode(nnode.val);
		tnode.left = convert(nnode.left, tnode);
		tnode.right = convert(nnode.right, tnode);
		return tnode;
	}

	/**
	 * Swap nodes in pairs: Given a linked list, swap every two adjacent nodes
	 * and return its head.
	 * 
	 * For example, Given 1->2->3->4, you should return the list as 2->1->4->3.
	 * 
	 * Your algorithm should use only constant space. You may not modify the
	 * values in the list, only nodes itself can be changed.
	 */
	public ListNode swapPairs(ListNode head) {
		if (head == null || head.next == null)
			return head;
		ListNode curr = head, succ = head.next;
		while (succ != null) {
			ListNode tmp = succ.next;
			curr.next = tmp;
			succ.next = curr;
			curr = tmp;
			succ = tmp == null ? tmp : tmp.next;
		}
		return head;
	}

	/**
	 * Sort a linked list in O(n log n) time using constant space complexity.
	 */
	public ListNode sortList(ListNode head) {
		if (head == null)
			return null;
		int length = length(head);
		sortList(head, 0, length);
		return head;
	}

	private int length(ListNode head) {
		if (head == null)
			return 0;
		int length = 0;
		ListNode node = head;
		while (node != null) {
			length += 1;
			node = node.next;
		}
		return length;
	}

	private ListNode pointerAt(ListNode head, int k) {
		if (head == null)
			return head;
		ListNode node = head;
		int i = 1;
		while (i < k) {
			i++;
			node = node.next;
		}
		return node;
	}

	private void sortList(ListNode head, int lb, int ub) {
		if (head == null || lb > ub)
			return;
		int mid = lb + (ub - lb) / 2;
		sortList(head, 0, mid);
		sortList(head, mid + 1, ub);
		merge(head, lb, mid);
	}

	private void merge(ListNode node, int lb, int mid) {
		if (node == null)
			return;
		ListNode midd = pointerAt(node, mid);
		ListNode left = node;
		ListNode curr = node;
		ListNode right = midd.next;
		while (left != midd && right != null) {
			if (left.val > right.val) {
				curr.val = left.val;
				left = left.next;
			} else {
				curr.val = right.val;
				right = right.next;
			}
			curr = curr.next;
		}
		while (left != midd) {
			curr.val = left.val;
			left = left.next;
			curr = curr.next;
		}
	}

	/**
	 * Given a linked list, reverse the nodes of a linked list k at a time and
	 * return its modified list. If the number of nodes is not a multiple of k
	 * then left-out nodes in the end should remain as it is. You may not alter
	 * the values in the nodes, only nodes itself may be changed. Only constant
	 * memory is allowed.
	 * 
	 * <pre>
	 * For example,
	 * Given this linked list: 1->2->3->4->5
	 * For k = 2, you should return: 2->1->4->3->5
	 * For k = 3, you should return: 3->2->1->4->5
	 * </pre>
	 */
	public ListNode reverseKGroup(ListNode head, int k) {
		int length = length(head);
		if (length == 0 || k > length)
			return head;
		return reverseKGroup(head, length, k);
	}

	private ListNode reverseKGroup(ListNode head, int n, int k) {
		if (head == null)
			return null;
		if (n < k)
			return head;
		ListNode node = pointerAt(head, k);
		if (node != null) {
			ListNode next = node.next;
			node.next = reverseKGroup(next, n - k, k);
			head = flip(head, node);
		}
		return head;
	}

	private ListNode flip(ListNode head, ListNode tail) {
		if (head == null)
			return null;
		tail = tail == null ? null : tail.next;
		ListNode node = head, flip = tail;
		while (node != tail) {
			ListNode tmp = new ListNode(node.val);
			tmp.next = flip;
			flip = tmp;
			node = node.next;
		}
		return flip;
	}

	/**
	 * Reverse a linked list from position m to n. Do it in-place and in
	 * one-pass.
	 * 
	 * <pre>
	 * For example:
	 * Given 1->2->3->4->5->NULL, m = 2 and n = 4,
	 * 
	 * return 1->4->3->2->5->NULL.
	 * Given m, n satisfy the following condition: 1 ≤ m ≤ n ≤ length
	 * of list.
	 * 1->2->3->4->5->6->7->8->9, m = 3, and n = 6
	 * prevS = 2, prevE = 5, start = 3, end = 6
	 * start.next = 7, prevE.next = 3, prevS.next = 6, end.next = 5
	 * 1->2->6->5->3->7->8->9
	 * 
	 * </pre>
	 */
	public ListNode reverseBetween(ListNode head, int m, int n) {
		if (head == null)
			return null;
		if (m > n)
			return head;
		int i = 1;
		ListNode node = head;

		ListNode start = head;
		ListNode end = head;

		ListNode prevS = head;
		ListNode prevE = head;

		while (true) {
			if (i == m) {
				start = node;
			} else if (i < m) {
				prevS = node;
			}
			if (i == n) {
				end = node;
				ListNode tmp = start.next;
				start.next = end.next;
				prevE.next = start;
				prevS.next = end;
				end.next = tmp;
				return head;
			} else if (i < n) {
				prevE = node;
			}
			i++;
			node = node.next;
		}
	}

	/**
	 * Copy List with Random Pointer. A linked list is given such that each node
	 * contains an additional random pointer which could point to any node in
	 * the list or null. Return a deep copy of the list.
	 */
	public RandomListNode copyRandomList(RandomListNode head) {
		if (head == null)
			return null;
		RandomListNode copy = null, tail = null, node = head;
		Map<RandomListNode, RandomListNode> rands;
		rands = new HashMap<RandomListNode, RandomListNode>();
		while (node != null) {
			RandomListNode tmp = new RandomListNode(node.label);
			if (copy == null) {
				copy = tmp;
				tail = copy;
			} else {
				tail.next = tmp;
				tail = tmp;
			}
			if (node.random != null) {
				tmp = new RandomListNode(node.random.label);
				rands.put(node, tmp);
			}
			node = node.next;
		}

		node = copy;
		RandomListNode curr = head;
		Map<RandomListNode, RandomListNode> next;
		next = new HashMap<RandomListNode, RandomListNode>();
		while (curr != null) {
			node.random = rands.get(curr);
			next.put(node, node.next);
			node = node.next;
			curr = curr.next;
		}

		node = copy;
		while (node != null) {
			if (node.random != null)
				node.random.next = next.get(node.random);
			node = node.next;
		}
		return copy;
	}

	/**
	 * Write a program to find the node at which the intersection of two singly
	 * linked lists begins. Notes:
	 * 
	 * If the two linked lists have no intersection at all, return null. The
	 * linked lists must retain their original structure after the function
	 * returns. You may assume there are no cycles anywhere in the entire linked
	 * structure. Your code should preferably run in O(n) time and use only O(1)
	 * memory.
	 * 
	 * <pre>
	 * For example, the following two linked lists:
	 * 
	 * A:          a1 → a2
	 *                    ↘
	 *                      c1 → c2 → c3
	 *                    ↗            
	 * B:     b1 → b2 → b3
	 * 
	 * </pre>
	 */
	public ListNode getIntersectionNode(ListNode headA, ListNode headB) { // FIXME
		if (headA == null || headB == null)
			return null;
		int a = length(headA);
		int b = length(headB);
		if (a > b) {
			return getIntersection(headB, headA);
		} else if (b > a) {
			return getIntersection(headA, headB);
		} else {
			ListNode nodeA = headA;
			ListNode nodeB = headB;
			while (nodeA != null) {
				if (nodeA == nodeB)
					return nodeA;
				nodeA = nodeA.next;
				nodeB = nodeB.next;
			}
			return null;
		}
	}

	private ListNode getIntersection(ListNode headA, ListNode headB) {
		if (isSubList(headA, headB))
			return headA;
		ListNode nodeA = headA;
		ListNode nodeB = headB;
		ListNode prevB = headB;
		ListNode prevA = headA;
		while (nodeA != null && nodeB != null) {
			if (prevB == nodeA)
				return prevB;
			if (prevA == nodeB)
				return prevA;
			if (nodeA == nodeB)
				return nodeA;
			prevA = nodeA;
			prevB = nodeB.next;
			nodeB = prevB == null ? null : prevB.next;
			nodeA = nodeA.next;
		}
		return null;
	}

	private boolean isSubList(ListNode headA, ListNode headB) {
		if (headA == null || headB == null)
			return false;
		ListNode nodeA = headA;
		ListNode nodeB = headB;
		while (nodeB != nodeA && nodeB != null) {
			nodeB = nodeB.next;
		}
		return nodeB == nodeA;
	}

	/**
	 * Given a linked list, return the node where the cycle begins. If there is
	 * no cycle, return null.
	 * 
	 * Follow up: Can you solve it without using extra space?
	 * 
	 * <pre>
	 * A->B->C->D->E->F->G->D
	 * 
	 * P1: A B C D E
	 *  
	 * P2: C E G C E
	 * </pre>
	 */
	public ListNode detectCycle(ListNode head) { // FIXME
		if (head == null || head.next == null)
			return null;
		ListNode next = head;
		ListNode nnext = head.next.next;
		ListNode pnnext = nnext;
		while (next != null && nnext != null) {
			if (next == nnext)
				return pnnext.next;
			pnnext = nnext;
			next = next.next;
			if (nnext.next == null)
				return null;
			nnext = nnext.next.next;
		}
		return null;
	}

	/**
	 * Given a sorted linked list, delete all duplicates such that each element
	 * appear only once.
	 * 
	 * <pre>
	 * For example, Given 1->1->2, return 1->2. 
	 * Given 1->1->2->3->3, return 1->2->3
	 * </pre>
	 */
	public ListNode deleteDuplicates1(ListNode head) {
		if (head == null)
			return null;
		ListNode node = null;
		ListNode curr = head;
		ListNode next = curr.next;
		while (next != null) {
			if (next.val != curr.val) {
				if (node == null) {
					head = curr;
					node = head;
				} else {
					node.next = curr;
					node = curr;
				}
			} else {
				if (next.next == null) {
					node.next = next;
					node = next;
				}
			}
			curr = next;
			next = next.next;
		}
		return head;
	}

	/**
	 * Remove Duplicates from Sorted List II: Given a sorted linked list, delete
	 * all nodes that have duplicate numbers, leaving only distinct numbers from
	 * the original list.
	 * 
	 * For example, Given 1->2->3->3->4->4->5, return 1->2->5. Given
	 * 1->1->1->2->3, return 2->3.
	 */
	public ListNode deleteDuplicates(ListNode head) {
		if (head == null)
			return null;

		ListNode node = null;
		ListNode pred = head;
		ListNode curr = head;
		ListNode next = curr.next;

		while (next != null) {
			if (curr.val != next.val) {
				if (node == null) {
					head = curr;
					node = head;
				} else {
					node.next = curr;
					node = curr;
				}
			}
			pred = curr;
			curr = next;
			if (pred.val == curr.val) {
				pred = curr;
				curr = curr.next;
			}
			next = curr.next;
			if (next == null) {
				node.next = curr;
				node = curr;
			}
		}
		return head;
	}

	/**
	 * Given a binary tree, check whether it is a mirror of itself (i.e,
	 * symmetric around its center).
	 * 
	 * <pre>
	 * For example, this binary tree is symmetric:
	 * 
	 *     1									  1
	 *    / \  									 / \
	 *   2   2    But the following is not:		2   2
	 *  / \ / \ 								 \   \
	 * 3  4 4  3								  3   3
	 * </pre>
	 * 
	 * Note: Bonus points if you could solve it both recursively and
	 * iteratively.
	 */
	public boolean isSymmetric(TreeNode root) {
		if (root == null) {
			return true;
		}
		return isSymmetric(root.left, root.right);
	}

	private boolean isSymmetric(TreeNode left, TreeNode right) {
		if (left == null && right == null) {
			return true;
		}
		if (left == null && right != null) {
			return false;
		}
		if (left != null && right == null) {
			return false;
		}
		int lval = left.val;
		int rval = right.val;
		return (lval == rval) && isSymmetric(left.left, right.left)
				&& isSymmetric(left.right, right.right);
	}

	public boolean isSymmetricIter(TreeNode root) {
		if (root == null) {
			return true;
		}
		Queue<TreeNode> lqueue = new Queue<TreeNode>();
		Queue<TreeNode> rqueue = new Queue<TreeNode>();

		if (root.left != null) {
			lqueue.enqueue(root.left);
		}
		if (root.right != null) {
			rqueue.enqueue(root.right);
		}

		while (!lqueue.isEmpty() && !rqueue.isEmpty()) {
			TreeNode lnode = lqueue.dequeue().val;
			TreeNode rnode = rqueue.dequeue().val;
			if (lnode.val != rnode.val) {
				return false;
			}
			if (lnode.left != null) {
				lqueue.enqueue(lnode.left);
			}
			if (lnode.right != null) {
				lqueue.enqueue(lnode.right);
			}

			if (rnode.left != null) {
				rqueue.enqueue(rnode.left);
			}
			if (rnode.right != null) {
				rqueue.enqueue(rnode.right);
			}
		}
		return lqueue.isEmpty() && rqueue.isEmpty();
	}

	/**
	 * Given a binary tree, populate each next pointer to point to its next
	 * right node. If there is no next right node, the next pointer should be
	 * set to NULL. Initially, all next pointers are set to NULL.
	 * 
	 * Note: You may only use constant extra space. You may assume that it is a
	 * perfect binary tree (i.e, all leaves are at the same level, and every
	 * parent has two children).
	 * 
	 * <pre>
	 * For example,
	 * Given the following perfect binary tree,
	 * 
	 *          1
	 *        /  \
	 *       2    3
	 *      / \  / \
	 *     4  5  6  7
	 * 
	 * After calling your function, the tree should look like:
	 * 
	 *          1 -> NULL
	 *        /  \
	 *       2 -> 3 -> NULL
	 *      / \  / \
	 *     4->5->6->7 -> NULL
	 * </pre>
	 */
	public TreeLinkNode connect(TreeLinkNode root) {
		if (root == null) {
			return null;
		}
		root.next = null;
		connect(root.left, root.right);
		return root;
	}

	private TreeLinkNode connect(TreeLinkNode left, TreeLinkNode right) {
		if (left == null && right == null) {
			return null;
		}
		if (left != null && left.left == null && left.right == null) {
			return left;
		}
		if (right != null && right.left == null && right.right == null) {
			return right;
		}
		left.next = right;
		TreeLinkNode node = connect(left.left, left.right);
		node.next = connect(right.left, right.right);
		return node;
	}

	/**
	 * Previous problem but for any binary tree.
	 */
	public TreeLinkNode connect2(TreeLinkNode left, TreeLinkNode right) {
		if (left == null && right == null) {
			return null;
		} else if (right == null) {
			left.next = null;
		} else if (left == null) {
			right.next = null;
		} else {
			left.next = right;
		}
		if (left.left == null && left.right == null) {
			return left;
		}
		if (right.left == null && right.right == null) {
			return right;
		}
		TreeLinkNode node = null;
		if (left.left == null || left.right == null) {
			node = left.left == null ? left.right : left.left;
			node.next = connect(right.left, right.right);
		}
		if (right.left == null || right.right != null) {
			node = right.left == null ? right.right : right.left;
			connect(left.left, left.right).next = node;
		}
		if ((left.left != null && left.right != null)
				&& (right.left != null && right.right != null)) {
			node = connect(left.left, left.right);
			node.next = connect(right.left, right.right);
		}
		return node;
	}

	/**
	 * Given a 2d grid map of '1's (land) and '0's (water), count the number of
	 * islands. An island is surrounded by water and is formed by connecting
	 * adjacent lands horizontally or vertically. You may assume all four edges
	 * of the grid are all surrounded by water.
	 * 
	 * <pre>
	 * Example 1: (Answer: 1)			Example 2: (Answer 3)
	 * 11110							11000
	 * 11010							11000
	 * 11000							00100
	 * 00000							00011
	 * </pre>
	 * 
	 * The basic idea of the following solution is merging adjacent lands, and
	 * the merging should be done recursively.
	 */
	public int numIslands(char[][] grid) { // TODO
		if (grid == null || grid.length == 0 || grid[0].length == 0)
			return 0;
		int count = 0;

		for (int i = 0; i < grid.length; i++) {
			for (int j = 0; j < grid[0].length; j++) {
				if (grid[i][j] == '1') {
					count++;
					merge(grid, i, j);
				}
			}
		}
		return count;
	}

	public void merge(char[][] grid, int i, int j) {
		// validity checking
		if (i < 0 || j < 0 || i > grid.length - 1 || j > grid[0].length - 1)
			return;

		// if current cell is water or visited
		if (grid[i][j] != '1')
			return;

		// set visited cell to '0'
		grid[i][j] = '0';

		// merge all adjacent land
		merge(grid, i - 1, j);
		merge(grid, i + 1, j);
		merge(grid, i, j - 1);
		merge(grid, i, j + 1);
	}

	/**
	 * Clone an undirected graph. Each node in the graph contains a label and a
	 * list of its neighbors.
	 * 
	 * OJ's undirected graph serialization:
	 * 
	 * Nodes are labeled uniquely. We use # as a separator for each node, and ,
	 * as a separator for node label and each neighbor of the node.
	 * 
	 * As an example, consider the serialized graph {0,1,2#1,2#2,2}.
	 * 
	 * The graph has a total of three nodes, and therefore contains three parts
	 * as separated by #.
	 */

	public UndirectedGraphNode cloneGraph(UndirectedGraphNode node) {
		if (node == null || node.neighbors == null)
			return node;
		List<UndirectedGraphNode> dup = new ArrayList<UndirectedGraphNode>();
		UndirectedGraphNode copy = new UndirectedGraphNode(node.label);
		for (UndirectedGraphNode nod : node.neighbors) {
			if (!node.visited) {
				copy.neighbors.add(visit(nod, dup));
			}
		}
		return copy;
	}

	private UndirectedGraphNode visit(UndirectedGraphNode node,
			List<UndirectedGraphNode> dup) {
		if (node.neighbors == null)
			return null;
		node.visited = true;
		UndirectedGraphNode copy = new UndirectedGraphNode(node.label);
		if (dup.indexOf(node) < 0) {
			for (UndirectedGraphNode nod : node.neighbors) {
				if (node.visited && node.label != nod.label)
					continue;
				if (node.label == nod.label)
					dup.add(nod);
				copy.neighbors.add(visit(nod, dup));
			}
		}
		return copy;
	}

	/**
	 * Given n non-negative integers representing an elevation map where the
	 * width of each bar is 1, compute how much water it is able to trap after
	 * raining. @see https://leetcode.com/problems/trapping-rain-water/
	 * 
	 * For example, Given [0, 1, 0, 2, 1, 0, 1, 3, 2, 1, 2, 1], return 6
	 */
	public int trap(int[] height) { // TODO
		int i = 0, j = 0, vol = 0;
		int length = height.length;
		while (i < length && (j + 2) < length) {
			if (height[j] == 0 || height[i] >= height[j + 1]) {
				if (existTrap(height, j + 1)) {
					int min = Math.min(height[j], height[j + 2]);
					vol += (j + 1 - i) * (min - height[j + 1]);
				}
				j++;
				continue;
			}
			for (int k = i + 1; k < j + 1; k++) {
				if (existTrap(height, k)) {
					vol -= height[k];
				}
			}
			i = j++;
		}
		return vol;
	}

	private boolean existTrap(int[] height, int i) {
		return (0 < i && i < height.length - 1) && (height[i] < height[i - 1])
				&& (height[i] < height[i + 1]);
	}

	public int trap1(int[] height) {
		int result = 0;

		if (height == null || height.length <= 2)
			return result;

		int left[] = new int[height.length];
		int right[] = new int[height.length];

		// scan from left to right
		int max = height[0];
		left[0] = height[0];
		for (int i = 1; i < height.length; i++) {
			if (height[i] < max) {
				left[i] = max;
			} else {
				left[i] = height[i];
				max = height[i];
			}
		}

		// scan from right to left
		max = height[height.length - 1];
		right[height.length - 1] = height[height.length - 1];
		for (int i = height.length - 2; i >= 0; i--) {
			if (height[i] < max) {
				right[i] = max;
			} else {
				right[i] = height[i];
				max = height[i];
			}
		}

		// calculate total
		for (int i = 0; i < height.length; i++) {
			result += Math.min(left[i], right[i]) - height[i];
		}

		return result;
	}

	/**
	 * Given a list, rotate the list to the right by k places, where k is
	 * non-negative.
	 * 
	 * For example: Given 1->2->3->4->5->NULL and k = 2, return
	 * 4->5->1->2->3->NULL.
	 */
	public ListNode rotateRight(ListNode head, int k) {
		int n = numNodes(head), i = 0;
		if (k > n) {
			head = rotateRight(head, n);
			return rotateRight(head, k - n);
		}
		ListNode node = head, prev = node;
		while (node.next != null) {
			if (i < n - k) {
				prev = node;
				i++;
			}
			node = node.next;
		}
		node.next = head;
		head = prev.next;
		prev.next = null;
		return head;
	}

	private int numNodes(ListNode head) {
		if (head == null)
			return 0;
		int count = 0;
		ListNode node = head;
		while (node != null) {
			count++;
			node = node.next;
		}
		return count;
	}

	/**
	 * Sort Colors. Given an array with n objects colored red, white or blue,
	 * sort them so that objects of the same color are adjacent, with the colors
	 * in the order red, white and blue.
	 * 
	 * Here, we will use the integers 0, 1, and 2 to represent the color red,
	 * white, and blue respectively. Note: You are not suppose to use the
	 * library's sort function for this problem.
	 * 
	 * 1.) Two passes on the array.
	 */
	public void sortColors(int[] nums) {
		int[] colors = new int[3];
		for (int i = 0; i < nums.length; i++) {
			if (nums[i] == 0) {
				colors[0]++;
			}
			if (nums[i] == 1) {
				colors[1]++;
			}
			if (nums[i] == 2) {
				colors[2]++;
			}
		}
		int k = 0, occ = 0;
		for (int i = 0; i < nums.length; i++) {
			if (i < colors[k] + occ) {
				nums[i] = k;
			}
			occ += colors[k++];
		}
	}

	/**
	 * Same as before but with one pass and constant space.
	 */
	public void sortColors2(int[] nums) { // TODO
		int j = 0, k = nums.length - 1;
		int tmp = 0;
		for (int i = 0; i <= k; ++i) {
			if (nums[i] == 0) {
				tmp = nums[j];
				nums[j] = nums[i];
				nums[i] = tmp;
				j++;
			} else if (nums[i] == 2) {
				tmp = nums[k];
				nums[k] = nums[i];
				nums[i] = tmp;
				--k;
				--i;
			}
		}
	}

	/**
	 * Partition List: Given a linked list and a value x, partition it such that
	 * all nodes less than x come before nodes greater than or equal to x.
	 * 
	 * You should preserve the original relative order of the nodes in each of
	 * the two partitions.
	 * 
	 * For example, Given 1->4->3->2->5->2 and x = 3, return 1->2->2->4->3->5
	 */
	public ListNode partition(ListNode head, int x) {
		if (head == null)
			return null;
		ListNode pred = head;
		ListNode node = head;
		ListNode prev = null;
		while (node != null && node.val != x) {
			prev = node;
			node = node.next;
		}
		if (node == null || prev == null)
			return null;
		ListNode succ = node.next;
		while (pred.next != node && succ != null) {
			prev = succ;
			if (x <= pred.val && succ.val < x) {
				int tmp = pred.val;
				pred.val = succ.val;
				succ.val = tmp;
				pred = pred.next;
				succ = succ.next;
			} else if (x <= pred.val) {
				succ = succ.next;
			} else {
				pred = pred.next;
			}
		}
		while (succ != null) {
			if (succ.val < x) {
				prev = succ;
				succ = succ.next;
			} else {
				prev.next = succ.next;
				succ.next = node;
				pred.next = succ;
				pred = succ;
				succ = prev.next;
			}
		}
		while (pred.next != node) {
			if (pred.val >= x) {
				int tmp = pred.val;
				node.val = pred.val;
				pred.val = tmp;
			}
			pred = pred.next;
		}
		return head;
	}

	/**
	 * Given a string S and a string T, find the minimum window in S which will
	 * contain all the characters in T in complexity O(n).
	 * 
	 * <pre>
	 * For example,
	 * S = "ADOBECODEBANC"
	 * T = "ABC"
	 * 
	 * Minimum window is "BANC".
	 * </pre>
	 * 
	 * Note: If there is no such window in S that covers all characters in T,
	 * return the empty string "".
	 * 
	 * If there are multiple such windows, you are guaranteed that there will
	 * always be only one unique minimum window in S.
	 */
	public String minWindow(String s, String t) { // TODO
		return null;
	}

	/**
	 * Given a binary tree containing digits from 0-9 only, each root-to-leaf
	 * path could represent a number.
	 * 
	 * An example is the root-to-leaf path 1->2->3 which represents the number
	 * 123.
	 * 
	 * Find the total sum of all root-to-leaf numbers.
	 * 
	 * <pre>
	 * For example,
	 * 
	 *     1
	 *    / \
	 *   2   3
	 * </pre>
	 * 
	 * The root-to-leaf path 1->2 represents the number 12. The root-to-leaf
	 * path 1->3 represents the number 13.
	 * 
	 * Return the sum = 12 + 13 = 25.
	 */
	public int sumNumbers(TreeNode root) {
		Queue<TreeNode> q = new Queue<TreeNode>();
		q.enqueue(root);
		int sum = 0;
		root.visited = true;
		List<Integer> list = new ArrayList<Integer>();
		while (!q.isEmpty()) {
			TreeNode node = q.dequeue().val;
			if (node.left != null && !node.left.visited) {
				sum = sumNumbers(node.left, list, sum);
				q.enqueue(node.left);
			}
			if (node.right != null && !node.right.visited) {
				sum = sumNumbers(node.right, list, sum);
				q.enqueue(node.right);
			}
		}
		return sum;
	}

	public int sumNumbers(TreeNode root, List<Integer> list, int sum) {
		if (root == null) {
			sum += getInt(list);
			int last = list.size() - 1;
			if (last >= 0) {
				list.remove(last);
			}
			return sum;
		}
		root.visited = true;
		list.add(root.val);
		if (root.left != null) {
			return sumNumbers(root.left, list, sum);
		}
		if (root.right != null) {
			return sumNumbers(root.right, list, sum);
		}
		return sum;
	}

	private int getInt(List<Integer> list) {
		if (list == null || list.isEmpty()) {
			return 0;
		}
		StringBuilder builder = new StringBuilder();
		for (int num : list) {
			builder.append(num);
		}
		return Integer.valueOf(builder.toString());
	}

	/**
	 * Given a binary tree, find the maximum path sum. The path may start and
	 * end at any node in the tree.
	 * 
	 * <pre>
	 * For example, given the below binary tree,
	 * 
	 *        1
	 *       / \
	 *      2   3
	 * 
	 * Return 6.
	 * </pre>
	 */
	public int maxPathSum(Node<Integer> root) {
		Queue<Node<Integer>> queue = new Queue<Node<Integer>>();
		queue.enqueue(root);
		Map<Node<Integer>, Boolean> map;
		map = new HashMap<Node<Integer>, Boolean>();
		int lsum = 0, rsum = 0, rval = root.val;
		while (!queue.isEmpty()) {
			Node<Integer> node = queue.dequeue().val;
			map.put(node, true);
			if (node.left != null && !visited(map, node.left)) {
				queue.enqueue(node.left);
				lsum = maxPathSum(node.left, lsum, map);
				System.out.println("Max sum left subtree: " + lsum);
			}
			if (node.right != null && !visited(map, node.right)) {
				queue.enqueue(node.right);
				rsum = maxPathSum(node.right, rsum, map);
				System.out.println("Max sum right subtree: " + rsum);
			}
		}
		return (rval < 0) ? Math.max(lsum, rsum) : lsum + rval + rsum;
	}

	private boolean visited(Map<Node<Integer>, Boolean> visited,
			Node<Integer> node) {
		return visited.containsKey(node);
	}

	private int maxPathSum(Node<Integer> node, int sum,
			Map<Node<Integer>, Boolean> visited) {
		visited.put(node, true);
		if (node == null) {
			return sum;
		}
		if (node.left == null) {
			return maxPathSum(node.right, sum + node.val, visited);
		} else if (node.right == null) {
			return maxPathSum(node.left, sum + node.val, visited);
		} else {
			return Math.max(maxPathSum(node.left, sum + node.val, visited),
					maxPathSum(node.right, sum + node.val, visited));
		}
	}

	/**
	 * Flatten Binary Tree to Linked List.
	 */
	public QNode<Integer> flatten(TreeNode root) {
		if (root == null) {
			return null;
		}
		QNode<Integer> head = new QNode<Integer>(root.val);
		QNode<Integer> list = head;
		flatten(root, list);
		return head;
	}

	public void flatten(TreeNode root, QNode<Integer> list) {
		if (root != null) {
			QNode<Integer> node = new QNode<Integer>(root.val);
			list.next = node;
			flatten(root.left, list);
			flatten(root.right, list);
		}
	}

	/**
	 * Given a binary tree, determine if it is a valid binary search tree (BST).
	 * 
	 * Assume a BST is defined as follows:
	 * 
	 * The left subtree of a node contains only nodes with keys less than the
	 * node's key. The right subtree of a node contains only nodes with keys
	 * greater than the node's key. Both the left and right subtrees must also
	 * be binary search trees.
	 */
	public boolean isValidBST(Node<Integer> root) {
		if (root == null) {
			return true;
		}
		if ((root.left == null) && (root.right == null)) {
			return true;
		}
		if ((root.left != null)
				&& ((Integer) root.left.val > (Integer) root.val)) {
			return false;
		}
		if ((root.right != null)
				&& ((Integer) root.right.val <= (Integer) root.val)) {
			return false;
		}
		int nval = (Integer) root.val;
		int lval = (Integer) root.left.val;
		int rval = (Integer) root.right.val;

		return ((lval <= nval) && (nval < rval)) && isValidBST(root.left)
				&& isValidBST(root.right);
	}

	/**
	 * Implement an iterator over a binary search tree (BST). Your iterator will
	 * be initialized with the root node of a BST.
	 * 
	 * Calling next() will return the next smallest number in the BST.
	 * 
	 * Note: next() and hasNext() should run in average O(1) time and uses O(h)
	 * memory, where h is the height of the tree.
	 */
	public class BSTIterator {

		Stack<TreeNode> stack = new Stack<TreeNode>();

		public BSTIterator(TreeNode root) {
			populate(root);
		}

		private void populate(TreeNode node) {
			while (node != null) {
				stack.push(node);
				node = node.left;
			}
		}

		public boolean hasNext() {
			return !stack.isEmpty();
		}

		public int next() {
			if (stack.isEmpty()) {
				throw new UnsupportedOperationException();
			}
			TreeNode toReturn = stack.peek();
			TreeNode node = stack.pop();
			node = node.right;
			while (node != null) {
				stack.push(node);
				node = node.left;
			}
			return toReturn.val;
		}
	}

	/**
	 * Given a binary tree, return iteratively the postorder traversal of its
	 * nodes' values. For example: Given binary tree {1,#,2,3}, return [3,2,1].
	 */
	public List<Integer> postorderTraversal(TreeNode root) {
		return null;
	}

	/**
	 * Invert a binary tree.
	 * 
	 * <pre>
	 *      4                     4
	 *    /   \        to		/   \
	 *   2     7			   7     2
	 *  / \   / \		      / \   / \
	 * 1   3 6   9			 9   6 3   1
	 * 
	 * </pre>
	 */
	public TreeNode invertTree(TreeNode root) {
		if (root == null) {
			return root;
		}
		TreeNode left = root.left;
		TreeNode right = root.right;
		root.right = invertTree(left);
		root.left = invertTree(right);
		return root;
	}

	/**
	 * Two elements of a binary search tree (BST) are swapped by mistake.
	 * 
	 * Recover the tree without changing its structure. Note: A solution using
	 * O(n) space is pretty straightforward.
	 * 
	 * Devise a O(n) space and a O(1) constant space solution.
	 */
	public void recoverTree1(TreeNode root) { // O(n) solution
		List<Integer> nums = new ArrayList<Integer>();
		recoverTree1(root, nums);
		int i = 0;
		for (i = 1; i < nums.size(); i++) {
			int prev = nums.get(i - 1);
			int curr = nums.get(i);
			if (curr >= prev)
				continue;
			break;
		}
		for (int j = nums.size() - 1; j > i; j++) {
			int prev = nums.get(i - 1);
			int curr = nums.get(i);
			if (curr >= prev)
				continue;
			nums.set(j, nums.get(i - 1));
			nums.set(i - 1, curr);
			break;
		}
	}

	private void recoverTree1(TreeNode root, List<Integer> list) {
		if (root == null)
			return;
		recoverTree1(root.left, list);
		list.add(root.val);
		recoverTree1(root.right, list);
	}

	/**
	 * Given a singly linked list where elements are sorted in ascending order,
	 * convert it to a height balanced BST.
	 * 
	 * [ 1, 2, 3, 4, 5, 6, 7, 8 ]
	 */
	public TreeNode sortedListToBST(ListNode head) {
		if (head == null) {
			return null;
		}
		TreeNode node = new TreeNode(head.val);
		if (head.next == null) {
			return node;
		}
		int n = 0;
		ListNode curr = head;
		while (curr != null) {
			n++;
			curr = curr.next;
		}
		return sortedListToBST(head, 0, n - 1);
	}

	private TreeNode sortedListToBST(ListNode node, int lb, int ub) {
		if (lb > ub) { // TODO
			return null;
		}
		int middle = lb + (ub - lb) / 2;
		TreeNode left = sortedListToBST(node, lb, middle - 1);
		TreeNode parent = new TreeNode(node.val);
		parent.left = left;
		node = node.next;
		parent.right = sortedListToBST(node, middle + 1, ub);
		return parent;
	}

	/**
	 * Two elements of a binary search tree (BST) are swapped by mistake.
	 * 
	 * Recover the tree without changing its structure. Devise an algorithm
	 * using O(1) constant space solution.
	 * 
	 * <pre>
	 *      4                     4
	 *    /   \        to		/   \
	 *   9     7			   2     7
	 *  / \   / \		      / \   / \
	 * 1   3 6   2			 1   3 6   9
	 * 
	 * </pre>
	 */
	public void recoverTree2(TreeNode root) {
		if (root == null)
			return;
		TreeNode first = null;
		TreeNode second = null;
		TreeNode pred = null;
		recoverTree2(root, first, second, pred);
	}

	private void recoverTree2(TreeNode root, TreeNode first, TreeNode second,
			TreeNode pred) {
		if (root == null)
			return;
		recoverTree2(root.left, first, second, pred);
		if (pred == null)
			pred = root;
		else {
			if (pred.val > root.val) {
				first = root;
				second = pred;
			}
			pred = root;
		}
		recoverTree2(root.right, first, second, pred);
	}

	public boolean isValidBST(TreeNode root) {
		if (root == null) {
			return true;
		}
		Node<Integer> node = new Node<Integer>(root.val);
		if (root.left != null) {
			node.left = new Node<Integer>(root.left.val);
		}
		if (root.right != null) {
			node.right = new Node<Integer>(root.right.val);
		}
		node.visited = root.visited;
		return isValidBST(node);
	}

	/**
	 * Given a binary tree, determine if it is height-balanced.
	 * 
	 * For this problem, a height-balanced binary tree is defined as a binary
	 * tree in which the depth of the two subtrees of every node never differ by
	 * more than 1.
	 */
	public boolean isBalanced(TreeNode root) {
		if (root == null) {
			return true;
		}
		int ldepth = maxDepth(root.left);
		int rdepth = maxDepth(root.right);
		return (Math.abs(rdepth - ldepth) <= 1) && isBalanced(root.left)
				&& isBalanced(root.right);
	}

	/**
	 * Given a binary tree, return iteratively the preorder traversal of its
	 * nodes' values. For example: Given binary tree {1,#,2,3}, return [1,2,3].
	 */
	public List<Integer> preorderTraversal(TreeNode root) {
		List<Integer> vals = new ArrayList<Integer>();
		Stack<TreeNode> stack = new Stack<TreeNode>();
		TreeNode node = root;
		while (!stack.isEmpty() || node != null) {
			if (node != null) {
				vals.add(node.val);
				if (node.right != null) {
					stack.push(node.right);
				}
				node = node.left;
			} else {
				node = stack.pop();
			}
		}
		return vals;
	}

	/**
	 * Given a binary tree, find its minimum depth.
	 * 
	 * The minimum depth is the number of nodes along the shortest path from the
	 * root node down to the nearest leaf node.
	 */
	public int minDepth(TreeNode root) {
		if (root == null) {
			return 0;
		}
		return 1 + Math.min(minDepth(root.left), minDepth(root.right));
	}

	public int maxDepth(TreeNode root) {
		if (root == null) {
			return 0;
		}
		return 1 + Math.max(maxDepth(root.left), maxDepth(root.right));
	}

	/**
	 * Binary Tree Right Side View: Given a binary tree, imagine yourself
	 * standing on the right side of it, return the values of the nodes you can
	 * see ordered from top to bottom. For example: Given the following binary
	 * tree,
	 */
	public List<Integer> rightSideView(TreeNode root) {
		List<Integer> list = new ArrayList<Integer>();
		while (root != null) {
			list.add(root.val);
			if (root.right != null) {
				root = root.right;
			} else {
				root = root.left;
			}
		}
		return list;
	}

	/**
	 * Remove element: Given an array and a value, remove all instances of that
	 * value in place and return the new length. The order of elements can be
	 * changed. It doesn't matter what you leave beyond the new length.
	 */
	public int removeInstance(int array[], int n) {
		return removeInstance(array, n, 0, array.length);
	}

	private int removeInstance(int array[], int n, int s, int l) {
		if (array == null)
			throw new IllegalArgumentException();
		if (s == array.length) {
			return l;
		}
		if (array[s] == n) {
			l--;
			int k = s, size = array.length - 1;
			while (k < size && array[k + 1] != n) {
				array[k] = array[k++];
			}
			if (k < size) {
				return removeInstance(array, n, k + 1, l);
			}
		}
		return removeInstance(array, n, s + 1, l);
	}

	public <T> void printInOrder(Node<T> root) {
		if (root == null) {
			return;
		}
		printInOrder(root.left);
		System.out.print(root.val + " ");
		printInOrder(root.right);
	}

	public void printInOrder(TreeNode root) {
		if (root == null) {
			return;
		}
		printInOrder(root.left);
		System.out.print(root.val + " ");
		printInOrder(root.right);
	}

	public <T> void printPreOrder(Node<T> root) {
		if (root == null) {
			return;
		}
		System.out.print(root.val + " ");
		printInOrder(root.left);
		printInOrder(root.right);
	}

	/**
	 * Given n non-negative integers representing the histogram's bar height
	 * where the width of each bar is 1, find the area of largest rectangle in
	 * the histogram.
	 */
	public int largestHistogram(double[] heights) {
		if (heights == null)
			throw new IllegalArgumentException();
		MergeSort sort = new MergeSort(heights);
		sort.sort(Order.ASC);
		double max = Integer.MIN_VALUE;
		int length = heights.length;
		for (int i = 0; i < length; i++) {
			double val = heights[i] * (length - i);
			if (val > max) {
				max = val;
			}
		}
		return (int) max;
	}

	/**
	 * Given pre-order and in-order traversal of a tree, construct the binary
	 * tree. Note: You may assume that duplicates do not exist in the tree.
	 * 
	 * <pre>
	 * Example: 
	 * Pre-order: F, B, A, D, C, E, G, I, H 
	 * In-order : A, B, C, D, E, F, G, H, I
	 * </pre>
	 */
	public Node<Integer> buildTree1(int[] preorder, int[] inorder) {
		return buildTree1(inorder, preorder, 0, inorder.length - 1);
	}

	public Node<Integer> buildTree1(int[] inorder, int[] preorder, int start,
			int end) { // TODO
		if (start > end) {
			return null;
		}
		/*
		 * Pick current node from Preorder traversal using preIndex and
		 * increment preIndex
		 */
		Node<Integer> node = new Node<Integer>(preorder[preIndex++]);
		if (start == end) {/* If this node has no children then return */
			return node;
		}
		/* Else find the index of this node in Inorder traversal */
		int inIndex = search(inorder, start, end, node.val);

		/*
		 * 77 144 76 68 Using index in Inorder traversal, construct left and
		 * right subtrees
		 */
		node.left = buildTree1(inorder, preorder, start, inIndex - 1);
		node.right = buildTree1(inorder, preorder, inIndex + 1, end);
		return node;
	}

	private int search(int array[], int start, int end, int value) {
		for (int i = 0; i <= end; i++) {
			if (array[i] == value) {
				return i;
			}
		}
		return -1;
	}

	/**
	 * Given in-order and post-order traversal of a tree, construct the binary
	 * tree. Note: You may assume that duplicates do not exist in the tree.
	 * 
	 * <pre>
	 * Example:
	 * In-order  : A, B, C, D, E, F, G, H, I
	 * Post-order: A, C, E, D, B, H, I, G, F
	 * Pre-order : F, B, A, D, C, E, G, I, H
	 * </pre>
	 */
	public Node<Integer> buildTree2(int[] inorder, int[] postorder) {
		if (postorder == null || inorder == null) {
			return null;
		}
		return buildTree2(inorder, 0, inorder.length - 1, postorder, 0,
				postorder.length - 1);
	}

	public Node<Integer> buildTree2(int[] inorder, int inStart, int inEnd,
			int[] postorder, int postStart, int postEnd) {
		if (inStart > inEnd || postStart > postEnd)
			return null;

		int rootValue = postorder[postEnd];
		Node<Integer> root = new Node<Integer>(rootValue);

		int i = 0;
		for (; i <= inorder.length; i++) {
			if (inorder[i] == rootValue)
				break;
		}
		root.left = buildTree2(inorder, inStart, i - 1, postorder, postStart,
				postStart + i - (inStart + 1));
		// Because k is not the length, it needs to -(inStart+1) to get the
		// length
		root.right = buildTree2(inorder, i + 1, inEnd, postorder, postStart + i
				- inStart, postEnd - 1);
		return root;
	}

	/**
	 * Find the contiguous subarray within an array (containing at least one
	 * number) which has the largest product.
	 * 
	 * For example, given the array [2, 3, -2, 4], the contiguous subarray [2,3]
	 * has the largest product = 6.
	 */
	public int[] contiguousArrayLargestProduct(int[] array) { // TODO
		if (array == null)
			return null;
		int i = 0, beg = 0, end = 0;
		int prod = 1, max = 1, min = 1;
		for (i = 0; i < array.length; i++) {
			if (array[i] > 0) {
				if (prod == 1) { // i.e. we start again
					beg = i;
				}
				prod *= array[i];
				min = Math.min(min * array[i], 1);
			} else if (array[i] == 0) {
				min = 1;
				prod = 1;
			} else {
				int temp = prod;
				prod = Math.max(min * array[i], 1);
				min = temp * array[i];
			}
			if (prod > max) { // we found a new max
				end = i;
				max = prod;
			}
		}
		return new int[] { beg, end };
	}

	/**
	 * Given n non-negative integers a1, a2, ..., an, where each represents a
	 * point at coordinate (i, ai). n vertical lines are drawn such that the two
	 * endpoints of line i is at (i, ai) and (i, 0). Find two lines, which
	 * together with x-axis forms a container, such that the container contains
	 * the most water.
	 * 
	 * Note: You may not slant the container.
	 */
	public int maxArea(int[] height) {
		if (height == null)
			throw new IllegalArgumentException();
		int max = Integer.MIN_VALUE;
		for (int i = 1; i < height.length - 1; i++) {
			int largest = Math.max(height[i - 1], height[i + 1]);
			int smallest = Math.min(largest, height[i]);
			if (max < smallest) {
				max = smallest;
			}
		}
		return max;
	}

	/**
	 * Given a set of candidate numbers (C) and a target number (T), find all
	 * unique combinations in C where the candidate numbers sum to T.
	 * 
	 * The same repeated number may be chosen from C unlimited number of times.
	 * 
	 * Note:
	 * 
	 * All numbers (including target) will be positive integers. Elements in a
	 * combination (a1, a2, … , ak) must be in non-descending order. (i.e, a1 ≤
	 * a2 ≤ … ≤ ak). The solution set must not contain duplicate combinations.
	 * 
	 * For example, given candidate set 2,3,6,7 and target 7, A solution set is:
	 * [7] [2, 2, 3]
	 * 
	 * Backtracking (think of DFS).
	 */
	public List<List<Double>> combinationSum1(double[] candidates, double target) {
		if (candidates == null)
			return null;
		new MergeSort(candidates).sort();
		List<List<Double>> allRslts = new ArrayList<List<Double>>();
		for (int i = 0; i < candidates.length; i++) {
			List<Double> list = new ArrayList<Double>();
			if (target == candidates[i]) {
				list.add(target);
				allRslts.add(list);
				continue;
			}
			combine1(candidates, i, target, list, allRslts);
		}
		return allRslts;
	}

	private void combine1(double[] candidates, int start, double target,
			List<Double> list, List<List<Double>> allRslts) {
		if (target >= 0) {
			list.add(candidates[start]);
			double diff = target - candidates[start];
			if (exists(candidates, start, diff)) {
				list.add(diff);
				if (!allRslts.contains(list)) {
					allRslts.add(new ArrayList<Double>(list));
				}
				int last = list.size() - 1;
				list.remove(last);
			}
			combine1(candidates, start, diff, list, allRslts);
		}
	}

	private boolean exists(double[] candidates, int start, double value) {
		if (candidates == null || candidates.length == 0) {
			return false;
		}
		for (int i = start; i < candidates.length; i++) {
			if (candidates[i] == value) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Find all possible combinations of k numbers that add up to a number n,
	 * given that only numbers from 1 to 9 can be used and each combination
	 * should be a unique set of numbers.
	 * 
	 * Ensure that numbers within the set are sorted in ascending order.
	 * 
	 * <pre>
	 * Example 1:
	 * 
	 * Input: k = 3, n = 7 Output: [[1,2,4]]
	 * 
	 * Example 2:
	 * 
	 * Input: k = 3, n = 9
	 * 
	 * Output: [[1,2,6], [1,3,5], [2,3,4]]
	 * 
	 * </pre>
	 */
	public List<List<Double>> combinationSum3(int k, int target) {
		double[] candidates = new double[9];
		for (int i = 1; i <= 9; i++) {
			candidates[i - 1] = i;
		}
		List<Double> list = new ArrayList<Double>();
		List<List<Double>> results = new ArrayList<List<Double>>();
		combine3(candidates, 0, k, target, list, results);
		return results;
	}

	private void combine3(double[] candidates, int start, int size,
			double target, List<Double> list, List<List<Double>> allRslts) {
		if (target >= 0 && start < candidates.length) {
			if (!list.contains(candidates[start])) {
				list.add(candidates[start]);
			} else {
				combine3(candidates, start + 1, size, target, list, allRslts);
			}
			double diff = target - candidates[start];
			if (exists(candidates, start, diff)) {
				list.add(diff);
				if (list.size() == size) {
					allRslts.add(new ArrayList<Double>(list));
				}
				int last = list.size() - 1;
				list.remove(last);
			}
			combine3(candidates, start, size, diff, list, allRslts);
		}
	}

	/**
	 * Given a collection of candidate numbers (C) and a target number (T), find
	 * all unique combinations in C where the candidate numbers sums to T.
	 * 
	 * Each number in C may only be used once in the combination.
	 * 
	 * Note:
	 * 
	 * All numbers (including target) will be positive integers. Elements in a
	 * combination (a1, a2, … , ak) must be in non-descending order. (i.e, a1 ≤
	 * a2 ≤ … ≤ ak). The solution set must not contain duplicate combinations.
	 * 
	 * For example, given candidate set 10,1,2,7,6,1,5 and target 8, A solution
	 * set is: [1, 7] [1, 2, 5] [2, 6] [1, 1, 6]
	 */
	public List<List<Double>> combinationSum2(double[] candidates, int target) {
		if (candidates == null || candidates.length == 0)
			return null;
		new MergeSort(candidates).sort();
		List<List<Double>> results = new ArrayList<List<Double>>();
		for (int i = 0; i < candidates.length; i++) {
			List<Double> list = new ArrayList<Double>();
			List<Double> visited = new ArrayList<Double>();
			combine2(candidates, i, target, visited, list, results);
		}
		return results;
	}

	private void combine2(double[] candidates, int start, double target,
			List<Double> visited, List<Double> list, List<List<Double>> allRslts) {
		if (target >= 0 && start < candidates.length) {
			if (!visited.contains(candidates[start])) {
				list.add(candidates[start]);
				visited.add(candidates[start]);
			} else {
				combine2(candidates, start + 1, target, visited, list, allRslts);
			}
			double diff = target - candidates[start];
			if (exists(candidates, start, diff, visited)) {
				list.add(diff);
				visited.add(diff);
				if (!allRslts.contains(list)) {
					allRslts.add(new ArrayList<Double>(list));
				}
				int last = list.size() - 1;
				list.remove(last);
			}
			combine2(candidates, start, diff, visited, list, allRslts);
		}
	}

	private boolean exists(double[] candidates, int start, double value,
			List<Double> visited) {
		if (candidates == null || candidates.length == 0) {
			return false;
		}
		for (int i = start; i < candidates.length; i++) {
			if (candidates[i] == value && !visited.contains(value)) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Say you have an array for which the i-th element is the price of a given
	 * stock on day i.
	 * 
	 * If you were only permitted to complete at most one transaction (i.e, buy
	 * one and sell one share of the stock), design an algorithm to find the
	 * maximum profit.
	 * 
	 * Quadratic time O(n*n). See below for linear time.
	 */
	public int maxProfit(int[] prices) {
		if (prices == null) {
			return 0;
		}
		int maxProfit = Integer.MIN_VALUE;
		for (int i = 0; i < prices.length; i++) {
			int max = Integer.MIN_VALUE;
			for (int j = 0; j < prices.length; j++) {
				int temp = prices[j] - prices[i];
				if (max < temp) {
					max = temp;
				}
			}
			if (max > maxProfit) {
				maxProfit = max;
			}
		}
		return maxProfit;
	}

	public int maxProfit1(int[] prices) { // TODO
		if (prices == null) {
			return 0;
		}
		int min = Integer.MAX_VALUE;
		int profit = 0;
		for (int i = 0; i < prices.length; i++) {
			profit = Math.max(profit, prices[i] - min);
			// The logic to making profit is to buy an
			// stock low and sell it at a high price
			min = Math.min(min, prices[i]);
		}
		return profit;
	}

	/**
	 * Say you have an array for which the i-th element is the price of a given
	 * stock on day i.
	 * 
	 * Design an algorithm to find the maximum profit. You may complete at most
	 * two transactions.
	 * 
	 * Note: You may not engage in multiple transactions at the same time (i.e,
	 * you must sell the stock before you buy again).
	 */
	public int maxProfit3(int[] prices) {
		if (prices == null) {
			return 0;
		}
		int profitT1 = 0, profitT2 = 0;
		int min = Integer.MAX_VALUE;
		for (int i = 0; i < prices.length; i++) {
			if (profitT1 < prices[i] - min) {
				profitT2 = profitT1;
			}
			profitT1 = Math.max(profitT1, prices[i] - min);
			min = Math.min(min, prices[i]);
		}
		return profitT1 + profitT2;
	}

	/**
	 * Say you have an array for which the i-th element is the price of a given
	 * stock on day i.
	 * 
	 * Design an algorithm to find the maximum profit. You may complete as many
	 * transactions as you like (i.e, buy one and sell one share of the stock
	 * multiple times). However, you may not engage in multiple transactions at
	 * the same time (i.e, you must sell the stock before you buy again).
	 */
	public int maxProfit2(int[] prices) { // TODO
		return 0;
	}

	public int maxProfit(int k, int[] prices) {
		if (prices.length < 2 || k <= 0)
			return 0;
		int[] local = new int[k + 1];
		int[] global = new int[k + 1];

		for (int i = 0; i < prices.length - 1; i++) {
			int diff = prices[i + 1] - prices[i];
			for (int j = k; j >= 1; j--) {
				local[j] = Math.max(global[j - 1] + Math.max(diff, 0), local[j]
						+ diff);
				global[j] = Math.max(local[j], global[j]);
			}
		}

		return global[k];
	}

	/**
	 * Implement pow(x, n). Complexity required O(log(n)).
	 */
	public double pow(double x, int n) {
		if (n == 0) {
			return 1;
		}
		if (n == 1) {
			return x;
		}
		if (n < 0) {
			return pow(1.0 / x, (-1 * n) / 2) * pow(1.0 / x, (-1 * n) / 2);
		}
		if (n % 2 == 0) {
			return pow(x, n / 2) * pow(x, n / 2);
		} else {
			return pow(x, n / 2) * pow(x, (n - 1) / 2);
		}
	}

	/**
	 * Suppose a sorted array is rotated at some pivot unknown to you
	 * beforehand. (i.e., 0 1 2 4 5 6 7 might become 4 5 6 7 0 1 2).
	 * 
	 * <pre>
	 * Find the minimum element in O(log(n)).
	 *  1. You may assume no duplicate exists in the array. 
	 *  2. The array may contain duplicates.
	 * </pre>
	 */
	public int findMin(int[] nums) {
		return findMin(nums, 0, nums.length - 1);
	}

	private int findMin(int[] nums, int start, int end) {
		if (start < end + 1) {
			int middle = start + (end - start) / 2;
			if (nums[middle] < Math.min(nums[middle + 1], nums[middle - 1])) {
				return nums[middle];
			}
			if (nums[middle] > Math.max(nums[middle - 1], nums[middle + 1])) {
				return findMin(nums, middle + 1, end);
			} else if (nums[middle - 1] < nums[middle]
					&& nums[middle] < nums[middle + 1]) {
				return findMin(nums, start, middle - 1);
			}
		}
		return -1;
	}

	/**
	 * There are two sorted arrays nums1 and nums2 of size m and n respectively.
	 * Find the median of the two sorted arrays. The overall run time complexity
	 * should be O(log (m+n)).
	 */
	public double findMedianSortedArrays(int[] A, int[] B) {
		int m = A.length;
		int n = B.length;

		if ((m + n) % 2 != 0) // odd
			return (double) findKth(A, B, (m + n) / 2, 0, m - 1, 0, n - 1);
		else { // even
			return (findKth(A, B, (m + n) / 2, 0, m - 1, 0, n - 1) + findKth(A,
					B, (m + n) / 2, 0, m - 1, 0, n - 1)) * 0.5;
		}
	}

	public static int findKth(int A[], int B[], int k, int aStart, int aEnd,
			int bStart, int bEnd) { // TODO

		int aLen = aEnd - aStart + 1;
		int bLen = bEnd - bStart + 1;

		// Handle special cases
		if (aLen == 0)
			return B[bStart + k];
		if (bLen == 0)
			return A[aStart + k];
		if (k == 0)
			return A[aStart] < B[bStart] ? A[aStart] : B[bStart];

		int aMid = aLen * k / (aLen + bLen); // a's middle count
		int bMid = k - aMid - 1; // b's middle count

		// make aMid and bMid to be array index
		aMid = aMid + aStart;
		bMid = bMid + bStart;

		if (A[aMid] > B[bMid]) {
			k = k - (bMid - bStart + 1);
			aEnd = aMid;
			bStart = bMid + 1;
		} else {
			k = k - (aMid - aStart + 1);
			bEnd = bMid;
			aStart = aMid + 1;
		}
		return findKth(A, B, k, aStart, aEnd, bStart, bEnd);
	}

	/**
	 * Given a complete binary tree, count the number of nodes. Definition of a
	 * complete binary tree from Wikipedia: In a complete binary tree every
	 * level, except possibly the last, is completely filled, and all nodes in
	 * the last level are as far left as possible. It can have between 1 and 2h
	 * nodes inclusive at the last level h.
	 */
	public int countNodes(TreeNode root) {
		if (root == null) {
			return 0;
		}
		return countNodes(root.left) + countNodes(root.right) + 1;
	}

	public TreeNode buildTree(int[] array) {
		TreeNode root = buildTree(array, 0, array.length - 1);
		return root;
	}

	public TreeNode buildTree(int[] array, int start, int end) {
		if (array == null || start > end) {
			return null;
		}
		int middle = start + (end - start) / 2;
		TreeNode node = new TreeNode(array[middle]);
		node.left = buildTree(array, start, middle - 1);
		node.right = buildTree(array, middle + 1, end);
		return node;
	}

	/**
	 * Implement int sqrt(int x). Compute and return the square root of x.
	 */
	public int mySqrt(int x) { // TODO
		int floor = 0;
		int ceiling = x;
		return mySqrt(x, floor, ceiling);
	}

	private int mySqrt(int x, int start, int end) {
		while (start != (end - 1)) {
			int middle = start + (end - start) / 2;
			if (middle * middle > x) {
				end = middle;
			} else if (middle * middle < x) {
				start = middle;
			} else {
				return middle;
			}
		}
		return start;
	}

	/**
	 * Divide two integers without using multiplication, division and mod
	 * operator.
	 * 
	 * If it is overflow, return MAX_INT.
	 */
	public int divide(int dividend, int divisor) { // TODO
		int count = 1;
		while (count * divisor <= dividend) {
			count *= 2;
		}
		int result = 0;
		while (count != 1) {
			count = count / 2;
			if ((result + count) * divisor > dividend) {
				continue;
			}
			result = result + count;
		}
		return result;
	}

	/**
	 * Given a set of non-overlapping intervals, insert a new interval into the
	 * intervals (merge if necessary).
	 * 
	 * You may assume that the intervals were initially sorted according to
	 * their start times.
	 * 
	 * Example 1: Given intervals [1,3],[6,9], insert and merge [2,5] in as
	 * [1,5],[6,9].
	 * 
	 * Example 2: Given [1,2],[3,5],[6,7],[8,10],[12,16], insert and merge [4,9]
	 * in as [1,2],[3,10],[12,16].
	 * 
	 * This is because the new interval [4,9] overlaps with [3,5],[6,7],[8,10].
	 */

	public List<int[]> insert(List<int[]> intervals, int[] interval) {
		if (intervals == null || interval == null) {
			return null;
		}
		int i = 0, j = 0;
		List<int[]> rslts = new ArrayList<int[]>();
		for (; i < intervals.size(); i++) {
			int[] intv = intervals.get(i);
			if (interval[0] > intv[1]) {
				rslts.add(intv);
			} else {
				break;
			}
		}
		if (i >= intervals.size()) {
			return rslts;
		}
		int ubound = intervals.get(i)[1];
		for (j = i; j < intervals.size(); j++) {
			int[] intv = intervals.get(j);
			if (interval[1] <= intv[0]) {
				ubound = Math.max(interval[1], intervals.get(j - 1)[1]);
				break;
			} else if (interval[1] <= intv[1]) {
				ubound = intv[1];
			} else {
				ubound = interval[1];
			}
		}
		rslts.add(new int[] { intervals.get(i)[0], ubound });
		if (j < intervals.size()) {
			for (int k = j; k < intervals.size(); k++) {
				rslts.add(intervals.get(k));
			}
		}
		return rslts;
	}

	/**
	 * Given an unsorted integer array, find the first missing positive integer.
	 * 
	 * For example, Given [1, 2, 0] return 3, and [3, 4, -1, 1] return 2.
	 * 
	 * Your algorithm should run in O(n) time and uses constant space.
	 */
	public int firstMissingPositive(int[] nums) { // TODO
		if (nums == null)
			throw new IllegalArgumentException();
		int i = 0;
		while (i < nums.length) {
			if (nums[i] != i + 1 && nums[i] >= 1 && nums[i] <= nums.length
					&& nums[nums[i] - 1] != nums[i]) {
				int tmp = nums[nums[i] - 1];
				nums[nums[i] - 1] = nums[i];
				nums[i] = tmp;
			} else {
				i++;
			}
		}
		for (i = 0; i < nums.length; i++) {
			if (nums[i] != i + 1) {
				return i + 1;
			}
		}
		return nums.length + 1;
	}

	/**
	 * Find peak element : A peak element is an element that is greater than its
	 * neighbors. Given an input array where num[i] ≠ num[i+1], find a peak
	 * element and return its index. The array may contain multiple peaks, in
	 * that case return the index to any one of the peaks is fine.
	 * 
	 * You may imagine that num[-1] = num[n] = -∞.
	 * 
	 * For example, in array [1, 2, 3, 1], 3 is a peak element and your function
	 * should return the index number 2.
	 * 
	 * You should do this in logarithmic time.
	 */
	public int getPeak(int[] nums) {
		if (nums == null)
			throw new IllegalArgumentException();
		return getPeak(nums, 0, nums.length - 1);
	}

	private int getPeak(int[] nums, int start, int end) {
		if (start < end) {
			int middle = start + (end - start) / 2;
			if (middle == 0) {
				if (nums[middle] > nums[middle + 1]) {
					return middle;
				}
			} else if (nums[middle] > Math.max(nums[middle - 1],
					nums[middle + 1])) {
				return middle;
			} else if (nums[middle] < nums[middle - 1]) {
				return getPeak(nums, start, middle - 1);
			} else {
				return getPeak(nums, middle + 1, end);
			}
		}
		return -1;
	}

	/**
	 * Maximum sub-array: Find the contiguous subarray within an array
	 * (containing at least one number) which has the largest sum.
	 * 
	 * For example, given the array [−2, 1, −3, 4, −1, 2, 1, −5, 4], the
	 * contiguous subarray [4, −1, 2, 1] has the largest sum = 6.
	 */
	public int[] contiguousArrayLargestSum(int[] array) { // TODO
		if (array == null)
			return null;
		int beg = 0, end = 0;
		int max = array[0], sum = array[0];
		for (int i = 1; i < array.length; i++) {
			int val = sum + array[i];
			if (val > 0) {
				if (sum == 0) {
					beg = i;
				}
				sum = val;
			} else {
				sum = 0;
			}
			if (sum > max) {
				end = i;
				max = sum;
			}
		}
		return new int[] { beg, end };
	}

	/**
	 * Median of two sorted arrays: There are two sorted arrays nums1 and nums2
	 * of size m and n respectively. Find the median of the two sorted arrays.
	 * The overall run time complexity should be O(log (m+n)).
	 */
	public int median(int[] lhs, int rhs[]) {
		if (lhs == null || rhs == null)
			throw new IllegalArgumentException();
		int i = 0, j = 0, k = 0, median = 0;
		int middle = (lhs.length + rhs.length - 2) / 2;
		while (i < middle) {
			if (lhs[j] < rhs[k]) {
				median = lhs[++j];
			} else if (lhs[j] > rhs[k]) {
				median = rhs[++k];
			} else {
				k++;
				j++;
				median = lhs[j];
			}
			i++;
		}
		return median;
	}

	/**
	 * Maximal Rectangle: Given a 2D binary matrix filled with 0's and 1's, find
	 * the largest rectangle containing all ones and return its area.
	 * 
	 * <pre>
	 * Example:
	 * [
	 *   [0, 0, 0, 0],
	 *   [0, 1, 1, 0],
	 *   [0, 1, 1, 0]
	 * ]
	 * </pre>
	 */
	public List<int[]> largestRectangle(int[][] matrix) {
		if (matrix == null)
			throw new IllegalArgumentException();
		int min[] = new int[matrix.length];
		int max[] = new int[matrix.length];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				if (matrix[i][j] == 1) {
					min[i] = j;
					break;
				}
				min[i] = -1;
			}
			if (min[i] < 0) {
				continue;
			}
			for (int j = min[i]; j < matrix[0].length; j++) {
				if (matrix[i][j] == 1) {
					max[i] = j;
				}
			}
		}
		int maxI = -1, maxJ = 0;
		int minI = -1, minJ = matrix[0].length;
		for (int i = 0; i < matrix.length; i++) {
			if (min[i] >= 0 && min[i] < minJ) {
				minJ = min[i];
			}
			if (min[i] >= 0 && minI < 0) {
				minI = i;
			}
			if (max[i] >= maxJ) {
				maxJ = max[i];
			}
			if (matrix[i][max[i]] == 1) {
				maxI = i;
			}
		}

		if (minI < 0 || maxI < 0)
			return null;

		List<int[]> list = new ArrayList<int[]>();
		list.add(new int[] { minI, minJ });
		list.add(new int[] { maxI, maxJ });
		return list;
	}

	/**
	 * Longest Consecutive Sequence: Given an unsorted array of integers, find
	 * the length of the longest consecutive elements sequence.
	 * 
	 * For example, Given [100, 4, 200, 1, 3, 2], The longest consecutive
	 * elements sequence is [1, 2, 3, 4]. Return its length: 4.
	 * 
	 * Your algorithm should run in O(n) complexity.
	 */
	public int longestConsecutiveSequence(int[] array) {
		if (array == null)
			throw new IllegalArgumentException();
		int pivot = pivot(array, 0, array.length - 1);
		int length = getLength(array, 0, pivot, array[pivot], -1)
				+ getLength(array, pivot, array.length - 1, array[pivot], 1);
		return length + 1;
	}

	private int getLength(int[] array, int lb, int ub, int x, int offset) {
		Map<Integer, Boolean> map = new HashMap<Integer, Boolean>();
		int length = 0;
		for (int i = 1; i <= (ub - lb); i++) {
			int val = x + (offset * i);
			map.put(val, false);
		}
		for (int i = lb; i <= ub; i++) {
			if (map.containsKey(array[i]) && !map.get(array[i])) {
				length++;
				map.put(array[i], true);
			}
		}
		return length;
	}

	private int pivot(int[] array, int lb, int ub) {
		if (array == null)
			throw new IllegalArgumentException();
		int x = array[ub];
		int i = lb - 1;
		for (int j = lb; j < ub; j++) {
			if (array[j] <= x) {
				i++;
				int temp = array[j];
				array[j] = array[i];
				array[i] = temp;
			}
		}
		int temp = array[i + 1];
		array[i + 1] = array[ub];
		array[ub] = temp;
		return i + 1;
	}

	/**
	 * Rotate an array of n elements to the right by k steps. For example, with
	 * n = 7 and k = 3, the array [1, 2, 3, 4, 5, 6, 7] is rotated to [5, 6, 7,
	 * 1, 2, 3, 4]. Could you do it in-place with O(1) extra space?
	 */
	public int[] rotate(int[] array, int n) {
		if (array == null)
			return null;
		int i = 0;
		int j = array.length - (n % array.length);
		while (true) {
			if (j >= array.length) {
				break;
			}
			int temp = array[i];
			array[i] = array[j];
			for (int k = j; k > i; k--) {
				array[k] = array[k - 1];
			}
			array[i + 1] = temp;
			i++;
			j++;
		}
		return array;
	}

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
	}

	private class QNode<T> {
		QNode<T> next;
		T val;

		public QNode(T val) {
			this.val = val;
		}
	}

	public class TreeNode {
		int val;
		TreeNode left;
		TreeNode right;
		boolean visited;

		TreeNode(int x) {
			val = x;
		}

		TreeNode(int x, TreeNode left, TreeNode right) {
			val = x;
			this.left = left;
			this.right = right;
		}
	}

	public class ListNode {
		int val;
		ListNode next;

		ListNode(int x) {
			val = x;
		}
	}

	public class Bucket {
		int low;
		int high;

		public Bucket() {
			low = -1;
			high = -1;
		}
	}

	public class TreeLinkNode {
		int val;
		TreeLinkNode left = null;
		TreeLinkNode right = null;
		TreeLinkNode next = null;

		TreeLinkNode(int val) {
			this.val = val;
		}
	}

	public class UndirectedGraphNode {
		int label;
		boolean visited = false;
		List<UndirectedGraphNode> neighbors;

		UndirectedGraphNode(int x) {
			label = x;
			neighbors = new ArrayList<UndirectedGraphNode>();
		}
	}

	class RandomListNode {
		int label;
		RandomListNode next = null;
		RandomListNode random = null;

		RandomListNode(int x) {
			this.label = x;
		}
	}

	public class Position {
		int i;
		int j;

		Position(int i, int j) {
			this.i = i;
			this.j = j;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + getOuterType().hashCode();
			result = prime * result + i;
			result = prime * result + j;
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			Position other = (Position) obj;
			if (!getOuterType().equals(other.getOuterType()))
				return false;
			if (i != other.i)
				return false;
			if (j != other.j)
				return false;
			return true;
		}

		private LeetCode getOuterType() {
			return LeetCode.this;
		}
	}
}
