package com.interviews.commons;

public class Node<T> {

	public T val;
	public boolean visited;
	public Node<T> left, right, parent;

	public Node(T val) {
		left = null;
		right = null;
		this.val = val;
		this.parent = null;
	}

}
