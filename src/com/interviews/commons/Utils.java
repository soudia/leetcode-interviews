package com.interviews.commons;

public class Utils<T> {

	public Node<T> successor(Node<T> node) {
		if (node.right != null) {
			return minimum(node.right);
		}
		Node<T> parent = node.parent;
		while (parent != null && parent.right == node) {
			node = parent;
			parent = parent.parent;
		}
		return null;
	}

	public Node<T> minimum(Node<T> right) {
		Node<T> node = right;
		while (node.left != null) {
			node = node.left;
		}
		return node;
	}
}
