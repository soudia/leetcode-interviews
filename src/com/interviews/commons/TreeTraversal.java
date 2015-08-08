package com.interviews.commons;

public class TreeTraversal {

	public static <U> void preOrder(Node<U> root) {
		if (root != null) {
			System.out.print(root.val);
			preOrder(root.left);
			preOrder(root.right);
		}
	}

	public static <U> void inOrder(Node<U> node) {
		if (node == null)
			return;
		inOrder(node.left);
		System.out.println(node.val);
		inOrder(node.right);
	}

	public static <U> void postOrder(Node<U> node) {
		if (node == null)
			return;
		preOrder(node.left);
		preOrder(node.right);
		System.out.print(node.val);
	}

	public static <E> void iterInOrder(Node<E> node) {
		java.util.Stack<Node<E>> stack = new java.util.Stack<Node<E>>();
		while (!stack.empty() || node != null) {
			if (node != null) {
				stack.push(node);
				node = node.left;
			} else {
				node = stack.pop();
				System.out.println(node.val);
				node = node.right;
			}
		}
	}

	public static <E> void iterPreOrder(Node<E> node) {
		java.util.Stack<Node<E>> stack = new java.util.Stack<Node<E>>();
		while (!stack.empty() || node != null) {
			if (node != null) {
				System.out.println(node.val);
				if (node.right != null) {
					stack.push(node.right);
				}
				node = node.left;
			} else {
				node = stack.pop();
			}
		}
	}

	public static <E> void iterPostOrder(Node<E> node) {
		java.util.Stack<Node<E>> stack = new java.util.Stack<Node<E>>();
		stack.push(node);
		Node<E> currNode, prevNode = null;
		while (!stack.isEmpty()) {
			currNode = stack.peek();
			if (prevNode == null || prevNode.left == currNode
					|| prevNode.right == currNode) {
				if (currNode.left != null) {
					stack.push(currNode.left);
				} else if (currNode.right != null) {
					stack.push(currNode.right);
				}
			} else if (currNode.left == prevNode) {
				stack.push(currNode.right);
			} else {
				System.out.println(node.val);
				stack.pop();
			}
			prevNode = currNode;
		}
	}
}
