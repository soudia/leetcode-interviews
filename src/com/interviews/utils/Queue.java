package com.interviews.utils;

public class Queue<E> {

	Node first;
	Node last;

	public void enqueue(E value) {
		Node node = new Node(value);
		if (first == null) {
			last = node;
			first = last;
		} else {
			last.next = node;
			last = node;
		}
	}

	public E poll() {
		if (first == null)
			return null;
		return first.value;
	}

	public E dequeue() {
		if (first == null)
			return null;
		Node node = first;
		first = node.next;
		return node.value;
	}

	public boolean isEmpty() {
		return first == null;
	}

	class Node {
		Node next = null;
		E value;

		public Node(E value) {
			this.value = value;
		}
	}
}
