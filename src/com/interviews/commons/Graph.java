package com.interviews.commons;

import java.util.ArrayList;
import java.util.List;

public class Graph<E> {

	public final List<Edge> edges;
	public final List<Vertex> vertices;

	public Graph() {
		this.edges = new ArrayList<Graph<E>.Edge>();
		this.vertices = new ArrayList<Graph<E>.Vertex>();
	}

	public Graph(final List<Edge> edges, final List<Vertex> vertices) {
		this.edges = edges;
		this.vertices = vertices;
	}

	private enum Direction {
		AB("A->B"), BA("B->"), ABBA("A->B->A");
		Direction(String description) {
			this.setDescription(description);
		}

		public void setDescription(String description) {
		}
	}

	public class Edge {
		public Vertex start = null, end = null;
		public Direction direction;
		public int distance = 0;

		Edge(Vertex start, Vertex end) {
			new Edge(start, end, Direction.AB);
		}

		Edge(Vertex start, Vertex end, Direction direction) {
			this.start = start;
			this.end = end;
			this.direction = direction;
		}

		public int hashCode() {
			int hStart = (start == null) ? 0 : start.hashCode();
			int hEnd = (end == null) ? 0 : end.hashCode();
			int hDir = direction.hashCode();
			return ((hStart + hEnd) * hEnd) + hDir;
		}

		public boolean equals(Object obj) {
			if (obj == null)
				return false;
			if (!(obj instanceof Graph.Edge))
				return false;
			@SuppressWarnings("unchecked")
			Edge edge = (Edge) obj;
			return this.start.equals(edge.start) && this.end.equals(edge.end)
					&& this.direction.equals(edge.direction);
		}
	}

	public class Vertex {
		E key;
		Vertex parent = null;
		String color = "WHITE";
		int time = 0;

		Vertex(E key) {
			this.key = key;
		}
	}
}
