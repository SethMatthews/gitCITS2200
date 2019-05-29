import java.io.*;
import java.util.*;

public class CITS2200ProjectTester {
	public static void loadGraph(MyCITS2200Project project, String path) {
		// The graph is in the following format:
		// Every pair of consecutive lines represent a directed edge.
		// The edge goes from the URL in the first line to the URL in the second line.
		try {
			BufferedReader reader = new BufferedReader(new FileReader(path));
			while (reader.ready()) {
				String from = reader.readLine();
				String to = reader.readLine();
				//System.out.println("Adding edge from " + from + " to " + to);
				project.addEdge(from, to);
			}
		} catch (Exception e) {
			System.out.println("There was a problem:");
			System.out.println(e.toString());
		}
	}

	public static void main(String[] args) {
		// Change this to be the path to the graph file.
		
		//String pathToGraphFile = "/Users/samuel/Documents/GitHub/gitCITS2200/example_graph.txt"; //SAM MAC PATH
		String pathToGraphFile = "C:/Users/User/Documents/GitHub/gitCITS2200/example_graph.txt"; //SAM PC PATH
		// Create an instance of your implementation.
		MyCITS2200Project proj = new MyCITS2200Project();
		// Load the graph into the project.
		loadGraph(proj, pathToGraphFile);
		for(int i=0; i<proj.graph.size(); i++) {
			for(int j=0; j<proj.graph.size(); j++) {
				System.out.println((Integer.toString(i)) + " to " + Integer.toString(j) + " shortestPath: " + Integer.toString(proj.getShortestPath(proj.mapA.get(i), proj.mapA.get(j))));
			}
		}
		String[] centers = proj.getCenters();
		for(int i=0; i<centers.length; i++) {
			System.out.println(centers[i]);
		}
		String[][] cc = proj.getStronglyConnectedComponents();
		for(int i=0; i<cc.length; i++) {
			for(int j=0; j<cc[i].length; j++) {
				System.out.print(cc[i][j]+" ");
			}
			System.out.println("");		
		}
		String[] path = proj.getHamiltonianPath();
		for(int i=0; i<path.length; i++) {
			System.out.println(path[i]);
		}

		// Write your own tests!
		//System.out.println(Integer.toString(proj.getShortestPath("/wiki/Flow_network", "/wiki/Ford%E2%80%93Fulkerson_algorithm")));
	}
}