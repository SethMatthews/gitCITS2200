//CITS2200 Project - 21973441 & 22244433


import java.util.Arrays;
import java.util.ArrayList;
import java.util.HashMap; // import the HashMap class
import java.util.LinkedList;
import java.util.Queue;
import java.util.Stack;

public class MyCITS2200Project implements CITS2200Project {

	public ArrayList<ArrayList<Integer>> graph;
	public ArrayList<ArrayList<Integer>> graphTranspose;
	public HashMap<Integer,String> mapA;
	public HashMap<String,Integer> mapB;
	public HashMap<Integer,Boolean> visited;
	//private ArrayList<Int> neighbours;

	public MyCITS2200Project () {
		this.mapA = new HashMap<Integer,String>();
		this.mapB = new HashMap<String,Integer>();
		this.visited = new HashMap<Integer,Boolean>();
		this.graph= new ArrayList<ArrayList<Integer>>();	//lists of each index's neighbours
		this.graphTranspose= new ArrayList<ArrayList<Integer>>();
		//this.neighbours = new ArrayList<Int>();
		}

	//Seth function//
	//adds vertice to both hashmaps
	//hashmaps link an interger with 
	public void addVert(String url){
		mapB.put(url,mapB.size());
		mapA.put(mapA.size(),url);

	}

	//Seth function//
	//returns true if vertice exists
	public boolean checkVert(String url){
		return mapB.containsKey(url);
	}

	//assign index for new url and updates both hashmaps
	//appends neighbour on the list located at the index respective to vertex number
	public void addEdge(String urlFrom, String urlTo) {
		if (!checkVert(urlFrom)) {	//if vertex is not in hashmap already
			addVert(urlFrom);
		}
		if (!checkVert(urlTo)) {	//if vertex is not in hashmap already
			addVert(urlTo);
		}

		int vertIndexTo = mapB.get(urlTo); //convert urlTo to an index through hashmap
		int vertIndexFrom = mapB.get(urlFrom); //convert urlFrom to an index through hashmap
		System.out.println(Integer.toString(vertIndexFrom) + ", " + Integer.toString(vertIndexTo));

		while(graph.size()<mapA.size()) {
			graph.add(null);	
		}
		if(graph.get(vertIndexFrom)==null) {
			ArrayList<Integer> newVert = new ArrayList<Integer>();
			newVert.add(vertIndexTo);
			graph.set(vertIndexFrom, newVert);
			System.out.println(graph.get(vertIndexFrom));
			
			//System.out.println("Add successful!");
		}
		else {
			graph.get(vertIndexFrom).add(vertIndexTo); //appends vertex(urlTo) in vertex(urlFrom) neighbour list
			System.out.println(graph.get(vertIndexFrom));
		}

		while(graphTranspose.size()<mapA.size()){
			graphTranspose.add(null);
		}
		if(graphTranspose.get(vertIndexTo)==null){
			ArrayList<Integer> newVertTranspose = new ArrayList<Integer>();
			newVertTranspose.add(vertIndexFrom);
			graphTranspose.set(vertIndexTo,newVertTranspose);
			System.out.println(graphTranspose.get(vertIndexTo));
		}
		else {
			graphTranspose.get(vertIndexTo).add(vertIndexFrom);
			System.out.println(graphTranspose.get(vertIndexTo));
		}
		//System.out.println(graph.get(vertIndexFrom));
		//System.out.println(graphTranspose.get(vertIndexTo));
		System.out.println("Graph: " + Integer.toString(graph.size()) + ", TransposeGraph: " + Integer.toString(graphTranspose.size()));
		System.out.println("");

	}


	public int getShortestPath(String urlFrom, String urlTo) {
		int rootNode = mapB.get(urlFrom);
		int finalNode = mapB.get(urlTo);
		HashMap<Integer,Integer> depth = new HashMap<Integer,Integer>();
		Queue<Integer> queue = new LinkedList<Integer>();

		queue.add(rootNode);
		depth.put(rootNode,0);

		while(!(queue.size()==0)) {
			int currentNode = queue.remove();
			if(currentNode == finalNode) {
				return depth.get(currentNode);
			}

			ArrayList<Integer> currentNodeChildren = graph.get(currentNode);
			if(currentNodeChildren!=null) {
				for(Integer i=0; i<currentNodeChildren.size(); i++) {
					int newVert = currentNodeChildren.get(i);
					if(!depth.containsKey(newVert)) {
						queue.add(newVert);
						depth.put(newVert,depth.get(currentNode)+1);
					}
				}
			}
		}
		return -1;

	}

	public void printMatrix(int[][] matrix) {
	    for (int row = 0; row < matrix.length; row++) {
	        for (int col = 0; col < matrix[row].length; col++) {
	            System.out.printf("%4d", matrix[row][col]);
	        }
	        System.out.println();
	    }
    }

    public void printArrayList(ArrayList<ArrayList<Integer>> a) {
    	for(int i=0; i<a.size(); i++) {
    		System.out.println(a.get(i));
    	}
    	System.out.println("");
    }

	public int[][] getAdjMatrix () {

		int[][] adjacency = new int[graph.size()][graph.size()];
		for(int i=0; i<graph.size(); i++) {
			ArrayList<Integer> adjacencyList = graph.get(i);

			//if adjacenyList.contins(j) set to 
			for(int j=0; j<graph.size(); j++) {
				if(i==j) {
					adjacency[i][j]=0;
				}
				else if(adjacencyList==null) {
					adjacency[i][j]=999;
				}
				else if(adjacencyList.contains(j)){
					adjacency[i][j]=1;
				}
				else {
					adjacency[i][j]=999;
				}
			}
		}
		return adjacency;
	}

	private static int[][] getFloydMatrix (int[][] a) {

		int size = a.length;
		for (int i=0; i<size; i++) {
			for (int j=0; j<size; j++) {
				for (int k=0; k<size; k++) {
					if(a[i][k] + a[k][j] < a[i][j]) {
						a[i][j] = Math.min(a[i][j], a[i][k] + a[k][j]);		
					}
				}
			}	
		}
		return a;
	}


	private String[] converter(ArrayList<String> s) {
		Object[] centers = s.toArray();
        String[] strCenters = Arrays.copyOf(centers, centers.length, String[].class);
        return strCenters;
	}

	public String[] getCenters() {
		
		int[][] adjacency = getAdjMatrix();
		int[][] floyd = getFloydMatrix(adjacency);
		printMatrix(floyd);
		int mostJumps = 999;
		ArrayList<String> potentialCenter = new ArrayList<String>();
		for (int i=0; i<floyd.length; i++) {
			for(int j=0; j<floyd.length; j++) {
				int potentialMJ = 0;
				if(floyd[i][j] == 999 || floyd[i][j] > mostJumps) {
					break;
				}
				else if(floyd[i][j]>potentialMJ) {
					potentialMJ = floyd[i][j];
				}
				if(j==floyd.length-1) {
					if(potentialMJ<mostJumps) {
						mostJumps = potentialMJ;
						potentialCenter.clear();
						potentialCenter.add(mapA.get(i));
					}
					else if(potentialMJ == mostJumps) {
						potentialCenter.add(mapA.get(i));
					}
				}
			}
		}
		String[] strCenters = converter(potentialCenter);
		return strCenters;
	}

	private Stack<Integer> depthFirstSearch(int vertex, Stack<Integer> vertStack, ArrayList<ArrayList<Integer>> searchGraph) {

		ArrayList<Integer> children = searchGraph.get(vertex);
		visited.put(vertex,true);
		if(children!=null) {
			for(int i=0; i<children.size(); i++) {
				if(!visited.containsKey(children.get(i))) {
					depthFirstSearch(children.get(i),vertStack,searchGraph);
				}
			}
		}
		vertStack.push(vertex);
		return vertStack;
    } 

    public String[][] sccConversion (ArrayList<Stack<Integer>> components) {
    	String[][] convertedComponents = new String[graph.size()][graph.size()];
    	for(int i=0; i<components.size(); i++) {
    		int j=0;
    		Stack<Integer> currStack = components.get(i);
    		while(!currStack.empty()) {
    			String url = mapA.get(currStack.pop());
    			convertedComponents[i][j] = url;
    			j++;
    		}
    	}
    	return convertedComponents;
    }

	public String[][] getStronglyConnectedComponents() {
		ArrayList<Stack<Integer>> alComponents = new ArrayList<Stack<Integer>>();
		Stack<Integer> vertStack = new Stack<Integer>();
		for(int i=0; i<graph.size(); i++) {
			if(visited.size()==graph.size()) {
				break;
			}
			if(!visited.containsKey(i)) {
				vertStack = depthFirstSearch(i,vertStack,graph);
			}
		}
		visited.clear();
		while(!vertStack.empty()) {
			Stack<Integer> tempStack = new Stack<Integer>();
			tempStack = depthFirstSearch(vertStack.pop(),tempStack,graphTranspose);
			if(tempStack.size()>1) {
				alComponents.add(tempStack);
			}
		}
		String[][] components = sccConversion(alComponents);
		return components;
	}

	public String[] getHamiltonianPath() {
		return null;

	}



}