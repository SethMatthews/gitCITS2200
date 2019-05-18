//CITS2200 Project - 21973441 & SETHNUMBER


import java.util.Arrays;
import java.util.ArrayList;
import java.util.HashMap; // import the HashMap class
import java.util.LinkedList;
import java.util.Queue;

public class MyCITS2200Project implements CITS2200Project {

	public ArrayList<ArrayList<Integer>> graph;
	public HashMap<Integer,String> mapA;
	public HashMap<String,Integer> mapB;
	public HashMap<Integer,Boolean> visited;
	//private ArrayList<Int> neighbours;

	public MyCITS2200Project () {
		this.mapA = new HashMap<Integer,String>();
		this.mapB = new HashMap<String,Integer>();
		this.visited = new HashMap<Integer,Boolean>();
		this.graph= new ArrayList<ArrayList<Integer>>();	//lists of each index's neighbours
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
		
		while(graph.size()-1<vertIndexFrom) {
			graph.add(null);
		}
		
		if(graph.get(vertIndexFrom)==null) {
			ArrayList<Integer> newVert = new ArrayList<Integer>();
			newVert.add(vertIndexTo);
			//System.out.println("Boutta add at index " + Integer.toString(vertIndexFrom) + " with graph size: " + Integer.toString(graph.size()));
			graph.set(vertIndexFrom, newVert);
			//System.out.println("Add successful!");
		}
		else {
			graph.get(vertIndexFrom).add(vertIndexTo); //appends vertex(urlTo) in vertex(urlFrom) neighbour list
		}
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

	private int[][] getAdjMatrix () {

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

	private Stack depthFirstSearch(int vertex, Stack vertStack) {

		ArrayList<Integer> children = graph.get(vertex);
		visited.put(vertex,true);
		if(children.size()!=null) {
			for(int i=0; i<children.size(); i++) {
				if(!visited.get(vertex)) {
					depthFirstSearch(children.get(i),vertStack);
				}
			}
			vertStack.push(i);
			return vertStack;
		}
		else {
			vertStack.push(i);
			return vertStack;
		}
    } 
  

    private ArrayList<ArrayList<Integer>> getTranspose() { 

    } 
  
    void fillOrder(int v, boolean visited[], Stack stack) { 

    } 

	public String[][] getStronglyConnectedComponents() {
		Stack<Integer> vertStack = new Stack<Integer>();
		for(int i=0; i<graph.(size); i++) {
			if(visited.size()==graph.size()) {
				break;
			}
			if(!visited.containsKey(i)) {
				vertStack = depthFirstSearch(i,vertStack);
			}
		}
		System.out.println(vertStack);
	}

	public String[] getHamiltonianPath() {
		return null;

	}



}