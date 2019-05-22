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
	public boolean[][] bitCheck;
	public int[][] adjacency;
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

	//gets shortest path from urlFrom to urlTo
	public int getShortestPath(String urlFrom, String urlTo) {
		int rootNode = mapB.get(urlFrom);
		int finalNode = mapB.get(urlTo);
		
		//HashMap used to store depth of each node from root node
		HashMap<Integer,Integer> depth = new HashMap<Integer,Integer>();
		//Queue to hold discovered vertexes that need to be explored in the depth first search
		Queue<Integer> queue = new LinkedList<Integer>();

		//adds root node to the queue and sets its depth to 0 in the HashMap
		queue.add(rootNode);
		depth.put(rootNode,0);

		//Breadth first search implementation to test form matches to the goal node
		while(!(queue.size()==0)) {
			int currentNode = queue.remove();
			if(currentNode == finalNode) {
				//returns target node depth when node is found
				return depth.get(currentNode);
			}

			//tests if the current node has children and if so, adds unexplored nodes to the queue to be searched
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
		//Returns -1 if no matching node to the target node can be found
		return -1;

	}


	//Function for printing matrixes (NOT FOR FINAL PROGRAM)
	public void printMatrix(int[][] matrix) {
	    for (int row = 0; row < matrix.length; row++) {
	        for (int col = 0; col < matrix[row].length; col++) {
	            System.out.printf("%4d", matrix[row][col]);
	        }
	        System.out.println();
	    }
    }

    //Function for printing ArrayLists (NOT FOR FINAL PROGRAM)
    public void printArrayList(ArrayList<ArrayList<Integer>> a) {
    	for(int i=0; i<a.size(); i++) {
    		System.out.println(a.get(i));
    	}
    	System.out.println("");
    }

    //Takes graph ArrayList and converts it into an adjacency list for use in the Floyd algorithm
	public int[][] getAdjMatrix () {

		int[][] adjacency = new int[graph.size()][graph.size()];
		
		//Loop through all nodes and their children to test for adjacency between nodes
		//0 indicates adjacency[x][x] (a vertexes adjacency to itself)
		//1 indicates that x is adjacent to y for adjacency[x][y]
		//999 indicates y is non adjacent to y for adjacency[x][y]
		for(int i=0; i<graph.size(); i++) {
			ArrayList<Integer> adjacencyList = graph.get(i);
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

	//Finds Floyd matrix given an adjacency matrix
	//Floyd matrix is the shortest path from any node to another
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

	//Simple function to convert an ArrayList of Strings to a String array
	private String[] converter(ArrayList<String> s) {
		Object[] centers = s.toArray();
        String[] strCenters = Arrays.copyOf(centers, centers.length, String[].class);
        return strCenters;
	}

	//Function to find the vertexes that are centers of the graph
	public String[] getCenters() {
		
		//call to the adjacency matrix function to get the matrix needed for the floyd function
		adjacency = getAdjMatrix();
		//call to the floyd function to return the matrix needed to calculate centers
		int[][] floyd = getFloydMatrix(adjacency);

		//set initial most jumps variable to set the number a potential center has to 'beat'
		int mostJumps = 999;

		//for each row of the floyd matrix (node of the graph), test whether it is a center
		ArrayList<String> potentialCenter = new ArrayList<String>();
		for (int i=0; i<floyd.length; i++) {
			int potentialMJ = 0;
			for(int j=0; j<floyd.length; j++) {
				//first, test if any values are 999 indicating a node is unreachable or if any node is larger than the recorded MJ. If true, break and process the next node
				if(floyd[i][j] == 999 || floyd[i][j] > mostJumps) {
					break;
				}
				//take note of the highest value in the node and store it in potentialMJ
				else if(floyd[i][j]>potentialMJ) {
					potentialMJ = floyd[i][j];
				}
				if(j==floyd.length-1) {
					//if PotentialMJ is less than mostJumps, clear the pottentialCenter ArrayList and add the current node as a better potential center
					if(potentialMJ<mostJumps) {
						mostJumps = potentialMJ;
						potentialCenter.clear();
						potentialCenter.add(mapA.get(i));
					}
					//if potentialMJ is equal to mostJumps, add the current node to the list of potential centers
					else if(potentialMJ == mostJumps) {
						potentialCenter.add(mapA.get(i));
					}
				}
			}
		}
		//convert and return the final list of centers
		String[] strCenters = converter(potentialCenter);
		return strCenters;
	}

	//implementation of a recursive depth first search that returns a stack for further processing
	private Stack<Integer> depthFirstSearch(int vertex, Stack<Integer> vertStack, ArrayList<ArrayList<Integer>> searchGraph) {

		//find all children of the current node
		ArrayList<Integer> children = searchGraph.get(vertex);
		//mark the current node as visited
		visited.put(vertex,true);
		//while the current node is not a leaf, recursively call another DFS on its children
		if(children!=null) {
			for(int i=0; i<children.size(); i++) {
				if(!visited.containsKey(children.get(i))) {
					depthFirstSearch(children.get(i),vertStack,searchGraph);
				}
			}
		}
		//when a leaf is found, push it onto the stack and return it
		vertStack.push(vertex);
		return vertStack;
    } 

    //converts the stack returned by the secondary DFS on the Transposed graph into the required 2D string array
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

    //finds strongly connected components in the graph
	public String[][] getStronglyConnectedComponents() {
		
		//create an arraylist that will eventually contain stacks of strongly connected components
		ArrayList<Stack<Integer>> alComponents = new ArrayList<Stack<Integer>>();
		//stack to hold results from DFS
		Stack<Integer> vertStack = new Stack<Integer>();
		//Perform DFS on graph nodes until all nodes are in the vertStack
		for(int i=0; i<graph.size(); i++) {
			//break if all nodes are in the stack
			if(visited.size()==graph.size()) {
				break;
			}
			//DFS from the current node if it's not in the stack
			if(!visited.containsKey(i)) {
				vertStack = depthFirstSearch(i,vertStack,graph);
			}
		}
		//clear the visited list in preparation for the second DFS with the Transposed graph
		visited.clear();
		//while there are still nodes to be processed, perform DFS on each node popped from the stack
		while(!vertStack.empty()) {
			int currNode = vertStack.pop();
			//check if currentNode has been visited and don't process if it has
			if(!visited.containsKey(currNode)) {
				Stack<Integer> tempStack = new Stack<Integer>();
				tempStack = depthFirstSearch(currNode,tempStack,graphTranspose);
				if(tempStack.size()>1) {
					alComponents.add(tempStack);
				}
			}
		}
		//convert components to the required form and return them
		String[][] components = sccConversion(alComponents);
		return components;
	} 

	public int flipNth (int binNumber, int n) {
		return (binNumber ^ (1 << n));
	}

	public int setNth1 (int binNumber, int n) {
		return (binNumber |= (1 << n));
	}

	public int setNth0 (int binNumber, int n) {
		return (binNumber &= ~(1 << n));
	}

	public boolean notIn (int n, int binNumber) {
		return (((1<<n) & binNumber)==0);
	}

	public int getNext (int n, int binNumber) {

		System.out.println("");
		System.out.println("Starting search with current node: " + Integer.toString(n));
		System.out.println("and binary string: " + Integer.toBinaryString(binNumber));
		System.out.println("");

		if(binNumber==-1 || binNumber==(1<<graph.size())-1){
			return binNumber;
		}
		for(int i=0; i<graph.size(); i++) {
			if(!notIn(i,binNumber)) {
				System.out.println("continued as " + Integer.toString(i) + " was already in the string");
				continue;
			}
			if(bitCheck[i][setNth1(binNumber,i)]) {
				System.out.println("very sad");
				return(-1);
			}
			if(adjacency[i][n]==1) {
				binNumber=setNth1(binNumber,i);
				bitCheck[i][binNumber]=true;
				System.out.println("boutta recurse");
				int tempSolution = getNext(i,binNumber);
				System.out.println(Integer.toBinaryString(tempSolution));
				if(tempSolution==(1<<graph.size())-1) {
					return tempSolution;
				}
				binNumber=setNth0(binNumber,i);
			}
		}
		System.out.println("");
		System.out.println("recursion complete");
		System.out.println(Integer.toBinaryString(binNumber));
		System.out.println("");
		return binNumber;
	}

	public String[] getHamiltonianPath() {
		
		final int END_STATE = (1<<graph.size())-1;
		System.out.println(Integer.toBinaryString(END_STATE));
		printMatrix(adjacency = getAdjMatrix());
		bitCheck = new boolean[graph.size()][1<<graph.size()];

		int result;
		for(int i=0; i<graph.size(); i++) {
			result = getNext(i,1<<i);
			if(result==END_STATE) {
				System.out.println("success");
			}
		}
		return null;

	}



}