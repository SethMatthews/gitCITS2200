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
	public int[][] bitCheck;

	public MyCITS2200Project () {
		this.mapA = new HashMap<Integer,String>();
		this.mapB = new HashMap<String,Integer>();
		this.visited = new HashMap<Integer,Boolean>();
		this.graph= new ArrayList<ArrayList<Integer>>();	//lists of each index's neighbours
		this.graphTranspose= new ArrayList<ArrayList<Integer>>();
		}

	/**
	 * adds a vertice to both hashmaps so it can be converted from a url to its corresponding node number and back
	 * 
	 * @param url the URL to be tested
	 */
	public void addVert(String url){
		mapB.put(url,mapB.size());
		mapA.put(mapA.size(),url);

	}

	/**
	 * checks if a vertice is present in mapB hashmap
	 * 
	 * @param url the URL to be tested
	 * @return true if the url is present else false
	 */
	public boolean checkVert(String url){
		return mapB.containsKey(url);
	}

	/**
	 * Adds an edge to the Wikipedia page graph. If the pages do not
	 * already exist in the graph, they will be added to the graph.
	 * 
	 * @param urlFrom the URL which has a link to urlTo.
	 * @param urlTo the URL which urlFrom has a link to.
	 */
	public void addEdge(String urlFrom, String urlTo) {
		if (!checkVert(urlFrom)) {	//if vertex is not in hashmap already
			addVert(urlFrom);
		}
		if (!checkVert(urlTo)) {	//if vertex is not in hashmap already
			addVert(urlTo);
		}

		int vertIndexTo = mapB.get(urlTo); //convert urlTo to an index through hashmap
		int vertIndexFrom = mapB.get(urlFrom); //convert urlFrom to an index through hashmap

		while(graph.size()<mapA.size()) {
			graph.add(null);	
		}
		if(graph.get(vertIndexFrom)==null) {
			ArrayList<Integer> newVert = new ArrayList<Integer>();
			newVert.add(vertIndexTo);
			graph.set(vertIndexFrom, newVert);
		}
		else {
			graph.get(vertIndexFrom).add(vertIndexTo); //appends vertex(urlTo) in vertex(urlFrom) neighbour list
		}

		while(graphTranspose.size()<mapA.size()){
			graphTranspose.add(null);
		}
		if(graphTranspose.get(vertIndexTo)==null){
			ArrayList<Integer> newVertTranspose = new ArrayList<Integer>();
			newVertTranspose.add(vertIndexFrom);
			graphTranspose.set(vertIndexTo,newVertTranspose);
		}
		else {
			graphTranspose.get(vertIndexTo).add(vertIndexFrom);
		}

	}

	/**
	 * Finds the shorest path in number of links between two pages.
	 * If there is no path, returns -1.
	 * 
	 * @param urlFrom the URL where the path should start.
	 * @param urlTo the URL where the path should end.
	 * @return the legnth of the shorest path in number of links followed.
	 */
	public int getShortestPath(String urlFrom, String urlTo) {
		int rootNode = mapB.get(urlFrom);
		int finalNode = mapB.get(urlTo);

		//return 0 if calculating the shortest path from a node to itself
		if(rootNode==finalNode) {
			return 0;
		}
		
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

	/**
	 * converts an ArrayList of Strings to a String array
	 *
	 * @param s the ArrayList of Strings corresponding to nodes that are centers in the graph
	 * @return all centers converted into strings and stored in a string array
	 */
	private String[] converter(ArrayList<String> s) {
		Object[] centers = s.toArray();
        String[] strCenters = Arrays.copyOf(centers, centers.length, String[].class);
        return strCenters;
	}

	/**
	 * performs a breadth first search from each node in the graph to find the shortest path from each node to each other node
	 *
	 * @return an array containing the shortest path from each node to each other node
	 */
	public int[][] BFS() {

		int[][] distanceList = new int[graph.size()][graph.size()];
		for(int i=0; i<graph.size(); i++) {
			
			//Queue to hold discovered vertexes that need to be explored in the depth first search
			Queue<Integer> queue = new LinkedList<Integer>();

			//adds root node to the queue and sets its depth to 0 in the HashMap
			queue.add(i);
			distanceList[i][i]=0;

			//Breadth first search implementation to test form matches to the goal node
			while(!(queue.size()==0)) {
				int currentNode = queue.remove();

				//tests if the current node has children and if so, adds unexplored nodes to the queue to be searched
				ArrayList<Integer> currentNodeChildren = graph.get(currentNode);
				if(currentNodeChildren!=null) {
					for(Integer j=0; j<currentNodeChildren.size(); j++) {
						int currentChild = currentNodeChildren.get(j);
						if(i==currentChild) {
							distanceList[i][currentChild]=0;
						}
						else if(distanceList[i][currentChild]==0) {
							queue.add(currentChild);
							distanceList[i][currentChild] = distanceList[i][currentNode]+1;
						}
					}
				}
			}
			for(int j=0; j<graph.size(); j++){
				if(distanceList[i][j]==0 && i!=j) {
					distanceList[i][j]=999;
				}
			}
		}
		return distanceList;
	}

	/**
	 * Finds all the centers of the page graph. The order of pages
	 * in the output does not matter. Any order is correct as long as
	 * all the centers are in the array, and no pages that aren't centers
	 * are in the array.
	 * 
	 * @return an array containing all the URLs that correspond to pages that are centers.
	 */
	public String[] getCenters() {
		
		//call to the distances function to return the matrix needed to calculate centers
		int[][] distances = BFS();

		//set initial most jumps variable to set the number a potential center has to 'beat'
		int mostJumps = 999;

		//for each row of the distances matrix (node of the graph), test whether it is a center
		ArrayList<String> potentialCenter = new ArrayList<String>();
		for (int i=0; i<distances.length; i++) {
			int potentialMJ = 0;
			for(int j=0; j<distances.length; j++) {
				//first, test if any values are 999 indicating a node is unreachable or if any node is larger than the recorded MJ. If true, break and process the next node
				if(distances[i][j] == 999 || distances[i][j] > mostJumps) {
					break;
				}
				//take note of the highest value in the node and store it in potentialMJ
				else if(distances[i][j]>potentialMJ) {
					potentialMJ = distances[i][j];
				}
				if(j==distances.length-1) {
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
		String[] strCenters;
		if(!potentialCenter.isEmpty()) {
			strCenters = converter(potentialCenter);
		}
		else{
			strCenters = new String[]{"No centers found"};
		}
		return strCenters;
	}

	/**
	 * implementation of a recursive depth first search that returns a stack for further processing
	 *
	 * @param vertex the current vertex being searched from
	 * @param vertStack current stack being added to by the depth first search
	 * @param searchGraph the graph that the DFS is performed upon
	 * @return a stack containing all objects found in the search
	 */
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

    /**
	 * converts the stack returned by the secondary DFS on the Transposed graph into the required 2D string array
	 *
	 * @param components the ArrayLists of Stacks of strongly connected components
	 * @return a jagged array where each row is an array of strongly connected components
	 */
    public String[][] sccConversion (ArrayList<Stack<Integer>> components) {
    	String[][] convertedComponents = new String[components.size()][];
    	for(int i=0; i<components.size(); i++) {
    		Stack<Integer> currStack = components.get(i);
    		convertedComponents[i] = new String[currStack.size()];
    		for(int j=0;!currStack.empty(); j++) {
    			convertedComponents[i][j] = mapA.get(currStack.pop());
    		}
    	}
    	return convertedComponents;
    }

    /**
	 * Finds all the strongly connected components of the page graph.
	 * Every strongly connected component can be represented as an array 
	 * containing the page URLs in the component. The return value is thus an array
	 * of strongly connected components. The order of elements in these arrays
	 * does not matter. Any output that contains all the strongly connected
	 * components is considered correct.
	 * 
	 * @return an array containing every strongly connected component.
	 */
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
				alComponents.add(tempStack);
			}
		}
		//convert components to the required form and return them
		String[][] components = sccConversion(alComponents);
		return components;
	} 

	/**
	 * changes the n'th node of the state to 0 if it was 1 or 1 if it was 0, counting from the right hand side of the state
	 *
	 * @param n the last node visited in the Hamiltonian path
	 * @param binNumber the current state of nodes visited
	 * @return the new state created when the n'th node is flipped
	 */
	public int flipNth (int binNumber, int n) {
		return (binNumber ^ (1 << n));
	}

	/**
	 * sets the n'th node of the state to 1, counting from the right hand side of the state
	 *
	 * @param n the last node visited in the Hamiltonian path
	 * @param binNumber the current state of nodes visited
	 * @return the new state created when the n'th node is set to 1
	 */
	public int setNth1 (int binNumber, int n) {
		return (binNumber |= (1 << n));
	}

	/**
	 * sets the n'th node of the state to 0, counting from the right hand side of the state
	 *
	 * @param n the last node visited in the Hamiltonian path
	 * @param binNumber the current state of nodes visited
	 * @return the new state created when the n'th node is set to 0
	 */
	public int setNth0 (int binNumber, int n) {
		return (binNumber &= ~(1 << n));
	}

	/**
	 * checks whether node n has been visited in the current state
	 *
	 * @param n the last node visited in the Hamiltonian path
	 * @param binNumber the current state of nodes visited
	 * @return True if n is not in binNumber of false if it is
	 */
	public boolean notIn (int n, int binNumber) {
		return (((1<<n) & binNumber)==0);
	}

	/**
	 * Recursive function to find the next node of a potential Hamiltonian path
	 *
	 * @param n the last node visited in the Hamiltonian path
	 * @param binNumber the current state of nodes visited
	 * @return the last node visited in the first index and the current state in the second index or -1 if the state lead to a dead end.
	 */
	public int[] getNext (int n, int binNumber) {

		//create a returnArray with current parent as n and current state as binNumber
		int[] returnArray = new int[]{n,binNumber};
		//check current state is failure or endstate
		if(binNumber==-1 || binNumber==(1<<graph.size())-1){
			return returnArray;
		}
		//attempt to branch one node further in current path
		ArrayList<Integer> children = graphTranspose.get(n);
		if(children!=null) {
			for(int child=0; child<children.size(); child++) {
				int i = children.get(child);
				//if node i already in path, skip processing
				if(!notIn(i,binNumber)) {
					continue;
				}
				//if this state and parent combination has been calculated before, break out of loop 
				if(bitCheck[i][setNth1(binNumber,i)]!=0) {
					returnArray[1]=-1;
					return returnArray;
				}
				//if node i not in path, check if node i is adjacent to node n
				//set state to state+parent
				binNumber=setNth1(binNumber,i);
				//set returnArray state to state+parent
				returnArray[1] = binNumber;
				//record parent and state+parent in bitcheck array with value equal to child
				bitCheck[i][binNumber]=n;
				//System.out.println("boutta recurse");
				//create an array to hold the result of recursive call
				int[] tempSolution = getNext(i,binNumber);
				//check if recursive call state is equal to end state and if so, return it
				if(tempSolution[1]==(1<<graph.size())-1) {
					return tempSolution;
				}
				//if recursive calls don't lead to the endstate, revert changes to binNumber
				binNumber = setNth0(binNumber,i);
				returnArray[1]= binNumber;
			}
		}
		return returnArray;
	}

	/**
	 * Finds the returned Hamiltonian path by backtracking from the last node visited
	 *
	 * @param parent the last node visited in the Hamiltonian path
	 * @param state the current state of nodes visited
	 * @return a Hamiltonian path of the page graph in string array form.
	 */
	public String[] pathFetcher(int parent, int state) {
		
		/*
		* Use the state and parent values to crawl back through the path the
		* getNext function took and add each node's url form to the final path
		*/
		String[] path = new String[graph.size()];
		for(int i=0; i<graph.size(); i++) {
			path[i] = mapA.get(parent);
			int tempParent = bitCheck[parent][state];
			state = setNth0(state,parent);
			parent = tempParent;
		}
		return path;
	}

	/**
	 * Finds a Hamiltonian path in the page graph. There may be many
	 * possible Hamiltonian paths. Any of these paths is a correct output.
	 * This method should never be called on a graph with more than 20
	 * vertices. If there is no Hamiltonian path, this method will
	 * return an empty array. The output array should contain the URLs of pages
	 * in a Hamiltonian path. The order matters, as the elements of the
	 * array represent this path in sequence. So the element [0] is the start
	 * of the path, and [1] is the next page, and so on.
	 * 
	 * @return a Hamiltonian path of the page graph.
	 */
	public String[] getHamiltonianPath() {
		
		//set a state integer that is x 1's where x is the number of nodes in the graph
		final int END_STATE = (1<<graph.size())-1;
		bitCheck = new int[graph.size()][1<<graph.size()];

		int[] result;
		String[] finalString = null;
		//make a call to the get next function for each node
		for(int i=0; i<graph.size(); i++) {
			result = getNext(i,1<<i);
			//if a result of one of the calls matches the end state, get the path and break out of the loop
			if(result[1]==END_STATE) {
				finalString = pathFetcher(result[0], result[1]);
				break;
			}
		}
		//if no path has been found, set the string array to reflect that
		if(finalString==null) {
			String[] s = new String[]{};
			finalString = s;
		}
		return finalString;
	}



}