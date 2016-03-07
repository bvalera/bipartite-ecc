require(igraph)

# Calculates the separatedness of clusters (See http://arxiv.org/pdf/0707.1616v3.pdf).
# This can be used as a criterion instead of a fixed number of communities.
bipartiteModularity <- function(g, s) {
  
  nodesType0 <- which(V(g)$type)
  nodesType1 <- which(!V(g)$type)
  adj <- get.adjacency(g)
  
  deg <- rowSums(adj)
  
  sum(sapply(nodesType0, function(i) {
    sapply(nodesType1, function(j) {
      if (s[i] == s[j]) {
        adj[i,j] - ((deg[i] * deg[j]) / length(E(g)))
      } else {
        0
      } 
    })  
  })) / length(E(g))  
}

getEdgeClusteringCoefficients <- function(graph) {
  
  # For all edges
  sapply(E(graph), function(edgeId) {
      
    # Compute the edge clustering coefficient
    edge <- get.edge(graph, edgeId)
    n1 <- edge[1] # First incident node of the edge
    n2 <- edge[2] # second incident node of the edge
    #####
    
    # The number of common neighbors of the neighbors of n1 and n2 
    # correspond to the number of squares the edge is contained in.
    q=0 # Qijmn
    v1 <- neighbors(graph, n1) #m
    v2 <- neighbors(graph, n2) #n
    
    nv1 <- unlist(neighborhood(graph, 1, v1))
      
    # Number of observed squares
    q <- length(which(nv1 %in% v2))
    # Number of possible squares
    possibleQ <- length(v1) * length(v2)
    # ECC
    q / possibleQ
  })
}

# Compute bipartite clusters using the edge clustering coefficient approach.
# graph: Bipartite graph
# maximumNumberOfClusters: Maximum number of clusters when the algorithm should stop decomposing the graph. 
#(Bipartite modularity could be used instead)
getBipartiteClusters <- function(graph, maximumNumberOfClusters) {
  
  # Calculate connected components
  components <- clusters(graph)
  
  while (components$no < maximumNumberOfClusters) {
    print(components$no)
    # Calculate edge clustering coefficients
    ecc <- getEdgeClusteringCoefficients(graph)
    
    # Remove edge with lowest clustering coefficient
    edgeToRemove <- which(ecc == min(ecc))
    graph <- delete.edges(graph, edgeToRemove)
    
    components <- clusters(graph)
  }
  
  # Extract and return cluster membership vector. 
  #Each connected component of the graph correspond to one cluster.
  # Each element corresponds to a node in the graph and the value corresponds to the cluster of the node. 
  components$membership
}


# Start the script call: start()
start <- function() {
  
  g <- read.graph("PATH_TO_THE_GRAPH_FILE/marvel_reduced.net", "pajek")
  print(getBipartiteClusters(g, 2))
  print(V(g)$id)
}