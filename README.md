# spingraphSDP
SDP solver for some relaxed problems in very large sparse graphs  
Feel free to use these lines as you wish.  
NOTE THAT THIS IS A WORK IN PROGRESS: not guaranteed to work.

## Goemans-Williamson algorithm using spins

### Info:
This is an implementation of the Goemans-Williamson algorithm using spins.  
It is inspired from this: http://web.stanford.edu/~montanar/SDPgraph/home.html

### To compile:
"gcc spinmaxcut.c -O9 -o spinmaxcut -lm".

### To execute:
"./spinmaxcut edgelist.txt k t embedding.txt".
- "edgelist.txt" should contain the undirected unweighted graph: one edge on each line (two unsigned long (nodes' ID)) separated by a space.
- k is the dimension of the embedding (e.g. k = 20)
- t is the number of iterations (e.g. t = 1000)
- embedding.txt: will contain an embeding of the graph on a sphere in dimension k.
- lab.txt: will contain the assignment for each node -1 and 1 (random hyperplane cut)
- it will print the size of the cut in the terminal

### Scalability:

Sparse graphs with 10M edges in few minutes on a commodity machine

### Initial contributors:
Maximilien Danisch  
May 2018  
http://bit.ly/danisch  
maximilien.danisch@gmail.com
