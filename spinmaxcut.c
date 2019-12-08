/*
Maximilien Danisch
May 2018
http://bit.ly/danisch
maximilien.danisch@gmail.com

Info:
Feel free to use these lines as you wish.
This is an implementation of the Goemans-Williamson algorithm using spins.
It is inspired from this: http://web.stanford.edu/~montanar/SDPgraph/home.html


To compile:
"gcc spinmaxcut.c -O9 -o spinmaxcut -lm".

To execute:
- "./spinmaxcut edgelist.txt k t embedding.txt".
- "edgelist.txt" should contain the undirected unweighted graph: one edge on each line (two unsigned long (nodes' ID)) separated by a space.
- k is the dimension of the embedding
- t is the number of iterations
- embedding.txt: will contain an embeding of the graph on a sphere in dimension k.
- lab.txt: will contain the assignment for each node -1 and 1 (random hyperplane cut)
- it will print the size of the cut in the terminal
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <strings.h>//to use "bzero"
#include <math.h>//to use "sqrt", "log" and "cos"
#include <time.h>//to estimate the runing time

#define NLINKS 100000000 //maximum number of edges for memory allocation, will increase if needed
#define PI 3.14159265358979323846

typedef struct {
	unsigned long s;
	unsigned long t;
} edge;

//edge list structure:
typedef struct {
	unsigned long n;//number of nodes
	unsigned long e;//number of edges
	edge *edges;//list of edges
	unsigned long *cd;//cumulative degree cd[0]=0 length=n+1
	unsigned long *adj;//concatenated lists of neighbors of all nodes
} adjlist;

//compute the maximum of three unsigned long
inline unsigned long max3(unsigned long a,unsigned long b,unsigned long c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

//reading the edgelist from file
adjlist* readedgelist(char* input){
	unsigned long e1=NLINKS;
	adjlist *g=malloc(sizeof(adjlist));
	FILE *file;

	g->n=0;
	g->e=0;
	file=fopen(input,"r");
	g->edges=malloc(e1*sizeof(edge));
	while (fscanf(file,"%lu %lu", &(g->edges[g->e].s), &(g->edges[g->e].t))==2) {
		g->n=max3(g->n,g->edges[g->e].s,g->edges[g->e].t);
		if (++(g->e)==e1) {
			e1+=NLINKS;
			g->edges=realloc(g->edges,e1*sizeof(edge));
		}
	}
	fclose(file);
	g->n++;

	g->edges=realloc(g->edges,g->e*sizeof(edge));

	return g;
}

//building the adjacency matrix
void mkadjlist(adjlist* g){
	unsigned long i,u,v;
	unsigned long *d=calloc(g->n,sizeof(unsigned long));

	for (i=0;i<g->e;i++) {
		d[g->edges[i].s]++;
		d[g->edges[i].t]++;
	}

	g->cd=malloc((g->n+1)*sizeof(unsigned long));
	g->cd[0]=0;
	for (i=1;i<g->n+1;i++) {
		g->cd[i]=g->cd[i-1]+d[i-1];
		d[i-1]=0;
	}

	g->adj=malloc(2*g->e*sizeof(unsigned long));

	for (i=0;i<g->e;i++) {
		u=g->edges[i].s;
		v=g->edges[i].t;
		g->adj[ g->cd[u] + d[u]++ ]=v;
		g->adj[ g->cd[v] + d[v]++ ]=u;
	}

	free(d);
	//free(g->edges);
}

//freeing memory
void free_adjlist(adjlist *g){
	free(g->edges);
	free(g->cd);
	free(g->adj);
	free(g);
}

//https://en.wikipedia.org/wiki/Fisher-Yates_shuffle
//not used yet. Maybe nodes should be shuffled between each iterations.
void shuff(unsigned long n, unsigned long *tab){
	unsigned long i,j,tmp;
	for (i=n;i>1;i--){
     	j=rand()%i;
		tmp=tab[i];
		tab[i]=tab[j];
		tab[j]=tmp;
	}
}

//generating number from standard normal distribution.
//https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
//will be used to generate points on a hypersphere
double gaussian() {
	double u1, u2;
	u1 = (rand()+1.) / (RAND_MAX+1.);
	u2 = (rand()+1.) / (RAND_MAX+1.);
	return sqrt(-2.*log(u1))*cos(2.*PI*u2);
}

//generating n points on a hypersphere of dimension k uniformly at random.
double *init_embedding(unsigned long n, unsigned k){
	unsigned long i;
	unsigned j;
	double s;
	double *emb=malloc(k*n*sizeof(double));//emb[k*i...k*(i+1)]=vector associated to node i
	for (i=0;i<n;i++){
		s=0.;
		for (j=0;j<k;j++){
			emb[i*k+j]=gaussian();
			s+=emb[i*k+j]*emb[i*k+j];
		}
		s=sqrt(s);
		for (j=0;j<k;j++){
			emb[i*k+j]/=s;
		}
	}
	return emb;
}


double *spinmaxcut(adjlist *g,unsigned k, unsigned t) {
	unsigned long n=g->n,u,v,i,tmp1,tmp2;
	unsigned t2,j;
	double s;
	double *emb=init_embedding(n,k);
	//double cut;//to print the values of the SDP after each iteration

	for (t2=0;t2<t;t2++) {
		//shuff(g->n,nodes);
		//cut=0;
		for (u=0;u<n;u++) {
			tmp1=k*u;
			bzero(emb+tmp1,k*sizeof(double));//the vector associated to node u is set to the null vector
			for (i=g->cd[u];i<g->cd[u+1];i++){
				v=g->adj[i];
				tmp2=k*v;
				for (j=0;j<k;j++){
					emb[tmp1+j]-=emb[tmp2+j];
				}
			}
			s=0;
			for (j=0;j<k;j++){
				s+=emb[tmp1+j]*emb[tmp1+j];
			}
			if (s>0){/////////////////////
				s=sqrt(s);
				//cut+=s;
				for (j=0;j<k;j++){
					emb[tmp1+j]/=s;
				}
			}
			else{
				for (j=0;j<k;j++){
					emb[tmp1+j]=gaussian();
					s+=emb[tmp1+j]*emb[tmp1+j];
				}
				s=sqrt(s);
				for (j=0;j<k;j++){
					emb[tmp1+j]/=s;
				}
			}
		}
		//cut=(2.*g->e+cut)/4.;
		//printf("%le\n",cut);
	}

	return emb;
}

char* hyperplanecut(adjlist* g,double* emb,unsigned k, unsigned long *cut, double *cutsdp){
	double* vect=malloc(k*sizeof(double));
	unsigned j,l;
	unsigned long i,u,tmp;
	double s;
	char *lab1=malloc(g->n*sizeof(char)),*lab2=malloc(g->n*sizeof(char)),*lab3;//lab[i] = label of node i
	unsigned long cut2;

	*cut=0;
	for (l=0;l<10;l++){//several random cut
		//random direction
		for (j=0;j<k;j++){
			vect[j]=gaussian();
		}

		//assigning label
		for (u=0;u<g->n;u++) {
			s=0;
			tmp=k*u;
			for (j=0;j<k;j++){
				s+=emb[tmp+j]*vect[j];
			}
			if (s<0)
				lab1[u]=-1;
			else
				lab1[u]=1;
		}

		//computing size of the cut
		cut2=0;
		for (i=0;i<g->e;i++) {
			if (lab1[g->edges[i].s] != lab1[g->edges[i].t]){
				cut2++;
			}
		}
		//printf("cut scal = %le\n",(g->e-cut3)/2);
		if (cut2>*cut){
			*cut=cut2;
			lab3=lab1;
			lab1=lab2;
			lab2=lab3;
		}
	}

	*cutsdp=0.;
	for (i=0;i<g->e;i++) {
		for (j=0;j<k;j++) {
			(*cutsdp)+=emb[g->edges[i].s*k+j]*emb[g->edges[i].t*k+j];
		}
	}
	*cutsdp=(g->e-(*cutsdp))/2.;

	return lab2;
}

int main(int argc,char** argv){
	adjlist* g;
	double *emb,cutsdp;
	unsigned long u,cut;
	unsigned t,k,i;
	char* lab;
	FILE* file;

	srand(time(NULL));

	time_t t1,t2;

	t1=time(NULL);

	printf("Number of dimensions of the embedding: %s\n",argv[2]);
	k=atoi(argv[2]);

	printf("Number of iterations: %s\n",argv[3]);
	t=atoi(argv[3]);

	printf("Reading edgelist from file %s\n",argv[1]);
	g=readedgelist(argv[1]);

	printf("Number of nodes: %lu\n",g->n);
	printf("Number of edges: %lu\n",g->e);

	printf("Building the adjacency list\n");
	mkadjlist(g);

	printf("Computing the Goemans-Williamson embedding using %u-dimensional spins\n",k);

	emb=spinmaxcut(g,k,t);

	printf("Printing resulting embedding in file %s\n",argv[4]);
	file=fopen(argv[4],"w");
	for (u=0;u<g->n;u++){
		for (i=0;i<k-1;i++){
			fprintf(file,"%lf ",emb[u*k+i]);
		}
		fprintf(file,"%lf\n",emb[u*k+k-1]);
	}
	fclose(file);

	printf("Random cut of the hypersphere\n");

	lab=hyperplanecut(g,emb,k,&cut,&cutsdp);

	printf("Printing labels in file %s\n",argv[5]);
	file=fopen(argv[5],"w");
	for (u=0;u<g->n;u++){
		fprintf(file,"%d\n",lab[u]);
	}
	fclose(file);

	printf("Objective of the sdp relaxation = %le\n",cutsdp);
	printf("Size of the cut = %lu\n",cut);

	t2=time(NULL);

	printf("- Overall time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));

	return 0;
}

