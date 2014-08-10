#ifndef _OPF_H_
#define _OPF_H_


#include "cimage.h"
#include "scene.h"
#include "adjacency.h"
#include "adjacency3.h"
#include "set.h"
#include "heap.h"
#include "realheap.h"
#include <time.h>

#define MAXDENS 1000   // Maximum value for pdf computation
#define MAXARCW 1000   // Maximum arc weight for supervised
		       // classification using complete graphs
// Status of the graph nodes for classification

#define PROTOTYPE   1

//variables used for precomputed distances
extern bool PrecomputedDistance; 
extern float  **Distance;

typedef struct _felem {
  float *feat;
} FElem;

typedef struct _features {
  FElem *elem;
  int  nfeats;
  int  nelems;
  int  nrows,ncols;
  int  Imax;
} Features;

typedef struct _features3 {
  FElem *elem;
  int  nfeats;
  int  nelems;
  int  xsize,ysize,zsize;
  int  Imax;
} Features3;

typedef struct _cnode {
  int    dens;    // density
  float *feat;    // feature vector
  int   *adj;     // list of adjacent nodes
  float *dist;    // distance from adjacent nodes  
  int    position; // position in the image/scene
} CNode;

typedef struct _snode {
  int   dens;    // probability density value
  int   pathval; // path value
  int   label;   // node label
  int   pred;    // predecessor node
  int truelabel; // true label if it is known
  int position;  // index in the image/scene
  float *feat;    // feature vector
  Set  *adj;     // list of adjacent nodes
  char  status;  // 0 - nothing, 1 - prototype 
} SNode;

typedef struct _subgraph {
  SNode *node;   // nodes of the image/scene subgraph
  int   nnodes;  // number of nodes
  int   nfeats;  // number of features
  int   bestk;   // number of adjacent nodes
  int   nlabels; // number of clusters
  float df;      // radius in the feature space for density computation
  float dm;      // radius for markov feature computation
  float di;      // radius in the image space for density computation
  int   mindens; // minimum density value
  int   maxdens; // maximum density value
  float K1,K2;  // Constants for PDF computation
  int   Imax;    // Maximum image value;
} Subgraph;

// ---- Component tree for pdf filtering ------------

typedef struct _sgctnode {
  int  level;   /* gray level */
  int  comp;    /* representative pixel of this node */
  int  dad;     /* dad node in the maxtree */
  int *son;     /* son nodes in the maxtree */
  int  numsons; /* number of sons in the maxtree */
  int  size;    /* number of pixels of the node */ 
} SgCTNode;

typedef struct _sgctree {
  SgCTNode *node;     /* nodes of the mtree */
  int      *cmap;     /* component map */
  int       root;     /* root node of the mtree */
  int       numnodes; /* number of nodes of the maxtree */
} SgCTree;


/*----------- Constructor and destructor ------------------------*/

Subgraph *CreateSubgraph(int nnodes); // Allocate nodes without
				      // feature space and adjacency list

void DestroySubgraph(Subgraph **sg); // Deallocate memory for subgraph


/*----------- Create subgraph (training set) with image/scene sampling
  ----------- and indicate the image/scene position of each node ----- */

// Get samples within a binary mask (optional) by skipping dx pixels
// along x and dy pixels along y. if mask is null, then the image
// domain is used. 

Subgraph *UnifSampl(Image *img, Image *mask, int dx, int dy);
Subgraph *UnifSampl3(Scene *scn, Scene *mask, int dx, int dy, int dz);

// Compute samples for image compression
Subgraph *ChessSampl(Image *img, int xinit);
Subgraph *VertSampl(Image *img);

// Get samples within a binary mask (optional) by choosing half below
// and half above a given threshold. 

Subgraph *ThresSampl3(Scene *scn, Scene *mask, int thres, int nnodes);
Subgraph *ThresNSampl3(Scene *scn, Scene *mask, int *thres, int N, int nnodes);

/*----------- Feature computation based on Markov properties -----------*/

// The intensity/color of a node and of its neighboors within a di
// radius in the image/scene are used to form its feature vector. 

void MarkovFeat3(Subgraph *sg, Scene *scn, Scene *mask, float dm);

// Compute features for the image/cimage 

Features *MSImageFeats(Image *img, int nscales);
Features *MSCImageFeats(CImage *cimg, int nscales);
Features *MarkovImageFeats(Image *img, float dm);
Features *KMarkovImageFeats(Image *img); // for image compression
Features *MarkovCImageFeats(CImage *cimg, float dm);
Features3 *MarkovSceneFeats(Scene *scn, Scene *mask, float dm);
void      DestroyFeatures(Features **f);
void      DestroyFeatures3(Features3 **f);

/*--------- Training functions compute an optimum path forest for
  --------- further classification ----------------------------- */

// Unsupervised training
void UnsupTrain(Subgraph *sg);
// Supervised training for complete graph
void SupTrainCompGraph(Subgraph *sg);
// Supervised training for knn-graph 
void SupTrainKnnGraph(Subgraph *sg);
// Semi-supervised training for knn-graph 
void SemiSupTrainKnnGraph(Subgraph *sg);

/*--- Classification functions ---------------------------------- */
Image *ImagePDF(Subgraph *sg, Features *f, int label);
Image *ImageLabel(Subgraph *sg, Image *img, int vol);
Image *ImageClassKnnGraph(Subgraph *sg, Features *f);
Scene *SceneClassKnnGraph(Subgraph *sg, Scene *mask, Features3 *f);
Image *ImageClassCompGraph(Subgraph *sg, Features *f);
Scene *SceneCluster(Subgraph *sg, Scene *scn, Scene *mask);

// Classify nodes of the evaluation/test set for KNNGraph
void ClassifyKnnGraph(Subgraph *sgtrain, Subgraph *sg);
// Classify nodes of evaluation/test set for CompGraph
void ClassifyCompGraph(Subgraph *sgtrain, Subgraph *sg);

/*------------ Auxiliary functions ------------------------------ */

// Compute normalized cut
float NormalizedCut( Subgraph *sg );
// Compute the influence zones of the pdf's maxima
void UnsupOPF(Subgraph *sg);
// Compute Euclidean istance between feature vectors
float EuclDist(float *f1, float *f2, int n);
// Compute  chi-squared distance between feature vectors
float ChiSquaredDist(float *f1, float *f2, int n);
// Create arcs for knn subgraph
void CreateArcs(Subgraph *sg, int knn);
// Create adjacent list between nodes with the same true lebale in subgraph: a knn graph 
void CreateArcsbyTrueLabel(Subgraph *sg, int knn);
// Create symmetric arcs for subgraph
void CreateSymArcs(Subgraph *sg);
// Destroy arcs in knn subgraph
void DestroyArcs(Subgraph *sg);
// Compute probability density function
void PDF(Subgraph *sg);
// Compute probability density function with space constraint
void SPDF(Subgraph *sg, Image *img);
// Create classification node
CNode *CreateCNode(Subgraph *sg);
// Destroy classification node
void DestroyCNode(CNode **node);
// Compute Markov features for classification node
void MarkovNodeFeat3(Scene *scn, Scene *mask, CNode *node, int p, float dm, int Imax);
// Compute the knn arcs for node
void NodeArcs(Subgraph *sg, CNode *node);
void NodeArcsByLabel(Subgraph *sg, CNode *node, int label);
void NodeArcsByTrueLabel(Subgraph *sg, CNode *node, int label);
// Compute node density for classification node 
void NodePD(Subgraph *sg, CNode *node);
// Compute node density for classification node by considering a desired label
void NodePDByLabel(Subgraph *sg, CNode *node, int label);
void NodePDByTrueLabel(Subgraph *sg, CNode *node, int label);
// Compute the predecessor in the best path to the classification node
int PredBestPath(Subgraph *sg, CNode *node);
// Compute the predecessor in the best path to the classification node using Bayes tie break
int PredBestPathBayesTie(Subgraph *sg, CNode *node);
// Find root node and check flatzone
int RootNode(Subgraph *sg, int p, char *flatzone);
void MSTPrototypes(Subgraph *sg);
// Compute the influnce zones of the maxima using their true labels
void SupOPF(Subgraph *sg);
void SupOPFwithErrors(Subgraph *sg); // allow errors in the training set
// Compute influence zones of prototypes
void SemiSupOPF(Subgraph *sg);
// Set features to the subgraph's nodes 
void SetSubgraphFeatures(Subgraph *sg, Features *f);
// Set 3D features to the subgraph's nodes 
void SetSubgraphFeatures3(Subgraph *sg, Features3 *f);
// Convert sample distances into integer arc weights
int ArcWeight(float dist);
// Copy subgraph (does not copy Arcs)
Subgraph *CopySubgraph(Subgraph *g);
// Split subgraph into two parts such that the size of the first part
// is given by a percentual of samples.
void SplitSubgraph(Subgraph *sg, Subgraph **sg1, Subgraph **sg2, float perc1);
// Takes the same amount of samples from subgraph sg and returns a subgraph composed of 
// them. The number of samples in given as a parameter.
Subgraph *SubSampleLabeledSubgraph( Subgraph *sg, int nsamples );
// Move errors above max_abs_error from the evaluation set to the
// training set, returning number of errors
int MoveErrorsToTrainingSet(Subgraph **sgtrain, Subgraph **sgeval, int max_abs_error);
// Remove nodes with distinct labels and distance=0
void RemoveInconsistentNodes(Subgraph **sg);
// Merge subgraphs with no arcs and feats
Subgraph *MergeSubgraphs(Subgraph *sg1, Subgraph *sg2);
//Compute accuracy of classification 
float Accuracy(Subgraph *g);
// Copy SNode
void CopySNode(SNode *dest, SNode *src, int nfeats);
//Swap nodes
void SwapSNode(SNode *a, SNode *b);
//Replace errors from evaluating set by non prototypes from training set
void SwapErrorsbyNonPrototypes(Subgraph **sgtrain, Subgraph **sgeval);
//Replace errors from evaluating set by randomly samples from training set
void SwapErrorsbySamples(Subgraph **sgtrain, Subgraph **sgeval);
//Resets subgraph fields (pred)
void ResetCompGraph(Subgraph *sg);
//Executes the learning procedure for CompGraph replacing the 
//missclassified samples in the evaluation set by non prototypes from
//training set
void LearningCompGraph(Subgraph **sgtrain, Subgraph **sgeval, int iterations);
//Executes the learning procedure for KnnGraph replacing the 
//missclassified samples in the evaluation set by non prototypes from
//training set
void LearningKnnGraph(Subgraph **sgtrain, Subgraph **sgeval, int iterations, int kmax);
// Compute the best k with minimum cut
void BestkMinCut(Subgraph *sg, int kmax);
// Compute the best k with minimum errors
void BestkMinError(Subgraph *sg, int kmax);
// Compute the best k with minimum cut for n clusters
void BestkMinCutNClusters(Subgraph *sg, int nclusters, int kmax);
// Eliminate maxima below H
void ElimMaxBelowH(Subgraph *sg, float H);
// Eliminate maxima below Area
void ElimMaxBelowArea(Subgraph *sg, int A);
// Eliminate maxima below Volume
void ElimMaxBelowVolume(Subgraph *sg, int V);
//read subgraph from opf format file
Subgraph *ReadSubgraph(char *file);
//write subgraph to disk
void WriteSubgraph(Subgraph *g, char *file);
//sort subgraph: order = INCREASING or DECREASING
void SortSubgraph(Subgraph **cg, char order);
/*read the precomputed distances between feature vectors*/
float **ReadDistances(char *fileName);
/*normalize features*/
void NormalizeFeatures(Subgraph *sg);

//---- Functions for pdf filtering

int       SgAncestor(int *dad, int *cmap, int rq);
int       SgRepresentative(int *cmap, int p);
SgCTree *CreateSgMaxTree(Subgraph *g);
void     DestroySgCTree(SgCTree **ctree);
int     *SgVolumeOpen(Subgraph *g, int thres);
int      SgVolumeLevel(SgCTree *ctree, int *level, int i, 
		       int thres, int cumvol);
int      *SgAreaOpen(Subgraph *g, int thres);
int       SgAreaLevel(SgCTree *ctree, int *level, int i, int thres);
void      SgCumSize(SgCTree *ctree, int i);



#endif
