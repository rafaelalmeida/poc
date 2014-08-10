#ifndef _SUBGRAPH_H_
#define _SUBGRAPH_H_


#include "cimage.h"
#include "scene.h"
#include "adjacency.h"
#include "adjacency3.h"
#include "set.h"
#include "heap.h"

#define MAXDENS 100000 // Maximum value for pdf computation
#define MAXARCW 100000 // Maximum arc weight for supervised
		       // classification using complete graphs

// Status of the graph nodes for classification

#define PROTOTYPE   1

typedef struct _felem {
  float *feat;
} FElem;

typedef struct _features {
  FElem *elem;
  int  nfeats;
  int  nelems;
  int  nrows,ncols;
} Features;

typedef struct _features3 {
  FElem *elem;
  int  nfeats;
  int  nelems;
  int  xsize,ysize,zsize;
} Features3;

typedef struct _cnode {
  int    dens;    // density
  float *feat;    // feature vector
  int   *adj;     // list of adjacent nodes
  float *dist;    // distance from adjacent nodes
  int    position; // position in the image/scene
} CNode;

typedef struct _snode {
  int dens;      // probability density value
  int pathval;   // path value
  int label;     // node label
  int pred;      // predecessor node
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
  double K1,K2;  // Constants for PDF computation
} Subgraph;

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
Subgraph *ChessSampl(Image *img);

// Get samples within a binary mask (optional) by choosing half below
// and half above a given threshold. 

Subgraph *ThresSampl3(Scene *scn, Scene *mask, int thres, int nnodes);
Subgraph *ThresNSampl3(Scene *scn, Scene *mask, int *thres, int N, int nnodes);

/*----------- Feature computation based on Markov properties -----------*/

// The intensity/color of a node and of its neighboors within a di
// radius in the image/scene are used to form its feature vector. 

void MarkovFeat3(Subgraph *sg, Scene *scn, float dm);

// Compute features for the image/cimage 

Features *MSImageFeats(Image *img, int nscales);
Features *MSCImageFeats(CImage *cimg, int nscales);
Features *MarkovImageFeats(Image *img, float dm);
Features *MarkovCImageFeats(CImage *cimg, float dm);
Features3 *MarkovSceneFeats(Scene *scn, float dm);
void      DestroyFeatures(Features **f);
void      DestroyFeatures3(Features3 **f);

/*--------- Training functions compute an optimum path forest for
  --------- further classification ----------------------------- */

// Compute arcs based on the best k, density function, and influence
// zones of the maxima
void UnsupTrain(Subgraph *sg, int kmax);
void UnsupTrainNClusters(Subgraph *sg, int nclusters, int kmax); 

// Compute supervised training for complete graph using optimum
// prototypes
void SupTrainCompGraph(Subgraph *sg);
// Compute supervised training for knn-graph 
void SupTrainKnnGraph(Subgraph *sg, int kmax);
// Compute semi-supervised training for knn-graph 
void SemiSupTrainKnnGraph(Subgraph *sg, int kmax);

/*--- Classification functions ---------------------------------- */

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
// Compute Euclidean distance
float EuclDist(float *f1, float *f2, int n);
// Create arcs for knn subgraph
void CreateArcs(Subgraph *sg, int knn);
// Create symmetric arcs for subgraph
void CreateSymArcs(Subgraph *sg);
// Destroy arcs in knn subgraph
void DestroyArcs(Subgraph *sg);
// Compute probability density function
void PDF(Subgraph *sg);
// Create classification node
CNode *CreateCNode(Subgraph *sg);
// Destroy classification node
void DestroyCNode(CNode **node);
// Compute Markov features for classification node
void MarkovNodeFeat3(Scene *scn, CNode *node, int p, float dm);
// Compute the knn arcs for node
void NodeArcs(Subgraph *sg, CNode *node);
// Compute node density for classification node 
void NodePD(Subgraph *sg, CNode *node);
// Compute the predecessor in the best path to the classification node
int PredBestPath(Subgraph *sg, CNode *node);
// Find root node and check flatzone
int RootNode(Subgraph *sg, int p, char *flatzone);
// Optimum prototype identification
void OptimumPrototypes(Subgraph *sg);
// Compute the influnce zones of the maxima using their true labels
void SupOPF(Subgraph *sg);
// Compute influence zones of prototypes
void SemiSupOPF(Subgraph *sg);
// Set features to the subgraph's nodes 
void SetSubgraphFeatures(Subgraph *sg, Features *f);
// Set 3D features to the subgraph's nodes 
void SetSubgraphFeatures3(Subgraph *sg, Features3 *f);
// Convert sample distances into integer arc weights
int ArcWeight(float dist);
// Split subgraph into two parts such that the size of the first part
// is given by a percentual of samples.
void SplitSubgraph(Subgraph *sg, Subgraph **sg1, Subgraph **sg2, float perc1);
// Move errors above max_abs_error from the evaluation set to the
// training set, returning number of errors
int MoveErrorsToTrainingSet(Subgraph **sgtrain, Subgraph **sgeval, int max_abs_error);
// Remove nodes with distinct labels and distance=0
void RemoveInconsistentNodes(Subgraph **sg);
// Merge subgraphs with no arcs and feats
Subgraph *MergeSubgraphs(Subgraph *sg1, Subgraph *sg2);
//Compute accuracy of classification for data sets
float Accuracy(Subgraph *g);
//Swap nodes
void SwapSNode(SNode *a, SNode *b);
//Replace errors from evaluating set by non prototypes from training set
void SwapErrorsbyNonPrototypes(Subgraph **sgtrain, Subgraph **sgeval);
//Resets subgraph fields (pred)
void ResetCompGraph(Subgraph *sg);
//Executes the learning procedure for CompGraph replacing the 
//missclassified samples in the evaluation set by non prototypes from
//training set
void LearningCompGraph(Subgraph **sgtrain, Subgraph **sgeval, int iterations);

#endif
