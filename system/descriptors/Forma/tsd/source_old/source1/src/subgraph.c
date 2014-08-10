#include "subgraph.h"

/*----------- Constructor and destructor ------------------------*/

// Allocate nodes without feature space and adjacency list 

Subgraph *CreateSubgraph(int nnodes)
{
  Subgraph *sg=(Subgraph *)calloc(1,sizeof(Subgraph));
  int i;

  sg->nnodes = nnodes;
  sg->node   = (SNode *)calloc(nnodes,sizeof(SNode));
  if (sg->node == NULL){
    Error("Cannot allocate nodes","CreateSubgraph");
  }

  for (i=0; i < sg->nnodes; i++) {
    sg->node[i].feat   = NULL;
    sg->node[i].adj    = NULL;
  }

  return(sg);
}

// Deallocate memory for subgraph

void DestroySubgraph(Subgraph **sg) 
{
  int i;

  if ((*sg)!=NULL){
    for (i=0; i < (*sg)->nnodes; i++) {
      if ((*sg)->node[i].feat != NULL)
	free((*sg)->node[i].feat);
      if ((*sg)->node[i].adj != NULL)
	DestroySet(&(*sg)->node[i].adj);
    }
    free((*sg)->node);
    free((*sg));
    *sg = NULL;
  }
}


// Create subgraph by uniform sampling 

Subgraph *UnifSampl(Image *img, Image *mask, int dx, int dy)
{
  Subgraph *sg=NULL;
  int x,y,nnodes,i,p;
  
  if (mask != NULL){
    nnodes=0;
    for (y=0; y < img->nrows; y=y+dy) 
      for (x=0; x < img->ncols; x=x+dx){
	p = x + img->tbrow[y];
	if (mask->val[p]){
	  nnodes++;
	}
      }
    sg = CreateSubgraph(nnodes);
    i=0;
    for (y=0; y < img->nrows; y=y+dy) 
      for (x=0; x < img->ncols; x=x+dx){
	p = x + img->tbrow[y];
	if (mask->val[p]){
	  sg->node[i].position = p;
	  i++;
	}
      }
  }else{
    nnodes=0;
    for (y=0; y < img->nrows; y=y+dy) 
      for (x=0; x < img->ncols; x=x+dx){
	p = x + img->tbrow[y];
	nnodes++;
      }  
    sg = CreateSubgraph(nnodes);
    i=0;
    for (y=0; y < img->nrows; y=y+dy) 
      for (x=0; x < img->ncols; x=x+dx){
	p = x + img->tbrow[y];
	sg->node[i].position = p;
	i++;
      }
  }

  return(sg);
}

Subgraph *UnifSampl3(Scene *scn, Scene *mask, int dx, int dy, int dz)
{
  Subgraph *sg=NULL;
  int x,y,z,nnodes,i,p;
  
  if (mask != NULL){
    nnodes=0;
    for (z=0; z < scn->zsize; z=z+dz) 
      for (y=0; y < scn->ysize; y=y+dy){
	for (x=0; x < scn->xsize; x=x+dx){
	  p = x + scn->tby[y] + scn->tbz[z];
	  if (mask->data[p]){
	    nnodes++;
	  }
	}
      }
    sg = CreateSubgraph(nnodes);
    i=0;
    for (z=0; z < scn->zsize; z=z+dz) 
      for (y=0; y < scn->ysize; y=y+dy){
	for (x=0; x < scn->xsize; x=x+dx){
	  p = x + scn->tby[y] + scn->tbz[z];
	  if (mask->data[p]){
	    sg->node[i].position = p;
	    i++;
	  }
	}
      }
  }else{
    nnodes=0;
    for (z=0; z < scn->zsize; z=z+dz) 
      for (y=0; y < scn->ysize; y=y+dy){
	for (x=0; x < scn->xsize; x=x+dx){
	  p = x + scn->tby[y] + scn->tbz[z];
	  nnodes++;
	}
      }
    sg = CreateSubgraph(nnodes);
    i=0;
    for (z=0; z < scn->zsize; z=z+dz) 
      for (y=0; y < scn->ysize; y=y+dy){
	for (x=0; x < scn->xsize; x=x+dx){
	  p = x + scn->tby[y] + scn->tbz[z];
	  sg->node[i].position = p;
	  i++;
	}
      }
  }

  return(sg);
}

// Compute samples for image compression
Subgraph *ChessSampl(Image *img)
{
  Subgraph *sg=NULL;
  int x,y,nnodes,i,p,xo;
  
  nnodes=0;
  xo = 1;
  for (y=0; y < img->nrows-1; y=y+1){ 
    for (x=xo; x < img->ncols-1; x=x+2){
      nnodes++;
    }
    xo = 1 - xo; 
  }

  sg = CreateSubgraph(nnodes);
  sg->nlabels=MaximumValue(img)+1;

  i=0;
  xo = 1;
  for (y=0; y < img->nrows-1; y=y+1){ 
    for (x=xo; x < img->ncols-1; x=x+2){
      p = x + img->tbrow[y];      
      sg->node[i].position  = p;
      sg->node[i].truelabel = img->val[p];      
      i++;
    }
    xo = 1 - xo; 
  }

  return(sg);
}

// Create subgraph by threshold sampling

Subgraph *ThresSampl3(Scene *scn, Scene *mask, int thres, int nnodes)
{
  Subgraph *sg=NULL;
  Scene *used=CreateScene(scn->xsize,scn->ysize,scn->zsize);
  int i,p,n;
  char side;

  sg   = CreateSubgraph(nnodes);
  i    = 0;
  n    = scn->xsize*scn->ysize*scn->zsize;
  side = 0; 

  if (mask != NULL){
    while (i < nnodes){
      p = RandomInteger(0,n-1);
      if (mask->data[p]){
	if (used->data[p]==0){
	  if ((side==0)&&(scn->data[p]<=thres)){
	    side = 1; // set side to above
	    sg->node[i].position = p;
	    used->data[p]=1;
	    i++;
	  }else{
	    if ((side==1)&&(scn->data[p]>thres)){
	      side = 0; // set side to below
	      sg->node[i].position = p;
	      used->data[p]=1;
	      i++;
	    }
	  }
	}
      }
    }
  }else{
    while (i < nnodes){
      p = RandomInteger(0,n-1);
      if (used->data[p]==0){
	if ((side==0)&&(scn->data[p]<=thres)){
	  side = 1; // set side to above
	  sg->node[i].position = p;
	  used->data[p]=1;
	  i++;
	}else{
	  if ((side==1)&&(scn->data[p]>thres)){
	    side = 0; // set side to below
	    sg->node[i].position = p;
	    used->data[p]=1;
	    i++;
	  }
	}
      }
    }
  }
  DestroyScene(&used);

  return(sg);
}

Subgraph *ThresNSampl3(Scene *scn, Scene *mask, int *thres, int N, int nnodes)
{
  Subgraph *sg=NULL;
  Scene *used=CreateScene(scn->xsize,scn->ysize,scn->zsize);
  int i,p,n;
  int side;

  sg   = CreateSubgraph(nnodes);
  i    = 0;
  n    = scn->xsize*scn->ysize*scn->zsize;
  side = 0;

  if (mask != NULL){
    while (i < nnodes){
      p = RandomInteger(0,n-1);
      if (mask->data[p]){
	if (used->data[p]==0){
	  if ((side==0)&&(scn->data[p]<=thres[0])){
	    side = 1; // set side to above
	    sg->node[i].position = p;
	    used->data[p]=1;
	    i++;
	  }else if ((side<N-1)&&(scn->data[p]>thres[side-1])&&(scn->data[p]<=thres[side])){
	    side++; // set side to above
	    sg->node[i].position = p;
	    used->data[p]=1;
	    i++;
	  }else{
	    if ((side==N-1)&&(scn->data[p]>thres[N-1])){
	      side = 0; // set side to below
	      sg->node[i].position = p;
	      used->data[p]=1;
	      i++;
	    }
	  }
	}
      }
    }
  }else{
    while (i < nnodes){
      p = RandomInteger(0,n-1);
      if (used->data[p]==0){
	if ((side==0)&&(scn->data[p]<=thres[0])){
	  side = 1; // set side to above
	  sg->node[i].position = p;
	  used->data[p]=1;
	  i++;
	}else if ((side<N-1)&&(scn->data[p]>thres[side-1])&&(scn->data[p]<=thres[side])){
	  side++; // set side to above
	  sg->node[i].position = p;
	  used->data[p]=1;
	  i++;
	}else{
	  if ((side==N-1)&&(scn->data[p]>thres[N-1])){
	    side = 0; // set side to below
	    sg->node[i].position = p;
	    used->data[p]=1;
	    i++;
	  }
	}
      }
    }
  }
  DestroyScene(&used);

  for (i=0; i < sg->nnodes; i++) 
    sg->node[i].status = 0; // training node


  return(sg);
}

/*----------- Feature computation based on Markov properties -----------*/

// The intensity/color of a node and of its neighboors within a dm
// radius in the image/scene are used to form its feature vector. It
// assumes that the sampling has been done.

void MarkovFeat3(Subgraph *sg, Scene *scn, float dm)
{
  AdjRel3 *A;
  int     p,q,i,j,xysize=scn->xsize*scn->ysize,aux;
  Voxel   u,v;

  if (dm >= 1.0){
    A=Spheric(dm);
    sg->dm=dm;
    sg->nfeats = A->n - 1;
    for (i=0; i < sg->nnodes; i++) 
      sg->node[i].feat = AllocFloatArray(sg->nfeats);

    for (i=0; i < sg->nnodes; i++) {
      p   = sg->node[i].position;
      aux = p%xysize;
      u.x = aux%scn->xsize;
      u.y = aux/scn->xsize;
      u.z = p/xysize;
      for (j=1; j < A->n; j++) {
	v.x = u.x + A->dx[j];
	v.y = u.y + A->dy[j];
	v.z = u.z + A->dz[j];
	if (ValidVoxel(scn,v.x,v.y,v.z)){
	  q = v.x + scn->tby[v.y] + scn->tbz[v.z];
	  sg->node[i].feat[j]=scn->data[q];
	}else
	  sg->node[i].feat[j]=scn->data[p];
      }
    }
    DestroyAdjRel3(&A);
  }else{
    sg->dm=0.0;
    sg->nfeats = 1;
    for (i=0; i < sg->nnodes; i++) 
      sg->node[i].feat = AllocFloatArray(sg->nfeats);

    for (i=0; i < sg->nnodes; i++) {
      p   = sg->node[i].position;
      sg->node[i].feat[0]=scn->data[p];
    }
  }
}

Features *MSImageFeats(Image *img, int nscales)
{
  Features *f=(Features *)calloc(1,sizeof(Features));
  AdjRel   *A=NULL;
  int       s,i,Imax=MaximumValue(img);
  Image    *img1,*img2;

  f->ncols = img->ncols;
  f->nrows = img->nrows;
  f->nelems = f->ncols*f->nrows;
  f->elem   = (FElem *)calloc(f->nelems,sizeof(FElem));
  for (i=0; i < f->nelems; i++) {
    f->elem[i].feat = AllocFloatArray(nscales);
    f->nfeats       = nscales;
  }
      
  img1  = CopyImage(img);

  for (s=1; s <= nscales; s=s+1) {
    A    = Circular(s);
    img2 = AsfOCRec(img1,A);
    for (i=0; i < f->nelems; i++) {
      f->elem[i].feat[s-1] = (float)img2->val[i]/(float)Imax;
    }
    DestroyImage(&img2);
    DestroyAdjRel(&A);
  }
  DestroyImage(&img1);
  return(f);
}

Features *MSCImageFeats(CImage *cimg, int nscales)
{
  Features *f=(Features *)calloc(1,sizeof(Features));
  AdjRel   *A=NULL;
  int       s,i,j;
  Image    *img;

  f->ncols  = cimg->C[0]->ncols;
  f->nrows  = cimg->C[1]->ncols;
  f->nelems = cimg->C[0]->ncols*cimg->C[0]->nrows;
  f->elem   = (FElem *)calloc(f->nelems,sizeof(FElem));
  for (i=0; i < f->nelems; i++) {
    f->elem[i].feat = AllocFloatArray(3*nscales);
    f->nfeats       = 3*nscales;
  }

  for (j=0; j < 3; j=j+1) {
    for (s=1; s <= nscales; s=s+1) {
      A   = Circular(s);
      img = AsfCORec(cimg->C[j],A);
      for (i=0; i < f->nelems; i++) {
	f->elem[i].feat[s-1+(j*nscales)] = (float)img->val[i]/255.0;
      }
      DestroyImage(&img);
      DestroyAdjRel(&A);
    }
  }
  return(f);
}

Features *MarkovImageFeats(Image *img, float dm)
{
  Features *f=(Features *)calloc(1,sizeof(Features));
  AdjRel   *A=NULL;
  int       p,q,i,j;
  Pixel     u,v;

  if (dm >= 1.0){
    A=Circular(dm);
    f->ncols  = img->ncols;
    f->nrows  = img->nrows;
    f->nelems = f->ncols*f->nrows;
    f->nfeats = A->n-1;
    f->elem   = (FElem *)calloc(f->nelems,sizeof(FElem));
    for (i=0; i < f->nelems; i++) {
      f->elem[i].feat = AllocFloatArray(f->nfeats);
    }

    for (p=0; p < f->nelems; p++) {
      u.x = p%img->ncols;
      u.y = p/img->ncols;
      for (j=1; j < A->n; j++) {
	v.x = u.x + A->dx[j];
	v.y = u.y + A->dy[j];
	if (ValidPixel(img,v.x,v.y)){
	  q = v.x + img->tbrow[v.y];
	  f->elem[p].feat[j-1]=img->val[q];
	}else
	  f->elem[p].feat[j-1]=img->val[p];
      }
    }
    DestroyAdjRel(&A);
  }else{
    f->ncols  = img->ncols;
    f->nrows  = img->nrows;
    f->nelems = f->ncols*f->nrows;
    f->nfeats = 1;
    f->elem   = (FElem *)calloc(f->nelems,sizeof(FElem));
    for (i=0; i < f->nelems; i++) {
      f->elem[i].feat = AllocFloatArray(f->nfeats);
    }
    
    for (p=0; p < f->nelems; p++) {
      f->elem[p].feat[0]=img->val[p];
    }
  }

  return(f);
}

Features *MarkovCImageFeats(CImage *cimg, float dm)
{
  Features *f=(Features *)calloc(1,sizeof(Features));
  AdjRel   *A=NULL;
  int       p,q,i,j,k;
  Pixel     u,v;
  Image    *img=cimg->C[0];

  if (dm >= 1.0){
    A=Circular(dm);
    f->ncols  = img->ncols;
    f->nrows  = img->nrows;
    f->nelems = f->ncols*f->nrows;
    f->nfeats = 3*A->n-3;
    f->elem   = (FElem *)calloc(f->nelems,sizeof(FElem));
    for (i=0; i < f->nelems; i++) {
      f->elem[i].feat = AllocFloatArray(f->nfeats);
    }

    for (p=0; p < f->nelems; p++) {
      u.x = p%img->ncols;
      u.y = p/img->ncols;
      for (j=1, k=0; j < A->n; j++, k=k+3) {
	v.x = u.x + A->dx[j];
	v.y = u.y + A->dy[j];
	if (ValidPixel(img,v.x,v.y)){
	  q = v.x + img->tbrow[v.y];
	  f->elem[p].feat[k]=cimg->C[0]->val[q];
	  f->elem[p].feat[k+1]=cimg->C[1]->val[q];
	  f->elem[p].feat[k+2]=cimg->C[2]->val[q];
	}else{
	  f->elem[p].feat[k]=cimg->C[0]->val[p];
	  f->elem[p].feat[k+1]=cimg->C[1]->val[p];
	  f->elem[p].feat[k+2]=cimg->C[2]->val[p];
	}
      }
    }
    DestroyAdjRel(&A);
  }else{
    f->ncols  = img->ncols;
    f->nrows  = img->nrows;
    f->nelems = f->ncols*f->nrows;
    f->nfeats = 3;
    f->elem   = (FElem *)calloc(f->nelems,sizeof(FElem));
    for (i=0; i < f->nelems; i++) {
      f->elem[i].feat = AllocFloatArray(f->nfeats);
    }
    
    for (p=0; p < f->nelems; p++) {
      f->elem[p].feat[0]=cimg->C[0]->val[p];
      f->elem[p].feat[1]=cimg->C[1]->val[p];
      f->elem[p].feat[2]=cimg->C[2]->val[p];
    }
  }
  return(f);
}

void DestroyFeatures(Features **f)
{
  int i;

  if ((*f)!=NULL){
    for (i=0; i < (*f)->nelems; i++) 
      free((*f)->elem[i].feat);
    free((*f)->elem);
    free(*f);
    (*f)=NULL;
  }
}

Features3 *MarkovSceneFeats(Scene *scn, float dm)
{
  Features3 *f=(Features3 *)calloc(1,sizeof(Features3));
  AdjRel3   *A=NULL;
  int       p,q,i,j,xysize=scn->xsize*scn->ysize,aux;
  Voxel     u,v;

  if (dm >= 1.0){
    A=Spheric(dm);
    f->xsize  = scn->xsize;
    f->ysize  = scn->ysize;
    f->zsize  = scn->zsize;
    f->nelems = f->xsize*f->ysize*f->zsize;
    f->nfeats = A->n-1;
    f->elem   = (FElem *)calloc(f->nelems,sizeof(FElem));
    for (i=0; i < f->nelems; i++) {
      f->elem[i].feat = AllocFloatArray(f->nfeats);
    }

    for (p=0; p < f->nelems; p++) {
      aux = p%xysize;
      u.x = aux%scn->xsize;
      u.y = aux/scn->xsize;
      u.z = p/xysize;
      for (j=1; j < A->n; j++) {
	v.x = u.x + A->dx[j];
	v.y = u.y + A->dy[j];
	v.z = u.z + A->dz[j];
	if (ValidVoxel(scn,v.x,v.y,v.z)){
	  q = v.x + scn->tby[v.y] + scn->tbz[v.z];
	  f->elem[p].feat[j-1]=scn->data[q];
	}else
	  f->elem[p].feat[j-1]=scn->data[p];
      }
    }
    DestroyAdjRel3(&A);
  }else{
    f->xsize  = scn->xsize;
    f->ysize  = scn->ysize;
    f->zsize  = scn->zsize;
    f->nelems = f->xsize*f->ysize*f->zsize;
    f->nfeats = 1;
    f->elem   = (FElem *)calloc(f->nelems,sizeof(FElem));
    for (i=0; i < f->nelems; i++) {
      f->elem[i].feat = AllocFloatArray(f->nfeats);
    }
    
    for (p=0; p < f->nelems; p++) {
      f->elem[p].feat[0]=scn->data[p];
    }
  }

  return(f);
}

void DestroyFeatures3(Features3 **f)
{
  int i;

  if ((*f)!=NULL){
    for (i=0; i < (*f)->nelems; i++) 
      free((*f)->elem[i].feat);
    free((*f)->elem);
    free(*f);
    (*f)=NULL;
  }
}


/*--------- Training functions compute an optimum path forest for
  --------- further classification ----------------------------- */

// Compute arcs based on the best k --- the one which minimizes the
// normalized cut. Then compute density function and influence zones
// of the maxima.

void UnsupTrain(Subgraph *sg, int kmax)
{
  int k;
  float mincut=REAL_MAX,nc;

  // Find the best k

  for (k=1; (k <= kmax)&&(mincut != 0.0); k++) {
    CreateArcs(sg,k);
    PDF(sg);
    UnsupOPF(sg);
    nc = NormalizedCut(sg);
    if (nc < mincut){
      mincut=nc;
      sg->bestk =k;
    }
    DestroyArcs(sg);
  }

  printf("bestk %d ncut %f \n",sg->bestk,mincut);
  CreateArcs(sg,sg->bestk);
  PDF(sg);
  UnsupOPF(sg);
}

void UnsupTrainNClusters(Subgraph *sg, int nclusters, int kmax)
{
  int k, *nlabels=AllocIntArray(kmax+1);
  float  mincut=REAL_MAX, *nclist=AllocFloatArray(kmax+1);

  // Compute cuts for all values of k

  for (k=1; (k <= kmax); k++) {
    CreateArcs(sg,k);
    PDF(sg);
    UnsupOPF(sg);
    nclist[k] = NormalizedCut(sg);
    nlabels[k] = sg->nlabels;
    if (nclist[k] < mincut){
      sg->bestk=k;
      mincut=nclist[k];
    }
    DestroyArcs(sg);
  }

  if (nlabels[sg->bestk]!=nclusters){ // mincut does not give the
				     // desired number of clusters
    
    mincut=REAL_MAX;
    sg->bestk=0;

    for (k=1; (k <= kmax); k++) { // select k which gives nclusters
				  // and lowest cut
      if (nlabels[k]==nclusters){
	if (nclist[k] < mincut){
	  sg->bestk=k;
	  mincut=nclist[k];
	}
      }
    }
    if (sg->bestk == 0){ // there is no such a k, then select k which
			// gives more than nclusters and lowest cut.

      
      for (k=1; (k <= kmax); k++) { 
	if (nlabels[k]>nclusters){
	  if (nclist[k] < mincut){
	    sg->bestk=k;
	    mincut=nclist[k];
	  }
	}
      }
      
      if (sg->bestk == 0) { // there is no such a k, then select k=1.
	sg->bestk = 1;
      }
    }
  }
  free(nclist);
  free(nlabels);
  printf("bestk is %d \n",sg->bestk);
  CreateArcs(sg,sg->bestk);
  PDF(sg);
  UnsupOPF(sg);
}

// Supervised training using a complete graph

void SupTrainCompGraph(Subgraph *sg)
{
  int p,q,tmp,weight;
  GQueue *Q = NULL;
  int *pathval = NULL;

  // compute optimum prototypes

  OptimumPrototypes(sg);

  // initialization

  pathval = AllocIntArray(sg->nnodes);
  Q = CreateGQueue(MAXARCW+2, sg->nnodes, pathval);

  for (p = 0; p < sg->nnodes; p++) {
    if (sg->node[p].status==PROTOTYPE){
      pathval[p]         = 0;
      sg->node[p].label  = sg->node[p].truelabel;
      sg->node[p].pred   = NIL;
      InsertGQueue( &Q, p );
    }else{ // non-prototypes
      pathval[p] = INT_MAX;
    }
  }

  // IFT with fmax

  while ( !EmptyGQueue( Q ) ) {
    p  = RemoveGQueue( Q );

    sg->node[p].pathval = pathval[p];

    for (q=0; q < sg->nnodes; q++){ 
      if (p!=q){
	if (pathval[p] < pathval[q]){
	  if (sg->node[p].truelabel!=sg->node[q].truelabel)
	    weight = MAXARCW+1;
	  else
	    weight = ArcWeight(EuclDist(sg->node[p].feat,sg->node[q].feat,sg->nfeats)/sg->df);
	  tmp  = MAX(pathval[p],weight);
	  if ( tmp < pathval[ q ] ) {
	    if (pathval[q]!=INT_MAX)
	      RemoveGQueueElem(Q, q);
	    pathval[q]        = tmp;
	    sg->node[q].pred  = p;
	    sg->node[q].label = sg->node[p].label;
	    InsertGQueue( &Q, q );
	  }
	}
      }
    }
  }

  DestroyGQueue( &Q );
  free( pathval );

  // Verify if error is really zero in training set and reset status
  // of the prototypes.

  for (p=0; p < sg->nnodes; p++){
    sg->node[p].status = 0;
    if (sg->node[p].label!=sg->node[p].truelabel)
      Warning("Class. error in the training set","SupTrainCompGraph");
  }
}

// Compute arcs based on the best k --- the one which minimizes errors
// in the training set. Then, compute density function and influence
// zones of the maxima.

void SupTrainKnnGraph(Subgraph *sg, int kmax)
{
  int k,p;
  int min_nerrors=INT_MAX,nerrors;

  // Find the best k
  
  for (k=1; (k <= kmax)&&(min_nerrors!=0); k++) {
    CreateArcs(sg,k);
    PDF(sg);
    SupOPF(sg);
    nerrors=0;
    for (p=0; p < sg->nnodes; p++)
      if (sg->node[p].label!=sg->node[p].truelabel)
	nerrors++;
    if (nerrors < min_nerrors){
      min_nerrors=nerrors;
      sg->bestk =k;
    }
    DestroyArcs(sg);
  }
  
  printf("bestk %d nerrors %d \n",sg->bestk,min_nerrors);
  
  CreateArcs(sg,sg->bestk);
  PDF(sg); 
  SupOPF(sg);
    
}

void SemiSupTrainKnnGraph(Subgraph *sg, int kmax)
{
  int k;
  float mincut=REAL_MAX,nc;


  // Find the best k

  for (k=1; (k <= kmax)&&(mincut != 0.0); k++) {
    CreateArcs(sg,k);
    PDF(sg);
    UnsupOPF(sg);
    nc = NormalizedCut(sg);
    if (nc < mincut){
      mincut=nc;
      sg->bestk =k;
    }
    DestroyArcs(sg);
  }

  printf("bestk %d mincut %f \n",sg->bestk,mincut);

  CreateArcs(sg,sg->bestk);
  PDF(sg); 
  SemiSupOPF(sg);

}

/*--------- Classification functions ----------------------------- */

Image *ImageClassKnnGraph(Subgraph *sg, Features *f)
{
  CNode  *node=CreateCNode(sg);
  int     i,p,pred;
  Image  *label=CreateImage(f->ncols,f->nrows);

  for (p=0; p < f->nelems; p++) {
    for (i=0; i < sg->nfeats; i++)
      node->feat[i]=f->elem[p].feat[i];
    NodeArcs(sg,node);
    NodePD(sg,node);
    pred=PredBestPath(sg,node);
    label->val[p]=sg->node[pred].label;
  }
  
  DestroyCNode(&node);

  return(label);
} 

Scene *SceneClassKnnGraph(Subgraph *sg, Scene *mask, Features3 *f)
{
  CNode  *node=CreateCNode(sg);
  int     i,p,pred;
  Scene  *label=CreateScene(f->xsize,f->ysize,f->zsize);

  if (mask!=NULL) {
    for (p=0; p < f->nelems; p++) {
      if (mask->data[p]) {
	for (i=0; i < sg->nfeats; i++)
	  node->feat[i]=f->elem[p].feat[i];
	NodeArcs(sg,node);
	NodePD(sg,node);
	pred=PredBestPath(sg,node);
	label->data[p]=sg->node[pred].label;
      }
    }
  }else{
    for (p=0; p < f->nelems; p++) {
      for (i=0; i < sg->nfeats; i++)
	node->feat[i]=f->elem[p].feat[i];
      NodeArcs(sg,node);
      NodePD(sg,node);
      pred=PredBestPath(sg,node);
      label->data[p]=sg->node[pred].label;
    }
  }
  
  DestroyCNode(&node);

  return(label);
}

// Classify nodes of evaluation/test set
void ClassifyKnnGraph(Subgraph *sgtrain, Subgraph *sg)
{
  CNode  *node=CreateCNode(sgtrain);
  int     i,p,pred;

  for (p=0; p < sg->nnodes; p++) {
    for (i=0; i < sg->nfeats; i++)
      node->feat[i]=sg->node[p].feat[i];
    NodeArcs(sgtrain,node);
    NodePD(sgtrain,node);
    pred=PredBestPath(sgtrain,node);
    sg->node[p].label=sgtrain->node[pred].label;
  }

  DestroyCNode(&node);
} 

Image *ImageClassCompGraph(Subgraph *sg, Features *f)
{
  Image *label=CreateImage(f->ncols,f->nrows);
  int p,q,n=f->nelems,tmp,minval,weight;


  for (q=0; q < n; q++) {
    if (sg->node[0].position!=q){
      weight = ArcWeight(EuclDist(f->elem[q].feat,sg->node[0].feat,sg->nfeats)/sg->df);
      minval = MAX(sg->node[0].pathval,weight);
      label->val[q] = sg->node[0].label;
      for (p=1; p < sg->nnodes; p++) {
	if (sg->node[p].position!=q){
	  weight = ArcWeight(EuclDist(f->elem[q].feat,sg->node[p].feat,sg->nfeats)/sg->df);
	  tmp  = MAX(sg->node[p].pathval,weight);
	  if (tmp < minval) {
	    label->val[q]=sg->node[p].label;
	  }
	}else{
	  label->val[q]=sg->node[p].label;
	  break;
	}
      }
    }else{
      label->val[q]=sg->node[0].label;
    }
  }
  return(label);
}


Scene *SceneCluster(Subgraph *sg, Scene *scn, Scene *mask)
{
  CNode  *node=CreateCNode(sg);
  int     p,n,pred;
  Scene  *label=CreateScene(scn->xsize, scn->ysize, scn->zsize);

  n = scn->xsize * scn->ysize * scn->zsize;
  
  for (p=0; p < n; p++) {
    if (mask->data[p]){
      MarkovNodeFeat3(scn,node,p,sg->dm);
      NodeArcs(sg,node);
      NodePD(sg,node);
      pred=PredBestPath(sg,node);
      label->data[p]=sg->node[pred].label;
    }
  }
  DestroyCNode(&node);

  return(label);
} 

/* ----------  Auxiliary functions ------------------------------- */

// Normalized cut

float NormalizedCut( Subgraph *sg ) {
  int l, p, q;
  Set *Saux;
  float ncut, dist;
  float *acumIC; //acumulate weights inside each class
  float *acumEC; //acumulate weights between the class and a distinct one

  ncut = 0.0;
  acumIC = AllocFloatArray( sg->nlabels );
  acumEC = AllocFloatArray( sg->nlabels );

  for ( p = 0; p < sg->nnodes; p++ ) {
    for ( Saux = sg->node[ p ].adj; Saux != NULL; Saux = Saux->next ) {
      q = Saux->elem;
      dist = (float) EuclDist(sg->node[p].feat,sg->node[q].feat,sg->nfeats);
      dist = sqrtf( dist );
      if ( dist != 0.0 ) {
	if ( sg->node[ p ].label == sg->node[ q ].label ) {
	  acumIC[ sg->node[ p ].label ] += 1.0 / dist; // intra-class weight
	} else { // inter - class weight
	  acumEC[ sg->node[ p ].label ] += 1.0 / dist; // inter-class weight
	}
      }
    }
  }

  for( l = 0; l < sg->nlabels; l++ ) {
    if ( ( acumIC[ l ] + acumEC[ l ] ) != 0.0 ) ncut += (float) acumEC[ l ] / ( acumIC[ l ] + acumEC[ l ] );
  }
  free( acumEC );
  free( acumIC );
  return( ncut );
}

// Influence zones of the maxima

void UnsupOPF(Subgraph *sg) {
  int p, q, l, tmp;
  GQueue *Q = NULL;
  Set *Saux=NULL;
  int *pathval = NULL, rp, rq;
  char flatzone_p,flatzone_q;

  pathval = AllocIntArray(sg->nnodes);
  Q = CreateGQueue(sg->maxdens + 1, sg->nnodes, pathval);
  SetRemovalPolicy(Q, MAXVALUE);

  for (p = 0; p < sg->nnodes; p++) {
    pathval[ p ] = sg->node[ p ].dens - 1; // 1-dome
    sg->node[ p ].pred  = NIL;
    InsertGQueue( &Q, p );
  }
  
  while ( !EmptyGQueue( Q ) ) {
    p  = RemoveGQueue( Q );
    rp = RootNode(sg,p,&flatzone_p);

    if ( sg->node[ p ].pred == NIL ) {
      pathval[ p ] = sg->node[ p ].dens;
    }

    sg->node[ p ].pathval = pathval[ p ];

    for ( Saux = sg->node[ p ].adj; Saux != NULL; Saux = Saux->next ) {
      q = Saux->elem;
      if ( pathval[ p ] > pathval[ q ] ) {
	tmp = MIN( pathval[ p ], sg->node[ q ].dens );
	if ( tmp > pathval[ q ] ) {
	  UpdateGQueue( &Q, q, tmp ); 
	  sg->node[ q ].pred  = p; 
	}
      } else {
	if ( pathval[ p ] == pathval[ q ] ) {
	  rq = RootNode(sg,q,&flatzone_q);
	  if ((rp != rq)&& 
	      (flatzone_p)&&
	      (flatzone_q)&& 
	      (sg->node[rp].dens==sg->node[rq].dens)){// merge maxima
	    sg->node[rq].pred = rp;
	  }
	}
      }
    }
  }
  DestroyGQueue( &Q );
  free( pathval );

  l = 0;
  for (p=0; p < sg->nnodes; p++) 
    if (sg->node[p].pred==NIL){
      sg->node[p].label = l;
      l++;
    }

  for (p=0; p < sg->nnodes; p++){ // label propagation
    rp = RootNode(sg,p,&flatzone_p);
    sg->node[p].label = sg->node[rp].label;
  }
  sg->nlabels = l;
}

// Influence zones of the maxima with their true labels

void SupOPF(Subgraph *sg) {
  int p, q, tmp;
  GQueue *Q = NULL;
  Set *Saux=NULL;
  int *pathval = NULL, rp, rq;
  char flatzone_p,flatzone_q;

  pathval = AllocIntArray(sg->nnodes);
  Q = CreateGQueue(sg->maxdens + 1, sg->nnodes, pathval);
  SetRemovalPolicy(Q, MAXVALUE);

  for (p = 0; p < sg->nnodes; p++) {
    pathval[ p ] = sg->node[ p ].dens - 1; // 1-dome
    sg->node[ p ].pred    = NIL;
    InsertGQueue( &Q, p );
  }
  
  while ( !EmptyGQueue( Q ) ) {
    p  = RemoveGQueue( Q );
    rp = RootNode(sg,p,&flatzone_p);

    if ( sg->node[ p ].pred == NIL ) {
      pathval[ p ] = sg->node[ p ].dens;
    }

    sg->node[ p ].pathval = pathval[ p ];

    for ( Saux = sg->node[ p ].adj; Saux != NULL; Saux = Saux->next ) {
      q = Saux->elem;
      if ( pathval[ p ] > pathval[ q ] ) {
	tmp = MIN( pathval[ p ], sg->node[ q ].dens );
	if ( tmp > pathval[ q ] ) {
	  UpdateGQueue( &Q, q, tmp ); 
	  sg->node[ q ].pred   = p; 
	}
      } else {
	if ( pathval[ p ] == pathval[ q ] ) {
	  rq = RootNode(sg,q,&flatzone_q);
	  if ((rp != rq)&& 
	      (sg->node[rp].truelabel == sg->node[rq].truelabel)&&
	      (flatzone_p)&&
	      (flatzone_q)&& 
	      (sg->node[rp].dens==sg->node[rq].dens)){// merge maxima
	    sg->node[rq].pred = rp;
	  }
	}
      }
    }
  }
  DestroyGQueue( &Q );
  free( pathval );

  for (p=0; p < sg->nnodes; p++) 
    if (sg->node[p].pred==NIL){
      sg->node[p].label = sg->node[p].truelabel;
    }

  for (p=0; p < sg->nnodes; p++){ // label propagation
    rp = RootNode(sg,p,&flatzone_p);
    sg->node[p].label = sg->node[rp].label;
  }
}

// Influence zones of the prototypes

void SemiSupOPF(Subgraph *sg) {
  int p, q, tmp;
  GQueue *Q = NULL;
  Set *Saux=NULL;
  int *pathval = NULL;

  pathval = AllocIntArray(sg->nnodes);
  Q = CreateGQueue(sg->maxdens + 1, sg->nnodes, pathval);
  SetRemovalPolicy(Q, MAXVALUE);

  for (p = 0; p < sg->nnodes; p++) {
    if (sg->node[p].status==PROTOTYPE){
      pathval[ p ] = sg->node[p].dens; 
      sg->node[ p ].pred    = NIL;
      sg->node[ p ].label   = sg->node[p].truelabel;
      InsertGQueue( &Q, p );
    }else
      sg->node[p].pathval = pathval[p]=INT_MIN;
  }

  while ( !EmptyGQueue( Q ) ) {
    p  = RemoveGQueue( Q );

    sg->node[ p ].pathval = pathval[ p ];

    for ( Saux = sg->node[ p ].adj; Saux != NULL; Saux = Saux->next ) {
      q = Saux->elem;
      if ( pathval[ p ] > pathval[ q ] ) {
	tmp = MIN( pathval[ p ], sg->node[ q ].dens );
	if ( tmp > pathval[ q ] ) {
	  if (pathval[q]!=INT_MIN)
	    RemoveGQueueElem( Q, q); 
	  sg->node[ q ].pred   = p;
	  sg->node[ q ].label  = sg->node[p].label;
	  pathval[q]           = tmp;
	  InsertGQueue( &Q, q );
	}
      }
    }
  }
  DestroyGQueue( &Q );
  free( pathval );

}

// Find root node and identify if the path is on a flatzone

int RootNode(Subgraph *sg, int p, char *flatzone)
{
  int density=sg->node[p].dens;

  *flatzone=1;
  while (sg->node[p].pred!=NIL){
    if (sg->node[p].dens!=density)
      *flatzone=0; 
    p = sg->node[p].pred;
  }
  return(p);
}

// Compute prototypes by setting their label as true label. Other
// nodes have their label set to NIL. It also computes maximum arc
// weight in sg->df for the subsequent IFT.

void OptimumPrototypes(Subgraph *sg) {
  int p,q,weight;
  GQueue *Q = NULL;
  int *pathval = NULL;
  float dist;

  // maximum arc weight

  sg->df=REAL_MIN;
  for (p=0; p < sg->nnodes; p++){
    for (q=0; q < sg->nnodes; q++){
      if (p!=q){
	dist = EuclDist(sg->node[p].feat,sg->node[q].feat,sg->nfeats);
	if (dist > sg->df)
	  sg->df=dist;
      }
    }
  }
  
  // initialization

  pathval = AllocIntArray(sg->nnodes);
  Q = CreateGQueue(MAXARCW + 1, sg->nnodes, pathval);

  for (p = 0; p < sg->nnodes; p++) {
    pathval[ p ] = INT_MAX;

  }
   
  pathval[0]  = 0;
  sg->node[0].pred = NIL;
  InsertGQueue( &Q, 0 );

  // Prim's algorithm for Minimum Spanning Tree

  while ( !EmptyGQueue( Q ) ) {
    p  = RemoveGQueue( Q );

    for (q=0; q < sg->nnodes; q++){
      if (Q->L.elem[q].color!=BLACK){
	if (p!=q){
	  weight = ArcWeight(EuclDist(sg->node[p].feat,sg->node[q].feat,sg->nfeats)/sg->df);
	  if ( weight < pathval[ q ] ) {
	    if (pathval[q]!=INT_MAX)
	      RemoveGQueueElem(Q, q);
	    pathval[q]=weight;
	    sg->node[q].pred = p;
	    InsertGQueue( &Q, q );
	  }
	}
      }
    }
  }

  DestroyGQueue( &Q );
  free( pathval );

  // Prototype identification

  for (p=0; p < sg->nnodes; p++)     
    if (sg->node[p].pred!=NIL){
      q = sg->node[p].pred;
      if (sg->node[p].truelabel!=sg->node[q].truelabel){
	sg->node[p].status=sg->node[q].status=PROTOTYPE;
	printf("Prototypes %d %d\n",p,q);
      }
    }

}


// Compute Euclidean distance between feature vectors
 
float EuclDist(float *f1, float *f2, int n)
{
  int i;
  float dist=0; 

  for (i=0; i < n; i++)
    dist += (f1[i]-f2[i])*(f1[i]-f2[i]); 

  return(dist);
}

// Create adjacent list in subgraph: a knn graph 

void CreateArcs(Subgraph *sg, int knn)
{
  int    i,j,l,k;
  float  dist;
  int   *nn=AllocIntArray(knn+1);
  float *d=AllocFloatArray(knn+1);

  /* Create graph with the knn-nearest neighbors */

  sg->df=0;
  
  for (i=0; i < sg->nnodes; i++) {
    for (l=0; l < knn; l++) 
      d[l]=REAL_MAX;
      
    for (j=0; j < sg->nnodes; j++) {
      if (j!=i){
	d[knn] = EuclDist(sg->node[i].feat,sg->node[j].feat,sg->nfeats);
	nn[knn]= j;
	k      = knn;
	while ((k > 0)&&(d[k]<d[k-1])){
	  dist    = d[k];
	  l       = nn[k];
	  d[k]    = d[k-1];
	  nn[k]   = nn[k-1];
	  d[k-1]  = dist;
	  nn[k-1] = l;
	  k--;
	}
      }
    }
      
    for (l=0; l < knn; l++){
      if (d[l] > sg->df) sg->df = d[l];
      InsertSet(&(sg->node[i].adj),nn[l]);	
    }
  }
  
  free(d);
  free(nn);  
}

void CreateSymArcs(Subgraph *sg)
{
  int    i,j;
  float  dist;

  /* Create graph with neighbors closer than sg->df */

  for (i=0; i < sg->nnodes; i++) {
    for (j=0; j < sg->nnodes; j++) {
      if (i!=j){
	dist = EuclDist(sg->node[i].feat,sg->node[j].feat,sg->nfeats);
	if (dist <= sg->df)
	  InsertSet(&(sg->node[i].adj),j);	
      }
    }
  }
}

void DestroyArcs(Subgraph *sg)
{
  int i;

  for (i=0; i < sg->nnodes; i++) 
    DestroySet(&(sg->node[i].adj));
}

void PDF(Subgraph *sg)
{
  int     i,nelems;
  float   dist;
  double  value;
  Set    *adj=NULL;

  sg->K2    = (2.0*(double)sg->df/9.0);
  sg->K1    = ((double)1.0/sqrt(2.0*PI*(double)sg->df/9.0));

  sg->mindens = INT_MAX; sg->maxdens = INT_MIN;
  for (i=0; i < sg->nnodes; i++) {
    nelems=1;    
    adj=sg->node[i].adj;
    value=1.0;
    while(adj != NULL) {
      dist = EuclDist(sg->node[i].feat,sg->node[adj->elem].feat,sg->nfeats); 
      nelems++;
      value += sg->K1*exp(-dist/sg->K2);
      adj = adj->next;
    }
    sg->node[i].dens = (int)(MAXDENS*value / nelems) + 1;
    if (sg->node[i].dens < sg->mindens)
      sg->mindens = sg->node[i].dens;
    if (sg->node[i].dens > sg->maxdens)
      sg->maxdens = sg->node[i].dens;
  }
}

CNode *CreateCNode(Subgraph *sg)
{
  CNode *node=(CNode *)calloc(1,sizeof(CNode));

  node->feat = AllocFloatArray(sg->nfeats);
  node->adj  = AllocIntArray(sg->bestk+1);
  node->dist = AllocFloatArray(sg->bestk+1);

  return(node);
}

void DestroyCNode(CNode **node)
{
  if ((*node)!=NULL){
    free((*node)->feat);
    free((*node)->adj);
    free((*node)->dist);
    free((*node));
    *node = NULL;
  }
}

void MarkovNodeFeat3(Scene *scn, CNode *node, int p, float dm)
{
  int     q,j,xysize=scn->xsize*scn->ysize,aux;
  Voxel   u,v;
  AdjRel3 *A;

  node->position = p;

  if (dm==0.0)
    node->feat[0]=scn->data[p];
  else{
    A = Spheric(dm);
    aux = p%xysize;
    u.x = aux%scn->xsize;
    u.y = aux/scn->xsize;
    u.z = p/xysize;
    for (j=0; j < A->n; j++) {
      v.x = u.x + A->dx[j];
      v.y = u.y + A->dy[j];
      v.z = u.z + A->dz[j];
      if (ValidVoxel(scn,v.x,v.y,v.z)){
	q = v.x + scn->tby[v.y] + scn->tbz[v.z];
	node->feat[j]=scn->data[q];
      }else
	node->feat[j]=scn->data[p];
    }
    DestroyAdjRel3(&A);
  }
}

void NodeArcs(Subgraph *sg, CNode *node)
{
  int    j,l,k;
  float  dist;

  for (l=0; l <= sg->bestk; l++) 
    node->dist[l]=REAL_MAX;

  for (j=0; j < sg->nnodes; j++) {
    if (sg->node[j].pathval != INT_MIN){ // only nodes achieved by
					 // optimum paths during training
      node->dist[sg->bestk] = EuclDist(node->feat,sg->node[j].feat,sg->nfeats);
      node->adj[sg->bestk]  = j;
      k                     = sg->bestk;
      while ((k > 0)&&(node->dist[k]<node->dist[k-1])){
	dist            = node->dist[k];
	l               = node->adj[k];
	node->dist[k]   = node->dist[k-1];
	node->adj[k]     = node->adj[k-1];
	node->dist[k-1] = dist;
	node->adj[k-1]   = l;
	k--;
      }
    }   
  }
}

void NodePD(Subgraph *sg, CNode *node)
{
  int     i,nelems;
  float   dist;
  double  value;

  nelems=1;
  value=1.0;
  for (i=0; i < sg->bestk; i++) {
    dist = EuclDist(node->feat,sg->node[node->adj[i]].feat,sg->nfeats); 
    nelems++;
    value += sg->K1*exp(-dist/sg->K2);
  }
  node->dens = (int)(MAXDENS*value / nelems) + 1;
}

int PredBestPath(Subgraph *sg, CNode *node)
{
  int maxpathval, i, pred, tmp;
  
  maxpathval=MIN(node->dens,sg->node[node->adj[0]].pathval);;
  pred      =node->adj[0];
  for (i=1; i < sg->bestk; i++){
    tmp=MIN(node->dens,sg->node[node->adj[i]].pathval);
    if (tmp > maxpathval){
      pred = node->adj[i];
      maxpathval = tmp;
    }
  }

  return(pred);
}

void SetSubgraphFeatures(Subgraph *sg, Features *f)
{
  int i,j;

  sg->nfeats = f->nfeats;  
  for (i=0; i < sg->nnodes; i++){
    sg->node[i].feat = AllocFloatArray(sg->nfeats);
    for (j=0; j < sg->nfeats; j++) 
      sg->node[i].feat[j] = f->elem[sg->node[i].position].feat[j];
  }
}

void SetSubgraphFeatures3(Subgraph *sg, Features3 *f)
{
  int i,j;

  sg->nfeats = f->nfeats;  
  for (i=0; i < sg->nnodes; i++){
    sg->node[i].feat = AllocFloatArray(sg->nfeats);
    for (j=0; j < sg->nfeats; j++) 
      sg->node[i].feat[j] = f->elem[sg->node[i].position].feat[j];
  }
}

int ArcWeight(float dist)
{
  return((int)(MAXARCW*dist));
}

// Split subgraph into two parts such that the size of the first part
// is given by a percentual of samples.

void SplitSubgraph(Subgraph *sg, Subgraph **sg1, Subgraph **sg2, float perc1)
{
  int *label=AllocIntArray(sg->nlabels),i,j,i1,i2;
  int *nelems=AllocIntArray(sg->nlabels),totelems;

  for (i=0; i < sg->nnodes; i++) {
    label[sg->node[i].truelabel]++;
  }

  for (i=0; i < sg->nnodes; i++) {
    nelems[sg->node[i].truelabel]=(int)MAX((int)(perc1*label[sg->node[i].truelabel]),1);
  }

  free(label);

  totelems=0;
  for (j=0; j < sg->nlabels; j++) 
    totelems += nelems[j];

  *sg1 = CreateSubgraph(totelems);
  *sg2 = CreateSubgraph(sg->nnodes-totelems);
  (*sg1)->nfeats = sg->nfeats;
  (*sg2)->nfeats = sg->nfeats;

  for (i1=0; i1 < (*sg1)->nnodes; i1++) 
    (*sg1)->node[i1].feat = AllocFloatArray((*sg1)->nfeats);
  for (i2=0; i2 < (*sg2)->nnodes; i2++) 
    (*sg2)->node[i2].feat = AllocFloatArray((*sg2)->nfeats);

  (*sg1)->nlabels = sg->nlabels;
  (*sg2)->nlabels = sg->nlabels;

  i1=0;
  while(totelems > 0){
    i = RandomInteger(0,sg->nnodes-1);
    if (sg->node[i].status!=NIL){
      if (nelems[sg->node[i].truelabel]>0){ // copy node to sg1
	(*sg1)->node[i1].position = sg->node[i].position;
	for (j=0; j < (*sg1)->nfeats; j++) 
	  (*sg1)->node[i1].feat[j]=sg->node[i].feat[j];
	(*sg1)->node[i1].truelabel = sg->node[i].truelabel;
	i1++;
	nelems[sg->node[i].truelabel] = nelems[sg->node[i].truelabel] - 1;
	sg->node[i].status = NIL;
	totelems--;
      }
    }
  }

  i2=0;
  for (i=0; i < sg->nnodes; i++) {
    if (sg->node[i].status!=NIL){
      (*sg2)->node[i2].position = sg->node[i].position;
      for (j=0; j < (*sg2)->nfeats; j++) 
	(*sg2)->node[i2].feat[j]=sg->node[i].feat[j];
      (*sg2)->node[i2].truelabel = sg->node[i].truelabel;
      i2++;
    }
  }

  free(nelems);
}

// Move errors above max_abs_error from the evaluation set to the
// training set

int MoveErrorsToTrainingSet(Subgraph **sgtrain, Subgraph **sgeval, int max_abs_error)
{
  int i,j,k,l,nerrors;
  Subgraph *sgtrain_aux, *sgeval_aux;


  nerrors=0;
  for (i=0; i < (*sgeval)->nnodes; i++)
    if (abs((*sgeval)->node[i].label-(*sgeval)->node[i].truelabel)>max_abs_error)
      nerrors++;

  if (nerrors==0)
    return(0);

  sgtrain_aux = CreateSubgraph((*sgtrain)->nnodes+nerrors);
  sgeval_aux  = CreateSubgraph((*sgeval)->nnodes-nerrors);

  sgtrain_aux->nfeats=(*sgtrain)->nfeats;
  sgeval_aux->nfeats=(*sgeval)->nfeats;


  sgtrain_aux->nlabels = (*sgtrain)->nlabels;
  sgeval_aux->nlabels = (*sgeval)->nlabels;

  for (i=0; i < (*sgtrain)->nnodes; i++){ // Copy training set    
    sgtrain_aux->node[i].feat = AllocFloatArray(sgtrain_aux->nfeats);
    for (j=0; j < sgtrain_aux->nfeats; j++) 
      sgtrain_aux->node[i].feat[j] = (*sgtrain)->node[i].feat[j]; 
    sgtrain_aux->node[i].truelabel = (*sgtrain)->node[i].truelabel;
    sgtrain_aux->node[i].position = (*sgtrain)->node[i].position;
  }

  k=i;
  l=0;
  for (i=0; i < (*sgeval)->nnodes; i++){ // Copy errors 
    if (abs((*sgeval)->node[i].label-(*sgeval)->node[i].truelabel)>max_abs_error){
      sgtrain_aux->node[k].feat = AllocFloatArray(sgtrain_aux->nfeats);
      for (j=0; j < sgtrain_aux->nfeats; j++) 
	sgtrain_aux->node[k].feat[j] = (*sgeval)->node[i].feat[j]; 
      sgtrain_aux->node[k].truelabel = (*sgeval)->node[i].truelabel;
      sgtrain_aux->node[k].position  = (*sgeval)->node[i].position;
      k++;
    }else{ // Copy non-errors
      sgeval_aux->node[l].feat = AllocFloatArray(sgeval_aux->nfeats);
      for (j=0; j < sgtrain_aux->nfeats; j++) 
	sgeval_aux->node[l].feat[j] = (*sgeval)->node[i].feat[j]; 
      sgeval_aux->node[l].truelabel = (*sgeval)->node[i].truelabel;
      sgeval_aux->node[l].position  = (*sgeval)->node[i].position;
      l++;
    }
  }

  DestroySubgraph(sgtrain);
  DestroySubgraph(sgeval);

  *sgtrain = sgtrain_aux;
  *sgeval  = sgeval_aux;

  printf("Training set size is %d\n",(*sgtrain)->nnodes);

  return(nerrors);
}

// Remove nodes with distinct labels and distance=0

void RemoveInconsistentNodes(Subgraph **sg)
{
  Subgraph *newsg;
  int i,j,k,num_of_inconsistents=0;

  for (i=0; i < (*sg)->nnodes; i++) 
    (*sg)->node[i].status=0;

  for (i=0; i < (*sg)->nnodes; i++) 
    if ((*sg)->node[i].status!=NIL){ // consistent node      
      for (j=0; j < (*sg)->nnodes; j++) 
	if ((i!=j)&&((*sg)->node[j].status!=NIL)){ // consistent node
	  if ((*sg)->node[i].truelabel!=(*sg)->node[j].truelabel)
	    if (EuclDist((*sg)->node[i].feat,(*sg)->node[j].feat,(*sg)->nfeats)==0.0){// i is not consistent with j, then mark j as inconsistent node
	      (*sg)->node[j].status=NIL;
	      num_of_inconsistents++;
	    }
	}
    }

  if (num_of_inconsistents>0){
    printf("There are %d inconsistent nodes among %d nodes\n",num_of_inconsistents,(*sg)->nnodes);
    newsg = CreateSubgraph((*sg)->nnodes - num_of_inconsistents);
    newsg->nfeats = (*sg)->nfeats;
    for (i=0; i < newsg->nnodes; i++) 
      newsg->node[i].feat = AllocFloatArray(newsg->nfeats);

    k=0;
    newsg->nlabels = (*sg)->nlabels;
    for (i=0; i < (*sg)->nnodes; i++) {
      if ((*sg)->node[i].status!=NIL){ // consistent node
	newsg->node[k].position  = (*sg)->node[i].position;
	newsg->node[k].truelabel = (*sg)->node[i].truelabel;
	for (j=0; j < newsg->nfeats; j++) 
	  newsg->node[k].feat[j] = (*sg)->node[i].feat[j];
	k++;
      }
    }
    newsg->nlabels=(*sg)->nlabels;
    DestroySubgraph(sg);
    *sg=newsg;
  }
}

// Merge subgraphs with no arcs and feats

Subgraph *MergeSubgraphs(Subgraph *sg1, Subgraph *sg2)
{
  Subgraph *sg3;
  int i1,i2,i3,nnodes;

  nnodes=sg1->nnodes+sg2->nnodes;
  for (i1=0; i1 < sg1->nnodes; i1++)
    for (i2=0; i2 < sg2->nnodes; i2++) 
      if (sg1->node[i1].position==sg2->node[i2].position){
	nnodes--;
        if (sg1->node[i1].status==PROTOTYPE){
	  sg2->node[i2].status=NIL; //marked for removal
	}else{ 
	  if (sg2->node[i2].status==PROTOTYPE){
	    sg1->node[i1].status=NIL; //marked for removal
	  }else{
	    sg2->node[i2].status=NIL; //marked for removal
	  }
	}
      }

  sg3 = CreateSubgraph(nnodes);
  i3  = 0;
  for (i1=0; i1 < sg1->nnodes; i1++) 
    if (sg1->node[i1].status!=NIL){
      sg3->node[i3].position  = sg1->node[i1].position; 
      sg3->node[i3].status    = sg1->node[i1].status;
      sg3->node[i3].truelabel = sg1->node[i1].truelabel;
      sg3->node[i3].adj=NULL;
      sg3->node[i3].feat=NULL;
      i3++;
    }
  for (i2=0; i2 < sg2->nnodes; i2++) 
    if (sg2->node[i2].status!=NIL){
      sg3->node[i3].position  = sg2->node[i2].position; 
      sg3->node[i3].status    = sg2->node[i2].status;
      sg3->node[i3].truelabel = sg2->node[i2].truelabel;
      sg3->node[i3].adj=NULL;
      sg3->node[i3].feat=NULL;
      i3++;
    }
  return(sg3);
}
