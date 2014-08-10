#include "complgraph.h"


//Todas funcoes abaixo assumem que os nós do
//grafo sao armazenados de forma ordenada pelo label. 


//Embaralha os nós de uma mesma classe de forma 
//aleatoria, preservando a ordenacao dos labels.
void RandomizeComplGraph(ComplGraph *cg){
  int i,j,c,offset=0;
  
  srand((int)time(NULL));
  for(c=1; c<=cg->nclasses; c++){
    for(i=0; i<cg->nclass[c]; i++){
      j = RandomInteger(0, cg->nclass[c]-1);
      GraphNodeSwap(&cg->node[offset+j], 
		    &cg->node[offset+i]);
    }
    offset += cg->nclass[c];
  }
}


//Retorna um subconjunto do grafo "cg" com amostras 
//de todas classes na proporcao dada por "rate".
//O subconjunto retornado eh removido do grafo "cg"
//original.
//A funcao "RandomizeComplGraph" deve ser chamada
//previamente para garantir amostras aleatorias.
ComplGraph *RemoveSamplesFromComplGraph(ComplGraph **cg, 
					float rate){
  ComplGraph *cg0,*cg1,*cg2;
  int p1,p2,j,i,c;
  int nsamples,ndata;
  int *tmp;

  cg0 = *cg;
  ndata = cg0->node[0].ndata;
  nsamples = 0;
  tmp = (int *)calloc(cg0->nclasses+1, sizeof(int));
  for(c=1; c<=cg0->nclasses; c++){
    tmp[c] = ceil(cg0->nclass[c]*rate);
    nsamples += tmp[c];
  }

  cg1 = CreateComplGraph(cg0->n-nsamples, cg0->nclasses, ndata);
  cg2 = CreateComplGraph(nsamples,        cg0->nclasses, ndata);

  for(c=1; c<=cg0->nclasses; c++){
    cg1->nclass[c] = cg0->nclass[c] - tmp[c];
    cg2->nclass[c] = tmp[c];
  }
  free(tmp);

  p1 = p2 = 0;
  for(c=1; c<=cg0->nclasses; c++){
    for(i=0; i<cg1->nclass[c]; i++){
      GraphNodeSwap(&cg0->node[p1+p2],&cg1->node[p1]);
      p1++;
    }
    for(j=0; j<cg2->nclass[c]; j++){
      GraphNodeSwap(&cg0->node[p1+p2],&cg2->node[p2]);
      p2++;
    }
  }

  //Troca cg e cg1.
  *cg = cg1;
  DestroyComplGraph(&cg0);

  return cg2;
}


//Junta dois complete graphs.
//A ordenacao dos labels eh preservada.
ComplGraph *MergeComplGraph(ComplGraph *cg1, 
			    ComplGraph *cg2){
  ComplGraph *cg;
  int p1,p2,j,i,c;
  int ndata,nclasses,nclass;

  if(cg1==NULL && cg2==NULL) return NULL;
  if(cg1==NULL)
    return CloneComplGraph(cg2);
  if(cg2==NULL)
    return CloneComplGraph(cg1);

  ndata = cg1->node[0].ndata;
  nclasses = MAX(cg1->nclasses, cg2->nclasses);
  cg = CreateComplGraph(cg1->n+cg2->n, nclasses, ndata);

  for(c=1; c<=nclasses; c++){
    nclass = 0;
    if(c<=cg1->nclasses)
      nclass += cg1->nclass[c];
    if(c<=cg2->nclasses)
      nclass += cg2->nclass[c];
    cg->nclass[c] = nclass;
  }

  p1 = p2 = 0;
  for(c=1; c<=nclasses; c++){
    if(c<=cg1->nclasses){
      for(i=0; i<cg1->nclass[c]; i++){
	GraphNodeCopy(&cg->node[p1+p2],
		      &cg1->node[p1]);
	p1++;
      }
    }
    if(c<=cg2->nclasses){
      for(j=0; j<cg2->nclass[c]; j++){
	GraphNodeCopy(&cg->node[p1+p2],
		      &cg2->node[p2]);
	p2++;
      }
    }
  }

  return cg;
}


//Testa a base de treinamento com a base de teste e 
//troca os nos classificados incorretamente com nos 
//nao semente da base de treinamento.
//A funcao "RandomizeComplGraph" deve ser chamada
//previamente para garantir trocas aleatorias.
// @return Taxa de acerto.
float ComplGraphRetraining(ComplGraph *cgTraining, 
			   ComplGraph *cgTesting){
  float rate;
  char msg[512];
  int i, nerrados = 0, nrelevant = 0;
  BMap *bm_error;

  if(cgTraining->nclasses != cgTesting->nclasses){
    sprintf(msg,"Incompatible number of classes");
    Error(msg,"ComplGraphRetraining");
  }

  bm_error = BMapNew(cgTesting->n);
  rate = ComplexClassifyGraph(cgTraining, cgTesting, bm_error);
  IdentifyNonRelevants(cgTraining);
  SwapErrorbyNotRelevant(cgTraining, cgTesting, bm_error);

  //verificando se sobraram elementos errados no teste 
  //pra serem classificados
  for(i=0; i<cgTesting->n; i++)
    if(_fast_BMapGet(bm_error,i))
      nerrados++;
  printf("\nelementos errados restantes: %d",nerrados);
  
  if(nerrados > 0) //se ainda sobrarem elementos errados 
    SwapErrorbyRelevant(cgTraining, cgTesting, bm_error);
	
  //verificando quantos nao relevantes sobraram
  for(i = 0; i < cgTraining->n; i++) 
    if(!cgTraining->node[i].is_relevant)
      nrelevant++;
  printf("\nsobraram %d nao relevantes",nrelevant);
  
  if(nrelevant > 0)
    SwapNotRelevantbyNewSamples(cgTraining, cgTesting, bm_error);
  
  nrelevant = 0;
  //verificando quantos nao relevantes sobraram
  for(i = 0; i < cgTraining->n; i++) 
    if(!cgTraining->node[i].is_relevant)
      nrelevant++;
  printf("\nsobraram %d nao relevantes",nrelevant);
  
  BMapDestroy(bm_error);
  
  for(i=0; i<cgTraining->n; i++){
    cgTraining->node[i].is_seed = 0;
    cgTraining->node[i].is_relevant = 1;
  }
  
  return rate;
}


//Grava o "ComplGraph" para disco em formato binario
//onde os vetores de caracteristicas sao gravados de 
//forma sequencial. O cabecalho armazena o tamanho
//dos vetores, bem como o numero de classes e o numero
//de registros por classe.
void  WriteComplGraph(ComplGraph *cg, char *filename){
  char msg[512];
  FILE *fp;
  int i,c,ndata;

  fp = fopen(filename,"wb");
  if(fp == NULL){
    sprintf(msg,"Cannot save %s",filename);
    Error(msg,"WriteComplGraph");
  }

  ndata = cg->node[0].ndata;
  fwrite(&cg->n,        sizeof(int), 1, fp);
  fwrite(&cg->nclasses, sizeof(int), 1, fp);
  fwrite(&ndata,        sizeof(int), 1, fp);

  for(c=1; c<=cg->nclasses; c++)
    fwrite(&(cg->nclass[c]), sizeof(int), 1, fp);
  
  for(i=0; i<cg->n; i++)
    fwrite(cg->node[i].data, sizeof(real), ndata, fp);

  fclose(fp);
}


//Le o "ComplGraph" a partir de um arquivo na forma 
//binaria.
ComplGraph *ReadComplGraph(char *filename){
  ComplGraph *cg;
  char msg[512];
  FILE *fp;
  int n,nclasses,ndata,i,j,c,m;

  fp = fopen(filename,"rb");
  if(fp == NULL){
    sprintf(msg,"Cannot open %s",filename);
    Error(msg,"ReadComplGraph");
  }

  fread(&n,sizeof(int), 1, fp);
  fread(&nclasses,sizeof(int), 1, fp);
  fread(&ndata,sizeof(int), 1, fp);

  cg = CreateComplGraph(n, nclasses, ndata);

  m = 0;
  for(c=1; c<=nclasses; c++){
    fread(&(cg->nclass[c]), sizeof(int), 1, fp);
    m += cg->nclass[c];
  }
  
  if(n!=m){
    sprintf(msg,"Bad or corrupted file");
    Error(msg,"ReadComplGraph");
  }
  
  c = 1;
  j = 0;
  for(i=0; i<n; i++){
    cg->node[i].label = c;
    fread(cg->node[i].data, sizeof(real), ndata, fp);
    j++;
    if(j>=cg->nclass[c]){
      c++;
      j = 0;
    }
  }
  fclose(fp);

  return cg;
}

inline void GraphNodeCopy(GraphNode *dest, GraphNode *src){
  memcpy(dest->data,
	 src->data, 
	 src->ndata*sizeof(real));
  dest->ndata = src->ndata;
  dest->label = src->label;
  dest->pred  = src->pred;
  dest->cost  = src->cost;
  dest->is_seed = src->is_seed;
  dest->relevance_value = src->relevance_value;
  dest->is_relevant = src->is_relevant;
}

inline real GraphNodeDistance(GraphNode *a, GraphNode *b){
  real sum = 0.0f;
  int i;

  for(i=0; i<a->ndata; i++)
    sum += (b->data[i]-a->data[i])*(b->data[i]-a->data[i]);
  return (real)sqrtf(sum);
}


inline void GraphNodeSwap(GraphNode *a, GraphNode *b){
  GraphNode tmp;
  
  tmp = *a;
  *a = *b;
  *b = tmp;
}


ComplGraph *CreateComplGraph(int n, 
			     int nclasses, 
			     int ndata){
  ComplGraph *cg;
  int i;

  cg = (ComplGraph *) malloc(sizeof(ComplGraph));
  if(cg == NULL)
    Error(MSG1,"CreateComplGraph");

  cg->n = n;
  cg->node = (GraphNode *) malloc(sizeof(GraphNode)*n);
  if(cg->node == NULL)
    Error(MSG1,"CreateComplGraph");

  cg->nclasses = nclasses;
  cg->nclass = (int *) calloc(nclasses+1, sizeof(int));
  if(cg->nclass == NULL)
    Error(MSG1,"CreateComplGraph");

  for(i=0; i<cg->n; i++){
    cg->node[i].data = (real *) malloc(sizeof(real)*ndata);
    cg->node[i].ndata = ndata;
    cg->node[i].pred = NIL;
    cg->node[i].cost = REAL_MAX;
    cg->node[i].is_seed = 0;
    cg->node[i].relevance_value = 0;
    cg->node[i].is_relevant = 1;
  }

  return cg;
}


void        DestroyComplGraph(ComplGraph **graph){
  ComplGraph *cg = *graph;
  int i;

  if(cg == NULL)
    return;
  for(i=0; i<cg->n; i++)
    free(cg->node[i].data);
  free(cg->node);
  free(cg->nclass);
  free(cg);
  *graph = NULL;
}


ComplGraph *CloneComplGraph(ComplGraph *cg){
  ComplGraph *clone;
  int ndata,i,c;

  ndata = cg->node[0].ndata;
  clone = CreateComplGraph(cg->n, cg->nclasses, ndata);

  for(i=0; i<cg->n; i++){
    GraphNodeCopy(&clone->node[i],
		  &cg->node[i]);
  }

  for(c=1; c<=cg->nclasses; c++)
    clone->nclass[c] = cg->nclass[c];

  return clone;
}


/* Calcula a árvore geradora mínima de um 
   Complete Graph.  */
void        ComplGraphMST(ComplGraph *cg){
  RealHeap *H;
  real *cost, dist;
  int i,p,q;
  
  /* vetor de real com os custos. */
  cost = (real *) malloc(sizeof(real)*cg->n);

  for(i=0; i<cg->n; i++)
    cost[i] = REAL_MAX;
  
  H = CreateRealHeap(cg->n+1, cost);

  cost[0] = 0.0; /* começa pelo nó 0 */
  InsertRealHeap(H, 0);

  while(!IsEmptyRealHeap(H)) {
    RemoveRealHeap(H, &p);
    /* todos os nós do grafo completo */
    for(q=0; q<cg->n; q++){ 
      if(p == q || H->color[q] == BLACK) continue;
      dist = GraphNodeDistance(&(cg->node[p]), 
			       &(cg->node[q]));
      if(dist < cost[q]){
	cg->node[q].pred = p;
	cost[q] = dist;
	if(H->color[q] == WHITE)
	  InsertRealHeap(H, q);
	else
	  GoUpRealHeap(H, H->pos[q]);
      }
    }
  }
  DestroyRealHeap(&H);

  /* atualiza os custos dos nós na estrutura 
     ComplGraph. */
  for(i=0; i<cg->n; i++)
    cg->node[i].cost = cost[i];
  
  free(cost);
}


/* Treina um ComplGraph. É calculada a árvore 
   geradora mínima, e os nós pai e filho com 
   labels diferentes viram sementes de um watershed. */
void        ComplGraphTraining(ComplGraph *cg){
  RealHeap *H;
  real *cost, dist, aux;
  int i,j,p,q;
  
  ComplGraphMST(cg);
  printf("MST generated.\n");
  
  /* vetor de real com os custos. */
  cost = (real *) malloc(sizeof(real)*cg->n);

  for(i=0; i<cg->n; i++){
    cost[i] = REAL_MAX;

    /* se a classe do nó i for diferente da classe 
       de seu predecessor na MST, marcar os dois nós 
       como sementes */
    j = cg->node[i].pred;
    if(j != NIL && 
       cg->node[i].label != cg->node[j].label){
      cg->node[i].is_seed = 1;
      cg->node[j].is_seed = 1;
    }
  }
  
  H = CreateRealHeap(cg->n+1, cost);

  /* insere as sementes na fila */
  j = 0;
  for(i=0; i<cg->n; i++){
    if(cg->node[i].is_seed){
      j++;
      cost[i] = 0.0;
      InsertRealHeap(H, i);
    }
  }
  printf("Número de sementes: %d\n", j);

  while(!IsEmptyRealHeap(H)){
    RemoveRealHeap(H, &p);
    /* todos os nós do grafo completo */
    for(q=0; q<cg->n; q++){
      if(p == q || cost[p] >= cost[q]) continue;
      aux = GraphNodeDistance(&(cg->node[p]),
			      &(cg->node[q]));
      dist = MAX(cost[p], aux);
      if(dist < cost[q]){
	cg->node[q].pred = p;
	cg->node[q].label = cg->node[p].label;
	cost[q] = dist;
	if (H->color[q] == WHITE)
	  InsertRealHeap(H, q);
	else
	  GoUpRealHeap(H, H->pos[q]);
      }
    }
  }
  DestroyRealHeap(&H);
  
  /* atualiza os custos dos nós na estrutura 
     ComplGraph. */
  for(i=0; i<cg->n; i++)
    cg->node[i].cost = cost[i];
  
  free(cost);
}



int         ComplGraphTestNode(ComplGraph *cg, 
			       GraphNode *node, int *p){
  int i, j, iMin, nclasses;
  int *label; /* numero de labels diferentes */
  int *index;
  real *cost, aux, minCost;
  
  nclasses = cg->nclasses;
  minCost = REAL_MAX;
  cost = (real *) malloc(sizeof(real)*cg->n);
  for(i=0; i<cg->n; i++){
    aux = GraphNodeDistance(&(cg->node[i]), 
			    node);
    cost[i] = MAX(cg->node[i].cost, aux);
    if(cost[i] < minCost)
      minCost = cost[i];
  }

  label = (int *) calloc(nclasses+1, sizeof(int));
  index = (int *) calloc(nclasses+1, sizeof(int));
  
  for(i=0; i<cg->n; i++){
    if(cost[i] == minCost){
      label[cg->node[i].label]++;
      index[cg->node[i].label] = i;
    }
  }
  
  free(cost);
  
  iMin = 0;
  j = 0;
  for(i=1; i<nclasses+1; i++){
    if(label[i] > j){
      j = label[i];
      iMin = i;
    }
  }

  *p = index[iMin];
  
  free(label);
  free(index);
  
  if(iMin == 0)
    iMin = 1;
  
  return iMin;
}


/* Calcula capacidade de treinamento de cada classe */
void        ComputeLearningCapacity(ComplGraph *cg){
  float *LC = NULL;
  int i,c;

  LC = (float *)calloc(cg->nclasses+1, sizeof(float));
  
  for(i=0; i<cg->n; i++){
    //storing seeds number
    if(cg->node[i].is_seed)
      LC[cg->node[i].label]++; 
  }

  //calculating LC for each class i
  for(c=1; c<=cg->nclasses; c++)
    LC[c] = 1.0-(LC[c]/cg->nclass[c]);

  for(c=1; c<=cg->nclasses; c++)
    printf("\nLC[%d]: %f",c,LC[c]);
  
  free(LC);
}


void SetRelevanceValuePath(ComplGraph *cg, 
			   int nodeID, int value){
  while(!cg->node[nodeID].is_seed){
    cg->node[nodeID].relevance_value += value; 
    nodeID = cg->node[nodeID].pred;
  }
  //setando relevance_value da raiz
  cg->node[nodeID].relevance_value += value; 
}


float ComplexClassifyGraph(ComplGraph *cgTraining, 
			   ComplGraph *cgTesting,
			   BMap *bm_error){
  int nclasses, result, c, i, nerrors, nodeID;
  float Acc, **error_matrix, error;

  nclasses = cgTraining->nclasses;
  for(i=0; i<cgTraining->n; i++)
    cgTraining->node[i].relevance_value = 0;

  //matriz de erros ei = ei1+ei2
  //ei1 = FP    ei2 = FN
  error_matrix = (float **)calloc(nclasses+1, sizeof(float *));
  for(c=1; c<=nclasses; c++)
    error_matrix[c] = (float *)calloc(2, sizeof(float));

  nerrors = 0;
  for(i=0; i<cgTesting->n; i++){
    result = ComplGraphTestNode(cgTraining,
				&(cgTesting->node[i]), 
				&nodeID);
    
    if(result != cgTesting->node[i].label){
      _fast_BMapSet1(bm_error, i);
      nerrors++;
      //atualizando FP
      error_matrix[result][0]++; 
      //atualizando FN
      error_matrix[cgTesting->node[i].label][1]++;
      //-1 no caminho
      SetRelevanceValuePath(cgTraining, nodeID, -1); 
    }
    else{
      //+1 no caminho
      SetRelevanceValuePath(cgTraining, nodeID, 1);
    }
  }

  for(c=1; c<=nclasses; c++){
    error_matrix[c][0] /= (float)(cgTesting->n - cgTesting->nclass[c]);
    error_matrix[c][1] /= (float)cgTesting->nclass[c];
  }
	
  error = 0.0;
  for(c=1; c<=nclasses; c++)
    error += (error_matrix[c][0]+error_matrix[c][1]);
  
  Acc = 1.0-(error/(2.0*nclasses));//Acurácia

  for(c=1; c<=nclasses; c++)
    free(error_matrix[c]);
  free(error_matrix);

  printf("\nnerros: %d",nerrors);
  return Acc;
}


float SimpleClassifyGraph(ComplGraph *cgTraining, 
			  ComplGraph *cgTesting){
  int i,j,c,numnseeds,result,nerrors,nclasses;
  int *offset=NULL,*end=NULL, nodeID;
  float rate,**error_matrix,error;
  char msg[512];
  bool flag;
  BMap *bm_error;

  if(cgTraining->nclasses != cgTesting->nclasses){
    sprintf(msg,"Incompatible number of classes");
    Error(msg,"ComplGraphRetraining");
  }
  nclasses = cgTraining->nclasses;
	
  //matriz de erros ei = ei1+ei2
  //ei1 = FP    ei2 = FN
  error_matrix = (float **)calloc(nclasses+1, sizeof(float *));
  for(c=1; c<=nclasses; c++)
    error_matrix[c] = (float *)calloc(2, sizeof(float));

  bm_error = BMapNew(cgTesting->n);
  
  nerrors = 0;
  for(i=0; i<cgTesting->n; i++){
    result = ComplGraphTestNode(cgTraining, 
				&(cgTesting->node[i]), &nodeID);
    
    if(result != cgTesting->node[i].label){
      _fast_BMapSet1(bm_error, i);
      nerrors++;
			
      error_matrix[result][0]++; //atualizando FP
      error_matrix[cgTesting->node[i].label][1]++; //atualizando FN
    }
  }
  
  for(c=1; c<=nclasses; c++){
    error_matrix[c][0] /= (float)(cgTesting->n -
				  cgTesting->nclass[c]);
    error_matrix[c][1] /= (float)cgTesting->nclass[c];
  }
	
  error = 0.0;
  for(c=1; c<=nclasses; c++)
    error += (error_matrix[c][0]+error_matrix[c][1]);
  
  printf("\nnerrors: %d\n",nerrors);


  rate = 1.0-(error/(2.0*nclasses));//Learning Rate
  ComputeLearningCapacity(cgTraining);
  
  offset = (int *)calloc(nclasses+1, sizeof(int));
  end    = (int *)calloc(nclasses+1, sizeof(int));
  for(c=2; c<=nclasses; c++)
    offset[c] = offset[c-1] + cgTraining->nclass[c-1];

  for(c=1; c<=nclasses; c++)
    end[c] = offset[c] + cgTraining->nclass[c];
  
  numnseeds = 0;
  for(i=0; i<cgTraining->n; i++)
    if(!cgTraining->node[i].is_seed)
      numnseeds++;

  for(i=0; i<cgTesting->n; i++){
    if(_fast_BMapGet(bm_error,i)){
      c = cgTesting->node[i].label;
      flag = false;
      for(j=offset[c]; j<end[c]; j++){
	if(!cgTraining->node[j].is_seed){
	  flag = true;
	  break;
	}
      }
      offset[c] = j;
      
      if(flag){
	// faz a troca SOMENTE ELEMENTOS DE MESMA CLASSE
	GraphNodeSwap(&cgTraining->node[j], 
		      &cgTesting->node[i]);
	cgTraining->node[j].is_seed = 1;
	numnseeds--;
	if(numnseeds == 0)
	  break;

	offset[c]++;
      }
    }
  }
  
  for(i=0; i<cgTraining->n; i++)
    cgTraining->node[i].is_seed = 0;
  
  BMapDestroy(bm_error);
  for(c=1; c<=nclasses; c++)
    free(error_matrix[c]);
  free(error_matrix);
  free(offset);
  free(end);
  
  return rate;
}

void IdentifyNonRelevants(ComplGraph *cg){
  int i;
  
  for(i = 0; i<cg->n; i++)
    if(cg->node[i].relevance_value <= 0)
      cg->node[i].is_relevant = 0;
}


void SwapErrorbyNotRelevant(ComplGraph *cgTraining, 
			    ComplGraph *cgTesting,
			    BMap *bm_error){
  int *offset=NULL,*end=NULL;
  int i,j,c,numnrelevants,nclasses, contador1 = 0, contador2 = 0;;
  bool flag;

  nclasses = cgTraining->nclasses;

  offset = (int *)calloc(nclasses+1, sizeof(int));
  end    = (int *)calloc(nclasses+1, sizeof(int));
  for(c=2; c<=nclasses; c++)
    offset[c] = offset[c-1] + cgTraining->nclass[c-1];
  
  for(c=1; c<=nclasses; c++)
    end[c] = offset[c] + cgTraining->nclass[c];

  numnrelevants = 0;
  for(i=0; i<cgTraining->n; i++)
    if(!cgTraining->node[i].is_relevant)
      numnrelevants++;

  contador2 = numnrelevants;
  for(i=0; i<cgTesting->n; i++){
    if(numnrelevants == 0)
      break;
    if(_fast_BMapGet(bm_error,i) &&
       cgTesting->node[i].is_relevant){
      c = cgTesting->node[i].label;
      flag = false;
      for(j=offset[c]; j<end[c]; j++){
	if(!cgTraining->node[j].is_relevant){
	  flag = true;
	  break;
	}
      }
      
      offset[c] = j;
      
      if(flag){
	contador1++;
	// faz a troca SOMENTE ELEMENTOS DE MESMA CLASSE
	GraphNodeSwap(&cgTraining->node[j], &cgTesting->node[i]);
	//marca o elemento i como nao errado, pois na proxima 
	//etapa iremos verificar se aindam existem elementos 
	//errados para serem trocados
	_fast_BMapSet0(bm_error, i); 
	//marcando elemento do teste como nao relevante, para 
	//ele nao voltar no treinamento na funcao 
	//SwapNotRelevantbyNewSamples
	cgTesting->node[i].is_relevant = 0; 

	cgTraining->node[j].is_relevant = 1;
	numnrelevants--;
	
	offset[c]++;
      }
    }
  }
  
  free(offset);
  free(end);
  
  printf("\nelementos errados trocados por nao relevantes: %d",contador1);
  printf("\nnao relevantes: %d",contador2);
}


void SwapErrorbyRelevant(ComplGraph *cgTraining, 
			 ComplGraph *cgTesting,
			 BMap *bm_error){
  int *offset=NULL,*end=NULL;
  int i,j,c,numrelevants,nclasses, contador1 = 0, contador2 = 0;;
  bool flag;

  nclasses = cgTraining->nclasses;
  
  offset = (int *)calloc(nclasses+1, sizeof(int));
  end    = (int *)calloc(nclasses+1, sizeof(int));
  for(c=2; c<=nclasses; c++)
    offset[c] = offset[c-1] + cgTraining->nclass[c-1];

  for(c=1; c<=nclasses; c++)
    end[c] = offset[c] + cgTraining->nclass[c];

  numrelevants = 0;
  for(i=0; i<cgTraining->n; i++)
    if(cgTraining->node[i].is_relevant &&
       !cgTraining->node[i].is_seed)
      numrelevants++;
  
  contador2 = numrelevants;
  for(i=0; i<cgTesting->n; i++){
    if(numrelevants == 0)
      break;
    if(_fast_BMapGet(bm_error,i) &&
       cgTesting->node[i].is_relevant){
      c = cgTesting->node[i].label;
      flag = false;
      for(j=offset[c]; j<end[c]; j++){
	if(cgTraining->node[j].is_relevant && 
	   !cgTraining->node[j].is_seed){
	  flag = true;
	  break;
	}
      }
      
      offset[c] = j;
      
      if(flag){
	contador1++;
	// faz a troca SOMENTE ELEMENTOS DE MESMA CLASSE
	GraphNodeSwap(&cgTraining->node[j], &cgTesting->node[i]);
	//marca o elemento i como nao errado, pois na proxima 
	//etapa iremos verificar se aindam existem elementos 
	//errados para serem trocados
	_fast_BMapSet0(bm_error, i); 

	cgTesting->node[i].is_relevant = 1; 
	cgTraining->node[j].is_relevant = 1;
	cgTraining->node[j].is_seed = 0;
	numrelevants--;
	
	offset[c]++;
      }
    }
  }
  
  free(offset);
  free(end);

  printf("\nelementos errados trocados por relevantes nao prototipos: %d",contador1);
  printf("\nrelevantes: %d",contador2);
}


void SwapNotRelevantbyNewSamples(ComplGraph *cgTraining,
				 ComplGraph * cgTesting,
				 BMap *bm_error){
  int *offset=NULL,*end=NULL;
  int i,j,c,numnrelevants,nclasses, contador1 = 0, contador2 = 0;
  bool flag;

  nclasses = cgTraining->nclasses;

  offset = (int *)calloc(nclasses+1, sizeof(int));
  end    = (int *)calloc(nclasses+1, sizeof(int));
  for(c=2; c<=nclasses; c++)
    offset[c] = offset[c-1] + cgTraining->nclass[c-1];

  for(c=1; c<=nclasses; c++)
    end[c] = offset[c] + cgTraining->nclass[c];

  numnrelevants = 0;
  for(i=0; i<cgTraining->n; i++)
    if(!cgTraining->node[i].is_relevant)
      numnrelevants++;

  contador2 = numnrelevants;
  for(i=0; i<cgTesting->n; i++){
    if(numnrelevants == 0)
      break;
    if(cgTesting->node[i].is_relevant){
      c = cgTesting->node[i].label;
      flag = false;
      for(j=offset[c]; j<end[c]; j++){
	if(!cgTraining->node[j].is_relevant){
	  flag = true;
	  break;
	}
      }
      
      offset[c] = j;
      
      if(flag){
	contador1++;
	// faz a troca SOMENTE ELEMENTOS DE MESMA CLASSE
	GraphNodeSwap(&cgTraining->node[j], &cgTesting->node[i]);
	//marca o elemento i como nao errado, pois na proxima 
	//etapa iremos verificar se aindam existem elementos 
	//errados para serem trocados
	_fast_BMapSet0(bm_error, i); 

	cgTesting->node[i].is_relevant = 0;
	cgTraining->node[j].is_relevant = 1;
	numnrelevants--;
	
	offset[c]++;
      }
    }
  }
  
  free(offset);
  free(end);

  printf("\nelementos trocados por nao relevantes: %d",contador1);
  printf("\nnao relevantes: %d",contador2);
}


ComplGraph *OPF(char *descriptorFileName, 
		char *accFileName, int iterations){
  ComplGraph *cgTraining = NULL, *cgTesting = NULL;
  float Acc = 0.0, totaltime, AccRejected = 0.0;
  int i;
  FILE *fp;
  timer tic, toc;
  ComplGraph *cgRejected = NULL, *cgMerged = NULL;

  fp = fopen(accFileName,"w");

  gettimeofday(&tic,NULL);

  /*rodando algoritmo de treinamento*/
  for (i = 1; i<= iterations; i++){
    if(i == 1){
      cgTraining = ReadComplGraph(descriptorFileName);
      RandomizeComplGraph(cgTraining);
      cgTesting = RemoveSamplesFromComplGraph(&cgTraining, 0.5);
    }
    ComplGraphTraining(cgTraining);
    Acc = ComplGraphRetraining(cgTraining, cgTesting);

    fprintf(fp,"%d	%f\n",i,Acc);
  }

  //reduzindo treinamento retirando nao relevantes
  cgRejected = RemoveNotRelevants(&cgTraining);

  //projetando classificador final
  ComplGraphTraining(cgTraining);

  //computando acurácia levando em consideracao os rejeitados
  if(cgRejected->n != 0){
    cgMerged = MergeComplGraph(cgRejected, cgTesting);
    AccRejected = SimpleClassifyGraph(cgTraining, cgMerged);
    printf("\nAccRejected = %f",AccRejected);
  }

  DestroyComplGraph(&cgTesting);
  fclose(fp);
  
  printf("\nAcc = %f",Acc);
  
  gettimeofday(&toc,NULL);
  totaltime = (toc.tv_sec-tic.tv_sec)*1000.0 + (toc.tv_usec-tic.tv_usec)*0.001;
  printf("\nTime: %f ms\n",totaltime);

  return cgTraining;
}


ComplGraph * RemoveNotRelevants(ComplGraph **cg){
  int numrelevant = 0, ndata, i, j, k;
  ComplGraph *cg0, *cg1, *cg2;
  cg0 = *cg;

  //verificando o numero de elementos relevantes, 
  //que sera o tamanho do novo grafo
  for(i = 0; i < cg0->n; i++)
    if(cg0->node[i].is_relevant)
      numrelevant++;

  ndata = cg0->node[0].ndata;
  cg1 = CreateComplGraph(numrelevant, cg0->nclasses, ndata);
  cg2 = CreateComplGraph(cg0->n-numrelevant, cg0->nclasses, ndata);

  //copiando elementos relevantes de cg0 para cg1
  j = 0;
  k = 0;
  for (i = 0; i < cg0->n; i++){
    if(cg0->node[i].is_relevant){
      GraphNodeCopy(&(cg1->node[j]), &(cg0->node[i]));
      cg1->nclass[cg1->node[j].label]++;
      j++;
    }
    else{ //copia para o grafo de rejeitados
      GraphNodeCopy(&(cg2->node[k]), &(cg0->node[i]));
      cg2->nclass[cg2->node[k].label]++;
      k++;
    }
  }

  DestroyComplGraph(&cg0);
  *cg = cg1;

  return cg2;
}


void ShowPred(ComplGraph *cgTraining, int nodeID){
  printf("\nnodeID: %d-----------------------------",nodeID);
  while(!cgTraining->node[nodeID].is_seed){
    printf("\nnodeId: %d",nodeID);
    nodeID = cgTraining->node[nodeID].pred;
  }
  printf("\nnodeId: %d (seed)",nodeID);
}

