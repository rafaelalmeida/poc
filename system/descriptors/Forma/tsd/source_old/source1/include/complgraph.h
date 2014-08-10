
#ifndef _COMPLGRAPH_H_
#define _COMPLGRAPH_H_

#include <values.h>
#include <time.h>
#include "common.h"
#include "curve.h"
#include "comptime.h"
#include "realheap.h"
#include "scene.h"
#include "set.h"

typedef struct _graphnode {
  real *data;
  int   ndata;

  int   label;
  int   pred;
  real  cost;
  int relevance_value;
  BMap *flag; 
  /*bitmap de flags (0 - FALSE 1 - TRUE)
  posicao	significado
		0			is_seed
		1			is_relevant
		2			is_error
		3			is_outlier
*/
} GraphNode;


typedef struct _complgraph {
  int n;
  GraphNode *node;

  int nclasses;
  int *nclass;

} ComplGraph;

//Todas funcoes abaixo assumem que os nós do
//grafo sao armazenados de forma ordenada pelo label. 

/*Reseta flags de um grafo*/
void ResetFlagComplGraph(ComplGraph *cg);

//Embaralha os nós de uma mesma classe de forma 
//aleatoria, preservando a ordenacao dos labels.
void RandomizeComplGraph(ComplGraph *cg);

//Retorna um subconjunto do grafo "cg" com amostras 
//de todas classes na proporcao dada por "rate".
//O subconjunto retornado eh removido do grafo "cg"
//original.
//A funcao "RandomizeComplGraph" deve ser chamada
//previamente para garantir amostras aleatorias.
ComplGraph *RemoveSamplesFromComplGraph(ComplGraph **cg, 
					float rate);

//Junta dois complete graphs.
//A ordenacao dos labels eh preservada.
ComplGraph *MergeComplGraph(ComplGraph *cg1, 
			    ComplGraph *cg2);


//Testa a base de treinamento com a base de teste e 
//troca os nos classificados incorretamente com nos 
//nao semente da base de treinamento.
//A funcao "RandomizeComplGraph" deve ser chamada
//previamente para garantir trocas aleatorias.
// @return Taxa de acerto.
float ComplGraphRetraining(ComplGraph *cgTraining, 
			   ComplGraph *cgTesting);


/* Classifica um conjunto de teste do modo complexo, 
somando +1 aos nos do caminho com acerto e -1 aos nos do caminho com erro.
@return Acurácia.*/
float ComplexClassifyGraph(ComplGraph *cgTraining, 
			   ComplGraph *cgTesting);

/* Classifica um conjunto de teste do modo simples, 
atribui a cada elemento uma label.
@return Acurácia.*/
float SimpleClassifyGraph(ComplGraph *cgTraining, 
			  ComplGraph *cgTesting);

/*Grava o "ComplGraph" para disco em formato binario
onde os vetores de caracteristicas sao gravados de 
forma sequencial. O cabecalho armazena o tamanho
dos vetores, bem como o numero de classes e o numero
de registros por classe.*/
void  WriteComplGraph(ComplGraph *cg, char *filename);


//Le o "ComplGraph" a partir de um arquivo na forma binaria.
ComplGraph *ReadComplGraph(char *filename);

/*Troca elementos classificados incorretamente no teste por 
elementos nao relevantes do treinamento, independente se 
forem prototipos ou nao*/
void	    SwapErrorbyNotRelevant(ComplGraph *cgTraining, 
				   ComplGraph *cgTesting);

/*Troca elementos classificados incorretamente no teste por 
elementos  relevantes não protótipos do treinamento, 
independente se forem prototipos ou nao*/
void	    SwapErrorbyRelevant(ComplGraph *cgTraining, 
				ComplGraph *cgTesting);

/*Troca elementos do teste por elementos não relevantes 
do treinamento*/
void	    SwapNotRelevantbyNewSamples(ComplGraph *cgTraining, 
					ComplGraph * cgTesting);

/*Troca elementos errados do teste por elementos  não  protótipos do conjunto de treinamento, mesmo sendo
de classes diferentes. Analisa a capacidade de treinamento das classes*/
void SwapErrorBySampleFromDifferentClass(ComplGraph *cgTraining, ComplGraph * cgTesting);

inline void GraphNodeCopy(GraphNode *dest, GraphNode *src);
inline void GraphNodeSwap(GraphNode *a, GraphNode *b);
inline real GraphNodeDistance(GraphNode *a, GraphNode *b);

ComplGraph *CreateComplGraph(int n, int nclasses, int ndata);
void        DestroyComplGraph(ComplGraph **graph);
ComplGraph *CloneComplGraph(ComplGraph *cg);

void        ComplGraphMST(ComplGraph *cg);
void        ComplGraphTraining(ComplGraph *cg);
int         ComplGraphTestNode(ComplGraph *cg, GraphNode *node,  int *p);

/*calcula LC de um grafo e retorna em um vetor de floats*/
float  *GetLearningCapacity(ComplGraph *cg);

/*Incrementa/Decrementa valor de relevancia de todos os 
nos do caminho de um no utilizado para classificacao
correta ate a sua raiz*/
void	    SetRelevanceValuePath(ComplGraph *cg, 
				  int nodeID, int value);

/* Identifica se os nos de um determinado grafo sao 
nao relevantes. Um dado nó é dito ser nao relevante 
caso relevance_value <= 0.*/
void	 IdentifyNonRelevants(ComplGraph *cg);

void	 ShowPred(ComplGraph *cgTraining, int nodeID);

/*ordena um grafo de acordo com suas labels*/
void ComplGraphSort(ComplGraph **cg, char order);

/*executa o algoritmo OPF e retorna o classificador 
projetado final*/
ComplGraph *OPF(char *descriptorFileName, 
		char *accFileName, 
		int iterations);

/*remove elementos outliers (prototipos nao relevantes) de 
um grafo e retorna um grafo de rejeitados*/
ComplGraph * RemoveOutliers(ComplGraph **cg);

/*remove elementos nao relevantes  (contagem <= 0) de 
um grafo e retorna um grafo reduzido*/
ComplGraph * RemoveNotRelevants(ComplGraph **cg);

#endif

