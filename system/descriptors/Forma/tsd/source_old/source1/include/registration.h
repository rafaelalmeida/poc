
#ifndef _REGISTRATION_H
#define	_REGISTRATION_H

#include<set.h>
#include<scene.h>
#include<curve.h>

#define n_rotZ 9
#define n_tX 5
#define n_tY 5
#define candidates 225



typedef struct _searchargs {
  float T[candidates][4][4];
  double err[candidates];
  
  
} SearchArgs;

typedef struct _regparameters{
	float RX,RY,RZ,TX,TY,TZ;
	float T[4][4];
} RegParameters;

typedef struct _search{
	RegParameters p;
	double score;
}Search;

Set *getSetBorder(Scene *scn);
Set *getSetFaixa(Scene *scn);
Scene *getFaixa(Scene *scn);

Scene *alignPMS(Scene *scn1,Scene *scn2,RegParameters *p);
Scene *alignPMS_bin(Scene *scn1,Scene *scn2,RegParameters *p);

Scene *doGridRegistrationBin(Scene *source, Scene *target, float T[4][4]);
Scene *doGridRegistration(Scene *source_b, Scene *target_b, Scene *source, Scene *target, float T[4][4]); // gera o mosaico de todo o conteudo das scenas
Scene *doGridRegistrationAND(Scene *source_b, Scene *target_b, Scene *source, Scene *target, float T[4][4]); // gera o mosaico apenas na area de interseccao entre as mascaras
Scene * RegMasks(Scene *scn1, Scene *scn2); // label 1: onde scn1 =1 e scn2 =0
																						// label 2:      scn1 =0 e scn2 =1
																						// label 3:      scn1 =1 e scn2 =1 								

double CalcErroMI(Scene *target,Set *setfaixa_t,Scene *source,Scene *faixa_s,Curve *h1,int L1, int L2,float T[4][4]);
Curve *NormHistogram3_masc(Scene *scn, Scene *masc);
double CalcErro_M3(Scene *target,Set *setfaixa_t,Scene *source,Scene *faixa_s,float T[4][4]);
double CalcErro_M2(Scene *target, Set *S_faixa, Scene *source, float t[4][4]);
double CalcErro_M1(Scene *target, Set *S, Scene *source_border, float t[4][4]);

void gerarTransformacoes(Scene *source, Scene *target, RegParameters P[candidates]);
void SearchRegistrationFunction(int method, Scene *source_bin, Scene *target_bin, Scene *source, Scene *target, RegParameters *p);
void createIdentity(float T[4][4]);


#endif	/* _REGISTRATION_H */

