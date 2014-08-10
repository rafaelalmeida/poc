#include<registration.h>
#include<adjacency3.h>
#include<morphology3.h>
#include<geometric3.h>
#include<algebra.h>
#include<matrix.h>

Set *getSetBorder(Scene *scn){
    
    AdjRel3 *adj=Spheric(1.0);
    Scene *borda=GetBorder3(scn, adj);
    Set *S=NULL;
    int i, n;
    n=scn->xsize*scn->ysize*scn->zsize;
    for(i=0;i<n;i++){
        if(borda->data[i]==1){
            InsertSet(&S, i);
        }
    }
    DestroyScene(&borda);
    DestroyAdjRel3(&adj);
    return S;
}

Set *getSetFaixa(Scene *scn){
    
    Scene *faixa=getFaixa(scn);
    Set *S=NULL;
    int i, n;
    n=scn->xsize*scn->ysize*scn->zsize;
    for(i=0;i<n;i++){
        if(faixa->data[i]==1){
            InsertSet(&S, i);
        }
    }
    DestroyScene(&faixa);
    return S;
}

Scene *getFaixa(Scene *scn){
	
	Set *set=NULL;
	Scene *d=DilateBin3(scn,&set,3.0);
	DestroySet(&set);
	Scene *e=ErodeBin3(scn,&set,3.0);
	Scene *faixa=Diff3(d,e);

	DestroyScene(&d);
	DestroyScene(&e);
	DestroySet(&set);
	return faixa;
}

Scene *alignPMS(Scene *scn1,Scene *scn2,RegParameters *p){
	Scene *align=CreateScene(scn1->xsize,scn1->ysize,scn1->zsize);
	float dz=(float)(-(int)scn2->zsize/2 + (int)align->zsize/2);
  (*p).RX=0.0; (*p).RY=0.0; (*p).RZ=0.0; (*p).TX=0.0; (*p).TY=0.0; (*p).TZ=dz;
  translacao((*p).T,0.0,0.0,dz);
	transformScene(scn2,(*p).T,align);
	return align;
}
Scene *alignPMS_bin(Scene *scn1,Scene *scn2,RegParameters *p){
	Scene *align=CreateScene(scn1->xsize,scn1->ysize,scn1->zsize);
	float dz=(float)(-(int)scn2->zsize/2 + (int)align->zsize/2);
  (*p).RX=0.0; (*p).RY=0.0; (*p).RZ=0.0; (*p).TX=0.0; (*p).TY=0.0; (*p).TZ=dz;
  translacao((*p).T,0.0,0.0,dz);
	transformScene_bin(scn2,(*p).T,align);
	return align;
}

Scene * RegMasks(Scene *scn1, Scene *scn2){
	Scene *result=CreateScene(scn1->xsize,scn1->ysize,scn1->zsize);
	int n=scn1->xsize*scn1->ysize*scn1->zsize;
	int i;
	
	for (i=0;i<n;i++){
	  if(scn1->data[i]!=0) scn1->data[i]=1;
		if(scn2->data[i]!=0) scn2->data[i]=2;
		result->data[i]=scn1->data[i]+scn2->data[i];
		
	}

  return result;

}

Scene * doGridRegistrationBin(Scene *source, Scene *target_align, float T[4][4]){ //gera imagem cinza quadriculada
    
    
    Scene *transform=CreateScene(target_align->xsize,target_align->ysize,target_align->zsize);
    transformScene_bin(target_align,T,transform);

    int n=source->xsize*source->ysize*source->zsize;
    int i;
    for(i=0;i<n;i++){
        if(source->data[i]==1)
            source->data[i]=150;
        if(transform->data[i]==1)
            transform->data[i]=80;
    }
    Scene *resultante=CreateScene(source->xsize,source->ysize,source->zsize);

    int delta=(source->xsize)/40;
    int Qd_x,Qd_y,Qd_z;
    int x,y,z;
    for(z=0;z<source->zsize;z++){
        Qd_z=z/delta;
        for(y=0;y<source->ysize;y++){
            Qd_y=y/delta;
            for(x=0;x<source->xsize;x++){
                Qd_x=x/delta;
                if(Qd_z%2==0){
                    if((Qd_y%2)==(Qd_x%2))
                        resultante->data[resultante->tbz[z]+resultante->tby[y]+x]=source->data[source->tbz[z]+source->tby[y]+x];
                    else
                        resultante->data[resultante->tbz[z]+resultante->tby[y]+x]=transform->data[transform->tbz[z]+transform->tby[y]+x];
                }
                else{
                    if((Qd_y%2)!=(Qd_x%2))
                        resultante->data[resultante->tbz[z]+resultante->tby[y]+x]=source->data[source->tbz[z]+source->tby[y]+x];
                    else
                        resultante->data[resultante->tbz[z]+resultante->tby[y]+x]=transform->data[transform->tbz[z]+transform->tby[y]+x];
                }
                
            }
        }
    }
        
    DestroyScene(&transform);
    return(resultante);

}


Scene * doGridRegistration(Scene *source_b, Scene *target_b_align, Scene *source, Scene *target_align, float T[4][4]){ //gera imagem cinza quadriculada
    Scene *transform_b=CreateScene(target_b_align->xsize,target_b_align->ysize,target_b_align->zsize);
    transformScene_bin(target_b_align,T,transform_b);
    
    
    Scene *transform=CreateScene(target_align->xsize,target_align->ysize,target_align->zsize);
    transformScene(target_align,T,transform);
    
    Scene *resultante=CreateScene(source->xsize,source->ysize,source->zsize);
    int delta=(source->xsize)/20;
    int Qd_x,Qd_y,Qd_z;
    int x,y,z;
    for(z=0;z<source->zsize;z++){
        Qd_z=z/delta;
        for(y=0;y<source->ysize;y++){
            Qd_y=y/delta;
            for(x=0;x<source->xsize;x++){
                Qd_x=x/delta;
                if(Qd_z%2==0){
                    if((Qd_y%2)==(Qd_x%2)&& source_b->data[source->tbz[z]+source->tby[y]+x]!=0)
                        resultante->data[resultante->tbz[z]+resultante->tby[y]+x]=source->data[source->tbz[z]+source->tby[y]+x];
                    else if((Qd_y%2)!=(Qd_x%2)&& transform_b->data[source->tbz[z]+source->tby[y]+x]!=0)
                        resultante->data[resultante->tbz[z]+resultante->tby[y]+x]=transform->data[transform->tbz[z]+transform->tby[y]+x];
                }
                else{
                    if((Qd_y%2)!=(Qd_x%2)&& source_b->data[source->tbz[z]+source->tby[y]+x]!=0)
                        resultante->data[resultante->tbz[z]+resultante->tby[y]+x]=source->data[source->tbz[z]+source->tby[y]+x];
                    else if((Qd_y%2)==(Qd_x%2)&& transform_b->data[source->tbz[z]+source->tby[y]+x]!=0)
                        resultante->data[resultante->tbz[z]+resultante->tby[y]+x]=transform->data[transform->tbz[z]+transform->tby[y]+x];
                }
                
            }
        }
    }
    

    DestroyScene(&transform_b);
    DestroyScene(&transform);
    return(resultante);

}
Scene * doGridRegistrationAND(Scene *source_b, Scene *target_b_align, Scene *source, Scene *target_align, float T[4][4]){ //gera imagem cinza quadriculada
    Scene *transform_b=CreateScene(target_b_align->xsize,target_b_align->ysize,target_b_align->zsize);
    transformScene_bin(target_b_align,T,transform_b);
    
    Scene *and=And3(source_b,transform_b);
    
    
    Scene *transform=CreateScene(target_align->xsize,target_align->ysize,target_align->zsize);
    transformScene(target_align,T,transform);
    
    Scene *resultante=CreateScene(source->xsize,source->ysize,source->zsize);
    int delta=(source->xsize)/20;
    int Qd_x,Qd_y,Qd_z;
    int x,y,z;
    for(z=0;z<source->zsize;z++){
        Qd_z=z/delta;
        for(y=0;y<source->ysize;y++){
            Qd_y=y/delta;
            for(x=0;x<source->xsize;x++){
                Qd_x=x/delta;
                if(Qd_z%2==0){
                    if((Qd_y%2)==(Qd_x%2)&& and->data[source->tbz[z]+source->tby[y]+x]!=0)
                        resultante->data[resultante->tbz[z]+resultante->tby[y]+x]=source->data[source->tbz[z]+source->tby[y]+x];
                    else if((Qd_y%2)!=(Qd_x%2)&& and->data[source->tbz[z]+source->tby[y]+x]!=0)
                        resultante->data[resultante->tbz[z]+resultante->tby[y]+x]=transform->data[transform->tbz[z]+transform->tby[y]+x];
                }
                else{
                    if((Qd_y%2)!=(Qd_x%2)&& and->data[source->tbz[z]+source->tby[y]+x]!=0)
                        resultante->data[resultante->tbz[z]+resultante->tby[y]+x]=source->data[source->tbz[z]+source->tby[y]+x];
                    else if((Qd_y%2)==(Qd_x%2)&& and->data[source->tbz[z]+source->tby[y]+x]!=0)
                        resultante->data[resultante->tbz[z]+resultante->tby[y]+x]=transform->data[transform->tbz[z]+transform->tby[y]+x];
                }
            }
        }
    }
    

    DestroyScene(&transform_b);
    DestroyScene(&transform);
    return(resultante);

}

void gerarTransformacoes(Scene *source, Scene *target, RegParameters P[candidates]){

    float RotZ[n_rotZ]={-5.,-1.,-0.5,-0.1,0.,0.1,0.5,1.,5.};
    float Dx[n_tX]={-5.,-1.,0.,1.,5.};
    float Dy[n_tY]={-5.,-1.,0.,1.,5.};
    
    int i,j,k,count;
    count =0;
    for(i=0;i<n_rotZ;i++){
    	for(j=0;j<n_tX;j++){
    		for(k=0;k<n_tY;k++){
    				RotTrans(0.0,0.0,RotZ[i],Dx[j],Dy[k],0.0,P[count].T,source,target);
    				P[count].RZ=RotZ[i];
        		P[count].TX=Dx[j];
        		P[count].TY=Dy[k];
        		count++;
    		}
    	}
    }
}
double CalcErroMI(Scene *target,Set *setfaixa_t,Scene *source,Scene *faixa_s,Curve *h1,int L1, int L2,float T[4][4]){ 

//h1: histograma da faixa do target
//h2: histograma da faixa do source
//h12: histograma conjunto 

double mi=0.0;
Curve *h2;
h2  = CreateCurve(L2+1); //histograma da faixa de source

/* 1. Calculo do histograma de Source e JoinHistograma (h12) */
Set *seed=setfaixa_t;
int p,i,j,l1,l2;
Voxel vox_t, vox_s;
Matrix *h12=CreateMatrix(L1+1,L2+1); //x:max_target y:max_source+1

int nelems=1;  //recebe 1 para evitar divisao por zero.
int h12_elements=1;
while(seed!=NULL){
    p=seed->elem;
    vox_t.z=p/(target->xsize*target->ysize);
    vox_t.y=(p-target->tbz[vox_t.z])/(target->xsize);
    vox_t.x=(p-target->tbz[vox_t.z])%(target->xsize);
    vox_s=Transform_Voxel(T, vox_t);
    if(ValidVoxel(source, vox_s.x, vox_s.y, vox_s.z) && faixa_s->data[faixa_s->tbz[vox_s.z]+faixa_s->tby[vox_s.y]+vox_s.x]==1){  // se o ponto tranformaco caiu na faixa de source
        l1=target->data[p];
        l2=source->data[source->tbz[vox_s.z]+source->tby[vox_s.y]+vox_s.x];
        h2->Y[l2]++;
        h12->val[h12->tbrow[l2]+l1]=(h12->val[h12->tbrow[l2]+l1])+1;		
        h12_elements++;
    }
    else{ // se o ponto caiu fora da faixa do source, o vai para o histograma de source o valor MaximumValue3(source)+1
        l1=target->data[p];
        l2=L2;
        h2->Y[l2]++;
        //h12->val[h12->tbrow[l2]+l1]=(h12->val[h12->tbrow[l2]+l1])+1;
    }
    seed=seed->next;
    nelems++;
}

/* 2. normalizar os histogramas */
for(i=0;i<h12->nrows;i++){
    h2->Y[i]=h2->Y[i]/nelems;
    for(j=0;j<h12->ncols;j++){
        h12->val[h12->tbrow[i]+j]=(h12->val[h12->tbrow[i]+j])/h12_elements;
    }
}

/* 3. calculo Informacao Mutua */
for(l1=0;l1<=L1;l1++){
    for(l2=0;l2<=L2;l2++){
        if(h1->Y[l1]* h2->Y[l2]!=0)
            mi=mi+( h12->val[h12->tbrow[l2]+l1] * log (1+ h12->val[h12->tbrow[l2]+l1]/( h1->Y[l1]* h2->Y[l2]  )));
    }
}


/* 4. destruir as variaveis  */
DestroyCurve(&h2);
DestroyMatrix(&h12);
DestroySet(&seed);
return (mi);
}

Curve *NormHistogram3_masc(Scene *scn, Scene *masc){
  int i, n, nbins,nelems=0;
  Curve *hist = NULL;

  nbins = MaximumValue3(scn)+1;

  hist  = CreateCurve(nbins);
  n = scn->xsize * scn->ysize * scn->zsize;
  
  for (i = 0; i < n; i++){
  	if(masc->data[i]==1){
  		hist->Y[scn->data[i]]++;
  		nelems++;
  	}
  }

  for (i = 0; i < nbins; i++){
    hist->X[i] = i;
    hist->Y[i] = hist->Y[i]/nelems;
  }
 
  return (hist);
}
double CalcErro_M3(Scene *target,Set *setfaixa_t,Scene *source,Scene *faixa_s,float T[4][4]){ 
	double err=0.0;
  Set *seed=setfaixa_t;
  int p,i_max,dif;
  Voxel vox_t, vox_s;
  i_max=MaximumValue3(source);
  int nelems=0;
  while(seed!=NULL){
      p=seed->elem;
      vox_t.z=p/(target->xsize*target->ysize);
      vox_t.y=(p-target->tbz[vox_t.z])/(target->xsize);
      vox_t.x=(p-target->tbz[vox_t.z])%(target->xsize);
      vox_s=Transform_Voxel(T, vox_t);
			if( ValidVoxel(source, vox_s.x, vox_s.y, vox_s.z) && faixa_s->data[faixa_s->tbz[vox_s.z]+faixa_s->tby[vox_s.y]+vox_s.x]==1){  // se o ponto tranformaco caiu na faixa de source
				dif=target->data[p] - source->data[source->tbz[vox_s.z]+source->tby[vox_s.y]+vox_s.x];
				err=err+pow(dif,2);	
			}
			else{ // se o ponto caiu fora da faixa do source, o vai para o histograma de source o valor MaximumValue3(source)+1
				err=err+pow(i_max,2);
			}
			seed=seed->next;
			nelems++;
	}
	
  DestroySet(&seed);
	return (err);
}

double CalcErro_M2(Scene *target, Set *S_faixa, Scene *source, float t[4][4]){

    double err=0.0;
    Set *seed=S_faixa;
    int p;
    Voxel vox_t, vox_s;
    while(seed!=NULL){
        p=seed->elem;
        vox_t.z=p/(target->xsize*target->ysize);
        vox_t.y=(p-target->tbz[vox_t.z])/(target->xsize);
        vox_t.x=(p-target->tbz[vox_t.z])%(target->xsize);
        vox_s=Transform_Voxel(t, vox_t);

        if(ValidVoxel(source, vox_s.x, vox_s.y, vox_s.z)&& VoxelValue(source, vox_s)!=VoxelValue(target, vox_t))
            err++; 
        if(!ValidVoxel(source, vox_s.x, vox_s.y, vox_s.z))
            err++;    
        seed=seed->next;
        
    }
    DestroySet(&seed);
    return err;
}

double CalcErro_M1(Scene *target, Set *S, Scene *source_border, float t[4][4]){

    double err=0.0;
    Set *seed=S;
    int p;
    Voxel vox_t, vox_s;
    while(seed!=NULL){
        p=seed->elem;
        vox_t.z=p/(target->xsize*target->ysize);
        vox_t.y=(p-target->tbz[vox_t.z])/(target->xsize);
        vox_t.x=(p-target->tbz[vox_t.z])%(target->xsize);
        vox_s=Transform_Voxel(t, vox_t);
        if(ValidVoxel(source_border, vox_s.x, vox_s.y, vox_s.z)&& VoxelValue(source_border, vox_s)==0) //nao coincidiu com a borda
            err++;
        if(!ValidVoxel(source_border, vox_s.x, vox_s.y, vox_s.z))
            err++;
        seed=seed->next;
    }
    DestroySet(&seed);
    
    return err;
}

void createIdentity(float T[4][4]){
	int a,b; 
    for(a=0;a<4;a++){
        for(b=0;b<4;b++){
            if(a==b)
                T[a][b]=1.0;
            else
                T[a][b]=0.0;
        }
    }
}
void SearchRegistrationFunction(int method, Scene *source_bin, Scene *target_bin, Scene *source, Scene *target, RegParameters *p){
    // inicializando p
    createIdentity((*p).T);
    (*p).RX=0.0; (*p).RY=0.0; (*p).RZ=0.0; (*p).TX=0.0; (*p).TY=0.0; (*p).TZ=0.0;
  
    RegParameters node; // nó i
    RegParameters P[candidates]; // search basis
    gerarTransformacoes(source_bin,target_bin,P);
    
   	int k,i,j,l=0;
    double d_M=0., d=0.;
    //declaracao das variaveis
    AdjRel3 *adj;
    Set *S;
    Scene *source_border,*grad_s,*grad_t,*faixa_s,*faixa_t;
    Curve *hist_faixaT;
    int L1=0,L2=0;
    //Search Algorithm
    switch(method){
    	case 1:{
    	 adj=Spheric(1.0);
  	   source_border=GetBorder3(source_bin, adj);
	     S=getSetBorder(target_bin);
	     d_M=CalcErro_M1(target_bin,S,source_border,(*p).T);
    	 break;}
    	case 2:{
    		S=getSetFaixa(target_bin);
    	  d_M=CalcErro_M2(target_bin,S,source_bin,(*p).T);
    	  break;}
    	case 3:{
    		adj=Spheric(1.0);
    		grad_s=MorphGrad3(source,adj);
    		grad_t=MorphGrad3(target,adj); 
		    faixa_s=getFaixa(source_bin);
   		  S=getSetFaixa(target_bin);
   		  d_M=CalcErro_M3(grad_t,S,grad_s,faixa_s,(*p).T);
   		  break;}
			case 4:{
				faixa_s=getFaixa(source_bin);
    		faixa_t=getFaixa(target_bin);
    		S=getSetFaixa(target_bin);
        hist_faixaT=NormHistogram3_masc(target,faixa_t);
    		DestroyScene(&faixa_t);
    		L1=MaximumValue3(target);
        L2=MaximumValue3(source)+1;
    		d_M=CalcErroMI(target,S,source,faixa_s,hist_faixaT,L1,L2,(*p).T);
    		break;}
    	case 5:{
    	  adj=Spheric(1.0);
    		grad_s=MorphGrad3(source,adj); 
		    grad_t=MorphGrad3(target,adj); 
    		faixa_s=getFaixa(source_bin);
    		faixa_t=getFaixa(target_bin);
    		S=getSetFaixa(target_bin);
    		hist_faixaT=NormHistogram3_masc(grad_t,faixa_t);
    		DestroyScene(&faixa_t);	
     		L1=MaximumValue3(grad_t);
		    L2=MaximumValue3(grad_s)+1;
				d_M=CalcErroMI(grad_t,S,grad_s,faixa_s,hist_faixaT,L1,L2,(*p).T);
    		break;}
    
    }
    
    do {
    		k=-1; 
        for(i=0;i<candidates;i++){
            MultMatrices(P[i].T, (*p).T, node.T); 
            node.RZ=P[i].RZ+p->RZ;  node.TX=P[i].TX+p->TX; node.TY=P[i].TY+p->TY;
       
            switch (method){
            	case 1:
            		
            		d=CalcErro_M1(target_bin,S,source_border,node.T); 
            		break;
            	case 2:
            		d=CalcErro_M2(target_bin,S,source_bin,node.T); 
            		break;
            	case 3:
            		d=CalcErro_M3(grad_t,S,grad_s,faixa_s,node.T);
            		break;
            	case 4:
            		d=CalcErroMI(target,S,source,faixa_s,hist_faixaT,L1,L2,node.T);
            		break;
            	case 5:
            		d=CalcErroMI(grad_t,S,grad_s,faixa_s,hist_faixaT,L1,L2,node.T);	
            		break;
            }
            if(method==1||method==2||method==3){
            	if(d<d_M){
            		d_M=d;
            		k=i; 
            	}
            }	
            else if(method==4||method==5){
            	 if(d>d_M){
            		d_M=d;
            		k=i;
            	}
            }
         }
         if(k!=-1){
         	
           	MultMatrices(P[k].T, (*p).T, node.T);
           	(*p).RZ+=P[k].RZ; (*p).TX+=P[k].TX;(*p).TY+=P[k].TY;
           	for(i=0;i<4;i++)
           		for(j=0;j<4;j++)
           			(*p).T[i][j]=node.T[i][j];
           
           	
           	printf("\tSearch iteration %d: RZ/%.1f TX/%.1f TY/%.1f   Score: %f\n",++l, P[k].RZ,P[k].TX,P[k].TY,d_M);
         }
         
    } while(k!=-1);
    
    printf("\tFinal Registration RegParameters: RZ%.1f TX%.1f TY%.1f\n",(*p).RZ,(*p).TX,(*p).TY);

    // Destruindo Variaveis  
    switch(method){
    	case 1:
    		DestroyAdjRel3(&adj);
      	DestroySet(&S);
      	DestroyScene(&source_border);
      	break;
      case 2:
      	DestroySet(&S);
      	break;
      case 3:
      	DestroyAdjRel3(&adj);
      	DestroySet(&S);	
      	DestroyScene(&grad_s);
      	DestroyScene(&grad_t);
      	DestroyScene(&faixa_s);
      	break;
      case 4:
      	DestroySet(&S);
      	DestroyScene(&faixa_s);	
      	DestroyCurve(&hist_faixaT);
      	break;
      case 5: 
        DestroyAdjRel3(&adj);	
      	DestroySet(&S);
      	DestroyScene(&faixa_s);	
      	DestroyCurve(&hist_faixaT);
      	DestroyScene(&grad_s);
      	DestroyScene(&grad_t);	
      	break;	
    }
   
}


