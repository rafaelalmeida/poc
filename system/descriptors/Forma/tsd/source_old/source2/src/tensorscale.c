#include "tensorscale.h" 

TensorScale *CreateTensorScale(Image *img, Image *mask, int m_pairs, int L, float stdDev){
    Point *epsilon;
    TensorScale *ts;
    Vector *tau;
    int i,j,k,p,v;
    float val_p,val1,val2,val,threshold,threshold2,difacc1,dif;
    float x,y, taux, tauy;
    float a2,b2,teta,l1E,l2E,uu,vv,aux;
    float gSxy, gSy2_x2;
    float Su2v2,Su4,Sv4,Su2,Sv2;
    float sin_teta, cos_teta;
    int L2;
    int x1,x2,y1,y2,a,b,c,d,pp,qq;
    float dx1,dx2,dy1,dy2;
    int ncols,nrows,destroymask;

    if (DEBUG) printf("CreateTensorScale:: Entrou!\n");

    L2 = L*L;
    ncols = img->ncols;
    nrows = img->nrows;
    if(mask==NULL){
        mask = CreateImage(ncols,nrows);
        SetImage(mask, 1);
        destroymask = 1;
    }
    else {
        destroymask = 0;
    }
    if (DEBUG) printf("CreateTensorScale:: Verificou mascara!\n");

    ts = (TensorScale *)malloc(sizeof(TensorScale));
    ts->orientation = CreateDImage(ncols,nrows);
    ts->anisotropy = CreateDImage(ncols,nrows);
    ts->thickness = CreateDImage(ncols,nrows);
    ts->m_pairs = m_pairs;
    //ts->L = L;

    if (DEBUG) printf("CreateTensorScale:: Criou imgs orientation, anisotropy e thickness!\n");

    tau = (Vector *)malloc(sizeof(Vector)*m_pairs);
    epsilon = (Point *)malloc(sizeof(Point)*m_pairs);

    teta = 0.0;
    for(i=0;i<m_pairs;i++){
        tau[i].x = cos(teta);
        tau[i].y = sin(teta);
        tau[i].z = 0.0;

        teta += ((float)PI/m_pairs); 
    }

    threshold = 3.0*stdDev;
    threshold2 = 2*threshold;

    if (DEBUG) printf("CreateTensorScale:: Vai percorrer a imagem...\n");
    for(i=1;i<nrows-1;i++){
        for(j=1;j<ncols-1;j++){
            p=img->tbrow[i]+j; //estranho isso!
            if (DEBUG) printf("\tCreateTensorScale:: p=%d\t", p);
            val_p = 2*img->val[p];
            if (DEBUG) printf("val_p recebeu 2*img->val[p]\t");
            val_p += img->val[p-ncols];
            if (DEBUG) printf("p-ncols=%d\t", p-ncols);
            val_p += img->val[p+ncols];
            if (DEBUG) printf("p+ncols=%d\t", p+ncols);
            val_p += img->val[p-1];
            if (DEBUG) printf("p-1=%d\t", p-1);
            val_p += img->val[p+1];
            if (DEBUG) printf("p+1=%d\t", p+1);
            val_p /= 6.0;
            if (DEBUG) printf("val_p / 6 !!\n");


            if(mask->val[p] == 0) {
                continue;
            }

            //--------- Sample lines --------------------------------------
            if (DEBUG) printf("\tCreateTensorScale:: INICIO - Sample Lines...\n");
            gSxy = gSy2_x2 = 0.0;
            for(k=0;k<m_pairs;k++){
                x = taux = tau[k].x;
                y = tauy = tau[k].y;
                v=1;
                difacc1 =  0.0;
                val1 = val2 = val_p;
                for(;v<L;v++){
                    x1 = (int)(x+j);
                    y1 = (int)(y+i);
                    if(x1>=ncols-1 || y1>=nrows-1 || x1<=0 || y1<=0)
                        break;
                    x2 = x1+1;
                    y2 = y1+1;

                    pp = img->tbrow[y1];
                    qq = img->tbrow[y2];

                    a = img->val[pp+x1];
                    b = img->val[pp+x2];
                    c = img->val[qq+x1];
                    d = img->val[qq+x2];

                    dx1 = x+j-x1; dx2 = x2-x-j;
                    dy1 = y+i-y1; dy2 = y2-y-i;

                    val = val1;
                    val1 = ((float)d*dx1+(float)c*dx2)*dy1+((float)b*dx1+(float)a*dx2)*dy2;
                    val = (val+val1)/2.0;
                    dif = ((val>val_p)?(val-val_p):(val_p-val));
                    difacc1 += dif;
                    if( dif > threshold )
                        break;

                    x1 = (int)(-x+j);
                    y1 = (int)(-y+i);
                    if(x1>=ncols-1 || y1>=nrows-1 || x1<=0 || y1<=0)
                        break;
                    x2 = x1+1;
                    y2 = y1+1;

                    pp = img->tbrow[y1];
                    qq = img->tbrow[y2];

                    a = img->val[pp+x1];
                    b = img->val[pp+x2];
                    c = img->val[qq+x1];
                    d = img->val[qq+x2];

                    val = val2;
                    val2 = ((float)d*dx2+(float)c*dx1)*dy2+((float)b*dx2+(float)a*dx1)*dy1;
                    val = (val+val2)/2.0;
                    dif = ((val>val_p)?(val-val_p):(val_p-val));
                    difacc1 += dif;
                    if( dif > threshold )
                        break;
                    if( difacc1 > threshold2 ){
                        v -=1;
                        x -= taux;
                        y -= tauy;
                        break;
                    }

                    x += taux;
                    y += tauy;
                }
                epsilon[k].x = x;
                epsilon[k].y = -y;

                gSxy += x*(-y);
                gSy2_x2 += (y+x)*(y-x); //(y*y-x*x);
            }
            if (DEBUG) printf("\tCreateTensorScale:: FIM - Sample Lines!\n");

            //-------------------- TETA -----------------------------------
            if (DEBUG) printf("\tCreateTensorScale:: INICIO - TETA...\n");
            if(gSy2_x2==0.0){ 
                if(gSxy>0.0) teta=PI/2.0;
                else teta=-PI/2.0;
            } else {
                teta = atanf((gSxy+gSxy)/gSy2_x2);

                if(gSxy<0.0 && teta<0.0) teta+=PI;
                else if(gSxy>0.0 && teta>0.0) teta-=PI;
                else if(teta==0.0 && gSy2_x2>0.0) teta=PI;
            }
            teta /= 2.0;
            if (DEBUG) printf("\tCreateTensorScale:: FIM - TETA!\n");

            //----------------- A & B ---------------------------------
            if (DEBUG) printf("\tCreateTensorScale:: INICIO - A & B...\n");
            sin_teta = sin(teta);
            cos_teta = cos(teta);
            Su2v2 = Su2 = Sv2 = Su4 = Sv4 = 0.0;
            for(k=0;k<m_pairs;k++) {
                uu = epsilon[k].x;
                vv = epsilon[k].y;
                aux = uu*cos_teta-vv*sin_teta;
                vv = vv*cos_teta+uu*sin_teta;
                uu = aux*aux;
                vv *= vv;
                Su2v2 += uu*vv;
                Su2 += uu;
                Sv2 += vv;
                Su4 += uu*uu;
                Sv4 += vv*vv;
            }
            aux = Su2v2*Su2v2-Su4*Sv4;
            a2 = (aux/(Su2v2*Sv2-Su2*Sv4));
            b2 = (aux/(Su2*Su2v2-Su4*Sv2));

            if(isnan(a2) || a2<0.0) a2=1.0;
            if(isnan(b2) || b2<0.0) b2=1.0;
            if(a2>L2) a2=L2;
            if(b2>L2) b2=L2;

            if(a2>=b2){
                l1E=a2;
                l2E=b2;
            } else {
                l1E=b2;
                l2E=a2;
            }
            if (DEBUG) printf("\tCreateTensorScale:: A & B - ts->X\t");
            ts->anisotropy->val[p]=ROUND((INT_MAX)*(1.0-sqrt(l2E/l1E)));
            ts->thickness->val[p]=ROUND((INT_MAX)*sqrt(sqrt(l2E)/L));
            if (DEBUG) printf("A & B - ts->X gravou!\n");

            if(teta<0) teta+=(float)PI;
            if(teta>PI) teta-=((int)(teta/PI))*PI;
            teta = PI-teta;

            //ts->orientation->val[p]=ROUND(teta*((INT_MAX)/PI));
            ts->orientation->val[p]=ROUND((teta/PI)*360); //AGORA O ANGULO FICA ENTRE 0 e 360
            //printf("orientation[%d]=%lf\t teta=%f\tINT_MAX=%d\tPI=%f\t**=%f\n", p, ts->orientation->val[p], teta, INT_MAX, PI, teta/PI*360);
            if (DEBUG) printf("\tCreateTensorScale:: FIM - A & B!\n");
        }
        if (DEBUG) printf("CreateTensorScale:: linha %d ja foi!\n", i);
    }
    if (DEBUG) printf("CreateTensorScale:: Percorreu a imagem!\n");

    free(tau);
    free(epsilon);

    if(destroymask) DestroyImage(&mask);

    return ts;
}

Image *TSEDistTrans(Image *bin) {
  Image *Dx=NULL,*Dy=NULL,*cost;
  Queue *Q=NULL;
  int i,p,q,n;
  Pixel u,v;
  int *sq=NULL,tmp=INT_MAX,dx,dy;
  AdjRel *A;

  n  = MAX(bin->ncols,bin->nrows);
  sq = AllocIntArray(n);
  for (i=0; i < n; i++) 
    sq[i]=i*i;

  A = Circular(1.5);
  cost = CreateImage(bin->ncols, bin->nrows);
  Dx = CreateImage(bin->ncols,bin->nrows);
  Dy = CreateImage(bin->ncols,bin->nrows);  
  n  = bin->ncols*bin->nrows;
  Q = CreateQueue(bin->ncols+bin->nrows,n);

  for(p = 0; p < n; p++){
    if(bin->val[p] > 0)
      cost->val[p] = INT_MAX;	  
    else{
      cost->val[p]=0;    
      InsertQueue(Q,cost->val[p]%Q->C.nbuckets,p);
    }
  }
  
  while(!EmptyQueue(Q)) {
    p   = RemoveQueue(Q);
    u.x = p%bin->ncols;
    u.y = p/bin->ncols;
    for (i=1; i < A->n; i++){
      v.x = u.x + A->dx[i];
      v.y = u.y + A->dy[i];
      if(ValidPixel(bin,v.x,v.y)){
	q = v.x + bin->tbrow[v.y];
	if (cost->val[p] < cost->val[q]){
	  dx  = Dx->val[p] + abs(v.x-u.x);
	  dy  = Dy->val[p] + abs(v.y-u.y);
	  tmp = sq[dx] + sq[dy];
	  if (tmp < cost->val[q]){
	    if (cost->val[q] == INT_MAX)
	      InsertQueue(Q,tmp%Q->C.nbuckets,q);
	    else
	      UpdateQueue(Q,q,cost->val[q]%Q->C.nbuckets,tmp%Q->C.nbuckets);
	    cost->val[q]  = tmp;
	    Dx->val[q] = dx;
	    Dy->val[q] = dy;
	  }
	}
      }
    }
  }

  free(sq);
  DestroyQueue(&Q);
  DestroyAdjRel(&A);
  DestroyImage(&Dx);
  DestroyImage(&Dy);
  
  return(cost);
}

AnnImg *TSEDistTransNew(Image *bin, Image *lambda){
  Image *Dx=NULL,*Dy=NULL;
  Queue *Q=NULL;
  int i,p,q,n;
  Pixel u,v;
  int *sq=NULL,tmp=INT_MAX,dx,dy;
  AdjRel *A;
  AnnImg *aimg;
  
  aimg = CreateAnnImg(bin);
  n  = MAX(bin->ncols,bin->nrows);
  sq = AllocIntArray(n);
  for (i=0; i < n; i++) 
    sq[i]=i*i;

  A = Circular(1.5);
  Dx = CreateImage(bin->ncols,bin->nrows);
  Dy = CreateImage(bin->ncols,bin->nrows);  
  n  = bin->ncols*bin->nrows;
  Q = CreateQueue(bin->ncols+bin->nrows,n);

  if (lambda != NULL) aimg->label = CopyImage(lambda);
  else aimg->label = CreateImage(bin->ncols,bin->nrows);

  for (p=0; p < n; p++) {
    aimg->cost->val[p]  = INT_MAX;
    aimg->pred->val[p]  = NIL;
    if ((aimg->label->val[p] != 0) || (bin->val[p] == 0)) {
      aimg->cost->val[p]  = 0;
      aimg->root->val[p]  = p;
      InsertQueue(Q,aimg->cost->val[p]%Q->C.nbuckets,p);
    }
  }

  while(!EmptyQueue(Q)) {
    p   = RemoveQueue(Q);
    u.x = p%bin->ncols;
    u.y = p/bin->ncols;
    for (i=1; i < A->n; i++){
      v.x = u.x + A->dx[i];
      v.y = u.y + A->dy[i];
      if(ValidPixel(bin,v.x,v.y)){
        q = v.x + bin->tbrow[v.y];
        if (aimg->cost->val[p] < aimg->cost->val[q]){
          dx  = Dx->val[p] + abs(v.x-u.x);
          dy  = Dy->val[p] + abs(v.y-u.y);
          tmp = sq[dx] + sq[dy];
          if (tmp < aimg->cost->val[q]){
            if (aimg->cost->val[q] == INT_MAX)
              InsertQueue(Q,tmp%Q->C.nbuckets,q);
            else
              UpdateQueue(Q,q,aimg->cost->val[q]%Q->C.nbuckets,tmp%Q->C.nbuckets);
            aimg->cost->val[q]  = tmp;
            aimg->pred->val[q]  = p;
            aimg->label->val[q] = aimg->label->val[p];
            aimg->root->val[q]  = aimg->root->val[p];
            Dx->val[q] = dx;
            Dy->val[q] = dy;
          }
        }
      }
    }
  }

  free(sq);
  DestroyQueue(&Q);
  DestroyAdjRel(&A);
  DestroyImage(&Dx);
  DestroyImage(&Dy);
  
  return(aimg);
}


TensorScale *CreateBinaryTensorScale(Image *bin, int m_pairs){
  AnnImg *aimg;
  Image *edt,*edt2,*bin2;
  Point *epsilon;
  TensorScale *ts = NULL;
  Vector *tau;
  int i,j,k,p,v,vi,n;
  float x,y,xc,yc,taux, tauy, aux;
  float a2,b2,b1,teta,u1,v1,u2,v2,aa,acc,wt,w;
  float gSxy, gSy2_x2;
  float sin_teta, cos_teta;
  float *lt_sqrt;
  int d,d1,d2,dmax;
  int ncols,nrows;

  ncols = bin->ncols;
  nrows = bin->nrows;
  n = ncols*nrows;

  bin2 = CopyImage(bin);
  p = bin2->tbrow[nrows-1];
  for(i=0; i<ncols; i++) {
    bin2->val[i] = 0;
    bin2->val[p+i] = 0;
  }
  for(i=0; i<nrows; i++) { 
    p = bin2->tbrow[i]; 
    bin2->val[p] = 0; 
    bin2->val[p+ncols-1] = 0;
  } 
  
  //------ Euclidean distance transform ----------------
  aimg = TSEDistTransNew(bin2, NULL);
  edt2 = CopyImage(aimg->cost);
  DestroyAnnImg(&aimg);

  dmax = MaximumValue(edt2);
  lt_sqrt = (float *)malloc((dmax+1)*sizeof(float));
  for(i=0;i<=dmax;i++)
    lt_sqrt[i] = sqrtf((float)i);

  edt = CreateImage(ncols,nrows);
  for(p=0; p<n; p++){
    d = edt2->val[p];
    d = ROUND(lt_sqrt[d]);
    edt->val[p] = d;
  }

  //---------------------------------------------------

  ts = (TensorScale *)malloc(sizeof(TensorScale));
  ts->orientation = CreateDImage(ncols,nrows);
  ts->anisotropy  = CreateDImage(ncols,nrows);
  ts->thickness   = CreateDImage(ncols,nrows);
  ts->m_pairs     = m_pairs;

  tau = (Vector *)malloc(sizeof(Vector)*m_pairs);
  epsilon = (Point *)malloc(sizeof(Point)*m_pairs);

  teta = 0.0;
  for(i=0;i<m_pairs;i++){
    tau[i].x = cosf(teta);
    tau[i].y = sinf(teta);
    tau[i].z = 0.0;

    teta += ((float)PI/m_pairs); 
  }

  for(i=1;i<nrows-1;i++){
    for(j=1;j<ncols-1;j++){
      p = bin2->tbrow[i]+j;
      if(bin2->val[p] == 0) continue;

      vi = edt->val[p];
      //--------- Sample lines --------------------------------------
      gSxy = gSy2_x2 = 0.0;
      xc = j+0.5; yc = i+0.5;
      for(k=0;k<m_pairs;k++){
        taux = tau[k].x;
        tauy = tau[k].y;
        v = vi;
        d1 = d2 = 0;

        while(1){
          x = v*taux;
          y = v*tauy;

          if(d1==0){
            d1 = edt->val[(int)(xc+x) + edt->tbrow[(int)(yc+y)]];
            if(d1 == 0)
              break;
          }
          
          if(d2==0){
            d2 = edt->val[(int)(xc-x) + edt->tbrow[(int)(yc-y)]];
            if(d2 == 0)
              break;
          }

          d = (d1<d2)?d1:d2;
          d1 -= d;
          d2 -= d;
          v += d;
        }

        epsilon[k].x = x;
        epsilon[k].y = -y;

        gSxy -= x*y;            //gSxy += x*(-y);
        gSy2_x2 += (y+x)*(y-x); //(y*y-x*x);
      }

      //-------------------- TETA -----------------------------------
      
      if(gSy2_x2==0.0){ 
        if(gSxy>0.0) teta=PI/2.0;
        else teta=-PI/2.0;
      }
      else{
        teta = atanf((gSxy+gSxy)/gSy2_x2);
        
        if(gSxy<0.0 && teta<0.0) teta+=PI;
        else if(gSxy>0.0 && teta>0.0) teta-=PI;
        else if(teta==0.0 && gSy2_x2>0.0) teta=PI;
      }
      teta /= 2.0;

      //----------------- A & B ---------------------------------
      b2   = (float)edt2->val[p];
      b1   = lt_sqrt[(int)b2];
      
      acc = wt = 0.0;
      sin_teta = sinf(teta);
      cos_teta = cosf(teta);
      for(k=0;k<m_pairs;k++){
        x = epsilon[k].x;
        y = epsilon[k].y;
        
        v1 = y*cos_teta + x*sin_teta;
        u1 = x*cos_teta - y*sin_teta;
        
        v2 = v1*v1;
        u2 = u1*u1;

        if(v2<b2){
          aa = b2*u2/(b2-v2);
          if(aa>=b2){
            w = (v1<0.0)?(b1+v1):(b1-v1);
            acc += w*aa;
            wt  += w;
          }
        }
      }

      if(wt>0.0)
        a2 = acc/wt;
      else
        a2 = b2;

      aux = 1.0-b2/a2;
      if(aux<0.0) aux = 0.0;
      ts->anisotropy->val[p] = sqrtf(aux);
      ts->thickness->val[p]  = b1; //sqrtf(b2);

      if(teta<0.0) teta+=(float)PI;
      if(teta>PI)  teta = PI;
      teta = PI-teta;

      ts->orientation->val[p] = teta;
      //printf("orientation=%lf\n", ts->orientation->val[p]);
    }
  }

  free(tau);
  free(epsilon);
  free(lt_sqrt);
  DestroyImage(&edt);
  DestroyImage(&edt2);
  DestroyImage(&bin2);

  return ts;
}


void DestroyTensorScale(TensorScale **ts){
  if(*ts != NULL){
    DestroyDImage(&((*ts)->orientation));
    DestroyDImage(&((*ts)->anisotropy));
    DestroyDImage(&((*ts)->thickness));
    free(*ts);
  }
  *ts = NULL;
}

CImage *ConvertTS2CImage(TensorScale *ts){
  CImage *cimg;
  int i,j,w,h,p;
  int RR,GG,BB,tRGB;
  int HH,SS,VV;
  double tmax;

  w = ts->anisotropy->ncols;
  h = ts->anisotropy->nrows;
  cimg = CreateCImage(w,h);

  tmax = MaximumDImageValue(ts->thickness);
  for(i=0; i<h; i++){
    for(j=0; j<w; j++){
      p=ts->anisotropy->tbrow[i]+j;
      SS=ROUND(ts->anisotropy->val[p]*255.0);
      VV=ROUND(255.0*sqrt(ts->thickness->val[p]/tmax));
      HH=ROUND(ts->orientation->val[p]*(255.0/PI));

      tRGB = HSV2RGB(triplet(HH,SS,VV));
      RR = t0(tRGB);
      GG = t1(tRGB);
      BB = t2(tRGB);
      
      cimg->C[0]->val[p]=RR;
      cimg->C[1]->val[p]=GG;
      cimg->C[2]->val[p]=BB;      
    }
  }

  return cimg;
}


void OutputTSColorSpace(char *filename, int size){
  CImage *cimg_colorspace;
  int i,j,x,y,p;
  float val,teta;
  int HH,SS,VV,tRGB,RR,GG,BB;

  cimg_colorspace = CreateCImage(size,size);
  VV=255;
  for(i=0;i<size;i++){
    for(j=0;j<size;j++){
      y = -(i-size/2);
      x = j-size/2;
      val = sqrt(pow(x,2)+pow(y,2));
      if(x==0 && y>0) teta=PI/2;
      else if(x==0 && y<=0) teta=-PI/2;
      else teta = atan((float)y/(float)x); 
      if(x<0) teta+=PI;
      if(teta<0) teta+=2*PI;
      if(val<size/2-5){
        HH = (ROUND(teta*256.0/PI)%256);
        SS = ROUND(val*(255.0/(size/2-5)));
        tRGB = HSV2RGB(triplet(HH,SS,VV));
        RR = t0(tRGB);
        GG = t1(tRGB);
        BB = t2(tRGB);
        p = j+cimg_colorspace->C[0]->tbrow[i];
        cimg_colorspace->C[0]->val[p]=RR;
        cimg_colorspace->C[1]->val[p]=GG;
        cimg_colorspace->C[2]->val[p]=BB;
      }
    }
  }

  WriteCImage(cimg_colorspace,filename);
  DestroyCImage(&cimg_colorspace);
}

/***********************************************************/
/* TSD - Tensor Scale Descriptor                           */
/***********************************************************/

/*****************TSD EXTRACTION ALGORITHM******************/
float *TSOrientationHistogram(TensorScale *ts){
    float *hist;
    float ratio,sum;
    double an,th;
    int w,h,i,j,p,bin=-1;

    if (DEBUG) printf("TSOrientationHistogram:: Iniciou\n");
    //ratio = (float)HISTOGRAMSIZE/PI;
    //ratio = (float)HISTOGRAMSIZE/255.0;
    ratio = (float) 1 / (360/HISTOGRAMSIZE); //igual a 1/(360/180)=0.5
    if (DEBUG) printf("TSOrientationHistogram:: ratio=%f\n", ratio);
    hist = (float *)calloc((HISTOGRAMSIZE+1), sizeof(float));
    //memset(hist, 0, sizeof(float)*(HISTOGRAMSIZE+1));
    w = ts->anisotropy->ncols;
    h = ts->anisotropy->nrows;

    if (DEBUG) printf("TSOrientationHistogram:: Vai percorrer a TS Image\n");

    for(i=0; i<h; i++){
        for(j=0; j<w; j++){
            p = ts->anisotropy->tbrow[i]+j;
            an = ts->anisotropy->val[p];
            th = ts->thickness->val[p];

            if(th>0.0){
                bin = ROUND(ts->orientation->val[p]*ratio);
                if (DEBUG) printf("\tORIENTATION=%lf\tbin=%d\tan=%lf\n", ts->orientation->val[p], bin, an);
                hist[bin] += an;
            }
        }
    }
    if (DEBUG) printf("TSOrientationHistogram:: Percorreu a img!\n");

    //copia os valores do bin HISTOGRAMSIZE para o bin zero, pois 360=0 graus
    hist[0] += hist[HISTOGRAMSIZE];
    hist[HISTOGRAMSIZE] = 0.0; //e limpa o ultimo bin

    //Normalizacao do histograma
    sum = 0.0;
    for(i=0;i<HISTOGRAMSIZE;i++)
        sum += hist[i];
    for(i=0;i<HISTOGRAMSIZE;i++)
        hist[i] /= sum;

    return hist;
}


CImage *TSShowHistograms(float *hist1, float *hist2, int offset){
  CImage *cimg;
  float max1,max2,max;
  int col,row,col2,row2,i;
  int xsize=HISTOGRAMSIZE*2, ysize=100*2;
  
  cimg = CreateCImage(10*2+xsize,10*2+ysize);
  SetImage(cimg->C[0], 255);
  SetImage(cimg->C[1], 255);
  SetImage(cimg->C[2], 255);
  DrawCImageLineDDA(cimg, 10, ysize+10, 10+xsize, ysize+10, triplet(0, 0, 0));

  max1=max2=0.0;
  for(i=0;i<HISTOGRAMSIZE;i++){
    if(hist1[i]>max1) max1=hist1[i];
    if(hist2[i]>max2) max2=hist2[i];
  }

  if(max1==0 || max2==0) exit(1);
  max = (max1>max2)?max1:max2;

  for(i=0;i<HISTOGRAMSIZE-1;i++){
    row = 10+ysize*(max-hist1[(i+offset)%HISTOGRAMSIZE])/max;
    col = xsize*(i/((float)HISTOGRAMSIZE-1.0))+10;
    row2 = 10+ysize*(max-hist1[(i+1+offset)%HISTOGRAMSIZE])/max;
    col2 = xsize*((i+1)/((float)HISTOGRAMSIZE-1.0))+10;
    DrawCImageLineDDA(cimg, col, row, col2, row2, triplet(255, 0, 0));

    row = 10+ysize*(max-hist2[i])/max;
    col = xsize*(i/((float)HISTOGRAMSIZE-1.0))+10;
    row2 = 10+ysize*(max-hist2[i+1])/max;
    col2 = xsize*((i+1)/((float)HISTOGRAMSIZE-1.0))+10;
    DrawCImageLineDDA(cimg, col, row, col2, row2, triplet(0, 0, 255));
  }

  return cimg;
}
/***********************************************************/

/******************TSD SIMILARITY ALGORITHM*****************/
float TSHistogramMatch(float *hist1, float *hist2, int *offset){
  float *newhist;
  float *newh1,*newh2;
  int i,j,p;
  float max,correlacao;
  float score;
  float dabs,aux;
  //int maxoffset;
  //!FILE *file1,*file2,*file3;

  newhist = (float *)malloc(sizeof(float)*2*HISTOGRAMSIZE);
  newh1=newhist;
  newh2=newhist+HISTOGRAMSIZE;

  //Ajuste no histograma
  newh1[0] = (2.0*hist1[0]+hist1[HISTOGRAMSIZE-1]+hist1[1])/4.0;
  newh2[0] = (2.0*hist2[0]+hist2[HISTOGRAMSIZE-1]+hist2[1])/4.0;
  for(i=1;i<HISTOGRAMSIZE-1;i++){
    newh1[i] = (2.0*hist1[i]+hist1[i-1]+hist1[i+1])/4.0;
    newh2[i] = (2.0*hist2[i]+hist2[i-1]+hist2[i+1])/4.0;
  }
  newh1[HISTOGRAMSIZE-1] = (2.0*hist1[HISTOGRAMSIZE-1]+hist1[HISTOGRAMSIZE-2]+hist1[0])/4.0;
  newh2[HISTOGRAMSIZE-1] = (2.0*hist2[HISTOGRAMSIZE-1]+hist2[HISTOGRAMSIZE-2]+hist2[0])/4.0;

  //Correlacao
  //maxoffset = ROUND((24.0*HISTOGRAMSIZE)/180.0); // 24 graus.
  *offset=0;
  max=0.0;
  for(i=0;i<HISTOGRAMSIZE;i++){
    correlacao=0.0;
    //if(i==maxoffset) i=HISTOGRAMSIZE-maxoffset; //angulo entre -24 e 24.
    for(p=i,j=0;j<HISTOGRAMSIZE;j++,p++){
      if(p==HISTOGRAMSIZE) p=0;
      correlacao+=(newh1[p]*newh2[j]);
    }

    if(correlacao>max){
      max=correlacao;
      *offset=i;
    }
  }

  //!file1 = fopen("histogram1d.txt","w");
  //!file2 = fopen("histogram2.txt","w");
  //!file3 = fopen("histogram1.txt","w");
  dabs = 0.0;
  for(p=*offset,j=0;j<HISTOGRAMSIZE;j++,p++){
    if(p==HISTOGRAMSIZE) p=0;

    //!fprintf(file1,"%d %f\n",j,newh1[p]);
    //!fprintf(file2,"%d %f\n",j,newh2[j]);
    //!fprintf(file3,"%d %f\n",j,newh1[j]);

    aux=(newh1[p]-newh2[j]);
    aux=(aux<0.0)?(-aux):(aux);
    dabs+=aux;
  }
  score = 1.0 - dabs;

  free(newhist);
  //!fclose(file1);
  //!fclose(file2);
  //!fclose(file3);

  return score;
}
/***********************************************************/



/***********************************************************/
/* TSDIZ - Tensor Scale Descriptor with Influence Zones    */
/***********************************************************/

/****************TSDIZ EXTRACTION ALGORITHM*****************/
/* 
   Calcula a média das orientações de ts, ponderada pelas
   respectivas anisotropias, em cada área definida por
   mask na imagem img (método de Mardia).
*/
double *TSMaskedMean(TensorScale *ts, Image *img, Image *mask) {
  Curve        *meanAux;
  double       *mean;
  double       meanX, meanY, angle;
  int          n, p, LMax;

  n         = img->ncols*img->nrows;
  LMax      = MaximumValue(mask);
  meanAux   = CreateCurve(LMax);
  mean      = (double *) calloc(LMax, sizeof(double));

  for (p = 0; p < n; p++) {
    if (img->val[p] != 0) {
      if (ts->orientation->val[p] == 0.0) ts->orientation->val[p] = PI;
      meanAux->X[mask->val[p]-1] += ts->anisotropy->val[p]*sin(2.0*ts->orientation->val[p]);
      meanAux->Y[mask->val[p]-1] += ts->anisotropy->val[p]*cos(2.0*ts->orientation->val[p]);
    }
  }
  
  angle = 0.0;
  meanX = meanY = 0.0;
  for (p = 0; p < LMax; p++) {
    meanX = meanAux->X[p];
    meanY = meanAux->Y[p];
    if (meanY == 0.0) {
      if (meanX > 0.0) angle = PI/2.0;
      else angle = 3.0*PI/2.0;
    }
    else {
      angle = atan(meanX/meanY);
      if (meanY < 0.0) angle += PI;
      if ((meanX < 0.0) && (meanY > 0.0)) angle += 2.0*PI;
    }
    angle  /= 2.0;
    mean[p] = angle;
  }

  DestroyCurve(&meanAux);
  
  return mean;
}

FeatureVector1D *TSDIZ_ExtractionAlgorithm(Image *bin, int nsamples) {
  FeatureVector1D *desc=NULL;
  Image           *bin2=NULL, *cont=NULL, *segments=NULL, *aux=NULL;
  TensorScale     *ts;
  AnnImg          *aimgCont;
  double          *mean;
  int             n, p, LMaxCont, ncols, nrows, m_pairs;
  
  ncols = bin->ncols;
  nrows = bin->nrows;
  n     = ncols*nrows;

  for (p = 0; p < n; p++) {
    if (bin->val[p] != 0) bin->val[p] = 1;
  }

  aux = MBB(bin);
  bin2 = AddFrame(aux,1,0);
  DestroyImage(&aux);

  ncols = bin2->ncols;
  nrows = bin2->nrows;
  n     = ncols*nrows;

  
  aux      = LabelContPixel(bin2);
  LMaxCont = MaximumValue(aux);
  cont     = CreateImage(ncols,nrows);

  /* Dividindo o contorno em segmentos */
  segments = CreateImage(ncols, nrows);
  for (p=0; p < n; p++) {
    if (aux->val[p] != 0){
      segments->val[p] =(int)(((aux->val[p]-1)*nsamples)/LMaxCont)+1;
      cont->val[p]     = 0;
      
    }
    else cont->val[p] = 1;
  }
  DestroyImage(&aux);
  aimgCont = TSEDistTransNew(cont, segments);

  /* 
     Calculando tensor scale e as medias das orientacoes
     em cada zona de influencia formada pelos segmentos de
     contorno.
  */
  m_pairs   = 24;
  ts        = CreateBinaryTensorScale(bin2, m_pairs);
  mean      = TSMaskedMean(ts, bin2, aimgCont->label);

  /*
    Convertendo de radianos para graus
  */
  desc = CreateFeatureVector1D(nsamples);
  for (p = 0; p < nsamples; p++) {
    desc->X[p] = 180.0*mean[p]/PI;
  }
  
  DestroyAnnImg(&aimgCont);
  DestroyImage(&bin2);
  DestroyImage(&cont);
  DestroyImage(&segments);
  DestroyTensorScale(&ts);
  free(mean);
  
  return desc;

}

/***********************************************************/

/*****************TSDIZ SIMILARITY ALGORITHM****************/
double TSDIZ_SimilarityAlgorithm(FeatureVector1D *desc1, FeatureVector1D *desc2) {

  double *d1x, *d2x;
  int cont, xmax, i, index, index2, mincont, esp, angle, minangle;
  double mindist, dist1, dist2, angle_aux1, angle_aux2, dist_aux1,dist_aux2, min;

  d1x = desc1->X;
  d2x = desc2->X;

  xmax = desc1->n;
  mindist = (double)INT_MAX;
  esp = 0;
  mincont = 0;
  minangle = 0;
  
  for (cont = 0; cont < xmax; cont++) {
    for (angle = 0; angle < 180; angle++) {
      dist1  = 0.0;
      dist2  = 0.0;
      index  = cont+1;
      index2 = cont-1;
      for (i = xmax-1; i >= 0; i--) {
        index--;
        index2++;
        
        if (index < 0) index += xmax;
        if (index2 > xmax-1) index2 -= xmax;
        
        angle_aux2 = (double)angle - d2x[index2];
        angle_aux1 = d2x[index] - (double)angle;
        
        if (angle_aux1 < 0.0) angle_aux1 += 180.0;
        if (angle_aux2 < 0.0) angle_aux2 += 180.0;
        
        dist_aux1 = fabs(angle_aux1 - d1x[i]);
        dist_aux2 = fabs(angle_aux2 - d1x[i]);
        if (dist_aux2 > 90.0) dist_aux2 = 180.0 - dist_aux2;
        if (dist_aux1 > 90.0) dist_aux1 = 180.0 - dist_aux1;
        
        dist1 += dist_aux1;
        dist2 += dist_aux2;

        if ((dist1 > mindist) && (dist2 > mindist)) break;
      }
      
      min = MIN(dist1,dist2);
      if (min < mindist) {
        mindist = min;
        minangle = angle;
        if (min == dist1) {
          mincont = cont+1;
          esp = 0;
        }
        else {
          mincont = cont-1;
          esp = 1;
        }
      }
    }
  }

  return mindist;
}

/***********************************************************/

/***********************************************************/
/* TSCS - Tensor Scale Contour Salience                    */
/***********************************************************/

/****************TSCS EXTRACTION ALGORITHM******************/
/*
  Input:  a => angle in any interval
  Return: angle in interval [0.0, 180.0]
*/
double FixAngle(double a){
  int n;
  
  if(a < 0.0) a += PI;
  if(a >= PI) {
    n  = (int)(a/PI);
    a -= PI*n;
  }
  return a;
}

/* Input:  a => angle in interval [0.0, 180.0]
           b => angle in interval [0.0, 180.0]
   Return: angle in interval [0.0, 90.0]
*/
double AngularDistance(double a, double b){
  double d;
  
  d = a - b;
  if(d < 0.0)  d = -d;
  if(d > PI/2) d = PI-d;
  return d;
}

/* 
   Input:  a => angle in interval [0.0, 180.0]
           b => angle in interval [0.0, 180.0]
   Return: angle in interval [0.0, 180.0]
*/
float AngularMean(float a, float b){
  float m,d,t1,t2;
  
  d = AngularDistance(a, b);
  
  t1 = MAX(a,b);
  t2 = MIN(a,b);
  
  if((t2+PI-t1) < (t1-t2))
    m = t1 + d/2.0;
  else m = t2 + d/2.0;
  
  m = FixAngle(m);
  
  return m;
}

FeatureVector2D *TSCS_ExtractionAlgorithm(Image *bin, double threshold) {
  TensorScale     *ts=NULL, *ts2=NULL;
  FeatureVector2D *desc=NULL;
  Image           *bin2=NULL, *cont=NULL;
  Curve3D         *saliences=NULL, *saliences2=NULL, *saliences3=NULL, *influence_area=NULL;
  AnnImg          *aimg=NULL;
  int             i, j, n, p, LMaxCont, ncols, nrows, m_pairs, prox, ant, aux, dist;
  double          sumz, ni, ne;

  n     = bin->ncols*bin->nrows;
  ncols = bin->ncols;
  nrows = bin->nrows;
  
  for (i = 0; i < n; i++) {
    if (bin->val[i] != 0) bin->val[i] = 1;
  }
  
  bin2 = CopyImage(bin);
  p = bin2->tbrow[nrows-1];
  for(i = 0; i < ncols; i++) {
    bin2->val[i]   = 0;
    bin2->val[p+i] = 0;
  }
  for(i = 0; i< nrows; i++){
    p = bin2->tbrow[i];
    bin2->val[p] = 0;
    bin2->val[p+ncols-1] = 0;
  }
  
  cont     = LabelContPixel(bin2);
  aimg     = TSEDistTransNew(bin2, cont);
  LMaxCont = MaximumValue(cont);
  
  m_pairs = 24;
  ts      = CreateBinaryTensorScale(bin2, m_pairs);
  
  ts2              = (TensorScale *)malloc(sizeof(TensorScale));
  ts2->orientation = CreateDImage(ncols,nrows);
  ts2->anisotropy  = CreateDImage(ncols,nrows);
  ts2->thickness   = CreateDImage(ncols,nrows);

  /*
    Procurar elipse de maior anisotropia em cada zona de influência.
  */
  for (i = 0; i < n; i++) {
    if (bin2->val[i] != 0)  {
      if (ts2->anisotropy->val[aimg->root->val[i]] < ts->anisotropy->val[i]) {
        ts2->anisotropy->val[aimg->root->val[i]]  = ts->anisotropy->val[i];
        ts2->orientation->val[aimg->root->val[i]] = ts->orientation->val[i];
        ts2->thickness->val[aimg->root->val[i]]   = ts->thickness->val[i];
      }
    }
  }
  DestroyAnnImg(&aimg);

  saliences  = CreateCurve3D(LMaxCont);
  saliences2 = CreateCurve3D(LMaxCont);
  saliences3 = CreateCurve3D(LMaxCont);

  for (i = 0; i < n; i++) {
    if (cont->val[i] != 0)  { 
      saliences->X[cont->val[i]-1] = i%ncols;
      saliences->Y[cont->val[i]-1] = i/ncols;
      saliences->Z[cont->val[i]-1] = ts2->orientation->val[i];
    }
  }

  /*
    Se a maior anisotropia for muito pequena, pegar como orientação a média das orientações dos vizinhos.
  */
  for (i = 0; i < n; i++) {
    if (cont->val[i] != 0) {
      if (ts2->anisotropy->val[i] < 0.001) {
        prox = cont->val[i]%LMaxCont;
        ant =  cont->val[i]-2;
        if (ant < 0) ant += LMaxCont;
        saliences->Z[cont->val[i]-1] = AngularMean(saliences->Z[prox],saliences->Z[ant]);
      }
    }
  }

  /*
    Fazer a diferença da orientação calculada de um pixel do contorno para a orientação de seu pixel vizinho, no sentido horário.
    Considerar somente saliências acima do limiar.
  */
  for (i = 0; i < LMaxCont; i++) {
    prox = (i+1)%LMaxCont;
    ant =  i-1;
    if (ant < 0) ant += LMaxCont;

    saliences2->X[i] = saliences3->X[i] = saliences->X[i];
    saliences2->Y[i] = saliences3->Y[i] = saliences->Y[i];
    saliences2->Z[i] = 180.0*AngularDistance(saliences->Z[prox], saliences->Z[ant])/PI;
    if (saliences2->Z[i] > threshold) saliences3->Z[i] = saliences2->Z[i];
    else saliences3->Z[i] = 0.0;
  }

  /*
    Calcular valor das saliências achadas.
  */
  influence_area = Saliences(bin2,10);
        
  for (i = 0; i < LMaxCont; i++) {
    if (saliences3->Z[i] > 0.0) {
      p = (int)saliences3->X[i] + cont->tbrow[(int)saliences3->Y[i]];
      saliences3->Z[i] = influence_area->Z[cont->val[p]-1];
    }
  }

  SortCurve3D(saliences3, 0, (LMaxCont - 1), DECREASING);


  /*
    Verificar se existem saliencias muito proximas.
  */
  for (i = 0; i < LMaxCont; i++) {
    sumz = 0.0;
    aux  = 0;
    if (fabs(saliences3->Z[i]) > 0.0) {
      for (j = 0; j < LMaxCont; j++) {
        dist = MAX (abs((int)saliences3->X[i] - (int)saliences3->X[j]),abs((int)saliences3->Y[i] - (int)saliences3->Y[j]));
        if ((dist <= 10) && (fabs(saliences3->Z[j]) > 0.0)) {
          sumz            += saliences3->Z[j];
          saliences3->Z[j] = 0.0;
          aux++;
        }
      }
    }
    if (aux > 0) saliences3->Z[i] = sumz;
  }

  ni  = 0.0;
  ne  = 0.0;
  aux = 0;
  for (i = 0; i < LMaxCont; i++) {
    if (saliences3->Z[i] > 0.0) {
      ni += saliences3->Z[i];
      aux++;
    }
    else if (saliences3->Z[i] < 0.0) {
      ne += fabs(saliences3->Z[i]);
      aux++;
    }
  }
  for (i = 0; i < LMaxCont; i++) {
    if (saliences3->Z[i] > 0.0) 
      saliences3->Z[i] /= ni;
    else if (saliences3->Z[i] < 0.0) 
      saliences3->Z[i] /= ne;
  }

  /*
    Criar e ordenar descritor pelas distancias.
  */
  desc = CreateFeatureVector2D(aux);
  aux = 0;
  for (i = 0; i < LMaxCont; i++) {
    if (fabs(saliences3->Z[i]) > 0.0) { 
      desc->X[aux] = saliences3->Z[i];
      j = (int)saliences3->X[i] + cont->tbrow[(int)saliences3->Y[i]];
      desc->Y[aux] = ((double)(cont->val[j]-1))/LMaxCont;
      aux++;
    }
  }

  SortFeatureVector2D(desc, 0, (desc->n - 1), INCREASING);
  DescInvertXY(desc);

  DestroyCurve3D(&influence_area);
  DestroyCurve3D(&saliences);
  DestroyCurve3D(&saliences2);
  DestroyCurve3D(&saliences3);
  DestroyTensorScale(&ts);
  DestroyTensorScale(&ts2);
  DestroyImage(&bin2);
  DestroyImage(&cont);
  DestroyAnnImg(&aimg);

  return(desc);

}
/***********************************************************/
/******* FUNCOES ADICIONADAS POR MIM ***********************/
//Encontra o minimo valor entre 3
float min(float x, float y, float z) {

    if ( (x<=y) && (x<=z) ) {
        return x;
    } else if ( (y<=x) && (y<=z) ) {
        return y;
    } else if ( (z<=x) && (z<=y) ) {
        return z;
    }
    return -1;
}

//Encontra o maximo valor entre 3
float max(float x, float y, float z) {

    if ( (x>=y) && (x>=z) ) {
        return x;
    } else if ( (y>=x) && (y>=z) ) {
        return y;
    } else if ( (z>=x) && (z>=y) ) {
        return z;
    }
    return -1;
}

void RGB2HSV_tsd(CImage *RGB, ImageHSV **HSV) {

    float r, g, b;
    float minVal, maxVal, delta;
    int i, j;

    for (i = 0; i < RGB->C[0]->nrows; i++) {
        for (j = 0; j < RGB->C[0]->ncols; j++) {

            //normaliza de 0 a 1
            r = (float) RGB->C[0]->val[j+RGB->C[0]->tbrow[i]]/255;
            g = (float) RGB->C[1]->val[j+RGB->C[1]->tbrow[i]]/255;
            b = (float) RGB->C[2]->val[j+RGB->C[2]->tbrow[i]]/255;

            minVal = min(r, g, b);
            maxVal = max(r, g, b);
            delta = maxVal - minVal;

            //H
            if (delta == 0) {
                HSV[i][j].H = 0;
            } else if (maxVal==r && g>=b) {
                HSV[i][j].H = 60*((g-b)/delta);
            } else if (maxVal==r && g<b) {
                HSV[i][j].H = 60*((g-b)/delta) + 360;
            } else if (maxVal==g) {
                HSV[i][j].H = 60*((b-r)/delta) + 120;
            } else if (maxVal==b) {
                HSV[i][j].H = 60* ((r-g)/delta) + 240;
            }

            //S
            if (maxVal==0) {
                HSV[i][j].S = 0;
            } else {
                HSV[i][j].S = delta/maxVal;
            }

            //V
            HSV[i][j].V = maxVal;

            //normalizando S e V entre 0 e 255
            HSV[i][j].S *= 255;
            HSV[i][j].V *= 255;

            //colocando H de 0 a 1
            //HSV[i][j].H /= 360;

        }
    }
}

float *ReadTSD(char *filename) {
    FILE *f;
    float *tsd=NULL;

    tsd = (float*) calloc(HISTOGRAMSIZE, sizeof(float));

    if ((f = fopen(filename, "rb")) == NULL) {
        printf("erro abrindo arquivo: %s\n", filename);
        exit(0);
    }
    fread(tsd, sizeof(float), HISTOGRAMSIZE, f);
    fclose(f);
    return tsd;
}

void WriteTSD(float *hist, char *filename) {
    FILE *fout;
    if ((fout = fopen(filename, "wb")) == NULL) {
        printf("ERRO CRIANDO ARQUIVO DE FV\n");
        exit(0);
    }
    fwrite(hist, sizeof(float), HISTOGRAMSIZE, fout);
    fclose(fout);
}
