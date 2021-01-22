#include "couplingTypes.h"
#include "dataprocess.h"
#include "interpoutil.h"
#include "importpart3mesh.h"
#include "commpart1.h"
#include "BoundaryDescr3D.h"

namespace coupler {

void DatasProc3D::mesh1dforDensityInterpo()
{
   LO nzb=bdesc->nzb;
   if(p1->periods[2]==1){
     if(p1->mype_z==0){
       for(LO l=0;l<nzb;l++){
	 mesh1ddens[l]=p1->pzcoords[p1->nz0-nzb+l]-2.0*cplPI;
	 mesh1ddens[p1->lk0+nzb+l]=p1->pzcoords[p1->lk2+l+1];          
       }
     } else if(p1->mype_z==p1->npz-1){
	 for(LO l=0;l<nzb;l++){
	   mesh1ddens[l]=p1->pzcoords[p1->lk1-nzb+l];
	   mesh1ddens[p1->lk0+nzb+l]=p1->pzcoords[l]+2.0*cplPI;
	 }        
     }else{
	 for(LO l=0;l<nzb;l++){
	   mesh1ddens[l]=p1->pzcoords[p1->lk1-nzb+l];
	   mesh1ddens[p1->lk0+nzb+l]=p1->pzcoords[p1->lk2+l+1];
	 }
     }
     for(LO k=0;k<p1->lk0;k++){  
       mesh1ddens[nzb+k]=p1->pzcoords[p1->lk1+k];
     }      
   }

}

void DatasProc3D::InterpoDensity3D()
{
  CV* yin;
  CV* yout;
  double* xout;
  LO nzb=bdesc->nzb;
  yin=new CV[p1->lk0+2*nzb];  
  if(preproc==true){
   for(LO i=0;i<p1->li0;i++){
      for(LO j=0;j<p1->lj0;j++){
        for(LO l=0;l<nzb;l++){
          yin[l]=bdesc->lowdenz[i][j][l];
          yin[p1->lk0+nzb+l]=bdesc->updenz[i][j][l];
        }
        for(LO k=0;k<p1->lk0;k++){  
          yin[nzb+k]=densin[i][j][k];
        }
        xout=p3->pzcoords[i];
        yout=densinterpo[i][j];
        Lag3dArray(yin,mesh1ddens,p1->lk0+2*nzb,yout,xout,p3->mylk0[i]);
      }
    }   
  }
  delete[] yin;
  yin=NULL;
}

void DatasProc3D::mesh1dforPotentialInterpo()
{
  LO nzb=bdesc->nzb;
  for(LO i=0;i<p3->li0;i++){
    for(LO l=0;l<nzb;l++){
      mesh1dpotent[i][l]=bdesc->lowzpart3[i][l];
      mesh1dpotent[i][p3->mylk0[i]+nzb+l]=bdesc->upzpart3[i][l];
    }

    for(LO j=0;j<p3->lj0/2;j++){
      for(LO k=0;k<p3->mylk0[i];k++){
	mesh1dpotent[i][k+nzb]=p3->pzcoords[i][k];
      }
    }
  }
}


void DatasProc3D::InterpoPotential3D()
{
  CV* yin;
  CV* yout;
  LO nzb=bdesc->nzb;
  if(preproc==true){
    for(LO i=0;i<p3->li0;i++){
      yin=new CV[p3->mylk0[i]+2*nzb];
      for(LO j=0;j<p3->lj0/2;j++){
        for(LO l=0;l<nzb;l++){
          yin[l]=bdesc->lowpotentz[i][j][l];
          yin[p3->mylk0[i]+nzb+l]=bdesc->uppotentz[i][j][l];
        }
        for(LO k=0;k<p3->mylk0[i];k++){
          yin[k+nzb]=potentinterpo[i][j][k];
        }
        yout=potentpart1[i][j];
        Lag3dArray(yin,mesh1dpotent[i],p3->mylk0[i]+2*nzb,yout,p1->pzp,p1->lk0);
     }

     delete[] yin; 
    }
  }
}


void gemXgcDatasProc3D::interpoDensityAlongZ(double*** box)
{
  bool debug;
  double* yin;
  double* yout;
  double* xout;
  LO nzb=bdesc->nzb;
  yin=new double[p1->lk0+2*nzb];
  if(preproc==true){
//   printf("reach interpo\n");
   MPI_Barrier(MPI_COMM_WORLD);  
// printf("after reach interpo\n");

   for(LO i=0;i<p1->li0;i++){
      for(LO j=0;j<p1->lj0;j++){
//     if(i == 0) printf("mype0: %d, j0: %d \n", p1->mype, j);
//        MPI_Barrier(MPI_COMM_WORLD);

	for(LO l=0;l<nzb;l++){
          yin[l]=bdesc->lowdenzgemxgc[i][j][l];
          yin[p1->lk0+nzb+l]=bdesc->updenzgemxgc[i][j][l];
        }

//	if(i == 0) printf("mype1: %d, j1: %d", p1->mype, j);	
//        MPI_Barrier(MPI_COMM_WORLD);
        for(LO k=0;k<p1->lk0;k++){
          yin[nzb+k]=densCpl[i][j][k];
        }

//	if(i == 0) printf("mype2: %d, j2: %d", p1->mype, j);
//  	MPI_Barrier(MPI_COMM_WORLD);
          xout=p3->theta_geo[i];

	  yout=box[i][j];

	  Lag3dArray(yin,bdesc->thetameshgem,p1->lk0+2*nzb,yout,xout,p3->mylk0[i]);
  
	  debug = false;
	  if (debug) {
            if (p1->mype == 1) {
	      for (LO k=0; k< p1->lk0+2*nzb; k++) printf("i: %d, j: %d, k: %d, xin: %f,  yin: %f \n",
		i, j, k, bdesc->thetameshgem[k], yout[k]);
              for (LO k=0; k< p3->mylk0[i]; k++) printf("i: %d, j: %d, k: %d, xout: %f,  yout: %f \n",
	        	i, j, k, xout[k], yout[k]);
	    }
	  }
//  	  if(i == 0) printf("mype: %d, j: %d \n", p1->mype, j);
//          MPI_Barrier(MPI_COMM_WORLD);

        }
//         if(i == 0) printf("first i=0, mype: %d \n", p1->mype);
//         MPI_Barrier(MPI_COMM_WORLD);
 
    }
  }
  delete[] yin;
  yin=NULL;
}

void gemXgcDatasProc3D::interpoDensityAlongY()
{
  bool debug = false;
  double* yin;
  double* xin;
  double* yout = new double[p1->nphi];
  double* xout;
  LO nzb=bdesc->nzb;
  yin=new double[p1->lj0+2*nzb];
  xin=new double[p1->lj0+2*nzb];
  for(LO i=0;i<p1->li0;i++){
    for(LO k=0;k<p3->mylk0[i];k++){
      for(LO j=0;j<p1->lj0;j++){
        yin[nzb+j]=densinterone[i][j][k]; 
      }
      for(LO b=0;b<nzb;b++){
        yin[b]=densinterone[i][p1->lj0+b-nzb][k];
        yin[p1->lj0+nzb-1+b]=densinterone[i][b][k];
      }
      xout=bdesc->ymeshxgc[i][k];
      Lag3dArray(yin,bdesc->ymeshgem,p1->lj0+2*nzb,yout,xout,p1->nphi); 
      for(LO j=0;j<p1->nphi;j++){
        densintertwo[i][j][k]=yout[j];
        if (debug) {
	  if (p1->mype == 1) {
	    printf("k:%d, j:%d, densintertwo: %f, densinterone: %f \n", k, j, densintertwo[i][j][k], densinterone[i][j][k]);
	  }
        }
      }      
    }
  }
}

void gemXgcDatasProc3D::InterpoPotential3DAlongZ(double*** boxin, double*** boxout)
{
  double* yin;
  double* yout;
  LO nzb=bdesc->nzb;
  if(preproc==true){
    for(LO i=0;i<p1->li0;i++){
      yin=new double[p3->mylk0[i]+2*nzb];
      for(LO j=0;j<p1->lj0;j++){
        for(LO l=0;l<nzb;l++){
          yin[l]=bdesc->lowpotentzgemxgc[i][j][l];
          yin[p3->mylk0[i]+nzb+l]=bdesc->uppotentzgemxgc[i][j][l];
        }
        for(LO k=0;k<p3->mylk0[i];k++){
          yin[k+nzb]=boxin[i][j][k];
        }
        yout=boxout[i][j];
        Lag3dArray(yin,bdesc->thflxmeshxgc[i],p3->mylk0[i]+2*nzb,yout,p1->thflx[i],p1->lk0);
     }
    }
  }

}


} /*interpo.cc*/
