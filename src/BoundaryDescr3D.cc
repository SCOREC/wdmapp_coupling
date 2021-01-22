#include "BoundaryDescr3D.h"
#include "couplingTypes.h"
#include "importpart3mesh.h"
#include "dataprocess.h"
#include "commpart1.h"
#include "testutilities.h"
#include "sendrecv_impl.h"
#include <mpi.h>
#include <math.h>

namespace coupler {

void BoundaryDescr3D::initGeneXgc(
    const Part3Mesh3D& p3m3d,
    const Part1ParalPar3D &p1pp3d)
{
  if(preproc==true){
    nzb=p1pp3d.nzb;
    updenz=new CV**[p1pp3d.li0];
    lowdenz=new CV**[p1pp3d.li0];
    for(LO i=0;i<p1pp3d.li0;i++){
      updenz[i]=new CV*[p1pp3d.lj0];
      lowdenz[i]=new CV*[p1pp3d.lj0];
      for(LO j=0;j<p1pp3d.lj0;j++){
	updenz[i][j]=new CV[nzb];
	lowdenz[i][j]=new CV[nzb];
      }
    }

    if(p1pp3d.mype_z==0){
      if(lowpbmat==NULL){
        lowpbmat=new CV*[p1pp3d.li0];
        for(LO i=0;i<p1pp3d.li0;i++){
          lowpbmat[i]=new CV[p1pp3d.lj0];      
        } 
      }
     }else if(p1pp3d.mype_z==p1pp3d.npz-1){
      if(uppbmat==NULL){
        uppbmat=new CV*[p1pp3d.li0];
        for(LO i=0;i<p1pp3d.li0;i++){
          uppbmat[i]=new CV[p1pp3d.lj0]; 
        }
       }
     }

    uppotentz=new CV**[p3m3d.xboxinds[p1pp3d.mype_x][0]];
    lowpotentz=new CV**[p3m3d.xboxinds[p1pp3d.mype_x][0]];
    upzpart3=new double*[p3m3d.xboxinds[p1pp3d.mype_x][0]];
    lowzpart3=new double*[p3m3d.xboxinds[p1pp3d.mype_x][0]];
    for(LO i=0;i<p3m3d.xboxinds[p1pp3d.mype_x][0];i++){
      uppotentz[i]=new CV*[p3m3d.lj0/2];
      lowpotentz[i]=new CV*[p3m3d.lj0/2]; 
      upzpart3[i]=new double[nzb];
      lowzpart3[i]=new double[nzb];
      for(LO j=0;j<p3m3d.lj0/2;j++){
	uppotentz[i][j]=new CV[nzb];
	lowpotentz[i][j]=new CV[nzb];
      }
    }    

    mesh1ddens = new double[p1pp3d.lk0+2*p1pp3d.nzb];
    mesh1dpotent=new double*[p3m3d.li0];
    for(LO i=0;i<p3m3d.li0;i++){
      mesh1dpotent[i]=new double[p3m3d.mylk0[i]+2*p1pp3d.nzb];
    }

    initpbmat(p1pp3d);
  } 
}

void BoundaryDescr3D::initpbmat(const Part1ParalPar3D &p1pp3d)
{
   LO num;
   if(p1pp3d.mype_z==0){
     for(LO i=0;i<p1pp3d.li0;i++){
       num=p1pp3d.li1+i;
       for(LO j=0;j<p1pp3d.lj0;j++){  
	 lowpbmat[i][j]=exp(CV(0.0,1.0)*2.0*cplPI*double(p1pp3d.n0_global)
                        *double(j+p1pp3d.ky0_ind)*p1pp3d.q_prof[num]);
       }
     } 
   } else if(p1pp3d.mype_z==p1pp3d.npz-1){
     for(LO i=0;i<p1pp3d.li0;i++){
       num=p1pp3d.li1+i;
       for(LO j=0;j<p1pp3d.lj0;j++){
         uppbmat[i][j]=exp(-CV(0.0,1.0)*2.0*cplPI*double(p1pp3d.n0_global)
                        *double(j+p1pp3d.ky0_ind)*p1pp3d.q_prof[num]);
       }
     }
   }
}

void BoundaryDescr3D::initGemXgc(const Part3Mesh3D& p3m3d,const Part1ParalPar3D &p1pp3d)
{
  nzb=p1pp3d.nzb;
  ymeshgem=new double[p1pp3d.lj0+4];

  for(LO j=0;j<p1pp3d.lj0;j++){
    ymeshgem[2+j]=p1pp3d.y_gem[j];
  }

  ymeshgem[1]=p1pp3d.y_gem[p1pp3d.lj0-1]-p1pp3d.ly;
  ymeshgem[0]=p1pp3d.y_gem[p1pp3d.lj0-2]-p1pp3d.ly;
  ymeshgem[p1pp3d.lj0+2]=p1pp3d.y_gem[0]+p1pp3d.ly;
  ymeshgem[p1pp3d.lj0+3]=p1pp3d.y_gem[1]+p1pp3d.ly;

  thetameshgem=new double[p1pp3d.lk0+4];
  for(LO k=0;k<p1pp3d.lk0;k++){
    thetameshgem[2+k]=p1pp3d.theta[p1pp3d.lk1+k];
  }
  if(p1pp3d.mype_z!=0 && p1pp3d.mype_z!=p1pp3d.npz-1){
    thetameshgem[1]=p1pp3d.theta[p1pp3d.lk1-1];
    thetameshgem[0]=p1pp3d.theta[p1pp3d.lk1-2];
    thetameshgem[p1pp3d.lk0+2]=p1pp3d.theta[p1pp3d.lk2+1];
    thetameshgem[p1pp3d.lk0+3]=p1pp3d.theta[p1pp3d.lk2+2];
  }else if(p1pp3d.mype_z==0){
    double delta = p1pp3d.theta[1] - p1pp3d.theta[0];
    thetameshgem[1]=-delta + p1pp3d.theta[0];
    thetameshgem[0]=-2.0*delta + p1pp3d.theta[0];
    thetameshgem[p1pp3d.lk0+2]=p1pp3d.theta[p1pp3d.lk2+1];
    thetameshgem[p1pp3d.lk0+3]=p1pp3d.theta[p1pp3d.lk2+2];     
  }else {
    double delta = p1pp3d.theta[p1pp3d.lk2] - p1pp3d.theta[p1pp3d.lk2-1];
    thetameshgem[1]=p1pp3d.theta[p1pp3d.lk1-1];
    thetameshgem[0]=p1pp3d.theta[p1pp3d.lk1-2];
    thetameshgem[p1pp3d.lk0+2]=p1pp3d.theta[p1pp3d.lk2] + delta;
    thetameshgem[p1pp3d.lk0+3]=p1pp3d.theta[p1pp3d.lk2] + 2.0*delta;
  }
  /*FIXME: is the boundary buffer required?*/
  ymeshxgc = new double**[p3m3d.li0];
  for(LO i=0;i<p3m3d.li0;i++){
    ymeshxgc[i]=new double*[p3m3d.mylk0[i]];
    for(LO k=0;k<p3m3d.mylk0[i];k++)
      ymeshxgc[i][k]=new double[p3m3d.nphi];
  }
  /*FIXME: there is problems*/
  for(LO i=0;i<p3m3d.li0;i++){
    for(LO k=0;k<p3m3d.mylk0[i];k++){
      for(LO j=0;j<p3m3d.nphi;j++)
        ymeshxgc[i][k][j]=p3m3d.y_xgc[i][j][k];
    }
  }
  /*3rd order central Lagrange interpolation*/
  thflxmeshxgc=new double*[p3m3d.li0];
  for(LO i=0;i<p3m3d.li0;i++){
    thflxmeshxgc[i]=new double[p3m3d.mylk0[i]+4];
    for(LO k=0;k<p3m3d.mylk0[i];k++){
      thflxmeshxgc[i][k+2]=p3m3d.theta_flx[i][p3m3d.mylk1[i]+k];
    }
    /*for lower boundary*/
    if(p3m3d.mylk1[i]==0){
      thflxmeshxgc[i][1]=p3m3d.theta_flx[i][p3m3d.versurf[i+p1pp3d.li1]-1]-2.0*cplPI; 
      thflxmeshxgc[i][0]=p3m3d.theta_flx[i][p3m3d.versurf[i+p1pp3d.li1]-2]-2.0*cplPI;
    }else if(p3m3d.mylk1[i]==1){
      thflxmeshxgc[i][1]=p3m3d.theta_flx[i][0];
      thflxmeshxgc[i][0]=p3m3d.theta_flx[i][p3m3d.versurf[i+p1pp3d.li1]-1]-2.0*cplPI;
    }else{
      thflxmeshxgc[i][1]=p3m3d.theta_flx[i][p3m3d.mylk1[i]-1];
      thflxmeshxgc[i][0]=p3m3d.theta_flx[i][p3m3d.mylk1[i]-2];
    }
    /*for upper boundary*/
    if(p3m3d.mylk2[i]==p3m3d.versurf[i+p1pp3d.li1]-1){
      thflxmeshxgc[i][p3m3d.mylk0[i]+2]=p3m3d.theta_flx[i][0]+2.0*cplPI;
      thflxmeshxgc[i][p3m3d.mylk0[i]+3]=p3m3d.theta_flx[i][1]+2.0*cplPI;
    }else if(p3m3d.mylk2[i]==p3m3d.versurf[i+p1pp3d.li1]-2){ 
      thflxmeshxgc[i][p3m3d.mylk0[i]+2]=p3m3d.theta_flx[i][p3m3d.versurf[i+p1pp3d.li1]-1];
      thflxmeshxgc[i][p3m3d.mylk0[i]+3]=p3m3d.theta_flx[i][0]+2.0*cplPI;      
    }else{
      thflxmeshxgc[i][p3m3d.mylk0[i]+2]=p3m3d.theta_flx[i][p3m3d.mylk2[i]+1];
      thflxmeshxgc[i][p3m3d.mylk0[i]+3]=p3m3d.theta_flx[i][p3m3d.mylk2[i]+2];
    }
  }
  
  updenzgemxgc=new double**[p1pp3d.li0];
  lowdenzgemxgc=new double**[p1pp3d.li0];
  for(LO i=0;i<p1pp3d.li0;i++){
    updenzgemxgc[i]=new double*[p3m3d.nphi];
    lowdenzgemxgc[i]=new double*[p3m3d.nphi];
    for(LO j=0;j<p1pp3d.lj0;j++){
      updenzgemxgc[i][j]=new double[nzb];
      lowdenzgemxgc[i][j]=new double[nzb];
    }
  }

  uppotentzgemxgc=new double**[p3m3d.li0];
  lowpotentzgemxgc=new double**[p3m3d.li0];
  upzpart3=new double*[p3m3d.li0];
  lowzpart3=new double*[p3m3d.li0];
  for(LO i=0;i<p3m3d.li0;i++){
    uppotentzgemxgc[i]=new double*[p3m3d.nphi];
    lowpotentzgemxgc[i]=new double*[p3m3d.nphi]; 
    upzpart3[i]=new double[nzb];
    lowzpart3[i]=new double[nzb];
    for(LO j=0;j<p3m3d.nphi;j++){
      uppotentzgemxgc[i][j]=new double[nzb];
      lowpotentzgemxgc[i][j]=new double[nzb];
    }
  }   

}

} // BoundaryDescr3D.cc 
 
