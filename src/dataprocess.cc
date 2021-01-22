#include "dataprocess.h"
#include "importpart3mesh.h"
#include "commpart1.h"
#include "BoundaryDescr3D.h"
#include "sendrecv_impl.h"
#include <cassert>
#include <cmath>

namespace coupler {

DatasProc3D::DatasProc3D(const Part1ParalPar3D* p1pp3d,
    const Part3Mesh3D* p3m3d,
    const BoundaryDescr3D* bdesc_,
    bool pproc,
    TestCase test_case,
    bool ypar,
    int nummode)
  : preproc(pproc),
    testcase(test_case),
    yparal(ypar)
 {
    p1=p1pp3d;
    p3=p3m3d;
    bdesc=bdesc_;
    init();
    AllocDensityArrays();
    AllocPotentArrays();
    AllocBufferForIntepo();
    AllocMatXYZtoPlane();
    mesh1dforDensityInterpo();
    Initmattoplane(); 
    Prepare_mats_from_planes();
    if(testcase==TestCase::t0) {
      TestInitPotentAlongz(p3m3d, p1->npy, nummode);
    }
  }

void DatasProc3D::init()
{
  if(preproc==true){
    if(yparal==true){
      if(p1->li0%p1->npy==0){
        part1li0=p1->li0/p1->npy;
        part3li0=part1li0;
      } else{
        if(p1->mype_y==p1->npy-1){
          part1li0=p1->li0%p1->npy;
          part3li0=part1li0;
        }  else{
          part1li0=p1->li0%p1->npy;
          part3li0=part1li0;
        }
      }
      part1lj0=2*p1->ny0;
      part3lj0=p1->ny0;   // here may need rethinking.
    } else{
      part1li0=p1->li0;
      part3li0=p1->li0;
      part1lj0=2*p1->ny0;
      part3lj0=p1->ny0;   // here may need rethinking.
    }
  }
  sum=0;
  for(LO i=0;i<p3->li0;i++)  sum+=p3->mylk0[i];
 
}

void DatasProc3D::AllocDensityArrays()
{
  if(yparal==false){
    densin=new CV**[p1->li0];
    for(LO i=0;i<p1->li0;i++){
      densin[i]=new CV*[p1->lj0];
      for(LO j=0;j<p1->lj0;j++)
        densin[i][j]=new CV[p1->lk0];
    }
    densinterpo=new CV**[p1->li0];
    for(LO i=0;i<p1->li0;i++){
      densinterpo[i]=new CV*[p1->lj0];
      for(GO j=0;j<p1->lj0;j++){
        densinterpo[i][j]=new CV[p3->mylk0[i]];
      }
    }

   densintmp=new CV[p1->lj0];
   densouttmp=new double[p1->lj0*2];

   denspart3=new double**[p3->li0];
   for(LO i=0;i<p3->li0;i++){
     denspart3[i]=new double*[p3->lj0];
     for(LO j=0; j<p3->lj0; j++)
        denspart3[i][j]=new double[p3->mylk0[i]];
   }

   densTOpart3=new double**[p3->li0];
   for(LO i=0;i<p3->li0;i++){
     densTOpart3[i]=new double*[p3->lj0];
     for(LO j=0; j<p3->lj0; j++)
        densTOpart3[i][j]=new double[p3->mylk0[i]];
   }
   denssend = new double[p3->blockcount*p3->lj0];   
 } 
}


void DatasProc3D::AllocPotentArrays()
{ 
  if(yparal==false){
    potentin=new double**[p3->li0];
    for(LO i=0;i<p3->li0;i++){
      potentin[i]=new double*[p1->y_res_back];
      for(LO j=0;j<p1->y_res_back;j++)
        potentin[i][j]=new double[p1->lk0];
    }

    potentintmp=new double[p1->y_res_back];
    potentouttmp=new CV[p1->y_res_back/2+1];
 
    potentinterpo=new CV**[p3->li0];
    for(LO i=0;i<p3->li0;i++){
      potentinterpo[i]=new CV*[p3->lj0/2];
      for(LO j=0;j<p3->lj0/2;j++)
        potentinterpo[i][j]=new CV[p3->mylk0[i]];
    }
  
    potentpart1=new CV**[p1->li0];
    for(LO i=0;i<p1->li0;i++){
      potentpart1[i]=new CV*[p1->lj0];
      for(LO j=0;j<p1->lj0;j++){
        potentpart1[i][j]=new CV[p1->lk0];
      }
    }
 
   potentsend = new CV[p1->blockcount*p1->lj0];
  
  }
}

void DatasProc3D::AllocBufferForIntepo()
{
    mesh1ddens = new double[p1->lk0+2*p1->nzb];
    mesh1dpotent=new double*[p3->li0];
    for(LO i=0;i<p3->li0;i++){
      mesh1dpotent[i]=new double[p3->mylk0[i]+2*p1->nzb];
    }
}

void DatasProc3D::AllocMatXYZtoPlane()
{
   mattoplane=new double***[p3->li0];
   for(LO i=0;i<p3->li0;i++){
     mattoplane[i] = new double**[p1->n_cuts];
     for(LO j=0;j<p1->n_cuts;j++){
       mattoplane[i][j]=new double*[p3->lj0];
       for(LO k=0;k<p3->lj0;k++){
	 mattoplane[i][j][k]=new double[p3->mylk0[i]];
       }    
     }
   }  

  mat_to_plane=new CV***[p3->li0];
   for(LO i=0;i<p3->li0;i++){
     mat_to_plane[i] = new CV**[p1->n_cuts];
     for(LO j=0;j<p1->n_cuts;j++){
       mat_to_plane[i][j]=new CV*[p3->lj0];
       for(LO k=0;k<p3->lj0;k++){
         mat_to_plane[i][j][k]=new CV[p3->mylk0[i]];
       }
     }
   }

  mat_from_weight=new double***[p1->li0];
  mat_from_ind_plane=new LO***[p1->li0];
  mat_from_ind_n=new LO***[p1->li0];
  for(LO i=0;i<p1->li0;i++){
    mat_from_weight[i]=new double**[p1->y_res_back];
    mat_from_ind_plane[i]=new LO**[p1->y_res_back];
    mat_from_ind_n[i]=new LO**[p1->y_res_back];
    for(LO j=0;j<p1->y_res_back;j++){
      mat_from_weight[i][j]=new double*[p1->lk0];
      mat_from_ind_plane[i][j]=new LO*[p1->lk0];
      mat_from_ind_n[i][j]=new LO*[p1->lk0];
      for(LO k=0;k<p1->lk0;k++){
        mat_from_weight[i][j][k]=new double[4];
        mat_from_ind_plane[i][j][k]=new LO[2];
        mat_from_ind_n[i][j][k]=new LO[4];
      }
    }
  }
}

void DatasProc3D::DistriPotentRecvfromPart3(const Array2d<double>* fieldfromXGC)
{ 
  double** tmp;
  tmp = new double*[p3->lj0];
  for(LO j=0;j<p3->lj0;j++){
    tmp[j] = new double[p3->blockcount]; 
  }
  double* array;
  array = fieldfromXGC->data();
  for(LO j=0;j<p3->lj0;j++){
    for(GO i=0;i<p3->blockcount;i++)
      tmp[p3->lj0-j-1][i]=array[j*p3->blockcount+i]/p1->norm_fact_field;
  }
  LO xl=0; 
  GO sumbegin=0;       
  GO numnode=0;  
  double sum_in;
  bool debug=false;
  for(LO i=0;i<p3->li0;i++){
    xl=p3->li1+i;
    double** datain= new double*[p3->lj0];   
    for(LO j=0;j<p3->lj0;j++){
      numnode=sumbegin;  
      datain[j]=new double[p3->versurf[xl]];    
      for(LO m=0;m<p3->versurf[xl];m++){
        datain[j][m]=tmp[j][numnode];
        numnode=numnode+1; 
      }
      reshuffleforward(datain[j],p3->nstart[xl],p3->versurf[xl]);
    }
    sumbegin+=p3->versurf[xl];
    
    if(debug){
      if(p1->mype==0){
        sum_in=0.0;
	printminmax2d(datain,p3->lj0,p3->versurf[xl],p1->mype,"datain",i);   
        printSumm2D(datain,p3->lj0,p3->versurf[xl],sum_in,p1->mype,"datain",i);
      }
    }
 
    for(LO j=0;j<p1->y_res_back;j++){
      for(LO k=0;k<p1->lk0;k++){
        potentin[i][j][k]=
        +datain[mat_from_ind_plane[i][j][k][0]][mat_from_ind_n[i][j][k][0]]*mat_from_weight[i][j][k][0]
        +datain[mat_from_ind_plane[i][j][k][0]][mat_from_ind_n[i][j][k][1]]*mat_from_weight[i][j][k][1]
        +datain[mat_from_ind_plane[i][j][k][1]][mat_from_ind_n[i][j][k][2]]*mat_from_weight[i][j][k][2]
        +datain[mat_from_ind_plane[i][j][k][1]][mat_from_ind_n[i][j][k][3]]*mat_from_weight[i][j][k][3];
      }
    }

    for(LO j=0;j<p3->lj0;j++){
      free(datain[j]);
    }
    free(datain);
    datain=NULL;
  }  
  assert(sumbegin==p3->blockcount); 

   if(debug){
     printminmax3d(potentin,p3->li0,p1->y_res_back,p1->lk0,p1->mype,"potentin",0);
   }
    

// Fourier transform  
  RealdataToCmplxdata3D();
  for(LO j=0;j<p3->lj0;j++)
    free(tmp[j]);
  free(tmp);
  tmp=NULL;
 } 



//Distribute the sub global potential  2d array received from part3 and reorder the sub 2darray.  
void DatasProc3D::oldDistriPotentRecvfromPart3(const Array2d<double>* fieldfromXGC)
{ 
  double** tmp;
  tmp = new double*[p3->lj0];
  for(LO j=0;j<p3->lj0;j++){
    tmp[j] = new double[p3->blockcount]; 
  }
  double* array;
  array = fieldfromXGC->data();
  for(LO j=0;j<p3->lj0;j++){
    for(GO i=0;i<p3->blockcount;i++)
      tmp[j][i]=array[j*p3->blockcount+i];
  }
  double** subtmp;
  for(LO j=0;j<p3->lj0;j++){
    GO sumbegin=0;       
    subtmp = new double*[p3->li0];
    LO xl=0;
    for(LO i=0;i<p3->li0;i++){
      xl=p3->li1+i;
      subtmp[i]=new double[p3->versurf[xl]];
      for(LO m=0;m<p3->versurf[xl];m++){
        subtmp[i][m]=tmp[j][sumbegin];
        sumbegin=sumbegin+1; 
      }
      reshuffleforward(subtmp[i],p3->nstart[xl],p3->versurf[xl]);
      for(LO k=0;k<p3->mylk0[i];k++){
        potentin[i][j][k]=subtmp[i][p3->mylk1[i]+k];
      }
      free(subtmp[i]);
    }
     assert(sumbegin==p3->blockcount); 
  }
  free(subtmp);
  for(LO j=0;j<p3->lj0;j++)
    free(tmp[j]);
  free(tmp);
 } 

// Assemble the potential sub2d array in each process into a bigger one, which is straightforwardly transferred by
// the adios2 API from coupler to Part1.
void DatasProc3D::AssemPotentSendtoPart1()
{
  LO* recvcount = new LO[p1->npz];
  LO* rdispls = new LO[p1->npz];

  MPI_Datatype mpitype = getMpiType(LO());      
  MPI_Allgather(&p1->lk0,1,mpitype,recvcount,1,mpitype,p1->comm_z); 
  rdispls[0]=0;
  for(LO i=1;i<p1->npz;i++){
    rdispls[i]=rdispls[i-1]+recvcount[i-1];
  }

  CV* tmp = new CV[p1->nz0]; 
  CV* blocktmp = new CV[p1->blockcount];

  GO sumbegin;
  for(LO j=0;j<p1->lj0;j++){ 
    for(GO h=0;h<p1->blockcount;h++){
      blocktmp[h] = CV({0.0,0.0});
    } 
    for(LO i=0;i<p1->li0;i++){
      sumbegin=0;
      for(LO h=0;h<i;h++){
        sumbegin+=(GO)p1->nz0;
      }     
      for(LO h=0;h<p1->nz0;h++){
        tmp[h]=CV({0.0,0.0});
      } 
      MPI_Allgatherv(potentpart1[i][j],p1->lk0,MPI_CXX_DOUBLE_COMPLEX,tmp,recvcount,rdispls,
                    MPI_CXX_DOUBLE_COMPLEX,p1->comm_z);    
      for(LO m=0;m<p1->nz0;m++){
        blocktmp[sumbegin+m]=tmp[m];
      }      
    }    
    assert(sumbegin+p1->nz0==p1->blockcount); 
    for(GO h=0;h<p1->blockcount;h++){
      potentsend[j*p1->blockcount+h] = blocktmp[h]; 	    
    }
  }
  delete[] blocktmp; 
  delete[] tmp,recvcount,rdispls;
  blocktmp=NULL;
  tmp=NULL;
  recvcount=NULL;
  rdispls=NULL;
}

////Distribute the subglobal density  2d array received from part1 to the processes.
void DatasProc3D::DistriDensiRecvfromPart1(const Array2d<CV>* densityfromGENE)
{
  CV** tmp; 
  tmp = new CV*[p1->lj0];
  for(LO j=0;j<p1->lj0; j++){
    tmp[j] = new CV[p1->blockcount];
  }
  CV* array = densityfromGENE->data();
  for(LO j=0;j<p1->lj0;j++){
    for(GO i=0;i<p1->blockcount;i++){
      tmp[j][i]=array[j*p1->blockcount+i];
    }
  }
  for(LO i=0;i<p1->li0;i++){
    for(LO j=0;j<p1->lj0;j++){
      for(LO k=0;k<p1->lk0;k++)  densin[i][j][k]=CV(0.0,0.0);
    }
  }
  
  CV** blocktmp= new CV*[p1->li0];
  for(LO i=0;i<p1->li0;i++)
    blocktmp[i]=new CV[p1->nz0]; 

  for(LO j=0;j<p1->lj0;j++){
    for(LO i=0;i<p1->li0;i++){
      GO sumbegin=0;
      for(LO h=0;h<i;h++){
        sumbegin+=(GO)p1->nz0;
      }
      for(LO m=0;m<p1->nz0;m++){
        blocktmp[i][m]=tmp[j][sumbegin+m];
      }
      for(LO k=0;k<p1->lk0;k++){
        densin[i][j][k]=blocktmp[i][p1->lk1+k];
      }
    }
  }
   for(LO i=0;i<p1->li0;i++){
     free(blocktmp[i]);
   }
   free(blocktmp);
   blocktmp=NULL;
   for(LO i=0;i<p1->lj0;i++){
     free(tmp[i]);
   }  
   free(tmp);
   tmp=NULL;
}

// Assemble the density sub2d array in each process into a global one, which is straightforwardly transferred by
// the adios2 API from coupler to Part3.
void DatasProc3D::oldAssemDensiSendtoPart3(BoundaryDescr3D& bdesc)
{
  zDensityBoundaryBufAssign(densin,bdesc);
  InterpoDensity3D();
  CmplxdataToRealdata3D();
  DensityToPart3();
  
  LO* recvcount = new LO[p1->npz];
  LO* rdispls = new LO[p1->npz];
  double* blocktmp = new double[p3->blockcount];

  for(LO j=0;j<p3->lj0;j++){
    LO xl=0;
    for(GO h=0;h<p3->blockcount;h++){
      blocktmp[h] = 0.0;
    }

    for(LO h=0;h<p1->npz;h++){
      recvcount[h]=0;
      rdispls[h]=0;
    }

    for(LO i=0;i<p1->li0;i++){
      MPI_Datatype mpitype = getMpiType(LO());      
      MPI_Allgather(&p3->mylk0[i],1,mpitype,recvcount,1,mpitype,p1->comm_z); 
      rdispls[0]=0;
      for(LO k=1;k<p1->npz;k++){
	rdispls[k]=rdispls[k-1]+recvcount[k-1];
      }

      xl=p1->li1+i;    
      double* tmp = new double[p3->versurf[xl]];
      double* tmp_one;
      tmp_one=denspart3[i][j];
      MPI_Allgatherv(tmp_one,p3->mylk0[i],MPI_DOUBLE,tmp,recvcount,rdispls,
                    MPI_DOUBLE,p1->comm_z);    
 
      reshufflebackward(tmp,p3->nstart[xl],p3->versurf[xl]);
      GO sumbegin=0;
      for(LO h=0;h<i;h++){
        sumbegin+=GO(p3->versurf[h+p3->li1]);
      } 
      for(LO m=0;m<p3->versurf[xl];m++){
        blocktmp[sumbegin+m]=tmp[m];
      }    
     if(i==p1->li0-1){
        assert((sumbegin+(GO)p3->versurf[xl]) == p3->blockcount);
      }

      free(tmp); 
    }
   
    for(GO h=0;h<p3->blockcount;h++){
        denssend[j*p3->blockcount+h] = blocktmp[h];
    } 
  }
  free(recvcount);
  free(rdispls);
  free(blocktmp);
}

// Assemble the density sub2d array in each process into a global one, which is straightforwardly transferred by
// the adios2 API from coupler to Part3.
void DatasProc3D::AssemDensiSendtoPart3(BoundaryDescr3D& bdesc)
{
  double*** tmpmat = new double**[p3->li0];
    for(LO i=0;i<p3->li0;i++){
      tmpmat[i]=new double*[p3->lj0];
      for(LO j=0;j<p3->lj0;j++){
        tmpmat[i][j]=new double[p3->mylk0[i]];
        for(LO k=0;k<p3->mylk0[i];k++){
          tmpmat[i][j][k]=0.0;
        }
      }
    }

  zDensityBoundaryBufAssign(densin,bdesc);

  InterpoDensity3D(); 

  bool debug=false;
  if(debug){
    printminmax(densinterpo,p1->li0,p1->lj0,p3->mylk0,p1->mype,"densinterpo",0);
    MPI_Barrier(MPI_COMM_WORLD);
    CV sum=CV(0.0,0.0);
    printSumm3D(densinterpo,p1->li0,p1->lj0,p3->mylk0,sum,
    MPI_COMM_WORLD,"densinterpo",0);
  }

// don't understand the following operation  
  
  CV*** loc_data=new CV**[p1->li0];
  for(LO i=0;i<p1->li0;i++){
    loc_data[i]=new CV*[p3->lj0];
    for(LO j=0;j<p3->lj0;j++) loc_data[i][j]=new CV[p3->mylk0[i]];
  }
  CV tmp1;

  for(LO i=0;i<p1->li0;i++){
    for(LO j=0;j<p1->lj0;j++){
      for(LO k=0;k<p3->mylk0[i];k++){
        loc_data[i][j][k]= densinterpo[i][j][k];  //FIXME: Here pointer is better      
      }
      if(j>0){
        for(LO k=0;k<p3->mylk0[i];k++){
          loc_data[i][p3->lj0-j][k]=std::conj(densinterpo[i][j][k]);
        }
      }
    } 
    for(LO j=0;j<p1->n_cuts;j++){
      for(LO k=0;k<p3->mylk0[i];k++){   
        tmp1=CV(0.0,0.0);
        for(LO h=0;h<p3->lj0;h++){
          tmp1+=mat_to_plane[i][j][h][k]*loc_data[i][h][k];
        }       
        tmpmat[i][j][k]+=tmp1.real();
      }
    }     
  }

  for(LO i=0;i<p1->li0;i++){
    for(LO j=0;j<p3->lj0;j++) free(loc_data[i][j]);
    free(loc_data[i]);
  }
  free(loc_data);

//don't understand the above operation   

  LO* recvcount = new LO[p1->npz];
  LO* rdispls = new LO[p1->npz];
  double* blocktmp = new double[p3->blockcount];

  double** tmp=new double*[p3->li0];
  for(LO i=0;i<p3->li0;i++){
    tmp[i]=new double[p3->versurf[p3->li1+i]];
  }

  GO sumbegin;
  LO xl;
  LO num;

  for(LO j=0;j<p3->lj0;j++){
    xl=0;
    for(GO h=0;h<p3->blockcount;h++){
      blocktmp[h] = 0.0;
    }
    for(LO i=0;i<p3->li0;i++){
      MPI_Datatype mpitype = getMpiType(LO());      
      for(LO h=0;h<p1->npz;h++){
        recvcount[h]=0;
        rdispls[h]=0;
      }
      MPI_Allgather(&p3->mylk0[i],1,mpitype,recvcount,1,mpitype,p1->comm_z); 
      rdispls[0]=0;
      for(LO k=1;k<p1->npz;k++){
	rdispls[k]=rdispls[k-1]+recvcount[k-1];
      }
 
      xl=p1->li1+i;   

      debug=false;
      if(debug){
	num=0;
	for(LO h=0;h<p1->npz;h++){
	  num+=recvcount[h];
	}	   
	std::cout<<"num versurf[xl]="<<num<<" "<<p3->versurf[xl]<<'\n';
	assert(num==p3->versurf[xl]);
      }     

      MPI_Allgatherv(tmpmat[i][j],p3->mylk0[i],MPI_DOUBLE,tmp[i],recvcount,rdispls,
                    MPI_DOUBLE,p1->comm_z);    
      reshufflebackward(tmp[i],p3->nstart[xl],p3->versurf[xl]);

      sumbegin=0;
      for(LO h=0;h<i;h++){
        sumbegin+=GO(p3->versurf[h+p3->li1]);
      } 
      for(LO m=0;m<p3->versurf[xl];m++){
        blocktmp[sumbegin+m]=tmp[i][m];
      }    
     if(i==p1->li0-1){
        assert((sumbegin+(GO)p3->versurf[xl]) == p3->blockcount);
      }
    }   
    for(GO h=0;h<p3->blockcount;h++){
        denssend[j*p3->blockcount+h] = blocktmp[h]*p1->norm_fact_dens;
    }
  }
 
  free(recvcount);
  free(rdispls);
  free(blocktmp);
  recvcount=NULL;
  rdispls=NULL;
  blocktmp=NULL;

  for(LO i=0;i<p3->li0;i++) free(tmp[i]);
  free(tmp);
  tmp=NULL;  

  for(LO i=0;i<p3->li0;i++){
    for(LO j=0;j<p3->lj0;j++){
      free(tmpmat[i][j]);
    }
    free(tmpmat[i]);
  }
  free(tmpmat);
  tmpmat=NULL;
}

//I dont's understand the function of the following matrix.
void DatasProc3D::oldInitmattoplane()
{
  double y_cut;
  LO tmp_ind;
  LO ind_l_tmp;
  LO ind_h_tmp;
  for(LO i=0;i<p3->li0;i++){
    for(LO k=0;k<p3->mylk0[i];k++){
      for(LO j=0;j<p3->lj0;j++){
        y_cut=p1->C_y[0]*(p1->q_prof[i+p3->li1]*p3->pzcoords[i][k]-p1->phi_cut[j])/p1->dy;
        y_cut=remainder(remainder(y_cut,double(p1->y_res))+double(p1->y_res),double(p1->y_res));
      
        tmp_ind=LO(y_cut);
        ind_l_tmp=remainder(remainder(tmp_ind,p1->y_res)+p1->y_res,p1->y_res);
        ind_h_tmp=remainder(remainder(tmp_ind+1,p1->y_res)+p1->y_res,p1->y_res);

        mattoplane[i][j][ind_h_tmp][k]=y_cut-double(tmp_ind);
        mattoplane[i][j][ind_l_tmp][k]=1.0-(y_cut-double(tmp_ind));
      }
    }
  }
}


void DatasProc3D::Initmattoplane()
{
  for(LO i=0;i<p3->li0;i++){
    for(LO o=0;o<p1->n_cuts;o++){
      for(LO j=0;j<p1->lj0;j++){
        for(LO k=0;k<p3->mylk0[i];k++){
         if(j==0){
            mat_to_plane[i][o][j][k]=CV(1.0,0.0);
          }else{
            mat_to_plane[i][o][j][k]=exp(CV(0,1.0)*double(j*p1->n0_global)*(p1->q_prof[i+p3->li1]
            *p3->pzcoords[i][k]-p1->L_tor*double(o+1)));
            mat_to_plane[i][o][p3->lj0-j][k]=exp(-CV(0,1.0)*double(j*p1->n0_global)*(p1->q_prof[i+p3->li1]
            *p3->pzcoords[i][k]-p1->L_tor*double(o+1)));
          } 
        }
     }
   }
 }
 bool debug=false;
 if(debug){
  CV sum=CV(0.0,0.0);
  for(LO i=0;i<p3->li0;i++){
    for(LO o=0;o<p1->n_cuts;o++){
      for(LO j=0;j<p1->lj0;j++){
        for(LO k=0;k<p3->mylk0[i];k++){
          sum=sum+mat_to_plane[i][o][j][k];
        }
      }
    }
  }
  std::cout<<"sum of mat_to_plane, mype_x="<<p1->mype_x<<" "<<sum<<'\n';
}

}


//The function of this routines is not clear so far.
void DatasProc3D::DensityToPart3()
{
  for(LO i=0;i<p3->li0;i++){
    for(LO k=0;k<p3->mylk0[i];k++){
      for(LO j=0;j<p3->lj0;j++){
        double tmp=0.0;
        for(LO l=0;l<p3->lj0;j++){
          tmp+=mattoplane[i][j][l][k]*denspart3[i][l][k];
        }
        densTOpart3[i][j][k]=tmp;
      }
    }
  }
}

// not very clear about the function of this routine
void DatasProc3D::Prepare_mats_from_planes()
{
  double dphi=2.0*cplPI/double(p1->n0_global*p1->n_cuts);
  double* phi_l=new double[p1->n_cuts];
  for(int i=0;i<p1->n_cuts;i++){
    phi_l[i]=double(i)*dphi;
  }
  double Ly=2.0*cplPI/(double(p1->n0_global)*p1->rhostar*p1->minor_r)*std::abs(p1->C_y[1]);
  double dy_inv=Ly/double(p1->y_res_back);
  
  double q,y,chi_red,phi,phi_red,chi_red_l,chi_red_r,phi_red_l,phi_red_r,dist_phi,
         dist_l,dist_r,w_plane_left,w_plane_right,chi_l, chi_u,dchi;
  LO count_l,count_r,ipl_l,ipl_r,ind_l,ind_u;

  for(int i=0;i<p3->li0;i++){
    q=p1->q_prof[i+p3->li1];
    double* tmp = new double[p3->versurf[i+p3->li1]]; 
    for(int j=0;j<p1->y_res_back;j++){
      y=double(j)*dy_inv;
      for(int k=0;k<p1->lk0;k++){
        chi_red=p1->pzcoords[p1->lk1+k];
        phi=q*chi_red-(y/p1->C_y[i+p3->li1])*(p1->rhostar*p1->minor_r);

        phi_red=remainder(phi,2.0*cplPI/double(p1->n0_global));
        if(phi_red<0) phi_red=2.0*cplPI/double(p1->n0_global)+phi_red;
        count_l=int((phi-phi_red)/(2.0*cplPI/double(p1->n0_global)));
        count_r=count_l;
        ipl_l=int(phi_red/dphi);
        ipl_r=ipl_l+1;       

        if(ipl_r==p1->n_cuts){
          ipl_r=0;
          count_r=count_r+1;
        }

        chi_red_l=(y/p1->C_y[i+p3->li1]*(p1->rhostar*p1->minor_r)+phi_l[ipl_l]
               +count_l*2.0*cplPI/double(p1->n0_global))/q;
        chi_red_r=(y/p1->C_y[i+p3->li1]*(p1->rhostar*p1->minor_r)+phi_l[ipl_r]
               +count_r*2.0*cplPI/double(p1->n0_global))/q;   

        chi_red_l=remainder(chi_red_l+cplPI,2.0*cplPI);
        if(chi_red_l<0)  chi_red_l=2.0*cplPI+chi_red_l;
        chi_red_l=chi_red_l-cplPI;
       
        chi_red_r=remainder(chi_red_r+cplPI,2.0*cplPI);
        if(chi_red_r<0) chi_red_r=2.0*cplPI+chi_red_r;
        chi_red_r=chi_red_r-cplPI;        

        phi_red_l=phi_l[ipl_l];
        phi_red_r=phi_l[ipl_r];

        mat_from_ind_plane[i][j][k][0]=ipl_l;
        mat_from_ind_plane[i][j][k][1]=ipl_r;

        dist_phi=sqrt(pow(phi_red_l-phi_red_r,2)+pow(chi_red_l-chi_red_r,2));
        dist_l=sqrt(pow(phi_red-phi_red_r,2)+pow(chi_red-chi_red_r,2));
        dist_r=sqrt(pow(phi_red-phi_red_l,2)+pow(chi_red-chi_red_l,2));        
        w_plane_left=dist_l/dist_phi;
        w_plane_right=dist_r/dist_phi;

        //left_plane 
        for(LO m=0;m<p3->versurf[i+p3->li1];m++){
          tmp[m]=std::abs(p3->zcoordsurf[i][m]-chi_red_l);
        } 

        ind_u=minloc(tmp,p3->versurf[i+p3->li1]);

        mat_from_ind_n[i][j][k][0]=ind_u;
        chi_u=p3->zcoordsurf[i][ind_u];
 
        if((chi_red_l-chi_u)>0){
          ind_l=ind_u+1;
          if(ind_l>p3->versurf[i+p3->li1]-1){
            ind_l=0;
            chi_l=p3->zcoordsurf[i][ind_l]+2.0*cplPI;
          }else{
            chi_l=p3->zcoordsurf[i][ind_l];
          }
          mat_from_ind_n[i][j][k][1]=ind_l;
        } else{
          ind_l=ind_u-1;
          if(ind_l<0){
            ind_l=p3->versurf[i+p3->li1]-1;
            chi_l=p3->zcoordsurf[i][ind_l]-2.0*cplPI;
          } else{
            chi_l=p3->zcoordsurf[i][ind_l];
          }
          mat_from_ind_n[i][j][k][1]=ind_l;      
        }  

        dchi=chi_l-chi_u;
      
        mat_from_weight[i][j][k][0]=(chi_l-chi_red_l)/dchi*w_plane_left;
        mat_from_weight[i][j][k][1]=(chi_red_l-chi_u)/dchi*w_plane_left;

        ////right plane
        for(LO m=0;m<p3->versurf[i+p3->li1];m++){
          tmp[m]=std::abs(p3->zcoordsurf[i][m]-chi_red_r);
        }

        ind_u=minloc(tmp,p3->versurf[i+p3->li1]);

        mat_from_ind_n[i][j][k][2]=ind_u;
        chi_u=p3->zcoordsurf[i][ind_u];       

        if((chi_red_r-chi_u)>0){
          ind_l=ind_u+1;
          if(ind_l>p3->versurf[p3->li1+i]-1){
            ind_l=0;
            chi_l=p3->zcoordsurf[i][ind_l]+2.0*cplPI;
          }else{
            chi_l=p3->zcoordsurf[i][ind_l];
          }
          mat_from_ind_n[i][j][k][3]=ind_l;
        }else{
          ind_l=ind_u-1;
          if(ind_l<0){
            ind_l=p3->versurf[p3->li1+i]-1;
            chi_l=p3->zcoordsurf[i][ind_l]-2.0*cplPI;
          }else{
            chi_l=p3->zcoordsurf[i][ind_l];
          }
          mat_from_ind_n[i][j][k][3]=ind_l;
        }       

        dchi=chi_l-chi_u;
        mat_from_weight[i][j][k][2]=(chi_l-chi_red_r)/dchi*w_plane_right;
        mat_from_weight[i][j][k][3]=(chi_red_r-chi_u)/dchi*w_plane_right;  
      }
    }
    free(tmp);
  }
  for(LO i=0;i<p3->li0;i++){
    free(p3->zcoordsurf[i]);
  }
  free(p3->zcoordsurf);

  bool debug=false;
  if(debug){
    LO* inds1=new LO[p1->li0];
    for(LO i=0;i<p1->li0;i++) inds1[i]=4;
    LO* inds2=new LO[p1->li0];
    for(LO i=0;i<p1->li0;i++) inds2[i]=2;

    double sum_weight=0.0;

    printSumm4D(mat_from_weight,p1->li0,p1->y_res_back,p1->lk0,inds1, sum_weight,
       MPI_COMM_WORLD,"mat_from_weight",0);

    printminmax4d(mat_from_ind_n,p1->li0,p1->y_res_back,p1->lk0,4,
     MPI_COMM_WORLD, "mat_from_ind_n",0);

    printminmax4d(mat_from_ind_plane,p1->li0,p1->y_res_back,p1->lk0,2,
     MPI_COMM_WORLD, "mat_from_ind_plane",0);

    printminmax4d(mat_from_weight,p1->li0,p1->y_res_back,p1->lk0,4,
     MPI_COMM_WORLD, "mat_from_weight",0);

  }
  
}

void DatasProc3D::TestInitPotentAlongz(const Part3Mesh3D* p3m3d,const LO npy, const LO n) 
{
  if(npy==1){
    LO li0,lj0,lk0;
    li0=p3->li0;
    lj0=p3->lj0;
    double ylen;
    double sum;
    double dy=2.0*cplPI/double(lj0);
    for(LO i=0;i<li0;i++){
      lk0=p3->mylk0[i];
      for(LO k=0;k<lk0;k++){
        ylen=0.0;
        for(LO j=0;j<lj0;j++){
          ylen=double(j)*dy;
          sum=0.0;
          for(LO h=0;h<n;h++){
            sum+=cos(double(h+1)*ylen-cplPI);
          }
          potentin[i][j][k]=sum;
        }
      }
    }
  }
}


}

