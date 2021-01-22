#include "dataprocess.h"
#include "importpart3mesh.h"
#include "commpart1.h"
#include "BoundaryDescr3D.h"
#include "sendrecv_impl.h"
#include "interpoutil.h"
#include <cassert>
#include <cmath>

namespace coupler {

  gemXgcDatasProc3D::gemXgcDatasProc3D(const Part1ParalPar3D* p1pp3d,
    const Part3Mesh3D* p3m3d,
    const BoundaryDescr3D* bdesc_,
    bool pproc,
    TestCase test_case, 
    bool ypar)
 :  preproc(pproc),
    testcase(test_case),
    yparal(ypar)  
 {  
   p1 = p1pp3d;
   p3 = p3m3d;
   bdesc = bdesc_;  

   allocDensityArrays(); 

   allocPotentArrays();

   allocSendRecvbuff();

 }


 void gemXgcDatasProc3D::allocDensityArrays()
 {
   densin=new double**[p1->imx+1];
   for(LO i=0;i<p1->imx+1;i++){
     densin[i]=new double*[p1->jmx+1];
     for(LO j=0;j<p1->jmx+1;j++)
       densin[i][j]=new double[2];       
   }

   densCpl=new double**[p1->li0];
   for(LO i=0;i<p1->li0;i++){ 
     densCpl[i]=new double*[p1->lj0];
       for(LO j=0;j<p1->lj0;j++)
         densCpl[i][j]=new double[p1->lk0];     
   }

   densXgc = new double**[p3->li0];
   for(LO i=0;i<p3->li0;i++){
     densXgc[i]=new double*[p3->lj0];
     for(LO j=0;j<p3->lj0;j++)
        densXgc[i][j]=new double[p3->mylk0[i]];
   } 

   densinterone = new double**[p1->li0];
   for (LO i=0; i<p1->li0; i++) {
     densinterone[i] = new double*[p1->lj0];
     for (LO j=0; j<p1->lj0; j++) {
       densinterone[i][j] = new double[p3->mylk0[i]];
     }
   }      

   densintertwo = new double**[p1->li0];
   for (LO i=0; i<p1->li0; i++) {
     densintertwo[i] = new double*[p1->nphi];
     for (LO j=0; j<p1->nphi; j++) {
       densintertwo[i][j] = new double[p3->mylk0[i]];
     }
   }
 
   denssend = new double[p3->blockcount*p3->nphi];
 }

 void gemXgcDatasProc3D::allocPotentArrays()
 {
   pot_gem_fl=new double***[p1->li0];
   potyCpl=new double**[p1->li0];
   for(LO i=0;i<p1->li0;i++){
     pot_gem_fl[i]=new double**[p1->lj0];
     potyCpl[i]=new double*[p1->lj0];
     for(LO j=0;j<p1->lj0;j++){
       pot_gem_fl[i][j]=new double*[p3->mylk0[i]];
       potyCpl[i][j]=new double[p3->mylk0[i]];
	 for(LO k=0;k<p3->mylk0[i];k++)
	   pot_gem_fl[i][j][k]=new double[4];
     }
   }
  
   potythCpl=new double**[p1->li0];
   for(LO i=0;i<p1->li0;i++){
     potythCpl[i]=new double*[p1->lj0];
     for(LO j=0;j<p1->lj0;j++) potythCpl[i][j]=new double[p1->lk0];
   }

   potGem=new double[p1->tli0*p1->lj0*p1->glk0];
 } 

 void gemXgcDatasProc3D::allocSendRecvbuff()
 {
   numsend=new LO[p1->NP];
   numrecv=new LO[p1->NP];
   sdispls=new LO[p1->NP];
   rdispls=new LO[p1->NP];

   // sending side is in tube_comm X grid_comm collective
   sendnum=0;
   recvnum=0;
   for(LO i=0;i<p1->NP;i++){
     numsend[i]=0;
     numrecv[i]=0;
     sdispls[i]=0;
     rdispls[i]=0;
   }
   
   LO rank;
   for (LO i=0; i<p1->sendOverlap_x.size(); i++) {
     for (LO j=0; j<p1->sendOverlap_th.size(); j++) {
       rank = p1->sendOverlap_x[i][0]*p1->npz + p1->sendOverlap_th[j][0];
       numsend[rank] = numsend[rank] + (p1->sendOverlap_x[i][2] - p1->sendOverlap_x[i][1]+1)
         *(p1->sendOverlap_th[j][2] - p1->sendOverlap_th[j][1]+1)*p1->lj0;
     }
   }
   sendnum = sendnum + numsend[0]; 
   for(LO i=1;i<p1->NP;i++){ 
     sendnum = sendnum + numsend[i];
     sdispls[i]=sdispls[i-1]+numsend[i-1];
   }

   for(LO i=0;i<p1->recvOverlap_x.size();i++){
     for(LO j=0;j<p1->recvOverlap_th.size();j++){
        rank=p1->recvOverlap_th[j][0]*p1->tnpx + p1->recvOverlap_x[i][0];
        numrecv[rank] = numrecv[rank] + (p1->recvOverlap_x[i][2]-p1->recvOverlap_x[i][1]+1)
         *(p1->recvOverlap_th[j][2]-p1->recvOverlap_th[j][1]+1)*p1->lj0;
     }
   }
 
   recvnum = recvnum + numrecv[0];
   for(LO i=1;i<p1->NP;i++){
     rdispls[i]=rdispls[i-1]+numrecv[i-1];
     recvnum = recvnum + numrecv[i];
   }

   bool debug = true;
   if (debug) {
     printf("sendnum: %d, recvnum: %d \n", sendnum, recvnum);
   }
   if (debug) {
     for (LO i=0; i<p1->NP; i++){
       printf("i: %d, numsend: %d, numrecv: %d, sdispls: %d, rdispls: %d\n", i, numsend[i], 
	 numrecv[i], sdispls[i],rdispls[i]);
     }
   

     LO sendtot;
     LO recvtot;
     MPI_Reduce(&sendnum, &sendtot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
     MPI_Reduce(&recvnum, &recvtot, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
     if (p1->mype == 0) {
       printf("sendtot: %d, recvtot: %d \n", sendtot, recvtot);
       assert(sendtot == recvtot);
     }
   }   

   debug = false;
   if (debug) {
     LO* sumsend = new LO[p1->NP];
     LO* recvtmp = new LO[p1->NP];
     MPI_Allreduce(numsend, sumsend, p1->NP, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

     for (LO i=0; i<p1->NP; i++)
     MPI_Gather(&numsend[i], 1, MPI_INT, recvtmp, 1, MPI_INT, i, MPI_COMM_WORLD);

     for (LO i=0; i<p1->NP; i++) {
       if (recvtmp[i] != numrecv[i]) {
	 printf("mype: %d, i: %d, recvtmp: %d, numrecv: %d \n", p1->mype, i, recvtmp[i], numrecv[i]);
       }
     }

     LO rnum = 0;
     for (LO i=0; i<p1->NP; i++) rnum = rnum + recvtmp[i];
     if(recvnum == rnum) printf("mype: %d, recvnum equals rnum \n", p1->mype);
   
     recvnum = recvnum + recvtmp[0];
     for(LO i=1;i<p1->NP;i++){
       rdispls[i]=rdispls[i-1]+recvtmp[i-1];
       recvnum = recvnum + recvtmp[i];
     }


  //   if (recvnum = sumsend[p1->mype]) printf("mype: %d, recvnum does equal to sendnum \n", p1->mype); 

     double* sendbuf=new double[sendnum];
     double* recvbuf=new double[recvnum];

     for(LO i=0; i<sendnum; i++) {
       sendbuf[i] = 0.0;
     }

     for (LO i=0; i<recvnum; i++) {
       recvbuf[i] = 0.0;
     }

     MPI_Alltoallv(sendbuf,numsend,sdispls,MPI_DOUBLE,recvbuf,numrecv,rdispls,MPI_DOUBLE,MPI_COMM_WORLD);
     printf("after mpye: %d\n",p1->mype);
   }
 }

 void gemXgcDatasProc3D::DistriDensiRecvfromGem(const Array3d<double>* densityfromGEM)
 {    
   bool debug = false;
//   double* array=densityfromGEM->data();
//   LO n=0;
  
   MPI_Barrier(MPI_COMM_WORLD); 
   densityFromGemToCoupler(densityfromGEM);  
// printf("1111 \n");
   zDensityBoundaryBufAssign(densCpl);
// printf("2222 \n");
   debug = false;
   if (debug) {	
   for(LO i=0;i<p1->li0;i++){ 
     for(LO k=0;k<p1->lk0;k++){
       for(LO j=0;j<p1->lj0;j++){
       if (p1->mype == 0) {
             printf("k:%d, j:%d, denCpl: %f \n", k, j, densCpl[i][j][k]);
	   }
	 }  
       }
     }
   }
   
   interpoDensityAlongZ(densinterone);
//   printf("3333 \n");

   printf("after interpo, mype: %d \n", p1->mype);
   interpoDensityAlongY();
   MPI_Barrier(MPI_COMM_WORLD);
  
   //fixme: distribute and assemble the datas to get the matrix to be sent to xgc by the adios2 routine     
   distriDataAmongCollective(p1, p3, densintertwo, denssend);   
}


void gemXgcDatasProc3D::DistriPotentRecvfromXGC(const Array2d<double>* fieldfromXGC)
{
  double*** tmp=new double**[p1->li0];
  for(LO i=0;i<p1->li0;i++){
    tmp[i]=new double*[p3->nphi];
    for(LO j=0;j<p3->nphi;j++){
      tmp[i][j]=new double[p3->versurf[i+p1->li1]]; //The boundary buffer is assigned. 
    }
  }    

  //FIXME: distribute potentfromXGC->datas() to tmp;
  double* potent=fieldfromXGC->data();  
  LO n=0;
  for(LO i=0;i<p1->li0;i++){
    for(LO j=0;j<p3->nphi;j++){
      for(LO k=0;k<p3->versurf[i+p1->li1];j++){
        n++;
        tmp[i][j][k]=potent[n];
      }
    }
  }  


  // First interpolation along the field line; 
  double* tmppotent;
  double* tmpflx=new double[4];
  double* tmplength;
  for(LO i=0;i<p1->li0;i++){
    for(LO j=0;j<p1->lj0;j++){ 
/*
      double* tmpthflx=new double[2.0*bdesc->nzb+p3->mylk0[i]];
      for(LO h=0;h<bdesc->nzb;h++) tmpthflx[h]=bdesc->lowpotentgemxgc[h];
      for(LO l=0;l<bdesc->) 
*/ 
      for(LO k=0;k<p3->mylk0[i];k++){ 
        for(LO h=0;h<4;h++){
          tmppotent=p3->thetaflx_pot[i][j][k][h];
          for(LO l=0;l<4;l++) tmpflx[l]=tmp[i][j][p3->thetaflx_ind_pot[i][j][k][h][l]];
          pot_gem_fl[i][j][k][h]=Lag3dInterpo1D(tmppotent,tmpflx,p3->thetaflx_pot[i][j][k][h][4]); 
        }
        tmppotent=pot_gem_fl[i][j][k];
        tmplength=p3->nodesdist_fl[i][j][k];
        potyCpl[i][j][k]=Lag3dInterpo1D(tmppotent,tmplength,p3->nodesdist_fl[i][j][k][4]);     
      }
    } 
  }  

  // The 2nd interpolation along theta
  InterpoPotential3DAlongZ(potyCpl,potythCpl);  
  potentFromCouplerToGem(fieldfromXGC);
}


void gemXgcDatasProc3D::densityFromGemToCoupler(const Array3d<double>* densityfromGEM)
{
   double* array=densityfromGEM->data();
// sending happens in MPI_COMM_WORLD, but the process in grid_comm and tube_comm is mappped to comm_x and comm_y    
   double* sendbuf=new double[sendnum];
   double* recvbuf=new double[recvnum];
   LO xnum=0;
   LO thnum=0;
   LO nrank=0;
   LO n, m;
   bool debug;

   for (LO i=0; i<sendnum; i++) {
     sendbuf[i] = 0.0;
   }

   for (LO i=0; i<recvnum; i++) {
     recvbuf[i] = 0.0;
   }

   //for initialize the sendbuf  
   for(LO i=0;i<p1->npx;i++){
     for(LO k=0;k<p1->npz;k++){
       nrank=i*p1->npz+k;
       xnum = 0;
       thnum = 0;
       if(numsend[nrank]!=0){
         while(i!=p1->sendOverlap_x[xnum][0]){
	   xnum+=1;
	   assert(xnum<=p1->sendOverlap_x.size());
	 }
	 while(k!=p1->sendOverlap_th[thnum][0]){
	   thnum+=1;
	   assert(thnum<=p1->sendOverlap_th.size());
	 }      
		
	 for (LO i1=p1->sendOverlap_x[xnum][1]; i1<p1->sendOverlap_x[xnum][2]+1; i1++) {
	   for (LO j1=0; j1<p1->lj0; j1++) {
	     for (LO k1=p1->sendOverlap_th[thnum][1]; k1<p1->sendOverlap_th[thnum][2]+1; k1++) {
//               nrank = p1->sendOverlap_x[xnum][0]*p1->npz + p1->sendOverlap_th[thnum][0]; 
     	       n = (i1-p1->tli1)*p1->lj0*p1->glk0+j1*p1->glk0+k1-p1->glk1;
               m = (i1-p1->sendOverlap_x[xnum][1])*p1->lj0
		 *p1->sendOverlap_th[thnum][3]+j1*p1->sendOverlap_th[thnum][3]
		 +k1-p1->sendOverlap_th[thnum][1];
	       sendbuf[sdispls[nrank]+m]=array[n];  //=densin[i1-p1->tli0][j1][k1-p1->glk0];
	       assert(m<numsend[nrank]);
	       assert(n<p1->tli0*p1->lj0*p1->glk0);
             }
	   }
	 }
    	 
       }
     }    
   }
  
   MPI_Barrier(MPI_COMM_WORLD);   
//   printf("before alltoallv mype: %d\n",p1->mype );

   debug = false;
   if (debug) {
   if (p1->mype == 1) {
     for (LO i=0; i<sendnum; i++)
     printf("i: %d, sendbuf: %f \n", i, sendbuf[i]);
   } 
   }
   MPI_Alltoallv(sendbuf,numsend,sdispls,MPI_DOUBLE,recvbuf,numrecv,rdispls,MPI_DOUBLE,MPI_COMM_WORLD);

   debug = false;
   if (debug) {
   if (p1->mype < 6) {
     for (LO i=0; i<recvnum; i++){
       if(isnan(recvbuf[i])) {
         printf("i: %d, recvbuf[i] is nan. \n", i);
	 exit(1);
       }
//       printf("i: %d, recvbuf: %f \n", i, recvbuf[i]);
     }
   } 
   }

   LO tubenum=0;
   LO gridnum=0;
   nrank=0;
   
   for(LO k=0; k<p1->gnpz; k++){
     for(LO i=0; i<p1->tnpx; i++){
       nrank = k*p1->tnpx + i;
       if(numrecv[nrank]!=0){
	 tubenum=0;
	 while(i!=p1->recvOverlap_x[tubenum][0]){
	   tubenum+=1;
	   assert(tubenum<=p1->recvOverlap_x.size());
	 }
	 gridnum=0;
	 while(k!=p1->recvOverlap_th[gridnum][0]){
	   gridnum+=1;
	   assert(gridnum<=p1->recvOverlap_th.size());
	 }              

	 for(LO i1=p1->recvOverlap_x[tubenum][1]; i1<p1->recvOverlap_x[tubenum][2] + 1; i1++){
	   for(LO j1=0; j1<p1->lj0; j1++){
	     for(LO k1=p1->recvOverlap_th[gridnum][1]; k1<p1->recvOverlap_th[gridnum][2] + 1; k1++){ 
               n = rdispls[nrank]+(i1-p1->recvOverlap_x[tubenum][1])
		 *p1->lj0*p1->recvOverlap_th[gridnum][3]+j1*p1->recvOverlap_th[gridnum][3]+k1-p1->recvOverlap_th[gridnum][1];
	       densCpl[i1-p1->li1][j1][k1-p1->lk1] = recvbuf[n];
	     }
	   }
	 }          
       }
     } 
   }  
}

// It's a reverse procedure of densityFromGemToCoupler
void gemXgcDatasProc3D::potentFromCouplerToGem(const Array2d<double>* fieldfromXGC)
{
   double* array = fieldfromXGC->data();
   GO scounts=recvnum;
   GO rcounts=sendnum;
   double* sendbuf=new double[scounts];
   double* recvbuf=new double[rcounts];

   //for initialize the sendbuf  
   LO tubenum=0;
   LO gridnum=0;
   LO nrank=0;
   LO n, m;

   for(LO k=0; k<p1->gnpz; k++){
     for(LO i=0; i<p1->tnpx; i++){
       nrank = k*p1->tnpx + i;
       if(numrecv[nrank] != 0){
	 tubenum=0;
	 while(i!=p1->recvOverlap_x[tubenum][0]){
	   tubenum+=1;
	   assert(tubenum<=p1->recvOverlap_x.size());
	 }
	 gridnum=0;
	 while(k!=p1->recvOverlap_th[gridnum][0]){
	   gridnum+=1;
	   assert(gridnum<=p1->recvOverlap_th.size());
	 }              

	 for(LO i1 = p1->recvOverlap_x[tubenum][1]; i1 < p1->recvOverlap_x[tubenum][2] + 1; i1++){
	   for(LO j1 = 0; j1<p1->lj0; j1++){
	     for(LO k1 = p1->recvOverlap_th[gridnum][1]; k1 < p1->recvOverlap_th[gridnum][2] + 1; k1++){
               n = rdispls[nrank]+(i1-p1->recvOverlap_x[tubenum][1])
		 *p1->lj0*p1->recvOverlap_th[gridnum][3]+j1*p1->recvOverlap_th[gridnum][3]+k1-p1->recvOverlap_th[gridnum][1];
	       sendbuf[n] = potythCpl[i1-p1->li1][j1][k1-p1->lk1];
	     }
	   }
	 }          
       }
     }
   } 

   MPI_Alltoallv(sendbuf,numrecv,rdispls,MPI_DOUBLE,recvbuf,numsend,sdispls,MPI_DOUBLE,MPI_COMM_WORLD);

   LO xnum=0;
   LO thnum=0;
   nrank=0;
   
   for (LO i=0; i<scounts; i++) {
     sendbuf[i] = 0.0;
   }

   for (LO i=0; i<rcounts; i++) {
     recvbuf[i] = 0.0;
   }

   //for initialize the sendbuf  
   for(LO i=0;i<p1->npx;i++){
     for(LO k=0;k<p1->npz;k++){
       nrank=i*p1->npz+k;
       xnum = 0;
       thnum = 0;
       if(numsend[nrank]!=0){
	 while(i!=p1->sendOverlap_x[xnum][0]){
	   xnum+=1;
	   assert(xnum<=p1->sendOverlap_x.size());
	 }
	 while(k!=p1->sendOverlap_th[thnum][0]){
	   thnum+=1;
	   assert(thnum<=p1->sendOverlap_th.size());
	 }      
		
	 for(LO i1=p1->sendOverlap_x[xnum][1]; i1<p1->sendOverlap_x[xnum][2] + 1; i1++){
	   for(LO j1=0; j1<p1->lj0; j1++){
	     for(LO k1=p1->sendOverlap_th[thnum][1]; k1<p1->sendOverlap_th[thnum][2]+1; k1++){
               n = (i1-p1->tli1)*p1->lj0*p1->glk0+j1*p1->glk0+k1-p1->glk1;
               m = (i1-p1->sendOverlap_x[xnum][1])*p1->lj0
                 *p1->sendOverlap_th[thnum][3]+j1*p1->sendOverlap_th[thnum][3]
                 +k1-p1->sendOverlap_th[thnum][1];
               potGem[n] = recvbuf[sdispls[nrank]+m];
               assert(m<numsend[nrank]);
               assert(n<p1->tli0*p1->lj0*p1->glk0);

/*
		potGem[(k1-p1->sendOverlap_th[thnum][1])*p1->lj0*(p1->imx+1)+j1*(p1->imx+1)+i1]
		=recvbuf[sdispls[n]+(i1-p1->sendOverlap_x[xnum][1])*p1->lj0
		 *p1->sendOverlap_th[thnum][3]+j1*p1->sendOverlap_th[thnum][3]
		 +k1-p1->sendOverlap_th[thnum][1]];
*/
	     }
	   }
	 }     
       }
     }
   } 
   MPI_Reduce(MPI_IN_PLACE,potGem,p1->tli0*p1->lj0*p1->glk0,MPI_DOUBLE,MPI_SUM,0,p1->tube_comm);

}

void gemXgcDatasProc3D::distriDataAmongCollective(const Part1ParalPar3D* p1, const Part3Mesh3D* p3, double*** inmatrix, double* outmatrix)
{
  bool debug;
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

  for(LO j=0; j<p3->nphi; j++){
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

      mpitype = getMpiType(double());
      MPI_Barrier(p1->comm_z);

      MPI_Allgatherv(inmatrix[i][j], p3->mylk0[i], mpitype, tmp[i], recvcount, rdispls,
                     mpitype, p1->comm_z);
      MPI_Barrier(p1->comm_z);

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
        outmatrix[j*p3->blockcount+h] = blocktmp[h];
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

}



} /*gemXgcDataproc.cc*/ 
