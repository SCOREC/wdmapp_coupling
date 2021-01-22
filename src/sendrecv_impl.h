#ifndef SENDRECV_IMPL_H
#define SENDRECV_IMPL_H

#include "couplingTypes.h"
#include "testutilities.h"
#include <mpi.h>
#include <type_traits> // is_same
#include <cassert>

namespace {
//anonymous namespace is not accessible outside this file
//Gerrett - better way to hide this function?

//TODO I think there is a cleaner way - ask Gerrett
template<class T>
MPI_Datatype getMpiType(T foo) {
  MPI_Datatype mpitype;
  //determine the type based on what is being sent
  if( std::is_same<T, double>::value ) {
    mpitype = MPI_DOUBLE;
  } else if ( std::is_same<T, coupler::CV>::value ) { 
    //https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node48.htm
    mpitype = MPI_CXX_DOUBLE_COMPLEX;
  } else if ( std::is_same<T, coupler::GO>::value ) { 
    mpitype = MPI_UNSIGNED_LONG;
  } else if ( std::is_same<T, coupler::LO>::value ) {
    mpitype = MPI_INT;
  } else {
    assert(false);
    fprintf(stderr, "Unknown type in %s... exiting\n", __func__);
    exit(EXIT_FAILURE);
  }
  return mpitype;
}

}

namespace coupler {

template<class T>
void mpisendrecv_aux2D(const MPI_Comm comm,const LO nzb,const LO lx,const LO ly,const LO lz,
     T*** lowbuf,T*** upbuf,T*** box)
{
      MPI_Datatype mpitype = getMpiType(T());
      T* sendbuf=new T[lx*ly*nzb];
      T* recvbuf=new T[lx*ly*nzb];
      MPI_Status status;
      int rank_source,rank_dest;
      MPI_Cart_shift(comm,0,1,&rank_source,&rank_dest);
      for(LO i=0;i<lx;i++){
       for(LO j=0;j<ly;j++){
         for(LO k=0;k<nzb;k++){
           sendbuf[i*ly*nzb+j*nzb+k]=box[i][j][k];
         }
       }
      }
      MPI_Sendrecv(sendbuf,lx*ly*nzb,mpitype,rank_dest,100,recvbuf,lx*ly*nzb,mpitype,rank_source,100,
                  comm,&status);
      for(LO i=0;i<lx;i++){
       for(LO j=0;j<ly;j++){
         for(LO k=0;k<nzb;k++){
              upbuf[i][j][k]=recvbuf[i*ly*nzb+j*nzb+k];
         }
       }
      }
      MPI_Cart_shift(comm,0,-1, &rank_source, &rank_dest);
      for(LO i=0;i<lx;i++){
       for(LO j=0;j<ly;j++){
         for(LO k=0;k<nzb;k++)
            sendbuf[i*ly*nzb+j*nzb+k]=box[i][j][lz-nzb+k];
       }
     }
      MPI_Sendrecv(sendbuf,lx*ly*nzb,mpitype,rank_dest,102,recvbuf,lx*ly*nzb,mpitype,rank_source,102,
                  comm,&status);
      for(LO i=0;i<lx;i++){
       for(LO j=0;j<ly;j++){
         for(LO k=0;k<nzb;k++)
              lowbuf[i][j][k]=recvbuf[i*ly*nzb+j*nzb+k];
       }
      }
  delete[] sendbuf;
  delete[] recvbuf;
}

template<class T>
void mpisendrecv_aux1D(const MPI_Comm comm,const LO nzb,const LO xind,const LO yind,const LO zind,
     T* lowbuf,T* upbuf,const T* box1d)
{

      MPI_Datatype mpitype = getMpiType(T());
      MPI_Status status;
      int rank_source, rank_dest;
      MPI_Cart_shift(comm,0,1,&rank_source,&rank_dest);
      MPI_Sendrecv(box1d,nzb,mpitype,rank_dest,100,upbuf,nzb,mpitype,rank_source,100,
            comm,&status);
      MPI_Cart_shift(comm,0,-1,&rank_source,&rank_dest);
      MPI_Sendrecv(&box1d[zind-nzb],nzb,mpitype,rank_dest,102,lowbuf,nzb,mpitype,rank_source,102,
            comm,&status);
}

/*
template<class T1, class T2, typename T3>
void istriDataAmongCollective(T1* p1, T2* p3, T3*** inmatrix, T3* outmatrix)
{
  bool debug;
  LO* recvcount = new LO[p1->npz];
  LO* rdispls = new LO[p1->npz];
  T3* blocktmp = new T3[p3->blockcount];

  T3** tmp=new T3*[p3->li0];
  for(LO i=0;i<p3->li0;i++){
    tmp[i]=new T3[p3->versurf[p3->li1+i]];
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

      mpitype = getMpiType(T3());
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
*/

}

#endif /*sendrecv_impl.h*/
