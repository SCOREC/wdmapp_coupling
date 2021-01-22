#include "adios2Routines.h"
#include "couplingTypes.h"
#include "dataprocess.h"
#include "importpart3mesh.h"
#include "commpart1.h"
#include "BoundaryDescr3D.h"
#include "testutilities.h"
#include <string>
#include <fstream>

int main(int argc, char **argv){
  int rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::string execut = filename(__FILE__);
  const std::string dir = "../coupling";
  const int time_step = atoi(argv[1]), RK_count = 4;

  fprintf(stderr, "%d number of time steps: %d\n", rank, time_step);

  const std::string xmlfile = "adios2cfg.xml";
  adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugON);
  adios2::Variable<double> senddensity;
  adios2::Variable<coupler::CV> sendfield;

  coupler::adios2_handler gDens(adios,"gem_density");
  coupler::adios2_handler cDens(adios,"density");
  coupler::adios2_handler xFld(adios,"field");
  coupler::adios2_handler cFld(adios,"cpl_field");
  coupler::adios2_handler gMesh(adios,"gem_mesh");
  coupler::adios2_handler gThf(adios,"gem_thfl");
  coupler::adios2_handler gGrd(adios,"gem_grid");
  coupler::adios2_handler gQprof(adios,"gem_qprf");
  coupler::adios2_handler xCouple(adios,"xgc_mesh");
  coupler::adios2_handler xRzcoords(adios,"xgc_rzcoords");  

  std::string model="global";
  int m =0;
  coupler::GO start[2] = {0,0};
  const coupler::Array1d<int>* gmesh=coupler::receive_pproc<int>(dir,gMesh,model);
  coupler::GO ntheta = (coupler::GO)gmesh->val(4);
  coupler::GO nr = (coupler::GO)gmesh->val(1);
  coupler::GO nnode = (coupler::GO)gmesh->val(5);
  printf("%d,%d,%d,%d,%d,%d \n", gmesh->val(0),gmesh->val(1),gmesh->val(2),
            gmesh->val(3), gmesh->val(4), gmesh->val(5));
  printf("mype: %d \n", rank);
  MPI_Barrier(MPI_COMM_WORLD);
/*
  coupler::GO count[2] = {2,(coupler::GO)gmesh->val(5)};
  coupler::GO count2[2] = {0,10};
*/
  const coupler::Array1d<double>* thfl_qprof=coupler::receive_pproc<double>(dir,gQprof,model);
  for(coupler::LO i=0; i<2; i++) if(!rank) fprintf(stderr,"array[%d]: %f\n",i,thfl_qprof->val(i));

  coupler::Array1d<double>* rzcoords = coupler::receive_pproc<double>(dir,gGrd,model);
  if(!rank) fprintf(stderr, "rzcoords[2*nnode-2]: %f, rzcoords[2*nnode-1]: %f \n", 
          rzcoords->val(2*nnode-2),rzcoords->val(2*nnode-1));

  coupler::Array1d<coupler::LO>* xcouple = coupler::receive_pproc<coupler::LO>(dir, xCouple, model);
  if(!rank) fprintf(stderr, "xcouple[5]: %d, xcouple[6]: %d \n", xcouple->val(5), xcouple->val(6));

  //intialize GEM class
  const bool preproc = true;
  const bool ypar = false;
  coupler::TestCase test_case = coupler::TestCase::off;
  coupler::CouplingCase ccase = coupler::CouplingCase::gemxgc; 
  
  coupler::Part1ParalPar3D p1pp3d(gmesh, thfl_qprof, test_case, preproc);
  MPI_Barrier(MPI_COMM_WORLD);

  coupler::Part1ParalPar3D* p1 = &p1pp3d;  
  coupler::Part3Mesh3D p3m3d(p1, xcouple, rzcoords, preproc, test_case);
  MPI_Barrier(MPI_COMM_WORLD);
  
  coupler::BoundaryDescr3D bdesc(p3m3d, p1pp3d, ccase, test_case, preproc);  
 
  coupler::Part3Mesh3D*     p3 = &p3m3d;
  coupler::BoundaryDescr3D* bound = &bdesc;

  coupler::gemXgcDatasProc3D  gxdp3d(p1,p3,bound,preproc,test_case,ypar);
  coupler::gemXgcDatasProc3D  gx3d_ = gxdp3d; 

  coupler::GO* start_adios = new coupler::GO[3];
  coupler::GO* count = new coupler::GO[3];
  coupler::GO dim0, dim1, dim2;
      
  coupler::GO INDS1d[2]={0,p3m3d.lj0*p3m3d.blockcount};

  coupler::Array3d<double>* densityfromGEM;
  coupler::Array2d<double>* densitytoXGC; 
  coupler::Array2d<double>* fieldfromXGC; 
  double* densitytmp;
  
  for (int i = 0; i < 1; i++) {
    for (int j = 0; j < 1; j++) {
      m = i*RK_count+j;
      start_adios[0] = p1->glk1;
      start_adios[1] = 0;
      start_adios[2] = p1->tli1;
      count[0] = p1->glk0;
      count[1] = p1->lj0;
      count[2] = p1->tli0;      
      densityfromGEM = coupler::receive_pproc_3d<double>(dir, gDens, start_adios, count, m, MPI_COMM_WORLD); 

      gxdp3d.DistriDensiRecvfromGem(densityfromGEM);
      MPI_Barrier(MPI_COMM_WORLD);

      densitytoXGC = new coupler::Array2d<double>(
                                  p3m3d.activenodes,p3m3d.lj0,p3m3d.blockcount,p3m3d.lj0,
                                  p3m3d.blockstart);
      densitytmp = densitytoXGC->data();
      for(int h=0; h<p3m3d.nphi*p3m3d.blockcount; h++){
        densitytmp[h] = gxdp3d.denssend[h];
       if(isnan(densitytmp[h])) {
         printf("h: %d, densitytmp[h] is nan. \n", h);
         exit(1);
       }
      }

//      realsum=0.0;

      coupler::send_from_coupler(adios,dir,densitytoXGC,cDens.IO,cDens.eng,cDens.name,senddensity,MPI_COMM_WORLD,m);

      coupler::GO start_1[2]={0, p3m3d.blockstart+p3m3d.cce_first_node-1};
      coupler::GO count_1[2]={coupler::GO(p3m3d.nphi), p3m3d.blockcount};
      fieldfromXGC = coupler::receive_field(dir, xFld,start_1, count_1, MPI_COMM_WORLD,m);
/*
      gxdp3d.DistriPotentRecvfromXGC(fieldfromXGC);
      coupler::destroy(fieldfromXGC);          
*/
    }
  }

  coupler::destroy(densitytoXGC); 
  coupler::destroy(densityfromGEM); 
  gDens.close();
  cDens.close();

  xFld.close();
/*
  cFld.close();
*/
  gMesh.close();
//  gThf.close();
  gGrd.close();
  gQprof.close();
  xCouple.close();
//  xRzcoords.close();

  MPI_Finalize();
  printf("MPI is finalized.\n");
  return 0;
}
/*gemxgc.cc*/
