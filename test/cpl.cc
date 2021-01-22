#include "adios2Routines.h"
#include "couplingTypes.h"
#include "dataprocess.h"
#include "importpart3mesh.h"
#include "commpart1.h"
#include "BoundaryDescr3D.h"
#include "testutilities.h"
#include <string>
#include <fstream> 
/*
void exParFor() {
  Kokkos::parallel_for(
      4, KOKKOS_LAMBDA(const int i) {
        printf("Hello from kokkos thread i = %i\n", i);
      });
}
*/

int main(int argc, char **argv){
  int rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
/*
  Kokkos::initialize(argc, argv);
  if(!rank) {
    printf("Hello World on Kokkos execution space %s\n",
         typeid(Kokkos::DefaultExecutionSpace).name());
    exParFor();
  }
*/
  if(argc != 2) {
    if(!rank) printf("Usage: %s <number of timesteps>\n", argv[1]);
    exit(EXIT_FAILURE);
  }
  const std::string dir = "../coupling";
  const int time_step = atoi(argv[1]), RK_count = 4;

  fprintf(stderr, "%d number of time steps: %d\n", rank, time_step);

  const std::string xmlfile = "adios2cfg.xml";
  adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugON);
  adios2::Variable<double> senddensity;
  adios2::Variable<coupler::CV> sendfield;

  coupler::adios2_handler gDens(adios,"gene_density");
  coupler::adios2_handler cDens(adios,"cpl_density");
  coupler::adios2_handler xFld(adios,"xgc_field");
  coupler::adios2_handler cFld(adios,"cpl_field");
  coupler::adios2_handler gQP(adios,"gene_pproc_qp");
  coupler::adios2_handler gRX(adios,"gene_pproc_rx");
  coupler::adios2_handler gCy(adios,"gene_cy_array");
  coupler::adios2_handler gInt(adios,"gene_pproc_i");
  coupler::adios2_handler xXcoord(adios,"xgc_x_coordss");
  coupler::adios2_handler xSurf(adios,"xgc_numsurfs");
  coupler::adios2_handler xZcoord(adios,"xgc_z_coordss");
  coupler::adios2_handler xVsurf(adios,"xgc_versurfs");
  coupler::adios2_handler xCce(adios,"xgc_cce_data");

  //receive GENE's preproc mesh discretization values
  std::string model = "local";
  coupler::Array1d<double>* q_prof = coupler::receive_pproc<double>(dir, gQP,model);
  coupler::Array1d<double>* gene_xval = coupler::receive_pproc<double>(dir, gRX,model);//matching gene's xval arr
  coupler::Array1d<int>* gene_parpar = coupler::receive_pproc<int>(dir, gInt,model);

  //intialize GENE class
  const bool preproc = true;
  const bool ypar = false;
  coupler::TestCase test_case = coupler::TestCase::off;
  coupler::CouplingCase ccase = coupler::CouplingCase::genexgc;
  coupler::Array1d<double>* gene_cy = coupler::receive_pproc<double>(dir, gCy,model);
  coupler::Part1ParalPar3D p1pp3d(gene_parpar->data(),gene_xval->data(),q_prof->data(), gene_cy->data(), preproc);

  //receive XGC's preproc mesh discretization values
  coupler::Array1d<double>* xgc_xcoords = coupler::receive_pproc<double>(dir, xXcoord,model);
  coupler::LO start1d=p1pp3d.mype_x*2;
  coupler::LO count1d=2;
  int* xgc_numsurf=coupler::receive_exact<int>(dir,xSurf,start1d,count1d);
  int numsurf = xgc_numsurf[0];
  int block_count = xgc_numsurf[1];
  double* xgc_zcoords = coupler::receive_exact<double>(dir,xZcoord, 0, block_count);
  int* xgc_versurf = coupler::receive_exact<int>(dir,xVsurf, 0, p1pp3d.nx0);
  coupler::Array1d<int>* xgc_cce = coupler::receive_pproc<int>(dir, xCce,model);
  coupler::Part3Mesh3D p3m3d(p1pp3d, numsurf, block_count, xgc_versurf, xgc_cce->data(), xgc_xcoords->data(), xgc_zcoords, preproc);
  const int nummode = 1;
  coupler::BoundaryDescr3D bdesc(p3m3d, p1pp3d, ccase, test_case, preproc);
  if(!p1pp3d.mype)std::cerr << "0.8"<< "\n"; 

  coupler::Part1ParalPar3D* mesh1=&p1pp3d;
  coupler::Part3Mesh3D*     mesh3=&p3m3d;
  coupler::BoundaryDescr3D* bound=&bdesc;
  coupler::DatasProc3D dp3d(mesh1, mesh3,bound, preproc, test_case, ypar, nummode);
  if(!p1pp3d.mype) std::cerr << "0.9"<< "\n";

  coupler::destroy(q_prof);
  coupler::destroy(gene_xval);
  coupler::destroy(gene_parpar);

  dp3d.InitFourierPlan3D();

  int m;
  double realsum;
  coupler::CV cplxsum;
  bool debug = false;
  coupler::LO* inds3d=new coupler::LO[p1pp3d.li0];
  for(coupler::LO h=0;h<p1pp3d.li0;h++) inds3d[h]=p1pp3d.lk0;
  for (int i = 0; i < time_step; i++) {
    for (int j = 0; j < RK_count; j++) {
      coupler::GO start[2]={0, p1pp3d.blockstart};
      coupler::GO count[2]={coupler::GO(p1pp3d.lj0), p1pp3d.blockcount};
      std::cout<<"mype, start count"<<p1pp3d.mype<<" "<<start[0]<<" "<<start[1]<<" "<<count[0]<<" "<<count[1]<<'\n';
      m=i*RK_count+j;
      coupler::Array2d<coupler::CV>* densityfromGENE = coupler::receive_density(dir, gDens,start,count,MPI_COMM_WORLD,m);
      MPI_Barrier(MPI_COMM_WORLD);

      dp3d.DistriDensiRecvfromPart1(densityfromGENE);
      cplxsum=coupler::CV(0.0,0.0);
      MPI_Barrier(MPI_COMM_WORLD);
      coupler::printSumm3D(dp3d.densin,p1pp3d.li0,p1pp3d.lj0,inds3d,cplxsum,
      MPI_COMM_WORLD,"densityfromGENE",m);


      dp3d.AssemDensiSendtoPart3(bdesc);
      coupler::Array2d<double>* densitytoXGC = new coupler::Array2d<double>(
                                                    p3m3d.activenodes,p3m3d.lj0,p3m3d.blockcount,p3m3d.lj0, 
                                                    p3m3d.blockstart);
      double* densitytmp = densitytoXGC->data();
      for(int h=0;h<p3m3d.lj0*p3m3d.blockcount;h++){
        densitytmp[h] = dp3d.denssend[h]; 
      }
      realsum=0.0;
      coupler::GO INDS1d[2]={0,p3m3d.lj0*p3m3d.blockcount};
      MPI_Barrier(MPI_COMM_WORLD);
      std::cout<<"p3m3d.lj0*p3m3d.blockcount="<<p3m3d.lj0*p3m3d.blockcount<<'\n';
      coupler::printSumm1D(dp3d.denssend,INDS1d,realsum,p1pp3d.comm_x,"densitytoXGC",m);

      if(!debug){
        coupler::send_from_coupler(adios,dir,densitytoXGC,cDens.IO,cDens.eng,cDens.name,senddensity,MPI_COMM_WORLD,m);    
        coupler::destroy(densitytoXGC);
      }

      if(!debug){
        coupler::GO start_1[2]={0,p3m3d.blockstart+p3m3d.cce_first_node-1};
        coupler::GO count_1[2]={coupler::GO(p3m3d.lj0),p3m3d.blockcount}; 
        coupler::Array2d<double>* fieldfromXGC = coupler::receive_field(dir, xFld,start_1,count_1,MPI_COMM_WORLD,m);
        dp3d.DistriPotentRecvfromPart3(fieldfromXGC);
        coupler::destroy(fieldfromXGC);
      }

      if(debug){
        coupler::Array2d<double>* fieldfromXGC = new coupler::Array2d<double>(
               coupler::GO(1), coupler::GO(1),coupler::GO(p3m3d.lj0), p3m3d.blockcount, 
               p3m3d.blockstart+p3m3d.cce_first_node-1);
        double* tmp=fieldfromXGC->data();
        for(coupler::GO h=0;h<p3m3d.lj0*p3m3d.blockcount-1;h++){
           tmp[h]=1.0;
        }
        dp3d.DistriPotentRecvfromPart3(fieldfromXGC);
        coupler::destroy(fieldfromXGC);
      }

      dp3d.AssemPotentSendtoPart1();
      coupler::Array2d<coupler::CV>* fieldtoGENE = new coupler::Array2d<coupler::CV>(
                                                   p1pp3d.totnodes,p1pp3d.lj0,p1pp3d.blockcount,
                                                   p1pp3d.lj0,p1pp3d.blockstart);
      coupler::CV* fieldtmp = fieldtoGENE->data(); 
      for(coupler::GO h=0;h<p1pp3d.lj0*p1pp3d.blockcount;h++){
        fieldtmp[h] = dp3d.potentsend[h];
      }         

  
      coupler::send_from_coupler(adios,dir,fieldtoGENE,cFld.IO,cFld.eng,cFld.name,sendfield,MPI_COMM_WORLD,m);
      coupler::printminmax(dp3d.potentpart1, p1pp3d.li0,p1pp3d.lj0,inds3d,p1pp3d.mype,"fieldtoGENE",m);
      cplxsum=coupler::CV(0.0,0.0);
      coupler::printSumm3D(dp3d.potentpart1, p1pp3d.li0,p1pp3d.lj0,inds3d,cplxsum,
      MPI_COMM_WORLD,"fieldtoGENE",m);

      coupler::destroy(fieldtoGENE);
      std::cerr << p1pp3d.mype << " done loop " << i << " " << j << "\n";
    }
  }

  gDens.close();
  cDens.close();
  xFld.close();
  cFld.close();
  gQP.close();
  gRX.close();
  gInt.close();
  gCy.close();
  xXcoord.close();
  xSurf.close();
  xZcoord.close();
  xVsurf.close();
  xCce.close();

  std::cerr << p1pp3d.mype << " before kokkos finalize\n";
//  Kokkos::finalize();
  std::cerr << p1pp3d.mype << " done kokkos finalize\n";
  MPI_Finalize();
  std::cout<<"MPI is finalized."<<'\n';
  return 0;
} 
/*cpl.cc*/
