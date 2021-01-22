#include "adios2Routines.h"
/*
void exParFor() {
  Kokkos::parallel_for(
      4, KOKKOS_LAMBDA(const int i) {
        printf("Hello from kokkos thread i = %i\n", i);
      });
}
*/
int main(int argc, char **argv){
  int rank, nprocs;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
/*
  Kokkos::initialize(argc, argv);
  if(!rank) {
    printf("Hello World on Kokkos execution space %s\n",
         typeid(Kokkos::DefaultExecutionSpace).name());
    exParFor();
  }
*/
  adios2::ADIOS adios(MPI_COMM_WORLD, adios2::DebugON);
  adios2::Variable<double> send_dbl_var;
  adios2::Variable<coupler::CV> send_cv_var;
  const std::string dir = "../coupling";
  const int time_step = 1, RK_count = 4;

  coupler::adios2_handler gDens(adios,"gene_density");
  coupler::adios2_handler cDens(adios,"cpl_density");
  coupler::adios2_handler xFld(adios,"xgc_field");
  coupler::adios2_handler cFld(adios,"cpl_field");

  int m;
  for (int i = 0; i < time_step; i++) {
    for (int j = 0; j < RK_count; j++) {
      m=i*time_step+RK_count;
      coupler::GO start[2]={0, 100}; //FIXME
      coupler::GO count[2]={10, 42}; //FIXME
      MPI_Comm subcomm = MPI_COMM_WORLD;
      coupler::Array2d<coupler::CV>* density = coupler::receive_density(dir, gDens, start, count, subcomm, m);
      coupler::Array2d<double> densitySend(2,2,1,1,0); //FIXME
      coupler::send_density(dir, &densitySend, cDens, send_dbl_var);
      coupler::destroy(density);
      coupler::destroy(&densitySend);

      coupler::Array2d<double>* field = coupler::receive_field(dir, xFld, start, count, subcomm, m);
      coupler::Array2d<coupler::CV> fieldSend(2,2,1,1,0); //FIXME
      coupler::send_field(dir, &fieldSend, cFld, send_cv_var);
      coupler::destroy(field);
      coupler::destroy(&fieldSend);
    }
  }

//  Kokkos::finalize();
  MPI_Finalize();
  return 0;
}

