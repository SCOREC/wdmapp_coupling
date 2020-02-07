#include<assert.h>
#include<adios2.h>
#include<iostream>
#include <mpi.h>
#include<vector>
#include<string.h>


void write_density_to_XGC(const std::vector<float> &dens, int rank, int size);
void read_field_from_XGC(std::vector<float> &field, int rank, int size);

void write_field_to_GENE(const std::vector<float> &field,int rank, int size);
void read_density_from_GENE(std::vector<float> &dens, int rank, int size);


int main(int argc, char **argv){
  int rank, size, step = 0;
  int field_step = 0, density_step = 0;
  std::vector<float> density, field;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  while(step < 5)
  {
  //get density data from GENE
  read_density_from_GENE(density, rank, size);

  //CPL-XGC Writer - send density data to XGC in the same format
  write_density_to_XGC(density, rank, size);
  density_step++;

  //get field data from XGC
  read_field_from_XGC(field, rank, size);
  
  //send field data to GENE
  write_field_to_GENE(field, rank, size);
  field_step++;
  std::cout << "This is for time step " << step <<std::endl;
  step++;
  }
  //  end loop
  //
  MPI_Barrier(MPI_COMM_WORLD);
  if( !rank )
  {
    assert (density[0] == (float)0);
    std::cout << "The asserted first density value is " << density[0] << "\n";

    assert (field[0] == (float)0);
    std::cout << "The asserted first field value is " << field[0] << "\n";

    std::cout << "cpl done\n";
  }

  MPI_Finalize();
  return 0;
}


void write_density_to_XGC(const std::vector<float> &dens, int rank, int size)
{
  const std::size_t Nx = dens.size();

  adios2::ADIOS cx_adios(MPI_COMM_WORLD, adios2::DebugON);
  adios2::IO cdensIO = cx_adios.DeclareIO("cxIO");
  cdensIO.SetEngine("Sst");

  auto bp_cdens = cdensIO.DefineVariable<float>("bp_cdens",{size * Nx}, {rank * Nx}, {Nx});
  adios2::Engine cdensWriter = cdensIO.Open("cdens.bp", adios2::Mode::Write);
  cdensWriter.BeginStep();
  cdensWriter.Put<float>(bp_cdens, dens.data());
  cdensWriter.EndStep();
  cdensWriter.Close();
}


void read_field_from_XGC(std::vector<float> &field, int rank, int size)
{
  std::cout << rank << " start reading from xgc\n";
  //  variable declaration
  adios2::Engine xfieldReader;
  adios2::IO xfieldIO;

  adios2::ADIOS xc_adios(MPI_COMM_WORLD, adios2::DebugON);
  xfieldIO = xc_adios.DeclareIO("xcIO");
  xfieldIO.SetEngine("Sst");
  xfieldReader = xfieldIO.Open("xfield.bp", adios2::Mode::Read);
  std::cout << rank << " XGC-to-coupling field engine created\n";

  xfieldReader.BeginStep();
  adios2::Variable<float> bp_xfield = xfieldIO.InquireVariable<float>("bp_xfield");
  std::cout << rank << " Incoming variable of size " << bp_xfield.Shape()[0] << "\n";

  const std::size_t t_size = bp_xfield.Shape()[0];
  const std::size_t my_start = (t_size / size) * rank;
  const std::size_t my_count = (t_size / size);

  const adios2::Dims start{my_start};
  const adios2::Dims count{my_count};
  const adios2::Box<adios2::Dims> sel(start, count);
  std::cout << " xField Reader of rank " << rank << " reading " << my_count
    << " floats starting at element " << my_start << "\n";

  field.resize(my_count);
  bp_xfield.SetSelection(sel);
  xfieldReader.Get(bp_xfield, field.data());
  xfieldReader.EndStep();
  xfieldReader.Close();
  std::cout << rank << " done reading from xgc\n";
}


void write_field_to_GENE(const std::vector<float> &field, int rank, int size)
{
  std::cout << rank << " start writing to gene\n";
 // local array size Nx = 10,   
  const std::size_t Nx = field.size();
  std::cout << rank << " Nx: " << Nx << "\n";

  adios2::ADIOS cg_adios(MPI_COMM_WORLD, adios2::DebugON);
  adios2::IO cfieldIO = cg_adios.DeclareIO("cgIO");
  cfieldIO.SetEngine("Sst");

  auto bp_cfield = cfieldIO.DefineVariable<float>("bp_cfield",{size * Nx}, {rank * Nx}, {Nx});
  adios2::Engine cfieldWriter = cfieldIO.Open("cfield.bp", adios2::Mode::Write);
/*  std::cout << " xField Writerer of rank " << rank << " writing " << my_count
    << " floats starting at element " << my_start << "\n";
*/

  cfieldWriter.BeginStep();
  cfieldWriter.Put<float>(bp_cfield, field.data());
  cfieldWriter.EndStep();
  cfieldWriter.Close();
  std::cout << rank << " done writing to gene\n";
}


void read_density_from_GENE(std::vector<float> &dens, int rank, int size)
{
  adios2::IO gdensIO;
  adios2::Engine gdensReader;

  adios2::ADIOS gc_adios(MPI_COMM_WORLD, adios2::DebugON);
  gdensIO = gc_adios.DeclareIO("gcIO");
  gdensIO.SetEngine("Sst");

  gdensReader = gdensIO.Open("gdens.bp", adios2::Mode::Read);
  std::cout << "GENE-to-coupling density engine created\n";

  gdensReader.BeginStep();
  adios2::Variable<float> bp_gdens = gdensIO.InquireVariable<float>("bp_gdens");
  std::cout << "Incoming variable of size " << bp_gdens.Shape()[0] << "\n";

  const std::size_t t_size = bp_gdens.Shape()[0];
  const std::size_t my_start = (t_size / size) *rank;
  const std::size_t my_count = (t_size / size);
  std::cout << "gDensity Reader of rank "<< rank <<" reading "<< my_count 
  << " floats starting from element "<< my_start << "\n";
  const adios2::Dims start{my_start};
  const adios2::Dims count{my_count};
  const adios2::Box<adios2::Dims> sel(start, count);
  dens.resize(my_count);

  bp_gdens.SetSelection(sel);
  gdensReader.Get(bp_gdens, dens.data());
  gdensReader.EndStep();
  gdensReader.Close();
}