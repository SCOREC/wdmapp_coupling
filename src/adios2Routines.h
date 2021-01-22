#ifndef ADIOS2ROUTINES_H_
#define ADIOS2ROUTINES_H_

#include <adios2.h>
#include <iostream>
#include <cassert>
//#include <Kokkos_Core.hpp> //not used
#include <typeinfo> //not used
#include <string>
#include <fstream>
#include "couplingConstants.h" //not used
#include "couplingTypes.h"

namespace coupler { 
      
  struct adios2_handler{
    adios2::IO IO;
    adios2::Engine eng;
    std::string name;
    adios2_handler(adios2::ADIOS &adios, const std::string name_):
            name(name_),IO(adios.DeclareIO(name_))  {}
    void close() {
      assert(eng);
      eng.Close();
    }
    std::string get_name() const { return name; };
  };

  /** Storage of double precision 2D array data
 *       */
  template<class T>
  class Array3d {
    public:
      Array3d(GO dim0, GO dim1, GO dim2, GO* start) :
        DIM0(dim0), DIM1(dim1), DIM2(dim2), START(start) {
          vals = new T[dim0*dim1*dim2];
      }
      ~Array3d() {
        DIM0 = DIM1 = DIM2 = 0;
        for(LO i=0;i<3;i++) START[i]=0;
        delete [] vals;
      }
      T val(long i) const {
        assert(i<(DIM0*DIM1*DIM2));
        return vals[i];
      }
      T* data() const { return vals; };
      GO dim0() const { return DIM0; };
      GO dim1() const { return DIM1; };
      GO dim2() const { return DIM2; };
      GO* start() const { return START; };
    private:
      T* vals;
      GO DIM0;
      GO DIM1;
      GO DIM2;
      GO* START;
  };

  /** Storage of double precision 2D array data
   *  and associated meta data
   */
  template<class T> 
  class Array2d {
    public:
      Array2d(GO gH, GO gW, GO lH, GO lW, GO start) :
        globH(gH), globW(gW), locH(lH), locW(lW), locFirstCol(start) {
          vals = new T[locH*locW];
      }   
      ~Array2d() {
        globH = globW = locH = locW = 0;
        delete [] vals;
      }   
      T val(long i) const {
        assert(i<(locH*locW));
        return vals[i];
      }   
      T* data() const { return vals; };
      GO globalH() const { return globH; };
      GO globalW() const { return globW; };
      GO localH() const { return locH; };
      GO localW() const { return locW; };
      GO start_col() const { return locFirstCol; };
    private:
      T* vals;
      GO globH;
      GO globW;
      GO locH;
      GO locW;
      GO locFirstCol;
  };

  template<class T>
  class Array1d {
    public:
      Array1d(GO gW, GO lW, GO start) :
        globW(gW), locW(lW), locFirstCol(start) {
        vals = new T[locW];
      }
      ~Array1d() {
        globW = locW = 0;
        delete [] vals;
      }
      T val(long i) const {
        assert(i<(locW));
        return vals[i];
      }
      T* data() const { return vals; };
      GO globalW() const { return globW; };
      GO localW() const { return locW; };
      GO start_col() const { return locFirstCol; };
    private:
      T* vals;
      GO globW;
      GO locW;
      GO locFirstCol;
  };
 
  /** Destroy Array1d and Array2d objects
   */
  template<typename T>
    void destroy(Array1d<T>* a) {
      delete a;
    }

  template<typename T>
    void destroy(Array2d<T>* a) {
      delete a;
    }

  template<typename T>
    void destroy(Array3d<T>* a) {
      delete a;
    }


  /** Sanity check values
   */
  template<typename T>
  void printSomeDensityVals(const Array2d<T>* density) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //asserting the density values received from GENE
    if(!rank)
    {
      for (int i = 0; i < 10; i++)
      {
        std::cerr << rank <<  ": first 10 density at "<< i
          << " is "<< density->val(i) <<"\n";
      }
      for (int i = 0; i < 10; i++)
      {
        std::cerr << rank << ": first 10 for rank 1 at: [67236]" << " + "<< i
          << " is " << density->val(67236 + i) << "\n";
      }
    }
  
    if(rank == 1)
    {
      for (int i = 0; i < 10; i++)
      {
        int offset = ((density->localW() - 1) * density->localH()) + 67235 - 9;
        std::cerr << rank << ": last 10 for rank 0 at: [67235 - 9]" << " + "<< i
          << " is " << density->val(offset  + i) << "\n";
      }
      int last_ten = (density->localH() * density->localW()) - 10;
      for (int i = 0; i < 10; i++)
      {
        std::cerr << rank <<  ": last 10 density at " << last_ten + i << " is "
          << density->val(last_ten + i) <<"\n";
      }
    }
  }
  
  
  template<typename T> 
  Array1d<T>* receive1d_from_ftn(const std::string dir, const std::string name,
      adios2::IO &read_io, adios2::Engine &eng,const std::string model) {
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  
    const std::string fname = dir + "/" + name + ".bp";
  
    if(!eng){
      read_io.SetEngine("Sst");
      read_io.SetParameters({
          {"DataTransport","RDMA"},
          {"OpenTimeoutSecs", "480"}
          });
      eng = read_io.Open(fname, adios2::Mode::Read);
      if(!rank) std::cerr << rank << ": " << name << " engine created\n";
    }
    else{
      std::cerr << rank << ": receive engine already exists \n";
    }
  
    eng.BeginStep();
    adios2::Variable<T> adios_var = read_io.InquireVariable<T>(name);

    const auto total_size = adios_var.Shape()[0];
    if(!rank) std::cout <<name<< "  total_size " <<total_size << "\n";
    GO my_start, my_count;
    if(model=="local"){
      my_start =(GO)(total_size / nprocs) * rank;
      my_count =(GO)(total_size / nprocs);
    }else if(model=="global"){
      my_start = 0;
      my_count = GO(total_size);     
    }else{
      std::cout<<"Error: 'model' is not correctly assigned."<<'\n';
      exit(1);
    }

    if(!rank)std::cout << " Reader of rank " << rank << " reading " << my_count
              << " floats starting at element " << my_start << "\n";
  
    const adios2::Dims start{my_start};
    const adios2::Dims count{my_count};
 
    const adios2::Box<adios2::Dims> sel(start, count);
    Array1d<T>* field = new Array1d<T>{total_size, my_count, my_start};
  
    adios_var.SetSelection(sel);
    eng.Get(adios_var, field->data());
    eng.EndStep();
    return field;
  }
  
  
  template<typename T> 
  T* receive1d_exact_ftn(const std::string dir, const std::string name,
      adios2::IO &read_io, adios2::Engine &eng, GO li1, GO li0, MPI_Comm &comm) {
    int rank, nprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);
  
    const std::string fname = dir + "/" + name + ".bp";
  
    if(!eng){
      read_io.SetEngine("Sst");
      read_io.SetParameters({
          {"DataTransport","RDMA"},
          {"OpenTimeoutSecs", "480"}
          });
      eng = read_io.Open(fname, adios2::Mode::Read);
      if(!rank) std::cerr << rank << ": " << name << " engine created\n";
    }
    else{
      std::cerr << rank << ": receive engine already exists \n";
    }
  
    eng.BeginStep();
    adios2::Variable<T> adios_var = read_io.InquireVariable<T>(name);
  
    const auto total_size = adios_var.Shape()[0];
    std::cerr << rank << ": total_size "<<total_size <<" \n";
    const auto my_start = li1;
    const auto my_count = li0;
    std::cout << " Reader of rank " << rank << " of "<<nprocs<<" ranks, reading " << my_count
              << " floats starting at element " << my_start << "\n";
  
    const adios2::Dims start{my_start};
    const adios2::Dims count{my_count};
  
    const adios2::Box<adios2::Dims> sel(start, count);
    T* val = new T[li0];
    adios_var.SetSelection(sel);
    eng.Get(adios_var, val);
    eng.EndStep();
    return val;
  } 
  
  template<typename T>
  Array2d<T>* receive2d_from_ftn(const std::string dir, const std::string name,
      adios2::IO &read_io, adios2::Engine &eng, GO start[2], GO count[2],MPI_Comm &comm,const int m) {
    int rank, nprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);
    std::cout<<"rank="<<rank<<'\n';
  
    const std::string fname = dir + "/" + name + ".bp";
    std::cout<<fname<<'\n'; 
    if(m==0){
      std::cout<<"creat engine for: "<<name<<'\n';
      read_io.SetEngine("Sst");
      read_io.SetParameters({
          {"DataTransport","RDMA"},
          {"OpenTimeoutSecs", "800"}
          });
      std::cout<<"engine parameters are set"<<'\n';
      eng = read_io.Open(fname, adios2::Mode::Read);
      if(!rank) std::cerr << rank << ": " << name << " engine created\n";
    } else{
      std::cerr << rank << ": receive engine already exists \n";
      assert(eng);
    }
 
    eng.BeginStep();
    adios2::Variable<T> adVar = read_io.InquireVariable<T>(name);
 
    const auto ftn_glob_height = adVar.Shape()[0]; 
    const auto ftn_glob_width = adVar.Shape()[1]; 
    std::cout<<"Shape 0 1="<<ftn_glob_width<<" "<<ftn_glob_height<<'\n';
    //fortran to C transpose
    const auto c_glob_height = ftn_glob_width;
    const auto c_glob_width = ftn_glob_height;

    Array2d<T>* a2d = new Array2d<T>(c_glob_height, c_glob_width,
        count[0], count[1], start[0]);

// Here, count,start take care of the fortran to C transpose. 
    const::adios2::Dims my_start({start[0], start[1]});
    const::adios2::Dims my_count({count[0], count[1]});
    const adios2::Box<adios2::Dims> sel(my_start, my_count);
  
    adVar.SetSelection(sel);
    eng.Get<T>(adVar, a2d->data());
    eng.EndStep();
    std::cerr << rank <<  ": receive " << name << " done \n";
    return a2d;
  }
  
  template<typename T>
  Array3d<T>* receive3d_from_ftn(const std::string dir, const std::string name,
      adios2::IO &read_io, adios2::Engine &eng, GO* start,GO* count,MPI_Comm &comm,const int m) {
    int rank, nprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);
 
    const std::string fname = dir + "/" + name + ".bp";
    if(m==0){
      std::cout<<"creat engine for: "<<name<<'\n';
      read_io.SetEngine("Sst");
      read_io.SetParameters({
          {"DataTransport","RDMA"},
          {"OpenTimeoutSecs", "800"}
          });
      std::cout<<"engine parameters are set"<<'\n';
      eng = read_io.Open(fname, adios2::Mode::Read);
      if(!rank) std::cerr << rank << ": " << name << " engine created\n";
    } else{
      std::cerr << rank << ": receive engine already exists \n";
      assert(eng);
    }
 
    eng.BeginStep();
    adios2::Variable<T> adVar = read_io.InquireVariable<T>(name);
 
    const auto dim0 = adVar.Shape()[0]; 
    const auto dim1 = adVar.Shape()[1]; 
    const auto dim2 = adVar.Shape()[2];

    std::cout<<"rank Shape 0 1 2="<<rank<<" "<<dim0<<" "<<dim1<<" "<<dim2<<'\n';
    //fortran to C transpose
    const auto DIM0 = dim2;
    const auto DIM1 = dim1;
    const auto DIM2 = dim0;

    Array3d<T>* a3d = new Array3d<T>(DIM0,DIM1,DIM2,start);

// Here, count,start take care of the fortran to C transpose. 
    const::adios2::Dims my_start({start[0], start[1], start[2]});
    const::adios2::Dims my_count({count[0], count[1], count[2]});
    const adios2::Box<adios2::Dims> sel(my_start, my_count);
  
    adVar.SetSelection(sel);
    eng.Get<T>(adVar, a3d->data());
    eng.EndStep();
    std::cerr << rank <<  ": receive " << name << " done \n";
    return a3d;
  }


  /* send columns (start_col) to (start_col + localW) */
  template<typename T>
  void send2d_from_C(const Array2d<T>* a2d, const std::string dir,
      const std::string name, adios2::IO coupling_io,adios2::Engine engine,
      adios2::Variable<T> &send_id) 
{
    const::adios2::Dims g_dims({a2d->globalW(), a2d->globalH()});
    const::adios2::Dims g_offset({0,a2d->start_col()});
    const::adios2::Dims l_dims({a2d->localW(), a2d->localH()});
    const std::string fname = dir + "/" + name + ".bp";
    if (!engine){
      coupling_io.SetEngine("Sst");
      send_id = coupling_io.DefineVariable<T>(name,
          g_dims, g_offset, l_dims);
      adios2::Engine  engine = coupling_io.Open(fname, adios2::Mode::Write);
    }
    else{
      std::cerr << ": field engine already created \n";
    }

    engine.BeginStep();
    engine.Put<T>(send_id, a2d->data());
    engine.EndStep();
    std::cout<<"The cpl_density was written"<<'\n';
 }
  
  /** Receive PreProc values from GENE
   */
  template<typename T>
  Array1d<T>* receive_pproc(const std::string cce_folder,
      adios2_handler &handler, const std::string model) { 
      std::string name = handler.get_name();
    return receive1d_from_ftn<T>(cce_folder,name, handler.IO, handler.eng, model);
  }

  template<typename T>
  Array2d<T>* receive_pproc_2d(const std::string cce_folder,
      adios2_handler &handler,GO my_start[2], GO my_count[2], const int m, MPI_Comm comm = MPI_COMM_WORLD) { 
      std::string name = handler.get_name();
    return receive2d_from_ftn<T>(cce_folder,name, handler.IO, handler.eng,my_start, my_count, comm, m);
  }

  template<typename T>
  Array3d<T>* receive_pproc_3d(const std::string cce_folder,
      adios2_handler &handler, GO* my_start, GO* my_count, const int m, MPI_Comm comm = MPI_COMM_WORLD) {
      std::string name = handler.get_name();
   return receive3d_from_ftn<T>(cce_folder, name, handler.IO, handler.eng, my_start, my_count, comm, m);
  } 

  template<typename T>
  T* receive_exact(const std::string cce_folder,
      adios2_handler &handler, GO my_start, GO my_count, MPI_Comm comm = MPI_COMM_WORLD) { 
      std::string name = handler.get_name();
    return receive1d_exact_ftn<T>(cce_folder,name, handler.IO, handler.eng, my_start, my_count, comm);
  }
  

template<typename T>
 void send_from_coupler(adios2::ADIOS &adios,const std::string cce_folder,
        const Array2d<T>* a2d, adios2::IO& sendIO,adios2::Engine& engine,const std::string fldname,                        
        adios2::Variable<T> &send_id, const MPI_Comm comm,const int m)  //     
 {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const std::string fld_name;

    const::adios2::Dims g_dims({a2d->globalW(), a2d->globalH()});
    const::adios2::Dims g_offset({0,a2d->start_col()});
    const::adios2::Dims l_dims({a2d->localW(), a2d->localH()});

    const std::string fname = cce_folder + "/" + fldname + ".bp";
  
    if(m==0){
      sendIO.SetEngine("Sst");
      sendIO.SetParameters({
      {"OpenTimeoutSecs", "480"}
          });
      send_id = sendIO.DefineVariable<T>(fldname,
          g_dims, g_offset, l_dims);
      printf("before cpl_density open \n"); 
      MPI_Barrier(MPI_COMM_WORLD);
      engine = sendIO.Open(fname, adios2::Mode::Write,comm);
    }
    if(engine) std::cout<<"sending Engine for "<<fldname<<" is created."<<'\n';   
    engine.BeginStep(adios2::StepMode::Append);
    engine.Put<T>(send_id, a2d->data());
    engine.EndStep();
    std::cout<<"rank="<<rank<<" "<<"The "<<fldname<<" was written"<<'\n';
  }
  

  Array2d<CV>* receive_density(const std::string cce_folder,
      adios2_handler &handler,GO my_start[2],GO my_count[2], MPI_Comm comm, const int m);
 
  void send_density(const std::string cce_folder, const Array2d<double>* density,
      const adios2_handler &handler, adios2::Variable<double> &send_id);

  Array2d<double>* receive_field(const std::string cce_folder,
      adios2_handler &handler, GO my_start[2], GO my_count[2], MPI_Comm comm,const int m);

  void send_field(const std::string cce_folder, const Array2d<CV>* field,
      const adios2_handler &handler, adios2::Variable<CV> &send_id); 

  void AdiosProTransFortrandCpp2D(LO rankout,LO mypxout,LO mypyout, const LO mypex,const LO mypey,
      const LO npx,const LO npy);

  void send_density_coupler(adios2::ADIOS &adios,const std::string cce_folder, const Array2d<double>* density, adios2::Variable<double> &send_id,const MPI_Comm comm);

}//end namespace coupler

#endif
