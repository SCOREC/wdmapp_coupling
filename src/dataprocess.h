#ifndef DATAPROCESS_H
#include "couplingTypes.h"
#include "testutilities.h"
#include <fftw3.h>
#include "importpart3mesh.h"
#include "commpart1.h"
#include "BoundaryDescr3D.h"

namespace coupler {

// forward declare
template<class T>
class Array2d;

template<class T>
class Array3d;

class BoundaryDescr3D;
class Part3Mesh3D;
class Part1ParalPar3D;

class DatasProc3D {
public:
  LO part1li0; // part1li0 is the element number of subdomain of x on each
               // process belonging to comm_y. This x subdomain belongs to 2d
               // box on the x-y plane storing the input charged ion density for
               // fourier transform.
  LO part3li0; // part3li0 is the element number of subdomain of x on each
               // process belonging to comm_y. This x subdomain belongs to 2d 
               // box on the x-y plane storing the input electrostatic 
               // potential.
  LO part1lj0; // the count of elements on y domain on each process after
               // backward Fourier transform
  LO part3lj0; ////the count of elements on y domain on each process after
               ///forward Fourier transform  
  LO sum;
  // here, pointers must be assigned a NULL;
  CV** densrecv = NULL; // the 2d array density recived by the coupler.
                           // and sent by part1. This is received in comm_x and comm_y communicators.
  CV*** densin = NULL; // input 3d density in complex number
  CV*** densinterpo = NULL; 
  CV* densintmp = NULL; // temporary 2d density array prepared for backward
                        // fourier transform
  double* densouttmp = NULL; // store the x-y 2d real density after backward 
                             // fourier transform
  double*** denspart3 = NULL; // storing the density being sent to the part3
  double*** densTOpart3 = NULL;

  double*  denssend = NULL; // the 1d array density  sent to part3
                             //  would be sent in comm_x and comm_y communicators


  double** potentrecv = NULL; // the 2d array potential received by the coupler 
                                 // and sent by part3. This is received in comm_x and comm_y communicators.
  double*** potentin = NULL; // the input real electrostatic potential in 3d xyz
  double* potentintmp = NULL;
  CV* potentouttmp = NULL; // 
  CV*** potentinterpo = NULL; // temporary xy 2d potential array for forward 
                              // fourier transform
  CV*** potentpart1 = NULL; // storing the electrostatic potential being sent
                            // to the part1.
  CV*  potentsend = NULL; // the 1d array complex potential sent  to part1, and would be sent 
                           // in comm_x and comm_y communicators. 

  // matrix for the transformation between planes and xyz
  double**** mattoplane=NULL;
 
  CV****     mat_to_plane=NULL; 
  double**** mat_from_weight=NULL;
  int****    mat_from_ind_plane=NULL;
  int****    mat_from_ind_n=NULL; 


  fftw_plan plan_forward = NULL, plan_backward = NULL;
  LO myli0;
  // For interpolation mesh
  double* mesh1ddens; // 1d mesh for interpolating density;
  double** mesh1dpotent; // the mesh for interpoating potential;

  /* constructor
   * optional argument supports setting
   * the prepoc and yparal modes
   */
  DatasProc3D(const Part1ParalPar3D* p1pp3d,
      const Part3Mesh3D* p3m3d,
      const BoundaryDescr3D* bdesc_,
      bool pproc = true,
      TestCase test_case = TestCase::off,
      bool ypar = false,
      int nummode = 1);
  ~DatasProc3D(){};
  void InterpoDensity3D();
  void InterpoPotential3D();
  //routines for Fourier transform
  void CmplxdataToRealdata3D();
  void RealdataToCmplxdata3D();
  void InitFourierPlan3D();
  void AllocMatXYZtoPlane();
  void Prepare_mats_from_planes();
  void Initmattoplane();

  void oldAssemDensiSendtoPart3(BoundaryDescr3D& bdesc);
  void AssemDensiSendtoPart3(BoundaryDescr3D& bdesc); 
  void DistriDensiRecvfromPart1(const Array2d<CV>* densityfromGENE);
  void AssemPotentSendtoPart1();
  void DensityToPart3();
  void oldDistriPotentRecvfromPart3(const Array2d<double>* fieldfromXGC);
  void DistriPotentRecvfromPart3(const Array2d<double>* fieldfromXGC); 
  void oldInitmattoplane();  

//boundary buffer
  void zPotentBoundaryBufAssign(const BoundaryDescr3D& bdesc);
  void zDensityBoundaryBufAssign(CV*** box,BoundaryDescr3D& bdesc);

//interpolation
  void mesh1dforDensityInterpo();
  void mesh1dforPotentialInterpo();

  LO getP1li0() { return p1->li0; };
  LO getP1ny0() { return p1->ny0; };
  LO getP1npy() { return p1->npy; };


private:
  const bool preproc;
  const TestCase testcase;
  const bool yparal;
  const Part1ParalPar3D* p1;
  const Part3Mesh3D* p3;
  const BoundaryDescr3D* bdesc;
  /* helper function for destructor */
  void FreeFourierPlan3D(); // called from the destructor - does that make sense?
  /* helper functions for constructor */
  void init();
  void AllocDensityArrays();
  void AllocPotentArrays();
  void AllocBufferForIntepo();
  void TestInitPotentAlongz(const Part3Mesh3D* p3m3d,const LO npy, const LO n);
  /* helper functions for CmplxdataToRealdata3D and RealdataToCmplxdata3D */
  void ExecuteRealToCmplx();
  void ExecuteCmplxToReal();

  };

class gemXgcDatasProc3D {
  public:
    double*** densin = NULL;  // Store the input density datas from GEM; 
    double*** densCpl = NULL;  // Store the density on the coupling subcommunicator of GEM
    double*** densinterone = NULL;  // Store the density interpolated along theta
    double*** densintertwo = NULL;  // Store the density interpolated along y
    double*** densXgc = NULL;
    double*  denssend = NULL;

    double**** pot_gem_fl = NULL;
    double*** potyCpl = NULL;
    double*** potythCpl = NULL;
    double* potGem = NULL;

    LO* numsend = NULL;
    LO* numrecv = NULL;
    LO* sdispls = NULL;
    LO* rdispls = NULL;
    LO sendnum;
    LO recvnum;

    gemXgcDatasProc3D(const Part1ParalPar3D* p1pp3d,
      const Part3Mesh3D* p3m3d,
      const BoundaryDescr3D* bdesc_,
      const bool pproc = true,
      const TestCase test_case = TestCase::off,
      const bool ypar = false);

    ~gemXgcDatasProc3D(){};
    
    void DistriDensiRecvfromGem(const Array3d<double>* densityfromGEM);
    void DistriPotentRecvfromXGC(const Array2d<double>* fieldfromXGC);
    void densityFromGemToCoupler(const Array3d<double>* densityfromGEM);  
 
  private:
    const bool preproc;
    const TestCase testcase;
    const bool yparal;
    const Part1ParalPar3D* p1;
    const Part3Mesh3D* p3;
    const BoundaryDescr3D* bdesc;

    void allocDensityArrays();
    void allocPotentArrays();
    void allocSendRecvbuff();
    void interpoDensityAlongZ(double*** box);
    void interpoDensityAlongY();
    void InterpoPotential3DAlongZ(double*** boxyin, double*** boxout); 
    void potentFromCouplerToGem(const Array2d<double>* fieldfromXGC);
    void zPotentBoundaryBufAssign(const double*** box,BoundaryDescr3D& bdesc);
    void zMeshPotentBoundaryBufAssign(BoundaryDescr3D& bdesc);
    void zDensityBoundaryBufAssign(double*** box);
    void distriDataAmongCollective(const Part1ParalPar3D* p1, const Part3Mesh3D* p3, double*** inmatrix, double* outmatrix);
};


 
void TransposeComplex(CV** InMatrix,CV** OutMatrix, DatasProc3D& dp3d,
     Part1ParalPar3D& p1pp3d);


} // namespace coupler

#endif
