/*
    AHMED -- Another software library on Hierarchical Matrices for
             Elliptic Differential equations

    Copyright (c) 2012 Mario Bebendorf

    You should have received a copy of the license along with 
    this software; if not, see AHMED's internet site.
*/


#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <assert.h>
#include <string.h>
#include "matgen_omp.h"
//#include "matgen_sqntl.h"
#include "basmod.h"
#include "panel.h"
#include "bemcluster.h"
#include "bemblcluster.h"
#include "solvers.h"
#include "matgen_laplace.h"
#include "H.h"
#include "blas.h"
#include "grid.h"
#include "sllist.h"

//include <fstream.h>
//using namespace std;

//#include "cmn_nrml.h"

template<class T1,class T2>
struct BEMMatrix : public Matrix<double>
{
  mblock<double>** blcks;
  bemblcluster<T1,T2>* blclTree;
  mblock<float>** U;     // the preconditioner in H-Matrix format
  blcluster *blU;    // the corresponding block cluster tree
  float* xf;
  
  BEMMatrix(unsigned N, bemblcluster<T1,T2>* bl) : Matrix<double>(N, N), blclTree(bl)
  { allocmbls(bl, blcks); xf = new float[N]; U=NULL; blU = NULL;}
  virtual ~BEMMatrix() {
    freembls(blclTree, blcks); if(U) freembls(blU, U);
    delete [] xf;
    if(blU) delete blU;
    delete blclTree;
  }

  void amux(double d, double* x, double* y) const {
    mltaHeHVec(d, blclTree, blcks, x, y);
  }

  void precond_apply(double* x) const {
    blas::copy(blU->getn1(), x, xf);
    HCholesky_solve(blU, U, xf);
    blas::copy(blU->getn1(), xf, x);
  }
};



extern void readINPUT(char*, double&, unsigned&, unsigned&, double&, bool&);
extern void mlta_massCL(unsigned, panel*,unsigned*,unsigned*, double, double*, double*);

double xcenter;

/*
double DirData2(vec3d v)
{
  return v[0] + v[1];
}

double NeumData2(vec3d v, vec3d n)
{
  return n[0] + n[1];
}

*/
double DirData2(vec3d v)
{
  v[0] -= 10.0*xcenter;
  return 1.0/sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

double NeumData2(vec3d v, vec3d n)
{
  v[0] -= 10.0*xcenter;
  double r = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  double r3 = -r*r*r;
  vec3d f = v / r3;
  return f*n;
}


//specify input grid(s)
//example: ./PrecondBEM /home/data/BEM/Trias/icosa*.vtk
int main(int argc, const char* argv[]) 
{
  
  if(argc<2) {
    std::cerr << "Please specify input file."<<std::endl;
    exit(1);
  }

  std::cout << "Solution of the inner Dirichlet problem for the Laplacian.";
  std::cout << std::endl;

  // Input-Werte lesen
  unsigned bmin, rankmax;
  double eta, eps, mem_orig, allmem;
  char fname[100];
  //  char fname2[100];
  bool wrtPart;

  readINPUT(fname, eps, rankmax, bmin, eta, wrtPart);
  std::cout << "Parameters read: " << "eps=" << eps
	    << ", bmin=" << bmin << ", eta=" << eta << std::endl;

  // Gitter lesen
  unsigned long nVertices, nPanels;
  vertex* Vertices;
  panel* Panels;
 

for(int iarg=1; iarg < argc; ++iarg) {
  //read a .vtk file 
  readVTKfile((char*) argv[iarg], nVertices, Vertices, nPanels, Panels);
  std::cout << "Read " << nPanels << " panels and " << nVertices << " points "
            << "from '" << argv[iarg] << "'." << std::endl;


  //double scale = 0.04;
  //eps = scale/nPanels;
  //eps = scale/(nPanels*sqrt(nPanels));
  //eps = scale/(nPanels*nPanels);
  std::cout <<"eps = " <<eps<<std::endl;

  /*
  BC* Pan_BC = new BC[nPanels];
  for (unsigned i=0; i<nPanels; ++i) Pan_BC[i].typ = 1;
  refinegrid(Panels, nPanels, Vertices, nVertices, Pan_BC);
  writegrid("tmp.grd", Panels, nPanels, Vertices, nVertices);
  */


  // maximale x-Koord. feststellen, um Testlsg. zu generieren
  xcenter = 0.0;
  for (unsigned i=0; i<nVertices; ++i) {
    double e = Vertices[i].point[0];
    if (xcenter<e) xcenter = e;
  }

  PRINT(xcenter);

  // --------------------------------------------------------------------------
  // Cluster-Baum der Panels erzeugen
  //

  std::cout << "Generating panel cluster tree ... " << std::flush;
  
  unsigned* op_perm_pan = new unsigned[nPanels];
  unsigned* po_perm_pan = new unsigned[nPanels];
  assert(op_perm_pan!=NULL && po_perm_pan!=NULL);
  for(unsigned i=0;i<nPanels;i++){
    op_perm_pan[i]=po_perm_pan[i]=i;
  }
  
  bemcluster<panel>* clTreePan = new bemcluster<panel>(Panels,op_perm_pan, 0, nPanels);
  clTreePan->createClusterTree(bmin, op_perm_pan, po_perm_pan);
  unsigned nClustersPan = clTreePan->getncl();
  std::cout << "done, " << nClustersPan << " clusters -- ";
  std::cout << inMB(nClustersPan*sizeof(bemcluster<panel>)) << " MB." << std::endl;

  /*  
      n_cluster<panel> *n_clTree = new n_cluster<panel>((panel**)clindxPan, clTreePan);
      std::cout << (int) n_clTree->have_cmn_nrml(clTreePan->getson(0), clTreePan->getson(1)) << std::endl;
      delete n_clTree;
      exit(0);
*/


  // --------------------------------------------------------------------------
  // Cluster-Baum der Vertices erzeugen
  //

  std::cout << "Generating vertex cluster tree ... " << std::flush;
  unsigned* op_perm_vtx = new unsigned[nVertices];
  unsigned* po_perm_vtx = new unsigned[nVertices];
  assert(op_perm_vtx!=NULL && po_perm_vtx!=NULL);
  for(unsigned i=0;i<nVertices;i++){
    op_perm_vtx[i]=po_perm_vtx[i]=i;
	Vertices[i].index=i;
  }

  /*

  // ------------------------------------------------------------------------
  // write mesh to NEUTRAL format file
  //
  
  std::ofstream osGRD("out.surf");
  osGRD << "surfacemesh" << std::endl;
  osGRD << nVertices << std::endl;
  for (unsigned j=0; j<nVertices; ++j) {
    osGRD << Vertices[j].point[0] << ' '
	  << Vertices[j].point[1] << ' '
	  << Vertices[j].point[2] << std::endl;
  }
  
  osGRD << nPanels << std::endl;
  for (unsigned i=0; i<nPanels; ++i)
    osGRD << (Panels[i].p1->idx)+1 << ' '
	  << (Panels[i].p2->idx)+1 << ' '
	  << (Panels[i].p3->idx)+1 << std::endl;
  
  osGRD.close();


  double *col = new double[nPanels];
  cluster* p = clTreePan->son1;
  for (unsigned i=p->nbeg; i<p->nend; ++i) col[clindxPan[i]->idx] = -1.0;
  p = clTreePan->son2->son1;
  for (unsigned i=p->nbeg; i<p->nend; ++i) col[clindxPan[i]->idx] = 0.0;
  p = clTreePan->son2->son2;
  for (unsigned i=p->nbeg; i<p->nend; ++i) col[clindxPan[i]->idx] = 1.0;

  std::ofstream osSOL("out.sol");
  
  osSOL << "solution Potential -size=" << nPanels
	<< " -type=surfaceelement" << std::endl;
  for (unsigned i=0; i<nPanels; ++i) osSOL << col[i] << std::endl;
  osSOL.close();

  exit(0);
  */

  /*
  double *v = new double[nVertices];

  for (unsigned ii=0; ii<nPanels; ++ii) {
    LDGCL3_bl(1, clindxPan+ii, nVertices, clindxVtx, v);
    double sum = 0.0;
    for (unsigned jj=0; jj<nVertices; ++jj) sum += v[jj];
    panel* p = (panel*) clindxPan[ii];
    double d = sum+0.5*p->area;
    std::cout << ii << ' ' << p->idx << ' ' << d << std::endl;
    if (fabs(d)>1e-3) {
      std::cout << p->p1->point << ' ' << p->p2->point << ' '
		<< p->p3->point << std::endl;
      exit(1);
    }
  }

  delete [] v;
  */

  bemcluster<vertex>* clTreeVtx = 
    new bemcluster<vertex>(Vertices,op_perm_vtx, 0, nVertices);

  clTreeVtx->createClusterTree(bmin, op_perm_vtx, po_perm_vtx);
  unsigned nClustersVtx = clTreeVtx->getncl();
  std::cout << "done, " << nClustersVtx << " clusters -- ";
  std::cout << inMB(nClustersVtx*sizeof(bemcluster<vertex>)) << " MB." << std::endl;

  // nun Panel-Vertex-Paare clustern
  std::cout << "Generating panel-vertex tree ... " << std::flush;

  bemblcluster<panel,vertex>* blclTreePanVtx = 
    new bemblcluster<panel,vertex>(0, 0, nPanels, nVertices);
  unsigned nblcksPanVtx = 0;
  blclTreePanVtx->subdivide(clTreePan, clTreeVtx, eta*eta, nblcksPanVtx);
  std::cout << "done, " << nblcksPanVtx << " blocks -- ";
  std::cout << inMB(blclTreePanVtx->size()) << " MB." << std::endl;

  // delete clTreeVtx;
  
  
  // --------------------------------------------------------------------------
  // Generate Double Layer Matrix
  //
  //

  double time0 = realtime(0.0);
  
  MATGENDLCL_DUF MatGenDL(Panels,op_perm_pan, Vertices,op_perm_vtx);
  mblock<double>** B;
  allocmbls(blclTreePanVtx, B);

  //matgenGeH_sqntl(MatGenDL, blclTreePanVtx, blclTreePanVtx, false, eps, rankmax, B);
  matgenGeH_omp(MatGenDL, nblcksPanVtx, blclTreePanVtx, eps, rankmax, B);

  // agglH(blclTreePanVtx, B, eps, rankmax);

  mem_orig = sizeof(double) * nPanels * nVertices;
  allmem = sizeH(blclTreePanVtx, B);
  std::cout << (char) 13 << "Needed storage " << inMB(allmem)
	    << " MB, without approximation " << inMB(mem_orig) << " MB.  "
	    << std::endl << "Compressed to " << 100.0*allmem/mem_orig << "%. "
	    << "Approximation took " << realtime(time0) << "s." << std::endl;

  if (wrtPart) {
    std::cout << "writing matrix partition ..." << std::flush;
    std::ofstream os("matrixB.ps");
    psoutputGeH(os, blclTreePanVtx, nPanels, B);
    os.close();
    std::cout << " done." << std::endl;
  }



  // --------------------------------------------------------------------------
  // generate Dirichlet Data
  //
 
  std::cout << "Generating Dirichlet data ... " << std::flush;
  double* DirData = new double[nVertices];
  genL_dirichlet_bc(DirData2, nVertices, nPanels, Panels,po_perm_vtx, DirData);
 
  // check error
  std::cout << "done. L2-error is "
	    << Dir_L2_Error(DirData2,nPanels, Panels,po_perm_vtx,DirData) << '.'
	    << std::endl;


  // --------------------------------------------------------------------------
  // generate right hand side
  //

  double* b = new double[nPanels];
  assert(b!=NULL);
  blas::setzero(nPanels, b);

  std::cout << "Multiplying Dirichlet data with double layer matrix ... "
	    << std::flush;
   
  mltaGeHVec(1.0, blclTreePanVtx, B, DirData, b);
  mlta_massCL(nPanels, Panels, op_perm_pan, po_perm_vtx, 0.5, DirData, b);
  std::cout << "done." << std::endl;

  delete [] DirData;

  // delete B
  freembls(blclTreePanVtx, B);
  delete blclTreePanVtx;


  
  // --------------------------------------------------------------------------
  // Generate Single Layer Matrix
  //

  // nun Panel-Paare clustern
  std::cout << "Generating panel-panel tree ... " << std::flush;
  bemblcluster<panel,panel>* blclTreePan = new bemblcluster<panel,panel>(0, 0, nPanels, nPanels);
  unsigned nblcksPan = 0;
  blclTreePan->subdivide_sym(clTreePan, eta*eta, nblcksPan);
  std::cout << "done, " << nblcksPan << " blocks -- ";
  std::cout << inMB(blclTreePan->size()) << " MB." << std::endl;

  //delete clTreePan;

  time0 = realtime(0.0);

  MATGENSLCC_DUF MatGenSL(Panels,op_perm_pan);

  // matrix A is a list of pointers mblocks
  BEMMatrix<panel,panel> A(nPanels, blclTreePan);
  //matgenHeH_sqntl(MatGenSL, A.blclTree, A.blclTree, false, eps, rankmax, A.blcks);
  matgenHeH_omp(MatGenSL, nblcksPan, A.blclTree, eps, rankmax, A.blcks);

  mem_orig = 0.5 * sizeof(double) * (nPanels+1) * nPanels;
  allmem = sizeH(A.blclTree, A.blcks);
  std::cout << (char) 13 << "Needed storage " << inMB(allmem)
	    << " MB, without approximation " << inMB(mem_orig) << " MB."
	    << std::endl << "Compressed to " << 100.0*allmem/mem_orig << "%. "
   	    << "Approximation took " << realtime(time0) << "s." << std::endl;

  if (wrtPart) {
    std::cout << "writing matrix partition ..." << std::flush;
    std::ofstream os("matrixA.ps");
    psoutputHeH(os, A.blclTree, nPanels, A.blcks);
    os.close();
    std::cout << " done." << std::endl;
  }

  /*
  std::cout << "Agglomerating matrix (eps=" << eps << ") ... " << std::flush;
  time0 = realtime(0.0);
  agglH(A.blclTree, A.blcks, eps, rankmax);
  std::cout << "done -- " << realtime(time0) << "s." << std::endl;

  if (wrtPart) {
    std::cout << "writing matrix partition ..." << std::flush;
    std::ofstream os("matrixAggl.ps");
    psoutputHeH(os, A.blclTree, nPanels, A.blcks);
    os.close();
    std::cout << " done." << std::endl;
  }

  allmem = sizeH(A.blclTree, A.blcks);
  std::cout << "Needed storage " << inMB(allmem)
	    << " MB, without approx. " << inMB(mem_orig) << " MB. "
	    << "Compressed to " << 100.0*allmem/mem_orig << "%. "
   	    << std::endl;
  */
  // --------------------------------------------------------------------------
  // Cholesky factorization
  //


  double delta = 0.1;
  PRINT(delta);

  std::cout << "Computing preconditioner ... " << std::flush;
  time0 = realtime(0.0);
  genCholprecond(A.blclTree, A.blcks, delta, rankmax, A.blU, A.U, true);
  std::cout << "done -- " << realtime(time0) << "s, "
	    << inMB(sizeH(A.blU, A.U, 'U')) << " MB." << std::endl;

  if (wrtPart) {
    std::cout << "writing matrix partitions ..." << std::flush;
    std::ofstream osU("matrixU.ps");
    psoutputUtH(osU, A.blU, nPanels, A.U);
    osU.close();
    std::cout << " done." << std::endl;
  }


  // --------------------------------------------------------------------------
  // Solve System
  //

  // generating right hand side
  time0 = realtime(0.0);

  double* x = new double[nPanels];
  assert(x!=NULL);
  blas::setzero(nPanels, x);

  double acc = 1E-14;
  unsigned steps = nPanels;

  // solve the system
  if (CG(A, b, x, acc, steps))
    std::cout << "CG_Cholesky: iteration did not converge.";
  else
    std::cout << "CG_Cholesky converged to " << acc << " in "
	      << steps << " steps.";

  std::cout << " Solution took " << realtime(time0) << "s." << std::endl;
 
  //Write vtk file
  //char colorMode= 'P';
  //writeVTKFile(nVertices, Vertices, nPanels, Panels, colorMode, x, po_perm_pan);


  // -------------------------------------------------------------------------
  // compare with exact solution
  //

  std::cout << "Validating the solution ..." << std::flush;
  
  double l2err = 0.0;
  double l2no = 0.0;
  double err = 0.0;
  double no = 0.0;

  for (unsigned i=0; i<nPanels; ++i) {
    panel* pan =&Panels[op_perm_pan[i]];
    vec3d v;
    v[0] = pan->getcenter(0);
    v[1] = pan->getcenter(1);
    v[2] = pan->getcenter(2);

    double nd = NeumData2(v, pan->normal);;

    err += SQR(nd-x[i]);
    no  += SQR(nd);

    l2err += pan->area * SQR(nd-x[i]);
    l2no  += pan->area * SQR(nd);
  }
  
  std::cout << " done, L2 rel. error " << sqrt(l2err/l2no)
	    << ", rel error " << sqrt(err/no) << std::endl;


  delete [] x;
  delete [] b;
  delete clTreeVtx;
  delete [] op_perm_pan;
  delete [] po_perm_pan;
  delete clTreePan;
  delete [] op_perm_vtx;
  delete [] po_perm_vtx;
  delete [] Panels;
  delete [] Vertices;
}
  return 0;
}
