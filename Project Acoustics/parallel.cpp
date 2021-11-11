#include "headers.hpp"

extern int myRank;
extern int nbTasks;

//================================================================================
// Build the local numbering and list of nodes for MPI communications
//================================================================================

void buildListsNodesMPI(Mesh& m)
{
  if(myRank == 0)
    printf("== build local numbering and list of nodes for MPI communications\n");

  //==== Build mask for nodes belonging to each MPI process (i.e. interior + interface)
  
  IntMatrix maskNodesEachProc(nbTasks, m.nbOfNodes);
  for(int nTask=0; nTask<nbTasks; nTask++){
    for(int nGlo=0; nGlo<m.nbOfNodes; nGlo++){
      maskNodesEachProc(nTask,nGlo) = 0;
    }
  }
  for(int iTriGlo=0; iTriGlo<m.nbOfElems; iTriGlo++){
    int nTask = m.elemNum(iTriGlo,2)-1;
    IntVector Elem = m.Elem2Nodes(iTriGlo);
    for(int i = 0; i < Elem.rows(); ++i){
      int nGlo = Elem(i)-1;
      maskNodesEachProc(nTask,nGlo) = 1;
    }
  }
  
  //==== Build local numbering for nodes belonging to the current MPI process (i.e. interior + interfaces)
  
  IntVector nodesGloToLoc(m.nbOfNodes);
  int nLoc = 1;
  for(int nGlo=0; nGlo<m.nbOfNodes; nGlo++){
    if(maskNodesEachProc(myRank, nGlo)){
      nodesGloToLoc(nGlo) = nLoc;  // this node belongs to the current MPI process
      nLoc++;
    }
    else{
      nodesGloToLoc(nGlo) = -1;  // this node does not belong to the current MPI process
    }
  }

  //==== Build list with nodes to exchange between the current MPI process and each neighboring MPI process (i.e. interfaces)
  
  m.numNodesToExch.resize(nbTasks);
  m.nodesToExch.resize(nbTasks,m.nbOfNodes);
  for(int nTask=0; nTask<nbTasks; nTask++){
    m.numNodesToExch(nTask) = 0;
    if(nTask != myRank){
      int count = 0;
      for(int nGlo=0; nGlo<m.nbOfNodes; nGlo++){
        if(maskNodesEachProc(myRank,nGlo) && maskNodesEachProc(nTask,nGlo)){
          m.nodesToExch(nTask,count) = nodesGloToLoc(nGlo);
          count++;
        }
      }
      m.numNodesToExch(nTask) = count;
      printf("   -> task %i send/recv %i nodes with task %i\n", myRank, m.numNodesToExch(nTask), nTask);
    }
  }
  
  if(m.numNodesToExch.maxCoeff() > 0){
    m.nodesToExch.conservativeResize(nbTasks, m.numNodesToExch.maxCoeff());
  }
  
  //==== Build local arrays for nodes/lines/triangles
  
  Matrix    coordsMyRank(m.nbOfNodes,2);
  IntMatrix elemNumMyRank(m.nbOfElems,3);
  SpIntMatrix elemNodesMyRank(m.nbOfElems,m.elemNodes.cols());
  IntVector BndNodesMyRank(m.nbOfBndNodes);
  
  nLoc = 0;
  int iEltLoc = 0;
  for(int nGlo=0; nGlo<m.nbOfNodes; nGlo++){
    if(nodesGloToLoc(nGlo) >= 0){
      coordsMyRank(nLoc,0) = m.coords(nGlo,0);
      coordsMyRank(nLoc,1) = m.coords(nGlo,1);
      nLoc++;
    }
  }
  for(int iElt=0; iElt<m.nbOfElems; iElt++){
    if(m.elemNum(iElt,2)-1 == myRank){
      elemNumMyRank(iEltLoc,0) = iEltLoc+1;
      for (int i = 1; i < 3; ++i)
        elemNumMyRank(iEltLoc,i) = m.elemNum(iElt,i);
      IntVector Elem = m.Elem2Nodes(iElt);
      for(int i = 0; i < Elem.rows(); ++i)
        elemNodesMyRank.coeffRef(iEltLoc,i) = nodesGloToLoc(Elem(i)-1); // Here
      iEltLoc++;
    }
  }

  int nbndLoc = 0;
  for(int nGlo=0; nGlo<m.nbOfBndNodes; nGlo++){
    if(nodesGloToLoc(m.BndNodes(nGlo)-1) >= 0){
      BndNodesMyRank(nbndLoc) = nodesGloToLoc(m.BndNodes(nGlo)-1);
      nbndLoc++;
    }
  }

  coordsMyRank.conservativeResize(nLoc,2);
  elemNumMyRank.conservativeResize(iEltLoc,3);
  elemNodesMyRank.conservativeResize(iEltLoc,m.elemNodes.cols());
  BndNodesMyRank.conservativeResize(nbndLoc);

  m.nbOfNodes = nLoc;
  m.nbOfElems = iEltLoc;
  m.nbOfBndNodes = nbndLoc;

  m.coords = coordsMyRank;
  m.elemNum = elemNumMyRank;
  m.elemNodes = elemNodesMyRank;
  m.BndNodes = BndNodesMyRank;
}

//================================================================================
// MPI-parallel exchange/add the interface terms
//================================================================================

void exchangeAddInterfMPI(Vector& vec, Mesh& m)
{
  MPI_Request *requestSnd;
  MPI_Request *requestRcv;
  MPI_Status status;
  requestSnd = new MPI_Request[nbTasks];
  requestRcv = new MPI_Request[nbTasks];
  
  double **bufferSnd;
  double **bufferRcv;
  bufferSnd = new double*[nbTasks];
  bufferRcv = new double*[nbTasks];
  
  for(int nTask=0; nTask<nbTasks; nTask++){
    int numToExch = m.numNodesToExch(nTask);
    if(numToExch > 0){
      bufferSnd[nTask] = new double[numToExch];
      bufferRcv[nTask] = new double[numToExch];
      for(int nExch=0; nExch<numToExch; nExch++)
        bufferSnd[nTask][nExch] = vec(m.nodesToExch(nTask,nExch)-1);
      MPI_Isend(bufferSnd[nTask], numToExch, MPI_DOUBLE, nTask, 0, MPI_COMM_WORLD, &requestSnd[nTask]);
      MPI_Irecv(bufferRcv[nTask], numToExch, MPI_DOUBLE, nTask, 0, MPI_COMM_WORLD, &requestRcv[nTask]);
    }
  }
  
  for(int nTask=0; nTask<nbTasks; nTask++){
    int numToExch = m.numNodesToExch(nTask);
    if(numToExch > 0){
      MPI_Wait(&requestRcv[nTask], &status);
      for(int nExch=0; nExch<numToExch; nExch++)
        vec(m.nodesToExch(nTask,nExch)-1) += bufferRcv[nTask][nExch];
      delete bufferRcv[nTask];
      MPI_Wait(&requestSnd[nTask], &status);
      delete bufferSnd[nTask];
    }
  }
  
  delete[] bufferSnd;
  delete[] bufferRcv;
  delete requestSnd;
  delete requestRcv;
}

//================================================================================
// MPI-parallel exchange/add (2) the interface terms
//================================================================================

void exchangeAddInterfMPI2(Vector& vec, Mesh& m)
{
  MPI_Request *requestSnd;
  MPI_Request *requestRcv;
  MPI_Status status;
  requestSnd = new MPI_Request[nbTasks];
  requestRcv = new MPI_Request[nbTasks];
  
  double **bufferSnd;
  double **bufferRcv;
  bufferSnd = new double*[nbTasks];
  bufferRcv = new double*[nbTasks];
  
  for(int nTask=0; nTask<nbTasks; nTask++){
    int numToExch = m.numNodesToExch(nTask);
    if(numToExch > 0){
      bufferSnd[nTask] = new double[numToExch];
      bufferRcv[nTask] = new double[numToExch];
      for(int nExch=0; nExch<numToExch; nExch++)
        bufferSnd[nTask][nExch] = vec(m.nodesToExch(nTask,nExch));
      MPI_Isend(bufferSnd[nTask], numToExch, MPI_DOUBLE, nTask, 0, MPI_COMM_WORLD, &requestSnd[nTask]);
      MPI_Irecv(bufferRcv[nTask], numToExch, MPI_DOUBLE, nTask, 0, MPI_COMM_WORLD, &requestRcv[nTask]);
    }
  }
  
  for(int nTask=0; nTask<nbTasks; nTask++){
    int numToExch = m.numNodesToExch(nTask);
    if(numToExch > 0){
      MPI_Wait(&requestRcv[nTask], &status);
      for(int nExch=0; nExch<numToExch; nExch++)
        vec(m.nodesToExch(nTask,nExch)-1) = bufferRcv[nTask][nExch];
      delete bufferRcv[nTask];
      MPI_Wait(&requestSnd[nTask], &status);
      delete bufferSnd[nTask];
    }
  }
  
  delete[] bufferSnd;
  delete[] bufferRcv;
  delete requestSnd;
  delete requestRcv;
}

