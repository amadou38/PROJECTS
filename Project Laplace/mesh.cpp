#include "headers.hpp"

extern int myRank;
extern int nbTasks;

//================================================================================
// Read the mesh from a gmsh-file (.msh) and store in a mesh-structure 'm'
//================================================================================

void readMsh(Mesh& m, string fileName)
{
  if(myRank == 0)
    printf("== read mesh file\n");
  
  // Open mesh file
  ifstream meshfile(fileName.c_str());
  if(!meshfile){
    printf("ERROR: Mesh '%s' not opened.\n", fileName.c_str());
    exit(EXIT_FAILURE);
  }
 
  // Read mesh file
  double dummy;
  while(meshfile.good()){
    
    string line;
    getline(meshfile,line);
    
    // Read NODES
    if(line.compare("$Nodes") == 0){
      
      // Read the number of nodes
      meshfile >> m.nbOfNodes;
      // Read the node coordinates
      m.coords.resize(m.nbOfNodes,2); // Review for 3D
      for(int nGlo=0; nGlo<m.nbOfNodes; ++nGlo)
        meshfile >> dummy >> m.coords(nGlo,0) >> m.coords(nGlo,1);// >> m.coords(nGlo,2); // Activate for 3D
      if(myRank == 0)
        printf("   -> %i nodes\n", m.nbOfNodes);
    }  
    // Read Interior NODES 
    if(line.compare("$IntNodes") == 0){
      
      // Read the number of nodes on the boundary
      meshfile >> m.nbOfIntNodes;
      
      // Read the nodes from the boundary
      m.IntNodes.resize(m.nbOfIntNodes); // Review for 3D
      for(int nGlo=0; nGlo<m.nbOfIntNodes; ++nGlo)
        meshfile >> m.IntNodes(nGlo);
    
      if(myRank == 0)
        printf("   -> %i nodes on the interior\n", m.nbOfIntNodes);
    }  
    // Read Boundary NODES 
    if(line.compare("$Boundary") == 0){
      
      // Read the number of nodes on the boundary
      meshfile >> m.nbOfBndNodes;
      
      // Read the nodes from the boundary
      m.BndNodes.resize(m.nbOfBndNodes); // Review for 3D
      for(int nGlo=0; nGlo<m.nbOfBndNodes; ++nGlo)
        meshfile >> m.BndNodes(nGlo);
    
      if(myRank == 0)
        printf("   -> %i nodes on the boundary\n", m.nbOfBndNodes);
    }
    // Read each edges from the boundary 
    if(line.compare("$Freeedg") == 0){
      
      // Read the number of edges from the boundary
      meshfile >> m.nbOfFreeedg;
      
      // Read the edges from the boundary
      m.Freeedg.resize(m.nbOfFreeedg); // Review for 3D
      for(int nGlo=0; nGlo<m.nbOfFreeedg; ++nGlo)
        meshfile >> m.Freeedg(nGlo);
    
      if(myRank == 0)
        printf("   -> %i of freeedges from the boundary\n", m.nbOfFreeedg);
    }
  }
  meshfile.close();
  // Read ELEMENTS
  IntMatrix elemNods;
  readInfo(m.nbOfElems, m.elemNum, m.elemNodes, fileName, "$Elements");
  if(myRank == 0)
    printf("   -> %i elements\n", m.nbOfElems);
  readInfo(m.nbOfEdges, m.edgNum, m.edgNodes, fileName, "$E2edg");
  if(myRank == 0)
    printf("   -> %i edges\n", m.nbOfEdges);
  readInfo(m.nbOfEneighs, m.eneighNum, m.eneighNodes, fileName, "$E2neigh");
  if(myRank == 0)
    printf("   -> %i eneighs\n", m.nbOfEneighs);
  readInfo(m.nbOfConns, m.connNum, m.connNodes, fileName, "$Connect");
  if(myRank == 0)
    printf("   -> %i connectics\n", m.nbOfConns);
}

//================================================================================
// Save a solution 'u' in a gmsh-file (.msh)
//================================================================================

void saveToMsh(Vector& u, Mesh& m, string name, string fileName)
{
  if(nbTasks > 1){
    ostringstream ss;
    ss << fileName << "_" << myRank;
    fileName = ss.str();
  }

  ofstream posFile(fileName.c_str());
  posFile << "$MeshFormat" << endl;
  posFile << "2.2 0 0" << endl;
  posFile << "$EndMeshFormat" << endl;
  
  posFile << "$ElementNodeData" << endl;
  posFile << "2" << endl;
  posFile << "\"" << name.c_str() << "\"" << endl;  // name of the view
  posFile << "" << endl;
  posFile << "1" << endl;
  posFile << "0" << endl; // ("Time")
  posFile << "4" << endl;
  posFile << "0" << endl; // ("timeStep")
  posFile << "1" << endl; // ("numComp")
  posFile << m.nbOfElems << endl;   // total number of elementNodeData in this file
  posFile << myRank << endl;
  for(int iEltLoc=0; iEltLoc<m.nbOfElems; iEltLoc++){
    // posFile << m.elemNum(iEltLoc,0)+m.nbOfBndNodes << " " << m.elemNum(iEltLoc,1);
    for(int n=0; n<m.elemNum(iEltLoc,1); n++){
      int nLoc = m.elemNodes.coeffRef(iEltLoc,n)-1;
      //cout<< "\nnloc: " << nLoc << endl;
      posFile << " " << u(nLoc);
    }
    posFile << endl;
  }
  posFile << "$EndElementNodeData" << endl;
  posFile.close();
}


void saveSol(Vector& u, string fileName)
{
  if(nbTasks > 1){
    ostringstream ss;
    ss << fileName << "_" << myRank;
    fileName = ss.str();
  }

  ofstream posFile(fileName.c_str());
  
  for(int i=0; i < u.rows(); i++)
    posFile << u(i) << endl;

  posFile.close();
}


void readInfo(int& nbOf, IntMatrix& Num, SpIntMatrix& Nodes, string fileName, string startName)
{
  ifstream meshfile(fileName.c_str());
  
  // Read mesh file
  double dummy;
  while(meshfile.good()){
    
  string line;
  getline(meshfile,line);
    
  if(line.compare(startName) == 0){
      
    // Read the total number of elements
    nbOf = 0;
    meshfile >> nbOf;
      
    // Temporary arrays for elements infos/nodes
    int nbOfmaxEdge = 12;        // maximum of edges per element
    Num.resize(nbOf,3);          // gmsh numerbing
    Nodes.resize(nbOf,nbOfmaxEdge);      // associated nodes
    
    // Read infos/nodes for all the elements
    IntVector nbOfEachElts(nbOfmaxEdge-1);
    for (int i = 0; i < nbOfEachElts.rows(); ++i)
      nbOfEachElts(i) = 0;
    for(int iGlo=0; iGlo<nbOf; iGlo++)
    {
      getline(meshfile, line);
      int infos;
      // Save element infos
      meshfile >> Num(iGlo,0);
      meshfile >> Num(iGlo,1);
      //meshfile >> infos;
      //meshfile >> dummy >> dummy;
      //if(infos > 2)
      //  meshfile >> dummy >> Num(iGlo,2); // partition tag
        meshfile >> Num(iGlo,2); // partition tag
      //else
      //  Num(iGlo,2) = 0;
      //for(int j=5; j<=infos; j++)
      //  meshfile >> dummy;                // useless infos
        
      // Check element type and save element nodes
      nbOfEachElts(Num(iGlo,1)-2)++;
      for (int i = 0; i < Num(iGlo,1); ++i)
        meshfile >> Nodes.coeffRef(iGlo,i);
    }
    int m = Num(0,1);
    for(int iGlo=1; iGlo<nbOf; iGlo++)
      if (m < Num(iGlo,1))
        m = Num(iGlo,1);
    Nodes.conservativeResize(nbOf,m);
    
    /*if(myRank == 0){
      for (int i = 0; i < nbOfEachElts.rows(); ++i)
        if (nbOfEachElts(i) != 0)
          printf(" Type with %d edges  -> %i elements\n", i+2, nbOfEachElts(i));
    }*/
  }
}
meshfile.close();   
}
