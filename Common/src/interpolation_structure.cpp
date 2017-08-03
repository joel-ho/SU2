/*!
 * \file interpolation_structure.cpp
 * \brief Main subroutines used by SU2_FSI
 * \author H. Kline
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/interpolation_structure.hpp"

CInterpolator::CInterpolator(void) {

  nZone = 0;
  Geometry = NULL;

  donor_geometry  = NULL;
  target_geometry = NULL;

  donorZone  = 0;
  targetZone = 0;

  Buffer_Receive_nVertex_Donor     = NULL;
  Buffer_Receive_nFace_Donor       = NULL;
  Buffer_Receive_nFaceNodes_Donor  = NULL;
  Buffer_Send_nVertex_Donor        = NULL;
  Buffer_Send_nFace_Donor          = NULL;
  Buffer_Send_nFaceNodes_Donor     = NULL;
  Buffer_Receive_GlobalPoint       = NULL;
  Buffer_Send_GlobalPoint          = NULL;
  Buffer_Send_FaceIndex            = NULL;
  Buffer_Receive_FaceIndex         = NULL;
  Buffer_Send_FaceNodes            = NULL;
  Buffer_Receive_FaceNodes         = NULL;
  Buffer_Send_FaceProc             = NULL;
  Buffer_Receive_FaceProc          = NULL;

  Buffer_Send_Coord                = NULL;
  Buffer_Send_Normal               = NULL;
  Buffer_Receive_Coord             = NULL;
  Buffer_Receive_Normal            = NULL;

}

CInterpolator::~CInterpolator(void) {

  //if (Buffer_Receive_nVertex_Donor!= NULL) delete[] Buffer_Receive_nVertex_Donor;
}


CInterpolator::CInterpolator(CGeometry ***geometry_container, CConfig **config, unsigned int iZone, unsigned int jZone) {

  /* Store pointers*/
  Geometry = geometry_container;

  donorZone  = iZone;
  targetZone = jZone;

  donor_geometry  = geometry_container[donorZone][MESH_0];
  target_geometry = geometry_container[targetZone][MESH_0];

  /*--- Initialize transfer coefficients between the zones ---*/
    /* Since this is a virtual function, call it in the child class constructor  */
  //Set_TransferCoeff(targetZone,donorZone,config);
  /*--- Initialize transfer coefficients between the zones ---*/
  //Set_TransferCoeff(Zones,config);

  //Buffer_Receive_nVertex_Donor = NULL;

}

inline void CInterpolator::Set_TransferCoeff(CConfig **config) { }

void CInterpolator::Determine_ArraySize(bool faces, int markDonor, int markTarget, unsigned long nVertexDonor, unsigned short nDim) {
  unsigned long nLocalVertex_Donor = 0, nLocalFaceNodes_Donor=0, nLocalFace_Donor=0;
  unsigned long iVertex, iPointDonor = 0;
  /* Only needed if face data is also collected */
  unsigned long inode;
  unsigned long donor_elem, jElem, jPoint;
  unsigned short iDonor;
  unsigned int nFaces=0, iFace, nNodes=0;
  bool face_on_marker = true;

#ifdef HAVE_MPI
  int rank = MASTER_NODE;
  int nProcessor = SINGLE_NODE;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif

  for (iVertex = 0; iVertex < nVertexDonor; iVertex++) {
    iPointDonor = donor_geometry->vertex[markDonor][iVertex]->GetNode();
    if (donor_geometry->node[iPointDonor]->GetDomain()) {
      nLocalVertex_Donor++;
      if (faces) {
        /*--- On Donor geometry also communicate face info ---*/
        if (nDim==3) {
          for (jElem=0; jElem<donor_geometry->node[iPointDonor]->GetnElem(); jElem++) {
            donor_elem = donor_geometry->node[iPointDonor]->GetElem(jElem);
            nFaces = donor_geometry->elem[donor_elem]->GetnFaces();
            for (iFace=0; iFace<nFaces; iFace++) {
              face_on_marker=true;
              nNodes = donor_geometry->elem[donor_elem]->GetnNodesFace(iFace);
              for (iDonor=0; iDonor<nNodes; iDonor++) {
                /*--- Local index of the node on face --*/
                inode = donor_geometry->elem[donor_elem]->GetFaces(iFace, iDonor);
                jPoint = donor_geometry->elem[donor_elem]->GetNode(inode);
                face_on_marker = (face_on_marker && (donor_geometry->node[jPoint]->GetVertex(markDonor) !=-1));
              }
              if (face_on_marker ) {
                nLocalFace_Donor++;
                nLocalFaceNodes_Donor+=nNodes;
              }
            }
          }
        }
        else {
          /*--- in 2D we use the edges ---*/
          nNodes=2;
          nFaces = donor_geometry->node[iPointDonor]->GetnPoint();
          for (iFace=0; iFace<nFaces; iFace++) {
            face_on_marker=true;
            for (iDonor=0; iDonor<nNodes; iDonor++) {
              inode = donor_geometry->node[iPointDonor]->GetEdge(iFace);
              jPoint = donor_geometry->edge[inode]->GetNode(iDonor);
              face_on_marker = (face_on_marker && (donor_geometry->node[jPoint]->GetVertex(markDonor) !=-1));
            }
            if (face_on_marker ) {
              nLocalFace_Donor++;
              nLocalFaceNodes_Donor+=nNodes;
            }
          }
        }
      }
    }
  }

  Buffer_Send_nVertex_Donor[0] = nLocalVertex_Donor;
  if (faces) {
    Buffer_Send_nFace_Donor[0] = nLocalFace_Donor;
    Buffer_Send_nFaceNodes_Donor[0] = nLocalFaceNodes_Donor;
  }

  /*--- Send Interface vertex information --*/
#ifdef HAVE_MPI
  SU2_MPI::Allreduce(&nLocalVertex_Donor, &MaxLocalVertex_Donor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
  SU2_MPI::Allgather(Buffer_Send_nVertex_Donor, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nVertex_Donor, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  if (faces) {
    SU2_MPI::Allreduce(&nLocalFace_Donor, &nGlobalFace_Donor, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nLocalFace_Donor, &MaxFace_Donor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nLocalFaceNodes_Donor, &nGlobalFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nLocalFaceNodes_Donor, &MaxFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_nFace_Donor, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nFace_Donor, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_nFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    MaxFace_Donor++;
  }
#else
  MaxLocalVertex_Donor    = nLocalVertex_Donor;
  Buffer_Receive_nVertex_Donor[0] = Buffer_Send_nVertex_Donor[0];
  if (faces) {
    nGlobalFace_Donor       = nLocalFace_Donor;
    nGlobalFaceNodes_Donor  = nLocalFaceNodes_Donor;
    MaxFaceNodes_Donor      = nLocalFaceNodes_Donor;
    MaxFace_Donor           = nLocalFace_Donor+1;
    Buffer_Receive_nFace_Donor[0] = Buffer_Send_nFace_Donor[0];
    Buffer_Receive_nFaceNodes_Donor[0] = Buffer_Send_nFaceNodes_Donor[0];
  }
#endif

}

void CInterpolator::Collect_VertexInfo(bool faces, int markDonor, int markTarget, unsigned long nVertexDonor, unsigned short nDim, unsigned long &nLocalVertex_Donor)
{
  unsigned long iVertex, iPointDonor = 0, iVertexDonor, nBuffer_Coord, nBuffer_Point;
  unsigned short iDim;
  
  nLocalVertex_Donor = 0;
  
  /* Only needed if face data is also collected */
  su2double  *Normal;

#ifdef HAVE_MPI
  int rank;
  int nProcessor = SINGLE_NODE;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
#endif


  for (iVertex = 0; iVertex < MaxLocalVertex_Donor; iVertex++) {
    Buffer_Send_GlobalPoint[iVertex] = 0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Buffer_Send_Coord[iVertex*nDim+iDim] = 0.0;
      if (faces)
        Buffer_Send_Normal[iVertex*nDim+iDim] = 0.0;
    }
  }

  /*--- Copy coordinates and point to the auxiliar vector --*/
  nLocalVertex_Donor = 0;

  for (iVertexDonor = 0; iVertexDonor < nVertexDonor; iVertexDonor++) {
    iPointDonor = donor_geometry->vertex[markDonor][iVertexDonor]->GetNode();
    if (donor_geometry->node[iPointDonor]->GetDomain()) {
      Buffer_Send_GlobalPoint[nLocalVertex_Donor] = donor_geometry->node[iPointDonor]->GetGlobalIndex();
      for (iDim = 0; iDim < nDim; iDim++)
        Buffer_Send_Coord[nLocalVertex_Donor*nDim+iDim] = donor_geometry->node[iPointDonor]->GetCoord(iDim);

      if (faces) {
        Normal =  donor_geometry->vertex[markDonor][iVertexDonor]->GetNormal();
        for (iDim = 0; iDim < nDim; iDim++)
          Buffer_Send_Normal[nLocalVertex_Donor*nDim+iDim] = Normal[iDim];
      }
      nLocalVertex_Donor++;
    }
  }
  nBuffer_Coord = MaxLocalVertex_Donor*nDim;
  nBuffer_Point = MaxLocalVertex_Donor;

#ifdef HAVE_MPI
  SU2_MPI::Allgather(Buffer_Send_Coord, nBuffer_Coord, MPI_DOUBLE, Buffer_Receive_Coord, nBuffer_Coord, MPI_DOUBLE, MPI_COMM_WORLD);
  SU2_MPI::Allgather(Buffer_Send_GlobalPoint, nBuffer_Point, MPI_UNSIGNED_LONG, Buffer_Receive_GlobalPoint, nBuffer_Point, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
  if (faces) {
    SU2_MPI::Allgather(Buffer_Send_Normal, nBuffer_Coord, MPI_DOUBLE, Buffer_Receive_Normal, nBuffer_Coord, MPI_DOUBLE, MPI_COMM_WORLD);
  }
#else
  for (iVertex = 0; iVertex < nBuffer_Coord; iVertex++)
    Buffer_Receive_Coord[iVertex] = Buffer_Send_Coord[iVertex];

  for (iVertex = 0; iVertex < nBuffer_Point; iVertex++)
    Buffer_Receive_GlobalPoint[iVertex] = Buffer_Send_GlobalPoint[iVertex];

  if (faces) {
    for (iVertex = 0; iVertex < nBuffer_Coord; iVertex++)
      Buffer_Receive_Normal[iVertex] = Buffer_Send_Normal[iVertex];
  }
#endif
}

void CInterpolator::Collect_VertexInfo(bool faces, int markDonor, int markTarget, unsigned long nVertexDonor, unsigned short nDim) {

	unsigned long  nLocalVertex_Donor = 0;
	
	Collect_VertexInfo(faces, markDonor, markTarget, nVertexDonor, nDim, nLocalVertex_Donor);
}

int CInterpolator::Find_InterfaceMarker(CConfig *config, unsigned short val_marker_interface) {
    
  unsigned short nMarker = config->GetnMarker_All();
  unsigned short iMarker;

  for (iMarker = 0; iMarker < nMarker; iMarker++) {

    /*--- If the tag GetMarker_All_FSIinterface(iMarker) equals the index we are looping at ---*/
    if (config->GetMarker_All_FSIinterface(iMarker) == val_marker_interface ) {

      /*--- We have identified the identifier for the interface marker ---*/
      return iMarker;
    }
  }

  return -1;
}


/* Nearest Neighbor Interpolator */
CNearestNeighbor::CNearestNeighbor(void):  CInterpolator() { }

CNearestNeighbor::CNearestNeighbor(CGeometry ***geometry_container, CConfig **config,  unsigned int iZone, unsigned int jZone) :  CInterpolator(geometry_container, config, iZone, jZone) {

  /*--- Initialize transfer coefficients between the zones ---*/
  Set_TransferCoeff(config);

}

CNearestNeighbor::~CNearestNeighbor() {}

void CNearestNeighbor::Set_TransferCoeff(CConfig **config) {

  int iProcessor, pProcessor, nProcessor;
  int markDonor, markTarget, Target_check, Donor_check;

  unsigned short iDim, nDim, iMarkerInt, nMarkerInt, iDonor;    

  unsigned long nVertexDonor, nVertexTarget, Point_Target, jVertex, iVertexTarget;
  unsigned long Global_Point_Donor, pGlobalPoint=0;

  su2double *Coord_i, Coord_j[3], dist, mindist, maxdist;

#ifdef HAVE_MPI

  int rank = MASTER_NODE;
  int *Buffer_Recv_mark=NULL, iRank;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  
  if (rank == MASTER_NODE) 
    Buffer_Recv_mark = new int[nProcessor];

#else

  nProcessor = SINGLE_NODE;

#endif

  /*--- Initialize variables --- */
  
  nMarkerInt = (int) ( config[donorZone]->GetMarker_n_FSIinterface() / 2 );
  
  nDim = donor_geometry->GetnDim();

  iDonor = 0;
  
  Buffer_Receive_nVertex_Donor = new unsigned long [nProcessor];


  /*--- Cycle over nMarkersInt interface to determine communication pattern ---*/

  for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; iMarkerInt++) {

    /*--- On the donor side: find the tag of the boundary sharing the interface ---*/
    markDonor  = Find_InterfaceMarker(config[donorZone],  iMarkerInt);
      
    /*--- On the target side: find the tag of the boundary sharing the interface ---*/
    markTarget = Find_InterfaceMarker(config[targetZone], iMarkerInt);

    #ifdef HAVE_MPI
    
    Donor_check  = -1;
    Target_check = -1;
        
    /*--- We gather a vector in MASTER_NODE to determines whether the boundary is not on the processor because of the partition or because the zone does not include it ---*/
    
    SU2_MPI::Gather(&markDonor , 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
    
    if (rank == MASTER_NODE)
      for (iRank = 0; iRank < nProcessor; iRank++)
        if( Buffer_Recv_mark[iRank] != -1 ) {
          Donor_check = Buffer_Recv_mark[iRank];
          break;
        }
    
    SU2_MPI::Bcast(&Donor_check , 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
    
    
    SU2_MPI::Gather(&markTarget, 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
    
    if (rank == MASTER_NODE)
      for (iRank = 0; iRank < nProcessor; iRank++)
        if( Buffer_Recv_mark[iRank] != -1 ) {
          Target_check = Buffer_Recv_mark[iRank];
          break;
        }

    SU2_MPI::Bcast(&Target_check, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
        
    #else
    Donor_check  = markDonor;
    Target_check = markTarget;  
    #endif
    
    /*--- Checks if the zone contains the interface, if not continue to the next step ---*/
    if(Target_check == -1 || Donor_check == -1)
      continue;

    if(markDonor != -1)
      nVertexDonor  = donor_geometry->GetnVertex( markDonor );
    else
      nVertexDonor  = 0;
    
    if(markTarget != -1)
      nVertexTarget = target_geometry->GetnVertex( markTarget );
    else
      nVertexTarget  = 0;
    
    Buffer_Send_nVertex_Donor  = new unsigned long [ 1 ];

    /* Sets MaxLocalVertex_Donor, Buffer_Receive_nVertex_Donor */
    Determine_ArraySize(false, markDonor, markTarget, nVertexDonor, nDim);

    Buffer_Send_Coord          = new su2double     [ MaxLocalVertex_Donor * nDim ];
    Buffer_Send_GlobalPoint    = new unsigned long [ MaxLocalVertex_Donor ];
    Buffer_Receive_Coord       = new su2double     [ nProcessor * MaxLocalVertex_Donor * nDim ];
    Buffer_Receive_GlobalPoint = new unsigned long [ nProcessor * MaxLocalVertex_Donor ];

    /*-- Collect coordinates, global points, and normal vectors ---*/
    Collect_VertexInfo( false, markDonor, markTarget, nVertexDonor, nDim );

    /*--- Compute the closest point to a Near-Field boundary point ---*/
    maxdist = 0.0;

    for (iVertexTarget = 0; iVertexTarget < nVertexTarget; iVertexTarget++) {

      Point_Target = target_geometry->vertex[markTarget][iVertexTarget]->GetNode();

      if ( target_geometry->node[Point_Target]->GetDomain() ) {

        target_geometry->vertex[markTarget][iVertexTarget]->SetnDonorPoints(1);
        target_geometry->vertex[markTarget][iVertexTarget]->Allocate_DonorInfo(); // Possible meme leak?

        /*--- Coordinates of the boundary point ---*/
        Coord_i = target_geometry->node[Point_Target]->GetCoord();

        mindist    = 1E6; 
        pProcessor = 0;

        /*--- Loop over all the boundaries to find the pair ---*/

        for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
          for (jVertex = 0; jVertex < MaxLocalVertex_Donor; jVertex++) {
            Global_Point_Donor = iProcessor*MaxLocalVertex_Donor+jVertex;

            /*--- Compute the dist ---*/
            dist = 0.0; 
            for (iDim = 0; iDim < nDim; iDim++) {
              Coord_j[iDim] = Buffer_Receive_Coord[ Global_Point_Donor*nDim+iDim];
              dist += pow(Coord_j[iDim] - Coord_i[iDim], 2.0);
            }

            if (dist < mindist) {
              mindist = dist; pProcessor = iProcessor; pGlobalPoint = Buffer_Receive_GlobalPoint[Global_Point_Donor];
            }

            if (dist == 0.0) break;
          }

        } 

        /*--- Store the value of the pair ---*/
        maxdist = max(maxdist, mindist);
        target_geometry->vertex[markTarget][iVertexTarget]->SetInterpDonorPoint(iDonor, pGlobalPoint);
        target_geometry->vertex[markTarget][iVertexTarget]->SetInterpDonorProcessor(iDonor, pProcessor);
        target_geometry->vertex[markTarget][iVertexTarget]->SetDonorCoeff(iDonor, 1.0);
      }
    }

    delete[] Buffer_Send_Coord;
    delete[] Buffer_Send_GlobalPoint;
    
    delete[] Buffer_Receive_Coord;
    delete[] Buffer_Receive_GlobalPoint;

    delete[] Buffer_Send_nVertex_Donor;

  }

  delete[] Buffer_Receive_nVertex_Donor;

  #ifdef HAVE_MPI
  if (rank == MASTER_NODE) 
    delete [] Buffer_Recv_mark;
  #endif
}



CIsoparametric::CIsoparametric(CGeometry ***geometry_container, CConfig **config, unsigned int iZone, unsigned int jZone)  :  CInterpolator(geometry_container, config, iZone, jZone) {

  /*--- Initialize transfer coefficients between the zones ---*/
  Set_TransferCoeff(config);

  /*--- For fluid-structure interaction data interpolated with have nDim dimensions ---*/
 // InitializeData(Zones,nDim);
}

CIsoparametric::~CIsoparametric() {}

void CIsoparametric::Set_TransferCoeff(CConfig **config) {
  unsigned long iVertex, jVertex;
  unsigned long  dPoint, inode, jElem, nElem;
  unsigned short iDim, iDonor=0, iFace;

  unsigned short nDim = donor_geometry->GetnDim();

  unsigned short nMarkerInt;
  unsigned short iMarkerInt;

  int markDonor=0, markTarget=0;
  int Target_check, Donor_check;

  long donor_elem=0, temp_donor=0;
  unsigned int nNodes=0;
  /*--- Restricted to 2-zone for now ---*/
  unsigned int nFaces=1; //For 2D cases, we want to look at edges, not faces, as the 'interface'
  bool face_on_marker=true;

  unsigned long nVertexDonor = 0, nVertexTarget= 0;
  unsigned long Point_Target = 0;

  unsigned long iVertexDonor, iPointDonor = 0;
  unsigned long jGlobalPoint = 0;
  int iProcessor;

  unsigned long nLocalFace_Donor = 0, nLocalFaceNodes_Donor=0;

  unsigned long faceindex;

  su2double dist = 0.0, mindist=1E6, *Coord, *Coord_i;
  su2double myCoeff[10]; // Maximum # of donor points
  su2double  *Normal;
  su2double *projected_point = new su2double[nDim];
  su2double tmp, tmp2;
  su2double storeCoeff[10];
  unsigned long storeGlobal[10];
  int storeProc[10];

  int rank = MASTER_NODE;
  int nProcessor = SINGLE_NODE;
  Coord = new su2double[nDim];
  Normal = new su2double[nDim];

#ifdef HAVE_MPI

  int *Buffer_Recv_mark=NULL, iRank;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  
  if (rank == MASTER_NODE) 
    Buffer_Recv_mark = new int[nProcessor];

#else

  nProcessor = SINGLE_NODE;

#endif

  /*--- Number of markers on the FSI interface ---*/
  nMarkerInt = (config[donorZone]->GetMarker_n_FSIinterface())/2;

  /*--- For the number of markers on the interface... ---*/
  for (iMarkerInt=1; iMarkerInt <= nMarkerInt; iMarkerInt++) {
    /*--- Procedure:
    * -Loop through vertices of the aero grid
    * -Find nearest element and allocate enough space in the aero grid donor point info
    *    -set the transfer coefficient values
    */

    /*--- On the donor side: find the tag of the boundary sharing the interface ---*/
    markDonor  = Find_InterfaceMarker(config[donorZone],  iMarkerInt);
      
    /*--- On the target side: find the tag of the boundary sharing the interface ---*/
    markTarget = Find_InterfaceMarker(config[targetZone], iMarkerInt);

    #ifdef HAVE_MPI
    
    Donor_check  = -1;
    Target_check = -1;
        
    /*--- We gather a vector in MASTER_NODE to determines whether the boundary is not on the processor because of the partition or because the zone does not include it ---*/
    
    SU2_MPI::Gather(&markDonor , 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
    
    if (rank == MASTER_NODE)
      for (iRank = 0; iRank < nProcessor; iRank++)
        if( Buffer_Recv_mark[iRank] != -1 ) {
          Donor_check = Buffer_Recv_mark[iRank];
          break;
        }
    
    SU2_MPI::Bcast(&Donor_check , 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
    
    SU2_MPI::Gather(&markTarget, 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
    
    if (rank == MASTER_NODE)
      for (iRank = 0; iRank < nProcessor; iRank++)
        if( Buffer_Recv_mark[iRank] != -1 ) {
          Target_check = Buffer_Recv_mark[iRank];
          break;
        }

    SU2_MPI::Bcast(&Target_check, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
        
    #else
    Donor_check  = markDonor;
    Target_check = markTarget;  
    #endif
    
    /*--- Checks if the zone contains the interface, if not continue to the next step ---*/
    if(Target_check == -1 || Donor_check == -1)
      continue;

    if(markDonor != -1)
      nVertexDonor  = donor_geometry->GetnVertex( markDonor );
    else
      nVertexDonor  = 0;

    if(markTarget != -1)
      nVertexTarget = target_geometry->GetnVertex( markTarget );
    else
      nVertexTarget  = 0;
    
    Buffer_Send_nVertex_Donor    = new unsigned long [1];
    Buffer_Send_nFace_Donor      = new unsigned long [1];
    Buffer_Send_nFaceNodes_Donor = new unsigned long [1];

    Buffer_Receive_nVertex_Donor    = new unsigned long [nProcessor];
    Buffer_Receive_nFace_Donor      = new unsigned long [nProcessor];
    Buffer_Receive_nFaceNodes_Donor = new unsigned long [nProcessor];

    /* Sets MaxLocalVertex_Donor, Buffer_Receive_nVertex_Donor */
    Determine_ArraySize(true, markDonor, markTarget, nVertexDonor, nDim);

    Buffer_Send_Coord       = new su2double [MaxLocalVertex_Donor*nDim];
    Buffer_Send_Normal      = new su2double [MaxLocalVertex_Donor*nDim];
    Buffer_Send_GlobalPoint = new unsigned long [MaxLocalVertex_Donor];

    Buffer_Receive_Coord       = new su2double [nProcessor*MaxLocalVertex_Donor*nDim];
    Buffer_Receive_Normal      = new su2double [nProcessor*MaxLocalVertex_Donor*nDim];
    Buffer_Receive_GlobalPoint = new unsigned long [nProcessor*MaxLocalVertex_Donor];

    /*-- Collect coordinates, global points, and normal vectors ---*/
    Collect_VertexInfo(true, markDonor,markTarget,nVertexDonor,nDim);

    Buffer_Send_FaceIndex    = new unsigned long[MaxFace_Donor];
    Buffer_Send_FaceNodes    = new unsigned long[MaxFaceNodes_Donor];
    Buffer_Send_FaceProc     = new unsigned long[MaxFaceNodes_Donor];

    Buffer_Receive_FaceIndex = new unsigned long[MaxFace_Donor*nProcessor];
    Buffer_Receive_FaceNodes = new unsigned long[MaxFaceNodes_Donor*nProcessor];
    Buffer_Receive_FaceProc  = new unsigned long[MaxFaceNodes_Donor*nProcessor];

    nLocalFace_Donor=0;
    nLocalFaceNodes_Donor=0;

    /*--- Collect Face info ---*/

    for (iVertex = 0; iVertex < MaxFace_Donor; iVertex++) {
      Buffer_Send_FaceIndex[iVertex] = 0;
    }
    for (iVertex=0; iVertex<MaxFaceNodes_Donor; iVertex++) {
      Buffer_Send_FaceNodes[iVertex] = 0;
      Buffer_Send_FaceProc[iVertex]  = 0;
    }

    Buffer_Send_FaceIndex[0] = rank * MaxFaceNodes_Donor;

    if (nDim==2) nNodes=2;

    for (iVertexDonor = 0; iVertexDonor < nVertexDonor; iVertexDonor++) {
      iPointDonor = donor_geometry->vertex[markDonor][iVertexDonor]->GetNode();

      if (donor_geometry->node[iPointDonor]->GetDomain()) {

    if (nDim==3)  nElem = donor_geometry->node[iPointDonor]->GetnElem();
    else          nElem =donor_geometry->node[iPointDonor]->GetnPoint();

    for (jElem=0; jElem < nElem; jElem++) {
      if (nDim==3) {
        temp_donor = donor_geometry->node[iPointDonor]->GetElem(jElem);
        nFaces = donor_geometry->elem[temp_donor]->GetnFaces();
        for (iFace=0; iFace<nFaces; iFace++) {
          /*-- Determine whether this face/edge is on the marker --*/
          face_on_marker=true;
          nNodes = donor_geometry->elem[temp_donor]->GetnNodesFace(iFace);
          for (iDonor=0; iDonor<nNodes; iDonor++) {
            inode = donor_geometry->elem[temp_donor]->GetFaces(iFace, iDonor);
            dPoint = donor_geometry->elem[temp_donor]->GetNode(inode);
            face_on_marker = (face_on_marker && (donor_geometry->node[dPoint]->GetVertex(markDonor) !=-1));
          }

          if (face_on_marker ) {
            for (iDonor=0; iDonor<nNodes; iDonor++) {
              inode = donor_geometry->elem[temp_donor]->GetFaces(iFace, iDonor);
              dPoint = donor_geometry->elem[temp_donor]->GetNode(inode);
              // Match node on the face to the correct global index
              jGlobalPoint=donor_geometry->node[dPoint]->GetGlobalIndex();
              for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
                for (jVertex = 0; jVertex < Buffer_Receive_nVertex_Donor[iProcessor]; jVertex++) {
                  if (jGlobalPoint ==Buffer_Receive_GlobalPoint[MaxLocalVertex_Donor*iProcessor+jVertex]) {
                    Buffer_Send_FaceNodes[nLocalFaceNodes_Donor]=MaxLocalVertex_Donor*iProcessor+jVertex;
                    Buffer_Send_FaceProc[nLocalFaceNodes_Donor]=iProcessor;
                  }
                }
              }
              nLocalFaceNodes_Donor++; // Increment total number of face-nodes / processor
            }
            /* Store the indices */
            Buffer_Send_FaceIndex[nLocalFace_Donor+1] = Buffer_Send_FaceIndex[nLocalFace_Donor]+nNodes;
            nLocalFace_Donor++; // Increment number of faces / processor
          }
        }
      }
      else {
        /*-- Determine whether this face/edge is on the marker --*/
        face_on_marker=true;
        for (iDonor=0; iDonor<nNodes; iDonor++) {
          inode = donor_geometry->node[iPointDonor]->GetEdge(jElem);
          dPoint = donor_geometry->edge[inode]->GetNode(iDonor);
          face_on_marker = (face_on_marker && (donor_geometry->node[dPoint]->GetVertex(markDonor) !=-1));
        }
        if (face_on_marker ) {
          for (iDonor=0; iDonor<nNodes; iDonor++) {
            inode = donor_geometry->node[iPointDonor]->GetEdge(jElem);
            dPoint = donor_geometry->edge[inode]->GetNode(iDonor);
            // Match node on the face to the correct global index
            jGlobalPoint=donor_geometry->node[dPoint]->GetGlobalIndex();
            for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
              for (jVertex = 0; jVertex < Buffer_Receive_nVertex_Donor[iProcessor]; jVertex++) {
                if (jGlobalPoint ==Buffer_Receive_GlobalPoint[MaxLocalVertex_Donor*iProcessor+jVertex]) {
                  Buffer_Send_FaceNodes[nLocalFaceNodes_Donor]=MaxLocalVertex_Donor*iProcessor+jVertex;
                  Buffer_Send_FaceProc[nLocalFaceNodes_Donor]=iProcessor;
                }
              }
            }
            nLocalFaceNodes_Donor++; // Increment total number of face-nodes / processor
          }
          /* Store the indices */
          Buffer_Send_FaceIndex[nLocalFace_Donor+1] = Buffer_Send_FaceIndex[nLocalFace_Donor]+nNodes;
          nLocalFace_Donor++; // Increment number of faces / processor
        }
      }
    }
      }
    }

    //Buffer_Send_FaceIndex[nLocalFace_Donor+1] = MaxFaceNodes_Donor*rank+nLocalFaceNodes_Donor;
#ifdef HAVE_MPI
    SU2_MPI::Allgather(Buffer_Send_FaceNodes, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG, Buffer_Receive_FaceNodes, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_FaceProc, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG, Buffer_Receive_FaceProc, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_FaceIndex, MaxFace_Donor, MPI_UNSIGNED_LONG, Buffer_Receive_FaceIndex, MaxFace_Donor, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
    for (iFace=0; iFace<MaxFace_Donor; iFace++) {
      Buffer_Receive_FaceIndex[iFace] = Buffer_Send_FaceIndex[iFace];
    }
    for (iVertex = 0; iVertex < MaxFaceNodes_Donor; iVertex++)
      Buffer_Receive_FaceNodes[iVertex] = Buffer_Send_FaceNodes[iVertex];
    for (iVertex = 0; iVertex < MaxFaceNodes_Donor; iVertex++)
      Buffer_Receive_FaceProc[iVertex] = Buffer_Send_FaceProc[iVertex];
#endif

    /*--- Loop over the vertices on the target Marker ---*/
    for (iVertex = 0; iVertex<nVertexTarget; iVertex++) {
      mindist=1E6;
      for (unsigned short iCoeff=0; iCoeff<10; iCoeff++) {
    storeCoeff[iCoeff]=0;
      }
      Point_Target = target_geometry->vertex[markTarget][iVertex]->GetNode();

      if (target_geometry->node[Point_Target]->GetDomain()) {

    Coord_i = target_geometry->node[Point_Target]->GetCoord();
    /*---Loop over the faces previously communicated/stored ---*/
    for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {

      nFaces = (unsigned int)Buffer_Receive_nFace_Donor[iProcessor];

      for (iFace = 0; iFace< nFaces; iFace++) {
        /*--- ---*/

        nNodes = (unsigned int)Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace+1] -
                (unsigned int)Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace];

        su2double *X = new su2double[nNodes*nDim];
        faceindex = Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace]; // first index of this face
        for (iDonor=0; iDonor<nNodes; iDonor++) {
          jVertex = Buffer_Receive_FaceNodes[iDonor+faceindex]; // index which points to the stored coordinates, global points
          for (iDim=0; iDim<nDim; iDim++) {
            X[iDim*nNodes+iDonor]=
                Buffer_Receive_Coord[jVertex*nDim+iDim];
          }
        }
        jVertex = Buffer_Receive_FaceNodes[faceindex];

        for (iDim=0; iDim<nDim; iDim++) {
          Normal[iDim] = Buffer_Receive_Normal[jVertex*nDim+iDim];
        }

        /* Project point used for case where surfaces are not exactly coincident, where
         * the point is assumed connected by a rigid rod normal to the surface.
         */
        tmp = 0;
        tmp2=0;
        for (iDim=0; iDim<nDim; iDim++) {
          tmp+=Normal[iDim]*Normal[iDim];
          tmp2+=Normal[iDim]*(Coord_i[iDim]-X[iDim*nNodes]);
        }
        tmp = 1/tmp;
        tmp2 = tmp2*sqrt(tmp);
        for (iDim=0; iDim<nDim; iDim++) {
          // projection of \vec{q} onto plane defined by \vec{n} and \vec{p}:
          // \vec{q} - \vec{n} ( (\vec{q}-\vec{p} ) \cdot \vec{n})
          // tmp2 = ( (\vec{q}-\vec{p} ) \cdot \vec{N})
          // \vec{n} = \vec{N}/(|N|), tmp = 1/|N|^2
          projected_point[iDim]=Coord_i[iDim] + Normal[iDim]*tmp2*tmp;
        }

        Isoparameters(nDim, nNodes, X, projected_point,myCoeff);

        /*--- Find distance to the interpolated point ---*/
        dist = 0.0;
        for (iDim=0; iDim<nDim; iDim++) {
          Coord[iDim] = Coord_i[iDim];
          for(iDonor=0; iDonor< nNodes; iDonor++) {
            Coord[iDim]-=myCoeff[iDonor]*X[iDim*nNodes+iDonor];
          }
          dist+=pow(Coord[iDim],2.0);
        }

        /*--- If the dist is shorter than last closest (and nonzero nodes are on the boundary), update ---*/
        if (dist<mindist ) {
          /*--- update last dist ---*/
          mindist = dist;
          /*--- Store info ---*/
          donor_elem = temp_donor;
          target_geometry->vertex[markTarget][iVertex]->SetDonorElem(donor_elem); // in 2D is nearest neighbor
          target_geometry->vertex[markTarget][iVertex]->SetnDonorPoints(nNodes);
          for (iDonor=0; iDonor<nNodes; iDonor++) {
            storeCoeff[iDonor] = myCoeff[iDonor];
            jVertex = Buffer_Receive_FaceNodes[faceindex+iDonor];
            storeGlobal[iDonor] =Buffer_Receive_GlobalPoint[jVertex];
            storeProc[iDonor] = (int)Buffer_Receive_FaceProc[faceindex+iDonor];
          }
        }
      
        delete [] X;
      }
    }
    /*--- Set the appropriate amount of memory and fill ---*/
    nNodes =target_geometry->vertex[markTarget][iVertex]->GetnDonorPoints();
    target_geometry->vertex[markTarget][iVertex]->Allocate_DonorInfo();

    for (iDonor=0; iDonor<nNodes; iDonor++) {
      target_geometry->vertex[markTarget][iVertex]->SetInterpDonorPoint(iDonor,storeGlobal[iDonor]);
      //cout <<rank << " Global Point " << Global_Point<<" iDonor " << iDonor <<" coeff " << coeff <<" gp " << pGlobalPoint << endl;
      target_geometry->vertex[markTarget][iVertex]->SetDonorCoeff(iDonor,storeCoeff[iDonor]);
      target_geometry->vertex[markTarget][iVertex]->SetInterpDonorProcessor(iDonor, storeProc[iDonor]);
    }
      }
    }

    delete[] Buffer_Send_nVertex_Donor;
    delete[] Buffer_Send_nFace_Donor;
    delete[] Buffer_Send_nFaceNodes_Donor;

    delete[] Buffer_Receive_nVertex_Donor;
    delete[] Buffer_Receive_nFace_Donor;
    delete[] Buffer_Receive_nFaceNodes_Donor;

    delete[] Buffer_Send_Coord;
    delete[] Buffer_Send_Normal;
    delete[] Buffer_Send_GlobalPoint;

    delete[] Buffer_Receive_Coord;
    delete[] Buffer_Receive_Normal;
    delete[] Buffer_Receive_GlobalPoint;

    delete[] Buffer_Send_FaceIndex;
    delete[] Buffer_Send_FaceNodes;
    delete[] Buffer_Send_FaceProc;

    delete[] Buffer_Receive_FaceIndex;
    delete[] Buffer_Receive_FaceNodes;
    delete[] Buffer_Receive_FaceProc;
  }
  delete [] Coord;
  delete [] Normal;
  
  delete [] projected_point;

  #ifdef HAVE_MPI
  if (rank == MASTER_NODE) 
    delete [] Buffer_Recv_mark;
  #endif
}

void CIsoparametric::Isoparameters(unsigned short nDim, unsigned short nDonor,
    su2double *X, su2double *xj, su2double *isoparams) {
  short iDonor,iDim,k; // indices
  su2double tmp, tmp2;
  
  su2double *x     = new su2double[nDim+1];
  su2double *x_tmp = new su2double[nDim+1];
  su2double *Q     = new su2double[nDonor*nDonor];
  su2double *R     = new su2double[nDonor*nDonor];
  su2double *A     = new su2double[nDim+1*nDonor];
  su2double *A2    = NULL;
  su2double *x2    = new su2double[nDim+1];
  
  bool *test  = new bool[nDim+1];
  bool *testi = new bool[nDim+1];
  
  su2double eps = 1E-10;
  
  short n = nDim+1;

  if (nDonor>2) {
    /*--- Create Matrix A: 1st row all 1's, 2nd row x coordinates, 3rd row y coordinates, etc ---*/
    /*--- Right hand side is [1, \vec{x}']'---*/
    for (iDonor=0; iDonor<nDonor; iDonor++) {
      isoparams[iDonor]=0;
      A[iDonor] = 1.0;
      for (iDim=0; iDim<n; iDim++)
        A[(iDim+1)*nDonor+iDonor]=X[iDim*nDonor+iDonor];
    }

    x[0] = 1.0;
    for (iDim=0; iDim<nDim; iDim++)
      x[iDim+1]=xj[iDim];

    /*--- Eliminate degenerate rows:
     * for example, if z constant including the z values will make the system degenerate
     * TODO: improve efficiency of this loop---*/
    test[0]=true; // always keep the 1st row
    for (iDim=1; iDim<nDim+1; iDim++) {
      // Test this row against all previous
      test[iDim]=true; // Assume that it is not degenerate
      for (k=0; k<iDim; k++) {
        tmp=0; tmp2=0;
        for (iDonor=0;iDonor<nDonor;iDonor++) {
          tmp+= A[iDim*nDonor+iDonor]*A[iDim*nDonor+iDonor];
          tmp2+=A[k*nDonor+iDonor]*A[k*nDonor+iDonor];
        }
        tmp  = pow(tmp,0.5);
        tmp2 = pow(tmp2,0.5);
        testi[k]=false;
        for (iDonor=0; iDonor<nDonor; iDonor++) {
          // If at least one ratio is non-matching row iDim is not degenerate w/ row k
          if (A[iDim*nDonor+iDonor]/tmp != A[k*nDonor+iDonor]/tmp2)
            testi[k]=true;
        }
        // If any of testi (k<iDim) are false, row iDim is degenerate
        test[iDim]=(test[iDim] && testi[k]);
      }
      if (!test[iDim]) n--;
    }

    /*--- Initialize A2 now that we might have a smaller system --*/
    A2 = new su2double[n*nDonor];
    iDim=0;
    /*--- Copy only the rows that are non-degenerate ---*/
    for (k=0; k<nDim+1; k++) {
      if (test[k]) {
        for (iDonor=0;iDonor<nDonor;iDonor++ ) {
          A2[nDonor*iDim+iDonor]=A[nDonor*k+iDonor];
        }
        x2[iDim]=x[k];
        iDim++;
      }
    }
    /*--- Initialize Q,R to 0 --*/
    for (k=0; k<nDonor*nDonor; k++) {
      Q[k]=0;
      R[k]=0;
    }
    /*--- TODO: make this loop more efficient ---*/
    /*--- Solve for rectangular Q1 R1 ---*/
    for (iDonor=0; iDonor<nDonor; iDonor++) {
      tmp=0;
      for (iDim=0; iDim<n; iDim++)
        tmp += (A2[iDim*nDonor+iDonor])*(A2[iDim*nDonor+iDonor]);

      R[iDonor*nDonor+iDonor]= pow(tmp,0.5);
      if (tmp>eps && iDonor<n) {
        for (iDim=0; iDim<n; iDim++)
          Q[iDim*nDonor+iDonor]=A2[iDim*nDonor+iDonor]/R[iDonor*nDonor+iDonor];
      }
      else if (tmp!=0) {
        for (iDim=0; iDim<n; iDim++)
          Q[iDim*nDonor+iDonor]=A2[iDim*nDonor+iDonor]/tmp;
      }
      for (iDim=iDonor+1; iDim<nDonor; iDim++) {
        tmp=0;
        for (k=0; k<n; k++)
          tmp+=A2[k*nDonor+iDim]*Q[k*nDonor+iDonor];

        R[iDonor*nDonor+iDim]=tmp;

        for (k=0; k<n; k++)
          A2[k*nDonor+iDim]=A2[k*nDonor+iDim]-Q[k*nDonor+iDonor]*R[iDonor*nDonor+iDim];
      }
    }
    /*--- x_tmp = Q^T * x2 ---*/
    for (iDonor=0; iDonor<nDonor; iDonor++)
      x_tmp[iDonor]=0.0;
    for (iDonor=0; iDonor<nDonor; iDonor++) {
      for (iDim=0; iDim<n; iDim++)
        x_tmp[iDonor]+=Q[iDim*nDonor+iDonor]*x2[iDim];
    }

    /*--- solve x_tmp = R*isoparams for isoparams: upper triangular system ---*/
    for (iDonor = n-1; iDonor>=0; iDonor--) {
      if (R[iDonor*nDonor+iDonor]>eps)
        isoparams[iDonor]=x_tmp[iDonor]/R[iDonor*nDonor+iDonor];
      else
        isoparams[iDonor]=0;
      for (k=0; k<iDonor; k++)
        x_tmp[k]=x_tmp[k]-R[k*nDonor+iDonor]*isoparams[iDonor];
    }
  }
  else {
    /*-- For 2-donors (lines) it is simpler: */
    tmp =  pow(X[0*nDonor+0]- X[0*nDonor+1],2.0);
    tmp += pow(X[1*nDonor+0]- X[1*nDonor+1],2.0);
    tmp = sqrt(tmp);

    tmp2 = pow(X[0*nDonor+0] - xj[0],2.0);
    tmp2 += pow(X[1*nDonor+0] - xj[1],2.0);
    tmp2 = sqrt(tmp2);
    isoparams[1] = tmp2/tmp;

    tmp2 = pow(X[0*nDonor+1] - xj[0],2.0);
    tmp2 += pow(X[1*nDonor+1] - xj[1],2.0);
    tmp2 = sqrt(tmp2);
    isoparams[0] = tmp2/tmp;
  }

  /*--- Isoparametric coefficients have been calculated. Run checks to eliminate outside-element issues ---*/
  if (nDonor==4) {
    /*-- Bilinear coordinates, bounded by [-1,1] ---*/
    su2double xi, eta;
    xi = (1.0-isoparams[0]/isoparams[1])/(1.0+isoparams[0]/isoparams[1]);
    eta = 1- isoparams[2]*4/(1+xi);
    if (xi>1.0) xi=1.0;
    if (xi<-1.0) xi=-1.0;
    if (eta>1.0) eta=1.0;
    if (eta<-1.0) eta=-1.0;
    isoparams[0]=0.25*(1-xi)*(1-eta);
    isoparams[1]=0.25*(1+xi)*(1-eta);
    isoparams[2]=0.25*(1+xi)*(1+eta);
    isoparams[3]=0.25*(1-xi)*(1+eta);

  }
  if (nDonor<4) {
    tmp = 0.0; // value for normalization
    tmp2=0; // check for maximum value, to be used to id nearest neighbor if necessary
    k=0; // index for maximum value
    for (iDonor=0; iDonor< nDonor; iDonor++) {
      if (isoparams[iDonor]>tmp2) {
        k=iDonor;
        tmp2=isoparams[iDonor];
      }
      // [0,1]
      if (isoparams[iDonor]<0) isoparams[iDonor]=0;
      if (isoparams[iDonor]>1) isoparams[iDonor] = 1;
      tmp +=isoparams[iDonor];
    }
    if (tmp>0)
      for (iDonor=0; iDonor< nDonor; iDonor++)
        isoparams[iDonor]=isoparams[iDonor]/tmp;
    else {
      isoparams[k] = 1.0;
    }
  }
  
  delete [] x;
  delete [] x_tmp;
  delete [] Q;
  delete [] R;
  delete [] A;
  if (A2 != NULL) delete [] A2;
  delete [] x2;
  
  delete [] test;
  delete [] testi;

}


/* Mirror Interpolator */
CMirror::CMirror(CGeometry ***geometry_container, CConfig **config,  unsigned int iZone, unsigned int jZone) :  CInterpolator(geometry_container, config, iZone, jZone) {

  /*--- Initialize transfer coefficients between the zones ---*/
  Set_TransferCoeff(config);

}

CMirror::~CMirror() {}

void CMirror::Set_TransferCoeff(CConfig **config) {
  unsigned long iVertex, jVertex;
  unsigned long iPoint;
  unsigned short iDonor=0, iFace=0, iTarget=0;

  unsigned short nMarkerInt;
  unsigned short iMarkerInt;

  int markDonor=0, markTarget=0;
  int Target_check, Donor_check;

  unsigned int nNodes=0, iNodes=0;
  unsigned long nVertexDonor = 0, nVertexTarget= 0;
  unsigned long Point_Donor = 0;
  unsigned long Global_Point = 0;
  unsigned long pGlobalPoint = 0;
  int iProcessor;

  unsigned long nLocalFace_Donor = 0, nLocalFaceNodes_Donor=0;

  unsigned long faceindex;

  int rank = MASTER_NODE;
  int nProcessor = SINGLE_NODE;

#ifdef HAVE_MPI
  int *Buffer_Recv_mark=NULL, iRank;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  
  if (rank == MASTER_NODE) 
    Buffer_Recv_mark = new int[nProcessor];

#else

  nProcessor = SINGLE_NODE;

#endif

  su2double *Buffer_Send_Coeff, *Buffer_Receive_Coeff;
  su2double coeff;

  /*--- Number of markers on the interface ---*/
  nMarkerInt = (config[targetZone]->GetMarker_n_FSIinterface())/2;

  /*--- For the number of markers on the interface... ---*/
  for (iMarkerInt=1; iMarkerInt <= nMarkerInt; iMarkerInt++) {
   /*--- Procedure:
    * -Loop through vertices of the aero grid
    * -Find nearest element and allocate enough space in the aero grid donor point info
    *    -set the transfer coefficient values
    */

    /*--- On the donor side: find the tag of the boundary sharing the interface ---*/
    markDonor  = Find_InterfaceMarker(config[donorZone],  iMarkerInt);

    /*--- On the target side: find the tag of the boundary sharing the interface ---*/
    markTarget = Find_InterfaceMarker(config[targetZone], iMarkerInt);

    #ifdef HAVE_MPI

    Donor_check  = -1;
    Target_check = -1;

    /*--- We gather a vector in MASTER_NODE to determines whether the boundary is not on the processor because of the partition or because the zone does not include it ---*/

    SU2_MPI::Gather(&markDonor , 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    if (rank == MASTER_NODE)
      for (iRank = 0; iRank < nProcessor; iRank++)
        if( Buffer_Recv_mark[iRank] != -1 ) {
          Donor_check = Buffer_Recv_mark[iRank];
          break;
        }

    SU2_MPI::Bcast(&Donor_check , 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    SU2_MPI::Gather(&markTarget, 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    if (rank == MASTER_NODE)
      for (iRank = 0; iRank < nProcessor; iRank++)
        if( Buffer_Recv_mark[iRank] != -1 ) {
          Target_check = Buffer_Recv_mark[iRank];
          break;
        }

    SU2_MPI::Bcast(&Target_check, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

    #else
    Donor_check  = markDonor;
    Target_check = markTarget;  
    #endif

    /*--- Checks if the zone contains the interface, if not continue to the next step ---*/
    if(Target_check == -1 || Donor_check == -1)
      continue;

    if(markDonor != -1)
      nVertexDonor  = donor_geometry->GetnVertex( markDonor );
    else
      nVertexDonor  = 0;

    if(markTarget != -1)
      nVertexTarget = target_geometry->GetnVertex( markTarget );
    else
      nVertexTarget  = 0;

    /*-- Collect the number of donor nodes: re-use 'Face' containers --*/
    nLocalFace_Donor=0;
    nLocalFaceNodes_Donor=0;
    for (jVertex = 0; jVertex<nVertexDonor; jVertex++) {
      Point_Donor =donor_geometry->vertex[markDonor][jVertex]->GetNode(); // Local index of jVertex

      if (donor_geometry->node[Point_Donor]->GetDomain()) {
        nNodes = donor_geometry->vertex[markDonor][jVertex]->GetnDonorPoints();
        nLocalFaceNodes_Donor+=nNodes;
        nLocalFace_Donor++;
      }
    }
    Buffer_Send_nFace_Donor= new unsigned long [1];
    Buffer_Send_nFaceNodes_Donor= new unsigned long [1];

    Buffer_Receive_nFace_Donor = new unsigned long [nProcessor];
    Buffer_Receive_nFaceNodes_Donor = new unsigned long [nProcessor];

    Buffer_Send_nFace_Donor[0] = nLocalFace_Donor;
    Buffer_Send_nFaceNodes_Donor[0] = nLocalFaceNodes_Donor;

    /*--- Send Interface vertex information --*/
#ifdef HAVE_MPI
    SU2_MPI::Allreduce(&nLocalFaceNodes_Donor, &MaxFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nLocalFace_Donor, &MaxFace_Donor, 1, MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_nFace_Donor, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nFace_Donor, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_nFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, Buffer_Receive_nFaceNodes_Donor, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    MaxFace_Donor++;
#else
    nGlobalFace_Donor       = nLocalFace_Donor;
    nGlobalFaceNodes_Donor  = nLocalFaceNodes_Donor;
    MaxFaceNodes_Donor      = nLocalFaceNodes_Donor;
    MaxFace_Donor           = nLocalFace_Donor+1;
    Buffer_Receive_nFace_Donor[0] = Buffer_Send_nFace_Donor[0];
    Buffer_Receive_nFaceNodes_Donor[0] = Buffer_Send_nFaceNodes_Donor[0];
#endif

    /*-- Send donor info --*/
    Buffer_Send_FaceIndex   = new unsigned long[MaxFace_Donor];
    Buffer_Send_FaceNodes   = new unsigned long[MaxFaceNodes_Donor];
    //Buffer_Send_FaceProc    = new unsigned long[MaxFaceNodes_Donor];
    Buffer_Send_GlobalPoint = new unsigned long[MaxFaceNodes_Donor];
    Buffer_Send_Coeff       = new su2double[MaxFaceNodes_Donor];

    Buffer_Receive_FaceIndex= new unsigned long[MaxFace_Donor*nProcessor];
    Buffer_Receive_FaceNodes= new unsigned long[MaxFaceNodes_Donor*nProcessor];
    //Buffer_Receive_FaceProc = new unsigned long[MaxFaceNodes_Donor*nProcessor];
    Buffer_Receive_GlobalPoint = new unsigned long[MaxFaceNodes_Donor*nProcessor];
    Buffer_Receive_Coeff    = new su2double[MaxFaceNodes_Donor*nProcessor];

    for (iVertex=0; iVertex<MaxFace_Donor; iVertex++) {
      Buffer_Send_FaceIndex[iVertex]=0;
    }
    for (iVertex=0; iVertex<MaxFaceNodes_Donor; iVertex++) {
      Buffer_Send_FaceNodes[iVertex]=0;
      //Buffer_Send_FaceProc[iVertex]=0;
      Buffer_Send_GlobalPoint[iVertex]=0;
      Buffer_Send_Coeff[iVertex]=0.0;
    }
    for (iVertex=0; iVertex<MaxFace_Donor; iVertex++) {
      Buffer_Send_FaceIndex[iVertex]=0;
    }

    Buffer_Send_FaceIndex[0]=rank*MaxFaceNodes_Donor;
    nLocalFace_Donor=0;
    nLocalFaceNodes_Donor=0;

    for (jVertex = 0; jVertex<nVertexDonor; jVertex++) {

      Point_Donor =donor_geometry->vertex[markDonor][jVertex]->GetNode(); // Local index of jVertex
      if (donor_geometry->node[Point_Donor]->GetDomain()) {
        nNodes = donor_geometry->vertex[markDonor][jVertex]->GetnDonorPoints();
        for (iDonor=0; iDonor<nNodes; iDonor++) {
          Buffer_Send_FaceNodes[nLocalFaceNodes_Donor] = donor_geometry->node[Point_Donor]->GetGlobalIndex();
          Buffer_Send_GlobalPoint[nLocalFaceNodes_Donor] =
              donor_geometry->vertex[markDonor][jVertex]->GetInterpDonorPoint(iDonor);
          Buffer_Send_Coeff[nLocalFaceNodes_Donor] =
              donor_geometry->vertex[markDonor][jVertex]->GetDonorCoeff(iDonor);
          nLocalFaceNodes_Donor++;
        }
        Buffer_Send_FaceIndex[nLocalFace_Donor+1] =Buffer_Send_FaceIndex[nLocalFace_Donor]+nNodes;
        nLocalFace_Donor++;
      }
    }

#ifdef HAVE_MPI
    SU2_MPI::Allgather(Buffer_Send_FaceNodes, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG, Buffer_Receive_FaceNodes, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_GlobalPoint, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG,Buffer_Receive_GlobalPoint, MaxFaceNodes_Donor, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_Coeff, MaxFaceNodes_Donor, MPI_DOUBLE,Buffer_Receive_Coeff, MaxFaceNodes_Donor, MPI_DOUBLE, MPI_COMM_WORLD);
    SU2_MPI::Allgather(Buffer_Send_FaceIndex, MaxFace_Donor, MPI_UNSIGNED_LONG, Buffer_Receive_FaceIndex, MaxFace_Donor, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
    for (iFace=0; iFace<MaxFace_Donor; iFace++) {
      Buffer_Receive_FaceIndex[iFace] = Buffer_Send_FaceIndex[iFace];
    }
    for (iVertex = 0; iVertex < MaxFaceNodes_Donor; iVertex++) {
      Buffer_Receive_FaceNodes[iVertex] = Buffer_Send_FaceNodes[iVertex];
      Buffer_Receive_GlobalPoint[iVertex] = Buffer_Send_GlobalPoint[iVertex];
      Buffer_Receive_Coeff[iVertex] = Buffer_Send_Coeff[iVertex];
    }
#endif
    /*--- Loop over the vertices on the target Marker ---*/
    for (iVertex = 0; iVertex<nVertexTarget; iVertex++) {

      iPoint = target_geometry->vertex[markTarget][iVertex]->GetNode();
      if (target_geometry->node[iPoint]->GetDomain()) {
        Global_Point = target_geometry->node[iPoint]->GetGlobalIndex();
        nNodes = 0;
        for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
          for (iFace = 0; iFace < Buffer_Receive_nFace_Donor[iProcessor]; iFace++) {
            faceindex = Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace]; // first index of this face
            iNodes = (unsigned int)Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace+1]- (unsigned int)faceindex;
            for (iTarget=0; iTarget<iNodes; iTarget++) {
              if (Global_Point == Buffer_Receive_GlobalPoint[faceindex+iTarget])
                nNodes++;
              //coeff =Buffer_Receive_Coeff[faceindex+iDonor];
            }
          }
        }

        target_geometry->vertex[markTarget][iVertex]->SetnDonorPoints(nNodes);
        target_geometry->vertex[markTarget][iVertex]->Allocate_DonorInfo();

        iDonor = 0;
        for (iProcessor = 0; iProcessor < nProcessor; iProcessor++) {
          for (iFace = 0; iFace < Buffer_Receive_nFace_Donor[iProcessor]; iFace++) {

            faceindex = Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace]; // first index of this face
            iNodes = (unsigned int)Buffer_Receive_FaceIndex[iProcessor*MaxFace_Donor+iFace+1]- (unsigned int)faceindex;
            for (iTarget=0; iTarget<iNodes; iTarget++) {
              if (Global_Point == Buffer_Receive_GlobalPoint[faceindex+iTarget]) {
                coeff =Buffer_Receive_Coeff[faceindex+iTarget];
                pGlobalPoint = Buffer_Receive_FaceNodes[faceindex+iTarget];
                target_geometry->vertex[markTarget][iVertex]->SetInterpDonorPoint(iDonor,pGlobalPoint);
                target_geometry->vertex[markTarget][iVertex]->SetDonorCoeff(iDonor,coeff);
                target_geometry->vertex[markTarget][iVertex]->SetInterpDonorProcessor(iDonor, iProcessor);
                //cout <<rank << " Global Point " << Global_Point<<" iDonor " << iDonor <<" coeff " << coeff <<" gp " << pGlobalPoint << endl;
                iDonor++;
              }
            }
          }
        }
      }
    }
    delete[] Buffer_Send_nFace_Donor;
    delete[] Buffer_Send_nFaceNodes_Donor;

    delete[] Buffer_Receive_nFace_Donor;
    delete[] Buffer_Receive_nFaceNodes_Donor;

    delete[] Buffer_Send_FaceIndex;
    delete[] Buffer_Send_FaceNodes;
    delete[] Buffer_Send_GlobalPoint;
    delete[] Buffer_Send_Coeff;

    delete[] Buffer_Receive_FaceIndex;
    delete[] Buffer_Receive_FaceNodes;
    delete[] Buffer_Receive_GlobalPoint;
    delete[] Buffer_Receive_Coeff;

  }

  #ifdef HAVE_MPI
  if (rank == MASTER_NODE) 
    delete [] Buffer_Recv_mark;
  #endif
}

/*--- Radial Basis Function Interpolator ---*/
CRadialBasisFunction::CRadialBasisFunction(void):  CInterpolator() { }

CRadialBasisFunction::CRadialBasisFunction(CGeometry ***geometry_container, CConfig **config,  unsigned int iZone, unsigned int jZone) :  CInterpolator(geometry_container, config, iZone, jZone) {

  /*--- Initialize transfer coefficients between the zones ---*/
  Set_TransferCoeff(config);

}

CRadialBasisFunction::~CRadialBasisFunction() {}

void CRadialBasisFunction::Set_TransferCoeff(CConfig **config) {

	int rank=MASTER_NODE, nProcessor=SINGLE_NODE;
  int iProcessor, pProcessor;
  int nPolynomial;
  int mark_donor, mark_target, target_check, donor_check;
  int *skip_row, *calc_polynomial_check;

  unsigned short iDim, nDim, iMarkerInt, nMarkerInt;    

  unsigned long iVertexDonor, jVertexDonor, iVertexTarget, iLocalM, iCount, jCount;
  unsigned long nVertexDonor, nVertexTarget, nVertexDonorInDomain;
  unsigned long nGlobalVertexDonor, iGlobalVertexDonor_end, nLocalM;
  unsigned long point_donor, point_target;
  unsigned long *nVertexDonorInDomain_arr, *nLocalM_arr;
  
  su2double val_i, val_j;
  su2double interface_coord_tol=1e6*numeric_limits<double>::epsilon();
  su2double *Coord_i, *Coord_j;
  su2double *local_M, *global_M_val_arr, *Buffer_recv_local_M;
  su2double *P;
  su2double *C_inv_trunc, *C_tmp;
  su2double *target_vec, *coeff_vec;
  
  su2double *P_limits, *P_tmp;
  
  CSymmetricMatrix *global_M, *Mp;

#ifdef HAVE_MPI

  int *Buffer_Recv_mark=NULL, iRank;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nProcessor);
  
  if (rank == MASTER_NODE) 
    Buffer_Recv_mark = new int[nProcessor];

#endif

  /*--- Initialize variables --- */
  
  nMarkerInt = (int) ( config[donorZone]->GetMarker_n_FSIinterface() / 2 );
  
  nDim = donor_geometry->GetnDim();

  
  Buffer_Receive_nVertex_Donor = new unsigned long [nProcessor];


  /*--- Cycle over nMarkersInt interface to determine communication pattern ---*/

  for (iMarkerInt = 1; iMarkerInt <= nMarkerInt; iMarkerInt++) {

    /*--- On the donor side: find the tag of the boundary sharing the interface ---*/
    mark_donor  = Find_InterfaceMarker(config[donorZone],  iMarkerInt);
      
    /*--- On the target side: find the tag of the boundary sharing the interface ---*/
    mark_target = Find_InterfaceMarker(config[targetZone], iMarkerInt);

#ifdef HAVE_MPI

    donor_check  = -1;
    target_check = -1;
        
    /*--- We gather a vector in MASTER_NODE to determines whether the boundary is not on the processor because of the partition or because the zone does not include it ---*/
    
    SU2_MPI::Gather(&mark_donor , 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
    
    if (rank == MASTER_NODE)
      for (iRank = 0; iRank < nProcessor; iRank++)
        if( Buffer_Recv_mark[iRank] != -1 ) {
          donor_check = Buffer_Recv_mark[iRank];
          break;
        }
    
    SU2_MPI::Bcast(&donor_check , 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
    
    
    SU2_MPI::Gather(&mark_target, 1, MPI_INT, Buffer_Recv_mark, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
    
    if (rank == MASTER_NODE)
      for (iRank = 0; iRank < nProcessor; iRank++)
        if( Buffer_Recv_mark[iRank] != -1 ) {
          target_check = Buffer_Recv_mark[iRank];
          break;
        }

    SU2_MPI::Bcast(&target_check, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
        
#else
    donor_check  = mark_donor;
    target_check = mark_target;  
#endif
    
    /*--- Checks if the zone contains the interface, if not continue to the next step ---*/
    if(target_check == -1 || donor_check == -1)
      continue;

    if(mark_donor != -1)
      nVertexDonor  = donor_geometry->GetnVertex( mark_donor );
    else
      nVertexDonor  = 0;
    
    if(mark_target != -1)
      nVertexTarget = target_geometry->GetnVertex( mark_target );
    else
      nVertexTarget  = 0;
    
    Buffer_Send_nVertex_Donor  = new unsigned long [ 1 ];

    /*--- Sets MaxLocalVertex_Donor, Buffer_Receive_nVertex_Donor ---*/
    Determine_ArraySize(false, mark_donor, mark_target, nVertexDonor, nDim);

    Buffer_Send_Coord          = new su2double     [ MaxLocalVertex_Donor * nDim ];
    Buffer_Send_GlobalPoint    = new unsigned long [ MaxLocalVertex_Donor ];
    Buffer_Receive_Coord       = new su2double     [ nProcessor * MaxLocalVertex_Donor * nDim ];
    Buffer_Receive_GlobalPoint = new unsigned long [ nProcessor * MaxLocalVertex_Donor ];

    /*-- Collect coordinates, global points, and normal vectors ---*/
    Collect_VertexInfo( false, mark_donor, mark_target, nVertexDonor, nDim, nVertexDonorInDomain);
    
    /*--- Collect information about number of donor vertices in domain ---*/
    /*--- Also calculate total number of donor vertices across all processors in domain on boundary ---*/
    nVertexDonorInDomain_arr = new unsigned long [nProcessor];
#ifdef HAVE_MPI
    SU2_MPI::Allgather(&nVertexDonorInDomain, 1, MPI_UNSIGNED_LONG, nVertexDonorInDomain_arr, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
    SU2_MPI::Allreduce(&nVertexDonorInDomain, &nGlobalVertexDonor, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
    nVertexDonorInDomain_arr[MASTER_NODE] = nVertexDonorInDomain;
    nGlobalVertexDonor = nVertexDonorInDomain;
#endif
    
		/*--- Calculate number of vertices on boundary prior to current processor ---*/
    iGlobalVertexDonor_end = 0;
    for (iProcessor=0; iProcessor<=rank; iProcessor++) {
    	iGlobalVertexDonor_end += nVertexDonorInDomain_arr[iProcessor];
    }
    
    /*--- Send information about size of local_M array ---*/
		nLocalM = nVertexDonorInDomain*(nVertexDonorInDomain+1)/2 \
		  + nVertexDonorInDomain*(nGlobalVertexDonor-iGlobalVertexDonor_end);
		
    nLocalM_arr = new unsigned long [nProcessor];
#ifdef HAVE_MPI
    SU2_MPI::Allgather(&nLocalM, 1, MPI_UNSIGNED_LONG, nLocalM_arr, 1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
#else
    nLocalM_arr[MASTER_NODE] = nLocalM;
#endif
    
    /*--- Initialize local M array and calculate values ---*/
    if (rank == MASTER_NODE) { cout << "Assembling donor zone radial basis function matrix... " << endl; }
    local_M = new su2double [nLocalM];  
    Coord_i = new su2double [nDim];
    Coord_j = new su2double [nDim];
    iCount=0;
    for (iVertexDonor=0; iVertexDonor<nVertexDonorInDomain; iVertexDonor++) {
    	for (iDim=0; iDim<nDim; iDim++) { Coord_i[iDim] = Buffer_Send_Coord[iVertexDonor*nDim + iDim]; }
    	
    	for (jVertexDonor=iVertexDonor; jVertexDonor<nVertexDonorInDomain; jVertexDonor++) { 
	  		for (iDim=0; iDim<nDim; iDim++) { Coord_j[iDim] = Buffer_Send_Coord[jVertexDonor*nDim + iDim]; }

	    		Get_Distance(Coord_i, Coord_j, nDim, val_i);
			    Get_RadialBasisValue(val_i, config[donorZone]);
			    local_M[iCount] = val_i;
			    iCount++;
	        
    	 }
			if (rank+1 < nProcessor ) {
		  	for (iProcessor=rank+1; iProcessor<nProcessor; iProcessor++) {
		  		for (jVertexDonor=0; jVertexDonor<nVertexDonorInDomain_arr[iProcessor]; jVertexDonor++) {
			  		for (iDim=0; iDim<nDim; iDim++) { Coord_j[iDim] = Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_Donor+jVertexDonor)*nDim + iDim]; }    		
							Get_Distance(Coord_i, Coord_j, nDim, val_i);
						  Get_RadialBasisValue(val_i, config[donorZone]);
						  local_M[iCount] = val_i;
						  iCount++;
		  		}
		  	}
    	}
    }
    
#ifdef HAVE_MPI
    if (rank != MASTER_NODE) {
    	SU2_MPI::Send(local_M, nLocalM, MPI_DOUBLE, MASTER_NODE, 0, MPI_COMM_WORLD);
    }
    
    /*--- Assemble global_M ---*/
    if (rank == MASTER_NODE) {
      global_M_val_arr = new su2double [nGlobalVertexDonor*(nGlobalVertexDonor+1)/2];
    	
    	/*--- Copy master node local_M to global_M ---*/
    	iCount = 0;
    	for (iLocalM=0; iLocalM<nLocalM; iLocalM++) {
    		global_M_val_arr[iCount] = local_M[iLocalM];
    		iCount++;
    	}
    	
    	/*--- Receive local_M from other processors ---*/
			if (nProcessor > SINGLE_NODE) {
		  	for (iProcessor=1; iProcessor<nProcessor; iProcessor++) {
			  	Buffer_recv_local_M = new su2double[nLocalM_arr[iProcessor]];
			  	SU2_MPI::Recv(Buffer_recv_local_M, nLocalM_arr[iProcessor], MPI_DOUBLE, iProcessor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			  	
			  	/*--- Copy processor's local_M to global_M ---*/
			  	for (iLocalM=0; iLocalM<nLocalM_arr[iProcessor]; iLocalM++) {
			  		global_M_val_arr[iCount] = Buffer_recv_local_M[iLocalM];
			  		iCount++;
			  	}
			  	
			  	delete [] Buffer_recv_local_M;
		  	}
    	}
    	
    	/*--- Initialize global_M ---*/
    	global_M = new CSymmetricMatrix;
    	global_M->Initialize(nGlobalVertexDonor, global_M_val_arr);
    }
    
#else
    global_M = new CSymmetricMatrix;
    global_M->Initialize(nVertexDonorInDomain, local_M);
#endif
    
    if (rank == MASTER_NODE) { cout << "\tDone." << endl; }
    
    /*--- Invert M matrix ---*/
    if (rank == MASTER_NODE)
    {
      cout << "Inverting donor zone radial basis function matrix... " << endl;;
    	switch (config[donorZone]->GetKindRadialBasisFunction())
    	{
    	  /*--- Cholesky decompose for basis functions giving positive definite matrix ---*/
				case WENDLAND_C2:
				case INV_MULTI_QUADRIC:
				case GAUSSIAN:

#ifdef HAVE_LAPACK
          cout << "\twith LAPACK... " << endl;
          global_M->CalcInv_dsptri();
#else
					global_M->CholeskyDecompose(true);
#endif

					break;
					
    	  /*--- LU decompose otherwise ---*/
				case THIN_PLATE_SPLINE:
				case MULTI_QUADRIC:
				
#ifdef HAVE_LAPACK
          cout << "\twith LAPACK... " << endl;
          global_M->CalcInv_dsptri();
#else
					global_M->LUDecompose();
#endif

					break;
			}
			
#ifndef HAVE_LAPACK
	    global_M->CalcInv(true);
#endif
      cout << "\tDone." << endl;
    }
    
    calc_polynomial_check = new int [nDim];
    
		/*--- Calculate C_inv_trunc ---*/
		if (rank == MASTER_NODE)
		{
		
		  if ( config[donorZone]->GetRadialBasisFunctionPolynomialOption() ) {
			  
			  /*--- Fill P matrix and get minimum and maximum values ---*/
//			  P_limits = new su2double [nDim*2];
			  P = new su2double [nGlobalVertexDonor*(nDim+1)];
			  iCount = 0;
			  for (iProcessor=MASTER_NODE; iProcessor<nProcessor; iProcessor++)
			  {
				  for (iVertexDonor=0; iVertexDonor<nVertexDonorInDomain_arr[iProcessor]; iVertexDonor++)
				  {
					  P[iCount*(nDim+1)] = 1;
					  for (iDim=0; iDim<nDim; iDim++)
					  {
						  P[iCount*(nDim+1)+iDim+1] = Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_Donor+iVertexDonor)*nDim + iDim];
					
//						  if (iCount == 0)
//						  {
//							  P_limits[iDim*2] = P[iCount*(nDim+1)+iDim+1];
//							  P_limits[iDim*2+1] = P[iCount*(nDim+1)+iDim+1];
//							  calc_polynomial_check[iDim] = 1;
//						  }
//						  else
//						  {
//							  /*--- Get minimum coordinate value ---*/
//							  P_limits[iDim*2] = \
//							  (P[iCount*(nDim+1)+iDim+1] < P_limits[iDim*2]) ? \
//							  P[iCount*(nDim+1)+iDim+1] : P_limits[iDim*2];

//							  /*--- Get maximum coordinate value ---*/
//							  P_limits[iDim*2+1] = \
//							  (P[iCount*(nDim+1)+iDim+1] > P_limits[iDim*2+1]) ? \
//							  P[iCount*(nDim+1)+iDim+1] : P_limits[iDim*2+1];
//						  }
					
					  }
					  iCount++;
				  }
			  }
			  
			  skip_row = new int [nDim+1];
			  skip_row[0] = 1;
			  for (int i=1; i<nDim+1; i++) {
			    skip_row[i] = 0;
			  }
			  cout << "Identifying surfaces lying on plane and discarding one polynomial term..." << endl;
			  Check_PolynomialTerms(nDim+1, nGlobalVertexDonor, P, skip_row, calc_polynomial_check, interface_coord_tol, true, nPolynomial);
			  cout << "\tDone." << endl;
//			  /*--- Check for any constant coordinates which will make RBF matrix singular ---*/
//			  nPolynomial = nDim;
//			  for (iDim=0; iDim<nDim; iDim++)
//			  {
//				  if ( (P_limits[iDim*2+1]-P_limits[iDim*2]) < interface_coord_tol )
//				  {
//					  calc_polynomial_check[iDim] = 0;
//					  nPolynomial--;
//				  }
//			  }
//			
//			  if (nPolynomial < nDim)
//			  {
//				  P_tmp = new su2double [nGlobalVertexDonor*(nPolynomial+1)];
//				
//				  iCount = 0;
//				  for (iVertexDonor=0; iVertexDonor<nGlobalVertexDonor; iVertexDonor++)
//				  {
//					  P_tmp[iCount] = 1;
//					  iCount++;
//					  for (iDim=0; iDim<nDim; iDim++)
//					  {
//						  if (calc_polynomial_check[iDim] == 1)
//						  {
//							  P_tmp[iCount] = P[iVertexDonor*(nDim+1)+iDim+1];
//							  iCount++;
//						  }
//					  }
//				  }
//				
//				  for (iVertexDonor=0; iVertexDonor<nGlobalVertexDonor*(nPolynomial+1); iVertexDonor++) {
//				    P[iVertexDonor] = P_tmp[iVertexDonor];
//				  }
//				  delete [] P_tmp;
//			  }
			
        /*--- Calculate Mp ---*/
        cout << "Calculating Mp..." << endl;
      	Mp = new CSymmetricMatrix;
      	Mp->Initialize(nPolynomial+1);
      	for (int m=0; m<nPolynomial+1; m++)
      	{
      		for (int n=m; n<nPolynomial+1; n++)
      		{
      		
      			val_i = 0;
      			for (iVertexDonor=0; iVertexDonor<nGlobalVertexDonor; iVertexDonor++)
      			{
      			
      				val_j = 0;
      				for (jVertexDonor=0; jVertexDonor<nGlobalVertexDonor; jVertexDonor++)
		    			{
		    				val_j += global_M->Read(iVertexDonor, jVertexDonor)*P[jVertexDonor*(nPolynomial+1)+n];
		    			}
		    			
		    			val_i += val_j*P[iVertexDonor*(nPolynomial+1)+m];
		    			
      			}
      			
      			Mp->Write(m, n, val_i);
      			
      		}
      	}
      	cout << "\tDone." << endl;
      	Mp->CalcInv(true);
      	
      	/*--- Calculate M_p*P*M_inv ---*/
      	cout << "Calculating M_p*P*M_inv..." << endl;
      	C_inv_trunc = new su2double [(nGlobalVertexDonor+nPolynomial+1)*nGlobalVertexDonor];
      	for (int m=0; m<nPolynomial+1; m++)
      	{
      		for (iVertexDonor=0; iVertexDonor<nGlobalVertexDonor; iVertexDonor++)
      		{
      			val_i = 0;
      			for (int n=0; n<nPolynomial+1; n++)
      			{
      				val_j = 0;
      				for (jVertexDonor=0; jVertexDonor<nGlobalVertexDonor; jVertexDonor++)
      				{
      					val_j += P[jVertexDonor*(nPolynomial+1)+n]*global_M->Read(jVertexDonor, iVertexDonor);
      				}
      				val_i += val_j*Mp->Read(m, n);
      			}
      			
      			/*--- Save in row major order ---*/
	      		C_inv_trunc[m*nGlobalVertexDonor+iVertexDonor] = val_i;
      		}  		
      	}
      	cout << "\tDone." << endl;
      	
      	/*--- Calculate (I - P'*M_p*P*M_inv) ---*/
      	cout << "Calculating I - P'*M_p*P*M_inv..." << endl;
      	C_tmp = new su2double [nGlobalVertexDonor*nGlobalVertexDonor];
      	for (iVertexDonor=0; iVertexDonor<nGlobalVertexDonor; iVertexDonor++) 
      	{
      		for (jVertexDonor=0; jVertexDonor<nGlobalVertexDonor; jVertexDonor++)
      		{
      			val_i = 0;
      			for (int m=0; m<nPolynomial+1; m++)
      			{
      				val_i += P[iVertexDonor*(nPolynomial+1)+m]*C_inv_trunc[m*nGlobalVertexDonor+jVertexDonor];
      			}
      		  
            /*--- Save in row major order ---*/
      			C_tmp[iVertexDonor*(nGlobalVertexDonor)+jVertexDonor] = -val_i;

      			
      			if (jVertexDonor==iVertexDonor) { C_tmp[iVertexDonor*(nGlobalVertexDonor)+jVertexDonor] += 1; }
      			
      		}
      	}
      	cout << "\tDone." << endl;
      	
      	/*--- Calculate M_inv*(I - P'*M_p*P*M_inv) ---*/
      	cout << "Calculating M_inv*(I - P'*M_p*P*M_inv)..." << endl;

        global_M->MatMatMult(true, C_tmp, nGlobalVertexDonor);
        
        cout << "\tDone." << endl;
      	
      	/*--- Write to C_inv_trunc matrix ---*/
      	cout << "Writing to C_inv_trunc matrix... ";
      	for (iVertexDonor=0; iVertexDonor<nGlobalVertexDonor; iVertexDonor++)
      	{
      		for (jVertexDonor=0; jVertexDonor<nGlobalVertexDonor; jVertexDonor++)
      		{
      			C_inv_trunc[(iVertexDonor+nPolynomial+1)*(nGlobalVertexDonor)+jVertexDonor] = C_tmp[iVertexDonor*(nGlobalVertexDonor)+jVertexDonor];
      		}
      	}
      	cout << "\tDone." << endl;
    	} // endif RadialBasisFunction_PolynomialOption
    	
    	else {
    	  C_inv_trunc = new su2double [nGlobalVertexDonor*nGlobalVertexDonor];
    	  for (iVertexDonor=0; iVertexDonor<nGlobalVertexDonor; iVertexDonor++)
      	{
      		for (jVertexDonor=0; jVertexDonor<nGlobalVertexDonor; jVertexDonor++)
      		{
      			C_inv_trunc[iVertexDonor*nGlobalVertexDonor+jVertexDonor] = global_M->Read(iVertexDonor, jVertexDonor);
      		}
      	}
    	}
    	
    } // endif (rank == MASTER_NODE)
    
#ifdef HAVE_MPI
	  SU2_MPI::Bcast(&nPolynomial, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
		SU2_MPI::Bcast(calc_polynomial_check, nDim, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);
    
    if (rank != MASTER_NODE)
    {
    	C_inv_trunc = new su2double [(nGlobalVertexDonor+nPolynomial+1)*nGlobalVertexDonor];
    }

  	SU2_MPI::Bcast(C_inv_trunc, (nGlobalVertexDonor+nPolynomial+1)*nGlobalVertexDonor, MPI_DOUBLE, MASTER_NODE, MPI_COMM_WORLD);
#endif
    
    /*--- Calculate H matrix ---*/
    cout << "Rank " << rank << ": Calculating H matrix... ";
    if (config[donorZone]->GetRadialBasisFunctionPolynomialOption()) {
      target_vec = new su2double [nGlobalVertexDonor+nPolynomial+1];
    }
    else {
      target_vec = new su2double [nGlobalVertexDonor];
    }
    
    coeff_vec = new su2double [nGlobalVertexDonor];
    
    for (iVertexTarget = 0; iVertexTarget < nVertexTarget; iVertexTarget++) 
    {
      point_target = target_geometry->vertex[mark_target][iVertexTarget]->GetNode();
    
		  if ( target_geometry->node[point_target]->GetDomain() ) 
		  {
		  
		    iCount = 0;
		    if (config[donorZone]->GetRadialBasisFunctionPolynomialOption())
		    {
		    	target_vec[iCount] = 1;
		    	iCount++;
        }
	    	for (iDim=0; iDim<nDim; iDim++) 
	    	{
	    		Coord_i[iDim] = target_geometry->node[point_target]->GetCoord(iDim);
	    		if (config[donorZone]->GetRadialBasisFunctionPolynomialOption())
	    		{
	      		if (calc_polynomial_check[iDim] == 1)
	      		{
	      			target_vec[iCount] = Coord_i[iDim];
	      			iCount++;
	      		}
	    		}
	    	}
				
		  	for (iProcessor=0; iProcessor<nProcessor; iProcessor++)
		  	{
		  		for (iVertexDonor=0; iVertexDonor<nVertexDonorInDomain_arr[iProcessor]; iVertexDonor++)
		  		{
							for (iDim=0; iDim<nDim; iDim++) { Coord_j[iDim] = Buffer_Receive_Coord[(iProcessor*MaxLocalVertex_Donor+iVertexDonor)*nDim + iDim]; }    		
							Get_Distance(Coord_i, Coord_j, nDim, val_i);
							Get_RadialBasisValue(val_i, config[donorZone]);
							target_vec[iCount] = val_i;
							iCount++;
		  		}
		  	}
		  	
		  	for (iVertexDonor=0; iVertexDonor<nGlobalVertexDonor; iVertexDonor++)
		  	{
		  		coeff_vec[iVertexDonor] = 0;
		  		for (jVertexDonor=0; jVertexDonor<iCount; jVertexDonor++) // May have error here, original loop ends at: jVertexDonor < nGlobalVertexDonor+nPolynomial+1
		  		{ 
		  			coeff_vec[iVertexDonor] += target_vec[jVertexDonor]*C_inv_trunc[jVertexDonor*nGlobalVertexDonor+iVertexDonor];
		  		}
		  	}
				
//				cout << "Rank: " << rank << ", iVertexTarget: " << iVertexTarget << ", coeff_vec: ";
		  	iCount = 0;
		  	for (iVertexDonor=0; iVertexDonor<nGlobalVertexDonor; iVertexDonor++)
		  	{
//		  	  cout << coeff_vec[iVertexDonor] << ", ";
		  		if ( coeff_vec[iVertexDonor] != 0 ) // ( abs(coeff_vec[iVertexDonor]) > interface_coord_tol )
		  		{
		  			iCount++;
		  		}
		  	}
//		  	cout << endl;
//		  	cout << "Rank: " << rank << ", iVertexTarget: " << iVertexTarget << ", nDonorPoints: " << iCount << endl;
		  	target_geometry->vertex[mark_target][iVertexTarget]->SetnDonorPoints(iCount);
		  	target_geometry->vertex[mark_target][iVertexTarget]->Allocate_DonorInfo();
		  	
		  	iCount = 0;
		  	jCount = 0;
		  	for (iProcessor=0; iProcessor<nProcessor; iProcessor++)
		  	{
		  		for (iVertexDonor=0; iVertexDonor<nVertexDonorInDomain_arr[iProcessor]; iVertexDonor++)
		  		{
		  			if ( coeff_vec[iCount] != 0 ) // ( abs(coeff_vec[iCount]) > interface_coord_tol )
						{
							point_donor = Buffer_Receive_GlobalPoint[iProcessor*MaxLocalVertex_Donor+iVertexDonor];
						  target_geometry->vertex[mark_target][iVertexTarget]->SetInterpDonorPoint(jCount, point_donor);
						  target_geometry->vertex[mark_target][iVertexTarget]->SetInterpDonorProcessor(jCount, iProcessor);	
				  		target_geometry->vertex[mark_target][iVertexTarget]->SetDonorCoeff(jCount, coeff_vec[iCount]);
				  		jCount++;
						}
						iCount++;
		  		}
		  	}
		  	
		  }
		}
		cout << "Done." << endl;
    
    /*--- Memory management ---*/
    delete [] nVertexDonorInDomain_arr;
    delete [] nLocalM_arr;
    delete [] local_M;
    delete [] Coord_i;
    delete [] Coord_j;
    delete [] calc_polynomial_check;
    delete [] C_inv_trunc;
    delete [] target_vec;
    delete [] coeff_vec;   
    
    if ( rank == MASTER_NODE )
    {
      delete global_M;
      
      if ( config[donorZone]->GetRadialBasisFunctionPolynomialOption() ) {
        delete [] skip_row;
//        delete [] P_limits;
        delete [] P;
        delete Mp;
        delete [] C_tmp;
      }
      
    }
    
    delete[] Buffer_Send_Coord;
    delete[] Buffer_Send_GlobalPoint;
    
    delete[] Buffer_Receive_Coord;
    delete[] Buffer_Receive_GlobalPoint;

    delete[] Buffer_Send_nVertex_Donor;
    
#ifdef HAVE_MPI
    if (rank == MASTER_NODE)
    {
      delete [] global_M_val_arr;
    }
#endif

  }

  delete[] Buffer_Receive_nVertex_Donor;

#ifdef HAVE_MPI
  if (rank == MASTER_NODE) 
    delete [] Buffer_Recv_mark;
#endif
}

void CRadialBasisFunction::Check_PolynomialTerms(int m, unsigned long n, double *P, int *skip_row, int* keep_row, double max_diff_tol, bool P_transposed, int &n_polynomial)
{
  int n_rows, *write_row;
  unsigned long iCount, jCount;
  double sum, max_diff, max_coeff, *coeff, *P_tmp;
  CSymmetricMatrix *PPT;
  
  n_rows = 0;
  for (int i=0; i<m; i++) {
    if (skip_row[i] == 0) { n_rows++; }
  }
  
  PPT = new CSymmetricMatrix;
  PPT->Initialize(n_rows);

  iCount = 0;
  for (int i = 0; i < m; i ++) {
    if (skip_row[i] == 0) {
    
      jCount = 0;
      for (int j = 0; j < m; j ++){
        if (skip_row[j] == 0) {
        
          sum = 0.0;
          for (unsigned long k = 0; k < n; k ++)
          {
            sum += (P_transposed)? P[k*m+i]*P[k*m+j] : P[i*n+k]*P[j*n+k];
          }
          PPT->Write(iCount, jCount, sum);
          
          jCount++;
        }
      }
      
      iCount++;
    }
  }
  
  PPT->CholeskyDecompose(true);
  PPT->CalcInv(true);

  coeff = new double [n_rows];
  iCount = 0;
  for (int i = 0; i < m; i ++) {
    if (skip_row[i] == 0) {
      coeff[iCount] = 0;
      for (unsigned long j = 0; j < n; j += 1)
      {
        coeff[iCount] += (P_transposed)? P[j*m+i] : P[i*n+j];
      }
      iCount++;
    }
  } 
  
  PPT->MatVecMult(coeff);
  
  max_diff = 0;
  for (unsigned long i = 0; i < n; i ++)
  {
    sum = 0;
    iCount = 0;
    for (int j = 0; j < m; j ++)
    {
      if (skip_row[j] == 0) {
        sum += (P_transposed)? coeff[iCount]*P[i*m+j] : coeff[iCount]*P[j*n+i];
        iCount++;
      }
    }
    max_diff = (abs(1-sum) > max_diff)? abs(1-sum):max_diff;
  }
  
  for (int i=0; i<n_rows; i++) {
    if (max_diff < max_diff_tol) { keep_row[i] = 0; }
    else {keep_row[i] = 1;}
  }
  
  cout << "max_diff: " << max_diff << endl;
  
  /*--- If points lie on plane ---*/
  if (max_diff < max_diff_tol)
  {
    iCount = 0;
    for (int i=0; i<n_rows; i++) {
      
      if (i == 0) {
        iCount = i;
      }
      else {
        iCount = ( abs(coeff[i]) > max_coeff )? i : iCount;
      }
      max_coeff = coeff[iCount];
    }

    for (int i=0; i<n_rows; i++) {
      if ( i != iCount ) {
        keep_row[i] = 1;
      }
    }
    
    n_polynomial = n_rows - 1;
    
    write_row = new int [m];
    iCount = 0;
    jCount = 0;
    for (int i=0; i<m; i++) {
      if (skip_row[i] == 1) {
        write_row[i] = 1;
        jCount++;
      }
      else if (keep_row[iCount] == 1) {
        write_row[i] = 1;
        iCount++;
        jCount++;
      }
      else {iCount++;}
    }
    
    P_tmp = new double [jCount*n];
    iCount = 0;
    if (P_transposed) {
      for (int i=0; i<n; i++) {
        for (int j=0; j<m; j++) {
          if (write_row[j] == 1) {
            P_tmp[iCount] = P[i*m+j];
            iCount++;
          }
        }
      }
    }
    else {
      for (int i=0; i<m; i++) {
        for (int j=0; j<n; j++) {
          if (write_row[i] == 1) {
            P_tmp[iCount] = P[i*n+j];
            iCount++;
          }
        }
      }
    }
    
    for (int i=0; i<jCount*n; i++) {
      P[i] = P_tmp[i];
    }
    
  }
  else {
    n_polynomial = n_rows;
  }
  
  cout << "Coefficients: ";
  for (int i = 0; i < n_rows; i ++)
  {
    cout << coeff[i] << ", ";
  }
  cout << endl;
  
  cout << "keep_row: ";
  for (int i = 0; i < n_rows; i ++)
  {
    cout << keep_row[i] << ", ";
  }
  cout << endl;
  
  delete PPT;
  delete [] coeff;
  if (max_diff<max_diff_tol) {
    delete [] write_row;
    delete [] P_tmp;    
  }
  
}

void CRadialBasisFunction::Get_Distance(su2double *Coord_i, su2double *Coord_j, unsigned short nDim, su2double &dist)
{
	dist = 0;
	for (unsigned short k=0; k<nDim; k++) {
		dist += pow((Coord_i[k] - Coord_j[k]), 2);
	}
	dist = sqrt(dist);
}

void CRadialBasisFunction::Get_RadialBasisValue(su2double &dist, CConfig *config)
{	
	switch (config->GetKindRadialBasisFunction()) {
	
		case WENDLAND_C2:
		  dist /= config->GetRadialBasisFunctionParameter();
		  dist = (dist<1) ? pow((1-dist), 4)*(4*dist+1) : 0;
		break;
			
		case INV_MULTI_QUADRIC:
	    dist = 1/sqrt(dist*dist + \
	    	(config->GetRadialBasisFunctionParameter())* \
	    	(config->GetRadialBasisFunctionParameter()));
		break;
		
		case GAUSSIAN:
		  dist /= config->GetRadialBasisFunctionParameter();
		  dist = exp(-dist*dist);
		break;
			
		case THIN_PLATE_SPLINE:
		  dist /= config->GetRadialBasisFunctionParameter();
		  dist = (dist>0) ? dist*dist*log(dist) : 0;
		break;
		
		case MULTI_QUADRIC:
	    dist = sqrt(dist*dist + \
	    	(config->GetRadialBasisFunctionParameter())* \
	    	(config->GetRadialBasisFunctionParameter()));
		break;
		
	}
}

/*--- Symmetric matrix class definitions ---*/
CSymmetricMatrix::CSymmetricMatrix()
{
	initialized = false;
	inversed = false;
	decomposed = none;
	
	val_vec = NULL;
	decompose_vec = NULL;
	inv_val_vec = NULL;
	perm_vec = NULL;
}

CSymmetricMatrix::~CSymmetricMatrix()
{
	if ( val_vec ) { delete [] val_vec; }
	if ( decompose_vec ) { delete [] decompose_vec; }
	if ( inv_val_vec ) { delete [] inv_val_vec; }
	if ( perm_vec ) { delete [] perm_vec; }
}

void CSymmetricMatrix::Initialize(int N)
{
	int i;
	
	sz = N;	
	num_val = sz*(sz+1)/2;
	val_vec = new double [num_val];
	for (i=0; i<num_val; i++){val_vec[i] = 0.0;}
	
	initialized = true;
}

void CSymmetricMatrix::Initialize(int N, double *formed_val_vec)
{
	int i;
	
	sz = N;	
	num_val = sz*(sz+1)/2;
	
	val_vec = new double [num_val];
	for (unsigned long i=0; i<num_val; i++) {val_vec[i] = formed_val_vec[i];}
	
	initialized = true;
}

inline int CSymmetricMatrix::CalcIdx(int i, int j)
{	
	return max(i, j) + (2*sz-min(i, j)-1)*min(i, j)/2;
}

inline int CSymmetricMatrix::CalcIdxFull(int i, int j)
{
	return i*sz + j;
}

inline int CSymmetricMatrix::GetSize()
{
  return sz;
}

void CSymmetricMatrix::Write(int i, int j, double val)
{
	if (! initialized) {
		throw invalid_argument("Matrix not initialized.");
	}
	else if (i<0 || i>=sz || j<0 || j>=sz) {
		throw out_of_range("Index to write to matrix out of bounds.");
	}
	val_vec[CalcIdx(i, j)] = val;
}

void CSymmetricMatrix::CholeskyDecompose(bool overwrite)
{
	int i, j, k;
	double *vec, sum;
	
	if (! initialized) {
		throw invalid_argument("Matrix not initialized.");
	}
	
	/*--- Point to correct vector ---*/
	if (overwrite) {
		vec = val_vec;
	}
	else {
		decompose_vec = new double [num_val];
		for (i=0; i<num_val; i++){decompose_vec[i] = val_vec[i];}
		vec = decompose_vec;
	}
	
	/*--- Decompose matrix ---*/
	for (j=0; j<sz; j++) {
		for (i=j; i<sz; i++) {
		
			sum = 0.0;
			for (k=0; k<j; k++) {
				if (k<j) {
					sum += vec[CalcIdx(i, k)]*vec[CalcIdx(j, k)];
				}
			}
			
			if (i==j) {
				vec[CalcIdx(i, i)] = sqrt(vec[CalcIdx(i, i)] - sum);
			}
			else {
				vec[CalcIdx(i, j)] = (vec[CalcIdx(i, j)] - sum)/vec[CalcIdx(j, j)];
			}
			
		}
	}
	
	decomposed = cholesky;
	
}

void CSymmetricMatrix::LUDecompose()
{
	bool interchange_row;
	int i, j, k, pivot_idx, tmp_perm_idx;
	double pivot, *tmp_row;
	
	if (! initialized) {
		throw invalid_argument("Matrix not initialized.");
	}
	
	/*--- Copy matrix values to LU matrix ---*/
	decompose_vec = new double [sz*sz];
	perm_vec = new int [sz];
	for (i=0; i<sz; i++) {
		for (j=i; j<sz; j++) {
			decompose_vec[CalcIdxFull(i, j)] = val_vec[CalcIdx(i, j)];
			decompose_vec[CalcIdxFull(j, i)] = decompose_vec[CalcIdxFull(i, j)];
		}
		perm_vec[i] = i;
	}
	
	/*--- Decompose LU matrix ---*/
	tmp_row = new double [sz];
	for (j=0; j<sz-1; j++) {
	
		/*--- Search for pivot and interchange rows ---*/
		interchange_row = false;
		pivot = decompose_vec[CalcIdxFull(j, j)];
		pivot_idx = j;
		for ( i=j+1; i<sz; i++ ) {
			if ( abs(decompose_vec[CalcIdxFull(i, j)]) > abs(pivot) ) {
				pivot = decompose_vec[CalcIdxFull(i, j)];
				pivot_idx = i;
				interchange_row = true;
			}
		}
		
		if ( interchange_row ) {
			
			for ( k=0; k<sz; k++ ) {
				tmp_row[k] = decompose_vec[CalcIdxFull(j, k)];
				decompose_vec[CalcIdxFull(j, k)] = \
					decompose_vec[CalcIdxFull(pivot_idx, k)];
				decompose_vec[CalcIdxFull(pivot_idx, k)] = tmp_row[k];
			}
			
			tmp_perm_idx = perm_vec[j];
			perm_vec[j] = perm_vec[pivot_idx];
			perm_vec[pivot_idx] = tmp_perm_idx;
			
		}
		
		/*--- Perform elimination ---*/
		for ( k=j+1; k<sz; k++ ) {
			decompose_vec[CalcIdxFull(k, j)] /= pivot;
		}
		
		for ( k=j+1; k<sz; k++ ) {
			for (i=j+1; i<sz; i++) {
				decompose_vec[CalcIdxFull(i, k)] -= \
					decompose_vec[CalcIdxFull(j, k)]*decompose_vec[CalcIdxFull(i, j)];
			}
		}
		
	}
	
	delete [] tmp_row;
	
	decomposed = lu;
	
}

void CSymmetricMatrix::CalcInv(bool overwrite)
{
	int i, j, k;
	double *vec, sum, *write_vec;
	
	if ( ! initialized ) {
		throw invalid_argument("Matrix not initialized.");
	}
	
	/*--- Decompose matrix if not already done ---*/
	if ( decomposed == none ) { LUDecompose(); }
	
	/*--- Calculate inverse from decomposed matrices ---*/
	switch ( decomposed ) {
	
		case cholesky:
		
			/*--- Point to correct vector ---*/
			if ( decompose_vec ) { vec = decompose_vec; }
			else { vec = val_vec; }
	
			/*--- Initialize inverse matrix ---*/
			inv_val_vec = new double [num_val];
			for (i=0; i<num_val; i++){inv_val_vec[i] = 0.0;}	
	
			/*---        Calculate L inverse       ---*/
			/*--- Solve smaller and smaller system ---*/
			for (j=0; j<sz; j++) { 
		
				inv_val_vec[CalcIdx(j, j)] = 1.0;
		
				/*--- Forward substitution ---*/
				for (i=j; i<sz; i++) {
		
					if (i==j) {
						inv_val_vec[CalcIdx(i, i)] = 1/ReadL(i, i);
					}
					else {
						sum = 0.0;
						for (k=j; k<i; k++) {
							if (k<i) {
								sum += vec[CalcIdx(i, k)]*inv_val_vec[CalcIdx(k, j)];
							}
						}
						inv_val_vec[CalcIdx(i, j)] = -sum/ReadL(i, i);
					}
			
				}
		
			} // L inverse in inv_val_vec
	
			/*--- Multiply inversed matrices ---*/
			for (j=0; j<sz; j++) {
				for (i=j; i<sz; i++) {
		
					sum = 0.0;
					for (k=i; k<sz; k++) {
						sum += inv_val_vec[CalcIdx(k, i)]*inv_val_vec[CalcIdx(k, j)];
					}
					vec[CalcIdx(i, j)] = sum;
			
				}
			} // Inverse values in vec
	
			/*--- Memory management ---*/
			delete [] inv_val_vec;
			inv_val_vec = NULL;
	
			if (decompose_vec && ! overwrite) {
				inv_val_vec = decompose_vec;
			}
			else if (decompose_vec && overwrite) {
				delete [] val_vec;
				val_vec = decompose_vec;
			}
			
			break;
			
		case lu:
			
			/*--- Point to correct vector ---*/
			vec = decompose_vec;
			
			/*--- Initialize inverse matrix ---*/
			inv_val_vec = new double [sz*sz];
			
			/*--- Invert L and U matrices in place ---*/
			for ( j=0; j<sz; j++ ) {
				inv_val_vec[CalcIdxFull(j, j)] = 1/vec[CalcIdxFull(j, j)];
				for ( i=j+1; i<sz; i++ ) {
				
					inv_val_vec[CalcIdxFull(i, j)] = -vec[CalcIdxFull(i, j)];
					inv_val_vec[CalcIdxFull(j, i)] = \
						-vec[CalcIdxFull(j, i)]*inv_val_vec[CalcIdxFull(j, j)];
					
					if (j+1 <= i) {
						for ( k=j+1; k<i; k++ ) {
						
							inv_val_vec[CalcIdxFull(i, j)] -= \
								vec[CalcIdxFull(i, k)]*inv_val_vec[CalcIdxFull(k, j)];
							
							inv_val_vec[CalcIdxFull(j, i)] -= \
								vec[CalcIdxFull(k, i)]*inv_val_vec[CalcIdxFull(j, k)];
						}
						inv_val_vec[CalcIdxFull(j, i)] /= vec[CalcIdxFull(i, i)];
					}
					
				}
			}
			
			/*--- Multiple U_inv with L_inv ---*/
			for ( i=0; i<sz; i++ ) {
				for ( j=0; j<sz; j++ ) {
					vec[CalcIdxFull(i, j)] = 0.0;
					for ( k=max(i, j); k<sz; k++ ) {
						vec[CalcIdxFull(i, j)] += inv_val_vec[CalcIdxFull(i, k)]* \
						( (k==j) ? 1 : inv_val_vec[CalcIdxFull(k, j)] );
					}
				}
			} // Permuted inverse matrix in vec (equal to decompose_vec, 
			  // which is a full matrix for LU decomposition)
			
			/*--- Get correct vector to write to ---*/
			delete [] inv_val_vec;
			if ( overwrite ) {
				write_vec = val_vec;
				inv_val_vec = NULL;
			}
			else { 
				inv_val_vec = new double [num_val];
				write_vec = inv_val_vec;
			}
			
			/*--- Permutate multiplied matrix to recover A_inv ---*/
			for (j=0; j<sz; j++) {
				k = perm_vec[j];
				for ( i=k; i<sz; i++ ) {
					write_vec[CalcIdx(i, k)] = vec[CalcIdxFull(i, j)];
				}
			}
			write_vec = NULL;

			/*--- Memory management ---*/
			vec = NULL;
			delete [] decompose_vec;
		
			break;
		
	}
	
	decompose_vec = NULL;
	decomposed = none;
	inversed = true;
	
}

void CSymmetricMatrix::CalcInv_dsptri()
{

#ifdef HAVE_LAPACK

  char *uplo;
  int *ipiv, *info;
  double *work;
  
  uplo = new char [1];
  uplo[0] = 'L';
  
  ipiv = new int [sz];
  info = new int [1];
  work = new double [sz];
  
  dsptrf_(uplo, &sz, val_vec, ipiv, info);
  dsptri_(uplo, &sz, val_vec, ipiv, work, info);

  delete [] uplo;
  delete [] ipiv;
  delete [] info;
  delete [] work;
  
  if ( decompose_vec ) { delete [] decompose_vec; }
  decompose_vec = NULL;
	decomposed = none;
	inversed = true;
	
#endif
  
}

void CSymmetricMatrix::CalcInv_dpotri() {}

double CSymmetricMatrix::Read(int i, int j)
{
	if (! initialized) {
		throw invalid_argument("Matrix not initialized.");
	}
	else if (i<0 || i>=sz || j<0 || j>=sz) {
		throw out_of_range("Index to read from matrix out of bounds.");
	}
	return val_vec[CalcIdx(i, j)];
}

double CSymmetricMatrix::ReadL(int i, int j)
{
	double *p;
	
	if (! initialized) {
		throw invalid_argument("Matrix not initialized.");
	}
	else if (decomposed == none) {
		throw invalid_argument("Matrix not decomposed yet or results have been deleted.");
	}
	else if (i<0 || i>=sz || j<0 || j>=sz) {
		throw out_of_range("Index to read from L matrix out of bounds.");
	}
	
	if (decompose_vec){ p = decompose_vec; }
	else {p = val_vec;}
	
	switch (decomposed) {

		case cholesky:
			if (i>=j){ return p[CalcIdx(i, j)];}
			else { return 0.0; }
			break;
			
		case lu:
			if (i>j){ return p[CalcIdxFull(i, j)];}
			else if (i==j) { return 1.0; }
			else { return 0.0; }
			break;
		
	}
	
}

double CSymmetricMatrix::ReadU(int i, int j)
{
	double *p;
	
	if (! initialized) {
		throw invalid_argument("Matrix not initialized.");
	}
	else if (decomposed == none) {
		throw invalid_argument("Matrix not decomposed yet or results have been deleted.");
	}
	else if (i<0 || i>=sz || j<0 || j>=sz) {
		throw out_of_range("Index to read from L matrix out of bounds.");
	}
	
	if (decompose_vec){ p = decompose_vec; }
	else {p = val_vec;}
	
	switch (decomposed) {

		case cholesky:
			return 0.0;
			break;
			
		case lu:
			if (j>=i){ return p[CalcIdxFull(j, i)];}
			else { return 0.0; }
			break;
		
	}
	
}

double CSymmetricMatrix::ReadInv(int i, int j)
{
	double *p;

	if (! initialized) {
		throw invalid_argument("Matrix not initialized.");
	}
	else if (i<0 || i>=sz || j<0 || j>=sz) {
		throw out_of_range("Index to read from matrix out of bounds."); 
	}

	if (inversed) {
		if (inv_val_vec) { p = inv_val_vec; }
		else { p = val_vec; }
		
		return p[CalcIdx(i, j)];
	}
	
	else { throw invalid_argument("Matrix inverse not calculated yet."); }
	
}

void CSymmetricMatrix::VecMatMult(double *v)
{
	double *tmp_res;
	
	tmp_res = new double [sz];
	for (int i=0; i<sz; i++) {
		tmp_res[i] = 0.0;
		for (int k=0; k<sz; k++) {
			tmp_res[i] += v[k]*Read(k, i);
		}
	}
		
	delete [] v;
	v = tmp_res;
	tmp_res = NULL;
	
}

void CSymmetricMatrix::VecMatMult(double *v, int N)
{
	double *tmp_res;
	
	tmp_res = new double [N];
	for (int i=0; i<N; i++) {
		tmp_res[i] = 0.0;
		for (int k=0; k<sz; k++) {
			tmp_res[i] += v[k]*Read(k, i);
		}
	}
	
	delete [] v;
	v = tmp_res;
	tmp_res = NULL;
}

void CSymmetricMatrix::VecMatMult(double *v, double *res, int N)
{	
	for (int i=0; i<N; i++) {
		res[i] = 0.0;
		for (int k=0; k<sz; k++) {
			res[i] += v[k]*Read(k, i);
		}
	}
}

void CSymmetricMatrix::MatVecMult(double *v)
{
	double *tmp_res;
	
	tmp_res = new double [sz];
	for (int i=0; i<sz; i++) {
		tmp_res[i] = 0.0;
		for (int k=0; k<sz; k++) {
			tmp_res[i] += v[k]*Read(i, k);
		}
	}
		
	for (int i=0; i<sz; i++) {
	  v[i] = tmp_res[i];
	}
		
	delete [] tmp_res;

}

void CSymmetricMatrix::MatMatMult(bool left_side, double *mat_vec, int N)
{
	double *tmp_res;
	
	tmp_res = new double [sz*N];
	
	if ( ! left_side ) {
	  throw invalid_argument("Matrix right multiply not implemented yet.");
	}
	
#ifdef HAVE_LAPACK
  
  char side[1]={'R'}, uplo[1]={'L'}; // Right side because mat_vec in row major order
  double *val_full, alpha=1, beta=0;
  
  /*--- Copy packed storage to full storage to use BLAS level 3 routine ---*/
  val_full = new double [sz*sz];
  for (unsigned long i=0; i<sz; i++) {
    for (unsigned long j=i; j<sz; j++) {
      val_full[i+sz*j] = val_vec[CalcIdx(i, j)]; // val_full in column major storage
      if (i != j) {
        val_full[j+sz*i] = val_full[i+sz*j];
      }
    }
  }

  dsymm_(side, uplo, &N, &sz, &alpha, val_full, &sz, mat_vec, &N, &beta, tmp_res, &N);
  for (int i=0; i<sz*N; i++) { mat_vec[i] = tmp_res[i]; }
  
  delete [] val_full;
  
#else

  for (unsigned long i=0; i<sz; i++) {
	  for (unsigned long j=0; j<N; j++) {
		  tmp_res[i*N+j] = 0;
		  for (unsigned long k=0; k<sz; k++) {
			  tmp_res[i*N+j] += val_vec[CalcIdx(i, k)]*mat_vec[k*N+j];
		  }
	  }
  }
	
	for (unsigned long i=0; i<sz*N; i++) {
		mat_vec[i] = tmp_res[i];
	}

#endif

	delete [] tmp_res;
	
}

void CSymmetricMatrix::Print()
{
	cout << "Matrix values:" << endl;
	cout << "["<< endl;
	for (int i=0; i<sz; i++) {
		cout << "[ ";
		for (int j=0; j<sz; j++) {
			cout << Read(i, j) << ", ";
		}
		cout << " ], " << endl;
	}
	cout << "]" << endl;
}

void CSymmetricMatrix::PrintInv()
{
	cout << "Matrix inverse values:" << endl;
	cout << "["<< endl;
	for (int i=0; i<sz; i++) {
		cout << "[ ";
		for (int j=0; j<sz; j++) {
			cout << ReadInv(i, j) << ", ";
		}
		cout << " ], " << endl;
	}
	cout << "]" << endl;
}

void CSymmetricMatrix::PrintLU()
{
	cout << "LU matrix:" << endl;
	cout << "["<< endl;
	for (int i=0; i<sz; i++) {
		cout << "[ ";
		for (int j=0; j<sz; j++) {
			cout << decompose_vec[CalcIdxFull(i, j)] << ", ";
		}
		cout << " ], " << endl;
	}
	cout << "]" << endl;
}

void CSymmetricMatrix::CheckInv()
{
	double sum;
	
	cout << "Multiplying matrix by its inverse:" << endl;
	cout << "["<< endl;
	for (int i=0; i<sz; i++) {
		cout << "[ ";
		for (int j=0; j<sz; j++) {
			sum = 0.0;
			for (int k=0; k<sz; k++) {sum += Read(i, k)*ReadInv(k, j);}
			cout << sum << ", ";
		}
		cout << " ], " << endl;
	}
	cout << "]" << endl;	
}
