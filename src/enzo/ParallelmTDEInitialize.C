
/***********************************************************************
/
/  INITIALIZE A COLLAPSE TEST
/
/  written by: Greg Bryan
/  date:       May, 1998
/  modified1:
/
/  PURPOSE:
/    Set up a number of spherical objects
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.
//
#ifdef USE_MPI
#include "mpi.h"
#ifdef USE_MPE
#include "mpe.h"
#endif /* USE_MPE */
#endif /* USE_MPI */

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "fortran.def"
#include "error.h"
#include "message.h"

/* function prototypes */
 
void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void WriteListOfInts(FILE *fptr, int N, int nums[]);
int CommunicationAllSumIntegerValues(int *Values, int Number);

static char *mTDEDensityName               = NULL;
static char *mTDEVelName1                  = NULL;
static char *mTDEVelName2                  = NULL;
static char *mTDEVelName3                  = NULL;
static int   mTDESubgridsAreStatic         = TRUE;
static int   mTDENumberOfInitialGrids      = 1;
static float mTDEBHMass                    = FLOAT_UNDEFINED;
#define MAX_INITIAL_GRIDS 10

int ParallelmTDEInitialize(FILE *fptr, FILE *Outfptr, 
			  HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  const char *DensName = "Density";
  const char *TEName   = "TotalEnergy";
  const char *GEName   = "GasEnergy";
  const char *Vel1Name = "x-velocity";
  const char *Vel2Name = "y-velocity";
  const char *Vel3Name = "z-velocity";
  const char *GravPotName = "GravPotential";

  char *DensityName = NULL, *VelxName = NULL, *VelyName = NULL, *VelzName = NULL;
  /* declarations */

  char line[MAX_LINE_LENGTH];
  int i, j, dim, gridnum, ret, SubgridsAreStatic, region;
  HierarchyEntry *Subgrid;
  gridnum = 0;
 
  /* Set default parameters and names */
 
  int   mTDEGridDimension[MAX_INITIAL_GRIDS][MAX_DIMENSION];
  int   mTDEGridLevel[MAX_INITIAL_GRIDS];
  FLOAT mTDEGridLeftEdge[MAX_INITIAL_GRIDS][MAX_DIMENSION];
  FLOAT mTDEGridRightEdge[MAX_INITIAL_GRIDS][MAX_DIMENSION];
  for (i = 0; i < MAX_INITIAL_GRIDS; i++)
    mTDEGridLevel[i] = 1;
  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
    mTDEGridLeftEdge[0][dim] = DomainLeftEdge[dim];
    mTDEGridRightEdge[0][dim] = DomainRightEdge[dim];
    mTDEGridDimension[0][dim] = MetaData.TopGridDims[dim];
  }

  mTDEGridLevel[0] = 0;

  /* read parameters */
  /* Read input from file. */
 
  char *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    ret = 0;
 
    /* Read parameters */
 
    if (sscanf(line, "mTDEDensityName = %s", dummy) == 1)
      mTDEDensityName = dummy;
    if (sscanf(line, "mTDExVelName = %s", dummy) == 1)
      mTDEVelName1 = dummy;
    if (sscanf(line, "mTDEyVelName = %s", dummy) == 1)
      mTDEVelName2 = dummy;
    if (sscanf(line, "mTDEzVelName = %s", dummy) == 1)
      mTDEVelName3 = dummy;

    ret += sscanf(line, "mTDEBHMass = %"FSYM, &mTDEBHMass);

        /* If the dummy char space was used, then make another. */
 
    if (dummy[0] != 0) {
      dummy = new char[MAX_LINE_LENGTH];
      ret++;
    }
      
    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "mTDE") && line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);
  } 

   /* -------------------------------------------------------------------- */
  /* Generate the root grid and set-up the hierarchy. */
 
  HierarchyEntry *GridsList;
  GridsList = &TopGrid;
 
  /* Initialize the root grid. */

  int TotalRefinement = nint(POW(FLOAT(RefineBy), mTDEGridLevel[gridnum]));

  if (GridsList->GridData->ParallelmTDEInitializeGrid(
			   mTDEDensityName,
			   mTDEVelName1,
			   mTDEVelName2,
			   mTDEVelName3,
               mTDEBHMass,
               mTDESubgridsAreStatic,
			   TotalRefinement) == FAIL) {
      ENZO_FAIL("Error in grid->ParallelmTDEInitializeGrid.\n");
  }

  /* set up field names and units */

  int count = 0;
  DataLabel[count++] = (char*) DensName;
  DataLabel[count++] = (char*) TEName;
  DataLabel[count++] = (char*) Vel1Name;
  DataLabel[count++] = (char*) Vel2Name;
  DataLabel[count++] = (char*) Vel3Name;

  if(WritePotential)
	DataLabel[count++] = (char*) GravPotName;

  for (i = 0; i < count; i++)
    DataUnits[i] = NULL;

    delete dummy;

  return SUCCESS;

}

int ParallelmTDEReInitialize(HierarchyEntry *TopGrid, TopGridData &MetaData)
{
  int dim, gridnum = 0;
   char *DensityName = NULL, *VelxName = NULL,*VelyName = NULL, *VelzName = NULL;

  if (MyProcessorNumber == ROOT_PROCESSOR)
  printf("mTDE: ReInitializing grid %"ISYM"\n", gridnum);

  if (mTDENumberOfInitialGrids > 1) {
    if (mTDEDensityName)
      sprintf(DensityName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
	      mTDEDensityName, gridnum);
    if (mTDEVelName1)
      sprintf(VelxName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
	      mTDEVelName1, gridnum);
    if (mTDEVelName2)
      sprintf(VelyName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
	      mTDEVelName2, gridnum);
   if (mTDEVelName3)
      sprintf(VelzName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
	      mTDEVelName3, gridnum);
   } else {
    DensityName           = mTDEDensityName;
    VelxName              = mTDEVelName1;
    VelyName              = mTDEVelName2;
    VelzName              = mTDEVelName3;
  }

  /* Call grid initializer.  Use TotalRefinement = -1 to flag real read. */
 
  int TotalRefinement = -1;
  HierarchyEntry *Temp = TopGrid;
  while (Temp != NULL) {
 
    if (Temp->GridData->ParallelmTDEInitializeGrid(
			   DensityName,
			   VelxName,
			   VelyName,
			   VelzName,
               mTDEBHMass,
               mTDESubgridsAreStatic,
			   TotalRefinement) == FAIL) {
      ENZO_FAIL("Error in grid->ParallelmTDEInitializeGrid.\n");

    }
 
    Temp = Temp->NextGridThisLevel;
  }
 
  return SUCCESS;
}

