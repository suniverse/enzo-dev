/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A COLLAPSE TEST)
/
/  written by: Xinyu Li
/  date:       May, 2019
/  modified1:
/
/  PURPOSE: set up problem to collapse a FDM halo
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
#ifdef USE_MPI
#include "mpi.h"
#ifdef USE_MPE
#include "mpe.h"
#endif /* USE_MPE */
#endif /* USE_MPI */

#include <hdf5.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "flowdefs.h"
#include "error.h"
#include "CosmologyParameters.h"

void my_exit(int status);
 
#ifdef PROTO /* Remove troublesome HDF PROTO declaration. */
#undef PROTO
#endif

/********************* PROTOTYPES *********************/
int ReadFile(char *name, int Rank, int Dims[], int StartIndex[],
       int EndIndex[], int BufferOffset[], float *buffer,
       inits_type **tempbuffer, int Part, int Npart);
 
int ReadIntFile(char *name, int Rank, int Dims[], int StartIndex[],
    int EndIndex[], int BufferOffset[], int *buffer,
    int **tempbuffer, int Part, int Npart);

void ReadAttribute(hid_t dset_id, int *Attribute, char *AttributeName, FILE *log_fptr, int io_log);
int ReadAttr(char *Fname, int *Rank, int Dims[], int *NSeg, int *LSeg, FILE *log_fptr);

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

/*******************************************************/
#define READFILE ReadFile

int grid::ParallelmTDEInitializeGrid(char *mTDEDensityName, 
                                    char *mTDEVelName1,
                                    char *mTDEVelName2,
                                    char *mTDEVelName3,
                                    float mTDEBHMass,
                                    int mTDESubgridsAreStatic,
                                    int TotalRefinement)
{
  /* declarations */
  int idim, dim, vel, ibx;
  int DeNum;
 
  int ExtraField[2];
 
  inits_type *tempbuffer = NULL;
 
  FILE *log_fptr;
 
#ifdef IO_LOG
  int         io_log = 1;
#else
  int         io_log = 0;
#endif

  char pid[MAX_TASK_TAG_SIZE];
  sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);
 
  char *logname = new char[MAX_NAME_LENGTH];
  strcpy(logname, "TSlog.");
  strcat(logname,pid);
 
  if (io_log) {
    log_fptr = fopen(logname, "a");
    fprintf(log_fptr, "\n");
    fprintf(log_fptr, "TSIG ParallelRootGridIO = %"ISYM"\n", ParallelRootGridIO);
    fprintf(log_fptr, "Processor %"ISYM", Target processor %"ISYM"\n",
        MyProcessorNumber, ProcessorNumber);
    fprintf(log_fptr, "TotalRefinement = %"ISYM"\n", TotalRefinement);
  }

  /* Determine if the data should be loaded in or not. */
 
  int ReadData = TRUE, Offset[] = {0,0,0};
 
  if (ParallelRootGridIO == TRUE && TotalRefinement == 1)
    ReadData = FALSE;
 
  if (io_log) fprintf(log_fptr, "ReadData = %"ISYM"\n", ReadData);

    /* Calculate buffer Offset (same as Grid unless doing ParallelRootGridIO
     (TotalRefinement = -1 if used as a signal that we should really load
     in the data regardless of the value of ParallelRootGridIO). */
 
  if (ParallelRootGridIO == TRUE && TotalRefinement == -1)
    for (dim = 0; dim < GridRank; dim++)
      Offset[dim] = nint((GridLeftEdge[dim] - DomainLeftEdge[dim])/CellWidth[dim][0]);

  int i, j, k, m, field, sphere, size, iden;
  float xdist,ydist,zdist;

  /* create fields */

  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;

  int ivel = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  if (GridRank > 1) 
    FieldType[NumberOfBaryonFields++] = Velocity2;
  if (GridRank > 2)
    FieldType[NumberOfBaryonFields++] = Velocity3;

  //printf("%d \n", NumberOfBaryonFields);
  if(WritePotential)
    FieldType[NumberOfBaryonFields++] = GravPotential;

  /* Set the subgrid static flag. */
 
  SubgridsAreStatic = mTDESubgridsAreStatic;

  if (ProcessorNumber == MyProcessorNumber) {
 
  /* Skip following if NumberOfBaryonFields == 0. */
 
  if (NumberOfBaryonFields > 0) {

      int DensNum = -1, GENum = -1, Vel1Num = -1, 
                         Vel2Num=-1, Vel3Num=-1, TENum=-1,
                         B1Num=-1, B2Num=-1, B3Num=-1, PhiNum=-1;
      IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, 
                           Vel2Num, Vel3Num, TENum,
                           B1Num, B2Num, B3Num, PhiNum);
      
    /* Determine the size of the fields. */
 
    int size = 1;
 
    for (dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];
 
    /* Allocate space for the fields. */
 
    if (ReadData == TRUE) {
      this->AllocateGrids();
    }
    
    if (mTDEDensityName !=NULL && ReadData){
      if (ReadFile(mTDEDensityName, GridRank, GridDimension,
              GridStartIndex, GridEndIndex, Offset, BaryonField[0],
              &tempbuffer, 0, 1) == FAIL) {
        ENZO_FAIL("Error reading density field.\n");
      }
    }

    if (mTDEVelName1 != NULL && ReadData){
      if (ReadFile(mTDEVelName1, GridRank, GridDimension,
              GridStartIndex, GridEndIndex, Offset, BaryonField[Vel1Num],
              &tempbuffer, 0, 1) == FAIL) {
        ENZO_FAIL("Error reading real part of wave function.\n");
      }
    }
  
    if (mTDEVelName2 != NULL && ReadData){
      if (ReadFile(mTDEVelName2, GridRank, GridDimension,
              GridStartIndex, GridEndIndex, Offset, BaryonField[Vel2Num],
              &tempbuffer, 0, 1) == FAIL) {
        ENZO_FAIL("Error reading imaginary part of wave function.\n");
      }
    }
 
    if (mTDEVelName3 != NULL && ReadData){
      if (ReadFile(mTDEVelName3, GridRank, GridDimension,
              GridStartIndex, GridEndIndex, Offset, BaryonField[Vel2Num],
              &tempbuffer, 0, 1) == FAIL) {
        ENZO_FAIL("Error reading imaginary part of wave function.\n");
      }
    }
  } // end: if (NumberOfBaryonFields > 0)
  } // end: if (ProcessorNumber == MyProcessorNumber)
  OldTime = Time;
 
  if (io_log) fclose(log_fptr);

  /* Set various units. */
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, 
    VelocityUnits, CriticalDensity = 1, BoxLength = 1, mu = 0.6;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, 
	   &VelocityUnits, Time);

  // If use particle, initial particles according to the FDM values and turn off QuantumPressure
  int CollapseTestParticleCount = 0;
  int SetupLoopCount, npart = 0;
  int ParticleCount = 0;
  int ind, indxp, indxn, indyp, indyn, indzp, indzn;
  int ip,in,jp,jn,kp,kn;
  double x,y,z,vx,vy,vz;

  if (MyProcessorNumber == ROOT_PROCESSOR){
    fprintf(stdout, "initialize BH on processor %d \n", MyProcessorNumber);
    for (SetupLoopCount = 0; SetupLoopCount < 1; SetupLoopCount++) {
     if (SetupLoopCount > 0) {
      /* If particles already exist (coarse particles), then delete. */
        if (NumberOfParticles > 0) this->DeleteParticles();
      /* Use count from previous loop to set particle number. */
        NumberOfParticles = npart;
        npart = 0;
      /* Allocate space. */
        this->AllocateNewParticles(NumberOfParticles);
      /* Particle values will be set below. */
      }
	// set some test particles
	ParticleCount = 1;

	// set many particles
	while (ParticleCount > 0) {
        if (SetupLoopCount > 0) {
		    ParticleMass[npart] = mTDEBHMass*2e33/(DensityUnits*POW(LengthUnits,3));
	        ParticleNumber[npart] = CollapseTestParticleCount++;
            ParticleType[npart] = PARTICLE_TYPE_BLACK_HOLE;
         // Set random position within cell.
		    ParticlePosition[0][npart] = 0.5;
		    ParticlePosition[1][npart] = 0.5;
		    ParticlePosition[2][npart] = 0.5;
			ParticleVelocity[0][npart] = 0.0;
			ParticleVelocity[1][npart] = 0.0;
			ParticleVelocity[2][npart] = 0.0;
		}
		ParticleCount--;
		npart++;

        }// end while
   } // end loop SetupLoopCount
   NumberOfParticles = npart;
   printf("mTDEInitialize: Number of Particles = %d on Processor %d\n", NumberOfParticles, MyProcessorNumber);


  }
  return SUCCESS;
}
