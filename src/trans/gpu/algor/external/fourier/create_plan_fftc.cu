#define cufftSafeCall(err) __cufftSafeCall(err, __FILE__, __LINE__)
#include <cufft.h>
#include <stdio.h>

static const char *_cudaGetErrorEnum(cufftResult error)
{
  switch (error)
    {
#define cr(x) case CUFFT_##x: return #x
      cr (SUCCESS);
      cr (INVALID_PLAN);
      cr (ALLOC_FAILED);
      cr (INVALID_TYPE);
      cr (INVALID_VALUE);
      cr (INTERNAL_ERROR);
      cr (EXEC_FAILED);
      cr (SETUP_FAILED);
      cr (INVALID_SIZE);
      cr (UNALIGNED_DATA);
#undef cr
    }
  return "UNKNOWN";
}

inline void __cufftSafeCall(cufftResult err, const char *file, const int line)
{
  if (CUFFT_SUCCESS != err) 
    {
      fprintf(stderr, "CUFFT at '%s:%d'\n", file, line);
      fprintf(stderr, "CUFFT error %d: %s\nterminating!\n", err, _cudaGetErrorEnum (err));
	  fflush(stderr);
      cudaDeviceReset(); 
    } /*else {
		fprintf(stderr, "CUFFT call at %s, %i returned code %s\n",file,line,_cudaGetErrorEnum(err));
		fflush(stderr);
	}*/

}


static int allocatedWorkspace=0;
static void* planWorkspace;
static int planWorkspaceSize=100*1024*1024; //100MB
 
extern "C"
void
//create_plan_fftc_(cufftHandle *PLANp, int *ISIGNp, int *Np, int *LOTp, int *STRIDEp, int *DISTp)
create_plan_fftc_(cufftHandle *PLANp, int *ISIGNp, int *Np, int *LOTp, int *ISTRIDEp, int *IDISTp, int *OSTRIDEp, int *ODISTp)
{

  int ISIGN = *ISIGNp;
  int N = *Np;
  int LOT = *LOTp;
  
  cufftHandle plan;
  
  if (cudaDeviceSynchronize() != cudaSuccess)
    {
      fprintf(stderr, "%s, %i :Cuda error: Failed to synchronize\n",__FILE__,__LINE__);
      return;	
    }
  
  
  // //create a single re-usable workspace
  // if(!allocatedWorkspace){
  //   allocatedWorkspace=1;
  //   //allocate plan workspace
  //   cudaMalloc(&planWorkspace,planWorkspaceSize);
  // }
  //
  // //disable auto allocation so we can re-use a single workspace (created above)
  //  cufftSetAutoAllocation(plan, false);
  
  int embed[1];
  int istride, ostride;
  int idist, odist;
  
  #ifdef TRANS_SINGLE
  cufftType cufft_1 = CUFFT_R2C;
  cufftType cufft_2 = CUFFT_C2R;
  #else
  cufftType cufft_1 = CUFFT_D2Z;
  cufftType cufft_2 = CUFFT_Z2D;
  #endif
  
  embed[0] = 1;
  istride   = *ISTRIDEp;
  idist     = *IDISTp;
  ostride   = *OSTRIDEp;
  odist     = *ODISTp;
  
  
  cufftSafeCall (cufftCreate (&plan));
  
/*
  if(0){
    fprintf(stderr,"CreatePlan cuFFT\n","N=",N);
    fprintf(stderr,"%s %d \n","plan=",plan);
    fprintf(stderr,"%s %d \n","LOT=",LOT);
    fprintf(stderr,"%s %d \n","ISIGN=",ISIGN);
    fprintf(stderr,"%s %d \n","Np=",*Np);
    fprintf(stderr,"%s %d \n","ISTRIDEp=",*ISTRIDEp);
    fprintf(stderr,"%s %d \n","IDISTp=",*IDISTp);
    fprintf(stderr,"%s %d \n","OSTRIDEp=",*OSTRIDEp);
    fprintf(stderr,"%s %d \n","ODISTp=",*ODISTp);
    fflush (stderr);
  }
*/ 
  
  cufftType type = ISIGN < 0 ? cufft_1 : cufft_2;
  
  cufftSafeCall (cufftPlanMany (&plan, 1, &N, embed, istride, idist, embed, ostride, odist, type, LOT));
  
  if (cudaDeviceSynchronize() != cudaSuccess)
    {
      fprintf(stderr, "%s, %i :Cuda error: Failed to synchronize\n",__FILE__,__LINE__);
	  fflush(stderr);
      return;	
    }
  
  *PLANp=plan;
  
/*
  fprintf(stderr,"cuFFT plan %i created succesfully\n",plan);
  fflush(stderr);
*/  
  // // get size used by this plan
  // size_t workSize;
  // cufftGetSize(plan,&workSize);
  //
  // // exit if we don't have enough space for the work area in the re-usable workspace
  // if(workSize > planWorkspaceSize){
  //   printf("create_plan_fftc: plan workspace size not large enough - exiting\n");
  // exit(1);
  // }
  
  
  return;

}

