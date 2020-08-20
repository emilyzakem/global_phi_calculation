#include "fintrf.h"
C======================================================================
C     CADVECT calls advection schemes to calculate scalar transports
C     [wflux,itertrac]=mit_difimpl_flux(timestep,rac,
C                               dzf,dzc,
C                               mask,
C                               diffkz,
C                               scalar);
C
C======================================================================
C     Gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

C     Declarations
      implicit none

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer*4 nlhs, nrhs

C     Function declarations:
      mwPointer mxCreateNumericArray
      integer*4 mxClassIDFromClassName
      integer mxIsDouble, mxIsNumeric
      mwPointer mxGetM, mxGetN
      mwSize mxGetNumberOfDimensions
      mwPointer mxGetDimensions
      mwPointer mxGetPr

C     Arguments for mxCreateNumericArray
      integer*4 classid, complexflag
      mwSize ndim, dims(4)
    
C     Pointers to input/output mxArrays:
      mwPointer deltat
      mwPointer dzf, dzc, rac
      mwPointer mask
      mwPointer diffkz
      mwPointer scalar 
      mwPointer wflux
      mwPointer itertrac

C     Array information:
      mwPointer mrows, ncols
      mwSize size

CC    For writing error to the screen
C     mexErrMsgIdAndTxt causes crash
      character*250 line
      integer*4 k
      integer*4 mexPrintf
      integer*4 errCode
      
C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 7) then
C         call mexErrMsgIdAndTxt ('MATLAB:cadvect:nInput     ',
C     &                           'Fourteen inputs required. ')
         write(line,*)'ERROR: Fourteen inputs required.'
         k=mexPrintf(line//achar(13)) 
         return
      elseif(nlhs .gt. 2) then
C         call mexErrMsgIdAndTxt ('MATLAB:cadvect:nOutput    ',
C     &                           'Too many output arguments.')
         write(line,*)'ERROR: Too many output arguments.'
         k=mexPrintf(line//achar(13)) 
         return
      endif

C     Validate inputs
C     Check that the input is a number.
      if(mxIsNumeric(prhs(1)) .eq. 0) then
C         call mexErrMsgIdAndTxt ('MATLAB:cadvect:NonNumeric',
C     &                           'Input must be a number.  ')
         write(line,*)'ERROR: Input must be numeric.'
         k=mexPrintf(line//achar(13)) 
         return     
      endif
      
C     Check data type of input argument.
      if (mxIsDouble(prhs(1)) .eq. 0) then
C         call mexErrMsgTxt('Input argument must be of type double.')
         write(line,*)'ERROR: Input must be double precision.'
         k=mexPrintf(line//achar(13)) 
         return     
      endif
      
C     Get the size of the input arrays using scalar as template.
C     Use mxGetDimensions to get dimensions into dims
      ndim = mxGetNumberOfDimensions(prhs(nrhs))
      call mxCopyPtrToInteger4(mxGetDimensions(prhs(nrhs)),  
     &           dims, mxGetNumberOfDimensions(prhs(nrhs)))

C      write(line,*)'ndims: ',ndim
C      k=mexPrintf(line//achar(13)) 
C      write(line,*)'dimensions: ',dims(1),'x',dims(2),
C     &                        'x',dims(3),'x',dims(4)
C      k=mexPrintf(line//achar(13)) 

C     Create Fortran arrays from the input arguments.
      deltat = mxGetPr(prhs(1))
      rac    = mxGetPr(prhs(2))
      dzf    = mxGetPr(prhs(3))
      dzc    = mxGetPr(prhs(4))
      mask   = mxGetPr(prhs(5))
      diffkz = mxGetPr(prhs(6))
      scalar = mxGetPr(prhs(7))
      
C     Output will be wflux
C     Create matrix for the return argument.
      classid = mxClassIDFromClassName('double')
      complexflag = 0

      plhs(1) = mxCreateNumericArray(ndim, dims, 
     &          classid, 
     &          complexflag)
     
      plhs(2) = mxCreateNumericArray(ndim, dims, 
     &          classid, 
     &          complexflag)


      wflux = mxGetPr(plhs(1))
      itertrac = mxGetPr(plhs(2))

CC     Call the computational subroutines.
         call mit_implicit_rmix(dims(1),dims(2),dims(3),dims(4),
     &           %val(deltat), %val(mask), %val(rac),
     &           %val(dzf), %val(dzc),%val(diffkz),%val(scalar),
     &           %val(wflux),%val(itertrac),errCode)
     
CC     Check data type of input argument.
      if (errCode .ne. 0) then
CC         call mexErrMsgTxt('Input argument must be of type double.')
         write(line,*)'ERROR:SOLVE_TRIDIAGONAL:Single matrix.'
         k=mexPrintf(line//achar(13)) 
         write(line,*)'errCode=',errCode
         k=mexPrintf(line//achar(13)) 
         return     
      endif
      
      return
      end