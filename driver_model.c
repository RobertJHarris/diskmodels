/* driver.c 
   rjh, aug 25, 2016
   code that calls a model creation function to generate a fits file
   that corresponds to an image of a disk with some parameterized 
   model 

   to incorporate in detail:
   1) polarization due to magnetic field, (mis)/aligned grains
   2) scattering processes 
   3) subpixellization?
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Python.h"
#include <time.h>

#include "fitsio.h"
#include "constants.h"


int main(int argc, char *argv[]){

  if(argc != 2){
    fprintf(stderr, "Usage %s:  [parameter filename]\n",argv[0]);
    exit(-1);
  }
  /* parameter file  that contains the parameterization*/
  char *paramfile;
  double *params; 
  params = (double *)malloc(sizeof(double)*NPARAMS);

  paramfile = argv[1];
  
  /* check to make sure this is what you wanted ... */
  fprintf(stderr,"paramfile = %s \n",paramfile);
  
  /* create function here to pass the parameters to the function that creates the model */
  get_parameters(params,paramfile);
  create_model(params);          
  return 0;
}

void get_parameters(double *params, char *paramfile){
  FILE *PARAMFILE;
  PARAMFILE = fopen(paramfile,"r");
  int k;

  if(PARAMFILE == NULL){
    fprintf(stderr,"Could not open parameter file. Dying.\n");
    exit(-1);
  }
  for(k=0; k < NPARAMS; k++){
    fscanf(PARAMFILE,"%lf",params+k);
  }
  return;
}
