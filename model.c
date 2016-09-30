   /* test code for coordinate transformation */
/* priority 1: adaptive integrator for radiative transfer eqn */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Python.h"
#include <time.h>
#include "image.h"
#include "fitsio.h"
#include "constants.h"


int sign(int x) {
  return (x > 0) - (x < 0);
}

 /* helper function to return the z coordinate of the z'=0 (i.e., disk)plane */

void define_coordinate_transform(double **combined_transform, double inc, double pa){
  double **inc_transform;
  double **pa_transform;  /* matrices that hold relvant coordinate transforms*/
  int i,j,k;

  inc_transform =(double **)malloc(sizeof(double*)*3);
  pa_transform =(double **)malloc(sizeof(double*)*3);
  
  for(i = 0;i<3;i++){
      inc_transform[i] = (double*) malloc(sizeof(double)*3); 
      pa_transform[i] = (double*) malloc(sizeof(double)*3); 
      combined_transform[i] = (double*) malloc(sizeof(double)*3); 
      for(j=0;j<3;j++){
	inc_transform[i][j]=0;
	pa_transform[i][j] =0;
	combined_transform[i][j]=0;
      }
  }

  inc_transform[0][0] = (pa_transform[2][2] = 1);
  inc_transform[1][1] = (inc_transform[2][2] = cos(inc*DTOR));
  inc_transform[1][2] = -1*(inc_transform[2][1] = -sin(inc*DTOR));
  pa_transform[0][0] = (pa_transform[1][1] = cos(pa*DTOR));
  pa_transform[0][1] = -1*(pa_transform[1][0] = -sin(pa*DTOR));
  
  for(i = 0; i<3;i++)
    for(j = 0; j<3;j++)
      for(k = 0;k<3;k++)
	combined_transform[i][j] += inc_transform[i][k]*pa_transform[k][j];

  for(i = 0; i<3;i++){
    free(inc_transform[i]);
    free(pa_transform[i]);
  }
  fprintf(stderr,"\n");
  if(DEBUG){
    for(i = 0; i<3;i++){
      for(j = 0; j<3;j++)
	fprintf(stderr,"%lf ",combined_transform[i][j]);
    } 
    fprintf(stderr,"\n");
  }

  free(inc_transform);
  free(pa_transform);
  return;
}

double pf(double tau){
  return (tau  < 1) ? 1 : 0;  /* only allow optically thin pol emission */
}

void compute_stokes(double *image, double *tau,  double *q, double *u, int pix){
  int i = 0;
  double angle =124.0;
  
  for(i = 0; i < pix*pix; i++){
    double pol_i = image[i] * pf(tau[i]);
    q[i]  = pol_i * cos(2*angle*DTOR);
    u[i]  = pol_i * sin(2*angle*DTOR);
  }  
  return;
}

void compute_even_odd(double *image, double *even, double *odd, int npix){
  int i;
  int j;
  for(i = 1; i < npix;i++){
    even[i] = (image[i] + image[npix -i])/2.0;
    odd[i] = (image[i] - image[npix-i])/2.0;
  }    
  for(j = 1; j < npix;j++){
    even[j*npix] = (image[j*npix] + image[npix * (npix - j)])/2.0;
    odd[j*npix] = (image[j*npix] - image[npix * (npix - j)])/2.0;
  }
  for(i = 1; i < npix; i++){
    for(j = 1; j < npix; j++){
      even[i+j*npix] = (image[i  + j *npix] + image[npix * (npix - j) + npix - i])/2.0;
      odd[i+j*npix] = (image[i  + j *npix] - image[npix * (npix - j) + npix - i])/2.0;
    }
  }
  return;
}


double return_zconst_plane(double x, double y, double zval, double inc, double pa){  
  double cpa  = cos(pa*DTOR), spa = sin(pa*DTOR), cinc = cos(inc*DTOR), sinc = sin(inc*DTOR);
  double det  = cpa*cpa*cinc + spa*spa*cinc;
  double yval = 1./det * (( -1*spa * cinc * (x - sinc*spa*zval)) + cpa * (y + sinc*cpa*zval)); 
  return sinc * yval + cinc * zval; 
}

void create_grid(double *x, double *y, double *z, int npix, double inc, double pa){
  int i;
  for(i = 0; i < npix*npix*NZ;i++){
    int zind = (NZ/2 - i /(npix*npix));
    double zmax = 100*AU;
    double zmin = 0.0001*AU;
    double pix_scale = IMSIZE/npix;
    double zval = zind > 0 ? pow(10,(log10(zmax) - log10(zmin))/(NZ/2)*zind + log10(zmin)): 
      -1*pow(10,(log10(zmax) - log10(zmin))/(NZ/2)*(abs(zind)) + log10(zmin));  
    x[i] = ((i % npix) - npix/2)*DIST*pix_scale;  
    // if((x[i] < (4.4*AU)) && (x[i] > (4.3*AU)))
    //  fprintf(stderr,"pix = %d \n",i%npix);
    y[i] = (((i / npix)%npix - npix/2))*DIST*pix_scale;  
    z[i] = return_zconst_plane(x[i],y[i],zval,inc,pa); /*currently only the disk plane */
  }

}


void initialize_image(double *im, int npix){
  int i;
  for(i = 0; i < npix*npix;i++)
    im[i] = 0;
  return;
}

void compute_temp_tau(double *x, double *y, double *z, double **combined_transform, 
		      double *params, double *temp, 
		      double *tau, int npix, double subpixels){
  int i,j,k,l;

  double model_x=0;
  double model_y=0; 
  double model_z=0;
  double model_r;
  double model_phi;
  double h, d_tau, alpha=0;

  double sigma_10   = params[0];
  double rc        = params[1]*AU;
  double gamma     = params[2];
  double q_disk    = params[3];

  double pix_scale = IMSIZE/npix;

  for(i = 0; i < npix*npix*NZ;i++){
    double dz = i < npix*npix*(NZ-1) ? (z[i] - z[i+ npix*npix]) : 0 ; /*stupid kludge for now */ 
    model_x = x[i]*combined_transform[0][0] +
      y[i]*combined_transform[0][1] +
      z[i]*combined_transform[0][2];
    
    model_y = x[i]*combined_transform[1][0] +
      y[i]*combined_transform[1][1] +
      z[i]*combined_transform[1][2];
    
    model_z = x[i]*combined_transform[2][0] +
      y[i]*combined_transform[2][1] +
      z[i]*combined_transform[2][2];
  
    /* convert to r,theta,z */ 
    model_r = sqrt(model_x*model_x + model_y*model_y);
    model_phi= atan2(model_y,model_x) ;
    model_phi = model_phi > 0 ? model_phi : model_phi + PI;
    
    if((model_r > 0 )){ // 2*SUBPIXREGION*AU)){
      temp[i]      = T10*pow(model_r/(10*AU),-1*q_disk);
      h         =  sqrt(KB*temp[i]*AD_GAMMA/MP) * (1/sqrt((GG * MSTAR/(pow(model_r,3)))));
      /* gas to dust  + opacity */
      alpha     = sigma_10*exp(-1.0*pow(model_r/rc,2-gamma))*
	pow((model_r /rc), -1*gamma) *pow((10*AU /rc), 1*gamma) 
	/(sqrt(2*PI)*h)*exp(-model_z*model_z/(2*h*h))
	* exp(1.0*pow(10*AU/rc,2-gamma));
      	
      	
      d_tau  = alpha*dz; 
    }
    else if((model_r == 0)){
      temp[i]=1400.0;
      alpha = 1e100;
      d_tau  = 1;
    }
    else {
      double avg_temp = 0;
      double avg_h    = 0;
      double avg_alpha= 0;
      double local_temp;
      double local_h;
      double local_alpha;
      double subpix = subpixels*subpixels;
      
      for(k = -(subpixels-1)/2; k <= (subpixels -1)/2; k++){
	for(l = -(subpixels-1)/2; l <= (subpixels -1)/2; l++){
	  model_x = (x[i] + k*pix_scale*DIST/subpixels)*combined_transform[0][0] +
	    (y[i] + l*pix_scale*DIST/subpixels)*combined_transform[0][1] +
	    z[i]*combined_transform[0][2];
	  
	  model_y = (x[i] + k*pix_scale*DIST/subpixels)*combined_transform[1][0] +
	    (y[i] + l*pix_scale*DIST/subpixels)*combined_transform[1][1] +
	    z[i]*combined_transform[1][2];
	  
	  model_z = (x[i] + k*pix_scale*DIST/subpixels)*combined_transform[2][0] +
	    (y[i] + l*pix_scale*DIST/subpixels)*combined_transform[2][1] +
	    z[i]*combined_transform[2][2];
	  
	  model_r = sqrt(model_x*model_x + model_y*model_y);	  
	  
	  if(model_r != 0) {
	    local_temp = T10*pow(model_r/(10*AU),-1*q_disk);	  
	    avg_temp += local_temp;
	  }
	  
	  if(model_r > SUBPIXREGION*AU) 
	    continue;
	  if(model_r == 0){
	    subpix -=1;
	    continue;
	  }
	  local_h    = sqrt(KB*local_temp*AD_GAMMA/MP) * (1/sqrt((GG * MSTAR/(pow(model_r,3)))));
	  local_alpha = 0.1 *sigma_10/(sqrt(2*PI)*local_h)*exp(-model_z*model_z/(2*local_h*local_h))*
	    pow(model_r / rc, -1*gamma)*exp(-1.0*pow(model_r/rc,2-gamma)); 
	  
	  
	  if((model_r != 0) &&  (local_temp < DUST_SUB_T))
	    avg_alpha +=  local_alpha;
	  
	}
      }
      
      temp[i] = (subpix != 0) ? avg_temp/subpix : 10;
      alpha = (subpix != 0) ? avg_alpha/subpix : 0;

    }

    if((i % npix == (int)(npix/2))  && ((i /(int)(npix) % npix) == (int)(npix/2)))
      {} //      fprintf(stderr,"dtau = %g\n",d_tau); 
    if(i / (npix*npix) == 0)
      tau[i] = d_tau;
    else
      tau[i] = tau[i - npix*npix] + d_tau;    
  }
  return;
}


double total(double *arr, int npix){
  int i= 0;
  double total = 0;
  for(i=0;i<npix*npix;i++)
    total = total+arr[i];
  return total;
}


void integrate_radiative_transfer(double *im, double *proj_tau, double *temp, double *tau, int npix){
  FILE *fp;
  int i,j;
  double pix_scale = IMSIZE/npix;
  double flux;
  double delta_tau_lim  = 0.1; 
  double delta_tau;
  for(i = npix*npix*NZ-1; i >= npix*npix;i--){ 
    delta_tau = (tau[i] - tau[i - npix*npix]);
    int steps = (int)(delta_tau / delta_tau_lim);
    
    double source_fn = temp[i] > 0 ? 2*PLANCK*NU*NU*NU/(CC*CC) * 
      1.0/(exp(PLANCK*NU/(KB*temp[i]))-1) : 0 ;
    
    if(steps <= 1)
      im[i % (npix*npix)] =im[i % (npix*npix)] 
	+ source_fn*exp(-1*tau[i])*(tau[i] - tau[i-npix*npix]);
    else
      {
	delta_tau = (tau[i] - tau[i-npix*npix])/steps;
	for(j=0;j<steps;j++)     
	  im[i % (npix*npix)] =im[i % (npix*npix)] 
	    + source_fn*exp(-1*(tau[i-npix*npix]+j*delta_tau))*delta_tau;
      }
  
    if(i > npix*npix*NZ - npix*npix){
      proj_tau[i-(npix*npix*NZ - npix*npix)] = tau[i];
    }
  }
  
  

    
  for(i = 0; i < npix*npix;i++)
    im[i] = im[i] * pow(pix_scale*DTOR/3600.0,2)*TOJY;
  
  return;
}


void populate_b_field(double *bfield, double *x, double *y, double *z, 
		      double **combined_transform, int npix){
  
  int i;
  double model_x, model_y, model_z,model_r, model_phi;
  double model_b_x, model_b_y, model_b_z;

  /* compute (again) the disk coordinate system coordinates for our grid */

  for(i = 0; i < npix*npix*NZ;i++){
    model_x = x[i]*combined_transform[0][0] +
      y[i]*combined_transform[0][1] +
      z[i]*combined_transform[0][2];
    
    model_y = x[i]*combined_transform[1][0] +
      y[i]*combined_transform[1][1] +
      z[i]*combined_transform[1][2];
    
    model_z = x[i]*combined_transform[2][0] +
      y[i]*combined_transform[2][1] +
      z[i]*combined_transform[2][2];
  
    /* convert to r,theta,z */ 
    model_r = sqrt(model_x*model_x + model_y*model_y);
    model_phi= atan2(model_y,model_x) ;
    //    model_phi = model_phi > 0 ? model_phi : model_phi + 2*PI;
    
    /* create toroidal bfield */
    model_b_x = -1*sin(model_phi);
    model_b_y = cos(model_phi);
    model_b_z = 0;

    
    bfield[3*i] = model_b_x*combined_transform[0][0] +
      model_b_y*combined_transform[1][0];
    bfield[3*i+1] = model_b_x*combined_transform[0][1] +
      model_b_y*combined_transform[1][1];
    bfield[3*i+2] = model_b_x*combined_transform[0][2] +
      model_b_y*combined_transform[1][2];
      
  }
}

void create_model(double *params){

  /* image / vis size declarations */
  int npix = (int)(NX);
  int size[2] = {(int)npix,(int)npix};
  int dim     = 2;
  int cube[3] = {(int)NZ,(int)npix,(int)npix};
  int cubed   = 3;
  double pix_scale = IMSIZE/npix;
  struct image *modelim;

  /* subgrid for high gradient regions */

  //double subpixels = atof(argv[4]); 
  double subpixels = 1;


  /* model parameterization */
  double disk_param[6] = {params[0],
			  params[1],
			  params[2],
			  params[3],
			  params[4],
			  params[5]};

  //  double env_param[3]  = {0E-18,0.5,1.0};
  double inc = disk_param[4];
  double pa  = disk_param[5];

  /* arrays to hold grids */

  double *xi,*yi,*zi;
  double *temp;       /* temperature */
  double *tau;        /* optical depth CUBE -- optical depth at given slice of cube. */
  double *tauplane;   /* the first plane ot the optical depth cube. tau through source */

  double *b_field;
  double *q_image;
  double *u_image; 

  /* counters */
  int i;
  
  /* coordinate transformation */

  double **combined_transform;

  modelim  = (struct image*) malloc(sizeof(struct image));
  modelim->im = (double *)malloc(sizeof(double)*npix*npix);
  modelim->npix = (int)(npix);
  modelim->pix_scale = pix_scale;

  /*setup image grid*/

  xi = (double *)malloc(sizeof(double)*npix*npix*NZ);
  yi = (double *)malloc(sizeof(double)*npix*npix*NZ);
  zi = (double *)malloc(sizeof(double)*npix*npix*NZ);

  temp = (double*) malloc(sizeof(double)*npix*npix*NZ);
  tau = (double*) malloc(sizeof(double)*npix*npix*NZ);
  tauplane = (double*)malloc(sizeof(double)*npix*npix);
  


  /* B-field */  
  /* this will hold the B-field in the coordinate frame of the observer*/
  
  b_field =  (double*) malloc(sizeof(double)*npix*npix*NZ*3);
  q_image = (double*)malloc(sizeof(double)*npix*npix);
  u_image = (double*)malloc(sizeof(double)*npix*npix);

  combined_transform =(double **)malloc(sizeof(double*)*3);


  /* check to make sure that we have enough memory. probably need to expand this... */

  
  if(modelim->im == NULL || combined_transform == NULL || xi == NULL || yi == NULL || zi==NULL || temp==NULL || tau==NULL ||  tauplane==NULL){
    fprintf(stderr, "Memory allocation failure. Run on a better machine, dummy.\n");
    exit(-1);
  }

  /* initialize and define geometry */
  
  if(DEBUG)
    fprintf(stderr,"Defining coordinate system...");

  define_coordinate_transform(combined_transform,inc,pa);
  create_grid(xi,yi,zi,npix,inc,pa);

  if(DEBUG)
    fprintf(stderr,"done.\n");
  initialize_image(modelim->im, npix);
  if(DEBUG)
    fprintf(stderr,"Computing temperature and optical depth from model...");
  compute_temp_tau(xi, yi, zi, combined_transform, disk_param, temp, tau, npix,subpixels);  
  if(DEBUG)
    fprintf(stderr,"done.\n");

  /* create b_field */
  initialize_python();
  //  make_fits_model_im("zi.fits",zi,cube,cubed,pix_scale);      
  //populate_b_field(b_field,xi,yi,zi,combined_transform, npix);

  /* free memory we don't need anymore */
  free(xi);
  free(yi); 
  free(zi);
  free(b_field); 
  free(u_image); 
  free(q_image);
  for(i = 0; i < 3; i++)
    free(combined_transform[i]);
  free(combined_transform);
  
  if(DEBUG)
    fprintf(stderr,"Integrating radiative transfer eqn...");  
  
  integrate_radiative_transfer(modelim->im,tauplane,temp,tau,npix);
  /*  compute_stokes(modelim->im,tauplane,q_image,u_image,npix);*/
  
  /* free memory we don't need anymore */
  free(temp);
  free(tau);
  free(tauplane);
  if(DEBUG)
       fprintf(stderr,"done.\n");

  /* deposit images and shiz into fits files */


  if(DEBUG)
    fprintf(stderr,"Exporting images to fits...."); 
  
  if(DEBUG){
    fprintf(stderr,"%d  %d\n",size[0],size[1]);
    fprintf(stderr,"%d\n",dim);
    fprintf(stderr,"%lf\n",pix_scale);
  }
  make_fits_model_im("disk.fits",modelim->im,size,dim,pix_scale);   
  //make_fits_model_im("tau.fits",tauplane,size,dim,pix_scale);   
  //make_fits_model_im("q.fits",q_image,size,dim,pix_scale);   
  //make_fits_model_im("u.fits",u_image,size,dim,pix_scale);   
  //make_fits_model_im("optical_depth.fits",tau,cube,cubed,pix_scale);   
  //make_fits_model_im("temp.fits",tau,cube,cubed,pix_scale);     
  
  //  make_fits_model_im("before.fits",modelim->im,size,dim,pix_scale);   
  fprintf(stderr,"The total flux is %lf\n",total(modelim->im,npix));

  /* output x component of bfield */
  //
  //{
  //  double *bx;
  //  bx = (double *)malloc(sizeof(double)*npix*npix*NZ);
  //  for(i = 0; i < npix*npix*NZ; i++)
  //  bx[i] = b_field[3*i];
  //  make_fits_model_im("bx.fits",bx,cube,cubed,pix_scale);   
  //  free(bx);
  //}
  
  free(modelim->im);    
  free(modelim);        
  return;
 
}




