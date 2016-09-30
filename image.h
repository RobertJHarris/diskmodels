#ifndef _MAKE_VIS_H
#define _MAKE_VIS_H 

struct image{
  double *im;
  double npix;
  double pix_scale;
};


/*
struct model_fft {
  double *real;
  double *imag;
  double *u;
  double *v;
  int udim;
  int vdim;
  double du;
  double dv;

};

struct visibility_datum{
  double real;
  double imag;
  double u;
  double v;
  double wgt;
};

struct data_vis{
  struct visibility_datum *visibility;
  int npts;
};

void read_data(char *filename, struct data_vis *vis);
void compute_fft(struct image *im, int *dims, int dim, struct model_fft *fft);
*/
#endif
