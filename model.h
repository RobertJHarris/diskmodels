#ifndef _MODEL_H
#define _MODEL_H

int check_params(double *par);
int sign(int x);
void define_coordinate_transform(double **combined_transform, double inc, double pa);
void compute_stokes(double *image, double *tau,  double *q, double *u, int pix);
void compute_even_odd(double *image, double *even, double *odd, int npix);
void pad_image(double *image, double *padded_image,int pix, int padfactor);
void resort_image(double *array, int npix);
void create_grid(double *x, double *y, double *z, int npix, double inc, double pa);
void initialize_image(double *im, int npix);
void compute_temp_tau(double *x, double *y, double *z, double **combined_transform,double *params, double *temp, double *tau, int npix, double subpixels);
void integrate_radiative_transfer(double *im, double *proj_tau, double *temp, double *tau, int npix);
void populate_b_field(double *bfield, double *x, double *y, double *z,double **combined_transform, int npix);
void interpolate_fft(struct model_fft *fft, struct data_vis *model);
void interpolate_fft_test(struct model_fft *fft, struct model_fft *interpol_fft);
void interpolate_fft_bilin(struct model_fft *fft, struct data_vis *model);
void image_fidelity(struct model_fft *vis, struct model_fft *interp_vis);
double pf(double tau);
double return_zconst_plane(double x, double y, double zval, double inc, double pa);
double compute_log_likelihood(struct data_vis *vis, struct data_vis *model);
double total(double *arr, int npix);
double compute_model_likelihood(struct data_vis *vis, struct walkers *mcmc);
double compute_model_likelihood_params(struct data_vis *vis, double *params);

double compute_external_sum_model(double *params);
#endif
