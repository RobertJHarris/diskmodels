#ifndef _FITSIO_H
#define _FITSIO_H 


void initialize_python();
void quit_python();
void init_numpy();
void  *set_path();
void make_fits_model_im(char *filename,double *im, int *size, int dim, double pix_scale);
long  return_n_elements_from_dict(PyObject* presult);
int return_python_array_from_dict(PyObject* dict, char *key, double *arr);
//int read_vis_file(char *filename, struct data_vis *vis);
//int read_vis_ascii_file(char *filename, struct data_vis *vis);

#endif
