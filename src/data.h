#ifndef _DATA_H
#define _DATA_H

#include <src/common.h>

typedef struct {
  char *name;
  char *interpolationvarsmap; 
  bool READDATAFILE,BINARYMPI,ASCII,INTERPOLATION, INTERPOLATIONMAP;
  #ifdef UNSTEADY
    char *name_m1;
    bool M1,MM1;
    #if _RESTIME_BW > 2
      char *name_m2;
      bool M2,MM2;
    #endif
    #if _RESTIME_BW > 3
      char *name_m3;
      bool M3,MM3;
    #endif
  #endif
} input_t;


long _ai_mpidatafile(gl_t *gl, long i, long j, long k);

void read_data_file(input_t input, np_t *np, gl_t *gl);

void write_data_file(np_t *np, gl_t *gl);

void read_data_file_binary(char *filename, np_t *np, gl_t *gl, long level);

void read_data_file_mpi(char *filename, np_t *np, gl_t *gl, long level);

void write_data_file_binary(char *filename, np_t *np, gl_t *gl);

void write_data_file_mpi(char *filename, np_t *np, gl_t *gl);

void read_data_file_ascii(char *filename, np_t *np, gl_t *gl, long level);

void write_data_file_ascii(char *filename, np_t *np, gl_t *gl);

void read_data_file_interpolation(input_t input, np_t *np, gl_t *gl);

void write_data_file_interpolation(char *filename, np_t *np, gl_t *gl);


#endif /* _DATA_H */



