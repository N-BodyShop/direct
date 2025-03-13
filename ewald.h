#ifndef EWALD_HINCLUDED
#define EWALD_HINCLUDED

#include "kd.h"

void diaEwald(PARTICLE *,int,int,double);
void blkEwald(PARTICLE *,int,PARTICLE *,int,int,double);
void umkEwald(PARTICLE *,int,PARTICLE *,int,int,double);

void diaEwald2D(PARTICLE *,int,int,float *);
void blkEwald2D(PARTICLE *,int,PARTICLE *,int,int,float *);
void umkEwald2D(PARTICLE *,int,PARTICLE *,int,int,float *);
#endif
