#ifndef KD_HINCLUDED
#define KD_HINCLUDED

#include <stdio.h>

#define DARK	1
#define GAS		2
#define STAR	4

#define SOFT_UNI	1
#define SOFT_PLUM	2
#define SOFT_SPLINE	3

typedef struct Particle {
	int iOrder;
	int iMark;
	float fMass;
	float fSoft;
	float r[3];
	double a[3];
	double dPot;
	} PARTICLE;

typedef struct kdContext {
	int nParticles;
	int nGas;
	int nDark;
	int nStar;
	int bGas;
	int bDark;
	int bStar;
	int nActive;
	int nMark;
	float fTime;
	double G;
	float fPeriod[3];
	float fCenter[3];
	int iChkptInterval;
	char *pszChkptName;
	int iSoftType;
	int bPeriodic;
	int iBlockSize;
	int uSecond;
	int uMicro;
	PARTICLE *p;
	} * KD;


struct chkptHeader {
	int nParticles;
	int nGas;
	int nDark;
	int nStar;
	int bGas;
	int bDark;
	int bStar;
	int nActive;
	int nMark;
	double G;
	float fPeriod[3];
	float fCenter[3];
	int iSoftType;
	int bPeriodic;
	int iBlockSize;
	int uSecond;
	int uMicro;
	int lBlock;
	int kBlock;
	};


void kdTime(KD,int *,int *);
int kdInit(KD *,double,float *,float *,int,char *);
int kdReadTipsy(KD,FILE *,int,int,int);
void kdSetSoft(KD,float);
void kdInMark(KD,char *);
void kdMarkOrder(KD);
void kdOrder(KD);
void kdGravSimple(KD,int,int);
void kdGrav(KD,int,int,int,int);
void kdRestart(KD,int);
void kdOutAccel(KD,char *);
void kdOutPot(KD,char *);
void kdFinish(KD);

#endif











