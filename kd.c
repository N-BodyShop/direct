#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <assert.h>
#include "kd.h"
#include "grav.h"
#include "ewald.h"
#include "tipsydefs.h"


void kdTime(KD kd,int *puSecond,int *puMicro)
{
	struct rusage ru;

	getrusage(0,&ru);
	*puMicro = ru.ru_utime.tv_usec - kd->uMicro;
	*puSecond = ru.ru_utime.tv_sec - kd->uSecond;
	if (*puMicro < 0) {
		*puMicro += 1000000;
		*puSecond -= 1;
		}
	kd->uSecond = ru.ru_utime.tv_sec;
	kd->uMicro = ru.ru_utime.tv_usec;
	}


int kdInit(KD *pkd,double G,float *fPeriod,float *fCenter,
		   int iChkptInterval,char *pszChkptName)
{
	KD kd;
	int j;

	kd = (KD)malloc(sizeof(struct kdContext));
	assert(kd != NULL);
	kd->G = G;
	kd->pszChkptName = (char *)malloc(strlen(pszChkptName)+1);
	assert(kd->pszChkptName != NULL);
	strcpy(kd->pszChkptName,pszChkptName);
	kd->iChkptInterval = iChkptInterval;
	for (j=0;j<3;++j) {
		kd->fPeriod[j] = fPeriod[j];
		kd->fCenter[j] = fCenter[j];
		}
	*pkd = kd;
	return(1);
	}


void kdFinish(KD kd)
{
	free(kd->pszChkptName);
	free(kd->p);
	free(kd);
	}


int kdReadTipsy(KD kd,FILE *fp,int bGas,int bDark,int bStar)
{
	int i,j,nCnt;
	struct dump h;
	struct gas_particle gp;
	struct dark_particle dp;
	struct star_particle sp;

	fread(&h,sizeof(struct dump),1,fp);
	kd->nParticles = h.nbodies;
	kd->nDark = h.ndark;
	kd->nGas = h.nsph;
	kd->nStar = h.nstar;
	kd->fTime = h.time;
	kd->nActive = 0;
	if (bDark) kd->nActive += kd->nDark;
	if (bGas) kd->nActive += kd->nGas;
	if (bStar) kd->nActive += kd->nStar;
	kd->nMark = kd->nActive;
	kd->bDark = bDark;
	kd->bGas = bGas;
	kd->bStar = bStar;
	/*
	 ** Allocate particles.
	 */
	kd->p = (PARTICLE *)malloc(kd->nActive*sizeof(PARTICLE));
	assert(kd->p != NULL);
	/*
	 ** Read Stuff!
	 */
	nCnt = 0;
	for (i=0;i<h.nsph;++i) {
		fread(&gp,sizeof(struct gas_particle),1,fp);
		if (bGas) {
			kd->p[nCnt].fMass = gp.mass;
			kd->p[nCnt].fSoft = gp.hsmooth;
			kd->p[nCnt].iOrder = nCnt;
			kd->p[nCnt].iMark = 1;
			for (j=0;j<3;++j) kd->p[nCnt].r[j] = gp.pos[j];
			++nCnt;
			}
		}
	for (i=0;i<h.ndark;++i) {
		fread(&dp,sizeof(struct dark_particle),1,fp);
		if (bDark) {
			kd->p[nCnt].fMass = dp.mass;
			kd->p[nCnt].fSoft = dp.eps;
			kd->p[nCnt].iOrder = nCnt;
			kd->p[nCnt].iMark = 1;
			for (j=0;j<3;++j) kd->p[nCnt].r[j] = dp.pos[j];
			++nCnt;
			}
		}
	for (i=0;i<h.nstar;++i) {
		fread(&sp,sizeof(struct star_particle),1,fp);
		if (bStar) {
			kd->p[nCnt].fMass = sp.mass;
			kd->p[nCnt].fSoft = sp.eps;
			kd->p[nCnt].iOrder = nCnt;
			kd->p[nCnt].iMark = 1;
			for (j=0;j<3;++j) kd->p[nCnt].r[j] = sp.pos[j];
			++nCnt;
			}
		}
	return(kd->nParticles);
	}


void kdSetSoft(KD kd,float fSoft)
{
	int i;
	
	for (i=0;i<kd->nActive;++i) {
		kd->p[i].fSoft = fSoft;
		}
	}


void kdInMark(KD kd,char *pszFile)
{
	FILE *fp;
	char ach[80];
	int i,iCnt,iDum;

	fp = fopen(pszFile,"r");
	if (!fp) {
		fprintf(stderr,"Could not open mark array, %s\n",pszFile);
		exit(1);
		}
	fgets(ach,80,fp);	/* ignore the array header! */
	iCnt = 0;
	for (i=0;i<kd->nGas;++i) {
		if (kd->bGas) {
			fscanf(fp,"%d",&kd->p[iCnt++].iMark);
			}
		else fscanf(fp,"%d",&iDum);
		}
	for (i=0;i<kd->nDark;++i) {
		if (kd->bDark) {
			fscanf(fp,"%d",&kd->p[iCnt++].iMark);
			}
		else fscanf(fp,"%d",&iDum);
		}
	for (i=0;i<kd->nStar;++i) {
		if (kd->bStar) {
			fscanf(fp,"%d",&kd->p[iCnt++].iMark);
			}
		else fscanf(fp,"%d",&iDum);
		}
	fclose(fp);
	}


void kdMarkOrder(KD kd)
{
	PARTICLE *p,t;
	int i,j;

	p = kd->p;
	i = 0;
	j = kd->nActive-1;
	while (i < j) {
		while (p[i].iMark) if (++i > j) goto done;
		while (!p[j].iMark) if (i > --j) goto done;
		t = p[i];
		p[i] = p[j];
		p[j] = t; 
		}
 done:
	kd->nMark = i;
	}


int cmpParticles(const void *v1,const void *v2)
{
	PARTICLE *p1=(PARTICLE *)v1,*p2=(PARTICLE *)v2;
	
	return(p1->iOrder - p2->iOrder);
	}


void kdOrder(KD kd)
{
	qsort(kd->p,kd->nActive,sizeof(PARTICLE),cmpParticles);
	}


/*
 ** A simple driver. No support for check-pointing here!
 */
void kdGravSimple(KD kd,int iSoftType,int bPeriodic)
{
	PARTICLE *p;
	int i,n;

	p = kd->p;
	n = kd->nActive;
	for (i=0;i<n;++i) {
		p[i].a[0] = 0.0;
		p[i].a[1] = 0.0;
		p[i].a[2] = 0.0;
		p[i].dPot = 0.0;
		}
	/*
	 ** Process one huge diagonal Block!
	 */
	if (bPeriodic) {
		diaEwald(p,n,iSoftType,kd->fPeriod[0]);
		}
	else {
		diaGrav(p,n,iSoftType);
		}
	for (i=0;i<n;++i) {
		p[i].a[0] *= kd->G;
		p[i].a[1] *= kd->G;
		p[i].a[2] *= kd->G;
		p[i].dPot *= kd->G;
		}
	}


void WriteChkpt(KD kd,int lBlock,int kBlock)
{
 	FILE *fp;
	struct rusage ru;
	struct chkptHeader h;
	int j;
	
	h.nParticles = kd->nParticles;
	h.nGas = kd->nGas;
	h.nDark = kd->nDark;
	h.nStar = kd->nStar;
	h.bGas = kd->bGas;
	h.bDark = kd->bDark;
	h.bStar = kd->bStar;
	h.nActive = kd->nActive;
	h.nMark = kd->nMark;
	h.G = kd->G;
	for (j=0;j<3;++j) {
		h.fPeriod[j] = kd->fPeriod[j];
		h.fCenter[j] = kd->fCenter[j];
		}
	h.iSoftType = kd->iSoftType;
	h.bPeriodic = kd->bPeriodic;
	h.iBlockSize = kd->iBlockSize;
	/*
	 ** Store time used up to now.
	 */
	getrusage(0,&ru);
	h.uMicro = ru.ru_utime.tv_usec - kd->uMicro;
	h.uSecond = ru.ru_utime.tv_sec - kd->uSecond;
	if (h.uMicro < 0) {
		h.uMicro += 1000000;
		h.uSecond -= 1;
		}
	h.lBlock = lBlock;
	h.kBlock = kBlock;
	fp = fopen(kd->pszChkptName,"wb");
	assert(fp != NULL);
	fwrite(&h,sizeof(struct chkptHeader),1,fp);
	fwrite(kd->p,sizeof(PARTICLE),kd->nActive,fp);
	fclose(fp);
	}


void ReadChkpt(KD kd,int *plBlock,int *pkBlock)
{
	FILE *fp;
	struct rusage ru;
	struct chkptHeader h;
	int j;

	fp = fopen(kd->pszChkptName,"rb");
	if (fp == NULL) {
		fprintf(stderr,"Sorry the check point file %s could not be opened\n",
				kd->pszChkptName);
		kdFinish(kd);
		exit(1);
		}
	fread(&h,sizeof(struct chkptHeader),1,fp);
	kd->nParticles = h.nParticles;
	kd->nGas = h.nGas;
	kd->nDark = h.nDark;
	kd->nStar = h.nStar;
	kd->bGas = h.bGas;
	kd->bDark = h.bDark;
	kd->bStar = h.bStar;
	kd->nActive = h.nActive;
	kd->nMark = h.nMark;
	kd->G = h.G;
	for (j=0;j<3;++j) {
		kd->fPeriod[j] = h.fPeriod[j];
		kd->fCenter[j] = h.fCenter[j];
		}
	kd->iSoftType = h.iSoftType;
	kd->bPeriodic = h.bPeriodic;
	kd->iBlockSize = h.iBlockSize;
	/*
	 ** Fix up CPU time.
	 */
	getrusage(0,&ru);
	kd->uMicro = ru.ru_utime.tv_usec - h.uMicro;
	kd->uSecond = ru.ru_utime.tv_sec - h.uSecond;
	if (kd->uMicro < 0) {
		kd->uMicro += 1000000;
		kd->uSecond -= 1;
		}
	/*
	 ** Allocate particles.
	 */
	kd->p = (PARTICLE *)malloc(kd->nActive*sizeof(PARTICLE));
	assert(kd->p != NULL);
	fread(kd->p,sizeof(PARTICLE),kd->nActive,fp);
	*plBlock = h.lBlock;
	*pkBlock = h.kBlock;
	fclose(fp);
	}


/*
 ** Driver routine for direct. With support for check points!
 */
void Grav(KD kd,int lBlock,int kBlock,int iSoftType,int bPeriodic,
		  int bVerbose)
{
	PARTICLE *p,*q;
	int k,l,n,nk,nb,nbk,bs,rs,rsk,iStride,bc,i;

	p = kd->p;
	n = kd->nMark;
	q = &kd->p[n];
	nk = kd->nActive - kd->nMark;
	bs = kd->iBlockSize;
	nb = n/bs;
	rs = n%bs;
	nbk = nk/bs;
	rsk = nk%bs;
	iStride = 16384/bs;
	if (bPeriodic) iStride = iStride/100;
	if (!iStride) iStride = 1;
	if (lBlock || kBlock) {
		l = lBlock;
		k = kBlock;
		bc = 1;
		if (k < nb) goto Restart1;
		else {
			k -= nb;
			goto Restart2;
			}
		}
	for (i=0;i<n;++i) {
		p[i].a[0] = 0.0;
		p[i].a[1] = 0.0;
		p[i].a[2] = 0.0;
		p[i].dPot = 0.0;
		}
	/*
	 ** First do all the diagonal blocks.
	 */
	bc = 1;
	for (k=0;k<nb;++k,++bc) {
		if (bVerbose) {
			if (!(bc%iStride)) {
				printf("Block:(%d,%d)\n",k,k);
				fflush(stdout);
				bc = 0;
				}
			}
		if (bPeriodic) {
			diaEwald(&p[k*bs],bs,iSoftType,kd->fPeriod[0]);
			}
		else {
			diaGrav(&p[k*bs],bs,iSoftType);
			}
		}
	if (bVerbose && rs) {
		printf("Block:(%d,%d)\n",k,k);
		fflush(stdout);
		}
	if (bPeriodic) {
		diaEwald(&p[k*bs],rs,iSoftType,kd->fPeriod[0]);
		}
	else {
		diaGrav(&p[k*bs],rs,iSoftType);
		}
	/*
	 ** Now do the off-diagonal blocks.
	 */
	bc = 1;
	for (l=1;l<nb;++l) {
		for (k=l;k<nb;++k,++bc) {
			if (bc == kd->iChkptInterval) {
				WriteChkpt(kd,l,k);
				if (bVerbose) {
					printf("Check point (%d,%d) written\n",l-1,k);
					fflush(stdout);
					}
				bc = 1;
				}
		Restart1:
			if (bVerbose) {
				if (!(bc%iStride)) {
					printf("Block:(%d,%d)\n",l-1,k);
					fflush(stdout);
					}
				}
			if (bPeriodic) {
				blkEwald(&p[(l-1)*bs],bs,&p[k*bs],bs,iSoftType,kd->fPeriod[0]);
				}
			else {
				blkGrav(&p[(l-1)*bs],bs,&p[k*bs],bs,iSoftType);
				}
			}
		}
	bc = 1;
	for (l=0;l<nb;++l,++bc) {
		if (!rs) break;
		if (bPeriodic) {
			blkEwald(&p[l*bs],bs,&p[nb*bs],rs,iSoftType,kd->fPeriod[0]);
			}	
		else {
			blkGrav(&p[l*bs],bs,&p[nb*bs],rs,iSoftType);
			}
		}
	/*
	 ** We are now done the mutual interactions now we need to do the
	 ** interactions due to the unmarked particles on the marked ones.
	 */
	bc = 1;
	if (nk) {
		if (bVerbose) {
			printf("Now doing unmarked particle interactions\n");
			fflush(stdout);
			}
		for (l=0;l<nb;++l) {
			for (k=0;k<nbk;++k,++bc) {
				if (bc == kd->iChkptInterval) {
					WriteChkpt(kd,l,k+nb);
					if (bVerbose) {
						printf("Check point (%d,%d) written\n",l,k+nb);
						fflush(stdout);
						}
					bc = 1;
					}
			Restart2:
				if (bVerbose) {
					if (!(bc%iStride)) {
						printf("Block:(%d,%d)\n",l,k+nb);
						fflush(stdout);
						}
					}
				if (bPeriodic) {
					umkEwald(&p[l*bs],bs,&q[k*bs],bs,iSoftType,kd->fPeriod[0]);
					}
				else {
					umkGrav(&p[l*bs],bs,&q[k*bs],bs,iSoftType);
					}
				}
			if (bPeriodic) {
				umkEwald(&p[l*bs],bs,&q[k*bs],rsk,iSoftType,kd->fPeriod[0]);
				}
			else {
				umkGrav(&p[l*bs],bs,&q[k*bs],rsk,iSoftType);
				}
			}
		for (k=0;k<nbk;++k,++bc) {
			if (bPeriodic) {
				umkEwald(&p[l*bs],rs,&q[k*bs],bs,iSoftType,kd->fPeriod[0]);
				}
			else {
				umkGrav(&p[l*bs],rs,&q[k*bs],bs,iSoftType);
				}
			}
		if (bPeriodic) {
			umkEwald(&p[l*bs],rs,&q[k*bs],rsk,iSoftType,kd->fPeriod[0]);
			}
		else {
			umkGrav(&p[l*bs],rs,&q[k*bs],rsk,iSoftType);
			}
		}
	for (i=0;i<n;++i) {
		p[i].a[0] *= kd->G;
		p[i].a[1] *= kd->G;
		p[i].a[2] *= kd->G;
		p[i].dPot *= kd->G;
		}
	}


void kdGrav(KD kd,int iBlockSize,int iSoftType,int bPeriodic,
			int bVerbose)
{
	kd->iSoftType = iSoftType;
	kd->bPeriodic = bPeriodic;
	kd->iBlockSize = iBlockSize;
	Grav(kd,0,0,iSoftType,bPeriodic,bVerbose);
	}


void kdRestart(KD kd,int bVerbose)
{
	int l,k;

	ReadChkpt(kd,&l,&k);
	Grav(kd,l,k,kd->iSoftType,kd->bPeriodic,bVerbose);
	}


void kdOutAccel(KD kd,char *pszFile)
{
	FILE *fp;
	int i,iCnt;

	fp = fopen(pszFile,"w");
	assert(fp != NULL);
	fprintf(fp,"%d\n",kd->nParticles);
	iCnt = 0;
	for (i=0;i<kd->nGas;++i) {
		if (kd->bGas) fprintf(fp,"%.17g\n",kd->p[iCnt++].a[0]);
		else fprintf(fp,"0\n");
		}
	for (i=0;i<kd->nDark;++i) {
		if (kd->bDark) fprintf(fp,"%.17g\n",kd->p[iCnt++].a[0]);
		else fprintf(fp,"0\n");
		}
	for (i=0;i<kd->nStar;++i) {
		if (kd->bStar) fprintf(fp,"%.17g\n",kd->p[iCnt++].a[0]);
		else fprintf(fp,"0\n");
		}
	iCnt = 0;
	for (i=0;i<kd->nGas;++i) {
		if (kd->bGas) fprintf(fp,"%.17g\n",kd->p[iCnt++].a[1]);
		else fprintf(fp,"0\n");
		}
	for (i=0;i<kd->nDark;++i) {
		if (kd->bDark) fprintf(fp,"%.17g\n",kd->p[iCnt++].a[1]);
		else fprintf(fp,"0\n");
		}
	for (i=0;i<kd->nStar;++i) {
		if (kd->bStar) fprintf(fp,"%.17g\n",kd->p[iCnt++].a[1]);
		else fprintf(fp,"0\n");
		}
	iCnt = 0;
	for (i=0;i<kd->nGas;++i) {
		if (kd->bGas) fprintf(fp,"%.17g\n",kd->p[iCnt++].a[2]);
		else fprintf(fp,"0\n");
		}
	for (i=0;i<kd->nDark;++i) {
		if (kd->bDark) fprintf(fp,"%.17g\n",kd->p[iCnt++].a[2]);
		else fprintf(fp,"0\n");
		}
	for (i=0;i<kd->nStar;++i) {
		if (kd->bStar) fprintf(fp,"%.17g\n",kd->p[iCnt++].a[2]);
		else fprintf(fp,"0\n");
		}
	fclose(fp);
	}


void kdOutPot(KD kd,char *pszFile)
{
	FILE *fp;
	int i,iCnt;

	fp = fopen(pszFile,"w");
	assert(fp != NULL);
	fprintf(fp,"%d\n",kd->nParticles);
	iCnt = 0;
	for (i=0;i<kd->nGas;++i) {
		if (kd->bGas) fprintf(fp,"%.17g\n",kd->p[iCnt++].dPot);
		else fprintf(fp,"0\n");
		}
	for (i=0;i<kd->nDark;++i) {
		if (kd->bDark) fprintf(fp,"%.17g\n",kd->p[iCnt++].dPot);
		else fprintf(fp,"0\n");
		}
	for (i=0;i<kd->nStar;++i) {
		if (kd->bStar) fprintf(fp,"%.17g\n",kd->p[iCnt++].dPot);
		else fprintf(fp,"0\n");
		}
	fclose(fp);
	}










