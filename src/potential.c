/***************************************
 * Threading potential access routines *
 * Copyright (C) 1994 David T. Jones   *
 ***************************************/

/* This code may not be distributed, copied or used without the permission of the Copyright holder(s) */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>

#include "potential.h"

#ifndef FALSE
#define FALSE 0
#define TRUE 1
#endif

#define SCALEFAC 1000

#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))

/*
 * This module provides a constant interface to the pairwise and solvation
 * threading potential functions 
 */

/* Local conformation tables */
short int             sr_de[4][4][21][21][INTERVALS][TOPOMAX];

/* MR Residue-Residue interaction matrices */
short int             mr_de[4][4][21][21][NRR];

/* LR Residue-Residue interaction matrices */
short int             lr_de[4][4][21][21][NRR];

/* Residue accessibility matrices */
short int             acc_de[21][NACC];

/* Residue OOI number matrices */
short int             ooi_de[21][(OOIMAX-OOIMIN) / OOIDIV + 1];

const float mindist[TOPOMAX][4][4] =
{
  {
    { 2.71, 2.28, 3.14, 1.91 },
    { 2.56, 2.73, 3.15, 2.05 },
    { 1.86, 1.73, 1.65, 1.86 },
    { 2.55, 1.39, 3.16, 2.29 },
  },
  {
    { 4.38, 3.62, 3.38, 3.32 },
    { 4.12, 3.68, 2.66, 3.26 },
    { 3.08, 2.61, 2.30, 2.37 },
    { 4.18, 3.55, 2.71, 3.47 },
  },
  {
    { 3.89, 2.90, 3.04, 3.84 },
    { 2.61, 2.32, 2.26, 3.03 },
    { 2.79, 2.06, 2.53, 2.63 },
    { 3.58, 2.69, 2.73, 3.51 },
  },
  {
    { 3.80, 2.77, 2.68, 3.61 },
    { 2.98, 2.49, 1.97, 3.53 },
    { 2.67, 1.73, 2.67, 2.49 },
    { 3.45, 2.33, 2.63, 3.46 },
  },
  {
    { 3.50, 2.51, 3.06, 3.22 },
    { 2.54, 2.50, 1.69, 2.79 },
    { 2.80, 2.03, 2.76, 2.53 },
    { 3.16, 2.94, 2.61, 3.37 },
  },
  {
    { 3.75, 2.78, 3.07, 3.66 },
    { 2.45, 2.22, 2.02, 2.61 },
    { 2.98, 1.83, 3.07, 2.66 },
    { 3.33, 2.90, 2.79, 3.72 },
  },
  {
    { 3.78, 2.45, 3.16, 3.51 },
    { 2.81, 2.83, 2.57, 2.89 },
    { 2.89, 1.83, 3.08, 2.72 },
    { 3.57, 3.25, 2.69, 3.80 },
  },
  {
    { 3.96, 3.12, 3.00, 3.71 },
    { 3.09, 2.51, 2.90, 2.87 },
    { 3.01, 1.86, 3.09, 2.66 },
    { 3.87, 3.41, 2.70, 3.75 },
  },
  {
    { 3.71, 3.30, 3.01, 3.80 },
    { 3.49, 2.64, 1.88, 3.36 },
    { 3.09, 2.62, 3.01, 2.69 },
    { 3.37, 2.52, 2.68, 3.56 },
  },
  {
    { 3.91, 2.63, 3.10, 3.81 },
    { 2.72, 2.62, 2.13, 2.84 },
    { 2.98, 1.93, 3.17, 2.65 },
    { 3.53, 2.34, 2.74, 3.73 },
  },
  {
    { 3.96, 2.86, 3.05, 3.69 },
    { 2.62, 2.80, 2.11, 2.90 },
    { 3.09, 1.85, 3.05, 2.72 },
    { 3.96, 3.24, 2.70, 3.60 },
  }
};

const float maxdist[TOPOMAX][4][4] =
{
  {
    { 4.11, 5.27, 6.18, 2.86 },
    { 5.11, 6.42, 7.40, 3.92 },
    { 3.68, 4.92, 5.81, 2.57 },
    { 5.25, 6.40, 7.24, 3.97 },
  },
  {
    { 7.71, 8.60, 9.55, 6.37 },
    { 8.63, 9.80, 10.78, 7.38 },
    { 7.06, 8.13, 9.03, 5.84 },
    { 8.87, 9.75, 10.64, 7.46 },
  },
  {
    { 11.00, 12.22, 12.86, 9.80 },
    { 12.15, 13.23, 13.86, 10.91 },
    { 10.46, 11.48, 11.98, 9.16 },
    { 12.10, 13.28, 14.00, 10.89 },
  },
  {
    { 14.52, 15.42, 16.19, 13.15 },
    { 15.52, 16.53, 17.22, 14.21 },
    { 13.65, 14.75, 15.29, 12.52 },
    { 15.60, 16.53, 17.41, 14.33 },
  },
  {
    { 17.75, 18.79, 19.68, 16.50 },
    { 18.96, 20.08, 20.66, 17.63 },
    { 16.98, 18.00, 18.86, 15.73 },
    { 19.00, 20.06, 20.78, 17.68 },
  },
  {
    { 21.15, 22.38, 23.02, 19.93 },
    { 22.09, 23.00, 23.59, 20.83 },
    { 20.18, 21.25, 21.92, 18.91 },
    { 22.32, 23.55, 24.07, 21.09 },
  },
  {
    { 24.59, 25.49, 26.34, 23.27 },
    { 25.54, 26.60, 27.03, 24.28 },
    { 23.59, 24.67, 25.39, 22.42 },
    { 25.69, 26.81, 27.46, 24.39 },
  },
  {
    { 27.98, 28.73, 29.73, 26.64 },
    { 28.76, 29.69, 30.34, 27.52 },
    { 26.89, 27.92, 28.55, 25.51 },
    { 28.97, 29.96, 30.88, 27.75 },
  },
  {
    { 31.29, 32.26, 33.10, 29.99 },
    { 32.35, 33.51, 33.44, 31.08 },
    { 30.22, 31.15, 32.11, 29.00 },
    { 32.38, 33.34, 34.30, 31.19 },
  },
  {
    { 34.39, 35.43, 35.84, 33.20 },
    { 34.98, 35.86, 36.51, 34.06 },
    { 33.60, 34.57, 34.84, 32.36 },
    { 35.76, 36.78, 36.96, 34.55 },
  },
  {
    { 37.75, 38.83, 39.41, 36.52 },
    { 38.74, 39.60, 40.25, 37.35 },
    { 36.54, 37.68, 38.52, 35.32 },
    { 38.50, 39.55, 39.50, 37.26 },
  }
};


float           distscale[TOPOMAX][4][4];

extern void fail(char *);

/* Return accessibility range */
int
                accrange(int acc)
{
    if (acc < 12)
	return 0;
    else if (acc < 36)
	return 1;
    else if (acc < 44)
	return 2;
    else if (acc < 87)
	return 3;
    return 4;
}

/* Byte-swap for little-endian systems */
void byteswap(void *f, register int n)
{
    register unsigned int *vp;

    vp = (unsigned int *)f;

    while (n--)
    {
        *vp = ((*vp >> 24) & 0xFF) | ((*vp >> 8) & 0x0000FF00) | ((*vp << 8) & 0x00FF0000) | ((*vp << 24) & 0xFF000000);
	*vp++;
    }
}

/* Read SR potentials into array and convert to fixed point integer */
static int      cread_sr(FILE * ifp)
{
    int        a, b, i, j, k, t;
    float      (* Fsr_de)[4][4][21][21][INTERVALS][TOPOMAX];

    if ((Fsr_de = malloc(sizeof(*Fsr_de))) == NULL)
	fail("Out of memory in cread_sr!");
    
    if (!fread(Fsr_de, 1, sizeof(*Fsr_de), ifp))
	return TRUE;

#if !defined(linux) && !defined(__alpha)
    byteswap(Fsr_de, sizeof(*Fsr_de)/4);
#endif

    for (i=0; i<21; i++)
	for (j=0; j<21; j++)
	    for (a=0; a<4; a++)
		for (b=0; b<4; b++)
		    for (t=0; t<TOPOMAX; t++)
			for (k=0; k<INTERVALS; k++)
			    sr_de[a][b][i][j][k][t] = (*Fsr_de)[a][b][i][j][k][t] * SCALEFAC;

    free(Fsr_de);

    return FALSE;
}

/* Read non-SR potentials into array and convert to fixed point integer */
static int      cread_other(FILE * ifp)
{
    int        a, b, i, j, k;
    float      (* Fmr_de)[4][4][21][21][NRR];
    float      (* Flr_de)[4][4][21][21][NRR];
    float      (* Facc_de)[21][NACC];
    float      (* Fooi_de)[21][(OOIMAX-OOIMIN) / OOIDIV + 1];

    Fmr_de = malloc(sizeof(* Fmr_de));
    Flr_de = malloc(sizeof(* Flr_de));
    Facc_de = malloc(sizeof(* Facc_de));
    Fooi_de = malloc(sizeof(* Fooi_de));

    if (Fmr_de == NULL || Flr_de == NULL || Facc_de == NULL || Fooi_de == NULL)
	fail("Out of memory in cread_other!");

    if (!fread(Fmr_de, 1, sizeof(*Fmr_de), ifp) || !fread(Flr_de, 1, sizeof(*Flr_de), ifp) || !fread(Facc_de, 1, sizeof(*Facc_de), ifp) || !fread(Fooi_de, 1, sizeof(*Fooi_de), ifp))
	return TRUE;

#if !defined(linux) && !defined(__alpha)
    byteswap(Fmr_de, sizeof(*Fmr_de)/4);
    byteswap(Flr_de, sizeof(*Flr_de)/4);
    byteswap(Facc_de, sizeof(*Facc_de)/4);
    byteswap(Fooi_de, sizeof(*Fooi_de)/4);
#endif

    for (i=0; i<21; i++)
	for (j=0; j<21; j++)
	    for (a=0; a<4; a++)
		for (b=0; b<4; b++)
		    for (k=0; k<NRR; k++)
		    {
			mr_de[a][b][i][j][k] = (*Fmr_de)[a][b][i][j][k] * SCALEFAC;
			lr_de[a][b][i][j][k] = (*Flr_de)[a][b][i][j][k] * SCALEFAC;
		    }
    
    for (i=0; i<21; i++)
	for (j=0; j<NACC; j++)
	    acc_de[i][j] = (*Facc_de)[i][j] * SCALEFAC;

    for (i=0; i<21; i++)
	for (j=0; j<(OOIMAX-OOIMIN) / OOIDIV + 1; j++)
	    ooi_de[i][j] = (*Fooi_de)[i][j] * SCALEFAC;

    free(Fmr_de);
    free(Flr_de);
    free(Facc_de);
    free(Fooi_de);

    return FALSE;
}

/* Return potential term for at1(a) -> at2(b) */
int             pairpot(
			   const int a, const int b,
			   const int at1, const int at2,
			   const int t,
			   const float d
)
{
    int    k;

    if (t >= 23)
    {
	k = (int) ((d-4.0F) / RRWIDTH);
	if (k < 0)
	    k = 0;
	if (k >= NRR)
	    k = NRR - 1;
	
	return lr_de[at1][at2][a][b][k];
    }
    else if (t > TOPOMAX)
    {
	k = (int) ((d-4.0F) / RRWIDTH);
	if (k < 0)
	    k = 0;
	if (k >= NRR)
	    k = NRR - 1;
	
	return mr_de[at1][at2][a][b][k];
    }
    else if (t >= 1)
    {
	k = (int) ((d - mindist[t-1][at1][at2]) * distscale[t-1][at1][at2]);
	if (k < 0)
	    k = 0;
	else if (k >= INTERVALS)
	    return 0;
	return sr_de[at1][at2][a][b][k][t-1];
    }

    return 0;
}

/* Return potential term for residue a with given % accessibility */
int             solvpot(const int a, const int acc)
{
    return acc_de[a][accrange(acc)];
}

/* Load potentials - return TRUE on failure */
int             potinit(void)
{
    int             a, b, i;
    char            potname[160];
    FILE           *ifp;

    /* Read potentials */

    if (getenv("THREAD_DIR"))
    {
	strcpy(potname, getenv("THREAD_DIR"));
	if (potname[strlen(potname)-1] != '/')
	    strcat(potname, "/potentials.dat");
    }
    else
	strcpy(potname, "./data/potentials.dat");

    ifp = fopen(potname, "rb");
    if (!ifp)
	return TRUE;
    if (cread_sr(ifp) || cread_other(ifp))
    {
	fclose(ifp);
	return TRUE;
    }

    fclose(ifp);

    for (i = 0; i < TOPOMAX; i++)
	for (a = 0; a < 4; a++)
	    for (b = 0; b < 4; b++)
		distscale[i][a][b] = INTERVALS / (maxdist[i][a][b] - mindist[i][a][b]);

    return FALSE;
}
