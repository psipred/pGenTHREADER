/****************************************************************************
 *
 * pGen(Dom)THREADER version 8.2 (Jul 2009)
 * written by David T. Jones,
 * Department of Computer Science,
 * University College London,
 * Gower Street,
 * London
 *
 * Email: d.jones@cs.ucl.ac.uk
 *
 ****************************************************************************/

/*
 * This program attempts to align a sequence and its profile with multiple profiles and evaluate
 * each match via knowledge-based potentials.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <sys/resource.h>
#include <unistd.h>

#include "potential.h"

#ifndef FALSE
#define FALSE 0
#define TRUE 1
#endif

#define BIG 1000000000
#define VBIG 1e32F

/* "Tweakable" Parameters */

/* Scaling factor for fixed-point arithmetic */
#define SCALEFAC 1000

/* Max-size-of-problem parameters */
#define MAXSEQLEN 10000
#define MAXATOMS 10000
#define NALN 10

#define SQR(x) ((x)*(x))
#define ABS(x) ((x)>=0?(x):-(x))
#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))
#define CH malloc_verify(), printf("Heap OK at line : %d.\n",__LINE__);
#define vecprod(a,b,c) (a[0]=b[1]*c[2]-b[2]*c[1],a[1]=b[2]*c[0]-b[0]*c[2],a[2]=b[0]*c[1]-b[1]*c[0])
#define dotprod(a,b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])


/* Constants for calculating CB coords (Prot. Eng. Vol.2 p.121) */
#define	TETH_ANG 0.9128
#define CACBDIST 1.538

const float     ZERO = 0.0F, ONE = 1.0F, TWO = 2.0F, THREE = 3.0F;

char            brkid[11], seq1[MAXSEQLEN], seq2[10000];
char            modseq[MAXSEQLEN], tsstruc[MAXSEQLEN], coreflg[MAXSEQLEN];
char            ssstruc[10000], ssstrel[10000], *resid[MAXSEQLEN];
short           modsdx[MAXSEQLEN], modsdx2[MAXSEQLEN], gaps[MAXSEQLEN];
short           sstidx[MAXSEQLEN], tooi[MAXSEQLEN];
int             maxgplen[MAXSEQLEN];
float           targf1[MAXSEQLEN][20], targf2[MAXSEQLEN][20];
int             tpltsc1[MAXSEQLEN][20], tpltsc2[MAXSEQLEN][20], firstid,
                seq1len, seq2len, nseqs;
float           rmax;

int             lc1, lc2;

int           **pmat;

int             psidata;

struct hash
{
    short           resid, acc;
    char            inscode, sstruc;
}
hashtbl[MAXSEQLEN];

struct SSTRUC
{
    short           start, length, type;
}
sstlist[100];

float           e_contrib[MAXSEQLEN], e_cav[MAXSEQLEN], contrib_min, contrib_max;

struct HIT
{
    char            brkid[11];
    float           epair, esolv, nwsc, perco_a, perco_b;
    int             laln, lena, lenb, sfrom, sto, ffrom, fto;
} hits[100000];

typedef struct
{
    short           length;
    short           posn_a[MAXSEQLEN], posn_b[MAXSEQLEN];
}
ALNS;

int             sstlen;
FILE           *efp;

const char     *rnames[] =
{
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
    "GAP", "UNK"
};

const char     *ssnames[] =
{
    "COIL", "HELIX", "STRAND", "TURN"
};

enum sscodes
{
    COIL, HELIX, STRAND, TURN
};

const char     *atmnames[] =
{
    "CA", "CB", "O ", "N ", "C "
};

#define NPAIRS 5

const short     atompair[][2] =
{
    {CBATOM, CBATOM},
    {CBATOM, NATOM},
    {CBATOM, OATOM},
    {NATOM, CBATOM},
    {OATOM, CBATOM},
};

const char     *rescodes = "ARNDCQEGHILKMFPSTWYV-X";
const char     *sscodes = "CHE";

/* Amino acid composition of SWISS-PROT 47.0 */
float dbaaf[20] =
{
    0.078558, 0.053457, 0.041873, 0.053116, 0.015533, 0.039327, 0.066123, 0.069598, 0.022845, 0.059217,
    0.096416, 0.059123, 0.023859, 0.040095, 0.048484, 0.068456, 0.054442, 0.011580, 0.030628, 0.067271
};

short           mutcode[22] =
{
    16, 11, 3, 6, 19, 8, 3, 0, 18, 12, 12, 1, 19, 18, 15, 16, 15, 18, 13, 9, 20, 20
};

short ssmat[3][3] = {
    {  135, -46, -472 },
    { -200,  90, -374 },
    { -162, -398, 215 }
};

/*  BLOSUM 50 */
short           aamat2[23][23] =
{
    {500, -200, -100, -200, -100, -100, -100, 000, -200, -100, -200, -100, -100, -300, -100, 100, 000, -300, -200, 000, -200, -100, -100},
    {-200, 700, -100, -200, -400, 100, 000, -300, 000, -400, -300, 300, -200, -300, -300, -100, -100, -300, -100, -300, -100, 000, -100},
    {-100, -100, 700, 200, -200, 000, 000, 000, 100, -300, -400, 000, -200, -400, -200, 100, 000, -400, -200, -300, 400, 000, -100},
    {-200, -200, 200, 800, -400, 000, 200, -100, -100, -400, -400, -100, -400, -500, -100, 000, -100, -500, -300, -400, 500, 100, -100},
    {-100, -400, -200, -400, 1300, -300, -300, -300, -300, -200, -200, -300, -200, -200, -400, -100, -100, -500, -300, -100, -300, -300, -200},
    {-100, 100, 000, 000, -300, 700, 200, -200, 100, -300, -200, 200, 000, -400, -100, 000, -100, -100, -100, -300, 000, 400, -100},
    {-100, 000, 000, 200, -300, 200, 600, -300, 000, -400, -300, 100, -200, -300, -100, -100, -100, -300, -200, -300, 100, 500, -100},
    {000, -300, 000, -100, -300, -200, -300, 800, -200, -400, -400, -200, -300, -400, -200, 000, -200, -300, -300, -400, -100, -200, -200},
    {-200, 000, 100, -100, -300, 100, 000, -200, 1000, -400, -300, 000, -100, -100, -200, -100, -200, -300, 200, -400, 000, 000, -100},
    {-100, -400, -300, -400, -200, -300, -400, -400, -400, 500, 200, -300, 200, 000, -300, -300, -100, -300, -100, 400, -400, -300, -100},
    {-200, -300, -400, -400, -200, -200, -300, -400, -300, 200, 500, -300, 300, 100, -400, -300, -100, -200, -100, 100, -400, -300, -100},
    {-100, 300, 000, -100, -300, 200, 100, -200, 000, -300, -300, 600, -200, -400, -100, 000, -100, -300, -200, -300, 000, 100, -100},
    {-100, -200, -200, -400, -200, 000, -200, -300, -100, 200, 300, -200, 700, 000, -300, -200, -100, -100, 000, 100, -300, -100, -100},
    {-300, -300, -400, -500, -200, -400, -300, -400, -100, 000, 100, -400, 000, 800, -400, -300, -200, 100, 400, -100, -400, -400, -200},
    {-100, -300, -200, -100, -400, -100, -100, -200, -200, -300, -400, -100, -300, -400, 1000, -100, -100, -400, -300, -300, -200, -100, -200},
    {100, -100, 100, 000, -100, 000, -100, 000, -100, -300, -300, 000, -200, -300, -100, 500, 200, -400, -200, -200, 000, 000, -100},
    {000, -100, 000, -100, -100, -100, -100, -200, -200, -100, -100, -100, -100, -200, -100, 200, 500, -300, -200, 000, 000, -100, 000},
    {-300, -300, -400, -500, -500, -100, -300, -300, -300, -300, -200, -300, -100, 100, -400, -400, -300, 1500, 200, -300, -500, -200, -300},
    {-200, -100, -200, -300, -300, -100, -200, -300, 200, -100, -100, -200, 000, 400, -300, -200, -200, 200, 800, -100, -300, -200, -100},
    {000, -300, -300, -400, -100, -300, -300, -400, -400, 400, 100, -300, 100, -100, -300, -200, 000, -300, -100, 500, -400, -300, -100},
    {-200, -100, 400, 500, -300, 000, 100, -100, 000, -400, -400, 000, -300, -400, -200, 000, 000, -500, -300, -400, 500, 200, -100},
    {-100, 000, 000, 100, -300, 400, 500, -200, 000, -300, -300, 100, -100, -400, -100, 000, -100, -200, -200, -300, 200, 500, -100},
    {-100, -100, -100, -100, -200, -100, -100, -200, -100, -100, -100, -100, -100, -200, -200, -100, 000, -300,
     -100, -100, -100, -100, -100}
};

/* BLOSUM 62 */
short           aamat[23][23] =
{
    {400, -100, -200, -200, 000, -100, -100, 000, -200, -100, -100, -100, -100, -200, -100, 100, 000, -300, -200, 000, -200, -100, -100},
    {-100, 500, 000, -200, -300, 100, 000, -200, 000, -300, -200, 200, -100, -300, -200, -100, -100, -300, -200, -300, -100, 000, -100},
    {-200, 000, 600, 100, -300, 000, 000, 000, 100, -300, -300, 000, -200, -300, -200, 100, 000, -400, -200, -300, 300, 000, -100},
    {-200, -200, 100, 600, -300, 000, 200, -100, -100, -300, -400, -100, -300, -300, -100, 000, -100, -400, -300, -300, 400, 100, -100},
    {000, -300, -300, -300,1000, -300, -400, -300, -300, -100, -100, -300, -100, -200, -300, -100, -100, -200, -200, -100, -300, -300, -100},
    {-100, 100, 000, 000, -300, 500, 200, -200, 000, -300, -200, 100, 000, -300, -100, 000, -100, -200, -100, -200, 000, 300, -100},
    {-100, 000, 000, 200, -400, 200, 500, -200, 000, -300, -300, 100, -200, -300, -100, 000, -100, -300, -200, -200, 100, 400, -100},
    {000, -200, 000, -100, -300, -200, -200, 600, -200, -400, -400, -200, -300, -300, -200, 000, -200, -200, -300, -300, -100, -200, -100},
    {-200, 000, 100, -100, -300, 000, 000, -200, 800, -300, -300, -100, -200, -100, -200, -100, -200, -200, 200, -300, 000, 000, -100},
    {-100, -300, -300, -300, -100, -300, -300, -400, -300, 400, 200, -300, 100, 000, -300, -200, -100, -300, -100, 300, -300, -300, -100},
    {-100, -200, -300, -400, -100, -200, -300, -400, -300, 200, 400, -200, 200, 000, -300, -200, -100, -200, -100, 100, -400, -300, -100},
    {-100, 200, 000, -100, -300, 100, 100, -200, -100, -300, -200, 500, -100, -300, -100, 000, -100, -300, -200, -200, 000, 100, -100},
    {-100, -100, -200, -300, -100, 000, -200, -300, -200, 100, 200, -100, 500, 000, -200, -100, -100, -100, -100, 100, -300, -100, -100},
    {-200, -300, -300, -300, -200, -300, -300, -300, -100, 000, 000, -300, 000, 600, -400, -200, -200, 100, 300, -100, -300, -300, -100},
    {-100, -200, -200, -100, -300, -100, -100, -200, -200, -300, -300, -100, -200, -400, 700, -100, -100, -400, -300, -200, -200, -100, -100},
    {100, -100, 100, 000, -100, 000, 000, 000, -100, -200, -200, 000, -100, -200, -100, 400, 100, -300, -200, -200, 000, 000, -100},
    {000, -100, 000, -100, -100, -100, -100, -200, -200, -100, -100, -100, -100, -200, -100, 100, 500, -200, -200, 000, -100, -100, -100},
    {-300, -300, -400, -400, -200, -200, -300, -200, -200, -300, -200, -300, -100, 100, -400, -300, -200, 1100, 200, -300, -400, -300, -100},
    {-200, -200, -200, -300, -200, -100, -200, -300, 200, -100, -100, -200, -100, 300, -300, -200, -200, 200, 700, -100, -300, -200, -100},
    {000, -300, -300, -300, -100, -200, -200, -300, -300, 300, 100, -200, 100, -100, -200, -200, 000, -300, -100, 400, -300, -200, -100},
    {-200, -100, 300, 400, -300, 000, 100, -100, 000, -300, -400, 000, -300, -300, -200, 000, -100, -400, -300, -300, 400, 100, -100},
    {-100, 000, 000, 100, -300, 300, 400, -200, 000, -300, -300, 100, -100, -300, -100, 000, -100, -300, -200, -200, 100, 400, -100},
    {-100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, 000}
};

/* Local conformation tables */
extern short int sr_de[4][4][21][21][INTERVALS][TOPOMAX];

/* MR Residue-Residue interaction matrices */
extern short int mr_de[4][4][21][21][NRR];

/* LR Residue-Residue interaction matrices */
extern short int lr_de[4][4][21][21][NRR];

/* Residue accessibility matrices */
extern short int acc_de[21][NACC];

/* Residue OOI number matrices */
extern short int ooi_de[21][(OOIMAX - OOIMIN) / OOIDIV + 1];

extern float    maxdist[TOPOMAX][4][4];
extern float    mindist[TOPOMAX][4][4];
extern float    distscale[TOPOMAX][4][4];

/* Structure generation arrays */
float           (**distmat)[4][4];
float           x[MAXSEQLEN][5], y[MAXSEQLEN][5], z[MAXSEQLEN][5];
float           phi[MAXSEQLEN], psi[MAXSEQLEN], omega[MAXSEQLEN];
short           trelacc[MAXSEQLEN];
float           tmat[MAXSEQLEN][22];

/* Vars for Gotoh-N-W alignment routine */
int           **pat;

int             nsst;

/* Control parameters */
int             SHUFCOUNT = 0;
int             MINSSTREL = 20;
float           SC_PARAM[20];	/* TEST parameters */
int             GP_OPEN[3] = { 748, 1263, 711 };
int             GP_EXT[3] = { 64, 162, 284 };
int             NWCUTOFF = 40;
float           DCUTOFF = 10.0F;
int             ssscore = 100;
int             prtalnflg = FALSE,  verbflg = FALSE, motifflg = FALSE, filtgapflg = FALSE, localmode = FALSE;
int             mod3flg = FALSE;
int             htmlflg = FALSE;
int             seqssflg = FALSE;
int             quickalnflg = FALSE;
int             ooisolv = FALSE;
int             ssstflg = FALSE;
char            *ssfname, *htmname;
char            motpref[160];

/* Dump a rude message to standard error and exit */
void
                fail(char *errstr)
{
    fprintf(stderr, "\n*** %s\n\n", errstr);
    exit(-1);
}

void           *allocmat(int rows, int columns, int size, int clrflg)
{
    int             i;
    void          **p;

    p = malloc(rows * sizeof(void *));

    if (p == NULL)
	fail("allocmat: malloc [] failed!");
    if (clrflg)
    {
	for (i = 0; i < rows; i++)
	    if ((p[i] = calloc(columns, size)) == NULL)
		fail("allocmat: calloc [][] failed!");
    }
    else
	for (i = 0; i < rows; i++)
	    if ((p[i] = malloc(columns * size)) == NULL)
		fail("allocmat: malloc [][] failed!");

    return p;
}

void            freemat(void *p, int rows)
{
    int             i;

    for (i = rows - 1; i >= 0; i--)
	free(((void **) p)[i]);
    free(p);
}

/* Convert AA letter to numeric code (0-22) */
int
                aanum(int ch)
{
    const static int aacvs[] =
    {
	999, 0, 20, 4, 3, 6, 13, 7, 8, 9, 22, 11, 10, 12, 2,
	22, 14, 5, 1, 15, 16, 22, 19, 17, 22, 18, 21
    };

    return (isalpha(ch) ? aacvs[ch & 31] : 22);
}

/* Convert string to lower case */
void
                lc_str(char *str)
{
    int             i;

    for (i = strlen(str) - 1; i >= 0; i--)
	if (isupper(str[i]))
	    str[i] = tolower(str[i]);
}


unsigned int rng_x=123456789, rng_y=362436069, rng_z=521288629, rng_w=88675123;

/* Fast Marsaglia XOR-shift RNG with period 2^128-1 */
unsigned int xor128(void)
{
    unsigned int t;

    t = (rng_x^(rng_x<<11));
    rng_x = rng_y;
    rng_y = rng_z;
    rng_z = rng_w;

    return rng_w = (rng_w^(rng_w>>19))^(t^(t>>8));
}

/* Generate random number 0<=x<1 */
#define ran0()  (xor128()*(1.0/4294967296.0))

/* randint(a,b) : return random integer a <= n <= b */
#define randint(low,high) ((low) + (int)(((high)-(low)+1) * ran0()))


/* Attempt to generate a random state unique to this process/host/time */
void
                randomise(void)
{
    int i;

    rng_x = (unsigned int)time(NULL);
    rng_y = (unsigned int)getppid();
    rng_z = (unsigned int)gethostid();
    rng_w = (unsigned int)getpid();

    if (verbflg)
	printf("Random seeds: %u %u %u %u\n", rng_x, rng_y, rng_z, rng_w);

    /* Warm up the generator */
    for (i=0; i<100; i++)
	xor128();
}


/* Perform random shuffle of sequence skipping GAPS/UNKs */
void            shufseqng(char *s, int len)
{
    int             i, r;
    char            temp;

    for (i = len - 1; i >= 1; i--)
	if (s[i] < 20)
	{
	    r = ran0() * (i + 1);
	    if (s[r] < 20)
	    {
		temp = s[i];
		s[i] = s[r];
		s[r] = temp;
	    }
	}
}

/* Perform random shuffle of sequence */
void            shufseq(char *s, int len)
{
    int             i, r;
    char            temp;

    for (i = len - 1; i >= 1; i--)
    {
	r = ran0() * (i + 1);
	temp = s[i];
	s[i] = s[r];
	s[r] = temp;
    }
}

/* Perform random rotation of sequence */
void            rotseq(char *s, int len)
{
    int             i, n;
    char            temp;

    n = randint(1, len - 1);

    while (n--)
    {
	temp = s[0];
	for (i = 0; i < len - 1; i++)
	    s[i] = s[i + 1];
	s[len - 1] = temp;
    }
}

/* Print alignment */
void
                prtalign(ALNS * newtem)
{
    int             i, b, id = 0, nb, p1, p2, seqflg;
    char            r1, r2, sim;

    nb = newtem->length / 60 + 1;

    for (b = 0; b < nb; b++)
    {
	printf("         ");
	for (i = 0; i < 58; i++)
	{
	    if (b * 60 + i + 3 > newtem->length)
		break;
	    p1 = newtem->posn_a[b * 60 + i + 3];
	    if (p1 && !(p1 % 10))
	    {
		printf("%3d", p1);
		i += 2;
	    }
	    else
		printf(" ");
	}
	putchar('\n');
	for (seqflg = 0; seqflg < 2; seqflg++)
	{
	    if (seqflg == 1)
		printf("%-8.8s ", brkid);
	    else
		printf("         ");
	    for (i = 0; i < 60; i++)
	    {
		if (b * 60 + i >= newtem->length)
		    break;
		p1 = newtem->posn_a[b * 60 + i + 1];

		switch (seqflg)
		{
		case 0:
		    r1 = (p1) ? sscodes[tsstruc[p1 - 1]] : '-';
		    break;
		case 1:
		    r1 = (p1) ? rescodes[seq1[p1 - 1]] : '-';
		    break;
		}

		putchar(r1);
	    }
	    putchar('\n');
	}
	printf("         ");
	for (i = 0; i < 60; i++)
	{
	    if (b * 60 + i >= newtem->length)
		break;
	    p1 = newtem->posn_a[b * 60 + i + 1];
	    p2 = newtem->posn_b[b * 60 + i + 1];
	    r1 = (p1) ? rescodes[seq1[p1 - 1]] : '-';
	    r2 = (p2) ? rescodes[seq2[p2 - 1]] : '-';
	    if (r1 == r2 && r1 != '-')
	    {
		id++;
		sim = '|';
	    }
	    else
		sim = ' ';
	    putchar(sim);
	}
	putchar('\n');
	printf("%-8.8s ", motpref[0] ? motpref : "Query");
	for (i = 0; i < 60; i++)
	{
	    if (b * 60 + i >= newtem->length)
		break;
	    p2 = newtem->posn_b[b * 60 + i + 1];
	    r2 = (p2) ? rescodes[seq2[p2 - 1]] : '-';
	    putchar(r2);
	}
	putchar('\n');
	if (ssstflg && seqssflg)
	{
	    printf("         ");
	    for (i = 0; i < 60; i++)
	    {
		if (b * 60 + i >= newtem->length)
		    break;
		p2 = newtem->posn_b[b * 60 + i + 1];
		r2 = (p2) ? sscodes[ssstruc[p2 - 1]] : '-';
		putchar(r2);
	    }
	    putchar('\n');
	}
	printf("         ");
	for (i = 0; i < 58; i++)
	{
	    if (b * 60 + i + 3 > newtem->length)
		break;
	    p2 = newtem->posn_b[b * 60 + i + 3];
	    if (p2 && !(p2 % 10))
	    {
		printf("%3d", p2);
		i += 2;
	    }
	    else
		printf(" ");
	}
	puts("\n\n");
    }
    printf("Percentage Identity = %3.1f.\n\n", 100.0F * id / MIN(seq1len, seq2len));
}

/* Print alignment in HTML format */
void
                htmlalign(ALNS * newtem)
{
    int             i, b, id = 0, nb, p1, p2, seqflg;
    char            r1, r2;
    FILE *hfp;

    hfp = fopen(htmname, "a");
    if (!hfp)
	fail("Cannot append to HTML output file!");

    nb = newtem->length / 60 + 1;

    fprintf(hfp, "<hr><body text=\"#ffffff\" bgcolor=\"#000000\"><tt><pre><br><a name=\"%s\"></a>", brkid);

    for (b = 0; b < nb; b++)
    {
	fprintf(hfp, "         ");
	for (i = 0; i < 58; i++)
	{
	    if (b * 60 + i + 3 > newtem->length)
		break;
	    p1 = newtem->posn_a[b * 60 + i + 3];
	    if (p1 && !(p1 % 10))
	    {
		fprintf(hfp, "<font color=\"#ffffff\">%3d</font>", p1);
		i += 2;
	    }
	    else
		fprintf(hfp, " ");
	}
	fputc('\n', hfp);
	for (seqflg = 0; seqflg < 2; seqflg++)
	{
	    if (seqflg == 1)
		fprintf(hfp, "<font color=\"#ffffff\">%-8.8s</font> ", brkid);
	    else
		fprintf(hfp, "         ");
	    for (i = 0; i < 60; i++)
	    {
		if (b * 60 + i >= newtem->length)
		    break;
		p1 = newtem->posn_a[b * 60 + i + 1];
		p2 = newtem->posn_b[b * 60 + i + 1];

		switch (seqflg)
		{
		case 0:
		    r1 = (p1) ? sscodes[tsstruc[p1 - 1]] : '-';
		    switch(r1)
		    {
		    case 'H':
			fprintf(hfp, "<font color=\"#ff00ff\">%c</font>", r1);
			break;
		    case 'E':
			fprintf(hfp, "<font color=\"#ffff00\">%c</font>", r1);
			break;
		    default:
			fprintf(hfp, "<font color=\"#ffffff\">%c</font>", r1);
			break;
		    }
		    break;
		case 1:
		    r1 = (p1) ? rescodes[seq1[p1 - 1]] : '-';
		    r2 = (p2) ? rescodes[seq2[p2 - 1]] : '-';
		    if (r1 == r2 && r1 != '-')
			fprintf(hfp, "<font color=\"#00ff00\">%c</font>", r1);
		    else if (p1 && p2 && aamat[seq1[p1 - 1]][seq2[p2 - 1]] > 100)
			fprintf(hfp, "<font color=\"#ff0000\">%c</font>", r1);
		    else
			fprintf(hfp, "<font color=\"#6060ff\">%c</font>", r1);
		    break;
		}
	    }
	    fputc('\n', hfp);
	}
	fprintf(hfp, "<font color=\"#ffffff\">%-8.8s</font> ", motpref[0] ? motpref : "Query");

	for (i = 0; i < 60; i++)
	{
	    if (b * 60 + i >= newtem->length)
		break;
	    p1 = newtem->posn_a[b * 60 + i + 1];
	    p2 = newtem->posn_b[b * 60 + i + 1];
	    r1 = (p1) ? rescodes[seq1[p1 - 1]] : '-';
	    r2 = (p2) ? rescodes[seq2[p2 - 1]] : '-';
	    if (r1 == r2 && r1 != '-')
	    {
		id++;
		fprintf(hfp, "<font color=\"#00ff00\">%c</font>", r2);
	    }
	    else if (p1 && p2 && aamat[seq1[p1 - 1]][seq2[p2 - 1]] > 100)
		fprintf(hfp, "<font color=\"#ff0000\">%c</font>", r2);
	    else
		fprintf(hfp, "<font color=\"#6060ff\">%c</font>", r2);
	}
	fputc('\n', hfp);
	if (ssstflg && seqssflg)
	{
	    fprintf(hfp, "         ");
	    for (i = 0; i < 60; i++)
	    {
		if (b * 60 + i >= newtem->length)
		    break;
		p2 = newtem->posn_b[b * 60 + i + 1];
		r2 = (p2) ? sscodes[ssstruc[p2 - 1]] : '-';
		switch(r2)
		{
		case 'H':
		    fprintf(hfp, "<font color=\"#ff00ff\">%c</font>", r2);
		    break;
		case 'E':
		    fprintf(hfp, "<font color=\"#ffff00\">%c</font>", r2);
		    break;
		default:
		    fprintf(hfp, "<font color=\"#ffffff\">%c</font>", r2);
		    break;
		}
	    }
	}
	fputc('\n', hfp);
	fprintf(hfp, "         ");
	for (i = 0; i < 58; i++)
	{
	    if (b * 60 + i + 3 > newtem->length)
		break;
	    p2 = newtem->posn_b[b * 60 + i + 3];
	    if (p2 && !(p2 % 10))
	    {
		fprintf(hfp, "<font color=\"#ffffff\">%3d</font>", p2);
		i += 2;
	    }
	    else
		fprintf(hfp, " ");
	}
	fprintf(hfp, "\n\n");
    }
    fprintf(hfp, "</pre><font color=\"#ffffff\">Percentage Identity = %3.1f%%</font></body><P>\n", 100.0F * id / MIN(seq1len, seq2len));

    fclose(hfp);
}


/* Compute energy sums for final model */
void
                e_final(float *pair, float *solv, float *perco_a, float *perco_b)
{
    int             i, j, n, t, npair = 0, nsolv = 0, ooisum = 0, nids = 0,
	nlr, nsr;
    short           atpair, ooi[MAXSEQLEN], resa, resb;
    float           e, sr = ZERO, lr = ZERO;


    *pair = *solv = ZERO;
    nlr = nsr = 0;

    for (i = 0; i < seq1len; i++)
    {
	if (modsdx[i] >= 0 && seq1[i] == seq2[modsdx[i]])
	    nids++;
	ooi[i] = 0;
	e_contrib[i] = ZERO;
    }
    contrib_min = VBIG;
    contrib_max = -VBIG;

    if (ooisolv)
	for (i = 0; i < seq1len; i++)
	    if (seq1[i] < GAP && modsdx[i] >= 0)
		for (j = i + 1; j < seq1len; j++)
		    if (seq1[j] < GAP && modsdx[j] >= 0 && distmat[i][j][CBATOM][CBATOM] <= 10.0F)
		    {
			ooi[i]++;
			ooi[j]++;
		    }

    for (i = 0; i < seq1len; i++)
	if (seq1[i] < GAP && modsdx[i] >= 0)
	{
	    for (j = i + 1; j < seq1len; j++)
		if (seq1[j] < GAP && modsdx[j] >= 0)
		{
		    t = modsdx[j] - modsdx[i];
		    for (atpair = 0; atpair < NPAIRS; atpair++)
		    {
			resa = seq2[modsdx[i]];
			resb = seq2[modsdx[j]];
			if (seq2[modsdx[i]] < GAP && seq2[modsdx[j]] < GAP && coreflg[i] && coreflg[j] && distmat[i][j][((resa == GLY && atompair[atpair][0] == CBATOM) ? CAATOM : atompair[atpair][0])][((resb == GLY && atompair[atpair][1] == CBATOM) ? CAATOM : atompair[atpair][1])] < DCUTOFF)
			{
			    e = pairpot(resa, resb, (resa == GLY && atompair[atpair][0] == CBATOM) ? CAATOM : atompair[atpair][0], (resb == GLY && atompair[atpair][1] == CBATOM) ? CAATOM : atompair[atpair][1], t, distmat[i][j][((resa == GLY && atompair[atpair][0] == CBATOM) ? CAATOM : atompair[atpair][0])][((resb == GLY && atompair[atpair][1] == CBATOM) ? CAATOM : atompair[atpair][1])]);
/*                          printf("%c %c %s %s %3d %5.2f %7.3f\n", rescodes[resa], rescodes[resb], atmnames[atompair[atpair][0]], atmnames[atompair[atpair][1]], t, distmat[i][j][((resa == GLY && atompair[atpair][0] == CBATOM) ? CAATOM : atompair[atpair][0])][((resb == GLY && atompair[atpair][1] == CBATOM) ? CAATOM : atompair[atpair][1])], e); */

			    if (t <= TOPOMAX)
			    {
				sr += e;
				nsr++;
			    }
			    else
			    {
				lr += e;
				nlr++;
			    }

			    if (t <= TOPOMAX)
			    {
				e_contrib[i] += e;
				e_contrib[j] += e;
			    }

			    npair++;
			}
		    }
		}

	    if (seq2[modsdx[i]] < GAP)
	    {
		if (!ooisolv)
		    *solv += solvpot(seq2[modsdx[i]], trelacc[i]);
		nsolv++;
	    }
	}


    if (ooisolv)
    {
	for (i = 0; i < seq1len; i++)
	    if (modsdx[i] >= 0 && seq2[modsdx[i]] < GAP)
		ooisum += ooi_de[seq2[modsdx[i]]][MAX(0, MIN(OOIMAX, ooi[i]) - OOIMIN) / OOIDIV];

/*    printf("sr = %f lr = %f ooi = %f\n", sr, lr, av_ooi); */

	*solv = ooisum;
    }

    for (i = 0; i < seq1len; i++)
    {
	for (e = n = 0, j = -7; j <= 7; j++)
	    if (i + j >= 0 && i + j < seq1len)
	    {
		e += e_contrib[i + j];
		n++;
	    }
	e_cav[i] = e / n;
    }

    for (i = 0; i < seq1len; i++)
    {
	if (e_cav[i] < contrib_min)
	    contrib_min = e_cav[i];
	if (e_cav[i] > contrib_max)
	    contrib_max = e_cav[i];
    }

    *pair = sr + lr;

    if (verbflg)
    {
	printf("SR = %f, LR = %f\n", -RTCONST * sr / SCALEFAC, -RTCONST * lr / SCALEFAC);
	printf("avSR = %f, avLR = %f, avTOT = %f\n", -RTCONST * sr / SCALEFAC / nsr, -RTCONST * lr / SCALEFAC / nlr, -RTCONST * (sr + lr) / SCALEFAC / (nlr + nsr));
    }

    *perco_a = 100.0 * nsolv / seq1len;
    *perco_b = 100.0 * nsolv / seq2len;
}


/*
 * Trace back highest scoring path
 */
void
                trace(short *posa, short *posb, int mati, int matj,
		      int pati, int patj, int lasti, int lastj, short *n)
{
    int             pij = pat[pati][patj], i, j;

    for (i = lasti + 1; i < mati; i++)
    {
	*(++posa) = i;
	*(++posb) = 0;
	(*n)++;
    }
    for (j = lastj + 1; j < matj; j++)
    {
	*(++posa) = 0;
	*(++posb) = j;
	(*n)++;
    }
    *(++posa) = mati;
    *(++posb) = matj;
    (*n)++;

    if (!pij)
	return;

    if (pij == 1)
	trace(posa, posb, mati + 1, matj + 1, pati + 1, patj + 1,
	      mati, matj, n);
    if (pij < 1)
	trace(posa, posb, mati + 1, matj - pij, pati + 1, patj - pij,
	      mati, matj, n);
    if (pij > 1)
	trace(posa, posb, mati + pij, matj + 1, pati + pij, patj + 1,
	      mati, matj, n);
}

int
                seqscore(const char *seq1, const char *seq2, ALNS * aln, const int alntype,
                         const int seq1len, const int seq2len)
{
    short          *posa, *posb;
    int             trace_back = (aln != NULL);
    int             now = 0, last = 1;
    int             pati, patj, mati, matj, i, j, k, l;
    int             toprows[MAXSEQLEN + 1], topcol, toprow;
    int             maxrows[MAXSEQLEN + 1], maxscore, maxflg = FALSE, maxcol,
	maxrow, diag, row, col, envclass, gap_open, gap_extend;
    int             mat[2][MAXSEQLEN + 1];


    if (trace_back)
	for (i = 1; i <= seq1len; i++)
	    pat[i][seq2len] = 0;

    for (j = seq2len; j > 0; j--)
    {

	if (trace_back)
	    pat[seq1len][j] = 0;

	for (i = seq1len; i > 0; i--)
	{
	    /* Get matrix element */

	    if (!alntype)
	    {
		/* Classic GenTHREADER sequence->profile */
		if (seq2[j-1] < 20)
		    mat[now][i] = tpltsc1[i - 1][seq2[j - 1]];
		else
		    mat[now][i] = -1;
	    }
	    else if (alntype == 1)
	    {
	        //profile-profile aln
		int score = 0, sumwt = 0;
		float profprof = 0, ppsum = 0;

		if (seq1[i - 1] < 20)
		{
		    score = tpltsc2[j - 1][seq1[i - 1]] * 501;
		    sumwt = 501;
		}

		if (seq2[j - 1] < 20)
		{
		    score += tpltsc1[i - 1][seq2[j - 1]] * 659;
		    sumwt += 659;
		}

		for (k=0; k<20; k++)
		{
		    int vecwt = (tpltsc1[i - 1][k] > 0);

		    profprof += targf1[i - 1][k] * tpltsc2[j - 1][k] * vecwt;
		    ppsum += targf1[i - 1][k] * vecwt;
		}

		for (k=0; k<20; k++)
		{
		    int vecwt = (tpltsc2[j - 1][k] > 0);

		    profprof += targf2[j - 1][k] * tpltsc1[i - 1][k] * vecwt;
		    ppsum += targf2[j - 1][k] * vecwt;
		}

		if (ppsum > 0.0F)
		    profprof /= ppsum;

		//		fprintf(stderr,"profprof = %f\n", profprof);
		if (sumwt > 0)
		    mat[now][i] = 0.154133F * profprof + 0.761082F * score / sumwt;
		else
		    mat[now][i] = -1;
	    }
	    else
		fail("Unknown alignment type!");


	    if (alntype && ssstflg && ssstrel[j-1] >= MINSSTREL)
	      mat[now][i] += ssmat[ssstruc[j-1]][tsstruc[i-1]];

	    if (seq1[i-1] < 20 && seq2[j-1] < 20)
		mat[now][i] += 0.087529F * solvpot(seq2[j-1], trelacc[i-1]);

	    if (j != seq2len && i != seq1len)
	    {
		diag = mat[last][i + 1];

		maxrow = maxrows[i];
		toprow = toprows[i];

		envclass = tsstruc[i - 1];

		gap_open = GP_OPEN[envclass];
		gap_extend = GP_EXT[envclass];

		if (toprow)
		    row = maxrow - gap_open - (toprow - j) * gap_extend + gap_extend;
		else
		    row = maxrow - gap_open;
		if (topcol)
		    col = maxcol - gap_open - (topcol - i) * gap_extend + gap_extend;
		else
		    col = maxcol - gap_open;

		if (diag > col && diag > row)
		{
		    mat[now][i] += diag;
		    if (trace_back)
			pat[i][j] = 1;
		}
		else
		{
		    if (row > col)
		    {
			mat[now][i] += row;
			if (trace_back)
			    pat[i][j] = -(toprow - j);
		    }
		    else
		    {
			mat[now][i] += col;
			if (trace_back)
			    pat[i][j] = topcol - i;
		    }
		}

		//fprintf(stderr,"i=%d j=%d toprow=%d maxrow=%d mat=%d\n", i, j, toprow, maxrow, mat[now][i]);

		if (diag > maxrows[i])
		{
		    maxrows[i] = diag;
		    toprows[i] = j + 1;
		}
		if (diag > maxcol)
		{
		    maxcol = diag;
		    topcol = i + 1;
		}

		if (i == 1 || j == 1)
		{
		    if (!maxflg || mat[now][i] > maxscore)
		    {
			maxflg = TRUE;
			maxscore = mat[now][i];
/*			printf("MAXSCORE = %d\n"); */
			if (trace_back)
			{
			    pati = i;
			    patj = j;
			    mati = matj = 1;
			    if (i == 1)
				matj = j;
			    if (j == 1)
				mati = i;
			}
		    }
		}
	    }
	    else if (j == seq2len)
	    {
		maxrows[i] = mat[now][i];
		toprows[i] = j;
	    }
	    else if (i == seq1len)
	    {
		maxcol = mat[now][i];
		topcol = i;
	    }
	}
	now = !now;
	last = !last;
    }

    if (!trace_back)
	return (maxscore/100);

    posa = aln->posn_a;
    posb = aln->posn_b;
    aln->length = 0;

    //call to trace
    trace(posa, posb, mati, matj, pati, patj, 0, 0, &aln->length);

    posa += aln->length;
    posb += aln->length;

    if (*posa == seq1len)
    {
	for (i = *(posb) + 1; i <= seq2len; i++)
	{
	    *(++posa) = 0;
	    *(++posb) = i;
	    (aln->length)++;
	}
	if (aln->length > MAXSEQLEN)
	    fail("score : max. align length exceeded!");
	return (maxscore/100);
    }
    if (*posb == seq2len)
	for (i = *(posa) + 1; i <= seq1len; i++)
	{
	    *(++posb) = 0;
	    *(++posa) = i;
	    (aln->length)++;
	}
    if (aln->length > MAXSEQLEN)
	fail("score : max. align length exceeded!");


    return (maxscore/100);


}


/* Align sequence to structural template */
void
                seqfit(float *epair, float *esolv, float *perco_a, float *perco_b, ALNS * tplt)
{
    int             i, j, k;

    /* Build structure gap penalty array */
    for (i = 0; i < nsst; i++)
	if (!i)
	    for (j = 0; j < sstlist[0].start; j++)
		maxgplen[j] = sstlist[0].start - j;
	else
	    for (j = sstlist[i - 1].start + sstlist[i - 1].length; j < sstlist[i].start; j++)
		maxgplen[j] = sstlist[i].start - j;
    for (j = sstlist[nsst - 1].start + sstlist[nsst - 1].length; j < seq1len; j++)
	maxgplen[j] = BIG;

    for (k = 0; k < seq1len; k++)
	modsdx[k] = -1;
    for (k = 0; k < seq2len; k++)
	modsdx2[k] = -1;
    for (k = 1; k <= tplt->length; k++)
	if (tplt->posn_a[k] > 0 && tplt->posn_b[k] > 0)
	{
	    modsdx[tplt->posn_a[k] - 1] = tplt->posn_b[k] - 1;
	    modsdx2[tplt->posn_b[k] - 1] = tplt->posn_a[k] - 1;
	}

    if (prtalnflg)
    {
	if (htmlflg)
	    htmlalign(tplt);

	prtalign(tplt);
    }

    e_final(epair, esolv, perco_a, perco_b);

    if (verbflg)
	printf("Epair = %f, Esolv = %f\n", *epair * -RTCONST / SCALEFAC, *esolv * -RTCONST / SCALEFAC);
}

/* Change some DSSP Gs to Hs */
void
                helix310(char *struc, int len)
{
    int             i, j, run;

    for (i = 0; i < len; i++)
	if (struc[i] == 'G')
	{
	    for (run = 1; run <= len - i && struc[i + run] == 'G'; run++);
	    if (run >= 4 || (run == 3 && ((i && struc[i - 1] == 'H') || (i < len-run-3 && struc[i + run] == 'H'))))
	    {
		for (j = i; j < i + run; j++)
		    struc[j] = 'H';
		i += run - 1;
	    }
	}
}

/* Smooth DSSP secondary structure definitions */
void
                smooth_sst(void)
{
    int             i;

    for (i = 1; i < seq1len - 1; i++)
	if (tsstruc[i - 1] == tsstruc[i + 1])
	    tsstruc[i] = tsstruc[i - 1];
}

/* Extend DSSP secondary structure definitions */
void
                ext_sst(void)
{
    int             i;

    for (i = 1; i < seq1len - 1; i++)
	if (!tsstruc[i] && !tsstruc[i + 1] && tsstruc[i - 1])
	{
	    tsstruc[i] = tsstruc[i - 1];
	    i++;
	}
	else if (!tsstruc[i - 1] && !tsstruc[i] && tsstruc[i + 1])
	    tsstruc[i] = tsstruc[i + 1];
}

/* Locate secondary structures */
void
                loc_sst(void)
{
    int             i, istart, l;

    nsst = 0;
    for (i = 0; i < seq1len; i++)
	if (tsstruc[i])
	{
	    istart = i;
	    while (i < seq1len && tsstruc[i] == tsstruc[istart])
		i++;
	    l = i - istart;
	    if (l < 2)
		continue;
	    sstlist[nsst].start = istart;
	    sstlist[nsst].length = l;
	    switch (tsstruc[istart])
	    {
	    case HELIX:
		sstlist[nsst].type = HELIX;
		if (verbflg)
		    printf("HELIX  at %3d, length %2d.\n", istart, l);
		break;
	    case STRAND:
		sstlist[nsst].type = STRAND;
		if (verbflg)
		    printf("STRAND at %3d, length %2d.\n", istart, l);
		break;
	    default:
		fail("Unknown secondary structure code!");
	    }
	    sstlen += l;
	    nsst++;
	}
    if (verbflg)
	putchar('\n');

    for (i = 0; i < seq1len; i++)
	sstidx[i] = -1;
    for (i = 0; i < nsst; i++)
	for (l = 0; l < sstlist[i].length; l++)
	    sstidx[sstlist[i].start + l] = i;
}

void
                readxyz(char *buf, float *x, float *y, float *z)
{
    char            temp[9];

    temp[8] = '\0';
    strncpy(temp, buf, 8);
    *x = atof(temp);
    strncpy(temp, buf + 8, 8);
    *y = atof(temp);
    strncpy(temp, buf + 16, 8);
    *z = atof(temp);
}


/* Build structural template */
void
                maketemplate(FILE * dfp, int distflg)
{
    char            buf[512], whichatm[MAXSEQLEN], *cp, ri[10];
    int             i, j, k, nres, natoms, n = 0,
                    consv;
    float           dv, cca[3], nca[3], xx[3], yy[3], sx, sy, rmean = ZERO, sum, targf[20];
    short           at1, at2;

    seq1len = i = natoms = 0;

    if (!fgets(buf, 512, dfp))
	return;

    while (!feof(dfp))
    {
	if (fread(buf, 1, 302, dfp) != 302)
	    break;

	if (distflg)
	    if (sscanf(buf + 12, "%f%f%f", phi + i, psi + i, omega + i) != 3 ||
		sscanf(buf + 39, "%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f",
		       &x[i][NATOM], &y[i][NATOM], &z[i][NATOM],
		       &x[i][CAATOM], &y[i][CAATOM], &z[i][CAATOM],
		       &x[i][CATOM], &y[i][CATOM], &z[i][CATOM],
		       &x[i][OATOM], &y[i][OATOM], &z[i][OATOM],
		       &x[i][CBATOM], &y[i][CBATOM], &z[i][CBATOM]) != 15)
		fail("Bad tdb file!");

	sscanf(buf+174, "%s", ri);

	if (resid[i] != NULL)
	    free(resid[i]);

	resid[i] = strdup(ri);

	cp = buf + 181;
	for (j = 0; j < 20; j++,cp+=6)
	{
	  consv = atoi(cp);
	  tpltsc1[i][j] = consv;
	}

	if (buf[5] == '-')
	{
	    seq1[i] = GAP;
	    whichatm[i] = 0;
	    tsstruc[i] = COIL;
	}
	else
	{
	    seq1[i] = aanum(buf[5]);
	    tsstruc[i] = buf[7];
	    sscanf(buf + 9, "%d", &j);
	    trelacc[i] = j;
	    whichatm[i] = 31;
	}
	++i;
    }

    seq1len = nres = i;

    helix310(tsstruc, nres);
    for (i = 0; i < nres; i++)
	switch (tsstruc[i])
	{
	case 'H':
	    tsstruc[i] = HELIX;
	    break;
	case 'E':
	case 'A':
	case 'P':
	    tsstruc[i] = STRAND;
	    break;
	default:
	    tsstruc[i] = COIL;
	    break;
	}

    for (i = 0; i < nres; i++)
    {
	for (k=0; k<20; k++)
	    targf[k] = dbaaf[k] * exp(0.00318 * tpltsc1[i][k]);

	for (sum=k=0; k<20; k++)
	    sum += targf[k];

	for (k=0; k<20; k++)
	    targf1[i][k] = targf[k] / sum;
    }

    if (!distflg)
	return;

    /* Check atoms */
    for (i = 0; i < nres; i++)
    {
	if (seq1[i] != GAP && !(whichatm[i] & (1 << NATOM)))
	{
	    printf("FATAL: Missing N atom in %d!\n", i + 1);
	    exit(1);
	}
	if (!(whichatm[i] & (1 << CAATOM)))
	{
	    if (verbflg)
		printf("WARNING: Missing CA atom in %d!\n", i + 1);
	    seq1[i] = GAP;
	}
	if (!(whichatm[i] & (1 << CBATOM)))
	{
	    if (!(whichatm[i] & (1 << CAATOM)) || !(whichatm[i] & (1 << CATOM)) || !(whichatm[i] & (1 << NATOM)))
	    {
		/* Not much left of this residue! */
		if (verbflg)
		    printf("WARNING: Missing main chain atom in %d!\n", i + 1);
		seq1[i] = GAP;
		continue;
	    }

	    /* Reconstruct CB atom */
	    nca[0] = x[i][CAATOM] - x[i][NATOM];
	    nca[1] = y[i][CAATOM] - y[i][NATOM];
	    nca[2] = z[i][CAATOM] - z[i][NATOM];
	    cca[0] = x[i][CAATOM] - x[i][CATOM];
	    cca[1] = y[i][CAATOM] - y[i][CATOM];
	    cca[2] = z[i][CAATOM] - z[i][CATOM];
	    for (k = 0; k < 3; k++)
		xx[k] = nca[k] + cca[k];
	    vecprod(yy, nca, cca);
	    sx = CACBDIST * cos(TETH_ANG) / sqrt(dotprod(xx, xx));
	    sy = CACBDIST * sin(TETH_ANG) / sqrt(dotprod(yy, yy));
	    x[i][CBATOM] = x[i][CAATOM] + xx[0] * sx + yy[0] * sy;
	    y[i][CBATOM] = y[i][CAATOM] + xx[1] * sx + yy[1] * sy;
	    z[i][CBATOM] = z[i][CAATOM] + xx[2] * sx + yy[2] * sy;
	    whichatm[i] |= 1 << CBATOM;
	    if (seq1[i] != GLY && verbflg)
		fprintf(stderr, "WARNING: dummy CB atom constructed for %d (%s)!\n", i + 1, rnames[seq1[i]]);
	}
    }

    if (!nres)
	return;

    distmat = allocmat(nres, nres, sizeof(**distmat), FALSE);

    /* Calculate interatomic distance template */
    for (i = 0; i < nres; i++)
	tooi[i] = 0;

    rmax = ZERO;

    for (i = 0; i < nres; i++)
	for (j = i; j < nres; j++)
	    for (at1 = CAATOM; at1 <= NATOM; at1++)
		for (at2 = CAATOM; at2 <= NATOM; at2++)
		    if (i != j || at1 != at2)
			if ((whichatm[i] & (1 << at1)) && (whichatm[j] & (1 << at2)))
			{
			    dv = sqrtf(SQR(x[i][at1] - x[j][at2]) + SQR(y[i][at1] - y[j][at2]) + SQR(z[i][at1] - z[j][at2]));
			    if (j - i > TOPOMAX && at1 == CBATOM && at2 == CBATOM)
			    {
				if (dv > rmax)
				    rmax = dv;
				rmean += dv;
				n++;
			    }
			    distmat[i][j][at1][at2] = distmat[j][i][at2][at1] = dv;
			    if (j - i > 1 && at1 == CBATOM && at2 == CBATOM && dv <= 9.5F)
			    {
				tooi[i]++;
				tooi[j]++;
			    }
			}
			else
			    distmat[i][j][at1][at2] = distmat[j][i][at2][at1] = 999.0;
		    else
			distmat[i][j][at1][at2] = ZERO;

    rmean /= n;

    for (i = 0; i < seq1len; i++)
	coreflg[i] = tsstruc[i] != COIL;

    for (i = 0; i < seq1len; i++)
	if (coreflg[i] && (!i || !coreflg[i - 1]) && (i == seq1len - 1 || !coreflg[i + 1]))
	    coreflg[i] = 0;
}

/* Write PDB file */
void            writepdb(char *fname, char *brkid)
{
    FILE           *ofp;
    int             atnum, i;
    float coorderr;

    ofp = fopen(fname, "w");
    if (ofp != NULL)
    {
	fprintf(ofp, "HEADER  %s BASED GenTHREADER MODEL\n", brkid);
	for (atnum = i = 0; i < seq2len; i++)
	    if (modsdx2[i] < 0 || seq1[modsdx2[i]] == GAP)
	    {
		fprintf(ofp, "REMARK  UNALIGNED ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  0.00  9.99\n",
		   ++atnum, " N  ", rnames[seq2[i]], i + 1, ZERO, ZERO, ZERO);
		fprintf(ofp, "REMARK  UNALIGNED ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  0.00  9.99\n",
		   ++atnum, " CA ", rnames[seq2[i]], i + 1, ZERO, ZERO, ZERO);
		fprintf(ofp, "REMARK  UNALIGNED ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  0.00  9.99\n",
		   ++atnum, " C  ", rnames[seq2[i]], i + 1, ZERO, ZERO, ZERO);
		fprintf(ofp, "REMARK  UNALIGNED ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  0.00  9.99\n",
		   ++atnum, " O  ", rnames[seq2[i]], i + 1, ZERO, ZERO, ZERO);
		if (seq2[i] != GLY)
		    fprintf(ofp, "REMARK  UNALIGNED ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  0.00  9.99\n",
			    ++atnum, " CB ", rnames[seq2[i]], i + 1, ZERO, ZERO, ZERO);
	    }
	    else
	    {
		if (seq2[i] < 20)
		    coorderr = -0.003852 * tpltsc1[modsdx2[i]][seq2[i]] + 0.000858 * ooi_de[seq2[i]][MAX(0, MIN(OOIMAX, tooi[modsdx2[i]]) - OOIMIN) / OOIDIV] + 0.000301 * e_contrib[modsdx2[i]] + 4.83;
		else
		    coorderr = 9.9;
		if (coorderr < 1.0)
		    coorderr = 1.0;
		if (coorderr > 10.0)
		    coorderr = 10.0;
		fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00 %5.2f\n",
			++atnum, " N  ", rnames[seq2[i]], i + 1, x[modsdx2[i]][NATOM], y[modsdx2[i]][NATOM], z[modsdx2[i]][NATOM], coorderr);
		fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00 %5.2f\n",
			++atnum, " CA ", rnames[seq2[i]], i + 1, x[modsdx2[i]][CAATOM], y[modsdx2[i]][CAATOM], z[modsdx2[i]][CAATOM], coorderr);
		fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00 %5.2f\n",
			++atnum, " C  ", rnames[seq2[i]], i + 1, x[modsdx2[i]][CATOM], y[modsdx2[i]][CATOM], z[modsdx2[i]][CATOM], coorderr);
		fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00 %5.2f\n",
			++atnum, " O  ", rnames[seq2[i]], i + 1, x[modsdx2[i]][OATOM], y[modsdx2[i]][OATOM], z[modsdx2[i]][OATOM], coorderr);
		if (seq2[i] != GLY)
		    fprintf(ofp, "ATOM   %4d %s %s  %4d    %8.3f%8.3f%8.3f  1.00 %5.2f\n",
			++atnum, " CB ", rnames[seq2[i]], i + 1, x[modsdx2[i]][CBATOM], y[modsdx2[i]][CBATOM], z[modsdx2[i]][CBATOM], coorderr);
	    }
	fprintf(ofp, "TER\nEND\n");
	fclose(ofp);
    }
}

/* Write profile comparison file */
void            writeprofile(char *fname)
{
    FILE           *ofp;
    int             i, k;

    ofp = fopen(fname, "w");
    if (ofp != NULL)
    {
	for (i = 0; i < seq2len; i++)
	    if (modsdx2[i] >= 0 && seq1[modsdx2[i]] != GAP)
	    {
		fprintf(ofp, "%2d", seq1[modsdx2[i]]);
		for (k=0; k<20; k++)
		    fprintf(ofp, " %d", tpltsc1[modsdx2[i]][k]);
		fprintf(ofp, "\n%2d", seq2[i]);
		for (k=0; k<20; k++)
		    fprintf(ofp, " %d", tpltsc2[i][k]);
		fprintf(ofp, "\n");
	    }

	fclose(ofp);
    }
}

/* Read PSI AA frequency data */
int             getpsi(FILE * lfil)
{
    int             i, j, k, naa, transtab[20];
    float           sum, targf[20];
    char            buf[256];

    if (!fgets(buf, 256, lfil))
	fail("Bad PSI format!");

    for (j = 0, i = 0; i < 80 && buf[i] && j < 20; i++)
	if (isalpha(buf[i]))
	    transtab[j++] = aanum(buf[i]);

    if (j != 20)
	fail("Bad PSI format!");

    naa = 0;
    while (!feof(lfil))
    {
	if (!fgets(buf, 256, lfil))
	    break;
	if (sscanf(buf, "%d", &i) != 1)
	    break;
	if (i != naa + 1)
	    fail("Bad PSI format!");
	seq2[naa] = buf[6];
	if (sscanf(buf + 7, "%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d", &tpltsc2[naa][transtab[0]], &tpltsc2[naa][transtab[1]], &tpltsc2[naa][transtab[2]], &tpltsc2[naa][transtab[3]], &tpltsc2[naa][transtab[4]], &tpltsc2[naa][transtab[5]], &tpltsc2[naa][transtab[6]], &tpltsc2[naa][transtab[7]], &tpltsc2[naa][transtab[8]], &tpltsc2[naa][transtab[9]], &tpltsc2[naa][transtab[10]], &tpltsc2[naa][transtab[11]], &tpltsc2[naa][transtab[12]], &tpltsc2[naa][transtab[13]], &tpltsc2[naa][transtab[14]], &tpltsc2[naa][transtab[15]], &tpltsc2[naa][transtab[16]], &tpltsc2[naa][transtab[17]], &tpltsc2[naa][transtab[18]], &tpltsc2[naa][transtab[19]]) != 20)
	    fail("Bad PSI format!");
	naa++;
    }

    for (i=0; i<naa; i++)
    {
	for (k=0; k<20; k++)
	    targf[k] = dbaaf[k] * exp(0.318 * tpltsc2[i][k]);

	for (sum=k=0; k<20; k++)
	    sum += targf[k];

	for (k=0; k<20; k++)
	    targf2[i][k] = targf[k] / sum;

	for (k=0; k<20; k++)
	    tpltsc2[i][k] *= 100;
    }

    psidata = TRUE;

    return naa;
}

/* Read PSI mtx PSSM data */
int             getmtx(FILE *lfil)
{
    int             i, j, k, naa;
    float           sum, targf[20];
    char            buf[256];

    rewind(lfil);

    if (fscanf(lfil, "%d", &naa) != 1)
	fail("Bad mtx file - sequence length not found!");

    if (naa > MAXSEQLEN)
	fail("Input sequence too long!");

    if (fscanf(lfil, "%s", seq2) != 1)
      fail("Bad mtx file - no sequence!");

    while (!feof(lfil))
    {
	if (!fgets(buf, 65536, lfil))
	    fail("Bad mtx file!");
	if (!strncmp(buf, "-32768 ", 7))
	{
	    for (j=0; j<naa; j++)
	    {
		if (sscanf(buf, "%*d%d%*d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%*d%d", &tpltsc2[j][ALA],  &tpltsc2[j][CYS], &tpltsc2[j][ASP],  &tpltsc2[j][GLU],  &tpltsc2[j][PHE],  &tpltsc2[j][GLY],  &tpltsc2[j][HIS],  &tpltsc2[j][ILE],  &tpltsc2[j][LYS],  &tpltsc2[j][LEU],  &tpltsc2[j][MET],  &tpltsc2[j][ASN],  &tpltsc2[j][PRO],  &tpltsc2[j][GLN],  &tpltsc2[j][ARG],  &tpltsc2[j][SER],  &tpltsc2[j][THR],  &tpltsc2[j][VAL],  &tpltsc2[j][TRP],  &tpltsc2[j][TYR]) != 20)
		    fail("Bad mtx format!");
		if (!fgets(buf, 65536, lfil))
		    break;
	    }
	}
    }

    for (i=0; i<naa; i++)
    {
	for (k=0; k<20; k++)
	    targf[k] = dbaaf[k] * exp(0.00318 * tpltsc2[i][k]);

	for (sum=k=0; k<20; k++)
	    sum += targf[k];

	for (k=0; k<20; k++)
	    targf2[i][k] = targf[k] / sum;
    }

    psidata = TRUE;

    return naa;
}


/* Read PSIPRED VFORMAT prediction data */
int             getpsipredv(FILE * lfil)
{
    int             naa;
    float confc, confh, confe;
    char            buf[256];

    if (!fgets(buf, 256, lfil))
	fail("Bad PSIPRED VFORMAT file!");

    naa = 0;
    while (!feof(lfil))
    {
	if (!fgets(buf, 256, lfil))
	    break;
	if (sscanf(buf+10, "%f%f%f", &confc, &confh, &confe) != 3)
	    break;
	seq2[naa] = buf[5];
	switch (buf[7])
	{
	case 'H':
	    ssstruc[naa] = HELIX;
	    break;
	case 'E':
	    ssstruc[naa] = STRAND;
	    break;
	default:
	    ssstruc[naa] = COIL;
	    break;
	}
	ssstrel[naa++] = 100*(2*MAX(MAX(confc, confh),confe)-(confc+confh+confe)-MIN(MIN(confc, confh),confe));
    }

    ssstflg = TRUE;

    return naa;
}

/* Read PSIPRED HFORMAT prediction data */
int             getpsipredh(FILE * lfil)
{
    int             naa, nsst, nrel;
    char            buf[256], *p;

    rewind(lfil);

    while (!feof(lfil))
    {
	if (!fgets(buf, 256, lfil))
	    fail("Bad PSIPRED format!");
	if (!strncmp(buf, "Conf:", 5))
	    break;
    }

    nsst = nrel = naa = 0;
    while (!feof(lfil))
    {
	p = buf+5;
	while (*++p)
	{
	    if (isdigit(*p))
		ssstrel[nrel++] = (*p - '0')*10;
	}
	if (!fgets(buf, 256, lfil))
	    fail("Bad PSIPRED format!");
	p = buf+5;
	while (*++p)
	    if (isalpha(*p))
		switch (*p)
		{
		case 'H':
		    ssstruc[nsst++] = HELIX;
		    break;
		case 'E':
		    ssstruc[nsst++] = STRAND;
		    break;
		default:
		    ssstruc[nsst++] = COIL;
		    break;
		}

	if (!fgets(buf, 256, lfil))
	    fail("Bad PSIPRED format!");
	p = buf+5;
	while (*++p)
	{
	    if (isalpha(*p))
		seq2[naa++] = *p;
	    else if (isdigit(*p))
		seq2[naa++] = -(*p - '0' + 1);
	}
	if (!fgets(buf, 256, lfil))
	    fail("Bad PSIPRED format!");
	while (!feof(lfil))
	{
	    if (!fgets(buf, 256, lfil))
		break;
	    if (!strncmp(buf, "Conf:", 5))
		break;
	}
    }

    if (nsst != nrel || naa != nsst)
    {
	fprintf(stderr, "Incorrect PSIPRED prediction file format!\n");
	return -1;
    }

    ssstflg = TRUE;

    if (verbflg)
	puts("Parsed PSIPRED output.");

    return naa;
}

/*
 * This routine will read in one sequence from a database file. The sequence
 * can be in any of the supported formats. Returns length of sequence.
 */
int
                getseq(char *dbname, char *dseq, FILE * lfil)
{
    int             i, j, len;
    short           badln, fformat = -1;
    enum
    {
	unknown, embl, genbank, staden, fastp, codata, owl, intelgen, gcg
    };
    char            temp[8192], split;
    int             offset;

    offset = j = 0;

    if (fgets(temp, 8192, lfil) && !strncmp(temp, "# PSIPRED H", 11))
	return getpsipredh(lfil);

    if (!strncmp(temp, "# PSIPRED V", 11))
	return getpsipredv(lfil);

    rewind(lfil);

    if (fgets(temp, 8192, lfil) && !strncmp(temp, "Conf:", 5))
	return getpsipredh(lfil);

    if (fgets(temp, 8192, lfil) && !strncmp(temp, "Conf:", 5))
	return getpsipredh(lfil);

    rewind(lfil);

    if (fgets(temp, 8192, lfil) && isdigit(temp[0]))
	return getmtx(lfil);

    if (fgets(temp, 8192, lfil) && !strncmp(temp, "Last position-specific", 22))
	return getpsi(lfil);

    rewind(lfil);

    if (fgets(temp, 8192, lfil) && fgets(temp, 8192, lfil))
    {
	if (!strncmp(temp + 14, "999.999", 7))
	{
	    ssstflg = TRUE;
	    do
	    {
		dseq[j] = temp[5];
		switch (temp[7])
		{
		case 'H':
		    ssstruc[j++] = HELIX;
		    break;
		case 'E':
		case 'A':
		case 'P':
		    ssstruc[j++] = STRAND;
		    break;
		default:
		    ssstruc[j++] = COIL;
		    break;
		}
	    }
	    while (fgets(temp, 8192, lfil));
	    dseq[j] = '\0';
	    return j;
	}
    }

    rewind(lfil);

    if (!fgets(temp, 8192, lfil))
	return (-1);

    /* Look for old-style PSI file */
    if (!strncmp(temp, "# PSI", 5))
	return getpsi(lfil);

    if (strstr(temp, "of:") != NULL && strstr(temp, "check:") != NULL)
	fformat = gcg;
    else if ((temp[0] == '<') && (temp[19] == '>'))
	fformat = staden;
    else if (strncmp(temp, "ID   ", 5) == 0)
	fformat = embl;
    else if (strncmp(temp, "LOCUS     ", 10) == 0)
	fformat = genbank;
    else if (strncmp(temp, "ENTRY", 5) == 0)
	fformat = codata;
    else if (temp[0] == ';')
	fformat = intelgen;
    else if (temp[0] == '>' && (temp[1] == '>' || temp[3] == ';'))
	fformat = owl;
    else if (temp[0] == '>')
	fformat = fastp;
    else
	fformat = unknown;

    switch (fformat)
    {
    case gcg:
	sscanf(strstr(temp, "of:") + 3, "%s", dbname);
	while (strstr(temp, "..") == NULL)
	    fgets(temp, 8192, lfil);
	fgets(temp, 8192, lfil);
	break;
    case embl:
	strncpy(dbname, temp + 5, 70);
	while (temp[0] != ' ')
	    fgets(temp, 8192, lfil);
	break;

    case genbank:
	while (strncmp(temp, "ORIGIN", 6) != 0)
	{
	    fgets(temp, 8192, lfil);
	    if (strncmp(temp, "DEFINITION", 10) == 0)
		strncpy(dbname, temp + 12, 70);
	}
	fgets(temp, 8192, lfil);
	break;

    case codata:
	strncpy(dbname, temp + 6, 70);
	while (strncmp(temp, "SEQUENCE", 8) != 0)
	    fgets(temp, 8192, lfil);
	fgets(temp, 8192, lfil);
	break;

    case owl:
	fgets(temp, 8192, lfil);
	strncpy(dbname, temp, 70);
	fgets(temp, 8192, lfil);
	break;

    case fastp:
	strncpy(dbname, temp + 1, 70);
	fgets(temp, 8192, lfil);
	break;

    case staden:
	strncpy(dbname, temp + 1, 18);
	offset = 20;
	break;

    case intelgen:
	while (*temp == ';')
	    fgets(temp, 8192, lfil);
	fgets(temp, 8192, lfil);
	break;

    default:
	do
	{
	    len = strlen(temp);
	    for (badln = i = 0; i < len; i++)
		if (islower(temp[i]) || temp[i] == 'J' || temp[i] == 'O' || temp[i] == 'U')
		{
		    badln = TRUE;
		    break;
		}
	    if (badln && !fgets(temp, 8192, lfil))
		return (-1);
	}
	while (badln);
	strcpy(dbname, "<NO NAME>");
	break;
    }

    if (dbname[(len = strlen(dbname)) - 1] == '\n')
	dbname[--len] = '\0';
    if (len >= 70)
	dbname[70] = '\0';

    for (;;)
    {
	if (!strncmp(temp, "//", 2))
	    break;
	len = strlen(temp);
	for (i = offset; i < len && j < MAXSEQLEN; i++)
	{
	    split = islower(temp[i]) ? toupper(temp[i]) : temp[i];
	    if (split == '@' || (fformat == owl && split == '*'))
	    {
		dseq[j] = '\0';
		while (fgets(temp, 8192, lfil));
		return (j);
	    }
	    if (isalpha(split))
		dseq[j++] = split;
	    else if (temp[i] == '\n')
		break;
	}
	if (staden)
	    offset = 0;
	if (!fgets(temp, 8192, lfil))
	    break;
    }

    if (j == MAXSEQLEN)
	fprintf(stderr, "\nWARNING: sequence %s over %d long; truncated!\n",
	       dbname, MAXSEQLEN);

    dseq[j] = '\0';
    return (j);
}

void
                usage(char *cmdname)
{
    printf("usage: %s {options} [sequence-file output-file [list-file]]\n", cmdname);
    puts("\nwhere options are as follows:");
    puts("-cnnn = set pair potential cutoff distance (default: 10.0)");
    puts("-rnnn = set number of sequence shuffles (default: 0)");
    puts("-ynnn = set secondary structure score weighting (default: 100)");
    puts("-hnnn = set minimum reliability for SST filter (default: 0.5)");
    puts("-Xa,b,c = set CHE gap opening penalties (default: 11)");
    puts("-Od,e,f = set CHE gap extension penalties (default: 1)");
    puts("-Ffnm = read secondary structure data from file fnm");
    puts("-Hfnm = output HTML alignment to file fnm");
    puts("-S = include predicted secondary structure in alignment output");
    puts("-u = use contact number solvation scoring");
    puts("-m{xxx} = generate PDB files (xxx optional prefix)");
    puts("-p{q} = print alignments (q selects quick align mode)");
    puts("-v = select Verbose Mode");
    puts("\nExample: threader -d100 -c12 -v -p globin.seq output.fil");
    exit(1);
}

int             main(int argc, char **argv)
{
    int             aa, i, j, k, l, nn = 0, domflg=FALSE, filtered, rfrom = 0, rto = 0, start1, start2, end1, end2, alnmode1, alnmode2, nid;
    char           *cmdstr, *cp;
    char            desc[512], wpname[80];
    float           tpair, tsolv, psum, psumsq, ssum, ssumsq, smax, pmax,
                    nwsc, nwsc1, nwsc2, nwsc3, nwsc4, dum, perid;
    float           pair_av, solv_av, tot_av, pair_sd, solv_sd, tot_sd;
    FILE           *ifp, *ofp, *tfp, *ssfp;
    char            tdbname[512], buf[MAXSEQLEN], templname[512];
    ALNS            tplt;

    struct rlimit limits;

    if (verbflg)
    {
	printf("GenTHREADER - Genomic Protein Sequence Threading Program\n");
	printf("Build date : %s\n",__DATE__);
	printf("Program written by David T. Jones\n");
	printf("Copyright (C) 2001 University College London\n");
	printf("Portions Copyright (C) 1997 University of Warwick\n");
	printf("Portions Copyright (C) 1989,1991 David T. Jones\n\n");
    }

    getrlimit(RLIMIT_STACK, &limits);

    if (limits.rlim_cur < 7340032)
        fail("Stack too small! (type 'limit stacksize 7m' at command prompt)");

    for (cmdstr = *argv++, argc--; argc && **argv == '-'; argv++, argc--)
	switch (*(*argv + 1))
	{
	case 'r':
	    SHUFCOUNT = atoi(*argv + 2);
	    break;
	case 'h':
	    MINSSTREL = 100.0*atof(*argv + 2);
	    if (MINSSTREL > 100)
	      fail("Max. secondary structure reliability should be <= 1.0!");
	    break;
	case 'B':
	    cp = *argv + 2;
	    for (i=0; i<20 && cp; i++)
		SC_PARAM[i] = atof(strsep(&cp, "/:;,"));
	    break;
	case 'O':
	    cp = *argv + 2;
	    for (i=0; i<3 && cp; i++)
		GP_OPEN[i] = 0.5 + 100.0 * atof(strsep(&cp, "/:;,"));
	    break;
	case 'X':
	    cp = *argv + 2;
	    for (i=0; i<3 && cp; i++)
		GP_EXT[i] = 0.5 + 100.0 * atof(strsep(&cp, "/:;,"));
	    break;
	case 'C':
	    NWCUTOFF = atoi(*argv + 2);
	    break;
	case 'S':
	    seqssflg = TRUE;
	    break;
	case 'g':
	    filtgapflg = TRUE;
	    break;
	case 'y':
	    ssscore = atoi(*argv + 2);
	    break;
	case 'c':
	    DCUTOFF = atof(*argv + 2);
	    break;
	case 'm':
	    motifflg = TRUE;
	    strcpy(motpref, *argv + 2);
	    break;
	case 'p':
	    prtalnflg = TRUE;
	    if (argv[0][2] == 'q')
		quickalnflg = TRUE;
	    if (argv[0][2] == 'm')
		mod3flg = TRUE;
	    break;
	case 'H':
	    htmname = strdup(*argv + 2);
	    htmlflg = TRUE;
	    break;
	case 'F':
	    ssfname = strdup(*argv + 2);
	    ssstflg = TRUE;
	    break;
	case 'u':
	    ooisolv = TRUE;
	    break;
        case 'l':
	    localmode = TRUE;
	    break;
	case 'v':
	    verbflg = TRUE;
	    break;
	default:
	    usage(cmdstr);
	}

    if (argc < 2)
	usage(cmdstr);

    randomise();

    ifp = fopen(argv[0], "r");
    if (!ifp)
	fail("Unable to open seq file!");
    seq2len = getseq(desc, seq2, ifp);

    if (ssstflg)
    {
	ssfp = fopen(ssfname, "r");
	if (!ssfp)
	    fail("Cannot open PSIPRED file!");
	ssstflg = FALSE;
	if (fgets(buf, 160, ssfp) && !strncmp(buf, "# PSIPRED H", 11))
	    if (getpsipredh(ssfp) != seq2len)
		fail("Length mismatch between sequence and PSIPRED file!");
	if (!strncmp(buf, "# PSIPRED V", 11))
	    if (getpsipredv(ssfp) != seq2len)
		fail("Length mismatch between sequence and PSIPRED file!");
	if (!strncmp(buf, "Conf:", 5))
	    if (getpsipredh(ssfp) != seq2len)
		fail("Length mismatch between sequence and PSIPRED file!");
	if (fgets(buf, 160, ssfp) && !strncmp(buf, "Conf:", 5))
	    if (getpsipredh(ssfp) != seq2len)
		fail("Length mismatch between sequence and PSIPRED file!");
	fclose(ssfp);
	if (!ssstflg)
	    fail("Cannot parse PSIPRED file!");
    }

    fclose(ifp);

    if (seq2len > MAXSEQLEN - 1)
	fail("Sequence too long!");

    if (seq2len < 1)
	fail("No sequence could be extracted! (the input file may not be correctly formatted)");

    if (verbflg)
	printf("\n%d residues read:\n%s\n\n", seq2len, seq2);

    for (i = 0; i < seq2len; i++)
    {
	aa = aanum(seq2[i]);
	switch (aa)
	{
	case 20:
	    /* ASP is more common than ASN! */
	    seq2[i] = ASP;
	    break;
	case 21:
	    /* GLU is more common than GLN! */
	    seq2[i] = GLU;
	    break;
	case 22:
	    seq2[i] = UNK;
	    break;
	default:
	    seq2[i] = aa;
	    break;
	}
    }

    for (filtered = i = 0; i < seq2len - 3; i++)
	if (seq2[i] == UNK && seq2[i + 1] == UNK && seq2[i + 2] == UNK)
	{
	    filtered = TRUE;
	    break;
	}

    if (!quickalnflg)
    {
	if (verbflg)
	    puts("\nReading mean force potential tables...");
	if (potinit())
	    fail("Unable to load potentials!");
    }

    ofp = fopen(argv[1], "w");
    if (!ofp)
	fail("Unable to open output file!");


    pat = allocmat(MAXSEQLEN + 1, MAXSEQLEN + 1, sizeof(int), FALSE);

    if (argc > 2)
	ifp = fopen(argv[2], "r");
    else if ((ifp = fopen("psichain.lst", "r")) == NULL && getenv("THREAD_DIR"))
    {
	strcpy(templname, getenv("THREAD_DIR"));
	if (templname[strlen(templname) - 1] != '/')
	    strcat(templname, "/");
	strcat(templname, "psichain.lst");
	ifp = fopen(templname, "r");
    }

    if (!ifp)
	fail("Cannot open protein list (default: psichain.lst)!");

    while (!feof(ifp))
    {
        if (!fgets(buf, 256, ifp))
            break;
        cp = buf+strlen(buf);


        while (--cp > buf)
            if (*cp == ' ')
                break;


	if (sscanf(cp, "%s", brkid) != 1)
	    break;

 	  while(*cp == ' ')
            --cp;

	  while(--cp > buf)
	    if (*cp == ' ')
	     break;
	  sscanf(cp, "%d", &rto);

	  while(*cp == ' ')
           --cp;

          while(--cp > buf)
            if (*cp == ' ')
	     break;
	  sscanf(cp, "%d", &rfrom);

//	  fprintf(stderr, "CHECK %s %d %d %s\n", brkid, rfrom, rto, cp);

	if (brkid[0] == '#')
	    continue;
	if (verbflg)
	    printf("Generating template for : %s ...\n", brkid);
	/* Read coords from PDB file */
	if (getenv("TDB_DIR"))
	{
	    strcpy(tdbname, getenv("TDB_DIR"));
	    if (tdbname[strlen(tdbname) - 1] != '/')
		strcat(tdbname, "/");
	}
	else if (getenv("THREAD_DIR"))
	{
	    strcpy(tdbname, getenv("THREAD_DIR"));
	    if (tdbname[strlen(tdbname) - 1] != '/')
		strcat(tdbname, "/");
	}
	else
	  strcpy(tdbname, "./tdb/");

	strcat(tdbname, brkid);
	strcat(tdbname, ".tdb");

	tfp = fopen(tdbname, "r");
	if (tfp == NULL)
	{
	    fprintf(stderr, "*** main: Cannot open TDB file (%s)!\n", tdbname);
	    fflush(stderr);
	    continue;
	}

	maketemplate(tfp, FALSE);

	if (!seq1len)
	{
	    fprintf(stderr, "*** main: %s - length error!\n", brkid);
	    fflush(stderr);
	    continue;
	}
        int last=seq2len;
        int naln=0;
	while( naln < NALN )
	{
          //truncate seq2 and proceed
          int x=0, s2from = 0, s2len=0;
          int s2to   = last;
          char frag[10000];

          i=s2from;

          while(i < s2to)
	  {
	    frag[x]=seq2[i];
            i++; x++;
	  }
          s2len=(s2to-s2from);

	  if (psidata)
	  {
 	    nwsc = seqscore(seq1, frag, NULL, 1, seq1len, s2len);
 	  }
          else
	  {
	    nwsc = seqscore(seq1, frag, NULL, 0, seq1len, s2len);
	  }

	  if (verbflg)
	      fprintf(stderr,"NWSC = %f\n", nwsc);

	  if (quickalnflg || nwsc >= NWCUTOFF)
	  {
	    if (psidata)
	    {
		nwsc = seqscore(seq1, frag, &tplt, 1, seq1len, s2len);
            }
	    else
		nwsc = seqscore(seq1, frag, &tplt, 0, seq1len, s2len);

	    if (quickalnflg)
	    {
		if (htmlflg)
		    htmlalign(&tplt);
  		   prtalign(&tplt);
		exit(0);
	    }

	    nid = l = 0;
	    start1 = start2 = BIG;
	    end1 = end2 = 0;
            int fst=1;


	    for (k = 1; k <= tplt.length; k++)
	    {


		if (tplt.posn_a[k] > 0 && tplt.posn_b[k] > 0)
		{
		    if (seq1[tplt.posn_a[k]-1] == seq2[tplt.posn_b[k]-1])
			nid++;

                    if (fst)
                    {
		     start1 = tplt.posn_a[k];
                     start2 = tplt.posn_b[k];
                     fst=0;
                    }
		    end1   = tplt.posn_a[k];
		    end2   = tplt.posn_b[k];
		    l++;
		}
	    }


	    perid = 100.0F * nid / MIN(end1-start1, end2-start2);

	    if (perid > 40.0F)
		nwsc += perid;

	    if (! htmlflg)
               strcpy(hits[nn].brkid,brkid);


	    // fprintf(stdout,"hits brkid %s start %d end %d %d\n" ,brkid,start2,end2,nn);

            if(! prtalnflg || ( start2 == rfrom && end2 == rto ))
	    {
	      //  fprintf(stdout,"hits brkid %s start %d end %d %d\n" ,brkid,start2,end2,nn);
	     hits[nn].laln = l;
	     hits[nn].lena = seq1len;
	     hits[nn].lenb = seq2len;
             hits[nn].sfrom= start2;
             hits[nn].sto  = end2;
	     hits[nn].nwsc = nwsc;

             domflg=TRUE;
	    }

	    if (verbflg)
	      printf("ALNLEN %s %4d %4d %4d %4d %4d %4d %4d\n", brkid, start1, end1, start2, end2, l,rfrom,rto);
	  }
          last=start2;


	  if ( ( !prtalnflg &&  nwsc >= NWCUTOFF && end2-start2 >=10 ) || ( start2==rfrom && end2==rto ) )
  	  {
	    //find another hit in same sequence
	    //fprintf(stdout,"REWIND hits brkid %s start %d end %d %d\n" ,brkid,start2,end2,nn);
             rewind(tfp);
             naln++;
          }
          else
	  {
	    if( ! prtalnflg )
	    {
	      //         fprintf(stdout,"FINISHED hits brkid %s start %d end %d %d\n" ,brkid,start2,end2,nn);
	         naln=NALN;
                 break;
            }
            else
	        continue;
          }
	  //fprintf(stdout,"TEMPLATE hits brkid %s start %d end %d %d\n" ,brkid,start2,end2,nn);

        maketemplate(tfp, TRUE);


	loc_sst();

	pmat = allocmat(seq2len, seq1len, sizeof(int), TRUE);

	if (prtalnflg)
	{
	    printf("\n\n>>> Alignment with %s:\n", brkid);
	}
	//where proper alignment happens and printing

	if (verbflg)
	    printf("Hits brkid %s start %d end %d %d\n" ,brkid,start2,end2,nn);

 	seqfit(&hits[nn].epair, &hits[nn].esolv, &hits[nn].perco_a, &hits[nn].perco_b, &tplt);

	if (mod3flg)
	{
	    printf(">P1;%c%c%c%c\nstructureX:p%c%c%c%c:%s:%c:%s:%c:PDB::0.00:0.00\n", brkid[0], brkid[1], brkid[2], brkid[3], brkid[0], brkid[1], brkid[2], brkid[3], resid[0], brkid[4] == '0' ? ' ' : brkid[4], resid[seq1len - 1], brkid[4] == '0' ? ' ' : brkid[4]);
	    for (i = 1; i <= tplt.length; i++)
		if (tplt.posn_a[i] <= 0)
		    putchar('-');
		else if (seq1[tplt.posn_a[i] - 1] < GAP)
		    putchar(rescodes[seq1[tplt.posn_a[i] - 1]]);
		else
		    putchar('X');
	    puts("*");
	    puts(">P1;SEQ\nsequence:SEQ:::::SEQ::0.00:0.00");
	    for (i = 1; i <= tplt.length; i++)
		if (tplt.posn_b[i] <= 0)
		    putchar('-');
		else if (seq2[tplt.posn_b[i] - 1] < GAP)
		    putchar(rescodes[seq2[tplt.posn_b[i] - 1]]);
		else
		    putchar('X');
	    puts("*");
	}

	if (SHUFCOUNT > 0)
	{
	    if (verbflg)
		puts("\nShuffling sequence...");
	    memcpy(buf, seq2, seq2len);
	    psum = psumsq = ssum = ssumsq = ZERO;
	    pmax = smax = -VBIG;
	    for (i = 0; i < SHUFCOUNT; i++)
	    {
		shufseq(seq2, seq2len);
		seqfit(&tpair, &tsolv, &dum, &dum, &tplt);
		pmax = MAX(tpair, pmax);
		smax = MAX(tsolv, smax);
		psum += tpair;
		psumsq += SQR(tpair);
		ssum += tsolv;
		ssumsq += SQR(tsolv);
	    }
	    memcpy(seq2, buf, seq2len);
	    pair_av = psum / SHUFCOUNT;
	    solv_av = ssum / SHUFCOUNT;
	    pair_sd = sqrt(psumsq / SHUFCOUNT - SQR(pair_av));
	    solv_sd = sqrt(ssumsq / SHUFCOUNT - SQR(solv_av));

	    if (verbflg)
	    {
		printf("Min. pair energy for shuffled sequences: %f\n", -(RTCONST / SCALEFAC) * pmax);
		printf("Mean pair energy for shuffled sequences: %f\n", -(RTCONST / SCALEFAC) * pair_av);
		printf("Match pair energy %f s.d. below mean\n", -(pair_av - hits[nn].epair) / pair_sd);
		printf("Min. solvation energy for shuffled sequences: %f\n", -(RTCONST / SCALEFAC) * smax);
		printf("Mean solvation energy for shuffled sequences: %f\n", -(RTCONST / SCALEFAC) * solv_av);
		printf("Match solvation energy %f s.d. below mean\n", -(solv_av - hits[nn].esolv) / solv_sd);
	    }

	    if (ofp != NULL)
		fprintf(ofp, "%7.1f %7.1f %6.2f %6.1f %6.1f %6.2f %6.2f %6.1f %5.1f %5.1f %d %d %s\n",
			-(RTCONST / SCALEFAC) * hits[nn].epair,
			-(RTCONST / SCALEFAC) * pmax,
			-(pair_av - hits[nn].epair) / pair_sd,
			-(RTCONST / SCALEFAC) * hits[nn].esolv,
			-(RTCONST / SCALEFAC) * smax,
			-(solv_av - hits[nn].esolv) / solv_sd,
			-(pair_av - hits[nn].epair) / pair_sd - (solv_av - hits[nn].esolv) / solv_sd,
			hits[nn].nwsc,
			hits[nn].perco_a,
			hits[nn].perco_b,
                        hits[nn].sfrom,
                        hits[nn].sto,
			hits[nn].brkid);
	    fflush(ofp);
	}

	if (motifflg)
	{
	    strcpy(wpname, motpref);
	    strcat(wpname, "_");
	    strcat(wpname, brkid);
	    strcat(wpname, ".model.pdb");
	    writepdb(wpname, brkid);
	}

	if ( nwsc >= NWCUTOFF  || prtalnflg )
	{
	    nn++;
        }

       	freemat(pmat, seq2len);
        freemat(distmat, seq1len);

	if (verbflg)
	    printf("%s %d %d %f %d\n",brkid,start2,end2,nwsc,nn);

        if(  ( ! prtalnflg && start2 < 10 && nwsc < NWCUTOFF ) || ( prtalnflg  &&  nwsc >= NWCUTOFF ) )
	{
          naln=NALN;
          continue;
	}


        }// while naln++;


	fclose(tfp);

	fflush(stdout);

        continue;
    }//while ifp

    if (SHUFCOUNT)
	return 0;

    pair_av = solv_av = pair_sd = solv_sd = ZERO;

    for (i = 0; i < nn; i++)
    {
	pair_av += hits[i].epair;
	solv_av += hits[i].esolv;
    }

    pair_av /= nn;
    solv_av /= nn;

    for (i = 0; i < nn; i++)
    {
	pair_sd += hits[i].epair - pair_av;
	solv_sd += hits[i].esolv - solv_av;
    }

    pair_sd = -SQR(pair_sd) / nn;
    solv_sd = -SQR(solv_sd) / nn;

    for (i = 0; i < nn; i++)
    {
	pair_sd += SQR(hits[i].epair - pair_av);
	solv_sd += SQR(hits[i].esolv - solv_av);
    }

    pair_sd = sqrt(pair_sd / (nn - 1));
    solv_sd = sqrt(solv_sd / (nn - 1));

    for (tot_av = ZERO, i = 0; i < nn; i++)
	tot_av += hits[i].epair + hits[i].esolv * pair_sd / solv_sd;
    tot_av /= nn;

    for (tot_sd = ZERO, i = 0; i < nn; i++)
	tot_sd += SQR(hits[i].epair + hits[i].esolv * pair_sd / solv_sd - tot_av);

    tot_sd = sqrt(tot_sd / (nn - 1));


    if (ofp != NULL)
	for (i = 0; i < nn; i++)
	{
	    fprintf(ofp, "%7.1f %5.2f %5.2f %7.1f %5.2f %5.2f %7.1f %5.2f %7.1f %5.2f %6.1f %4d %4d %4d %4d %4d %s\n",
		    -(RTCONST / SCALEFAC) * hits[i].epair,
		    -(pair_av - hits[i].epair) / pair_sd,
		    (hits[i].perco_a >= 50.0 && hits[i].perco_b >= 50.0) ? -(pair_av - hits[i].epair) / pair_sd : -9.99,
		    -(RTCONST / SCALEFAC) * hits[i].esolv,
		    -(solv_av - hits[i].esolv) / solv_sd,
		    (hits[i].perco_a >= 50.0 && hits[i].perco_b >= 50.0) ? -(solv_av - hits[i].esolv) / solv_sd : -9.99,
		    -(RTCONST / SCALEFAC) * (hits[i].epair + hits[i].esolv * (pair_sd / solv_sd)),
		    -(tot_av - hits[i].epair - hits[i].esolv * pair_sd / solv_sd) / tot_sd,
		    (hits[i].perco_a >= 50.0 && hits[i].perco_b >= 50.0) ? -(RTCONST / SCALEFAC) * (hits[i].epair + hits[i].esolv * (pair_sd / solv_sd)) : 9999.9999,
		    (hits[i].perco_a >= 50.0 && hits[i].perco_b >= 50.0) ? -(tot_av - hits[i].epair - hits[i].esolv * pair_sd / solv_sd) / tot_sd : -9.99,
		    hits[i].nwsc,
		    hits[i].laln,
		    hits[i].lena,
		    hits[i].lenb,
                    hits[i].sfrom,
                    hits[i].sto,
		    hits[i].brkid);

	    fflush(ofp);
	}

    freemat(pat, MAXSEQLEN+1);

    return 0;
}
