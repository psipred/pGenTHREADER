/* PRINTCONS - Print Residue Conservation */

/* Copyright (C) 2010 David T. Jones */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

#define MAXSEQLEN 5000

int             profile[MAXSEQLEN][20];

int             seqlen;

char seq[MAXSEQLEN];

enum aacodes
{
    ALA, ARG, ASN, ASP, CYS,
    GLN, GLU, GLY, HIS, ILE,
    LEU, LYS, MET, PHE, PRO,
    SER, THR, TRP, TYR, VAL,
    UNK
};

void
err(char *s)
{
    fprintf(stderr, "%s\n", s);
}

void
fail(char *s)
{
    err(s);
    exit(1);
}

/* Convert AA letter to numeric code (0-20) */
int
aanum(ch)
    int             ch;
{
    static const int      aacvs[] =
    {
        999, 0, 20, 4, 3, 6, 13, 7, 8, 9, 20, 11, 10, 12, 2,
        20, 14, 5, 1, 15, 16, 20, 19, 17, 20, 18, 20
    };

    return (isalpha(ch) ? aacvs[ch & 31] : 20);
}

/* Make 1st level prediction averaged over specified weight sets */
void
predict()
{
    int             aa, i, j, cv;

    for (i=0; i<seqlen; i++)
    {
        aa = seq[i];
        if (aa >= 20)
            cv = 0;
        else
            cv = (int)(profile[i][aa]/100.0);
        if (cv > 9)
            cv = 9;
        if (cv < 0)
            cv = 0;
        putchar('0'+cv);
    }

    putchar('\n');
}

/* Read PSI AA frequency data */
int             getmtx(FILE *lfil)
{
    int             aa, i, j, naa;
    char            buf[256], *p;

    if (fscanf(lfil, "%d", &naa) != 1)
        fail("Bad mtx file - no sequence length!");

    if (naa > MAXSEQLEN)
        fail("Input sequence too long!");

    if (fscanf(lfil, "%s", seq) != 1)
        fail("Bad mtx file - no sequence!");

    while (!feof(lfil))
    {
        if (!fgets(buf, 65536, lfil))
            fail("Bad mtx file!");
        if (!strncmp(buf, "-32768 ", 7))
        {
            for (j=0; j<naa; j++)
            {
                if (sscanf(buf, "%*d%d%*d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%*d%d", &profile[j][ALA],  &profile[j][CYS], &profile[j][ASP],  &profile[j][GLU],  &profile[j][PHE],  &profile[j][GLY],  &profile[j][HIS],  &profile[j][ILE],  &profile[j][LYS],  &profile[j][LEU],  &profile[j][MET],  &profile[j][ASN],  &profile[j][PRO],  &profile[j][GLN],  &profile[j][ARG],  &profile[j][SER],  &profile[j][THR],  &profile[j][VAL],  &profile[j][TRP],  &profile[j][TYR]) != 20)
                    fail("Bad mtx format!");
                if (!fgets(buf, 65536, lfil))
                    break;
            }
        }
    }

    return naa;
}

main(int argc, char **argv)
{
    int             i, niters;
    FILE *ifp;

    /* malloc_debug(3); */
    if (argc != 2)
        fail("usage : mtx2cons mtx-file");
    ifp = fopen(argv[1], "r");
    if (!ifp)
        exit(1);
    seqlen = getmtx(ifp);
    fclose(ifp);

    for (i=0; i<seqlen; i++)
    {
        putchar(seq[i]);
        seq[i] = aanum(seq[i]);
    }

    putchar('\n');

    predict();

    return 0;
}
