/***************************************
 *   Threading potential definitions   *
 * Copyright (C) 1994 David T. Jones   *
 ***************************************/

#ifdef __alpha
#define long int
#endif

#define TOPOMAX 11
#define INTERVALS 20
#define NACC 5
#define OOIMIN 6
#define OOIMAX 30
#define OOIDIV 1
#define NRR 40
#define RRWIDTH (1.0F)
#define RTCONST (0.582F)
#define SIGMA (0.02F)
#define unSCALED

enum atmcodes
{
    CAATOM, CBATOM, OATOM, NATOM, CATOM
};

enum aacodes
{
    ALA, ARG, ASN, ASP, CYS,
    GLN, GLU, GLY, HIS, ILE,
    LEU, LYS, MET, PHE, PRO,
    SER, THR, TRP, TYR, VAL,
    GAP, UNK
};


int             pairpot(const int, const int, const int, const int, const int, const float);
int             potinit(void);
int             solvpot(const int, const int);
