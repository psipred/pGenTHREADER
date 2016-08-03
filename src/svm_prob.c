/* Estimate 3Dhit score for mGenTHREADER hit & calculate p-value */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

#define TRUE 1
#define FALSE 0

#define SQR(x) ((x)*(x))
#define MIN(x,y) ((x)<(y)?(x):(y))
#define MAX(x,y) ((x)>(y)?(x):(y))

/* logistic 'squashing' function (+/- 1.0) */
#define logistic(x) (1.0 / (1.0 + exp(-(x))))


void
err(s)
    char           *s;
{
    fprintf(stderr, "%s\n", s);
}

void
fail(s)
    char           *s;
{
    fprintf(stderr, "%s\n", s);
    exit(1);
}


main(argc, argv)
    int             argc;
    char          **argv;
{
    int             i, niters, laln, lena, lenb, margflg;
    FILE *ifp, *cfp;
    char brkid[40], guess[160], buf[256];
    float nwsc, pair, zpair, solv, zsolv, svm_score, tprob, prob, maxprob, wt[9];

    /* malloc_debug(3); */
    if (argc != 2 && argc != 11)
	fail("usage : psi_prob input-file {svm weights}");

    if (argc == 11)
	for (i=0; i<9; i++)
	    wt[i] = atof(argv[2+i]);
    
    ifp = fopen(argv[1], "r");
    if (!ifp)
	fail("Cannot open file!");

    maxprob = margflg = 0;
    guess[0] = '\0';

    while (!feof(ifp))
    {
	if (fscanf(ifp, "%f%f%*s%f%f%*s%*s%*s%*s%*s%f%d%d%d%s", &pair, &zpair, &solv, &zsolv, &nwsc, &laln, &lena, &lenb, brkid) != 9)
	    break;

	/* Use linear SVM to estimate 3Dhit score */

	if (argc != 11)
	{
	    svm_score =  35.85;
	    svm_score += -4.000026    * (pair + 117.865) / 98.625;
	    svm_score += -0.355853    * zpair;
	    svm_score += -4.515154    * (solv + 1.467) / 4.031;
	    svm_score += -2.688433    * zsolv;
	    svm_score += 21.434526    * (nwsc - 124.79) / 148.838;
	    svm_score += 1.352219     * (laln - 112.307) / 58.384;
	    svm_score += 0.000000     * (lena - 261.948) / 181.058;
	    svm_score += 0.081323     * (lenb - 160.821) / 79.312;
	}
	else
	{
	    svm_score = wt[0];
	    svm_score += wt[1] * (pair + 117.865) / 98.625;
	    svm_score += wt[2] * zpair;
	    svm_score += wt[3] * (solv + 1.467) / 4.031;
	    svm_score += wt[4] * zsolv;
	    svm_score += wt[5] * (nwsc - 124.79) / 148.838;
	    svm_score += wt[6] * (laln - 112.307) / 58.384;
	    svm_score += wt[7] * (lena - 261.948) / 181.058;
	    svm_score += wt[8] * (lenb - 160.821) / 79.312;
	}

	if (laln < 20)
	    svm_score = laln;
	
	prob = 50.807224 * exp(-0.232529 * svm_score);

	if (prob > 1.0)
	  prob = 1.0;
	if (prob < 0.0)
	  prob = 0.0;

	if (prob < 0.1)
	{
	    if (prob < 0.0001)
		printf("CERT   ");
	    else if (prob < 0.001)
		printf("HIGH   ");
	    else if (prob < 0.01)
		printf("MEDIUM ");
	    else
		printf("LOW    ");
	}
	else
	    printf("GUESS  ");
	if (prob < 0.001)
	    printf("%7.3f %6.0e", svm_score, prob);
	else
	    printf("%7.3f %6.3f", svm_score, prob);

	printf(" %7.1f %5.1f %6.1f %4d %4d %4d %s\n", pair, solv, nwsc, laln, lena, lenb, brkid);
    }
}
