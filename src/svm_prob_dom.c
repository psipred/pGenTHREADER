/* Estimate 3Dhit score for pDomTHREADER hit & calculate p-value */

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
    int             i, niters, laln, lena, lenb, margflg, sfrom,sto;
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
	if (fscanf(ifp, "%f%f%*s%f%f%*s%*s%*s%*s%*s%f%d%d%d%d%d%s", &pair, &zpair, &solv, &zsolv, &nwsc, &laln, &lena, &lenb, &sfrom,&sto,brkid) != 11)
	    break;

	/* Use linear SVM to estimate 3Dhit score */

            svm_score = -10.627101;
	    svm_score += -0.287522   * (1/(1+exp(-0.1*(pair+50))));
	    svm_score += 0.18282     * (1/(1+exp(-1*(solv+1))));
	    svm_score += 21.70844    * (1/(1+exp(-0.01*(nwsc-50))));
	    svm_score += 0.644233    * (laln/lena);
            svm_score += -0.19508    * (1/(1+exp(-0.05*(solv-50))));
	
	prob = 1/(1+ exp(1.16808 * (svm_score+3.66866)));

	if (prob > 1.0)
	  prob = 1.0;
	if (prob < 0.0)
	  prob = 0.0;

        if(svm_score > 0 )
        {
	if (prob < 0.01)
	{
	    if (prob < 0.00001)
		printf("CERT   ");
	    else if (prob < 0.0001)
		printf("HIGH   ");
	    else if (prob < 0.001)
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
	printf(" %7.1f %5.1f %6.1f %4d %4d %4d %4d %4d %s\n", pair, solv, nwsc, laln, lena, lenb, sfrom, sto, brkid);
        }
    }

    return 0;
}
