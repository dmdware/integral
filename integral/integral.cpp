// integral.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define DR		0.0001f
#define SD		0.01f

double P(double i)
{
	return 1.0f / (i*i);
}

double pull(double r, double I)
{
	double dr = DR;
	double i;
	double k = 0;

	//printf("i=%0.19lf i<=%0.19lf i+=%0.19lf\r\n", 1 - DR, r, dr);

	for (i = 1 - DR; i > r; i -= dr)
	{
		k += P(i) * dr * I;
		//printf("pull%0.19lf dr%0.19lf\r\n I%0.19lf", P(i), dr, I);
	}

	dr = i - r;

	if (dr <= 0.0f)
		return k;

	k += P(i) * dr * I;

	return k;
}

double push(double R, double I)
{
	double dR = DR;
	double i;
	double k = 0;

	for (i = 1 - DR; i > R; i -= dR)
	{
		k += pull(i, I) * dR;
		//printf("pull%0.19lf dR%0.19lf\r\n", pull(i, I), dR);
	}

	dR = i - R;

	if (dR <= 0.0f)
		return k;

	k += pull(i, I) * dR;
	//printf("pull%0.19lf dR%0.19lf\r\n", pull(i, I), dR);

	return k;
}

#define NI	12

#define SR	0.15f
#define II	0.001f

int main()
{
	int i = 0;
	double r = SR;
	double rs[NI];
	for (; i < NI; ++i)
	{
		double k = push(r, II);
		double j = pull(r, II);
		double r2 = r - k - j
			;
		printf("t=%d r=%0.19lf->%0.19lf\r\n", i, r, r2);
		rs[i] = r;
		r = r2;
	}

	double vs[NI];
	for (i = 1; i < NI; ++i)
	{
		vs[i] = rs[i] - rs[i - 1];
		printf("t=%d v+Ir+a=%0.19lf\r\n", i, vs[i]);
	}

	double as[NI];
	for (i = 2; i < NI; ++i)
	{
		as[i] = vs[i] - vs[i - 1];
		printf("t=%d v(i)+Ir(i)+a(i) - v(i-1)+Ir(i-1)+a(i-1)=%0.19lf\r\n", i, as[i]);
	}

	double rs2[NI];
	double minI = 0.0f;
	double maxI = 0.2f;
	double midI = 0.0785f;
	double minGM = 0.00000001f;
	double maxGM = 7.5f;
	double midGM = (minGM + maxGM) / 2.0f;
	double r21 = SR;
	double r22;
	double vs2[NI];
	double as2[NI];
	double dr2[NI];

	goto reset;

	for (i = 0; i < NI; ++i)
	{
		as2[i] = midGM / (r21*r21);
		vs2[i] = (i==0? pow(2.0f * midGM / (r21*r21), 0.5f) :vs2[i - 1]) + as2[i];
		r22 = r21 - vs2[i] + r21 * midI;
		rs2[i] = r21;
		dr2[i] = rs2[i] - rs[i];
		r21 = r22;

		//if (r21 < 0.00000001)
			//goto dec;

		//printf("vs2[%d]=%0.19lf\r\n", i, vs2[i]);
		//printf("rs2[%d]=%0.19lf\r\n", i, rs2[i]);

			continue;

		reset:
			r21 = SR;
			vs2[0] = pow(2.0f * midGM / (r21*r21), 0.5f);
		i = -1;
		printf("reset\r\n");
		;
	}

narrow:
	if (rs2[i - 1] > rs[i - 1] - 0.00000001f)
	{
		inc:
		printf("rs2[%d] = %0.19lf   > rs[%d]=%0.19lf ++  GM%0.19lf\r\n", i-1, rs2[i-1], i-1, rs[i-1], midGM);
		for (i = 0; i < NI; ++i)
		{
			printf("t=%d   %0.19lf - %0.19lf = %0.19lf   a = %0.19lf   v = %0.19lf \r\n", i, rs2[i], rs[i], rs2[i] - rs[i], as2[i], vs2[i]);
		}
		printf("GM=%0.19lf    I=%0.19lf\r\n", midGM, midI);
		minGM = midGM;
		midGM = (midGM + maxGM + maxGM) / 3.0f;
		printf("GM=%0.19lf    I=%0.19lf\r\n", midGM, midI);
		//maxGM = (maxGM * 1.01f);
		system("pause");
		goto reset;
	}
	else if(rs2[i - 1] < rs[i - 1] + 0.00000001f)
	{
		dec:
		printf("rs2[%d] = %0.19lf   < rs[%d]=%0.19lf  -- GM%0.19lf\r\n", i-1, rs2[i-1], i - 1, rs[i - 1], midGM);
		for (i = 0; i < NI; ++i)
		{
			printf("t=%d   %0.19lf - %0.19lf = %0.19lf   a = %0.19lf  v = %0.19lf  \r\n", i, rs2[i], rs[i], rs2[i] - rs[i], as2[i], vs2[i]);
		}
		printf("GM=%0.19lf    I=%0.19lf\r\n", midGM, midI);
		maxGM = midGM;
		midGM = (midGM + minGM + minGM) / 3.0f;
		printf("GM=%0.19lf    I=%0.19lf\r\n", midGM, midI);
		//minGM = (minGM * 0.99f);
		system("pause");
		goto reset;
	}
	else
	{
		for (i = 0; i < NI; ++i)
		{
			printf("rs2[%d] = %0.19lf\r\n", i, rs2[i]);
		}
		printf("GM=%0.19lf    I=%0.19lf\r\n", midGM, midI);
		system("pause");
	}
	printf("narrow\r\n");
	goto reset;


	//printf("\r\npush=%0.19lf, pull=%0.19lf\r\n", k, j);
	system("pause");

    return 0;
}

