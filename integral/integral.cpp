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


#define DR2		0.00035f
#define DR3		0.00001f

double Q2(double r, double I)
{
	double dr = DR3;
	double i;
	double k = 0;


	//printf("pull%0.19lf dR%0.19lf\r\n", r, dr);


	//i = r;
	i = 1;
	//for (i = r; i + dr <= 1; i += dr)
	//if(i + dr <= 1)
	for (; i - dr >= r; i -= dr)
	{
		double q2 = k;// Q(i + dr, I);
		double mp = q2;// fmax(q2, I);
		if (i >= 1)
			mp = I;
		//double mp = fmax(q2, (i*i)/((i-I)*(i-I)) );
		//double s = (1 + I);
		//double mp = (-(1.0f - i) * s + (1.0f - (i-dr)) * s) + q2;
		double qq = (i - mp);
		//double qq = - dr + mp;
		//if (qq <= 0)
		//return 0;
		//if (i != qq)
		{
			//printf("mp%0.19lf       s%0.19lf      si%0.19lf    q2   %0.19lf       dr/r=%0.19lf\r\n", mp, s, (-(1.0f - i) * s + (1.0f - (i - dr)) * s), q2,    dr/r);
			//printf("i%0.19lf qq%0.19lf r%0.19lf dr%0.19lf mp%0.19lf\r\n", i, qq, r, dr, mp);
		}
		//k += (-1.0f + (i*i) / (qq*qq)) * dr;
		//k = (-1.0f + (i*i) / (qq*qq)) * dr + mp;
		k = (-1.0f + (i*i) / (qq*qq)) * dr + mp;
		//printf("pull%0.19lf dR%0.19lf    i%0.19lf    qq%0.19lf    dd%0.19lf    iq%0.19lf\r\n", i, dr, i, qq, (-1.0f + (i*i) / (qq*qq)) * dr, (i*i) / (qq*qq));
		//printf("mp->k  = %0.19lf    ->    %0.19lf\r\n", mp, k);
		//system("pause");
	}

	//dr = 1 - r;

	//if (dr <= 0.0f)
	//return k;
	//return k;

	dr = i - r;

	if (dr <= 0.0f)
		return k;

	//return k;

	{
		double q2 = k;// Q(i + dr, I);
		double mp = fmax(q2, I);
		if (i >= 1)
			mp = I;
		//double mp = fmax(q2, (i*i) / ((i - I)*(i - I)));
		//double s = (1 + I);
		//double mp = (-(1.0f - i) * s + (1.0f - (i - dr)) * s) + q2;
		double qq = (i - mp);
		//double qq = -dr + mp;
		//if (qq <= 0)
		//return 0;
		//k += (-1.0f + (i*i) / (qq*qq)) * dr;
		//k = (-1.0f + (i*i) / (qq*qq)) * dr + mp;
		k = (-1.0f + (i*i) / (qq*qq)) * dr + mp;
	}

	return k;
}

double Q(double r, double I)
{
	double dr = DR2;
	double i;
	double k = 0;


	//printf("pull%0.19lf dR%0.19lf\r\n", r, dr);

	i = r;
	//for (i = r; i + dr <= 1; i += dr)
	if(i + dr <= 1)
	{
		double q2 = Q(i + dr, I);
		double mp = fmax(q2, I);
		double qq = (i - mp);
		//if (qq <= 0)
			//return 0;
		//if (i != qq)
		{
			//printf("mp%0.19lf\r\n", mp);
			//printf("i%0.19lf qq%0.19lf r%0.19lf dr%0.19lf mp%0.19lf\r\n", i, qq, r, dr, mp);
		}
		k += (-1.0f + (i*i) / (qq*qq)) * dr + mp;
		//printf("pull%0.19lf dR%0.19lf\r\n", i, dr);
	}

	dr = 1 - r;

	if (dr <= 0.0f)
		return k;
	
	{
		double q2 = Q(i + dr, I);
		double mp = fmax(q2, I);
		double qq = (i - mp);
		//if (qq <= 0)
			//return 0;
		k += (-1.0f + (i*i) / (qq*qq)) * dr + mp;
	}

	return k;
}

#define NI	13
#define SR	0.2f
#define II	0.001f
#define III	II

void tQ()
{


	int i = 0;
	double r = SR;// 0.3f;
	double rs[NI];
	double rs2[NI];
	double r21 = SR;// 0.3f;
	for (; i < NI; ++i)
	{
		double k = push(r, III);
		double j = pull(r, III);
		double r2 = r - k - j
			;
#define	I4	0.000159f
		//double r22 = r21 - Q2(r + 0.00134f, 0.000014f) - 0.00134f;
		double r22 = r21 - Q2(r - 0.02f, 0.000134f);
		printf("t=%d    r=%0.19lf   d=%0.19lf -     %0.19lf       =      %0.19lf       %% %0.19lf     d%% %0.19lf \r\n", i, r22,  r22-r21, r2, r22-r2, (r22 - r2)/r2*100.0f, fabs(r22-r21)/(r22-r21)*(r22 - r21)/(r2-r)*100.0f);
		rs[i] = r;
		r = r2;
		r21 = r22;
	}
	system("pause");
}


int main()
{
	tQ();

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
	double midI = 0.001;
	double minGM = 0.00000001f;
	double maxGM = 7.5f;
	double midGM = II;	// (minGM + maxGM) / 2.0f;
	double r21 = SR;
	double r22;
	double vs2[NI];
	double as2[NI];
	double dr2[NI];

	goto reset;

	for (i = 0; i < NI; ++i)
	{
		as2[i] = midGM / (r21*r21);
		vs2[i] = (i==0? pow(2.0f * midGM / (r21), 0.5f) :vs2[i - 1]) + as2[i];
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
			vs2[0] = pow(2.0f * midGM / (r21), 0.5f);
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

