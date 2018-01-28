// integral.cpp : Defines the entry point for the console application.
//

#include "sys/includes.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "sys/texture.h"

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

//#define NI	128
#define NI	12
//#define SR	0.9f
#define II	(0.001f/1.0f)
#define III	II


#define SR	0.15f //0.42f

void tQ()
{
	return;

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

void line(double x1, double y1, double x2, double y2, int wx, int wy, int channels, unsigned char v, unsigned char usec, unsigned char *data, double skipi, double skn)
{
	double i;
	double d = sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
	int x;
	int y;

	for (i = skipi; i <= d; i += skn)
	{
		x = x1 + (x2 - x1)*i/d;
		y = y1 + (y2 - y1)*i/d;

		if (x < 0)
			continue;
		if (y < 0)
			continue;
		if (x >= wx)
			continue;
		if (y >= wy)
			continue;

		data[channels*wx*y + channels*x + usec] = v;
	}
}

void draw(double *rs, double *rs2, double *vs, double *vs2, double midI, double midGM)
{
	texdata t;
	int i;
	int j;
	ltexinit(&t);
	double maxv = -99999999999;
	double minv = 9999999999999;
	double maxp = -99999999999;
	double minp = 9999999999999;
	
#define WX	640
#define WY	920

	t.sizex = WX;
	t.sizey = WY;
	t.channels = 3;
	t.data = (unsigned char*)malloc(sizeof(unsigned char)*WX*WY * 3);

	memset(t.data, -1, sizeof(unsigned char) *WX*WY * 3);

	for (i = 0; i < NI; ++i)
	{
		if (vs[i] > maxv)
			maxv = vs[i];
		if (vs2[i] > maxv)
			maxv = vs2[i];
		if (vs[i] < minv)
			minv = vs[i];
		if (vs2[i] < minv)
			minv = vs2[i];

		if (rs[i] > maxp)
			maxp = rs[i];
		if (rs2[i] > maxp)
			maxp = rs2[i];
		if (rs[i] < minp)
			minp = rs[i];
		if (rs2[i] < minp)
			minp = rs2[i];
	}

	line(0, WY - WY*((-minp)) / (maxp - minp), WX - 1,WY -  WY*(0 - minp) / (maxp - minp), WX, WY, 3, 0, 0, t.data, 0, 3);

	for (i = 0; i < NI - 1; ++i)
	{
		line(WX*(double)i / (double)(NI-1), WY-  WY*(rs[i] - minp) / (maxp - minp), WX*(double)(i + 1) / (double)(NI-1), WY - WY*(rs[i + 1] - minp) / (maxp - minp), WX, WY, 3, 0, 1, t.data, 1, 3);
		line(WX*(double)i / (double)(NI-1), WY -WY*(rs2[i] - minp) / (maxp - minp), WX*(double)(i + 1) / (double)(NI-1), WY - WY*(rs2[i + 1] - minp) / (maxp - minp), WX, WY, 3, 0, 2, t.data, 2, 3);
	}

	savepng("out.png", t.data, t.sizex, t.sizey, t.channels);
	ltexfree(&t);
}

void upp(double r0, double *v0, double midI, double midGM, double v2)
{
	double r = 1;
	double v = 0;
	double a;

	while (r > r0 && r==r && fabs(r)!=INFINITY)
	{
		a = midGM / (r*r);
		v = v + a;
		//r = r - v - (1 - r) * midI;
		r = r - v + r * midI;
	}

	*v0 = v;

	*v0 = pow(2.0f * midGM / (r0*r0), 0.5f);
	//*v0 = v2 - (1 - r) * midI;
	//*v0 = 0;
}


void render()
{

	int i = 0;
	double r = SR;
	double rs[NI];
	for (; i < NI; ++i)
	{
		double k = push(r, II);
		double j = pull(r, II);
		double r2 = r - k - j
		//double r2 = r - Q2(r, 0.000134f);
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
	double midI = (0.0785000026226043701/1.0f);// 0.00000019;
	double minGM = 0.0000000000000000000001f;
	double maxGM = 70.5f;
	double midGM = (0.0000034402408530694/1.0f);// II;	// (minGM + maxGM) / 2.0f;
	double r21 = SR;
	double r22;
	double vs2[NI];
	double as2[NI];
	double dr2[NI];
	double absmax;
	double maxs;
	int maxi;
	double lastmidGM = 0;
	double lastmidI = 0;
	double smin;
	double smax;

#define CHECK	0.0000000000000000001f
#define CHECK2	0.0000000000000000001f

	goto reset;
	while (true)
	{

		for (i = 0; i < NI; ++i)
		{
			as2[i] = midGM / (r21*r21);
			vs2[i] = (i == 0 ? (pow(2.0f * midGM / (r21), 0.5f)) : vs2[i - 1]) + as2[i];

			vs2[i] = (i == 0 ? pow(2.0f * midGM / (r21*r21), 0.5f) : vs2[i - 1]) + as2[i];
			//if (i == 0)
			//	upp(r21, &vs2[i], midI, midGM, vs[i]);
			//vs2[i] = (i == 0 ? 0 : vs2[i - 1]) + as2[i];
			//r22 = r21 - vs2[i] - (1 - r21) * midI;
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

			if (fabs(lastmidGM - midGM) < CHECK)
			{
				for (i = 0; i < NI; ++i)
				{
					printf("t=%d   %0.19lf - %0.19lf = %0.19lf   a = %0.19lf  v = %0.19lf  \r\n", i, rs2[i], rs[i], rs2[i] - rs[i], as2[i], vs2[i]);
				}
				printf("GM=%0.19lf    I=%0.19lf\r\n", midGM, midI);
				printf("GM=%0.19lf    I=%0.19lf\r\n", midGM, midI);

				goto advi;
			}

			r21 = SR;
			//upp(r21, &vs2[i], midI, midGM, vs[i]);
			vs2[0] = ( pow(2.0f * midGM / (r21), 0.5f));
			i = -1;
			//printf("reset\r\n");
			lastmidGM = midGM;
			;
		}

		goto done;

	narrow:
		if (rs2[i - 1] > rs[i - 1] //- 0.00000001f
			)
		{
		inc:
			//printf("rs2[%d] = %0.19lf   > rs[%d]=%0.19lf ++  GM%0.19lf\r\n", i - 1, rs2[i - 1], i - 1, rs[i - 1], midGM);
			for (i = 0; i < NI; ++i)
			{
				//printf("t=%d   %0.19lf - %0.19lf = %0.19lf   a = %0.19lf   v = %0.19lf \r\n", i, rs2[i], rs[i], rs2[i] - rs[i], as2[i], vs2[i]);
			}
			//printf("GM=%0.19lf    I=%0.19lf\r\n", midGM, midI);
			minGM = midGM;
			midGM = (midGM + maxGM + maxGM) / 3.0f;
			//printf("GM=%0.19lf    I=%0.19lf\r\n", midGM, midI);
			//maxGM = (maxGM * 1.01f);
			//system("pause");
			goto reset;
		}
		else if (rs2[i - 1] < rs[i - 1] //+ 0.00000001f
			)
		{
		dec:
			//printf("rs2[%d] = %0.19lf   < rs[%d]=%0.19lf  -- GM%0.19lf\r\n", i - 1, rs2[i - 1], i - 1, rs[i - 1], midGM);
			for (i = 0; i < NI; ++i)
			{
				//printf("t=%d   %0.19lf - %0.19lf = %0.19lf   a = %0.19lf  v = %0.19lf  \r\n", i, rs2[i], rs[i], rs2[i] - rs[i], as2[i], vs2[i]);
			}
			//printf("GM=%0.19lf    I=%0.19lf\r\n", midGM, midI);
			maxGM = midGM;
			midGM = (midGM + minGM + minGM) / 3.0f;
			//printf("GM=%0.19lf    I=%0.19lf\r\n", midGM, midI);
			//minGM = (minGM * 0.99f);
			//system("pause");
			goto reset;
		}
		else
		{
			//for (i = 0; i < NI; ++i)
			{
			//	printf("rs2[%d] = %0.19lf\r\n", i, rs2[i]);
			}
			//printf("GM=%0.19lf    I=%0.19lf\r\n", midGM, midI);
			//system("pause");
			printf("narrow\r\n");
			goto advi;
		}
		printf("narrow\r\n");
		goto reset;
	}

advi:

	printf("%0.19lf - %0.19lf\r\n", lastmidI, midI);
	if (fabs(lastmidI - midI) < CHECK2)
	{
	done:
		printf("%0.19lf - %0.19lf\r\n", lastmidI, midI);
		draw(rs, rs2, vs, vs2, midI, midGM);
		system("pause");
		exit(0);
	}

	goto done;

	lastmidGM = 0;

	//minI = 0.0f;
	//maxI = 0.2f;
	//midI = 0.019;
	minGM = 0.0000000000000000000001f;
	maxGM = 70.5f;
	midGM = II;	// (minGM + maxGM) / 2.0f;
	absmax = 0;
	maxs = 0;
	smin = 0;
	smax = 0;

	for (i = 0; i < NI; ++i)
	{
		if (fabs(rs2[i] - rs[i]) > absmax)
		{
			maxs = rs2[i] - rs[i];
			absmax = fabs(maxs);
		}

		if (rs2[i] - rs[i] < smin)
			smin = rs2[i] - rs[i];
		if (rs2[i] - rs[i] > smax)
			smax = rs2[i] - rs[i];
	}

#if 0
	if (smin < 0 && smax > 0)
	{

		lastmidI = midI;
		//printf("midI %0.19lf -> %0.19lf\r\n", midI, (midI + maxI) / 2.0);
		minI = midI;
		midI = (midI + maxI) / 2.0;
	}
	else 
#endif
		if (maxs < 0)
	{
		lastmidI = midI;
		//printf("midI %0.19lf -> %0.19lf\r\n", midI, (midI + minI) / 2.0);
		maxI = midI;
		midI = (midI + minI) / 2.0;
	}
	else if (maxs > 0)
	{
		lastmidI = midI;
		//printf("midI %0.19lf -> %0.19lf\r\n", midI, (midI + maxI) / 2.0);
		minI = midI;
		midI = (midI + maxI) / 2.0;
	}
	else
	{
		goto done;
	}

	goto reset;

	//printf("\r\npush=%0.19lf, pull=%0.19lf\r\n", k, j);
	system("pause");
}


int main()
{

#if 0
#define INCALS	#include "alskdjasd.txt"
	char *tt =
		//#include "begraw.h"
#include "begraw.h"
#include "alskdjasd.txt"
#include "endraw.h"
		
		R"###(
#include "alskdjasd.txt"
//#include "endraw.h"
)###";
#endif

	render();

	tQ();

	int i = 0;
	double r = SR;
	double rs[NI];
	for (; i < NI; ++i)
	{
		double k = push(r, II);
		double j = pull(r, II);
		//double r2 = r - k - j
		double r2 = r - Q2(r, 0.000134f);
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
	double midI = 0.019;
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

