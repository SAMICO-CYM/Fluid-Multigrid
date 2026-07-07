#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int Nx, n_level, n_relax, k;
double dt, h, h2, xleft, xright, Cahn, gam;

double* dvector(long i_start, long i_end)
{
	double* v;

	v = (double*)malloc((i_end - i_start + 1 + 1) * sizeof(double));
	return v - i_start + 1; // 인덱스 1 추가
}

void free_dvector(double* v, long i_start, long i_end)
{
	free(v + i_start - 1);
}

void vec_copy(double* a, double* b, int i_start, int i_end) // 복사본 생성
{
	int i;

	for (i = i_start; i <= i_end;i++)
		a[i] = b[i];
}

void vec_add(double* a, double* b, double* c, int i_start, int i_end) // 벡터 덧셈
{
	int i;

	for (i = i_start; i <= i_end;i++)
		a[i] = b[i] + c[i];
}

void vec_sub(double* a, double* b, double* c, int i_start, int i_end) // 벡터 뺄셈
{
	int i;

	for (i = i_start; i <= i_end;i++)
		a[i] = b[i] - c[i];
}

double norm1D(double* a, int i_start, int i_end) // norm 계산
{
	int i;
	double value = 0.0;

	for (i = i_start; i <= i_end; i++)
		value += a[i] * a[i];
	return sqrt(value / (i_end - i_start + 1.0)); // i_start = 1, i_end = Nx
}

double diff_norm1D(double* a, double* b, int Nx) // a-b 벡터의 norm
{
	double* d, value;

	d = dvector(1, Nx);
	vec_sub(d, a, b, 1, Nx);
	value = norm1D(d, 1, Nx);
	free_dvector(d, 1, Nx);
	return value;
}

void init_data(double* oc)
{
	int i;

	for (i = 1;i <= Nx;i++)
	{
		oc[i] = 0.2 * (1.0 - 2.0 * rand() / (double)RAND_MAX);
	}
}

void source(double* oc, double* f)
{
	int i;

	for (i = 1;i <= Nx;i++)
	{
		f[i] = oc[i] / dt + oc[i];
	}
}

void restrict1D(double* cf, double* cc, int Nxc) // restriction; k -> k-1
{
	int i;

	for (i = 1; i <= Nxc; i++)
	{
		cc[i] = 0.5 * (cf[2 * i - 1] + cf[2 * i]);
	}
}

void prolong1D(double* cc, double* cf, int Nxc) // interpolation; k-1 -> k
{
	int i;

	for (i = 1;i <= Nxc;i++)
	{
		cf[2 * i - 1] = cc[i];
		cf[2 * i] = cc[i];
	}
}

double Lap(double* c, int i, int Nxt) // L 연산자에서 공간과 관련된 항 계산
{
	double ht, dcdx_L = 0.0, dcdx_R = 0.0;

	ht = (xright - xleft) / (double)Nxt;
	if (i > 1)
		dcdx_L = (c[i] - c[i - 1]) / ht;
	if (i < Nxt)
		dcdx_R = (c[i + 1] - c[i]) / ht;
	return (dcdx_R - dcdx_L) / ht;
}

void nonL(double* NSO, double* nc, int Nxt)
{
	int i;

	for (i = 1;i <= Nxt;i++)
	{
		NSO[i] = nc[i] / dt + pow(nc[i], 3) - Cahn * Lap(nc, i, Nxt); // N 연산자
	}
}

void csource(double* fc, double* cf_new, double* ff, double* cc_new, int Nxc)
{
	double* NSO, * NSOC, * defectf, * defectc;

	NSO = dvector(1, 2 * Nxc);
	NSOC = dvector(1, Nxc);
	defectf = dvector(1, 2 * Nxc);
	defectc = dvector(1, Nxc);
	nonL(NSO, cf_new, 2 * Nxc);
	nonL(NSOC, cc_new, Nxc);
	vec_sub(defectf, ff, NSO, 1, 2 * Nxc); // alpha_k
	restrict1D(defectf, defectc, Nxc);
	vec_add(fc, defectc, NSOC, 1, Nxc); // f^n_{k-1}
	free_dvector(NSO, 1, 2 * Nxc);
	free_dvector(NSOC, 1, Nxc);
	free_dvector(defectf, 1, 2 * Nxc);
	free_dvector(defectc, 1, Nxc);
}

void relaxGS(double* nc, double* f, int Nxt) // Gauss Seiden
{
	int i, iter;
	double ht2, coeff, sor;


	ht2 = pow((xright - xleft) / (double)Nxt, 2);
	for (iter = 1; iter <= n_relax; iter++)
	{
		for (i = 1; i <= Nxt; i++)
		{
			coeff = 1.0 / dt + 3.0*pow(nc[i], 2);
			sor = f[i] + 2.0*pow(nc[i], 3);
			if (i > 1)
			{
				coeff += Cahn / ht2;
				sor += Cahn*nc[i - 1] / ht2;
			}
			if (i < Nxt)
			{
				coeff += Cahn / ht2;
				sor += Cahn * nc[i + 1] / ht2;
			}
			nc[i] = sor / coeff;
		}
	}
}

void print_vector(FILE* fpt, double* a, int nl, int nh) // matlab 시각화하기 위해 텍스트 파일로 기록
{
	int i;

	for (i = nl; i <= nh; i++)
		fprintf(fpt, "%16.14f \n", a[i]);
}

void vcycle(double* cf_new, double* ff, int Nxf, int ilevel)
{
	int i;

	relaxGS(cf_new, ff, Nxf); // pre-smoothing

	if (ilevel > 0)
	{
		int Nxc = Nxf / 2;
		double* cc_new, * fc, * wcc_new, * correct_c;

		cc_new = dvector(1, Nxc);
		fc = dvector(1, Nxc);
		wcc_new = dvector(1, Nxf);
		correct_c = dvector(1, Nxf);
		restrict1D(cf_new, cc_new, Nxc); // k > k-1 restrict
		csource(fc, cf_new, ff, cc_new, Nxc);
		vec_copy(wcc_new, cc_new, 1, Nxc);
		vcycle(wcc_new, fc, Nxc, ilevel - 1); // 다시 V cycle
		vec_sub(wcc_new, wcc_new, cc_new, 1, Nxc);
		prolong1D(wcc_new, correct_c, Nxc); // Interpolation
		vec_add(cf_new, cf_new, correct_c, 1, Nxf); //  uf_new = uf_new + defect
		relaxGS(cf_new, ff, Nxf); // post-smoothing
		free_dvector(cc_new, 1, Nxc);
		free_dvector(fc, 1, Nxc);
		free_dvector(wcc_new, 1, Nxc);
		free_dvector(correct_c, 1, Nxf);
	}
}

int main()
{
	extern int Nx, n_level, n_relax, k;
	extern double dt, h, h2, xleft, xright, Cahn, gam;
	char buffer[20], buffermu[20];
	int n, ns, max_it, it_MG, max_it_MG;
	double T, tol, error, * nc, * oc, * c_tmp, * f;
	FILE* fpt;

	xleft = 0.0; xright = 1.0; Nx = 128;
	h = (xright - xleft) / Nx; h2 = pow(h, 2);
	gam = 4.0 * h / (2 * sqrt(2.0) * atanh(0.9)); Cahn = pow(gam, 2);
	T = 10.0; max_it = 100; dt = T / max_it; ns = max_it / 5;
	n_relax = 3; max_it_MG = 100; tol = 1.0e-7;
	n_level = (int)(log(Nx) / log(2.0) + 0.1) - 1;
	oc = dvector(1, Nx);
	nc = dvector(1, Nx);
	c_tmp = dvector(1, Nx);
	f = dvector(1, Nx);
	init_data(nc);
	sprintf(buffer, "AC1D.m");
	fpt = fopen(buffer, "w");
	print_vector(fpt, nc, 1, Nx);

	for (n = 1;n <= max_it; n++)
	{
		vec_copy(oc, nc, 1, Nx);
		source(oc, f);
		it_MG = 0; error = 2.0 * tol;
		vec_copy(c_tmp, nc, 1, Nx);
		while (it_MG <= max_it_MG && error > tol)
		{
			vcycle(nc, f, Nx, n_level);
			error = diff_norm1D(c_tmp, nc, Nx);
			vec_copy(c_tmp, nc, 1, Nx); it_MG++;
		}
		printf("n = %3d MG iterations = %d error = %16.14f\n", n, it_MG, error);

	}
	if (n % ns == 0)
	{
		print_vector(fpt, nc, 1, Nx);
	}
	fclose(fpt);
	return 0;
}