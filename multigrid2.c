#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int Nx, n_level, n_relax, k;
double dt, h, h2, xleft, xright;

void initialization(double* u_init) // 초기값 u_init
{
	int i;
	double x, pi = 4.0 * atan(1.0); // 4.0*pi/4 = pi

	for (i = 1; i <= Nx; i++)
	{
		x = xleft + ((double)i - 0.5) * h; // 각 단위구간에서의 중점값
		u_init[i] = cos(2.0 * pi * x); // 각 x값에서의 초기조건 데이터 조정
	}
}

void source(double* u_n, double* f)
{
	int i;

	for (i = 1; i <= Nx; i++)
	{
		f[i] = u_n[i] / dt; // f_i 값 조정
	}
}

void relaxGS(double* u_sp1, double* f, int Nxt) // Gauss Seiden
{
	int i, iter;
	double ht2, coeff, sor;


	ht2 = pow((xright - xleft) / (double)Nxt, 2);
	for (iter = 1; iter <= n_relax; iter++)
	{
		for (i = 1; i <= Nxt; i++)
		{
			coeff = 1.0 / dt;
			sor = f[i];
			if (i > 1)
			{
				coeff += 1.0 / ht2;
				sor += u_sp1[i - 1] / ht2;
			}
			if (i < Nxt)
			{
				coeff += 1.0 / ht2;
				sor += u_sp1[i + 1] / ht2;
			}
			u_sp1[i] = sor / coeff;
		}
	}
}

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

double Lap1D(double* u, int i, int Nxt) // L 연산자에서 공간과 관련된 항 계산
{
	double ht, dudx_L = 0.0, dudx_R = 0.0;

	ht = (xright - xleft) / (double)Nxt;
	if (i > 1)
		dudx_L = (u[i] - u[i - 1]) / ht;
	if (i < Nxt)
		dudx_R = (u[i + 1] - u[i]) / ht;
	return (dudx_R - dudx_L) / ht;
}

void restrict1D(double* uf, double* uc, int Nxc) // restriction; k -> k-1
{
	int i;

	for (i = 1; i <= Nxc; i++)
	{
		uc[i] = 0.5 * (uf[2 * i - 1] + uf[2 * i]);
	}
}

void prolong1D(double* uc, double* uf, int Nxc) // interpolation; k-1 -> k
{
	int i;

	for (i = 1;i <= Nxc;i++)
	{
		uf[2 * i - 1] = uc[i];
		uf[2 * i] = uc[i];
	}
}

void print_data(FILE* fpt, double* a, int i_start, int i_end) // matlab 시각화하기 위해 텍스트 파일로 기록
{
	int i;

	for (i = i_start; i <= i_end; i++)
		fprintf(fpt, "%16.15f \n", a[i]);
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

void vcycle(double* uf_new, double* ff, int Nxf, int ilevel)
{
	int i;

	relaxGS(uf_new, ff, Nxf); // pre-smoothing

	if (ilevel > 0)
	{
		int Nxc = Nxf / 2;
		double* fc, * defectc, * defect;

		fc = dvector(1, Nxc);
		defectc = dvector(1, Nxc);
		defect = dvector(1, Nxf);
		for (i = 1;i <= Nxf;i++)
		{
			defect[i] = ff[i] - uf_new[i] / dt + Lap1D(uf_new, i, Nxf); // 뒷 두 항 = L 연산자 적용된 값; defect 계산
		}
		restrict1D(defect, fc, Nxc); // k > k-1 restrict
		for (i = 1;i <= Nxc;i++)
			defectc[i] = 0.0; // 다시 V cycle 돌리기 위해 초기값 제로 함수로 두기.
		vcycle(defectc, fc, Nxc, ilevel - 1); // 다시 V cycle
		prolong1D(defectc, defect, Nxc); // Interpolation
		vec_add(uf_new, uf_new, defect, 1, Nxf); //  uf_new = uf_new + defect
		relaxGS(uf_new, ff, Nxf); // post-smoothing
		free_dvector(fc, 1, Nxc);
		free_dvector(defectc, 1, Nxc);
		free_dvector(defect, 1, Nxf);
	}
}

int main()
{
	char buffer[20];
	int it, ns, max_iteration, it_MG, max_it_MG, count = 1;
	double tol, error, * u_new, * u_old, * u_tmp, * f;
	FILE* fpt;

	xleft = 0.0; xright = 1.0; Nx = 64; max_iteration = 200;
	ns = max_iteration / 5; h = (xright - xleft) / (double)Nx; // ns = 출력되는 데이터 개수: 40단위로 출력
	h2 = pow(h, 2); dt = 0.5 * h2; n_relax = 2; // n_relax = \nu 
	n_level = (int)(log(Nx) / log(2.0) + 0.1) - 1; // k값; 여기서는 5
	u_old = dvector(1, Nx);
	u_new = dvector(1, Nx);
	u_tmp = dvector(1, Nx);
	f = dvector(1, Nx);
	initialization(u_new);
	sprintf(buffer, "Heat1D_MG.m");
	fpt = fopen(buffer, "w");
	print_data(fpt, u_new, 1, Nx);
	max_it_MG = 1000; tol = 1.0e-7; // 최대 반복 횟수, 허용 오차

	for (it = 1; it <= max_iteration; it++)
	{
		vec_copy(u_old, u_new, 1, Nx);
		source(u_old, f);
		it_MG = 0; error = 2.0 * tol;
		vec_copy(u_tmp, u_new, 1, Nx);
		while (it_MG <= max_it_MG && error > tol) // 최대 반복 횟수를 넘지 않고 주어진 한계 오차보다 작아질 때까지 Vcycle 반복
		{
			vcycle(u_new, f, Nx, n_level);
			error = diff_norm1D(u_tmp, u_new, Nx);
			vec_copy(u_tmp, u_new, 1, Nx);
			it_MG++;
		}
		printf("it=%3d : iteration %d error = %16.14f\n \n", it, it_MG, error);
		if (it % ns == 0)
		{
			print_data(fpt, u_new, 1, Nx); count++;
			printf("\n count=%d \n", count);
		}
	}
	fclose(fpt);
	printf("Nx      = %d\n", Nx);
	printf("max_it  = %d\n", max_iteration);
	printf("ns      = %d\n", ns);
	printf("dt      = %f\n", dt);
	printf("n_level = %d\n", n_level);
	return 0;
}