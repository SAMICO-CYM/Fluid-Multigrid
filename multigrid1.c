#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int Nx;
double dt, h, h2, xleft, xright;

void initialization(double *u_init) // 초기값 u_init
{
	int i;
	double x, pi = 4.0 * atan(1.0); // 4.0*pi/4 = pi

	for (i = 1; i <= Nx; i++)
	{
		x = xleft + ((double)i - 0.5) * h; // 각 단위구간에서의 중점값
		u_init[i] = cos(2.0 * pi * x); // 각 x값에서의 초기조건 데이터 조정
	}
}

void source(double *u_n, double *f)
{
	int i;

	for (i = 1; i <= Nx; i++)
	{
		f[i] = u_n[i] / dt; // f_i 값 조정
	}
}

void GaussSeidel(double *u_sp1, double *f, int Nx) // 가우스-세이델 구현, u_sp1 = u^{s+1}
{
	int i;
	double coeff, sor;

	for (i = 1; i <= Nx; i++)
	{
		coeff = 1.0 / dt;
		sor = f[i];
		if (i > 1)
		{
			coeff += 1.0 / h2;
			sor += u_sp1[i - 1] / h2;
		}
		if (i < Nx)
		{
			coeff += 1.0 / h2;
			sor += u_sp1[i + 1] / h2;
		}
		u_sp1[i] = sor / coeff;
	}
}

double *dvector(long i_start, long i_end)
{
	double *v;

	v = (double*)malloc((i_end - i_start + 1 + 1) * sizeof(double));
	return v - i_start + 1;
}

void free_dvector(double *v, long i_start, long i_end)
{
	free(v + i_start - 1);
}

void print_data(FILE* fpt, double *a, int i_start, int i_end)
{
	int i;

	for (i = i_start; i <= i_end; i++)
		fprintf(fpt, "%16.15f \n", a[i]);
}

void vec_copy(double* a, double* b, int i_start, int i_end)
{
	int i;

	for (i = i_start; i <= i_end;i++)
		a[i] = b[i];
}

void vec_sub(double* a, double* b, double* c, int i_start, int i_end)
{
	int i;

	for (i = i_start; i <= i_end;i++)
		a[i] = b[i] - c[i];
}

double norm1D(double* a, int i_start, int i_end)
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

int main(void)
{
	char buffer[20];
	int n, print_interval, Nt, it_GS, max_iteration, count = 1;
	double tol, error, * u_np1, * u_n, * f;
	FILE* fpt;
	
	xleft = 0.0; xright = 1.0;
	Nx = 1024; Nt = 100; print_interval = Nt / 5;
	h = (xright - xleft) / (double)Nx; h2 = pow(h, 2);
	dt = 2 * h2;
	max_iteration = 1000; tol = 1.0e-7;

	u_n = dvector(1, Nx);
	u_np1 = dvector(1, Nx);
	f = dvector(1, Nx);
	initialization(u_n);
	sprintf(buffer, "Heat1D_Gs.m");
	fpt = fopen(buffer, "w");
	print_data(fpt, u_n, 1, Nx);
	printf("\nThe nubmer of prints is %d \n", count);
	vec_copy(u_np1, u_n, 1, Nx);

	for (n = 1;n <= Nt;n++)
	{
		source(u_n, f);
		it_GS = 0; error = 1.0;
		while (it_GS <= max_iteration && error > tol) 
		{
			GaussSeidel(u_np1, f, Nx); // GS로 다음 n step의 u값 구하기
			error = diff_norm1D(u_n, u_np1, Nx); // n+1과 n과의 차이
			vec_copy(u_n, u_np1, 1, Nx); // n+1을 다시 n에 덮어쓰기
			it_GS++; // 반복수 증가
		}

		printf("n = %3d GS Iterations %d error = %16.14f\n", n, it_GS, error);
		if (n % print_interval == 0) // 20 반복단위마다 기록
		{
			print_data(fpt, u_n, 1, Nx); count++;
			printf("\nThe number of prints is %d \n\n", count);
		}
	}

	fclose(fpt);
	printf("Nx = %d\n", Nx);
	printf("Nx = %d\n", Nt);
	printf("dt = %f\n", dt);
	printf("print_interval = %d\n", print_interval);
	return 0;
}