#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int op = 0;

void imprimeM(int n, float **m, float *b);
void leMatriz(float **m, int n, float dA, float dB, float dP, float dC, float dD);
float **criaMatriz(int n);
void liberaMatriz(int n, float **m);
float *criaB(int n, float **m);
void triangulariza(int n, float **m, float *b);
float *substituicaoRegressiva(int n, float **m, float *b);
void imprimeVetor(int n, float *a);
void gauss(int n, float **m);
void seidel(int n, float **m);
float *criaX(int n);
float *criaXSeidel(int n, float *b, float **m);

int main()
{
	int n;
	float **m;
	float *b, *x;
	float dA, dB, dP, dC, dD;

	printf("Insira o valor de N: ");
	scanf("%d", &n);
	printf("Digite o valor de dA: ");
	scanf("%f", &dA);
	printf("Digite o valor de dB: ");
	scanf("%f", &dB);
	printf("Digite o valor de dP: ");
	scanf("%f", &dP);
	printf("Digite o valor de dC: ");
	scanf("%f", &dC);
	printf("Digite o valor de dD: ");
	scanf("%f", &dD);

	m = criaMatriz(n);

	leMatriz(m, n, dA, dB, dP, dC, dD);
	gauss(n, m);

	leMatriz(m, n, dA, dB, dP, dC, dD);
	seidel(n, m);

	liberaMatriz(n, m);
	return 0;
}

//
void gauss(int n, float **m)
{

	float *b = criaB(n, m);
	
	printf("\n\t-- Matriz Inicial --\n");
	imprimeM(n, m, b);

	triangulariza(n, m, b);
	printf("\n\t-- Matriz Triangularizada --\n");
	imprimeM(n, m, b);

	float *x = substituicaoRegressiva(n, m, b);
	printf("\n\t-- Vetor Resposta --\n");
	imprimeVetor(n, x);

	free(b);
	free(x);
	printf("\n\t-- Numero de operações: %d\n", op);
}

void seidel(int n, float **m)
{
	float *b = criaB(n, m);
	double tol = 0.0000000001, distRel = 1.0, soma;
	double d[n];
	double dmax, xmax;
	int i, j;
	float *xa = criaXSeidel(n, b, m);

	float *x = criaXSeidel(n, b, m);
	printf("\n\t-- Matriz Inicial --\n");
	imprimeM(n,m,b);

	while (distRel > tol)
	{
		dmax = 0;
		xmax = 0;
		for (i = 0; i < n; i++)
		{
			soma = 0;
			//esquerda do pivor
			j = i-2;
			if(j >= 0){
				soma = soma + m[i][j] * x[j];
			}
			j = i-1;
			if(j >= 0){
				soma = soma + m[i][j] * x[j];
			}
			
			//direita do pivor
			j = i+1;
			if(j < n){
				soma = soma + m[i][j] * xa[j];
			}
			j=i+2;
			if(j < n){
				soma = soma + m[i][j] * xa[j];
			}

			x[i] = (b[i] - soma) / m[i][i];
			//printf("\nx[i]: %f",x[i]);
			if (x[i] > xmax){
				xmax = x[i];
			}
			
			d[i] = fabs(xa[i] - x[i]);
			//printf("\nd[i]: %f",d[i]);
			
			if (d[i] > dmax){
				dmax = d[i];
			}

		} // fim for i
		distRel = dmax / xmax;
		printf("%lf ", distRel);
	}
	
	printf("\n\t-- Vetor Resposta --\n");
	imprimeVetor(n, x);
	free(b);
	free(xa);
	free(x);
}

//
void triangulariza(int n, float **m, float *b)
{
	op = 0;
	int indicePivo, auxInt;
	float aux, maior, mult;

	for (int k = 0; k < (n - 1); k++)
	{
		indicePivo = k;
		maior = fabs(m[k][k]);
		for (int i = k + 1; i <= k + 2; i++)
		{
			if (k >= n - 2)
			{
				if (k == n - 2)
				{
					if (fabs(m[i][k]) > maior)
					{
						maior = fabs(m[i][k]);
						indicePivo = i;
					}
					i++;
				}
			}

			else
			{
				if (fabs(m[i][k]) > maior)
				{
					maior = fabs(m[i][k]);
					indicePivo = i;
				}
			}
		}
		//troca linha pela maior
		if (indicePivo != k)
		{
			for (int j = k; j <= n - 1; j++)
			{
				aux = m[k][j];
				m[k][j] = m[indicePivo][j];
				m[indicePivo][j] = aux;
			}
			aux = b[k];
			b[k] = b[indicePivo];
			b[indicePivo] = aux;
		}
		for (int i = (k + 1); i <= k + 2; i++)
		{
			auxInt = i;
			if (k >= n - 2)
			{
				if (k == n - 2)
				{
					mult = m[i][k] / m[k][k];
					i += 2;
				}
			}
			else
			{
				mult = m[i][k] / m[k][k];
			}
			op++;
			m[auxInt][k] = 0; // para visualização da matriz triangularizada
			printf("m linha %d = %f; ",i, mult);

			//for(j=(k+1);j<n;j++)
			for (int j = k + 1; j < n; j++)
			{
				m[auxInt][j] = m[auxInt][j] - mult * m[k][j];
				op = op + 2;
			} // fim  j
			b[auxInt] = b[auxInt] - mult * b[k];
			op = op + 2;
		} // fim linha i

		// Mostrando a matriz intermediaria
		printf("\n\n\t-- Matriz -> Etapa %d --\n", k);
		imprimeM(n, m, b);
		printf("\n");
	}
}

float *substituicaoRegressiva(int n, float **m, float *b)
{
	float *x = criaX(n);
	float soma;

	x[n - 1] = b[n - 1] / m[n - 1][n - 1];
	for (int i = (n - 2); i >= 0; i--)
	{
		soma = b[i];
		for (int j = i + 1; (j < n) && (j < (i + 5)); j++)
		{
			soma = soma - m[i][j] * x[j];
		}
		x[i] = soma / m[i][i];
	}

	return x;
}

// ----------------------------------------   Auxiliares  ---------------------------------------------

void leMatriz(float **m, int n, float dA, float dB, float dP, float dC, float dD)
{
	int i, j, sub;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			sub = i - j;
			switch (sub)
			{
			case -2:
				m[i][j] = dD;
				break;
			case -1:
				m[i][j] = dC;
				break;
			case 0:
				m[i][j] = dP;
				break;
			case 1:
				m[i][j] = dB;
				break;
			case 2:
				m[i][j] = dA;
				break;
			default:
				m[i][j] = 0;
			}
		}
	}
}

float **criaMatriz(int n)
{
	float **m = (float **)malloc(n * sizeof(float *));
	for (int i = 0; i < n; i++)
	{
		m[i] = (float *)malloc(n * sizeof(float));
	}
	return m;
}

float *criaB(int n, float **m)
{
	float *b = (float *)malloc(n * sizeof(float));

	for (int i = 0; i < n; i++)
	{
		b[i] = 0;
		for (int j = 0; j < n; j++)
		{
			b[i] = b[i] + m[i][j];
		}
	}
	return b;
}

float *criaX(int n)
{
	float *b = (float *)malloc(n * sizeof(float));

	for (int i = 0; i < n; i++)
	{
		b[i] = 0;
	}
	return b;
}

float *criaXSeidel(int n, float *b, float **m)
{
	float *x = (float *)malloc(n * sizeof(float));

	for (int i = 0; i < n; i++)
	{
		x[i] = b[i] / m[i][i];
	}

	return x;
}
void imprimeM(int n, float **m, float *b)
{
	printf("\n");
	for (int i = 0; i < n; i++)
	{
		printf(" ");
		for (int j = 0; j < n; j++)
		{
			printf(" %.3f ", m[i][j]);
		}
		printf("| %.3f ", b[i]);
		printf("\n");
	}
}

void imprimeVetor(int n, float *a)
{
	printf("\n--- Vetor ---\n");
	for (int i = 0; i < n; i++)
	{
		printf("| %f |", a[i]);
		printf("\n");
	}
	printf("\n");
}

void liberaMatriz(int n, float **m)
{
	for (int i = 0; i < n; i++)
	{
		free(m[i]);
	}
	free(m);
}
