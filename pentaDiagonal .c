#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int op = 0;

void imprimeM(int n, double **m, double *b);
void leMatriz(double **m, int n, double dA, double dB, double dP, double dC, double dD);
double **criaMatriz(int n);
void liberaMatriz(int n, double **m);
double *criaB(int n, double **m);
void triangulariza(int n, double **m, double *b);
double *substituicaoRegressiva(int n, double **m, double *b);
void imprimeVetor(int n, double *a);
void gauss(int n, double **m);
void seidel(int n, double **m);
double *criaX(int n);
double *criaXSeidel(int n, double *b, double **m);
void erro(int n, double*a);

int main(){
	int n;
	double **m;
	double *b, *x;
	double dA, dB, dP, dC, dD;

	printf("Insira o valor de N: ");
	scanf("%d", &n);
	printf("Digite o valor de dA: ");
	scanf("%lf", &dA);
	printf("Digite o valor de dB: ");
	scanf("%lf", &dB);
	printf("Digite o valor de dP: ");
	scanf("%lf", &dP);
	printf("Digite o valor de dC: ");
	scanf("%lf", &dC);
	printf("Digite o valor de dD: ");
	scanf("%lf", &dD);

	m = criaMatriz(n);

	leMatriz(m, n, dA, dB, dP, dC, dD);
	gauss(n, m);

	leMatriz(m, n, dA, dB, dP, dC, dD);
	seidel(n, m);

	liberaMatriz(n, m);
	return 0;
}

//
void gauss(int n, double **m){

	double *b = criaB(n, m);
	double *erro = (double*)malloc(sizeof(double));
	printf("\n\t01 - Iniciando Eliminacao de Gauss\n");
	/*printf("\n\t-- Matriz Inicial --\n");
		
	imprimeM(n, m, b);*/

	triangulariza(n, m, b);
	/*printf("\n\t-- Matriz Triangularizada --\n");
	imprimeM(n, m, b);*/

	double *x = substituicaoRegressiva(n, m, b);
	printf("\n\t01.1 - Gravando Vetor Resposta\n");
	imprimeVetor(n, x);
	
	printf("\n\t01.1 - Gravando Vetor Resposta\n");
	erro(n, x);

	free(b);
	free(x);
	printf("\n\t01.2 - Numero de operações em Eliminacao de Gauss: %d\n", op);
}

void seidel(int n, double **m){
	double *b = criaB(n, m);
	double tol = 0.0000000001, distRel = 1.0, soma;
	double d[n];
	double dmax, xmax;
	int i, j, it=0;
	
	double *xa = criaXSeidel(n, b, m);
	double *x = criaXSeidel(n, b, m);

	op = 0;

	printf("\n\t02 - Iniciando Eliminicao de Gauss-Seidel\n");
	
	while ((distRel > tol)){
		it++;
		dmax = 0;
		xmax = 0;
		for (i = 0; i < n; i++){
			soma = 0;
			//esquerda do pivo
			for(j = i-2; j <= i-1; j++){
				if(j >= 0){
					soma = soma + (m[i][j] * x[j]);
					op += 2;
				}
			}
			//direita do pivo
			for(j = i+1; j <= i+2; j++){
				if(j >= 0){
					soma = soma + (m[i][j] * xa[j]);
					op += 2;				
				}
			}
			x[i] = (b[i] - soma) / m[i][i];
			op += 2;
			if (fabs(x[i]) > xmax){
				xmax = fabs(x[i]);

			}
			
			d[i] = fabs(xa[i] - x[i]);
			op++;
			
			if (d[i] > dmax){
				dmax = d[i];
			}

		} // fim for i
		distRel = dmax / xmax;
		op++;
		for(i = 0; i < n; i++){
			xa[i]=x[i];
		}
	}

	printf("\n\t02.1 - Vetor Resposta\n");
	imprimeVetor(n, x);
	free(b);
	free(xa);
	free(x);
	
	printf("\n\t02.2 - Numero de operações em Gauss-Seidel: %d\n", op);
}

//
void triangulariza(int n, double **m, double *b)
{
	op = 0;
	int indicePivo, auxInt;
	double aux, maior, mult;

	for (int k = 0; k < (n - 1); k++){
		indicePivo = k;
		maior = fabs(m[k][k]);
		for (int i = k + 1; i <= k + 2; i++){
			if (k >= n - 2){
				if (k == n - 2){
					if (fabs(m[i][k]) > maior){
						maior = fabs(m[i][k]);
						indicePivo = i;
					}
					i++;
				}
			}

			else{
				if (fabs(m[i][k]) > maior){
					maior = fabs(m[i][k]);
					indicePivo = i;
				}
			}
		}
		//troca linha pela maior
		if (indicePivo != k){
			for (int j = k; j <= n - 1; j++){
				aux = m[k][j];
				m[k][j] = m[indicePivo][j];
				m[indicePivo][j] = aux;
			}
			aux = b[k];
			b[k] = b[indicePivo];
			b[indicePivo] = aux;
		}
		for (int i = (k + 1); i <= k + 2; i++){
			auxInt = i;
			if (k >= n - 2)
			{
				if (k == n - 2){
					mult = m[i][k] / m[k][k];
					i += 2;
				}
			}
			else{
				mult = m[i][k] / m[k][k];
			}
			op++;
			m[auxInt][k] = 0; // para visualização da matriz triangularizada
			//printf("m linha %d = %lf; ",i, mult);

			
			for (int j = k + 1; j < n; j++){
				m[auxInt][j] = m[auxInt][j] - mult * m[k][j];
				op = op + 2;
			} // fim  j
			b[auxInt] = b[auxInt] - mult * b[k];
			op = op + 2;
		} // fim linha i

		// Mostrando a matriz intermediaria
		//printf("\n\n\t-- Matriz -> Etapa %d --\n", k);
		//imprimeM(n, m, b);
		//printf("\n");
	}
}

double *substituicaoRegressiva(int n, double **m, double *b)
{
	double *x = criaX(n);
	double soma;

	x[n - 1] = b[n - 1] / m[n - 1][n - 1];
	for (int i = (n - 2); i >= 0; i--){
		soma = b[i];
		for (int j = i + 1; (j < n) && (j < (i + 5)); j++){
			soma = soma - m[i][j] * x[j];
		}
		x[i] = soma / m[i][i];
	}

	return x;
}

void erro(int n, double* a){
	double emax = 0, emedio = 0, e;
	for(int i; i<n;i++){
		e = fabs(a[i] - 1);
		if(e > emax){
			emax = e;
		}
		emedio += e;
	}
	emedio /= n;
	printf("\nErro Max: %lf\nErro Medio: %lf\n\n",emax,emedio);
}

// ----------------------------------------   Auxiliares  ---------------------------------------------

void leMatriz(double **m, int n, double dA, double dB, double dP, double dC, double dD)
{
	int i, j, sub;
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			sub = i - j;
			switch (sub){
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

double **criaMatriz(int n)
{
	double **m = (double **)malloc(n * sizeof(double *));
	for (int i = 0; i < n; i++){
		m[i] = (double *)malloc(n * sizeof(double));
	}
	return m;
}

double *criaB(int n, double **m)
{
	double *b = (double *)malloc(n * sizeof(double));

	for (int i = 0; i < n; i++){
		b[i] = 0;
		for (int j = 0; j < n; j++){
			b[i] = b[i] + m[i][j];
		}
	}
	return b;
}

double *criaX(int n)
{
	double *b = (double *)malloc(n * sizeof(double));

	for (int i = 0; i < n; i++){
		b[i] = 0;
	}
	return b;
}

double *criaXSeidel(int n, double *b, double **m)
{
	double *x = (double *)malloc(n * sizeof(double));

	for (int i = 0; i < n; i++){
		x[i] = b[i] / m[i][i];
	}

	return x;
}
void imprimeM(int n, double **m, double *b)
{
	printf("\n");
	for (int i = 0; i < n; i++){
		printf(" ");
		for (int j = 0; j < n; j++){
			printf(" %.3lf ", m[i][j]);
		}
		printf("| %.3lf ", b[i]);
		printf("\n");
	}
}

void imprimeVetor(int n, double *a)
{
	printf("\n--- Vetor ---\n");
	for (int i = 0; i < n; i++){
		printf("| %.3lf |", a[i]);
		printf("\n");
	}
	printf("\n");
}

void liberaMatriz(int n, double **m)
{
	for (int i = 0; i < n; i++){
		free(m[i]);
	}
	free(m);
}
