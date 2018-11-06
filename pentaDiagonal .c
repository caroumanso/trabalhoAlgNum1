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
void erro(int n, float*a);

int main(){
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
void gauss(int n, float **m){
	float *b = criaB(n, m);
	
	printf("\n\n\t01 - Iniciando Eliminacao de Gauss\n");
	printf("\n\t-- Matriz Inicial --\n");
	if(n<=20){
		imprimeM(n, m, b);
	}else{
		printf("\t\t - Matriz muito grande.");
	}

	triangulariza(n, m, b);
	float *x = substituicaoRegressiva(n, m, b);
	
	printf("\n\t-- Matriz Final --\n");
	if(n<=20){
		imprimeM(n, m, b);
	}else{
		printf("\t\t - Matriz muito grande.");
	}

	printf("\n\t01.1 - Gravando Vetor Resposta\n");
	if(n<=20){
		imprimeVetor(n, x);
	}else{
		printf("\t\t - Matriz muito grande.");
	}
	
	printf("\n\t01.2 - Numero de operações em Eliminacao de Gauss: %d\n", op);
	printf("\n\t01.3 - Analizando Erro Gauss: \n");
	erro(n, x);
	free(b);
	free(x);
}

void seidel(int n, float **m){
	float *b = criaB(n, m);
	float tol = 0.0000000001, distRel = 1.0, soma;
	float d[n];
	float dmax, xmax, max;
	int i, j, it=0;
	
	float *xa = criaXSeidel(n, b, m);
	float *x = criaXSeidel(n, b, m);

	op = 0;

	printf("\n\t02 - Iniciando Metodo de Gauss-Seidel\n");
	printf("\n\t-- Matriz Inicial --\n");
	if(n<=20){
		imprimeM(n, m, b);
	}else{
		printf("\t\t - Matriz muito grande.");
	}
	
	while (tol < distRel){
		it++;
		dmax = 0;
		xmax = 0;
		float xamax = 0;
		for (i = 0; i < n; i++){
			soma = 0;
			d[i]=0;
			//esquerda do pivo
			for(j = i-2; j < i; j++){
				if(j >= 0){
					soma = soma + (m[i][j] * x[j]);
					op += 2;
				}
			}
			//direita do pivo
			for(j = i+1; j <= i+2; j++){
				if(j < n){
					soma = soma + (m[i][j] * xa[j]);
					op += 2;				
				}
			}
			x[i] = (b[i] - soma) / m[i][i];
			op += 2;

			max = fabs(x[i]);
			if ( max > xmax ){
				xmax = max;
			}

			d[i] = fabs(x[i] - xa[i]);
			op++;
			if (d[i] > dmax){
				dmax = d[i];
			}
		}
		distRel = dmax / xmax;
		op++;
		for(i = 0; i < n; i++){
			xa[i]=x[i];
		}
	}

	printf("\n\t-- Matriz Final --\n");
	if(n<=20){
		imprimeM(n, m, b);
	}else{
		printf("\t\t - Matriz muito grande.");
	}
	

	printf("\n\t02.1 - Gravando Vetor Resposta\n");
	if(n<=20){
		imprimeVetor(n, x);
	}else{
		printf("\t\t - Matriz muito grande.");
	}
	
	printf("\n\t02.2 - Numero de operações em Gauss-Seidel: %d   - Numero de Iterações: %d\n", op,it);
	printf("\n\t02.3 - Analizando Erro de Seidel: \n");
	erro(n, x);
	
	free(b);
	free(xa);
	free(x);
}

//
void triangulariza(int n, float **m, float *b){
	op = 0;
	int indicePivo, auxInt;
	float aux, maior, mult;

	for (int k = 0; k < (n - 1); k++){
		indicePivo = k;
		maior = fabs(m[k][k]);
		//encontra o pivo
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
			

			
			for (int j = k + 1; j < n; j++){
				m[auxInt][j] = m[auxInt][j] - mult * m[k][j];
				op += 2;
			} 
			b[auxInt] = b[auxInt] - mult * b[k];
			op += 2;
		} 
	}
}

float *substituicaoRegressiva(int n, float **m, float *b){
	float *x = criaX(n);
	float soma;

	x[n - 1] = b[n - 1] / m[n - 1][n - 1];
	op++;
	for (int i = (n - 2); i >= 0; i--){
		soma = b[i];
		for (int j = i + 1; (j < n) && (j < (i + 5)); j++){
			soma = soma - (m[i][j] * x[j]);
			op+=2;
		}
		x[i] = soma / m[i][i];
		op++;
	}

	return x;
}

void erro(int n, float* a){
	float emax = 0, emedio = 0, e;
	for(int i=0; i<n;i++){
		e = fabs(a[i] - 1);
		if(e > emax){
			emax = e;
		}
		emedio += e;
	}
	emedio /= n;
	printf("\n\tErro Max: %.15f\n\tErro Medio: %.15f\n\n",emax,emedio);
}

// ----------------------------------------   Auxiliares  ---------------------------------------------

void leMatriz(float **m, int n, float dA, float dB, float dP, float dC, float dD){
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

float **criaMatriz(int n){
	float **m = (float **)malloc(n * sizeof(float *));
	for (int i = 0; i < n; i++){
		m[i] = (float *)malloc(n * sizeof(float));
	}
	return m;
}

float *criaB(int n, float **m){
	float *b = (float *)malloc(n * sizeof(float));

	for (int i = 0; i < n; i++){
		b[i] = 0;
		for (int j = 0; j < n; j++){
			b[i] = b[i] + m[i][j];
		}
	}
	return b;
}

float *criaX(int n){
	float *b = (float *)malloc(n * sizeof(float));

	for (int i = 0; i < n; i++){
		b[i] = 0;
	}
	return b;
}

float *criaXSeidel(int n, float *b, float **m){
	float *x = (float *)malloc(n * sizeof(float));

	for (int i = 0; i < n; i++){
		x[i] = b[i] / m[i][i];
	}

	return x;
}
void imprimeM(int n, float **m, float *b){
	printf("\n");
	for (int i = 0; i < n; i++){
		printf(" ");
		for (int j = 0; j < n; j++){
			printf(" %.3f ", m[i][j]);
		}
		printf("| %.3f ", b[i]);
		printf("\n");
	}
}

void imprimeVetor(int n, float *a){
	printf("\n\t---    Vetor    ---\n\t");
	for (int i = 0; i < n; i++){
		if(i<10){
			printf("| x[0%d] = %.7f |", i, a[i]);
		}else{
			printf("| x[%d] = %.7f |", i, a[i]);
		}
		printf("\n\t");
	}
	printf("\n");
}

void liberaMatriz(int n, float **m){
	for (int i = 0; i < n; i++){
		free(m[i]);
	}
	free(m);
}
