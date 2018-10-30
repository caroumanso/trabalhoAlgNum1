#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int op = 0;

void imprimeM(int n, float** m, float* b);
void leMatriz(float**m,int n,float dA,float dB,float dP,float dC,float dD);
float** criaMatriz(int n);
void liberaMatriz(int n, float** m);
float* criaB(int n, float**m);
void triangulariza(int n, float** m, float* b);
float* substituicaoRegressiva(int n, float** m , float* b);
void imprimeVetor(int n, float* a);
void gauss(int n, float** m);
void seidel(int n,)

int main(){
	int n;
	float** m;
	float *b, *x;
	float dA,dB, dP, dC, dD;	
	
	
	printf("Insira o valor de N: ");
	scanf("%d",&n);
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
	leMatriz(m,n,dA,dB,dP,dC,dD);
	
	//
	gauss(n,m);
	
	
	leMatriz(m,n,dA,dB,dP,dC,dD);
	//
	seidel(n,m);
	
	liberaMatriz(n,m);
	return 0; 
}



//
void gauss(int n, float** m){
	
	b = criaB(n, m);
	
	imprimeM(n, m, b);
	
	triangulariza(n, m, b);
	imprimeM(n, m, b);
	
	x = substituicaoRegressiva(n, m, b);
	imprimeVetor(n,x);
	
	free(b);
	free(x);
	printf("\nNumero de operações: %d\n",op);
}



/*
void seidel(int n, float** m){
	float tol = 0.000001;
    float distRel = 1.0;
    float d[N];
    int dmax,xmax;
    
    while(distRel > tol){
      iter = 0;
      iter++;
      dmax = 0;
      xmax = 0;
      for(i=0;i< n;i++){
		  soma=0;
		  for(j=0;j<i ;j++){
		      soma = soma + A[i][j]*X[j];
		  }
		  printf("\n Vetor1");
		  for(j=(i+1);j<n;j++){
		      soma = soma + A[i][j]*XA[j];
		  }
		  X[i]= (b[i]-soma )/A[i][i];
		  if(X[i] > xmax){
		    xmax = X[i];
		  }
		  printf("\n Vetor2");
		  d[i]=fabs(XA[i]-X[i]);        
		  if(d[i]>dmax){
		  	dmax = d[i];
		  }
     } // fim for i
  	distRel = dmax/xmax;
}*/

//
void triangulariza(int n, float** m, float* b){
	op = 0;	
	int indicePivo, auxInt;
	float aux, maior, mult;

	for(int k=0;k<(n-1);k++)
    {
       indicePivo = k;
       maior = fabs(m[k][k]);
       for(int i=k+1; i<= k+2;i++){
			if(k >= n-2){
				if(k == n-2){
					if(fabs(m[i][k]) > maior){
						maior = fabs(m[i][k]);
						indicePivo = i;
					}
					i++;
				}
					
			}
				
			else{
		   		if(fabs(m[i][k]) > maior){
		   			maior = fabs(m[i][k]);
		   			indicePivo = i;
		   		}
			}
       } 
       //troca linha pela maior
       if(indicePivo != k){
       		for(int j = k; j <= n-1; j++){
       			aux = m[k][j];
       			m[k][j] = m[indicePivo][j];
       			m[indicePivo][j] = aux;
       		}
       		aux = b[k];
       		b[k] = b[indicePivo];
       		b[indicePivo] = aux;
       }
	
   	 	for(int i=(k+1);i<= k+2;i++){
   	 		auxInt = i;
			if(k >= n-2){
				if(k == n-2){
					 mult = m[i][k]/m[k][k];
					i += 2;
				}
			}
			else{
           		mult = m[i][k]/m[k][k];
			}
			op++;
        	m[auxInt][k]=0; // para visualização da matriz triangularizada
        	printf("   m = %f;", mult);
         
           //for(j=(k+1);j<n;j++)
        	for(int j=k+1;j<n;j++)
        	{
        		m[auxInt][j]= m[auxInt][j] - mult*m[k][j];
        		op = op +2;
        	} // fim  j
           b[auxInt]= b[auxInt]- mult*b[k];
				op = op + 2;
        } // fim linha i
	
      // Mostrando a matriz intermediaria
      printf("\n - Matriz apos a etapa %d ------\n", k);
    	imprimeM(n, m, b);
		printf("\n");
	

    }
}

float* substituicaoRegressiva(int n, float** m , float* b){
	float* x = criaX(n);
	float soma;
	
	x[n-1]= b[n-1]/m[n-1][n-1]; 
    for(int i=(n-2);i>=0;i--)  
    {
      soma=b[i];
      for(int j=i+1;(j<n)&&(j<(i+5));j++)
      {
        soma = soma - m[i][j]*x[j];
      }
      x[i]= soma/m[i][i];
      printf(" %f\n",soma);
    }
    
    return x;
}

// ----------------------------------------   Auxiliares  ---------------------------------------------

void leMatriz(float**m,int n,float dA,float dB,float dP,float dC,float dD){
	int i,j, sub;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			sub = i-j;
			switch (sub){
				case -2 : m[i][j] = dD; break;
				case -1 : m[i][j] = dC; break;
				case 0 : m[i][j] = dP; break;
				case 1 : m[i][j] = dB; break;
				case 2 : m[i][j] = dA; break;
				default : m[i][j] = 0;
			}
		}
	}
}

float** criaMatriz(int n){
	float** m = (float**) malloc (n*sizeof(float*));
	for(int i = 0; i<n;i++){
		m[i] = (float*) malloc (n*sizeof(float));
	}
	return m;
}

float* criaB(int n, float**m){
	float* b = (float*) malloc (n*sizeof(float));
	
	for(int i=0;i<n;i++)
    {
      b[i]=0;
      for(int j=0;j< n;j++)
      {
           b[i]= b[i] + m[i][j];
      }
      //printf("%.1f",b[i]);
	}
	//imprimeM(n,m,b);
	return b;
}

float* criaX(int n){
	float* b = (float*) malloc (n*sizeof(float));
	
	for(int i=0;i<n;i++)
    {
      b[i]=0;
	}
	return b;
}

void imprimeM(int n, float** m, float* b){
	printf("\n---------------Matriz---------------\n"); 
	for(int i = 0;i<n;i++){
		printf(" "); 
		for(int j=0;j<n;j++){
			printf(" %.3f ",m[i][j]);
		}
		printf("| %.3f ",b[i]);
		printf("\n");
		 
	}
}

void imprimeVetor(int n, float* a){
	for(int i=0;i<n;i++){
		printf("| %f |",a[i]);
		printf("\n");
	}
	printf("\n");
}

void liberaMatriz(int n, float** m){
	for(int i = 0; i<n;i++){
		free(m[i]);
	}
	free(m);
}


