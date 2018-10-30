#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int op = 0;

void imprimeM(int n, float** m);
void leMatriz(float**m,int n,float dA,float dB,float dP,float dC,float dD);
float** criaMatriz(int n);
void liberaMatriz(int n, float** m);
float* criaB(int n, float**m);
void triangulariza(int n, float** m, float* b);

int main(){
	int n;
	float** m;
	float* b;
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
	b = criaB(n,m);
	leMatriz(m,n,dA,dB,dP,dC,dD);
	imprimeM(n,m);
	triangulariza(n, m, b);
	imprimeM(n, m);
	liberaMatriz(n,m);
	printf("\nNumero de operações: %d\n",op);
	
	return 0; 
}

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
	}
	return b;
}

void imprimeM(int n, float** m){
	printf("\n---------------Matriz---------------\n"); 
	for(int i = 0;i<n;i++){
		printf(" "); 
		for(int j=0;j<n;j++){
			printf("  %.3f ",m[i][j]);
		}
		printf("\n"); 
	}
	//printf("\n--------------- Fim  ---------------\n"); 
}

void liberaMatriz(int n, float** m){
	for(int i = 0; i<n;i++){
		free(m[i]);
	}
	free(m);
}

/*void jacobi(int n, float** m){

}

void seidel(int n, float** m){

}*/

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
						op++; // fabs é uma operação contavel?
					}
					i++;
				}
					
			}
				
			else{
		   		if(fabs(m[i][k]) > maior){
		   			maior = fabs(m[i][k]);
		   			indicePivo = i;
		   			op++;
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
/*	 Impressão da matriz intermediaria
		imprimeM(n, m);
		printf("\n");*/
       

	
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
    	imprimeM(n,m);
		printf("\n");


    }
}
