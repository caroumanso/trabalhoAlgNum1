#include <stdio.h>
#include <stdlib.h>
 
void imprimeM(int n, float** m);
void leMatriz(float**m,int n,float dA,float dB,float dP,float dC,float dD);
float** criaMatriz(int n);
void liberaMatriz(int n, float** m);
float* criaB(int n, float**m);

int main(){
	int n;
	float** m;
	float* b;
	float dA,dB, dP, dC, dD;
	printf("insira um m e os bagulhos...");
	scanf("%d",&n);
	m = criaMatriz(n);
	b = criaB(n,m);
	scanf(" %f %f %f %f %f", &dA, &dB, &dP,&dC, &dD);
	leMatriz(m,n,dA,dB,dP,dC,dD);
	imprimeM(n,m);
	liberaMatriz(n,m);
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
	
	for(i=0;i<n;i++)
    {
      b[i]=0;
      for(j=0;j< n;j++)
      {
           b[i]= b[i] + m[i][j];
      }
	}
	return b;
}

void imprimeM(int n, float** m){
	
	for(int i = 0;i<n;i++){
		printf(" "); 
		for(int j=0;j<n;j++){
			printf("%.2f ",m[i][j]);
		}
		printf("\n"); 
	}
}

void liberaMatriz(int n, float** m){
	for(int i = 0; i<n;i++){
		free(m[i]);
	}
	free(m);
}

void jacobi(int n, float** m){

}

void seidel(int n, float** m){

}

void triangulariza(int n, float** m){
	
	int indicePivo;
	float aux;
	
	for(int k=0;k<(n-1);k++)
    {
       printf("\n--------- Etapa %d -------------\n", k);
       
       //acha maior
       indicePivo = k;
       maior = fabs(A[k][k]);
       for(i=k+1; i<=(n-1);i++){
       		if(fabs(A[i][k]) > maior){
       			maior = fabs(A[i][k]);
       			indicePivo = i;
       		}
       }
       
       //troca linha pela maior
       if(indicePivo != k){
       		for(j=k; j<=n-1;j++){
       			aux = A[k][j];
       			A[k][j] = A[indicePivo][j];
       			A[indicePivo][j] = aux;
       		}
       		aux = b[k];
       		b[k] = b[indicePivo];
       		b[indicePivo] = aux;
       }
        
       for(i=(k+1);i<n;i++)
       {
           m = A[i][k]/A[k][k];
           //A[i][k]=0; // para visualização da matriz triangularizada
           printf("   m = %f;", m);
         
           //for(j=(k+1);j<n;j++)
           for(j=(k);j<n;j++)
           {
              A[i][j]= A[i][j]- m*A[k][j];
           } // fim  j
           b[i]= b[i]- m*b[k];

        } // fim linha i

      // Mostrando a matriz intermediaria
      printf("\n - Matriz apos a etapa %d ------\n", k);
      for(i=0;i< n;i++)
     {
       for(j=0;j< n;j++)
       {
        printf("  %15.6f ", A[i][j]);
        }
        printf(" | b[%d]: %10.3f\n", i, b[i]);
      }


    }
}
