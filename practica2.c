#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define N 512
float Mat[N][N];
float MatDD[N][N];
float V1[N];
float V2[N];
float V3[N];
float V4[N];


void InitData(){
	int i,j;
	srand(4422543);
	for( i = 0; i < N; i++ )
 		for( j = 0; j < N; j++ ){Mat[i][j]=(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX)));
 	if ( (abs(i - j) <= 3) && (i != j))
 		MatDD[i][j] = (((i*j)%3) ? -1 : 1)*(rand()/(1.0*RAND_MAX));
 	else if ( i == j )
 		MatDD[i][j]=(((i*j)%3)?-1:1)*(10000.0*(rand()/(1.0*RAND_MAX)));
 	else MatDD[i][j] = 0.0;
 	}
for( i = 0; i < N; i++ ){
	V1[i]=(i<N/2)?(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX))):0.0;
 	V2[i]=(i>=N/2)?(((i*j)%3)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX))):0.0;
 	V3[i]=(((i*j)%5)?-1:1)*(100.0*(rand()/(1.0*RAND_MAX)));
	}
}
void PrintVect( float vect[N], int from, int numel ){
	int i;
	for( i = from; i < from + numel; i++){
        	printf("%f ",vect[i]);
	}
	printf("\n");
}
void PrintRow( float mat[N][N], int row, int from, int numel ){
	int i,j;
	i = row;
	for(j = from; j < from + numel; j++){
        	printf("%f ",mat[i][j]);
        }
	printf("\n");
}
void MultEscalar( float vect[N], float vectres[N], float alfa ){
	int i;
	for( i=0; i<N; i++){
		vectres[i]=vect[i]*alfa;
	}
}
float Scalar( float vect1[N], float vect2[N] ){
        int i;
        float Escalar;
        Escalar = 0;
        for( i=0; i < N; i++){
                Escalar+=vect1[i]*vect2[i];
        }
        printf("\n");
        return Escalar;
}
float Magnitude( float vect[N] ){
	return sqrtf(Scalar(vect,vect));
}
int Ortogonal(float vect1[N],float vect2[N]){
	if (Scalar(vect1,vect2) == 0 )
		return 1;
	else return 0;
}
void Projection( float vect1[N], float vect2[N], float vectres[N] ){
	float val1=Scalar(vect1,vect2);
	float val2=Magnitude(vect2);
	float val3=val1/val2;
	MultEscalar(vect2,vectres,val3);
}
float Infininorm( float M[N][N] ){
	int i,j;
	float suma=0;
	float max=0;
	float abs=0;
	for ( i=0;i<N;i++ ){
		for ( j=0;j<N;j++ ){
			abs=fabs(M[i][j]);
			suma=suma+abs;
		}
		if (suma>max)
			max=suma;
		suma=0;
	}
	return max;
}
float Onenorm( float M[N][N] ){
	int i,j;
        float suma=0;
        float max=0;
        float abs=0;
        for ( i=0;i<N;i++ ){
                for ( j=0;j<N;j++ ){
                        abs=fabs(M[j][i]);
                        suma=suma+abs;
                }
               	if (suma>max)
                        max=suma;
                suma=0;
        }
        return max;
}
float NormFrobenius( float M[N][N] ){
	int i,j;
        float suma=0;
        float max=0;
        float abs=0;
        for ( i=0;i<N;i++ ){
                for ( j=0;j<N;j++ ){
                        abs=fabs(M[i][j])*fabs(M[i][j]);
                        suma=suma+abs;
                }
	}
	return sqrtf(suma);
}
int DiagonalDom( float M[N][N] ){
	int i,j;
        float suma=0;
        float abs=0;
	int var=1;
        for ( i=0;i<N;i++ ){
                suma=0;
		for ( j=0;j<N;j++ ){
			if ( j!=i )
                        	abs=fabs(M[i][j]);
                        	suma=suma+abs;
                }
		if ( fabs(M[i][i])<suma )
			var=0;	
        }
	return var;
}
int Jacobi( float M[N][N] , float vect[N], float sol[N], unsigned iter ){
	int i,j,k;
	float temp[N];
	if (DiagonalDom(M)==0.0) 
		return 0;
	else
		for (j=0;j<N;j++){
			sol[j]=0;
		}
		for (i=0;i<iter;i++){
			for (j=0;j<N;j++){
				temp[j]=vect[j];
				for (k=0;k<N;k++ ){
					if ( j!=k )
						temp[j]-=M[k][j]*sol[k];		
				}
				temp[j]/=M[j][j];
			}
		
			for (j=0;j<N;j++){
				sol[j]=temp[j];
			}
		return 1;		
		}	
	

}
int main(){
        InitData();
        printf("V1 del 0 al 9 i del 256 al 265:\n");
        PrintVect(V1,0,10);
        PrintVect(V1,256,10);
        printf("\nV2 del 0 al 9 i del 256 al 265:\n");
        PrintVect(V2,0,10);
        PrintVect(V2,256,10);
        printf("\nV3 del 0 al 9 i del 256 al 265:\n");
        PrintVect(V3,0,10);
        PrintVect(V3,256,10);
        printf("\nMat fila 0 i fila 100 del 0 al 9:\n");
        PrintRow(Mat,0,0,10);
        PrintRow(Mat,100,0,10);
        printf("\nMatDD fila 0 del 0 al 9 i fila 100 del 95 al 104:\n");       PrintRow(MatDD,0,0,10);
        PrintRow(MatDD,100,95,10);
        printf("\nInfininorma de Mat = %f",Infininorm(Mat));
        printf("\nNorma ú de Mat = %f",Onenorm(Mat));
        printf("\nNorma de Frobenius de Mat = %f", NormFrobenius(Mat));
        if(DiagonalDom(Mat)) printf("\nLa matriu Mat és diagonal dominant\n");
        else printf("\nLa matriu Mat no és diagonal dominant\n");
        printf("\nInfininorma de MatDD = %f",Infininorm(MatDD));
        printf("\nNorma ú de MatDD = %f",Onenorm(MatDD));
        printf("\nNorma de Frobenius de MatDD = %f", NormFrobenius(MatDD));
        if(DiagonalDom(MatDD)) printf("\nLa matriu MatDD és diagonal dominant\n");
        else printf("\nLa matriu MatDD no és diagonal dominant\n");
        printf("\nEscalar <V1,V2> = %f Escalar <V1,V3> = %f Escalar <V2,V3> = %f\n",Scalar(V1,V2),Scalar(V1,V3),Scalar(V2,V3));
        printf("\nMagnitud V1,V2 i V3 = %f %f %f\n", Magnitude(V1),Magnitude(V2),Magnitude(V3));
        if(Ortogonal(V1,V2)) printf("\nV1 i V2 són ortogonals");
        else printf("\nV1 i V2 no són ortogonals");
        if(Ortogonal(V1,V3)) printf("\nV1 i V3 són ortogonals");
        else printf("\nV1 i V3 no són ortogonals");
        if(Ortogonal(V2,V3)) printf("\nV2 i V3 són ortogonals");
        else printf("\nV2 i V3 no són ortogonals");
        printf("\n\nEls elements 0 al 9 i 256 al 265 del resultat de multiplicar V3x2.0 són:\n");
        MultEscalar(V3,V4,2);
        int i;
        for( i = 0; i<10; i++){
                printf("%f ",V4[i]);
        }
        printf("\n");
        for( i = 256; i<266; i++){
                printf("%f ",V4[i]);
        }
        printf("\nEls elements 0 a 9 del resultat de la projecció de V2 sobre V3 són:\n");
        Projection (V2, V3, V4);
        for( i = 0; i<10; i++){
                printf("%f ",V4[i]);
        }
        printf("\nEls elements 0 a 9 del resultat de la projecció de V1 sobre V2 són:\n");
        Projection (V1, V2, V4);
	for( i = 0; i<10; i++){
                printf("%f ",V4[i]);
        }
        int trobat = Jacobi(MatDD, V3, V4, 1);
        if( trobat==0 ) printf("\nLa matriu MatDD no és diagonal dominant, no es pot aplicar Jacobi\n");
        else{
                printf("\n\nEls elements 0 a 9 de la solució (1 iter) del sistema d'equacions són:\n");
                for( i = 0; i<10; i++){
                        printf("%f ",V4[i]);
                }
                printf("\nEls elements 0 a 9 de la solució (1000 iters) del sistema d'equacions són:\n");
                trobat = Jacobi(MatDD, V3, V4, 1000);
                for( i = 0; i<10; i++){
                        printf("%f ",V4[i]);
                }
        }
        trobat = Jacobi(Mat, V3, V4, 1);
        if( trobat==0 ) printf("\nLa matriu Mat no és diagonal dominant, no es pot aplicar Jacobi\n");
        else{
                printf("\n\nEls elements 0 a 9 de la solució (1 iter) del sistema d'equacions són:\n");
                for( i = 0; i<10; i++){
                        printf("%f ",V4[i]);
                }
                printf("\nEls elements 0 a 9 de la solució (1000 iters) del sistema d'equacions són:\n");
                trobat = Jacobi(Mat, V3, V4, 1000);
                for( i = 0; i<10; i++){
                        printf("%f ",V4[i]);
                }
        }
        printf("\n");
}
         
