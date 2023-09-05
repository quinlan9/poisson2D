#include<stdio.h>
#include"poisson.h"

int main(){
	FILE *fp1,*fp2;
	double a[maxn][maxn],b[maxn][maxn],c[maxn][maxn];
	double sum,diff=0.0;
	fp1=fopen("grid_outcomes_2d.txt","r");
	//fp2=fopen("grid_outcomes_pscw.txt","r");
	fp2=fopen("grid_outcomes_fence.txt","r");
	for(int i=0;i<maxn;i++){
		for(int j=0;j<maxn;j++){
			fscanf(fp1,"%lf",&a[i][j]);
			fscanf(fp2,"%lf",&b[i][j]);

		}
	}
    	if (fp1 == NULL||fp2==NULL) {
        printf("failed to open file！\n");
        return 1;
    }

	for(int i=0;i<maxn;i++){
                for(int j=0;j<maxn;j++){
      			c[i][j]=a[i][j] - b[i][j];
			diff=c[i][j];
			sum=sum+diff*diff;           
                }
        }
	printf("difference:%lf\n",sum);
	fclose(fp1);
	fclose(fp2);
}
