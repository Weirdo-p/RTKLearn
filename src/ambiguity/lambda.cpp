#include "navigation/ambiguity/lambda.h"

int matinv(double *A, int n);


/* new matrix -------------------------------------------------------------*/
extern double *mat(int n, int m)
{
	double *p;

	if (n <= 0 || m <= 0) return NULL;
	if (!(p = (double *)malloc(sizeof(double)*n*m))) 
	{
		printf("matrix memory allocation error: n=%d,m=%d\n", n, m);
	}
	return p;
}

/* new integer matrix ------------------------------------------------------*/
int *imat(int n, int m)
{
	int *p;

	if (n <= 0 || m <= 0) return NULL;
	if (!(p = (int *)malloc(sizeof(int)*n*m))) {
		printf("integer matrix memory allocation error: n=%d,m=%d\n", n, m);
	}
	return p;
}

/* zero matrix -------------------------------------------------------------*/
double *zeros(int n, int m)
{
	double *p;

	if (n <= 0 || m <= 0) return NULL;
	if (!(p = (double *)calloc(sizeof(double), n*m))) 
	{
		printf("matrix memory allocation error: n=%d,m=%d\n", n, m);
	}

	return p;
}

/* identity matrix --------------------------------------------------------*/
double *eye(int n)
{
	double *p;
	int i;

	if ((p = zeros(n, n))) for (i = 0; i < n; i++) p[i + i * n] = 1.0;
	return p;
}

/* multiply matrix -----------------------------------------------------------*/
void matmul(const char *tr, int n, int k, int m, double alpha, const double *A, const double *B, double beta, double *C)
{
	double d;
	int i, j, x, f = tr[0] == 'N' ? (tr[1] == 'N' ? 1 : 2) : (tr[1] == 'N' ? 3 : 4);

	for (i = 0; i < n; i++) for (j = 0; j < k; j++) {
		d = 0.0;
		switch (f) {
		case 1: for (x = 0; x < m; x++) d += A[i + x * n] * B[x + j * m]; break;
		case 2: for (x = 0; x < m; x++) d += A[i + x * n] * B[j + x * k]; break;
		case 3: for (x = 0; x < m; x++) d += A[x + i * m] * B[x + j * m]; break;
		case 4: for (x = 0; x < m; x++) d += A[x + i * m] * B[j + x * k]; break;
		}
		if (beta == 0.0) C[i + j * n] = alpha * d; else C[i + j * n] = alpha * d + beta * C[i + j * n];
	}
}

/* copy matrix --------------------------------------------------------------*/
void matcpy(double *A, const double *B, int n, int m)
{
	memcpy(A, B, sizeof(double)*n*m);
}

/* solve linear equation -----------------------------------------------------*/
int solve(const char *tr, const double *A, const double *Y, int n, int m, double *X)
{
	double *B = mat(n, n);
	int info;

	matcpy(B, A, n, n);
	if (!(info = matinv(B, n))) matmul(tr[0] == 'N' ? "NN" : "TN", n, m, n, 1.0, B, Y, 0.0, X);   //X=B'*Y=E/Z'
	free(B);
	return info;
}

/* LU decomposition ----------------------------------------------------------*/
static int ludcmp(double *A, int n, int *indx, double *d)
{
	double big, s, tmp, *vv = mat(n, 1);
	int i, imax = 0, j, k;

	*d = 1.0;
	for (i = 0; i < n; i++) {
		big = 0.0; for (j = 0; j < n; j++) if ((tmp = fabs(A[i + j * n])) > big) big = tmp;
		if (big > 0.0) vv[i] = 1.0 / big; else { free(vv); return -1; }
	}
	for (j = 0; j < n; j++) {
		for (i = 0; i < j; i++) {
			s = A[i + j * n]; for (k = 0; k < i; k++) s -= A[i + k * n] * A[k + j * n]; A[i + j * n] = s;
		}
		big = 0.0;
		for (i = j; i < n; i++) {
			s = A[i + j * n]; for (k = 0; k < j; k++) s -= A[i + k * n] * A[k + j * n]; A[i + j * n] = s;
			if ((tmp = vv[i] * fabs(s)) >= big) { big = tmp; imax = i; }
		}
		if (j != imax) {
			for (k = 0; k < n; k++) {
				tmp = A[imax + k * n]; A[imax + k * n] = A[j + k * n]; A[j + k * n] = tmp;
			}
			*d = -(*d); vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (A[j + j * n] == 0.0) { free(vv); return -1; }
		if (j != n - 1) {
			tmp = 1.0 / A[j + j * n]; for (i = j + 1; i < n; i++) A[i + j * n] *= tmp;
		}
	}
	free(vv);
	return 0;
}
/* LU back-substitution ------------------------------------------------------*/
static void lubksb(const double *A, int n, const int *indx, double *b)
{
	double s;
	int i, ii = -1, ip, j;

	for (i = 0; i < n; i++) {
		ip = indx[i]; s = b[ip]; b[ip] = b[i];
		if (ii >= 0) for (j = ii; j < i; j++) s -= A[i + j * n] * b[j]; else if (s) ii = i;
		b[i] = s;
	}
	for (i = n - 1; i >= 0; i--) {
		s = b[i]; for (j = i + 1; j < n; j++) s -= A[i + j * n] * b[j]; b[i] = s / A[i + i * n];
	}
}

/* inverse of matrix ---------------------------------------------------------*/
int matinv(double *A, int n) {
	double d, *B;
	int i, j, *indx;

	indx = imat(n, 1); B = mat(n, n); matcpy(B, A, n, n);
	if (ludcmp(B, n, indx, &d)) { free(indx); free(B); return -1; }
	for (j = 0; j < n; j++) {
		for (i = 0; i < n; i++) A[i + j * n] = 0.0; A[j + j * n] = 1.0;
		lubksb(B, n, indx, A + j * n);
	}
	free(indx); free(B);
	return 0;
}


/* LD factorization (Q=L'*diag(D)*L) -----------------------------------------*/
static int LD(int n, const double *Q, double *L, double *D)   //ŽËŽŠ·ÖœâËã·šÍ¬ FMFAC5£š²ÎŒû"Ä£ºý¶ÈŸ­µä°æ.pdf"£© 
{
    int i,j,k,info=0;
    double a,*A=mat(n,n); //mat():žøŸØÕó·ÖÅäÄÚŽæ
    
    memcpy(A,Q,sizeof(double)*n*n); //œ«QÕóžŽÖÆµœAÖÐ
    for (i=n-1;i>=0;i--) {          //FMFAC5ÖÐÊÇŽÓi=n¿ªÊŒµÄ£¬ÎªºÎÓÐÇø±ð£¿ÓÉÓÚŸØÕóÔªËØÊÇŽÓ0µœn*n-1£¬¶øFMFAC5ÖÐÊÇ1µœn*n£¬¹ÊŽËŽŠÐèÒª¶àŒõ1
        if ((D[i]=A[i+i*n])<=0.0) {info=-1; break;}  //D(i,i)=Q(i,i),×¢ÒâŽËŽŠDÊÇÒ»žöÁÐÏòÁ¿¶ø·ÇŸØÕó
        a=sqrt(D[i]);                                //a=sqrt(Q(i,i))
        for (j=0;j<=i;j++) L[i+j*n]=A[i+j*n]/a;      //L(i,1:i)=Q(i,1:i)/sqrt(Q(i,i))
        for (j=0;j<=i-1;j++) for (k=0;k<=j;k++) A[j+k*n]-=L[i+k*n]*L[i+j*n];  //Q(j,1:j)=Q(j,1:j)-L(i,1:j)L(i,j) Í¬ÑùµÄÓÉÓÚŸØÕóÔ­Òò£¬jŽÓ0¿ªÊŒ   
        for (j=0;j<=i;j++) L[i+j*n]/=L[i+i*n];       //L(i,1:i)=L(i,1:i)/L(i,i)
    }
    free(A);
    if (info) printf("%s : LD factorization error\n");
    return info;
}

/* LDLT factorization  ( Q=L*diag(D)*L' ) ------2017.5.9----------------------------------------- */
static int LDLT(int n, const double *Q, double *L, double *D)
{
	int i, j, k, info = 0;
	double  *A = mat(n, n); 

	memcpy(A, Q, sizeof(double)*n*n); 
	for (i = 0; i < n; i++)   L[i + i*n] = 1;
	for (j = 0; j < n - 1; j++)
	{
		if ((D[j] = A[j + j*n]) <= 0.0) { info = -1; break; }
		for (i = 1; i < n; i++) 
			L[i + j*n] = A[i + j*n] / D[j];
		for (k = j + 1; k < n; k++)
			for (i = k; i < n; i++)  
				A[i + k*n] = A[i + k*n] - A[k + j*n] * L[i + j*n];
	}
	D[n - 1] = A[n*n - 1];
	free(A);
	return info;
}

/* integer gauss transformation ----------------------------------------------*/
static void gauss(int n, double *L, double *Z, int i, int j)  //ÀàËÆZTRANËã·š£¬µ«iºÍjµ÷»»ÁË,j±íÊŸÁÐ
{
    int k,mu;
    
    if ((mu=(int)ROUND(L[i+j*n]))!=0) {   //ËÄÉáÎåÈë
        for (k=i;k<n;k++) L[k+n*j]-=(double)mu*L[k+i*n];   //žßË¹ÏûÔª
        for (k=0;k<n;k++) Z[k+n*j]-=(double)mu*Z[k+i*n];   //ZŸØÕó×öÏàÍ¬µÄ±ä»¯
    }
}

/* permutationsœ»»» --------------------------------------------------------------*/
static void perm(int n, double *L, double *D, int j, double del, double *Z)    //ÌõŒþ·œ²îÅÅÐò
{                                                                              //ÓëSRC1Ëã·šÖÐµÄÊµÏÖÏàÍ¬
    int k;
    double eta,lam,a0,a1;
    
    eta=D[j]/del;
    lam=D[j+1]*L[j+1+j*n]/del;   //lamda3=...
    D[j]=eta*D[j+1]; D[j+1]=del;
    for (k=0;k<=j-1;k++) {     //ÓÉÓÚŸØÕóÔ­Òò£¬³õÖµ±äÎª0
        a0=L[j+k*n]; a1=L[j+1+k*n];   //lamda1=...  lamda2=...
        L[j+k*n]=-L[j+1+j*n]*a0+a1;
        L[j+1+k*n]=eta*a0+lam*a1;
    }
    L[j+1+j*n]=lam;
    for (k=j+2;k<n;k++) SWAP(L[k+j*n],L[k+(j+1)*n]);
    for (k=0;k<n;k++) SWAP(Z[k+j*n],Z[k+(j+1)*n]);
}

/* lambda reduction (z=Z'*a, Qz=Z'*Q*Z=L'*diag(D)*L) (ref.[1]) ---------------*/
static void reduction(int n, double *L, double *D, double *Z)  //ÀàËÆSRC1Ëã·š£¬°ÑžßË¹ÏûÔªºÍÌõŒþ·œ²îÅÅÐò·ÅÔÚÒ»ÆðœøÐÐ£šÄ£ºý¶ÈŸ­µä°æ.pdf£©
{                                                              //ÎÒÈÏÎªžÃËã·šÓëSRC1Ö»ÊÇ¿ŽÆðÀŽ²»Í¬£¬µ«ÔËÐÐ¹ý³ÌÊÇÏàÍ¬µÄ
    int i,j,k;
    double del;
    
    j=n-2; k=n-2;    //ŽËŽŠkÖži1,ÓÉÓÚŸØÕóÔ­Òò±äÎªn-2   jÖži,Í¬ÑùµÄÓÉÓÚŸØÕóÔ­Òò±äÎªn-2
    while (j>=0) {   //while i>1
        if (j<=k) for (i=j+1;i<n;i++) gauss(n,L,Z,i,j);   //ZTRANËã·š
        del=D[j]+L[j+1+j*n]*L[j+1+j*n]*D[j+1];            //ŒÆËãdelta
        if (del+1E-6<D[j+1]) { /* compared considering numerical error */
            perm(n,L,D,j,del,Z);
            k=j; j=n-2;   //ŒŽi1=i£¬iž³³õÖµÒÔ±ãœøÐÐÏÂÒ»ŽÎÑ­»·¡££šÏà±ÈÓÚSRC1Ëã·šÉÙÁËsw²ÎÊý²¢ÇÒÉÙÁËswÄÇÒ»²ãÑ­»·£©
        }
        else j--;
    }
}

/* modified lambda (mlambda) search (ref. [2]) -------------------------------*/
static int search(int n, int m, const double *L, const double *D,       
                  const double *zs, double *zn, double *s)
{                                                            //n:float parameters(nb=16)   m: fixed solutions(=2)   zs:ŒŽz=Z'*a
    int i,j,k,c,nn=0,imax=0;
    double newdist,maxdist=1E99,y;   
    double *S=zeros(n,n),*dist=mat(n,1),*zb=mat(n,1),*z=mat(n,1),*step=mat(n,1);
	//²Î¿ŒLDLT_MSearch.m
	//*zbÎªC_z£¬ŒŽÌõŒþ¹ÀÖµ   *zÎª¶Ô*zsÈ¡Õûºóœá¹û    yÎªzb[k]-z[k]    *znÖÐŽæŽ¢ËÑË÷µœµÄºòÑ¡×é     *sŽæŽ¢ºòÑ¡×é¶ÔÓŠµÄ¿ÕŒäŽóÐ¡
    
    k=n-1; dist[k]=0.0;
    zb[k]=zs[k];  //µÚÒ»žöÌõŒþ¹ÀÖµ³õÊŒ»¯Îªž¡µãœâ
    z[k]=ROUND(zb[k]); y=zb[k]-z[k]; step[k]=SGN(y);     //SGN(x): ((x)<=0.0?-1.0:1.0)   ROUND(x): (floor((x)+0.5))    
    for (c=0;c<LOOPMAX;c++) {         //LOOPMAX: maximum count of search loop  =10000
        newdist=dist[k]+y*y/D[k];
        if (newdist<maxdist) {
            if (k!=0) {                 //ŒÌÐøÏòÏÂÒ»²ãËÑ
                dist[--k]=newdist;
                for (i=0;i<=k;i++)
                    S[k+i*n]=S[k+1+i*n]+(z[k+1]-zb[k+1])*L[k+1+i*n];
                zb[k]=zs[k]+S[k+k*n];
                z[k]=ROUND(zb[k]); y=zb[k]-z[k]; step[k]=SGN(y);
            }
            else {                      //ŽæŽ¢ºòÑ¡×é
                if (nn<m) {
                    if (nn==0||newdist>s[imax]) imax=nn;
                    for (i=0;i<n;i++) zn[i+nn*n]=z[i];
                    s[nn++]=newdist;
                }
                else {
                    if (newdist<s[imax]) {
                        for (i=0;i<n;i++) zn[i+imax*n]=z[i];
                        s[imax]=newdist;
                        for (i=imax=0;i<m;i++) if (s[imax]<s[i]) imax=i;
                    }
                    maxdist=s[imax];
                }
                z[0]+=step[0]; y=zb[0]-z[0]; step[0]=-step[0]-SGN(step[0]);
            }
        }
        else {                           //·µ»ØÉÏÒ»²ã
            if (k==n-1) break;
            else {
                k++;
                z[k]+=step[k]; y=zb[k]-z[k]; step[k]=-step[k]-SGN(step[k]);
            }
        }
    }
    for (i=0;i<m-1;i++) { /* sort by s */      //¶ÔºòÑ¡×éÅÅÐò£¬°ŽÕÕ¶ÔÓŠ¿ÕŒäÓÉÐ¡µœŽóµÄË³Ðò
        for (j=i+1;j<m;j++) {
            if (s[i]<s[j]) continue;
            SWAP(s[i],s[j]);
            for (k=0;k<n;k++) SWAP(zn[k+i*n],zn[k+j*n]); 
        }
    }
    free(S); free(dist); free(zb); free(z); free(step);
    
    if (c>=LOOPMAX) {
        printf("%s : search loop count overflow\n");
        return -1;
    }
    return 0;
}
/* lambda/mlambda integer least-square estimation ------------------------------
* integer least-square estimation. reduction is performed by lambda (ref.[1]),
* and search by mlambda (ref.[2]).
* args   : int    n      I  number of float parameters ,n=nb=32
*          int    m      I  number of fixed solutions ,m=2
*          double *a     I  float parameters (n x 1)
*          double *Q     I  covariance matrix of float parameters (n x n)  ŒŽQb£š²»ÊÇ»ùÏßµÄ£¬¶øÊÇrtklibÀï¶Ô¹Ì¶šœâµÄÃüÃû£©
*          double *F     O  fixed solutions (n x m)  ŒŽb
*          double *s     O  sum of squared residuals of fixed solutions (1 x m)
* return : status (0:ok,other:error)
* notes  : matrix stored by column-major order (fortran convension)
*-----------------------------------------------------------------------------*/
extern int lambda(int n, int m, const double *a, const double *Q, double *F, double *s)
{
    int info;
    double *L,*D,*Z,*z,*E;   
    
    if (n<=0||m<=0) return -1;
    L=zeros(n,n); D=mat(n,1); Z=eye(n); z=mat(n,1),E=mat(n,m);
    
    /* LD factorization */
    if (!(info=LD(n,Q,L,D))) {
        
        /* lambda reduction */
        reduction(n,L,D,Z);   //žßË¹ÏûÔª&ÌõŒþ·œ²îÅÅÐò£¬µÃµœÐÂµÄLDL·ÖœâºÍÕûÊý±ä»»ŸØÕóZ
        matmul("TN",n,1,n,1.0,Z,a,0.0,z); /* z=Z'*a */
        
        /* mlambda search */
        if (!(info=search(n,m,L,D,z,E,s))) {       //ŽËŽŠzÎªÕûÊý±ä»»ºóµÄž¡µãœâzs    EÎªºòÑ¡×ézn      sÎªºòÑ¡×é¶ÔÓŠµÄ¿ÕŒäŽóÐ¡
            
            info=solve("T",Z,E,n,m,F); /* F=E/Z' */      //¶ÔzÄæ±ä»»µÃµœa
        }
    }
    free(L); free(D); free(Z); free(z); free(E);
    return info;    //info=0
}