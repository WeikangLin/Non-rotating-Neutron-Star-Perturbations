# include<math.h>
# include<stdio.h>
# define pi 3.14159265358979
# define G 6.673e-11
# define c 2.99792458e8
# define M 1.99e30
# define SatD 2.68e17
# define l 2
# define Cri 1E-7

static double a[3000],b[3000],n;
static int j=0;

double che(e)
double e;
{
double le;
double p,lp;
int i=0;
le=log10(e);
i=j-1;
for(;i>0;i=i-1)
{
if(le<=a[i]&&le>a[i-1])
  {
  lp=b[i-1]*(le-a[i])/(a[i-1]-a[i])+b[i]*(le-a[i-1])/(a[i]-a[i-1]);
  break;
  }
}
p=pow(10,lp);
return(p);
}

double ch(p)
double p;
{
double lp;
double e,le;
int i=0;
lp=log10(p);
i=j-1;
for(;i>0;i=i-1)
{
if(lp<=b[i]&&lp>b[i-1])
  {
  le=a[i-1]*(lp-b[i])/(b[i-1]-b[i])+a[i]*(lp-b[i-1])/(b[i]-b[i-1]);
  break;
  }
}
e=pow(10,le);
return(e);
}

double f(r,p,e,m)
double r,p,e,m;
{
double g1, g2;
g1=-G*(e+p/c/c)*(m+4*pi*r*r*r*p/c/c);
g2=(r*r-2*G*r*m/c/c);
return(g1/g2);
}

double fm(r,e)
double r,e;
{
double mf;
mf=4*pi*r*r*e;
return(mf);
}

double Bf(r,p,m,B)
double r,p,m,B;
{
double dB, A;
A=1/(1-2*G*m/r/c/c);
dB=2*G/c/c/r/r*(m+4*pi*r*r*r*p/c/c)*A*B;
return(dB);
}

double DH1(r,m,A,p,e,H1,H0,K,V)
double r,m,A,p,e,H1,H0,K,V;
{
double dH1;
dH1=-1/r*(l+1+2*m*A/r+4*pi*r*r*A*(p-e))*H1+1/r*A*(H0+K-16*pi*(e+p)*V);
return(dH1);
}

double DK(r,H0,H1,Dv,K,e,p,A,W)
double r,H0,H1,Dv,K,e,p,A,W;
{
double dK;
dK=1/r*H0+0.5*l*(l+1)/r*H1-((l+1)/r-0.5*Dv)*K-8*pi*(e+p)*sqrt(A)/r*W;
return(dK);
}

double DW(r,W,A,gamma,p,B,X,V,H0,K)
double r,W,A,gamma,p,B,X,V,H0,K;
{
double dW;
dW=-(l+1)/r*W+r*sqrt(A)*(1/gamma/p/sqrt(B)*X-l*(l+1)/r/r*V+0.5*H0+K);
return(dW);
}

double DX(r,X,e,p,B,Dv,H0,w,H1,K,V,A,F,W)
double r,X,e,p,B,Dv,H0,w,H1,K,V,A,F,W;
{
double dX;
dX=-l/r*X+(e+p)*sqrt(B)*(0.5*(1/r-0.5*Dv)*H0+0.5*(r*w*w/B+0.5*l*(l+1)/r)*H1+0.5*(1.5*Dv-1/r)*K-0.5*l*(l+1)*Dv/r/r*V-1/r*(4*pi*(e+p)*sqrt(A)+w*w*sqrt(A)/B-0.5*r*r*F)*W);
return(dX);
}

double H0f(r,B,X,m,p,A,H1,K,w)
double r,B,X,m,p,A,H1,K,w;
{
double H0;
H0=8*pi*r*r*r/sqrt(B)*X-(0.5*l*(l+1)*(m+4*pi*r*r*r*p)-w*w*r*r*r/A/B)*H1+(0.5*(l+2)*(l-1)*r-w*w*r*r*r/B-1/r*A*(m+4*pi*r*r*r*p)*(3*m-r+4*pi*r*r*r*p))*K;
H0=H0/(3*m+0.5*(l+2)*(l-1)*r+4*pi*r*r*r*p);
return(H0);
}

double Vf(r,w,e,p,B,A,Dp,W,H0,X)
double r,w,e,p,B,A,Dp,W,H0,X;
{
double V;
V=1/w/w/(e+p)*B*(1/sqrt(B)*X+1/r*Dp/sqrt(A)*W-0.5*(e+p)*H0);
return(V);
}

double gammaf(e)
double e;
{
double le, p, gamma;
int i=j-1;
e=e/G*c*c;
le=log10(e);
for(;i>0;i=i-1)
  {
  if(le<=a[i]&&le>a[i-1])
    {
    p=b[i-1]*(le-a[i])/(a[i-1]-a[i])+b[i]*(le-a[i-1])/(a[i]-a[i-1]);
    p=pow(10,p);
    gamma=(e+p/c/c)/e*(b[i]-b[i-1])/(a[i]-a[i-1]);
    break;
    }
  }
return(gamma);
}

double Ff(r,A,e,p,m,Dp)
double r,A,e,p,m,Dp;
{
double F;
F=8*pi/r/r*sqrt(A);
F=F*(e+3*p+r*Dp-(m/4/pi/r/r/r+p)*(4+A*(m/r-4*pi*r*r*e)));
return(F);
}



main()
{
FILE *inf, *outf, *mfile, *pfile, *rhofile, *Bfile, *Bcor, *Wfile, *Wfile1, *Wfile2, *Vfile, *Vfile1, *Vfile2;
char filename[20],opfile[20];
double dr=0.05, copy;
double r, r0=1.0, R, RR, drx, rx;
double Gc2=G/c/c, Gc4=G/c/c/c/c;
double p0, p, e, e0, m, mR, dm, A, B=1.0, dB, BR, Bfactor, m1, m2, m3, m4, p1, p2, p3, p4, B1, B2, B3, B4, power;
double H1, H0, K, W, X, F, V, Dv, gamma, Dp, N;
double DH11,DK1,DW1,DX1,H01,H02,K1,K2,x,X1,X2,Xp1,Xp2,W1,W2,V1,V2,V01,V02,V0;
double w,wcheck, o[2],wi;
double aR,bR,gR,hR,kR,n,Y1,Y2,Z,DZ,DDZ,VZ,Ar1,Ar2,Ai1,Ai2,ar,ai,Ar[2],Ai[2],Br[2],Bi[2];
int t, q;
double times;

printf("Please enter the data file's name: ");
scanf("%s",filename);
if((inf=fopen(filename,"r"))==NULL)
  {printf("file not found!\n");exit(0);}
while(fscanf(inf,"%lf",&a[j])==1){fscanf(inf,"%lf",&b[j]);j++;}
printf("Please enter a output file's name: ");
scanf("%s",opfile);
outf=fopen(opfile,"w");
power=pow(10.0,b[0]);

for(times=3.0;times>=3.0;times=times-0.5)
{
  mfile=fopen("minside.txt","w");
  rhofile=fopen("einside.txt","w");
  pfile=fopen("pinside.txt","w");
  Bfile=fopen("Binside.txt","w");
  r=r0;
  printf("Please enter the central density (1.0E18 for example): ");
  scanf("%lf",&e0);
  printf("Please enter the trial range of the part of complex frequency (3.0E-5 for example):\n");
  printf("omega1= "); scanf("%lf",&o[0]); printf("omega2= "); scanf("%lf",&o[1]);
  printf("No erro occurs. Processing...\n");
  e=e0;
  e0=e0*Gc2;
  p=che(e);
  p0=p*Gc4;
  m=4*r*r*pi*e*r/3;
  for(;p>power; r=r+dr)
    {
    fprintf(rhofile,"%.15e",e*Gc2);   /*-- e,p,m in G=c=1 --*/
    fprintf(rhofile,"\n");
    fprintf(pfile,"%.15e",p*Gc4);
    fprintf(pfile,"\n");
    fprintf(Bfile,"%f",B);
    fprintf(Bfile,"\n");
    fprintf(mfile,"%.15e",m*Gc2);
    fprintf(mfile,"\n");
    p1=f(r,p,e,m);
    m1=fm(r,e);
    if((p+dr*p1/2)>power) e=ch(p+dr*p1/2); else break;
    B1=Bf(r,p,m,B);
    p2=f(r+dr/2,p+dr*p1/2,e,m+dr*m1/2);
    m2=fm(r+dr/2,e);
    if((p+dr*p2/2)>power) e=ch(p+dr*p2/2); else break;
    B2=Bf(r+dr/2,p+dr*p1/2,m+dr*m1/2,B+dr*B1/2);
    p3=f(r+dr/2,p+dr*p2/2,e,m+dr*m2/2);
    m3=fm(r+dr/2,e);
    if((p+dr*p3)>power) e=ch(p+dr*p3); else break;
    B3=Bf(r+dr/2,p+dr*p2/2,m+dr*m2/2,B+dr*B2/2);
    p4=f(r+dr,p+dr*p3,e,m+dr*m3);
    m4=fm(r+dr,e);
    B4=Bf(r+dr,p+dr*p3,m+dr*m3,B+dr*B3);
    p=p+dr*(p1+2*p2+2*p3+p4)/6;
    e=ch(p);
    dm=dr*(m1+2*m2+2*m3+m4)/6;
    m=m+dm;
    B=B+dr*(B1+2*B2+2*B3+B4)/6;
    }
  R=r;
  mR=m*Gc2;
  BR=1-2*G*m/r/c/c;
  Bfactor=BR/B;
  gamma=(b[1]-b[0])/(a[1]-a[0]);
  N=1/(gamma-1);
  RR=R-(N+1)*(p-dr*(p1+2*p2+2*p3+p4)/6)/(p1+2*p2+2*p3+p4)*6;
  fclose(rhofile);
  fclose(pfile);
  fclose(Bfile);
  fclose(mfile);
  Bfile=fopen("Binside.txt","r");
  Bcor=fopen("Bcor.txt","w");
  while(fscanf(Bfile,"%lf",&copy)==1){fprintf(Bcor,"%f",copy*Bfactor);fprintf(Bcor,"\n");}
  fclose(Bfile);
  fclose(Bcor);

q=1;
wcheck=0;
for(t=0;;t++)
{
if(t==0)
w=o[t];
else
w=o[q];

/*-- inner solutions --*/

   /*-- Solution 1--*/
Wfile=fopen("W1.txt","w");
Vfile=fopen("V1.txt","w");
mfile=fopen("minside.txt","r");
rhofile=fopen("einside.txt","r");
pfile=fopen("pinside.txt","r");
fscanf(pfile,"%lf",&p);
fscanf(rhofile,"%lf",&e);
Bcor=fopen("Bcor.txt","r");
fscanf(Bcor,"%lf",&B);
fclose(Bcor);
fclose(rhofile);
fclose(pfile);

W=1.0; K=(e+p);
X=(e+p)*sqrt(B)*((4*pi/3*(e+3*p)-w*w/B/l)*W+0.5*K); 
H1=(2*l*K+16*pi*(e+p)*W)/l/(l+1);

rhofile=fopen("einside.txt","r");
pfile=fopen("pinside.txt","r");
Bcor=fopen("Bcor.txt","r");


r=r0; 
while(r<=R)
{
fscanf(pfile,"%lf",&p);
fscanf(rhofile,"%lf",&e);
fscanf(Bcor,"%lf",&B);
fscanf(mfile,"%lf",&m);
if((wcheck-w)/w>-Cri&&(wcheck-w)/w<Cri)
{
fprintf(Wfile,"%.15e",sqrt(1-2*m/r)*W);
fprintf(Wfile,"\n");
}

Dp=-(e+p)*(m+4*pi*r*r*r*p)/r/r/(1-2*m/r);
Dv=-2*Dp/(e+p);


A=1/(1-2*m/r);
gamma=gammaf(e);
H0=H0f(r,B,X,m,p,A,H1,K,w);
V=Vf(r,w,e,p,B,A,Dp,W,H0,X);
if(r==r0)V01=V;else;
if((wcheck-w)/w>-Cri&&(wcheck-w)/w<Cri)
{
fprintf(Vfile,"%.15e",V);
fprintf(Vfile,"\n");
}
F=Ff(r,A,e,p,m,Dp);

DH11=DH1(r,m,A,p,e,H1,H0,K,V);
DK1=DK(r,H0,H1,Dv,K,e,p,A,W);
DW1=DW(r,W,A,gamma,p,B,X,V,H0,K);
DX1=DX(r,X,e,p,B,Dv,H0,w,H1,K,V,A,F,W);
H1=H1+DH11*dr;
K=K+DK1*dr;
W=W+DW1*dr;
X=X+DX1*dr;
r=r+dr;
}
fclose(Wfile);
fclose(Vfile);
fclose(Bcor);
fclose(rhofile);
fclose(pfile);
fclose(mfile);
X1=X; Xp1=DX1; K1=K; H01=H0f(r,B,X,m,p,A,H1,K,w); W1=W; V1=V;

/*-- 2nd solution --*/
Wfile=fopen("W2.txt","w");
Vfile=fopen("V2.txt","w");
mfile=fopen("minside.txt","r");
rhofile=fopen("einside.txt","r");
pfile=fopen("pinside.txt","r");
fscanf(pfile,"%lf",&p);
fscanf(rhofile,"%lf",&e);
Bcor=fopen("Bcor.txt","r");
fscanf(Bcor,"%lf",&B);
fclose(Bcor);
fclose(rhofile);
fclose(pfile);

W=1.0; K=-(e+p);
X=(e+p)*sqrt(B)*((4*pi/3*(e+3*p)-w*w/B/l)*W+0.5*K); 
H1=(2*l*K+16*pi*(e+p)*W)/l/(l+1);

rhofile=fopen("einside.txt","r");
pfile=fopen("pinside.txt","r");
Bcor=fopen("Bcor.txt","r");


r=r0; 
while(r<=R)
{
fscanf(pfile,"%lf",&p);
fscanf(rhofile,"%lf",&e);
fscanf(Bcor,"%lf",&B);
fscanf(mfile,"%lf",&m);
if((wcheck-w)/w>-Cri&&(wcheck-w)/w<Cri)
{
fprintf(Wfile,"%.15e",sqrt(1-2*m/r)*W);
fprintf(Wfile,"\n");
}

Dp=-(e+p)*(m+4*pi*r*r*r*p)/r/r/(1-2*m/r);
Dv=-2*Dp/(e+p);

A=1/(1-2*m/r);
gamma=gammaf(e);
H0=H0f(r,B,X,m,p,A,H1,K,w);
V=Vf(r,w,e,p,B,A,Dp,W,H0,X);
if(r==r0)V02=V;else;
if((wcheck-w)/w>-Cri&&(wcheck-w)/w<Cri)
{
fprintf(Vfile,"%.15e",V);
fprintf(Vfile,"\n");
}
if(r==r0)V02=V;else;
F=Ff(r,A,e,p,m,Dp);

DH11=DH1(r,m,A,p,e,H1,H0,K,V);
DK1=DK(r,H0,H1,Dv,K,e,p,A,W);
DW1=DW(r,W,A,gamma,p,B,X,V,H0,K);
DX1=DX(r,X,e,p,B,Dv,H0,w,H1,K,V,A,F,W);
H1=H1+DH11*dr;
K=K+DK1*dr;
W=W+DW1*dr;
X=X+DX1*dr;
r=r+dr;
}
fclose(Wfile);
fclose(Vfile);
fclose(Bcor);
fclose(rhofile);
fclose(pfile);
fclose(mfile);
X2=X; Xp2=DX1;K2=K; H02=H0f(r,B,X,m,p,A,H1,K,w); W2=W; V2=V;

/*-- determing the unique solution --*/
x=-(X1-(RR-R)/(N+1)*Xp1)/(X2-(RR-R)/(N+1)*Xp2);
H0=H01+x*H02; K=K1+x*K2; W=W1+x*W2; V=V1+x*V2;V0=V01+x*V02;

/*-- Solution Review, May be skipped  --*/
if((wcheck-w)/w>-Cri&&(wcheck-w)/w<Cri)
{
Wfile1=fopen("W1.txt","r");
Wfile2=fopen("W2.txt","r");
Wfile=fopen("W.txt","w");
Vfile1=fopen("V1.txt","r");
Vfile2=fopen("V2.txt","r");
Vfile=fopen("V.txt","w");
r=r0;
while(r<=R)
{
fscanf(Wfile1,"%lf",&W1);
fscanf(Wfile2,"%lf",&W2);
fscanf(Wfile2,"\n");
W=W1+x*W2;
fprintf(Wfile,"%e",W/(1+x));
fprintf(Wfile,"\n");
fscanf(Vfile1,"%lf",&V1);
fscanf(Vfile2,"%lf",&V2);
V=V1+x*V2;
fprintf(Vfile,"%e",V/V0);
fprintf(Vfile,"\n");
r=r+dr;
}
fclose(Wfile1);
fclose(Wfile2);
fclose(Wfile);
fclose(Vfile1);
fclose(Vfile2);
fclose(Vfile);
}
if((wcheck-w)/w>-Cri&&(wcheck-w)/w<Cri)break;
wcheck=w;
/*-- BC at the surface --*/
n=0.5*(l-1)*(l+2);
aR=-(n*R+3*mR)/(w*w*R*R-(n+1)*mR/R);
bR=(n*R*(R-2*mR)-w*w*R*R*R*R+mR*(R-3*mR));
bR=bR/(R-2*mR)/(w*w*R*R-(n+1)*mR/R);
gR=n*(n+1)*R*R+3*n*mR*R+6*mR*mR;
gR=gR/R/R/(n*R+3*mR);
hR=-n*R*R+3*n*mR*R+3*mR*mR;
hR=hR/(R-2*mR)/(n*R+3*mR);
kR=-R*R/(R-2*mR);
Y1=K;
Y2=aR*H0+bR*K;
Z=(kR*Y1-Y2)/(kR*gR-hR);
DZ=(gR*Y2-hR*Y1)/(gR*kR-hR);

  /*-- check=gR*Y2-hR*Y1; --*/

/*-- Integrating the wave equation outward --*/

for(r=R;r<25.0/w;r=r+dr)
{
drx=dr/(1-2*mR/r);
VZ=(1-2*mR/r)/r/r/r/(n*r+3*mR)/(n*r+3*mR);
VZ=VZ*(2*n*n*(n+1)*r*r*r+6*n*n*mR*r*r+18*n*mR*mR*r+18*mR*mR*mR);
DDZ=(VZ-w*w)*Z;
Z=Z+DZ*drx;
DZ=DZ+DDZ*drx;
}
r=r-dr;
rx=r+2*mR*log(r/2/mR-1);

/*-- Matching the solution at infinity --*/

Ar1=2*cos(w*rx)-2*(n+1)/w/r*sin(w*rx)+1/w/w/r/r*(1.5*mR*w*(1+2/n)*sin(w*rx)-n*(n+1)*cos(w*rx));
Ai1=2*sin(w*rx)+2*(n+1)/w/r*cos(w*rx)-1/w/w/r/r*(1.5*mR*w*(1+2/n)*cos(w*rx)+n*(n+1)*sin(w*rx));
Ar2=-2*w*sin(w*rx)-2*(n+1)*cos(w*rx)/r+1/w/r/r*(1.5*mR*w*(1+2/n)*cos(w*rx)+n*(n+1)*sin(w*rx))+(1-2*mR/r)*2*(n+1)/w/r/r*sin(w*rx);
 Ai2=2*w*cos(w*rx)-2*(n+1)*sin(w*rx)/r+1/w/r/r*(1.5*mR*w*(1+2/n)*sin(w*rx)-n*(n+1)*cos(w*rx))-(1-2*mR/r)*2*(n+1)/w/r/r*cos(w*rx);

ar=(Ai2*Z-Ai1*DZ)/(Ar1*Ai2-Ar2*Ai1);
ai=-(Ar1*DZ-Ar2*Z)/(Ar1*Ai2-Ar2*Ai1); /*-- a minus sign in front because we need the complex conjugate --*/
if(t==0)
{
Ar[t]=ar;
Ai[t]=ai;
}
else
{
Ar[q]=ar;
Ai[q]=ai;
Br[0]=(o[0]*Ar[1]-o[1]*Ar[0])/(o[0]-o[1]);
Br[1]=(Ar[0]-Ar[1])/(o[0]-o[1]);
Bi[0]=(o[0]*Ai[1]-o[1]*Ai[0])/(o[0]-o[1]);
Bi[1]=(Ai[0]-Ai[1])/(o[0]-o[1]);
w=-(Br[0]*Br[1]+Bi[0]*Bi[1])/(Br[1]*Br[1]+Bi[1]*Bi[1]);
if (w<=o[0]){o[1]=o[0];o[0]=w;Ar[1]=Ar[0];Ai[1]=Ai[0];q=0;}
else if(w>=o[1]){o[0]=o[1];o[1]=w;Ar[0]=Ar[1];Ai[0]=Ai[1];q=1;}
else if((o[1]-w)>(w-o[0])){o[1]=w;q=1;}
else {o[0]=w;q=0;}
}
}
printf("%d ",t);
wi=(Br[0]*Bi[1]-Bi[0]*Br[1])/(Br[1]*Br[1]+Bi[1]*Bi[1]);
fprintf(outf,"%f    %f     %f    %f\n", mR/G*c*c/M,RR/1000,(w*3E5/2/pi),1/wi/c);
/*-- fprintf(outf,"%f       %e\n", sqrt(mR/R/R/R*1E6),w/2/pi*3E5); --*/  /*-- f in msec-1 --*/
/*-- fprintf(outf,"%f       %e\n", mR/G*c*c/M,w/2/pi*3E5); --*/
}
fclose(outf);
  printf("Results have been saved. \nPress any key to contiune\n",opfile);
  getch();
}
