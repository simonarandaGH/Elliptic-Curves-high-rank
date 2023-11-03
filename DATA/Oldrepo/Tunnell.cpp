

// (C)2023 Simon Aranda.

/******* 
This small C++ program calculates in a few seconds all congruent numbers below 4*10^6. It is nothing more than the Tunnell's T. criterion, applied in a broad way. It is trivial to extend it to 10^9, etc. As well as modifying it to print the A_[] vectors.
*******/

  
typedef std::int_fast64_t I64;
typedef std::int_fast32_t I32;

#define CLIM 2000
#define ASIZE (CLIM*CLIM)
static array<I32,ASIZE> A1,A2,A3,A4;

// n= core*(mxq^2)
I32 I32maxsquare(const I32 p)
{
double a= p;
double b= sqrt(a);
I32 c=floor(b)+1;
for(I32 s= c;s>0;s--) {
  if(p%s>0)continue;
	I32 ss= (s*s);
  if((p%ss)==0)return s;
	}
return 1;
}

// core()
I32 I32core(I32 p)
{
I32 s= I32maxsquare(p);
I32 ss= (s*s);
I32 cr= (p/ss);
return cr;
}

int main(int argc,char* argv[])
{
A1.fill(0);A2.fill(0);A3.fill(0);A4.fill(0);
for(I32 x= 0;x< CLIM;x++)  {
  for(I32 y= 0;y< CLIM;y++)  {
    for(I32 z= 0;z< CLIM;z++)  {
      int pval=1;       // pound value
      if(x!=0)pval*=2;
      if(y!=0)pval*=2;
      if(z!=0)pval*=2;
      I64 n1= +2*x*x +1*y*y + 8*z*z; // n odd a
      I64 n2= +2*x*x +1*y*y +32*z*z; // n odd b
      if(n1<ASIZE) A1[n1]+= pval;
      if(n2<ASIZE) A2[n2]+= pval;
      //
      I64 n3= +4*x*x +1*y*y + 8*z*z; // n even a
      I64 n4= +4*x*x +1*y*y +32*z*z; // n even b
      if(n3<ASIZE) A3[n3]+= pval;
      if(n4<ASIZE) A4[n4]+= pval;
      }
    }
  }

#define KMOD 8
array<I32,KMOD> VRM;
I32 r0,r1,r2;r0=r1=r2=0;
VRM.fill(0);

for(I32 br= 1;br< ASIZE;br++)  // brute number
  {
  I32 n= I32core(br);
  bool even= (n%2==0); // n = core number.
  I32 a,b,i;
  if(even) {
    i= n/2;a =A3[i];b =A4[i]; // even core
    }
  else {
    i=n;a =A1[i];b =A2[i]; // Odd core
    }
  int r=0;        // r() function
  if(a==0 && b==0)r=1;
  if(a>0 && b>0 && a==2*b)r=2;

  if(r==0)r0++;
  if(r==1)r1++;
  if(r==2)r2++;
  if(r>0) {
    int m= br%8;VRM[m]++;
    }

  if(false) {
    cout<<br<<", "<<n<<", "<<r;
    cout<<";  ("<<i<<": "<<a<<"/"<<b<<") ";
    if(r>0)cout<<" CN.";
    cout<<endl;
    }
  } // for

cout<<"r-distribution: 0,1,2,02,12 :"<<endl;
cout<<r0<<endl;
cout<<r1<<endl;
cout<<r2<<endl;
cout<<(r0+r2)<<endl;
cout<<(r1+r2)<<endl;

cout<<"VRM:"<<endl;
I32 st=0;
for(unsigned i=0;i<VRM.size();i++) {
  cout<<i<<"%8: "<<VRM[i]<<endl;
  st+=VRM[i];
  }
cout<<"t= "<<st<<endl;

return 0;
}

