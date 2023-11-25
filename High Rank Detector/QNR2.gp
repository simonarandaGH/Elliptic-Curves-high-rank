\\\\\\\\\
\\ QNR2.gp
\\ QNR version 2.0
\\ (C)Nov2023.Simon Aranda
\\ PARIGP script
\\ 2-phase closed loop high rank detector
\\

\\ for any fam [-,-].
\\ q rational points counter
\\ n number
\\ r rank of the curve

pqnr(q,n,r)=
{
print1(" (" q "/" n "/" r ") ");
}

{
\\ init vars
default(logfile,"temp_log000.txt");
qex=0; \\ same rank cntr
nx=0;  \\ n with max rank
rx=0;  \\ max rank
gqmin=oo;  \\ global cntr minimal

\\ tunning
BSZ=10^4; \\ block size bytes

for(bnum= 0 ,oo, \\ work in blocks
  bbas=bnum*BSZ; \\ block offset

  print;
  printsep1(" ","(-2,-1,0,-1,n^2)  #blk=",bnum,bbas);
  pqnr(gqmin,nx,rx);print(" ." qex);

  Q1=vector(BSZ);   \\ first Q
  E=vector(BSZ);    \\ from init

  \\ adjustable first h
  naiveh= 80*(bnum+1);  \\ common h
  for(m= 1,BSZ,         \\ the points
    n= bbas+m;

    \\ the curve coeff form
    coef=[-2,-1,0,-1,n^2];

    te=ellinit(coef,D=1);if(#te<1,next;);
    E[m]=te;
    tl=ellratpoints(te,naiveh); \\ not stored
    Q1[m]=#tl;
    );
  \\ first min max
  pin1=vecmin(Q1);pax1=vecmax(Q1);

  \\ tunning first cut
  cut1= (pax1+pin1)>>1;
  printsep(" ","Cut1: ",pin1,"/",pax1,"-->",cut1);

  \\ init second phase
  Q =vector(BSZ,i,0);P=vector(BSZ,i,[]);
  qrph=0; \\ call cntr 
  qmin=oo; \\ qmin because Q sparse

  \\ tunning second h
  naiveh*= 10; \\ common
  for(m= 1,BSZ,    \\ the points
    if(Q1[m]<cut1,next;);
    n= bbas+m;qrph++;

    \\ optional h*Q1 or only h.
    lp=ellratpoints(E[m],naiveh*Q1[m]);
    tt=#lp;Q[m]=tt;P[m]=lp;
    if(tt<qmin,qmin=tt;); \\ qmin because..
    );
  print("@ellratpoints: " qrph);

  \\ second min max
  qmax=vecmax(Q);
  if(gqmin==oo,gqmin=qmax;);

  \\ Design options
  \\ low old gqmin or not. 
  \\ gqmin=gqmin>>1;if(gqmin<1,gqmin=1;);
  
  \\ skip block or not
  \\ skip= qmax < gqmin;
  \\ skip=0;
  \\ if(skip,    printsep(" ","blk skip: ",qmax,"<",gqmin);
  \\ next;        );

  \\ tunning second cut
  if(qmax<gqmin,
    qcut2=qmax-(qmin>>1);
    gqmin=gqmin -qmin;if(gqmin<1,gqmin=1;);
    ,
    qcut2= (qmax+qmin)>>1;
    );

  printsep(" ","Cut2",qmin,"/",qmax,"-->",qcut2);

  qrnk=0; \\ call cntr
  for(m=1,BSZ,
    if(Q[m]<qcut2,next;);
    eri=ellrankinit(E[m]);
    [rn,tb,tc,tg]=ellrank(eri,1,P[m]);qrnk++;

    \\ process rank result
    n=bbas+m;
    if(rn<rx,
      \\ Option aprox Q or not
      gqmin=max(gqmin,Q[m]+(rx-rn));
      
      print1("butnot: ");pqnr(Q[m],n,rn);print;
      next;
      );
    if(rn==rx,
      nx=n;qex++;gqmin=Q[m];
      print1("EQUAL: ");pqnr(Q[m],n,rn);
      print("....................");
      next;
      );
    if(rn>rx,
      rx=rn;nx=n;gqmin=Q[m];qex=1;
      print1("HIGH:  ");pqnr(Q[m],n,rn);
      print("....................");
      next;
      );
    ); \\ end for

  print("@ellrank: " qrnk );
  ); \\ end for all blk

}
