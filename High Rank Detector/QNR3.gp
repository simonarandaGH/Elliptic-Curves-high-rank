\\\\\\\\\
\\ QNR3.gp
\\ QNR version 3.0
\\ (C)Nov2023.Simon Aranda
\\ PARIGP script
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
default(logfile,"temp_log777.txt");


qex=0; \\ same rank cntr
nx=0;  \\ n with max rank
rx=0;  \\ max rank
qx=0;

BSZ=10^3; 
for(bnum= 0,oo,
  bbas=bnum*BSZ; \\ block offset
  bcen=bbas+(BSZ>>1);

  print;printsep1(" "
  ,"(xxx)  #blk=",bnum,bbas+1);
  pqnr(qx,nx,rx);
  printsep(" ",".",qex);

  E=vector(BSZ,i,0);    \\ from init
  for(m= 1,BSZ,     \\ the points
    n= bbas+m;
    \\ the curve coeff form
    \\ coef= [ 0 , n+1 , 0 , n ,0];
    \\ coef= [ 0 , -10*n+1 , 0 , -10*n ,0];
    \\ coef= [0,4*n^2+12*n-3,0,32*(n+3),0];

    coef= [0, 4*n^2 -12*n-3 ,0, 32*(-n+3)  ,0];
    te=ellinit(coef,1);
    if(#te<1,
      print("error ellinit " m);E[m]=0;
      ,
      \\ mm=ellminimalmodel(te);
      E[m]=te;
      );
    );
  print("inits done");
  \\ examine block
  Q=vector(BSZ,i,0);P=vector(BSZ,i,[]);  \\ 
  h= bcen;
  for(m= 1,BSZ,
    if(E[m]==0,next;);
    tl= ellratpoints(E[m],h,0);
    tt=#tl;Q[m]=tt;P[m]=tl;
    );
  \\ the best
  qqxx=vecmax(Q); qqll=vecmin(Q);
  if(qqxx==qqll,
    lim=oo;
    ,
    half=(qqxx+qqll)>>1;lim=max(qx,half);
    );
  printsep(" ",qqll,qqxx,"lim=",lim);
  qrnk=0;
  for(xm= 1,BSZ,
    xq=Q[xm];if(xq< lim,next;);
    xn=bbas+xm;
    xri=ellrankinit(E[xm]);
    [xr,tb,tc,tg]=ellrank(xri,,P[xm]);
    qrnk++;
    if(xr==rx,
      nx=xn;qx=xq>>1;qex++;
      pqnr(xq,xn,xr);
      print(" <-- EQUAL.................");
      );
    if(xr>rx,
      rx=xr;nx=xn;qx=xq>>1;qex=1;
      pqnr(xq,xn,xr);
      print(" <-- HIGH..................");
      );
    );
  print(qrnk);
  );
}
