\\\\\\\\\
\\ QNR.gp
\\ (C)2023.Simon Aranda
\\ fam [-n^2,+1].
\\ q number of points counter
\\ n number
\\ r rank of the curve C[n].

pqnr(q,n,r)=
{
print1(" (" q "/" n "/" r ") ");
}

{
default(logfile,"temp_qnrnn1.txt");

rx=0;BSZ=10^4;gqmin=oo;
for(bnum= 0,oo,
  bbas=bnum*BSZ; \\ blk n offset
  print;
  printsep(" ","#blk=",bnum,bbas," gqmin="
    ,gqmin," rx=",rx);

  \\ cache mem
  Q1=vector(BSZ);E=vector(BSZ);
  Q =vector(BSZ,i,0);P=vector(BSZ,i,[]);

  h= 200*(bnum+1); \\ common h
  for(m= 1,BSZ,               \\ the points
    n= bbas+m;coef=[-n^2,+1];
    te=ellinit(coef,D=1);if(#te<1,next;);
    tl=ellratpoints(te,h);
    Q1[m]=#tl;E[m]=te;
    );
  pin1=vecmin(Q1);pax1=vecmax(Q1);
  tt1= (pax1-pin1)>>1;cut1= pin1+tt1;
  printsep(" ","Cut1: ",pin1,pax1,cut1);

  qrph=0;
  h= (bnum+1)*BSZ; \\ common
  for(m= 1,BSZ,               \\ the points
    if(Q1[m]<cut1,next;);
    n= bbas+m;
    lp=ellratpoints(E[m],h);qrph++;
    Q[m]=#lp;P[m]=lp;
    );
  print("ellratpoints calls= " qrph);
  qmax=vecmax(Q);
  if(gqmin==oo,gqmin=qmax;);
  if(qmax<gqmin,gqmin--;); \\ time decr
  if((qmax+4)<gqmin,
    printsep(" ","blck skip: ",qmax,gqmin);
    next;
    );

  \\ second cut
  if(qmax<gqmin,qcut=qmax-1;,qcut=gqmin-1;);
  printsep(" ","gqmin=",gqmin,"qmax=",qmax,"cut=",qcut);

  qrnk=0;
  for(m=1,BSZ,
    if(Q[m]<qcut,next;);
    eri=ellrankinit(E[m]);
    [rn,tb,tc,tg]=ellrank(eri,1,P[m]);qrnk++;

    n=bbas+m;
    if(rn<rx,
      gqmin=max(gqmin,Q[m]);
      print1("butnot: ");pqnr(Q[m],n,rn);print;
      next;
      );
    if(rn==rx,
      print1("EQUAL: ");pqnr(Q[m],n,rn);print;
      gqmin=max(gqmin,Q[m]);
      next;
      );
    if(rn>rx,
      rx=rn;
      print1("HIGH:  ");pqnr(Q[m],n,rn);print;
      gqmin=Q[m];
      next;
      );
    ); \\ end 2 for

  print("ellrank calls=" qrnk );
  ); \\ end for all blk

}
