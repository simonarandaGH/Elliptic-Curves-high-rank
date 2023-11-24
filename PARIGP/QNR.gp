\\\\\\\\\
\\ QNR.gp
\\ (C)2023.Simon Aranda
\\ for any fam [-,-].
\\ q rat points counter
\\ n number
\\ r rank of the curve.



pqnr(q,n,r)=
{
print1(" (" q "/" n "/" r ") ");
}

{
default(logfile,"temp_qnr.txt");

qex=0;nx=0;rx=0;BSZ=10^4; \\ block size bytes
gqmin=oo;

for(bnum= 0 ,oo,
  bbas=bnum*BSZ; \\ blk n offset
  print;
  printsep(" ","(*,*)  #blk=",bnum,bbas," gqmin="
    ,gqmin," qex=",qex," nx=",nx," rx=",rx);

  Q1=vector(BSZ);E=vector(BSZ);
  h= 200*(bnum+1); \\ common h
  for(m= 1,BSZ,    \\ the points
    n= bbas+m;

    coef=[0,-n]; \\ any

    te=ellinit(coef,D=1);if(#te<1,next;);
    tl=ellratpoints(te,h);
    Q1[m]=#tl;E[m]=te;
    );
  pin1=vecmin(Q1);pax1=vecmax(Q1);
  cut1= (pax1+pin1)>>1;
  printsep(" ","Cut1: ",pin1,"/",pax1,"-->",cut1);

  Q =vector(BSZ,i,0);P=vector(BSZ,i,[]);
  qrph=0;
  h= BSZ*(bnum+1); \\ common
  for(m= 1,BSZ,    \\ the points
    if(Q1[m]<cut1,next;);
    n= bbas+m;
    lp=ellratpoints(E[m],h);qrph++;
    Q[m]=#lp;P[m]=lp;
    );
  print("@ellratpoints: " qrph);
  qmax=vecmax(Q);

  if(gqmin==oo,gqmin=qmax;);
  if(qmax<gqmin,
    printsep(" ","blk skip: ",qmax,"<",gqmin);
    gqmin=gqmin-qmax;if(gqmin<1,gqmin=1;);
    next;
    );

  \\ second cut
  qcut2=qmax>>1;
  printsep(" ","Cut2",gqmin,"/",qmax,"-->",qcut2);

  qrnk=0;
  for(m=1,BSZ,
    if(Q[m]<qcut2,next;);
    eri=ellrankinit(E[m]);
    [rn,tb,tc,tg]=ellrank(eri,1,P[m]);qrnk++;
    n=bbas+m;
    if(rn<rx,
      gqmin=max(gqmin,Q[m]+(rx-rn));
      print1("butnot: ");pqnr(Q[m],n,rn);print;
      next;
      );
    if(rn==rx,
      nx=n;
      print1("EQUAL: ");pqnr(Q[m],n,rn);
      print("....................");
      gqmin=Q[m];qex++;
      next;
      );
    if(rn>rx,
      rx=rn;nx=n;
      print1("HIGH:  ");pqnr(Q[m],n,rn);
      print("....................");
      gqmin=Q[m];qex=1;
      next;
      );
    ); \\ end for

  print("@ellrank: " qrnk );
  ); \\ end for all blk

}
