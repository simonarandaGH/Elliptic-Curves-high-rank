
\\ nn1.gp
\\ (C)2023.Simon Aranda

{
default(logfile,"temp_nn1seq.txt");

rx=0;qrpmin=0;
BLK=10000;print("BlkSz=" BLK);
\\ start number
for(bnum=  75  ,oo,
  bbeg=bnum*BLK;bend=(bnum+1)*BLK-1;
  print;
  printsep(" ","new blk=",bnum
  ,"  n:",bbeg,"..",bend);

  qover=0;
  for(tim=1,oo,
    if(tim>1,
      qrpmin=floor(0.7*qrpmin);
      );
    print;
    printsep(" ","srch blk:",bnum,".",tim);

    for(n= bbeg,bend,
      if(n%1000==0,
        printsep(" ","tick  n=",n
        ,"min=",qrpmin,"rx=",rx;);
        );
      
      coef=[-n^2,+1];
      EE=ellinit(coef);if(#EE<1,next;);
      h=100+ ceil(n/4);
      rpl=ellratpoints(EE,h);
      qrp=#rpl;
      if(qrp< qrpmin,next;);

      ERI=ellrankinit(EE);
      [rn,tb,tc,gpl]=ellrank(ERI,,rpl);
      qover++;
  
      if(rn<rx,
        qrpmin=qrp+1;
        printsep(" "
        ,"butnot:",qrp,n,rn
        ,"min->",qrpmin);
        next;
        );
      if(rn==rx,
        qrpmin=qrp;
        print;printsep(" "
        ,qrp,";",n,rn
        ," <---EQUAL: min->",qrpmin
        ," <---EQUALRANK");
        print;
        next;
        );
      rx=rn;
      qrpmin=qrp;
      print;
      printsep(
        " ",qrp,";",n,rn
        ," <---HIGH: min->",qrpmin
        ," <---HIGHRANK");
      print;
      );
    printsep(" ","doend: qover=",qover);
    if(qover>0,break;);
    );
  printsep(" ","block done. rx=",rx);
  );
}

