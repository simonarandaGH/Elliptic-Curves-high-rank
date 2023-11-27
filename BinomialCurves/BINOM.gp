\\ binom.gp
\\(C)2023.SimonAranda


{
default(logfile,"temp_log.txt");

for(k=1,oo,
  n=k^2;

  print1(n " : ");
  for(m=0,n,
    coef= [0 ,0 ,0, -binomial(n,m)^2 , +n^3];
    ei=ellinit(coef,1);
    mm=ellminimalmodel(ei);
    eri=ellrankinit(mm);
    r=ellrank(eri);rn=r[1];
    print1(" " rn);
    );
  print;
  );

}

