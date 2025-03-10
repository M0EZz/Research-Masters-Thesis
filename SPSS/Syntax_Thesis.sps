* Encoding: UTF-8.

CROSSTABS
  /TABLES=UpdatedAIAct BY UpdatedRisk
  /FORMAT=AVALUE TABLES
  /STATISTICS=CHISQ GAMMA D CTAU 
  /CELLS=COUNT TOTAL 
  /COUNT ROUND CELL.

LOGISTIC REGRESSION VARIABLES Q5RiskPresence
  /METHOD=ENTER AIAct_Continuous Year Topic_1 Topic_2 Topic_3 Topic_4 Topic_5 Topic_6 Topic_7 
    Topic_8 Topic_9 Outlet_Broadsheet Outlet_Tabloid Topic_0 
  /CONTRAST (Topic_1)=Indicator(1)
  /CONTRAST (Topic_2)=Indicator(1)
  /CONTRAST (Topic_3)=Indicator(1)
  /CONTRAST (Topic_4)=Indicator(1)
  /CONTRAST (Topic_5)=Indicator(1)
  /CONTRAST (Topic_6)=Indicator(1)
  /CONTRAST (Topic_7)=Indicator(1)
  /CONTRAST (Topic_8)=Indicator(1)
  /CONTRAST (Topic_9)=Indicator(1)
  /CONTRAST (Outlet_Broadsheet)=Indicator(1)
  /CONTRAST (Outlet_Tabloid)=Indicator(1)
  /SAVE=PRED PGROUP
  /PRINT=GOODFIT CI(95)
  /CRITERIA=PIN(0.05) POUT(0.10) ITERATE(20) CUT(0.5).

DATASET ACTIVATE DataSet3.

/* This macro computes Krippendorff's alpha reliability estimate for judgments */.
/* made at any level of measurement, any number of judges, with or */.
/* without missing data.  The macro assumes the data file is set up */.
/* in a SPSS data file with judges as the variables and the units being */.
/* judged in the rows.  The entries in the data matrix should be */.
/* the coding (quantified or numerically coded for nominal judgments) given */.
/* to the unit in that row by the judge in that column.  Once the macro is */.
/* activated (by running the command set below), the syntax is */.
/* */.
/* KALPHA judges = judgelist/level = a/detail = b/boot = z.\n/* */.
/* where 'judgelist' is a list of variable names holding the names of the */.
/* judges, 'a' is the level of measurement (1 = nominal, 2 = ordinal, */.
/* 3 = interval, 4 = ratio), 'b' is set to 1 if you desire SPSS to print */.
/* the coincidence and delta matrices, and 'z' is the number of bootstrap */.
/* samples desired for inference;  z must be at least 1000 and is truncated to the */.
/* lowest 1000 entered (for example, 2300 is truncated to 2000) */.
/* The '/level' and '/detail' and '/boot' subcommands are */.
/* optional and default to 1,0, and 0, respectively, if omitted */.
/* */.
/* Missing data should be represented with a 'period' character */.
/* Units that are not coded by at least one judge are excluded from */.
/* the analysis */.
/* */.
/* This macro is version 4.0, updated with a new bootstrapping algorithm on Dec 18, 2018  */.
/* */.
/* */.
/* Written by Andrew F. Hayes */.
/* http://www.afhayes.com */.


DEFINE kalpha (judges = !charend ('/')/level = !charend('/') !default(1)/detail = !charend('/') 
    !default(0)/boot = !charend('/') !default(0)/seed=!charend('/') !default(random)).
PRESERVE.
SET MXLOOP = 900000000.
SET LENGTH = NONE.
SET SEED = !seed.
SET PRINTBACK = OFF.
MATRIX.
get dat/variables = !judges/file = */names = vn/missing = -9999999.
compute btn = !boot.
do if (!boot > 0).
  compute btn = trunc(!boot/1000)*1000.
end if.
do if (!boot > 0 and btn = 0).
  print/title = "Number of bootstraps must be at least 1000.".
end if.
compute btprob = 0.

/* FIRST WE CREATE THE DATA FILE EXCLUDING OBJECTS WITH ONLY ONE JUDGMENT */.
/* THAT DATA FILE IS HELD IN DAT AND DAT3 */.

compute rw = 1.
loop i = 1 to nrow(dat).
  compute good = 0.
  loop j = 1 to ncol(dat).
    do if (dat(i,j) <> -9999999).
      compute good = good + 1.
    end if.
  end loop.
  do if (good > 1).
    compute dat(rw,:) = dat(i,:).
    compute rw = rw+1.
  end if.
end loop.
compute dat = dat(1:(rw-1),:).
compute nj = ncol(dat).
compute nobj = nrow(dat).
compute dat3 = dat.

/* NOW WE CREATE A SINGLE COLUMN OF DATA TO FIGURE OUT HOW MANY */.
/* UNIQUE JUDGMENTS ARE MADE, AND WE SORT IT */.

compute m = reshape(t(dat),(nobj*nj),1).
compute allm = nobj*nj.
compute j = 0.
loop i = 1 to nrow(m).
  do if m(i,1) <> -9999999.
    compute j = j + 1.
    compute m(j,:)=m(i,:).
  end if.
end loop.
compute m = m(1:j,1).
compute mss = nrow(m).
compute mss = allm-mss.
compute mtmp = m.
compute mtmp(GRADE(m)) = m.
compute m = mtmp.
compute m2 = make(nrow(m),1,m(1,1)).
compute yass = csum((m = m2))/nrow(m).

do if (yass <> 1).
  compute des = design(m).
  compute uniq = ncol(des).
  compute coinc = make(uniq,uniq,0).
  compute delta = coinc.
  compute map = make(uniq,1,0).
  loop i = 1 to nrow(m).
    loop j = 1 to uniq.
      do if (des(i,j) = 1).
        compute map(j,1) = m(i,1).
      end if.
    end loop.
  end loop.
  loop i = 1 to nobj.
    loop j = 1 to nj.
      do if dat(i,j) <> -9999999.
        loop k = 1 to uniq.
          do if dat(i,j) = map(k,1).
            compute dat(i,j) = k.
            BREAK.
          end if.
        end loop.
      end if.
    end loop.
  end loop.
  compute datms = (dat <> -9999999).
  compute mu = rsum(datms).
  compute nprs = csum(mu&*(mu-1))*.5.
  compute btalp = make((btn+1),1,-999).



/* THIS CONSTRUCTS THE COINCIDENCE MATRIX FROM THE MATRIX DATA */.

  loop k = 1 to nobj.
    compute temp = make(uniq, uniq, 0).
    loop i = 1 to nj.
      loop j = 1 to nj.
        do if (dat(k,i) <> -9999999 AND dat(k,j) <> -9999999 AND i <> j).
          compute temp(dat(k,i),dat(k,j)) = temp(dat(k,i),dat(k,j)) + (1/(mu(k,1)-1)).
        end if.
      end loop.
    end loop.
    compute coinc = coinc + temp.
  end loop.
  compute q = reshape(coinc, (nrow(coinc)*ncol(coinc)), 1).
  compute q = csum(q > 0).
  compute nc = rsum(coinc).
  compute n = csum(nc).
  compute coinct = coinc.
  compute dmat = diag(coinc).
  compute nzero = csum(dmat > 0).
  compute bootm = nprs.
  compute nx = (dmat/n)&**bootm.
  compute nx=rnd(btn*csum(nx)).
  compute numone = 0.

/* THIS CONSTRUCTS THE EXPECTED MATRIX */.

  compute expect = coinc.
  loop i = 1 to uniq.
    loop j = 1 to uniq.
      do if (i = j).
        compute expect(i,j)=nc(i,1)*(nc(j,1)-1)/(n-1).
      else if (i <> j).
        compute expect(i,j)=nc(i,1)*nc(j,1)/(n-1).
      end if.
    end loop.
  end loop.



  loop z = 1 to (btn + 1).

/* HERE IS WHERE WE START DOING THE BOOTSTRAPPING */.
    do if (z > 1).
      compute btalp(z,1)=1.
      compute rchfirst=-1.
      loop u = 1 to nobj.
        compute muloop=(mu(u,1)*(mu(u,1)-1))/2.
        loop ppp= 1 to muloop.
          compute rchoose=trunc(uniform(1,1)*nprs)+1.
          do if (ppp = 2 and rchfirst=rchoose).
            compute rchoose=trunc(uniform(1,1)*nprs)+1.          
          end if.
          compute rchfirst=rchoose.
          compute btalp(z,1)=btalp(z,1)-(er(rchoose,1)/(mu(u,1)-1)).
        end loop.
      end loop.
      do if (btalp(z,1) <= -1).
        compute btalp(z,1)=-1.
      end if. 
    end if.

    do if (z = 1).
      do if (!level = 2).
        compute delta = make(uniq,uniq,0).
          loop i = 1 to uniq.
            loop j = i to uniq.
              do if (i <> j).
                compute delta(i,j) = (csum(nc(i:j,1))-(nc(i,1)/2)-(nc(j,1)/2))**2.
                compute delta(j,i) = delta(i,j).
              end if.
            end loop.
          end loop.
          compute v = {"Ordinal"}.
          do if (z = 1).
            compute deltat = delta.
          end if.
        end if.
      do if (!level = 1).
        compute delta = 1-ident(uniq).
        compute v = {"Nominal"}.
        compute deltat = delta.
      end if.
      do if (!level = 3).
        loop i = 1 to uniq.
          loop j = i to uniq.
            do if (i <> j).
              compute delta(i,j) = (map(i,1)-map(j,1))**2.
              compute delta(j,i) = delta(i,j).
            end if.
          end loop.
        end loop.
        compute v = {"Interval"}.
        compute deltat = delta.
      end if.
      do if (!level = 4).
        loop i = 1 to uniq.
          loop j = i to uniq.
            do if (i <> j).
              compute delta(i,j) = ((map(i,1)-map(j,1))/(map(i,1)+map(j,1)))**2.
              compute delta(j,i) = delta(i,j).
            end if.
          end loop.
        end loop.
        compute v = {"Ratio"}.
        compute deltat = delta.
      end if.
      compute num = csum(rsum(delta&*coinc)).
      compute den = csum(rsum(delta&*expect)).
      do if (den > 0).
        compute alp = 1-(num/den).
        compute btalp(1,1)=alp.
        compute expdis=csum(rsum((expect&*delta)))/n.
      end if.


  /* this is new */.
      compute er=make(nprs,3,0).
      compute cnt=0.
      loop k = 1 to nrow(dat).
        loop i = 1 to (ncol(dat)-1).
          loop j = (i+1) to ncol(dat).
            compute v1=dat(k,i).
            compute v2=dat(k,j).
            do if (v1 <> -9999999 and v2 <> -9999999).
              compute cnt=cnt+1.
              compute er(cnt,1:2)={v1,v2}.
              compute er(cnt,3)=delta(v1,v2).
            end if.
          end loop.
        end loop.
      end loop.
      compute er=er(:,3).
      loop i = 1 to nprs.
        compute er(i,1)=(2*er(i,1))/(expdis*csum(mu)).
      end loop.
    end if.

  end loop.
  compute alpfirst = btalp(1,1).

/* NOW WE CALCULATE CI AND P(Q) FROM BOOTSTRAPPING */.
  do if (btn > 0).
    compute btalp=btalp(2:nrow(btalp),1).

/* NOW WE SORT THE BOOTSTRAP ESTIMATES */.

    compute btalptmp = btalp.
    compute btalptmp(GRADE(btalp)) = btalp.
    compute btalp = btalptmp.
    compute btalp = btalp(1:nrow(btalp),1).
    compute mn = csum(btalp)/btn.
    compute low95 = trunc(.025*btn).
    compute high95 = trunc(.975*btn)+1.
    compute low95 = btalp(low95,1).
    compute high95 = btalp(high95,1).
    compute median = btalp(0.50*btn).
    compute q = {.9, 0; .8, 0; .7, 0; 0.67, 0; .6, 0; .5, 0}.
    loop i = 1 to 6.
      compute qcomp = (btalp < q(i,1)).
      compute qcomp = csum(qcomp)/btn.
      compute q(i,2)=qcomp.
    end loop.
  end if.
  do if (btalp(1,1) = -999).
    compute btprob = 1.
  end if.

  print/title = "Krippendorff's Alpha Reliability Estimate".
  do if (btn = 0 or btprob = 1).
    compute res = {alpfirst, nobj, nj, nprs}.
    compute lab = {"Alpha", "Units", "Obsrvrs", "Pairs"}.
  end if.
  do if (btn > 0 and btprob = 0).
    compute res = {alpfirst, low95, high95, nobj, nj, nprs}.
    compute lab = {"Alpha", "LL95%CI", "UL95%CI", "Units", "Observrs", "Pairs"}.
  end if.
  print res/title = " "/rnames = v/cnames = lab/format = F10.4.
  do if (btn > 0 and btprob = 0).
    print q/title = "Probability (q) of failure to achieve an alpha of at least alphamin:"/clabels 
    = "alphamin" "q"/format = F10.4.
    print btn/title = "Number of bootstrap samples:".
  end if.
  print vn/title = "Judges used in these computations:"/format = a8.
  do if (!detail = 1).
    print/title = "====================================================".
    print coinct/title = "Observed Coincidence Matrix"/format = F9.2.
    print expect/title = "Expected Coincidence Matrix"/format = F9.2.
    print deltat/title = "Delta Matrix"/format F9.2.
    compute tmap = t(map).
    print tmap/title "Rows and columns correspond to following unit values"/format = F9.2.
  end if.
else.
  print/title = "ERROR: Input Reliability Data Matrix Exhibits No Variation.".
end if.
do if (btprob = 1).
  print/title = "A problem was encountered when bootstrapping, so these results are not printed".
end if.
print/title = "Examine output for SPSS errors and do not interpret if any are found".
END MATRIX.
RESTORE.
!ENDDEFINE.

kalpha judges = Coder1 Coder2/level = 2/boot = 0/detail = 0.

DATASET ACTIVATE DataSet2.
PLUM Risk BY Topic_1 Topic_2 Topic_3 Topic_4 Topic_5 Topic_6 Topic_7 Topic_8 Topic_9 
    Outlet_Broadsheet Outlet_Tabloid WITH Year AIAct_Continuous
  /CRITERIA=CIN(95) DELTA(0) LCONVERGE(0) MXITER(100) MXSTEP(5) PCONVERGE(1.0E-6) SINGULAR(1.0E-8)
  /LINK=LOGIT
  /PRINT=FIT PARAMETER SUMMARY.

REGRESSION
  /DESCRIPTIVES MEAN STDDEV CORR SIG N
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA COLLIN TOL CHANGE
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT absRes
  /METHOD=ENTER Topic_1 Topic_2 Topic_3 Topic_4 Topic_5 Topic_6 Topic_7 Topic_8 Topic_9 
    Outlet_Broadsheet Outlet_Tabloid Year Topic_0 AIAct_Continuous
  /SCATTERPLOT=(*ZRESID ,*ZPRED)
  /RESIDUALS NORMPROB(ZRESID).

REGRESSION
  /DESCRIPTIVES MEAN STDDEV CORR SIG N
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS CI(95) R ANOVA COLLIN TOL CHANGE
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT Difference_Score
  /METHOD=ENTER Topic_1 Topic_2 Topic_3 Topic_4 Topic_5 Topic_6 Topic_7 Topic_8 Topic_9 
    Outlet_Broadsheet Outlet_Tabloid Year Topic_0 AIAct_Continuous.
