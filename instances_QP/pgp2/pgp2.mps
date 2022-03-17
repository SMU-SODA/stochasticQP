* This Power Generation Problem comes from J.L. Higle and S. Sen,
* Stochastic Decomposition, Kluwer Academic Publishers 1996, who credit
* F.V. Louveaux and Y. Smeers, �Optimal investment for electricity generation:
* A stochastic model and a test problem�, in: Yu. Ermoliev and R.J-B Wets (eds.)
* Numerical Techniques for Stochastic Optimization, Springer Verlag, 1988.
*
*23412345678  12345678  123456789012   12345678  123456789012
NAME          PGP2
ROWS
 N  FOBJ
 G  MXDEMD
 L  BUDGET
 L  CAPEQ1
 L  CAPEQ2
 L  CAPEQ3
 L  CAPEQ4
 G  DNODE1
 G  DNODE2
 G  DNODE3
COLUMNS
*23412345678  12345678  123456789012   12345678  123456789012
    INVEQ1    FOBJ         10.0        MXDEMD       1.0
    INVEQ1    BUDGET       10.0        CAPEQ1      -1.0
    INVEQ2    FOBJ          7.0        MXDEMD       1.0
    INVEQ2    BUDGET        7.0        CAPEQ2      -1.0
    INVEQ3    FOBJ         16.0        MXDEMD       1.0
    INVEQ3    BUDGET       16.0        CAPEQ3      -1.0
    INVEQ4    FOBJ          6.0        MXDEMD       1.0
    INVEQ4    BUDGET        6.0        CAPEQ4      -1.0
    EQ1ND1    FOBJ         40.0        CAPEQ1       1.0
    EQ1ND1    DNODE1        1.0
    EQ1ND2    FOBJ         24.0        CAPEQ1       1.0
    EQ1ND2    DNODE2        1.0
    EQ1ND3    FOBJ          4.0        CAPEQ1       1.0
    EQ1ND3    DNODE3        1.0
    EQ2ND1    FOBJ         45.0        CAPEQ2       1.0
    EQ2ND1    DNODE1        1.0
    EQ2ND2    FOBJ         27.0        CAPEQ2       1.0
    EQ2ND2    DNODE2        1.0
    EQ2ND3    FOBJ          4.5        CAPEQ2       1.0
    EQ2ND3    DNODE3        1.0
    EQ3ND1    FOBJ         32.0        CAPEQ3       1.0
    EQ3ND1    DNODE1        1.0
    EQ3ND2    FOBJ         19.2        CAPEQ3       1.0
    EQ3ND2    DNODE2        1.0
    EQ3ND3    FOBJ          3.2        CAPEQ3       1.0
    EQ3ND3    DNODE3        1.0
    EQ4ND1    FOBJ         55.0        CAPEQ4       1.0
    EQ4ND1    DNODE1        1.0
    EQ4ND2    FOBJ         33.0        CAPEQ4       1.0
    EQ4ND2    DNODE2        1.0
    EQ4ND3    FOBJ          5.5        CAPEQ4       1.0
    EQ4ND3    DNODE3        1.0
    PEN1      FOBJ       1000.0        CAPEQ1      -1.0
    PEN2      FOBJ       1000.0        CAPEQ2      -1.0
    PEN3      FOBJ       1000.0        CAPEQ3      -1.0
    PEN4      FOBJ       1000.0        CAPEQ4      -1.0
RHS
    RHS       MXDEMD       15.0
    RHS       BUDGET      220.0
    RHS       DNODE1        5.0
    RHS       DNODE2        4.0
    RHS       DNODE3        3.0
BOUNDS
  UP bnd       EQ1ND1       11.2
  UP bnd       EQ1ND2       10.6
  UP bnd       EQ1ND3       11.2
  UP bnd       EQ2ND1       10.2
  UP bnd       EQ2ND2       11.2
  UP bnd       EQ2ND3       11.2
  UP bnd       EQ3ND1       11.2
  UP bnd       EQ3ND2       11.2
  UP bnd       EQ3ND3       11.2
  UP bnd       EQ4ND1       11.2
  UP bnd       EQ4ND2       11.2
  UP bnd       EQ4ND3       11.2
  LO bnd       EQ1ND1       0.1
  LO bnd       EQ1ND2       0.1
  LO bnd       EQ1ND3       0.1
  LO bnd       EQ2ND1       0.1
  LO bnd       EQ2ND2       0.1
  LO bnd       EQ2ND3       0.1
  LO bnd       EQ3ND1       0.1
  LO bnd       EQ3ND2       0.1
  LO bnd       EQ3ND3       0.1
  LO bnd       EQ4ND1       0.1
  LO bnd       EQ4ND2       0.1
  LO bnd       EQ4ND3       0.1
QUADOBJ
    EQ1ND1       EQ1ND1             11.  
    EQ1ND1       EQ1ND2             11.
    EQ1ND1       EQ1ND3             11. 
    EQ1ND2       EQ1ND2             11. 
    EQ1ND2       EQ1ND3             11.  
    EQ1ND3       EQ1ND3             11.
    EQ2ND1       EQ2ND1             11.  
    EQ2ND1       EQ2ND2             11.
    EQ2ND1       EQ2ND3             11. 
    EQ2ND2       EQ2ND2             11. 
    EQ2ND2       EQ2ND3             11.  
    EQ2ND3       EQ2ND3             11.
    EQ3ND1       EQ3ND1             11.  
    EQ3ND1       EQ3ND2             11.
    EQ3ND1       EQ3ND3             11. 
    EQ3ND2       EQ3ND2             11. 
    EQ3ND2       EQ3ND3             11.  
    EQ3ND3       EQ3ND3             11.
    EQ4ND1       EQ4ND1             11.  
    EQ4ND1       EQ4ND2             11.
    EQ4ND1       EQ4ND3             11. 
    EQ4ND2       EQ4ND2             11. 
    EQ4ND2       EQ4ND3             11.  
    EQ4ND3       EQ4ND3             11.
ENDATA
