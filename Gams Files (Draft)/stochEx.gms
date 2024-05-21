sets
i 'factories' /f1*f3/
j 'distribution centers' /d1*d5/
s 'scenarios' /lo,mid,hi/
;
parameter capacity(i) /f1 500, f2 450, f3 650/;
table demand(j,s) ’possible outcomes for demand’
    lo mid hi
d1 150 160 170
d2 100 120 135
d3 250 270 300
d4 300 325 350
d5 600 700 800
;
parameter prob(s) ’probabilities’ /
lo 0.25
mid 0.5
hi 0.25 /;
table transcost(i,j) ’unit transportation cost’
    d1  d2   d3   d4   d5
f1 2.49 5.21 3.76 4.85 2.07
f2 1.46 2.54 1.83 1.86 4.76
f3 3.26 3.08 2.60 3.76 4.45
 ;
scalar prodcost ’unit production cost’ /14/;
scalar price ’sales price’ /24/;
scalar wastecost ’cost of removal of overstocked products’ /4/;
*----------------------------------------------------------------------- * Form the Benders master problem *-----------------------------------------------------------------------
set
iter ’max Benders iterations’ /iter1*iter25/
dyniter(iter) ’dynamic subset’
;

positive variables
ship(i,j) ’shipments’
product(i) ’production’
slackproduct(i) ’slack’
received(j) ’quantity sent to market’
;


free variables
zmaster ’objective variable of master problem’
theta ’extra term in master obj’
;

equations
    masterobj ’master objective function’
    production(i) ’calculate production in each factory’
    receive(j) ’calculate quantity to be send to markets’
    prodcap(i) ’production capacity’
    optcut(iter) ’Benders optimality cuts’
;

parameter
    cutconst(iter) ’constants in optimality cuts’
    cutcoeff(iter,j) ’coefficients in optimality cuts’
;


masterobj..
    zmaster =e= sum((i,j), transcost(i,j)*ship(i,j))
    +sum(i,prodcost*product(i)) + theta;
    
receive(j)..        received(j) =e= sum(i, ship(i,j));
production(i)..     product(i) =e= sum(j, ship(i,j));
prodcap(i)..        product(i) + slackproduct(i) =e= capacity(i); 
optcut(dyniter)..   theta =g= cutconst(dyniter) +
                        sum(j, cutcoeff(dyniter,j)*received(j));

model masterproblem /masterobj, receive, production, prodcap, optcut/;


*-----------------------------------------------------------------------
* Form the Benders’ subproblem
* Notice in equation selling we use the level value received.l, i.e.
* this is a constant
*-----------------------------------------------------------------------
positive variables
sales(j) ’sales (actually sold)’
waste(j) ’overstocked products’
slacksales(j) ’slack’
;

free variables
    zsub ’objective variable of sub problem’
;

equations
    subobj ’subproblem objective function’
    selling(j) ’part of received is sold’
    selmax(j) ’upperbound on sales’
;

parameter demnd(j) ’demand’;


subobj..
    zsub =e= -sum(j, price*sales(j)) + sum(j, wastecost*waste(j));
selling(j).. sales(j) + waste(j) =e= received.l(j);
selmax(j).. sales(j) + slacksales(j) =e= demnd(j);

model subproblem /subobj,selling,selmax/;
*-------------------------------------------------------------------
* Benders algorithm
*-------------------------------------------------------------------

*
* step 1: solve master without cuts
*
dyniter(iter) = NO;
cutconst(iter) = 0;
cutcoeff(iter,j) = 0;
theta.fx = 0;
solve masterproblem minimizing zmaster using lp;

*
* repair bounds
*

theta.lo = -INF;
theta.up = INF;

scalar lowerbound /-INF/;
scalar upperbound /INF/;
parameter objsub(s);
scalar objmaster;
objmaster = zmaster.l;

option limrow = 0;
option limcol = 0;
subproblem.solprint = 0;
masterproblem.solprint = 0;

loop(iter,
*
* solve subproblems
*
    dyniter(iter) = yes;
    loop(s,
        demnd(j) = demand(j,s);
        solve subproblem minimizing zsub using lp;
        objsub(s) = zsub.l;
        cutconst(iter) = cutconst(iter) - prob(s)*sum(j,(-selmax.m(j))*demand(j,s));
        cutcoeff(iter,j) = cutcoeff(iter,j) - prob(s)*(-selling.m(j));
    );
    
    upperbound = min(upperbound, objmaster + sum(s, prob(s)*objsub(s)));

*
* convergence test
*
    display lowerbound,upperbound;
    display lowerbound;
    abort$( (upperbound-lowerbound) < 0.01*(1+abs(lowerbound)) ) "Converged";

*
* solve masterproblem
*
solve masterproblem minimizing zmaster using lp;
display ship.l;
lowerbound = zmaster.l;
objmaster = zmaster.l-theta.l;
display objmaster;
);