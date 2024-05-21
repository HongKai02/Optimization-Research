sets
crops /'wheat', 'corn', 'sugarbeets'/
scenarios /lo,mid,hi/
;

table yield(crops,scenarios) ’possible outcomes for demand’
            lo  mid hi
wheat       2   2.5 3
corn        2.4 3   3.6
sugarbeets  16  20  24
;

parameter alpha 'confidence level for CVaR';
alpha = [2/3];


parameter probability(scenarios) ’probability of each scenario materializing’
/
lo [1/3]
mid [1/3]
hi [1/3]
/;

parameter cost(crops) 'cost to plant per acre'
/
wheat 150,
corn 230,
sugarbeets 260
/;

*----------------------------------------------------------------------- * Form the Benders master problem *-----------------------------------------------------------------------
************************** Sets and Variables Needed to Perform Benders **************************
sets iter ’max Benders iterations’ /iter1*iter2500/,
     dyniter(iter) ’dynamic subset’,
     dynscen(scenarios) 'dynamic subset of scenarios',
     curiter(iter) 'dynamic subset of iter',
     curscen(scenarios) 'scenario that the problem should look at'
;

parameters oldgamma(iter),
          oldX(crops, iter),
          oldU(scenarios, iter),
          pi(scenarios, crops, iter),
          phi(scenarios, iter);

************************** Initial Master Problem **************************
positive variables x(crops) '# of acres of land to allocate to wheat, corn, and sugarbeets';

free variables zmaster ’objective variable of master problem’,
               theta ’extra term in master obj’,
               gamma 'VaR'
;

equations masterobj ’master objective function’,
          plantamount 'constraint on first stage variables'
;

masterobj..
    zmaster =e= sum(crops, cost(crops) * x(crops)) + gamma + [1/(1-alpha)] * theta;
    
plantamount..
    sum(crops, x(crops)) =l= 500;

gamma.lo = -50000;
theta.lo = -50000;

model initialmasterproblem /masterobj, plantamount/;

************************** Subproblem **************************
positive variables p(crops, scenarios) '# of tons of each type of crops to sell in each scenario',
                   y(crops, scenarios) '# of tons of wheat/corn to buy in each scenario',
                   u(scenarios) 'auxiliar variables for reformulation of CVaR',
                   subX(crops)
;

free variables zsub,
               subgamma;

equations subobj(scenarios),
          CVaRConstraint(scenarios),
          wheatfeeding(scenarios),
          cornfeeding(scenarios),
          sellingBeets(scenarios),
          setx(crops, iter),
          setgamma,
          CVaRConstraintDebug(scenarios)
          ;

* Lots of equations here are indexed by scenarios, but there should only be one active scenario at a time.
subobj(curscen)..
    zsub =e= probability(curscen) * u(curscen);
    
CVaRConstraint(curscen)..
    u(curscen) =g= (238 * y('wheat', curscen) + 210 * y('corn', curscen) - 170 * p('wheat', curscen) - 150 * p('corn', curscen) - 36 * p('sugarbeets', curscen)) - subgamma;
    
CVaRConstraintDebug(curscen)..
    u(curscen) =g= (238 * y('wheat', curscen) + 210 * y('corn', curscen) - 170 * p('wheat', curscen) - 150 * p('corn', curscen) - 36 * p('sugarbeets', curscen));

wheatfeeding(curscen)..
    yield('wheat', curscen) * subX('wheat') + y('wheat', curscen) - p('wheat', curscen) =g= 200;

cornfeeding(curscen)..
    yield('corn', curscen) * subX('corn') + y('corn', curscen) - p('corn', curscen) =g= 240;

sellingBeets(curscen)..
    p('sugarbeets', curscen) =l= yield('sugarbeets', curscen) * subX('sugarbeets');

setx(crops, curiter)..
    subX(crops) =e= oldX(crops, curiter);

setgamma(curiter)..
    subgamma =e= oldgamma(curiter);
    
model subproblem /subobj, CVaRConstraint, wheatfeeding, cornfeeding, sellingBeets, setx, setgamma/;

************************** Master Problem **************************
equation cutcon(iter);

cutcon(dyniter)..
    theta =g= sum(scenarios, probability(scenarios) * oldU(scenarios, dyniter)) +
              sum(scenarios, sum(crops, pi(scenarios, crops, dyniter) * (x(crops) - oldx(crops, dyniter)))
              + phi(scenarios, dyniter) * (gamma - oldgamma(dyniter)))

model masterproblem /masterobj, plantamount, cutcon/;

************************** Benders Algorithm **************************
*
* Set initial values for bounds
*

*theta.lo = -INF;
*theta.up = INF;
*
* Step 1: Solve master problem without cuts
*
dyniter(iter) = NO;
curiter(iter) = NO;
curscen(scenarios) = NO;
*theta.fx = 0;
solve initialmasterproblem minimizing zmaster using lp;
display x.l;

scalar lowerbound /-INF/;
scalar upperbound /INF/;
parameter objsub(scenarios);
scalar objmaster;
objmaster = zmaster.l;

* Save the x and gamma values
oldX(crops, 'iter1') = x.l(crops);
oldgamma('iter1') = gamma.l;

*option limrow = 0;
*option limcol = 0;
*subproblem.solprint = 0;
*masterproblem.solprint = 0;

*
* Step 2: Solve subproblems
*
scalar counter /0/;
loop(iter,
    dyniter(iter) = yes;
    curiter(iter) = yes;
    loop(scenarios,
        curscen(scenarios) = yes;
        solve subproblem minimizing zsub using lp;
* Save values from this subproblem
        display scenarios, u.l;
        oldU(curscen, iter) = u.l(curscen);
        display "WHATDUPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP";
        display oldU, subgamma.l, y.l, p.l;
        pi(scenarios, crops, iter) = setx.m(crops, iter);
        phi(scenarios, iter) = setgamma.m(iter);
        objsub(scenarios) = zsub.l;
*        abort "Ran one subproblem";
        curscen(scenarios) = no;
    );
    
upperbound = min(upperbound, objmaster + [1/(1-alpha)] * sum(scenarios, probability(scenarios) * oldU(scenarios, iter)) + oldgamma(iter));

*
* Step 3: Check for convergence
*
    display lowerbound,upperbound;
    abort$(abs(upperbound-lowerbound) < 1000) "Converged";
*    abort$( (upperbound-lowerbound) < 0.01*(1+(lowerbound)) ) "Converged";

*
* Step 4: Solve masterproblem, and go back to step 2
*
    solve masterproblem minimizing zmaster using lp;
* Save values from master problem
    oldX(crops, iter) = x.l(crops);
    oldgamma(iter) = gamma.l;
    lowerbound = zmaster.l;
    objmaster = zmaster.l-theta.l;
    display objmaster;
* Update curiter
    curiter(iter) = no;
    counter = counter + 1;
    display lowerbound, upperbound;
    abort$(counter = 100) "Done with 10 iter";
);