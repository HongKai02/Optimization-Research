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
alpha = [1/3];


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

parameter weight;
weight = 0.01;

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
          phi(scenarios, iter),
          oldY(crops, scenarios, iter),
          oldP(crops, scenarios, iter)
          ;
          

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
    zmaster =e= sum(crops, cost(crops) * x(crops)) + weight * gamma + theta;
    
plantamount..
    sum(crops, x(crops)) =l= 500;

* Set lower bounds for gamma and theta (or else they'd go to -infinity)
gamma.lo = -500000000000;
theta.lo = -500000000000;

model initialmasterproblem /masterobj, plantamount/;

************************** Subproblem **************************
positive variables p(crops, scenarios) '# of tons of each type of crops to sell in each scenario',
                   y(crops, scenarios) '# of tons of wheat/corn to buy in each scenario',
                   u(scenarios) 'auxiliary variables for reformulation of CVaR',
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
          setgamma
;

* Lots of equations here are indexed by scenarios, but there should only be one active scenario at a time.
subobj(curscen)..
    zsub =e= (weight) * ([1/(1-alpha)] * (probability(curscen) * u(curscen))) + (1-weight) * (probability(curscen) * (238 * y('wheat', curscen) + 210 * y('corn', curscen) - 170 * p('wheat', curscen) - 150 * p('corn', curscen) - 36 * p('sugarbeets', curscen)));
    
CVaRConstraint(curscen)..
    u(curscen) =g= (238 * y('wheat', curscen) + 210 * y('corn', curscen) - 170 * p('wheat', curscen) - 150 * p('corn', curscen) - 36 * p('sugarbeets', curscen)) - subgamma;
    
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
    theta =g= sum(scenarios, weight * [1/(1-alpha)] * (probability(scenarios) * oldU(scenarios, dyniter))) +
              sum(scenarios, (1-weight) * probability(scenarios) * (238 * oldY('wheat', scenarios, dyniter) + 210 * oldY('corn', scenarios, dyniter) - 170 * oldP('wheat', scenarios, dyniter) - 150 * oldP('corn', scenarios, dyniter) - 36 * oldP('sugarbeets', scenarios, dyniter))) + 
              sum(crops, sum(scenarios, pi(scenarios, crops, dyniter) * (x(crops) - oldx(crops, dyniter)))) +
              sum(scenarios, phi(scenarios, dyniter) * (gamma - oldgamma(dyniter)));

model masterproblem /masterobj, plantamount, cutcon/;

************************** Benders Algorithm **************************
*
* Step 1: Solve master problem without cuts
*
dyniter(iter) = NO;
curiter(iter) = NO;
curscen(scenarios) = NO;
solve initialmasterproblem minimizing zmaster using lp;
display zmaster.l;

scalar lowerbound /-INF/;
scalar upperbound /INF/;
parameter objsub(scenarios);
scalar objmaster;
objmaster = zmaster.l - theta.l;
display objmaster;

* Save the x and gamma values
oldX(crops, 'iter1') = x.l(crops);
oldgamma('iter1') = gamma.l;


oldU(scenarios, 'iter1') = 0;
pi(scenarios, crops, 'iter1') = 0;
phi(scenarios, 'iter1') = 0;
oldY(crops, scenarios, 'iter1') = 0;
oldP(crops, scenarios, 'iter1') = 0;

scalar counter /1/;
    loop(iter,
    if (ord(iter) <> 1,    
        solve masterproblem minimizing zmaster using lp;
* Save values from master problem
        oldX(crops, iter) = x.l(crops);
        oldgamma(iter) = gamma.l;
        lowerbound = zmaster.l;
        objmaster = zmaster.l-theta.l;
        display lowerbound, x.l, oldY, oldP, pi, phi, gamma.l, theta.l;
* Update curiter
        counter = counter + 1;
    );
    
*
* Step 2: Solve subproblems
*
    dyniter(iter) = yes;
    curiter(iter) = yes;
    loop(scenarios,
        curscen(scenarios) = yes;
        solve subproblem minimizing zsub using lp;
        display x.l;
        display p.l, y.l, gamma.l, theta.l;
* Save values from this subproblem
        oldU(curscen, iter) = u.l(curscen);
        pi(scenarios, crops, iter) = setx.m(crops, iter);
        phi(scenarios, iter) = setgamma.m(iter);
        oldY(crops, scenarios, iter) = y.l(crops, scenarios);
        oldP(crops, scenarios, iter) = p.l(crops, scenarios);
        objsub(scenarios) = zsub.l;
        curscen(scenarios) = no;
    );
    
upperbound = min(upperbound, objmaster + (weight) * ([1/(1-alpha)] * sum(scenarios, probability(scenarios) * oldU(scenarios, iter))) + (1-weight) * (sum(scenarios, probability(scenarios) * (238 * oldY('wheat', scenarios, iter) + 210 * oldY('corn', scenarios, iter) - 170 * oldP('wheat', scenarios, iter) - 150 * oldP('corn', scenarios, iter) - 36 * oldP('sugarbeets', scenarios, iter)))));

*
* Step 3: Check for convergence
*
    display lowerbound,upperbound, counter, oldU, oldY, oldP;
    abort$(abs(upperbound-lowerbound) < 10) "Converged";
    curiter(iter) = no;
);