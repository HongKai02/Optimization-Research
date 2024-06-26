option limrow = 0, limcol = 0;
****************************** Data ******************************

set scenarios /'good', 'average', 'bad'/;
set sellablecrops /'wheat', 'corn', 'sugarbeetsexp', 'sugarbeetscheap'/;
set crops /'wheat', 'corn', 'sugarbeets'/;
set cattlefeed(crops) /'wheat', 'corn' /;

parameter alpha 'confidence level for CVaR';
alpha = [0/3];
* When alpha is 2/3, we look at the average of the worst case after that, which is just the 'bad' case. When we only look at this case, 599950 is our optimal solution, so this makes sense
* If we set aplha to be 1/3, it should look at the 2/3 tail, and the average of that should be the average of the bad and average case, but this does not seem to be the case

parameter probability(scenarios) 'probability of each scenario materializing'
/
good [1/3],
bad [1/3],
average [1/3]
/;

table yield(crops, scenarios)
                good    average     bad
wheat           3       2.5         2
corn            3.6     3           2.4
sugarbeets      24      20          16
;

parameter cost(crops) 'cost to plant per acre'
/
wheat 150,
corn 230,
sugarbeets 260
/;

*parameter price(crops)

positive variable x(crops) '# of acres of land to allocate to wheat, corn, and sugarbeets',
                  p(sellablecrops, scenarios) '# of tons of each type of crops to sell in each scenario',
                  y(cattlefeed, scenarios) '# of tons of wheat/corn to buy in each scenario',
                  u(scenarios) 'auxiliar variables for reformulation of CVaR'
                  ;
                  
variable obj, gamma;


****************************** Model ******************************


equations objective, firstStage, wheatfeeding, cornfeeding, sellingBeets, expBeets, CVaR;

objective..
    obj =e= sum(crops, cost(crops) * x(crops)) + gamma + [1/(1-alpha)] * sum(scenarios, probability(scenarios) * u(scenarios));
*    sum(scenarios, probability(scenarios) * (238 * y('wheat', scenarios) + 210 * y('corn', scenarios) - 170 * p('wheat', scenarios) - 150 * p('corn', scenarios) - 36 * p('sugarbeetsexp', scenarios) - 10 * p('sugarbeetscheap', scenarios)));
    
* Ax = b
firstStage..
    sum(crops, x(crops)) =l= 500;
    
*** There are equations binding the first stage variable with second stage variables ***
wheatfeeding(scenarios)..
    yield('wheat', scenarios) * x('wheat') + y('wheat', scenarios) - p('wheat', scenarios) =g= 200;

cornfeeding(scenarios)..
    yield('corn', scenarios) * x('corn') + y('corn', scenarios) - p('corn', scenarios) =g= 240;

sellingBeets(scenarios)..
    p('sugarbeetsexp', scenarios) + p('sugarbeetscheap', scenarios) =l= yield('sugarbeets', scenarios) * x('sugarbeets');
    
expBeets(scenarios)..
    p('sugarbeetsexp', scenarios) =l= 6000;

*** u(Omega) >= d(omega)^T y(omega) - gamma ***
CVaR(scenarios)..
    u(scenarios) =g= (238 * y('wheat', scenarios) + 210 * y('corn', scenarios) - 170 * p('wheat', scenarios) - 150 * p('corn', scenarios) - 36 * p('sugarbeetsexp', scenarios) - 10 * p('sugarbeetscheap', scenarios)) - gamma

model statisticalFarming /all/;
    
solve statisticalFarming using lp min obj;

display x.l, p.l, y.l;


