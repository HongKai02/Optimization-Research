option limrow = 0, limcol = 0;
****************************** Data ******************************

set scenarios /'good', 'average', 'bad'/;
set crops /'wheat', 'corn', 'sugarbeets'/;


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

positive variable x(crops) '# of acres of land to allocate to wheat, corn, and sugarbeets',
                  p(crops, scenarios) '# of tons of each type of crops to sell in each scenario',
                  y(crops, scenarios) '# of tons of wheat/corn to buy in each scenario'
                  ;
                  
variable obj;

****************************** Model ******************************

equations objective, firstStage, wheatfeeding, cornfeeding, sellingBeets;

objective..
    obj =e= sum(crops, cost(crops) * x(crops)) + sum(scenarios, probability(scenarios) * (238 * y('wheat', scenarios) + 210 * y('corn', scenarios) - 170 * p('wheat', scenarios) - 150 * p('corn', scenarios) - 36 * p('sugarbeets', scenarios) ));
    
firstStage..
    sum(crops, x(crops)) =l= 500;
    
wheatfeeding(scenarios)..
    yield('wheat', scenarios) * x('wheat') + y('wheat', scenarios) - p('wheat', scenarios) =g= 200;

cornfeeding(scenarios)..
    yield('corn', scenarios) * x('corn') + y('corn', scenarios) - p('corn', scenarios) =g= 240;

sellingBeets(scenarios)..
    p('sugarbeets', scenarios) =l= yield('sugarbeets', scenarios) * x('sugarbeets');


model statisticalFarming /all/;
    
solve statisticalFarming using lp min obj;

display x.l, p.l, y.l;


