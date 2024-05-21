* free variable: Unconcstrained decision variables, the objective is always the free variable

option limrow = 0, limcol = 0;

positive variable x1, x2, x3;

set scenarios /bad, average, good/;

parameters delta(scenarios)

positive variable u(scenarios);

variable obj;

equations objective, constraint1, constraint2, con3, con4, obj2;

objective..
    obj =e= 150*x1 + 230*x2 + 260 * x3 + gamma + [1/(1-alpha)] * sum(scenarios, u(scenarios));
    
* Ax = b
constraint1..
    x1 + x2 + x3 =l= 500;

* u(omega) >= d(omega)^T y(omega) - gamma
constraint2(scenarios)..
    u(scenarios) =g= 

constraint2..
    4*x1 + 7*x2 =l= 802;

model hw1_1 /objective, constraint1, constraint2/;
    
solve hw1_1 using lp max obj;

display obj.l, x1.l, x2.l;

