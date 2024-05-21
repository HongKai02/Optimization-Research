* free variable: Unconcstrained decision variables, the objective is always the free variable

option limrow = 0, limcol = 0;

positive variable x1, x2;

variable obj;

equations objective, constraint1, constraint2, con3, con4, obj2;

objective..
    obj =e= 70*x1 + 90*x2;
    
constraint1..
    4*x1 + 3*x2 =l= 40;

constraint2..
    4*x1 + 7*x2 =l= 802;

model hw1_1 /objective, constraint1, constraint2/;
    
solve hw1_1 using lp max obj;

display obj.l, x1.l, x2.l;


