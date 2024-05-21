option limrow = 0, limcol = 0;

positive variable x, y;

variable obj;

equations objective, constraint1, constraint2, constraint3, constraint4;

objective..
    obj =e= -y - [x/4];
    
constraint1..
    y - x =l= 5;

constraint2..
    y - [1/2]*x =l= [15/2];

constraint3..
    y + [1/2]*x =l= [35/2];
        
constraint4..
    -y + x =l= 10;
    
x.up = 16;

model original /objective, constraint1, constraint2, constraint3, constraint4/;
    
solve original using lp minimizing obj;


