% bx=[1.0, 1.5, 3.7, 4.3, 5.6, 6.1];
% by=[0.0, 1.3, 0.5, -1.0, 2.1, 3.5 ];
bx=[0,1,1,0.5,0,-1];
by=[0,0,1,2,4,6 ];

ex=load('hehe1.csv');
ey=load('hehe2.csv');

plot(bx,by);
plot(bx,by,'ro');
hold;
plot(ex,ey);