*****GAMS*****
parameter
*cost parameters
pm production cost for manufactuer/150/
k integration cost/18000/
c0 cost depende to delivery time/50/
c1 cost depende to delivery time/15/
c2 cost depende to delivery time/50/
c3 cost depende to delivery time/15/
s  shortage cost /30/
h  holding cost /3/
*sensivity parameters
beta1 price sensivity for traditional channel/0.75/
beta2 price sensivity for traditional channel/0.8/
beta3 price sensivity for traditional channel/0.75/
teta return sensivity for BODH/0.01/
delta1 the customer switch between traditional and BODH channel/0.2/
delta2 the customer switch between traditional and omni-channel/0.2/
delta3 the customer switch between BODH and omni-channel/0.2/
gama2 delivery sensivity for BODH/0.9/
gama3 delivery sensivity for omni-channel/0.5/
tu1 the customer who switch between online channels to offline channels/0.8/
tu2 the customer who switch between BODH to omni-channels/0.6/
phi2 return quantity independent of refund price for BODH/0.001/
phi3 return quantity independent of refund price for omni-channel/0.001/
phii2 sensivity refund price to return quantity for BODH/0.002/
phii3 sensivity refund price to return quantity for omni-channel/0.002/
*R return quantity of BODH
mu     Stochastic demand factor (mean)/0/
Sigma  Stochastic demand factor (Standard deviation)/3/
*o1 percentage of online customers who will to pay BOPS
*demand parameters
y1 deterministic demand for traditional retailing
y2 deterministic demand for e-tailer
y3 deterministic demand for omni-channel retailing
AL Market potential of the product/420/
n1 indexes of consumer’s channel preferences towards traditional services
n2 indexes of consumer’s channel preferences towards BODH services
n3 indexes of consumer’s channel preferences towards BOPS services/0.1/
AL1 Number of consumers who prefer the traditional respectively/200/
AL2 Number of consumers who prefer the BODH respectively/190/
AL3 Number of consumers who prefer the omni-channel respectively/420/
ALbops Number of consumers who prefer the BOPS respectively/42/
f      Probability density function
ff     Cumulative distribution function
pi /3.1415/
positive variables
p1 sale price of traditional retailing
p2 sale price of BODH retailing
p3 sale of omni-channel retailing
q1 quantity amount for traditional retailing
q2 quantity amount for e-tailer
q3 quantity amount for omni-channel
l2 delivery lead-time of BODH
l3 delivery lead-time of omni-channel
r refund price for online channels
x1
x2
x3
variable
Z;
equations
of the profit of whole supply chain
co1
co2
co3;
*y1=e=AL1-beta*p1+(tu1/2)*(gama2*l2)+delta1*(p2-p1)+delta2*(p3-p1);
*y2=e=AL2-beta*p2-gama2*l2+tu2*(l3-l2)-delta1*(p2-p1)-delta3*(p3-p2)+teta*r;
*y3=e=AL3-beta*p3-gama3*l3+(tu1/2)*(gama2*l2)-tu2*(l3-l2)-delta2*(p3-p1)-delta3*(p3-p2) + teta*r;
*R=e=phi2+phii2*pr;
*ALoc=n3*AL;
*AL=ALoc=(1-o1)ALe+ALs+ALbops
*ALbops=ALoc*n4
of..Z=e=(p1+s-pm)*q1-s*(AL1-beta1*p1+(tu1/2)*(gama2*l2)+delta1*(p2-p1)+delta2*(p3-p1))-   (p1+s+h) * (  (q1-(AL1-beta1*p1+(tu1/2)*(gama2*l2)+delta1*(p2-p1)+delta2*(p3-p1)))*(1/2) * errorf(  (    (q1 - (AL1-beta1*p1+(tu1/2)*(gama2*l2)+delta1*(p2-p1)+delta2*(p3-p1)) - x1 )/(3*sqrt(2)) ) )  -  9  * (  exp( (-sqr(q1-(AL1-beta1*p1+(tu1/2)*(gama2*l2)+delta1*(p2-p1)+delta2*(p3-p1))))/18 )  - exp ( (-sqr (x1))/18)   ) )   +     (p2+s-pm) *q2 - s*(AL2-beta2*p2-gama2*l2+tu2*(l3-l2)-delta1*(p2-p1)-delta3*(p3-p2)+teta*r)-(p2+s+h)*  (  (q2-(AL2-beta2*p2-gama2*l2+tu2*(l3-l2)-delta1*(p2-p1)-delta3*(p3-p2)+teta*r))*(1/2) * errorf(  (    (q2 - (AL2-beta2*p2-gama2*l2+tu2*(l3-l2)-delta1*(p2-p1)-delta3*(p3-p2)+teta*r) - x2 )/(3*sqrt(2)) ) )  -  9  * (  exp( (-sqr(q2-(AL2-beta2*p2-gama2*l2+tu2*(l3-l2)-delta1*(p2-p1)-delta3*(p3-p2)+teta*r)))/18 )  - exp ( (-sqr (x2))/18)   ) )  - (phi2+phii2*r)*(r+h) - sqr(c0-c1*l2) +(p3+s-pm)*q3 - s*(AL3-beta3*p3-gama3*l3+(tu1/2)*(gama2*l2)-tu2*(l3-l2)-delta2*(p3-p1)-delta3*(p3-p2) + teta*r)-(p3+s+h) *  (  (q3-(AL3-beta3*p3-gama3*l3+(tu1/2)*(gama2*l2)-tu2*(l3-l2)-delta2*(p3-p1)-delta3*(p3-p2) + teta*r))*(1/2) * errorf(  (    (q3 - (AL3-beta3*p3-gama3*l3+(tu1/2)*(gama2*l2)-tu2*(l3-l2)-delta2*(p3-p1)-delta3*(p3-p2) + teta*r) - x3 )/(3*sqrt(2)) ) )  -  9  * (  exp( (-sqr(q3-(AL3-beta3*p3-gama3*l3+(tu1/2)*(gama2*l2)-tu2*(l3-l2)-delta2*(p3-p1)-delta3*(p3-p2) + teta*r)))/18 )  - exp ( (-sqr (x3))/18)   ) ) - (phi3+phii3*r)*(r+h) - sqr(c2-c3*l3) - k;
co1..r=l=p2;
co2..p2=l=p1;
co3..p3=l=p1;
model example/all/;
solve example using NLP max Z;
display p1.l,p2.l,p3.l,q1.l,q2.l,q3.l,r.l,l2.l,l3.l,x1.l,x2.l,x3.l;


