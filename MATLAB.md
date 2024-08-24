%case (5) 


clc;
clear;
close all;
tic,

%%Sets
        I  = 1;       %Index set of KKT Cases
        J  = 2;       %Delivery time
        H  = 3;       %Index set of channels

        
%%Parameters
%%%cost parameters     
e= exp(1);
pm=150;
k=18000;
c0=50;
c1=15;
c2=50;
c3=15;
s =30;
h=3;

  %%%sensevity parameters
beta   = 0.75;    %price sensivity for channels
teta   = 0.11;    %return sensivity for BODH
delta1 = 0.2 ;    %the customer switch between traditional and BODH channel
delta2 = 0.2 ;    %the customer switch between traditional and omni-channel
delta3 = 0.2 ;    %the customer switch between BODH and omni-channel
gama2  = 0.7 ;    %delivery sensivity for BODH
gama3  = 0.5 ;    %delivery sensivity for omnichannel
tu1    = 0.8 ;    %the customer who switch between online channels to offline channels
tu2    = 0.5 ;    %the customer who switch between BODH to omni-channels
phi    = 0.01;    %return quantity independent of refund price for BODH 
phii   = 0.02;    %sensivity refund price to return quantity for BODH

  %%%demand parameters
Alpha1 = 300;     %Number of consumers who prefer the traditional respectively
Alpha2 = 250;     %Number of consumers who prefer the BODH respectively
Alpha3 = 360;     %Number of consumers who prefer the omnichannel respectively
mu     = 0  ;     %Stochastic demand factor (mean)
Sigma  = 10 ;     %Stochastic demand factor (Standard deviation) 
  
        
  %%Variables        
         syms r  ;              %Refund price for online market
         syms p1 ;              %The price for traditional channel 
         syms p2 ;              %The price for BODH channel
         syms p3 ;              %The price for omni-channel
         syms l1 ;              %delivery time for BODH channel
         syms l2 ;              %delivery time for omni-channel
   
    %%Equations
    y = zeros(H,I); y = poly2sym(y);        %Price dependent deterministic demand in channel h   
    ru= zeros(H,I); ru = poly2sym(ru);      %Only for simplification 
    
    l1=((p1-pm)*(tu1/2)*gama2+(p2-pm)*(-gama2-tu2)+(p3-pm)*((tu1/2)*gama2+tu2)+(2*c1*c0))/(2*(c1)^2);  %sensivity for delivery lead time for BODH channel
    l2=((p2-pm)*(tu2)+((p3-pm)*(-gama3-tu2))+(2*c3*c2))/(2*(c3)^2);                    %sensivity for delivery lead time for omni-channel    

    
    
    for i = 1 : I
      
       y(1,i)=Alpha1 - beta*p1 + ((tu1/2)*gama2*l1)+delta1*(p2-p1) + delta2*(p3-p1);
       y(2,i)=Alpha2 - beta*p2 - gama2*l1 + tu2*(l2-l1) - delta1*(p2-p1) + delta3*(p3-p2) + teta*r;
       y(3,i)=Alpha3 - beta*p3 - gama3*l2 + ((tu1/2)*(gama2*l1)) - tu2*(l2-l1) - delta2*(p3-p1) - delta3*(p3-p2);
        
        
        ru(1,i)= (p1+s-pm)/(p1+s+h);
        ru(2,i)= (p2+s-pm)/(p2+s+h);
        ru(3,i)= (p3+s-pm)/(p3+s+h);
    end    
        
        
  %% Solution method

  %%Save results at:
    empty_individual.Price=[];
    empty_individual.Pricer=[];
    empty_individual.Quantity=[];
    empty_individual.deliverytime=[];
    empty_individual.TRevenue=[];
    Results = repmat(empty_individual,I,1);             % Results for each I  
    
  %%The optimal order quantities in each I
    F = zeros(H,I); F = poly2sym(F);    
    q = zeros(H,I); q = poly2sym(q);    
   
    for i = 1 : I
        for h = 1 : H
            F(h,i)=.5*(1+erf((ru(h,i)-mu)/(Sigma*(2^.5))));  % Cumulative probability distribution function (Normal) for price-dependent stochastic demand in market segment h from zone i.
            q(h,i)=y(h,i)+inv(F(h,i));                       % Order quantity for market segment h in zone i.
        end
    end
    
                                        

%% Pricing model
    syms x;
    Zi =.5*(1+erf((x-mu)/(Sigma*(2^.5))));   % Stochastic demand factor (Normal distribution) for channel i
     
    TR  =zeros(1,I)  ;        % Total Revenue for each Case  
    p   =zeros(H,I)  ;        % The price for channel h for each Case
    Q   =zeros(H,I)  ;        % Total order quantity for each Case
    l   =zeros(J,I)  ;        % delivery lead-time for channele h for each Case
    r   =zeros(1,I)  ;        % Returned quantity for channele h for each Case
    Landa =zeros(3,I) ;       % Lagrangian multiper for each Case

 
 
    %%%%%% Case (5)
    
       i= 1;
       syms p1; syms p2;  syms p3; syms r;   
       EQ1=y(1,i)+ inv(F(1,i))+ y(3,i)+ inv(F(3,i))- int((Zi),-Sigma*sqrt(6),inv(F(1,i))) - int((Zi),-Sigma*sqrt(6),inv(F(3,i))) + delta3*(p2-p1)+ beta*((2*pm)-p1-p1)+ delta1*(p2-p1);
       EQ2=y(2,i)+inv(F(2,i))+ delta1*(p1-p2) + delta3*(p1-p2) + beta*(pm-p2)- int((Zi),-Sigma*sqrt(6),inv(F(2,i)));
       EQ3= ((p2-pm)*teta)-(phii*((2*r)+h))-phi;
       EQ4= p1-p3;     
       EQ5= abs(p2-r) -(p2-r); 
       EQ6= abs(p1-p2)-(p1-p2);
       EQ7= abs(p1-p3)-(p1-p3);
       
       Gsym = (EQ1).^2+(EQ2).^2+(EQ3).^2+(EQ4).^2+(EQ5).^2+(EQ6).^2+(EQ7).^2;       
       Gstr = sym2str(Gsym);
       Gstr = replace(Gstr,{'p1','p2','p3'},{'x(1)','x(2)','x(3)'});
       Gstr = strrep(Gstr,'erf','EEE');
       Gstr = strrep(Gstr,'r','x(4)');
       Gstr = strrep(Gstr,'EEE','erf');     
       G = str2func(strcat('@(x)',Gstr));
       
       lb = [  0,   0,   0,   0];
       ub = [400, 400, 400, 400];
       options = optimoptions('particleswarm','display','off');
       Solution = particleswarm(G,4,lb,ub,options);
       
       if (G(Solution)>1e-2)
           disp('  This case is infeasible.');
           fprintf('\n');
       end
       disp(Solution);
       fprintf('  Error = %s\n',G(Solution));
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       p(1,i) = Solution(1); 
       p(2,i) = Solution(2); 
       p(3,i) = Solution(3);
       r(1,i) = Solution(4);
            
       y(1,i)=max(0,(Alpha1 - beta*p(1,i) + ((tu1/2)*(gama2*l(1,i)))+delta1*(p(2,i)-p(1,i)) + delta2*(p(3,i)-p(1,i))));
       y(2,i)=max(0,(Alpha2 - beta*p(2,i) - gama2*l(1,i) + tu2*(l(2,i)-l(1,i)) - delta1*(p(2,i)-p(1,i)) + delta3*(p(3,i)-p(2,i)) + teta*r(1,i)));
       y(3,i)=max(0,(Alpha3 - beta*p(3,i) - gama3*l(2,i) + ((tu1/2)*(gama2*l(1,i))) - tu2*(l(2,i)-l(1,i)) - delta2*(p(3,i)-p(1,i)) - delta3*(p(3,i)-p(2,i))));
        
       ru(1,i)= (p(1,i)+s-pm)/(p(1,i)+s+h);
       ru(2,i)= (p(2,i)+s-pm)/(p(2,i)+s+h);
       ru(3,i)= (p(3,i)+s-pm)/(p(3,i)+s+h);  
        
     
       for h = 1 : H
            F(h,i)=.5*(1+erf((ru(h,i)-mu)/(Sigma*(2^.5))));  % Cumulative probability distribution function for price-dependent stochastic demand in market segment h from zone i.
            y(h,i)= sym2poly(y(h,i));
            q(h,i)=y(h,i)+inv(F(h,i));        Q(h,i)= sym2poly(q(h,i));      % Order quantity for market segment h in zone i.
       end
                                               
    
    
    l(1,i)=((p(1,i)-pm)*(tu1/2)*gama2+(p(2,i)-pm)*(-gama2-tu2)+(p(3,i)-pm)*((tu1/2)*gama2+tu2)+(2*c1*c0))/(2*(c1)^2);  %sensivity for delivery lead time for BODH channel
    l(2,i)=((p(2,i)-pm)*(tu2)+((p(3,i)-pm)*(-gama3-tu2))+(2*c3*c2))/(2*(c3)^2);                    %sensivity for delivery lead time for omni-channel    
      
     
    
      TR(1,i)=(p(1,i)+s-pm)*(y(1,i)+inv(.5*(1+erf((ru(1,i)-mu)/(Sigma*(2^.5))))))-s*y(1,i)-(p(1,i)+s+h)*(int((Zi),-Sigma*sqrt(6),inv(.5*(1+erf((ru(1,i)-mu)/(Sigma*(2^.5)))))))  + (p(2,i)+s-pm)*(y(2,i)+inv(.5*(1+erf((ru(2,i)-mu)/(Sigma*(2^.5))))))-s*y(2,i)-(p(2,i)+s+h)*(int((Zi),-Sigma*sqrt(6),inv(.5*(1+erf((ru(2,i)-mu)/(Sigma*(2^.5)))))))  + (p(3,i)+s-pm)*(y(3,i)+inv(.5*(1+erf((ru(3,i)-mu)/(Sigma*(2^.5))))))-s*y(3,i)-(p(3,i)+s+h)*(int((Zi),-Sigma*sqrt(6),inv(.5*(1+erf((ru(3,i)-mu)/(Sigma*(2^.5))))))) - (phi+phii*r(1,i))*(r(1,i)+h) - (c0-c1*l(1,i))^2  - (c2-c3*l(2,i))^2 - k;
    
      Landa(3,i) = (y(3,i)+inv(.5*(1+erf((ru(3,i)-mu)/(Sigma*(2^.5)))))) + (p(1,i)-pm)*delta2 + (p(2,i)-pm)*delta3 - (p(3,i)-pm)*(beta+delta2+delta3)- int((Zi),-Sigma*sqrt(6),inv(.5*(1+erf((ru(3,i)-mu)/(Sigma*(2^.5)))))) ; 

        
      Results(i).Price = [p(1,i),p(2,i),p(3,i)];
      Results(i).Pricer = r(1,i);
      Results(i).Quantity = [Q(1,i), Q(2,i), Q(3,i)];
      Results(i).deliverytime = [l(1,i), l(2,i)];
      Results(i).TRevenue = TR(1,i);
      Results(i).Landa = [Landa(1,i),Landa(2,i),Landa(3,i)];
      
