clear all;

close all ;

alpha = 1
gamma = 0:0.005:1 ;

beta = gamma.*alpha ;

theta = 0:0.01:1 ;
theta = theta .* (pi/4) ;

for m = 1 : length(beta)
    for n = 1 : length(theta)
        A1 = alpha ;
        A2 = sqrt(beta(m)) ;
        
        a = [ -A1 + A2/cos(pi/4)*cos(pi/4-theta(n)) ; +A1 - A2/cos(pi/4)*sin(pi/4-theta(n)) ] ;
        b = [ +A1 - A2/cos(pi/4)*cos(pi/4+theta(n)) ; +A1 - A2/cos(pi/4)*sin(pi/4+theta(n)) ] ;
        c = [ +A1 - A2/cos(pi/4)*cos(pi/4-theta(n)) ; +A1 + A2/cos(pi/4)*sin(pi/4-theta(n)) ] ;
        
        D(m,n) = min( norm(a-b), norm(a-c) ) ;
    
    end


end

D

figure(1)
mesh( theta./(pi/4)*45, gamma, D )

figure(2);
plot(  theta./(pi/4)*45, D( 41, : ), theta./(pi/4)*45, D( 61, : ), theta./(pi/4)*45, D( 81, : ), theta./(pi/4)*45, D( 101, : ), theta./(pi/4)*45, D( 121, : ), theta./(pi/4)*45, D( 141, : )   )

grid

legend( 'P_1/P_2= 0.2' , 'P_1/P_2= 0.3' , 'P_1/P_2= 0.4' , 'P_1/P_2= 0.5' , 'P_1/P_2= 0.6' , 'P_1/P_2= 0.7' ) ;