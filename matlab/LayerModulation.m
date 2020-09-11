%% Hierarchical Modulation
close all ; clear all ;

global ER N1 N2 Ns

ER = [ 3.25 3.5 3.75 4 4.25 4.5 4.75 5.00 ] ;

TYPE = 'QAM' ;

N1 = 4 ;
N2 = 4 ;
Ns = N1*N2 ;

switch lower( TYPE )
    case 'am'
        LayerModulation_AM
    otherwise
        LayerModulation_QAM
end