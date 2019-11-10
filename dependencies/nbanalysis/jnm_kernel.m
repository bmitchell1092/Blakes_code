function [ jnmKrnl ] = jnm_kernel( kType, kWidth )
%jnm_kernel.m
%Jake Westerberg
%Vanderbilt University
%August 4, 2016
%v1.0.0
%   Adapted from code written by Wolf Zinke. Supplies kernel for 
%   convolution of spike train.

if exist( 'kType', 'var' ) == false || isempty( kType ) == true
    
    kType = 'gauss';
    
end

if exist( 'kWidth', 'var' ) == false || isempty( kWidth ) == true
    
    kWidth = 20;
    
end

switch lower( kType )

    case 'gauss'
        
        Half_BW = round( 4 * kWidth );
        x = -Half_BW : Half_BW;
        jnmKrnl = ( 1 / ( sqrt( 2 * pi ) * kWidth ) ) * ...
            exp( -1 * ( ( x.^2 ) / ( 2 * kWidth^2 ) ) );
        
    case 'psp'
        
        if length( kWidth ) == 2
          
            tg = kWidth( 2 );
            kWidth( 2 ) = [];
            
        else
            
            tg = 1;
            
        end 
      
        Half_BW = ceil( kWidth * 8 );
        x = 0 : Half_BW;
        jnmKrnl = [ zeros( 1, Half_BW ), ...
            ( 1 - ( exp( -( x ./ tg ) ) ) ) .* ( exp( -( x ./ kWidth) ) ) ];
      
end

jnmKrnl = jnmKrnl - jnmKrnl( end );
jnmKrnl( jnmKrnl < 0 ) = 0;
jnmKrnl = jnmKrnl / sum( jnmKrnl );

end

