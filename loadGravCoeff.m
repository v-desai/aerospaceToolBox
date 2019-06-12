function [C,S] = loadGravCoeff(model,gravOrder)
if model == 'EGM96'
    % To Download: https://cddis.nasa.gov/926/egm96/getit.html
    % Directly: ftp://cddis.gsfc.nasa.gov/pub/egm96/general_info/egm96_to360.ascii
    gravCoeff = load('egm96_to360_ascii.txt');
end
maxTerms = max([1,min([gravOrder,360])]);
C = zeros(maxTerms,maxTerms);
S = zeros(maxTerms,maxTerms);
if maxTerms >= 2
    index = 1;
    for i = 3:maxTerms+1
        for j = 1:i
            C(i,j) = gravCoeff(index,3);
            S(i,j) = gravCoeff(index,4);
            index = index+1;
        end  %  for j = 1:i
    end  %  for i = 3:maxTerms+1,
end  %  if ( maxTerms >= 2 )
