clc;

NAMES = { 'FresnelCS', ...
          'GeneralizedFresnelCS', ...
          'buildClothoid', ...
          'evalClothoid', ...
          'pointsOnClothoid' } ;

disp('---------------------------------------------------------');
for k=1:length(NAMES)
  N=NAMES{k} ;
  fprintf(1,'Compiling: %s\n',N) ;

  CMD = ['mex -output ../G1fitting/',N,' -largeArrayDims mex_',N,'.cc Clothoid.cc'] ;
  if isunix
    if ismac
      CMD = [CMD, ' -lstdc++ CXXFLAGS="\$CXXFLAGS -Wall -O2 -g0"'] ;
    else
      CMD = [CMD, ' -lstdc++ CXXFLAGS="\$CXXFLAGS -Wall -O2 -g0"'] ;
    end
  elseif ispc
  end
  disp(CMD);
  eval(CMD);
end
disp('----------------------- DONE ----------------------------');
