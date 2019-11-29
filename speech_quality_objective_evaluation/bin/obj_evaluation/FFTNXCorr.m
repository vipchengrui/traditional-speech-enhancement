function Y= FFTNXCorr( ref_VAD, startr, nr, deg_VAD, startd, nd)
% this function has other simple implementations, current implementation is
% consistent with the C version

% one way to do this (in time domain) =====
x1= ref_VAD( startr: startr+ nr- 1);
x2= deg_VAD( startd: startd+ nd- 1);
x1= fliplr( x1);
Y= conv( x2, x1);
% done =====

% % the other way to do this (in freq domain)===
% Nx= 2^ (ceil( log2( max( nr, nd))));
% x1= zeros( 1, 2* Nx);
% x2= zeros( 1, 2* Nx);
% x1( 1: nr)= fliplr( ref_VAD( startr: startr+ nr- 1));
% x2( 1: nd)= deg_VAD( startd: startd+ nd- 1);
% 
% if (nr== 491)
%     fid= fopen( 'mat_debug.txt', 'wt');
%     fprintf( fid, '%f\n', x1);
%     fclose( fid);
% end
% 
% x1_fft= fft( x1, 2* Nx);
% x2_fft= fft( x2, 2* Nx);
% 
% tmp1= ifft( x1_fft.* x2_fft, 2* Nx);
% 
% Ny= nr+ nd- 1;
% Y= tmp1( 1: Ny);
% % done ===========





















