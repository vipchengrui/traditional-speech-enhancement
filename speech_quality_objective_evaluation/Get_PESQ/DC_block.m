function mod_data= DC_block( data, Nsamples)

global Downsample DATAPADDING_MSECS SEARCHBUFFER

ofs= SEARCHBUFFER* Downsample;
mod_data= data;

%compute dc component, it is a little weird
facc= sum( data( ofs+ 1: Nsamples- ofs))/ Nsamples; 
mod_data( ofs+ 1: Nsamples- ofs)= data( ofs+ 1: Nsamples- ofs)- facc;

mod_data( ofs+ 1: ofs+ Downsample)= mod_data( ofs+ 1: ofs+ Downsample).* ...
    ( 0.5+ (0: Downsample- 1))/ Downsample;

mod_data( Nsamples- ofs: -1: Nsamples- ofs-Downsample+ 1)= ...
    mod_data( Nsamples- ofs: -1: Nsamples- ofs-Downsample+ 1).* ...
    ( 0.5+ (0: Downsample- 1))/ Downsample;


     
    





