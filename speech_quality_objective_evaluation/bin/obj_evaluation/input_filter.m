function [mod_ref_data, mod_deg_data]= input_filter( ref_data, ref_Nsamples, ...
    deg_data, deg_Nsamples)

mod_ref_data= DC_block( ref_data, ref_Nsamples);
mod_deg_data= DC_block( deg_data, deg_Nsamples);

mod_ref_data= apply_filters( mod_ref_data, ref_Nsamples);
mod_deg_data= apply_filters( mod_deg_data, deg_Nsamples);

