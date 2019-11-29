clear
cleanFile = 'DataSet/TIMIT_10_TEST-White+5db.wav';
enhancedFile_1 = 'Enhanced/out_TIMIT_10_TEST_IRMDNN_0308_10.wav';
[snr_mean, segsnr_mean]= comp_snr(cleanFile, enhancedFile_1);
