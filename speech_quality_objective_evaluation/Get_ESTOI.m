clear all

cleanFile = 'cleanSpeech\TIMIT_201_TEST.wav';
[x,FS0] = audioread(cleanFile);

x = x(1:4820096);

enhancedFile_1 = 'Weiner+UnwrappedPhase_201\out_Weiner_UP_201_TEST-White-5db.wav';
[y1,FS1] = audioread(enhancedFile_1);
d1 = estoi(x, y1, 8000);

enhancedFile_2 = 'Weiner+UnwrappedPhase_201\out_Weiner_UP_201_TEST-White+0db.wav';
[y2,FS2] = audioread(enhancedFile_2);
d2 = estoi(x, y2, 8000);

enhancedFile_3 = 'Weiner+UnwrappedPhase_201\out_Weiner_UP_201_TEST-White+5db.wav';
[y3,FS3] = audioread(enhancedFile_3);
d3 = estoi(x, y3, 8000);

enhancedFile_4 = 'Weiner+UnwrappedPhase_201\out_Weiner_UP_201_TEST-White+10db.wav';
[y4,FS4] = audioread(enhancedFile_4);
d4 = estoi(x, y4, 8000);