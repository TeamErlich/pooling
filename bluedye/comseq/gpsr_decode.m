% y = absorbance reading RS_pooling_10specimens

y_pre = [0.0612	0.0422	0.0416	0.0415	0.041	0.0794	0.0611	0.0542	0.0672	0.0443	0.0436	0.0446
0.0409	0.0873	0.0625	0.0423	0.0421	0.0655	0.0424	0.0419	0.044	0.0622	0.087	0.0639
0.0415	0.0621	0.0636	0.0404	0.0418	0.0562	0.0423	0.0423	0.0439	0.0442	0.0445	0.0451
0.0603	0.0561	0.0679	0.0534	0.0682	0.0551	0.0439	0.0784	0.0835	0.0537	0.0631	0.0521
0.0448	0.0643	0.0825	0.0568	0.0416	0.0422	0.0625	0.058	0.0428	0.0454	0.0537	0.0452
0.0438	0.0412	0.0633	0.0432	0.0626	0.063	0.0686	0.069	0.0689	0.0445	0.0473	0.0434
0.059	0.0564	0.0414	0.0598	0.0417	0.0417	0.0659	0.0418	0.074	0.0605	0.0437	0.0613
0.065	0.0403	0.0418	0.0686	0.0863	0.0414	0.0426	0.0627	0.0649	0.0442	0.1032	0.0432];

% reshape matrix to vector (y_pre, rows, cols)
y=reshape(y_pre,96,1);

A = load('RS_matrix.tab');

imagesc(A)

% plate reader output
imagesc(y_pre)
colorbar

A0 = normalizeRows(A,2)/2;

tau = 0.005*max(abs(A0'*y));

fractionalOutput = applyGPSR(y, A0, tau);
bar(fractionalOutput)
xlim([0 length(fractionalOutput)+1]);

discreteOutput = modifySolution(fractionalOutput, A0, y);
bar(discreteOutput);

x=discreteOutput;

%find pools
C.loc=find(x==1);

%which well in which plate each specimen is in, saved to Dinas_wells.txt
C.Plate_Well = mod(C.loc,96);
C.Plate_No   = ceil(C.loc/96);