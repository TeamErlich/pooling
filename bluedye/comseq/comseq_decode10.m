% 10 blue specimens

% calculated diff between measurement (ex = 630nm) and ref (em = 470nm)

y_pre = [0.0176	-0.0004	-0.0006	-0.0014	-0.0014	0.0349	0.0169	0.0113	0.0233	0.0008	0.0008	0.0018
-0.0038	0.0416	0.0171	-0.0016	-0.0015	0.0201	-0.0032	-0.0023	0.0006	0.016	0.0399	0.0198
-0.0016	0.0181	0.0192	-0.0028	-0.0016	0.0122	-0.0014	-0.001	0.0005	0.0005	0.0002	0.0019
0.0137	0.0116	0.0223	0.0089	0.0224	0.0092	-0.0006	0.0318	0.035	0.0082	0.0163	0.0083
-0.0001	0.0185	0.0374	0.0132	-0.0016	-0.001	0.018	0.014	-0.001	0.0015	0.0089	0.0023
-0.002	-0.0022	0.0184	-0.0014	0.0168	0.0165	0.0216	0.0234	0.0215	-0.0012	0.0015	0.0008
0.0154	0.0125	-0.001	0.0154	-0.0014	-0.0017	0.0215	-0.0015	0.0279	0.016	-0.0004	0.0169
0.0194	-0.0022	-0.0009	0.0248	0.0414	-0.0019	-0.0001	0.0183	0.0205	-0.0007	0.0552	0.0015];

imagesc(y_pre)
colorbar

% reshape matrix to 96x1 vector
y=reshape(y_pre,96,1);

% pooling matrix
A = load('RS_matrix.tab');

imagesc(A)

% normalize pooling matrix
A0 = normalizeRows(A,2)/2;

tau = 0.005*max(abs(A0'*y));

fractionalOutput = applyGPSR(y, A0, tau);
bar(fractionalOutput)

% xlim([length(fractionalOutput)+1]);

discreteOutput = modifySolution(fractionalOutput, A0, y);
bar(discreteOutput)

x=discreteOutput;

%find pools containing positive specimens
C.loc=find(x==1);

% well/plate info
C.Plate_Well=mod(C.loc,96);
C.Plate_No=ceil(C.loc/96);


