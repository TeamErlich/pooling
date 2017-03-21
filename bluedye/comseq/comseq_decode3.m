% 3 blue wells

% calculated diff between measurement (ex = 630nm) and ref (em = 470nm)

y_pre = [0.0087	-0.0032	-0.0032	-0.003	-0.0025	0.0088	-0.0026	-0.0021	-0.0018	-0.0018	0.0118	-0.0028
-0.002	0.0138	-0.0034	-0.0032	-0.0033	-0.0027	-0.0031	-0.0031	-0.0043	-0.002	-0.004	-0.0016
-0.0028	-0.0023	-0.003	-0.0019	0.0326	-0.002	-0.0028	-0.0019	-0.0024	-0.0013	-0.004	-0.0017
-0.0024	-0.0025	-0.0037	0.0193	-0.0032	-0.003	-0.0014	-0.0037	0.0072	-0.0018	-0.0033	0.0177
-0.0021	-0.0018	0.0081	-0.003	-0.002	-0.0032	-0.0013	0.007	-0.0024	0.0101	-0.0011	-0.0021
-0.0018	-0.004	0.0148	-0.0012	-0.0028	-0.0028	-0.0018	-0.002	-0.001	-0.0043	-0.003	-0.0045
-0.0021	0.0199	-0.0021	-0.0018	-0.002	-0.0028	-0.0024	0.0125	-0.0026	0.0152	-0.0023	-0.0017
-0.003	-0.0034	-0.0024	-0.0023	-0.0022	-0.003	0.0161	-0.0029	-0.0029	-0.0037	0.0089	-0.0018];

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


