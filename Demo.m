clear all
addpath('Datasets');
%% 1. Load Datesets
load Fdataset_ms
Wrr1 = drug_ChemS;
Wrr2 = drug_AtcS;
Wrr3 = drug_SideS;
Wrr4 = drug_DDIS;
Wrr5 = drug_TargetS;
Wrr = [Wrr1, Wrr2, Wrr3, Wrr4, Wrr5];
Wdd1 = disease_PhS;
Wdd2 = disease_DoS;
Wdd = [Wdd1, Wdd2];
Wdr = didr;
Wrd = Wdr';
[dn, dr] = size(Wdr);
min_mn = min(dn, dr);
%% 2. MSBMF algorithm
lambda1 = 0.1;
lambda2 = 0.01;
lambda3 = lambda2;
k = floor(min_mn * 0.7);
maxiter = 300;
tol1 = 2*1e-3;
tol2 = 1*1e-4;
[U, V, iter] = MSBMF(Wdr, Wdd, Wrr, lambda1, lambda2, lambda3, k, tol1, tol2, maxiter);
M_recovery = U * V';
