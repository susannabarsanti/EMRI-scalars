(* ::Package:: *)

SetSystemOptions["ParallelOptions"->"MKLThreadNumber"->21];
SetSystemOptions["ParallelOptions"->"ParallelThreadNumber"->21];
LaunchKernels[21];


Exit[];


<<SpinWeightedSpheroidalHarmonics`
<<KerrGeodesics`
<<Teukolsky`
ParallelNeeds["SpinWeightedSpheroidalHarmonics`"];
ParallelNeeds["KerrGeodesics`"];
ParallelNeeds["Teukolsky`"];


(*master equation*) 


\[CapitalDelta][r_] = r^2+a^2- 2 M r ;
\[CapitalOmega][r_,a_]=M^(1/2)/(r^(3/2)+ a M^(1/2)) ; 
pmk[r_,a_]= m \[CapitalOmega][r,a] - (m a)/(2 M (M+ Sqrt[M^2-a^2]))  ;


ddRH[r_]=Solve[\[CapitalDelta][r]^2 D[1/\[CapitalDelta][r] D[RH[r],r],r]-V[r]RH[r]==0,RH''[r]][[1,1,2]]//Simplify;
ddR\[Infinity][r_]=Solve[\[CapitalDelta][r]^2 D[1/\[CapitalDelta][r] D[R\[Infinity][r],r],r]-V[r]R\[Infinity][r]==0,R\[Infinity]''[r]][[1,1,2]]//Simplify;


KK[r_] = (r^2 + a^2)\[Omega]- m a;
\[Lambda]= \[ScriptCapitalE] ;
V[r_] = - (((KK[r])^2+4 I (r-M) KK[r] )/\[CapitalDelta][r])+ 8 I \[Omega] r  + \[Lambda]; 
cc = {-12 I \[Omega] M + \[Lambda] ( \[Lambda]+2) -12 a \[Omega] ( a \[Omega] - m ), 8 I a (3 a \[Omega] - \[Lambda] (a \[Omega] -m)), -24 I a M (a \[Omega] -m) + 12 a^2  (1 -2 (a \[Omega] -m)^2), 24 I a^3 (a \[Omega] -m) -24 M a^2 , 12 a^4};
\[Eta][r_]= Sum[cc[[i+1]]/r^i,{i,0,4}]; 
F[r_]=D[\[Eta][r],r]/\[Eta][r] \[CapitalDelta][r]/(r^2+a^2);
G[r_]= -((2 (r-M))/(r^2+a^2))+(r \[CapitalDelta][r])/(r^2+a^2)^2;
\[Alpha][r_]= -I KK[r] B[r]/\[CapitalDelta][r]^2+ 3 I D[KK[r],r ]+\[Lambda] + (6 \[CapitalDelta][r])/r^2;
B[r_]= 2 \[CapitalDelta][r] (- I KK[r]+r -M - (2 \[CapitalDelta][r])/r);
U1[r_]= V[r]+\[CapitalDelta][r]^2/B[r] (D[2 \[Alpha][r] + D[B[r],r]/\[CapitalDelta][r],r]-D[\[Eta][r],r]/\[Eta][r] (\[Alpha][r]+D[B[r],r]/\[CapitalDelta][r]));

U[r]=  (\[CapitalDelta][r] U1[r])/(r^2+a^2 )^2+ (G[r])^2+ (\[CapitalDelta][r] D[G[r],r])/(r^2+a^2)-F[r] G[r]; 
drs=(r^2+a^2)/\[CapitalDelta][r]; 
master=1/drs^2 Derivative[2][u][r]+1/drs D[1/drs,r]Derivative[1][u][r]-F[r] Derivative[1][u][r]/drs - U[r] u[r];

master=Collect[master,{Derivative[2][u][r],Derivative[1][u][r],u[r]},Simplify];


(*source term*) 


en[a_,p_,e_]:= KerrGeoEnergy[a,p,e,1];
L[a_,p_,e_]:=KerrGeoAngularMomentum[a,p,e,1];

\[Phi][a_,p_,e_,\[Chi]_?NumberQ]:= NIntegrate[V\[Phi][a,p,e,chi]/(J[a,p,e,chi] (Vr[a,p,e,chi])^(1/2)),{chi, 0, \[Chi]},WorkingPrecision->50];
t[a_,p_,e_,\[Chi]_?NumberQ]:= NIntegrate[Vt[a,p,e,chi]/(J[a,p,e,chi] (Vr[a,p,e,chi])^(1/2)),{chi, 0, \[Chi]},WorkingPrecision->50];
x[a_,p_,e_]:=L[a,p,e] -a en[a,p,e];

Vr[a_,p_,e_,\[Chi]_]:= x[a,p,e]^2+a^2+2 a x[a,p,e] en[a,p,e]-(2 M x[a,p,e]^2)/p (3 + e Cos[\[Chi]]) ;
V\[Phi][a_,p_,e_,\[Chi]_]:=x[a,p,e] + a en[a,p,e] - (2 M x[a,p,e])/p (1+ e Cos[\[Chi]]);

Vt[a_,p_,e_,\[Chi]_]:= a^2 en[a,p,e]-(2 a M x[a,p,e] )/p (1+ e Cos[\[Chi]]) + (en[a,p,e] p^2)/(1+ e Cos[\[Chi]])^2;

J[a_,p_,e_,\[Chi]_] := 1 - (2 M)/p (1+ e Cos[\[Chi]])+ (a/p)^2 (1+ e Cos[\[Chi]])^2;
erre[e_,p_,\[Chi]_]:= p/(1 + e Cos[\[Chi]]);


source[p_,e_,a_,l_,m_,k_]:=Block[{r,\[Theta],mp=1,M=1,\[Chi],\[Omega],S,dS,d2S,\[CapitalOmega]r,\[CapitalOmega]\[Phi],\[Omega]mk,Cnn ,Cmn,Cmm,Amn0,Amm0,Amn1,Amm1,Amm2,Ann0,vr,vt,v\[Phi],j,ur,y,rep,fr},
\[CapitalSigma] = (r^2+a^2 (Cos[\[Theta]])^2)//.{\[Theta]-> \[Pi]/2};

dt =  en[a,p,e]/\[CapitalSigma] ((r^2+a^2)^2/\[CapitalDelta][r]-a^2 (Sin[\[Theta]])^2) +( a L[a,p,e])/\[CapitalSigma] (1-(r^2+a^2)/\[CapitalDelta][r]) //.{\[Theta]-> \[Pi]/2} ; 

fr = KerrGeoFrequencies[a,p,e,1];

\[CapitalOmega]r= fr[[1]]; 
\[CapitalOmega]\[Phi]= fr[[3]]; 
\[Omega]mk = m \[CapitalOmega]\[Phi] + k \[CapitalOmega]r; 
\[Omega]=\[Omega]mk;

vr = Vr[a,p,e,\[Chi]];
v\[Phi] = V\[Phi][a,p,e,\[Chi]];
vt = Vt[a,p,e,\[Chi]];
j = J[a,p,e,\[Chi]];

ur = erre[e,p,\[Chi]]^-1;

Cnn = j/(4 p^4 vt) (p^2 en[a,p,e] -a x[a,p,e] (1+e Cos[\[Chi]])^2+e p Sin[\[Chi]] (vr)^(1/2) )^2; 
Cmn = (I x[a,p,e] j)/(2 Sqrt[2] p^3 vt) (1+e Cos[\[Chi]]) (p^2 en[a,p,e] -a x[a,p,e] (1+e Cos[\[Chi]])^2+e p Sin[\[Chi]] (vr)^(1/2) ); 
Cmm= -(( x[a,p,e]^2 j)/(2 p^2 vt)) (1+e Cos[\[Chi]])^2; 

Amn0 = 2/Sqrt[ \[Pi]] Cmn /(ur (1 - 2 M ur + a^2  ur^2)^2) (2 a^2 ur^3+ (I a (a \[Omega]-m)-4 M) ur^2+ 2 ur + I \[Omega]) (dS+ (a \[Omega] - m) S );
Amm0 = 1/Sqrt[2 \[Pi]] (Cmm S)/(ur^2 (1 - 2 M ur + a^2  ur^2)^2) (-2 I a^3 (a \[Omega] - m)ur^5+a (a \[Omega] - m)(6 I M + a (a \[Omega] - m ))ur^4-4 I a (a \[Omega] - m)ur^3+2 \[Omega] (I M + a (a \[Omega] - m))ur^2-2 I \[Omega] ur + \[Omega]^2);
Amn1 = 2/Sqrt[\[Pi]] Cmn/(ur(1 - 2 M ur + a^2  ur^2)) (dS+ (a \[Omega] - m) S );
Amm1 = -Sqrt[(2/\[Pi])] (Cmm S)/(ur^2 (1 - 2 M ur + a^2  ur^2)) (a^2 ur^3+ (I a (a \[Omega] -m)-2 M )ur^2+ ur + I \[Omega] );
Amm2 = -Sqrt[(1/(2 \[Pi]))] (Cmm S)/ur^2;
Ann0 = -Sqrt[(2/\[Pi])] Cnn /(1 - 2 M ur + a^2  ur^2)^2 (- 2 I a (dS + (a \[Omega] -m) S )ur + d2S + 2 (a \[Omega] -m) dS + ( (a \[Omega] -m)^2-2 )S );

S=Sqrt[2 \[Pi]]SpinWeightedSpheroidalHarmonicS[-2,l,m,a \[Omega],\[Pi]/2,0];
dS = D[Sqrt[2 \[Pi]] SpinWeightedSpheroidalHarmonicS[-2,l,m,a \[Omega],y,0],y]//.y-> \[Pi]/2;
d2S = D[Sqrt[2 \[Pi]] SpinWeightedSpheroidalHarmonicS[-2,l,m,a \[Omega],y,0],{y,2}]//.y-> \[Pi]/2;

Return[{Ann0+Amn0+Amm0,Amn1+Amm1,Amm2,((Ann0+Amn0+Amm0)//.{Sin[\[Chi]]->rep})//.{rep->-Sin[\[Chi]]},((Amn1+Amm1)//.{Sin[\[Chi]]->rep})//.{rep->-Sin[\[Chi]]},((Amm2)//.{Sin[\[Chi]]->rep})//.{rep->-Sin[\[Chi]]}}]
];


(*integrator*) 


Int[latus_,ecc_,\[ScriptL]_,mu_,kappa_, Aa_, rinf_]:=Block[{ fr,rpiu, l=\[ScriptL],dEd\[Sigma],Aout,m=mu,a = Aa,p = latus,e=ecc, k=kappa,r,M=1,mp=1,r1,\[Epsilon]=10^-16,eq,qSh, XHI,dXHI, X\[Infinity]I,dX\[Infinity]I,X,W,Y\[Infinity],YH, RH,RH0, RH1,RH2,\[ScriptCapitalE],Ain, Bin ,r2,R\[Infinity],Dtran,\[Theta]=\[Pi]/2,IH,I\[Infinity],S,dS,d2S,uHI,duHI, u\[Infinity]I, du\[Infinity]I,d0,b0,eqShp,eqShm,eqSh,eqS\[Infinity],eqS\[Infinity]m,eqS\[Infinity]p,f,\[CapitalDelta],df,d2f,\[Chi],ip,\[CapitalOmega]r,\[CapitalOmega]\[Phi],\[Omega]mk, \[Omega],pi,rs, y,vr,vt,v\[Phi],j,in,out,din,dout,Aout\[Infinity],Aouth,dEd\[Sigma]\[Infinity],dEd\[Sigma]h,dlm,alpha,eps,ci,im,dLd\[Sigma]\[Infinity],dLd\[Sigma]h},

fr = KerrGeoFrequencies[a,p,e,1];
\[CapitalOmega]r= fr[[1]]; 
\[CapitalOmega]\[Phi]= fr[[3]]; 
\[Omega]mk = m \[CapitalOmega]\[Phi] + k \[CapitalOmega]r; 
\[Omega]=\[Omega]mk;

pi= \[Omega]- (m a)/(2 M (M+ Sqrt[M^2-a^2]));

rpiu= M+Sqrt[M^2-a^2];
r1=M+Sqrt[M^2-a^2]+\[Epsilon];
r2=rinf/Abs[\[Omega]];

\[CapitalDelta][r_]=r^2-2M r+a^2;

eq=master;

\[ScriptCapitalE]= SpinWeightedSpheroidalEigenvalue[-2,l,m,a \[Omega]]; 

S=SpinWeightedSpheroidalHarmonicS[-2,l,m,a \[Omega],\[Pi]/2,0];
dS = D[SpinWeightedSpheroidalHarmonicS[-2,l,m,a \[Omega],y,0],y]//.y-> \[Pi]/2;
d2S = D[SpinWeightedSpheroidalHarmonicS[-2,l,m,a \[Omega],y,0],{y,2}]//.y-> \[Pi]/2;

rs[r] = r+(2 M (M + Sqrt[M^2 - a^2]))/(2  Sqrt[M^2 - a^2]) Log[(r-M - Sqrt[M^2 - a^2] )/(2 M)]- (2 M (M - Sqrt[M^2 - a^2]))/(2  Sqrt[M^2 - a^2]) Log[(r-M + Sqrt[M^2 - a^2] )/(2 M)];


dlm= Sqrt[2 M rpiu]((8-24 I M \[Omega] -16 M^2 \[Omega]^2)rpiu^2+(12 I a m - 16 M + 16 a m M \[Omega] + 24 I M^2 \[Omega] )rpiu - 4 a^2  m^2 - 12 I a m M + 8 M^2); 
alpha = (256 (2 M rpiu)^5 pi (pi^2+4 eps^2)(pi^2+16 eps^2) \[Omega]^3)/ci ; 
eps = Sqrt[M^2-a^2]/(4 M rpiu); 
ci= ((\[Lambda]+2)^2+4 a m \[Omega] - 4 a^2  \[Omega]^2)(\[Lambda]^2+36 a m \[Omega] -36 a^2  \[Omega]^2)+(2\[Lambda]+3)(96 a^2 \[Omega]^2-48 a m \[Omega])+144 \[Omega]^2 (M^2-a^2); 

Dtran=-((4\[Omega]^2)/cc[[1]]);

(*een if precision of Teuk is Machine Precision, when you evaluate it in a r`prec it becomes prec (a bit less)*)
RH[r_]=TeukolskyRadial[-2,l,m,a,\[Omega]]["In"][r]/dlm;
R\[Infinity][r_]=TeukolskyRadial[-2,l,m,a,\[Omega]]["Up"][r]*Dtran;

W[r_]=(RH[r]R\[Infinity]'[r]-R\[Infinity][r]RH'[r])/\[CapitalDelta][r];

Bin=W[p]/(2I \[Omega]  Dtran);

ip =  f source[p,e,a,l,m,k][[1]] -df source[p,e,a,l,m,k][[2]]+d2f source[p,e,a,l,m,k][[3]];
im =  f source[p,e,a,l,m,k][[4]] -df source[p,e,a,l,m,k][[5]]+d2f source[p,e,a,l,m,k][[6]]; 

vr = Vr[a,p,e,\[Chi]];
v\[Phi] = V\[Phi][a,p,e,\[Chi]];
vt = Vt[a,p,e,\[Chi]];
j = J[a,p,e,\[Chi]];


eqShp = ((mp \[CapitalOmega]r  vt/(j (vr)^(1/2)) ip Exp[+I (\[Omega]mk t[a,p,e,\[Chi]]-m \[Phi][a,p,e,\[Chi]])])//.{f->RH[r], df->RH'[r] , d2f->ddRH[r]})//.r-> erre[e,p,\[Chi]];
eqShm = ((mp \[CapitalOmega]r  vt/(j (vr)^(1/2)) im Exp[-I (\[Omega]mk t[a,p,e,\[Chi]]-m \[Phi][a,p,e,\[Chi]])])//.{f->RH[r], df->RH'[r] , d2f->ddRH[r]})//.r-> erre[e,p,\[Chi]]; 
eqSh[\[Chi]_] = eqShp+eqShm;
 
  
eqS\[Infinity]p = ((mp \[CapitalOmega]r  vt/(j (vr)^(1/2)) ip Exp[+I (\[Omega]mk t[a,p,e,\[Chi]]-m \[Phi][a,p,e,\[Chi]])])//.{f->R\[Infinity][r], df->R\[Infinity]'[r] , d2f->ddR\[Infinity][r]})//.r-> erre[e,p,\[Chi]]; 
eqS\[Infinity]m = ((mp \[CapitalOmega]r  vt/(j (vr)^(1/2)) im Exp[-I (\[Omega]mk t[a,p,e,\[Chi]]-m \[Phi][a,p,e,\[Chi]])])//.{f->R\[Infinity][r], df->R\[Infinity]'[r] , d2f->ddR\[Infinity][r]})//.r-> erre[e,p,\[Chi]]; 
eqS\[Infinity][\[Chi]_] = eqS\[Infinity]p+eqS\[Infinity]m;

Print[Precision[eqSh[0]-2]];
Aout = 1/(2 I \[Omega] Bin ) NIntegrate[eqSh[\[Chi]],{\[Chi],0,\[Pi]},MaxRecursion->14,AccuracyGoal->4,WorkingPrecision->Precision[eqSh[0]]-4];
dEd\[Sigma]\[Infinity]= Abs[Aout]^2/(4 \[Pi]  \[Omega]^2);
Print[Aout//Precision];


Aout\[Infinity] = -( cc[[1]]/(8 I \[Omega]^3 Bin  dlm))NIntegrate[eqS\[Infinity][\[Chi]],{\[Chi],0,\[Pi]},MaxRecursion->25,AccuracyGoal->4,WorkingPrecision->Precision[eqS\[Infinity][0]]-4];
Print[Aout\[Infinity]//Precision];


dEd\[Sigma]h=  alpha Abs[Aout\[Infinity]]^2/(4 \[Pi] (\[Omega]^2) );

dLd\[Sigma]\[Infinity] = (m/\[Omega])*dEd\[Sigma]\[Infinity]; 
dLd\[Sigma]h = (m/\[Omega])*dEd\[Sigma]h; 

Return[{dEd\[Sigma]h,dEd\[Sigma]\[Infinity],dLd\[Sigma]h,dLd\[Sigma]\[Infinity]}]
];


(*BH toolkit for reference*) 


IntTool[R_,\[ScriptL]_,mu_,A_]:=Block[{l=\[ScriptL],dEd\[Sigma],dEd\[Sigma]\[Infinity],dEd\[Sigma]\[Infinity]2,dEd\[Sigma]h2,dEd\[Sigma]h,  m=mu, a=A,r0=R,orbit, mode,s,n,k,Rin,\[Omega],M=1,Rup},

(* KerrGeoOrbit[a,r0,e,x]
a: the black hole spin
p: the semi-latus rectum
e: the eccentricity
x: cos\[Theta]inc = the orbital inclination.*)
orbit=KerrGeoOrbit[a,r0,0,1];

s=-2;n=0;k=0;
mode=TeukolskyPointParticleMode[s,l,m,n,k,orbit];

dEd\[Sigma]h = mode["Fluxes"]["Energy"]["\[ScriptCapitalH]"];
dEd\[Sigma]\[Infinity] =   mode["Fluxes"]["Energy"]["\[ScriptCapitalI]"]; 
dEd\[Sigma]h2 = dEd\[Sigma]h*2;
dEd\[Sigma]\[Infinity]2 = dEd\[Sigma]\[Infinity]*2;
dEd\[Sigma] = dEd\[Sigma]h2 + dEd\[Sigma]\[Infinity]2;
\[Omega] = (m M^(1/2))/(r0^(3/2)+ a M^(1/2)); 
Rin=TeukolskyRadial[s,l,m,a,\[Omega]]["In"];
Rup=TeukolskyRadial[s,l,m,a,\[Omega]]["Up"];
Return[{dEd\[Sigma]h,dEd\[Sigma]\[Infinity]}]
]


(*check with zero eccentricity*) 


IntTool[7`60,2,2,0.9`60]
Int[7`180,0,2,2,0,0.9`180,1000]


IntTool[7`60,2,-2,0.9`60]
Int[7`180,0,2,-2,0,0.9`180,1000]


(*Export grav_e_p_lm.dat. In each file, there are all the k for the lm of that file. 
Each k is computed until the fractional change to the accumulated sum is smaller then 10^-4 
for three consecutives value of k, when I stop. Otherwise I stop to (input) nmax.

The total flux is obtained by summing over all l, all m(even negative), n from 0 to nmax, 
where the non-zero-n fluxes are multiplied by a factor 2

In principle my idea was not only to print the .dat file but also to sum over k, returning 
the sum as the output of the function, but there is some mistake in the sum in the code so 
I never used the output.
*) 


tablesum[lmin_,lman_,kmax_,p_,e_,a_]:=Block[{M,i,j,m,flux,k,ics,array,tab,kfin,ind},
For[j=lmin,j<=lman,j++,
For[m=-j,m<=j,m++,
i=0;
ics={10^-35,10^-35,10^-35,10^-35};
array=ConstantArray[0,kmax+1];
If[m!=0,
For[k=0,k<=kmax,k++,
Print["(l,m,n) = ","(",j,",",m,",",k,")"];
Print["p = ","(",SetPrecision[p,7],")"];
If[k==0,
flux=Int[p,e,j,m,k,a,1000];,
flux=2*Int[p,e,j,m,k,a,1000];
];
array[[k+1]]=flux;
kfin=k;
If[flux[[2]]/ics[[2]] > 10^-4, ics = ics + flux,
ics=ics+flux;
i=i+1;
If[i==3, 
Break[];
]
]
];
tab=Table[Flatten[{ind-1,array[[ind]]}],{ind,1,kfin+1,1}];
Export["grav_e"<>ToString[IntegerPart[e]]<>StringDrop[ToString[FractionalPart[SetPrecision[e,2]]],2]<>"_p"<>ToString[IntegerPart[p]]<>StringDrop[ToString[FractionalPart[SetPrecision[p,6]]],2]<>"_lm"<>ToString[j]<>ToString[m]<>".dat",tab];
, If[j<4,For[k=1,k<=kmax,k++,
Print["(l,m,n) = ","(",j,",",m,",",k,")"];
Print["p = ","(",SetPrecision[p,7],")"];
flux=2*Int[p,e,j,m,k,a,1000];
array[[k]]=flux;
kfin=k;
If[flux[[2]]/ics[[2]] > 10^-4, ics = ics + flux,
ics=ics+flux;
i=i+1;
If[i==3, 
Break[];
]
]
]];
tab=Table[Flatten[{ind,array[[ind]]}],{ind,1,kfin,1}];
Export["grav_e"<>ToString[IntegerPart[e]]<>StringDrop[ToString[FractionalPart[SetPrecision[e,2]]],2]<>"_p"<>ToString[IntegerPart[p]]<>StringDrop[ToString[FractionalPart[SetPrecision[p,6]]],2]<>"_lm"<>ToString[j]<>ToString[m]<>".dat",tab];]; 
]];
Return[ics];
];


(*array p gives the grid in p, uniformily spaced in the quantity arrayup*)


arrayup[ee_,a_]:= Range[1/Sqrt[0.1`150*KerrGeoSeparatrix[a,ee,1]+10.02`150],1/Sqrt[0.1`150*KerrGeoSeparatrix[a,ee,1]+0.02`150],(1/Sqrt[0.1`150*KerrGeoSeparatrix[a,ee,1]+0.02`150]-1/Sqrt[0.1`150*KerrGeoSeparatrix[a,ee,1]+10.02`150])/40];


arrayp[ee_,a_]:= Table[(1/arrayup[ee,a][[i]])^2+0.9`150*KerrGeoSeparatrix[a,ee,1],{i,1,Length[arrayup[ee,a]]}];


(*final function to compute for different eccentricity values- see https://arxiv.org/abs/2102.02713*)


tablep[emin_,emax_,a_]:=Block[{M,i,pp,ee,m,flusso,up},
flusso = ConstantArray[0,IntegerPart[(emax-emin)*10.01+1]];
For[ee=emin,ee<=emax, ee=ee+0.1`150,flusso[[IntegerPart[(ee-emin)*10.01+1]]]=Partition[Flatten[ParallelTable[{arrayp[ee,a][[i]],ee,tablesum[2,2,1,arrayp[ee,a][[i]],ee,a]},{i,1,2}]],6];
Export["grav_e"<>ToString[IntegerPart[ee]]<>StringDrop[ToString[FractionalPart[SetPrecision[ee,2]]],2]<>"_a"<>ToString[IntegerPart[a*10.1]]<>".dat",Partition[Flatten[flusso[[IntegerPart[(ee-emin)*10.01+1]]]],6]];
];
];


Exit[];
