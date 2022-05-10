(* ::Package:: *)

SetSystemOptions["ParallelOptions"->"MKLThreadNumber"->41];
SetSystemOptions["ParallelOptions"->"ParallelThreadNumber"->41];
LaunchKernels[41];


Exit[];


<<SpinWeightedSpheroidalHarmonics`
<<KerrGeodesics`
<<Teukolsky`
ParallelNeeds["SpinWeightedSpheroidalHarmonics`"];
ParallelNeeds["KerrGeodesics`"];
ParallelNeeds["Teukolsky`"];


\[CapitalDelta]  = r^2+a^2- 2 M r ;
master =   \[CapitalDelta] D[\[CapitalDelta]  D[R[r],r],r] + V[r] R[r];
V[r] = -4 a m M r \[Omega]+r (2 M \[Lambda]-r \[Lambda]+r^3 \[Omega]^2)+a^2 (m^2-\[Lambda]+r (2 M+r) \[Omega]^2) ; 
master = Simplify[master//.{R[r]-> u[r]/Sqrt[a^2+r^2],R'[r]-> D[u[r]/Sqrt[a^2+r^2],r],  R''[r]-> D[u[r]/Sqrt[a^2+r^2],{r,2}]}];
drs = (a^2+r^2)/\[CapitalDelta];
master = Simplify[master//.{u''[r]-> u''[\[ScriptR]], u'[r]-> u'[\[ScriptR]] }];
master=Simplify[master//.{u''[\[ScriptR]]-> u''[r] drs^2+u'[r] D[drs,r],u'[\[ScriptR]]-> u'[r] drs}];
master=Simplify[master * (a^2+r^2)^(5/2) / (a^2+r^2)^4 ]//.{u''[r]-> u''[\[ScriptR]] };
master=Collect[master//.{u''[\[ScriptR]]->1/drs^2 u''[r]+1/drs D[1/drs,r]u'[r]},{u''[r],u'[r],u[r]},FullSimplify];


\[CapitalSigma] = (r^2+a^2 (Cos[\[Theta]])^2)//.{\[Theta]-> \[Pi]/2};
pi= \[Omega]- (m a)/(2 M (M+ Sqrt[M^2-a^2]))  ;
en[a_,p_,e_]:= KerrGeoEnergy[a,p,e,1];
L[a_,p_,e_]:=KerrGeoAngularMomentum[a,p,e,1];

dt =  en[a,p,e]/\[CapitalSigma] ((r^2+a^2)^2/\[CapitalDelta]-a^2 (Sin[\[Theta]])^2) +( a L[a,p,e])/\[CapitalSigma] (1-(r^2+a^2)/\[CapitalDelta]) //.{\[Theta]-> \[Pi]/2} ; 
source[\[Chi]_]:= - 4 \[Pi] di mp \[CapitalOmega]r/(2 \[Pi] (a^2+r^2)^(1/2))  Conjugate[S]/dt Vt[\[Chi]]/(J[\[Chi]] (Vr[\[Chi]])^(1/2))  (Exp[I (\[Omega]mk t[\[Chi]]-m \[Phi][\[Chi]])]+ Exp[-I (\[Omega]mk t[\[Chi]]-m \[Phi][\[Chi]])]) ;
x=L[a,p,e] -a en[a,p,e] ; 
Vr[\[Chi]_]:= x^2+a^2+2 a x en[a,p,e]-(2 M x^2)/p (3 + e Cos[\[Chi]]) ;
V\[Phi][\[Chi]_]:= x + a en[a,p,e] - (2 M x)/p (1+ e Cos[\[Chi]]);
Vt[\[Chi]_]:= a^2 en[a,p,e]-(2 a M x )/p (1+ e Cos[\[Chi]]) + (en[a,p,e] p^2)/(1+ e Cos[\[Chi]])^2;
J[\[Chi]_] := 1 - (2 M)/p (1+ e Cos[\[Chi]])+ (a/p)^2 (1+ e Cos[\[Chi]])^2;
\[Phi][\[Chi]_?NumberQ]:= NIntegrate[V\[Phi][chi]/(J[chi] (Vr[chi])^(1/2)),{chi, 0, \[Chi]},WorkingPrecision->150];
t[\[Chi]_?NumberQ]:= NIntegrate[Vt[chi]/(J[chi] (Vr[chi])^(1/2)),{chi, 0, \[Chi]},WorkingPrecision->150];
erre[\[Chi]_]:= p/(1 + e Cos[\[Chi]]);

\[CapitalOmega]r= KerrGeoFrequencies[a,p,e,1][[1]]; 
\[CapitalOmega]\[Phi]= KerrGeoFrequencies[a,p,e,1][[3]]; 
\[Omega]mk = m \[CapitalOmega]\[Phi] + k \[CapitalOmega]r; 


SpinWeightedSpheroidalHarmonicS[0,0,0,0,\[Pi]/2,0]


IntScal[latus_,ecc_,\[ScriptL]_,mu_,kappa_,dii_,A_,rinf_]:=Block[{l=\[ScriptL],rp,dEd\[Sigma],dEd\[Sigma]\[Infinity],dEd\[Sigma]h,dEd\[Sigma]\[Infinity]2,dEd\[Sigma]h2,  Aout,Aouth, m=mu, k=kappa,\[Omega]=\[Omega]mk//.{p-> latus, e-> ecc},r,M=1,mp=1,p=latus,e=ecc,r1,\[Epsilon]=10^-3,eq,eqS,eqSh,uHI,duHI,d0,b0,u\[Infinity]I,du\[Infinity]I,\[Psi],S,\[Lambda], W,drt,b,d,Y\[Infinity],YH, di=dii, a=A, normS,r2,\[Chi],dLd\[Sigma]\[Infinity],dLd\[Sigma]h},

PrintTemporary["prec t = ",Precision[t[\[Pi]/2]]];
eq=master;

r1=(M + Sqrt[M^2 - a^2])+2 \[Epsilon];
rp=(M + Sqrt[M^2 - a^2]);
r2=rinf/Abs[\[Omega]];

S=SpinWeightedSpheroidalHarmonicS[0,l,m,a \[Omega],\[Pi]/2,0];

YH[r_]=TeukolskyRadial[0,l,m,a,\[Omega]]["In"][r]*Sqrt[r^2+a^2]/Sqrt[rp^2+a^2];
Y\[Infinity][r_]=TeukolskyRadial[0,l,m,a,\[Omega]]["Up"][r]*Sqrt[r^2+a^2];

drt=(r^2+a^2)/\[CapitalDelta];

W=(YH[r]Y\[Infinity]'[r]-Y\[Infinity][r]YH'[r])/drt//.r->p;

eqS=(YH[r])//.r-> erre[\[Chi]];
eqS= eqS * source[\[Chi]]//.r-> erre[\[Chi]];

Print["p = ",Precision[p]];
Aout=W^-1 NIntegrate[eqS,{\[Chi],0,\[Pi]},MaxRecursion->12,AccuracyGoal->5,WorkingPrecision->Precision[eqS//.\[Chi]->0]-5];
Print["prec Aout = ",Precision[Aout]];
dEd\[Sigma]\[Infinity]=1/(4\[Pi]) \[Omega]^2 Abs[Aout]^2;

eqSh=(Y\[Infinity][r])//.r-> erre[\[Chi]];
eqSh= eqSh * source[\[Chi]]//.r-> erre[\[Chi]];

Aouth=W^-1 NIntegrate[eqSh,{\[Chi],0,\[Pi]},MaxRecursion->12,AccuracyGoal->5,WorkingPrecision->Precision[eqSh//.\[Chi]->0]-5];
Print["prec Aouth = ",Precision[Aouth]];
dEd\[Sigma]h=1/(4\[Pi]) \[Omega] pi  Abs[Aouth]^2;
dLd\[Sigma]\[Infinity] = (m/\[Omega])*dEd\[Sigma]\[Infinity]; 
dLd\[Sigma]h = (m/\[Omega])*dEd\[Sigma]h;

dEd\[Sigma]=dEd\[Sigma]h+ dEd\[Sigma]\[Infinity];

Return[{dEd\[Sigma]h,dEd\[Sigma]\[Infinity],dLd\[Sigma]h,dLd\[Sigma]\[Infinity]}]
];


tablesum[lmin_,lman_,kmax_,p_,e_,d_,a_]:=Block[{M,i,j,m,flux,k,ics,array,tab,kfin,ind},
For[j=lmin,j<=lman,j++,
For[m=-j,m<=j,m++,
If[EvenQ[j+Abs[m]],
i=0;
ics={10^-35,10^-35,10^-35,10^-35};
array=ConstantArray[0,kmax+1];
If[m!=0,
For[k=0,k<=kmax,k++,
Print["(l,m,n) = ","(",j,",",m,",",k,")"];
If[k==0,
flux=IntScal[p,e,j,m,k,d,a,1000];,
flux=2*IntScal[p,e,j,m,k,d,a,1000];
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
Export["scal_e"<>ToString[IntegerPart[e]]<>StringDrop[ToString[FractionalPart[SetPrecision[e,2]]],2]<>"_p"<>ToString[IntegerPart[p]]<>StringDrop[ToString[FractionalPart[SetPrecision[p,6]]],2]<>"_lm"<>ToString[j]<>ToString[m]<>".dat",tab];
, For[k=1,k<=kmax,k++,
Print["(l,m,n) = ","(",j,",",m,",",k,")"];
If[k==0,
flux=IntScal[p,e,j,m,k,d,a,1000];,
flux=2*IntScal[p,e,j,m,k,d,a,1000];
];
array[[k]]=flux;
kfin=k;
If[flux[[2]]/ics[[2]] > 10^-4, ics = ics + flux,
ics=ics+flux;
i=i+1;
If[i==3, 
Break[];
]
]
];
tab=Table[Flatten[{ind,array[[ind]]}],{ind,1,kfin,1}];
Export["scal_e"<>ToString[IntegerPart[e]]<>StringDrop[ToString[FractionalPart[SetPrecision[e,2]]],2]<>"_e"<>ToString[IntegerPart[e]]<>StringDrop[ToString[FractionalPart[SetPrecision[p,6]]],2]<>"_lm"<>ToString[j]<>ToString[m]<>".dat",tab];]; 
]]];
Return[ics];
];


arrayup[ee_,a_]:= Range[1/Sqrt[0.1`180*KerrGeoSeparatrix[a,ee,1]+10.02`180],1/Sqrt[0.1`180*KerrGeoSeparatrix[a,ee,1]+0.02`180],(1/Sqrt[0.1`180*KerrGeoSeparatrix[a,ee,1]+0.02`180]-1/Sqrt[0.1`180*KerrGeoSeparatrix[a,ee,1]+10.02`180])/40];


arrayp[ee_,a_]:= Table[(1/arrayup[ee,a][[i]])^2+0.9`180*KerrGeoSeparatrix[a,ee,1],{i,1,Length[arrayup[ee,a]]}];


tablep[emin_,emax_,dd_,a_]:=Block[{M,i,pp,ee,m,flusso,up},
flusso = ConstantArray[0,IntegerPart[(emax-emin)*10.01+1]];
For[ee=emin,ee<=emax, ee=ee+0.1`180,flusso[[IntegerPart[ee*10.01]]]=Partition[Flatten[ParallelTable[{arrayp[ee,a][[i]],ee,tablesum[1,5,90,arrayp[ee,a][[i]],ee,dd,a]},{i,1,Length[arrayp[ee,a]]}]],6];
Export["scal_e"<>ToString[IntegerPart[ee]]<>StringDrop[ToString[FractionalPart[SetPrecision[ee,2]]],2]<>"_d01a"<>ToString[IntegerPart[a*10.1]]<>".dat",Partition[Flatten[flusso[[IntegerPart[ee*10.01]]]],6]];
];
Export["scal_d01a"<>ToString[IntegerPart[a*10.1]]<>".dat",Partition[Flatten[flusso],6]];
];


tablep[0.2`60,0.2`60, 1`60,0.2`60];


Exit[];
