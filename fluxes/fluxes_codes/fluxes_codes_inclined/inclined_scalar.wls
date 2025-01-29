#!/usr/bin/env wolframscript
(* ::Package:: *)

(*Exit[];*)
(*CODICE FLUSSI SCALARI*)
(*Questa versione tiene conto di tutti gli erratum relativi al paper di Hughes del 2000. 
La working precision dipende dalla precisione della funzione integranda.
In questo codice viene realizzata una griglia equispaziata in u e x. Dove x \[EGrave] definito 
nel paper di Hughes https://arxiv.org/pdf/2102.02713.pdf. La maggior parte delle formule vengono
da l\[IGrave].
*)
(*In questo codice, la griglia \[EGrave] fatta rispetto a u=1/sqrt[r-Subscript[r, lso]], vengono stampati i modi in l, m e in k
e vengono salvati sul file solo i flussi e non tutto. Per la variazione dei parametri orbitale, vedi vecchio codice*)

SetDirectory[NotebookDirectory[]];
SetSystemOptions["ParallelOptions"->"MKLThreadNumber"->55];
SetSystemOptions["ParallelOptions"->"ParallelThreadNumber"->55];
LaunchKernels[55];
<<SpinWeightedSpheroidalHarmonics`
<<KerrGeodesics`
<<Teukolsky`
ParallelNeeds["SpinWeightedSpheroidalHarmonics`"];
ParallelNeeds["KerrGeodesics`"];
ParallelNeeds["Teukolsky`"];


<<SpinWeightedSpheroidalHarmonics`
<<KerrGeodesics`
<<Teukolsky`


InputParameter=300;


imax=20; (*number of points of the array, -1*)
eps= 2/100;
jmax=6;
a0=0.95`100;
arrayu[i_,a_]:=Range[1/Sqrt[1/10*KerrGeoISSO[a,i]+10+eps],1/Sqrt[1/10*KerrGeoISSO[a,i]+eps]+10^-6,(1/Sqrt[1/10*KerrGeoISSO[a,i]+eps]-1/Sqrt[1/10*KerrGeoISSO[a,i]+10+eps])/imax];
arrayp[i_,a_]:=Table[{(1/arrayu[i,a][[x]])^2+9/10*KerrGeoISSO[a,i],(1/arrayu[i,a][[x]])^2+9/10*KerrGeoISSO[a,i]-KerrGeoISSO[a,i]},{x,1,Length[arrayu[i,a]]}];
xgrid=SetPrecision[Range[-0.9,0.9,2/jmax],InputParameter]; (*gli estremi dovrebbero essere -1 e 1*)


Rgrid={};
InputParameter=300;
Table[
Rarray=arrayp[xgrid[[j]],a0];
AppendTo[Rgrid,Rarray];,{j,1,jmax}] (*dovrebbe essere jmax+1*)


Headerk={{"l" ,"m","k","R0","R0-RLSO","E","L","K","x","\[CapitalDelta]E\[Infinity]","\[CapitalDelta]EH","\[CapitalDelta]L\[Infinity]","\[CapitalDelta]LH","\[CapitalDelta]E","\[CapitalDelta]L"}};
Headerl={{"l" ,"R0","R0-RLSO","E","L","K","x","\[CapitalDelta]E\[Infinity]","\[CapitalDelta]EH","\[CapitalDelta]L\[Infinity]","\[CapitalDelta]LH","\[CapitalDelta]E","\[CapitalDelta]L"}};
a0=0.95`100;
FILEk="output_"<>ToString[DecimalForm[a0,3]]<>"_kmodes.tsv";
FILEl="output_"<>ToString[DecimalForm[a0,3]]<>"_lmodes.tsv";

Export["../results/S_"<>FILEk,Headerk];
Export["../results/S_"<>FILEl,Headerl];


SCALARINTEGRATION[l_,m_,k_,omega_,CP_]:=Block[{eq,str,data,precision,\[Epsilon]=10^(-8),rH,YH,\[Lambda],r,Y\[Infinity],drt,pW,Z\[Infinity],ZH,zh,z\[Infinity],dE\[Infinity],dEH,rh,A\[Infinity],AH,Omega},
If[m>l,Print["Error: m cannot be larger than l"];
Return[Null],
If[CP!= 0,Omega=SetPrecision[omega,OmegaPrecision],
If[L<0,Omega=m \[Omega]pmin,Omega=m \[Omega]p]];
\[Lambda]=SpinWeightedSpheroidalEigenvalue[0,l,m,a Omega]-(a Omega)^2+2 m  (a Omega) ;
rH=(M + \[Sqrt](M^2 - a^2))+2 \[Epsilon];
rh=(M + \[Sqrt](M^2 - a^2));
YH[r_]=TeukolskyRadial[0,l,m,a,Omega]["In"][r]*Sqrt[r^2+a^2]/Sqrt[rh^2+a^2];
Y\[Infinity][r_]=TeukolskyRadial[0,l,m,a,Omega]["Up"][r]*Sqrt[r^2+a^2];
drt=(R0^2+a^2)/\[CapitalDelta]0;
pW=(YH[R0]Y\[Infinity]'[R0]-Y\[Infinity][R0]YH'[R0])/drt;
precision=Max[Precision[((\[Gamma]+a^2 e \[Mu]minus Cos[\[Theta]New]^2)/(\[Sqrt](\[Beta](\[Mu]plus-\[Mu]minus Cos[\[Theta]New]^2))))*(Conjugate[SpinWeightedSpheroidalHarmonicS[0,l,m,a Omega][ArcCos[Sqrt[\[Mu]minus]Cos [\[Theta]New]],0]]/dt)E^(I(Omega tFunction-m \[Phi]Function))//.\[Theta]New->\[Pi]/3]-5,40];
If[CP!=0,
(*Print["Integrand precision",precision];*)
SOURCE=-4\[Pi] ((mp \[CapitalDelta]0 ScalarCharge)/(TTheta (R0^2+a^2)^(3/2)))NIntegrate[((\[Gamma]+a^2 e \[Mu]minus Cos[\[Theta]New]^2)/(\[Sqrt](\[Beta](\[Mu]plus-\[Mu]minus Cos[\[Theta]New]^2))))*(Conjugate[SpinWeightedSpheroidalHarmonicS[0,l,m,a Omega][ArcCos[Sqrt[\[Mu]minus]Cos [\[Theta]New]],0]]/dt)E^(I(Omega tFunction-m \[Phi]Function)),{\[Theta]New,0,2\[Pi]},AccuracyGoal->8,WorkingPrecision->precision],
S0=SpinWeightedSpheroidalHarmonicS[0,l,m,a Omega,\[Pi]/2,0];SOURCE=-4\[Pi] mp ScalarCharge (\[CapitalDelta]0/(a^2+R0^2)^((3/2)))(Conjugate[S0]/dt0);];
zh=(Y\[Infinity][R0])*SOURCE;
z\[Infinity]=(YH[R0]) *SOURCE;
Z\[Infinity]=(1/pW) drt z\[Infinity] ;
ZH=(1/pW)drt zh;
dE\[Infinity]=(1/(16\[Pi]))Omega^2 Abs[Z\[Infinity]]^2;
dEH=(1/(16\[Pi]))Omega( Omega- (m a)/(2 M (M+ Sqrt[M^2-a^2]))) Abs[ZH]^2;
If[NumberQ[dEH],dEH=dEH,dEH=0];
If[NumberQ[dE\[Infinity]],dE\[Infinity]=dE\[Infinity],dE\[Infinity]=0];
(*Print["dE\[Infinity]",Precision["dE\[Infinity]"]];(*Test precision*)
PrintTemporary["Precision dEH= ",Precision[dEH]," Precision dE\[Infinity]= ",Precision[dE\[Infinity]]];(*Test precision*)*)
A\[Infinity]=m (1/(16\[Pi]))Omega Abs[Z\[Infinity]]^2;
AH=m (1/(16\[Pi]))( Omega- (m a)/(2 M (M+ Sqrt[M^2-a^2]))) Abs[ZH]^2;
If[NumberQ[AH],AH=AH,AH=0];
If[NumberQ[A\[Infinity]],A\[Infinity]=A\[Infinity],A\[Infinity]=0];
data={{l,m,k,R0,\[CapitalDelta]R,e,L,Q2,x,dE\[Infinity],dEH,A\[Infinity],AH,dE\[Infinity]+dEH,A\[Infinity]+AH}};
str=OpenAppend["S_"<>FILEk];
WriteString[str,ExportString[data,"TSV"]];
Close[str];
Return[{dE\[Infinity],dEH,dE\[Infinity]+dEH,A\[Infinity],AH,A\[Infinity]+AH}];
]
]
MYINTEGRATIONGRID[kstart_,kend_,l_,m_,\[Epsilon]_,CP_]:=Block[{k,F,Ftot,dAH,dA\[Infinity],AH,A\[Infinity],dEH,FH,dE\[Infinity],F\[Infinity],i,j,Fmax,Flux,AFlux,Atot,A,AHN,A\[Infinity]N,FHN,F\[Infinity]N},
Ftot=0;
Atot=0;
i=0;
j=0;
Fmax=0;
dAH=0;
dA\[Infinity]=0;
dEH=0;
dE\[Infinity]=0;
For[k=kstart,k<kend,k=k+2,
Clear[F,A];
If[i<3 && j<2,
{F\[Infinity],FH,F,A\[Infinity],AH,A}=SCALARINTEGRATION[l,m,k,\[Omega],CP];
If[Abs[F]>Fmax,Fmax=Abs[F]];
If[Abs[F]<\[Epsilon] Fmax || Abs[F]<10^-25,i=i+1,i=0];
If[k==0,
Flux=F;AFlux=A;AHN=AH;A\[Infinity]N=A\[Infinity];FHN=FH;F\[Infinity]N=F\[Infinity],
Flux=2*F;AFlux=2*A;AHN=2*AH;A\[Infinity]N=2*A\[Infinity];FHN=2*FH;F\[Infinity]N=2*F\[Infinity];
];
Ftot=Ftot+Flux;
Atot=Atot+AFlux;
dAH=dAH+AHN;
dA\[Infinity]=dA\[Infinity]+A\[Infinity]N;
dEH=dEH+FHN;
dE\[Infinity]=dE\[Infinity]+F\[Infinity]N;
PrintTemporary[{l,m,k,Flux,Ftot,\[Omega]}];,
Break[];
];
];
Return[{dE\[Infinity],dEH,Ftot,dA\[Infinity],dAH,Atot}];
]
lResults[l_,\[Epsilon]_,CP_]:=Block[{str,data,dE\[Infinity]l,dEHl,dA\[Infinity]l,dAHl,EnergyFlux,kend,lEnergy,kmin,kmax,m,lAngularMomentum,AngularMomentumFlux,dA\[Infinity],dAH,dE\[Infinity],dEH},
lEnergy=0;
lAngularMomentum=0;
{dE\[Infinity]l,dEHl,dA\[Infinity]l,dAHl}={0,0,0,0};
For[m=-l,m<=l,m++,
{EnergyFlux,AngularMomentumFlux}={0,0};
If[CP!=0,

If[EvenQ[l+m]==True,kmin=0;kmax=40,kmin=1;kmax=41];
{dE\[Infinity],dEH,EnergyFlux,dA\[Infinity],dAH,AngularMomentumFlux}=MYINTEGRATIONGRID[kmin,kmax,l,m,\[Epsilon],CP],

If[EvenQ[m+l]==True && m!=0,
{dE\[Infinity],dEH,EnergyFlux,dA\[Infinity],dAH,AngularMomentumFlux}=SCALARINTEGRATION[l,m,0,0,0]]];

(*PrintTemporary[{m,EnergyFlux}];*)
lEnergy=lEnergy+EnergyFlux;
lAngularMomentum=lAngularMomentum+AngularMomentumFlux;
{dE\[Infinity]l,dEHl,dA\[Infinity]l,dAHl}={dE\[Infinity]l,dEHl,dA\[Infinity]l,dAHl}+{dE\[Infinity],dEH,dA\[Infinity],dAH};
];
data={{l,R0,\[CapitalDelta]R,e,L,Q2,x,dE\[Infinity]l,dEHl,dA\[Infinity]l,dAHl,lEnergy,lAngularMomentum}};
str=OpenAppend["S_"<>FILEl];
WriteString[str,ExportString[data,"TSV"]];
Close[str];
Return[{dE\[Infinity]l,dEHl,lEnergy,dA\[Infinity]l,dAHl,lAngularMomentum}];
]
ScalarResults[lend_,\[Epsilon]_,CP_]:=Block[{l,lFlux,i,Fmax,LFlux,aFlux,AFlux,dE\[Infinity],dEH,dA\[Infinity],dAH,dE\[Infinity]l,dEHl,dA\[Infinity]l,dAHl},
lFlux=0;
LFlux=0;
AFlux=0;
aFlux=0;
dE\[Infinity]=0;
dEH=0;
dA\[Infinity]=0;
dAH=0;
i=0;
Fmax=0;
For[l=0,l<lend,l++,
(*PrintTemporary["Scalar mode, l=",l];*)
If[i<3,
{dE\[Infinity]l,dEHl,lFlux,dA\[Infinity]l,dAHl,aFlux}=lResults[l,\[Epsilon]/10,CP];
LFlux=LFlux+lFlux;
dE\[Infinity]=dE\[Infinity]+dE\[Infinity]l;
dEH=dEH+dEHl;
dA\[Infinity]=dA\[Infinity]+dA\[Infinity]l;
dAH=dAH+dAHl;
AFlux=AFlux+aFlux;
If[Abs[lFlux]>Fmax,Fmax=Abs[lFlux]];
If[Abs[lFlux]<\[Epsilon] Fmax || Abs[lFlux]<10^-20,i=i+1,i=0],
Break[];
];
];
Return[{dE\[Infinity]l,dEHl,LFlux,dA\[Infinity]l,dAHl,AFlux}];
]


\[CapitalSigma]=SetPrecision[r^2+a^2 Cos[\[Theta]]^2,InputParameter];
\[CapitalDelta]=SetPrecision[r^2+a^2-2M r,InputParameter];
rH=SetPrecision[M+Sqrt[M^2-a^2]+10^-16,InputParameter];
\[CapitalOmega]h=SetPrecision[a/(2 M rH),InputParameter];
\[Omega]=\[Omega]mk;


\[CapitalDelta]0=SetPrecision[\[CapitalDelta]//.r->R0,InputParameter];


\[Omega]pmin=-Sqrt[M]/(R0^(3/2)-a Sqrt[M]);
qmin=1-3 M/R0-2a Sqrt[M/R0^3];
e0min=(1-2 M/R0-a Sqrt[M/R0^3])/(Sqrt[qmin]);
Lmin=-(Sqrt[M]/Sqrt[(qmin/R0)])(a^2/R0^2+1+2a Sqrt[M/R0^3]);
\[Omega]p=Sqrt[M]/(R0^(3/2)+a Sqrt[M]);
q=1-3 M/R0+2a Sqrt[M/R0^3];
e0=(1-2 M/R0+a Sqrt[M/R0^3])/(Sqrt[q]);
L0=Sqrt[M R0]/Sqrt[q] (a^2/R0^2+1-2a Sqrt[M/R0^3]);
dt0=1/\[CapitalDelta]0 ((R0^2+a^2+2 (M a^2)/R0)e0-(2M a L0)/R0);
dtmin=1/\[CapitalDelta]0 ((R0^2+a^2+2 (M a^2)/R0)e0min-(2M a Lmin)/R0);


\[Beta]=SetPrecision[a^2 (1-e^2),InputParameter];
Q2=SetPrecision[Q+(L-a e)^2,InputParameter];

(* the following function was built with an old version of the BHPT, which was giving a positive angular momentum for x=-1: 

Lfunc[a_,R0_,inc_]:=Block[{LL},
LL=KerrGeoConstantsOfMotion[a,R0,0,inc][["\[ScriptCapitalL]"]];
If[inc*LL<0,LL=-LL,LL=LL];
Return[LL]];*)

Lfunc[a_,R0_,inc_]:=KerrGeoConstantsOfMotion[a,R0,0,inc][["\[ScriptCapitalL]"]]
efunc[a_,R0_,inc_]:=KerrGeoConstantsOfMotion[a,R0,0,inc][["\[ScriptCapitalE]"]]
Qfunc[a_,R0_,inc_]:=KerrGeoConstantsOfMotion[a,R0,0,inc][["\[ScriptCapitalQ]"]] 
(**)
\[CapitalRho]=SetPrecision[(e(r^2+a^2)-a L)^2-\[CapitalDelta](r^2+(L-a e)^2+Q),InputParameter];
D\[Rho]=D[\[CapitalRho],r];
DD\[Rho]=D[\[CapitalRho],{r,2}];

c\[Iota]=SetPrecision[L/Sqrt[L^2+Q],InputParameter]; (*nota che in questo codice x non corrisponde a questa variabile!!!*)
(*x=SetPrecision[L/Sqrt[Q2+2a e L-(a^2) (e^2) ],InputParameter];*) (*vecchia, nel nuovo codice non serve definirla*)
\[Mu]minus=SetPrecision[(L^2+Q+\[Beta]-Sqrt[(L^2+Q+\[Beta])^2-4\[Beta] Q])/(2\[Beta]),InputParameter];
(*\[Mu]minus=0;*)
\[Mu]plus=SetPrecision[(L^2+Q+\[Beta]+Sqrt[(L^2+Q+\[Beta])^2-4\[Beta] Q])/(2\[Beta]),InputParameter];
\[CapitalSigma]0=SetPrecision[R0^2+a^2  \[Mu]minus Cos[\[Theta]New]^2,InputParameter] ;
\[Gamma]=SetPrecision[e((R0^2+a^2)^2/\[CapitalDelta]0-a^2)+a L (1-(R0^2+a^2)/\[CapitalDelta]0),InputParameter];
\[Delta]=SetPrecision[a e((R0^2+a^2)/\[CapitalDelta]0-1)-a^2 L/\[CapitalDelta]0,InputParameter];
dt=SetPrecision[1/\[CapitalDelta]0 ((R0^2+a^2+(2M R0 a^2 (1-\[Mu]minus Cos[\[Theta]New]^2))/\[CapitalSigma]0)e-(2M R0 a L)/\[CapitalSigma]0 ),InputParameter];
\[Mu]ratio=SetPrecision[\[Mu]minus/\[Mu]plus,InputParameter];
tFunction=SetPrecision[\[Gamma]/Sqrt[\[Beta] \[Mu]plus] (EllipticK[\[Mu]ratio]-EllipticF[\[Pi]/2-\[Theta]New,\[Mu]ratio])+a^2 e Sqrt[\[Mu]plus/\[Beta]](EllipticE[\[Pi]/2-\[Theta]New,\[Mu]ratio]-EllipticE[\[Pi]/2,\[Mu]ratio]+EllipticK[\[Mu]ratio]-EllipticF[\[Pi]/2-\[Theta]New,\[Mu]ratio]),InputParameter];
TTheta=SetPrecision[4*tFunction//.\[Theta]New->\[Pi]/2,InputParameter];
\[Phi]Function=SetPrecision[1/Sqrt[\[Beta] \[Mu]plus] (L(EllipticPi[\[Mu]minus,\[Pi]/2,\[Mu]ratio]-EllipticPi[\[Mu]minus,\[Pi]/2-\[Theta]New,\[Mu]ratio])+\[Delta](EllipticK[\[Mu]ratio]-EllipticF[\[Pi]/2-\[Theta]New,\[Mu]ratio])),InputParameter];
\[CapitalOmega]\[Theta]=SetPrecision[2 \[Pi]/TTheta,InputParameter];
\[CapitalOmega]\[Phi]=SetPrecision[(4*\[Phi]Function//.\[Theta]New->\[Pi]/2)/TTheta,InputParameter];
\[Omega]mk=SetPrecision[m \[CapitalOmega]\[Phi]+k \[CapitalOmega]\[Theta],InputParameter];


OmegaPrecision=300;
a=a0;
eps=10^-2;
M=1;
mp=1;
ScalarCharge=1;
If[Im[Lmin]==0,Lmin=Lmin,Lmin=0];
EpsScalar=0.001;
EpsTensorial=0.001;


Header={{"j" ,"i","R0","R0-RLSO","E","L","K","x","\[CapitalDelta]E\[Infinity]","\[CapitalDelta]EH","\[CapitalDelta]L\[Infinity]","\[CapitalDelta]LH","\[CapitalDelta]E","\[CapitalDelta]L"}};
FILE="../results/output_S_"<>ToString[DecimalForm[a,3]]<>"_fast.tsv";
Export[FILE,Header];
lMax=20; (*20*)
kMax=20;(*20*)


ParallelTable[
i=Mod[ij,imax+1];
j=Floor[ij/(imax+1)];
R0=Rgrid[[j+1,i+1,1]];
\[CapitalDelta]R=Rgrid[[j+1,i+1,2]];
x=xgrid[[j+1]];
Print[i,j,R0,x];
e=efunc[a,R0,x];
L=Lfunc[a,R0,x];
Q=Qfunc[a,R0,x];
{dE\[Infinity],dEH,\[CapitalDelta]E,dA\[Infinity],dAH,\[CapitalDelta]L}=ScalarResults[lMax,EpsScalar,-1];
Print["{S\[CapitalDelta]E,S\[CapitalDelta]L}//Precision= ",{\[CapitalDelta]E//Precision,\[CapitalDelta]L//Precision}];
data={{j,i,\[CapitalDelta]R,R0,e,L,Q2,x,dE\[Infinity],dEH,dA\[Infinity],dAH,\[CapitalDelta]E,\[CapitalDelta]L}};
Print[data];
str=OpenAppend[FILE];
WriteString[str,ExportString[data,"TSV"]];
Close[str];
,{ij,0,(jmax)*(imax+1)-1}](*se sopra \[EGrave] jmax+1 metti jmax+1 anche qui, altrimenti va bene jmax*)



