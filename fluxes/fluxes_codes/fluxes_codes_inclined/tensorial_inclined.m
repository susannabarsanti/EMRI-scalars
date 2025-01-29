(* ::Package:: *)

(*QUESTO CODICE E' UNA COPIA DEL CODICE Fast_Ev_Hughes_Tensorial.wls. E' STATO INVIATO A SUSANNA PER RUNNARLO.*)
SetDirectory[NotebookDirectory[]]
SetSystemOptions["ParallelOptions"->"MKLThreadNumber"->55];
SetSystemOptions["ParallelOptions"->"ParallelThreadNumber"->55];
LaunchKernels[55];
<<SpinWeightedSpheroidalHarmonics`
<<KerrGeodesics`
<<Teukolsky`
ParallelNeeds["SpinWeightedSpheroidalHarmonics`"];
ParallelNeeds["KerrGeodesics`"];
ParallelNeeds["Teukolsky`"];



InputParameter=140;


imax=20; (*number of points of the array, -1*)
eps= 2/100;
jmax=6;
a0=0.95`140;
arrayu[i_,a_]:=Range[1/Sqrt[1/10*KerrGeoISSO[a,i]+10+eps],1/Sqrt[1/10*KerrGeoISSO[a,i]+eps]+10^-6,(1/Sqrt[1/10*KerrGeoISSO[a,i]+eps]-1/Sqrt[1/10*KerrGeoISSO[a,i]+10+eps])/imax];
arrayp[i_,a_]:=Table[{(1/arrayu[i,a][[x]])^2+9/10*KerrGeoISSO[a,i],(1/arrayu[i,a][[x]])^2+9/10*KerrGeoISSO[a,i]-KerrGeoISSO[a,i]},{x,1,Length[arrayu[i,a]]}];
xgrid=SetPrecision[Range[-91/100,91/100,2/jmax],InputParameter];(*la griglia dovrebbe essere da -1 a 1, ho messo questi valori perch\[EGrave] fare il test set*)


Rgrid={};
InputParameter=140;
Table[
Rarray=arrayp[xgrid[[j]],a0];
AppendTo[Rgrid,Rarray];,{j,1,jmax}](*ho tolto il + 1*)


Headerk={{"l" ,"m","k","R0","R0-RLSO","E","L","K","x","\[CapitalDelta]E\[Infinity]","\[CapitalDelta]EH","\[CapitalDelta]L\[Infinity]","\[CapitalDelta]LH","\[CapitalDelta]E","\[CapitalDelta]L","\!\(\*SubscriptBox[\(Z\), \(H\)]\)","\!\(\*SubscriptBox[\(Z\), \(\[Infinity]\)]\)","\[Omega]"}};
Headerl={{"l" ,"R0","R0-RLSO","E","L","K","x","\[CapitalDelta]E\[Infinity]","\[CapitalDelta]EH","\[CapitalDelta]L\[Infinity]","\[CapitalDelta]LH","\[CapitalDelta]E","\[CapitalDelta]L"}};
a0=0.95`140;

FILEk="output_"<>ToString[DecimalForm[a0,3]]<>"_kmodes.tsv";
FILEl="output_"<>ToString[DecimalForm[a0,3]]<>"_lmodes.tsv";
Export["../results/GR_"<>FILEk,Headerk];
Export["../results/GR_"<>FILEl,Headerl];


(*controlla condizione su m, \[EGrave] m+k che deve essere diverso da zero no m solo!!*)
TENSORIALINTEGRATE[l_,m_,k_,omega_]:=Block[{Rin,str,data,Rup,\[Lambda],ZH,Z\[Infinity],dEH,dE\[Infinity],A\[Infinity],AH,Omega,precision,argomento1,argomento2,argomento3, argomento4, armonica},
If[m==0 && k==0,Return[{0,0,0,0,0,0}]];
Omega=SetPrecision[omega,OmegaPrecision];

(*Rin[r_]:=TeukolskyRadial[s,l,m,a,Omega,WorkingPrecision->(InputParameter-5)]["In"][r]/dlm;
Rup[r_]:=TeukolskyRadial[s,l,m,a,Omega,WorkingPrecision->(InputParameter-5)]["Up"][r]*Dinf;*)

Rin[r_]:=TeukolskyRadial[s,l,m,a,Omega,WorkingPrecision->135]["In"][r]/dlm;
Rup[r_]:=TeukolskyRadial[s,l,m,a,Omega,WorkingPrecision->135]["Up"][r]*Dinf;


\[Lambda]=SpinWeightedSpheroidalEigenvalue[-2,l,m,a Omega];
(*precision=Precision[(\[Pi]/(I Omega TTheta Bin)/.r->R0)((\[Gamma]+a^2 e \[Mu]minus Cos[\[Theta]New]^2)/Sqrt[\[Beta](\[Mu]plus-\[Mu]minus Cos[\[Theta]New]^2)] E^(I(Omega tFunction-m \[Phi]Function)) (Iplus/.R->Rin/.S->SpinWeightedSpheroidalHarmonicS[s,l,m,a Omega]//.\[Theta]->ArcCos[Sqrt[\[Mu]minus]Cos[\[Theta]New]]/.\[Phi]->0/.r->R0)/.\[Theta]New->0)]-5;*)
precision=20;

armonica= SpinWeightedSpheroidalHarmonicS[s,l,m,a Omega];
argomento1= (\[Gamma]+a^2 e \[Mu]minus Cos[\[Theta]New]^2)/Sqrt[\[Beta](\[Mu]plus-\[Mu]minus Cos[\[Theta]New]^2)] E^(I(Omega tFunction-m \[Phi]Function)) (Iplus/.R->Rin/.S->armonica//.\[Theta]->ArcCos[Sqrt[\[Mu]minus]Cos[\[Theta]New]]/.\[Phi]->0/.r->R0);
argomento2= (\[Gamma]+a^2 e \[Mu]minus Cos[\[Theta]New]^2)/Sqrt[\[Beta](\[Mu]plus-\[Mu]minus Cos[\[Theta]New]^2)] E^(-I(Omega tFunction-m \[Phi]Function)) (Iminus/.R->Rin/.S->armonica//.\[Theta]->ArcCos[Sqrt[\[Mu]minus]Cos[\[Theta]New]]/.\[Phi]->0/.r->R0);
argomento3= (\[Gamma]+a^2 e \[Mu]minus Cos[\[Theta]New]^2)/Sqrt[\[Beta](\[Mu]plus-\[Mu]minus Cos[\[Theta]New]^2)] E^(I(Omega tFunction-m \[Phi]Function)) (Iplus/.R-> Rup/.S->armonica//.\[Theta]->ArcCos[Sqrt[\[Mu]minus]Cos[\[Theta]New]]/.\[Phi]->0/.r->R0);
argomento4= (\[Gamma]+a^2 e \[Mu]minus Cos[\[Theta]New]^2)/Sqrt[\[Beta](\[Mu]plus-\[Mu]minus Cos[\[Theta]New]^2)] E^(-I(Omega tFunction-m \[Phi]Function)) (Iminus/.R->Rup/.S->armonica//.\[Theta]->ArcCos[Sqrt[\[Mu]minus]Cos[\[Theta]New]]/.\[Phi]->0/.r->R0);

(*Print["prec2 = ",Precision[(\[Pi]/(I Omega TTheta Bin)/.r->R0)(argomento1/.\[Theta]New->0)]-5];

Print["prec2 = ",Precision[(argomento1/.\[Theta]New->0)]-5];*)

(*ZH=(\[Pi]/(I Omega TTheta Bin)/.r->R0)(NIntegrate[argomento1 + argomento2,{\[Theta]New,0,\[Pi]},WorkingPrecision->precision,AccuracyGoal->8]+NIntegrate[(\[Gamma]+a^2 e \[Mu]minus Cos[\[Theta]New]^2)/Sqrt[\[Beta](\[Mu]plus-\[Mu]minus Cos[\[Theta]New]^2)] E^(-I(Omega tFunction-m \[Phi]Function)) (Iminus/.R->Rin/.S->SpinWeightedSpheroidalHarmonicS[s,l,m,a Omega]//.\[Theta]->ArcCos[Sqrt[\[Mu]minus]Cos[\[Theta]New]]/.\[Phi]->0/.r->R0),{\[Theta]New,0,\[Pi]},WorkingPrecision->precision,AccuracyGoal->8]);
Z\[Infinity]=((-\[Pi] c0)/(4 I (Omega^3) dlm TTheta Bin )/.r->R0)(NIntegrate[(\[Gamma]+a^2 e \[Mu]minus Cos[\[Theta]New]^2)/Sqrt[\[Beta](\[Mu]plus-\[Mu]minus Cos[\[Theta]New]^2)] E^(I(Omega tFunction-m \[Phi]Function)) (Iplus/.R-> Rup/.S->SpinWeightedSpheroidalHarmonicS[s,l,m,a Omega]//.\[Theta]->ArcCos[Sqrt[\[Mu]minus]Cos[\[Theta]New]]/.\[Phi]->0/.r->R0),{\[Theta]New,0,\[Pi]},WorkingPrecision->precision,AccuracyGoal->8]+NIntegrate[(\[Gamma]+a^2 e \[Mu]minus Cos[\[Theta]New]^2)/Sqrt[\[Beta](\[Mu]plus-\[Mu]minus Cos[\[Theta]New]^2)] E^(-I(Omega tFunction-m \[Phi]Function)) (Iminus/.R->Rup/.S->SpinWeightedSpheroidalHarmonicS[s,l,m,a Omega]//.\[Theta]->ArcCos[Sqrt[\[Mu]minus]Cos[\[Theta]New]]/.\[Phi]->0/.r->R0),{\[Theta]New,0,\[Pi]},WorkingPrecision->precision,AccuracyGoal->8]);*)

ZH=(\[Pi]/(I Omega TTheta Bin)/.r->R0)(NIntegrate[argomento1 + argomento2,{\[Theta]New,0,\[Pi]},WorkingPrecision->precision,AccuracyGoal->8]);
Z\[Infinity]=((-\[Pi] c0)/(4 I (Omega^3) dlm TTheta Bin )/.r->R0)(NIntegrate[argomento3 + argomento4,{\[Theta]New,0,\[Pi]},WorkingPrecision->precision,AccuracyGoal->8]);


dE\[Infinity]=Abs[ZH]^2/(4\[Pi] Omega^2);
dEH=\[Alpha]lm Abs[Z\[Infinity]]^2/(4\[Pi] Omega^2);
If[NumberQ[dEH],dEH=dEH,dEH=0];
If[NumberQ[dE\[Infinity]],dE\[Infinity]=dE\[Infinity],dE\[Infinity]=0];
A\[Infinity]=m Abs[ZH]^2/(4\[Pi] Omega^3);
AH=m \[Alpha]lm Abs[Z\[Infinity]]^2/(4\[Pi] Omega^3);
(*Print["Precision dEH=",Precision[dEH]];
PrintTemporary["Precision dEH= ",Precision[dEH]," Precision dE\[Infinity]= ",Precision[dE\[Infinity]]];
Print[{dE\[Infinity],dEH,dE\[Infinity]+dEH}];*)
If[NumberQ[AH],AH=AH,AH=0];
If[NumberQ[A\[Infinity]],A\[Infinity]=A\[Infinity],A\[Infinity]=0];
data={{l,m,k,R0,\[CapitalDelta]R,e,L,Q2,x,dE\[Infinity],dEH,A\[Infinity],AH,dE\[Infinity]+dEH,A\[Infinity]+AH,ZH,Z\[Infinity],Omega}};
str=OpenAppend["GR_"<>FILEk];
WriteString[str,ExportString[data,"TSV"]];
Close[str];
Return[{dE\[Infinity],dEH,dE\[Infinity]+dEH,A\[Infinity],AH,A\[Infinity]+AH}];
]


kINTEGRATION[kend_,\[ScriptL]_,\[ScriptM]_,\[CapitalEpsilon]_]:=Block[{Flux,str,data,Fmax,FH,F\[Infinity],dEH,dE\[Infinity],dAH,dA\[Infinity],F,k,i,l=\[ScriptL],m=\[ScriptM],Flist,Fk,Ak,A,AH,A\[Infinity],AFlux,AHN,A\[Infinity]N,FHN,F\[Infinity]N},
Flux=0;
AFlux=0;
Fmax=0;
dE\[Infinity]=0;
dEH=0;
dA\[Infinity]=0;
dAH=0;
i=0;
For[k=0,k<kend,k++,
{F\[Infinity],FH,Fk,A\[Infinity],AH,Ak}=TENSORIALINTEGRATE[l,m,k,\[Omega]];
If[Abs[Fk]>Fmax,Fmax=Abs[Fk],Fmax=Fmax];
If[Abs[Fk]<\[CapitalEpsilon] Fmax|| Abs[Fk]<10^-25,i=i+1,i=0];
If[k==0,
F=Fk;A=Ak;AHN=AH;A\[Infinity]N=A\[Infinity];FHN=FH;F\[Infinity]N=F\[Infinity],
F=2*Fk;A=2*Ak;AHN=2*AH;A\[Infinity]N=2*A\[Infinity];FHN=2*FH;F\[Infinity]N=2*F\[Infinity];
]; 
Flux=Flux+F; 
AFlux=AFlux+A;
dAH=dAH+AHN;
dA\[Infinity]=dA\[Infinity]+A\[Infinity]N;
dEH=dEH+FHN;
dE\[Infinity]=dE\[Infinity]+F\[Infinity]N;
PrintTemporary[{l,m,k,Fk,Flux,\[Omega]}];
If[i>= 3,
Break[];
];
];
Return[{dE\[Infinity],dEH,Flux,dA\[Infinity],dAH,AFlux}];
]


lINTEGRATION[l_,\[CapitalEpsilon]_]:=Block[{str,data,dE\[Infinity]l,dEHl,dA\[Infinity]l,dAHl,EnergyFlux,KEND,lEnergy,kmax,m,lAngularMom,AFlux,dA\[Infinity],dAH,dE\[Infinity],dEH},
lEnergy=0;
{dE\[Infinity]l,dEHl,dA\[Infinity]l,dAHl}={0,0,0,0};
lAngularMom=0;
kmax=20;
For[m=-l,m<=l,m++, (*La condizione su m \[EGrave] inserita in tensorial integrate, dove viene messo subito a zero il contributo m+k\[Equal]0 ??*)
{dE\[Infinity],dEH,EnergyFlux,dA\[Infinity],dAH,AFlux}=kINTEGRATION[kmax,l,m,\[CapitalEpsilon]];
lEnergy=lEnergy+EnergyFlux;
lAngularMom=lAngularMom+AFlux;
{dE\[Infinity]l,dEHl,dA\[Infinity]l,dAHl}={dE\[Infinity]l,dEHl,dA\[Infinity]l,dAHl}+{dE\[Infinity],dEH,dA\[Infinity],dAH};
(*PrintTemporary[{m,EnergyFlux}];*)
];
data={{l,R0,\[CapitalDelta]R,e,L,Q2,x,dE\[Infinity]l,dEHl,dA\[Infinity]l,dAHl,lEnergy,lAngularMom}};
str=OpenAppend["GR_"<>FILEl];
WriteString[str,ExportString[data,"TSV"]];
Close[str];
Return[{dE\[Infinity]l,dEHl,lEnergy,dA\[Infinity]l,dAHl,lAngularMom}];
];



TensorialResults[lend_,\[CapitalEpsilon]_]:=Block[{l,lFlux,i,Fmax,LFlux,AFlux,lA,dE\[Infinity],dEH,dA\[Infinity],dAH,dE\[Infinity]l,dEHl,dA\[Infinity]l,dAHl},
lFlux=0;
LFlux=0;
AFlux=0;
dE\[Infinity]=0;
dEH=0;
dA\[Infinity]=0;
dAH=0;
i=0;
Fmax=0;
For[l=2,l<lend,l++,
(*PrintTemporary["Tensorial modes, l=",l];*)
{dE\[Infinity]l,dEHl,lFlux,dA\[Infinity]l,dAHl,lA}=lINTEGRATION[l,\[CapitalEpsilon]/10];
LFlux=LFlux+lFlux;
dE\[Infinity]=dE\[Infinity]+dE\[Infinity]l;
dEH=dEH+dEHl;
dA\[Infinity]=dA\[Infinity]+dA\[Infinity]l;
dAH=dAH+dAHl;
AFlux=AFlux+lA;
If[Abs[lFlux]>Fmax,Fmax=Abs[lFlux],Fmax=Fmax];
If[Abs[lFlux]<\[CapitalEpsilon] Fmax,i=i+1,i=0];
If[i>= 3,
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
\[ScriptCapitalK]=SetPrecision[(r^2+a^2)\[Omega]- m a,InputParameter];
\[Rho]=SetPrecision[-(1/(r-I a Cos[\[Theta]])),InputParameter];
\[Rho]b=SetPrecision[-(1/(r+I a Cos[\[Theta]])),InputParameter];
\[CapitalTheta]=SetPrecision[Sqrt[Q-Cot[\[Theta]]^2 L^2-a^2 Cos[\[Theta]]^2 (1-e^2)],InputParameter];
Lp[scalar_,s_]:=D[scalar,\[Theta]]-m /Sin[\[Theta]] scalar+a \[Omega] Sin[\[Theta]]scalar+s Cot[\[Theta]]scalar;
\[Alpha]lm = SetPrecision[(256 (2 M rH)^5 (\[Omega]-m \[CapitalOmega]h) ((\[Omega]-m \[CapitalOmega]h)^2+4 \[Epsilon]^2)((\[Omega]-m \[CapitalOmega]h)^2+16 \[Epsilon]^2) \[Omega]^3)/clm ,InputParameter]; 
\[Epsilon]= SetPrecision[Sqrt[M^2-a^2]/(4 M rH),InputParameter]; 
clm= SetPrecision[((\[Lambda]+2)^2+4 a m \[Omega] - 4 a^2  \[Omega]^2)(\[Lambda]^2+36 a m \[Omega] -36 a^2 \[Omega]^2)+(2\[Lambda]+3)(96 a^2 \[Omega]^2-48 a m \[Omega])+144 \[Omega]^2 (M^2-a^2),InputParameter];
dlm= SetPrecision[Sqrt[2M rH]((8-24 I M \[Omega] -16 M^2 \[Omega]^2)rH^2+(12 I a m - 16 M + 16 a m M \[Omega] + 24 I M^2 \[Omega] )rH - 4 a^2  m^2 - 12 I a m M + 8 M^2),InputParameter]; 
c0=SetPrecision[-12 I \[Omega] M +\[Lambda](\[Lambda]+2)-12 a \[Omega] (a \[Omega]-m ),InputParameter];
Dinf=SetPrecision[-4 \[Omega]^2/c0,InputParameter];
Bin=SetPrecision[W/(\[CapitalDelta] 2I \[Omega] Dinf),InputParameter];


Cnn=SetPrecision[Simplify[mp/(4\[CapitalSigma]^3 dt) (e(r^2+a^2)-a L)^2],InputParameter];
Cnmbplus=SetPrecision[Simplify[mp \[Rho]/(2Sqrt[2] \[CapitalSigma]^2 dt) (e(r^2+a^2)-a L)(I Sin[\[Theta]](a e - L/Sin[\[Theta]]^2)+\[CapitalTheta])],InputParameter];
Cmbmbplus=SetPrecision[Simplify[mp \[Rho]^2/(2\[CapitalSigma] dt) (I Sin[\[Theta]](a e - L/Sin[\[Theta]]^2)+\[CapitalTheta])^2],(*La precisione di tutti i parametri viene settata*)
InputParameter];
Cnmbminus=SetPrecision[Simplify[mp \[Rho]/(2Sqrt[2] \[CapitalSigma]^2 dt) (e(r^2+a^2)-a L)(I Sin[\[Theta]](a e - L/Sin[\[Theta]]^2)-\[CapitalTheta])],InputParameter];
Cmbmbminus=SetPrecision[Simplify[mp \[Rho]^2/(2\[CapitalSigma] dt) (I Sin[\[Theta]](a e - L/Sin[\[Theta]]^2)-\[CapitalTheta])^2],InputParameter];


L1=Lp[S[\[Theta],\[Phi]],1];
L2=Lp[S[\[Theta],\[Phi]],2];
L1L2=Lp[\[Rho]^-4 Lp[\[Rho]^3  S[\[Theta],\[Phi]],2],1];


Ann0=SetPrecision[Simplify[-2 Cnn/(\[Rho]^2 \[Rho]b \[CapitalDelta]^2) L1L2],InputParameter];(*La precisione di tutti i parametri viene settata*)
Anmb0plus=SetPrecision[Simplify[-2 Sqrt[2] Cnmbplus/(\[CapitalDelta] \[Rho]^3) ((I \[ScriptCapitalK]/\[CapitalDelta]-\[Rho]-\[Rho]b)L2+ (\[ScriptCapitalK]/\[CapitalDelta] ) a Sin[\[Theta]](\[Rho]b-\[Rho])S[\[Theta],\[Phi]])],InputParameter];(*La precisione di tutti i parametri viene settata*)
Anmb0minus=SetPrecision[Simplify[-2 Sqrt[2] Cnmbminus/(\[CapitalDelta] \[Rho]^3) ((I \[ScriptCapitalK]/\[CapitalDelta]-\[Rho]-\[Rho]b)L2+(\[ScriptCapitalK]/\[CapitalDelta])  a Sin[\[Theta]](\[Rho]b-\[Rho])S[\[Theta],\[Phi]])],InputParameter];(*La precisione di tutti i parametri viene settata*)
Ambmb0plus=SetPrecision[Simplify[S[\[Theta],\[Phi]] \[Rho]b/\[Rho]^3 Cmbmbplus((\[ScriptCapitalK]/\[CapitalDelta])^2+2I \[Rho] \[ScriptCapitalK]/\[CapitalDelta]+I D[\[ScriptCapitalK]/\[CapitalDelta],r])],InputParameter];
Ambmb0minus=SetPrecision[Simplify[S[\[Theta],\[Phi]] \[Rho]b/\[Rho]^3 Cmbmbminus((\[ScriptCapitalK]/\[CapitalDelta])^2+2I \[Rho] \[ScriptCapitalK]/\[CapitalDelta]+I D[\[ScriptCapitalK]/\[CapitalDelta],r])],InputParameter];
Anmb1plus=SetPrecision[Simplify[-2 Sqrt[2] Cnmbplus/(\[Rho]^3 \[CapitalDelta]) (L2+I a Sin[\[Theta]] S[\[Theta],\[Phi]](\[Rho]-\[Rho]b))],InputParameter]; (*Attenzione nel paper originale di Hughes c'\[EGrave] un fattore \[Rho] in pi\[UGrave]*)
Anmb1minus=SetPrecision[Simplify[-2 Sqrt[2] Cnmbminus/(\[Rho]^3 \[CapitalDelta]) (L2+I a Sin[\[Theta]] S[\[Theta],\[Phi]](\[Rho]-\[Rho]b))],InputParameter];
Ambmb1plus=SetPrecision[Simplify[2 (Cmbmbplus \[Rho]b S[\[Theta],\[Phi]] )/\[Rho]^3 (\[Rho]-I \[ScriptCapitalK]/\[CapitalDelta])],InputParameter];
Ambmb1minus=SetPrecision[Simplify[2 (Cmbmbminus \[Rho]b S [\[Theta],\[Phi]])/\[Rho]^3 (\[Rho]-I \[ScriptCapitalK]/\[CapitalDelta])],InputParameter];
Ambmb2plus=SetPrecision[Simplify[-S [\[Theta],\[Phi]] \[Rho]b/\[Rho]^3 Cmbmbplus ],InputParameter] ;
Ambmb2minus=SetPrecision[Simplify[-S[\[Theta],\[Phi]] \[Rho]b/\[Rho]^3 Cmbmbminus ],InputParameter]; 


A0plus=SetPrecision[Simplify[Ann0+Anmb0plus+Ambmb0plus],InputParameter];
A0minus=SetPrecision[Simplify[Ann0+Anmb0minus+Ambmb0minus],InputParameter];
A1plus=SetPrecision[Simplify[Ambmb1plus+Anmb1plus],InputParameter];
A1minus=SetPrecision[Simplify[Ambmb1minus+Anmb1minus],InputParameter];
A2plus=SetPrecision[Ambmb2plus,InputParameter];
A2minus=SetPrecision[Ambmb2minus,InputParameter];


Iplus=SetPrecision[R[r] A0plus- R'[r] A1plus+R''[r]A2plus,InputParameter];
Iminus=SetPrecision[R[r]A0minus-R'[r]A1minus+R''[r]A2minus,InputParameter];
W=SetPrecision[Rin[r]Rup'[r]-Rup[r]Rin'[r],InputParameter];


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

(* the following function has been used for a bug in the toolkit for x=-1, which has now been corrected 

Lfunc[a_,R0_,inc_]:=Block[{LL},
LL=KerrGeoConstantsOfMotion[a,R0,0,inc][["\[ScriptCapitalL]"]];
If[inc*LL<0,LL=-LL,LL=LL];
Return[LL]];
*)

Lfunc[a_,R0_,inc_]:=KerrGeoConstantsOfMotion[a,R0,0,inc][["\[ScriptCapitalL]"]]
efunc[a_,R0_,inc_]:=KerrGeoConstantsOfMotion[a,R0,0,inc][["\[ScriptCapitalE]"]]
Qfunc[a_,R0_,inc_]:=KerrGeoConstantsOfMotion[a,R0,0,inc][["\[ScriptCapitalQ]"]] 
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


s=-2;
OmegaPrecision=140;
a=a0;
eps=0.02;
M=1;
mp=1;
ScalarCharge=1;
If[Im[Lmin]==0,Lmin=Lmin,Lmin=0];
EpsScalar=0.001;
EpsTensorial=0.001;


Header={{"j" ,"i","R0-RLSO","R0","E","L","K","x","\[CapitalDelta]E\[Infinity]","\[CapitalDelta]EH","\[CapitalDelta]L\[Infinity]","\[CapitalDelta]LH","\[CapitalDelta]E","\[CapitalDelta]L"}};
FILEGR="../results/output_GR_"<>ToString[DecimalForm[a,3]]<>".tsv";
Export[FILEGR,Header];
lMax=20; (*20*)
kMax=20;(*20*)


ParallelTable[
i=Mod[ij,imax+1];
j=Floor[ij/(imax+1)];
Print[i," ",j];
R0=Rgrid[[j+1,i+1,1]];
\[CapitalDelta]R=Rgrid[[j+1,i+1,2]];
x=xgrid[[j+1]];
Print[R0," ",x];
e=efunc[a,R0,x];
L=Lfunc[a,R0,x];
Q=Qfunc[a,R0,x];
{dE\[Infinity],dEH,\[CapitalDelta]E,dA\[Infinity],dAH,\[CapitalDelta]L}=TensorialResults[lMax,EpsTensorial];
Print["{T\[CapitalDelta]E,T\[CapitalDelta]L}//Precision= ",{\[CapitalDelta]E//Precision,\[CapitalDelta]L//Precision}];
data={{j,i,\[CapitalDelta]R,R0,e,L,Q2,x,dE\[Infinity],dEH,dA\[Infinity],dAH,\[CapitalDelta]E,\[CapitalDelta]L}};
Print[data];
str=OpenAppend[FILEGR];
WriteString[str,ExportString[data,"TSV"]];
Close[str];
,{ij,0,(jmax)*(imax+1)-1}];(*ho tolto il +1 da jmax*)



Exit[]; 
