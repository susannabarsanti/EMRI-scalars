(* ::Package:: *)

(* ::Title:: *)
(*Code for the massive scalar fluxes computation*)


(* ::CodeText:: *)
(*If you make use of this code, please cite https://arxiv.org/abs/2212.03888, published in: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.131.051401.*)
(**)
(*Important Notes: *)
(*- Toolkit vs massless : factor \[Pi]/4*)
(**)
(*- THINGS TO IMPLEMENT: A FUNCTION FO R2 THAT VARIES WITH THE MASS OF THE SCALAR*)
(**)
(*If you implement a ParallelTable with the IntScal function, you do not obtain the right results (mixed parameteres in the parallelizations). Better a ForCycle. DistributeDefinitions[functions] when you parellalize is a good habit. ParallelDo better then ParallelTable? Parallelize? *)


SetSystemOptions["ParallelOptions"->"MKLThreadNumber"->10];
SetSystemOptions["ParallelOptions"->"ParallelThreadNumber"->10];
LaunchKernels[10];


Exit[];


<<SpinWeightedSpheroidalHarmonics`
<<KerrGeodesics`
<<Teukolsky`
ParallelNeeds["SpinWeightedSpheroidalHarmonics`"];
ParallelNeeds["KerrGeodesics`"];
ParallelNeeds["Teukolsky`"];


(* Master radial equation. *)
\[CapitalDelta]  = r^2+a^2- 2 M r ;
V[r] =  (\[Omega]-(a m)/\[Rho]^2)^2-\[CapitalDelta]/\[Rho]^8 (\[Lambda] \[Rho]^4+2 M r^3 +a^2 (r^2-4M r + a^2)) - (\[CapitalDelta] \[Mu]^2)/\[Rho]^2 ; 
\[Rho] =Sqrt[r^2+a^2 ]; 
drs = (a^2+r^2)/\[CapitalDelta]; (*drstar/dr*)
master=1/drs^2 u''[r]+1/drs D[1/drs,r]u'[r]+V[r]*u[r];


(* Series expansion at the horizon. *)
order=4;

p= \[Omega]- (m a)/(2 M (M+ Sqrt[M^2-a^2]))  ;
u[r_]:=Sum[Exp[-I p rt[r]]b[i](r-(M + Sqrt[M^2 - a^2]))^i,{i,0,order}];
seriesH=Collect[Series[master//.{rt[r]->0,rt'[r]->(r^2+a^2)/\[CapitalDelta],rt''[r]->D[(r^2+a^2)/\[CapitalDelta],r]},{r,M+Sqrt[M^2-a^2],order+2}],r-(M+Sqrt[M^2-a^2])];
Table[b[i]=Solve[Coefficient[seriesH,r-(M+Sqrt[M^2-a^2]),i]==0,b[i]][[1,1,2]],{i,1,order}];
uH[r_]=Sum[Exp[-I p rt[r]]b[i](r-(M+Sqrt[M^2-a^2]))^i//.rt[r]->r+(2 M (M + Sqrt[M^2 - a^2]))/(2  Sqrt[M^2 - a^2]) Log[(r-M - Sqrt[M^2 - a^2] )/(2 M)]- (2 M (M - Sqrt[M^2 - a^2]))/(2  Sqrt[M^2 - a^2]) Log[(r-M + Sqrt[M^2 - a^2] )/(2 M)],{i,0,order}];
Clear[u]


(* Series expansion at infinity. *)
order=4;

u[r_]:= Exp[I Sqrt[\[Omega]^2-\[Mu]^2] rt[r]] r^((I M \[Mu]^2)/Sqrt[-\[Mu]^2+\[Omega]^2]) Sum[d[i] (r)^{-i},{i,0,order}];
series\[Infinity]=Collect[Series[master//.{rt[r]->0,rt'[r]->(r^2+a^2)/\[CapitalDelta],rt''[r]->D[(r^2+a^2)/\[CapitalDelta],r]},{r,\[Infinity],order+2}],r];
Table[d[i+1]=Solve[Coefficient[series\[Infinity] * r^(-((I M \[Mu]^2)/Sqrt[-\[Mu]^2+\[Omega]^2]))//Simplify//Normal,r,-(i+2)]==0,d[i+1]][[1,1,2]],{i,0,order-1}];
u\[Infinity][r_]=(Exp[I Sqrt[\[Omega]^2-\[Mu]^2] rt[r]  ]r^((I M \[Mu]^2)/Sqrt[-\[Mu]^2+\[Omega]^2])  Sum[ d[i]r^{-i},{i,0,order}])//.rt[r]->r+(2 M (M + Sqrt[M^2 - a^2]))/(2  Sqrt[M^2 - a^2]) Log[(r-M - Sqrt[M^2 - a^2] )/(2 M)]- (2 M (M - Sqrt[M^2 - a^2]))/(2  Sqrt[M^2 - a^2]) Log[(r-M + Sqrt[M^2 - a^2] )/(2 M)];
Clear[u]


(* Source term for circular orbits. *)
\[Omega]p =M^(1/2)/(R0^(3/2)+ a M^(1/2)); 
\[CapitalSigma] = (r^2+a^2 (Cos[\[Theta]])^2)//.{\[Theta]-> \[Pi]/2};

v = Sqrt[M/r]; 
q = a/M; 
e= (1 - 2 v^2 + q v^3)/Sqrt[1 - 3 v^2+ 2 q v^3]; 
L = r v (1- 2 q v^3+ q^2 v^4)/Sqrt[1 - 3 v^2+ 2 q v^3];

dt =  e/\[CapitalSigma] ((r^2+a^2)^2/\[CapitalDelta]-a^2 (Sin[\[Theta]])^2) +( a L)/\[CapitalSigma] (1-(r^2+a^2)/\[CapitalDelta]) //.{\[Theta]-> \[Pi]/2} ; 
source= - 4 \[Pi]   di mp 1/(a^2+r^2)^(1/2)  Conjugate[S]/dt;
di=1; (*scalar charge set to 1. The fluxes scales with the charge as E\[Alpha]d^2*)


(*Function for the single mode flux computation.
- R: orbital radius;
- (l,m): modes indexes; 
- A: MBH spin
- ms: mass of the scalar*)

IntScal[R_,\[ScriptL]_,mu_,A_,ms_]:=Block[{
l=\[ScriptL], m=mu, \[Omega]=mu*\[Omega]p//.{r->R}, M=1,mp=1,R0=R,a=A,\[Mu]=ms,r,
r1,r2,\[Epsilon]=10^-5,
S,\[Lambda],W,drt,u\[Infinity]I,du\[Infinity]I,uHI,duHI,d0,b0,b,d,Y\[Infinity],YH,
eq,eqS,eqSh,rules,
Aout,Aouth,
dEd\[Sigma],dEd\[Sigma]\[Infinity],dEd\[Sigma]h},

eq=master;

r1=(M + Sqrt[M^2 - a^2])+2 \[Epsilon];

(*Check on r2. The values changes depending on the mass of the scalar and the radius the flux is computed. 
Conditions to impose so that u\[Infinity] don't go below 10^-50 when \[Mu]>\[Omega]. r2 furthest possible but if you put it too large then u\[Infinity] goes below 10^-50 and the flux at the horizon is not right.  
E.g. for \[Mu]=0.1, r=3, r2=700 works well, but not for \[Mu]=0.3, must be changed to r2=300.*)
If[\[Omega]>\[Mu], 
r2= 1200/Abs[\[Omega]], 
r2 =300;
PrintTemporary["r2 = ",N[r2,4]];
];

Print["\[Omega]= ",N[\[Omega],4]];
Print["Sqrt[\!\(\*SuperscriptBox[\(\[Omega]\), \(2\)]\)-\!\(\*SuperscriptBox[\(\[Mu]\), \(2\)]\)]= ",N[Sqrt[\[Omega]^2-\[Mu]^2],4]];

S=SpinWeightedSpheroidalHarmonicS[0,l,m,a Sqrt[\[Omega]^2-\[Mu]^2]][\[Pi]/2,0];
\[Lambda]=SpinWeightedSpheroidalEigenvalue[0,l,m,a Sqrt[\[Omega]^2-\[Mu]^2]]+2 m a Sqrt[\[Omega]^2-\[Mu]^2] -2 m a \[Omega]; 

Print["r=",N[R,3]];
Print["l=",l];
Print["m=",m];

uHI[b0_]=uH[r1]//.b[0]->b0;
duHI[b0_]=D[uH[r],r]//.{b[0]->b0,r->r1};

u\[Infinity]I[d0_]=u\[Infinity][r2][[1]]//.d[0]->d0;
du\[Infinity]I[d0_]=D[u\[Infinity][r][[1]],r]//.{d[0]->d0,r->r2};

(*u\[Infinity]I=xinf[r2];
du\[Infinity]I=D[xinf[r],r]//.{r->r2};*)

Print["N[u\[Infinity]I[1],4]= ", N[u\[Infinity]I[1],4]];
Print["N[du\[Infinity]I[1],4]=", N[du\[Infinity]I[1],4]];

rules={MaxSteps->Infinity,Method->"StiffnessSwitching",WorkingPrecision->65,AccuracyGoal->65,PrecisionGoal->65};

(* Forward integration from the horizon: invece che r1-r2 fai r1-R *)
YH= NDSolve[{eq==0,u[r1]==uHI[1],u'[r1]==duHI[1]},u,{r,r1,R},rules][[1]];

(* Inward integartion from the horizon: invece che r2-r1 fai r2-R*)
Y\[Infinity]= NDSolve[{eq==0,u[r2]==u\[Infinity]I[1],u'[r2]==du\[Infinity]I[1]},u,{r,r2,R},rules][[1]];

drt=(r^2+a^2)/\[CapitalDelta];

(* Define the wronskian *)
W=((u[r]//.YH)(u'[r]//.Y\[Infinity])-(u[r]//.Y\[Infinity])(u'[r]//.YH))/drt//.r->R;

eqS=((u[r]//.YH)*source)//.{r->R};

Aout=1/4 (W^-1 eqS)//.{r->R};

(*dEd\[Sigma]\[Infinity]=\[Omega]  Sqrt[\[Omega]^2-\[Mu]^2] Abs[Aout]^2 HeavisideTheta[\[Omega]^2-\[Mu]^2]*Abs[xinf[r2]]^2;*)

dEd\[Sigma]\[Infinity]=\[Omega]  Sqrt[\[Omega]^2-\[Mu]^2] Abs[Aout]^2 HeavisideTheta[\[Omega]^2-\[Mu]^2];

(*Print[N[Abs[xinf[r2]]^2,4]];(*conta solo quando il flusso all'infinito non \[EGrave] nullo ofc*)*)

eqSh=((u[r]//.Y\[Infinity])*source)//.{r->R};

Aouth= 1/4 (W^-1 eqSh)//.{r->R};

dEd\[Sigma]h=\[Omega] p Abs[Aouth]^2 ;

Return[{dEd\[Sigma]h, dEd\[Sigma]\[Infinity]}]
]


(*The toolkit fluxes differ from the massless IntScal output of \[Pi]/4, i.e  IntTool*\[Pi]/4 = IntScal.
In our formalism, the massless scalar fluxes differ from the toolkit by a factor 4, i.e. OUR = TOOL/4.
Finally, the correct fluxes in our work should be IntScal*1/\[Pi] for a single mode, and IntScal*2/\[Pi] if only the m>0 are computed.*)

IntTool[R_,\[ScriptL]_,mu_,A_]:=Block[{l=\[ScriptL],dEd\[Sigma],dEd\[Sigma]\[Infinity],dEd\[Sigma]h, m=mu,a,r0=R,orbit,mode,s,n,k},

If[A<0, a=-A; orbit=KerrGeoOrbit[a,r0,0,-1], a=A;orbit=KerrGeoOrbit[a,r0,0,1]] ;

s=0;n=0;k=0;
mode=TeukolskyPointParticleMode[s,l,m,n,k,orbit];

dEd\[Sigma]h = mode["Fluxes"]["Energy"]["\[ScriptCapitalH]"];
dEd\[Sigma]\[Infinity] = mode["Fluxes"]["Energy"]["\[ScriptCapitalI]"]; 

Return[{dEd\[Sigma]h,dEd\[Sigma]\[Infinity]}]
]


(*Radial grid. 
- nn: number of points; 
- eps: shift from the isco, the minimum radius will be ISCO+eps; 
- rfin: the maximum radius will be ISCO+rifn+eps *)

nn=3; 
eps = 2/100;
rfin = 15; 

(*Function for the radial spacing. In this way the grid is equally spaced in the variable u = 1/Sqrt[r-9/10 r_isco]*)
arrayu[a_]:= If[a>0, 
Range[1/Sqrt[1/10*KerrGeoISCO[a,1]+rfin+eps],1/Sqrt[1/10*KerrGeoISCO[a,1]+eps],
(1/Sqrt[1/10*KerrGeoISCO[a,1]+eps]-1/Sqrt[1/10*KerrGeoISCO[a,1]+rfin+eps])/nn], 
Range[1/Sqrt[1/10*KerrGeoISCO[-a,-1]+rfin+eps],1/Sqrt[1/10*KerrGeoISCO[-a,-1]+eps],
(1/Sqrt[1/10*KerrGeoISCO[-a,-1]+eps]-1/Sqrt[1/10*KerrGeoISCO[-a,-1]+rfin+2/100])/nn]];

arrayp[a_]:= If[a>0, 
Table[(1/arrayu[a][[i]])^2+9/10*KerrGeoISCO[a,1],{i,1,Length[arrayu[a]]}],
Table[(1/arrayu[a][[i]])^2+9/10*KerrGeoISCO[-a,-1],{i,1,Length[arrayu[a]]}]] ;


(*Function for the sum over the different modes. For the scalar case, only the l+m=even modes are present. 
- lmin: goes from 1 for the scalar circular case. It can be 0 for non-circular orbits; 
- lmax: used 10 for the circular massive case, it provides enough accuracy;
- a: MBH spin*)

(* IF YOU USE THIS FUNCTION YOU DON'T OBTAIN THE RIGHT RESULTS. 
tab[lmin_,lmax_,a_,mss_]:=Block[{nmult,multipoles,j,m,kk=1,x},

nmult=Length[Flatten[Table[Range[1,i],{i,lmin,lmax}]]];
multipoles=ConstantArray[0,nmult];

For[j=lmin,j<=lmax,j++,
	For[m=1,m<=j,m++,
	
multipoles[[kk]]={j,m};
kk = kk+1;
]
];

(*Select even multipoles for the scalar case*)
multipoles=Select[multipoles,EvenQ[#[[1]]+#[[2]]]&];
Print[Length[multipoles]];
Print[multipoles];

(*Factor 2/\[Pi] that accounts for the m<0 modes and the \[Pi] is due to difference in the normalization between the massive scalar integrator and the toolkit/massless one. The factor \[Pi] is not present in the data in the repository*)
x = Partition[Flatten[ParallelTable[{arrayp[a][[ics]],Sum[(2/\[Pi])*IntScal[arrayp[a][[ics]],multipoles[[l,1]],multipoles[[l,2]],a,mss],{l,1,Length[multipoles]}]},{ics,1,Length[arrayp[a]]}]],5];
Return[x];
];#THIS WAY OF PARELLALIZE IS WRONG. TRY PARALLELDO OPPURE PARELLALIZE.  

(*Exporting the data in a .dat file*)
ellmin=1;
ellmax=10; 
spin=0.9`100;
mass= 0.01`100 ; 
Export["massive_aXX.dat",tab[ellmin,ellmax,spin,mass]];*)



(*rag[a_]:={2.5`300};*)

tab[a_,mss_]:= ParallelTable[{arrayp[a][[i]],
IntScal[arrayp[a][[i]],1,1,1,a, mss]+
IntScal[arrayp[a][[i]],2,2,1,a, mss]+
IntScal[arrayp[a][[i]],3,1,1,a, mss]+
IntScal[arrayp[a][[i]],3,3,1,a, mss]+
IntScal[arrayp[a][[i]],4,2,1,a, mss]+
IntScal[arrayp[a][[i]],4,4,1,a, mss]+
IntScal[arrayp[a][[i]],5,1,1,a, mss]+
IntScal[arrayp[a][[i]],5,3,1,a, mss]+
IntScal[arrayp[a][[i]],5,5,1,a, mss]+
IntScal[arrayp[a][[i]],6,2,1,a, mss]+
IntScal[arrayp[a][[i]],6,4,1,a, mss]+
IntScal[arrayp[a][[i]],6,6,1,a, mss]+
IntScal[arrayp[a][[i]],7,1,1,a, mss]+
IntScal[arrayp[a][[i]],7,3,1,a, mss]+
IntScal[arrayp[a][[i]],7,5,1,a, mss]+
IntScal[arrayp[a][[i]],7,7,1,a, mss]+
IntScal[arrayp[a][[i]],8,2,1,a, mss]+
IntScal[arrayp[a][[i]],8,4,1,a, mss]+
IntScal[arrayp[a][[i]],8,6,1,a, mss]+
IntScal[arrayp[a][[i]],8,8,1,a, mss]+
IntScal[arrayp[a][[i]],9,1,1,a, mss]+
IntScal[arrayp[a][[i]],9,3,1,a, mss]+
IntScal[arrayp[a][[i]],9,5,1,a, mss]+
IntScal[arrayp[a][[i]],9,7,1,a, mss]+
IntScal[arrayp[a][[i]],9,9,1,a, mss]+
IntScal[arrayp[a][[i]],10,2,1,a, mss]+
IntScal[arrayp[a][[i]],10,4,1,a, mss]+
IntScal[arrayp[a][[i]],10,6,1,a, mss]+
IntScal[arrayp[a][[i]],10,8,1,a, mss]+
IntScal[arrayp[a][[i]],10,10,1,a,mss]},{i,1,Length[arrayp[a]]}]; 

spin=0.9`100;
mass= 0.01`100 ; 
Export["massive_aXX.dat",tab[spin,mass]];


Exit[];
