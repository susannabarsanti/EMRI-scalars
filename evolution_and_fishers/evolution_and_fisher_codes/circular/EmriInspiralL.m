(* ::Package:: *)

BeginPackage["EmriInspiralL`"];


(* ::Subsection:: *)
(*Function usage*)


KerrIsco::usage =
		"Isco for Kerr";
		
KerrOrbit::usage =
		"Equatorial orbital frequency and orbital energy for Kerr geodesics";
		
OrbitalEvoT::usage = 
		"Integrate the equations for the EMRI adiabatic evolution for a given (r0,\[Phi]0) and time T";
		
fftwave::usage = 
		"Fourier Transform of the signal in time";

overlap::usage = 
		"Overlap Maximised over time and phase between two waveforms";					
		
pattern::usage = 
		"Pattern function for LISA depending on time";
		
emritemplate::usage = 
		"Generate the full waveform containing +/x polarizations";

htime::usage = 
		"Sample the waveform in the time domain building the I-II channels for LISA";

findr0::usage = 
		"Function to find r0 for a given observing time";

window::usage = 
		"Function create a window to smooth the h(t) in time";
		
snrrwave::usage = 
		"Signal to noise ratio of a given signal";

derwave::usage = 
		"Driver for numerical derivatives of the waveform with respect to binary parameters";

finitederiv::usage = 
		"Derivative of the waveform in the fourier domain";

scalarN::usage = 
		"Compute the weighted inner product between two waveforms for a given PSD"

fisherN::usage = 
		"Compute the covariance matrix for the waveform h with derivative dh, for a given PSD";
		
invermatrix::usage = 
		"Invert matrix with singular value decomposition and choiche of # of pivot to remove";

integralN::usage = 
		"Numerical integral with different methods";
		
SnLISA::usage = 
		"Analytical fit of the LISA psd";

CoeffDer::usage =
		"Coefficients of finite derivatives";

transition::usage =
		"Thorne - Ori transition radius";

dopplerf::usage =
		"Doppler frequnecy";


dEGRdr = Import["./dEdrGR09HPI.mx"];


Begin["`Private`"];


precision = 200;
year = 365;
day = 3600 24;
GN = 6674 10^-14;
cc = 299792458 10^-3;
msun = 1477 10^-3;
Rs = 1495978707/10;
pc = 96939420213600/\[Pi];


(* ::Subsection:: *)
(*EMRI dynamics*)


(* ::Subsubsection::Closed:: *)
(*Find the transition region*)


transition[m1_,m2_]:=Block[{solX,X,T,\[Alpha],\[Beta],\[Kappa],R0,\[Tau]0,\[Eta],T0=-3,Tf=3,rI=10,\[Phi]I=0.78},

(* parameters calibrated for Schwarzschild *)
(*{\[Alpha],\[Beta],\[Kappa]}={7716/10^7,1604/10^5,1955/10^5};*)
(* parameters calibrated for Kerr a=0.9 *)
{\[Alpha],\[Beta],\[Kappa]}={3447/10^5,9039/10^5,4214/10^4};
{R0,\[Tau]0}={(\[Beta] \[Kappa])^(2/5) \[Alpha]^(-3/5),(\[Alpha] \[Beta] \[Kappa])^(-1/5)};

\[Eta]=m2/m1;

solX=NDSolve[{X''[T]==-X[T]^2-T,X[T0]==Sqrt[-T0],X'[T0]==-(1/(2 Sqrt[-T0]))},X,{T,T0,Tf},MaxSteps->Infinity,
		Method->"StiffnessSwitching",WorkingPrecision->precision];

Return[\[Eta]^(2/5) R0 X[0]//.solX[[1]]];

]


(* ::Subsubsection::Closed:: *)
(*Kerr ISCO*)


KerrIsco[\[Chi]_]:=Block[{Z1,Z2},

Z1=1+(1-\[Chi]^2)^(1/3) ((1+\[Chi])^(1/3)+(1-\[Chi])^(1/3));
Z2=Sqrt[3\[Chi]^2+Z1^2];

Return[3+Z2-Sign[\[Chi]]Sqrt[(3-Z1)(3+Z1+2Z2)]];

]


(* ::Subsubsection::Closed:: *)
(*Kerr orbital parameters*)


KerrOrbit[a_,r_]:=Block[{\[Omega],energy,denergy},

\[Omega] = 1/(r^(3/2)+ a);

energy = (r^2-2 r+ a Sqrt[r])/(r Sqrt[r^2-3r+ 2a Sqrt[r]]);
denergy = (-3 a^2+8 a Sqrt[r]+(-6+r) r)/(2 (2 a+(-3+r) Sqrt[r]) r^(3/2) Sqrt[2 a Sqrt[r]+(-3+r) r]);

Return[{\[Omega],energy,denergy}]

]


(* ::Subsubsection::Closed:: *)
(*Find starting radius for given source and T*)


findr0[m1_,m2_,spin_,rend_,T_]:=Block[{q=m2/m1,factsec,factday,fact,tend,rini,r,\[ScriptCapitalE],eqns,t},

\[ScriptCapitalE][r_] = dEGRdr[r,spin];
 		
factsec = cc/(m1 msun);
factday = factsec day ; (* change time in days in dimensionless time *)
fact = factsec;

tend = T day;
(* Stop the integration when it's at the transition T = 0 *)

eqns = r'[t]==-q fact \[ScriptCapitalE][r[t]]/KerrOrbit[spin,r[t]][[3]]; 

rini = NDSolveValue[{eqns,r[tend]==rend},r,{t,tend,0},
				    MaxSteps->10^6,Method->"StiffnessSwitching",WorkingPrecision->precision+1][0];

Return[rini];

]


(* ::Subsubsection::Closed:: *)
(*Adiabatic evolution for a given T*)


OrbitalEvoT[m1_,m2_,spin_,r0_,\[Phi]0_,T_]:=Block[{solr,t,sol,r,\[Psi],factday,factsec,tend,q=m2/m1,fact,eqns,IC,\[ScriptCapitalE]},

\[ScriptCapitalE][r_] = dEGRdr[r,spin];

factsec = cc/(m1 msun);
factday = factsec day ; (* change time in days in dimensionless time *)
fact = factsec;

tend = T day;

eqns = {r'[t]==-q fact \[ScriptCapitalE][r[t]]/KerrOrbit[spin,r[t]][[3]],
			\[Psi]'[t]==fact KerrOrbit[spin,r[t]][[1]]};

IC = {r[0]==r0,\[Psi][0]==\[Phi]0};

sol=NDSolveValue[{eqns,IC},{\[Psi],r},{t,0,tend},
				    Method->"StiffnessSwitching",WorkingPrecision->precision,MaxSteps->10^6];

Return[{sol[[1]],sol[[2]]}];

]


(* ::Subsection:: *)
(*Waveforms*)


(* ::Subsubsection::Closed:: *)
(*Pattern functions*)


pattern[pol_,\[Theta]s_,\[Phi]s_,\[Theta]L_,\[Phi]L_,t_]:=Block[{\[Gamma]=Sqrt[3]/2,T,\[Phi]T,Lz,LN,NLz,\[Psi],\[Theta],\[Phi],Fplus,Fcros},

T = year day;
\[Phi]T = 2\[Pi] t/T;

\[Theta] = ArcCos[1/2 Cos[\[Theta]s]-\[Gamma] Sin[\[Theta]s]Cos[\[Phi]T-\[Phi]s]];
(*\[Phi] = \[Phi]T + ArcTan[(Sqrt[3]Cos[\[Theta]s]+Sin[\[Theta]s]Cos[\[Phi]T-\[Phi]s])/(2Sin[\[Theta]s]Sin[\[Phi]T-\[Phi]s])]+If[pol==1,0,-\[Pi]/4];*)

\[Phi] = \[Phi]T + ArcTan[2Sin[\[Theta]s]Sin[\[Phi]T-\[Phi]s],Sqrt[3]Cos[\[Theta]s]+Sin[\[Theta]s]Cos[\[Phi]T-\[Phi]s]]+Which[pol==1,0,pol==2,-\[Pi]/4,pol!=1\[Or]pol!= 2,Print["BadPatternFunctionOption"];Return[]];

Lz = 1/2 Cos[\[Theta]L]-\[Gamma] Sin[\[Theta]L]Cos[\[Phi]T-\[Phi]L];
LN = Cos[\[Theta]L]Cos[\[Theta]s]+Sin[\[Theta]L]Sin[\[Theta]s]Cos[\[Phi]L-\[Phi]s];
NLz = 1/2 Sin[\[Theta]L]Sin[\[Theta]s]Sin[\[Phi]L-\[Phi]s]
	-\[Gamma] Cos[\[Phi]T](Cos[\[Theta]L]Sin[\[Theta]s]Sin[\[Phi]s]-Cos[\[Theta]s]Sin[\[Theta]L]Sin[\[Phi]L])
	-\[Gamma] Sin[\[Phi]T]((Cos[\[Theta]s]Sin[\[Theta]L]Cos[\[Phi]L]-Cos[\[Theta]L]Sin[\[Theta]s]Cos[\[Phi]s]));

(*\[Psi] = ArcTan[(Lz-LN Cos[\[Theta]s])/NLz];*)
\[Psi] = ArcTan[NLz,Lz-LN Cos[\[Theta]]];

Fplus = (1/2 (1+Cos[\[Theta]]^2)Cos[2\[Phi]]Cos[2\[Psi]]-Cos[\[Theta]]Sin[2\[Phi]]Sin[2\[Psi]]);
Fcros = (1/2 (1+Cos[\[Theta]]^2)Cos[2\[Phi]]Sin[2\[Psi]]+Cos[\[Theta]]Sin[2\[Phi]]Cos[2\[Psi]]);

Return[{Fplus,Fcros}];

]


(* ::Subsubsection::Closed:: *)
(*Window function to temper the FFT*)


window[n_,nn_,\[Alpha]_]:=Block[{w1,w2},

If[\[Alpha] == 0, 1,

	w1=1/2 (1+Cos[\[Pi]((2n)/(\[Alpha](nn-1))-1)]);
	w2=1/2 (1+Cos[\[Pi]((2n)/(\[Alpha](nn-1))-2/\[Alpha]+1)]);

	Piecewise[{{w1,0<=n<=\[Alpha] (nn-1)/2},{1,\[Alpha] (nn-1)/2<=n<=(nn-1)(1-\[Alpha]/2)},
	{w2,(nn-1)(1-\[Alpha]/2)<=n<=nn-1}}]

		]
]


(* ::Subsubsection::Closed:: *)
(*Doppler frequency*)


dopplerf[m1_,m2_,spin_,r0_,\[Phi]0_,\[Theta]s_,\[Phi]s_,T_,t_]:=Block[{M,\[Mu],fdopp,\[Omega],x,evolution,\[Omega]D,fD},

M = m1 msun; (* total mass *)
\[Mu] = m2 msun; (* reduced mass *)

evolution = OrbitalEvoT[m1,m2,spin,r0,\[Phi]0,T];

\[Omega]D = D[\[Omega][x] Rs/cc Sin[\[Theta]s]Cos[2\[Pi] x/(day year)-\[Phi]s],x]; (* doppler angular frequency *)

fD = (\[Omega]D/\[Pi])//.{\[Omega][x]->evolution[[1]]'[t],\[Omega]'[x]->evolution[[1]]''[t]};

Return[fD//.x->t];

]


(* ::Subsubsection:: *)
(*h+x polarizations*)


emritemplate[m1_,m2_,spin_,r0_,\[Phi]0_,\[Theta]s_,\[Phi]s_,\[Theta]k_,\[Phi]k_,T_,t_]:=Block[{evolution,\[Phi],\[Mu],M,\[Phi]D,r,\[Omega],LN},

M = m1 msun; (* primary mass *)
\[Mu] = m2 msun; (* secondary mass *)

evolution = OrbitalEvoT[m1,m2,spin,r0,\[Phi]0,T];

\[Omega] = evolution[[1]]'[t]; (* angular frequency *)

\[Phi]D = \[Omega] Rs/cc Sin[\[Theta]s]Cos[2\[Pi] t/(day year)-\[Phi]s]; (* doppler phase *)

\[Phi] = evolution[[1]][t]+\[Phi]D; (* total phase *)
r = evolution[[2]][t];

LN = Cos[\[Theta]k]Cos[\[Theta]s]+Sin[\[Theta]k]Sin[\[Theta]s]Cos[\[Phi]k-\[Phi]s];

Return[2 (M \[Omega]/cc)^(2/3) \[Mu] {Cos[2 \[Phi]](1+LN^2),-Sin[2 \[Phi]]2 LN}];

]


(* ::Subsubsection:: *)
(*Sampling in time*)


htime[m1_,m2_,spin_,r0_,\[Theta]s_,\[Phi]s_,\[Phi]0_,\[Theta]k_,\[Phi]k_,observingtime_]:=Block[{win,nn,timevec,T,\[Delta]t,t,
i,fmax,h,fdopp,hp,hx,hplus,hcross,rmax,fisco,FplusI,
FcrossI,FplusII,FcrossII,risco,fsampling},

T = day observingtime; (* total time in seconds *)

rmax = OrbitalEvoT[m1,m2,spin,r0,\[Phi]0,observingtime][[2]][T]; (* final radius after T *)	
risco = KerrIsco[spin]; (* isco radius/M *)
fdopp = dopplerf[m1,m2,spin,r0,\[Phi]0,\[Theta]s,\[Phi]s,observingtime,T]; (* max fdop @ rmax *)

fmax = cc KerrOrbit[spin,rmax][[1]]/(m1 msun \[Pi])+fdopp; (* f @ rmax + fdop @ rmax *)
fisco = cc KerrOrbit[spin,risco][[1]]/(m1 msun \[Pi]); (* Isco Kerr *)

Print["\!\(\*SubscriptBox[\(f\), \(isco\)]\) = ",N[fisco,4],"Hz"];
Print["\!\(\*SubscriptBox[\(f\), \(dopp\)]\) = ",N[fdopp,4],"Hz"];
Print["\!\(\*SubscriptBox[\(f\), \(max\)]\) = ",N[fmax,4],"Hz"];

fsampling = fisco;

nn = 2^(Ceiling[Log[2 fsampling T]/Log[2]]);
\[Delta]t = T/(nn-1); (* time step *)

If[fsampling < fmax, Print["Samplign freq smaller than fmax"];Return[]];
If[1/2/\[Delta]t < fsampling, Print["Nymquist criterion failed"];Return[]]; (*fNy = fs/2 > fmax *)
(* sampling is alwasy up to the folding freq fs/2 = N/T/2 = 1/2/\[Delta]t= fNy*)

timevec =  N[Table[i \[Delta]t,{i,0,nn-1}],precision];

If[timevec[[-1]]/day != observingtime, Print["end time \[NotEqual] T"];Return[]];

{hplus[t_],hcross[t_]} = emritemplate[m1,m2,spin,r0,\[Phi]0,\[Theta]s,\[Phi]s,\[Theta]k,\[Phi]k,observingtime,t]; (* analytical waveform *)

win = Table[window[i,nn,5/100],{i,0,nn-1}];

{hp,hx} = {win*hplus[timevec],win*hcross[timevec]};

If[hp[[1]]!=hp[[-1]] \[Or] hx[[1]]!=hx[[-1]],Print["Initial/final points of h are different"]];

{FplusI,FcrossI} = pattern[1,\[Theta]s,\[Phi]s,\[Theta]k,\[Phi]k,timevec];
{FplusII,FcrossII} = pattern[2,\[Theta]s,\[Phi]s,\[Theta]k,\[Phi]k,timevec];

Return[{timevec, FplusI*hp+FcrossI*hx,FplusII*hp+FcrossII*hx}];

]


(* ::Subsubsection::Closed:: *)
(*DFT of the signal*)


fftwave[hIt_,hIIt_,time_]:=Block[{hIf,hIIf,nsample,\[Delta]f,\[Delta]t,freqVec,i},

(* Make the # of points even *)

If[EvenQ[hIt//Length]==True,hIf=Fourier[hIt,FourierParameters->{1, -1}],
hIf=Fourier[Drop[hIt,1],FourierParameters->{1, -1}]];

If[EvenQ[hIIt//Length]==True,hIIf=Fourier[hIIt,FourierParameters->{1, -1}],
hIIf=Fourier[Drop[hIIt,1],FourierParameters->{1, -1}]];
If[Length[hIf]!=Length[hIIf],Print["DIFFERENT # OF SAMPLES"];Return[]];

nsample=Length[hIf];(* total number of points *)

\[Delta]t = time[[2]]-time[[1]]; 
\[Delta]f = 1/(\[Delta]t nsample);(* freqeuncy step *)

Print["\[Delta]t = ",N[\[Delta]t,4],"s - ","\[Delta]f = ",N[\[Delta]f,4],"Hz - ","# samples = ",nsample];

freqVec = Table[i \[Delta]f,{i,0,nsample/2}]; (* array of frequencies *)

If[1/2/\[Delta]t != freqVec[[-1]], Print["fNy \[NotEqual] fmax"];Return[]];

(* I am excluding the Nyquist frequency [nsample/2 + 1 \[Rule] nsample/2] The DC will be cut out by fmin*)

Return[{Thread[{freqVec[[1;;-2]], \[Delta]t hIf[[1;;nsample/2]]}],
			Thread[{freqVec[[1;;-2]], \[Delta]t hIIf[[1;;nsample/2]]}]}];
]


(* ::Subsection:: *)
(*Fisher matrix derivatives and data analysis tools*)


(* ::Subsubsection::Closed:: *)
(*SNR*)


snrrwave[m1_,m2_,spin_,r0_,\[Phi]0_,\[Theta]s_,\[Phi]s_,\[Theta]k_,\[Phi]k_,dL_,T_,WDT_]:=Block[{fmin=10^-4,rmax,fdopp,fmax,d,\[Rho],hI,hII,timevec,fact,hp,hx,FpI,FxI,FpII,FxII,hIt,hIIt},

d=10^6 dL pc; (* Luminosity distance *)

{timevec,hIt,hIIt}=htime[m1,m2,spin,r0,\[Theta]s,\[Phi]s,\[Phi]0,\[Theta]k,\[Phi]k,T];

{hI,hII} = fftwave[hIt,hIIt,timevec];

rmax = OrbitalEvoT[m1,m2,spin,r0,\[Phi]0,T][[2]][T day]; (* final radius after T *)	
fdopp = dopplerf[m1,m2,spin,r0,\[Phi]0,\[Theta]s,\[Phi]s,T,T day]; (* max fdop @ rmax *)
fmax = cc KerrOrbit[spin,rmax][[1]]/(m1 msun \[Pi])+fdopp; (* f @ rmax + fdop @ rmax *)

Print["frequency interval = [",N[fmin,4],N[fmax,4],"]Hz"];

fmax = -1;

\[Rho] = 1/d Sqrt[scalarN[hI,hI,"S",WDT,fmin,fmax]+scalarN[hII,hII,"S",WDT,fmin,fmax]];

Return[\[Rho]]

]


(* ::Subsubsection::Closed:: *)
(*Coefficients finite derivatives*)


CoeffDer[order_,derorder_,mode_]:=Block[{A,V,coeff},

A=Which[mode=="C",
Join[{ConstantArray[1,order]},Table[Range[-(order-1)/2,(order-1)/2]^j,{j,1,order-1}]],mode=="F",Join[{ConstantArray[1,order]},Table[Range[0,order-1]^j,{j,1,order-1}]]];

V = derorder! Table[KroneckerDelta[j,derorder],{j,0,order-1}];

coeff = DeleteCases[LinearSolve[A,V],0];

Return[coeff]

]


(* ::Subsubsection::Closed:: *)
(*Finite derivatives for the intrinsic parameters*)


finitederiv[distance_,i_,parc_,\[Epsilon]_,T_,order_]:=Block[{\[Delta],fact,dhI,dhII,incr,par\[Delta]\[Lambda],h\[Delta]\[Lambda],coeff,dh\[Lambda]I,dh\[Lambda]II,derorder=1,Npoints},

fact=If[i<=2,parc[[i]],1]; (* take the log derivative for the masses *)

Print[i];

Npoints = order-1;

If[parc[[i]]!=0,
	h\[Delta]\[Lambda] = ConstantArray[0,Npoints];
		par\[Delta]\[Lambda] = ConstantArray[parc,Npoints];
			\[Delta] = DeleteCases[Table[(j-Npoints/2) \[Epsilon][[i]],{j,0,Npoints}],0];
				Table[par\[Delta]\[Lambda][[j,i]] = parc[[i]](1+\[Delta][[j]]),{j,1,Npoints}],
		
	h\[Delta]\[Lambda] = ConstantArray[0,order];
		par\[Delta]\[Lambda] = ConstantArray[parc,order];
			\[Delta] = Table[j \[Epsilon][[i]],{j,0,Npoints}];
				Table[par\[Delta]\[Lambda][[j,i]] = \[Delta][[j]],{j,1,order}]
];

Table[h\[Delta]\[Lambda][[j]] = htime[Sequence@@Join[par\[Delta]\[Lambda][[j]],{T}]],{j,1,Length[par\[Delta]\[Lambda]]}];

(* Compute coefficients for finite differences *)
coeff = If[parc[[i]]!=0,CoeffDer[order,derorder,"C"],CoeffDer[order,derorder,"F"]];

(* Pre-factor for the derivative *)
incr = fact/distance/If[parc[[i]]!=0,parc[[i]] \[Epsilon][[i]],\[Epsilon][[i]]];

dh\[Lambda]I = incr {coeff} . h\[Delta]\[Lambda][[All,2]];
dh\[Lambda]II = incr {coeff} . h\[Delta]\[Lambda][[All,3]];

(*Return[{Thread[{h\[Delta]\[Lambda][[1,1]],dh\[Lambda]I[[1]]}],Thread[{h\[Delta]\[Lambda][[1,1]],dh\[Lambda]II[[1]]}]}];*)
{dhI,dhII} = fftwave[dh\[Lambda]I[[1]],dh\[Lambda]II[[1]],h\[Delta]\[Lambda][[1,1]]];

Return[{dhI,dhII}];

]


derwave[\[Epsilon]_,m1_,m2_,spin_,r0_,\[Phi]0_,\[Theta]s_,\[Phi]s_,\[Theta]k_,\[Phi]k_,dL_,observingtime_]:=Block[{timevec,vecder,
i,d,\[Rho],hI,hII,dhI,dhII,npar,hplus,hcross,hp,hx,FplusI,FcrossI,FplusII,FcrossII,t,T,fisco,
nn,\[Delta]t,order,win,risco,pos,LN,fdopp,fsampling,rmax,fmax},

d=10^6 dL pc; (* Luminosity distance *)

npar = Length[{m1,m2,spin,r0,\[Phi]0,\[Theta]s,\[Phi]s,\[Theta]k,\[Phi]k,d}];

T = day observingtime; (* total time in seconds *)

rmax = OrbitalEvoT[m1,m2,spin,r0,\[Phi]0,observingtime][[2]][T]; (* final radius after T *)	
risco = KerrIsco[spin];
fisco = cc KerrOrbit[spin,risco][[1]]/(m1 msun \[Pi]); (* Isco Kerr *)

fmax = cc KerrOrbit[spin,rmax][[1]]/(m1 msun \[Pi])+fdopp; (* f @ rmax + fdop @ rmax *)
fdopp = dopplerf[m1,m2,spin,r0,\[Phi]0,\[Theta]s,\[Phi]s,observingtime,T]; (* max fdop @ rmax *)

fsampling = fisco;

nn = 2^(Ceiling[Log[2 fsampling T]/Log[2]]);
\[Delta]t = T/(nn-1); (* time step *)

If[fsampling < fmax, Print["Samplign freq smaller than fmax"];Return[]];
If[1/2/\[Delta]t < fsampling, Print["Nymquist criterion failed"];Return[]]; (*fNy = fs/2 > fmax *)

timevec = N[Table[i \[Delta]t,{i,0,nn-1}],precision];
If[timevec[[-1]]/day != observingtime, Print["end time \[NotEqual] T"];Return[]];

(* Build the signal *)
{hplus[t_],hcross[t_]} = emritemplate[m1,m2,spin,r0,\[Phi]0,\[Theta]s,\[Phi]s,\[Theta]k,\[Phi]k,observingtime,t]; 

LN = Cos[\[Theta]k]Cos[\[Theta]s]+Sin[\[Theta]k]Sin[\[Theta]s]Cos[\[Phi]k-\[Phi]s];

win = Table[window[i,nn,5/100],{i,0,nn-1}];

{hp,hx} = {win*hplus[timevec],win*hcross[timevec]};

{FplusI,FcrossI} = pattern[1,\[Theta]s,\[Phi]s,\[Theta]k,\[Phi]k,timevec];
{FplusII,FcrossII} = pattern[2,\[Theta]s,\[Phi]s,\[Theta]k,\[Phi]k,timevec];

{hI,hII} = fftwave[FplusI hp+FcrossI hx,FplusII hp+FcrossII hx,timevec];

\[Rho] = snrrwave[m1,m2,spin,r0,\[Phi]0,\[Theta]s,\[Phi]s,\[Theta]k,\[Phi]k,dL,observingtime,12];

Print["SNR = ",\[Rho]];

dhI = dhII = ConstantArray[0,npar]; (* all paramse + dL *)

order = {11,11,11,11,11,9,9,9,9};
pos = {1,2,3,4,5,6,7,9,10};

If[Length[order]!=Length[pos],Print["Derivatives don't match parameters"];Return[],Continue];

vecder=Table[finitederiv[d,pos[[i]],
{m1,m2,spin,\[Alpha], r0,\[Theta]s,\[Phi]s,\[Phi]0,\[Theta]k,\[Phi]k},\[Epsilon],observingtime,order[[i]]],{i,1,Length[pos]}];(* all derivatives but distance and \[Phi]0 *)

Table[{dhI[[pos[[i]]]],dhII[[pos[[i]]]]}=vecder[[i]],{i,1,Length[pos]}];

PrintTemporary["distance derivative"]; (* log derivative *)
dhI[[-1]]=Thread[{hI[[All,1]],-(1/d)hI[[All,2]]}];
dhII[[-1]]=Thread[{hI[[All,1]],-(1/d)hII[[All,2]]}]; 

PrintTemporary["\[Phi]0 derivative"];
{hI,hII} = fftwave[-4FcrossI hp/(1+LN^2)LN+FplusI/LN hx(1+LN^2),-4FcrossII hp/(1+LN^2)LN+FplusII/LN hx(1+LN^2),timevec];
dhI[[7]]=Thread[{hI[[All,1]],1/d hI[[All,2]]}];
dhII[[7]]=Thread[{hI[[All,1]],1/d hII[[All,2]]}]; 

Return[{dhI,dhII}];

]


(* ::Subsubsection::Closed:: *)
(*LISA fit PSD*)


SnLISA[f_,T_]:=Block[{A= 9 10^-45,\[Beta],\[Alpha],\[Gamma],k,fk,fstar=1909 10^-5,L=25 10^8,Poms,Pacc,WDN},


{\[Beta],\[Alpha],\[Gamma],k,fk}=Which[
	T==0,{0,0,0,0,0},
	T==12,{292,171*10^-3,1680,1020,215*10^-5},
	T==6,{243,133*10^-3,917,482,258*10^-5}
];

Poms=(15 10^-12)^2 (1+((2 10^-3)/f)^4);
Pacc=(3 10^-15)^2 (1+((4 10^-4)/f)^2)(1+(f/(8 10^-3))^4);

WDN = If[T==0,0,A f^(-7/3) Exp[-f^\[Alpha]+\[Beta] f Sin[k f]](1+Tanh[\[Gamma](fk-f)])];

Return[10/(3L^2) (Poms+2(1+Cos[f/fstar]^2) Pacc/(2\[Pi] f)^4)(1+6/10 (f/fstar)^2)+WDN];

]


(* ::Subsubsection::Closed:: *)
(*Numerical scalar product*)


scalarN[wave1_,wave2_,mode_,WDT_,fmin_,fmax_]:=Block[{func,df,nini,nfin,h1,h2,\[Delta]f,freqVec,psd,ntot,noise,prec},

If[Length[wave1]!=Length[wave2],Print["Different # of samples"];Return[]];

nini=First@Position[wave1[[All,1]],_?(#>=fmin&)][[1]];
nfin= If[fmax == -1, -1,First@Position[wave1[[All,1]],_?(#<=fmax&)][[-1]]];

h1 = wave1[[nini;;nfin,2]];
h2 = wave2[[nini;;nfin,2]];

df = wave2[[2,1]]-wave2[[1,1]];
freqVec = wave1[[nini;;nfin,1]];

psd = SnLISA[freqVec,WDT];

func = h1 Conjugate[h2]/psd;

Return[4*Re[integralN[func,df,mode]]];

]


(* ::Subsubsection::Closed:: *)
(*Numerical integral*)


integralN[data_,dx_,mode_]:=Block[{int,i,nn},

nn = Length[data];

nn = If[mode=="S",If[EvenQ[nn],nn-1,nn],nn];

int = Which[

	mode=="R",Total[data,Method->"CompensatedSummation"],
	mode=="T",(data[[1]]/2 + data[[-1]]/2 + Sum[data[[i]], {i, 2, nn - 1}]),
	mode=="S",1/3*Sum[data[[2i-1]]+4data[[2i]]+data[[2i+1]], {i, 1, (nn-1)/2}],
	mode=="B",2/45 (7(data[[1]]+data[[nn]])+32Total[data[[Range[2,nn-1,2]]],Method->"CompensatedSummation"]
			+12Total[data[[Range[3,nn-2,4]]],Method->"CompensatedSummation"]+14Total[data[[Range[5,nn-2,4]]],Method->"CompensatedSummation"])];

Return[dx*int];

]


(* ::Subsubsection::Closed:: *)
(*Covariance with SVD*)


invermatrix[matrix_,thresh_,normalization_]:=Block[{dm,norm,U,S,Vh,matrixnorm,Kval,u,Inv},

dm=Sqrt[Diagonal[matrix]];
norm=If[TrueQ[normalization],dm\[TensorProduct]dm,1];

matrixnorm=matrix/norm;
{U,S,Vh}=SingularValueDecomposition[matrixnorm];

Kval = If[10^thresh<1,Count[Diagonal[S],u_/;u>10^thresh],Length[Diagonal[S]]-thresh];

Inv=Vh[[All,1;;Kval]] . DiagonalMatrix[1/Diagonal[S][[1;;Kval]]] . ConjugateTranspose[U][[1;;Kval,All]];

Return[Inv/norm]

]


(* ::Subsubsection::Closed:: *)
(*Fisher Matrix Fully Numerical*)


fisherN[dh_,mode_,WDT_,fmin_,fmax_]:=Block[{\[CapitalGamma],f,dim,i,j,k},

dim = Length[dh];

\[CapitalGamma] = ConstantArray[0,{dim,dim}];

(* ---- Compute the Fisher and Covariance Matrix ---- *)

Monitor[For[i=1,i<=dim,i++,
	For[j=1,j<=i,j++,
		
		\[CapitalGamma][[i,j]] = If[i==j,scalarN[dh[[i]],dh[[j]],mode,WDT,fmin,fmax],
								2scalarN[dh[[i]],dh[[j]],mode,WDT,fmin,fmax]];
		]
	],\[CapitalGamma]//MatrixForm];

\[CapitalGamma] = Normal[Symmetrize[\[CapitalGamma]]];

Return[\[CapitalGamma]]

]


(* ::Subsubsection::Closed:: *)
(*Overlap Maximised over time and phase*)


overlap[wave1_,wave2_,WDT_]:=Block[{cc=299792,nsample1=Length[wave1],
nsample2=Length[wave2],h1,h2,\[Delta]f,psd,\[ScriptCapitalO]1,\[ScriptCapitalO]2,\[ScriptCapitalO]12,nini,prec,noise},

If[nsample1!=nsample2,Print["Different # of samples"]];

nini=2;(*First@Position[wave1[[All,1]],_?(#>=10^-4&)][[1]];*)

h1 = wave1[[nini;;-1,2]];
h2 = wave2[[nini;;-1,2]];
\[Delta]f = wave1[[2,1]]-wave1[[1,1]];

psd = SnLISA[wave1[[nini;;-1,1]],WDT];
Print["Precision templates = ",{Precision[h1],Precision[h2]}];

\[ScriptCapitalO]1 = Re[Total[h1 Conjugate[h1]/psd]];
\[ScriptCapitalO]2 = Re[Total[h2 Conjugate[h2]/psd]];

\[ScriptCapitalO]12 = Max[Abs[InverseFourier[Union[{0},(h1 Conjugate[h2])/psd],FourierParameters->{-1, -1}]]];

Return[{\[ScriptCapitalO]12,\[ScriptCapitalO]1,\[ScriptCapitalO]2}];

]


(* ::Subsection::Closed:: *)
(*End*)


End[];

EndPackage[];
