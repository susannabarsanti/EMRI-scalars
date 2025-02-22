(* ::Package:: *)

BeginPackage["EmriInspiralLECC`"];


<<SpinWeightedSpheroidalHarmonics`
<<KerrGeodesics`
<<Teukolsky`


(* ::Subsection::Closed:: *)
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

integralN::usage = 
		"Numerical integral with different methods";
		
SnLISA::usage = 
		"Analytical fit of the LISA psd";

transition::usage =
		"Thorne - Ori transition radius";

dopplerf::usage =
		"Doppler frequnecy";	
		
dEddrcirc::usage = 
		"circular scalar flux";			
		
dEGRdrcirc::usage = 
		"circular gravitational flux";

dEddr::usage = 
		"scalar flux - [e,p]";
		
dEGRdr::usage = 
		"grav eccentric flux - [e,p]";

		
dLddr::usage = 
		"scalar angular momentum flux - [e,p]";
		
dLGRdr::usage = 
		"grav angular momentum flux - [e,p]";			


s0=Import["/Users/susannabarsanti/Documents/EMRI_ECC/DATA/scalar/a09/circ_scal_a09.dat"];
s1=Import["/Users/susannabarsanti/Documents/EMRI_ECC/DATA/scalar/a09/e01.dat"];
s2=Import["/Users/susannabarsanti/Documents/EMRI_ECC/DATA/scalar/a09/e02.dat"];
s3=Import["/Users/susannabarsanti/Documents/EMRI_ECC/DATA/scalar/a09/e03.dat"];
s4=Import["/Users/susannabarsanti/Documents/EMRI_ECC/DATA/scalar/a09/e04.dat"];
s5=Import["/Users/susannabarsanti/Documents/EMRI_ECC/DATA/scalar/a09/e05.dat"];


t0=Import["/Users/susannabarsanti/Documents/EMRI_ECC/DATA/tensorial/a09/circ_grav_a09.dat"];
t1=Import["/Users/susannabarsanti/Documents/EMRI_ECC/DATA/tensorial/a09/e01.dat"];
t2=Import["/Users/susannabarsanti/Documents/EMRI_ECC/DATA/tensorial/a09/e02.dat"];
t3=Import["/Users/susannabarsanti/Documents/EMRI_ECC/DATA/tensorial/a09/e03.dat"];
t4=Import["/Users/susannabarsanti/Documents/EMRI_ECC/DATA/tensorial/a09/e04.dat"];
t5=Import["/Users/susannabarsanti/Documents/EMRI_ECC/DATA/tensorial/a09/e05.dat"];


arrayup[ee_,a_]:=Range[1/Sqrt[0.1`200*KerrGeoSeparatrix[a,ee,1]+10.02`200],1/Sqrt[0.1`200*KerrGeoSeparatrix[a,ee,1]+0.02`200],(1/Sqrt[0.1`200*KerrGeoSeparatrix[a,ee,1]+0.02`200]-1/Sqrt[0.1`200*KerrGeoSeparatrix[a,ee,1]+10.02`200])/40];
arrayp[ee_,a_]:=Table[(1/arrayup[ee,a][[i]])^2+0.9`180*KerrGeoSeparatrix[a,ee,1],{i,1,Length[arrayup[ee,a]]}];


setab0=Table[{0,x-KerrGeoSeparatrix[0.9`200,0,1],Interpolation[Thread[{s0[[All,1]],s0[[All,2]]+s0[[All,3]]}]][x]},{x,arrayp[0,0.9`120][[41]],arrayp[0,0.9`150][[1]],(arrayp[0,0.9`91][[1]]- arrayp[0,0.9`190][[41]])/200}];
setab01=Table[{0.1`100,x-KerrGeoSeparatrix[0.9`200,0.1`200,1],Interpolation[Thread[{s1[[All,1]],s1[[All,2]]+s1[[All,3]]}]][x]},{x,arrayp[0.1`120,0.9`120][[41]],arrayp[0.1`120,0.9`120][[1]],(arrayp[0.1`185,0.9`185][[1]]- arrayp[0.1`185,0.9`185][[41]])/200}];
setab02=Table[{0.2`100,x-KerrGeoSeparatrix[0.9`200,0.2`200,1],Interpolation[Thread[{s2[[All,1]],s2[[All,2]]+s2[[All,3]]}]][x]},{x,arrayp[0.2`165,0.9`165][[41]],arrayp[0.2`180,0.9`180][[1]],(arrayp[0.2`185,0.9`185][[1]]- arrayp[0.2`185,0.9`185][[41]])/200}];
setab03=Table[{0.3`100,x-KerrGeoSeparatrix[0.9`200,0.3`200,1],Interpolation[Thread[{s3[[All,1]],s3[[All,2]]+s3[[All,3]]}]][x]},{x,arrayp[0.3`150,0.9`150][[41]],arrayp[0.3`100,0.9`100][[1]],(arrayp[0.3`185,0.9`185][[1]]- arrayp[0.3`185,0.9`185][[41]])/200}];
setab04=Table[{0.4`100,x-KerrGeoSeparatrix[0.9`200,0.4`200,1],Interpolation[Thread[{s4[[All,1]],s4[[All,2]]+s4[[All,3]]}]][x]},{x,arrayp[0.4`150,0.9`150][[41]],arrayp[0.4`180,0.9`180][[1]],(arrayp[0.4`185,0.9`185][[1]]- arrayp[0.4`185,0.9`185][[41]])/200}];
setab05=Table[{0.5`100,x-KerrGeoSeparatrix[0.9`200,0.5`200,1],Interpolation[Thread[{s5[[All,1]],s5[[All,2]]+s5[[All,3]]}]][x]},{x,arrayp[0.5`91,0.9`90][[41]],arrayp[0.5`110,0.9`110][[1]],(arrayp[0.5`105,0.9`105][[1]]- arrayp[0.5`92,0.9`92][[41]])/200}];


sltab0=Table[{0,x-KerrGeoSeparatrix[0.9`200,0,1],Interpolation[Thread[{s0[[All,1]],s0[[All,4]]+s0[[All,5]]}]][x]},{x,arrayp[0,0.9`120][[41]],arrayp[0,0.9`150][[1]],(arrayp[0,0.9`91][[1]]- arrayp[0,0.9`190][[41]])/200}];
sltab01=Table[{0.1`100,x-KerrGeoSeparatrix[0.9`200,0.1`200,1],Interpolation[Thread[{s1[[All,1]],s1[[All,4]]+s1[[All,5]]}]][x]},{x,arrayp[0.1`120,0.9`120][[41]],arrayp[0.1`120,0.9`120][[1]],(arrayp[0.1`185,0.9`185][[1]]- arrayp[0.1`185,0.9`185][[41]])/200}];
sltab02=Table[{0.2`100,x-KerrGeoSeparatrix[0.9`200,0.2`200,1],Interpolation[Thread[{s2[[All,1]],s2[[All,4]]+s2[[All,5]]}]][x]},{x,arrayp[0.2`165,0.9`165][[41]],arrayp[0.2`180,0.9`180][[1]],(arrayp[0.2`185,0.9`185][[1]]- arrayp[0.2`185,0.9`185][[41]])/200}];
sltab03=Table[{0.3`100,x-KerrGeoSeparatrix[0.9`200,0.3`200,1],Interpolation[Thread[{s3[[All,1]],s3[[All,4]]+s3[[All,5]]}]][x]},{x,arrayp[0.3`150,0.9`150][[41]],arrayp[0.3`100,0.9`100][[1]],(arrayp[0.3`185,0.9`185][[1]]- arrayp[0.3`185,0.9`185][[41]])/200}];
sltab04=Table[{0.4`100,x-KerrGeoSeparatrix[0.9`200,0.4`200,1],Interpolation[Thread[{s4[[All,1]],s4[[All,4]]+s4[[All,5]]}]][x]},{x,arrayp[0.4`150,0.9`150][[41]],arrayp[0.4`180,0.9`180][[1]],(arrayp[0.4`185,0.9`185][[1]]- arrayp[0.4`185,0.9`185][[41]])/200}];
sltab05=Table[{0.5`100,x-KerrGeoSeparatrix[0.9`200,0.5`200,1],Interpolation[Thread[{s5[[All,1]],s5[[All,4]]+s5[[All,5]]}]][x]},{x,arrayp[0.5`91,0.9`90][[41]],arrayp[0.5`110,0.9`110][[1]],(arrayp[0.5`105,0.9`105][[1]]- arrayp[0.5`92,0.9`92][[41]])/200}];


tetab0=Table[{0,x-KerrGeoSeparatrix[0.9`200,0,1],Interpolation[Thread[{t0[[All,1]],t0[[All,2]]+t0[[All,3]]}]][x]},{x,arrayp[0,0.9`120][[41]],arrayp[0,0.9`150][[1]],(arrayp[0,0.9`91][[1]]- arrayp[0,0.9`190][[41]])/200}];
tetab01=Table[{0.1`100,x-KerrGeoSeparatrix[0.9`200,0.1`200,1],Interpolation[Thread[{t1[[All,1]],t1[[All,2]]+t1[[All,3]]}]][x]},{x,arrayp[0.1`120,0.9`120][[41]],arrayp[0.1`120,0.9`120][[1]],(arrayp[0.1`185,0.9`185][[1]]- arrayp[0.1`185,0.9`185][[41]])/200}];
tetab02=Table[{0.2`100,x-KerrGeoSeparatrix[0.9`200,0.2`200,1],Interpolation[Thread[{t2[[All,1]],t2[[All,2]]+t2[[All,3]]}]][x]},{x,arrayp[0.2`165,0.9`165][[41]],arrayp[0.2`180,0.9`180][[1]],(arrayp[0.2`185,0.9`185][[1]]- arrayp[0.2`185,0.9`185][[41]])/200}];
tetab03=Table[{0.3`100,x-KerrGeoSeparatrix[0.9`200,0.3`200,1],Interpolation[Thread[{t3[[All,1]],t3[[All,2]]+t3[[All,3]]}]][x]},{x,arrayp[0.3`150,0.9`150][[41]],arrayp[0.3`100,0.9`100][[1]],(arrayp[0.3`185,0.9`185][[1]]- arrayp[0.3`185,0.9`185][[41]])/200}];
tetab04=Table[{0.4`100,x-KerrGeoSeparatrix[0.9`200,0.4`200,1],Interpolation[Thread[{t4[[All,1]],t4[[All,2]]+t4[[All,3]]}]][x]},{x,arrayp[0.4`150,0.9`150][[41]],arrayp[0.4`180,0.9`180][[1]],(arrayp[0.4`185,0.9`185][[1]]- arrayp[0.4`185,0.9`185][[41]])/200}];
tetab05=Table[{0.5`100,x-KerrGeoSeparatrix[0.9`200,0.5`200,1],Interpolation[Thread[{t5[[All,1]],t5[[All,2]]+t5[[All,3]]}]][x]},{x,arrayp[0.5`91,0.9`90][[41]],arrayp[0.5`110,0.9`110][[1]],(arrayp[0.5`105,0.9`105][[1]]- arrayp[0.5`92,0.9`92][[41]])/200}];


tltab0=Table[{0,x-KerrGeoSeparatrix[0.9`200,0,1],Interpolation[Thread[{t0[[All,1]],t0[[All,4]]+t0[[All,5]]}]][x]},{x,arrayp[0,0.9`120][[41]],arrayp[0,0.9`150][[1]],(arrayp[0,0.9`91][[1]]- arrayp[0,0.9`190][[41]])/200}];
tltab01=Table[{0.1`100,x-KerrGeoSeparatrix[0.9`200,0.1`200,1],Interpolation[Thread[{t1[[All,1]],t1[[All,4]]+t1[[All,5]]}]][x]},{x,arrayp[0.1`120,0.9`120][[41]],arrayp[0.1`120,0.9`120][[1]],(arrayp[0.1`185,0.9`185][[1]]- arrayp[0.1`185,0.9`185][[41]])/200}];
tltab02=Table[{0.2`100,x-KerrGeoSeparatrix[0.9`200,0.2`200,1],Interpolation[Thread[{t2[[All,1]],t2[[All,4]]+t2[[All,5]]}]][x]},{x,arrayp[0.2`165,0.9`165][[41]],arrayp[0.2`180,0.9`180][[1]],(arrayp[0.2`185,0.9`185][[1]]- arrayp[0.2`185,0.9`185][[41]])/200}];
tltab03=Table[{0.3`100,x-KerrGeoSeparatrix[0.9`200,0.3`200,1],Interpolation[Thread[{t3[[All,1]],t3[[All,4]]+t3[[All,5]]}]][x]},{x,arrayp[0.3`150,0.9`150][[41]],arrayp[0.3`100,0.9`100][[1]],(arrayp[0.3`185,0.9`185][[1]]- arrayp[0.3`185,0.9`185][[41]])/200}];
tltab04=Table[{0.4`100,x-KerrGeoSeparatrix[0.9`200,0.4`200,1],Interpolation[Thread[{t4[[All,1]],t4[[All,4]]+t4[[All,5]]}]][x]},{x,arrayp[0.4`150,0.9`150][[41]],arrayp[0.4`180,0.9`180][[1]],(arrayp[0.4`185,0.9`185][[1]]- arrayp[0.4`185,0.9`185][[41]])/200}];
tltab05=Table[{0.5`100,x-KerrGeoSeparatrix[0.9`200,0.5`200,1],Interpolation[Thread[{t5[[All,1]],t5[[All,4]]+t5[[All,5]]}]][x]},{x,arrayp[0.5`91,0.9`90][[41]],arrayp[0.5`110,0.9`110][[1]],(arrayp[0.5`105,0.9`105][[1]]- arrayp[0.5`92,0.9`92][[41]])/200}];


settab0=Table[{{setab0[[i,1]],setab0[[i,2]]},setab0[[i,3]]},{i,1,Length[setab0[[All,1]]]}];
settab1=Table[{{setab01[[i,1]],setab0[[i,2]]},setab01[[i,3]]},{i,1,Length[setab01[[All,1]]]}];
settab2=Table[{{setab02[[i,1]],setab0[[i,2]]},setab02[[i,3]]},{i,1,Length[setab02[[All,1]]]}];
settab3=Table[{{setab03[[i,1]],setab0[[i,2]]},setab03[[i,3]]},{i,1,Length[setab03[[All,1]]]}];
settab4=Table[{{setab04[[i,1]],setab0[[i,2]]},setab04[[i,3]]},{i,1,Length[setab04[[All,1]]]}];
settab5=Table[{{setab05[[i,1]],setab0[[i,2]]},setab05[[i,3]]},{i,1,Length[setab05[[All,1]]]}];


slttab0=Table[{{sltab0[[i,1]],sltab0[[i,2]]},sltab0[[i,3]]},{i,1,Length[sltab0[[All,1]]]}];
slttab1=Table[{{sltab01[[i,1]],sltab0[[i,2]]},sltab01[[i,3]]},{i,1,Length[sltab01[[All,1]]]}];
slttab2=Table[{{sltab02[[i,1]],sltab0[[i,2]]},sltab02[[i,3]]},{i,1,Length[sltab02[[All,1]]]}];
slttab3=Table[{{sltab03[[i,1]],sltab0[[i,2]]},sltab03[[i,3]]},{i,1,Length[sltab03[[All,1]]]}];
slttab4=Table[{{sltab04[[i,1]],sltab0[[i,2]]},sltab04[[i,3]]},{i,1,Length[sltab04[[All,1]]]}];
slttab5=Table[{{sltab05[[i,1]],sltab0[[i,2]]},sltab05[[i,3]]},{i,1,Length[sltab05[[All,1]]]}];


tettab0=Table[{{tetab0[[i,1]],tetab0[[i,2]]},tetab0[[i,3]]},{i,1,Length[tetab0[[All,1]]]}];
tettab1=Table[{{tetab01[[i,1]],tetab0[[i,2]]},tetab01[[i,3]]},{i,1,Length[tetab01[[All,1]]]}];
tettab2=Table[{{tetab02[[i,1]],tetab0[[i,2]]},tetab02[[i,3]]},{i,1,Length[tetab02[[All,1]]]}];
tettab3=Table[{{tetab03[[i,1]],tetab0[[i,2]]},tetab03[[i,3]]},{i,1,Length[tetab03[[All,1]]]}];
tettab4=Table[{{tetab04[[i,1]],tetab0[[i,2]]},tetab04[[i,3]]},{i,1,Length[tetab04[[All,1]]]}];
tettab5=Table[{{tetab05[[i,1]],tetab0[[i,2]]},tetab05[[i,3]]},{i,1,Length[tetab05[[All,1]]]}];


tlttab0=Table[{{tltab0[[i,1]],tltab0[[i,2]]},tltab0[[i,3]]},{i,1,Length[tltab0[[All,1]]]}];
tlttab1=Table[{{tltab01[[i,1]],tltab0[[i,2]]},tltab01[[i,3]]},{i,1,Length[tltab01[[All,1]]]}];
tlttab2=Table[{{tltab02[[i,1]],tltab0[[i,2]]},tltab02[[i,3]]},{i,1,Length[tltab02[[All,1]]]}];
tlttab3=Table[{{tltab03[[i,1]],tltab0[[i,2]]},tltab03[[i,3]]},{i,1,Length[tltab03[[All,1]]]}];
tlttab4=Table[{{tltab04[[i,1]],tltab0[[i,2]]},tltab04[[i,3]]},{i,1,Length[tltab04[[All,1]]]}];
tlttab5=Table[{{tltab05[[i,1]],tltab0[[i,2]]},tltab05[[i,3]]},{i,1,Length[tltab05[[All,1]]]}];


sse=Join[settab0,settab1,settab2,settab3,settab4,settab5];
ssl=Join[slttab0,slttab1,slttab2,slttab3,slttab4,slttab5];
tte=Join[tettab0,tettab1,tettab2,tettab3,tettab4,tettab5];
ttl=Join[tlttab0,tlttab1,tlttab2,tlttab3,tlttab4,tlttab5];


dEddr[e_,p_]= Interpolation[sse][e,p-KerrGeoSeparatrix[0.9`300,e,1]];
dLddr[e_,p_]= Interpolation[ssl][e,p-KerrGeoSeparatrix[0.9`300,e,1]];
dEGRdr[e_,p_]= Interpolation[tte][e,p-KerrGeoSeparatrix[0.9`300,e,1]];
dLGRdr[e_,p_]= Interpolation[ttl][e,p-KerrGeoSeparatrix[0.9`300,e,1]];


rplcirc = SetPrecision[KerrGeoSeparatrix[0.9`150,0,1]+0.1`150,100];
dEGRdrcirc[r_]:= Interpolation[Table[{x,Interpolation[Thread[{t0[[All,1]],t0[[All,2]]+t0[[All,3]]}]][x]},{x,rplcirc,12.34`100,(12.34`100- rplcirc)/200}]][r];
dEddrcirc[r_]:= Interpolation[Table[{x,Interpolation[Thread[{s0[[All,1]],s0[[All,2]]+s0[[All,3]]}]][x]},{x,rplcirc,12.34`100,(12.34`100- rplcirc)/200}]][r];


Begin["`Private`"];


precision = 45;
year = 365;
day = 3600 24;
GN = 6674 10^-14;
cc = 299792458 10^-3;
msun = 1477 10^-3;
Rs = 1495978707/10;
pc = 96939420213600/\[Pi];


(* ::Subsection::Closed:: *)
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


findr0[m1_,m2_,spin_,\[Alpha]_,eend_,T_]:=Block[{q=m2/m1,factsec,pend,factday,fact,IC,tend,ini,e,p,\[ScriptCapitalE],\[ScriptCapitalL],eqns,t,dEde,dEdp,H,dLde,dLdp},

\[ScriptCapitalE][p_,e_] = dEGRdr[e,p]+\[Alpha]^2 (dEddr[e,p])/4;
\[ScriptCapitalL][p_,e_] = dLGRdr[e,p]+\[Alpha]^2 (dLddr[e,p])/4;
 		
factsec = cc/(m1 msun);
factday = factsec day ; (* change time in days in dimensionless time *)
fact = factsec;

tend = T day;
(* Stop the integration when it's at the transition T = 0 *)

dEde[p_,e_]= D[KerrGeoEnergy[spin,p,e,1],e];
dEdp[p_,e_]= D[KerrGeoEnergy[spin,p,e,1],p]; 
dLde[p_,e_]= D[KerrGeoAngularMomentum[spin,p,e,1],e]; 
dLdp[p_,e_]= D[KerrGeoAngularMomentum[spin,p,e,1],p]; 

H = dEdp[p[t],e[t]] dLde[p[t],e[t]] -dEde[p[t],e[t]] dLdp[p[t],e[t]]; 

eqns = {p'[t]== -q fact H^-1 (- dEde[p[t],e[t]] \[ScriptCapitalL][p[t],e[t]] + dLde[p[t],e[t]] \[ScriptCapitalE][p[t],e[t]]),
e'[t]== -q fact H^-1 (dEdp[p[t],e[t]] \[ScriptCapitalL][p[t],e[t]] - dLdp[p[t],e[t]] \[ScriptCapitalE][p[t],e[t]])};

pend = KerrGeoSeparatrix[0.9`350,eend,1]+0.11`350;


IC = {p[tend]==pend, e[tend]== eend};

ini = NDSolveValue[{eqns,IC},{p,e},{t,tend,0},
				    Method->"StiffnessSwitching",WorkingPrecision->45,MaxSteps-> \[Infinity]];

Return[ini];

]


(* ::Subsubsection::Closed:: *)
(*Adiabatic evolution for a given T*)


OrbitalEvoT[m1_,m2_,spin_,\[Alpha]_,p0_,e0_,\[Phi]0_,T_]:=Block[{phiR,en,pr,L,r,r0,Vr,V\[Phi],Vt,J,d\[Chi]dt,x,\[Chi],solr,t,sol,p,e,\[Psi],factday,factsec,tend,q=m2/m1,fact,eqns,IC,\[ScriptCapitalE],\[ScriptCapitalL],dEde,dEdp,dLde,dLdp,H,d\[Phi]dt},

\[ScriptCapitalE][p_,e_] = dEGRdr[e,p]+\[Alpha]^2 (dEddr[e,p])/4;
\[ScriptCapitalL][p_,e_] = dLGRdr[e,p]+\[Alpha]^2 (dLddr[e,p])/4;

factsec = cc/(m1 msun);
factday = factsec day ; (* change time in days in dimensionless time *)
fact = factsec;

tend = T day;

dEde[p_,e_]= D[KerrGeoEnergy[spin,p,e,1],e];
dEdp[p_,e_]= D[KerrGeoEnergy[spin,p,e,1],p]; 
dLde[p_,e_]= D[KerrGeoAngularMomentum[spin,p,e,1],e]; 
dLdp[p_,e_]= D[KerrGeoAngularMomentum[spin,p,e,1],p]; 

H = dEdp[p[t],e[t]] dLde[p[t],e[t]] - dEde[p[t],e[t]] dLdp[p[t],e[t]]; 

en[a_,p_,e_]:= KerrGeoEnergy[a,p,e,1];
L[a_,p_,e_]:= KerrGeoAngularMomentum[a,p,e,1];

eqns = {p'[t]==   -q fact H^-1 (- dEde[p[t],e[t]] \[ScriptCapitalL][p[t],e[t]] + dLde[p[t],e[t]] \[ScriptCapitalE][p[t],e[t]]),
e'[t]== - q fact H^-1 (dEdp[p[t],e[t]] \[ScriptCapitalL][p[t],e[t]] - dLdp[p[t],e[t]] \[ScriptCapitalE][p[t],e[t]]),
			\[Psi]'[t]== fact KerrGeoFrequencies[spin,p[t],e[t],1][[3]], phiR'[t]== fact KerrGeoFrequencies[spin,p[t],e[t],1][[1]]};

(*p0 = (2 rp ra)/(ra+rp) ; 
Print["p0=",N[p0,6]];
e0 = ( ra - rp)/(ra+rp) ; 
Print["e0=",N[e0,6]];*)

IC = {p[0]==p0 ,e[0]==e0, \[Psi][0]==0, phiR[0]==0};


sol=NDSolveValue[{eqns,IC},{\[Psi],p,e,phiR},{t,0,tend},
				    Method->"StiffnessSwitching",WorkingPrecision->40,MaxSteps-> \[Infinity]];

Return[{sol[[1]],sol[[2]],sol[[3]],sol[[4]]}];

]


(* ::Subsection::Closed:: *)
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


dopplerf[m1_,m2_,spin_,p0_,e0_,\[Phi]0_,\[Theta]s_,\[Phi]s_,T_,t_,d_]:=Block[{M,\[Mu],fdopp,\[Omega],x,evolution,\[Omega]D,fD},

M = m1 msun; (* total mass *)
\[Mu] = m2 msun; (* reduced mass *)

evolution = OrbitalEvoT[m1,m2,spin,d,p0,e0,\[Phi]0,T];


\[Omega]D = D[\[Omega][x] Rs/cc Sin[\[Theta]s]Cos[2\[Pi] x/(day year)-\[Phi]s],x]; (* doppler angular frequency *)

fD = (\[Omega]D/\[Pi])//.{\[Omega][x]->evolution[[1]]'[t],\[Omega]'[x]->evolution[[1]]''[t]};

Return[fD//.x->t];

]


(* ::Subsubsection::Closed:: *)
(*h+x polarizations*)


emritemplate[m1_,m2_,spin_,p0_,e0_,\[Phi]0_,\[Theta]s_,\[Phi]s_,\[Theta]k_,\[Phi]k_,T_,t_,d_,nmax_]:=Block[{hPPp,hCCc,hPlus,hCross,A,n,a,b,c,evolution,\[Phi],\[Phi]R,\[Mu],M,\[Phi]D,p,e,\[Omega],LN},

M = m1 msun; (* primary mass *)
\[Mu] = m2 msun; (* secondary mass *)

evolution = OrbitalEvoT[m1,m2,spin,d,p0,e0,\[Phi]0,T];
Print["evo_temp"];

\[Omega] = evolution[[1]]'[t]; (* angular frequency *)

\[Phi]D = \[Omega] Rs/cc Sin[\[Theta]s]Cos[2\[Pi] t/(day year)-\[Phi]s]; (* doppler phase *)

\[Phi] = evolution[[1]][t]+\[Phi]D; (* total phase *)
p = evolution[[2]][t];
e = evolution[[3]][t];
\[Phi]R = evolution[[4]][t];

LN = Cos[\[Theta]k]Cos[\[Theta]s]+Sin[\[Theta]k]Sin[\[Theta]s]Cos[\[Phi]k-\[Phi]s];

a = ConstantArray[0,nmax];
b = ConstantArray[0,nmax];
c = ConstantArray[0,nmax];

A = (M \[Omega]/cc)^(2/3) \[Mu];

For[n=1,n<=nmax,n++,
a[[n]]= -n A (BesselJ[n-2,n e] - 2 e BesselJ[n-1,n e]+2/n BesselJ[n,n e]+ 2 e  BesselJ[n+1,n e] - BesselJ[n+2,n e] ) Cos[n (\[Phi]0+\[Phi]D+\[Phi])];

b[[n]]= -n A (1-e^2)^(1/2) (BesselJ[n-2,n e] - 2  BesselJ[n,n e]+BesselJ[n+2,n e]) Sin[n (\[Phi]0+\[Phi]D+\[Phi])] ;

c[[n]]= 2 A BesselJ[n,n e] Cos[n (\[Phi]0+\[Phi]D+\[Phi])]  ;
];


hPlus = ConstantArray[0,nmax];
hCross = ConstantArray[0,nmax];
Clear[n];

Print["hCross",hCross];

For[n=1,n<=nmax,n++,
hPlus[[n]]=  ((1+LN^2)(a[[n]] Cos[ArcCos[Sqrt[2]Cos[\[Phi]R]]]-b[[n]] Sin[ArcCos[Sqrt[2]Cos[\[Phi]R]]])+(1-LN^2) c[[n]]);

hCross[[n]]=-   2 LN (b[[n]] Cos[ArcCos[Sqrt[2]Cos[\[Phi]R]]]+ a[[n]] Sin[ArcCos[Sqrt[2]Cos[\[Phi]R]]]);
];


hPPp = Sum[hPlus[[n]],{n,1,nmax}];
hCCc = Sum[hCross[[n]],{n,1,nmax}];

Print["hplus"];
Print["hCC"];

Return[{hPPp,hCCc}];

]


(* ::Subsubsection::Closed:: *)
(*Sampling in time*)


htime[m1_,m2_,spin_,p0_,e0_,\[Theta]s_,\[Phi]s_,\[Phi]0_,\[Theta]k_,\[Phi]k_,observingtime_,d_,nmax_]:=Block[{win,nn,timevec,T,\[Delta]t,t,
i,fmax,h,fdopp,hp,hx,hplus,hcross,pfin,efin,psep,FplusI,
FcrossI,FplusII,FcrossII,risco,fsampling,fisco,evolution},

T = day observingtime; (* total time in seconds *)

evolution = OrbitalEvoT[m1,m2,spin,d,p0,e0,\[Phi]0,observingtime];
Print["evo_htim"];

pfin = evolution[[2]][T]; (* final radius after T *)	
efin = evolution[[3]][T];
psep = KerrGeoSeparatrix[spin,efin,1]; (* isco radius/M *)
fdopp= - dopplerf[m1,m2,spin,p0,e0,\[Phi]0,\[Theta]s,\[Phi]s,observingtime,T,d]; (* max fdop @ rmax *)



fmax = cc KerrGeoFrequencies[0.9`300,pfin,efin,1][[3]]/(m1 msun \[Pi])+fdopp; (* f @ rmax + fdop @ rmax *)
fisco = cc KerrGeoFrequencies[0.9`300,pfin,efin,1][[3]]/(m1 msun \[Pi]); (* Isco Kerr *)

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

{hplus[t_],hcross[t_]} = emritemplate[m1,m2,spin,p0,e0,\[Phi]0,\[Theta]s,\[Phi]s,\[Theta]k,\[Phi]k,observingtime,t,d,nmax]; (* analytical waveform *)

win = Table[window[i,nn,5/100],{i,0,nn-1}];

{hp,hx} = {win*hplus[timevec],win*hcross[timevec]};

If[hp[[1]]!=hp[[-1]] \[Or] hx[[1]]!=hx[[-1]],Print["Initial/final points of h are different"]];

{FplusI,FcrossI} = pattern[1,\[Theta]s,\[Phi]s,\[Theta]k,\[Phi]k,timevec];
{FplusII,FcrossII} = pattern[2,\[Theta]s,\[Phi]s,\[Theta]k,\[Phi]k,timevec];

Print["F+/Fx"];

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
