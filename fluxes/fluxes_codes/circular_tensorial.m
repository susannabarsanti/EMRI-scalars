(* ::Package:: *)

Exit[];


<<SpinWeightedSpheroidalHarmonics`
<<KerrGeodesics`
<<Teukolsky`
ParallelNeeds["SpinWeightedSpheroidalHarmonics`"];
ParallelNeeds["KerrGeodesics`"];
ParallelNeeds["Teukolsky`"];


(*Function for the single mode flux computation.
- R: orbital radius;
- (l,m): modes indexes; 
- A: MBH spin*)


IntGrav[R_,\[ScriptL]_,mu_,a_]:=Block[{l=\[ScriptL],dLd\[Sigma]h,dLd\[Sigma]\[Infinity],dEd\[Sigma]\[Infinity],dEd\[Sigma]h,m=mu,orbit,mode,s,n,k,\[Omega]},

orbit= If[a>0, KerrGeoOrbit[a,R,0,1],  KerrGeoOrbit[-a,R,0,-1]];

s=-2;n=0;k=0;
mode=TeukolskyPointParticleMode[s,l,m,n,k,orbit];

Print["(l,m) = ","(",l,",",m,")"];
Print["R = ","(",SetPrecision[R,7],")"];

dEd\[Sigma]h = mode["Fluxes"]["Energy"]["\[ScriptCapitalH]"];
dEd\[Sigma]\[Infinity] = mode["Fluxes"]["Energy"]["\[ScriptCapitalI]"]; 

\[Omega] = If[a>0, KerrGeoFrequencies[a,R,0,1][[3]], KerrGeoFrequencies[-a,R,0,-1][[3]]];

dLd\[Sigma]h = dEd\[Sigma]h* m / \[Omega];
dLd\[Sigma]\[Infinity] = dEd\[Sigma]\[Infinity]* m / \[Omega];
Return[{dEd\[Sigma]h,dEd\[Sigma]\[Infinity],dLd\[Sigma]h,dLd\[Sigma]\[Infinity]}]
];


(*Radial grid. 
- nn: number of points; 
- eps: shift from the isco, the minimum radius will be ISCO+eps; 
- rfin: the maximum radius will be ISCO+rifn+eps *)

nn=3; 
eps = 2/100;
rfin = 30; 

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
- lmin: goes from 2 for the tensorial circular case. It can be 0 for non-circular orbits; 
- lmax: usually 15 for the circular case provides enough accuracy;
- a: MBH spin*)

tab[lmin_,lmax_,a_]:=Block[{nmult,multipoles,j,m,kk=1,x},

nmult=Length[Flatten[Table[Range[1,i],{i,lmin,lmax}]]];
multipoles=ConstantArray[0,nmult];

For[j=lmin,j<=lmax,j++,
	For[m=1,m<=j,m++,
	
multipoles[[kk]]={j,m};
kk = kk+1;
]
];

(*Select even multipoles for the scalar case: 
multipoles=Select[multipoles,EvenQ[#[[1]]+#[[2]]]&];*)

Print[Length[multipoles]];
Print[multipoles];

(*Factor 2 that accounts for the m<0 modes*)
x = Partition[Flatten[ParallelTable[{arrayp[a][[ics]],Sum[2*IntGrav[arrayp[a][[ics]],multipoles[[l,1]],multipoles[[l,2]],a],{l,1,Length[multipoles]}]},{ics,1,Length[arrayp[a]]}]],5];
Return[x];
];


(*Exporting the data in a .dat file*)
ellmin=2;
ellmax=15;
spin=0.9`30;
Export["circ_scal_aXX.dat",tab[ellmin,ellmax,spin]];


Exit[];
