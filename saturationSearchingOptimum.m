(* ::Package:: *)

#/home/song/bin/WolframScript -script

ClearAll["Global`*"]

SeedRandom[];

AbsoluteTiming[

ParallelTable[
Block[{totP,totS,kd,stepNum,sampleSize,totT,x,t,des,init,totK,pars,vars,dvars,k11},

SetDirectory["/home/song/workspace/networkMotifs/scaffoldTitration"];
kd=10;
totP=0.1;totS=1;
stepNum=5;
sampleSize=10000;

totK=0.0001;
init={totK,totP,totS,0,0,0,totT,0,0,1.*10^-4};
des={-k[1]* x[1][t] *x[3][t]+k[2]* x[5][t]+k[3] *x[5][t]-k[7]* x[1][t]* x[7][t]+k[8]*x[8][t]+k11[t]*kd-kd*x[1][t],-k[4]* x[2][t] *x[4][t]+k[5]* x[6][t]+k[6] *x[6][t]-k[9]* x[2][t]* x[7][t]+k[10]* x[9][t],-k[1]* x[1][t] *x[3][t]+k[2] *x[5][t]+k[6]*x[6][t],
-k[4]*x[2][t]* x[4][t]+k[3]*x[5][t]+k[5]* x[6][t],
k[1]* x[1][t] *x[3][t]-k[2] *x[5][t]-k[3]* x[5][t],
k[4] *x[2][t]* x[4][t]-k[5] *x[6][t]-k[6] *x[6][t],-k[7] *x[1][t]* x[7][t]-k[9] *x[2][t]* x[7][t]+k[8]* x[8][t]+k[10]*x[9][t],
k[7]* x[1][t] *x[7][t]-k[8]* x[8][t],
k[9] *x[2][t] *x[7][t]-k[10] *x[9][t],0};

vars=Array[x,9];AppendTo[vars,k11];
dvars=Thread[Derivative[1][vars]];

Block[{T,ssthreshold,sampleSizeT,count,allAd,allUs,maxUs,maxAd,ks,num},

pars={};
For[num=1,num<=sampleSize,num++,
(*tot[n_]:=tot[n]=10^(RandomReal[]*4-3);*)
(*ksTest1=Array[k,10];*)
(*totT=1.*^-3;*)
allAd={};
allUs={};
sampleSizeT=100;

Block[{k,tPer,step,ad,us,ts,x4,i,sol},k[n_]:=k[n]=10^(RandomReal[]*6-3);
step=0;
tPer={};
ssthreshold=1.*^-5;
For[count=1,count<=sampleSizeT,count++,
totT=1.*10^(RandomReal[]*4-3);
(* Print[des]; *){sol}=NDSolve[{Through[dvars[t]]==des,Through[vars[0]]==init,With[{df=Through[dvars[t]]},WhenEvent[Norm[df]<ssthreshold,{AppendTo[tPer,t],step=step+1,If[step>stepNum,"StopIntegration"],k11[t]->10*k11[t]}]]},vars,{t,0,200000},MaxSteps->10000];
ts=tPer;
If[Length[ts]==stepNum+1 &&AllTrue[ts,Positive],
x4=Evaluate[x[4][ts-0.001]/.sol];xT=Evaluate[(x[7][ts-0.001]+x[8][ts-0.001]+x[9][ts-0.001])/.sol];

us=Sqrt[((Abs[(x4[[4]]-x4[[3]])]/totS)*Min[((Abs[(x4[[4]]-x4[[3]])]/Max[Abs[(x4[[3]]-x4[[1]])],0.001]+Abs[(x4[[4]]-x4[[3]])]/Max[Abs[(x4[[stepNum+1]]-x4[[4]])],0.001])/2)/10.0,1.0])];

ad=0.0001;
For[i=1,i<=stepNum,i++,
ad=ad*Sqrt[(Min[(Max[Abs[Evaluate[x[4][Range[ts[[i]],ts[[i+1]],1]]/.sol]-Evaluate[x[4][ts[[i]]]/.sol]]]/(0.2*totS)),1.0]*((0.01)/(Max[Abs[(x4[[i+1]]-x4[[i]])/totS],0.01])))];
];
ad=(ad/0.0001)^(1/(stepNum));
AppendTo[allUs,us];
AppendTo[allAd,ad];
];
];
maxUs=Max[allUs];
maxAd=Max[allAd];
ks=Array[k,10];
AppendTo[pars,Join[ks,{totT,totK,totP,totS,maxUs,maxAd,num,(ks[[2]]+ks[[3]])/ks[[1]],(ks[[5]]+ks[[6]])/ks[[4]],ks[[8]]/ks[[7]],ks[[10]]/ks[[9]]}]];
];
];
Export["saturationSearchingOptimum"<>ToString[parallelNum]<>".csv",pars];
];
],{parallelNum,15}]
]

Exit[];



