>
## //Parameters  
verbosity = 0;  
real D = 0.1;  
real H = 0.41;  
real cx0 = 0.2;  
real cy0 = 0.2; //center of cylinder  
real xa = 0.15;  
real ya = 0.2;  
real xe = 0.25;  
real ye = 0.2;  
int nn = 20;  
//TODO  
real Um = 1.5; //max velocity (Rey 100)  
real nu = 1e-3;  
//  
func U1 = 4.*Um*y*(H-y)/(H*H); //Boundary condition  
func U2 = 0.;  
real T=5;  
real dt = D/nn/Um; //CFL = 1  
real epspq = 1e-10;  
real eps = 1e-6;  
// Variables  
func Ub = Um*2./3.;  
real alpha = 1/dt;  
real Rey = Ub*D/nu;  
real t = 0.;  
//  
## // Mesh  
//border fr1(t=0, 2.2){y=0.1*sin(10*pi*t/2.2); x=t;   label=1;} // #swaving border  
//border fr2(t=0, H){x=2.2; y=t; label=2;}  
//border fr3(t=2.2, 0){x=t; y=H; label=1;}  
//border fr4(t=H, 0){x=0; y=t; label=22;}  
////border fr5(t=2*pi, 0){x=cx0+D*sin(t)/2;   y=cy0+D*cos(t)/2; label=3;}  
//mesh Th = buildmesh(fr1(5*nn) + fr2(nn) + fr3(5*nn) + fr4(nn));// + fr5(-nn*3)); // #mesh generation  
// part 
real d1=0.25; //(1) thickness
real d2=0.5;  //(3) thickness
real d3=0.25;  //(4) thickness
### //(A)  
border A1(t=0, 2.2){y=H; x=t; label=11;} 
border A2(t=0, d1){y=H+t; x=2.2; label=12;} 
border A3(t=2.2,0){y=H+d1; x=t; label=13;} // top  
border A4(t=d1,0){y=H+t; x=0; label=14;}   
mesh Th1 = buildmesh(A1(30)+A2(10)+A3(30)+A4(10));  
//plot(Th1,wait=1)  
### //(fr)    
border fr1(t=0, 2.2){y=0.1*sin(10*pi*t/2.2); x=t;label=1;} // #swaving border  
border fr2(t=0, H){x=2.2; y=t; label=2;}  
border fr3(t=2.2, 0){x=t; y=H; label=1;}  
border fr4(t=H, 0){x=0; y=t; label=22;}   
//border fr5(t=2*pi, 0){x=cx0+D*sin(t)/2;    y=cy0+D*cos(t)/2; label=3;}    
mesh Th2 = buildmesh(fr1(5*nn) + fr2(nn) + fr3(5*nn) + fr4(nn));// + fr5(-nn*3)); // #mesh generation  
//plot(Th2,wait=1);  
### //(B)  
border B1(t=2.2,0){y=0.1*sin(10*pi*t/2.2); x=t; label=31;}   
border B2(t=0, -d2){y=t; x=0; label=32;}  
border B3(t=0,2.2){y=-d2; x=t; label=33;}  
border B4(t=-d2,0){y=t; x=2.2; label=34;}  
mesh Th3 = buildmesh(B1(40)+B2(20)+B3(40)+B4(20));  
### //(C)  
border C1(t=2.2,0){y=-d2; x=t; label=41;}  
border C2(t=0, -d3){y=-d2+t; x=0; label=42;}  
border C3(t=0,2.2){y=-d2-d3; x=t; label=43;} // bottom  
border C4(t=d3,0){y=-d2-t; x=2.2; label=44;}  
mesh Th4 = buildmesh(C1(30)+C2(10)+C3(30)+C4(10));  
//  
mesh Th=Th1+Th2+Th3+Th4;  
plot(Th,wait=1);  
### // region labels
### // Th1=A,Th2=fr,Th3=B,Th4=C
int regTh1   = Th1(1.1, H+d1/2.0).region;        
int regTh2 = Th2(1.1, H/2.).region;   
int regTh3  = Th3(1.1, -d2/2.0).region;       
int regTh4  = Th4(1.1, -d2-d3/2.).region;  
## // thermal conductivity kappa     
fespace UUh(Th,P0);  
UUh kappa=1.;  
int n = 0;  
//real kappa=1;
### // A  
    for (int i = 0; i < Th1.nt; i++) {  
    int reg = Th[i].region;  
     if (reg == regTh1) {  
        kappa[][i] = 1.1;  
    } 
   }
### // fr  
    for (int i = Th1.nt; i < Th2.nt+Th1.nt; i++) {  
    int reg = Th[i].region;  
     if (reg == regTh2) {  
        kappa[][i] = 1.;  
    }  
   }  
### // B  
    for (int i = Th1.nt+Th2.nt; i < Th1.nt+Th2.nt+Th3.nt; i++) {  
    int reg = Th[i].region;  
     if (reg == regTh3) {  
        kappa[][i] = 1.2;  
    }  
   }  
### // C 
    for (int i = Th1.nt+Th2.nt+Th3.nt; i < Th1.nt+Th2.nt+Th3.nt+Th4.nt; i++) {  
    int reg = Th[i].region;  
     if (reg == regTh4) {  
        kappa[][i] = 1.1;  
    }  
   }  
plot(kappa,fill=true,wait=1);  
//  
### // Navier-Stokes eq
// Fespace  
fespace Mh(Th2, [P1]);  
Mh p;  
//
fespace Xh(Th2, [P2]);  
Xh u1, u2;  
//
fespace Wh(Th2, [P1dc]);  
Wh w; //vorticity  
//
// Macro  
macro grad(u) [dx(u), dy(u)] //  
macro div(u1, u2) (dx(u1) + dy(u2)) //  
//
// Problem  
varf von1 ([u1, u2, p], [v1, v2, q])  
    = //on(3, u1=0, u2=0)  // #comment out  
     on(1, u1=0, u2=0)  
     +on(2,22,u1=1,u2=0)  
    ;  
//
//remark : the value 100 in next varf is manualy   fitted, because free outlet.
### // Solver by variational form
varf vA (p, q) =  
    int2d(Th)(  
          grad(p)' * grad(q)  
    )  
    + int1d(Th, 2)(  
          100*p*q  
    )  
    ;  
//
varf vM (p, q)  
    = int2d(Th, qft=qf2pT)(  
          p*q  
    )  
    + on(2, p=0)  
    ;  
//
varf vu ([u1], [v1])  
    = int2d(Th)(  
          alpha*(u1*v1)  
        + nu*(grad(u1)' * grad(v1))  
    )  
    + on(1, u1=0)  
    ;  
//
varf vu1 ([p], [v1]) = int2d(Th)(p*dx(v1));  
varf vu2 ([p], [v1]) = int2d(Th)(p*dy(v1));   
//
varf vonu1 ([u1], [v1]) = on(1, u1=U1) + on(3, u1=0);  
varf vonu2 ([u1], [v1]) = on(1, u1=U2) + on(3, u1=0);  
//
matrix pAM = vM(Mh, Mh, solver=UMFPACK);  
matrix pAA = vA(Mh, Mh, solver=UMFPACK);  
matrix AU = vu(Xh, Xh, solver=UMFPACK);  
matrix vB1 = vu1(Mh, Xh);  
matrix vB2 = vu2(Mh, Xh);  
//
real[int] brhs1 = vonu1(0, Xh);
real[int] brhs2 = vonu2(0, Xh);
//
varf vrhs1(uu, vv) = int2d(Th)(convect([u1, u2], -dt, u1)*vv*alpha) + vonu1;  
varf vrhs2(v2, v1) = int2d(Th)(convect([u1, u2], -dt, u2)*v1*alpha) + vonu2;  
//
### // Uzawa function  
func real[int] JUzawa (real[int] & pp){  
    real[int] b1 = brhs1; b1 += vB1*pp;  
    real[int] b2 = brhs2; b2 += vB2*pp;  
    u1[] = AU^-1 * b1;  
    u2[] = AU^-1 * b2;  
    pp = vB1'*u1[];  
    pp += vB2'*u2[];  
    pp = -pp;  
    return pp;  
}  
//
// Preconditioner function  
func real[int] Precon (real[int] & p){  
    real[int] pa = pAA^-1*p;  
    real[int] pm = pAM^-1*p;  
    real[int] pp = alpha*pa + nu*pm;  
    return pp;  
}  
//
// Initialization  
p = 0;  
//
### // Time loop of solving UNOMI problem 
### // solving fluid
int ndt = T/dt;  
for(int i = 0; i < ndt; ++i){  
    // Update
    brhs1 = vrhs1(0, Xh);
    brhs2 = vrhs2(0, Xh);
//
// Solve
int res = LinearCG(JUzawa, p[], precon=Precon, nbiter=100, verbosity=10, veps=eps);
//    assert(res==1);
    eps = -abs(eps);
//
// Vorticity  
    w = -dy(u1) + dx(u2);  
//    plot(Th,[u1,u2],coef=0.05,cmm="time="+(dt*i),wait=0, WindowIndex = 1);
//
// Update  
    dt = min(dt, T-t);  
    t += dt;  
    if(dt < 1e-10*T) break;  
### // solving heat equation  
fespace VVh(Th,P1);  
    VVh temp,tt;  
//
solve Heat (temp,tt,solver=Crout,init=n)  
      = int2d(Th)(
      u1*dx(temp)*tt+u2*dy(temp)*tt+kappa*(dx(temp)*dx(tt)+dy(temp)*dy(tt)))  
      //-int2d(Th)(fQ*tt) 
#### // boundary conditions 
    + on(22, temp=500)//+on(3,temp=510)  
    + on(13,temp=500)  
    + on(43,temp=450)  
    ;
### // plot solutions: fluid & temperture
plot([u1,u2],temp,wait=0,fill=false,cmm="time="+(dt*i),coef=0.05);
//
//cout << "temp max = " << temp[].linfty;  
//    << ", u2 max = " << u2[].linfty
//    << ", p max = " << p[].max << endl;
}
### // loop end

![](https://github.com/user-attachments/assets/6a2ca028-86f2-4911-be15-8d0771d90992)