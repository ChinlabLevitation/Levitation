import { useState, useEffect, useRef } from "react";
import { LineChart, Line, XAxis, YAxis, Tooltip, Legend, ResponsiveContainer, ReferenceLine } from "recharts";

// ═══════════════════════════════════════════════════════════════════════════════
//  SYMMETRY UTILITIES
//  All five tensors are stored as 6-element upper-triangular arrays [t11, t12,
//  t13, t22, t23, t33], guaranteeing symmetry by construction.
//  The 6-vector encodes the independent DOF of a symmetric 3×3 matrix.
// ═══════════════════════════════════════════════════════════════════════════════

// Pack a symmetric 3×3 matrix (row-major 9-array) → 6-vector [t11,t12,t13,t22,t23,t33]
const pack   = m => [m[0], m[1], m[2], m[4], m[5], m[8]];

// Unpack 6-vector → symmetric row-major 9-array
const unpack = t => [
  t[0], t[1], t[2],
  t[1], t[3], t[4],
  t[2], t[4], t[5]
];

// Symmetrise an arbitrary 9-array: M → (M+Mᵀ)/2, then pack
const symmetrise = m => pack([
  m[0],           (m[1]+m[3])/2,  (m[2]+m[6])/2,
  (m[1]+m[3])/2,  m[4],           (m[5]+m[7])/2,
  (m[2]+m[6])/2,  (m[5]+m[7])/2, m[8]
]);

// Convenience constructors (return 6-vectors)
const diagT  = (a,b,c) => [a, 0, 0, b, 0, c];   // diagonal symmetric
const scalarT = s       => diagT(s, s, s);         // isotropic scalar × I
const zeroT  = ()       => diagT(0, 0, 0);

// Matrix-vector product using unpacked form
const mv3    = (t, v) => { const m = unpack(t);
  return [m[0]*v[0]+m[1]*v[1]+m[2]*v[2],
          m[3]*v[0]+m[4]*v[1]+m[5]*v[2],
          m[6]*v[0]+m[7]*v[1]+m[8]*v[2]]; };

// ═══════════════════════════════════════════════════════════════════════════════
//  MATH HELPERS
// ═══════════════════════════════════════════════════════════════════════════════
const norm4  = q => { const n = Math.hypot(...q); return q.map(x => x/n); };
const cross3 = (a,b) => [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
const mvT3   = (t,v) => { const m = unpack(t);   // mᵀ v = m v  (symmetric!)
  return [m[0]*v[0]+m[3]*v[1]+m[6]*v[2],
          m[1]*v[0]+m[4]*v[1]+m[7]*v[2],
          m[2]*v[0]+m[5]*v[1]+m[8]*v[2]]; };

function quatToR(q) {
  const [q0,q1,q2,q3] = q;
  return [1-2*(q2*q2+q3*q3), 2*(q1*q2-q0*q3),   2*(q1*q3+q0*q2),
          2*(q1*q2+q0*q3),   1-2*(q1*q1+q3*q3), 2*(q2*q3-q0*q1),
          2*(q1*q3-q0*q2),   2*(q2*q3+q0*q1),   1-2*(q1*q1+q2*q2)];
}
// R (3×3 row-major array) times vector
const Rv3 = (R,v) => [R[0]*v[0]+R[1]*v[1]+R[2]*v[2],
                      R[3]*v[0]+R[4]*v[1]+R[5]*v[2],
                      R[6]*v[0]+R[7]*v[1]+R[8]*v[2]];
// Rᵀ times vector
const RtV = (R,v) => [R[0]*v[0]+R[3]*v[1]+R[6]*v[2],
                      R[1]*v[0]+R[4]*v[1]+R[7]*v[2],
                      R[2]*v[0]+R[5]*v[1]+R[8]*v[2]];

const qdot = (q,w) => {
  const [q0,q1,q2,q3] = q;
  return [.5*(-q1*w[0]-q2*w[1]-q3*w[2]), .5*(q0*w[0]-q3*w[1]+q2*w[2]),
          .5*( q3*w[0]+q0*w[1]-q1*w[2]), .5*(-q2*w[0]+q1*w[1]+q0*w[2])];
};
const axQ = (ax,deg) => {
  const a = deg*Math.PI/180, n = Math.hypot(...ax), k = ax.map(x=>x/n);
  return [Math.cos(a/2), k[0]*Math.sin(a/2), k[1]*Math.sin(a/2), k[2]*Math.sin(a/2)];
};

// ═══════════════════════════════════════════════════════════════════════════════
//  PHYSICS ODE
//  All tensors passed as 6-vectors (symmetric by construction).
// ═══════════════════════════════════════════════════════════════════════════════
function makeRHS(M, Idiag, Xi, Gamma, Cr, Lambda, ar, az, Mg_g, g) {
  return s => {
    const x=s.slice(0,3), v=s.slice(3,6), q=norm4(s.slice(6,10)), wb=s.slice(10,13);
    const R = quatToR(q);

    // Temperature gradient in body frame
    const gTlab = [ar*x[0], ar*x[1], az*x[2]-Mg_g];
    const gTb   = RtV(R, gTlab);

    // Translational equation: Mv̇ = Mg - R[Γ(∇lnT)_b + Ξ v_b]
    const vb  = RtV(R, v);
    const Fb  = mv3(Gamma, gTb).map((f,i) => -(f + mv3(Xi, vb)[i]));
    const Fl  = Rv3(R, Fb);
    const vd  = [Fl[0]/M, Fl[1]/M, Fl[2]/M - g];

    // Quaternion kinematics
    const qd = qdot(q, wb);

    // Rotational equation: I ω̇_b = −ω_b×(I ω_b) − Λ(∇lnT)_b − Cr ω_b
    const LgT  = mv3(Lambda, gTb);
    const Cw   = mv3(Cr, wb);
    const tau  = [-(LgT[0]+Cw[0]), -(LgT[1]+Cw[1]), -(LgT[2]+Cw[2])];
    const Iw   = [Idiag[0]*wb[0], Idiag[1]*wb[1], Idiag[2]*wb[2]];
    const gyro = cross3(wb, Iw);
    const wbd  = [(tau[0]-gyro[0])/Idiag[0],
                  (tau[1]-gyro[1])/Idiag[1],
                  (tau[2]-gyro[2])/Idiag[2]];
    return [...v, ...vd, ...qd, ...wbd];
  };
}

function rk4(rhs, s, dt) {
  const k1=rhs(s), k2=rhs(s.map((v,i)=>v+.5*dt*k1[i]));
  const k3=rhs(s.map((v,i)=>v+.5*dt*k2[i])), k4=rhs(s.map((v,i)=>v+dt*k3[i]));
  const sn=s.map((v,i)=>v+dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6);
  const qn=norm4(sn.slice(6,10)); for(let i=0;i<4;i++) sn[6+i]=qn[i];
  return sn;
}

function simulate(p) {
  const {M,g,Idiag,Xi,Gamma,Cr,Lambda,ar,az,x0,v0,wb0,tilt,tEnd,nSteps} = p;
  // Use Gamma[3] = γ₂₂ as proxy for levitation offset (diagonal approx)
  const Mg_g = M*g/(Gamma[3]||1);
  const rhs  = makeRHS(M,Idiag,Xi,Gamma,Cr,Lambda,ar,az,Mg_g,g);
  const q0   = axQ([1,0,0], tilt);
  let s = [...x0,...v0,...q0,...wb0];
  const dt=tEnd/nSteps, out=[];
  for (let i=0;i<=nSteps;i++) {
    const R = quatToR(norm4(s.slice(6,10)));
    out.push({
      t:     +(i*dt).toFixed(3),
      x:     +s[0].toFixed(4), y:+s[1].toFixed(4), z:+s[2].toFixed(4),
      w1:    +s[10].toFixed(3), w2:+s[11].toFixed(3), w3:+s[12].toFixed(3),
      speed: +Math.hypot(s[3],s[4],s[5]).toFixed(4),
      tilt:  +(Math.acos(Math.max(-1,Math.min(1,R[8])))*180/Math.PI).toFixed(1),
      R:     [...R],
    });
    if (i<nSteps) s=rk4(rhs,s,dt);
  }
  return out;
}

// ═══════════════════════════════════════════════════════════════════════════════
//  3-D RENDERING
// ═══════════════════════════════════════════════════════════════════════════════
function proj3(pt, cx, cy, sc, elev, azim) {
  const ea=elev*Math.PI/180, aa=azim*Math.PI/180, [x0,y0,z0]=pt;
  const x1=x0*Math.cos(aa)+z0*Math.sin(aa);
  const y1=-x0*Math.sin(aa)*Math.sin(ea)+y0*Math.cos(ea)+z0*Math.cos(aa)*Math.sin(ea);
  const z1=-x0*Math.sin(aa)*Math.cos(ea)-y0*Math.sin(ea)+z0*Math.cos(aa)*Math.cos(ea);
  return [cx+x1*sc, cy-y1*sc, z1];
}
function buildCyl(hl=0.10, r=0.034, n=14) {
  const V=[], F=[];
  for(let i=0;i<n;i++){const t=2*Math.PI*i/n; V.push([r*Math.cos(t),r*Math.sin(t), hl]);}
  for(let i=0;i<n;i++){const t=2*Math.PI*i/n; V.push([r*Math.cos(t),r*Math.sin(t),-hl]);}
  V.push([0,0,hl]); V.push([0,0,-hl]);
  for(let i=0;i<n;i++){
    const j=(i+1)%n;
    F.push({v:[i,j,j+n,i+n],k:'s'});
    F.push({v:[2*n,i,j],k:'t'});
    F.push({v:[2*n+1,j+n,i+n],k:'b'});
  }
  return {V,F};
}
const CYL=buildCyl();

// ═══════════════════════════════════════════════════════════════════════════════
//  DEFAULTS  — all tensors as symmetric 6-vectors
// ═══════════════════════════════════════════════════════════════════════════════
const DEF = () => ({
  M:1, g:9.81,
  Idiag: [0.02, 0.02, 0.008],
  Xi:     scalarT(2.0),   // Ξ = 2·I  (isotropic)
  Gamma:  scalarT(5.0),   // Γ = 5·I
  Cr:     scalarT(0.5),   // Cᵣ = 0.5·I
  Lambda: zeroT(),        // Λ = 0  (symmetric body)
  ar:0.5, az:0.5,
  x0:[0.5,0.3,0.2], v0:[0,0,0], wb0:[0,0.5,2],
  tilt:30, tEnd:20, nSteps:1200,
});

// ═══════════════════════════════════════════════════════════════════════════════
//  SYMMETRIC TENSOR EDITOR
//  Displays 3×3 matrix but keeps only the 6 upper-triangle entries editable.
//  Lower triangle is shown as read-only mirrors to make symmetry visible.
// ═══════════════════════════════════════════════════════════════════════════════
function SymTensorEditor({name, desc, symbol, value, onChange, color}) {
  const [open,setOpen] = useState(false);
  // value is 6-vector [t11,t12,t13,t22,t23,t33]
  // Display indices: row i, col j → 6-vec index
  const idx = (i,j) => i<=j ? [0,1,2,4,5,6][i*3+j-i*(i+1)/2] :
                               [0,1,2,4,5,6][j*3+i-j*(j+1)/2]; // mirror
  // Mapping from (i,j) to 6-vector position
  const sixIdx = (i,j) => {
    if(i===0&&j===0) return 0; if(i===0&&j===1) return 1; if(i===0&&j===2) return 2;
    if(i===1&&j===0) return 1; if(i===1&&j===1) return 3; if(i===1&&j===2) return 4;
    if(i===2&&j===0) return 2; if(i===2&&j===1) return 4; if(i===2&&j===2) return 5;
    return 0;
  };
  const isUpper = (i,j) => i<=j;
  const isDiag  = [1,2,4].every(i=>value[i]===0);
  const diagVals = [value[0],value[3],value[5]];

  const ni = (upper) => ({
    background: upper ? '#060e1a' : '#0a1525',
    border: `1px solid ${upper ? '#1e3048' : '#0f1e30'}`,
    color: upper ? '#e2e8f0' : '#4a6080',
    borderRadius:2, padding:'1px 3px', fontSize:10, fontFamily:'monospace',
    width:46, textAlign:'center',
    cursor: upper ? 'text' : 'not-allowed',
  });

  return (
    <div style={{marginBottom:4}}>
      <div onClick={()=>setOpen(o=>!o)} style={{display:'flex',justifyContent:'space-between',
        alignItems:'center',cursor:'pointer',background:'#0d1828',borderRadius:4,padding:'3px 6px',
        border:`1px solid ${open?color:'#1e3048'}`}}>
        <span style={{color:color,fontSize:10,fontWeight:700}}>{symbol} — {name}</span>
        <span style={{color:'#475569',fontSize:9}}>
          {isDiag?`diag(${diagVals.map(v=>v.toFixed(2)).join(',')})`:' symm'}{open?'▲':'▼'}
        </span>
      </div>
      {open&&(
        <div style={{background:'#0a1220',border:'1px solid #1e3048',borderRadius:4,padding:6,marginTop:2}}>
          <div style={{color:'#64748b',fontSize:9,marginBottom:4,lineHeight:1.4}}>{desc}</div>
          <div style={{color:'#334155',fontSize:8,marginBottom:3}}>
            Upper triangle = independent DOF &nbsp;|&nbsp; Lower triangle = mirror (read-only)
          </div>
          <div style={{display:'flex',flexDirection:'column',gap:2,alignItems:'center'}}>
            {[0,1,2].map(ri=>(
              <div key={ri} style={{display:'flex',gap:2}}>
                {[0,1,2].map(ci=>{
                  const upper = isUpper(ri,ci);
                  const si    = sixIdx(ri,ci);
                  return (
                    <input key={ci} type="number" step={0.01} value={value[si]}
                      style={{...ni(upper), border:`1px solid ${ri===ci?color:upper?'#1e3048':'#0f1e30'}`}}
                      readOnly={!upper}
                      onChange={upper ? e=>{
                        const v=[...value]; v[si]=parseFloat(e.target.value)||0; onChange(v);
                      } : undefined}/>
                  );
                })}
              </div>
            ))}
          </div>
          <div style={{display:'flex',gap:4,marginTop:4,flexWrap:'wrap'}}>
            {[
              ['Isotropic', v=>scalarT(v[0])],
              ['Diag only', v=>diagT(v[0],v[3],v[5])],
              ['Zero',      ()=>zeroT()],
              ['Symmetrise input', v=>{
                // Allow user to type arbitrary values then auto-symmetrise
                const m=unpack(v); return symmetrise(m);
              }],
            ].map(([lbl,fn])=>(
              <button key={lbl} onClick={()=>onChange(fn(value))}
                style={{background:'#1e3048',color:'#94a3b8',border:'none',borderRadius:3,
                        padding:'2px 6px',cursor:'pointer',fontSize:9}}>{lbl}</button>
            ))}
          </div>
        </div>
      )}
    </div>
  );
}

// ═══════════════════════════════════════════════════════════════════════════════
//  SECTION WRAPPER
// ═══════════════════════════════════════════════════════════════════════════════
function Section({title,children,defaultOpen=false}){
  const [open,setOpen]=useState(defaultOpen);
  return(
    <div style={{marginBottom:4}}>
      <div onClick={()=>setOpen(o=>!o)} style={{display:'flex',justifyContent:'space-between',
        alignItems:'center',cursor:'pointer',background:'#0f1e30',borderRadius:4,
        padding:'4px 8px',border:'1px solid #1e3048',marginBottom:open?4:0}}>
        <span style={{color:'#7dd3fc',fontSize:10,fontWeight:700}}>{title}</span>
        <span style={{color:'#475569',fontSize:10}}>{open?'▲':'▼'}</span>
      </div>
      {open&&<div>{children}</div>}
    </div>
  );
}

// ═══════════════════════════════════════════════════════════════════════════════
//  MAIN APP
// ═══════════════════════════════════════════════════════════════════════════════
export default function App(){
  const [cfg,setCfg]       = useState(DEF);
  const [frames,setFrames] = useState(null);
  const [fi,setFi]         = useState(0);
  const [play,setPlay]     = useState(false);
  const [elev,setElev]     = useState(22);
  const [azim,setAzim]     = useState(38);
  const [busy,setBusy]     = useState(false);
  const [msg,setMsg]       = useState('Configure parameters, then click ▶ Run');
  const playRef=useRef(false), fiRef=useRef(0), rafRef=useRef(null);

  const sc=(k,v)=>setCfg(c=>({...c,[k]:v}));
  const sv=(k,i,v)=>setCfg(c=>{const a=[...c[k]];a[i]=v;return{...c,[k]:a};});

  function run(){
    setBusy(true);setPlay(false);playRef.current=false;setMsg('Simulating…');
    setTimeout(()=>{
      try{
        const t0=Date.now(), f=simulate(cfg);
        setFrames(f);setFi(0);fiRef.current=0;
        setMsg(`✓ Done in ${Date.now()-t0}ms · ${f.length} frames · vmax=${Math.max(...f.map(x=>x.speed)).toFixed(3)} m/s`);
      }catch(e){setMsg('Error: '+e.message);}
      setBusy(false);
    },10);
  }
  function resetDef(){setCfg(DEF());setFrames(null);setMsg('Reset to defaults.');}

  useEffect(()=>{
    if(!play||!frames)return;
    playRef.current=true;
    let last=null;
    const step=ts=>{
      if(!playRef.current)return;
      if(last===null)last=ts;
      if(ts-last>38){last=ts;fiRef.current=(fiRef.current+1)%frames.length;setFi(fiRef.current);}
      rafRef.current=requestAnimationFrame(step);
    };
    rafRef.current=requestAnimationFrame(step);
    return()=>{playRef.current=false;cancelAnimationFrame(rafRef.current);};
  },[play,frames]);

  const togglePlay=()=>{if(!frames)return;setPlay(v=>!v);playRef.current=!playRef.current;};
  const goReset=()=>{setPlay(false);playRef.current=false;setFi(0);fiRef.current=0;};
  const goEnd=()=>{setPlay(false);playRef.current=false;const f=frames.length-1;setFi(f);fiRef.current=f;};

  const W=360,H=300,CX=W/2,CY=H/2,SC=200;
  const fr=frames&&frames[fi];
  const pj=pt=>proj3(pt,CX,CY,SC,elev,azim);

  function scene(){
    if(!fr) return(
      <text x={W/2} y={H/2} textAnchor="middle" fill="#2a4a6a" fontSize="12" fontFamily="monospace">
        Click ▶ Run Simulation
      </text>
    );
    const R=fr.R, xc=[fr.x,fr.y,fr.z];
    const labAxes=[['X',[.35,0,0],'#4fc3f7'],['Y',[0,.35,0],'#81c784'],['Z',[0,0,.35],'#ffb74d']].map(([l,e,c])=>{
      const [ox,oy]=pj([0,0,0]), [ex,ey]=pj(e);
      return(<g key={l}>
        <line x1={ox} y1={oy} x2={ex} y2={ey} stroke={c} strokeWidth="1" strokeDasharray="4,3" opacity="0.5"/>
        <text x={ex+4} y={ey+4} fill={c} fontSize="11" fontFamily="monospace">{l}</text>
      </g>);
    });
    const pathFull=frames.filter((_,i)=>i%6===0).map((f,i)=>{
      const [px,py]=pj([f.x,f.y,f.z]); return`${i===0?'M':'L'}${px.toFixed(1)},${py.toFixed(1)}`;
    }).join(' ');
    const pathSoFar=frames.slice(0,fi+1).filter((_,i)=>i%3===0).map((f,i)=>{
      const [px,py]=pj([f.x,f.y,f.z]); return`${i===0?'M':'L'}${px.toFixed(1)},${py.toFixed(1)}`;
    }).join(' ');
    const tv=CYL.V.map(v=>{
      const wx=R[0]*v[0]+R[1]*v[1]+R[2]*v[2]+xc[0];
      const wy=R[3]*v[0]+R[4]*v[1]+R[5]*v[2]+xc[1];
      const wz=R[6]*v[0]+R[7]*v[1]+R[8]*v[2]+xc[2];
      return pj([wx,wy,wz]);
    });
    const cylFaces=[...CYL.F]
      .map(f=>({...f,depth:f.v.reduce((s,i)=>s+tv[i][2],0)/f.v.length}))
      .sort((a,b)=>a.depth-b.depth)
      .map((f,i)=>{
        const pts=f.v.map(vi=>`${tv[vi][0].toFixed(1)},${tv[vi][1].toFixed(1)}`).join(' ');
        const col=f.k==='t'?'#7986cb':f.k==='b'?'#5c6bc0':'#546e7a';
        return <polygon key={i} points={pts} fill={col} stroke="#78909c" strokeWidth="0.3" opacity="0.9"/>;
      });
    const bodyAxes=[[0,'#ef5350'],[1,'#66bb6a'],[2,'#42a5f5']].map(([ci,c])=>{
      const ei=[R[ci],R[3+ci],R[6+ci]];
      const tip=[xc[0]+ei[0]*0.22,xc[1]+ei[1]*0.22,xc[2]+ei[2]*0.22];
      const [bx,by]=pj(xc), [tx,ty]=pj(tip);
      return(<g key={ci}>
        <line x1={bx.toFixed(1)} y1={by.toFixed(1)} x2={tx.toFixed(1)} y2={ty.toFixed(1)}
              stroke={c} strokeWidth="2.5" strokeLinecap="round"/>
        <circle cx={tx.toFixed(1)} cy={ty.toFixed(1)} r="4" fill={c}/>
      </g>);
    });
    const [dx,dy]=pj(xc);
    return(<>
      {labAxes}
      <path d={pathFull}  fill="none" stroke="#4fc3f7" strokeWidth="0.5" opacity="0.15"/>
      <path d={pathSoFar} fill="none" stroke="#4fc3f7" strokeWidth="1.5" opacity="0.7"/>
      {cylFaces}{bodyAxes}
      <circle cx={dx.toFixed(1)} cy={dy.toFixed(1)} r="4" fill="white" opacity="0.95"/>
    </>);
  }

  const stride=Math.max(1,Math.floor((frames?.length||1)/300));
  const cd=frames?frames.filter((_,i)=>i%stride===0):[];

  const BG='#0e1117',P='#141e2e',BD='#1e3048';
  const lb={color:'#7c8fa6',fontSize:10,fontFamily:'monospace'};
  const ni={background:'#0a1220',border:'1px solid #1e3048',color:'#e2e8f0',
            borderRadius:3,padding:'2px 4px',fontSize:11,fontFamily:'monospace',width:56};
  const bt=(bg='#1d4ed8')=>({background:bg,color:'white',border:'none',borderRadius:4,
            padding:'4px 10px',cursor:'pointer',fontSize:11,fontWeight:700});
  const card={background:P,border:`1px solid ${BD}`,borderRadius:7,padding:8};

  return(
    <div style={{background:BG,minHeight:'100vh',color:'#e2e8f0',fontFamily:'monospace',padding:8}}>
      <div style={{display:'flex',justifyContent:'space-between',alignItems:'center',marginBottom:6}}>
        <div>
          <div style={{color:'#7dd3fc',fontSize:14,fontWeight:700}}>Rigid-Body Thermophoresis Simulator</div>
          <div style={{color:'#334155',fontSize:9}}>
            RK4 · 13-DOF quaternion ODE · Harmonic trap · All tensors symmetric (Onsager + dissipation) · 24 DOF total
          </div>
        </div>
        <button style={bt('#374151')} onClick={resetDef}>↺ Reset</button>
      </div>

      {/* DOF summary bar */}
      <div style={{background:'#0a1525',border:'1px solid #1e3048',borderRadius:5,padding:'4px 10px',
                   marginBottom:6,fontSize:9,color:'#64748b',display:'flex',gap:16,flexWrap:'wrap'}}>
        <span style={{color:'#7dd3fc',fontWeight:700}}>Symmetry constraints:</span>
        <span><span style={{color:'#4fc3f7'}}>Ξ</span> symmetric (dissipation) · 6 DOF</span>
        <span><span style={{color:'#81c784'}}>Γ</span> symmetric (Onsager) · 6 DOF</span>
        <span><span style={{color:'#f48fb1'}}>Λ</span> symmetric (Onsager) · 6 DOF</span>
        <span><span style={{color:'#ffb74d'}}>Cᵣ</span> symmetric (dissipation) · 6 DOF</span>
        <span style={{color:'#4ade80',fontWeight:700}}>Total: 24 DOF</span>
      </div>

      <div style={{display:'flex',gap:8,flexWrap:'wrap',alignItems:'flex-start'}}>

        {/* LEFT: params */}
        <div style={{...card,width:248,flexShrink:0,maxHeight:'92vh',overflowY:'auto'}}>

          <Section title="🔵 Particle — Mass & Inertia" defaultOpen={true}>
            <div style={{padding:'0 2px 4px'}}>
              {[['Mass M (kg)','M',0.1],['Gravity g (m/s²)','g',0.01]].map(([l,k,st])=>(
                <div key={k} style={{display:'flex',justifyContent:'space-between',alignItems:'center',marginBottom:2}}>
                  <span style={lb}>{l}</span>
                  <input type="number" step={st} value={cfg[k]} style={ni} onChange={e=>sc(k,parseFloat(e.target.value)||0)}/>
                </div>
              ))}
              <div style={{...lb,marginBottom:1}}>I = diag(I₁,I₂,I₃) kg·m²</div>
              <div style={{display:'flex',gap:2}}>
                {cfg.Idiag.map((v,i)=>(
                  <input key={i} type="number" step={0.001} value={v} style={{...ni,width:46}}
                    onChange={e=>sv('Idiag',i,parseFloat(e.target.value)||0)}/>
                ))}
              </div>
            </div>
          </Section>

          <Section title="🌡️ Temperature Field — Harmonic Trap" defaultOpen={true}>
            <div style={{padding:'0 2px 4px'}}>
              <div style={{color:'#64748b',fontSize:9,marginBottom:4,lineHeight:1.5}}>
                ∇lnT = (αr·x, αr·y, αz·z − Mg/γ)
              </div>
              {[['αr — radial','ar',0.05],['αz — axial','az',0.05]].map(([l,k,st])=>(
                <div key={k} style={{display:'flex',justifyContent:'space-between',alignItems:'center',marginBottom:2}}>
                  <span style={lb}>{l}</span>
                  <input type="number" step={st} value={cfg[k]} style={ni} onChange={e=>sc(k,parseFloat(e.target.value)||0)}/>
                </div>
              ))}
            </div>
          </Section>

          <Section title="⚙️ Response Tensors (symmetric, Onsager)" defaultOpen={true}>
            <div style={{padding:'0 2px 4px'}}>
              <div style={{color:'#475569',fontSize:9,marginBottom:4,lineHeight:1.4}}>
                All tensors stored as 6-vectors [t₁₁,t₁₂,t₁₃,t₂₂,t₂₃,t₃₃].<br/>
                Lower triangle is read-only (mirrors upper triangle).
              </div>
              <SymTensorEditor name="Resistance" symbol="Ξ"
                desc="F_drag = −Ξ v_b  |  symmetric by energy dissipation  |  6 DOF"
                color="#4fc3f7" value={cfg.Xi} onChange={v=>sc('Xi',v)}/>
              <SymTensorEditor name="Thermo Force Γ" symbol="Γ"
                desc="F_T = −Γ (∇lnT)_b  |  symmetric by Onsager reciprocity  |  6 DOF"
                color="#81c784" value={cfg.Gamma} onChange={v=>sc('Gamma',v)}/>
              <SymTensorEditor name="Thermo Torque Λ" symbol="Λ"
                desc="τ_T = −Λ (∇lnT)_b  |  symmetric by Onsager reciprocity  |  6 DOF  |  zero for mirror-symmetric bodies"
                color="#f48fb1" value={cfg.Lambda} onChange={v=>sc('Lambda',v)}/>
              <SymTensorEditor name="Rotational Drag" symbol="Cᵣ"
                desc="τ_drag = −Cᵣ ω_b  |  symmetric by energy dissipation  |  6 DOF"
                color="#ffb74d" value={cfg.Cr} onChange={v=>sc('Cr',v)}/>
            </div>
          </Section>

          <Section title="🚀 Initial Conditions" defaultOpen={true}>
            <div style={{padding:'0 2px 4px'}}>
              <div style={{display:'flex',justifyContent:'space-between',alignItems:'center',marginBottom:4}}>
                <span style={lb}>Tilt θ₀ (°)</span>
                <input type="number" step={1} value={cfg.tilt} style={ni} onChange={e=>sc('tilt',parseFloat(e.target.value)||0)}/>
              </div>
              {[['x₀ (m)','x0',0.05],['v₀ (m/s)','v0',0.1],['ω₀ (rad/s)','wb0',0.1]].map(([l,k,st])=>(
                <div key={k} style={{marginBottom:4}}>
                  <div style={{...lb,marginBottom:1}}>{l}</div>
                  <div style={{display:'flex',gap:2}}>
                    {cfg[k].map((v,i)=>(
                      <input key={i} type="number" step={st} value={v} style={{...ni,width:46}}
                        onChange={e=>sv(k,i,parseFloat(e.target.value)||0)}/>
                    ))}
                  </div>
                </div>
              ))}
            </div>
          </Section>

          <Section title="⏱ Simulation">
            <div style={{padding:'0 2px 4px'}}>
              {[['t_end (s)','tEnd',1],['Steps','nSteps',100]].map(([l,k,st])=>(
                <div key={k} style={{display:'flex',justifyContent:'space-between',alignItems:'center',marginBottom:2}}>
                  <span style={lb}>{l}</span>
                  <input type="number" step={st} value={cfg[k]} style={ni} onChange={e=>sc(k,parseFloat(e.target.value)||0)}/>
                </div>
              ))}
            </div>
          </Section>

          <button style={{...bt('#0f766e'),width:'100%',marginTop:4,padding:'7px'}} onClick={run} disabled={busy}>
            {busy?'⏳ Simulating…':'▶ Run Simulation'}
          </button>
          <div style={{marginTop:4,fontSize:9,lineHeight:1.5,
            color:msg.startsWith('✓')?'#4ade80':msg.startsWith('Error')?'#f87171':'#64748b'}}>
            {msg}
          </div>
          {frames&&(
            <div style={{marginTop:4,fontSize:9,color:'#475569',lineHeight:1.6}}>
              x:[{Math.min(...frames.map(f=>f.x)).toFixed(2)},{Math.max(...frames.map(f=>f.x)).toFixed(2)}]<br/>
              z:[{Math.min(...frames.map(f=>f.z)).toFixed(2)},{Math.max(...frames.map(f=>f.z)).toFixed(2)}]<br/>
              vmax={Math.max(...frames.map(f=>f.speed)).toFixed(3)} m/s
            </div>
          )}
        </div>

        {/* MIDDLE: 3D view */}
        <div style={{flex:'1 1 350px',minWidth:300}}>
          <div style={{...card,marginBottom:6}}>
            <div style={{display:'flex',justifyContent:'space-between',alignItems:'center',marginBottom:4}}>
              <span style={{color:'#7dd3fc',fontSize:10,fontWeight:700}}>
                3-D VIEW &nbsp;
                <span style={{color:'#ef5350'}}>■</span> e₁ &nbsp;
                <span style={{color:'#66bb6a'}}>■</span> e₂ &nbsp;
                <span style={{color:'#42a5f5'}}>■</span> e₃
              </span>
              {fr&&<span style={{color:'#475569',fontSize:9}}>t={fr.t.toFixed(2)}s  tilt={fr.tilt.toFixed(1)}°</span>}
            </div>
            <svg width={W} height={H}
              style={{display:'block',margin:'0 auto',background:'#060d1a',borderRadius:5,border:'1px solid #1e3048'}}>
              {scene()}
            </svg>
            <div style={{display:'flex',gap:16,marginTop:6,flexWrap:'wrap'}}>
              <label style={{...lb,display:'flex',alignItems:'center',gap:4}}>Elev {elev}°
                <input type="range" min={-80} max={80} value={elev} onChange={e=>setElev(+e.target.value)} style={{width:80}}/>
              </label>
              <label style={{...lb,display:'flex',alignItems:'center',gap:4}}>Azim {azim}°
                <input type="range" min={0} max={360} value={azim} onChange={e=>setAzim(+e.target.value)} style={{width:80}}/>
              </label>
            </div>
          </div>
          {frames&&(
            <div style={card}>
              <div style={{display:'flex',gap:5,alignItems:'center',marginBottom:4}}>
                <button style={bt()} onClick={goReset}>⏮</button>
                <button style={bt(play?'#92400e':'#1d4ed8')} onClick={togglePlay}>{play?'⏸ Pause':'▶ Play'}</button>
                <button style={bt()} onClick={goEnd}>⏭</button>
                <span style={{color:'#475569',fontSize:9,marginLeft:4}}>
                  {fr&&`t=${fr.t.toFixed(2)}s  frame ${fi}/${frames.length-1}`}
                </span>
              </div>
              <input type="range" min={0} max={frames.length-1} value={fi} style={{width:'100%'}}
                onChange={e=>{const f=+e.target.value;setFi(f);fiRef.current=f;}}/>
            </div>
          )}
        </div>

        {/* RIGHT: charts */}
        <div style={{flex:'1 1 250px',minWidth:240,display:'flex',flexDirection:'column',gap:6}}>
          {[
            {title:'CoM Position (m)',        lines:[{k:'x',c:'#4fc3f7',n:'x'},{k:'y',c:'#81c784',n:'y'},{k:'z',c:'#ffb74d',n:'z'}]},
            {title:'Angular Velocity (rad/s)', lines:[{k:'w1',c:'#ef9a9a',n:'ω₁'},{k:'w2',c:'#ce93d8',n:'ω₂'},{k:'w3',c:'#80cbc4',n:'ω₃'}]},
            {title:'Tilt (°) & Speed (m/s)',   lines:[{k:'tilt',c:'#fff176',n:'tilt°'},{k:'speed',c:'#f48fb1',n:'speed'}]},
          ].map(({title,lines},ti)=>(
            <div key={ti} style={card}>
              <div style={{color:'#7dd3fc',fontSize:9,fontWeight:700,marginBottom:3}}>{title}</div>
              <ResponsiveContainer width="100%" height={120}>
                <LineChart data={cd} margin={{top:2,right:4,bottom:0,left:-20}}>
                  <XAxis dataKey="t" stroke="#0f1f33" tick={{fill:'#334155',fontSize:8}}/>
                  <YAxis stroke="#0f1f33" tick={{fill:'#334155',fontSize:8}}/>
                  <Tooltip contentStyle={{background:'#0a1628',border:'none',fontSize:9}} itemStyle={{color:'#94a3b8'}}/>
                  <Legend wrapperStyle={{fontSize:9,color:'#64748b'}}/>
                  {lines.map(({k,c,n})=><Line key={k} type="monotone" dataKey={k} name={n} stroke={c} dot={false} strokeWidth={1.2}/>)}
                  {fr&&<ReferenceLine x={fr.t} stroke="#ff6e40" strokeWidth={1} strokeDasharray="3,2"/>}
                </LineChart>
              </ResponsiveContainer>
            </div>
          ))}
        </div>

      </div>
    </div>
  );
}
