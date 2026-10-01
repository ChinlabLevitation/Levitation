(* ::Package:: *)
(* ShapeLevitation.wl -- shared machinery for the shape notebooks (Shape1_Ellipsoid, Shape2_EllipsoidOffset, Shape3_Harmonics).

   Physics: free-molecule gas-surface stress of a convex body (notes v22, Eq. 50), integrated over a star-shaped surface
   r = ell R(n) n, gives the complete linear response about the centre of mass (CoM), in the body frame:
        (F, tau) = - Rsp . ( grad ln T , v , omega )          (6 x 9)
   Blocks: f_th = Rsp[[1;;3,1;;3]], f_drag = Rsp[[1;;3,4;;6]], B = Rsp[[1;;3,7;;9]],
           T_th = Rsp[[4;;6,1;;3]], B^T = Rsp[[4;;6,4;;6]], T_drag = Rsp[[4;;6,7;;9]].

   Two independent engines:
   * symbolicShape: exact perturbation series in eps for R = polynomial-in-n perturbations of the unit sphere.
     Everything reduces to moments  Int x^a y^b z^c dOmega  over the unit sphere, so it is fast and exact.
     Natural units: lengths ell, f_th in gamma0, f_drag in xi0, B and B^T in xi0 ell, T_th in gamma0 ell, T_drag in c0,
     mass in rho ell^3, inertia in rho ell^5.
   * makeParticle: SI numbers at finite eps by Gauss-Legendre (theta) x trapezoid (phi) quadrature -- used for simulations. *)

BeginPackage["ShapeLevitation`"];

sphereMoment::usage = "sphereMoment[poly, {x,y,z}] integrates a polynomial in the unit normal over the unit sphere.";
symbolicShape::usage = "symbolicShape[R, {x,y,z}, eps, order, density] -> mass, CoM, inertia and the 6x9 response about the CoM as series in eps (natural units).";
blocks::usage = "blocks[Rsp] splits a 6x9 response into f_th, f_drag, B, T_th, Bt, T_drag.";
gasModel::usage = "gasModel[<|\"PTorr\" -> 5, \"T\" -> 300, \"kappa\" -> 0.026, \"a\" -> 1|>] returns the gas parameters and the sphere coefficients.";
makeParticle::usage = "makeParticle[Rfun, ell, rhoP, gas, density] -> SI particle (mass, CoM, inertia, response about the CoM). Rfun[{x,y,z}] is the radius in units of ell.";
steadyMotions::usage = "steadyMotions[pt] finds every steady spin/orbit state at levitation and its linear stability.";
simulate::usage = "simulate[pt, G, q0, w0, tEnd] integrates the nonlinear rigid-body equations in a uniform gradient grad ln T = -G Z.";
orbitFromRun::usage = "orbitFromRun[run, t1, t2] measures radius, period, sense and drift of the horizontal CoM motion.";
convexityMargin::usage = "convexityMargin[Rfun] is the smallest principal curvature (units 1/ell) of the surface on a grid; positive means convex.";
shapeDiagram::usage = "shapeDiagram[Rfun, com, opts] draws the body (colour = R - 1) with body axes, geometric centre (black) and CoM (magenta).";
hoverDrift::usage = "hoverDrift[pt, u] gives the levitation gradient and terminal (glide) velocity for a non-rotating body with lab Z = u in the body frame.";
bodyMesh::usage = "bodyMesh[Rfun, n] precomputes a coloured surface mesh (units of ell) for animations.";
drawBody::usage = "drawBody[mesh, x, R, scale] places the mesh at x with orientation R, lengths multiplied by scale.";
hat::usage = "hat[a] is the skew matrix [a]x.";
quatR::usage = "quatR[q] is the rotation matrix of the unit quaternion q (body -> lab).";
qFromTo::usage = "qFromTo[a, b] is a unit quaternion rotating unit vector a onto b.";
eulerZXZ::usage = "eulerZXZ[a, b, c] is the rotation matrix Rz[a].Rx[b].Rz[c].";
$gEarth::usage = "g in m/s^2.";

Begin["`Private`"];

$gEarth = 9.81; $kB = 1.380649*^-23; $mAir = 4.81*^-26; $torr = 133.322;
hat[{a1_, a2_, a3_}] := {{0, -a3, a2}, {a3, 0, -a1}, {-a2, a1, 0}};
quatR[{q0_, q1_, q2_, q3_}] := (q0^2 - {q1, q2, q3}.{q1, q2, q3}) IdentityMatrix[3] + 2 Outer[Times, {q1, q2, q3}, {q1, q2, q3}] + 2 q0 hat[{q1, q2, q3}];
Qmat[{q0_, q1_, q2_, q3_}] := {{-q1, -q2, -q3}, {q0, -q3, q2}, {q3, q0, -q1}, {-q2, q1, q0}};
qFromTo[a_, b_] := Module[{ax = Cross[a, b], c = a.b},
   If[Norm[ax] < 10^-12, If[c > 0, {1., 0., 0., 0.}, Join[{0.}, Normalize[Cross[a, If[Abs[a[[1]]] < 0.9, {1, 0, 0}, {0, 1, 0}]]]]],
    N@Normalize[Join[{1 + c}, ax]]]];
eulerZXZ[a_, b_, c_] := RotationMatrix[a, {0, 0, 1}].RotationMatrix[b, {1, 0, 0}].RotationMatrix[c, {0, 0, 1}];
blocks[Rsp_] := <|"f_th" -> Rsp[[1 ;; 3, 1 ;; 3]], "f_drag" -> Rsp[[1 ;; 3, 4 ;; 6]], "B" -> Rsp[[1 ;; 3, 7 ;; 9]],
   "T_th" -> Rsp[[4 ;; 6, 1 ;; 3]], "Bt" -> Rsp[[4 ;; 6, 4 ;; 6]], "T_drag" -> Rsp[[4 ;; 6, 7 ;; 9]]|>;

(* ---------------- symbolic engine ---------------- *)
monoInt[a_, b_, c_] := monoInt[a, b, c] = If[OddQ[a] || OddQ[b] || OddQ[c], 0,
    2 Gamma[(a + 1)/2] Gamma[(b + 1)/2] Gamma[(c + 1)/2]/Gamma[(a + b + c + 3)/2]];
sphereMoment[expr_, vars_List] := Module[{e = Expand[expr]},
   If[e === 0, 0, Total[(#[[2]] monoInt @@ #[[1]]) & /@ CoefficientRules[e, vars]]]];

(* truncate a polynomial in eps *)
trunc[e_, eps_, n_] := Expand[e] /. eps^k_ /; k > n :> 0;
truncM[m_, eps_, n_] := Map[trunc[#, eps, n] &, m, {ArrayDepth[m]}];

(* homogeneous parts of a density polynomial rho(X,Y,Z): list of {degree, f_d(n)} *)
homParts[dens_, vars_] := Module[{t, e},
   e = Expand[dens /. Thread[vars -> t vars]];
   DeleteCases[Table[{d, Coefficient[e, t, d]}, {d, 0, Exponent[e, t]}], {_, 0}]];

symbolicShape[R0_, vars : {x_, y_, z_}, eps_, order_, density_ : 1, OptionsPattern[{"Accommodation" -> 1}]] := Module[
   {n = vars, R, grad, gs, invR2, xS, S, absS, invS, P, Id = IdentityMatrix[3], aa = OptionValue["Accommodation"], be, kt, kv,
    parts, mass, mom1, mom2, com, Icom, r, rx, Kq, KV, blk, Rsp},
   R = trunc[Normal@Series[R0, {eps, 0, order}], eps, order];
   (* tangential gradient on the unit sphere *)
   grad = D[R, {n}]; gs = trunc[grad - n (n.grad), eps, order];
   (* S = (surface normal area element per dOmega) = R^2 n - R grad_s R ; |S| = R^2 Sqrt[1 + |gs|^2/R^2] *)
   invR2 = trunc[1 + Sum[(k + 1) (1 - R)^k, {k, 1, order}], eps, order];
   xS = trunc[(gs.gs) invR2, eps, order];
   S = truncM[R^2 n - R gs, eps, order];
   absS = trunc[R^2 (1 + Sum[Binomial[1/2, k] xS^k, {k, 1, Floor[order/2]}]), eps, order];
   invS = trunc[invR2 (1 + Sum[Binomial[-1/2, k] xS^k, {k, 1, Floor[order/2]}]), eps, order];
   P = truncM[Outer[Times, S, S] invS, eps, order];
   (* mass properties with density rho (1 + ...), radial integrals done analytically *)
   parts = homParts[density, vars];
   mass = Sum[sphereMoment[trunc[pp[[2]] R^(pp[[1]] + 3)/(pp[[1]] + 3), eps, order], n], {pp, parts}];
   mom1 = Table[Sum[sphereMoment[trunc[pp[[2]] n[[i]] R^(pp[[1]] + 4)/(pp[[1]] + 4), eps, order], n], {pp, parts}], {i, 3}];
   mom2 = Table[Sum[sphereMoment[trunc[pp[[2]] (Boole[i == j] - n[[i]] n[[j]]) R^(pp[[1]] + 5)/(pp[[1]] + 5), eps, order], n], {pp, parts}], {i, 3}, {j, 3}];
   com = trunc[Normal@Series[mom1/mass, {eps, 0, order}], eps, order];
   Icom = truncM[mom2 - mass ((com.com) Id - Outer[Times, com, com]), eps, order];
   (* normalised kernels: Int Kq dOmega = Id (gamma0), -Int KV dOmega = Id (xi0) for the unit sphere *)
   be = 2 - aa + Pi aa/4;
   kt = 3/(16 Pi); kv = 3/(4 Pi (aa + be));
   Kq = kt (aa (absS Id - P) + 2 (2 - aa) P);
   KV = -kv ((aa/2) (absS Id - P) + be P);
   r = truncM[R n - com, eps, order]; rx = hat[r];
   blk = {Kq, -KV, truncM[KV.rx, eps, order], truncM[rx.Kq, eps, order], truncM[-rx.KV, eps, order], truncM[rx.KV.rx, eps, order] (aa + be)/aa};
   Rsp = Map[sphereMoment[#, n] &, blk, {3}];
   Rsp = ArrayFlatten[{{Rsp[[1]], Rsp[[2]], Rsp[[3]]}, {Rsp[[4]], Rsp[[5]], Rsp[[6]]}}];
   <|"R" -> R, "M" -> Expand[mass], "com" -> Expand[com], "I" -> Expand[Icom], "Rsp" -> Expand[Rsp], "blocks" -> Map[Expand, blocks[Rsp], {3}]|>];

(* ---------------- numerical engine (SI) ---------------- *)
gasModel[opts_Association : <||>] := Module[{p = Join[<|"PTorr" -> 5., "T" -> 300., "kappa" -> 0.026, "a" -> 1.|>, opts], N0, h, al, be},
   N0 = p["PTorr"] $torr/($kB p["T"]); h = $mAir/(2 $kB p["T"]); al = $mAir N0/Sqrt[Pi h]; be = 2 - p["a"] (1 - Pi/4);
   Join[p, <|"N" -> N0, "h" -> h, "alpha" -> al, "beta" -> be, "Ct" -> p["a"] Sqrt[h]/(5 Sqrt[Pi]), "Cn" -> 2 (2 - p["a"]) Sqrt[h]/(5 Sqrt[Pi]),
     "gamma0/ell^2" -> 16 Sqrt[Pi]/15 p["kappa"] p["T"] Sqrt[h], "xi0/ell^2" -> 4 Pi/3 al (p["a"] + be), "c0/ell^4" -> 4 Pi/3 p["a"] al|>]];

gaussLegendre[n_] := gaussLegendre[n] = Module[{r = NIntegrate`GaussRuleData[n, MachinePrecision]}, {2 r[[1]] - 1, 2 r[[2]]}];
quadGrid[nth_, nph_] := quadGrid[nth, nph] = Module[{x, w},
    {x, w} = gaussLegendre[nth];
    Flatten[Table[{ArcCos[-x[[i]]], 2. Pi (j - 1)/nph, w[[i]] 2. Pi/nph}, {i, nth}, {j, nph}], 1]];   (* weights for d(cos th) dphi *)

makeParticle[Rfun_, ell_, rhoP_, gas_, density_ : (1 &), nth_ : 48, nph_ : 96] := Module[
   {grid = quadGrid[nth, nph], xr, wr, m0 = 0., m1 = {0., 0., 0.}, m2 = ConstantArray[0., {3, 3}], com, I0, Rsp = ConstantArray[0., {6, 9}],
    al = gas["alpha"], aa = gas["a"], be = gas["beta"], Ct = gas["Ct"], Cn = gas["Cn"], kT = gas["kappa"] gas["T"], Id = IdentityMatrix[3],
    dR, geo, X, Y, Z},
   {xr, wr} = gaussLegendre[20];
   (* radius and its Cartesian gradient (analytic), projected onto the sphere *)
   dR = Function[Evaluate[{X, Y, Z}], Evaluate[D[Rfun[{X, Y, Z}], {{X, Y, Z}}]]];
   geo = Table[Module[{th = g[[1]], ph = g[[2]], n, Rv, gr, gs},
       n = {Sin[th] Cos[ph], Sin[th] Sin[ph], Cos[th]}; Rv = Rfun[n]; gr = dR @@ n; gs = gr - n (n.gr);
       {n, Rv, ell^2 (Rv^2 n - Rv gs), g[[3]]}], {g, grid}];
   (* mass properties *)
   Do[Module[{n = gg[[1]], Rv = gg[[2]], wa = gg[[4]]},
      Do[Module[{s = Rv (xr[[k]] + 1)/2, rv, dm}, rv = ell s n; dm = rhoP density[rv/ell] s^2 (Rv/2 wr[[k]]) wa ell^3;
        m0 += dm; m1 += dm rv; m2 += dm ((rv.rv) Id - Outer[Times, rv, rv])], {k, Length[xr]}]], {gg, geo}];
   com = m1/m0; I0 = m2 - m0 ((com.com) Id - Outer[Times, com, com]);
   (* response about the CoM *)
   Do[Module[{S = gg[[3]], Sm, P, KV, Kq, rx, blk},
      Sm = Norm[S]; P = Outer[Times, S, S]/Sm; rx = hat[ell gg[[2]] gg[[1]] - com];
      KV = -al ((aa/2) (Sm Id - P) + be P);
      Kq = Ct (Sm Id - P) + Cn P;
      blk = Join[kT Kq, -KV, KV.rx, 2];
      Rsp += gg[[4]] Join[blk, rx.blk]], {gg, geo}];
   <|"Rfun" -> Rfun, "ell" -> ell, "rhoP" -> rhoP, "gas" -> gas, "M" -> m0, "com" -> com, "I" -> (I0 + Transpose[I0])/2, "Rsp" -> Rsp,
    "blocks" -> blocks[Rsp], "gamma0" -> gas["gamma0/ell^2"] ell^2, "xi0" -> gas["xi0/ell^2"] ell^2, "c0" -> gas["c0/ell^4"] ell^4,
    "G0" -> m0 $gEarth/(gas["gamma0/ell^2"] ell^2)|>];

(* body-frame equations of motion at gradient G: state v, w, u (u = lab Z in the body frame) *)
bodyRHS[pt_, G_, {v_, w_, u_}] := Module[{b = pt["blocks"], M = pt["M"], I0 = pt["I"], F, tau},
   F = G b["f_th"].u - b["f_drag"].v - b["B"].w - M $gEarth u;
   tau = G b["T_th"].u - b["Bt"].v - b["T_drag"].w;
   {F/M - Cross[w, v], Inverse[I0].(tau - Cross[w, I0.w]), -Cross[w, u]}];

(* steady states: w = Om u, v constant in the body frame, levitation u.v = 0; exact within the linear-response model *)
steadyMotions[pt_] := Module[{b = pt["blocks"], ev, vecs, cands, sols = {}},
   {ev, vecs} = Eigensystem[Inverse[b["T_drag"]].b["T_th"]];
   cands = Select[Transpose[{ev, vecs}], Abs[Im[#[[1]]]] <= 10^-8 Max[Abs[ev], 10^-300] && Norm[Im[#[[2]]]] < 10^-8 &];
   Do[Module[{u0 = sg Normalize[Re[c[[2]]]], G0 = pt["G0"], vv, uu, Om, G, eqs, sol, st, J, vars, ee, rates, Rorb},
      vv = Array[Unique["v"] &, 3]; uu = Array[Unique["u"] &, 3]; Om = Unique["Om"]; G = Unique["G"];
      eqs = Join[Sequence @@ bodyRHS[pt, G, {vv, Om uu, uu}][[1 ;; 2]], {uu.vv, uu.uu - 1}];
      sol = Quiet@Check[FindRoot[eqs, Join[Transpose[{vv, ConstantArray[0., 3]}], Transpose[{uu, u0}], {{Om, G0 Re[c[[1]]]}, {G, G0}}],
          MaxIterations -> 200], $Failed];
      If[sol =!= $Failed,
       st = {vv, (Om /. sol) uu, uu} /. sol;
       vars = Array[Unique["s"] &, 9];
       J = D[Flatten[bodyRHS[pt, G /. sol, Partition[vars, 3]]], {vars}] /. Thread[vars -> Flatten[st]];
       ee = SortBy[Eigenvalues[J], Abs];
       rates = Rest[ee];                                                  (* drop the exact zero: |u| is conserved *)
       If[Abs[Om /. sol] < 10^-9, sol = sol /. (Om -> _) :> (Om -> 0.); st[[2]] = {0., 0., 0.}];
       Rorb = If[(Om /. sol) == 0, Infinity, Norm[st[[1]] - (st[[1]].st[[3]]) st[[3]]]/Abs[Om /. sol]];
       AppendTo[sols, <|"u" -> st[[3]], "Omega" -> Om /. sol, "G" -> G /. sol, "gradT (K/cm)" -> (G /. sol) pt["gas"]["T"]/100,
         "v_body" -> st[[1]], "vperp" -> Norm[st[[1]] - (st[[1]].st[[3]]) st[[3]]], "R" -> Rorb, "period" -> If[(Om /. sol) == 0, Infinity, 2 Pi/Abs[Om /. sol]],
         "residual" -> Norm[eqs /. sol], "rates" -> rates, "stable" -> Max[Re[rates]] < 0, "growth" -> Max[Re[rates]]|>]]], {c, cands}, {sg, {1, -1}}];
   DeleteDuplicatesBy[sols, Round[#["u"], 10^-6] &]];

hoverDrift[pt_, u_] := Module[{b = pt["blocks"], Mg = pt["M"] $gEarth, w, z, G, v},
   w = LinearSolve[b["f_drag"], u]; z = LinearSolve[b["f_drag"], b["f_th"].u];
   G = Mg (u.w)/(u.z); v = G z - Mg w;
   <|"G" -> G, "v_body" -> v, "vperp" -> v - (u.v) u, "speed" -> Norm[v - (u.v) u], "torque" -> G b["T_th"].u - b["Bt"].v|>];

(* compiled right-hand side, lengths in ell and times in ms *)
simulate[pt_, G_, q0_, w0_, tEnd_, opts___] := Module[{M0 = pt["M"], I0 = pt["I"], Iinv, Rm = pt["Rsp"], l0 = pt["ell"], t0 = 10.^-3, Zh = {0, 0, 1},
    ys, x, v, q, w, R, wr, rhs, sc, cf, f, y, t, sol},
   Iinv = Inverse[I0];
   ys = Array[Unique["y"] &, 13]; sc = Join[ConstantArray[l0, 3], ConstantArray[l0/t0, 3], {1, 1, 1, 1}, ConstantArray[1/t0, 3]];
   {x, v, q, w} = {sc[[1 ;; 3]] ys[[1 ;; 3]], sc[[4 ;; 6]] ys[[4 ;; 6]], ys[[7 ;; 10]], sc[[11 ;; 13]] ys[[11 ;; 13]]};
   R = quatR[q/Sqrt[q.q]];
   wr = -Rm.Join[Transpose[R].(-G Zh), Transpose[R].v, w];
   rhs = t0 Join[v, R.wr[[1 ;; 3]]/M0 - $gEarth Zh, 1/2 Qmat[q].w, Iinv.(wr[[4 ;; 6]] - Cross[w, I0.w])]/sc;
   cf = Compile[{{yy, _Real, 1}}, Evaluate[rhs /. Thread[ys -> Table[Compile`GetElement[yy, i], {i, 13}]]]];
   f[yy_?(VectorQ[#, NumericQ] &)] := cf[yy];
   sol = NDSolveValue[{y'[t] == f[y[t]], y[0] == Join[{0, 0, 0}, {0, 0, 0}, N[Normalize[q0]], N[w0] t0]}, y, {t, 0, 1000 tEnd},
     opts, MaxSteps -> 10^7, PrecisionGoal -> 8, AccuracyGoal -> 9];
   <|"pt" -> pt, "G" -> G, "tEnd" -> tEnd,
    "state" -> Function[{ts}, Module[{yy = sol[1000 ts]}, Join[l0 yy[[1 ;; 3]], l0/t0 yy[[4 ;; 6]], Normalize[yy[[7 ;; 10]]], yy[[11 ;; 13]]/t0]]]|>];

unwrap[a_List] := FoldList[#1 + Mod[#2 - #1 + Pi, 2 Pi] - Pi &, First[a], Rest[a]];
orbitFromRun[run_, t1_, t2_, n_ : 3000] := Module[{ts = N@Subdivide[t1, t2, n], xy, ctr, rad, ang, slope, tt},
   xy = run["state"][#][[1 ;; 2]] & /@ ts;
   ctr = Mean[xy]; rad = Norm[# - ctr] & /@ xy;
   ang = unwrap[ArcTan @@@ (# - ctr & /@ xy)];
   slope = Coefficient[Fit[Transpose[{ts, ang}], {1, tt}, tt], tt];
   <|"centre" -> ctr, "R" -> Mean[rad], "R rel. spread" -> StandardDeviation[rad]/Mean[rad], "Omega" -> slope, "period" -> 2 Pi/Abs[slope],
    "sense" -> If[slope > 0, "counter-clockwise seen from above", "clockwise seen from above"],
    "vertical drift (m/s)" -> (run["state"][t2][[3]] - run["state"][t1][[3]])/(t2 - t1)|>];

(* smallest principal curvature of the level set F = |r| - R(r/|r|) = 0 *)
convexityMargin[Rfun_, npts_ : 40] := Module[{X, Y, Z, F, gF, hF, g, h, vals},
   F = Sqrt[X^2 + Y^2 + Z^2] - Rfun[{X, Y, Z}/Sqrt[X^2 + Y^2 + Z^2]];
   gF = D[F, {{X, Y, Z}}]; hF = D[F, {{X, Y, Z}, 2}];
   g = Function[Evaluate[{X, Y, Z}], Evaluate[gF]]; h = Function[Evaluate[{X, Y, Z}], Evaluate[hF]];
   vals = Flatten@Table[Module[{n = N[{Sin[t] Cos[p], Sin[t] Sin[p], Cos[t]}], r, gv, P, K},
       r = Rfun[n] n; gv = g @@ r; P = IdentityMatrix[3] - Outer[Times, gv, gv]/(gv.gv);
       K = P.(h @@ r).P/Norm[gv];
       Sort[Eigenvalues[(K + Transpose[K])/2]][[2]]], {t, Pi/(2 npts), Pi, Pi/npts}, {p, 0, 2 Pi, Pi/npts}];
   Min[vals]];

shapeDiagram[Rfun_, com_ : {0, 0, 0}, opts___] := Module[{surf, ax, rng},
   rng = MinMax[Table[Rfun[{Sin[t] Cos[p], Sin[t] Sin[p], Cos[t]}] - 1, {t, 0.05, Pi, Pi/30}, {p, 0, 2 Pi, Pi/30}]];
   rng = Max[Abs[rng], 10^-6];
   surf = ParametricPlot3D[Rfun[{Sin[t] Cos[p], Sin[t] Sin[p], Cos[t]}] {Sin[t] Cos[p], Sin[t] Sin[p], Cos[t]}, {t, 0, Pi}, {p, 0, 2 Pi},
     Mesh -> None, PlotPoints -> 70, ColorFunctionScaling -> False,
     ColorFunction -> Function[{X, Y, Z, t, p}, ColorData["ThermometerColors"][0.5 + 0.5 (Rfun[{Sin[t] Cos[p], Sin[t] Sin[p], Cos[t]}] - 1)/rng]],
     PlotStyle -> Opacity[0.8]];
   ax = {{Red, Arrow[Tube[{{0, 0, 0}, {1.7, 0, 0}}, 0.015]]}, {Darker[Green], Arrow[Tube[{{0, 0, 0}, {0, 1.7, 0}}, 0.015]]},
     {Blue, Arrow[Tube[{{0, 0, 0}, {0, 0, 1.7}}, 0.015]]},
     Text[Style["x", 16, Red], {1.9, 0, 0}], Text[Style["y", 16, Darker[Green]], {0, 1.9, 0}], Text[Style["z", 16, Blue], {0, 0, 1.9}],
     {Black, Sphere[{0, 0, 0}, 0.045]}, {Magenta, Sphere[com, 0.06]}};
   Show[surf, Graphics3D[ax], opts, Boxed -> False, Axes -> False, Lighting -> "Neutral", ImageSize -> 300, ViewPoint -> {2.4, 1.3, 1.1},
    PlotRange -> {{-1.95, 1.95}, {-1.95, 1.95}, {-1.95, 1.95}}]];

bodyMesh[Rfun_, np_ : 30] := Module[{ts = N@Subdivide[0, Pi, np], ps = N@Most@Subdivide[0, 2 Pi, 2 np], pts, rv, polys, idx, rng},
   idx[i_, j_] := (i - 1) (2 np) + Mod[j - 1, 2 np] + 1;
   pts = Flatten[Table[Module[{n = {Sin[t] Cos[p], Sin[t] Sin[p], Cos[t]}}, {Rfun[n] n, Rfun[n] - 1}], {t, ts}, {p, ps}], 1];
   rng = Max[Abs[pts[[All, 2]]], 10^-6];
   polys = Flatten[Table[{idx[i, j], idx[i, j + 1], idx[i + 1, j + 1], idx[i + 1, j]}, {i, np}, {j, 2 np}], 1];
   <|"pts" -> pts[[All, 1]], "polys" -> polys, "cols" -> (ColorData["ThermometerColors"][0.5 + 0.5 #/rng] & /@ pts[[All, 2]])|>];
drawBody[mesh_, x_, R_, scale_] := GraphicsComplex[(x + scale R.#) & /@ mesh["pts"], {EdgeForm[], Polygon[mesh["polys"]]}, VertexColors -> mesh["cols"]];

End[];
EndPackage[];
