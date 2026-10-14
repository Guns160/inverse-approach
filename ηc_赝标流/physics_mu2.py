#!/usr/bin/env python3
"""Save all 2001-alpha diagnostics for one row of the 2D platform scan.

For every (Lambda, M2min) point, rebuild the standard-form SVD once and save
residual norm, H1 penalty norm, pole residues, QOC, singular values, projected
coefficients, and the stable orthogonal residual for all alpha-selection studies.
No noise budget or alpha-selection rule is applied by this worker.
"""
import numpy as np, math, argparse, hashlib, os, tempfile
from pathlib import Path

MC=1.133471090883; G2=0.038; G3=0.013; S_H=14.99
# Only mc and alpha_s are run; G2/G3 remain fixed inputs, not a claim of G3 scale independence.
ALPHA_S=0.302610568157
ALPHAS=np.logspace(-20,0,2001)
M2MIN_VALS=1.5+0.1*np.arange(136)

def sha256_file(path):
    digest=hashlib.sha256()
    with open(path,"rb") as stream:
        for block in iter(lambda:stream.read(1024*1024),b""):
            digest.update(block)
    return digest.hexdigest()

def publish_npz(target,arrays):
    target=Path(target); target.parent.mkdir(parents=True,exist_ok=True)
    fd,name=tempfile.mkstemp(prefix=target.name+f".{os.getpid()}.",suffix=".tmp",dir=target.parent)
    os.close(fd); temporary=Path(name)
    try:
        with temporary.open("wb") as stream:
            np.savez_compressed(stream,**arrays)
            stream.flush(); os.fsync(stream.fileno())
        try:
            os.link(temporary,target)
            return True
        except FileExistsError:
            return False
    finally:
        temporary.unlink(missing_ok=True)

def trapw(x):
    d=np.diff(x); w=np.empty_like(x); w[0]=d[0]/2; w[-1]=d[-1]/2; w[1:-1]=(d[:-1]+d[1:])/2; return w

def h_tridiag(s):
    h=np.diff(s); diag=np.zeros(len(s)); off=h/6-1/h; elem=h/3+1/h; diag[:-1]+=elem; diag[1:]+=elem; return diag,off
def chol_tridiag(d,off):
    r=np.empty_like(d); e=np.empty(len(d)-1); r[0]=math.sqrt(d[0])
    for i in range(len(e)): e[i]=off[i]/r[i]; r[i+1]=math.sqrt(d[i+1]-e[i]**2)
    return r,e
def solve_R(r,e,B):
    X=np.empty_like(B); X[-1]=B[-1]/r[-1]
    for i in range(len(r)-2,-1,-1): X[i]=(B[i]-e[i]*X[i+1])/r[i]
    return X
def solve_RT(r,e,B):
    X=np.empty_like(B); X[0]=B[0]/r[0]
    for i in range(1,len(r)): X[i]=(B[i]-e[i-1]*X[i-1])/r[i]
    return X
def solve_H(r,e,b): return solve_R(r,e,solve_RT(r,e,b))

# ---- A-A channel ----
def aa_rho_pert(s):
    s=np.asarray(s,float); out=np.zeros_like(s); q=4*MC*MC; m=s>q
    v=np.sqrt(1-q/s[m]); out[m]=3*v*(1-v*v/3)/(8*np.pi**2); return out
def aa_ope(L,m2):
    q=4*MC*MC; vel=np.linspace(0,np.sqrt(1-q/L),16001); sg=q/(1-vel*vel); ws=trapw(sg); rw=aa_rho_pert(sg)
    pert=np.empty_like(m2)
    for a in range(0,len(m2),256):
        mm=m2[a:a+256]; pert[a:a+len(mm)]=(np.exp(-sg[None,:]/mm[:,None])/mm[:,None]*rw)@ws
    x=np.linspace(1e-8,1-1e-8,4001); xb=1-x; wx=trapw(x)
    a3=-(45/8)*x*xb*(x**3+xb**3)-(69/72)*x**2*xb**2
    b3=-(3/4)*MC**2*x*xb*(x**4+xb**4)+MC**2*x**2*xb**2*(-(23/12)*(x**2+xb**2)+x*xb/3+2)
    c3=(2/5)*MC**4*x*xb*(x**5+xb**5)
    g2=np.empty_like(m2); g3=np.empty_like(m2)
    for a in range(0,len(m2),128):
        mm=m2[a:a+128,None]; E=np.exp(-MC**2/(x*xb*mm))
        i2=E*(1/(2*mm**2)-MC**2*(x**3+xb**3)/(2*mm**3*x**2*xb**2))
        i3=E*(a3/(2*mm**3*x**3*xb**3)+b3/(6*mm**4*x**4*xb**4)+c3/(24*mm**5*x**5*xb**5))
        g2[a:a+len(mm)]=G2*(i2@wx)/(6*np.pi); g3[a:a+len(mm)]=G3*(i3@wx)/(4*np.pi)**2
    return pert+g2+g3

# ---- P-P channel ----
def pp_li2_series(z):
    z=np.asarray(z,dtype=float); term=z.copy(); total=term.copy()
    for k in range(2,10001):
        term*=z; add=term/(k*k); updated=total+add
        if np.all(np.abs(add)<=4*np.spacing(np.maximum(1.0,np.abs(updated)))): return updated
        total=updated
    raise RuntimeError("Li2 series did not converge")
def pp_li2(x):
    values=np.asarray(x,dtype=float); flat=values.reshape(-1); out=np.empty_like(flat)
    az=flat==0; ao=flat==1; low=(flat>0)&(flat<=0.5); high=(flat>0.5)&(flat<1)
    out[az]=0; out[ao]=np.pi**2/6; out[low]=pp_li2_series(flat[low])
    z=1-flat[high]; out[high]=np.pi**2/6-np.log(flat[high])*np.log1p(-flat[high])-pp_li2_series(z)
    return out.reshape(values.shape)
def pp_finite_function(v):
    lp=np.log1p(v)-np.log1p(-v); lh=np.log1p(v)-np.log(2.0)
    bracket=(np.pi**2/6+lp*lh+2*pp_li2((1-v)/(1+v))-4*pp_li2(v)+pp_li2(v*v)+2*pp_li2((1+v)/2)-2*pp_li2((1-v)/2))
    return ((1+v*v)*bracket+(19/16-3*v+v*v/8+3*v**4/16)*lp+29*v/8-3*v*v/8+6*v*lh-4*v*np.log(v))
def pp_rho_pert_components(s, mc=MC, alpha_s=ALPHA_S):
    values=np.asarray(s,dtype=float); lo=np.zeros_like(values); nlo=np.zeros_like(values)
    threshold=4*mc*mc; above=values>threshold; at=values==threshold
    if np.any(above):
        sa=values[above]; v=np.sqrt(1-threshold/sa)
        lo[above]=3*sa*v/(8*np.pi**2); nlo[above]=sa*alpha_s*pp_finite_function(v)/(2*np.pi**3)
    nlo[at]=threshold*alpha_s/(4*np.pi)
    return lo,nlo
def pp_rho_pert(s):
    lo,nlo=pp_rho_pert_components(s); return lo+nlo
def pp_condensates(m2, mc=MC, x_nodes=4001):
    x=np.linspace(1e-8,1-1e-8,x_nodes); xb=1-x; q=x*xb; wx=trapw(x)
    cubic=mc*mc*(x**3+xb**3); a=mc*x+mc*xb; x3=x**3+xb**3
    b4=64-(416/3)*q-45*x3/q
    b6=((45*a*a/q-(2/3)*mc*mc)*x3-6*mc*mc*(x**4+xb**4)/q-12*mc*mc*(x**3+xb**3)
        -(478/9)*mc*mc*q+(287/9)*a*a-(320/9)*q*mc*mc)
    b8=((6*a*a/q+28*mc*mc)*mc*mc*(x**4+xb**4)-4*mc**4*(x**4+xb**4)
        +(16/5)*mc**4*(x**5+xb**5)/q+(30*a*a+88*mc*mc*q)*mc*mc*(x*x+xb*xb)
        -(112/3)*a**4+(64/3)*mc**4*q+128*mc**4*q*q)
    g2=np.empty_like(m2); g3=np.empty_like(m2)
    for start in range(0,len(m2),128):
        mm=m2[start:start+128,None]; expo=np.exp(-mc*mc/(q*mm))
        term4_g2=(-6*mc*mc-2*cubic/q)/(mm**2*x**2*xb**2)
        term6_g2=a*a*cubic/(mm**3*x**4*xb**4)
        i2=expo*(term4_g2+term6_g2)
        i3=expo*(b4/(mm**2*x**2*xb**2)+b6/(2*mm**3*x**3*xb**3)+b8/(6*mm**4*x**4*xb**4)
            -(2/15)*a*a*mc**4*(x**5+xb**5)/(mm**5*x**6*xb**6))
        g2[start:start+len(mm)]=G2*(i2@wx)/(48*np.pi)
        g3[start:start+len(mm)]=G3*(i3@wx)/(32*(4*np.pi)**2)
    return g2,g3
def pp_ope(upper,m2):
    threshold=4*MC*MC
    velocity=np.linspace(0,np.sqrt(1-threshold/upper),16001)
    s=threshold/(1-velocity*velocity); ws=trapw(s)
    rho_lo,rho_nlo=pp_rho_pert_components(s)
    lo=np.empty_like(m2); nlo=np.empty_like(m2)
    for start in range(0,len(m2),256):
        mm=m2[start:start+256]; weight=np.exp(-s[None,:]/mm[:,None])/mm[:,None]
        lo[start:start+len(mm)]=(weight*rho_lo)@ws
        nlo[start:start+len(mm)]=(weight*rho_nlo)@ws
    g2,g3=pp_condensates(m2)
    return lo+nlo+g2+g3

def sgrid_aa(L):
    n=int(np.floor((L-S_H)/.01+1e-10)); s=S_H+.01*np.arange(n+1)
    if L-s[-1]>1e-10: return np.append(s,L)
    s[-1]=L; return s
def sgrid_pp(upper):
    n=int(np.floor((upper-S_H)/0.01+1e-10)); s=S_H+0.01*np.arange(n+1)
    if upper-s[-1]>1e-10:
        s=np.append(s,upper)
    else:
        s[-1]=upper
    if np.min(np.diff(s))<0.009999:
        raise AssertionError(f"P-P s grid contains a spurious small cell at Lambda={upper}")
    return s

def hats(m2,s):
    K=np.zeros((len(m2),len(s)))
    for j,h in enumerate(np.diff(s)):
        q=h/m2; st=np.exp(-s[j]/m2); small=np.abs(q)<1e-4
        asc=np.empty_like(q); des=np.empty_like(q)
        z=q[small]; asc[small]=.5-z/3+z*z/8-z**3/30+z**4/144; des[small]=.5-z/6+z*z/24-z**3/120+z**4/720
        z=q[~small]; e=np.exp(-z); asc[~small]=(1-(1+z)*e)/z**2; des[~small]=(z-1+e)/z**2
        fac=h/m2; K[:,j]+=st*fac*des; K[:,j+1]+=st*fac*asc
    return K

def run_lambda(channel, Lam, outdir, m2min_values=None, run_id="", point_id="", worker_sha256="", manifest_sha256="", attempt_id="", scheduler_job_id="", submission_id=""):
    if worker_sha256:
        actual_worker_sha256=sha256_file(__file__)
        if actual_worker_sha256 != worker_sha256:
            raise ValueError(f"worker SHA-256 mismatch: expected {worker_sha256}, got {actual_worker_sha256}")
    masses={"1delta":(2.9841,),"2delta":(2.9841,3.51067),"3delta":(2.9841,3.51067,3.6378),"1pole":(2.9841,),"2pole":(2.9841,3.6378)}[channel]
    selected_m2min=M2MIN_VALS if m2min_values is None else np.asarray(m2min_values,dtype=float)
    raw_indices=(selected_m2min-1.5)/0.1
    indices=np.rint(raw_indices).astype(int)
    valid=np.isclose(raw_indices,indices,rtol=0.0,atol=1e-10)
    inside=(indices>=0)&(indices<136)
    if np.any(~np.isfinite(selected_m2min)) or np.any(~valid) or np.any(~inside):
        raise ValueError("M2min values must lie on the 1.5:0.1:15.0 grid")
    selected_m2min=1.5+0.1*indices
    is_pp = channel in ("1pole","2pole")
    rho_pert = pp_rho_pert if is_pp else aa_rho_pert
    ope = pp_ope if is_pp else aa_ope
    sgrid = sgrid_pp if is_pp else sgrid_aa
    m2all=np.arange(1.5,80.0001,0.01)
    s=sgrid(Lam); G=ope(Lam,m2all); HAT=hats(m2all,s)
    pcols=np.column_stack([np.exp(-m*m/m2all)/m2all for m in masses]); KF=np.column_stack((pcols,HAT))
    pab=len(masses); N=len(s)
    boundary=np.array([0.,rho_pert(np.array([s[-1]]))[0]])
    hd,ho=h_tridiag(s); rd,re=chol_tridiag(hd[1:-1],ho[1:-1])
    bif=np.zeros((N-2,2)); bif[0,0]=ho[0]; bif[-1,1]=ho[-1]; q=solve_H(rd,re,bif@boundary)
    fixed=np.array([pab,KF.shape[1]-1]); free=np.r_[np.arange(pab),np.arange(pab+1,KF.shape[1]-1)]
    iw=trapw(s)[1:-1]
    boundary_penalty=float(hd[0]*boundary[0]**2+hd[-1]*boundary[1]**2-q@(bif@boundary))
    nm=len(selected_m2min); na=len(ALPHAS)
    residues=np.empty((nm,pab,na)); residual=np.empty((nm,na)); penalty=np.empty((nm,na))
    qoc=np.empty((nm,na)); qoc_full_euclidean=np.empty((nm,na)); qoc_full_l2=np.empty((nm,na))
    solution_euclidean_norm_sq=np.empty((nm,na)); solution_l2_norm_sq=np.empty((nm,na))
    effective_dof=np.empty((nm,na)); data_norm=np.empty(nm); n_data=np.empty(nm,dtype=np.int64)
    orthogonal_residual_sq=np.empty(nm); direct_relative_error=np.empty((nm,6))
    stationarity_relative_error=np.empty((nm,6)); pole_orthogonality_backward_error=np.empty((nm,6))
    penalty_relative_error=np.empty((nm,6)); qoc_fd_relative_error=np.empty((nm,4)); pole_condition=np.empty(nm)
    singular_parts=[]; beta_parts=[]; offsets=[0]
    sample_indices=np.array([0,400,800,1200,1600,2000])
    qoc_fd_indices=np.array([400,800,1200,1600])
    for j,M2min in enumerate(selected_m2min):
        i1=round((M2min-1.5)/0.01)
        K=KF[i1:]; g=G[i1:]; m2=m2all[i1:]; w_m2=trapw(m2)
        sw=np.sqrt(w_m2); k0=sw[:,None]*K[:,free]; b0=sw*(g-K[:,fixed]@boundary)+k0[:,pab:]@q
        C=k0[:,:pab]; D=k0[:,pab:]; Q,R=np.linalg.qr(C,mode='reduced'); bp=b0-Q@(Q.T@b0); Dp=D-Q@(Q.T@D)
        A=solve_RT(rd,re,Dp.T).T; U,sv,Vh=np.linalg.svd(A,full_matrices=False); beta=U.T@bp
        den=sv[:,None]**2+ALPHAS
        orthogonal_remainder=bp-U@beta
        nullpart=float(orthogonal_remainder@orthogonal_remainder)
        residual[j]=np.sqrt(np.sum(((ALPHAS[None,:]*beta[:,None])/den)**2,axis=0)+nullpart)
        coefficients=(sv*beta)[:,None]/den
        y=Vh.T@coefficients
        t=solve_R(rd,re,y)
        penalty[j]=np.sqrt(np.maximum(np.sum(coefficients*coefficients,axis=0)+boundary_penalty,0.0))
        derivative_coefficients=-(sv*beta)[:,None]/den**2
        dy=Vh.T@derivative_coefficients
        dt=solve_R(rd,re,dy)
        derivative_log_factor=np.log(10)*ALPHAS
        qoc[j]=np.sum(iw[:,None]*(dt*derivative_log_factor)**2,axis=0)
        effective_dof[j]=np.sum(sv[:,None]**2/den,axis=0)+pab
        qt_d=Q.T@D
        residues[j]=np.linalg.solve(R,Q.T@b0[:,None]-qt_d@t)
        residue_derivative=np.linalg.solve(R,-qt_d@dt)
        continuum_internal=t-q[:,None]
        continuum_euclidean_sq=np.sum(continuum_internal**2,axis=0)+boundary@boundary
        continuum_l2_sq=np.sum(iw[:,None]*continuum_internal**2,axis=0)+trapw(s)[[0,-1]]@(boundary**2)
        solution_euclidean_norm_sq[j]=np.sum(residues[j]**2,axis=0)+continuum_euclidean_sq
        solution_l2_norm_sq[j]=np.sum(residues[j]**2,axis=0)+continuum_l2_sq
        scaled_residue_derivative=residue_derivative*derivative_log_factor
        scaled_continuum_derivative=dt*derivative_log_factor
        qoc_full_euclidean[j]=np.sum(scaled_residue_derivative**2,axis=0)+np.sum(scaled_continuum_derivative**2,axis=0)
        qoc_full_l2[j]=np.sum(scaled_residue_derivative**2,axis=0)+qoc[j]
        data_norm[j]=math.sqrt(float(np.sum(w_m2*g*g))); n_data[j]=len(g)
        orthogonal_residual_sq[j]=nullpart; pole_condition[j]=np.linalg.cond(C)
        singular_parts.append(sv); beta_parts.append(beta); offsets.append(offsets[-1]+len(sv))
        for k_index,k in enumerate(sample_indices):
            cont=np.r_[boundary[0],continuum_internal[:,k],boundary[1]]
            zsol=np.r_[residues[j,:,k],cont]
            transformed_residual=C@residues[j,:,k]+D@t[:,k]-b0
            direct=math.sqrt(float(np.sum(w_m2*(K@zsol-g)**2)))
            direct_relative_error[j,k_index]=abs(direct-residual[j,k])/max(direct,np.finfo(float).tiny)
            pole_backward_scale=np.linalg.norm(C)*(np.linalg.norm(C)*np.linalg.norm(residues[j,:,k])+np.linalg.norm(D)*np.linalg.norm(t[:,k])+np.linalg.norm(b0))
            pole_orthogonality_backward_error[j,k_index]=np.linalg.norm(C.T@transformed_residual)/max(pole_backward_scale,np.finfo(float).tiny)
            stationarity=A.T@(A@y[:,k]-bp)+ALPHAS[k]*y[:,k]
            stationarity_relative_error[j,k_index]=np.linalg.norm(stationarity)/max(np.linalg.norm(A.T@bp),np.finfo(float).tiny)
            direct_penalty=math.sqrt(max(float(np.sum(hd*cont*cont)+2*np.sum(ho*cont[:-1]*cont[1:])),0.0))
            penalty_relative_error[j,k_index]=abs(direct_penalty-penalty[j,k])/max(direct_penalty,np.finfo(float).tiny)
        log_alpha=np.log10(ALPHAS)
        for k_index,k in enumerate(qoc_fd_indices):
            fd_residue=(residues[j,:,k+1]-residues[j,:,k-1])/(log_alpha[k+1]-log_alpha[k-1])
            fd_continuum=(continuum_internal[:,k+1]-continuum_internal[:,k-1])/(log_alpha[k+1]-log_alpha[k-1])
            fd_qoc=np.sum(fd_residue**2)+np.sum(iw*fd_continuum**2)
            qoc_fd_relative_error[j,k_index]=abs(fd_qoc-qoc_full_l2[j,k])/max(qoc_full_l2[j,k],np.finfo(float).tiny)
        arrays=(residual[j],penalty[j],qoc[j],qoc_full_euclidean[j],qoc_full_l2[j],solution_euclidean_norm_sq[j],solution_l2_norm_sq[j],effective_dof[j],residues[j],sv,beta)
        if any(not np.all(np.isfinite(values)) for values in arrays):
            raise AssertionError(f"non-finite output at {channel} Lambda={Lam} M2min={M2min}")
        if np.any(np.diff(residual[j])<-1e-12*np.maximum(residual[j][1:],1.0)) or np.any(np.diff(penalty[j])>1e-12*np.maximum(penalty[j][:-1],1.0)):
            raise AssertionError(f"non-monotone residual or penalty at {channel} Lambda={Lam} M2min={M2min}")
        diagnostic_maxima={
            "direct_residual":float(np.max(direct_relative_error[j])),
            "stationarity":float(np.max(stationarity_relative_error[j])),
            "pole_orthogonality_backward_error":float(np.max(pole_orthogonality_backward_error[j])),
            "penalty":float(np.max(penalty_relative_error[j])),
            "qoc_finite_difference":float(np.max(qoc_fd_relative_error[j])),
        }
        diagnostic_limits={"direct_residual":1e-6,"stationarity":1e-8,"pole_orthogonality_backward_error":1e-12,"penalty":1e-8,"qoc_finite_difference":1e-3}
        failed={name:value for name,value in diagnostic_maxima.items() if value>diagnostic_limits[name]}
        if failed:
            raise AssertionError(f"solver diagnostic failed at {channel} Lambda={Lam} M2min={M2min}: {failed}")
    out=Path(outdir); out.mkdir(parents=True,exist_ok=True)
    target=out/f"{channel}_L{Lam:.1f}.npz"
    arrays=dict(
        schema_version=np.array("full-alpha-v2"),run_id=np.array(run_id),point_id=np.array(point_id),worker_sha256=np.array(worker_sha256),
        manifest_sha256=np.array(manifest_sha256),attempt_id=np.array(attempt_id),scheduler_job_id=np.array(scheduler_job_id),submission_id=np.array(submission_id),
        channel=np.array(channel),Lambda=np.array(Lam),M2max=np.array(80.0),M2step=np.array(0.01),
        M2min=selected_m2min,alpha=ALPHAS,masses=np.array(masses),s_grid=s,
        residues=residues,residual=residual,penalty=penalty,qoc=qoc,
        qoc_full_euclidean=qoc_full_euclidean,qoc_full_l2=qoc_full_l2,
        solution_euclidean_norm_sq=solution_euclidean_norm_sq,solution_l2_norm_sq=solution_l2_norm_sq,
        effective_dof=effective_dof,data_norm=data_norm,n_data=n_data,
        orthogonal_residual_sq=orthogonal_residual_sq,pole_condition=pole_condition,
        singular_values=np.concatenate(singular_parts),beta=np.concatenate(beta_parts),
        svd_offsets=np.array(offsets,dtype=np.int64),
        direct_sample_indices=sample_indices,direct_residual_relative_error=direct_relative_error,
        stationarity_relative_error=stationarity_relative_error,
        pole_orthogonality_backward_error=pole_orthogonality_backward_error,
        penalty_relative_error=penalty_relative_error,
        qoc_fd_indices=qoc_fd_indices,qoc_fd_relative_error=qoc_fd_relative_error)
    published=publish_npz(target,arrays)
    if not published:
        with np.load(target,allow_pickle=False) as existing:
            immutable_provenance=tuple(str(existing[key]) for key in ("run_id","point_id","worker_sha256","manifest_sha256"))
        expected_immutable=(run_id,point_id,worker_sha256,manifest_sha256)
        if immutable_provenance != expected_immutable:
            raise RuntimeError(f"existing result has different provenance: {target}")
    max_diagnostic=max(float(np.max(values)) for values in (direct_relative_error,stationarity_relative_error,pole_orthogonality_backward_error,penalty_relative_error,qoc_fd_relative_error))
    disposition="DONE" if published else "REUSED"
    print(f"{disposition} {channel} L={Lam} points={nm} alphas={na} max_diagnostic={max_diagnostic:.3e}",flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("channel",choices=("1delta","2delta","3delta","1pole","2pole")); ap.add_argument("Lambda",type=float); ap.add_argument("outdir")
    ap.add_argument("--m2min",type=float,nargs="+",default=None,help="optional smoke/pilot subset of the production M2min grid")
    ap.add_argument("--run-id",default=""); ap.add_argument("--point-id",default=""); ap.add_argument("--worker-sha256",default="")
    ap.add_argument("--manifest-sha256",default=""); ap.add_argument("--attempt-id",default=""); ap.add_argument("--scheduler-job-id",default=""); ap.add_argument("--submission-id",default="")
    a=ap.parse_args(); run_lambda(a.channel,a.Lambda,a.outdir,a.m2min,a.run_id,a.point_id,a.worker_sha256,a.manifest_sha256,a.attempt_id,a.scheduler_job_id,a.submission_id)
