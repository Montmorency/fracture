import sys

import numpy as np

import ase.io

from ase.lattice.cubic import BodyCenteredCubic
from quippy import Potential, BOHR

class TightBindingPot(object):
  """
  Interface to the lmto-tb code. Contains the input
  file and a routine for writing the ctrl.{atom} file.
  relax= 5 for FLETCHER/G
  """
  def __init__(self, alat=1.00, nk=1, n_sup=1):
    self.alat = alat   # lattice constant, angstrom
    self.nk   = nk     # number of k-points along cubic axes
    self.n    = n_sup  # we use an n x n x n supercell of primitive cell
#We've removed carbon from this input file.
    self.template = """HEADER auto-generated control file
%% const fp=0 cpl=1 xcf=4 gga=3
%% const pureFe=0 PRBmodel=0 CModel2=0 scale=1.0
%% const epsilon=1 sign=epsilon?-1:1 xFe=1/3
%% const nbas=%(NBAS)d nspec=%(NSPEC)d
%% const sd=1 ovlp=sd ul=1 u1=0 tbu=0 io=0 nitq=50
%% const verb=31 so=0 nsp=2 tetra=0 metal={fp?3:1} width=0.002 N=1
%% const au=0.529177 NdFe=6.8
%% const beta=0.06 kmix=300 nav=0
%% const show=0 mpol=0 nx=5
%% const dyn=0 temp=300 taup=10 taub=100 time=100000 tstep=5
%% const hess=F relax=0 nit=60 xtol=1d-3 gtol=5d-4 step=0.01 nkill=100 nitf=50
%% const fs=0.048377 K=1/0.6333328d-5 amass=1.09716d-3
VERS    TB=10 LM=7 FP=7 ASA=7
IO      SHOW={show} HELP=F VERBOS={verb} WKP=F
CONST   nit=100 conv=1d-4 qtol=1d-3 pair=F
        V0at=79.76508 V0bcc=V0at*2 V0fcc=V0at*4
        vfrac0=1 vfrac=vfrac0 Vbcc=vfrac*V0bcc Vfcc=vfrac*V0fcc
        q0=sqrt(8/3) q=q0 qJang=1.57693 qexp=1.582 qTB=1.614
        ahcp=(Vfcc/(sqrt(3)*q))^(1/3)        
        abcc=Vbcc^(1/3) afcc=Vfcc^(1/3) aeps=sqrt(3)*ahcp
        alat=%(ALAT)f 
        nk=%(NK)d mull=-1 bzj=1 ewtol=1d-6 ef0=0.643
        R=2.2 RC=0.8
%% const NdFe=7 V0at=79.76508 V0bcc=V0at*2 a0bcc=V0bcc^(1/3) 
        V0bcc={V0bcc} a0bcc={a0bcc}
# cut-offs in a.u.
%% ifdef scale
        fixscale=abcc
%% else
        fixscale=a0bcc
%% endif
        r1HFau=0.8*fixscale rcHFau=2*fixscale
        r1HFppau=0.75*fixscale rcHFppau=0.95*fixscale
        r1ssau=1.1*fixscale rcssau=2*fixscale 
        r1sdau=1.1*fixscale rcsdau=2*fixscale 
        r1ddau={PRBModel?1.1:0.9}*fixscale rcddau=1.4*fixscale 
        r1ppau={PRBModel?1.1:0.9}*fixscale rcppau=1.4*fixscale
        cutmod={PRBModel?1:2} 
%% ifdef sd
        rmaxhau=3*fixscale
%% else
        rmaxhau=1.4*fixscale
%% endif
        r1HHau=1.1*fixscale rcHHau=1.4*fixscale
# on-site terms
        q0s={sd?1:0} q0p=0 q0dFe={sd?7:NdFe}
        esFe=0.15 epFe=0.45 edFe=0 momFe={nsp==1?0:2}
        U=1 Us=U Up=U UdFe=U stniFe={sd?0.055:0.05}
        q0sC=2 q0pC=2 
        esc=-0.467663945 epc=0.08275667052 UC=1.238348985 
        q0H=1 esH=-0.085 UH=1.2 momH={nsp==1?0:1}
        spp=0 ppd=0 sdd=0 pdp=0 ddd=0 pdf=0 ddg=0
# hopping integrals
        qsp=1 qpp=1 qpd=1
        fsp=0 fpp=0 fpd=0
        odds=0 oddp=0 oddd=0 opp=0 osp=0 opd=0
# Fe-Fe
        r0ff=0.5*sqrt(3)*V0bcc^(1/3)
        qdds0=1 qddp0=1 qddd0=1 qss0=0.3 qsd0=0.57
        fdd0=0.65  fss0=-0.35 fsd0=-0.5
        fdds0=-fdd0*6 fddp0=fdd0*4.914539385 fddd0=fdd0*-2.232504465
        qdds=0.9 qddp=0.9 qddd=0.9 qss=qss0 qsd=0.3
        hddsr0=fdds0*exp(-qdds0*r0ff)
        hddpr0=fddp0*exp(-qddp0*r0ff)
        hdddr0=fddd0*exp(-qddd0*r0ff)
        hssr0=fss0*exp(-qss0*r0ff)
        hsdr0=fsd0*exp(-qsd0*r0ff)
        fdds=hddsr0*exp(qdds*r0ff)
        fddp=hddpr0*exp(qddp*r0ff)
        fddd=hdddr0*exp(qddd*r0ff)
        fss=hssr0*exp(qss*r0ff)
        fsd=hsdr0*exp(qsd*r0ff)
        qoss0=qss qosd0=qsd
        oss0=0.45 osd0=0.5
        ossr0=oss0*exp(-qoss0*r0ff)
        osdr0=osd0*exp(-qosd0*r0ff)
        qoss=qoss0 qosd=qosd0
        oss=ossr0*exp(qoss*r0ff)
        osd=osdr0*exp(qosd*r0ff)
# Fe-H
        r0HF=1.453500953
        qHFss0=0.592 qHFsd0=0.601
        fHFss0=-0.8365709269 fHFsd0=-0.5041736305
        hHFssr0=fHFss0*exp(-qHFss0*r0HF)
        hHFsdr0=fHFsd0*exp(-qHFsd0*r0HF)
        qHFss=0.7762840122 qHFsd=0.4544987809
        fHFss=hHFssr0*exp(qHFss*r0HF)
        fHFsd=hHFsdr0*exp(qHFsd*r0HF)

        ofacHFss=0.4676030053 ofacHFsd=0.399106628
        qoHFss0=0.552 qoHFsd0=0.412
        oHFss0=-ofacHFss*fHFss0 oHFsd0=-ofacHFsd*fHFsd0
        oHFssr0=oHFss0*exp(-qoHFss0*r0HF)
        oHFsdr0=oHFsd0*exp(-qoHFsd0*r0HF)
        qoHFss=0.2863260142 qoHFsd=0.473014452
        oHFss=oHFssr0*exp(qoHFss*r0HF)
        oHFsd=oHFsdr0*exp(qoHFsd*r0HF)

        fHHss=0 qHHss=0.5 
        fHFsp=0 qHFsp=0
        oHHss=0 oHFsp=0 qoHFsp=0
# Fe-Fe pair potential
%% ifdef sd
%%  ifdef scale
%%   ifdef PRBModel
         b0=536 m0=0 p0=1.49 b1=-371.2 m1=0 p1=1.413111 
%%   else
         b0=665.6 m0=0 p0=1.408429 b1=-536.8 m1=0 p1=1.362971
%%   endif
%%  else
         b0=698.666667 m0=0 p0=1.52 b1=-517.466667 m1=0 p1=1.4576
%%  endif
%% else
%%  ifdef scale
         b0=682.8 m0=0 p0=1.5165 b1=-466.8 m1=0 p1=1.435
%%  else
         b0=683.1 m0=0 p0=1.5376 b1=-459.5 m1=0 p1=1.4544
%%  endif
%% endif
# Fe-H pair potential
        qHF=2.69224661 nHF=-1 bHFfac=0.2995633136
        bHF0=1000 bHF=bHF0*bHFfac
%% ifdef scale
        ascale = alat
%% else
        ascale = 1
%% endif
        rmaxh=rmaxhau/ascale
        r1HF=r1HFau/ascale rcHF=rcHFau/ascale
        r1ss=r1ssau/ascale rcss=rcssau/ascale r1sd=r1sdau/ascale rcsd=rcsdau/ascale 
        r1dd=r1ddau/ascale rcdd=rcddau/ascale r1pp=r1ppau/ascale rcpp=rcppau/ascale
        r1HFpp=r1HFppau/ascale rcHFpp=rcHFppau/ascale
%% ifdef mpol
        spp=0 ppd=0 sdd=0 pdp=1 ddd=3 pdf=0 ddg=6
%% else
        spp=0 ppd=0 sdd=0 pdp=0 ddd=0 pdf=0 ddg=0
%% endif
        force=1 pv=1 mol=0
ITER    CONV=conv CONVC=qtol NIT={nitq} MIX=A{nx},b={beta}
DYN
%% if dyn==1|dyn==2|dyn==3
        MD[MODE={dyn} TSTEP={tstep/fs} TEMP={temp/K} TAUP={taup/fs}
           TIME={time/fs} TAUB={taub/fs}]
%% elseif relax>0
        MSTAT[MODE={relax} HESS={hess} XTOL={xtol} GTOL={gtol}
              STEP={step} NKILL={nkill}] NIT={nitf}
%% endif
STRUC   NBAS={nbas} NSPEC={nspec} NL={fp?5:3} ALAT=alat
        PLAT= %(PLAT_STR)s
SITE
        %(SITE_STR)s
BZ      NKABC=nk TETRA={tetra} METAL={metal}
        EFMAX=2 EF0=ef0 DELEF=0.01 N={N} W={width}
        NPTS=5001 BZJOB=bzj SAVDOS=T NOINV=F
        INVIT=F MULL=mull DOS=-4.5 1 EFMAX=2
HAM     NSPIN={nsp} ELIND=-0.8 GMAX=gmax REL=T SO={so}
        XCFUN={xcf} GGA={gga} FORCES=12
        PWMODE=pwmode PWEMIN=1 PWEMAX=pwemax OVEPS=oveps
SPEC    
        ATOM=Fe Z=26 R=R I=stniFe A=0.025 AMASS=55.845/{amass}
        IDU= 0 0 0 0 UH= 0 0 0 0  JH=stniFe stniFe stniFe stniFe 
        COLOUR=0.1 0.1 0.1  RADIUS=0.5
%%  ifdef fp
        LMX=2 LMXA=4 KMXA=4 LFOCA=1
        RSMH=0.95 0.95 0.95 0 EH=-0.1 -0.1 -0.1 -0.1
        RSMH2=0.95 0.95 0.95 EH2=-1.1 -1.1 -1.1
        Q=2 0 6 MMOM=0 0 2 PZ=0 {cpl?3.9:0}
%%  else
        IDXDN={sd?1:3} 3 1 QPOL= spp ppd sdd pdp ddd pdf ddg 0 0 0
%%  endif
        ATOM=H Z=1 R=RH I=stniH A=0.025 AMASS=1.00794/{amass}
        LMX=2 LMXL=2 LMXA=2
        IDU= 0 0 0 0 UH= 0 0 0 0  JH=stniH stniH stniH stniH
        RSMH=RH/1.5 RH/1.5 RH/1.5 0 EH=-0.1 -0.1 -0.1 -0.1
        MMOM=1 0 0
        COLOUR=0.9 0.2 0.2 RADIUS=0.2
%% ifndef fp
        IDXDN=1 3 3 
%% endif
START   CNTROL=T
        ATOM=Fe   P= 4 4 3 4 4 3
                  Q= q0s/{nsp}            esFe   Us
                     q0p/{nsp}            epFe   Up
                     (q0dFe+momFe)/{nsp}  edFe  UdFe
                     q0s/{nsp}            esFe   Us
                     q0p/{nsp}            epFe   Up
                     (q0dFe-momFe)/{nsp}  edFe  UdFe
        ATOM=H    P= 1 2 3 1 2 3
                  Q= (q0H+momH)/{nsp}     esH    UH
                     0                    0      0
                     0                    0      0
                     (q0H-momH)/{nsp}     esH   UH
                     0                    0      0
                     0                    0      0
OPTIONS ASA[ADNF[0] NSPH[0] TWOC[0] CCOR[1]]
ME
        2
        Fe Fe MEMODE=2 PPMODE=10 POLY=5 CUTMOD=cutmod CUTPP=r1pp rcpp
            | fss fsp fpp -fpp/2 fsd fpd -fpd/sqrt(3) 
                                                 fdds fddp fddd
        DECAY=qss qsp qpp qpp    qsd qpd  qpd    qdds qddp qddd
        CUT=r1ss rcss 0 0 0 0 0 0 r1sd rcsd 0 0 0 0 r1dd rcdd r1dd rcdd r1dd rcdd
            @ oss osp opp -opp/2 osd opd -opd/sqrt(3) 
                                                 odds oddp oddd
        DECAY=qoss qsp qpp qpp    qosd qpd  qpd    qdds qddp qddd
        CUT=r1ss rcss 0 0 0 0 0 0 r1sd rcsd 0 0 0 0 r1dd rcdd r1dd rcdd r1dd rcdd
            ! b0 m0 p0   b1 m1 p1   0 0 0
        Fe H MEMODE=2 PPMODE=10 POLY=5 CUTMOD=cutmod CUTPP=r1HFpp rcHFpp
            | fHFss fHFsp  0  0 fHFsd 0 0 0 0 0
        DECAY=qHFss qHFsp  0  0 qHFsd 0 0 0 0 0
        CUT= r1HF rcHF r1HF rcHF 0 0 0 0 r1HF rcHF 0 0 0 0 0 0 0 0 0 0
            @ oHFss oHFsp  0  0 oHFsd 0 0 0 0 0
        DECAY=qoHFss qoHFsp  0  0 qoHFsd 0 0 0 0 0
        CUT= r1HF rcHF r1HF rcHF 0 0 0 0 r1HF rcHF 0 0 0 0 0 0 0 0 0 0
            ! bHF nHF qHF  0 0 0   0 0 0
        H Fe MEMODE=2 PPMODE=10 POLY=5 CUTMOD=cutmod CUTPP=r1HFpp rcHFpp
            | fHFss fHFsp  0  0 fHFsd 0 0 0 0 0
        DECAY=qHFss qHFsp  0  0 qHFsd 0 0 0 0 0
        CUT= r1HF rcHF r1HF rcHF 0 0 0 0 r1HF rcHF 0 0 0 0 0 0 0 0 0 0
            @ oHFss oHFsp  0  0 oHFsd 0 0 0 0 0
        DECAY=qoHFss qoHFsp  0  0 qoHFsd 0 0 0 0 0
        CUT= r1HF rcHF r1HF rcHF 0 0 0 0 r1HF rcHF 0 0 0 0 0 0 0 0 0 0
            ! bHF nHF qHF  0 0 0  0 0 0
        H H MEMODE=2 PPMODE=0
            | 0 0 0 0 0 0 0 0 0 0 
        DECAY=0 0 0 0 0 0 0 0 0 0
        CUT=  0 0  0 0  0 0  0 0  0 0  0 0  
              0 0  0 0  0 0  0 0  0 0  0 0  
            @ 0 0 0 0 0 0 0 0 0 0
        DECAY=0 0 0 0 0 0 0 0 0 0
        CUT=  0 0  0 0  0 0  0 0  0 0  0 0  
              0 0  0 0  0 0  0 0  0 0  0 0  
            ! 0 0 0  0 0 0    0 0 0
TB      FORCES=force EVDISC=T RMAXH=rmaxh TRH=T RHO=T 3PV=pv
        MOL=mol GAMMA=F PAIR=pair SCALE={scale}
        UL={ul} IODEL={io} OVLP={ovlp} TBU={tbu} NOUAVG={nav} U1={u1}
EWALD   TOL=ewtol NKDMX=1999 NKRMX=1999
OPTIONS ASA[ADNF[0] NSPH[0] TWOC[0] CCOR[1]]"""

  def write_control_file(self, filename, atoms, nspec=2):
    """
    Generate control file for TBE code from template and Atoms object
    """
    plat_str = ''
    for i in range(3):
        plat_str += ("%12.8f"*3) % tuple(atoms.cell[i, :]/self.alat) + '\n      '
    plat_str = plat_str.rstrip()
    site_str = ''
    for sym, pos in zip(atoms.get_chemical_symbols(),
                        atoms.get_positions()):
        site_str += ('ATOM=%-2s POS=%12.8f %12.8f %12.8f\n        ' %
                     ((sym,) + tuple(pos/self.alat)))
    site_str = site_str.rstrip()
    dict = {
        'ALAT': self.alat/BOHR,
        'NBAS': len(atoms),
        'NK': self.nk,
        'PLAT_STR': plat_str,
        'SITE_STR': site_str,
        'NSPEC': nspec
        }    
    ctrl_file = open(filename, 'w')
    ctrl_file.write(self.template % dict)
    ctrl_file.close()
