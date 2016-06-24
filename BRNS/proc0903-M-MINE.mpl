# 
restart ;
iswitch:=0;
# 
if searchtext("Maple 6", kernelopts(version)) = 1 then
  libname := "C:\\Program Files\\BRNSPackage\\acglib-Maple6", "C:\\Program Files\\Maple 6\\lib", libname ;
else
  libname := "/home/dh12423/Documents/5_BRNS/BRNS_full_for_OMEN_15062016/acglib-Maple10", "/opt/maple17/lib","/home/dh12423/Documents/5_BRNS/BRNS_full_for_OMEN_15062016/share/share", libname ;
end if; 
with(acg) ; with(student) : with(linalg) : with(share) :
try
 with(macrofor) : 
catch : ;
end try;
# 
currentdir("/home/dh12423/Documents/5_BRNS/BRNS_full_for_OMEN_15062016"):
read "spread.m" ;
# 
dir_f:="/home/dh12423/Documents/5_BRNS/BRNS_full_for_OMEN_15062016/GeneratedFortranFiles":
# If creating a DLL, assign dummy values to all undefined variables
# 

if tot_time = 'tot_time' then tot_time := 1.0 fi ;
if depth_max = 'depth_max' then depth_max := 1.0 fi ;
if val = 'val' then val := 1.0 fi ;
if viq = 'viq' then viq := 0.0 fi ;
if vq0 = 'vq0' then vq0 := 1.0 fi ;
if viw = 'viw' then viw := 0.0 fi ;
if vw0 = 'vw0' then vw0 := 0.0 fi ;
if viDb = 'viDb' then viDb := 0.0 fi ;
if vDb0 = 'vDb0' then vDb0 := 0.0 fi ;
if vipor = 'vipor' then vipor := 0.0 fi ;
if vpor0 = 'vpor0' then vpor0 := 0.5 fi ;
if viarea = 'viarea' then viarea := 0.0 fi ;
if varea0 = 'varea0' then varea0 := 1.0 fi ;
if T_C = 'T_C' then T_C := 20.0 fi ;
if S = 'S' then S := 35.0 fi ;
if T = 'T' then T := T_C + 273.15 fi ;
if Dt = 'Dt' then Dt := 0.1 fi ;
if nnodes = 'nnodes' then nnodes := 3 fi ;
if vigrid = 'vigrid' then vigrid := 0 fi ;
if nsolids = 'nsolids' then nsolids := 0 fi ;
if ndissolved = 'ndissolved' then ndissolved := ncompo fi ;
if listsolids = 'listsolids' then listsolids := [] fi ;
if diffdata = 'diffdata' then diffdata := [seq(1.0,i=1..ncompo)] fi ;
if alphadata = 'alphadata' then alphadata := [seq(1.0,i=1..ncompo)] fi ;
if type_up = 'type_up' then type_up := [seq(2,i=1..ncompo)] fi ;
if bnddata_up = 'bnddata_up' then bnddata_up := [seq(1.0,i=1..ncompo)] fi ;
if type_down = 'type_down' then type_down := [seq(1,i=1..ncompo)] fi ;
if bnddata_down = 'bnddata_down' then bnddata_down := [seq(0.0,i=1..ncompo)] fi ;
if nrnotransp = 'nrnotransp' then nrnotransp := 0 fi ;
if listnotransp = 'listnotransp' then listnotransp := [] fi ;
if vic = 'vic' then vic := 2 fi ;
if iniconc = 'iniconc' then iniconc := [seq(1.0,i=1..ncompo)] fi ;
if listinput = 'listinput' then listinput := [seq(i,i=1..ncompo)] fi ;
if file_in_names = 'file_in_names' then file_in_names := [seq(dummy,i=1..ncompo)] fi ;
if noutput = 'noutput' then noutput := 0 fi ;
if nroutput = 'nroutput' then nroutput := 0 fi ;
if listoutput = 'listoutput' then listoutput := [] fi ;
if listroutput = 'listroutput' then listroutput := [] fi ;
if file_names = 'file_names' then file_names := [] fi ;
if file_rnames = 'file_rnames' then file_rnames := [] fi ;
if time_iniout = 'time_iniout' then time_iniout := 0.01 fi ;
if time_intvout = 'time_intvout' then time_intvout := 0.02 fi ;
if nopt_v = 'nopt_v' then nopt_v := 0 fi ;
if ntopt_v = 'ntopt_v' then ntopt_v := 0 fi ;
if nparam_opt = 'nparam_opt' then nparam_opt := 0 fi ;
if maxxmeas_v = 'maxxmeas_v' then maxxmeas_v := 0 fi ;
if maxspmeas_v = 'maxspmeas_v' then maxspmeas_v := 0 fi ;
if opt_name = 'opt_name' then opt_name := [] fi ;
if idpar_v = 'idpar_v' then idpar_v := [] fi ;
if filemeas_name = 'filemeas_name' then filemeas_name := [] fi ;
# Macrofor Extension
# [OPENFAPPEND,UNIT,FILE,STATUS]
`macrofort/openfappend` := proc(level,unit,file,status)
`macrofort/lprint`(cat(`macrofort/space`(level+6),`open(unit=`,unit,`,file='`,file,`',`)):
`macrofort/lprint`(cat(`macrofort/space`(level),` +      status='`,status,`',access='append')`))
end:
p0(ncompo,variables) ;
p1(nreactions,ncompo,variables,variables_old) ;
p2(variables) ;
p3(eqrxnId);
p4(ncompo): # pivoting with rationals - old version is p4old(ncompo)
p5(ncompo,nreactions) ;
p6(ncompo,nreactions) ;
p7(nreactions): # full dR/dC - p7old(nreactions) ;
p8(neqrxns) :
p9(ncompo) :
p10() :
nparphys := 11 :# real*8
phys_name := [al,q0,w0,Db0,por0,area0,t_celsius,salin,delt,depthmax,endt] :
phys_val := [val,vq0,vw0,vDb0,vpor0,varea0,T_C,S,Dt,depth_max,tot_time];
nparphys2 := 7 :# integer
phys_name2 := [iq,iw,iDb,ipor,igrid,iarea,ic] :
phys_val2 := [viq,viw,viDb,vipor,vigrid,viarea,vic];
acg8(nparphys,phys_name,phys_val,nparphys2,phys_name2,phys_val2,dir_f) ;
acg0(nsolids,ndissolved,ncompo,nreactions,bio_name,phys_name,phys_val,phys_name2,phys_val2,dir_f,nnodes):
acg1(type_up,bnddata_up,type_down,bnddata_down,dir_f) ;
acg2(ncompo,diffdata,alphadata,dir_f) ;
acg3(nparam,bio_name,bio_val,dir_f) ;
acg4(ncompo,func,dir_f) ;
acg5(pd,dir_f) ;
#acg6(type_down,nsolids,listsolids,dir_f) ; # openbound
acg7(ncompo,noutput,nroutput,listoutput,listroutput,file_names,file_rnames,time_iniout,time_intvout,dir_f) ;
acg0a(nopt_v,ntopt_v,nparam_opt);
acg0b(maxxmeas_v,maxspmeas_v);
acg9(nparam_opt,opt_name,dir_f) ;
acg10(nparam_opt,opt_name,dir_f) ;
acg11(nopt_v,idpar_v,ntopt_v,filemeas_name) ;
acg12(vic,iniconc,ncompo,listinput,file_in_names,dir_f);
acg13(ncompo) ;
acg14(nrnotransp,listnotransp,dir_f) ;
acg15(nreactions);
acg16(nsolids,listsolids,dir_f);
if(iswitch=1) then acg17(variables,dir_f) fi; # switches.f
try
  if (nswitches < 0) then nswitches := 0 fi;
catch:
  nswitches := 0;
end try;
acg17a(nswitches,dir_f);
try
  if (nparameters < 0) then nparameters := 0 fi;
catch:
  nparameters := 0;
end try;
acg17b(nparameters,dir_f);
try
  if (varporosity < 0) then varporosity := 0 fi;
catch:
  varporosity := 0;
end try;
acg17c(varporosity,dir_f);
#acg18(rjac,dir_f):
"now recompile fortran code. check additional option settings in drivervalues.f";

