X=['input.lammps'];

fid = fopen(X,'wt');

fprintf(fid, '# lammps input file\n');
fprintf(fid,'\n');
X=['# lammps input file\n'];fprintf(fid,X);
fprintf(fid,'\n');
X=['variable	k equal 5 # kcal/mol/A^2\n'];fprintf(fid,X);
X=['variable	zdes equal ',num2str(z0),'\n'];fprintf(fid,X);
fprintf(fid,'\n');
X=['processors	* * 1\n'];fprintf(fid,X);
fprintf(fid,'\n');
X=['boundary	p p p\n'];fprintf(fid,X);
X=['units		real\n'];fprintf(fid,X);
fprintf(fid,'\n');
X=['atom_style	full\n'];fprintf(fid,X);
X=['bond_style	harmonic\n'];fprintf(fid,X);
X=['angle_style	harmonic\n'];fprintf(fid,X);
X=['dihedral_style	opls\n'];fprintf(fid,X);
X=['improper_style	harmonic\n'];fprintf(fid,X);
fprintf(fid,'\n');
X=['pair_style	lj/cut/coul/long 12.0\n'];fprintf(fid,X);
X=['kspace_style	pppm 1.0e-4\n'];fprintf(fid,X);
X=['pair_modify 	shift yes mix arithmetic\n'];fprintf(fid,X);
X=['special_bonds	lj/coul 0.0 0.0 0.5 angle yes dihedral yes\n'];fprintf(fid,X);
fprintf(fid,'\n');
X=['read_data	./data.lammps\n'];fprintf(fid,X);
X=['include 	./PARM.lammps\n'];fprintf(fid,X);
fprintf(fid,'\n');
X=['group		H2O0 type 1\n'];fprintf(fid,X);
X=['group		H2O0 include molecule\n'];fprintf(fid,X);
X=['variable	tzlo equal bound(H2O0,zmin)+3\n'];fprintf(fid,X);
X=['region		btm block INF INF INF INF INF ${tzlo}\n'];fprintf(fid,X);
X=['group		btm region btm\n'];fprintf(fid,X);
X=['group		H2O subtract H2O0 btm\n'];fprintf(fid,X);
X=['group		GO type 3\n'];fprintf(fid,X);
X=['group		GO include molecule\n'];fprintf(fid,X);
X=['group		Car type 3 5 13\n'];fprintf(fid,X);
X=['group		OHH subtract GO Car\n'];fprintf(fid,X);
X=['region		cld cylinder x 0.1 ',num2str(z0),' 0.8 INF INF\n'];fprintf(fid,X);
X=['group		cld region cld\n'];fprintf(fid,X);
X=['group		thr intersect Car cld\n'];fprintf(fid,X);
fprintf(fid,'\n');
X=['compute	twat0 H2O0 temp\n'];fprintf(fid,X);
X=['compute	tOHH OHH temp\n'];fprintf(fid,X);
fprintf(fid,'\n');
X=['fix		mynve1 H2O0 nve\n'];fprintf(fid,X);
X=['fix		mynve2 OHH nve\n'];fprintf(fid,X);
X=['fix		myber1 H2O0 temp/berendsen 10 10 1\n'];fprintf(fid,X);
X=['fix_modify	myber1 temp twat0\n'];fprintf(fid,X);
X=['fix		myber2 OHH temp/berendsen 10 10 1\n'];fprintf(fid,X);
X=['fix_modify	myber2 temp tOHH\n'];fprintf(fid,X);
X=['fix		s1 H2O0 shake 1.0e-4 200 0 b 1 a 1\n'];fprintf(fid,X);
fprintf(fid,'\n');
X=['thermo 	10\n'];fprintf(fid,X);
X=['timestep 	0.005\n'];fprintf(fid,X);
X=['run		500\n'];fprintf(fid,X);
fprintf(fid,'\n');
X=['fix		myber1 H2O0 temp/berendsen 10 100 10\n'];fprintf(fid,X);
X=['fix_modify	myber1 temp twat0\n'];fprintf(fid,X);
X=['fix		myber2 OHH temp/berendsen 10 100 10\n'];fprintf(fid,X);
X=['fix_modify	myber2 temp tOHH\n'];fprintf(fid,X);
fprintf(fid,'\n');
X=['timestep 	0.1\n'];fprintf(fid,X);
X=['run		4500\n'];fprintf(fid,X);
fprintf(fid,'\n');
X=['unfix		myber1\n'];fprintf(fid,X);
X=['unfix		mynve1\n'];fprintf(fid,X);
fprintf(fid,'\n');
X=['fix		mynve1 H2O nve\n'];fprintf(fid,X);
X=['compute	twat H2O temp\n'];fprintf(fid,X);
X=['fix		myber1 H2O temp/berendsen 300 300 100\n'];fprintf(fid,X);
X=['fix_modify	myber1 temp twat\n'];fprintf(fid,X);
X=['fix		myber2 OHH temp/berendsen 300 300 100\n'];fprintf(fid,X);
X=['fix_modify	myber2 temp tOHH\n'];fprintf(fid,X);
fprintf(fid,'\n');
X=['timestep	1.0\n'];fprintf(fid,X);
X=['run		20000\n'];fprintf(fid,X);
X=['reset_timestep	0\n'];fprintf(fid,X);
fprintf(fid,'\n');
X=['dump 		dp1 all atom 100000 dump.lammpstrj\n'];fprintf(fid,X);
X=['unfix		mynve2\n'];fprintf(fid,X);
X=['unfix		myber2\n'];fprintf(fid,X);
X=['fix		mynve2 GO nve\n'];fprintf(fid,X);
X=['fix		myber1 H2O temp/berendsen 300 300 100\n'];fprintf(fid,X);
X=['fix_modify	myber1 temp twat\n'];fprintf(fid,X);
X=['compute	tGO GO temp\n'];fprintf(fid,X);
X=['fix		myber2 GO temp/berendsen 300 300 100\n'];fprintf(fid,X);
X=['fix_modify	myber2 temp tGO\n'];fprintf(fid,X);
fprintf(fid,'\n');
X=['variable	zave equal xcm(thr,z)\n'];fprintf(fid,X);
X=['fix		mytth thr spring tether ${k} NULL NULL ${zdes} 0\n'];fprintf(fid,X);
X=['fix		myat1 all ave/time 10 100 1000 v_zave v_zdes f_mytth[3] f_mytth file Result.dat\n'];fprintf(fid,X);
X=['compute cc1 H2O chunk/atom bin/1d z 0.0 1.0\n'];fprintf(fid,X);
X=['fix myac H2O ave/chunk 100 10000 1000000 cc1 density/number file waterprofile.dat\n'];fprintf(fid,X);
fprintf(fid,'\n');
X=['run 		1000000\n'];fprintf(fid,X);
X=['write_data	data.equilibrium\n'];fprintf(fid,X);
fprintf(fid,'\n');
fclose(fid);