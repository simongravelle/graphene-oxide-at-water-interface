% this matlab script generates 29 folders with the GO nanoparticle in between -28A to 28A from the interface

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the LAMMPS parameter file %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the water molecule properties
MASSwater=load('./WaterMolecule/file.mass'); 
PAIRCOEFFwater=load('./WaterMolecule/file.paircoeff');
BONDwater=load('./WaterMolecule/file.bond');
ANGLEwater=load('./WaterMolecule/file.angle');

% Load the GO nanoparticle properties
Positions=load(['GOnanoparticle/Positions.dat']);
Bonds=load(['GOnanoparticle/Bonds.dat']);
Angles=load(['GOnanoparticle/Angles.dat']);
Dihedrals=load(['GOnanoparticle/Dihedrals.dat']);
Impropers=load(['GOnanoparticle/Impropers.dat']);
MASSgraphitics=load(['GOnanoparticle/ffmass.dat']);
MASSgraphitics(:,1)=MASSgraphitics(:,1)+length(MASSwater(:,1));
PAIRCOEFFgraphitics=load(['GOnanoparticle/ffpair.dat']);
PAIRCOEFFgraphitics(:,4)=PAIRCOEFFgraphitics(:,3); 
PAIRCOEFFgraphitics(:,3)=PAIRCOEFFgraphitics(:,2); 
PAIRCOEFFgraphitics(:,2)=PAIRCOEFFgraphitics(:,1);
PAIRCOEFFgraphitics(:,[1:2])=PAIRCOEFFgraphitics(:,[1:2])+length(MASSwater(:,1));
BONDgraphitics=load(['GOnanoparticle/ffbonds.dat']);
BONDgraphitics(:,1)=BONDgraphitics(:,1)+length(BONDwater(:,1));
ANGLEgraphitics=load(['GOnanoparticle/ffangles.dat']);
ANGLEgraphitics(:,1)=ANGLEgraphitics(:,1)+length(ANGLEwater(:,1));
DIHEDRALgraphitics=load(['GOnanoparticle/ffdihedrals.dat']);
IMPROPERgraphitics=load(['GOnanoparticle/ffimpropers.dat']);

% Combine the water and GO particle properties
MASS=[MASSwater; MASSgraphitics];
PAIRCOEFF=[PAIRCOEFFwater; PAIRCOEFFgraphitics];
BOND=[BONDwater; BONDgraphitics];
ANGLE=[ANGLEwater; ANGLEgraphitics];
DIHEDRAL=[DIHEDRALgraphitics];
IMPROPER=[IMPROPERgraphitics];
Natomtypes=length(MASS(:,1)); 
Nbondtypes=length(BOND(:,1)); 
Nangletypes=length(ANGLE(:,1)); 
Ndihedraltypes=length(DIHEDRAL(:,1));  
Nimpropertypes=length(IMPROPER(:,1)); 

%%%%%%%%%%%%%%%%%%%
% Write PARM file %
%%%%%%%%%%%%%%%%%%%

disp('Start writing the data file')
fid = fopen('PARM.lammps','wt');
fprintf(fid, '#Mass\n\n');
for ii=1:length(MASS(:,1))
		fprintf(fid, 'mass ');
		fprintf(fid, num2str(MASS(ii,1)));
		fprintf(fid, ' ');
		fprintf(fid, num2str(MASS(ii,2)));

	fprintf(fid, '\n');
end
fprintf(fid, '\n');
fprintf(fid, '#Pair Coeff\n\n');
for ii=1:length(PAIRCOEFF(:,1))
		fprintf(fid, 'pair_coeff ');
		fprintf(fid, num2str(PAIRCOEFF(ii,1)));
		fprintf(fid, ' ');
		fprintf(fid, num2str(PAIRCOEFF(ii,2)));
		%fprintf(fid, ' lj/cut/coul/long ');
		fprintf(fid, ' ');
		fprintf(fid, num2str(PAIRCOEFF(ii,3)));
		fprintf(fid, ' ');
		fprintf(fid, num2str(PAIRCOEFF(ii,4)));

	fprintf(fid, '\n');
end

fprintf(fid, '\n');
fprintf(fid, '#Bond\n\n');
for ii=1:length(BOND(:,1))
		fprintf(fid, 'bond_coeff ');
		fprintf(fid, num2str(BOND(ii,1)));
		fprintf(fid, ' ');
		fprintf(fid, num2str(BOND(ii,2)));
		fprintf(fid, ' ');
		fprintf(fid, num2str(BOND(ii,3)));

	fprintf(fid, '\n');
end
fprintf(fid, '\n');
fprintf(fid, '#Angle\n\n');
for ii=1:length(ANGLE(:,1))
		fprintf(fid, 'angle_coeff ');
		fprintf(fid, num2str(ANGLE(ii,1)));
		fprintf(fid, ' ');
		fprintf(fid, num2str(ANGLE(ii,2)));
		fprintf(fid, ' ');
		fprintf(fid, num2str(ANGLE(ii,3)));

	fprintf(fid, '\n');
end
fprintf(fid, '\n');
fprintf(fid, '#Dihedral\n\n');
for ii=1:length(DIHEDRAL(:,1))
		fprintf(fid, 'dihedral_coeff ');
		fprintf(fid, num2str(DIHEDRAL(ii,1)));
		%fprintf(fid, ' opls ');
		fprintf(fid, ' ');
		fprintf(fid, num2str(DIHEDRAL(ii,2)));
		fprintf(fid, ' ');
		fprintf(fid, num2str(DIHEDRAL(ii,3)));
		fprintf(fid, ' ');
		fprintf(fid, num2str(DIHEDRAL(ii,4)));
		fprintf(fid, ' ');
		fprintf(fid, num2str(DIHEDRAL(ii,5)));
	fprintf(fid, '\n');
end
fprintf(fid, '\n');
fprintf(fid, '#Improper\n\n');
for ii=1:length(IMPROPER(:,1))
		fprintf(fid, 'improper_coeff ');
		fprintf(fid, num2str(IMPROPER(ii,1)));
		fprintf(fid, ' ');
		fprintf(fid, num2str(IMPROPER(ii,2)));
		fprintf(fid, ' ');
		fprintf(fid, num2str(IMPROPER(ii,3)));

	fprintf(fid, '\n');
end
disp('Done writing the data file')
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate all 29 folders and input files %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% box dimensions
txlo=-3.4*floor(1.4*20/3.4); txhi=-txlo;
tylo=-3.4*floor(1.4*20/3.4); tyhi=-tylo;
tzlo=-50; tzhi=-tzlo;

% GO properties
PGO=Positions;
PGO(:,3)=PGO(:,3)+2; % shift		
PnmpC=PGO(PGO(:,3)==3 | PGO(:,3)==5 | PGO(:,3)==13,:);
PGO(:,5)=PGO(:,5)-mean(PnmpC(:,5));
PGO(:,6)=PGO(:,6)-mean(PnmpC(:,6));
PGO(:,7)=PGO(:,7)-mean(PnmpC(:,7));
BGO=Bonds;
BGO(:,2)=BGO(:,2)+1; % shift OH water
AGO=Angles;
AGO(:,2)=AGO(:,2)+1; % shift HOH water
DGO=Dihedrals;
IGO=Impropers;

nbH2O=4500; % number of water molecule

pos=[-28:2:28];
for iij=1:length(pos)

	z0=pos(iij);

	cptatom=0;
	cptbond=0;
	cptangle=0;
	cptdihedrals=0;
	cptimpropers=0;
	cptmol=1;

	clear B Ag D A Ig;

	x=0; y=0; z=z0;
	for ii=1:length(BGO(:,1))
		cptbond=cptbond+1;
		B(cptbond,:)=[cptbond BGO(ii,2) BGO(ii,3)+cptatom BGO(ii,4)+cptatom];
	end
	for ii=1:length(AGO(:,1))
		cptangle=cptangle+1;
		Ag(cptangle,:)=[cptangle AGO(ii,2) AGO(ii,3)+cptatom AGO(ii,4)+cptatom AGO(ii,5)+cptatom];
	end
	for ii=1:length(DGO(:,1))
		cptdihedrals=cptdihedrals+1;
		D(cptdihedrals,:)=[cptdihedrals DGO(ii,2) DGO(ii,3)+cptatom DGO(ii,4)+cptatom DGO(ii,5)+cptatom DGO(ii,6)+cptatom];
	end
	for ii=1:length(IGO(:,1))
		cptimpropers=cptimpropers+1;
		Ig(cptimpropers,:)=[cptimpropers IGO(ii,2) IGO(ii,3)+cptatom IGO(ii,4)+cptatom IGO(ii,5)+cptatom IGO(ii,6)+cptatom];
	end
	cpt0=cptatom;
	for ii=1:length(PGO(:,1))
		cptatom=cptatom+1;
		A(cptatom,:)=[PGO(ii,1)+cpt0 cptmol PGO(ii,3) PGO(ii,4) PGO(ii,7)+x PGO(ii,6)+y PGO(ii,5)+z];
	end

	%%%%%%%%%%%%%%%%%%%
	% water molecules %
	%%%%%%%%%%%%%%%%%%%

	lengthAbeforewater=length(A(:,1)); % number of atmo before starting to insert water molecule (for overlap detection)

	dx=3.4; dy=3.4; dz=3.4;

	PH2O=load('./WaterMolecule/Position.dat');
	BH2O=load('./WaterMolecule/Bond.dat');
	AH2O=load('./WaterMolecule/Angle.dat');

	x=txlo+dx/4;
	y=tylo+dy/4;
	z=-46;
	cptwat=0;
	while z < tzhi
		while y < tyhi-dy/2
			while x < txhi-dx/2
				
			
				xeff=x+(rand-0.5)/10;
				yeff=y+(rand-0.5)/10;
				zeff=z+(rand-0.5)/10;
			
				for ii=1:length(PH2O(:,1))
					watmol(ii,:)=[cptatom cptmol PH2O(ii,3) PH2O(ii,4) PH2O(ii,5)+xeff PH2O(ii,6)+yeff PH2O(ii,7)+zeff];
				end
				mind=1e5;
				for ii=1:length(watmol(:,1))
					for jj=1:lengthAbeforewater
						d=sqrt((watmol(ii,5)-A(jj,5))^2+(watmol(ii,6)-A(jj,6))^2+(watmol(ii,7)-A(jj,7))^2);
						if d<mind
							mind=d;
						end
					end
				end
				if mind>3.6 && cptwat<nbH2O
					cptmol=cptmol+1; cptwat=cptwat+1;
					for ii=1:length(BH2O(:,1))
						cptbond=cptbond+1;
						B(cptbond,:)=[cptbond BH2O(ii,2) BH2O(ii,3)+cptatom BH2O(ii,4)+cptatom];
					end
					for ii=1:length(AH2O(:,1))
						cptangle=cptangle+1;
						Ag(cptangle,:)=[cptangle AH2O(ii,2) AH2O(ii,3)+cptatom AH2O(ii,4)+cptatom AH2O(ii,5)+cptatom];
					end
					for ii=1:length(PH2O(:,1))
						cptatom=cptatom+1;
						A(cptatom,:)=[cptatom cptmol PH2O(ii,3) PH2O(ii,4) PH2O(ii,5)+xeff PH2O(ii,6)+yeff PH2O(ii,7)+zeff];
					end
				end

				x=x+dx;
			end
			x=txlo+dx/4;
			y=y+dy;
		end
		y=tylo+dy/4;
		z=z+dz;
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%
	% write data.lammps file %
	%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	fid = fopen('data.lammps','wt');
	fprintf(fid, '# System\n\n');
	fprintf(fid, num2str(cptatom));
	fprintf(fid, ' atoms\n');
	fprintf(fid, num2str(cptbond));
	fprintf(fid, ' bonds\n');
	fprintf(fid, num2str(cptangle));
	fprintf(fid, ' angles\n');
	fprintf(fid, num2str(cptdihedrals));
	fprintf(fid, ' dihedrals\n');
	fprintf(fid, num2str(cptimpropers));
	fprintf(fid, ' impropers\n\n');
	fprintf(fid, num2str(Natomtypes));
	fprintf(fid, ' atom types\n');
	fprintf(fid, num2str(Nbondtypes));
	fprintf(fid, ' bond types\n');
	fprintf(fid, num2str(Nangletypes));
	fprintf(fid, ' angle types\n');
	fprintf(fid, num2str(Ndihedraltypes));
	fprintf(fid, ' dihedral types\n');
	fprintf(fid, num2str(Nimpropertypes));
	fprintf(fid, ' improper types\n');
	fprintf(fid, 'extra bond per atom ');
	fprintf(fid, num2str(2));
	fprintf(fid, '\n');
	fprintf(fid, 'extra angle per atom ');
	fprintf(fid, num2str(1));
	fprintf(fid, '\n');
	fprintf(fid, 'extra special per atom ');
	fprintf(fid, num2str(2));
	fprintf(fid, '\n');
	fprintf(fid, num2str([txlo txhi]));
	fprintf(fid, ' xlo xhi\n');
	fprintf(fid, num2str([tylo tyhi]));
	fprintf(fid, ' ylo yhi\n');
	fprintf(fid, num2str([tzlo tzhi]));
	fprintf(fid, ' zlo zhi\n\n');
	fprintf(fid, 'Atoms\n\n');
	for ii=1:length(A(:,1))
		for jj=1:7
			fprintf(fid, num2str(A(ii,jj)));
			fprintf(fid, '	');
		end
		fprintf(fid, '\n');
	end
	fprintf(fid, '\n');
	fprintf(fid, 'Bonds\n\n');
	for ii=1:length(B(:,1))
		for jj=1:4
			fprintf(fid, num2str(B(ii,jj)));
			fprintf(fid, '	');
		end
		fprintf(fid, '\n');
	end
	fprintf(fid, '\n');
	fprintf(fid, 'Angles\n\n');
	for ii=1:length(Ag(:,1))
		for jj=1:5
			fprintf(fid, num2str(Ag(ii,jj)));
			fprintf(fid, '	');
		end
		fprintf(fid, '\n');
	end
	fprintf(fid, '\n');
	fprintf(fid, 'Dihedrals\n\n');
	for ii=1:length(D(:,1))
		for jj=1:6
			fprintf(fid, num2str(D(ii,jj)));
			fprintf(fid, '	');
		end
		fprintf(fid, '\n');
	end
	fprintf(fid, '\n');
	fprintf(fid, 'Impropers\n\n');
	for ii=1:length(Ig(:,1))
		for jj=1:6
			fprintf(fid, num2str(Ig(ii,jj)));
			fprintf(fid, ' ');
		end
		fprintf(fid, '\n');
	end
	fclose(fid);

	if z0<0
		a=-z0;
		X=['../Datam',num2str(abs(a))];
		if ~exist(X, 'dir')
			mkdir(X);
		end	
		file=sprintf('../Datam%d',a);
		
		GenerateLAMMPSInput
		
		copyfile('./data.lammps',file);
		copyfile('./PARM.lammps',file);
		copyfile('./input.*',file);
		copyfile('./SurfEnergy.sh',file);
	else
		a=z0;
		X=['../Data',num2str(abs(a))];
		if ~exist(X, 'dir')
			mkdir(X);
		end	
		file=sprintf('../Data%d',a);
		
		GenerateLAMMPSInput
		
		copyfile('./data.lammps',file);
		copyfile('./PARM.lammps',file);
		copyfile('./input.*',file);
		copyfile('./SurfEnergy.sh',file);
	end
	
end

























